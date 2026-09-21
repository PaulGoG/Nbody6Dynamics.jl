# =============================================================================
# GPU validation driver: host record, GPU-gated suite, GPU and CPU pipelines,
# scaling benchmark, every stage logged into one directory under runs/
# =============================================================================

const _VALIDATION_STAGES = (:suite, :gpu, :cpu, :bench)

"""
    run_gpu_validation(; base_dir = _PROJECT_ROOT, stages = [:suite, :gpu, :cpu, :bench],
                       bench_n = [20000, 50000, 100000], bench_threads = [4, 8],
                       bench_gpu_lists = nothing, bench_tcrit = 0.25,
                       dry_run = false) -> String

Run the acceptance sequence of a CUDA host and collect everything under
`<base_dir>/runs/gpu_validation_<machine>_<timestamp>/`, `<machine>` being
the hostname with compact CPU and GPU tags so that
machines sharing a hostname stay distinguishable: `HOST_INFO.toml`
(hardware fingerprint with the GPU query, compute capabilities, CUDA path
and `nvcc` release, `gcc`, `gfortran` and glibc versions, the host
compilers `nvcc` can be offered and the verdict of the host-compiler probe,
Julia, package commit),
one `<stage>.log` per stage, the benchmark CSV, the `RUN_INFO.toml` and
`telemetry.csv` of every benchmark run, and `VALIDATION.toml` (command,
status, exit code and duration per stage, the run directories the pipelines
created). Stages, in order:

- `:suite` — the test suite with `NBODY6_GPU_TESTS=1` (CUDA build in a
  temporary tree, N = 1000 on one device, then on two when present);
- `:gpu` — `input_files/gpu/gpu_pipeline.toml`: CUDA build into
  `backend/Nbody6PPGPU-beijing-gpu` and the 2 × 25 000-star merger on device 0;
- `:cpu` — `input_files/gpu/cpu_pipeline.toml`: the AVX reference build and
  the same merger;
- `:bench` — `bench/gpu_scaling.jl` over `bench_n` × `bench_threads` ×
  `bench_gpu_lists` (default `[[0]]`, plus `[0, 1]` when two devices are
  visible) for `bench_tcrit` N-body time units; needs both binaries.

Each stage is a separate Julia process whose output is echoed and written
to its log with terminal escape sequences removed; a failing stage is
recorded and the later stages still run. `VALIDATION.toml` is rewritten
after every stage, so an interrupted sequence leaves a readable record.

A stage's verdict rests on its artefacts as well as its exit code: a `:gpu`
or `:cpu` stage that exits zero without leaving a run directory marked
`[pipeline] completed` is recorded `incomplete` with the reason,
because the run summary is written when the
engine exits and a process killed during post-processing or plotting would
otherwise pass.
With `dry_run = true` the host record and the planned commands are written
and nothing is executed. Returns the validation directory.

# Example
```julia
dir = run_gpu_validation(; stages = [:suite, :gpu])
TOML.parsefile(joinpath(dir, "VALIDATION.toml"))["results"]["suite"]["status"]
```
"""
function run_gpu_validation(;
    base_dir::AbstractString = _PROJECT_ROOT,
    stages::AbstractVector{Symbol} = collect(_VALIDATION_STAGES),
    bench_n::AbstractVector{<:Integer} = [20000, 50000, 100000],
    bench_threads::AbstractVector{<:Integer} = [4, 8],
    bench_gpu_lists::Union{Nothing,AbstractVector{<:AbstractVector{<:Integer}}} = nothing,
    bench_tcrit::Real = 0.25,
    dry_run::Bool = false,
)::String
    base_dir = abspath(base_dir)
    for s in stages
        s in _VALIDATION_STAGES || throw(
            ArgumentError(
                "unknown validation stage :$s; choose from " *
                join(string.(":", _VALIDATION_STAGES), ", "),
            ),
        )
    end
    isempty(stages) && throw(ArgumentError("stages must not be empty"))
    isempty(bench_n) && throw(ArgumentError("bench_n must not be empty"))
    isempty(bench_threads) && throw(ArgumentError("bench_threads must not be empty"))
    bench_tcrit > 0 || throw(ArgumentError("bench_tcrit must be > 0; got $bench_tcrit"))

    t_start = time()
    stamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    # The host record is taken first: its machine identity names the directory,
    # because the hostname alone does not distinguish the machines of a fleet.
    host = _validation_host_record(base_dir)
    out_dir = joinpath(base_dir, "runs", "gpu_validation_$(_safe_name(host["machine"]))_$stamp")
    mkpath(out_dir)

    open(joinpath(out_dir, "HOST_INFO.toml"), "w") do io
        TOML.print(io, host)
    end
    @info "GPU validation directory: $out_dir"
    @info "Machine: $(host["machine"])"
    @info "Host: $(host["host"]); GPU: $(host["gpu"]); nvcc: $(host["nvcc_release"])"
    @info "nvcc host-compiler probe: $(first(eachline(IOBuffer(host["nvcc_probe"]))))" flags =
        host["nvcc_host_flags"]

    n_devices = length(host["compute_capabilities"])
    gpu_lists =
        bench_gpu_lists === nothing ? (n_devices ≥ 2 ? [[0], [0, 1]] : [[0]]) :
        [Int.(l) for l in bench_gpu_lists]
    commands = _validation_commands(base_dir, bench_n, bench_threads, gpu_lists, bench_tcrit)

    summary = Dict{String,Any}(
        "started" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        "stages" => String.(stages),
        "dry_run" => dry_run,
        "results" => Dict{String,Any}(),
    )
    _write_validation_summary(out_dir, summary)

    for s in stages
        cmd = commands[s]
        entry = Dict{String,Any}("command" => _command_string(cmd))
        unmet = _stage_prerequisite(s, base_dir)
        if dry_run
            entry["status"] = "planned"
        elseif unmet !== nothing
            entry["status"] = "skipped"
            entry["reason"] = unmet
            @warn "Stage :$s skipped: $unmet"
        else
            log = joinpath(out_dir, "$s.log")
            @info "Stage :$s → $(basename(log))"
            t0 = time()
            code, attempts = _run_stage_with_retry(cmd, log, s)
            entry["exit_code"] = code
            attempts > 1 && (entry["attempts"] = attempts)
            code ≥ 128 && (entry["signal"] = code - 128)
            entry["seconds"] = round(time() - t0; digits = 1)
            entry["log"] = basename(log)
            # A zero exit code is necessary but not sufficient: a pipeline
            # stage killed after its engine finished leaves a run directory
            # that looks complete. Demand the completion marker as well.
            incomplete = code == 0 ? _stage_incomplete(s, base_dir, t0) : nothing
            entry["status"] =
                code != 0 ? "failed" : incomplete === nothing ? "passed" : "incomplete"
            incomplete === nothing || (entry["reason"] = incomplete)
            if code != 0
                @warn "Stage :$s failed with exit code $code; see $(basename(log))"
            elseif incomplete !== nothing
                @warn "Stage :$s exited 0 but did not finish: $incomplete; see $(basename(log))"
            else
                @info "Stage :$s passed ($(_format_elapsed(time() - t0)))"
            end
        end
        summary["results"][String(s)] = entry
        _write_validation_summary(out_dir, summary)
    end

    if !dry_run
        summary["run_dirs"] = _new_entries(joinpath(base_dir, "runs"), t_start, out_dir)
        summary["bench_results"] =
            _collect_bench_artefacts(joinpath(base_dir, "bench"), out_dir, t_start)
        summary["finished"] = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
        summary["seconds"] = round(time() - t_start; digits = 1)
        _write_validation_summary(out_dir, summary)
    end
    return out_dir
end

"""
    _stage_prerequisite(stage, base_dir) -> Union{Nothing,String}

`nothing` when `stage` can run, otherwise the reason it cannot: the
benchmark needs the CPU and the GPU build trees that the `:cpu` and `:gpu`
stages produce under `base_dir/backend`.
"""
function _stage_prerequisite(stage::Symbol, base_dir::AbstractString)::Union{Nothing,String}
    stage == :bench || return nothing
    missing_trees = String[]
    for tree in ("Nbody6PPGPU-beijing", "Nbody6PPGPU-beijing-gpu")
        isdir(joinpath(base_dir, "backend", tree, "build")) || push!(missing_trees, tree)
    end
    isempty(missing_trees) && return nothing
    return "no build under backend/" *
           join(missing_trees, " and backend/") *
           " (the :cpu and :gpu stages produce them)"
end

"""
    _stage_incomplete(stage, base_dir, t0) -> Union{Nothing,String}

`nothing` when a stage that exited zero also left the artefacts it is
supposed to produce, otherwise the reason it did not. The `:gpu` and `:cpu`
stages must each leave one run directory under `base_dir/runs` started at or
after `t0`, carrying the `[pipeline] completed` marker and an engine segment
that reached END RUN; `:suite` and `:bench` are judged by their exit code
alone, the benchmark because [`_collect_bench_artefacts`](@ref) reports what
it gathered.

An exit code is not enough on its own: the run summary is written when the
engine exits, so a stage killed during integration, post-processing or
plotting leaves a summary and, under libuv, may still report exit code 0.
"""
function _stage_incomplete(stage::Symbol, base_dir::AbstractString, t0::Real)::Union{Nothing,String}
    stage in (:gpu, :cpu) || return nothing
    runs_dir = joinpath(base_dir, "runs")
    isdir(runs_dir) || return "no runs directory under $(basename(base_dir))/runs"
    candidates = filter(readdir(runs_dir; join = true)) do d
        isdir(d) && !startswith(basename(d), "gpu_validation_") && mtime(d) ≥ t0 - 1
    end
    isempty(candidates) && return "the stage left no run directory under runs/"
    finished = filter(_pipeline_completed, candidates)
    if isempty(finished)
        newest = basename(argmax(mtime, candidates))
        return "runs/$newest has no [pipeline] completed marker: the pipeline was " *
               "interrupted after the engine exited (post-processing or plotting)"
    end
    # A pipeline also completes on partial output: the engine must have
    # reached END RUN for the stage to count.
    any(d -> _engine_completed(d) !== false, finished) && return nothing
    newest = basename(argmax(mtime, finished))
    return "runs/$newest completed its pipeline on partial output: the engine " *
           "ended without END RUN (killed, or halted on its energy check)"
end

"""A name with the characters outside `[A-Za-z0-9_-]` replaced by `_`, for a directory name."""
_safe_name(name::AbstractString) = replace(String(name), r"[^A-Za-z0-9_-]" => "_")

"""
    _validation_host_record(base_dir) -> Dict{String,Any}

The hardware fingerprint with the GPU query, plus the compute capabilities,
CUDA path and `nvcc` release, `gcc`, `gfortran` and glibc banner lines, the
banner of every host compiler the build could offer `nvcc` through `-ccbin`
([`_host_compiler_candidates`](@ref)), the outcome of the host-compiler
probe ([`_nvcc_probe_record`](@ref)), and the package commit of `base_dir`.
"""
function _validation_host_record(base_dir::AbstractString)::Dict{String,Any}
    d = _hardware_fingerprint(; gpu_probe = true)
    cuda_path = detect_cuda_path()
    d["date"] = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
    d["package_commit"] = _source_stamp(base_dir)
    d["compute_capabilities"] = detect_compute_capabilities()
    d["cuda_path"] = cuda_path
    d["nvcc_release"] = nvcc_release(cuda_path)
    d["gcc"] = _tool_banner(`gcc --version`)
    d["gfortran"] = _tool_banner(`gfortran --version`)
    d["glibc"] = _tool_banner(`getconf GNU_LIBC_VERSION`)
    d["host_compilers"] = Dict{String,String}(
        cc => _tool_banner(`$cc --version`) for cc in _host_compiler_candidates()
    )
    flags, verdict = _nvcc_probe_record(cuda_path, d["nvcc_release"])
    d["nvcc_host_flags"] = flags
    d["nvcc_probe"] = verdict
    return d
end

"""
    _nvcc_probe_record(cuda_path, release) -> (flags, verdict)

The host-compiler probe ([`_nvcc_host_compiler_flags`](@ref)) run for the
host record: the options it settled on and `"passed"`, or no options and
the probe's message — the `nvcc` output of every attempt — when no host
compiler works. Without an `nvcc` release the probe is skipped.
"""
function _nvcc_probe_record(
    cuda_path::AbstractString,
    release::AbstractString,
)::Tuple{Vector{String},String}
    isempty(release) && return (String[], "skipped: nvcc unavailable")
    flags = try
        _nvcc_host_compiler_flags(cuda_path)
    catch e
        e isa ErrorException || rethrow()
        return (String[], sprint(showerror, e))
    end
    return (flags, "passed")
end

"""First output line of `cmd`, or `"unavailable"` when the tool is absent or fails."""
function _tool_banner(cmd::Cmd)::String
    output = try
        read(pipeline(cmd; stderr = devnull), String)
    catch e
        e isa Union{ProcessFailedException,Base.IOError} || rethrow()
        return "unavailable"
    end
    first_line = something(iterate(eachline(IOBuffer(output))), ("",))[1]
    return isempty(strip(first_line)) ? "unavailable" : String(strip(first_line))
end

"""
    _validation_commands(base_dir, bench_n, bench_threads, gpu_lists, tcrit) -> Dict{Symbol,Cmd}

The four stage commands, each a Julia child process running in `base_dir`
with the environment it needs (`NBODY6_GPU_TESTS` for the suite, the CPU
and GPU backend trees for the benchmark).
"""
function _validation_commands(
    base_dir::AbstractString,
    bench_n::AbstractVector{<:Integer},
    bench_threads::AbstractVector{<:Integer},
    gpu_lists::AbstractVector{<:AbstractVector{<:Integer}},
    tcrit::Real,
)::Dict{Symbol,Cmd}
    julia = Base.julia_cmd()
    setup = joinpath(base_dir, "scripts", "run_setup.jl")
    gpu_cfg = joinpath(base_dir, "input_files", "gpu", "gpu_pipeline.toml")
    cpu_cfg = joinpath(base_dir, "input_files", "gpu", "cpu_pipeline.toml")
    bench = joinpath(base_dir, "bench", "gpu_scaling.jl")
    cpu_tree = joinpath(base_dir, "backend", "Nbody6PPGPU-beijing")
    gpu_tree = joinpath(base_dir, "backend", "Nbody6PPGPU-beijing-gpu")
    gpu_arg = join((join(string.(l), ",") for l in gpu_lists), ";")
    in_base(cmd) = Cmd(cmd; dir = base_dir)
    return Dict{Symbol,Cmd}(
        :suite => in_base(
            addenv(
                `$julia --project=$base_dir -e 'using Pkg; Pkg.test()'`,
                "NBODY6_GPU_TESTS" => "1",
            ),
        ),
        :gpu => in_base(`$julia $setup $gpu_cfg`),
        :cpu => in_base(`$julia $setup $cpu_cfg`),
        :bench => in_base(
            addenv(
                `$julia $bench $(join(string.(bench_n), ",")) $(join(string.(bench_threads), ",")) $gpu_arg $(string(tcrit))`,
                "NBODY6_CPU_BACKEND" => cpu_tree,
                "NBODY6_GPU_BACKEND" => gpu_tree,
            ),
        ),
    )
end

"""The command line of `cmd` as one string, with its added environment variables in front."""
function _command_string(cmd::Cmd)::String
    env = cmd.env === nothing ? String[] : String[e for e in cmd.env if startswith(e, "NBODY6_")]
    words = join(cmd.exec, " ")
    return isempty(env) ? words : join(env, " ") * " " * words
end

const _ANSI_ESCAPE = r"\e\[[0-9;?]*[A-Za-z]|\r"

"""
    _tee_run(cmd, logfile) -> Int

Run `cmd` with stdout and stderr merged, echoing every line to this
process's stdout and writing it to `logfile` with terminal escape
sequences and carriage returns removed. Returns the exit code, and
`128 + signal` for a process killed by a signal — a signal-killed process
reports `exitcode == 0` through libuv, so returning that alone would record
a crash as a success. A command that cannot be spawned throws.
"""
function _tee_run(cmd::Base.AbstractCmd, logfile::AbstractString)::Int
    out = Pipe()
    proc = run(pipeline(cmd; stdout = out, stderr = out); wait = false)
    close(out.in)
    open(logfile, "w") do log
        for line in eachline(out)
            println(stdout, line)
            println(log, replace(line, _ANSI_ESCAPE => ""))
            flush(log)
        end
    end
    wait(proc)
    return proc.termsignal == 0 ? Int(proc.exitcode) : 128 + Int(proc.termsignal)
end

"""
    _run_stage_with_retry(cmd, log, stage) -> (code, attempts)

Run `cmd` through [`_tee_run`](@ref), writing to `log`. A stage killed by a
signal is an upstream crash rather than a verdict on the package — Julia
1.13's optimiser segfaults intermittently during JIT compilation, and it has
cost whole validation runs — so that output is kept as
`<stage>.signal<N>.log` and the command is run once more. Returns the final
exit code and the number of attempts.
"""
function _run_stage_with_retry(
    cmd::Base.AbstractCmd,
    log::AbstractString,
    stage::Symbol,
)::Tuple{Int,Int}
    code = _tee_run(cmd, log)
    code < 128 && return (code, 1)
    signal = code - 128
    crash_log = string(splitext(log)[1], ".signal", signal, ".log")
    mv(log, crash_log; force = true)
    @warn "Stage :$stage was killed by signal $signal; its output is kept in " *
          "$(basename(crash_log)) and the stage is being run once more."
    return (_tee_run(cmd, log), 2)
end

function _write_validation_summary(out_dir::AbstractString, summary::Dict{String,Any})
    open(joinpath(out_dir, "VALIDATION.toml"), "w") do io
        TOML.print(io, summary)
    end
    return nothing
end

"""Names of the directories in `dir` modified after `t_start` (Unix seconds), except `exclude`."""
function _new_entries(dir::AbstractString, t_start::Real, exclude::AbstractString)::Vector{String}
    isdir(dir) || return String[]
    names = String[]
    for name in sort(readdir(dir))
        path = joinpath(dir, name)
        (isdir(path) && path != exclude && mtime(path) ≥ t_start - 1) || continue
        push!(names, name)
    end
    return names
end

"""
    _collect_bench_artefacts(bench_dir, out_dir, t_start) -> Vector{String}

Copy the benchmark CSVs written after `t_start` into `out_dir` and the
`RUN_INFO.toml` and `telemetry.csv` of every benchmark run started after
it into `out_dir/bench_runs/<run_id>/`. Returns the copied CSV names.
"""
function _collect_bench_artefacts(
    bench_dir::AbstractString,
    out_dir::AbstractString,
    t_start::Real,
)::Vector{String}
    copied = String[]
    results = joinpath(bench_dir, "results")
    if isdir(results)
        for name in sort(readdir(results))
            path = joinpath(results, name)
            (endswith(name, ".csv") && mtime(path) ≥ t_start - 1) || continue
            cp(path, joinpath(out_dir, name); force = true)
            push!(copied, name)
        end
    end
    runs = joinpath(bench_dir, "runs")
    for run_id in _new_entries(runs, t_start, out_dir)
        for file in ("RUN_INFO.toml", "telemetry.csv")
            src = joinpath(runs, run_id, file)
            isfile(src) || continue
            dest = joinpath(out_dir, "bench_runs", run_id)
            mkpath(dest)
            cp(src, joinpath(dest, file); force = true)
        end
    end
    return copied
end
