# =============================================================================
# Simulation execution with run ID tracking and real-time monitoring
# =============================================================================

"""
    generate_run_id(prefix::String = "run") -> String

Generate a unique run identifier: `{prefix}_YYYYMMDD_HHMMSS_{4hex}`.
"""
function generate_run_id(prefix::String = "run")::String
    ts = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    hex = bytes2hex(rand(UInt8, 2))
    return "$(prefix)_$(ts)_$(hex)"
end

# ---------------------------------------------------------------------------
# Elapsed-time formatting
# ---------------------------------------------------------------------------

"""
    _format_elapsed(seconds::Float64) -> String

Human-readable elapsed time: `"1.2 s"`, `"3m 42s"`, `"1h 05m 12s"`.
"""
function _format_elapsed(seconds::Float64)::String
    s = round(Int, seconds)
    if s < 60
        return @sprintf("%.1f s", seconds)
    elseif s < 3600
        m, sec = divrem(s, 60)
        return @sprintf("%dm %02ds", m, sec)
    else
        h, rem = divrem(s, 3600)
        m, sec = divrem(rem, 60)
        return @sprintf("%dh %02dm %02ds", h, m, sec)
    end
end

# ---------------------------------------------------------------------------
# Spinner characters for heartbeat
# ---------------------------------------------------------------------------

const _SPINNER = ['⠋', '⠙', '⠹', '⠸', '⠼', '⠴', '⠦', '⠧', '⠇', '⠏']

"""
    run_simulation(cfg::Nbody6Config; base_dir = _PROJECT_ROOT) -> String

Execute the Nbody6++ simulation with:
- Unique run ID and isolated output directory
- Bash wrapper for `ulimit -s unlimited` and CUDA env
- Real-time ADJUST line monitoring with wall-clock timer and heartbeat spinner
- Frozen config snapshot saved alongside output

Returns the path to the run directory.
"""
function run_simulation(cfg::Nbody6Config; base_dir::AbstractString = _PROJECT_ROOT)::String
    sim = cfg.simulation

    # --- Resolve input file (relative to the backend source tree) ---
    src_dir = joinpath(base_dir, cfg.install.install_dir)
    input_path = abspath(joinpath(src_dir, sim.input_file))
    isfile(input_path) || error("Input file not found: $input_path")

    # --- Generate run ID and create directory structure ---
    #   runs/<run_id>/
    #       output/      — simulation output files
    #       plots/       — post-processing plots (created later)
    #       config.toml  — frozen config
    run_id = generate_run_id(sim.run_id_prefix)
    run_dir = abspath(joinpath(base_dir, sim.runs_dir, run_id))
    out_dir = joinpath(run_dir, "output")
    mkpath(out_dir)

    @info "Run ID:  $run_id"
    @info "Run dir: $run_dir"

    return _execute_simulation(cfg, run_dir, out_dir, input_path; base_dir = base_dir)
end

"""
    _execute_simulation(cfg, run_dir, out_dir, input_path;
                        base_dir = _PROJECT_ROOT, label = "simulation") -> String

Shared execution core for [`run_simulation`](@ref) and the merger pipeline:
locates the binary, freezes the config into `run_dir`, copies the binary
into `out_dir` for reproducibility, writes the launch script, and runs it
under the teed run log (§9) with the opt-in live monitor. `input_path` must
be absolute (the launch script executes from `out_dir`). Returns `run_dir`.
"""
function _execute_simulation(
    cfg::Nbody6Config,
    run_dir::AbstractString,
    out_dir::AbstractString,
    input_path::AbstractString;
    base_dir::AbstractString = _PROJECT_ROOT,
    label::AbstractString = "simulation",
)::String
    sim = cfg.simulation
    run_id = basename(run_dir)

    # --- Locate binary ---
    src_dir = joinpath(base_dir, cfg.install.install_dir)
    binary = _find_binary(src_dir, sim.binary_name)

    # --- Save frozen config ---
    save_config(cfg, joinpath(run_dir, "config.toml"))

    # --- Copy binary for reproducibility ---
    local_binary = joinpath(out_dir, basename(binary))
    cp(binary, local_binary; force = true)
    chmod(local_binary, 0o755)

    # --- Build launch script ---
    stdout_path = joinpath(out_dir, cfg.postprocess.stdout_file)
    stderr_path = joinpath(out_dir, "err1000")
    launch_script =
        _write_launch_script(out_dir, local_binary, input_path, stdout_path, stderr_path, cfg)

    # --- Execute, teeing pipeline logs to the run directory (§9) ---
    _with_run_log(run_dir) do
        omp_threads = _effective_omp_threads(sim)
        threads_total = omp_threads * sim.mpi_ranks
        if threads_total > Sys.CPU_THREADS
            @warn "CPU oversubscription: $(sim.mpi_ranks) rank(s) × $omp_threads OpenMP threads " *
                  "exceed the $(Sys.CPU_THREADS) logical CPUs of this host"
        end
        @info "Backend threads: $omp_threads OpenMP × $(sim.mpi_ranks) MPI rank(s)"

        cpu_before = _children_cpu_times()
        t_start = time()
        @info "Starting $label..."
        process = cd(out_dir) do
            run(`bash $launch_script`; wait = false)
        end

        # The launch script execs the binary, so its PID is the process PID;
        # an mpirun launcher is handled through the process-tree scan.
        monitor = if sim.telemetry_interval > 0
            _start_telemetry(
                run_dir,
                getpid(process),
                t_start;
                interval = sim.telemetry_interval,
                gpu_probe = cfg.build.enable_gpu,
            )
        else
            nothing
        end

        # Live ticker is opt-in and interactive-only; the Fortran stdout is
        # captured to out1000 regardless.
        if sim.monitor && stderr isa Base.TTY
            _monitor_stdout_file(stdout_path, process, t_start)
        end
        wait(process)

        elapsed = time() - t_start
        telemetry =
            _finish_telemetry(monitor, cpu_before, _children_cpu_times(), elapsed, threads_total)
        merge!(telemetry, _backend_performance(stdout_path, stderr_path))

        if !success(process)
            @warn "Simulation exited with non-zero status ($(process.exitcode)) after $(_format_elapsed(elapsed))"
        end

        # --- Write run summary ---
        _write_run_summary(
            cfg,
            run_dir,
            run_id,
            stdout_path,
            out_dir,
            elapsed;
            telemetry = telemetry,
        )
        if haskey(telemetry, "cpu_efficiency")
            @info @sprintf(
                "CPU: %.1f s user + %.1f s system on %d thread(s); efficiency %.2f",
                telemetry["cpu_user_s"],
                telemetry["cpu_system_s"],
                threads_total,
                telemetry["cpu_efficiency"]
            )
        end

        @info "Simulation complete. Run: $run_id  ($(_format_elapsed(elapsed)))"
    end
    return run_dir
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""
Write a self-contained bash launch script that sets `ulimit`, environment
variables, and runs the simulation with proper I/O redirection.

A runtime-generated shell script is the one sanctioned bash use in this
package: `ulimit -s unlimited` (mandatory for the stack-heavy Fortran) can
only be applied to the child process from a wrapping shell — Julia's `run`
cannot set resource limits on the spawned binary.
"""
function _write_launch_script(
    run_dir::AbstractString,
    binary::AbstractString,
    input::AbstractString,
    stdout_file::AbstractString,
    stderr_file::AbstractString,
    cfg::Nbody6Config,
)::String
    path = joinpath(run_dir, "_launch.sh")
    sim = cfg.simulation
    build = cfg.build

    open(path, "w") do io
        println(io, "#!/bin/bash")
        println(io, "# Generated by Nbody6Dynamics.jl — DO NOT EDIT")
        println(io, "set -o pipefail")
        println(io, "ulimit -s unlimited")
        println(io, "export OMP_STACKSIZE=4096M")
        # OpenMP thread cap; 0 leaves the runtime default (inherited
        # OMP_NUM_THREADS or every logical CPU).
        sim.omp_threads > 0 && println(io, "export OMP_NUM_THREADS=$(sim.omp_threads)")

        # CUDA environment
        if build.enable_gpu
            cuda = isempty(build.cuda_path) ? detect_cuda_path() : build.cuda_path
            if !isempty(cuda)
                println(io, "export CUDA_HOME=\"$cuda\"")
                println(io, "export PATH=\"$cuda/bin:\$PATH\"")
                println(io, "export LD_LIBRARY_PATH=\"$cuda/lib64:\$LD_LIBRARY_PATH\"")
            end
        end

        # The simulation command (use stdbuf for line-buffered output if available)
        stdbuf_prefix = "command -v stdbuf >/dev/null 2>&1 && STDBUF='stdbuf -oL' || STDBUF=''"
        println(io, stdbuf_prefix)

        if build.enable_mpi && sim.mpi_ranks > 1
            println(
                io,
                "exec mpirun --bind-to none -np $(sim.mpi_ranks) " *
                "\$STDBUF \"$binary\" < \"$input\" > \"$stdout_file\" 2> \"$stderr_file\"",
            )
        else
            println(
                io,
                "exec \$STDBUF \"$binary\" < \"$input\" > \"$stdout_file\" 2> \"$stderr_file\"",
            )
        end
    end
    chmod(path, 0o755)
    return path
end

"""
Poll the stdout file while the simulation runs and print ADJUST summaries
with wall-clock elapsed time and a heartbeat spinner.

The spinner updates every 2 seconds on the same line, showing the user that
the process is alive even during quiet periods. When a diagnostics line
(ADJUST or TIME[NB]) is detected, it prints as a full log line and resets
the spinner.
"""
function _monitor_stdout_file(path::String, process::Base.Process, t_start::Float64 = time())
    last_pos = 0
    spin_idx = 1
    last_event_t = t_start   # wall-clock of last printed diagnostics line

    while process_running(process)
        new_output = false

        if isfile(path)
            sz = filesize(path)
            if sz > last_pos
                open(path, "r") do io
                    seek(io, last_pos)
                    while !eof(io)
                        line = readline(io; keep = false)
                        if _print_monitor_line(line, t_start)
                            new_output = true
                            last_event_t = time()
                        end
                    end
                    last_pos = position(io)
                end
            end
        end

        # Heartbeat spinner (overwritten in-place) when no new diagnostics
        if !new_output
            elapsed = time() - t_start
            idle = time() - last_event_t
            spin_ch = _SPINNER[mod1(spin_idx, length(_SPINNER))]
            spin_idx += 1
            # \r overwrites the line; \e[K clears to end of line
            print(stderr, "\r\e[K  $spin_ch  Running... $(_format_elapsed(elapsed))")
        end

        sleep(2.0)
    end

    # Clear spinner line
    print(stderr, "\r\e[K")

    # Final flush — read any remaining output
    if isfile(path)
        open(path, "r") do io
            seek(io, last_pos)
            while !eof(io)
                line = readline(io; keep = false)
                _print_monitor_line(line, t_start)
            end
        end
    end
end

"""
Parse and print an ADJUST or TIME[NB] line for real-time monitoring.

ADJUST lines carry energy/virial info; TIME[NB] lines carry particle counts.
Both are printed with the wall-clock elapsed time prefix.

Returns `true` if a diagnostics line was printed, `false` otherwise.
"""
function _print_monitor_line(line::AbstractString, t_start::Float64 = time())::Bool
    stripped = lstrip(line)
    elapsed_str = _format_elapsed(time() - t_start)

    # ── ADJUST line (key-value or positional) ──
    if startswith(stripped, "ADJUST:")
        # Clear any spinner residue
        print(stderr, "\r\e[K")
        tokens = split(replace(stripped, r"^ADJUST:\s*" => ""))
        length(tokens) >= 4 || return false
        try
            if uppercase(tokens[1]) == "TIME"
                kv = Dict{String,String}()
                i = 1
                while i < length(tokens)
                    kv[uppercase(tokens[i])] = tokens[i + 1]
                    i += 2
                end
                t_nb = parse(Float64, get(kv, "TIME", "0"))
                t_myr = parse(Float64, get(kv, "T[MYR]", "0"))
                qvir = parse(Float64, get(kv, "Q", "0"))
                de = parse(Float64, get(kv, "DE", "0"))
                @info @sprintf(
                    "[%s]  t_NB=%.4f  t_Myr=%.1f  |ΔE/E|=%.2e  Q_vir=%.3f",
                    elapsed_str,
                    t_nb,
                    t_myr,
                    abs(de),
                    qvir
                )
            else
                length(tokens) >= 8 || return false
                t_nb = parse(Float64, tokens[1])
                t_myr = parse(Float64, tokens[2])
                qvir = parse(Float64, tokens[3])
                de = parse(Float64, tokens[4])
                n = parse(Int, tokens[6])
                @info @sprintf(
                    "[%s]  t_NB=%.4f  t_Myr=%.1f  N=%d  |ΔE/E|=%.2e  Q_vir=%.3f",
                    elapsed_str,
                    t_nb,
                    t_myr,
                    n,
                    abs(de),
                    qvir
                )
            end
        catch
        end
        return true
    end

    # ── TIME[NB] line: particle counts ──
    if startswith(stripped, "TIME[NB]")
        print(stderr, "\r\e[K")
        m_n = match(r"\bN\s+(\d+)", stripped)
        m_np = match(r"NPAIRS\s+(\d+)", stripped)
        if !isnothing(m_n)
            n = parse(Int, m_n.captures[1])
            np = isnothing(m_np) ? 0 : parse(Int, m_np.captures[1])
            @info @sprintf("[%s]  N=%d  N_pairs=%d", elapsed_str, n, np)
            return true
        end
    end

    return false
end

"""
    _effective_omp_threads(sim::SimulationConfig) -> Int

OpenMP thread count the backend will run with: `sim.omp_threads` when set,
otherwise an inherited positive `OMP_NUM_THREADS`, otherwise every logical
CPU (`Sys.CPU_THREADS`, the OpenMP runtime default).
"""
function _effective_omp_threads(sim::SimulationConfig)::Int
    sim.omp_threads > 0 && return sim.omp_threads
    inherited = tryparse(Int, get(ENV, "OMP_NUM_THREADS", ""))
    return (inherited === nothing || inherited < 1) ? Sys.CPU_THREADS : inherited
end

"""
    _reported_omp_threads(stdout_path) -> Union{Nothing,Int}

OpenMP thread count the backend itself reports at start-up (`nbody6.F`
prints `RANK: <r> OpenMP Number of Threads: <n>` from the OpenMP runtime).
`nothing` when the stdout capture is absent or carries no such line.
"""
function _reported_omp_threads(stdout_path::AbstractString)::Union{Nothing,Int}
    isfile(stdout_path) || return nothing
    for line in eachline(stdout_path)
        m = match(r"OpenMP Number of Threads:\s*(\d+)", line)
        m === nothing || return parse(Int, m.captures[1])
    end
    return nothing
end

"""
    _write_run_summary(cfg, run_dir, run_id, stdout_path, out_dir, elapsed; telemetry = nothing)

Write the machine-readable run summary `RUN_INFO.toml` into `run_dir`:
run identity, wall-clock time, and the backend thread layout (effective
OpenMP threads, MPI ranks); provenance (package and backend commits); the
hardware fingerprint (§6; GPU probed when `build.enable_gpu`); the
`[telemetry]` table from [`_finish_telemetry`](@ref) when given; and the
output file inventory.
"""
function _write_run_summary(
    cfg::Nbody6Config,
    run_dir::String,
    run_id::String,
    stdout_path::String,
    out_dir::String = run_dir,
    elapsed::Float64 = 0.0;
    telemetry::Union{Nothing,Dict{String,Any}} = nothing,
)
    run_table = Dict{String,Any}(
        "id" => run_id,
        "date" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        "elapsed_seconds" => round(elapsed; digits = 1),
        "elapsed" => _format_elapsed(elapsed),
        "omp_threads" => _effective_omp_threads(cfg.simulation),
        "mpi_ranks" => cfg.simulation.mpi_ranks,
    )
    isfile(stdout_path) && (run_table["stdout_lines"] = countlines(stdout_path))
    reported = _reported_omp_threads(stdout_path)
    reported === nothing || (run_table["omp_threads_reported"] = reported)

    d = Dict{String,Any}(
        "run" => run_table,
        "provenance" => Dict{String,Any}(
            "package_commit" => _git_commit(_PROJECT_ROOT),
            "backend_commit" =>
                _git_commit(joinpath(_PROJECT_ROOT, "backend", "Nbody6PPGPU-beijing")),
        ),
        "hardware" => _hardware_fingerprint(; gpu_probe = cfg.build.enable_gpu),
    )
    telemetry === nothing || (d["telemetry"] = telemetry)
    if isdir(out_dir)
        files = sort(readdir(out_dir))
        d["output"] = Dict{String,Any}(
            "files" => files,
            "total_bytes" => sum(f -> filesize(joinpath(out_dir, f)), files; init = 0),
        )
    end

    open(joinpath(run_dir, "RUN_INFO.toml"), "w") do io
        TOML.print(io, d)
    end
    return nothing
end
