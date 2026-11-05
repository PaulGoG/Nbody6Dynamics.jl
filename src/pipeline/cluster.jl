# =============================================================================
# Batch-scheduler layer: one Slurm job per engine segment
# =============================================================================
#
# A run longer than the queue's time limit is a chain of jobs, each running
# one engine segment: the first job starts the pipeline, every later one
# resumes it, and the one that reaches the end time post-processes. The job
# script is generated into `<run_dir>/cluster/` at submission; every
# submission and every job start and end is appended to a ledger.

"""
Columns of the job ledger (`cluster/jobs.csv` of a run, `cluster_jobs.csv` of
the runs directory), in file order.
"""
const _LEDGER_COLUMNS =
    (:time, :run_id, :event, :job_id, :depends_on, :segment, :status, :host, :note)

"""Row of the job ledger as [`read_job_ledger`](@ref) returns it: the nine columns as strings."""
const _LedgerRow = NamedTuple{_LEDGER_COLUMNS,NTuple{9,String}}

"""Run statuses after which no further job is submitted (complete, or a decision is needed)."""
const _FINAL_STATUSES = ("completed", "halted", "watchdog", "stalled")

"""Segment statuses of an engine that ended without a stop request or `END RUN`."""
const _FAILED_STATUSES = ("failed", "killed")

"""
    _chain_stop_reason(info_path) -> String

Why a `resubmit` chain must not submit another job after the last segment
recorded in `info_path` (`RUN_INFO.toml`), or `""` when it may. A segment
that ended `failed` or `killed` stops the chain when it made no progress
(its `t_end` is missing or not beyond its `t_start`), or when the segment
before it ended `failed` or `killed` too: a deterministic failure would
otherwise be resubmitted until `max_segments`, or forever.
"""
function _chain_stop_reason(info_path::AbstractString)::String
    isfile(info_path) || return ""
    segments = get(TOML.parsefile(info_path), "segments", Any[])
    isempty(segments) && return ""
    last_segment = segments[end]
    last_segment isa AbstractDict || return ""
    status = string(get(last_segment, "status", ""))
    status in _FAILED_STATUSES || return ""
    index = get(last_segment, "index", length(segments))
    t_start = get(last_segment, "t_start", 0.0)
    t_end = get(last_segment, "t_end", nothing)
    if !(t_end isa Real) || !(t_start isa Real) || t_end ≤ t_start
        return "segment $index ended $status without progress " *
               "(t_end = $(something(t_end, "none")), t_start = $t_start)"
    end
    if length(segments) ≥ 2
        previous = segments[end - 1]
        if previous isa AbstractDict && string(get(previous, "status", "")) in _FAILED_STATUSES
            return "segments $(get(previous, "index", length(segments) - 1)) and $index both " *
                   "ended failed or killed (latest: $status)"
        end
    end
    return ""
end

# ---------------------------------------------------------------------------
# Job script
# ---------------------------------------------------------------------------

"""
    _run_path(run_dir, base_dir) -> String

`run_dir` resolved against `base_dir` and made absolute, without a trailing
path separator, so that its `basename` is the run id. Raises an
`ArgumentError` for an empty `run_dir`.
"""
function _run_path(run_dir::AbstractString, base_dir::AbstractString)::String
    isempty(run_dir) && throw(ArgumentError("the run directory must be a nonempty path"))
    path = abspath(_resolve_path(base_dir, run_dir))
    return isdirpath(path) && path != "/" ? dirname(path) : path
end

"""
    _sbatch_path(path) -> String

`path` as the argument of an `#SBATCH` option: unchanged, or double-quoted
when it holds whitespace (Slurm's parser of `#SBATCH` lines honours double
quotes; shell single quotes are not read there).
"""
_sbatch_path(path::AbstractString)::String =
    occursin(r"\s", path) ? "\"" * String(path) * "\"" : String(path)

"""
    job_script(cfg, run_dir; base_dir = cfg.config_dir) -> String

Text of the `sbatch` script that runs one segment of the run in `run_dir`
(resolved against `base_dir` when relative) under the `[cluster]` settings of
`cfg`. Nothing is written; [`submit_run`](@ref) stores the text as
`<run_dir>/cluster/job.sbatch`.

The `#SBATCH` header carries the job name (`cluster.job_name`, else the run
id), the partition, account, QOS, memory and generic resources when set, the
node, task and CPU counts, `--exclusive` when requested, the time limit,
`--signal=B:TERM@<stop_margin>` (the signal time clamped to 1–65535 s, the
range Slurm accepts) and the output file `<run_dir>/cluster/slurm_%j.out`,
followed by the `cluster.sbatch_options` verbatim. The body runs the
`cluster.prologue` lines, sets `JULIA_NUM_THREADS` (`cluster.julia_threads`,
or the job's `SLURM_CPUS_PER_TASK` with `cluster.cpus_per_task` as the
fallback outside a job) and `exec`s Julia on the package's
`scripts/run_cluster_segment.jl`. The `exec` makes the Julia driver the batch
shell itself, so the signal Slurm sends ahead of the time limit reaches the
driver, whose exit hook stops the engine with a dump.

# Arguments
- `cfg`: configuration whose `[cluster]` and `simulation.stop_margin` are rendered
- `run_dir`: run directory (`<runs_dir>/<run_id>`)
- `base_dir`: directory a relative `run_dir` is resolved against

# Returns
The script text, ending in a newline.
"""
function job_script(
    cfg::Nbody6Config,
    run_dir::AbstractString;
    base_dir::AbstractString = cfg.config_dir,
)::String
    cl = cfg.cluster
    dir = _run_path(run_dir, base_dir)
    julia = isempty(cl.julia) ? joinpath(Sys.BINDIR, "julia") : cl.julia
    signal_time = clamp(round(Int, cfg.simulation.stop_margin), 1, 65535)
    threads =
        cl.julia_threads > 0 ? string(cl.julia_threads) :
        "\"\${SLURM_CPUS_PER_TASK:-$(cl.cpus_per_task)}\""
    segment_script = joinpath(_PACKAGE_ROOT, "scripts", "run_cluster_segment.jl")

    io = IOBuffer()
    println(io, "#!/bin/bash")
    println(io, "#SBATCH --job-name=", isempty(cl.job_name) ? basename(dir) : cl.job_name)
    for (option, value) in (
        ("partition", cl.partition),
        ("account", cl.account),
        ("qos", cl.qos),
        ("mem", cl.mem),
        ("gres", cl.gres),
    )
        isempty(value) || println(io, "#SBATCH --$option=$value")
    end
    println(io, "#SBATCH --nodes=$(cl.nodes)")
    println(io, "#SBATCH --ntasks=$(cl.ntasks)")
    println(io, "#SBATCH --cpus-per-task=$(cl.cpus_per_task)")
    cl.exclusive && println(io, "#SBATCH --exclusive")
    println(io, "#SBATCH --time=$(cl.time_limit)")
    println(io, "#SBATCH --signal=B:TERM@$(signal_time)")
    println(io, "#SBATCH --output=", _sbatch_path(joinpath(dir, "cluster", "slurm_%j.out")))
    for option in cl.sbatch_options
        println(io, "#SBATCH ", option)
    end
    println(io, "# Generated by Nbody6Dynamics.jl — DO NOT EDIT")
    println(io, "set -o pipefail")
    for line in cl.prologue
        println(io, line)
    end
    println(io, "export JULIA_NUM_THREADS=$threads")
    println(
        io,
        "exec $(_sh_quote(julia)) --startup-file=no $(_sh_quote(segment_script)) $(_sh_quote(dir))",
    )
    return String(take!(io))
end

"""
    _write_job_script(cfg, run_dir) -> String

Write [`job_script`](@ref) of the run in `run_dir` (absolute) to
`<run_dir>/cluster/job.sbatch` with mode 0o755; an existing script is kept as
a `#k` sibling. Returns the path.
"""
function _write_job_script(cfg::Nbody6Config, run_dir::AbstractString)::String
    path = joinpath(run_dir, "cluster", "job.sbatch")
    mkpath(dirname(path))
    _backup_existing(path)
    write(path, job_script(cfg, run_dir; base_dir = run_dir))
    chmod(path, 0o755)
    return path
end

# ---------------------------------------------------------------------------
# Submission
# ---------------------------------------------------------------------------

"""
    _sbatch_cmd(cl, script; after = "", kind = "afterok", kill_on_invalid = false) -> Cmd

Submission command `<sbatch> --parsable [--kill-on-invalid-dep=yes]
[--dependency=<kind>:<after>] <script>` for the `[cluster]` settings `cl`; the
dependency is omitted when `after` is empty.
"""
function _sbatch_cmd(
    cl::ClusterConfig,
    script::AbstractString;
    after::AbstractString = "",
    kind::AbstractString = "afterok",
    kill_on_invalid::Bool = false,
)::Cmd
    args = String[cl.sbatch, "--parsable"]
    kill_on_invalid && push!(args, "--kill-on-invalid-dep=yes")
    isempty(after) || push!(args, "--dependency=$(kind):$(after)")
    push!(args, String(script))
    return Cmd(args)
end

"""
    _parse_job_id(reply) -> Union{Nothing,String}

Job id in the standard output of `sbatch`: `Submitted batch job <id>`, or the
`<id>` / `<id>;<cluster>` line `--parsable` prints (the last nonempty line is
read). `nothing` when neither form is found.
"""
function _parse_job_id(reply::AbstractString)::Union{Nothing,String}
    m = match(r"Submitted batch job (\d+)", reply)
    m === nothing || return String(something(m.captures[1]))
    lines = filter(!isempty, strip.(split(reply, '\n')))
    isempty(lines) && return nothing
    m = match(r"^(\d+)(?:;\S+)?$", last(lines))
    return m === nothing ? nothing : String(something(m.captures[1]))
end

"""
    _submit(cmd) -> String

Run the submission command `cmd` and return the job id it reports
([`_parse_job_id`](@ref)). Raises an error naming the command when it cannot
be started, exits with a nonzero status (its standard error is quoted), or
replies with no job id.
"""
function _submit(cmd::Cmd)::String
    line = Base.shell_escape(cmd)
    out = IOBuffer()
    err = IOBuffer()
    proc = try
        run(pipeline(ignorestatus(cmd); stdout = out, stderr = err))
    catch e
        e isa Base.IOError || rethrow()
        error("job submission failed: could not start `$line`: $(e.msg)")
    end
    reply = String(take!(out))
    success(proc) || error(
        "job submission failed: `$line` exited with status $(proc.exitcode): " *
        strip(String(take!(err))),
    )
    id = _parse_job_id(reply)
    id === nothing &&
        error("job submission failed: `$line` replied without a job id: \"$(strip(reply))\"")
    return id
end

"""
    submit_run(cfg; base_dir = cfg.config_dir, run_id = "", dry_run = false)
        -> @NamedTuple{run_dir::String, script::String, commands::Vector{String}, job_ids::Vector{String}}

Create a run and submit its chain of scheduler jobs, one engine segment per
job, under the `[cluster]` settings of `cfg` (`cluster.enabled` must be set;
the `[cluster]` checks of the configuration are applied again, so a
configuration built in memory is held to them too).

The run directory is `<runs_dir>/<run_id>` (`run_id` generated as by
[`run_pipeline`](@ref) when empty); a directory that already holds a
`RUN_INFO.toml` is refused, since a run in progress is continued by its own
chain or by [`resume_pipeline`](@ref). The configuration is frozen with
absolute paths into `<run_dir>/cluster/config.toml`, which the first job
reads, and the job script ([`job_script`](@ref)) is written to
`<run_dir>/cluster/job.sbatch`.

With `chain = "dependency"`, `max_segments` jobs are submitted at once, each
after the first with `--dependency=afterok:<previous job>` and
`--kill-on-invalid-dep=yes`, so a failed segment cancels the rest of the
chain. With `chain = "resubmit"` one job is submitted, and each job submits
its successor while the run is short of its end time
([`cluster_segment`](@ref)). Every submission is appended to the job ledgers
(`<run_dir>/cluster/jobs.csv` and `<runs_dir>/cluster_jobs.csv`, see
[`read_job_ledger`](@ref)) as soon as it is made, so a failure part-way
through the chain leaves the submissions already made on record.

With `dry_run`, nothing is submitted and no ledger row is written; the
commands are returned with `<job k>` placeholders for the dependency ids.

# Arguments
- `cfg`: the run's configuration, with `cluster.enabled = true`
- `base_dir`: directory the relative paths of `cfg` resolve against
- `run_id`: name of the run directory; empty = generated
- `dry_run`: render and write the run's files without submitting

# Returns
A named tuple: `run_dir`, `script` (path of `job.sbatch`), `commands` (the
submission command lines) and `job_ids` (empty for a dry run).

# Example
```julia
cfg = load_config("run.toml")    # [cluster] enabled = true
sub = submit_run(cfg; dry_run = true)
foreach(println, sub.commands)
sub = submit_run(cfg)
sub.job_ids                           # e.g. ["4711"]
```
"""
function submit_run(
    cfg::Nbody6Config;
    base_dir::AbstractString = cfg.config_dir,
    run_id::AbstractString = "",
    dry_run::Bool = false,
)
    cl = cfg.cluster
    cl.enabled || throw(
        ArgumentError(
            "submit_run: cluster.enabled is false; set `[cluster] enabled = true` to submit " *
            "the run as scheduler jobs, or call run_pipeline to run it in this process",
        ),
    )
    _validate_cluster(cfg)
    isempty(run_id) && (
        run_id = generate_run_id(
            (cfg.merger.enabled ? "merger_" : "") * cfg.simulation.run_id_prefix,
        )
    )
    runs_dir = _resolve_path(base_dir, cfg.simulation.runs_dir)
    run_dir = joinpath(runs_dir, run_id)
    isfile(joinpath(run_dir, "RUN_INFO.toml")) && throw(
        ArgumentError(
            "submit_run: $run_dir already holds a run (RUN_INFO.toml); a run in progress is " *
            "continued by its own job chain or by resume_pipeline, not submitted again",
        ),
    )

    mkpath(joinpath(run_dir, "cluster"))
    frozen = _with_absolute_paths(cfg, base_dir)
    save_config(frozen, joinpath(run_dir, "cluster", "config.toml"))
    script = _write_job_script(frozen, run_dir)

    n_jobs = cl.chain == "dependency" ? cl.max_segments : 1
    commands = String[]
    job_ids = String[]
    previous = ""
    for k in 1:n_jobs
        after = k == 1 ? "" : dry_run ? "<job $(k - 1)>" : previous
        cmd = _sbatch_cmd(cl, script; after = after, kill_on_invalid = k > 1)
        line = Base.shell_escape(cmd)
        push!(commands, line)
        dry_run && continue
        id = _submit(cmd)
        push!(job_ids, id)
        _record_job(
            run_dir,
            runs_dir;
            event = "submitted",
            job_id = id,
            depends_on = previous,
            note = line,
        )
        @info "Submitted job $id for run $run_id" * (isempty(previous) ? "" : " (after $previous)")
        previous = id
    end
    return (run_dir = String(run_dir), script = script, commands = commands, job_ids = job_ids)
end

"""Request file of a restart job (`<run_dir>/cluster/restart.toml`)."""
const _RESTART_REQUEST = joinpath("cluster", "restart.toml")

"""
    _read_restart_request(run_dir) -> Union{Nothing,Dict{String,Any}}

The restart request of `run_dir` ([`submit_restart`](@ref)) when one is pending, else
`nothing`.
"""
function _read_restart_request(run_dir::AbstractString)::Union{Nothing,Dict{String,Any}}
    path = joinpath(run_dir, _RESTART_REQUEST)
    isfile(path) || return nothing
    req = TOML.parsefile(path)
    return get(req, "status", "") == "pending" ? req : nothing
end

"""
    _consume_restart_request(run_dir, request, job_id)

Mark the restart request of `run_dir` consumed by the job `job_id` before its restart
starts: a job ended during the restart segment is continued by a resume of that segment,
not by the same restart again.
"""
function _consume_restart_request(
    run_dir::AbstractString,
    request::AbstractDict,
    job_id::AbstractString,
)
    consumed = merge(
        Dict{String,Any}(request),
        Dict{String,Any}(
            "status" => "consumed",
            "consumed" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
            "job_id" => String(job_id),
        ),
    )
    _atomic_write_toml(joinpath(run_dir, _RESTART_REQUEST), consumed)
    return nothing
end

"""
    submit_restart(run_dir; t_end = nothing, t_end_myr = nothing, tcrit_extra = nothing,
                   dump = nothing, kz = nothing, engine = :run, dry_run = false)
        -> @NamedTuple{run_dir::String, script::String, command::String, job_id::String,
                       request::Dict{String,Any}}

Submit a job chain whose first segment restarts the run in `run_dir` from a dump
([`restart_pipeline`](@ref)): the continuation of a completed run beyond its end time
(an extension), or of a run halted by the energy check or stopped by the stall watchdog
from an earlier dump, possibly with changed options (`kz`). The run's frozen
`config.toml` must have `[cluster] enabled = true`; its `[cluster]` settings govern the
chain, which continues with [`cluster_segment`](@ref) as for a run from
[`submit_run`](@ref).

The end of the restart is `t_end` [NB], `t_end_myr` [Myr] (converted with the run's time
unit from `global.30`), or `tcrit_extra` [NB] after the dump; with none of them, the
run's end time (its input's `TCRIT`, or a later end given to an earlier restart). The
dump is `dump` (a file name in `output/`) or the one a resume would take (the newest
complete restart dump before any halt). Both are resolved and checked here — the dump
complete and outside a chain subsystem (`NSUB = 0`), the span positive, no engine of the
run alive, no other request pending — and written to `<run_dir>/cluster/restart.toml`
with status `pending`; the first job of the chain consumes the request before it starts
the engine, so a job ended mid-segment is continued by a resume, not by a second restart.
A consumed request of an earlier restart is kept as `cluster/restart#k.toml`.

With `dry_run` the request is written to `cluster/restart_preview.toml` and the restart
input to `cluster/restart_preview.inp` instead, and nothing is submitted.

# Example
```julia
# continue a completed run to 680 Myr on the cluster
submit_restart("runs/merger_i8a_seed1"; t_end_myr = 680.0)
```
"""
function submit_restart(
    run_dir::AbstractString;
    t_end::Union{Nothing,Real} = nothing,
    t_end_myr::Union{Nothing,Real} = nothing,
    tcrit_extra::Union{Nothing,Real} = nothing,
    dump::Union{Nothing,AbstractString} = nothing,
    kz::Union{Nothing,AbstractDict{<:Integer,<:Integer}} = nothing,
    engine::Symbol = :run,
    dry_run::Bool = false,
)
    run_dir = _run_path(run_dir, pwd())
    count(!isnothing, (t_end, t_end_myr, tcrit_extra)) ≤ 1 ||
        throw(ArgumentError("submit_restart: give at most one of t_end, t_end_myr, tcrit_extra"))
    engine in (:run, :backend) ||
        throw(ArgumentError("submit_restart: engine must be :run or :backend; got :$engine"))
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    isfile(info_path) || throw(
        ArgumentError("submit_restart: $info_path not found; only pipeline runs can be restarted"),
    )
    cfg = load_config(joinpath(run_dir, "config.toml"))
    cfg.cluster.enabled || throw(
        ArgumentError(
            "submit_restart: cluster.enabled is false in $(joinpath(run_dir, "config.toml")); " *
            "set the run's [cluster] section first, or restart in this process with " *
            "restart_pipeline",
        ),
    )
    _validate_cluster(cfg)
    _read_restart_request(run_dir) === nothing || throw(
        ArgumentError(
            "submit_restart: $(joinpath(run_dir, _RESTART_REQUEST)) holds a pending request; " *
            "a job of the run has not consumed it yet",
        ),
    )
    info = TOML.parsefile(info_path)
    segments = get(info, "segments", Any[])
    out_dir = joinpath(run_dir, "output")
    if !isempty(segments) &&
       get(segments[end], "status", "") == "running" &&
       _segment_alive(segments[end], out_dir)
        throw(ArgumentError("submit_restart: the engine of $(basename(run_dir)) is running"))
    end
    stdout_path = joinpath(out_dir, cfg.postprocess.stdout_file)
    chosen = if dump === nothing
        marker = _resume_dump(out_dir, stdout_path)
        marker === nothing &&
            throw(ArgumentError("submit_restart: no complete restart dump in $out_dir"))
        marker.file
    else
        String(dump)
    end
    dump_path = joinpath(out_dir, chosen)
    isfile(dump_path) || throw(ArgumentError("submit_restart: dump not found: $dump_path"))
    _dump_complete(dump_path) ||
        throw(ArgumentError("submit_restart: $chosen is incomplete (truncated write)"))
    _restartable(dump_path) || throw(
        ArgumentError(
            "submit_restart: $chosen was written inside a chain subsystem (NSUB > 0); " *
            "choose an earlier dump",
        ),
    )
    markers = _dump_markers(stdout_path)
    k = findlast(m -> m.file == chosen, markers)
    t_dump = k === nothing ? _dump_time(chosen) : markers[k].time_nb
    original_inp = joinpath(out_dir, get(get(info, "run", Dict{String,Any}()), "input_file", ""))
    isfile(original_inp) || throw(
        ArgumentError("submit_restart: the original input (run.input_file) is missing in $out_dir"),
    )
    target = if t_end !== nothing
        Float64(t_end)
    elseif t_end_myr !== nothing
        unit = _time_unit_myr(out_dir)
        isfinite(unit) && unit > 0 || throw(
            ArgumentError("submit_restart: no time unit in $(joinpath(out_dir, "global.30"))"),
        )
        Float64(t_end_myr) / unit
    elseif tcrit_extra !== nothing
        t_dump + Float64(tcrit_extra)
    else
        _run_end_time(info, original_inp)
    end
    span = target - t_dump
    span > 0 || throw(
        ArgumentError(
            "submit_restart: the end time $(target) NB is not after the dump $(chosen) " *
            "(t = $(t_dump) NB)",
        ),
    )
    request = Dict{String,Any}(
        "status" => "pending",
        "requested" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        "dump" => chosen,
        "t_start" => t_dump,
        "tcrit_extra" => span,
        "t_end" => target,
        "engine" => String(engine),
        "kz" =>
            kz === nothing ? Dict{String,Any}() :
            Dict{String,Any}(string(i) => Int(v) for (i, v) in kz),
    )
    cluster_dir = joinpath(run_dir, "cluster")
    mkpath(cluster_dir)
    script = _write_job_script(cfg, run_dir)
    cmd = _sbatch_cmd(cfg.cluster, script)
    line = Base.shell_escape(cmd)
    if dry_run
        _atomic_write_toml(joinpath(cluster_dir, "restart_preview.toml"), request)
        _write_restart_inp(
            joinpath(cluster_dir, "restart_preview.inp"),
            original_inp,
            span;
            kz = kz,
        )
        return (run_dir = run_dir, script = script, command = line, job_id = "", request = request)
    end
    # The consumed request of an earlier restart is kept as cluster/restart#k.toml.
    _backup_existing(joinpath(run_dir, _RESTART_REQUEST))
    _atomic_write_toml(joinpath(run_dir, _RESTART_REQUEST), request)
    runs_dir = _resolve_path(cfg.config_dir, cfg.simulation.runs_dir)
    id = _submit(cmd)
    _record_job(
        run_dir,
        runs_dir;
        event = "submitted",
        job_id = id,
        note = "restart from $(chosen) to t = $(target): " * line,
    )
    @info "Submitted job $id: restart of $(basename(run_dir)) from $chosen (t = $t_dump) to t = $target"
    return (run_dir = run_dir, script = script, command = line, job_id = id, request = request)
end

# ---------------------------------------------------------------------------
# The job's entry
# ---------------------------------------------------------------------------

"""
    _save_job_record(run_dir, job_id)

Keep the scheduler's description of the job (`scontrol show jobid -dd`) as
`<run_dir>/cluster/job_<job_id>.txt`. Does nothing outside a Slurm job (empty
`job_id`) or without `scontrol` on the `PATH`; a failing `scontrol` is
reported as a warning and does not end the job.
"""
function _save_job_record(run_dir::AbstractString, job_id::AbstractString)
    (isempty(job_id) || Sys.which("scontrol") === nothing) && return nothing
    path = joinpath(run_dir, "cluster", "job_$(job_id).txt")
    try
        mkpath(dirname(path))
        write(path, read(`scontrol show jobid -dd $job_id`))
    catch e
        e isa Union{ProcessFailedException,Base.IOError,SystemError} || rethrow()
        @warn "scontrol show jobid $job_id failed; the job record is not kept" exception = e
    end
    return nothing
end

"""
    cluster_segment(run_dir) -> Dict{Symbol,Any}

Entry of a scheduler job (`scripts/run_cluster_segment.jl`): run one engine
segment of the run in `run_dir` and, with `chain = "resubmit"`, submit the
job of the next one.

The configuration is read from `<run_dir>/config.toml` when it exists, else
from `<run_dir>/cluster/config.toml`. The latter is the copy frozen at
submission and serves only the first job; from the first segment on,
`<run_dir>/config.toml` holds the whole configuration, `[cluster]` included,
and is the single source of truth for both the engine segment and the chain
(`chain`, `max_segments`, `sbatch`, the job script). Changes between jobs,
such as a different `simulation.wall_budget`, are therefore made in
`<run_dir>/config.toml`.

A pending restart request ([`submit_restart`](@ref)) is taken first: it is
marked consumed with the job id, then the restart segment runs
([`restart_pipeline`](@ref)). Otherwise, without a `RUN_INFO.toml` the
pipeline is started ([`run_pipeline`](@ref) with the run id
`basename(run_dir)`), and with one it is continued
([`resume_pipeline`](@ref)). The job id is taken from `SLURM_JOB_ID` (empty
outside a job); the ledgers receive a `started` row before the segment and an
`ended` row with the run status and the segment count after it, or with
status `error` and the exception's message before the exception propagates.
The scheduler's job description is kept as `cluster/job_<job_id>.txt` when
`scontrol` is available.

A successor `<sbatch> --parsable --dependency=afterany:<job_id>
cluster/job.sbatch` is submitted when `chain = "resubmit"`, the run status is
none of `completed`, `halted`, `watchdog`, `stalled`, the pipeline has not
finished, and the segment count is below `max_segments` (`0` = no bound); the
script is regenerated from the configuration if it is missing. A segment that
ended `failed` or `killed` without progress, or the second such segment in a
row, stops the chain instead (`_chain_stop_reason`): a warning names
the run, status and segment, and the ledgers receive a `chain_stopped` row
with the reason.

# Returns
The result of the pipeline call: empty while the run is short of its end
time, and once its pipeline has finished.

# Example
```julia
# what the generated job script runs, within a Slurm job
results = cluster_segment("/scratch/runs/merger_run_20261013_120000_ab12")
```
"""
function cluster_segment(run_dir::AbstractString)::Dict{Symbol,Any}
    run_dir = _run_path(run_dir, pwd())
    run_config = joinpath(run_dir, "config.toml")
    submit_config = joinpath(run_dir, "cluster", "config.toml")
    cfg_path = isfile(run_config) ? run_config : submit_config
    isfile(cfg_path) || throw(
        ArgumentError(
            "cluster_segment: neither $run_config nor $submit_config exists; the run directory " *
            "was not created by submit_run",
        ),
    )
    cfg = load_config(cfg_path)
    cl = cfg.cluster
    job_id = get(ENV, "SLURM_JOB_ID", "")
    runs_dir = _resolve_path(cfg.config_dir, cfg.simulation.runs_dir)
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    host = gethostname()

    _record_job(
        run_dir,
        runs_dir;
        event = "started",
        job_id = job_id,
        segment = string(_segment_count(info_path)),
        host = host,
    )
    _save_job_record(run_dir, job_id)

    request = _read_restart_request(run_dir)
    result = try
        if request !== nothing
            _consume_restart_request(run_dir, request, job_id)
            kz = Dict{Int,Int}(parse(Int, k) => Int(v) for (k, v) in get(request, "kz", Dict()))
            restart_pipeline(
                run_dir;
                tcrit_extra = Float64(request["tcrit_extra"]),
                dump = String(request["dump"]),
                kz = isempty(kz) ? nothing : kz,
                engine = Symbol(get(request, "engine", "run")),
                base_dir = cfg.config_dir,
            )
        elseif isfile(info_path)
            resume_pipeline(run_dir; base_dir = cfg.config_dir)
        else
            run_pipeline(cfg; base_dir = cfg.config_dir, run_id = basename(run_dir))
        end
    catch e
        _record_job(
            run_dir,
            runs_dir;
            event = "ended",
            job_id = job_id,
            segment = string(_segment_count(info_path)),
            status = "error",
            host = host,
            note = sprint(showerror, e),
        )
        rethrow()
    end

    status = _run_status(run_dir)
    segments = _segment_count(info_path)
    _record_job(
        run_dir,
        runs_dir;
        event = "ended",
        job_id = job_id,
        segment = string(segments),
        status = status,
        host = host,
    )

    continuing =
        cl.chain == "resubmit" && !(status in _FINAL_STATUSES) && !_pipeline_completed(run_dir)
    stop_reason = continuing ? _chain_stop_reason(info_path) : ""
    if !isempty(stop_reason)
        _record_job(
            run_dir,
            runs_dir;
            event = "chain_stopped",
            job_id = job_id,
            segment = string(segments),
            status = status,
            host = host,
            note = stop_reason,
        )
        @warn "Run $(basename(run_dir)) is at status $(status) after segment $(segments): " *
              "$(stop_reason); the chain is stopped and no further job is submitted " *
              "(inspect the run, then resume it with resume_pipeline or a new job)"
    elseif continuing && (cl.max_segments == 0 || segments < cl.max_segments)
        script = joinpath(run_dir, "cluster", "job.sbatch")
        isfile(script) || _write_job_script(cfg, run_dir)
        cmd = _sbatch_cmd(cl, script; after = job_id, kind = "afterany")
        line = Base.shell_escape(cmd)
        id = _submit(cmd)
        _record_job(
            run_dir,
            runs_dir;
            event = "submitted",
            job_id = id,
            depends_on = job_id,
            host = host,
            note = (isempty(job_id) ? "resubmitted: " : "resubmitted by job $(job_id): ") * line,
        )
        @info "Run $(basename(run_dir)) is at status $(status) after $(segments) segment(s); " *
              "job $(id) submitted for the next segment"
    elseif continuing
        @info "Run $(basename(run_dir)) is at status $(status) after $(segments) segment(s); " *
              "cluster.max_segments = $(cl.max_segments) reached, no further job submitted"
    end
    return result
end

# ---------------------------------------------------------------------------
# Ledger
# ---------------------------------------------------------------------------

"""
    _csv_field(value) -> String

`value` as one CSV field: unchanged, or enclosed in double quotes with inner
quotes doubled when it holds a comma, a double quote or a line break.
"""
function _csv_field(value)::String
    s = string(value)
    occursin(r"[,\"\n\r]", s) || return s
    return "\"" * replace(s, "\"" => "\"\"") * "\""
end

"""
    _append_ledger(path, row)

Append `row` (a [`_LedgerRow`](@ref)) to the ledger at `path`, writing the
header first when the file is created. The row goes out in one `write`, so
rows of concurrent jobs do not interleave.
"""
function _append_ledger(path::AbstractString, row::_LedgerRow)
    mkpath(dirname(path))
    fresh = !isfile(path) || filesize(path) == 0
    line = join((_csv_field(getfield(row, c)) for c in _LEDGER_COLUMNS), ',') * "\n"
    open(path, "a") do io
        fresh && write(io, join(_LEDGER_COLUMNS, ',') * "\n")
        write(io, line)
    end
    return nothing
end

"""
    _record_job(run_dir, runs_dir; event, job_id = "", depends_on = "", segment = "",
                status = "", host = gethostname(), note = "")

Append one row, stamped with the current time and the run id
`basename(run_dir)`, to both job ledgers: `<run_dir>/cluster/jobs.csv` and
`<runs_dir>/cluster_jobs.csv`. `event` is `submitted`, `started`, `ended` or
`chain_stopped`.
"""
function _record_job(
    run_dir::AbstractString,
    runs_dir::AbstractString;
    event::AbstractString,
    job_id::AbstractString = "",
    depends_on::AbstractString = "",
    segment::AbstractString = "",
    status::AbstractString = "",
    host::AbstractString = gethostname(),
    note::AbstractString = "",
)
    row = _LedgerRow((
        Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        String(basename(run_dir)),
        String(event),
        String(job_id),
        String(depends_on),
        String(segment),
        String(status),
        String(host),
        String(note),
    ))
    _append_ledger(joinpath(run_dir, "cluster", "jobs.csv"), row)
    _append_ledger(joinpath(runs_dir, "cluster_jobs.csv"), row)
    return nothing
end

"""
    _csv_records(text, path) -> Vector{Vector{String}}

Records of the CSV `text`: fields separated by commas, records by line breaks
(`\\n` or `\\r\\n`); a field enclosed in double quotes may hold commas, line
breaks and doubled quotes. Empty lines are skipped. Raises an
`ArgumentError` naming `path` for a quoted field that is not closed.
"""
function _csv_records(text::AbstractString, path::AbstractString)::Vector{Vector{String}}
    records = Vector{String}[]
    fields = String[]
    buf = IOBuffer()
    quoted = false
    content = false     # the current record holds a character or a separator
    chars = collect(text)
    n = length(chars)
    k = 1
    while k ≤ n
        c = chars[k]
        if quoted
            if c == '"' && k < n && chars[k + 1] == '"'
                write(buf, '"')
                k += 1
            elseif c == '"'
                quoted = false
            else
                write(buf, c)
            end
        elseif c == '"'
            quoted = true
            content = true
        elseif c == ','
            push!(fields, String(take!(buf)))
            content = true
        elseif c == '\n' || c == '\r'
            c == '\r' && k < n && chars[k + 1] == '\n' && (k += 1)
            if content
                push!(fields, String(take!(buf)))
                push!(records, fields)
                fields = String[]
            end
            content = false
        else
            write(buf, c)
            content = true
        end
        k += 1
    end
    quoted && throw(ArgumentError("job ledger $path ends inside a quoted field"))
    if content
        push!(fields, String(take!(buf)))
        push!(records, fields)
    end
    return records
end

"""
    read_job_ledger(path) -> Vector{NamedTuple}

Rows of a job ledger written by [`submit_run`](@ref) and
[`cluster_segment`](@ref) (`<run_dir>/cluster/jobs.csv`, or
`<runs_dir>/cluster_jobs.csv` for every run of a runs directory), as named
tuples (a `Vector` of `NamedTuple`s, not a table) with the nine `String`
fields `time`, `run_id`, `event` (`submitted`, `started`, `ended`,
`chain_stopped`), `job_id`, `depends_on`, `segment`, `status`, `host` and
`note`; empty cells read as `""`. An absent file gives an empty vector. Raises
an `ArgumentError` when the header is not the nine column names in this
order or a row has another number of cells.

# Example
```julia
rows = read_job_ledger(joinpath(run_dir, "cluster", "jobs.csv"))
[r.job_id for r in rows if r.event == "submitted"]
```
"""
function read_job_ledger(path::AbstractString)::Vector{_LedgerRow}
    rows = _LedgerRow[]
    isfile(path) || return rows
    records = _csv_records(read(path, String), path)
    isempty(records) && return rows
    header = Symbol.(strip.(records[1]))
    Tuple(header) == _LEDGER_COLUMNS || throw(
        ArgumentError(
            "job ledger $path has the header $(join(header, ',')); expected " *
            join(_LEDGER_COLUMNS, ','),
        ),
    )
    for (k, cells) in enumerate(records[2:end])
        length(cells) == length(_LEDGER_COLUMNS) || throw(
            ArgumentError(
                "job ledger $path: row $k has $(length(cells)) cells for " *
                "$(length(_LEDGER_COLUMNS)) columns",
            ),
        )
        push!(rows, _LedgerRow(Tuple(cells)))
    end
    return rows
end
