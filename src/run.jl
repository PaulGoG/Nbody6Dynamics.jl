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
    run_simulation(cfg::Nbody6Config; base_dir = cfg.config_dir, run_id = "") -> String

Execute the Nbody6++ simulation with:
- Unique run ID and isolated output directory
- Bash wrapper for `ulimit -s unlimited` and CUDA env
- Real-time ADJUST line monitoring with wall-clock timer and heartbeat spinner
- Frozen config snapshot saved alongside output

`simulation.input_file`, `install.install_dir` and `simulation.runs_dir`
resolve against `base_dir` when relative. Returns the path to the run
directory.
"""
function run_simulation(
    cfg::Nbody6Config;
    base_dir::AbstractString = cfg.config_dir,
    run_id::AbstractString = "",
)::String
    sim = cfg.simulation

    # --- Resolve the input file against the project directory ---
    input_path = _resolve_path(base_dir, sim.input_file)
    isfile(input_path) || error("Input file not found: $input_path")

    # --- Generate run ID and create directory structure ---
    #   runs/<run_id>/
    #       output/      — simulation output files
    #       plots/       — post-processing plots (created later)
    #       config.toml  — frozen config
    isempty(run_id) && (run_id = generate_run_id(sim.run_id_prefix))
    run_dir = joinpath(_resolve_path(base_dir, sim.runs_dir), run_id)
    out_dir = joinpath(run_dir, "output")
    mkpath(out_dir)

    @info "Run ID:  $run_id"
    @info "Run dir: $run_dir"

    return _execute_simulation(cfg, run_dir, out_dir, input_path; base_dir = base_dir)
end

# ---------------------------------------------------------------------------
# Restarts from the engine's COMMON dumps
# ---------------------------------------------------------------------------

"""
    _dump_time(name) -> Float64

Time suffix of a COMMON dump file name (`comm.1_12.5` → 12.5); `NaN` for
names that are not dumps.
"""
function _dump_time(name::AbstractString)::Float64
    m = match(r"^comm\.[12]_([0-9.eE+-]+)$", name)
    m === nothing && return NaN
    return something(tryparse(Float64, m.captures[1]), NaN)
end

"""
    _latest_dump(out_dir) -> Union{Nothing,String}

File name of the COMMON dump with the largest time suffix in `out_dir`
(`comm.1_<t>` at the end of a run or at a stop request, `comm.2_<t>` at
every adjustment in checkpoint mode), or `nothing` when none exists. The
time in a name is rounded to the digits of `DTADJ`.
"""
function _latest_dump(out_dir::AbstractString)::Union{Nothing,String}
    isdir(out_dir) || return nothing
    best = nothing
    best_t = -Inf
    for f in readdir(out_dir)
        t = _dump_time(f)
        isnan(t) && continue
        if t > best_t
            best_t = t
            best = f
        end
    end
    return best
end

"""
    _write_restart_inp(path, original_inp, tcrit_extra; tcrtp0 = nothing)

Write the restart input for `KSTART = 2` from the original input file: the
`&INNBODY6` block with `KSTART=2` (and `TCRTP0` replaced when given) and
the original `&ININPUT` block with `TCRIT` set to `tcrit_extra`, which the
engine adds to the saved time (`modify.F`: `TCRIT = TTOT + TCRIT`). Later
namelists are not read on restart and are omitted.

The increment is written with round-trip precision: the engine ends a run at
the first adjustment beyond `TCRIT − 20 DTMIN`, and four decimals moved the
end time by more than that at large N.
"""
function _write_restart_inp(
    path::AbstractString,
    original_inp::AbstractString,
    tcrit_extra::Real;
    tcrtp0::Union{Nothing,Real} = nothing,
)
    tcrit_extra > 0 || throw(ArgumentError("tcrit_extra must be > 0; got $tcrit_extra"))
    text = read(original_inp, String)
    m6 = match(r"&INNBODY6\s*\n(.*?)/", text)
    mi = match(r"&ININPUT\s*\n(.*?Level='[^']*')\s*/"ms, text)
    (m6 === nothing || mi === nothing) &&
        error("restart: could not locate &INNBODY6 and &ININPUT blocks in $original_inp")
    b6 = replace(m6.captures[1], r"KSTART\s*=\s*\d+" => "KSTART=2")
    if tcrtp0 !== nothing
        b6 = replace(b6, r"TCRTP0\s*=\s*[0-9.eE+-]+" => @sprintf("TCRTP0=%.6G", tcrtp0))
    end
    bi =
        replace(mi.captures[1], r"TCRIT\s*=\s*[0-9.eE+-]+" => "TCRIT=" * repr(Float64(tcrit_extra)))
    occursin("TCRIT=", bi) ||
        error("restart: no TCRIT entry found in the &ININPUT block of $original_inp")
    open(path, "w") do io
        println(io, "&INNBODY6")
        print(io, rstrip(b6), " /\n\n")
        println(io, "&ININPUT")
        print(io, rstrip(bi), " /\n")
    end
    return path
end

"""
    restart_simulation(run_dir; tcrit_extra, dump = nothing, tcrtp0 = nothing,
                       base_dir = nothing) -> String

Continue a finished or interrupted run from one of the engine's COMMON
dumps for `tcrit_extra` further N-body time units. The chosen dump
(default: the latest `comm.[12]_<t>` in `output/`) is copied to
`output/comm.1`, which is what `KSTART = 2` reads; a restart input is
written from the run's original input file; the engine runs in the same
output directory with stdout and stderr appended, so `out1000`, `lagr.7`,
`esc.11`, and the time-stamped snapshot and stellar-evolution files
continue. `RUN_INFO.toml` gains one entry in its `segments` list per
launch and the telemetry of each segment goes to its own CSV. The frozen
`config.toml` of the run carries absolute paths, so the engine tree is found
without a `base_dir`; one may still be given to override it. Returns
`run_dir`.
"""
function restart_simulation(
    run_dir::AbstractString;
    tcrit_extra::Real,
    dump::Union{Nothing,AbstractString} = nothing,
    tcrtp0::Union{Nothing,Real} = nothing,
    base_dir::Union{Nothing,AbstractString} = nothing,
)::String
    run_dir = abspath(run_dir)
    out_dir = joinpath(run_dir, "output")
    isdir(out_dir) || error("restart: $out_dir does not exist")
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    isfile(info_path) || error("restart: $info_path not found; only pipeline runs can be restarted")
    info = TOML.parsefile(info_path)
    original_inp = joinpath(out_dir, get(info["run"], "input_file", ""))
    isfile(original_inp) || error(
        "restart: original input file not recorded or missing (run.input_file in RUN_INFO.toml)",
    )
    cfg = load_config(joinpath(run_dir, "config.toml"))
    base_dir === nothing && (base_dir = cfg.config_dir)

    chosen = dump === nothing ? _latest_dump(out_dir) : String(dump)
    chosen === nothing && error("restart: no COMMON dump (comm.[12]_<t>) in $out_dir")
    dump_path = joinpath(out_dir, chosen)
    isfile(dump_path) || error("restart: dump not found: $dump_path")
    target = joinpath(out_dir, "comm.1")
    _backup_existing(target)
    cp(dump_path, target; force = true)
    @info "Restart from $chosen (t = $(_dump_time(chosen))) for $tcrit_extra more N-body time units"

    restart_inp = joinpath(out_dir, "restart.inp")
    _backup_existing(restart_inp)
    _write_restart_inp(restart_inp, original_inp, tcrit_extra; tcrtp0 = tcrtp0)

    return _execute_simulation(
        cfg,
        run_dir,
        out_dir,
        restart_inp;
        base_dir = base_dir,
        label = "restart",
        restart = (dump = chosen, tcrit_extra = Float64(tcrit_extra)),
    )
end

"""
    _execute_simulation(cfg, run_dir, out_dir, input_path;
                        base_dir = cfg.config_dir, label = "simulation",
                        restart = nothing) -> String

Shared execution core for [`run_simulation`](@ref), the merger pipeline,
and [`restart_simulation`](@ref): locates the binary, freezes the config
into `run_dir` (with its paths made absolute against `base_dir`, see
[`_with_absolute_paths`](@ref)), copies the binary and the input file into
`out_dir` for reproducibility, writes the launch script, and runs it under
the teed run log with the opt-in live monitor. `input_path` must be
absolute (the launch script executes from `out_dir`). With `restart = (;
dump, tcrit_extra)` the binary copy is reused, stdout and stderr are
appended, and the run summary records a further segment. Returns `run_dir`.
"""
function _execute_simulation(
    cfg::Nbody6Config,
    run_dir::AbstractString,
    out_dir::AbstractString,
    input_path::AbstractString;
    base_dir::AbstractString = cfg.config_dir,
    label::AbstractString = "simulation",
    restart::Union{Nothing,NamedTuple} = nothing,
)::String
    sim = cfg.simulation
    run_id = basename(run_dir)
    is_restart = restart !== nothing

    # --- Locate binary ---
    src_dir = _resolve_path(base_dir, cfg.install.install_dir)
    binary = _find_binary(src_dir, sim.binary_name, cfg.build)
    # A restart takes its sizes from the dump, which the same build wrote.
    is_restart || _check_engine_limits(input_path, _engine_limits(src_dir))

    # --- Provenance at launch: the record names the code that ran, not the
    # tree found when the record is written hours later (a `git pull` during
    # a long run used to change the recorded commit).
    provenance = Dict{String,Any}(
        "package_commit" => _source_stamp(_PACKAGE_ROOT),
        "backend_commit" => isempty(src_dir) ? "unknown" : _git_commit(src_dir),
    )

    # --- Save frozen config (absolute paths: loadable from anywhere) ---
    is_restart || save_config(_with_absolute_paths(cfg, base_dir), joinpath(run_dir, "config.toml"))

    # --- Copy binary, build record and input for reproducibility (kept on restart) ---
    local_binary = joinpath(out_dir, basename(binary))
    if !(is_restart && isfile(local_binary))
        cp(binary, local_binary; force = true)
        chmod(local_binary, 0o755)
        build_info = joinpath(dirname(binary), "BUILD_INFO.toml")
        isfile(build_info) && cp(build_info, joinpath(out_dir, "BUILD_INFO.toml"); force = true)
    end
    input_copy = joinpath(out_dir, basename(input_path))
    abspath(input_path) == abspath(input_copy) || cp(input_path, input_copy; force = true)

    # --- A stop request of an earlier segment would end this one at once ---
    stale_stop = joinpath(out_dir, _STOP_FILE)
    if isfile(stale_stop)
        rm(stale_stop)
        @info "Removed a stop request left in $(out_dir) by an earlier segment"
    end

    # --- Build launch script ---
    stdout_path = joinpath(out_dir, cfg.postprocess.stdout_file)
    stderr_path = joinpath(out_dir, "err1000")
    # A restart appends to the capture: this segment's stdout starts here.
    stdout_offset = is_restart && isfile(stdout_path) ? filesize(stdout_path) : 0
    launch_script = _write_launch_script(
        out_dir,
        local_binary,
        input_path,
        stdout_path,
        stderr_path,
        cfg;
        append = is_restart,
    )

    # --- Execute, teeing pipeline logs to the run directory ---
    _with_run_log(run_dir) do
        omp_threads = _effective_omp_threads(sim)
        threads_total = omp_threads * sim.mpi_ranks
        if threads_total > Sys.CPU_THREADS
            @warn "CPU oversubscription: $(sim.mpi_ranks) rank(s) × $omp_threads OpenMP threads " *
                  "exceed the $(Sys.CPU_THREADS) logical CPUs of this host"
        end
        @info "Backend threads: $omp_threads OpenMP × $(sim.mpi_ranks) MPI rank(s)"

        # Segment index: 1 for the initial launch, one more per restart.
        segment = 1 + _segment_count(joinpath(run_dir, "RUN_INFO.toml"))
        # N-body time the segment starts from: the dump's exact time when
        # the caller knows it, else the rounded one of the file name.
        t_start_nb = if !is_restart
            0.0
        elseif haskey(restart, :t_start)
            Float64(restart.t_start)
        else
            _dump_time(restart.dump)
        end

        cpu_before = _children_cpu_times()
        t_start = time()
        launch_date = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")
        @info "Starting $label..."
        # Own session (detach): the hangup of a terminal closing above the
        # pipeline never reaches the engine. An operator interrupt is turned
        # into `_terminate` below, since the terminal's Ctrl-C no longer is.
        process = cd(out_dir) do
            run(detach(`bash $launch_script`); wait = false)
        end

        # The record of the segment exists while the engine runs, so a
        # pipeline killed before the engine exits leaves it on file.
        opening = Dict{String,Any}(
            "index" => segment,
            "kind" => is_restart ? "restart" : "initial",
            "status" => "running",
            "date" => launch_date,
            "input" => basename(input_copy),
            "dump" => is_restart ? restart.dump : "",
            "tcrit_extra" => is_restart ? restart.tcrit_extra : 0.0,
            "t_start" => t_start_nb,
            "host" => gethostname(),
            "pid" => Int(getpid(process)),
            "package_commit" => provenance["package_commit"],
            "stdout_offset" => stdout_offset,
            "stop_requested" => "",
        )
        isempty(sim.gpu_list) || (opening["gpu_list"] = copy(sim.gpu_list))
        slurm_job_id = get(ENV, "SLURM_JOB_ID", "")
        isempty(slurm_job_id) || (opening["slurm_job_id"] = slurm_job_id)
        _open_segment(run_dir, run_id, opening; input_file = basename(input_copy))
        # Known to the exit hook until reaped, so a terminated driver stops it.
        engine_pid = Int32(getpid(process))
        _register_engine(
            _ActiveEngine(
                engine_pid,
                String(run_dir),
                String(out_dir),
                String(stdout_path),
                stdout_offset,
                segment,
                t_start,
                sim.stop_margin,
            ),
        )

        # The helper tasks below hand an operator interrupt to this task, the
        # one waiting on the engine: SIGINT lands in whichever task is current.
        waiting = current_task()

        # The launch script execs the binary, so its PID is the process PID;
        # an mpirun launcher is handled through the process-tree scan.
        monitor = if sim.telemetry_interval > 0
            _start_telemetry(
                run_dir,
                getpid(process),
                t_start;
                interval = sim.telemetry_interval,
                gpu_probe = cfg.build.enable_gpu,
                csv_name = segment == 1 ? "telemetry.csv" : "telemetry_$(segment).csv",
                interrupt_to = waiting,
            )
        else
            nothing
        end

        # Start-up watchdog: terminate a run that never advances past t = 0.
        watchdog =
            sim.startup_timeout > 0 ?
            _start_startup_watchdog(
                stdout_path,
                process,
                sim.startup_timeout;
                interrupt_to = waiting,
                offset = stdout_offset,
            ) : nothing
        # Completion monitor: an engine that printed END RUN but never exits
        # is terminated after the grace period and recorded as completed.
        completion =
            sim.exit_grace > 0 ?
            _start_completion_monitor(
                stdout_path,
                process,
                sim.exit_grace;
                interrupt_to = waiting,
                offset = stdout_offset,
            ) : nothing
        # Wall budget: a stop request `stop_margin` seconds before it expires,
        # termination at the budget.
        stopper =
            sim.wall_budget > 0 ?
            _start_stop_timer(
                out_dir,
                process,
                sim.wall_budget,
                sim.stop_margin;
                interrupt_to = waiting,
            ) : nothing
        # Checkpoint retention: only the newest periodic dumps are kept.
        pruner =
            sim.checkpoint_keep > 0 ?
            _start_dump_pruner(
                out_dir,
                process,
                sim.checkpoint_keep;
                protect = is_restart ? String[restart.dump] : String[],
                interrupt_to = waiting,
            ) : nothing

        # Live ticker is opt-in and interactive-only; the Fortran stdout is
        # captured to out1000 regardless.
        try
            if sim.monitor && stderr isa Base.TTY
                _monitor_stdout_file(
                    stdout_path,
                    process,
                    t_start;
                    live = sim.live_diagnostics,
                    live_interval = sim.live_interval,
                )
            end
            wait(process)
        catch e
            e isa InterruptException || rethrow()
            # Kill first, report after: the console may be gone with whoever
            # sent the interrupt, and a report that fails must not leave the
            # engine running.
            _terminate(process)
            _close_segment(
                run_dir,
                segment;
                status = "killed",
                exit_status = _exit_status(process),
                elapsed = time() - t_start,
            )
            watchdog === nothing || (watchdog.stop[] = true)
            completion === nothing || (completion.stop[] = true)
            stopper === nothing || (stopper.stop[] = true)
            pruner === nothing || (pruner.stop[] = true)
            monitor === nothing || (monitor.stop_requested = true)
            @warn "Interrupted; the engine was terminated (SIGTERM, SIGKILL after $(_KILL_GRACE_SECONDS) s if needed)"
            rethrow()
        finally
            _deregister_engine(engine_pid)
        end
        watchdog === nothing || (watchdog.stop[] = true)
        completion === nothing || (completion.stop[] = true)
        stopper === nothing || (stopper.stop[] = true)
        pruner === nothing || (pruner.stop[] = true)
        pruner === nothing || wait(pruner.task)

        elapsed = time() - t_start
        telemetry =
            _finish_telemetry(monitor, cpu_before, _children_cpu_times(), elapsed, threads_total)
        merge!(telemetry, _backend_performance(stdout_path, stderr_path))

        hung = watchdog !== nothing && watchdog.fired[]
        completed = _run_completed(stdout_path; offset = stdout_offset)
        killed_after_completion = completion !== nothing && completion.fired[]
        status =
            killed_after_completion ? "completed" :
            _segment_status(
                stdout_path;
                offset = stdout_offset,
                exit_status = _exit_status(process),
                watchdog = hung,
            )
        stop_dump = _stop_dump(stdout_path; offset = stdout_offset)
        t_end = if status == "stopped" && !isempty(stop_dump)
            last(_dump_markers(stdout_path; offset = stdout_offset)).time_nb
        else
            t_adjust = _last_adjust_time(stdout_path; offset = stdout_offset)
            isnan(t_adjust) ? t_start_nb : t_adjust
        end
        if killed_after_completion
            @warn "The engine printed END RUN but had not exited $(sim.exit_grace) s later; " *
                  "terminated and recorded as completed (exit status $(_exit_status(process)))"
        elseif !success(process) && !hung
            @warn "Simulation exited with non-zero status ($(_exit_status(process))) after $(_format_elapsed(elapsed))"
        end

        # --- Write run summary (before raising, so a watchdog kill is on record) ---
        _write_run_summary(
            cfg,
            run_dir,
            run_id,
            stdout_path,
            out_dir,
            elapsed;
            telemetry = telemetry,
            input_file = basename(input_copy),
            src_dir = src_dir,
            provenance = provenance,
            segment = merge(
                opening,
                Dict{String,Any}(
                    "status" => status,
                    "elapsed_seconds" => round(elapsed; digits = 1),
                    "exit_status" => _exit_status(process),
                    "completed" => status == "completed",
                    "watchdog" => hung,
                    "terminated_after_completion" => killed_after_completion,
                    "stop_dump" => stop_dump,
                    "t_end" => t_end,
                    "stop_requested" =>
                        (stopper !== nothing && stopper.requested[]) ? "wall_budget" : "",
                    "dumps_pruned" => pruner === nothing ? 0 : pruner.removed[],
                ),
            ),
        )
        if pruner !== nothing && pruner.removed[] > 0
            @info "Checkpoint retention: $(pruner.removed[]) periodic dump(s) removed " *
                  "($(Base.format_bytes(pruner.bytes[]))); the newest $(sim.checkpoint_keep) kept"
        end
        if hung
            error(
                "Simulation terminated by the start-up watchdog: no adjustment beyond t = 0 within " *
                "$(sim.startup_timeout) s. The engine hung after initialisation (recorded in " *
                "RUN_INFO.toml); check the interval values of the input file (DTADJ, DELTAT, DTPLOT " *
                "must have a short exact decimal expansion) or raise simulation.startup_timeout for large N",
            )
        end
        if haskey(telemetry, "cpu_efficiency")
            @info @sprintf(
                "CPU: %.1f s user + %.1f s system on %d thread(s); efficiency %.2f",
                telemetry["cpu_user_s"],
                telemetry["cpu_system_s"],
                threads_total,
                telemetry["cpu_efficiency"]
            )
        end

        if completed
            @info "Simulation complete. Run: $run_id  ($(_format_elapsed(elapsed)))"
        elseif status == "stopped"
            if isempty(stop_dump)
                @warn "Simulation stopped on a stop request without writing a dump; a resume " *
                      "starts from the last periodic dump. Run: $run_id  " *
                      "($(_format_elapsed(elapsed)))"
            else
                @info "Simulation stopped at t = $(t_end) N-body units on a stop request " *
                      "(dump $(stop_dump)); the run can be resumed. Run: $run_id  " *
                      "($(_format_elapsed(elapsed)))"
            end
        elseif status == "halted"
            @warn "Simulation halted by the engine's energy check at t = $(t_end) N-body units " *
                  "(exit status $(_exit_status(process))). Run: $run_id  " *
                  "($(_format_elapsed(elapsed))); the output is partial"
        elseif status == "killed" || status == "failed"
            @warn "Simulation ended without END RUN (status $(status), exit status " *
                  "$(_exit_status(process))). Run: $run_id  ($(_format_elapsed(elapsed))); " *
                  "the output is partial"
        end
    end
    return run_dir
end

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

"""`s` as one POSIX-shell word: single-quoted, embedded single quotes escaped."""
_sh_quote(s::AbstractString) = "'" * replace(String(s), "'" => "'\\''") * "'"

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
    cfg::Nbody6Config;
    append::Bool = false,
)::String
    redir_out = append ? ">>" : ">"
    redir_err = append ? "2>>" : "2>"
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
        # CUDA devices the engine may use (its own GPU_LIST variable; unset = all).
        isempty(sim.gpu_list) ||
            println(io, "export GPU_LIST=$(_sh_quote(join(sim.gpu_list, ' ')))")

        # CUDA environment
        if build.enable_gpu
            cuda = isempty(build.cuda_path) ? detect_cuda_path() : build.cuda_path
            if !isempty(cuda)
                println(io, "export CUDA_HOME=$(_sh_quote(cuda))")
                println(io, "export PATH=$(_sh_quote(joinpath(cuda, "bin"))):\"\$PATH\"")
                println(
                    io,
                    "export LD_LIBRARY_PATH=$(_sh_quote(joinpath(cuda, "lib64"))):\"\$LD_LIBRARY_PATH\"",
                )
            end
        end

        # The simulation command (use stdbuf for line-buffered output if available)
        stdbuf_prefix = "command -v stdbuf >/dev/null 2>&1 && STDBUF='stdbuf -oL' || STDBUF=''"
        println(io, stdbuf_prefix)

        engine_call =
            "$(_sh_quote(binary)) < $(_sh_quote(input)) " *
            "$redir_out $(_sh_quote(stdout_file)) $redir_err $(_sh_quote(stderr_file))"
        if build.enable_mpi && sim.mpi_ranks > 1
            println(io, "exec mpirun --bind-to none -np $(sim.mpi_ranks) \$STDBUF " * engine_call)
        else
            println(io, "exec \$STDBUF " * engine_call)
        end
    end
    chmod(path, 0o755)
    return path
end

"""
    _input_sizes(inp_path) -> Union{Nothing,@NamedTuple{n::Int,nnbopt::Int,nbin0::Int}}

Particle number `N`, neighbour number `NNBOPT` and primordial-binary number
`NBIN0` of an engine input file; `NBIN0` is 0 when the file does not set
it. `nothing` when `N` or `NNBOPT` is not found. Comment lines (first
non-blank character `!`) are not read: the headers of the shipped inputs
describe the run in words such as `N = 100 000`.
"""
function _input_sizes(
    inp_path::AbstractString,
)::Union{Nothing,@NamedTuple{n::Int,nnbopt::Int,nbin0::Int}}
    text = join((l for l in eachline(inp_path) if !startswith(lstrip(l), '!')), '\n')
    mn = match(r"\bN\s*=\s*(\d+)", text)
    mo = match(r"\bNNBOPT\s*=\s*(\d+)", text)
    (mn === nothing || mo === nothing) && return nothing
    mb = match(r"\bNBIN0\s*=\s*(\d+)", text)
    return (
        n = parse(Int, mn.captures[1]),
        nnbopt = parse(Int, mo.captures[1]),
        nbin0 = mb === nothing ? 0 : parse(Int, mb.captures[1]),
    )
end

"""
    _check_engine_limits(inp_path, limits)

Refuse an input the engine build cannot hold, with the three conditions the
engine tests at start-up (`data.F`, `input.F`, `verify.f`): `N + NBIN0 <
NMAX − 2`, `NBIN0 < KMAX − 2` and `NNBOPT ≤ min(N/2, LMAX − 50)`. Throws an
`ArgumentError` naming the quantity, the limit and the remedy; returns
`nothing` when the input fits, or when `limits` is `nothing` or the sizes
cannot be read from the input.
"""
function _check_engine_limits(inp_path::AbstractString, limits::Union{Nothing,NamedTuple})
    limits === nothing && return nothing
    sizes = _input_sizes(inp_path)
    sizes === nothing && return nothing
    name = basename(inp_path)
    rebuild = "configure a larger --with-par preset in build.configure_flags and rebuild the engine"
    if sizes.n + sizes.nbin0 ≥ limits.nmax - 2
        throw(
            ArgumentError(
                "N + NBIN0 = $(sizes.n + sizes.nbin0) of $name does not fit the engine build: " *
                "NMAX = $(limits.nmax), and the engine requires N + NBIN0 < NMAX − 2; $rebuild",
            ),
        )
    end
    if sizes.nbin0 ≥ limits.kmax - 2
        throw(
            ArgumentError(
                "NBIN0 = $(sizes.nbin0) of $name does not fit the engine build: " *
                "KMAX = $(limits.kmax), and the engine requires NBIN0 < KMAX − 2; $rebuild",
            ),
        )
    end
    nnbmax = min(sizes.n ÷ 2, limits.lmax - 50)
    if sizes.nnbopt > nnbmax
        throw(
            ArgumentError(
                "NNBOPT = $(sizes.nnbopt) of $name exceeds the engine's limit min(N/2, LMAX − 50) = " *
                "$nnbmax (N = $(sizes.n), LMAX = $(limits.lmax)); lower NNBOPT (merger.nbody6.nnbopt for a generated input)",
            ),
        )
    end
    return nothing
end

"""
    _adjust_advanced(stdout_path; offset = 0) -> Bool

Whether the captured stdout carries an `ADJUST:` line with `TIME > 0`, the
sign that the integration has advanced past initialisation. Only the lines
from byte `offset` on are read; `false` when the file is shorter.
"""
function _adjust_advanced(stdout_path::AbstractString; offset::Integer = 0)::Bool
    isfile(stdout_path) || return false
    filesize(stdout_path) < offset && return false
    return open(stdout_path) do io
        seek(io, offset)
        for line in eachline(io)
            m = match(r"^\s*ADJUST:\s+TIME\s+([0-9.E+-]+)", line)
            m === nothing && continue
            t = tryparse(Float64, m.captures[1])
            t !== nothing && t > 0 && return true
        end
        return false
    end
end

"""Seconds a terminated engine gets to exit after SIGTERM before SIGKILL."""
const _KILL_GRACE_SECONDS = 15.0

"""
    _terminate(process; grace = _KILL_GRACE_SECONDS)

Send SIGTERM to `process`, wait up to `grace` seconds for it to exit, and
send SIGKILL if it is still running.
"""
function _terminate(process::Base.Process; grace::Real = _KILL_GRACE_SECONDS)
    kill(process)
    deadline = time() + grace
    while process_running(process) && time() < deadline
        sleep(0.2)
    end
    if process_running(process)
        # The signal goes before the report: a report may fail (closed console
        # pipe) and the kill must not depend on it.
        kill(process, Base.SIGKILL)
        @warn "Engine still running $(grace) s after SIGTERM; SIGKILL sent"
    end
    return nothing
end

"""
    _forwarding_interrupts(body, target) -> Task

`@async body()`, except that an `InterruptException` raised inside the task
is re-thrown in `target` (the task waiting on the engine) instead of ending
the helper. SIGINT is delivered to whichever task the main thread is running
at that moment — under a long `wait` that is usually one of the helper
tasks polling the engine — while the engine's termination lives in the
waiting task's handler; without the hand-over the helper died and the
engine ran on. `target === nothing` keeps the plain behaviour.
"""
function _forwarding_interrupts(body::Function, target::Union{Nothing,Task})
    return @async begin
        try
            body()
        catch e
            (e isa InterruptException && target !== nothing && !istaskdone(target)) || rethrow()
            schedule(target, e; error = true)
        end
        nothing
    end
end

"""
    _start_startup_watchdog(stdout_path, process, timeout; interrupt_to = nothing,
                            offset = 0) -> (; stop, fired, task)

Asynchronous watchdog: unless the stdout from byte `offset` on (the
segment's own part of an appended capture) shows an adjustment beyond t = 0
within `timeout` seconds, the process is terminated (SIGTERM, then SIGKILL
after `_KILL_GRACE_SECONDS`) and `fired` is set. The kill precedes the log
record of it: on 2026-09-25 an orphaned pipeline's watchdog fired, its
warning failed on the closed console pipe, and the engine it had not yet
signalled ran for thirteen hours. Setting `stop` ends the watchdog quietly;
an `InterruptException` landing in the task goes to `interrupt_to`
([`_forwarding_interrupts`](@ref)).
"""
function _start_startup_watchdog(
    stdout_path::AbstractString,
    process::Base.Process,
    timeout::Real;
    interrupt_to::Union{Nothing,Task} = nothing,
    offset::Integer = 0,
)
    stop = Ref(false)
    fired = Ref(false)
    task = errormonitor(
        _forwarding_interrupts(interrupt_to) do
            deadline = time() + timeout
            while !stop[] && process_running(process)
                _adjust_advanced(stdout_path; offset = offset) && return nothing
                if time() > deadline
                    fired[] = true
                    _terminate(process)
                    @warn "Start-up watchdog: no adjustment beyond t = 0 after $(timeout) s; the engine was terminated"
                    return nothing
                end
                sleep(min(2.0, timeout))
            end
        end,
    )
    return (stop = stop, fired = fired, task = task)
end

"""Name of the file whose presence in the engine's working directory makes it save its state and exit (`intgrt.F`)."""
const _STOP_FILE = "STOP"

"""
    _request_stop(out_dir) -> String

Ask the engine running in `out_dir` to stop: create the `STOP` file it polls
for. Synchronous, so it can be called from a process-exit hook. Returns the
path.
"""
_request_stop(out_dir::AbstractString)::String = touch(joinpath(out_dir, _STOP_FILE))

"""
    _start_stop_timer(out_dir, process, budget, margin; interrupt_to = nothing)
        -> (; stop, requested, fired, task)

Asynchronous wall-budget timer, counted from the call: `budget − margin`
seconds on, the engine running in `out_dir` is asked to stop
([`_request_stop`](@ref)) and `requested` is set; an engine still running
`budget` seconds on is terminated (SIGTERM, then SIGKILL after
`_KILL_GRACE_SECONDS`) and `fired` is set. Each signal precedes its log
record, as in [`_start_startup_watchdog`](@ref). Setting `stop` ends the
timer quietly; an `InterruptException` landing in the task goes to
`interrupt_to` ([`_forwarding_interrupts`](@ref)).
"""
function _start_stop_timer(
    out_dir::AbstractString,
    process::Base.Process,
    budget::Real,
    margin::Real;
    interrupt_to::Union{Nothing,Task} = nothing,
)
    t0 = time()
    stop = Ref(false)
    requested = Ref(false)
    fired = Ref(false)
    task = errormonitor(
        _forwarding_interrupts(interrupt_to) do
            while !stop[] && process_running(process)
                if !requested[] && time() - t0 ≥ budget - margin
                    _request_stop(out_dir)
                    requested[] = true
                    @info "Wall budget: stop requested after $(round(time() - t0; digits = 1)) s " *
                          "(budget $(budget) s, margin $(margin) s)"
                end
                if time() - t0 ≥ budget
                    fired[] = true
                    _terminate(process)
                    @warn "Wall budget of $(budget) s reached with the engine still running; " *
                          "it was terminated"
                    return nothing
                end
                sleep(0.2)
            end
        end,
    )
    return (stop = stop, requested = requested, fired = fired, task = task)
end

"""
    _prune_dumps(out_dir, keep; protect = String[]) -> @NamedTuple{removed::Int,bytes::Int}

Delete all but the newest `keep` periodic restart dumps (`comm.2_<t>`, regular
files) in `out_dir`, ordered by the time in the name ([`_dump_time`](@ref)),
not lexicographically. Names in `protect` are never deleted, nor are the
`comm.1_<t>` dumps or any other file. Returns the number of files removed
and their total size in bytes; `(removed = 0, bytes = 0)` without touching
anything for `keep ≤ 0` or a missing directory.
"""
function _prune_dumps(
    out_dir::AbstractString,
    keep::Integer;
    protect::AbstractVector{<:AbstractString} = String[],
)::@NamedTuple{removed::Int, bytes::Int}
    (keep ≤ 0 || !isdir(out_dir)) && return (removed = 0, bytes = 0)
    candidates = Tuple{Float64,String}[]
    for f in readdir(out_dir)
        startswith(f, "comm.2_") || continue
        t = _dump_time(f)
        isnan(t) && continue
        isfile(joinpath(out_dir, f)) && push!(candidates, (t, f))
    end
    # Newest first; equal times by name.
    sort!(candidates; rev = true)
    removed = 0
    bytes = 0
    for (_, f) in Iterators.drop(candidates, keep)
        f in protect && continue
        path = joinpath(out_dir, f)
        isfile(path) || continue
        nbytes = filesize(path)
        rm(path)
        removed += 1
        bytes += nbytes
    end
    return (removed = removed, bytes = bytes)
end

"""
    _start_dump_pruner(out_dir, process, keep; protect = String[], interval = 10.0,
                       interrupt_to = nothing) -> (; stop, removed, bytes, task)

Asynchronous retention of the periodic restart dumps in `out_dir`: every
`interval` seconds while `process` runs, all but the newest `keep` are
deleted ([`_prune_dumps`](@ref), names in `protect` excepted); one more
pass follows the end of the loop, for the dumps written since the last.
`removed` and `bytes` count the files deleted and their size. Setting
`stop` ends the pruner within half a second; an `InterruptException`
landing in the task goes to `interrupt_to` ([`_forwarding_interrupts`](@ref)).
"""
function _start_dump_pruner(
    out_dir::AbstractString,
    process::Base.Process,
    keep::Integer;
    protect::AbstractVector{<:AbstractString} = String[],
    interval::Real = 10.0,
    interrupt_to::Union{Nothing,Task} = nothing,
)
    stop = Ref(false)
    removed = Ref(0)
    bytes = Ref(0)
    task = errormonitor(_forwarding_interrupts(interrupt_to) do
        while !stop[] && process_running(process)
            r = _prune_dumps(out_dir, keep; protect = protect)
            removed[] += r.removed
            bytes[] += r.bytes
            t_next = time() + interval
            while !stop[] && process_running(process) && time() < t_next
                sleep(min(0.5, interval))
            end
        end
        r = _prune_dumps(out_dir, keep; protect = protect)
        removed[] += r.removed
        bytes[] += r.bytes
        return nothing
    end)
    return (stop = stop, removed = removed, bytes = bytes, task = task)
end

"""
    _run_completed(stdout_path; offset = 0) -> Bool

Whether the captured stdout carries the engine's `END RUN` line
(`adjust.F`, printed when the termination criterion is met). Only the last
64 KiB are scanned, so the check stays cheap on long runs. Nothing before
byte `offset` is read, so a restart does not see the `END RUN` of an
earlier segment.
"""
function _run_completed(stdout_path::AbstractString; offset::Integer = 0)::Bool
    return occursin("END RUN", _stdout_tail(stdout_path; offset = offset))
end

"""
    _stdout_tail(stdout_path; offset = 0, window = 65536) -> String

Text of the stdout capture from byte `max(offset, filesize - window)` to the
end; `""` when the file is absent or `offset` is not below its size. The
lower bound keeps every check of a segment off what earlier segments of an
appended capture wrote.
"""
function _stdout_tail(
    stdout_path::AbstractString;
    offset::Integer = 0,
    window::Integer = 65536,
)::String
    isfile(stdout_path) || return ""
    size = filesize(stdout_path)
    offset ≥ size && return ""
    return open(stdout_path) do io
        seek(io, max(offset, size - window, 0))
        read(io, String)
    end
end

"""Element type of [`_dump_markers`](@ref)."""
const _DumpMarker = @NamedTuple{file::String, time_nb::Float64, line_end::Int}

"""
    _dump_markers(stdout_path; offset = 0) -> Vector{@NamedTuple{file, time_nb, line_end}}

One entry per restart dump the engine reports from byte `offset` on, in file
order: the dump's file name (`comm.1_<t>` / `comm.2_<t>`), its exact time
`TTOT` from the ` MYDUMP` line (the file name carries a rounded one), and the
byte position just after that line. The ` W MYDUMP` and ` R MYDUMP` lines do
not match. The whole file from `offset` on is read, not a tail window. Empty
when the file is absent.
"""
function _dump_markers(stdout_path::AbstractString; offset::Integer = 0)::Vector{_DumpMarker}
    markers = _DumpMarker[]
    isfile(stdout_path) || return markers
    filesize(stdout_path) ≤ offset && return markers
    open(stdout_path) do io
        seek(io, offset)
        while !eof(io)
            line = readline(io)
            m = match(r"^\s*MYDUMP\s+(\S+)\s+\S+\s+\d+\s+(comm\.[12]_\S+)\s*$", line)
            m === nothing && continue
            t = tryparse(Float64, m.captures[1])
            t === nothing && continue
            push!(
                markers,
                (file = String(m.captures[2]), time_nb = t, line_end = Int(position(io))),
            )
        end
    end
    return markers
end

"""
    _last_adjust_time(stdout_path; offset = 0) -> Float64

`TIME` of the last `ADJUST:` line from byte `offset` on, in N-body units;
`NaN` when there is none.
"""
function _last_adjust_time(stdout_path::AbstractString; offset::Integer = 0)::Float64
    isfile(stdout_path) || return NaN
    filesize(stdout_path) < offset && return NaN
    return open(stdout_path) do io
        seek(io, offset)
        t_last = NaN
        for line in eachline(io)
            m = match(r"^\s*ADJUST:\s+TIME\s+([0-9.E+-]+)", line)
            m === nothing && continue
            t = tryparse(Float64, m.captures[1])
            t === nothing || (t_last = t)
        end
        return t_last
    end
end

"""
    _segment_status(stdout_path; offset = 0, exit_status = nothing, watchdog = false)
        -> String

Status of a finished segment from the stdout it wrote (from byte `offset`
on), its exit status ([`_exit_status`](@ref)) and whether the start-up
watchdog fired. A segment record carries one of seven values:

  - `running`: the engine has been launched and has not exited;
  - `completed`: `END RUN` printed;
  - `halted`: the engine's energy check stopped the run;
  - `stopped`: the engine ended at a stop request and can be resumed;
  - `watchdog`: terminated by the start-up watchdog;
  - `killed`: ended by a signal, or its end was not witnessed;
  - `failed`: exited by itself without any of the above (runtime error).

This function returns one of the last six; the first matching of
`watchdog`, `completed`, `halted`, `stopped`, `killed`, `failed` wins.
"""
function _segment_status(
    stdout_path::AbstractString;
    offset::Integer = 0,
    exit_status::Union{Nothing,Integer} = nothing,
    watchdog::Bool = false,
)::String
    watchdog && return "watchdog"
    tail = _stdout_tail(stdout_path; offset = offset)
    occursin("END RUN", tail) && return "completed"
    occursin("CALCULATIONS HALTED", tail) && return "halted"
    if occursin("TERMINATION BY MANUAL INTERVENTION", tail) || occursin("COMMON SAVED AT", tail)
        return "stopped"
    end
    (exit_status === nothing || exit_status < 0) && return "killed"
    return "failed"
end

"""
    _stop_dump(stdout_path; offset = 0) -> String

File name of the `comm.1_<t>` dump the engine wrote at a stop request: the
last dump line from byte `offset` on that follows the last
`TERMINATION BY MANUAL INTERVENTION` line. `""` when there is no stop line,
or no such dump line after it (some inputs stop without writing a dump).
"""
function _stop_dump(stdout_path::AbstractString; offset::Integer = 0)::String
    isfile(stdout_path) || return ""
    filesize(stdout_path) ≤ offset && return ""
    return open(stdout_path) do io
        seek(io, offset)
        stopped = false
        name = ""
        for line in eachline(io)
            if occursin("TERMINATION BY MANUAL INTERVENTION", line)
                stopped = true
                name = ""
            elseif stopped
                m = match(r"^\s*MYDUMP\s+(\S+)\s+\S+\s+\d+\s+(comm\.[12]_\S+)\s*$", line)
                m === nothing && continue
                tryparse(Float64, m.captures[1]) === nothing && continue
                startswith(m.captures[2], "comm.1_") && (name = String(m.captures[2]))
            end
        end
        return name
    end
end

"""
    _open_segment(run_dir, run_id, segment; input_file = "")

Append the record of a segment that is starting to the `segments` list of
`run_dir/RUN_INFO.toml` (created when absent) and set `run.id`,
`run.segments` and `run.status = "running"`; `run.input_file` is set only
when not yet recorded, so a restart keeps the original input. Every other
key of the file is kept. Written before the engine's exit so that a run
whose pipeline dies leaves a record of the segment.
"""
function _open_segment(
    run_dir::AbstractString,
    run_id::AbstractString,
    segment::AbstractDict;
    input_file::AbstractString = "",
)
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    d = isfile(info_path) ? TOML.parsefile(info_path) : Dict{String,Any}()
    segments = Vector{Any}(get(d, "segments", Any[]))
    push!(segments, Dict{String,Any}(segment))
    d["segments"] = segments
    run_table = get!(d, "run", Dict{String,Any}())
    run_table["id"] = String(run_id)
    run_table["segments"] = length(segments)
    run_table["status"] = "running"
    if isempty(get(run_table, "input_file", "")) && !isempty(input_file)
        run_table["input_file"] = String(input_file)
    end
    _atomic_write_toml(info_path, d)
    return nothing
end

"""
    _close_segment(run_dir, index; status, exit_status = nothing, elapsed = nothing,
                   t_end = nothing, fields = Dict{String,Any}()) -> Bool

Close the open record of segment `index` in `run_dir/RUN_INFO.toml`: set its
`status` and `completed`, and `exit_status`, `elapsed_seconds` and `t_end`
when given, then the entries of `fields`; set `run.status`. No other key is
touched. Returns `false` without writing when the file or the entry is
missing. Uses only TOML parsing and file writes, since it also runs from a
process-exit hook.
"""
function _close_segment(
    run_dir::AbstractString,
    index::Integer;
    status::AbstractString,
    exit_status::Union{Nothing,Integer} = nothing,
    elapsed::Union{Nothing,Real} = nothing,
    t_end::Union{Nothing,Real} = nothing,
    fields::AbstractDict = Dict{String,Any}(),
)::Bool
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    isfile(info_path) || return false
    d = TOML.parsefile(info_path)
    segments = get(d, "segments", Any[])
    k = findfirst(s -> s isa AbstractDict && get(s, "index", nothing) == index, segments)
    k === nothing && return false
    entry = segments[k]
    entry["status"] = String(status)
    entry["completed"] = status == "completed"
    exit_status === nothing || (entry["exit_status"] = Int(exit_status))
    elapsed === nothing || (entry["elapsed_seconds"] = round(Float64(elapsed); digits = 1))
    t_end === nothing || (entry["t_end"] = Float64(t_end))
    for (key, value) in fields
        entry[String(key)] = value
    end
    run_table = get!(d, "run", Dict{String,Any}())
    run_table["status"] = String(status)
    _atomic_write_toml(info_path, d)
    return true
end

"""
An engine process of this session, as the exit hook needs it: where it runs,
which segment record is open for it, and how long it may take to stop.
"""
struct _ActiveEngine
    pid::Int32
    run_dir::String
    out_dir::String
    stdout_path::String
    stdout_offset::Int
    segment::Int
    t_launch::Float64
    stop_margin::Float64
end

"""
Engines launched by this process and not yet reaped. The vector is replaced,
never mutated: the exit hook may interrupt any task, and must not find a
half-updated container.
"""
const _ACTIVE_ENGINES = Ref{Vector{_ActiveEngine}}(_ActiveEngine[])

"""Add an engine to [`_ACTIVE_ENGINES`](@ref) (by replacing the vector)."""
_register_engine(e::_ActiveEngine) = (_ACTIVE_ENGINES[] = vcat(_ACTIVE_ENGINES[], [e]); nothing)

"""Remove the engine `pid` from [`_ACTIVE_ENGINES`](@ref) (by replacing the vector)."""
_deregister_engine(pid::Integer) =
    (_ACTIVE_ENGINES[] = filter(e -> e.pid != pid, _ACTIVE_ENGINES[]); nothing)

"""
    _child_state(pid) -> Tuple{Bool,Union{Nothing,Int}}

State of the child process `pid` by a non-blocking `waitpid`: `(true,
nothing)` while it runs; `(false, code)` once it has exited, with `code` its
exit code or the negated number of the signal that ended it; `(false,
nothing)` when `pid` is not (or no longer) a child of this process. A call
that finds the child exited reaps it. For the exit hook only
([`_stop_engines_at_exit`](@ref)); elsewhere the event loop reaps children.
"""
function _child_state(pid::Integer)::Tuple{Bool,Union{Nothing,Int}}
    status = Ref{Cint}(0)
    r = ccall(:waitpid, Cint, (Cint, Ptr{Cint}, Cint), pid, status, 1)   # 1 = WNOHANG
    r == 0 && return (true, nothing)
    if r == pid
        sig = status[] & 0x7f
        code = sig == 0 ? Int((status[] >> 8) & 0xff) : -Int(sig)
        return (false, code)
    end
    return (false, nothing)
end

"""
    _stop_engines_at_exit()

Process-exit hook: ask every engine of [`_ACTIVE_ENGINES`](@ref) still
running to stop ([`_request_stop`](@ref)), wait up to the largest
`stop_margin` for them to exit, terminate the rest (SIGTERM, SIGKILL after
`_KILL_GRACE_SECONDS`), and close their segment records with the status the
stdout shows and `stop_requested = "signal"`. A batch scheduler ends a job
with SIGTERM; the engine runs in its own session, so without the hook it
would outlive the driver, or die with it without writing a dump. Julia runs
the hook outside the event loop: only synchronous operations are used (file
writes, `waitpid` and `kill` by `ccall`, `Libc.systemsleep`), no `sleep`,
`wait`, `run`, tasks or logging macros.
"""
function _stop_engines_at_exit()
    engines = _ACTIVE_ENGINES[]
    isempty(engines) && return nothing
    t0 = time()
    alive = Dict{Int32,Bool}()
    codes = Dict{Int32,Union{Nothing,Int}}()
    for e in engines
        running, code = _child_state(e.pid)
        alive[e.pid] = running
        codes[e.pid] = code
        running && _request_stop(e.out_dir)
    end
    # Poll the engines still marked alive until none is or `deadline` passes.
    function poll(deadline::Float64)
        while any(e -> alive[e.pid], engines) && time() ≤ deadline
            Libc.systemsleep(0.2)
            for e in engines
                alive[e.pid] || continue
                running, code = _child_state(e.pid)
                running && continue
                alive[e.pid] = false
                codes[e.pid] = code
            end
        end
        return nothing
    end
    poll(t0 + maximum(e.stop_margin for e in engines))
    if any(e -> alive[e.pid], engines)
        for e in engines
            alive[e.pid] && ccall(:kill, Cint, (Cint, Cint), e.pid, 15)
        end
        poll(time() + _KILL_GRACE_SECONDS)
        for e in engines
            alive[e.pid] && ccall(:kill, Cint, (Cint, Cint), e.pid, 9)
        end
        Libc.systemsleep(0.2)
        for e in engines
            alive[e.pid] || continue
            running, code = _child_state(e.pid)
            running && continue
            alive[e.pid] = false
            codes[e.pid] = code
        end
    end
    for e in engines
        # One record that cannot be closed must not keep the others open.
        try
            status =
                _segment_status(e.stdout_path; offset = e.stdout_offset, exit_status = codes[e.pid])
            stop_dump = _stop_dump(e.stdout_path; offset = e.stdout_offset)
            markers = _dump_markers(e.stdout_path; offset = e.stdout_offset)
            t_end = if status == "stopped" && !isempty(stop_dump)
                last(markers).time_nb
            else
                t_adjust = _last_adjust_time(e.stdout_path; offset = e.stdout_offset)
                isnan(t_adjust) ? nothing : t_adjust
            end
            _close_segment(
                e.run_dir,
                e.segment;
                status = status,
                exit_status = codes[e.pid],
                elapsed = time() - e.t_launch,
                t_end = t_end,
                fields = Dict{String,Any}("stop_requested" => "signal", "stop_dump" => stop_dump),
            )
        catch err
            err isa Union{SystemError,Base.IOError,TOML.ParserError,ArgumentError} || rethrow()
        end
    end
    _ACTIVE_ENGINES[] = _ActiveEngine[]
    return nothing
end

"""
    _latest_mtime(dir) -> Float64

Most recent modification time (Unix seconds) of the files directly in
`dir`; `0.0` for an absent or empty directory.
"""
function _latest_mtime(dir::AbstractString)::Float64
    isdir(dir) || return 0.0
    latest = 0.0
    for f in readdir(dir; join = true)
        isfile(f) && (latest = max(latest, mtime(f)))
    end
    return latest
end

"""
    _start_completion_monitor(stdout_path, process, grace; interrupt_to = nothing,
                              offset = 0) -> (; stop, fired, task)

Asynchronous monitor: once the stdout from byte `offset` on (the segment's
own part of an appended capture) shows `END RUN`, the process is
terminated (SIGTERM, then SIGKILL after `_KILL_GRACE_SECONDS`) and `fired`
is set as soon as no file in the output
directory (the one holding `stdout_path`) has been modified for `grace`
seconds while the process is still alive. The engine writes its final
COMMON dump after `END RUN` when `KZ(1) > 0`; that write keeps a file
changing and therefore keeps the engine alive, whatever its duration. What
the monitor ends is an engine that has finished writing and does not exit,
as observed in tidal-field runs, which otherwise blocks the pipeline until
an external timeout. Setting `stop` ends the monitor quietly; an
`InterruptException` landing in the task goes to `interrupt_to`
([`_forwarding_interrupts`](@ref)). The kill precedes its log record, as in
[`_start_startup_watchdog`](@ref).
"""
function _start_completion_monitor(
    stdout_path::AbstractString,
    process::Base.Process,
    grace::Real;
    interrupt_to::Union{Nothing,Task} = nothing,
    offset::Integer = 0,
)
    out_dir = dirname(abspath(stdout_path))
    stop = Ref(false)
    fired = Ref(false)
    task = errormonitor(
        _forwarding_interrupts(interrupt_to) do
            completed = false
            while !stop[] && process_running(process)
                completed || (completed = _run_completed(stdout_path; offset = offset))
                if completed && time() - _latest_mtime(out_dir) ≥ grace
                    fired[] = true
                    _terminate(process)
                    @warn "Completion monitor: END RUN printed and the output directory idle for $(grace) s, " *
                          "but the engine had not exited; it was terminated"
                    return nothing
                end
                sleep(min(5.0, max(grace, 1.0)))
            end
        end,
    )
    return (stop = stop, fired = fired, task = task)
end

"""
Poll the stdout file while the simulation runs and print ADJUST summaries
with wall-clock elapsed time and a heartbeat spinner.

The spinner updates every 2 seconds on the same line, showing the user that
the process is alive even during quiet periods. When a diagnostics line
(ADJUST or TIME[NB]) is detected, it prints as a full log line and resets
the spinner.
"""
function _monitor_stdout_file(
    path::String,
    process::Base.Process,
    t_start::Float64 = time();
    live::Bool = false,
    live_interval::Real = 30.0,
)
    last_pos = 0
    spin_idx = 1
    last_event_t = t_start   # wall-clock of last printed diagnostics line
    last_live = t_start      # wall-clock of the last sparkline panel

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

        # Opt-in in-terminal sparklines of the diagnostics so far: printed
        # as a block below the log lines; the spinner resumes underneath.
        if live && time() - last_live ≥ live_interval
            panel = _live_diagnostics_panel(path)
            if panel !== nothing
                print(stderr, "\r\e[K")
                println(stderr, panel)
            end
            last_live = time()
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
    _live_diagnostics_panel(stdout_path; width = 60, height = 6) -> Union{Nothing,String}

In-terminal sparklines of the run so far, from the ADJUST records of the
stdout capture: the virial ratio `Q = T/|W|` and `log10 |ΔE/E|` against
time (Myr when the scaling is known, N-body units otherwise), drawn with
UnicodePlots at `width × height` characters each and returned as one
string. `nothing` with fewer than two adjustments, or when the file cannot
be parsed yet (a line may be half written).
"""
function _live_diagnostics_panel(
    stdout_path::AbstractString;
    width::Int = 60,
    height::Int = 6,
)::Union{Nothing,String}
    isfile(stdout_path) || return nothing
    diag = try
        read_diagnostics(stdout_path)
    catch e
        e isa Union{ArgumentError,ErrorException,Base.IOError,BoundsError} || rethrow()
        return nothing
    end
    adj = diag.adjust
    length(adj) ≥ 2 || return nothing
    physical = any(r -> r.time_myr > 0, adj)
    t = physical ? [r.time_myr for r in adj] : [r.time_nb for r in adj]
    xlabel = physical ? "t [Myr]" : "t [NB]"
    q = [r.qvir for r in adj]
    plots = String[]
    push!(
        plots,
        sprint(
            show,
            UnicodePlots.lineplot(
                t,
                q;
                xlabel = xlabel,
                ylabel = "Q",
                width = width,
                height = height,
                name = "Q = T/|W|",
            ),
        ),
    )
    nz = [i for i in eachindex(adj) if adj[i].de_rel != 0]
    if length(nz) ≥ 2
        push!(
            plots,
            sprint(
                show,
                UnicodePlots.lineplot(
                    t[nz],
                    log10.(abs.([adj[i].de_rel for i in nz]));
                    xlabel = xlabel,
                    ylabel = "log10|dE/E|",
                    width = width,
                    height = height,
                    name = "energy error",
                ),
            ),
        )
    end
    return join(plots, "\n")
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
        catch e
            # A malformed ADJUST line only costs one ticker update.
            e isa ArgumentError || e isa BoundsError || rethrow()
            @debug "Skipping unparseable ADJUST line in the live monitor" line = stripped exception =
                e
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
    _reported_gpu_devices(stderr_path) -> Vector{String}

Devices the engine's GPU library reports at initialisation
(`gpunb.velocity.cu` prints `# GPU initialization - rank: r; HOST h;
NGPU n; device: i <name>` to stderr, one line per device and rank), as
`"rank r: device i <name>"`, distinct. Empty without the file or the lines
(CPU builds print none).
"""
function _reported_gpu_devices(stderr_path::AbstractString)::Vector{String}
    isfile(stderr_path) || return String[]
    devices = String[]
    for line in eachline(stderr_path)
        m = match(r"# GPU initialization - rank:\s*(\d+);.*device:\s*(\d+)\s+(.*)$", line)
        m === nothing && continue
        entry = "rank $(m.captures[1]): device $(m.captures[2]) $(strip(m.captures[3]))"
        entry in devices || push!(devices, entry)
    end
    return devices
end

"""Exit status of a finished process: its exit code, or the negated signal number when it was terminated by a signal (Julia reports exit code 0 in that case)."""
_exit_status(p::Base.Process)::Int = p.termsignal != 0 ? -Int(p.termsignal) : Int(p.exitcode)

"""Number of entries in the `segments` list of an existing `RUN_INFO.toml` (0 without the file)."""
function _segment_count(info_path::AbstractString)::Int
    isfile(info_path) || return 0
    return length(get(TOML.parsefile(info_path), "segments", Any[]))
end

"""
    _write_run_summary(cfg, run_dir, run_id, stdout_path, out_dir, elapsed;
                       telemetry = nothing, input_file = "", segment = nothing)

Write the machine-readable run summary `RUN_INFO.toml` into `run_dir`:
run identity, wall-clock time, and the backend thread layout (effective
OpenMP threads, MPI ranks, the configured `gpu_list` and the devices the
engine reported); provenance (package and backend commits); the
`[build]` table copied from the binary's `BUILD_INFO.toml` when present;
the hardware fingerprint ([`_hardware_fingerprint`](@ref)); the
`[telemetry]` table from [`_finish_telemetry`](@ref) when given; the
output file inventory; and the `segments` list, one entry per launch
(initial run and restarts). On a restart the previous segments are kept,
`run.elapsed_seconds` accumulates, and the latest telemetry replaces the
table (per-segment CSVs remain). A last entry still `running` with the index
of `segment` (the record opened at launch) is replaced by `segment`, and
`run.status` takes the status of `segment`, or keeps the recorded one.
"""
function _write_run_summary(
    cfg::Nbody6Config,
    run_dir::String,
    run_id::String,
    stdout_path::String,
    out_dir::String = run_dir,
    elapsed::Float64 = 0.0;
    telemetry::Union{Nothing,Dict{String,Any}} = nothing,
    input_file::AbstractString = "",
    segment::Union{Nothing,Dict{String,Any}} = nothing,
    src_dir::AbstractString = "",
    provenance::Union{Nothing,Dict{String,Any}} = nothing,
)
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    previous = isfile(info_path) ? TOML.parsefile(info_path) : Dict{String,Any}()
    segments = Vector{Any}(get(previous, "segments", Any[]))
    if segment !== nothing
        # The record opened at launch (`_open_segment`) is completed in place.
        open_entry = isempty(segments) ? nothing : last(segments)
        if open_entry isa AbstractDict &&
           get(open_entry, "status", "") == "running" &&
           haskey(segment, "index") &&
           get(open_entry, "index", nothing) == segment["index"]
            segments[end] = segment
        else
            push!(segments, segment)
        end
    end
    elapsed_total =
        elapsed + sum(
            Float64(get(s, "elapsed_seconds", 0.0)) for
            s in segments[1:(end - (segment === nothing ? 0 : 1))];
            init = 0.0,
        )
    run_table = Dict{String,Any}(
        "id" => run_id,
        "date" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        "elapsed_seconds" => round(elapsed_total; digits = 1),
        "elapsed" => _format_elapsed(elapsed_total),
        "omp_threads" => _effective_omp_threads(cfg.simulation),
        "mpi_ranks" => cfg.simulation.mpi_ranks,
        "segments" => length(segments),
    )
    previous_run = get(previous, "run", Dict{String,Any}())
    if segment !== nothing && haskey(segment, "status")
        run_table["status"] = segment["status"]
    elseif haskey(previous_run, "status")
        run_table["status"] = previous_run["status"]
    end
    # The original input is recorded once; restarts must not replace it.
    previous_input = get(previous_run, "input_file", "")
    recorded_input = isempty(previous_input) ? String(input_file) : String(previous_input)
    isempty(recorded_input) || (run_table["input_file"] = recorded_input)
    isfile(stdout_path) && (run_table["stdout_lines"] = countlines(stdout_path))
    reported = _reported_omp_threads(stdout_path)
    reported === nothing || (run_table["omp_threads_reported"] = reported)
    isempty(cfg.simulation.gpu_list) || (run_table["gpu_list"] = copy(cfg.simulation.gpu_list))
    devices = _reported_gpu_devices(joinpath(dirname(stdout_path), "err1000"))
    isempty(devices) || (run_table["gpu_devices"] = devices)

    # Commits captured at launch when given (`_execute_simulation`); the
    # fallback reads the trees now, which is right only for a record written
    # in the same session as the launch.
    stamps =
        provenance === nothing ?
        Dict{String,Any}(
            "package_commit" => _source_stamp(_PACKAGE_ROOT),
            "backend_commit" => isempty(src_dir) ? "unknown" : _git_commit(src_dir),
        ) : Dict{String,Any}(provenance)
    d = Dict{String,Any}(
        "run" => run_table,
        "provenance" => stamps,
        "hardware" => _hardware_fingerprint(),
    )
    manifest = _snapshot_manifest(run_dir)
    manifest === nothing || (d["provenance"]["environment_manifest"] = manifest)
    telemetry === nothing || (d["telemetry"] = telemetry)
    isempty(segments) || (d["segments"] = segments)
    build_info = joinpath(out_dir, "BUILD_INFO.toml")
    isfile(build_info) && (d["build"] = TOML.parsefile(build_info))
    if isdir(out_dir)
        files = sort(readdir(out_dir))
        d["output"] = Dict{String,Any}(
            "files" => files,
            "total_bytes" => sum(f -> filesize(joinpath(out_dir, f)), files; init = 0),
        )
    end

    _atomic_write_toml(info_path, d)
    return nothing
end

"""
    _stamp_pipeline_completion(run_dir, phases, elapsed) -> Bool

Record in `run_dir/RUN_INFO.toml` that the pipeline ran to its end, as
`[pipeline]` with `completed = true`, the phases performed and the total
wall time. Returns `false`, without raising, when there is no run directory
or no summary to amend.

A run summary is written as soon as the engine exits, so its presence says
only that the simulation finished — post-processing and plotting come
afterwards and a run killed in between leaves a summary that looks
complete. This marker is what distinguishes a finished pipeline from an
interrupted one, and [`_pipeline_completed`](@ref) is what reads it back.
`engine_completed` states separately whether the last engine segment reached
END RUN, since the pipeline also completes on partial output.
"""
function _stamp_pipeline_completion(
    run_dir::AbstractString,
    phases::AbstractVector{<:AbstractString},
    elapsed::Real,
)::Bool
    isempty(run_dir) && return false
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    isfile(info_path) || return false
    d = TOML.parsefile(info_path)
    d["pipeline"] = Dict{String,Any}(
        "completed" => true,
        "finished" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        "phases" => String.(phases),
        "elapsed_seconds" => round(Float64(elapsed); digits = 1),
    )
    # The pipeline post-processes whatever the engine left, so its own
    # completion says nothing about the integration; state that separately.
    segments = get(d, "segments", Any[])
    isempty(segments) ||
        (d["pipeline"]["engine_completed"] = get(last(segments), "completed", false) === true)
    _atomic_write_toml(info_path, d)
    return true
end

"""
    _pipeline_completed(run_dir) -> Bool

Whether `run_dir/RUN_INFO.toml` carries the `[pipeline] completed` marker of
[`_stamp_pipeline_completion`](@ref). `false` for a missing directory,
a missing summary, an unreadable one, or a run that predates the marker.
"""
function _pipeline_completed(run_dir::AbstractString)::Bool
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    isfile(info_path) || return false
    return try
        pipeline_table = get(TOML.parsefile(info_path), "pipeline", Dict{String,Any}())
        get(pipeline_table, "completed", false) === true
    catch e
        e isa TOML.ParserError || rethrow()
        false
    end
end

"""
    _engine_completed(run_dir) -> Union{Nothing,Bool}

Whether the last engine segment recorded in `run_dir/RUN_INFO.toml` reached
`END RUN`. `nothing` when the summary is missing or unreadable, or lists no
segment (a post-processing-only pipeline).
"""
function _engine_completed(run_dir::AbstractString)::Union{Nothing,Bool}
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    isfile(info_path) || return nothing
    segments = try
        get(TOML.parsefile(info_path), "segments", Any[])
    catch e
        e isa TOML.ParserError || rethrow()
        return nothing
    end
    isempty(segments) && return nothing
    return get(last(segments), "completed", false) === true
end
