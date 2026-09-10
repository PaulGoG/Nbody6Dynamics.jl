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
    run_simulation(cfg::Nbody6Config; base_dir = _PROJECT_ROOT, run_id = "") -> String

Execute the Nbody6++ simulation with:
- Unique run ID and isolated output directory
- Bash wrapper for `ulimit -s unlimited` and CUDA env
- Real-time ADJUST line monitoring with wall-clock timer and heartbeat spinner
- Frozen config snapshot saved alongside output

Returns the path to the run directory.
"""
function run_simulation(
    cfg::Nbody6Config;
    base_dir::AbstractString = _PROJECT_ROOT,
    run_id::AbstractString = "",
)::String
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
    isempty(run_id) && (run_id = generate_run_id(sim.run_id_prefix))
    run_dir = abspath(joinpath(base_dir, sim.runs_dir, run_id))
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
(`comm.1_<t>` / `comm.2_<t>`, written every `NCOMM × DELTAT`), or
`nothing` when none exists.
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
    bi = replace(mi.captures[1], r"TCRIT\s*=\s*[0-9.eE+-]+" => @sprintf("TCRIT=%.4f", tcrit_extra))
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
                       base_dir = _PROJECT_ROOT) -> String

Continue a finished or interrupted run from one of the engine's COMMON
dumps for `tcrit_extra` further N-body time units. The chosen dump
(default: the latest `comm.[12]_<t>` in `output/`) is copied to
`output/comm.1`, which is what `KSTART = 2` reads; a restart input is
written from the run's original input file; the engine runs in the same
output directory with stdout and stderr appended, so `out1000`, `lagr.7`,
`esc.11`, and the time-stamped snapshot and stellar-evolution files
continue. `RUN_INFO.toml` gains one entry in its `segments` list per
launch and the telemetry of each segment goes to its own CSV. Returns
`run_dir`.
"""
function restart_simulation(
    run_dir::AbstractString;
    tcrit_extra::Real,
    dump::Union{Nothing,AbstractString} = nothing,
    tcrtp0::Union{Nothing,Real} = nothing,
    base_dir::AbstractString = _PROJECT_ROOT,
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
                        base_dir = _PROJECT_ROOT, label = "simulation",
                        restart = nothing) -> String

Shared execution core for [`run_simulation`](@ref), the merger pipeline,
and [`restart_simulation`](@ref): locates the binary, freezes the config
into `run_dir`, copies the binary and the input file into `out_dir` for
reproducibility, writes the launch script, and runs it under the teed run
log (§9) with the opt-in live monitor. `input_path` must be absolute (the
launch script executes from `out_dir`). With `restart = (; dump,
tcrit_extra)` the binary copy is reused, stdout and stderr are appended,
and the run summary records a further segment. Returns `run_dir`.
"""
function _execute_simulation(
    cfg::Nbody6Config,
    run_dir::AbstractString,
    out_dir::AbstractString,
    input_path::AbstractString;
    base_dir::AbstractString = _PROJECT_ROOT,
    label::AbstractString = "simulation",
    restart::Union{Nothing,NamedTuple} = nothing,
)::String
    sim = cfg.simulation
    run_id = basename(run_dir)
    is_restart = restart !== nothing

    # --- Locate binary ---
    src_dir = joinpath(base_dir, cfg.install.install_dir)
    binary = _find_binary(src_dir, sim.binary_name, cfg.build)

    # --- Save frozen config ---
    is_restart || save_config(cfg, joinpath(run_dir, "config.toml"))

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

    # --- Build launch script ---
    stdout_path = joinpath(out_dir, cfg.postprocess.stdout_file)
    stderr_path = joinpath(out_dir, "err1000")
    launch_script = _write_launch_script(
        out_dir,
        local_binary,
        input_path,
        stdout_path,
        stderr_path,
        cfg;
        append = is_restart,
    )

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

        # Segment index: 1 for the initial launch, one more per restart.
        segment = 1 + _segment_count(joinpath(run_dir, "RUN_INFO.toml"))

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
            )
        else
            nothing
        end

        # Start-up watchdog: terminate a run that never advances past t = 0.
        watchdog =
            sim.startup_timeout > 0 ?
            _start_startup_watchdog(stdout_path, process, sim.startup_timeout) : nothing
        # Completion monitor: an engine that printed END RUN but never exits
        # is terminated after the grace period and recorded as completed.
        completion =
            sim.exit_grace > 0 ? _start_completion_monitor(stdout_path, process, sim.exit_grace) :
            nothing

        # Live ticker is opt-in and interactive-only; the Fortran stdout is
        # captured to out1000 regardless.
        if sim.monitor && stderr isa Base.TTY
            _monitor_stdout_file(stdout_path, process, t_start)
        end
        wait(process)
        watchdog === nothing || (watchdog.stop[] = true)
        completion === nothing || (completion.stop[] = true)

        elapsed = time() - t_start
        telemetry =
            _finish_telemetry(monitor, cpu_before, _children_cpu_times(), elapsed, threads_total)
        merge!(telemetry, _backend_performance(stdout_path, stderr_path))

        hung = watchdog !== nothing && watchdog.fired[]
        completed = _run_completed(stdout_path)
        killed_after_completion = completion !== nothing && completion.fired[]
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
            segment = Dict{String,Any}(
                "index" => segment,
                "kind" => is_restart ? "restart" : "initial",
                "date" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
                "elapsed_seconds" => round(elapsed; digits = 1),
                "exit_status" => _exit_status(process),
                "completed" => completed,
                "watchdog" => hung,
                "terminated_after_completion" => killed_after_completion,
                "input" => basename(input_copy),
                "dump" => is_restart ? restart.dump : "",
                "tcrit_extra" => is_restart ? restart.tcrit_extra : 0.0,
            ),
        )
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
        isempty(sim.gpu_list) || println(io, "export GPU_LIST=\"$(join(sim.gpu_list, ' '))\"")

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
                "\$STDBUF \"$binary\" < \"$input\" $redir_out \"$stdout_file\" $redir_err \"$stderr_file\"",
            )
        else
            println(
                io,
                "exec \$STDBUF \"$binary\" < \"$input\" $redir_out \"$stdout_file\" $redir_err \"$stderr_file\"",
            )
        end
    end
    chmod(path, 0o755)
    return path
end

"""
    _adjust_advanced(stdout_path) -> Bool

Whether the captured stdout carries an `ADJUST:` line with `TIME > 0`, the
sign that the integration has advanced past initialisation.
"""
function _adjust_advanced(stdout_path::AbstractString)::Bool
    isfile(stdout_path) || return false
    for line in eachline(stdout_path)
        m = match(r"^\s*ADJUST:\s+TIME\s+([0-9.E+-]+)", line)
        m === nothing && continue
        t = tryparse(Float64, m.captures[1])
        t !== nothing && t > 0 && return true
    end
    return false
end

"""
    _start_startup_watchdog(stdout_path, process, timeout) -> (; stop, fired, task)

Asynchronous watchdog: unless the stdout shows an adjustment beyond t = 0
within `timeout` seconds, the process is terminated (SIGTERM) and `fired`
is set. Setting `stop` ends the watchdog quietly.
"""
function _start_startup_watchdog(stdout_path::AbstractString, process::Base.Process, timeout::Real)
    stop = Ref(false)
    fired = Ref(false)
    task = @async begin
        deadline = time() + timeout
        while !stop[] && process_running(process)
            _adjust_advanced(stdout_path) && return nothing
            if time() > deadline
                fired[] = true
                @warn "Start-up watchdog: no adjustment beyond t = 0 after $(timeout) s; terminating the engine"
                kill(process)
                return nothing
            end
            sleep(min(2.0, timeout))
        end
    end
    return (stop = stop, fired = fired, task = task)
end

"""
    _run_completed(stdout_path) -> Bool

Whether the captured stdout carries the engine's `END RUN` line
(`adjust.F`, printed when the termination criterion is met). Only the last
64 KiB are scanned, so the check stays cheap on long runs.
"""
function _run_completed(stdout_path::AbstractString)::Bool
    isfile(stdout_path) || return false
    tail = open(stdout_path) do io
        size = filesize(stdout_path)
        seek(io, max(0, size - 65536))
        read(io, String)
    end
    return occursin("END RUN", tail)
end

"""
    _start_completion_monitor(stdout_path, process, grace) -> (; stop, fired, task)

Asynchronous monitor: once the stdout shows `END RUN`, the process is given
`grace` seconds to exit; if it is still running afterwards it is terminated
(SIGTERM) and `fired` is set. The engine has been observed to finish its
integration, print its final tables, and never exit (tidal-field runs), which
otherwise blocks the pipeline until an external timeout. Setting `stop`
ends the monitor quietly.
"""
function _start_completion_monitor(stdout_path::AbstractString, process::Base.Process, grace::Real)
    stop = Ref(false)
    fired = Ref(false)
    task = @async begin
        deadline = Inf
        while !stop[] && process_running(process)
            if deadline == Inf && _run_completed(stdout_path)
                deadline = time() + grace
            end
            if time() > deadline
                fired[] = true
                @warn "Completion monitor: END RUN printed but the engine did not exit within $(grace) s; terminating it"
                kill(process)
                return nothing
            end
            sleep(min(5.0, max(grace, 1.0)))
        end
    end
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
the hardware fingerprint (§6; GPU probed when `build.enable_gpu`); the
`[telemetry]` table from [`_finish_telemetry`](@ref) when given; the
output file inventory; and the `segments` list, one entry per launch
(initial run and restarts). On a restart the previous segments are kept,
`run.elapsed_seconds` accumulates, and the latest telemetry replaces the
table (per-segment CSVs remain).
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
)
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    previous = isfile(info_path) ? TOML.parsefile(info_path) : Dict{String,Any}()
    segments = Vector{Any}(get(previous, "segments", Any[]))
    segment === nothing || push!(segments, segment)
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
    # The original input is recorded once; restarts must not replace it.
    previous_input = get(get(previous, "run", Dict{String,Any}()), "input_file", "")
    recorded_input = isempty(previous_input) ? String(input_file) : String(previous_input)
    isempty(recorded_input) || (run_table["input_file"] = recorded_input)
    isfile(stdout_path) && (run_table["stdout_lines"] = countlines(stdout_path))
    reported = _reported_omp_threads(stdout_path)
    reported === nothing || (run_table["omp_threads_reported"] = reported)
    isempty(cfg.simulation.gpu_list) || (run_table["gpu_list"] = copy(cfg.simulation.gpu_list))
    devices = _reported_gpu_devices(joinpath(dirname(stdout_path), "err1000"))
    isempty(devices) || (run_table["gpu_devices"] = devices)

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

    open(info_path, "w") do io
        TOML.print(io, d)
    end
    return nothing
end
