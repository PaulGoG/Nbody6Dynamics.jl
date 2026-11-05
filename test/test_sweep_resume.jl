# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Re-entry of a sweep after its driver was ended (roadmap 8.4): the state of
# every point read from its own run record, the relaunch of the unfinished
# points with the finished ones kept, the continuation of an interrupted point
# by the worker entry (fake engine, in process), and the hand-over of a
# termination signal from the driver to its workers.

@testset "Sweep re-entry" begin
    nd = Nbody6Dynamics
    base = mktempdir()
    # Stand-in for the engine, as in test_cluster.jl: integer time, one step per
    # 0.1 s, STOP honoured with a dump, KSTART=2 continues from the dump's time.
    fake_engine = raw"""
    #!/bin/bash
    input=$(cat)
    tcrit=$(printf '%s' "$input" | grep -o 'TCRIT=[0-9.]*' | head -n 1 | cut -d= -f2)
    tcrit=${tcrit%.*}
    dump() {
        printf '\x04\x00\x00\x00%4d\x04\x00\x00\x00\x04\x00\x00\x00abcd\x04\x00\x00\x00' "$2" > "$1"
    }
    if printf '%s' "$input" | grep -q 'KSTART=2'; then
        t0=$(dd if=comm.1 bs=1 skip=4 count=4 2>/dev/null | tr -d ' ')
        dump "comm.1_$t0.0" "$t0"
        printf '  MYDUMP    %d.0000000000000000      %d.0                         201 comm.1_%d.0                \n' "$t0" "$t0" "$t0"
        printf '  W MYDUMP J,II,NPARTMP=         201\n  NA-NS=          85\n'
        t=$((t0 + 1))
        end=$((t0 + tcrit))
    else
        printf '## header\nTIME a\n' > lagr.7
        t=0
        end=$tcrit
    fi
    while [ "$t" -le "$end" ]; do
        printf ' ADJUST:  TIME    %d.00000E+00  T[Myr]   0.100E+01  Q   0.500E+00  DE   0.000E+00 DELTA   0.000E+00 DETOT   0.000E+00 E  -2.500E-01\n' "$t"
        printf '%d.0 1\n' "$t" >> lagr.7
        : > "conf.3_$t"
        [ "$t" -eq "$end" ] && break
        if [ -f STOP ]; then
            printf '\n         TERMINATION BY MANUAL INTERVENTION\n'
            dump "comm.1_$t.0" "$t"
            printf '  MYDUMP    %d.0000000000000000      %d.0                         203 comm.1_%d.0                \n' "$t" "$t" "$t"
            printf '  W MYDUMP J,II,NPARTMP=         203\n  NA-NS=          85\n'
            printf '\n\n         COMMON SAVED AT TOFF/TIME/TTOT =  0.00000000E+00  TCOMP =  0.0\n'
            exit 0
        fi
        sleep 0.1
        t=$((t + 1))
    done
    printf '\n         END RUN    TIME[Myr] =    1.00  TOFF/TIME/TTOT=      0.00000000      %d.00000000      %d.00000000  CPUTOT =    0.0  ERRTOT = 0.00000D+00  DETOT = 0.00000D+00\n' "$end" "$end"
    dump "comm.1_$end.0" "$end"
    printf '  MYDUMP    %d.0000000000000000      %d.0                         204 comm.1_%d.0                \n' "$end" "$end" "$end"
    printf '  W MYDUMP J,II,NPARTMP=         204\n  NA-NS=          85\n'
    """
    engine_path = joinpath(base, "engine", "build", "nbody6++")
    mkpath(dirname(engine_path))
    write(engine_path, replace(fake_engine, r"^ {4}"m => ""))
    chmod(engine_path, 0o755)
    write(
        joinpath(base, "config.toml"),
        """
        [install]
        enabled = false
        install_dir = "engine"

        [build]
        enable_mpi = false
        enable_gpu = false
        enable_hdf5 = false

        [simulation]
        input_file = "unused.inp"
        runs_dir = "runs"
        binary_name = "nbody6++"
        omp_threads = 1
        telemetry_interval = 0.0
        startup_timeout = 0.0
        exit_grace = 0.0
        wall_budget = 3.0
        stop_margin = 1.0

        [postprocess]
        enabled = false

        [visualization]
        enabled = false

        [merger]
        enabled = false
        config_file = ""
        """,
    )
    # TCRIT = 30: the first segment is stopped part-way by the wall budget
    # (2 s of 0.1 s steps), a continuation completes it.
    write(
        joinpath(base, "merger.toml"),
        """
        [merger]
        seed = 5
        n_clusters = 1
        orbit_mode = "single"

        [merger.cluster1]
        N = 100
        rbar = 1.0
        model = "plummer"

        [merger.output]
        tcrit = 30.0
        dtadj = 1.0
        deltat = 1.0

        [merger.stellar]
        kz19 = 0
        """,
    )
    sweep_toml(seeds) = """
        [sweep]
        name = "reentry"
        pipeline_config = "config.toml"
        merger_config = "merger.toml"
        seeds = $(seeds)
        concurrency = 4
        omp_threads = 1
        poll_interval = 0.2
        """
    write(joinpath(base, "sweep.toml"), sweep_toml("[1, 2, 3, 4]"))
    cfg = load_sweep_config(joinpath(base, "sweep.toml"))

    record(dir; status, pipeline = false) = begin
        mkpath(dir)
        text = "[run]\nstatus = \"$status\"\n\n[[segments]]\nindex = 1\nstatus = \"$status\"\n"
        pipeline && (text *= "\n[pipeline]\ncompleted = true\n")
        write(joinpath(dir, "RUN_INFO.toml"), text)
    end

    @testset "state of a point at re-entry" begin
        dir = mktempdir()
        @test nd._sweep_resume_state(joinpath(dir, "absent")) == ("continue", "")
        record(joinpath(dir, "a"); status = "completed", pipeline = true)
        @test nd._sweep_resume_state(joinpath(dir, "a")) == ("done", "completed")
        record(joinpath(dir, "b"); status = "completed")    # not post-processed
        @test nd._sweep_resume_state(joinpath(dir, "b")) == ("continue", "completed")
        for s in ("halted", "watchdog", "stalled")
            record(joinpath(dir, s); status = s)
            @test nd._sweep_resume_state(joinpath(dir, s)) == ("final", s)
        end
        for s in ("stopped", "killed", "failed", "running")
            record(joinpath(dir, s); status = s)
            @test nd._sweep_resume_state(joinpath(dir, s)) == ("continue", s)
        end
    end

    @testset "re-entry relaunches only the unfinished points" begin
        sdir, pts = prepare_sweep(cfg; sweep_dir = joinpath(base, "sweep_stub"))
        run_of(p) = joinpath(sdir, p.id, "run")
        record(run_of(pts[1]); status = "completed", pipeline = true)
        record(run_of(pts[2]); status = "halted")
        record(run_of(pts[3]); status = "stopped")
        # pts[4] never started; the cut driver left every point `running`
        cut =
            Dict(p.index => Dict{String,Any}("status" => "running", "started" => "x") for p in pts)
        write_sweep_index(sdir, cfg, pts; status = cut)
        stub_log = joinpath(base, "stub.log")
        stub_script = """
            echo "\$(basename "\$(pwd -P)")" >> "\$1"
            mkdir -p run
            printf '[run]\\nstatus = "completed"\\n\\n[[segments]]\\nindex = 1\\nstatus = "completed"\\n\\n[pipeline]\\ncompleted = true\\n' > run/RUN_INFO.toml
            """
        stub = pdir -> Cmd(`bash -c $stub_script stub $stub_log`; dir = pdir)
        @test resume_sweep(cfg, sdir; worker_command = stub) == sdir
        @test sort(readlines(stub_log)) == sort([pts[3].id, pts[4].id])
        idx = read_sweep_index(sdir)
        @test idx["sweep"]["entries"] == 2
        @test [p["status"] for p in idx["points"]] == ["done", "halted", "done", "done"]
        @test idx["points"][3]["entry"] == 2 && idx["points"][4]["entry"] == 2
        @test isfile(joinpath(sdir, "sweep_summary.csv"))
        # a second re-entry has nothing left to run
        rm(stub_log)
        resume_sweep(cfg, sdir; worker_command = stub)
        @test !isfile(stub_log) && read_sweep_index(sdir)["sweep"]["entries"] == 3
        # a configuration that no longer describes the indexed points is refused
        write(joinpath(base, "sweep_other.toml"), sweep_toml("[1, 2, 3, 5]"))
        other = load_sweep_config(joinpath(base, "sweep_other.toml"))
        @test_throws ArgumentError resume_sweep(other, sdir; worker_command = stub)
    end

    @testset "a stopped point whose worker fails is not launched again" begin
        write(joinpath(base, "sweep_fail.toml"), sweep_toml("[1]") * "max_segments = 3\n")
        cfgf = load_sweep_config(joinpath(base, "sweep_fail.toml"))
        sdir, pts = prepare_sweep(cfgf; sweep_dir = joinpath(base, "sweep_fail"))
        rdir = joinpath(sdir, pts[1].id, "run")
        mkpath(rdir)
        write(
            joinpath(rdir, "RUN_INFO.toml"),
            "[run]\nstatus = \"stopped\"\n\n[[segments]]\nindex = 1\nstatus = \"stopped\"\n" *
            "t_start = 0.0\nt_end = 4.0\n",
        )
        cut =
            Dict(p.index => Dict{String,Any}("status" => "running", "started" => "x") for p in pts)
        write_sweep_index(sdir, cfgf, pts; status = cut)
        fail_log = joinpath(base, "fail.log")
        fail_script = """
            echo "\$(basename "\$(pwd -P)")" >> "\$1"
            exit 1
            """
        failing = pdir -> Cmd(`bash -c $fail_script stub $fail_log`; dir = pdir)
        @test resume_sweep(cfgf, sdir; worker_command = failing) == sdir
        @test readlines(fail_log) == [pts[1].id]    # once: the record still allows two segments
        idx = read_sweep_index(sdir)
        @test idx["points"][1]["status"] == "stopped" && idx["points"][1]["exit_status"] == 1
    end

    @testset "an interrupted point continues from its dump" begin
        write(joinpath(base, "sweep_one.toml"), sweep_toml("[1]"))
        one = load_sweep_config(joinpath(base, "sweep_one.toml"))
        sdir, pts = prepare_sweep(one; sweep_dir = joinpath(base, "sweep_fake"))
        pdir = joinpath(sdir, only(pts).id)
        run_dir = joinpath(pdir, "run")
        nd.continue_sweep_point(pdir)               # fresh start, stopped by the wall budget
        @test nd._run_status(run_dir) == "stopped"
        n = 0
        while nd._sweep_resume_state(run_dir)[1] == "continue" && n < 6
            nd.continue_sweep_point(pdir)
            n += 1
        end
        @test nd._sweep_resume_state(run_dir) == ("done", "completed")
        info = nd.TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))
        @test length(info["segments"]) ≥ 2 && info["segments"][1]["status"] == "stopped"
        times = [
            parse(Float64, first(split(l))) for
            l in readlines(joinpath(run_dir, "output", "lagr.7"))[3:end]
        ]
        @test times == collect(0.0:30.0)            # one record per step across the joins
        nd.continue_sweep_point(pdir)               # idempotent on a finished point
        @test length(nd.TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))["segments"]) ==
              length(info["segments"])
    end

    @testset "the driver hands a termination on to its workers" begin
        flag = joinpath(base, "worker_terminated")
        proc = run(
            `bash -c "trap 'kill %1 2>/dev/null; echo term > $flag; exit 0' TERM; sleep 30 & wait"`;
            wait = false,
        )
        sleep(0.5)
        nd._register_worker(nd._ActiveWorker(Int32(getpid(proc)), 1.0))
        t0 = time()
        @test nd._terminate_workers()
        @test isfile(flag) && time() - t0 < 10
        @test isempty(nd._ACTIVE_WORKERS[])
        @test !nd._terminate_workers()                       # nothing registered any more
        # No `wait(proc)`: the hand-over reaps the child by `waitpid`, as the exit
        # hook must, so libuv never learns of its exit and a wait would block.
    end
end
