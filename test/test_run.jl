# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Engine interval digit counter" begin
    term = Nbody6Dynamics._engine_digit_counter_terminates
    # The intervals of the runs that hung (2026-09-08 sweep controls, 2026-09-10 reproduction)
    for bad in (0.6302, 0.1576, 0.6085, 2.434)
        @test !term(bad)
    end
    # The intervals of the runs that completed, and dyadic values
    for good in (0.5604, 0.2802, 2.241, 0.5, 0.25, 0.1, 1.0, 0.15625, 0.0009765625, 2.521)
        @test term(good)
    end
    @test engine_interval(0.6302) == 161 / 256
    @test engine_interval(0.1576) == 81 / 512     # max_digits = 9 binds: ten digits break the engine's format
    @test engine_interval(0.1576; max_digits = 10) == 161 / 1024
    @test engine_interval(2.521) == 161 / 64
    @test engine_interval(0.5) == 0.5 &&
          engine_interval(1.0) == 1.0 &&
          engine_interval(0.25) == 0.25
    @test engine_interval(0.0) == 0.0 && engine_interval(-1.0) == -1.0
    @test engine_interval(0.01) == 5 / 512        # nearest 2⁻⁹ multiple
    @test engine_interval(1e-9) == 1 / 512        # never rounds to zero
    for dt in (0.6302, 0.1576, 2.521, 0.0636, 0.3043, 12.7, 100.0, 0.05)
        v = engine_interval(dt)
        @test term(v)
        @test v * 512 == round(v * 512)            # dyadic with ≤ 9 binary digits
        s = Nbody6Dynamics._decimal_string(v)
        @test !occursin('.', s) || length(split(s, '.')[2]) ≤ 9
        dt ≥ 0.25 && @test abs(v - dt) / dt ≤ 1 / 256
    end
    @test Nbody6Dynamics._decimal_string(161 / 256) == "0.62890625"
    @test Nbody6Dynamics._decimal_string(0.5) == "0.5"
    @test Nbody6Dynamics._decimal_string(1.0) == "1"
    @test Nbody6Dynamics._decimal_string(1 / 1024) == "0.0009765625"
    @test Nbody6Dynamics._decimal_string(2.515625) == "2.515625"
    # A logged rounding, and silence when the value is already dyadic
    @test (@test_logs (:info, r"DTADJ = 0.62890625 NB") Nbody6Dynamics._engine_interval_logged(
        "DTADJ",
        0.6302,
    )) == 161 / 256
    @test (@test_logs Nbody6Dynamics._engine_interval_logged("DELTAT", 0.5)) == 0.5

    # End to end: Myr intervals through the generator come out dyadic and exact
    ic_dir = mktempdir()
    mt = joinpath(ic_dir, "merger.toml")
    write(
        mt,
        """
        [merger]
        seed = 5
        n_clusters = 2
        orbit_mode = "kepler"

        [merger.cluster1]
        model = "plummer"
        N = 150
        rbar = 1.0
        imf = "kroupa"

        [merger.cluster2]
        model = "plummer"
        N = 150
        rbar = 1.0
        imf = "kroupa"

        [merger.orbit]
        apocentre = 5.0
        eccentricity = 0.5

        [merger.output]
        format = "nbody"
        truncate_jacobi = true
        output_dir = "$(ic_dir)"
        tcrit_myr = 2.0
        dtadj_myr = 0.25
        deltat_myr = 1.0

        [merger.stellar]
        dtplot_myr = 4.0
        """,
    )
    generate_merger_ic(load_merger_config(mt))
    inp = read(joinpath(ic_dir, "merger.inp"), String)
    for key in ("DTADJ", "DELTAT", "DTPLOT")
        m = match(Regex(key * "=([0-9.Ee+-]+)"), inp)
        @test m !== nothing
        v = parse(Float64, m.captures[1])
        @test term(v) && v * 1024 == round(v * 1024)
    end
end

# =====================================================================

@testset "Run completion monitor and exit grace" begin
    dir = mktempdir()
    out = joinpath(dir, "out1000")
    @test !Nbody6Dynamics._run_completed(out)
    write(out, " ADJUST:  TIME 1.0\n")
    @test !Nbody6Dynamics._run_completed(out)
    write(
        out,
        " ADJUST:  TIME 1.0\n\n         END RUN    TIME[Myr] =   20.50  TOFF/TIME/TTOT=  0.0\n rank PE N Total\n",
    )
    @test Nbody6Dynamics._run_completed(out)
    # The scan covers only the tail of a long file
    long = joinpath(dir, "long")
    open(long, "w") do io
        println(io, "END RUN")
        for _ in 1:20000
            println(io, "ADJUST: TIME 1.0 filler line to push the marker out of the window")
        end
    end
    @test !Nbody6Dynamics._run_completed(long)

    # A process that keeps running after END RUN is terminated after the grace period
    p = run(`sleep 30`; wait = false)
    mon = Nbody6Dynamics._start_completion_monitor(out, p, 1.0)
    wait(p)
    @test mon.fired[] && !process_running(p) && Nbody6Dynamics._exit_status(p) == -15
    # A process that exits by itself leaves the monitor silent
    q = run(`sleep 1`; wait = false)
    mon_q = Nbody6Dynamics._start_completion_monitor(out, q, 30.0)
    wait(q)
    mon_q.stop[] = true
    @test !mon_q.fired[] && Nbody6Dynamics._exit_status(q) == 0
    # Without END RUN nothing fires within the grace period
    r = run(`sleep 3`; wait = false)
    mon_r = Nbody6Dynamics._start_completion_monitor(joinpath(dir, "absent"), r, 1.0)
    wait(r)
    mon_r.stop[] = true
    @test !mon_r.fired[]
    # A file still being written (the final COMMON dump) keeps the engine alive
    busy = mktempdir()
    busy_out = joinpath(busy, "out1000")
    write(busy_out, "END RUN\n")
    @test Nbody6Dynamics._latest_mtime(joinpath(busy, "none")) == 0.0
    @test time() - Nbody6Dynamics._latest_mtime(busy) < 5
    b = run(`sleep 7`; wait = false)
    mon_b = Nbody6Dynamics._start_completion_monitor(busy_out, b, 2.0)
    writer = @async for _ in 1:12
        open(joinpath(busy, "comm.1"), "a") do io
            write(io, "x")
        end
        sleep(0.5)
    end
    wait(writer)
    @test process_running(b) && !mon_b.fired[]      # 6 s of writes, never idle for 2 s
    wait(b)
    mon_b.stop[] = true
    @test !mon_b.fired[] && Nbody6Dynamics._exit_status(b) == 0

    # Config key
    @test SimulationConfig().exit_grace == 120.0
    cfg_path = joinpath(dir, "c.toml")
    write(cfg_path, "[simulation]\nexit_grace = 0\n")
    @test load_config(cfg_path).simulation.exit_grace == 0.0
    write(cfg_path, "[simulation]\nexit_grace = -1\n")
    @test_throws ErrorException load_config(cfg_path)

    # Sweep outcome: a completed-then-terminated segment counts as completed
    run_dir = joinpath(dir, "run")
    mkpath(run_dir)
    write(
        joinpath(run_dir, "RUN_INFO.toml"),
        "[run]\nelapsed_seconds = 7.5\n[[segments]]\nexit_status = -15\ncompleted = true\n",
    )
    oc = Nbody6Dynamics._sweep_point_outcome(run_dir)
    @test oc["exit_status"] == -15 && oc["completed"] == true && oc["elapsed_seconds"] == 7.5
    write(
        joinpath(run_dir, "RUN_INFO.toml"),
        "[run]\nelapsed_seconds = 1.0\n[[segments]]\nexit_status = 0\n",
    )
    @test Nbody6Dynamics._sweep_point_outcome(run_dir)["completed"] == true
    write(
        joinpath(run_dir, "RUN_INFO.toml"),
        "[run]\nelapsed_seconds = 1.0\n[[segments]]\nexit_status = 1\n",
    )
    @test Nbody6Dynamics._sweep_point_outcome(run_dir)["completed"] == false

    # The sweep's pre-launch binary check honours the build variant
    sw = mktempdir()
    bk = joinpath(sw, "backend", "build")
    mkpath(bk)
    touch(joinpath(bk, "nbody6++.avx"))
    write(
        joinpath(sw, "pipeline.toml"),
        "[install]\ninstall_dir = \"backend\"\n[simulation]\ninput_file = \"x.inp\"\n",
    )
    cp(joinpath(@__DIR__, "..", "input_files", "merger_demo_small.toml"), joinpath(sw, "m.toml"))
    st = joinpath(sw, "sweep.toml")
    write(
        st,
        "[sweep]\nname = \"t\"\npipeline_config = \"pipeline.toml\"\nmerger_config = \"m.toml\"\nseeds = [1]\n",
    )
    @test Nbody6Dynamics._check_sweep_binary(load_sweep_config(st)) == joinpath(bk, "nbody6++.avx")
    write(
        joinpath(sw, "pipeline.toml"),
        "[install]\ninstall_dir = \"backend\"\n[build]\nenable_gpu = true\n[simulation]\ninput_file = \"x.inp\"\n",
    )
    @test_throws ErrorException Nbody6Dynamics._check_sweep_binary(load_sweep_config(st))
end

# =====================================================================

@testset "Run ID generation" begin
    id1 = generate_run_id()
    @test startswith(id1, "run_")
    @test length(id1) > 15

    id2 = generate_run_id("bench")
    @test startswith(id2, "bench_")

    # Two successive IDs should differ (hex suffix)
    id3 = generate_run_id()
    @test id1 != id3
end

# =====================================================================

@testset "Machine identity" begin
    # Model strings are stripped of vendor and marketing boilerplate.
    @test Nbody6Dynamics._hardware_tag("13th Gen Intel(R) Core(TM) i9-13900KS") == "i9-13900KS"
    @test Nbody6Dynamics._hardware_tag("AMD Ryzen 9 9950X 16-Core Processor") == "Ryzen-9-9950X"
    @test Nbody6Dynamics._hardware_tag("NVIDIA GeForce RTX 5090, 32607 MiB, 610.57.04, 12.0") ==
          "RTX-5090"
    @test Nbody6Dynamics._hardware_tag("Tesla T4, 15360 MiB, 575.57.08, 7.5") == "Tesla-T4"
    @test Nbody6Dynamics._hardware_tag("") == ""
    @test Nbody6Dynamics._hardware_tag("unavailable") == ""
    @test length(Nbody6Dynamics._hardware_tag("A"^40)) ≤ 16

    # The whole point: one hostname, three machines, three identities.
    shared = "workstation-01"
    ids = map((
        ("AMD Ryzen 9 9950X 16-Core Processor", "NVIDIA GeForce RTX 5090, 32607 MiB"),
        ("13th Gen Intel(R) Core(TM) i9-13900KS", "NVIDIA GeForce RTX 5070 Ti, 16303 MiB"),
        ("13th Gen Intel(R) Core(TM) i9-13900KS", "AMD Radeon Pro W7900, 46068 MiB"),
    )) do (cpu, gpu)
        Nbody6Dynamics._machine_id(Dict("host" => shared, "cpu_model" => cpu, "gpu" => gpu))
    end
    @test length(unique(ids)) == 3
    @test all(startswith(id, shared) for id in ids)
    @test ids[1] == "workstation-01-Ryzen-9-9950X-RTX-5090"

    # Only the first device of a multi-GPU record contributes.
    @test Nbody6Dynamics._machine_id(
        Dict(
            "host" => "h",
            "cpu_model" => "unknown",
            "gpu" => "Tesla T4, 15360 MiB; Tesla T4, 15360 MiB",
        ),
    ) == "h-Tesla-T4"

    # Without usable hardware strings the identity degrades to the hostname.
    @test Nbody6Dynamics._machine_id(Dict("host" => "plain")) == "plain"
    @test Nbody6Dynamics._machine_id(
        Dict("host" => "plain", "cpu_model" => "unknown", "gpu" => "unavailable"),
    ) == "plain"

    # The fingerprint always carries one.
    fp = Nbody6Dynamics._hardware_fingerprint()
    @test haskey(fp, "machine")
    @test startswith(fp["machine"], fp["host"])
end

# =====================================================================

@testset "Pipeline completion marker" begin
    mktempdir() do dir
        # No summary to amend.
        @test Nbody6Dynamics._stamp_pipeline_completion(dir, ["simulation"], 1.0) == false
        @test Nbody6Dynamics._pipeline_completed(dir) == false
        @test Nbody6Dynamics._stamp_pipeline_completion("", ["simulation"], 1.0) == false

        # A summary written by the engine phase alone is not "completed":
        # this is the state the three killed stages of 2026-09-11 left.
        info = joinpath(dir, "RUN_INFO.toml")
        open(info, "w") do io
            Nbody6Dynamics.TOML.print(
                io,
                Dict("run" => Dict("id" => "merger_cpu_x", "segments" => 1)),
            )
        end
        @test Nbody6Dynamics._pipeline_completed(dir) == false

        @test Nbody6Dynamics._stamp_pipeline_completion(
            dir,
            ["merger_ic", "simulation", "postprocess", "plots"],
            12.25,
        )
        @test Nbody6Dynamics._pipeline_completed(dir)

        # Amending preserves what the run summary already held.
        parsed = Nbody6Dynamics.TOML.parsefile(info)
        @test parsed["run"]["id"] == "merger_cpu_x"
        @test parsed["pipeline"]["phases"] == ["merger_ic", "simulation", "postprocess", "plots"]
        @test parsed["pipeline"]["elapsed_seconds"] ≈ 12.2 atol = 0.1

        # Unreadable summary: false, not an exception.
        write(info, "this is not TOML {{{")
        @test Nbody6Dynamics._pipeline_completed(dir) == false
    end
end

# =====================================================================

@testset "Commit stamp only at a repository root" begin
    quiet(cmd) = run(pipeline(cmd; stdout = devnull, stderr = devnull))
    mktempdir() do repo
        quiet(`git -C $repo init --quiet`)
        write(joinpath(repo, "f.txt"), "x")
        quiet(`git -C $repo add f.txt`)
        quiet(`git -C $repo -c user.email=t@t -c user.name=t commit --quiet -m init`)
        @test Nbody6Dynamics._git_commit(repo) != "unknown"
        # A plain subdirectory must not inherit the enclosing commit.
        sub = joinpath(repo, "backend")
        mkpath(sub)
        @test Nbody6Dynamics._git_commit(sub) == "unknown"
    end
end

# =====================================================================

@testset "Source provenance stamp" begin
    # A tree deployed by file copy has no .git: the version of its
    # Project.toml stands in for the commit rather than "unknown".
    mktempdir() do dir
        @test Nbody6Dynamics._source_stamp(dir) == "unknown"

        write(
            joinpath(dir, "Project.toml"),
            """
            name = "Deployed"
            uuid = "11111111-2222-3333-4444-555555555555"
            version = "0.3.1"
            """,
        )
        @test Nbody6Dynamics._source_stamp(dir) == "v0.3.1+nogit"

        # A Project.toml without a version stays unknown.
        write(joinpath(dir, "Project.toml"), "name = \"Deployed\"\n")
        @test Nbody6Dynamics._source_stamp(dir) == "unknown"
    end

    # This repository is a git checkout, so the commit wins.
    stamp = Nbody6Dynamics._source_stamp(Nbody6Dynamics._PROJECT_ROOT)
    @test stamp == "unknown" || !occursin("+nogit", stamp)
end

# =====================================================================

@testset "Elapsed time formatting" begin
    fmt = Nbody6Dynamics._format_elapsed

    # Sub-minute
    @test fmt(0.0) == "0.0 s"
    @test fmt(1.5) == "1.5 s"
    @test fmt(59.4) == "59.4 s"

    # Minutes
    @test fmt(60.0) == "1m 00s"
    @test fmt(122.0) == "2m 02s"
    @test fmt(3599.0) == "59m 59s"

    # Hours
    @test fmt(3600.0) == "1h 00m 00s"
    @test fmt(3661.0) == "1h 01m 01s"
    @test fmt(7384.0) == "2h 03m 04s"
end

# =====================================================================

@testset "Run summary and hardware provenance" begin
    hw = Nbody6Dynamics._hardware_fingerprint()
    for key in (
        "host",
        "os",
        "cpu_model",
        "cpu_threads",
        "total_memory_gib",
        "julia_version",
        "julia_threads",
        "blas_threads",
        "versioninfo",
    )
        @test haskey(hw, key)
    end
    @test occursin("Julia Version", hw["versioninfo"])

    # The resolved manifest travels with the run; an identical snapshot
    # is kept, a different one backed up rather than overwritten.
    mktempdir() do run_dir
        name = Nbody6Dynamics._snapshot_manifest(run_dir)
        @test name == "environment_manifest.toml"
        snapshot = joinpath(run_dir, name)
        @test occursin("julia_version", read(snapshot, String))
        @test Nbody6Dynamics._snapshot_manifest(run_dir) == name
        @test !isfile(joinpath(run_dir, "environment_manifest#1.toml"))
        write(snapshot, "stale = true\n")
        @test Nbody6Dynamics._snapshot_manifest(run_dir) == name
        @test read(joinpath(run_dir, "environment_manifest#1.toml"), String) == "stale = true\n"
    end
    @test hw["cpu_threads"] ≥ 1
    @test hw["total_memory_gib"] > 0
    @test hw["blas_threads"] ≥ 1
    @test !haskey(hw, "gpu")   # probed only on request
    hw_gpu = Nbody6Dynamics._hardware_fingerprint(; gpu_probe = true)
    @test haskey(hw_gpu, "gpu") && !isempty(hw_gpu["gpu"])

    # RUN_INFO.toml writer: structure and round-trip
    run_dir = mktempdir()
    out_dir = joinpath(run_dir, "output")
    mkpath(out_dir)
    write(joinpath(out_dir, "out1000"), "line1\nline2\n")
    cfg_info = Nbody6Config(
        InstallConfig(),
        BuildConfig(),
        SimulationConfig(),
        PostprocessConfig(),
        VisualizationConfig(),
        MergerPipelineConfig(),
    )
    Nbody6Dynamics._write_run_summary(
        cfg_info,
        run_dir,
        "testrun_a1b2",
        joinpath(out_dir, "out1000"),
        out_dir,
        42.5,
    )
    info = Nbody6Dynamics.TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))
    @test info["run"]["id"] == "testrun_a1b2"
    @test info["run"]["elapsed_seconds"] == 42.5
    @test info["run"]["stdout_lines"] == 2
    @test haskey(info["provenance"], "package_commit")
    @test haskey(info["provenance"], "backend_commit")
    @test info["hardware"]["cpu_threads"] ≥ 1
    @test "out1000" in info["output"]["files"]
    @test info["output"]["total_bytes"] > 0
end

# =====================================================================

@testset "export_for_paper provenance" begin
    src_dir = mktempdir()
    dest = mktempdir()
    fig_path = joinpath(src_dir, "plots", "energy.pdf")
    mkpath(dirname(fig_path))
    write(fig_path, "pdfbytes")
    write(
        joinpath(src_dir, "RUN_INFO.toml"),
        """
        [provenance]
        package_commit = "abc1234"
        backend_commit = "def5678-dirty"
        """,
    )
    out = export_for_paper([fig_path], dest; run_dir = src_dir)
    @test length(out) == 1
    @test isfile(out[1])
    @test startswith(basename(out[1]), basename(src_dir) * "__")
    sidecar = out[1] * ".provenance.toml"
    @test isfile(sidecar)
    prov = Nbody6Dynamics.TOML.parsefile(sidecar)
    @test prov["package_commit"] == "abc1234"
    @test prov["backend_commit"] == "def5678-dirty"
    # Never-overwrite: second export backs up, both exist
    export_for_paper([fig_path], dest; run_dir = src_dir)
    @test length(filter(f -> endswith(f, ".pdf"), readdir(dest))) ≥ 2
end

# =====================================================================

@testset "Backend thread control" begin
    sim0 = SimulationConfig()
    @test sim0.omp_threads == 0
    @test sim0.telemetry_interval == 5.0
    @test Nbody6Dynamics._effective_omp_threads(SimulationConfig(; omp_threads = 3)) == 3
    withenv("OMP_NUM_THREADS" => nothing) do
        @test Nbody6Dynamics._effective_omp_threads(sim0) == Sys.CPU_THREADS
    end
    withenv("OMP_NUM_THREADS" => "4") do
        @test Nbody6Dynamics._effective_omp_threads(sim0) == 4
        @test Nbody6Dynamics._effective_omp_threads(SimulationConfig(; omp_threads = 2)) == 2
    end
    withenv("OMP_NUM_THREADS" => "junk") do
        @test Nbody6Dynamics._effective_omp_threads(sim0) == Sys.CPU_THREADS
    end

    # Launch script: OMP_NUM_THREADS exported only for an explicit cap
    launch_dir = mktempdir()
    make_cfg(n) = Nbody6Config(
        InstallConfig(),
        BuildConfig(),
        SimulationConfig(; omp_threads = n),
        PostprocessConfig(),
        VisualizationConfig(),
        MergerPipelineConfig(),
    )
    args = (joinpath(launch_dir, "nbody6++"), "in.inp", "out1000", "err1000")
    script0 = read(Nbody6Dynamics._write_launch_script(launch_dir, args..., make_cfg(0)), String)
    @test !occursin("OMP_NUM_THREADS", script0)
    @test occursin("OMP_STACKSIZE", script0)
    script4 = read(Nbody6Dynamics._write_launch_script(launch_dir, args..., make_cfg(4)), String)
    @test occursin("export OMP_NUM_THREADS=4\n", script4)

    # Reported thread count from the backend's start-up banner
    out_path = joinpath(launch_dir, "out1000")
    write(out_path, "header\n RANK:  0  OpenMP Number of Threads:  12\n ADJUST: ...\n")
    @test Nbody6Dynamics._reported_omp_threads(out_path) == 12
    write(out_path, "no banner here\n")
    @test Nbody6Dynamics._reported_omp_threads(out_path) === nothing
    @test Nbody6Dynamics._reported_omp_threads(joinpath(launch_dir, "absent")) === nothing
end

# =====================================================================

@testset "Restart bookkeeping" begin
    # Dump selection by time suffix
    d = mktempdir()
    for f in ("comm.1_0.0", "comm.2_2.0", "comm.1_4.0", "comm.1", "conf.3_1", "comm.2_x")
        touch(joinpath(d, f))
    end
    @test Nbody6Dynamics._dump_time("comm.1_12.5") == 12.5
    @test isnan(Nbody6Dynamics._dump_time("comm.1")) && isnan(Nbody6Dynamics._dump_time("conf.3_1"))
    @test Nbody6Dynamics._latest_dump(d) == "comm.1_4.0"
    @test Nbody6Dynamics._latest_dump(mktempdir()) === nothing

    # Restart input from an original merger input: KSTART=2, TCRIT as increment,
    # only the first two namelists
    orig = joinpath(d, "merger.inp")
    write(
        orig,
        "&INNBODY6\nKSTART=1,TCOMP=1E+08,TCRTP0=3600,isernb=40,iserreg=40,iserks=0 /\n\n" *
        "&ININPUT\nN=1957,NFIX=1,NCRIT=10,NRAND=7,NNBOPT=44,NRUN=1,NCOMM=10,\n" *
        "ETAI=0.02,ETAR=0.02,RS0=0.1415,DTADJ=0.5000,DELTAT=1.0000,TCRIT=5.00,QE=2.000E-04,RBAR=6.3,ZMBAR=0.6,\n" *
        "KZ(1:10)=1 -1 2 0 0 0 3 0 0 0\nKZ(11:20)=0 1 0 0 0 0 0 0 3 0\n" *
        "DTMIN=2.698E-05,RMIN=6.128E-04,ETAU=0.1,ECLOSE=1,GMIN=1.000E-06,GMAX=0.01,SMAX=1,\n" *
        "Level='C' /\n\n&INSSE /\n\n&INDATA\nALPHAS=2.35 /\n",
    )
    rin = joinpath(d, "restart.inp")
    Nbody6Dynamics._write_restart_inp(rin, orig, 3.0; tcrtp0 = 900.0)
    txt = read(rin, String)
    @test occursin("KSTART=2,TCOMP=1E+08,TCRTP0=900,isernb=40", txt)
    @test occursin("TCRIT=3.0000,QE=2.000E-04", txt)
    @test occursin("KZ(11:20)=0 1 0 0 0 0 0 0 3 0", txt) && occursin("Level='C' /", txt)
    @test !occursin("&INSSE", txt) && !occursin("&INDATA", txt)
    @test_throws ArgumentError Nbody6Dynamics._write_restart_inp(rin, orig, 0.0)

    # RUN_INFO segments accumulate across launches
    run_dir = mktempdir()
    out_dir = joinpath(run_dir, "output")
    mkpath(out_dir)
    write(joinpath(out_dir, "out1000"), "x\n")
    cfg_r = Nbody6Config(
        InstallConfig(),
        BuildConfig(),
        SimulationConfig(),
        PostprocessConfig(),
        VisualizationConfig(),
        MergerPipelineConfig(),
    )
    seg(i, kind, el) = Dict{String,Any}(
        "index" => i,
        "kind" => kind,
        "elapsed_seconds" => el,
        "dump" => kind == "restart" ? "comm.1_4.0" : "",
    )
    Nbody6Dynamics._write_run_summary(
        cfg_r,
        run_dir,
        "r1",
        joinpath(out_dir, "out1000"),
        out_dir,
        10.0;
        input_file = "merger.inp",
        segment = seg(1, "initial", 10.0),
    )
    @test Nbody6Dynamics._segment_count(joinpath(run_dir, "RUN_INFO.toml")) == 1
    Nbody6Dynamics._write_run_summary(
        cfg_r,
        run_dir,
        "r1",
        joinpath(out_dir, "out1000"),
        out_dir,
        5.0;
        input_file = "restart.inp",
        segment = seg(2, "restart", 5.0),
    )
    info = Nbody6Dynamics.TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))
    @test info["run"]["segments"] == 2
    @test info["run"]["input_file"] == "merger.inp"   # original kept across restarts
    @test info["run"]["elapsed_seconds"] == 15.0
    @test length(info["segments"]) == 2 && info["segments"][2]["dump"] == "comm.1_4.0"

    # Launch script append mode for restarts
    args = (joinpath(run_dir, "nbody6++"), "in.inp", "out1000", "err1000")
    s_app =
        read(Nbody6Dynamics._write_launch_script(run_dir, args..., cfg_r; append = true), String)
    @test occursin(">> \"out1000\" 2>> \"err1000\"", s_app)
    s_new = read(Nbody6Dynamics._write_launch_script(run_dir, args..., cfg_r), String)
    @test occursin("> \"out1000\" 2> \"err1000\"", s_new) && !occursin(">>", s_new)

    # Signal terminations are recorded as negative signal numbers
    p_sig = run(`sleep 5`; wait = false)
    kill(p_sig)
    wait(p_sig)
    @test Nbody6Dynamics._exit_status(p_sig) == -15
    p_ok = run(`true`)
    @test Nbody6Dynamics._exit_status(p_ok) == 0

    # Start-up watchdog: stdout progress detection and termination of a stalled process
    wd_out = joinpath(run_dir, "wd_out1000")
    write(wd_out, " ADJUST:  TIME    0.00000E+00  T[Myr]   0.000E+00  Q   0.450E+00\n")
    @test !Nbody6Dynamics._adjust_advanced(wd_out)
    write(
        wd_out,
        read(wd_out, String) * " ADJUST:  TIME    5.00000E-01  T[Myr]   0.338E+01  Q   0.457E+00\n",
    )
    @test Nbody6Dynamics._adjust_advanced(wd_out)
    @test !Nbody6Dynamics._adjust_advanced(joinpath(run_dir, "absent"))
    stalled = run(`sleep 30`; wait = false)
    wd = Nbody6Dynamics._start_startup_watchdog(joinpath(run_dir, "absent"), stalled, 1.0)
    wait(stalled)
    wait(wd.task)
    @test wd.fired[] && !process_running(stalled)
    @test Nbody6Dynamics._exit_status(stalled) == -15
    fine = run(`sleep 2`; wait = false)
    wd2 = Nbody6Dynamics._start_startup_watchdog(wd_out, fine, 30.0)
    wait(wd2.task)
    @test !wd2.fired[]
    wd2.stop[] = true
    wait(fine)

    # restart_simulation refuses runs without the bookkeeping it needs
    @test_throws ErrorException restart_simulation(mktempdir(); tcrit_extra = 1.0)
end

# =====================================================================

@testset "Live diagnostics panel" begin
    fixture = joinpath(@__DIR__, "fixtures", "out1000")
    panel = Nbody6Dynamics._live_diagnostics_panel(fixture)
    @test panel isa String
    @test occursin("Q", panel) && occursin("t [", panel)
    @test count(==('\n'), panel) ≥ 8
    dir = mktempdir()
    one = joinpath(dir, "out1000")
    write(one, " ADJUST:  TIME    0.00000E+00  T[Myr]   0.000E+00  Q   0.450E+00  DE   0.000E+00\n")
    @test Nbody6Dynamics._live_diagnostics_panel(one) === nothing
    @test Nbody6Dynamics._live_diagnostics_panel(joinpath(dir, "absent")) === nothing
    # Config keys
    @test SimulationConfig().live_diagnostics == false && SimulationConfig().live_interval == 30.0
    cfg_path = joinpath(dir, "c.toml")
    write(cfg_path, "[simulation]\nlive_diagnostics = true\nlive_interval = 5\n")
    c = load_config(cfg_path)
    @test c.simulation.live_diagnostics && c.simulation.live_interval == 5.0
    write(cfg_path, "[simulation]\nlive_interval = 0.5\n")
    @test_throws ErrorException load_config(cfg_path)
end

# =====================================================================
