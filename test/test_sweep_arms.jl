# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Arms of a sweep (roadmap 4.3): named bundles of merger and pipeline
# overrides as the outermost axis, paired with identical initial conditions
# across the arms, the perturbed twin of the generator, CPU points without a
# device slot, and the next segment of a stopped point within one entry.

@testset "Sweep arms and the perturbed twin" begin
    nd = Nbody6Dynamics
    base = mktempdir()
    write(
        joinpath(base, "config.toml"),
        """
        [install]
        enabled = false
        install_dir = "engine"

        [build]
        enable_gpu = true

        [simulation]
        input_file = "unused.inp"
        runs_dir = "runs"
        omp_threads = 8

        [postprocess]
        enabled = false

        [visualization]
        enabled = false

        [merger]
        enabled = false
        config_file = ""
        """,
    )
    write(
        joinpath(base, "merger.toml"),
        """
        [merger]
        seed = 5
        n_clusters = 1
        orbit_mode = "single"

        [merger.cluster1]
        N = 400
        rbar = 1.0
        model = "plummer"

        [merger.cluster1.binaries]
        fraction = 0.1

        [merger.output]
        tcrit = 10.0
        dtadj = 1.0
        deltat = 1.0

        [merger.stellar]
        kz19 = 0
        """,
    )
    sweep_toml(arms; extra = "") = """
        [sweep]
        name = "arms"
        pipeline_config = "config.toml"
        merger_config = "merger.toml"
        seeds = [21, 22]
        devices = [0, 1]
        engines_per_device = 2
        $(extra)

        [sweep.grid]
        "merger.cluster1.rbar" = [1.0]

        $(arms)
        """
    arms_text = """
        [sweep.arms.base]

        [sweep.arms.twin]
        "merger.output.twin_amplitude" = 1e-10

        [sweep.arms.threads16]
        "pipeline.simulation.omp_threads" = 16

        [sweep.arms.cpu]
        "pipeline.build.enable_gpu" = false
        "pipeline.install.install_dir" = "engine_avx"
        "pipeline.simulation.wall_budget" = 30.0
        "pipeline.simulation.stop_margin" = 5.0
        """
    path = joinpath(base, "sweep.toml")
    write(path, sweep_toml(arms_text; extra = "max_segments = 4"))
    cfg = load_sweep_config(path)
    @test first.(cfg.arms) == ["base", "cpu", "threads16", "twin"]       # name order
    @test isempty(last(cfg.arms[1])) && cfg.max_segments == 4

    @testset "points, configurations and records" begin
        pts = sweep_points(cfg)
        @test length(pts) == 8                                          # 4 arms × 1 grid point × 2 seeds
        @test [p.arm for p in pts] == repeat(["base", "cpu", "threads16", "twin"]; inner = 2)
        @test pts[1].id == "001_arm=base_cluster1-rbar=1_seed=21"
        @test pts[8].id == "008_arm=twin_cluster1-rbar=1_seed=22"
        sdir, _ = prepare_sweep(cfg; sweep_dir = joinpath(base, "sweep"))
        conf(p) = nd.TOML.parsefile(joinpath(sdir, p.id, "config.toml"))
        merg(p) = nd.TOML.parsefile(joinpath(sdir, p.id, "merger.toml"))
        @test conf(pts[5])["simulation"]["omp_threads"] == 16         # after the sweep's own override
        @test conf(pts[1])["simulation"]["omp_threads"] == cfg.omp_threads
        cpu = conf(pts[3])
        @test cpu["build"]["enable_gpu"] == false &&
              cpu["install"]["install_dir"] == joinpath(base, "engine_avx")
        @test cpu["simulation"]["wall_budget"] == 30.0
        @test merg(pts[7])["merger"]["output"]["twin_amplitude"] == 1e-10
        @test !haskey(merg(pts[1])["merger"]["output"], "twin_amplitude")
        @test nd._point_uses_gpu(joinpath(sdir, pts[1].id)) &&
              !nd._point_uses_gpu(joinpath(sdir, pts[3].id))
        idx = read_sweep_index(sdir)
        @test idx["sweep"]["arms"] == ["base", "cpu", "threads16", "twin"]
        @test [p["arm"] for p in idx["points"]] == [p.arm for p in pts]
        cols, rows = sweep_summary(sdir)
        @test cols[5] == "arm" && [r["arm"] for r in rows] == [p.arm for p in pts]

        # Identical initial conditions across the arms of a seed; the twin moves one body.
        gen(p) = (
            out = joinpath(base, "ic", p.id);
            nd.generate_merger_ic(
                load_merger_config(joinpath(sdir, p.id, "merger.toml"));
                output_dir = out,
            );
            out
        )
        dirs = Dict(p.arm => gen(p) for p in pts if p.seed == 21)
        dat(arm) = readlines(joinpath(dirs[arm], "dat.10"))
        sha(arm) = nd.TOML.parsefile(joinpath(dirs[arm], "merger_ic.toml"))["meta"]["dat10_sha256"]
        @test dat("base") == dat("threads16") == dat("cpu")
        @test sha("base") == sha("threads16") == sha("cpu") && length(sha("base")) == 64
        changed = findall(dat("base") .!= dat("twin"))
        twin = nd.TOML.parsefile(joinpath(dirs["twin"], "merger_ic.toml"))["twin"]
        @test length(changed) == 1 && changed[1] == twin["body"]
        @test twin["amplitude"] == 1e-10 && twin["displacement"] > 0
        @test sha("twin") != sha("base")
        nb = parse(
            Int,
            match(r"NBIN0=(\d+)", read(joinpath(dirs["base"], "merger.inp"), String)).captures[1],
        )
        @test twin["body"] == 2nb + 1                                   # the first single star
        xb = parse.(Float64, split(dat("base")[twin["body"]]))[2:4]
        xt = parse.(Float64, split(dat("twin")[twin["body"]]))[2:4]
        @test isapprox(xt, xb .* (1 + 1e-10); rtol = 1e-13)
    end

    @testset "refusals" begin
        refused(arms, key) =
            try
                write(path, sweep_toml(arms))
                load_sweep_config(path)
                false
            catch e
                e isa ArgumentError && occursin(key, e.msg)
            end
        @test refused(
            "[sweep.arms.a]\n\"merger.cluster1.rbar\" = 2.0\n",
            "sweep.arms.a.merger.cluster1.rbar",
        )
        @test refused("[sweep.arms.a]\n\"merger.seed\" = 3\n", "sweep.arms.a.merger.seed")
        # keys the generator reads change the initial conditions; structural pipeline keys
        for key in (
            "merger.cluster1.W0",
            "merger.output.truncate_jacobi",
            "merger.tidal.kz14",
            "merger.orbit.eccentricity",
            "merger.output.format",
        )
            @test refused("[sweep.arms.a]\n\"$key\" = 1\n", "initial conditions")
        end
        @test refused(
            "[sweep.arms.a]\n\"pipeline.simulation.runs_dir\" = \"x\"\n",
            "structural key",
        )
        @test refused("[sweep.arms.a]\n\"pipeline.merger.config_file\" = \"x\"\n", "structural key")
        accepted =
            "[sweep.arms.a]\n\"merger.nbody6.etai\" = 0.01\n\"merger.stellar.level\" = \"B\"\n" *
            "\"merger.output.dtadj\" = 0.25\n\"pipeline.simulation.omp_threads\" = 3\n"
        write(path, sweep_toml(accepted))
        @test length(last(only(load_sweep_config(path).arms))) == 4
        @test refused("[sweep.arms.a]\n\"simulation.omp_threads\" = 3\n", "sweep.arms.a.simulation")
        @test refused(
            "[sweep.arms.a]\nmerger.output.twin_amplitude = 1e-10\n",
            "sweep.arms.a.merger",
        )
        @test refused("[sweep.arms.\"b c\"]\n", "sweep.arms.b c")
        write(path, sweep_toml(arms_text; extra = "max_segments = 0"))
        @test_throws ArgumentError load_sweep_config(path)
        bad = joinpath(base, "merger_bad.toml")
        write(
            bad,
            replace(
                read(joinpath(base, "merger.toml"), String),
                "deltat = 1.0\n" => "deltat = 1.0\ntwin_amplitude = 1e-3\n",
            ),
        )
        @test_throws ArgumentError load_merger_config(bad)                # beyond 1e-6
    end

    @testset "exit status of the driver hook" begin
        @test Nbody6Dynamics._driver_exit_status() == 143
        Nbody6Dynamics._DRIVER_FAILED[] = true
        @test Nbody6Dynamics._driver_exit_status() == 1
        Nbody6Dynamics._DRIVER_FAILED[] = false
    end

    @testset "the next segment within an entry" begin
        dir = mktempdir()
        seg(t0, t1) = "[[segments]]\nindex = 1\nstatus = \"stopped\"\nt_start = $t0\nt_end = $t1\n"
        write(joinpath(dir, "RUN_INFO.toml"), seg(0.0, 4.0))
        @test nd._continue_in_entry(dir, 4)
        @test !nd._continue_in_entry(dir, 1)                         # the default: one segment per entry
        write(joinpath(dir, "RUN_INFO.toml"), seg(4.0, 4.0))
        @test !nd._continue_in_entry(dir, 4)                         # no progress
        @test !nd._continue_in_entry(joinpath(dir, "absent"), 4)
    end
end
