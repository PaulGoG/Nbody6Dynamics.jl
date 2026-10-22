# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope. Run-infrastructure safeguards that need
# no engine: hand-written RUN_INFO.toml records and in-memory configurations.

@testset "Run safety" begin
    TOMLlib = Nbody6Dynamics.TOML
    in_memory(; install = InstallConfig(), build = BuildConfig(), sim = SimulationConfig()) =
        Nbody6Config(
            install,
            build,
            sim,
            PostprocessConfig(),
            VisualizationConfig(),
            MergerPipelineConfig(),
        )

    @testset "Resubmit chain stops after a deterministic failure" begin
        dir = mktempdir()
        info = joinpath(dir, "RUN_INFO.toml")
        record(segments...) = write(
            info,
            "[run]\nstatus = \"$(last(segments).status)\"\n\n" * join(
                (
                    "[[segments]]\nindex = $k\nstatus = \"$(s.status)\"\nt_start = $(s.t_start)\n" *
                    (s.t_end === nothing ? "" : "t_end = $(s.t_end)\n") for
                    (k, s) in enumerate(segments)
                ),
                "\n",
            ),
        )
        seg(status, t_start, t_end) = (status = status, t_start = t_start, t_end = t_end)
        reason = Nbody6Dynamics._chain_stop_reason

        @test reason(joinpath(dir, "absent.toml")) == ""
        # No progress: t_end not beyond t_start, or no t_end at all
        record(seg("failed", 5.0, 5.0))
        @test occursin("without progress", reason(info)) && occursin("segment 1", reason(info))
        record(seg("stopped", 0.0, 4.0), seg("killed", 4.0, nothing))
        @test occursin("without progress", reason(info)) && occursin("killed", reason(info))
        # Progress after a segment that ended normally: the chain continues
        record(seg("stopped", 0.0, 4.0), seg("failed", 4.0, 6.0))
        @test reason(info) == ""
        # Two failed or killed segments in a row stop it even with progress
        record(seg("killed", 0.0, 4.0), seg("failed", 4.0, 6.0))
        @test occursin("both", reason(info))
        # Statuses that are not failures never stop it here
        record(seg("stopped", 0.0, 0.0))
        @test reason(info) == ""

        # The ledger event the chain records when it stops
        Nbody6Dynamics._record_job(
            dir,
            dir;
            event = "chain_stopped",
            job_id = "7",
            segment = "2",
            status = "failed",
            note = "segment 2 ended failed without progress",
        )
        rows = Nbody6Dynamics.read_job_ledger(joinpath(dir, "cluster", "jobs.csv"))
        @test length(rows) == 1 && rows[1].event == "chain_stopped" && rows[1].status == "failed"
        @test "stalled" in Nbody6Dynamics._FINAL_STATUSES
    end

    @testset "A fresh run never reuses a run directory" begin
        dir = mktempdir()
        input = joinpath(dir, "case.inp")
        write(input, "1 100.0\n")
        run_dir = joinpath(dir, "runs", "r1")
        mkpath(run_dir)
        write(joinpath(run_dir, "RUN_INFO.toml"), "[run]\nid = \"r1\"\nstatus = \"stopped\"\n")
        cfg = in_memory(;
            install = InstallConfig(; enabled = false, install_dir = joinpath(dir, "engine")),
            sim = SimulationConfig(; input_file = input, runs_dir = joinpath(dir, "runs")),
        )
        err = try
            Nbody6Dynamics.run_simulation(cfg; base_dir = dir, run_id = "r1")
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin(run_dir, err.msg) &&
              occursin("resume_run", err.msg) &&
              occursin("restart_simulation", err.msg)
        @test !isdir(joinpath(run_dir, "output"))   # refused before anything is created
    end

    @testset "Sweep points are done only when their run completed" begin
        sdir = mktempdir()
        write_status(id, status) = begin
            run = joinpath(sdir, id, "run")
            mkpath(run)
            isempty(status) ||
                write(joinpath(run, "RUN_INFO.toml"), "[run]\nstatus = \"$status\"\n")
            joinpath(sdir, id)
        end
        point_status = Nbody6Dynamics._sweep_point_status
        @test point_status(0, joinpath(write_status("a", "completed"), "run")) == "done"
        @test point_status(0, joinpath(write_status("b", "halted"), "run")) == "halted"
        @test point_status(1, joinpath(write_status("c", "failed"), "run")) == "failed"
        @test point_status(1, joinpath(sdir, "a", "run")) == "failed"   # post-processing failed
        @test point_status(0, joinpath(write_status("d", ""), "run")) == "done"   # no record
        @test point_status(2, joinpath(sdir, "d", "run")) == "failed"

        # An index whose points were all marked done before the run status counted
        points = join(
            (
                "[[points]]\nindex = $k\nid = \"$id\"\nkind = \"merger\"\nstatus = \"done\"\n" *
                "dir = \"$(joinpath(sdir, id))\"\n" for
                (k, id) in enumerate(("a", "b", "c", "d"))
            ),
            "\n",
        )
        write(joinpath(sdir, "sweep_index.toml"), "[sweep]\nname = \"safety\"\n\n" * points)
        _, done =
            @test_logs (:warn, r"2 point\(s\) marked done") match_mode = :any Nbody6Dynamics._sweep_done_points(
                sdir,
            )
        @test [p["id"] for (p, _) in done] == ["a", "d"]
        @test all(r -> basename(r) == "run", last.(done))
    end

    @testset "Configuration combinations and defaults" begin
        dir = mktempdir()
        refused(body, key) = begin
            path = joinpath(dir, "case.toml")
            write(path, body)
            try
                load_config(path)
                false
            catch e
                e isa ArgumentError && occursin(key, e.msg)
            end
        end
        @test refused(
            "[simulation]\nrun_test = true\n\n[postprocess]\ndata_dir = \"elsewhere\"\n",
            "postprocess.data_dir",
        )
        @test !refused(
            "[simulation]\nrun_test = false\n\n[postprocess]\ndata_dir = \"elsewhere\"\n",
            "postprocess",
        )
        @test refused(
            "[postprocess]\nbinary_spans = true\nparse_stdout = false\n",
            "postprocess.binary_spans",
        )
        @test !refused("[postprocess]\nbinary_spans = true\n", "postprocess")
        @test refused("[simulation]\nbinary_name = \"\"\n", "simulation.binary_name")
        @test refused("[visualization.style]\nmarker_min = 0.0\n", "visualization.style.marker_min")
        @test refused(
            "[visualization.style]\nmarker_min = 10.0\nmarker_max = 4.0\n",
            "visualization.style.marker_max",
        )
        pp = PostprocessConfig()
        @test pp.population_m1_min == 0.47 && pp.population_m1_max == 0.76
    end

    @testset "CUDA toolkit of the launch" begin
        out = mktempdir()
        cuda_of = Nbody6Dynamics._launch_cuda_path
        @test cuda_of(BuildConfig(; enable_gpu = true, cuda_path = "/opt/cuda-x"), out) ==
              "/opt/cuda-x"
        gpu = BuildConfig(; enable_gpu = true)
        @test cuda_of(gpu, out) == detect_cuda_path()      # no build record
        write(joinpath(out, "BUILD_INFO.toml"), "enable_gpu = true\ncuda_path = \"\"\n")
        @test cuda_of(gpu, out) == detect_cuda_path()      # record without a path
        write(
            joinpath(out, "BUILD_INFO.toml"),
            "enable_gpu = true\ncuda_path = \"/usr/local/cuda-13.1\"\n",
        )
        @test cuda_of(gpu, out) == "/usr/local/cuda-13.1"
        # The configured path still wins over the record
        @test cuda_of(BuildConfig(; enable_gpu = true, cuda_path = "/opt/cuda-x"), out) ==
              "/opt/cuda-x"
        script = Nbody6Dynamics._write_launch_script(
            out,
            joinpath(out, "nbody6++.avx.gpu"),
            joinpath(out, "case.inp"),
            joinpath(out, "out1000"),
            joinpath(out, "err1000"),
            in_memory(; build = gpu),
        )
        @test occursin("export CUDA_HOME='/usr/local/cuda-13.1'\n", read(script, String))
    end

    @testset "Engine amendments carried by the fork are not rewritten" begin
        tree = mktempdir()
        amendments = TOMLlib.parsefile(Nbody6Dynamics._ENGINE_AMENDMENTS_FILE)["amendment"]
        ids = sort([String(a["id"]) for a in amendments])
        # A fork checkout: every file already holds the replacement text
        for a in amendments
            path = joinpath(tree, a["file"])
            mkpath(dirname(path))
            isfile(path) || write(path, "      SUBROUTINE X\n")
            open(path, "a") do io
                write(io, a["replace"], "      END\n")
            end
        end
        before = Dict(a["file"] => read(joinpath(tree, a["file"]), String) for a in amendments)
        stamps = Dict(a["file"] => mtime(joinpath(tree, a["file"])) for a in amendments)
        result = Nbody6Dynamics._apply_engine_amendments(tree)
        @test sort(result.applied) == ids && isempty(result.skipped)
        for a in amendments
            path = joinpath(tree, a["file"])
            @test read(path, String) == before[a["file"]]
            @test mtime(path) == stamps[a["file"]]
            @test !isfile(path * ".tmp")
        end
        # An upstream checkout: the search text is replaced once
        upstream = mktempdir()
        a = first(amendments)
        path = joinpath(upstream, a["file"])
        mkpath(dirname(path))
        write(path, "      SUBROUTINE X\n" * a["search"] * "\n      END\n")
        # the other amendments' files are absent from this tree: skipped with a warning
        result =
            @test_logs (:warn, r"not in the source tree") match_mode = :any Nbody6Dynamics._apply_engine_amendments(
                upstream,
            )
        @test String(a["id"]) in result.applied
        @test occursin(a["replace"], read(path, String)) &&
              !occursin(a["search"], read(path, String))
    end
end
