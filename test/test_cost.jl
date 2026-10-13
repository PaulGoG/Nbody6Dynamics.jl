@testset "Cost model" begin
    @testset "Rate from the table" begin
        @test cost_rate(1e5) == 64.0
        @test cost_rate(5e5) == 281.0
        @test cost_rate(1e6) == 1087.0
        # Log–log interpolation: the midpoint in log N gives the geometric mean
        @test cost_rate(sqrt(5e5 * 1e6)) ≈ sqrt(281.0 * 1087.0)
        # Extrapolation with the exponent of the nearest segment
        @test cost_rate(1e4) ≈ 64 * 10^(-log10(281 / 64) / log10(5))
        @test cost_rate(2e6) ≈ 1087.0^2 / 281.0
        @test_throws ArgumentError cost_rate(1e5; table = [(1e5, 64.0)])
        @test_throws ArgumentError cost_rate(1e5; table = [(1e6, 1087.0), (1e5, 64.0)])
        @test_throws ArgumentError cost_rate(0)
        # The configuration's form of a table
        table = [[1e5, 60.0], [1e6, 1000.0]]
        @test cost_rate(1e6; table = table) == 1000.0
        @test cost_rate(sqrt(1e5 * 1e6); table = table) ≈ sqrt(60.0 * 1000.0)
    end

    function cost_cfg(;
        wall_budget = 0.0,
        stop_margin = 120.0,
        checkpoint_keep = 0,
        budget_hours = 0.0,
    )
        return Nbody6Config(
            InstallConfig(),
            BuildConfig(),
            SimulationConfig(; wall_budget, stop_margin, checkpoint_keep),
            PostprocessConfig(),
            VisualizationConfig(),
            MergerPipelineConfig(),
            ClusterConfig(),
            CostConfig(; budget_hours),
        )
    end
    dump_size(n, nnbopt) = Nbody6Dynamics._checkpoint_dump_bytes(n, nnbopt)

    @testset "Estimate arithmetic" begin
        args = (; n = 5e5, span_nb = 100, nnbopt = 300, dtadj = 1.0, deltat = 10.0)
        e = estimate_cost(cost_cfg(; wall_budget = 36 * 3600 - 1800, stop_margin = 600); args...)
        @test e isa CostEstimate
        @test e.n == 500_000
        @test e.rate_s_per_nb == 281.0
        @test e.wall_s == 28100.0
        @test e.segments == 1
        @test startswith(e.basis, "interpolated between N = ")
        e = estimate_cost(cost_cfg(; wall_budget = 7200, stop_margin = 600); args...)
        @test e.segments == 5
        @test estimate_cost(cost_cfg(); args...).segments == 1

        # Dumps: one per adjustment, or the number kept
        e = estimate_cost(cost_cfg(; checkpoint_keep = 0); args...)
        @test e.dump_bytes == 101 * dump_size(500_000, 300)
        e = estimate_cost(cost_cfg(; checkpoint_keep = 2); args...)
        @test e.dump_bytes == 2 * dump_size(500_000, 300)
        @test e.snapshot_bytes == 40 * 500_000 * 11
        @test e.text_bytes == 130 * 500_000 * 11
        @test e.disk_bytes == e.dump_bytes + e.snapshot_bytes + e.text_bytes
        # Without intervals: one terminal dump, no snapshots
        e = estimate_cost(cost_cfg(); n = 5e5, span_nb = 100, nnbopt = 300)
        @test e.dump_bytes == dump_size(500_000, 300)
        @test e.snapshot_bytes == 0 && e.text_bytes == 0

        # Two-cluster configuration
        e = estimate_cost(cost_cfg(); args..., merger = true)
        @test e.rate_s_per_nb == 562.0
        @test e.wall_s == 56200.0
        @test occursin("merger factor", e.basis)

        # Budget
        @test estimate_cost(cost_cfg(; budget_hours = 1); args...).exceeds
        @test !estimate_cost(cost_cfg(; budget_hours = 10); args...).exceeds
        @test !estimate_cost(cost_cfg(); args...).exceeds

        @test_throws ArgumentError estimate_cost(cost_cfg(); n = 0, span_nb = 1)
        @test_throws ArgumentError estimate_cost(cost_cfg(); n = 1000, span_nb = 0)

        # Four-line summary
        e = estimate_cost(cost_cfg(; budget_hours = 1); args...)
        text = sprint(show, MIME"text/plain"(), e)
        @test length(split(text, '\n')) == 4
        @test occursin("exceeds the budget", text)

        fields = Nbody6Dynamics._cost_fields(e)
        @test fields["cost_wall_s"] == 28100.0
        @test fields["cost_segments"] == 1
        @test fields["cost_disk_bytes"] == e.disk_bytes
        @test fields["cost_rate_s_per_nb"] == 281.0
        @test fields["cost_basis"] == e.basis
    end

    @testset "Estimate from an input file" begin
        dir = mktempdir()
        inp = joinpath(dir, "fixture.inp")
        write(
            inp,
            """
            ! Fixture: N = 1 000 000 bodies, DTADJ = 7.0, TCRIT = 99
            &INNBODY6
            KSTART=1, TCOMP=100.0, TCRTP0=1.0E5, isernb=40, iserreg=40, iserks=640 /

            &ININPUT
            N=20000, NFIX=1, NCRIT=10, NRAND=10000, NNBOPT=300, NRUN=1, NCOMM=10,
            ETAI=0.02, ETAR=0.02, RS0=0.19, DTADJ=0.5, DELTAT=1.0, TCRIT=5.0D0,
            QE=1.0E-5, Level='C' /
            """,
        )
        @test Nbody6Dynamics._input_value(inp, "DTADJ") == 0.5
        @test Nbody6Dynamics._input_value(inp, "TCRIT") == 5.0
        @test Nbody6Dynamics._input_value(inp, "DTPLOT") === nothing

        e = estimate_cost(cost_cfg(), inp)
        @test e.n == 20_000
        @test e.span_nb == 5.0
        @test e.rate_s_per_nb == cost_rate(20_000)
        @test e.dump_bytes == 11 * dump_size(20_000, 300)
        @test e.snapshot_bytes == 40 * 20_000 * 6
        e = estimate_cost(cost_cfg(), inp; span_nb = 2.5)
        @test e.span_nb == 2.5
        @test e.dump_bytes == 6 * dump_size(20_000, 300)
        @test e.snapshot_bytes == 40 * 20_000 * 3

        bad = joinpath(dir, "no_nnbopt.inp")
        write(bad, "&ININPUT\nN=20000, TCRIT=5.0 /\n")
        @test_throws ArgumentError estimate_cost(cost_cfg(), bad)
    end

    @testset "Launch-time check" begin
        # Stand-in for the engine, as in the wall-budget tests of test_run.jl.
        fake_engine = raw"""
        #!/bin/bash
        tcrit=${FAKE_TCRIT:-1000}
        t=0
        while [ "$t" -le "$tcrit" ]; do
            printf ' ADJUST:  TIME    %d.00000E+00  T[Myr]   0.100E+01  Q   0.500E+00  DE   0.000E+00 DELTA   0.000E+00 DETOT   0.000E+00 E  -2.500E-01\n' "$t"
            if [ -z "$FAKE_IGNORE_STOP" ] && [ -f STOP ]; then
                printf '\n         TERMINATION BY MANUAL INTERVENTION\n'
                printf '  MYDUMP    %d.0000000000000000      %d.0                         201 comm.1_%d.0                \n' "$t" "$t" "$t"
                : > "comm.1_$t.0"
                printf '\n\n         COMMON SAVED AT TOFF/TIME/TTOT =  0.00000000E+00  TCOMP =  0.0\n'
                exit 0
            fi
            sleep 0.2
            t=$((t + 1))
        done
        printf '\n         END RUN    TIME[Myr] =    1.00  TOFF/TIME/TTOT=      0.00000000      %d.00000000      %d.00000000  CPUTOT =    0.0  ERRTOT = 0.00000D+00  DETOT = 0.00000D+00\n' "$tcrit" "$tcrit"
        """
        base = mktempdir()
        engine_path = joinpath(base, "engine", "build", "nbody6++")
        mkpath(dirname(engine_path))
        # `#!/bin/bash` at column 1, whether or not the literal was dedented.
        write(engine_path, replace(fake_engine, r"^ {4}"m => ""))
        chmod(engine_path, 0o755)
        write(joinpath(base, "in.inp"), "N=1000, NNBOPT=50, TCRIT=2.0, DTADJ=1.0, DELTAT=1.0\n")
        cfg_path = joinpath(base, "c.toml")
        function fake_cfg(; enabled = true, budget_hours = 0.0, on_exceed = "warn")
            write(
                cfg_path,
                """
                [install]
                enabled = false
                install_dir = "engine"

                [build]
                enable_mpi = false
                enable_gpu = false
                enable_hdf5 = false

                [simulation]
                input_file = "in.inp"
                runs_dir = "runs"
                binary_name = "nbody6++"
                omp_threads = 1
                monitor = false
                telemetry_interval = 0.0
                startup_timeout = 0.0
                exit_grace = 0.0

                [postprocess]
                enabled = false

                [visualization]
                enabled = false

                [cost]
                enabled = $enabled
                budget_hours = $(Float64(budget_hours))
                on_exceed = "$on_exceed"
                """,
            )
            return load_config(cfg_path)
        end
        launch(cfg, id) = withenv("FAKE_TCRIT" => "2", "FAKE_IGNORE_STOP" => nothing) do
            run_simulation(cfg; base_dir = base, run_id = id)
        end
        last_segment(run_dir) =
            Nbody6Dynamics.TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))["segments"][end]
        cost_keys =
            ("cost_wall_s", "cost_segments", "cost_disk_bytes", "cost_rate_s_per_nb", "cost_basis")

        # Refused before anything is launched
        refused = try
            launch(fake_cfg(; budget_hours = 1e-9, on_exceed = "refuse"), "refuse")
        catch err
            err
        end
        @test refused isa ArgumentError
        @test occursin("cost.budget_hours", sprint(showerror, refused))
        @test !isfile(joinpath(base, "runs", "refuse", "output", "_launch.sh"))

        # Warned, run to its end, estimate in the segment record
        warn_dir = @test_logs (:warn, r"exceeds the budget") match_mode = :any launch(
            fake_cfg(; budget_hours = 1e-9, on_exceed = "warn"),
            "warn",
        )
        seg = last_segment(warn_dir)
        @test seg["status"] == "completed"
        @test all(k -> haskey(seg, k), cost_keys)
        @test seg["cost_segments"] == 1
        @test startswith(seg["cost_basis"], "extrapolated below the table")
        @test seg["cost_rate_s_per_nb"] == round(cost_rate(1000); digits = 2)

        # Disabled: no estimate, no keys
        off_dir = launch(fake_cfg(; enabled = false, budget_hours = 1e-9), "off")
        seg = last_segment(off_dir)
        @test seg["status"] == "completed"
        @test !any(k -> haskey(seg, k), cost_keys)
    end
end
