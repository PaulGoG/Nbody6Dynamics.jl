# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Hardware telemetry" begin
    # /proc stat parser: command names may contain spaces and parentheses
    rec =
        "4242 (my (odd) proc) S 17 4242 4242 0 -1 4194560 100 0 0 0 350 25 0 0 20 0 8 0 100 " *
        "12345678 2048 18446744073709551615 1 1 0 0 0 0 0 0 0 0 0 0 0 0 17 3 0 0 0 0 0 0 0 0 0 0 0 0 0"
    @test Nbody6Dynamics._parse_proc_stat(rec) == (17, 350, 25)
    @test_throws ArgumentError Nbody6Dynamics._parse_proc_stat("garbage without parens")
    @test_throws ArgumentError Nbody6Dynamics._parse_proc_stat("1 (x) S 2 3")

    # nvidia-smi aggregation: means, sums, max; [N/A] → NaN
    two = "45, 10, 1234, 120.50, 60\n55, 30, 766, 79.50, 70\n"
    u, mu, mem, pw, tc = Nbody6Dynamics._parse_nvidia_smi(two)
    @test u ≈ 50.0 && mu ≈ 20.0 && mem ≈ 2000.0 && pw ≈ 200.0 && tc == 70.0
    na = Nbody6Dynamics._parse_nvidia_smi("12, 3, 500, [N/A], 41\n")
    @test na[1] == 12.0 && isnan(na[4]) && na[5] == 41.0
    @test all(isnan, Nbody6Dynamics._parse_nvidia_smi(""))
    @test_throws ArgumentError Nbody6Dynamics._parse_nvidia_smi("1, 2, 3\n")

    # Exact child accounting is monotone across a child launch
    cpu0 = Nbody6Dynamics._children_cpu_times()
    wait(run(`sleep 0.1`; wait = false))
    cpu1 = Nbody6Dynamics._children_cpu_times()
    if Sys.islinux()
        @test cpu1[1] ≥ cpu0[1] && cpu1[2] ≥ cpu0[2]
    else
        @test all(isnan, cpu1)
    end

    # Sampler on a live child process
    tele_dir = mktempdir()
    t0 = time()
    child = run(`sleep 1.2`; wait = false)
    @test_throws ArgumentError Nbody6Dynamics._start_telemetry(
        tele_dir,
        getpid(child),
        t0;
        interval = 0.0,
        gpu_probe = false,
    )
    mon = Nbody6Dynamics._start_telemetry(
        tele_dir,
        getpid(child),
        t0;
        interval = 0.2,
        gpu_probe = false,
    )
    wait(child)
    summary = Nbody6Dynamics._finish_telemetry(
        mon,
        cpu0,
        Nbody6Dynamics._children_cpu_times(),
        time() - t0,
        2,
    )
    @test summary["samples"] == length(mon.samples) ≥ 3
    @test summary["sampling_interval_s"] == 0.2
    @test summary["threads_total"] == 2
    csv_lines = readlines(joinpath(tele_dir, "telemetry.csv"))
    @test csv_lines[1] == join(string.(fieldnames(Nbody6Dynamics.TelemetrySample)), ",")
    @test length(csv_lines) == summary["samples"] + 1
    @test all(l -> count(==(','), l) == fieldcount(Nbody6Dynamics.TelemetrySample) - 1, csv_lines)
    if Sys.islinux()
        @test mon.samples[1].n_processes ≥ 1
        @test summary["peak_rss_mib"] > 0
        @test haskey(summary, "cpu_user_s") && haskey(summary, "cpu_efficiency")
        @test summary["cpu_efficiency"] ≥ 0
        @test haskey(summary, "peak_load_1min")
    end
    @test !haskey(summary, "mean_gpu_util_pct")   # no GPU probe requested

    # Backend-reported performance: stdout timing table and stderr Gflops
    @test Nbody6Dynamics._timing_key("Reg.GPU.S") == "reg_gpu_s"
    @test Nbody6Dynamics._timing_key("KS.Init.B") == "ks_init_b"
    @test Nbody6Dynamics._timing_key("Total") == "total"
    perf_dir = mktempdir()
    out_perf = joinpath(perf_dir, "out1000")
    err_perf = joinpath(perf_dir, "err1000")
    header = "  rank   PE   N        Total     Init.    Intgrt      Reg.      Irr.        KS    Adjust       OUT Reg.GPU.S    xtsub1   itides3"
    write(
        out_perf,
        "noise\n" *
        header *
        "\n" *
        "   0  0    1000      1.00000      0.10      0.80      0.30      0.40      0.05      0.02      0.01      0.00  0.00000E+00         0\n" *
        "ADJUST: ...\n" *
        header *
        "\n" *
        "   0  0     998      3.12553      0.19      2.86      0.83      1.14      0.15      0.11      0.00      0.83  0.00000E+00         2\n" *
        "   short line\n",
    )
    timing = Nbody6Dynamics._backend_timing_table(out_perf)
    @test timing["n_particles"] == 998        # last table wins
    @test timing["total"] ≈ 3.12553
    @test timing["reg"] ≈ 0.83 && timing["irr"] ≈ 1.14 && timing["ks"] ≈ 0.15
    @test timing["reg_gpu_s"] ≈ 0.83 && timing["xtsub1"] == 0.0 && timing["itides3"] == 2
    @test !haskey(timing, "rank") && !haskey(timing, "pe")
    write(out_perf, "no table\n")
    @test Nbody6Dynamics._backend_timing_table(out_perf) === nothing
    @test Nbody6Dynamics._backend_timing_table(joinpath(perf_dir, "absent")) === nothing
    write(
        err_perf,
        "[R.0 AVX Pot.A] Ni 1000  NTOT 1000  pot(s) 0.000796\n" *
        "[R.0 AVX Reg.F ] Nsend 1  Ngrav 1  <Ni> 1000   send(s) 0.000033 grav(s) 0.002927  Perf.(Gflops) 20.500020\n" *
        "# Open AVX regular force - rank: 0; threads: 4\n" *
        "[R.0 AVX Reg.F ] Nsend 3  Ngrav 3  <Ni> 900   send(s) 0.000040 grav(s) 0.003000  Perf.(Gflops) 30.5\n",
    )
    gf = Nbody6Dynamics._force_kernel_gflops(err_perf)
    @test gf["samples"] == 2 && gf["mean"] ≈ 25.5 && gf["peak"] ≈ 30.5
    write(err_perf, "nothing here\n")
    @test Nbody6Dynamics._force_kernel_gflops(err_perf) === nothing
    perf = Nbody6Dynamics._backend_performance(joinpath(perf_dir, "absent"), err_perf)
    @test isempty(perf)

    # Disabled sampler still carries the exact accounting
    none = Nbody6Dynamics._finish_telemetry(nothing, cpu0, cpu1, 10.0, 4)
    @test none["samples"] == 0 && none["sampling_interval_s"] == 0.0
    @test none["threads_total"] == 4

    # Run summary carries the thread layout and the telemetry table
    run_dir = mktempdir()
    out_dir = joinpath(run_dir, "output")
    mkpath(out_dir)
    write(joinpath(out_dir, "out1000"), " RANK:  0  OpenMP Number of Threads:  8\n")
    cfg_t = Nbody6Config(
        InstallConfig(),
        BuildConfig(),
        SimulationConfig(; omp_threads = 8, mpi_ranks = 1),
        PostprocessConfig(),
        VisualizationConfig(),
        MergerPipelineConfig(),
    )
    Nbody6Dynamics._write_run_summary(
        cfg_t,
        run_dir,
        "testrun_tele",
        joinpath(out_dir, "out1000"),
        out_dir,
        3.0;
        telemetry = summary,
    )
    info = Nbody6Dynamics.TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))
    @test info["run"]["omp_threads"] == 8
    @test info["run"]["omp_threads_reported"] == 8
    @test info["run"]["mpi_ranks"] == 1
    @test info["telemetry"]["samples"] == summary["samples"]
    @test info["telemetry"]["threads_total"] == 2
end

# =====================================================================

@testset "Telemetry readers and figure" begin
    dir = mktempdir()
    header = join(string.(fieldnames(TelemetrySample)), ",")
    row(
        i;
        gpu = "NaN",
    ) = "$(5.0 * i),3,$(200 + i),$(210 + i),$(10.0 * i),$(2.0 + 0.1 * i),1.5,$gpu,$gpu,$gpu,$gpu,$gpu"
    write(joinpath(dir, "telemetry.csv"), header * "\n" * join([row(i) for i in 1:6], "\n") * "\n")
    s = read_telemetry(joinpath(dir, "telemetry.csv"))
    @test length(s) == 6 && s[1].elapsed_s == 5.0 && s[6].rss_mib == 206.0 && s[1].n_processes == 3
    @test isnan(s[1].gpu_util_pct) && s[3].cores_busy == 2.3
    # Columns matched by name, in any order; a missing column or a ragged row is an error
    cols = split(header, ',')
    perm = reverse(cols)
    write(
        joinpath(dir, "perm.csv"),
        join(perm, ",") * "\n" * join(reverse(split(row(2), ',')), ",") * "\n",
    )
    @test read_telemetry(joinpath(dir, "perm.csv"))[1].elapsed_s == 10.0
    write(joinpath(dir, "short.csv"), "elapsed_s,rss_mib\n1,2\n")
    @test_throws ArgumentError read_telemetry(joinpath(dir, "short.csv"))
    write(joinpath(dir, "ragged.csv"), header * "\n1,2,3\n")
    @test_throws ArgumentError read_telemetry(joinpath(dir, "ragged.csv"))
    @test read_telemetry(joinpath(dir, "empty.csv") |> p -> (write(p, ""); p)) == TelemetrySample[]
    # Segments are concatenated with cumulative offsets
    write(
        joinpath(dir, "telemetry_2.csv"),
        header * "\n" * join([row(i) for i in 1:2], "\n") * "\n",
    )
    all_s = read_run_telemetry(dir)
    @test length(all_s) == 8 && all_s[7].elapsed_s == 30.0 + 5.0 && all_s[8].elapsed_s == 40.0
    @test read_run_telemetry(joinpath(dir, "nowhere")) == TelemetrySample[]
    # Figures: two panels without GPU samples, three with them
    vis = VisualizationConfig(; enabled = true, format = "png", dpi = 100, output_dir = dir)
    p2 = plot_telemetry(all_s, vis; filename = "tel_cpu")
    @test p2 !== nothing && isfile(p2)
    gpu_rows = join([row(i; gpu = "$(50 + i)") for i in 1:6], "\n")
    write(joinpath(dir, "gpu.csv"), header * "\n" * gpu_rows * "\n")
    p3 = plot_telemetry(read_telemetry(joinpath(dir, "gpu.csv")), vis; filename = "tel_gpu")
    @test p3 !== nothing && isfile(p3)
    @test plot_telemetry(s[1:1], vis; filename = "tel_one") === nothing
    # Time axis units follow the span
    @test MakieExt._telemetry_time_axis([0.0, 60.0])[1] == [0.0, 60.0]
    @test MakieExt._telemetry_time_axis([0.0, 600.0])[1] == [0.0, 10.0]
    @test MakieExt._telemetry_time_axis([0.0, 7200.0])[1] == [0.0, 2.0]
end

# =====================================================================
