# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Shipped GPU-host configurations" begin
    gdir = joinpath(@__DIR__, "..", "input_files", "gpu")
    gcfg = load_config(joinpath(gdir, "gpu_pipeline.toml"))
    @test gcfg.install.enabled && endswith(gcfg.install.install_dir, "Nbody6PPGPU-beijing-gpu")
    @test gcfg.build.enable_gpu && !gcfg.build.enable_hdf5 && isempty(gcfg.build.cuda_arch)
    @test gcfg.simulation.gpu_list == [0] && gcfg.simulation.omp_threads == 8
    @test gcfg.merger.enabled && gcfg.merger.config_file == "merger_50k.toml"
    ccfg = load_config(joinpath(gdir, "cpu_pipeline.toml"))
    @test ccfg.install.enabled && !ccfg.build.enable_gpu && isempty(ccfg.simulation.gpu_list)
    @test endswith(ccfg.install.install_dir, "Nbody6PPGPU-beijing")
    m = load_merger_config(joinpath(gdir, "merger_50k.toml"))
    @test length(m.clusters) == 2 && all(c -> c.N == 25000, m.clusters)
    @test m.nbody6.qe == 0.01 && m.output.tcrit_myr == 50.0
    @test m.orbit.apocentre == 10.0
end

# =====================================================================
# Engine-dependent tests: opt in with NBODY6_BINARY_TESTS=1. Without
# NBODY6_BACKEND_ROOT the backend is cloned and built in a temporary
# directory (what the Backend workflow does); with it, an existing
# build under <root>/backend/Nbody6PPGPU-beijing is used.
if get(ENV, "NBODY6_BINARY_TESTS", "0") == "1"
    @testset "Backend build and run (binary-gated)" begin
        work = mktempdir()
        root = get(ENV, "NBODY6_BACKEND_ROOT", "")
        artifacts = get(ENV, "NBODY6_BINARY_TEST_ARTIFACTS", "")
        isempty(artifacts) || mkpath(artifacts)
        function _keep(run_dir, tag)
            isempty(artifacts) && return
            for f in ("RUN_INFO.toml", "telemetry.csv", "telemetry_2.csv", "nbody6dynamics.log")
                src = joinpath(run_dir, f)
                isfile(src) && cp(src, joinpath(artifacts, "$(tag)_$(f)"); force = true)
            end
        end
        base_toml(runs_dir, extra_sim = "", merger = "enabled = false\nconfig_file = \"\"") = """
                                                                                  [install]
                                                                                  enabled = $(isempty(root))
                                                                                  install_dir = "backend/Nbody6PPGPU-beijing"
                                                                                  reinstall = false
                                                                                  clean_build = false

                                                                                  [build]
                                                                                  configure_flags = ["--enable-mcmodel=large", "--with-par=b1m"]
                                                                                  enable_mpi = false
                                                                                  enable_hdf5 = false
                                                                                  enable_gpu = false
                                                                                  nproc = 0

                                                                                  [simulation]
                                                                                  run_test = true
                                                                                  input_file = "../../input_files/N1k_quick.inp"
                                                                                  runs_dir = "$(runs_dir)"
                                                                                  binary_name = "nbody6++"
                                                                                  mpi_ranks = 1
                                                                                  omp_threads = 2
                                                                                  run_id_prefix = "bt"
                                                                                  monitor = false
                                                                                  telemetry_interval = 1.0
                                                                                  startup_timeout = 300.0
                                                                                  $(extra_sim)

                                                                                  [postprocess]
                                                                                  enabled = true

                                                                                  [visualization]
                                                                                  enabled = false

                                                                                  [merger]
                                                                                  $(merger)
                                                                                  """
        base = isempty(root) ? work : root
        if isempty(root)
            cp(joinpath(@__DIR__, "..", "input_files"), joinpath(work, "input_files"))
        end
        runs_dir = joinpath(work, "runs")

        # 1. Build (temporary tree only) and a single-cluster run
        cfg_path = joinpath(work, "single.toml")
        write(cfg_path, base_toml(runs_dir))
        cfg = load_config(cfg_path)
        t0 = time()
        results = run_pipeline(cfg; base_dir = base)
        @info "binary test: single run + build in $(round(time() - t0; digits = 1)) s"
        run_dir = Nbody6Dynamics._find_latest_run(cfg, base)
        @test isdir(run_dir)
        info = Nbody6Dynamics.TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))
        @test info["run"]["omp_threads"] == 2 && info["run"]["omp_threads_reported"] == 2
        @test info["run"]["segments"] == 1 && info["segments"][1]["exit_status"] == 0
        @test info["telemetry"]["cpu_user_s"] > 0 && info["telemetry"]["samples"] ≥ 1
        @test haskey(info["telemetry"], "backend_timing") &&
              info["telemetry"]["backend_timing"]["total"] > 0
        @test isfile(joinpath(run_dir, "telemetry.csv"))
        @test haskey(results, :snapshots) && length(results[:snapshots]) ≥ 2
        @test haskey(results, :diagnostics) && length(results[:diagnostics].adjust) ≥ 2
        _keep(run_dir, "single")

        # 2. Merger demo with frequent dumps, then a restart
        mtoml = read(joinpath(@__DIR__, "..", "input_files", "merger_demo_small.toml"), String)
        # Fixed seed for reproducibility of the engine-dependent runs
        mtoml = replace(mtoml, r"tcrit = [0-9.]+" => "tcrit = 2.0")
        mtoml =
            replace(mtoml, "[merger]\n" => "[merger]\nseed = 11\n"; count = 1) *
            "\n[merger.nbody6]\nncomm = 2\n"
        mpath = joinpath(work, "merger_demo.toml")
        write(mpath, mtoml)
        cfg_m_path = joinpath(work, "merger.toml")
        write(cfg_m_path, base_toml(runs_dir, "", "enabled = true\nconfig_file = \"$(mpath)\""))
        cfg_m = load_config(cfg_m_path)
        # the backend is built now: never rebuild
        cfg_m = Nbody6Config(
            InstallConfig(; enabled = false, install_dir = cfg_m.install.install_dir),
            cfg_m.build,
            cfg_m.simulation,
            cfg_m.postprocess,
            cfg_m.visualization,
            cfg_m.merger,
        )
        run_pipeline(cfg_m; base_dir = base)
        mrun = Nbody6Dynamics._find_latest_run(cfg_m, base)
        out = joinpath(mrun, "output")
        adjust_times(p) = [
            parse(Float64, match(r"TIME\s+([0-9.E+-]+)", l).captures[1]) for
            l in eachline(p) if startswith(lstrip(l), "ADJUST:")
        ]
        t_first = adjust_times(joinpath(out, "out1000"))
        @test maximum(t_first) ≈ 2.0
        @test Nbody6Dynamics._latest_dump(out) !== nothing
        restart_simulation(mrun; tcrit_extra = 1.0, base_dir = base)
        t_all = adjust_times(joinpath(out, "out1000"))
        @test maximum(t_all) ≈ 3.0 && length(t_all) > length(t_first)
        info_m = Nbody6Dynamics.TOML.parsefile(joinpath(mrun, "RUN_INFO.toml"))
        @test info_m["run"]["segments"] == 2 && info_m["segments"][2]["kind"] == "restart"
        @test info_m["run"]["input_file"] == "merger.inp"
        @test isfile(joinpath(mrun, "telemetry_2.csv"))
        @test any(startswith("conf.3_3"), readdir(out))
        _keep(mrun, "merger_restart")

        # 3. Point-mass tidal field
        ttoml = replace(
            mtoml,
            "[merger.nbody6]\nncomm = 2\n" => "[merger.nbody6]\nqe = 0.05\n\n[merger.tidal]\nkz14 = 2\ngmg = 1.0e11\nrg0 = 8.5\n",
        )
        ttoml = replace(ttoml, r"tcrit = [0-9.]+" => "tcrit = 1.0")
        tpath = joinpath(work, "merger_tidal.toml")
        write(tpath, ttoml)
        cfg_t_path = joinpath(work, "tidal.toml")
        write(cfg_t_path, base_toml(runs_dir, "", "enabled = true\nconfig_file = \"$(tpath)\""))
        cfg_t = load_config(cfg_t_path)
        cfg_t = Nbody6Config(
            InstallConfig(; enabled = false, install_dir = cfg_t.install.install_dir),
            cfg_t.build,
            cfg_t.simulation,
            cfg_t.postprocess,
            cfg_t.visualization,
            cfg_t.merger,
        )
        run_pipeline(cfg_t; base_dir = base)
        trun = Nbody6Dynamics._find_latest_run(cfg_t, base)
        tout = read(joinpath(trun, "output", "out1000"), String)
        @test occursin("POINT-MASS MODEL", tout)
        @test occursin(
            "KZ(11:20)=0 1 0 2 0 0 0 0 3 0",
            read(joinpath(trun, "output", "merger.inp"), String),
        )
        @test Nbody6Dynamics.TOML.parsefile(joinpath(trun, "RUN_INFO.toml"))["segments"][1]["exit_status"] ==
              0
        _keep(trun, "merger_tidal")

        # 4. Two-point sweep through the worker processes
        swork = joinpath(work, "sweep")
        mkpath(swork)
        stoml = replace(mtoml, r"tcrit = [0-9.]+" => "tcrit = 0.5")
        stoml = replace(stoml, r"N = 1000" => "N = 300")
        write(joinpath(swork, "merger.toml"), stoml)
        write(
            joinpath(swork, "config.toml"),
            replace(
                base_toml(runs_dir),
                "install_dir = \"backend/Nbody6PPGPU-beijing\"" => "install_dir = \"$(joinpath(base, "backend", "Nbody6PPGPU-beijing"))\"",
            ),
        )
        write(
            joinpath(swork, "sweep.toml"),
            """
[sweep]
name = "gated"
pipeline_config = "config.toml"
merger_config = "merger.toml"
seeds = [11]
concurrency = 2
omp_threads = 2
runs_dir = "runs"
poll_interval = 1.0

[sweep.grid]
"merger.orbit.eccentricity" = [0.0, 0.6]
""",
        )
        scfg = load_sweep_config(joinpath(swork, "sweep.toml"))
        t0 = time()
        sdir = run_sweep(scfg)
        @info "binary test: two-point sweep in $(round(time() - t0; digits = 1)) s"
        sidx = read_sweep_index(sdir)
        @test length(sidx["points"]) == 2
        @test all(p -> p["status"] == "done" && p["exit_status"] == 0, sidx["points"])
        for p in sidx["points"]
            @test isfile(joinpath(p["dir"], "run", "RUN_INFO.toml"))
            @test isfile(joinpath(p["dir"], "run", "output", "lagr.7"))
            @test isfile(joinpath(p["dir"], "sweep_point.log"))
        end
        scols, srows = sweep_summary(sdir)
        @test all(r -> r["n_final"] > 0 && r["t_final_myr"] > 0, srows)
        @test isfile(joinpath(sdir, "sweep_summary.csv"))
        sfigs = sweep_figures(sdir, sweep_visualization(scfg, sdir))
        @test all(isfile, sfigs)
    end
end

# GPU-dependent tests: opt in with NBODY6_GPU_TESTS=1 on a host with an
# NVIDIA device, nvidia-smi and nvcc. The backend is cloned and built
# with CUDA in a temporary tree (a GPU build in the working tree would
# replace the CPU binary's objects), then the 1k input runs on the first
# device and, with two or more devices, on the first two.
if get(ENV, "NBODY6_GPU_TESTS", "0") == "1"
    @testset "Backend GPU build and run (GPU-gated)" begin
        caps = detect_compute_capabilities()
        @test !isempty(caps)
        archs = cuda_arch_from_compute_cap.(caps)
        supported = nvcc_supported_archs(detect_cuda_path())
        @test !isempty(supported) && all(a -> a in supported, archs)
        n_dev = count(
            !isempty,
            strip.(split(read(`nvidia-smi --query-gpu=index --format=csv,noheader`, String), '\n')),
        )
        work = mktempdir()
        cp(joinpath(@__DIR__, "..", "input_files"), joinpath(work, "input_files"))
        runs_dir = joinpath(work, "runs")
        gpu_toml(gpu_list, install_enabled) = """
[install]
enabled = $(install_enabled)
install_dir = "backend/Nbody6PPGPU-beijing"
reinstall = false
clean_build = false

[build]
configure_flags = ["--enable-mcmodel=large", "--with-par=b1m"]
enable_mpi = false
enable_hdf5 = false
enable_gpu = true
nproc = 0

[simulation]
run_test = true
input_file = "../../input_files/N1k_quick.inp"
runs_dir = "$(runs_dir)"
binary_name = "nbody6++"
mpi_ranks = 1
omp_threads = 4
gpu_list = $(gpu_list)
run_id_prefix = "gt"
monitor = false
telemetry_interval = 1.0
startup_timeout = 600.0

[postprocess]
enabled = true

[visualization]
enabled = false

[merger]
enabled = false
config_file = ""
"""
        # 1. Build with CUDA for the detected capabilities and run on device 0
        cfg_path = joinpath(work, "gpu1.toml")
        write(cfg_path, gpu_toml("[0]", true))
        cfg = load_config(cfg_path)
        t0 = time()
        results = run_pipeline(cfg; base_dir = work)
        @info "GPU test: build + single-device run in $(round(time() - t0; digits = 1)) s"
        run_dir = Nbody6Dynamics._find_latest_run(cfg, work)
        info = Nbody6Dynamics.TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))
        @test info["segments"][1]["exit_status"] == 0
        @test endswith(info["build"]["binary"], ".gpu")
        @test info["build"]["cuda_arch"] == archs
        @test !isempty(info["build"]["nvcc_release"])
        @test info["run"]["gpu_list"] == [0]
        @test length(info["run"]["gpu_devices"]) == 1
        @test occursin("GPU Reg.F", info["telemetry"]["force_kernel_gflops"]["kernel"])
        @test info["telemetry"]["force_kernel_gflops"]["mean"] > 0
        @test info["hardware"]["gpu"] != "unavailable"
        @test haskey(info["telemetry"], "mean_gpu_util_pct")
        diag = results[:diagnostics]
        @test length(diag.adjust) ≥ 2
        @test maximum(abs(a.de_rel) for a in diag.adjust) < 1e-2
        @test isfile(joinpath(run_dir, "output", "BUILD_INFO.toml"))

        # 2. Two devices in one process, same binary
        if n_dev ≥ 2
            cfg2_path = joinpath(work, "gpu2.toml")
            write(cfg2_path, gpu_toml("[0, 1]", false))
            cfg2 = load_config(cfg2_path)
            run_pipeline(cfg2; base_dir = work)
            run_dir2 = Nbody6Dynamics._find_latest_run(cfg2, work)
            info2 = Nbody6Dynamics.TOML.parsefile(joinpath(run_dir2, "RUN_INFO.toml"))
            @test info2["segments"][1]["exit_status"] == 0
            @test info2["run"]["gpu_list"] == [0, 1]
            @test length(info2["run"]["gpu_devices"]) == 2
        else
            @info "GPU test: one device visible; the two-device case is skipped"
        end
    end
end

# =====================================================================
