using Test
using Nbody6Dynamics

# Temporary directory for test artifacts
const TESTDIR = mktempdir()

# ---------------------------------------------------------------------------
# Helper: write a Fortran binary record (for test data generation)
# ---------------------------------------------------------------------------
function _write_fortran_record(io::IO, data::Vector{T}) where {T}
    bytes = reinterpret(UInt8, data)
    marker = Int32(length(bytes))
    write(io, marker)
    write(io, bytes)
    write(io, marker)
end

function _write_fortran_record(io::IO, data::Vector{UInt8})
    marker = Int32(length(data))
    write(io, marker)
    write(io, data)
    write(io, marker)
end

@testset "Nbody6Dynamics.jl" begin

    # =====================================================================
    @testset "Configuration" begin
        # Write a minimal config
        cfg_path = joinpath(TESTDIR, "test_config.toml")
        write(
            cfg_path,
            """
[install]
enabled = false
install_dir = "test-nbody"

[build]
enable_mpi = false
enable_hdf5 = false
enable_gpu = false
cuda_path = "/usr/local/cuda"

[simulation]
run_test = false
runs_dir = "test-runs"
run_id_prefix = "test"

[postprocess]
enabled = false
read_escapers = true
escapers_file = "esc.11"
read_stellar_evo = true
stellar_evo_pattern = "sev.83_*"

[visualization]
enabled = false
dpi = 150
figsize = [10, 8]
""",
        )

        cfg = load_config(cfg_path)

        @test cfg.install.enabled == false
        @test cfg.install.install_dir == "test-nbody"
        @test cfg.build.enable_mpi == false
        @test cfg.build.cuda_path == "/usr/local/cuda"
        @test cfg.simulation.run_test == false
        @test cfg.simulation.run_id_prefix == "test"
        @test cfg.postprocess.read_escapers == true
        @test cfg.postprocess.escapers_file == "esc.11"
        @test cfg.postprocess.read_stellar_evo == true
        @test cfg.postprocess.stellar_evo_pattern == "sev.83_*"
        @test cfg.postprocess.read_binary_evo == true
        @test cfg.postprocess.binary_evo_pattern == "bev.82_*"
        @test cfg.visualization.dpi == 150
        @test cfg.visualization.figsize == (10, 8)
        @test cfg.install.source_url == "https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing.git"
    end

    # =====================================================================
    @testset "Config round-trip (save/load)" begin
        cfg_path = joinpath(TESTDIR, "test_config_rt.toml")
        write(
            cfg_path,
            """
[install]
enabled = true
ref = "v2026.07"

[build]
enable_gpu = true
cuda_path = "/opt/cuda"

[simulation]
run_id_prefix = "bench"
omp_threads = 6
telemetry_interval = 2.5

[postprocess]
enabled = true

[visualization]
dpi = 300
figsize = [12, 9]
""",
        )

        cfg = load_config(cfg_path)
        save_path = joinpath(TESTDIR, "test_config_saved.toml")
        save_config(cfg, save_path)

        cfg2 = load_config(save_path)
        @test cfg2.install.ref == "v2026.07"
        @test InstallConfig().ref == "618d7a4"
        @test cfg2.build.enable_gpu == cfg.build.enable_gpu
        @test cfg2.build.cuda_path == cfg.build.cuda_path
        @test cfg2.simulation.run_id_prefix == cfg.simulation.run_id_prefix
        @test cfg2.simulation.omp_threads == 6
        @test cfg2.simulation.telemetry_interval == 2.5
        @test cfg2.visualization.dpi == cfg.visualization.dpi
        @test cfg2.visualization.figsize == cfg.visualization.figsize
    end

    # =====================================================================
    @testset "Config validation" begin
        # Valid config (defaults only) passes validation — smoke test
        ok_path = joinpath(TESTDIR, "val_ok.toml")
        write(
            ok_path,
            """
[simulation]
mpi_ranks = 1

[visualization]
format = "pdf"
""",
        )
        @test load_config(ok_path) isa Nbody6Config

        # Each entry: (label, TOML violating exactly one rule, expected
        # substring of the error message naming the offending section.key).
        bad_cases = [
            ("install_dir", "[install]\ninstall_dir = \"\"\n", "install.install_dir"),
            ("nproc", "[build]\nnproc = -1\n", "build.nproc"),
            ("mpi_ranks", "[simulation]\nmpi_ranks = 0\n", "simulation.mpi_ranks"),
            ("omp_threads", "[simulation]\nomp_threads = -1\n", "simulation.omp_threads"),
            (
                "telemetry_interval",
                "[simulation]\ntelemetry_interval = -0.5\n",
                "simulation.telemetry_interval",
            ),
            (
                "startup_timeout",
                "[simulation]\nstartup_timeout = -1.0\n",
                "simulation.startup_timeout",
            ),
            (
                "mpi_no_mpi_build",
                "[build]\nenable_mpi = false\n\n[simulation]\nmpi_ranks = 4\n",
                "requires build.enable_mpi = true",
            ),
            ("runs_dir", "[simulation]\nruns_dir = \"\"\n", "simulation.runs_dir"),
            ("run_id_prefix", "[simulation]\nrun_id_prefix = \"\"\n", "simulation.run_id_prefix"),
            ("input_file", "[simulation]\ninput_file = \"\"\n", "simulation.input_file"),
            (
                "snapshot_format_hdf5",
                "[postprocess]\nsnapshot_format = \"hdf5\"\n",
                "postprocess.snapshot_format",
            ),
            (
                "snapshot_format_other",
                "[postprocess]\nsnapshot_format = \"parquet\"\n",
                "postprocess.snapshot_format",
            ),
            (
                "snapshot_pattern",
                "[postprocess]\nsnapshot_pattern = \"\"\n",
                "postprocess.snapshot_pattern",
            ),
            (
                "stdout_file",
                "[postprocess]\nparse_stdout = true\nstdout_file = \"\"\n",
                "postprocess.stdout_file",
            ),
            (
                "lagr_file",
                "[postprocess]\nread_lagr = true\nlagr_file = \"\"\n",
                "postprocess.lagr_file",
            ),
            (
                "escapers_file",
                "[postprocess]\nread_escapers = true\nescapers_file = \"\"\n",
                "postprocess.escapers_file",
            ),
            (
                "stellar_evo_pattern",
                "[postprocess]\nread_stellar_evo = true\nstellar_evo_pattern = \"\"\n",
                "postprocess.stellar_evo_pattern",
            ),
            ("format", "[visualization]\nformat = \"gif\"\n", "visualization.format"),
            ("column", "[visualization]\ncolumn = \"triple\"\n", "visualization.column"),
            ("units", "[visualization]\nunits = \"cgs\"\n", "visualization.units"),
            ("dpi", "[visualization]\ndpi = 50\n", "visualization.dpi"),
            ("figsize", "[visualization]\nfigsize = [0.0, 6.0]\n", "visualization.figsize"),
            (
                "marker_budget",
                "[visualization.style]\nmarker_budget = 0.0\n",
                "visualization.style.marker_budget",
            ),
            (
                "marker_bounds",
                "[visualization.style]\nmarker_min = 10.0\nmarker_max = 4.0\n",
                "0 < marker_min ≤ marker_max",
            ),
            (
                "q_log_threshold",
                "[visualization.style]\nq_log_threshold = -1.0\n",
                "visualization.style.q_log_threshold",
            ),
            ("q_floor", "[visualization.style]\nq_floor = 1.5\n", "visualization.style.q_floor"),
            (
                "zoom_frac",
                "[visualization.style]\nzoom_frac = 0.0\n",
                "visualization.style.zoom_frac",
            ),
            ("anim_fps", "[visualization.style]\nanim_fps = -5\n", "visualization.style.anim_fps"),
            (
                "anim_target_seconds",
                "[visualization.style]\nanim_target_seconds = 0.0\n",
                "visualization.style.anim_target_seconds",
            ),
            (
                "merger_config_file",
                "[merger]\nenabled = true\nconfig_file = \"\"\n",
                "merger.config_file",
            ),
        ]
        for (label, body, expected) in bad_cases
            path = joinpath(TESTDIR, "val_bad_$label.toml")
            write(path, body)
            @test_throws expected load_config(path)
        end
    end

    # =====================================================================
    @testset "Platform detection" begin
        platform = detect_platform()
        @test platform isa Symbol
        @test platform in [:fedora, :ubuntu, :debian, :unknown]

        # check_command should find basic tools
        @test Nbody6Dynamics.check_command("ls") == true
        @test Nbody6Dynamics.check_command("nonexistent_tool_xyz") == false
    end

    # =====================================================================
    @testset "CUDA detection" begin
        # detect_cuda_path returns a string (may be empty if no CUDA)
        cuda = detect_cuda_path()
        @test cuda isa String

        # cuda_env_vars returns a dict
        env = Nbody6Dynamics.cuda_env_vars("/usr/local/cuda")
        @test haskey(env, "CUDA_HOME")
        @test env["CUDA_HOME"] == "/usr/local/cuda"
        @test occursin("/usr/local/cuda/bin", env["PATH"])

        # Empty path returns empty dict
        env_empty = Nbody6Dynamics.cuda_env_vars("")
        @test isempty(env_empty)
    end

    # =====================================================================
    @testset "Engine source tree validation" begin
        # A local repository standing in for the upstream engine: `configure`
        # is tracked there, so a checkout without it is incomplete.
        quiet(cmd) = run(pipeline(cmd; stdout = devnull, stderr = devnull))
        origin = mktempdir()
        write(joinpath(origin, "configure"), "#!/bin/sh\nexit 0\n")
        quiet(`git -C $origin init --quiet`)
        quiet(`git -C $origin add configure`)
        quiet(`git -C $origin -c user.email=t@t -c user.name=t commit --quiet -m init`)
        inst = Nbody6Dynamics.InstallConfig(; source_url = origin, ref = "")
        src = joinpath(mktempdir(), "engine")

        # Absent → cloned
        Nbody6Dynamics._ensure_source_tree(src, inst)
        @test Nbody6Dynamics._source_tree_ready(src)

        # Complete → left as it is, untracked work untouched
        marker = joinpath(src, "untracked.txt")
        write(marker, "keep me")
        Nbody6Dynamics._ensure_source_tree(src, inst)
        @test isfile(marker) && read(marker, String) == "keep me"

        # Truncated checkout (an interrupted clone) → restored from git,
        # untracked work kept
        rm(joinpath(src, "configure"))
        @test !Nbody6Dynamics._source_tree_ready(src)
        @test_logs (:warn, r"checkout is incomplete") match_mode = :any Nbody6Dynamics._ensure_source_tree(
            src,
            inst,
        )
        @test Nbody6Dynamics._source_tree_ready(src)
        @test isfile(marker)

        # No usable git state left → removed and cloned again
        rm(joinpath(src, "configure"))
        rm(joinpath(src, ".git"); recursive = true)
        @test_logs (:warn, r"cannot be restored") match_mode = :any Nbody6Dynamics._ensure_source_tree(
            src,
            inst,
        )
        @test Nbody6Dynamics._source_tree_ready(src)
        @test !isfile(marker)
    end

    @testset "GPU build target" begin
        # configure arguments: the engine's configure finds nvcc only on PATH
        # (its --with-cuda fallback reuses the cached PATH check), so the
        # build exports the toolkit's bin and names the prefix explicitly.
        build_of(body) =
            (p = joinpath(mktempdir(), "cfg.toml"); write(p, body); load_config(p).build)
        cpu_args = Nbody6Dynamics._configure_args(
            build_of("[build]\nenable_gpu = false\n"),
            "/usr/local/cuda",
        )
        @test "--disable-gpu" in cpu_args && "--disable-mpi" in cpu_args
        @test !any(startswith("--with-cuda"), cpu_args)
        @test Nbody6Dynamics._configure_args(
            build_of("[build]\nenable_gpu = true\nconfigure_flags = [\"--with-par=b1m\"]\n"),
            "/opt/cuda-13.1",
        ) == ["--with-par=b1m", "--disable-mpi", "--with-cuda=/opt/cuda-13.1"]
        @test Nbody6Dynamics._configure_args(
            build_of("[build]\nenable_gpu = true\nconfigure_flags = [\"--with-cuda=/x\"]\n"),
            "/opt/cuda-13.1",
        ) == ["--with-cuda=/x", "--disable-mpi"]
        @test !any(
            startswith("--with-cuda"),
            Nbody6Dynamics._configure_args(build_of("[build]\nenable_gpu = true\n"), ""),
        )

        # Architecture names and code-generation flags
        @test cuda_arch_from_compute_cap("9.0") == "sm_90"
        @test cuda_arch_from_compute_cap("12.0") == "sm_120"
        @test cuda_arch_from_compute_cap(" 8.6\n") == "sm_86"
        @test_throws ArgumentError cuda_arch_from_compute_cap("sm_90")
        @test_throws ArgumentError cuda_arch_from_compute_cap("")
        @test Nbody6Dynamics._parse_compute_caps("9.0\n9.0\n12.0\n[N/A]\n\n") == ["9.0", "12.0"]
        @test Nbody6Dynamics._parse_compute_caps("") == String[]
        @test detect_compute_capabilities() isa Vector{String}
        @test cuda_gencode_flags(String[]) == ""
        @test cuda_gencode_flags(["sm_90"]) ==
              "-gencode arch=compute_90,code=sm_90 -gencode arch=compute_90,code=compute_90"
        flags = cuda_gencode_flags(["sm_120", "sm_90"])
        @test occursin("code=sm_120", flags) && occursin("code=sm_90", flags)
        @test endswith(flags, "-gencode arch=compute_120,code=compute_120")
        @test_throws ArgumentError cuda_gencode_flags(["sm90"])
        @test Nbody6Dynamics._parse_nvcc_release(
            "Cuda compilation tools, release 12.8, V12.8.93",
        ) == "12.8"
        @test Nbody6Dynamics._parse_nvcc_release("no banner") == ""
        @test nvcc_release("/nonexistent/cuda") == ""

        # Toolkit support: nvcc lists what it can compile, nvidia-smi what is installed
        @test Nbody6Dynamics._parse_nvcc_arch_list(
            "compute_50\ncompute_90\ncompute_100\ncompute_120\ncompute_90\n",
        ) == ["sm_50", "sm_90", "sm_100", "sm_120"]
        @test Nbody6Dynamics._parse_nvcc_arch_list("") == String[]
        @test nvcc_supported_archs("/nonexistent/cuda") == String[]
        @test nvcc_supported_archs() isa Vector{String}
        @test Nbody6Dynamics._check_cuda_arch_support(["sm_120"], String[], "") === nothing
        @test Nbody6Dynamics._check_cuda_arch_support(["sm_90"], ["sm_80", "sm_90"], "12.4") ===
              nothing
        unsupported = try
            Nbody6Dynamics._check_cuda_arch_support(["sm_90", "sm_120"], ["sm_80", "sm_90"], "12.4")
        catch err
            err
        end
        @test unsupported isa ErrorException &&
              occursin("sm_120", unsupported.msg) &&
              occursin("12.4", unsupported.msg) &&
              !occursin("for sm_90", unsupported.msg)

        # Resolution: an explicit list wins; without a device the list stays empty
        @test resolve_cuda_arch(
            BuildConfig(; enable_gpu = true, cuda_arch = ["sm_90", "sm_120"]),
        ) == ["sm_90", "sm_120"]
        if isempty(detect_compute_capabilities())
            resolved = @test_logs (:warn, r"nvcc default target") resolve_cuda_arch(
                BuildConfig(; enable_gpu = true),
            )
            @test resolved == String[]
        end

        # Configuration: parsing, validation, round-trip
        dir = mktempdir()
        write_cfg(body) = (p = joinpath(dir, "c.toml"); write(p, body); p)
        cfg = load_config(
            write_cfg(
                "[build]\nenable_gpu = true\ncuda_arch = [\"sm_90\", \"sm_120\"]\n" *
                "[simulation]\ngpu_list = [0, 1]\n",
            ),
        )
        @test cfg.build.cuda_arch == ["sm_90", "sm_120"]
        @test cfg.simulation.gpu_list == [0, 1]
        @test BuildConfig().cuda_arch == String[] && SimulationConfig().gpu_list == Int[]
        rt = joinpath(dir, "rt.toml")
        save_config(cfg, rt)
        cfg2 = load_config(rt)
        @test cfg2.build.cuda_arch == ["sm_90", "sm_120"] && cfg2.simulation.gpu_list == [0, 1]
        bad_arch = try
            load_config(write_cfg("[build]\ncuda_arch = [\"sm90\"]\n"))
        catch err
            err
        end
        @test bad_arch isa ErrorException && occursin("build.cuda_arch", bad_arch.msg)
        @test_throws ErrorException load_config(
            write_cfg("[build]\ncuda_arch = [\"compute_90\"]\n"),
        )
        no_gpu = try
            load_config(write_cfg("[simulation]\ngpu_list = [0]\n"))
        catch err
            err
        end
        @test no_gpu isa ErrorException && occursin("build.enable_gpu", no_gpu.msg)
        gpu_on = "[build]\nenable_gpu = true\n[simulation]\n"
        @test_throws ErrorException load_config(write_cfg(gpu_on * "gpu_list = [-1]\n"))
        @test_throws ErrorException load_config(write_cfg(gpu_on * "gpu_list = [0, 0]\n"))
        @test_throws ErrorException load_config(write_cfg(gpu_on * "gpu_list = [0, 1, 2, 3, 4]\n"))
        @test load_config(write_cfg(gpu_on * "gpu_list = [3, 1, 2, 0]\n")).simulation.gpu_list ==
              [3, 1, 2, 0]

        # CUFLAGS override built from the configured Makefile
        mk = joinpath(dir, "Makefile")
        write(mk, "NVCC = nvcc\nCUFLAGS =  -O3 -D CUDA_5 -I ../extra_inc/cuda\nLIBS = -lm\n")
        helper = Nbody6Dynamics._CUDA_HELPER_DIR
        @test Nbody6Dynamics._cuflags_with_arch(mk, ["sm_90"]) ==
              "-I $helper -O3 -D CUDA_5 -I ../extra_inc/cuda " * cuda_gencode_flags(["sm_90"])
        @test Nbody6Dynamics._cuflags_with_arch(mk, String[], ["-allow-unsupported-compiler"]) ==
              "-I $helper -O3 -D CUDA_5 -I ../extra_inc/cuda -allow-unsupported-compiler"
        @test Nbody6Dynamics._cuflags_with_arch(
            mk,
            ["sm_90"],
            ["-ccbin", "gcc-14"];
            helper_dir = "/h",
        ) ==
              "-I /h -O3 -D CUDA_5 -I ../extra_inc/cuda " *
              cuda_gencode_flags(["sm_90"]) *
              " -ccbin gcc-14"
        # The shipped helper headers replace the engine's 2012 copy, whose inline
        # functions read cudaDeviceProp fields that CUDA 13.0 removed.
        for name in ("helper_cuda.h", "helper_string.h", "LICENSE")
            @test isfile(joinpath(helper, name))
        end
        helper_src = read(joinpath(helper, "helper_cuda.h"), String)
        @test !occursin(r"deviceProp\.(clockRate|computeMode)", helper_src)
        @test occursin("cudaDevAttrClockRate", helper_src) &&
              occursin("checkCudaErrors", helper_src)
        # nvcc_flags: parsed, round-tripped, validated
        flags_cfg = load_config(
            write_cfg("[build]\nenable_gpu = true\nnvcc_flags = [\"-ccbin\", \"gcc-14\"]\n"),
        )
        @test flags_cfg.build.nvcc_flags == ["-ccbin", "gcc-14"]
        save_config(flags_cfg, rt)
        @test load_config(rt).build.nvcc_flags == ["-ccbin", "gcc-14"]
        @test BuildConfig().nvcc_flags == String[]
        @test_throws ErrorException load_config(write_cfg("[build]\nnvcc_flags = [\" \"]\n"))
        @test Nbody6Dynamics._cuflags_with_arch(joinpath(dir, "absent"), ["sm_90"]) ==
              "-I $helper -O3 " * cuda_gencode_flags(["sm_90"])
        @test Nbody6Dynamics._cuflags_with_arch(mk, String[]) ==
              "-I $helper -O3 -D CUDA_5 -I ../extra_inc/cuda"

        # Binary selection by suffix tags
        src = joinpath(dir, "backend")
        bdir = joinpath(src, "build")
        mkpath(bdir)
        @test_throws ErrorException Nbody6Dynamics._find_binary(src, "nbody6++", BuildConfig())
        for f in ("nbody6++.avx", "nbody6++.avx.gpu", "nbody6++.avx.mpi", "nbody6++.avx.mpi.gpu")
            touch(joinpath(bdir, f))
        end
        pick(; kw...) = basename(Nbody6Dynamics._find_binary(src, "nbody6++", BuildConfig(; kw...)))
        @test pick() == "nbody6++.avx"
        @test pick(; enable_gpu = true) == "nbody6++.avx.gpu"
        @test pick(; enable_mpi = true) == "nbody6++.avx.mpi"
        @test pick(; enable_mpi = true, enable_gpu = true) == "nbody6++.avx.mpi.gpu"
        rm(joinpath(bdir, "nbody6++.avx.gpu"))
        missing_gpu = try
            Nbody6Dynamics._find_binary(src, "nbody6++", BuildConfig(; enable_gpu = true))
        catch err
            err
        end
        @test missing_gpu isa ErrorException &&
              occursin("gpu = true", missing_gpu.msg) &&
              occursin("nbody6++.avx.mpi.gpu", missing_gpu.msg)
        @test_throws ErrorException Nbody6Dynamics._find_binary(
            joinpath(dir, "nowhere"),
            "nbody6++",
            BuildConfig(),
        )

        # Build record
        make_full(b, s) = Nbody6Config(
            InstallConfig(),
            b,
            s,
            PostprocessConfig(),
            VisualizationConfig(),
            MergerPipelineConfig(),
        )
        cfg_gpu = make_full(
            BuildConfig(; enable_gpu = true, cuda_arch = ["sm_90"]),
            SimulationConfig(; gpu_list = [0, 1]),
        )
        cfg_cpu = make_full(BuildConfig(), SimulationConfig())
        info_path = Nbody6Dynamics._write_build_info(
            src,
            cfg_gpu,
            ["--with-par=b1m", "--enable-mcmodel=large"],
            "/opt/cuda",
            ["sm_90"],
            joinpath(bdir, "nbody6++.avx.mpi.gpu"),
        )
        @test info_path == joinpath(bdir, "BUILD_INFO.toml")
        binfo = Nbody6Dynamics.TOML.parsefile(info_path)
        @test binfo["binary"] == "nbody6++.avx.mpi.gpu" && binfo["enable_gpu"] == true
        @test binfo["cuda_arch"] == ["sm_90"] && binfo["cuda_path"] == "/opt/cuda"
        @test binfo["configure_args"] == ["--with-par=b1m", "--enable-mcmodel=large"]
        @test haskey(binfo, "nvcc_release") && haskey(binfo, "backend_commit")
        # A CPU record in its own build tree (the record sits next to the binary)
        bdir_cpu = joinpath(dir, "backend_cpu", "build")
        mkpath(bdir_cpu)
        cpu_info = Nbody6Dynamics.TOML.parsefile(
            Nbody6Dynamics._write_build_info(
                dirname(bdir_cpu),
                cfg_cpu,
                String[],
                "",
                String[],
                joinpath(bdir_cpu, "nbody6++.avx"),
            ),
        )
        @test cpu_info["enable_gpu"] == false && !haskey(cpu_info, "cuda_arch")
        @test Nbody6Dynamics.TOML.parsefile(info_path)["cuda_arch"] == ["sm_90"]

        # Launch script: GPU_LIST exported only when configured
        args = (joinpath(dir, "nbody6++"), "in.inp", "out1000", "err1000")
        s_gpu = read(Nbody6Dynamics._write_launch_script(dir, args..., cfg_gpu), String)
        @test occursin("export GPU_LIST=\"0 1\"\n", s_gpu)
        s_cpu = read(Nbody6Dynamics._write_launch_script(dir, args..., cfg_cpu), String)
        @test !occursin("GPU_LIST", s_cpu)

        # Devices and kernel label the engine reports on stderr
        err = joinpath(dir, "err1000")
        write(
            err,
            "# GPU initialization - rank: 0; HOST node1; NGPU 2; device: 0 NVIDIA H200\n" *
            "# GPU initialization - rank: 0; HOST node1; NGPU 2; device: 1 NVIDIA H200\n" *
            "[R.0 GPU Reg.F ] Nsend 12  Ngrav 12 <Ni> 40   send(s) 0.01 grav(s) 0.02  nb(s) 0.0  out(s) 0.0  Perf.(Gflops) 1500.0\n" *
            "[R.0 GPU Reg.F ] Nsend 12  Ngrav 12 <Ni> 40   send(s) 0.01 grav(s) 0.02  nb(s) 0.0  out(s) 0.0  Perf.(Gflops) 2500.0\n" *
            "# GPU initialization - rank: 0; HOST node1; NGPU 2; device: 0 NVIDIA H200\n",
        )
        @test Nbody6Dynamics._reported_gpu_devices(err) ==
              ["rank 0: device 0 NVIDIA H200", "rank 0: device 1 NVIDIA H200"]
        @test Nbody6Dynamics._reported_gpu_devices(joinpath(dir, "absent")) == String[]
        gf = Nbody6Dynamics._force_kernel_gflops(err)
        @test gf["samples"] == 2 && gf["mean"] == 2000.0 && gf["peak"] == 2500.0
        @test gf["kernel"] == "GPU Reg.F"
        write(
            err,
            "[R.0 AVX Reg.F ] x Perf.(Gflops) 100.0\n[R.0 GPU Reg.F ] x Perf.(Gflops) 300.0\n",
        )
        @test Nbody6Dynamics._force_kernel_gflops(err)["kernel"] == "AVX Reg.F, GPU Reg.F"

        # Run summary: gpu_list, reported devices, and the build table
        run_dir = joinpath(dir, "run")
        out_dir = joinpath(run_dir, "output")
        mkpath(out_dir)
        write(joinpath(out_dir, "out1000"), "header\n")
        write(
            joinpath(out_dir, "err1000"),
            "# GPU initialization - rank: 0; HOST n; NGPU 1; device: 0 NVIDIA GeForce RTX 5090\n",
        )
        cp(info_path, joinpath(out_dir, "BUILD_INFO.toml"))
        Nbody6Dynamics._write_run_summary(
            cfg_gpu,
            run_dir,
            "gpu_a1b2",
            joinpath(out_dir, "out1000"),
            out_dir,
            1.0,
        )
        rinfo = Nbody6Dynamics.TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))
        @test rinfo["run"]["gpu_list"] == [0, 1]
        @test rinfo["run"]["gpu_devices"] == ["rank 0: device 0 NVIDIA GeForce RTX 5090"]
        @test rinfo["build"]["cuda_arch"] == ["sm_90"]
        @test haskey(rinfo["hardware"], "gpu")
        # A CPU run carries neither the GPU keys nor a build table without the record
        cpu_run = joinpath(dir, "cpu_run")
        mkpath(joinpath(cpu_run, "output"))
        Nbody6Dynamics._write_run_summary(
            cfg_cpu,
            cpu_run,
            "cpu_a1b2",
            joinpath(cpu_run, "output", "out1000"),
            joinpath(cpu_run, "output"),
            1.0,
        )
        cinfo = Nbody6Dynamics.TOML.parsefile(joinpath(cpu_run, "RUN_INFO.toml"))
        @test !haskey(cinfo["run"], "gpu_list") && !haskey(cinfo["run"], "gpu_devices")
        @test !haskey(cinfo, "build")
    end

    # =====================================================================
    @testset "GPU validation driver" begin
        # The dependency check accepts an nvcc that sits under the toolkit rather
        # than on PATH (the build exports <cuda_path>/bin itself).
        cuda_dir = mktempdir()
        mkpath(joinpath(cuda_dir, "bin"))
        touch(joinpath(cuda_dir, "bin", "nvcc"))
        dep_cfg(body) = (p = joinpath(cuda_dir, "cfg.toml"); write(p, body); load_config(p))
        @test !(
            "nvcc" in
            check_dependencies(dep_cfg("[build]\nenable_gpu = true\ncuda_path = \"$cuda_dir\"\n"))
        )
        @test !("nvcc" in check_dependencies(dep_cfg("[build]\nenable_gpu = false\n")))
        if !Nbody6Dynamics.check_command("nvcc")
            absent = joinpath(cuda_dir, "absent")
            @test "nvcc" in check_dependencies(
                dep_cfg("[build]\nenable_gpu = true\ncuda_path = \"$absent\"\n"),
            )
        end

        # nvcc host-compiler probe: classification of the compiler message
        @test Nbody6Dynamics._unsupported_host_compiler(
            "#error -- unsupported GNU version! gcc versions later than 14 are not supported!",
        )
        @test Nbody6Dynamics._unsupported_host_compiler("unsupported clang version")
        @test !Nbody6Dynamics._unsupported_host_compiler("probe.cu(1): error: expected a \";\"")
        # No nvcc under the path: the probe stands aside and the build reports it
        @test Nbody6Dynamics._nvcc_host_compiler_flags("/nonexistent/cuda") == String[]

        # Logged stage runner: merged output, exit code, escape sequences stripped
        vdir = mktempdir()
        log = joinpath(vdir, "stage.log")
        julia = Base.julia_cmd()
        @test Nbody6Dynamics._tee_run(
            `$julia -e 'print("\e[31mred\e[0m\n"); println(stderr, "err line")'`,
            log,
        ) == 0
        lines = readlines(log)
        @test "red" in lines && "err line" in lines
        @test Nbody6Dynamics._tee_run(`$julia -e 'exit(3)'`, log) == 3
        @test Nbody6Dynamics._tool_banner(
            `$julia -e 'println("banner line"); println("second")'`,
        ) == "banner line"
        @test Nbody6Dynamics._tool_banner(`/nonexistent/tool --version`) == "unavailable"

        # Probe against a stand-in nvcc: a script that mimics the toolkit's
        # verdicts (rejects the default compiler, accepts an override or a
        # -ccbin), so the search logic runs without a CUDA toolkit.
        function fake_nvcc(dir, body)
            mkpath(joinpath(dir, "bin"))
            path = joinpath(dir, "bin", "nvcc")
            write(path, "#!/bin/sh\n" * body)
            chmod(path, 0o755)
            return dir
        end
        accepts_override = fake_nvcc(
            mktempdir(),
            "case \" \$* \" in *' -allow-unsupported-compiler '*) exit 0;; esac\n" *
            "echo '#error -- unsupported GNU version! gcc versions later than 15 are not supported!' >&2\nexit 1\n",
        )
        flags =
            @test_logs (:warn, r"building with -allow-unsupported-compiler") Nbody6Dynamics._nvcc_host_compiler_flags(
                accepts_override;
                candidates = String[],
            )
        @test flags == ["-allow-unsupported-compiler"]
        needs_ccbin = fake_nvcc(
            mktempdir(),
            "case \" \$* \" in *' -ccbin /usr/bin/true '*) exit 0;; esac\n" *
            "echo 'unsupported GNU version! gcc versions later than 15 are not supported!' >&2\n" *
            "for i in \$(seq 1 60); do echo \"type_traits(\$i): error: identifier char8_t is undefined\" >&2; done\nexit 1\n",
        )
        flags =
            @test_logs (:warn, r"building with -ccbin /usr/bin/true") Nbody6Dynamics._nvcc_host_compiler_flags(
                needs_ccbin;
                candidates = ["/nonexistent/g++-99", "/usr/bin/true"],
            )
        @test flags == ["-ccbin", "/usr/bin/true"]
        hopeless = try
            Nbody6Dynamics._nvcc_host_compiler_flags(needs_ccbin; candidates = String[])
        catch err
            err
        end
        @test hopeless isa ErrorException
        @test occursin("default host compiler; -allow-unsupported-compiler", hopeless.msg)
        @test occursin("gcc15-c++", hopeless.msg) && occursin("lines omitted", hopeless.msg)
        @test occursin("--- default host compiler ---", hopeless.msg) &&
              occursin("--- -allow-unsupported-compiler ---", hopeless.msg)
        # Every attempt's own output is reported, so a rejected -ccbin candidate
        # can be diagnosed from the message alone.
        echoing = fake_nvcc(
            mktempdir(),
            "echo \"nvcc args: \$*\" >&2\n" *
            "echo 'unsupported GNU version! gcc versions later than 15 are not supported!' >&2\nexit 1\n",
        )
        per_attempt = try
            Nbody6Dynamics._nvcc_host_compiler_flags(echoing; candidates = ["/nonexistent/g++-99"])
        catch err
            err
        end
        @test per_attempt isa ErrorException
        @test occursin(
            r"--- -ccbin /nonexistent/g\+\+-99 ---\nnvcc args: [^\n]* -ccbin /nonexistent/g\+\+-99\n",
            per_attempt.msg,
        )
        @test occursin(
            r"--- -allow-unsupported-compiler ---\nnvcc args: [^\n]* -allow-unsupported-compiler\n",
            per_attempt.msg,
        )
        # glibc ≥ 2.42 declares rsqrt/rsqrtf with an exception specification
        # the CUDA headers lack: the probe retries with the feature-macro
        # override, alone or on top of the host-compiler choice.
        glibc_line =
            "/usr/include/bits/mathcalls.h(206): error: exception specification is incompatible " *
            "with that of previous function \"rsqrt\" (declared at line 629 of crt/math_functions.h)"
        @test Nbody6Dynamics._glibc_c2y_conflict(glibc_line)
        @test !Nbody6Dynamics._glibc_c2y_conflict(
            "unsupported GNU version! gcc versions later than 15",
        )
        glibc_msg = "echo '$glibc_line' >&2\nexit 1\n"
        glibc_only = fake_nvcc(
            mktempdir(),
            "case \" \$* \" in *' -U_GNU_SOURCE -D_DEFAULT_SOURCE '*) exit 0;; esac\n" * glibc_msg,
        )
        flags =
            @test_logs (:warn, r"glibc declares rsqrt and rsqrtf") Nbody6Dynamics._nvcc_host_compiler_flags(
                glibc_only;
                candidates = String[],
            )
        @test flags == ["-U_GNU_SOURCE", "-D_DEFAULT_SOURCE"]
        # A configured -ccbin stays in nvcc_flags; only the override is added
        flags =
            @test_logs (:warn, r"glibc declares rsqrt and rsqrtf") Nbody6Dynamics._nvcc_host_compiler_flags(
                glibc_only;
                nvcc_flags = ["-ccbin", "/usr/bin/true"],
                candidates = String[],
            )
        @test flags == ["-U_GNU_SOURCE", "-D_DEFAULT_SOURCE"]
        # Host compiler rejected AND the glibc conflict (Fedora 44 with CUDA
        # 13.1): the first accepted -ccbin candidate, with the override
        both = fake_nvcc(
            mktempdir(),
            "case \" \$* \" in\n" *
            "  *' -ccbin /usr/bin/true -U_GNU_SOURCE -D_DEFAULT_SOURCE '*) exit 0;;\n" *
            "  *' -ccbin /usr/bin/true '*) " *
            glibc_msg *
            ";;\n" *
            "esac\n" *
            "echo 'unsupported GNU version! gcc versions later than 15 are not supported!' >&2\nexit 1\n",
        )
        flags = @test_logs (:warn, r"building with -ccbin /usr/bin/true") (
            :warn,
            r"glibc declares rsqrt and rsqrtf",
        ) Nbody6Dynamics._nvcc_host_compiler_flags(
            both;
            candidates = ["/nonexistent/g++-99", "/usr/bin/true"],
        )
        @test flags == ["-ccbin", "/usr/bin/true", "-U_GNU_SOURCE", "-D_DEFAULT_SOURCE"]
        # The override does not help: every attempt is reported, with the
        # header-patch advice
        glibc_stuck = fake_nvcc(mktempdir(), glibc_msg)
        stuck = try
            Nbody6Dynamics._nvcc_host_compiler_flags(glibc_stuck; candidates = String[])
        catch err
            err
        end
        @test stuck isa ErrorException
        @test occursin("default host compiler with -U_GNU_SOURCE -D_DEFAULT_SOURCE", stuck.msg)
        @test occursin("math_functions.h", stuck.msg) && occursin("noexcept(true)", stuck.msg)
        stuck_pinned = try
            Nbody6Dynamics._nvcc_host_compiler_flags(
                glibc_stuck;
                nvcc_flags = ["-ccbin", "/usr/bin/false"],
                candidates = String[],
            )
        catch err
            err
        end
        @test stuck_pinned isa ErrorException &&
              occursin("configured host compiler", stuck_pinned.msg)
        @test occursin("math_functions.h", stuck_pinned.msg)
        # Host-record verdict of the probe
        @test Nbody6Dynamics._nvcc_probe_record("/nonexistent/cuda", "") ==
              (String[], "skipped: nvcc unavailable")
        flags, verdict = Nbody6Dynamics._nvcc_probe_record(echoing, "13.1")
        @test flags == String[] && occursin("--- default host compiler ---", verdict)
        flags, verdict =
            @test_logs (:warn, r"building with -allow-unsupported-compiler") Nbody6Dynamics._nvcc_probe_record(
                accepts_override,
                "13.1",
            )
        @test flags == ["-allow-unsupported-compiler"] && verdict == "passed"
        # A configured -ccbin is final: no search, the failure is reported as is
        pinned_fail = try
            Nbody6Dynamics._nvcc_host_compiler_flags(
                needs_ccbin;
                nvcc_flags = ["-ccbin", "/usr/bin/false"],
                candidates = ["/usr/bin/true"],
            )
        catch err
            err
        end
        @test pinned_fail isa ErrorException &&
              occursin("configured host compiler", pinned_fail.msg)
        @test Nbody6Dynamics._host_compiler_candidates() isa Vector{String}
        @test Nbody6Dynamics._output_excerpt("a\nb\nc", 1, 1) == "a\n… (1 lines omitted)\nc"
        @test Nbody6Dynamics._output_excerpt("a\nb", 5, 5) == "a\nb"
        @test Nbody6Dynamics._run_capture(
            `$(Base.julia_cmd()) -e 'println(stderr, "e"); print("o")'`,
            joinpath(vdir, "cap.log"),
        ) == (true, "e\no")

        # Dry run: host record and planned commands, nothing executed
        @test_throws ArgumentError run_gpu_validation(;
            base_dir = vdir,
            stages = [:nope],
            dry_run = true,
        )
        @test_throws ArgumentError run_gpu_validation(;
            base_dir = vdir,
            stages = Symbol[],
            dry_run = true,
        )
        @test_throws ArgumentError run_gpu_validation(;
            base_dir = vdir,
            bench_tcrit = 0.0,
            dry_run = true,
        )
        out = run_gpu_validation(;
            base_dir = vdir,
            stages = [:suite, :bench],
            bench_n = [1000],
            bench_threads = [2],
            bench_gpu_lists = [[0], [0, 1]],
            bench_tcrit = 0.5,
            dry_run = true,
        )
        @test startswith(basename(out), "gpu_validation_") && dirname(out) == joinpath(vdir, "runs")
        host = Nbody6Dynamics.TOML.parsefile(joinpath(out, "HOST_INFO.toml"))
        for key in (
            "host",
            "gpu",
            "compute_capabilities",
            "nvcc_release",
            "gcc",
            "gfortran",
            "package_commit",
            "glibc",
            "host_compilers",
            "nvcc_host_flags",
            "nvcc_probe",
        )
            @test haskey(host, key)
        end
        @test host["nvcc_probe"] isa String && host["nvcc_host_flags"] isa Vector
        summary = Nbody6Dynamics.TOML.parsefile(joinpath(out, "VALIDATION.toml"))
        @test summary["dry_run"] && summary["stages"] == ["suite", "bench"]
        @test summary["results"]["suite"]["status"] == "planned"
        @test occursin("NBODY6_GPU_TESTS=1", summary["results"]["suite"]["command"])
        bench_cmd = summary["results"]["bench"]["command"]
        @test occursin("gpu_scaling.jl 1000 2 0;0,1 0.5", bench_cmd)
        @test occursin("NBODY6_GPU_BACKEND=", bench_cmd) &&
              occursin("Nbody6PPGPU-beijing-gpu", bench_cmd)
        @test !haskey(summary["results"], "gpu") && isempty(filter(endswith(".log"), readdir(out)))
        # A stage whose prerequisites are missing is skipped with the reason, no process spawned
        bare = mktempdir()
        skipped = run_gpu_validation(; base_dir = bare, stages = [:bench])
        skipped_summary = Nbody6Dynamics.TOML.parsefile(joinpath(skipped, "VALIDATION.toml"))
        @test skipped_summary["results"]["bench"]["status"] == "skipped"
        @test occursin("Nbody6PPGPU-beijing-gpu", skipped_summary["results"]["bench"]["reason"])
        @test !isfile(joinpath(skipped, "bench.log")) && haskey(skipped_summary, "finished")
        @test Nbody6Dynamics._stage_prerequisite(:suite, bare) === nothing
        # Benchmark artefact collection on an empty bench tree is a no-op
        @test Nbody6Dynamics._collect_bench_artefacts(joinpath(vdir, "bench"), out, 0.0) == String[]
    end

    # =====================================================================
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
        cp(
            joinpath(@__DIR__, "..", "input_files", "merger_demo_small.toml"),
            joinpath(sw, "m.toml"),
        )
        st = joinpath(sw, "sweep.toml")
        write(
            st,
            "[sweep]\nname = \"t\"\npipeline_config = \"pipeline.toml\"\nmerger_config = \"m.toml\"\nseeds = [1]\n",
        )
        @test Nbody6Dynamics._check_sweep_binary(load_sweep_config(st)) ==
              joinpath(bk, "nbody6++.avx")
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
    @testset "Fortran binary I/O" begin
        # Create a synthetic Fortran binary file
        fpath = joinpath(TESTDIR, "test_fortran.bin")
        open(fpath, "w") do io
            # Write a record: 3 Float32 values
            data = Float32[1.0, 2.0, 3.0]
            marker = Int32(sizeof(data))
            write(io, marker)
            write(io, data)
            write(io, marker)

            # Write a second record: 2 Int32 values
            idata = Int32[42, 99]
            marker2 = Int32(sizeof(idata))
            write(io, marker2)
            write(io, idata)
            write(io, marker2)
        end

        open(fpath, "r") do io
            rec1 = Nbody6Dynamics.read_fortran_record(io, Float32, 3)
            @test rec1 ≈ Float32[1.0, 2.0, 3.0]

            rec2 = Nbody6Dynamics.read_fortran_record(io, Int32, 2)
            @test rec2 == Int32[42, 99]
        end

        # peek_record_size
        open(fpath, "r") do io
            sz = Nbody6Dynamics.peek_record_size(io)
            @test sz == Int32(12)  # 3 × 4 bytes
        end
    end

    # =====================================================================
    @testset "conf.3 reader" begin
        # Create a synthetic conf.3 file in standard format
        fpath = joinpath(TESTDIR, "conf.3_test")
        n_particles = 5

        open(fpath, "w") do io
            # Record 1: header integers
            hdr = Int32[n_particles, 1, 1, 20]
            _write_fortran_record(io, hdr)

            # Record 2: AS(1:20) parameters
            params = zeros(Float32, 20)
            params[1] = 0.5f0    # time
            params[3] = 1.0f0    # rbar
            params[4] = 0.5f0    # zmbar
            params[11] = 10.0f0  # tscale
            params[12] = 5.0f0   # vstar
            params[18] = 2.0f0   # rscale
            _write_fortran_record(io, params)

            # Particle records (standard: M, X, Y, Z, VX, VY, VZ, NAME)
            for i in 1:n_particles
                m = Float32(1.0 / n_particles)
                x, y, z = Float32.(randn(3))
                vx, vy, vz = Float32.(0.1 .* randn(3))
                name = Int32(i)
                buf = IOBuffer()
                write(buf, m, x, y, z, vx, vy, vz, name)
                _write_fortran_record(io, take!(buf))
            end
        end

        snap = read_conf3(fpath)
        @test nparticles(snap) == n_particles
        @test time_nb(snap.header) ≈ 0.5
        @test tscale(snap.header) ≈ 10.0
        @test rscale(snap.header) ≈ 2.0
        @test length(snap.mass) == n_particles
        @test size(snap.pos) == (3, n_particles)
        @test size(snap.vel) == (3, n_particles)
        @test snap.name == Int32.(1:n_particles)
        @test isempty(snap.rho)  # standard format has no density

        # Test read_all_conf3
        snaps = read_all_conf3(TESTDIR, "conf.3_")
        @test length(snaps) >= 1
    end

    # =====================================================================
    @testset "conf.3 reader — extended format" begin
        fpath = joinpath(TESTDIR, "conf.3_ext")
        n_particles = 3

        open(fpath, "w") do io
            hdr = Int32[n_particles, 1, 1, 20]
            _write_fortran_record(io, hdr)

            params = zeros(Float32, 20)
            params[1] = 1.0f0
            _write_fortran_record(io, params)

            # Extended: M, RHO, XNS, X, Y, Z, VX, VY, VZ, PHI, NAME
            for i in 1:n_particles
                buf = IOBuffer()
                write(buf, Float32(0.1 * i))   # mass
                write(buf, Float32(100.0))      # rho
                write(buf, Float32(0.01))       # xns
                write(buf, Float32.(randn(3))...)  # pos
                write(buf, Float32.(0.1 .* randn(3))...)  # vel
                write(buf, Float32(-0.5))       # phi
                write(buf, Int32(i))            # name
                _write_fortran_record(io, take!(buf))
            end
        end

        snap = read_conf3(fpath)
        @test nparticles(snap) == n_particles
        @test length(snap.rho) == n_particles
        @test length(snap.phi) == n_particles
        @test all(snap.rho .≈ 100.0f0)
        @test all(snap.phi .≈ -0.5f0)
    end

    # =====================================================================
    @testset "Diagnostics parser" begin
        diag_path = joinpath(TESTDIR, "out1000")
        write(
            diag_path,
            """
 Some header text...

 PHYSICAL SCALING:  R* = 1.234  M* = 5678.0  V* = 3.456  T* = 7.890
                    <M> = 0.567  SU = 1.0  AU = 1.0

 ADJUST:    0.0000      0.00  1.000  0.00E+00 -0.2500  10000    500   1.234
 ADJUST:    0.5000     50.00  0.987  1.23E-06 -0.2499   9998    498   1.245
 ADJUST:    1.0000    100.00  0.995  2.45E-06 -0.2498   9990    495   1.267

 END RUN
""",
        )

        diag = read_diagnostics(diag_path)

        @test length(diag.adjust) == 3
        @test diag.adjust[1].time_nb ≈ 0.0
        @test diag.adjust[2].time_myr ≈ 50.0
        @test diag.adjust[3].n == 9990
        @test diag.adjust[2].npairs == 498

        # Physical scaling
        @test haskey(diag.physical_scaling, "R*")
        @test diag.physical_scaling["R*"] ≈ 1.234
        @test diag.physical_scaling["V*"] ≈ 3.456

        # Unit extraction
        units = extract_scaling(diag)
        @test units.rbar ≈ 1.234
        @test units.tscale ≈ 7.890
    end

    # =====================================================================
    @testset "Lagrangian radii reader" begin
        lagr_path = joinpath(TESTDIR, "lagr.7")

        # Synthetic lagr.7 in the upstream lagr.f layout: ##/TIME header,
        # then one row per epoch (3 epochs, 5 radii columns, simplified)
        write(
            lagr_path,
            """
## Time, [Column number | Label]:    2 R_{lagr}   (synthetic)
TIME   0.01  0.05  0.20  0.50  1.00
0.0  0.01  0.05  0.20  0.50  1.00
0.5  0.012 0.055 0.22  0.52  1.05
1.0  0.015 0.060 0.25  0.55  1.10
""",
        )

        lagr = read_lagr(lagr_path)

        @test length(lagr.time) == 3
        @test lagr.time ≈ [0.0, 0.5, 1.0]
        @test size(lagr.radii) == (5, 3)
        @test lagr.radii[1, 1] ≈ 0.01
        @test lagr.radii[5, 3] ≈ 1.10
    end

    # =====================================================================
    @testset "Escaper reader" begin
        # Real esc.11 shape (escape.F): 5 NB-unit diagnostics, then the
        # physical-unit block T[Myr] M[M*] EESC VI[km/s] K* NAME, then the
        # direction angles ANGLE PHI / ANGLE THETA (tokens 12–13).  The last
        # line is truncated to 11 tokens (no angles) and must be skipped.
        esc_path = joinpath(TESTDIR, "esc.11")
        write(
            esc_path,
            """
         TTOT         BODY         RI           VI           STEP         T[Myr]       M[M*]      EESC      VI[km/s]     K*  NAME      ANGLE PHI     ANGLE THETA
   1.00000E+00  5.00000E-04  2.00000E+01  4.00000E+01  1.95312E-03  1.23400E+00  5.00000E-01 -1.23000E-01  1.56000E+01   0       101  1.20000E+01  2.50000E+01
   2.00000E+00  3.00000E-04  2.50000E+01  4.20000E+01  1.95312E-03  2.56700E+00  3.00000E-01  4.56000E-01  2.23000E+01   1       202  1.85000E+02 -3.00000E+01
   3.00000E+00  1.20000E-03  3.00000E+01  4.40000E+01  1.95312E-03  3.89000E+00  1.20000E+00 -7.89000E-01  1.01000E+01  14       303  3.40000E+02  6.00000E+01
   4.00000E+00  1.00000E-03  3.10000E+01  4.50000E+01  1.95312E-03  4.50000E+00  1.00000E+00  1.00000E-01  1.20000E+01   1       404
""",
        )

        escs = read_escapers(esc_path)
        @test length(escs) == 3   # the 11-token line is skipped

        @test escs[1].time_myr ≈ 1.234
        @test escs[1].mass_solar ≈ 0.500
        @test escs[1].escape_energy ≈ -0.123
        @test escs[1].velocity_kms ≈ 15.6
        @test escs[1].stellar_type == 0
        @test escs[1].name == 101
        @test escs[1].phi_deg ≈ 12.0
        @test escs[1].theta_deg ≈ 25.0

        @test escs[2].phi_deg ≈ 185.0
        @test escs[2].theta_deg ≈ -30.0

        @test escs[3].stellar_type == 14   # BH (Hurley convention)
        @test escs[3].name == 303
        @test escs[3].phi_deg ≈ 340.0
        @test escs[3].theta_deg ≈ 60.0

        # Empty file
        empty_path = joinpath(TESTDIR, "esc_empty.11")
        write(empty_path, "")
        @test isempty(read_escapers(empty_path))
    end

    # =====================================================================
    @testset "Stellar evolution reader" begin
        # Header time is TPHYS [Myr]; data-line token 1 is TTOT [NB] —
        # different clocks, both kept. 15-token v2026.07+ layout.
        sev_path = joinpath(TESTDIR, "sev.83_0")
        write(
            sev_path,
            """
  3  0.5000
   0.4315   1   101  0  1.20  0.800   0.123  -0.456  3.750  0.0  0.0  90.0  0.0  0.0  1.0e-10
   0.4315   2   202  1  0.80  1.200   1.500   0.200  4.100  0.0  0.0  5.20  0.0  0.0  1.0e-10
   0.4315   3   303  13 2.50  10.00   5.000   1.500  4.500  0.0  0.0  4.10  1.4  1.0e-5  1.0e-10
""",
        )

        sev = read_stellar_evolution(sev_path)
        @test sev.n_stars == 3
        @test sev.time_myr ≈ 0.5
        @test length(sev.records) == 3

        r1 = sev.records[1]
        @test r1.time_nb ≈ 0.4315   # per-line TTOT [NB], not the Myr header
        @test r1.index == Int32(1)
        @test r1.name == Int32(101)
        @test r1.stellar_type == Int32(0)  # MS
        @test r1.ri ≈ 1.20
        @test r1.mass_solar ≈ 0.800
        @test r1.log_luminosity ≈ 0.123
        @test r1.log_radius ≈ -0.456
        @test r1.log_teff ≈ 3.750
        @test r1.ms_lifetime_myr ≈ 90.0    # TM

        r3 = sev.records[3]
        @test r3.stellar_type == Int32(13)  # NS
        @test r3.mass_core ≈ 1.4

        # Test read_all_stellar_evolution
        sev2_path = joinpath(TESTDIR, "sev.83_1")
        write(
            sev2_path,
            """
  2  1.0000
   1.0000   1   101  2  1.50  0.750   0.500  -0.200  3.600  0.0  0.0  8.0  0.2  0.01  2.0
   1.0000   2   202  4  0.90  1.100   2.000   0.400  3.900  0.0  0.0  6.5  0.4  0.02  4.0
""",
        )

        sevs = read_all_stellar_evolution(TESTDIR, "sev.83_*")
        @test length(sevs) == 2
        @test sevs[1].time_myr < sevs[2].time_myr
    end

    # =====================================================================
    @testset "Binary evolution reader" begin
        # hrplot.F FORMAT 5: 32 tokens per line; header NPAIRS TPHYS [Myr];
        # token 1 of every data line is TTOT [NB]. Record 1 is a verbatim
        # engine line, record 2 carries a recognisable 1..12 tail.
        bev_path = joinpath(TESTDIR, "bev.82_0")
        write(
            bev_path,
            """
       2      0.0
   0.00000E+00       1       2       1       2  0  1   0  1.34446E+00  9.21880E-01  1.20647E+00  1.30277E+00  2.21134E-01  1.98278E-01 -1.99897E+00 -2.08678E+00 -6.28737E-01 -6.64486E-01  3.57640E+00  3.57233E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  4.98658E+05  6.20908E+05  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  1.52819E-01  1.40743E-01
   0.00000E+00       3       4       3       4  1  0   0  5.14787E-01  6.97120E-01  5.58338E+00  4.42220E+00  1.11885E+00  5.68179E-01  3.81907E-01 -9.85920E-01 -1.23430E-02 -2.97930E-01  3.86343E+00  3.70000E+00  1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0  9.0  10.0  11.0  12.0
   0.00000E+00       5       6
""",
        )
        bev = read_binary_evolution(bev_path)
        @test bev.n_pairs == 2
        @test bev.time_myr == 0.0
        @test length(bev.records) == 2        # the short line is skipped

        r1 = bev.records[1]
        @test r1.time_nb == 0.0
        @test (r1.index1, r1.index2, r1.name1, r1.name2) == (1, 2, 1, 2)
        @test (r1.stellar_type1, r1.stellar_type2, r1.stellar_type_cm) == (0, 1, 0)
        @test r1.ri ≈ 1.34446
        @test r1.eccentricity ≈ 0.92188
        @test r1.log_period_days ≈ 1.20647
        @test r1.log_semi_major_axis_rsun ≈ 1.30277
        @test r1.mass1 ≈ 0.221134
        @test r1.mass2 ≈ 0.198278
        @test r1.log_luminosity1 ≈ -1.99897
        @test r1.log_radius2 ≈ -0.664486
        @test r1.log_teff1 ≈ 3.5764
        @test r1.ms_lifetime1_myr ≈ 4.98658e5     # TM of a 0.22 M☉ dwarf
        @test r1.radius_envelope2 ≈ 0.140743

        r2 = bev.records[2]
        @test (r2.age1_myr, r2.age2_myr) == (1.0, 2.0)
        @test (r2.epoch1_myr, r2.epoch2_myr) == (3.0, 4.0)
        @test (r2.ms_lifetime1_myr, r2.ms_lifetime2_myr) == (5.0, 6.0)
        @test (r2.mass_core1, r2.mass_core2) == (7.0, 8.0)
        @test (r2.radius_core1, r2.radius_core2) == (9.0, 10.0)
        @test (r2.radius_envelope1, r2.radius_envelope2) == (11.0, 12.0)

        # The engine's period and semi-major axis columns obey Kepler III,
        # a³/P² = m₁ + m₂ in au, yr, M☉ — a check of the column map and of
        # the R☉ → pc → au conversions.
        P_yr = 10^r1.log_period_days / 365.25
        a_au = semi_major_axis_pc(r1) / Nbody6Dynamics._AU_IN_PC
        @test a_au^3 / P_yr^2 ≈ r1.mass1 + r1.mass2 rtol = 1e-2
        @test binding_energy(r1) ≈
              Nbody6Dynamics._G_PC_KMS2_MSUN * r1.mass1 * r1.mass2 / (2 * semi_major_axis_pc(r1))

        write(joinpath(TESTDIR, "bev.82_1"), "       0      1.5\n")
        bevs = read_all_binary_evolution(TESTDIR, "bev.82_*")
        @test length(bevs) == 2
        @test bevs[1].time_myr < bevs[2].time_myr
        @test bevs[2].n_pairs == 0 && isempty(bevs[2].records)
        @test isempty(read_all_binary_evolution(joinpath(TESTDIR, "absent"), "bev.82_*"))
        @test_throws ErrorException read_binary_evolution(joinpath(TESTDIR, "missing.82"))
    end

    # =====================================================================
    @testset "Binary population diagnostics" begin
        # Synthetic snapshot: four singles and one pair with known velocities.
        # Scalings: 10 M☉ per NB mass unit, 2 km/s per NB velocity unit.
        params = zeros(Float32, 20)
        params[3] = 1.0f0
        params[4] = 10.0f0
        params[11] = 1.0f0
        params[12] = 2.0f0
        names = Int32.(1:6)
        mass = Float32[0.1, 0.1, 0.1, 0.1, 0.3, 0.3]
        pos = zeros(Float32, 3, 6)
        vel = zeros(Float32, 3, 6)
        vel[1, 1:4] .= Float32[1, -1, 1, -1]     # singles: ±1 along x
        vel[1, 5], vel[1, 6] = 3.0f0, -3.0f0     # pair (5, 6): equal masses, c.m. at rest in x
        vel[2, 5] = vel[2, 6] = 0.5f0            # c.m. moving at 0.5 along y
        hdr = SnapshotHeader(Int32(6), Int32(1), Int32(1), Int32(20), params)
        snap = Snapshot(hdr, names, mass, pos, vel, Float32[], Float32[])

        _bin_rec(n1, n2, m1, m2, loga, e) = BinaryRecord(
            0.0,
            Int32(n1),
            Int32(n2),
            Int32(n1),
            Int32(n2),
            Int32(0),
            Int32(0),
            Int32(0),
            0.0,
            e,
            0.0,
            loga,
            m1,
            m2,
            zeros(18)...,
        )
        hard_rec = _bin_rec(5, 6, 3.0, 3.0, 3.0, 0.1)   # a = 10³ R☉
        soft_rec = _bin_rec(5, 6, 3.0, 3.0, 6.0, 0.9)   # a = 10⁶ R☉
        bev0 = BinaryEvolutionSnapshot(0.0, 2, [hard_rec, soft_rec])
        bev1 = BinaryEvolutionSnapshot(1.0, 1, [hard_rec])

        # Systems: four singles (m = 0.1, v = ±x̂) and the pair's c.m.
        # (m = 0.6, v = 0.5 ŷ). M = 1, v_c = 0.3 ŷ,
        # σ² = [4 × 0.1 × (1 + 0.09) + 0.6 × 0.04] / (3 × 1) = 0.46/3.
        sc = hardness_scale(snap, bev0)
        @test sc.n_systems == 5
        @test sc.m_mean ≈ 2.0 rtol = 1e-6            # Float32 particle masses
        @test sc.sigma_kms ≈ 2 * sqrt(0.46 / 3) rtol = 1e-6
        # A pair with a component missing from the snapshot is not reduced.
        sc_missing = hardness_scale(
            snap,
            BinaryEvolutionSnapshot(0.0, 1, [_bin_rec(5, 99, 3.0, 3.0, 3.0, 0.0)]),
        )
        @test sc_missing.n_systems == 6
        @test_throws ArgumentError hardness_scale(
            Snapshot(
                hdr,
                Int32[],
                Float32[],
                zeros(Float32, 3, 0),
                zeros(Float32, 3, 0),
                Float32[],
                Float32[],
            ),
            bev0,
        )

        x = binary_hardness(bev0, sc.m_mean, sc.sigma_kms)
        scale = sc.m_mean * sc.sigma_kms^2
        @test x[1] ≈ binding_energy(hard_rec) / scale
        @test x[1] > 1 && x[2] < 1
        @test x[1] / x[2] ≈ 1e3            # same masses, a ratio 10³
        @test semi_major_axis_pc(hard_rec) ≈ 1e3 * Nbody6Dynamics._RSUN_PC
        @test_throws ArgumentError binary_hardness(bev0, 0.0, 1.0)
        @test_throws ArgumentError binary_hardness(bev0, 1.0, -1.0)

        pop = binary_population(
            [bev0, bev1];
            n_stars = 6,
            m_mean = sc.m_mean,
            sigma_kms = sc.sigma_kms,
        )
        @test pop.classified
        @test pop.time_myr == [0.0, 1.0]
        @test pop.n_pairs == [2, 1]
        @test pop.n_hard == [1, 1]
        @test pop.n_soft == [1, 0]
        @test pop.binary_fraction ≈ [2 / 4, 1 / 5]   # N_b / (N_stars − N_b)
        pop_u = binary_population([bev0, bev1])
        @test !pop_u.classified
        @test all(iszero, pop_u.n_hard) && pop_u.n_soft == pop_u.n_pairs
        @test all(isnan, pop_u.binary_fraction)
        pop_v = binary_population(
            [bev0, bev1];
            n_stars = [6, 4],
            m_mean = [2.0, 2.0],
            sigma_kms = [sc.sigma_kms, 1e-3],
        )
        @test pop_v.binary_fraction ≈ [0.5, 1 / 3]
        @test pop_v.n_hard == [1, 1]
        @test_throws DimensionMismatch binary_population([bev0, bev1]; n_stars = [6, 6, 6])
        @test_throws ArgumentError binary_population(BinaryEvolutionSnapshot[])

        scales = binary_scales([bev0, bev1], [snap])
        @test scales.n_stars == [6, 6]
        @test scales.m_mean ≈ [2.0, 2.0] rtol = 1e-6
        @test scales.sigma_kms ≈ fill(sc.sigma_kms, 2) rtol = 1e-6
        @test binary_scales([bev0], Snapshot[]) === nothing

        # Figures: time series, single epoch, classified and unclassified
        # elements, period histograms, and the no-data paths.
        vis_b = VisualizationConfig(; format = "png", column = "single", output_dir = mktempdir())
        @test isfile(plot_binary_population(pop, vis_b; filename = "pop_series"))
        @test isfile(
            plot_binary_population(
                binary_population([bev0]; n_stars = 6),
                vis_b;
                filename = "pop_single",
            ),
        )
        @test isfile(
            plot_binary_orbital_elements(
                bev0,
                vis_b;
                m_mean = sc.m_mean,
                sigma_kms = sc.sigma_kms,
                filename = "ae_classified",
            ),
        )
        @test isfile(plot_binary_orbital_elements(bev0, vis_b; filename = "ae_plain"))
        @test isfile(plot_binary_period_distribution([bev0, bev1], vis_b; filename = "period_two"))
        @test isfile(plot_binary_period_distribution([bev0], vis_b; filename = "period_one"))
        empty_bev = BinaryEvolutionSnapshot(2.0, 0, BinaryRecord[])
        @test isfile(plot_binary_orbital_elements(empty_bev, vis_b; filename = "ae_empty"))
        @test isfile(plot_binary_period_distribution([empty_bev], vis_b; filename = "period_empty"))
        @test_throws ErrorException plot_binary_period_distribution(
            BinaryEvolutionSnapshot[],
            vis_b,
        )
        rm(vis_b.output_dir; recursive = true, force = true)
    end

    # =====================================================================
    @testset "Parameter sweeps" begin
        work = mktempdir()
        fixtures = joinpath(@__DIR__, "fixtures")
        base_pipeline = joinpath(work, "config.toml")
        write(
            base_pipeline,
            """
[install]
enabled = false
install_dir = "backend/Nbody6PPGPU-beijing"

[simulation]
run_test = true
input_file = "../../input_files/N1k_quick.inp"
runs_dir = "runs"
omp_threads = 8
run_id_prefix = "sw"

[postprocess]
enabled = true

[visualization]
enabled = false
format = "png"
column = "double"
output_dir = "plots"

[merger]
enabled = false
config_file = ""
""",
        )
        base_merger = joinpath(work, "merger.toml")
        write(
            base_merger,
            replace(
                read(joinpath(@__DIR__, "..", "input_files", "merger_demo_small.toml"), String),
                "[merger]\n" => "[merger]\nseed = 7\n";
                count = 1,
            ),
        )
        sweep_toml(extra = "") = """
[sweep]
name = "unit"
pipeline_config = "config.toml"
merger_config = "merger.toml"
seeds = [1, 2]
concurrency = 2
omp_threads = 3
$(extra)
[sweep.grid]
"merger.orbit.eccentricity" = [0.0, 0.5]
"merger.cluster2.N" = [300, 400]
"""
        spath = joinpath(work, "sweep.toml")
        write(spath, sweep_toml())
        scfg = load_sweep_config(spath)
        @test scfg.name == "unit"
        @test scfg.pipeline_config == base_pipeline && scfg.merger_config == base_merger
        @test first.(scfg.grid) == ["merger.cluster2.N", "merger.orbit.eccentricity"]   # key order
        @test last.(scfg.grid) == [[300, 400], [0.0, 0.5]]
        @test scfg.seeds == [1, 2] && scfg.concurrency == 2 && scfg.omp_threads == 3
        @test scfg.runs_dir == joinpath(work, "runs") && scfg.poll_interval == 2.0

        pts = sweep_points(scfg)
        @test length(pts) == 8                       # 2 × 2 grid × 2 seeds
        @test [p.index for p in pts] == 1:8
        @test pts[1].id == "001_cluster2-N=300_orbit-eccentricity=0_seed=1"
        @test pts[2].id == "002_cluster2-N=300_orbit-eccentricity=0_seed=2"
        @test pts[3].id == "003_cluster2-N=400_orbit-eccentricity=0_seed=1"   # first axis fastest
        @test pts[5].id == "005_cluster2-N=300_orbit-eccentricity=0.5_seed=1"
        @test pts[8].values == ["merger.cluster2.N" => 400, "merger.orbit.eccentricity" => 0.5]
        @test pts[8].seed == 2
        @test Nbody6Dynamics._axis_short("merger.cluster1.binaries.fraction") ==
              "cluster1-binaries-fraction"
        @test Nbody6Dynamics._format_axis_value(0.25) == "0.25"
        @test Nbody6Dynamics._format_axis_value(1000) == "1000"
        @test Nbody6Dynamics._format_axis_value(true) == "true"
        @test Nbody6Dynamics._format_axis_value("king W0=6") == "king-W0-6"

        d = Dict{String,Any}(
            "merger" => Dict{String,Any}("orbit" => Dict{String,Any}("apocentre" => 1.0)),
        )
        Nbody6Dynamics._set_nested!(d, "merger.orbit.eccentricity", 0.3)
        @test d["merger"]["orbit"]["eccentricity"] == 0.3
        @test_throws ArgumentError Nbody6Dynamics._set_nested!(d, "merger.tidal.kz14", 2)

        # Validation
        bad(extra_or_text) =
            (write(spath, extra_or_text); @test_throws ErrorException load_sweep_config(spath))
        bad(replace(sweep_toml(), "name = \"unit\"" => "name = \"unit sweep\""))
        bad(replace(sweep_toml(), "seeds = [1, 2]" => "seeds = [1, 1]"))
        bad(replace(sweep_toml(), "seeds = [1, 2]" => "seeds = []"))
        bad(replace(sweep_toml(), "concurrency = 2" => "concurrency = 0"))
        bad(replace(sweep_toml(), "omp_threads = 3" => "omp_threads = 0"))
        bad(
            replace(
                sweep_toml(),
                "\"merger.cluster2.N\" = [300, 400]\n" => "\"merger.seed\" = [3]\n",
            ),
        )
        bad(
            replace(
                sweep_toml(),
                "\"merger.cluster2.N\" = [300, 400]\n" => "\"orbit.apocentre\" = [3.0]\n",
            ),
        )
        bad(
            replace(
                sweep_toml(),
                "\"merger.cluster2.N\" = [300, 400]\n" => "\"merger.cluster2.N\" = []\n",
            ),
        )
        bad(
            replace(
                sweep_toml(),
                "merger_config = \"merger.toml\"" => "merger_config = \"absent.toml\"",
            ),
        )
        # An absent grid is a seed ensemble of the base configuration, not an error
        write(spath, replace(sweep_toml(), r"\[sweep\.grid\][\s\S]*" => ""))
        @test isempty(load_sweep_config(spath).grid)
        write(spath, sweep_toml("poll_interval = 0.5\n"))
        scfg = load_sweep_config(spath)
        @test scfg.poll_interval == 0.5

        # Preparation (dry run): derived configs, index, summary, figures
        sdir = joinpath(work, "sweep_unit")
        @test run_sweep(scfg; dry_run = true, sweep_dir = sdir) == sdir
        @test_throws ErrorException prepare_sweep(scfg; sweep_dir = sdir)   # never reuse a sweep dir
        @test isfile(joinpath(sdir, "base_config.toml")) &&
              isfile(joinpath(sdir, "base_merger.toml"))
        m3 = load_merger_config(joinpath(sdir, pts[3].id, "merger.toml"))
        @test m3.clusters[2].N == 400 && m3.orbit.eccentricity == 0.0 && m3.seed == 1
        m8 = load_merger_config(joinpath(sdir, pts[8].id, "merger.toml"))
        @test m8.clusters[2].N == 400 && m8.orbit.eccentricity == 0.5 && m8.seed == 2
        c3 = load_config(joinpath(sdir, pts[3].id, "config.toml"))
        @test !c3.install.enabled
        @test c3.install.install_dir == joinpath(work, "backend", "Nbody6PPGPU-beijing")
        @test c3.simulation.runs_dir == joinpath(sdir, pts[3].id)
        @test c3.simulation.omp_threads == 3 && c3.simulation.run_test && !c3.simulation.monitor
        @test c3.merger.enabled && c3.merger.config_file == joinpath(sdir, pts[3].id, "merger.toml")
        @test c3.visualization.column == "double"     # untouched base settings survive
        idx = read_sweep_index(sdir)
        @test idx["sweep"]["name"] == "unit" && idx["sweep"]["n_points"] == 8
        @test idx["sweep"]["axes"] == ["merger.cluster2.N", "merger.orbit.eccentricity"]
        @test all(p -> p["status"] == "pending", idx["points"])
        @test idx["points"][5]["values"]["merger.orbit.eccentricity"] == 0.5
        @test idx["points"][5]["dir"] == joinpath(sdir, pts[5].id)

        # Mark two points done with real output excerpts and summarise
        status = Dict{Int,Dict{String,Any}}()
        for i in (1, 3)
            run_out = joinpath(sdir, pts[i].id, "run", "output")
            mkpath(run_out)
            cp(joinpath(fixtures, "out1000"), joinpath(run_out, "out1000"))
            cp(joinpath(fixtures, "lagr.7"), joinpath(run_out, "lagr.7"))
            write(
                joinpath(sdir, pts[i].id, "run", "RUN_INFO.toml"),
                "[run]\nelapsed_seconds = 12.5\n\n[[segments]]\nexit_status = 0\n",
            )
            status[i] =
                Dict{String,Any}("status" => "done", "exit_status" => 0, "elapsed_seconds" => 12.5)
        end
        status[2] =
            Dict{String,Any}("status" => "failed", "exit_status" => 1, "elapsed_seconds" => 3.0)
        write_sweep_index(sdir, scfg, pts; status = status)
        columns, rows = sweep_summary(sdir)
        @test columns[1:5] == ["index", "id", "kind", "control_of", "seed"]
        @test columns[6:7] == idx["sweep"]["axes"]
        @test length(rows) == 8
        @test rows[1]["status"] == "done" &&
              rows[1]["elapsed_seconds"] == 12.5 &&
              rows[1]["exit_status"] == 0
        @test rows[1]["n_final"] > 0 &&
              isfinite(rows[1]["de_final"]) &&
              isfinite(rows[1]["q_final"])
        @test rows[1]["t_final_myr"] > 0
        @test rows[2]["status"] == "failed" && rows[2]["exit_status"] == 1
        @test rows[4]["status"] == "pending" &&
              isnan(rows[4]["de_final"]) &&
              rows[4]["n_final"] == -1
        csv_path = write_sweep_summary(sdir)
        lines = readlines(csv_path)
        @test length(lines) == 9
        @test lines[1] == join(columns, ",")
        @test startswith(
            lines[2],
            "1,001_cluster2-N=300_orbit-eccentricity=0_seed=1,merger,,1,300,0,done,0,true,12.5,",
        )
        @test occursin(",pending,-1,false,,,-1,-1,,", lines[5])   # NaN → empty cells

        vis_sw = sweep_visualization(scfg, sdir)
        @test vis_sw.output_dir == joinpath(sdir, "plots") &&
              vis_sw.format == "png" &&
              vis_sw.column == "double"
        figs = sweep_figures(sdir, vis_sw)
        @test length(figs) == 4 && all(isfile, figs)   # two seeds → comparison and ensemble figures
        @test isfile(
            plot_sweep_lagrangian(
                sdir,
                vis_sw;
                axis = "merger.orbit.eccentricity",
                fraction = 0.1,
                filename = "r10",
            ),
        )
        @test_throws ArgumentError plot_sweep_energy(sdir, vis_sw; axis = "merger.orbit.apocentre")
        rm(work; recursive = true, force = true)

        # Annotation corner: the least occupied of the data's bounding box
        xs = collect(0.0:0.1:1.0)
        @test Nbody6Dynamics._emptiest_corner(xs, xs) == :tl          # rising series frees the top-left
        @test Nbody6Dynamics._emptiest_corner(xs, 1 .- xs) == :tr     # falling series frees the top-right
        @test Nbody6Dynamics._emptiest_corner(Float64[], Float64[]) == :tl
    end

    # =====================================================================
    @testset "Seeded ensembles" begin
        # Type-7 quantiles and linear interpolation
        q = Nbody6Dynamics._quantile_sorted
        v = [1.0, 2.0, 3.0, 4.0, 5.0]
        @test q(v, 0.5) == 3.0 && q(v, 0.0) == 1.0 && q(v, 1.0) == 5.0
        @test q(v, 0.16) ≈ 1.64 && q(v, 0.84) ≈ 4.36
        @test q(v, 0.025) ≈ 1.1 && q(v, 0.975) ≈ 4.9
        @test q([7.0], 0.3) == 7.0
        @test_throws ArgumentError q(Float64[], 0.5)
        @test_throws ArgumentError q(v, 1.5)
        interp = Nbody6Dynamics._interpolate_linear
        @test interp([0.0, 1.0, 2.0], [0.0, 10.0, 0.0], [0.0, 0.5, 1.0, 1.5, 2.0]) ==
              [0.0, 5.0, 10.0, 5.0, 0.0]
        @test interp([1.0], [3.0], [1.0]) == [3.0]
        @test_throws ArgumentError interp([0.0, 1.0], [0.0, 1.0], [1.5])
        @test_throws DimensionMismatch interp([0.0, 1.0], [0.0], [0.5])

        # Members y_i = i t on different spans: common grid = intersection
        members = [(collect(0.0:0.5:10.0), i .* collect(0.0:0.5:10.0)) for i in 1:5]
        members[2] = (collect(1.0:0.5:12.0), 2 .* collect(1.0:0.5:12.0))   # starts later, ends later
        st = ensemble_statistics(members; n_grid = 19)
        @test st.n == 5
        @test st.time[1] == 1.0 && st.time[end] == 10.0 && length(st.time) == 19
        @test st.median ≈ 3 .* st.time
        @test st.q16 ≈ 1.64 .* st.time && st.q84 ≈ 4.36 .* st.time
        @test st.q025 ≈ 1.1 .* st.time && st.q975 ≈ 4.9 .* st.time
        single = ensemble_statistics(members[1:1]; n_grid = 5)
        @test single.n == 1 &&
              single.median == single.q025 == single.q975 ≈ collect(range(0.0, 10.0; length = 5))
        @test_throws ArgumentError ensemble_statistics(typeof(members[1])[])
        @test_throws ArgumentError ensemble_statistics([
            members[1],
            (collect(20.0:1.0:30.0), zeros(11)),
        ])
        @test_throws ArgumentError ensemble_statistics([([1.0, 0.5, 2.0], [0.0, 1.0, 2.0])])
        @test_throws ArgumentError ensemble_statistics(members; n_grid = 1)

        # Series extractors on a run directory built from the fixtures
        fixtures = joinpath(@__DIR__, "fixtures")
        rdir = mktempdir()
        @test Nbody6Dynamics._run_series(rdir, :lagrangian) === nothing
        @test Nbody6Dynamics._run_series(rdir, :energy) === nothing
        mkpath(joinpath(rdir, "output"))
        cp(joinpath(fixtures, "out1000"), joinpath(rdir, "output", "out1000"))
        cp(joinpath(fixtures, "lagr.7"), joinpath(rdir, "output", "lagr.7"))
        t_l, r_l = Nbody6Dynamics._run_series(rdir, :lagrangian; fraction = 0.5)
        @test length(t_l) == length(r_l) > 1 && issorted(t_l) && all(>(0), r_l)
        t_e, de = Nbody6Dynamics._run_series(rdir, :energy)
        @test all(>(0), de) && length(t_e) == length(de)
        t_n, n = Nbody6Dynamics._run_series(rdir, :n_stars)
        @test all(>(0), n)
        t_p, np = Nbody6Dynamics._run_series(rdir, :n_pairs)
        @test length(np) == length(t_n) && all(≥(0), np)
        @test_throws ArgumentError Nbody6Dynamics._run_series(rdir, :unknown)
        @test Nbody6Dynamics._series_label(:energy, 0.5) isa AbstractString
        @test_throws ArgumentError Nbody6Dynamics._series_label(:unknown, 0.5)

        # A seeds-only sweep (no grid axes) and a one-axis sweep, marked done with fixture output
        work = mktempdir()
        write(
            joinpath(work, "config.toml"),
            "[install]\nenabled = false\ninstall_dir = \"backend\"\n\n[simulation]\nrun_test = true\nruns_dir = \"runs\"\n\n[visualization]\nformat = \"png\"\n",
        )
        write(
            joinpath(work, "merger.toml"),
            read(joinpath(@__DIR__, "..", "input_files", "merger_demo_small.toml"), String),
        )
        function _done_sweep(name, grid_text, seeds_text)
            spath = joinpath(work, "sweep_$(name).toml")
            write(
                spath,
                "[sweep]\nname = \"$(name)\"\npipeline_config = \"config.toml\"\nmerger_config = \"merger.toml\"\nseeds = $(seeds_text)\n$(grid_text)",
            )
            scfg = load_sweep_config(spath)
            sdir = joinpath(work, "sweep_$(name)")
            prepare_sweep(scfg; sweep_dir = sdir)
            pts = sweep_points(scfg)
            status = Dict{Int,Dict{String,Any}}()
            for p in pts
                out = joinpath(sdir, p.id, "run", "output")
                mkpath(out)
                cp(joinpath(fixtures, "out1000"), joinpath(out, "out1000"))
                cp(joinpath(fixtures, "lagr.7"), joinpath(out, "lagr.7"))
                status[p.index] = Dict{String,Any}(
                    "status" => "done",
                    "exit_status" => 0,
                    "elapsed_seconds" => 1.0,
                )
            end
            write_sweep_index(sdir, scfg, pts; status = status)
            return scfg, sdir, pts
        end
        scfg0, sdir0, pts0 = _done_sweep("seeds", "", "[1, 2, 3]")
        @test isempty(scfg0.grid)
        @test [p.id for p in pts0] == ["001_seed=1", "002_seed=2", "003_seed=3"]
        ens0 = sweep_ensembles(sdir0, :lagrangian)
        @test length(ens0) == 1 && ens0[1].stats.n == 3 && isempty(ens0[1].values)
        @test ens0[1].stats.median ≈ ens0[1].stats.q975    # identical members → zero-width bands
        vis_e = VisualizationConfig(;
            format = "png",
            column = "single",
            output_dir = joinpath(work, "plots"),
        )
        @test isfile(plot_sweep_ensemble(sdir0, vis_e))
        @test isfile(plot_sweep_ensemble(sdir0, vis_e; quantity = :energy))
        @test isfile(plot_sweep_ensemble(sdir0, vis_e; quantity = :n_stars, filename = "ens_n"))
        @test isfile(plot_sweep_lagrangian(sdir0, vis_e; filename = "seeds_lagr"))   # no axes: one colour, no legend
        @test_throws ArgumentError plot_sweep_lagrangian(
            sdir0,
            vis_e;
            axis = "merger.orbit.eccentricity",
        )
        figs0 = sweep_figures(sdir0, vis_e)
        @test length(figs0) == 4 && all(isfile, figs0)

        scfg1, sdir1, pts1 = _done_sweep(
            "axes",
            "[sweep.grid]\n\"merger.orbit.eccentricity\" = [0.0, 0.5]\n\"merger.cluster2.N\" = [300, 400]\n",
            "[1, 2]",
        )
        ens1 = sweep_ensembles(sdir1, :energy)
        @test length(ens1) == 4 && all(e -> e.stats.n == 2, ens1)
        @test isfile(
            plot_sweep_ensemble(sdir1, vis_e; quantity = :lagrangian, filename = "ens_axes"),
        )
        @test isfile(
            plot_sweep_ensemble(
                sdir1,
                vis_e;
                axis = "merger.orbit.eccentricity",
                fixed = Dict("merger.cluster2.N" => 400),
                filename = "ens_fixed",
            ),
        )
        @test_throws ErrorException plot_sweep_ensemble(
            sdir1,
            vis_e;
            axis = "merger.orbit.eccentricity",
            fixed = Dict("merger.cluster2.N" => 999),
        )
        @test_throws ArgumentError plot_sweep_ensemble(
            sdir1,
            vis_e;
            fixed = Dict("merger.cluster2.N" => 400),
        )   # the axis itself
        @test_throws ArgumentError plot_sweep_ensemble(
            sdir1,
            vis_e;
            fixed = Dict("merger.orbit.apocentre" => 1.0),
        )
        @test length(sweep_figures(sdir1, vis_e)) == 4
        rm(work; recursive = true, force = true)
        rm(rdir; recursive = true, force = true)
    end

    # =====================================================================
    @testset "Remnant diagnostics" begin
        using StableRNGs
        rng_r = StableRNG(7)
        N = 2000
        pos, vel = sample_plummer(N, 1.0; rng = rng_r)
        mass = fill(1.0 / N, N)
        virialise!(mass, pos, vel)

        # Core radius: with exact Plummer densities the particle-weighted
        # Casertano–Hut estimator gives √0.3 a; the 6-neighbour estimate
        # sits a few per cent below.
        r_all = vec(sqrt.(sum(pos .^ 2; dims = 1)))
        ρ_exact = (1 .+ r_all .^ 2) .^ (-2.5)
        cr = core_radius(pos, ρ_exact)
        @test cr.r_core ≈ sqrt(0.3) rtol = 0.06
        @test all(abs.(cr.centre) .< 0.1)
        ρ6 = Nbody6Dynamics._local_density(pos, mass)
        @test length(ρ6) == N && all(≥(0), ρ6)
        @test 0.42 < core_radius(pos, ρ6).r_core < 0.62
        @test_throws ArgumentError Nbody6Dynamics._local_density(pos[:, 1:5], mass[1:5])
        @test_throws DimensionMismatch core_radius(pos, ρ6[1:10])
        @test_throws ArgumentError core_radius(pos, zeros(N))

        # Rotation: an isotropic sphere has λ_R at the noise level; adding
        # solid-body rotation about z recovers the axis and raises λ_R.
        rot0 = rotation_analysis(pos, vel, mass)
        @test rot0.lambda_r < 0.15
        @test rot0.lambda_peebles < 0.05
        @test length(rot0.profile.radius) == 8 && issorted(rot0.profile.radius)
        vrot = copy(vel)
        Ω = 0.5
        for i in 1:N
            vrot[1, i] -= Ω * pos[2, i]
            vrot[2, i] += Ω * pos[1, i]
        end
        rot1 = rotation_analysis(pos, vrot, mass)
        @test rot1.lambda_r > 0.6
        @test abs(rot1.axis[3]) > 0.99
        @test rot1.axis[3] > 0                          # L along +z for counter-clockwise motion
        @test all(>(0), rot1.profile.v_rot)              # v_φ > 0 in every shell
        @test rot1.profile.v_rot_over_sigma[end - 1] > rot1.profile.v_rot_over_sigma[1]   # Ω R grows outward
        @test rotation_analysis(pos, vrot, mass; nbins = 3).profile.radius |> length == 3
        @test_throws ArgumentError rotation_analysis(pos[:, 1:5], vrot[:, 1:5], mass[1:5])

        # Mass segregation: random masses give Λ ≈ 1; the heaviest stars
        # placed innermost give Λ ≫ 1 and a small half-mass ratio.
        m_rand = mass .* (1 .+ rand(rng_r, N))
        ms0 = mass_segregation(pos, m_rand; seed = 1)
        @test abs(ms0.lambda_msr - 1) < 3 * ms0.lambda_err + 0.3
        @test 0.7 < ms0.r_half_ratio < 1.3
        m_seg = zeros(N)
        m_seg[sortperm(r_all)] = sort(1 .+ 9 .* rand(rng_r, N); rev = true)
        ms1 = mass_segregation(pos, m_seg; seed = 1)
        @test ms1.lambda_msr > 5
        @test ms1.r_half_ratio < 0.4
        @test mass_segregation(pos, m_seg; seed = 1).lambda_msr == ms1.lambda_msr   # reproducible
        @test_throws ArgumentError mass_segregation(pos[:, 1:20], m_seg[1:20])       # ≤ n_massive
        @test_throws ArgumentError mass_segregation(pos, m_seg; n_massive = 1)
        @test Nbody6Dynamics._mst_length(pos, Int[]) == 0.0
        @test Nbody6Dynamics._mst_length([0.0 3.0 3.0; 0.0 0.0 4.0; 0.0 0.0 0.0], [1, 2, 3]) ≈ 7.0   # 3 + 4

        # Coalescence and the full time series on two synthetic snapshots:
        # separated clusters, then superposed with their orbital velocities.
        p1, v1 = sample_plummer(600, 0.3; rng = rng_r)
        m1 = fill(1.0 / 1200, 600)
        virialise!(m1, p1, v1)
        p2, v2 = sample_plummer(600, 0.3; rng = rng_r)
        m2 = fill(1.0 / 1200, 600)
        virialise!(m2, p2, v2)
        pos_o, vel_o, mass_o, ranges_o, _ = Nbody6Dynamics.setup_two_cluster_orbit(
            p1,
            v1,
            m1,
            p2,
            v2,
            m2,
            6.0,
            0.0;
            truncate_jacobi_flag = false,
        )
        function _rsnap(pos, vel, mass, t)
            params = zeros(Float32, 20)
            params[1] = t
            params[3] = 2.0f0     # rbar: 1 NB = 2 pc
            params[4] = 1.0f0
            params[11] = 3.0f0    # tscale: 1 NB = 3 Myr
            params[12] = 1.0f0
            Snapshot(
                SnapshotHeader(Int32(length(mass)), Int32(1), Int32(1), Int32(20), params),
                Int32.(1:length(mass)),
                Float32.(mass),
                Float32.(pos),
                Float32.(vel),
                Float32[],
                Float32[],
            )
        end
        s0 = _rsnap(pos_o, vel_o, mass_o, 0.0)
        # "Merged": the clusters overlap (0.5 apart, radii ≈ 0.35) and keep
        # their orbital velocities, so the pair's angular momentum survives
        # with the orbital sense; superposing them exactly would cancel it.
        pos_m = copy(pos_o)
        pos_m[:, ranges_o[1]] .= p1 .- [0.25, 0.0, 0.0]
        pos_m[:, ranges_o[2]] .= p2 .+ [0.25, 0.0, 0.0]
        s1 = _rsnap(pos_m, vel_o, mass_o, 1.0)
        coal = coalescence_time([s0, s1], ranges_o)
        @test coal.index == 2 && coal.time_nb == 1.0 && coal.time_myr == 3.0
        @test coal.n_distinct == [2, 1]
        far = coalescence_time([s0, s0], ranges_o)
        @test far.index === nothing && isnan(far.time_nb) && far.n_distinct == [2, 2]
        @test_throws ArgumentError coalescence_time([s0], ranges_o[1:1])
        L_orb = Nbody6Dynamics._orbital_angular_momentum(s0, ranges_o)
        @test abs(L_orb[3]) > 0.5 && abs(L_orb[1]) < 1e-6 && abs(L_orb[2]) < 1e-6   # orbit in the x–y plane

        diag = remnant_diagnostics([s0, s1], ranges_o; n_massive = 20, n_random = 20)
        @test diag.time == [0.0, 1.0] && diag.time_myr == [0.0, 3.0] && diag.rbar == [2.0, 2.0]
        @test diag.coalescence_time == 1.0 && diag.coalescence_time_myr == 3.0
        @test all(diag.n_bound .≥ 1100)
        @test all(0.9 .< diag.bound_mass_fraction .≤ 1.0)
        @test diag.r_core[1] > diag.r_core[2] > 0           # two separated clusters → one compact remnant
        @test diag.r_half[1] > diag.r_half[2] > 0
        @test diag.spin_alignment[1] ≈ 1.0 atol = 1e-3       # spin of the pair = orbital L at t = 0
        @test diag.spin_alignment[2] > 0.9                   # remnant keeps the orbital sense
        @test all(isfinite, diag.lambda_r) && all(isfinite, diag.lambda_peebles)
        @test all(isfinite, diag.lambda_msr) && all(isfinite, diag.segregation_ratio)
        @test isnan(diag.segregation_time)
        @test !isempty(diag.profile.radius)
        @test_throws ArgumentError remnant_diagnostics(Snapshot[], ranges_o)

        csv_path = joinpath(mktempdir(), "remnant.csv")
        @test write_remnant_diagnostics(csv_path, diag) == csv_path
        lines = readlines(csv_path)
        @test length(lines) == 4
        @test startswith(lines[1], "# coalescence_time_nb=1.0 coalescence_time_myr=3.0")
        @test lines[2] ==
              "time_nb,time_myr,rbar_pc,n_bound,bound_mass_fraction,r_core,r_half,lambda_r,lambda_peebles,spin_x,spin_y,spin_z,spin_alignment,lambda_msr,lambda_msr_err,segregation_ratio"
        @test startswith(lines[4], "1,3,2,")

        vis_r = VisualizationConfig(; format = "png", column = "single", output_dir = mktempdir())
        figs = remnant_figures(diag, vis_r)
        @test length(figs) == 4 && all(isfile, figs)
        vis_nb = VisualizationConfig(;
            format = "png",
            column = "single",
            units = "nbody",
            output_dir = vis_r.output_dir,
        )
        @test isfile(plot_remnant_structure(diag, vis_nb; filename = "structure_nb"))
        @test isfile(
            plot_mass_segregation_evolution(
                diag,
                vis_nb;
                lambda_threshold = 1.5,
                filename = "seg_nb",
            ),
        )
        @test isfile(
            plot_rotation_profile(
                diag.profile,
                vis_r;
                rbar = 2.0,
                lambda_r = 0.3,
                filename = "profile",
            ),
        )
        @test_throws ErrorException plot_rotation_profile(
            RotationProfile(Float64[], Float64[], Float64[], Float64[]),
            vis_r,
        )
        rm(vis_r.output_dir; recursive = true, force = true)
    end

    # =====================================================================
    @testset "Control configurations" begin
        demo = joinpath(@__DIR__, "..", "input_files", "merger_demo_small.toml")
        raw = Nbody6Dynamics.TOML.parsefile(demo)
        raw["merger"]["cluster2"]["N"] = 500
        raw["merger"]["cluster2"]["rbar"] = 4.0
        raw["merger"]["cluster1"]["binaries"] = Dict{String,Any}("fraction" => 0.1)
        raw["merger"]["nbody6"] = Dict{String,Any}("qe" => 0.01)
        c = control_merger_dict(raw)
        m = c["merger"]
        @test m["n_clusters"] == 1 && m["orbit_mode"] == "explicit"
        @test !haskey(m, "orbit") && !haskey(m, "cluster2")
        @test m["cluster1"]["N"] == 1500                                  # 1000 + 500
        @test m["cluster1"]["rbar"] ≈ (1000 * 2.0 + 500 * 4.0) / 1500     # N-weighted r_h
        @test m["cluster1"]["position"] == [0.0, 0.0, 0.0] &&
              m["cluster1"]["velocity"] == [0.0, 0.0, 0.0]
        @test m["cluster1"]["model"] == "king" && m["cluster1"]["W0"] == 6.0
        @test m["cluster1"]["binaries"]["fraction"] == 0.1                # cluster 1's population kept
        @test m["nbody6"]["qe"] == 0.01 && haskey(m, "output")             # other sections verbatim
        # NB intervals scale by (RBAR_est / r_h,c)^{3/2}: apocentre 12 × 500/1500 + N-weighted r_h
        rbar_c = m["cluster1"]["rbar"]
        factor = ((12.0 * 500 / 1500 + rbar_c) / rbar_c)^1.5
        @test m["output"]["tcrit"] ≈ raw["merger"]["output"]["tcrit"] * factor
        @test m["output"]["deltat"] ≈ raw["merger"]["output"]["deltat"] * factor
        raw_myr = deepcopy(raw)
        delete!(raw_myr["merger"]["output"], "tcrit")
        raw_myr["merger"]["output"]["tcrit_myr"] = 25.0
        raw_myr["merger"]["stellar"] = Dict{String,Any}("dtplot_myr" => 5.0)
        m_myr = control_merger_dict(raw_myr)["merger"]
        @test m_myr["output"]["tcrit_myr"] == 25.0 && !haskey(m_myr["output"], "tcrit")   # physical: verbatim
        @test m_myr["output"]["deltat"] ≈ raw["merger"]["output"]["deltat"] * factor     # NB: scaled
        @test m_myr["stellar"]["dtplot_myr"] == 5.0
        raw_ex = deepcopy(raw)
        raw_ex["merger"]["orbit_mode"] = "explicit"
        delete!(raw_ex["merger"], "orbit")
        raw_ex["merger"]["cluster1"]["position"] = [-6.0, 0.0, 0.0]
        raw_ex["merger"]["cluster2"]["position"] = [12.0, 0.0, 0.0]
        f_ex =
            Nbody6Dynamics._control_time_factor(raw_ex["merger"], [1000, 500], [2.0, 4.0], rbar_c)
        centre = (1000 * -6.0 + 500 * 12.0) / 1500
        rms = sqrt((1000 * (-6.0 - centre)^2 + 500 * (12.0 - centre)^2) / 1500)
        @test f_ex ≈ ((rms + rbar_c) / rbar_c)^1.5
        @test raw["merger"]["n_clusters"] == 2                            # input untouched
        @test_throws ArgumentError control_merger_dict(Dict{String,Any}())
        bad = deepcopy(raw)
        delete!(bad["merger"], "cluster2")
        @test_throws ArgumentError control_merger_dict(bad)

        work = mktempdir()
        dst = joinpath(work, "control.toml")
        @test write_control_merger_config(demo, dst) == dst
        ccfg = load_merger_config(dst)
        @test length(ccfg.clusters) == 1 &&
              ccfg.clusters[1].N == 2000 &&
              ccfg.orbit_mode == "explicit"
        @test_throws ErrorException write_control_merger_config(joinpath(work, "absent.toml"), dst)

        # The single-cluster generator path: one member block, no truncation, at rest
        small = replace(read(dst, String), r"N = \d+" => "N = 300")
        write(joinpath(work, "control_small.toml"), small)
        res = generate_merger_ic(
            load_merger_config(joinpath(work, "control_small.toml"));
            output_dir = joinpath(work, "ic"),
        )
        @test res.N_total == 300 &&
              length(res.cluster_ranges) == 1 &&
              length(res.cluster_ranges[1]) == 300
        @test isfile(joinpath(work, "ic", "dat.10")) && isfile(joinpath(work, "ic", "merger.inp"))
        @test parse_merger_summary(joinpath(work, "ic", "merger_summary.txt")) == [collect(1:300)]
        # Physical-time intervals: converted with the realised T* at generation
        phys = replace(small, r"tcrit = [0-9.]+" => "tcrit_myr = 20.0")
        phys = replace(phys, r"deltat = [0-9.]+" => "deltat_myr = 4.0")
        write(joinpath(work, "control_phys.toml"), phys)
        cfg_p = load_merger_config(joinpath(work, "control_phys.toml"))
        @test cfg_p.output.tcrit_myr == 20.0 && cfg_p.output.deltat_myr == 4.0
        res_p = generate_merger_ic(cfg_p; output_dir = joinpath(work, "ic_phys"))
        meta_p = Nbody6Dynamics.TOML.parsefile(joinpath(work, "ic_phys", "merger_ic.toml"))
        t_star = meta_p["meta"]["t_star_myr"]
        @test t_star > 0
        @test meta_p["output"]["tcrit"] ≈ 20.0 / t_star && meta_p["output"]["tcrit_myr"] ≈ 20.0
        # The interval is the dyadic rounding of the converted value (engine digit counter)
        @test meta_p["output"]["deltat"] == engine_interval(4.0 / t_star)
        @test abs(meta_p["output"]["deltat"] - 4.0 / t_star) ≤ 4.0 / t_star / 256
        inp_p = read(joinpath(work, "ic_phys", "merger.inp"), String)
        @test occursin(Nbody6Dynamics.Printf.@sprintf("TCRIT=%.2f", 20.0 / t_star), inp_p)
        @test occursin(
            "Time unit: T* =",
            read(joinpath(work, "ic_phys", "merger_summary.txt"), String),
        )
        both = replace(small, r"tcrit = [0-9.]+" => "tcrit = 1.0\ntcrit_myr = 20.0")
        write(joinpath(work, "both.toml"), both)
        @test_throws ErrorException load_merger_config(joinpath(work, "both.toml"))
        neg = replace(small, r"tcrit = [0-9.]+" => "tcrit_myr = -1.0")
        write(joinpath(work, "neg.toml"), neg)
        @test_throws ErrorException load_merger_config(joinpath(work, "neg.toml"))
        # dtplot below deltat is caught at generation once T* is known
        # (the control file carries a [merger.stellar] table with the scaled dtplot)
        @test occursin(r"dtplot = [0-9.]+", phys)
        late = replace(phys, r"dtplot = [0-9.]+" => "dtplot_myr = 1.0")
        write(joinpath(work, "late.toml"), late)
        @test_throws ErrorException load_merger_config(joinpath(work, "late.toml"))   # both physical: config time
        late_nb = replace(phys, r"dtplot = [0-9.]+" => "dtplot = 0.01")
        write(joinpath(work, "late_nb.toml"), late_nb)
        @test_throws ErrorException generate_merger_ic(
            load_merger_config(joinpath(work, "late_nb.toml"));
            output_dir = joinpath(work, "ic_late"),
        )
        kepler_one = replace(small, "orbit_mode = \"explicit\"" => "orbit_mode = \"kepler\"")
        write(joinpath(work, "kepler_one.toml"), kepler_one)
        @test_throws ErrorException load_merger_config(joinpath(work, "kepler_one.toml"))

        # Sweep with controls: companions interleaved, derived at the point's values
        write(
            joinpath(work, "config.toml"),
            "[install]\nenabled = false\ninstall_dir = \"backend\"\n\n[simulation]\nrun_test = true\nruns_dir = \"runs\"\n\n[visualization]\nformat = \"png\"\n",
        )
        cp(demo, joinpath(work, "merger.toml"))
        write(
            joinpath(work, "sweep.toml"),
            "[sweep]\nname = \"ctrl\"\npipeline_config = \"config.toml\"\nmerger_config = \"merger.toml\"\nseeds = [1]\ncontrols = true\n\n[sweep.grid]\n\"merger.cluster2.N\" = [300, 400]\n",
        )
        scfg = load_sweep_config(joinpath(work, "sweep.toml"))
        @test scfg.controls
        pts = sweep_points(scfg)
        @test length(pts) == 4
        @test [p.kind for p in pts] == ["merger", "control", "merger", "control"]
        @test pts[2].id == "002_cluster2-N=300_seed=1_control" && pts[2].control_of == pts[1].id
        @test pts[2].values == pts[1].values && pts[2].seed == pts[1].seed
        sdir = joinpath(work, "sweep_ctrl")
        prepare_sweep(scfg; sweep_dir = sdir)
        m3 = load_merger_config(joinpath(sdir, pts[3].id, "merger.toml"))
        m4 = load_merger_config(joinpath(sdir, pts[4].id, "merger.toml"))
        @test length(m3.clusters) == 2 && m3.clusters[2].N == 400
        @test length(m4.clusters) == 1 && m4.clusters[1].N == 1400 && m4.seed == 1
        idx = read_sweep_index(sdir)
        @test idx["sweep"]["controls"] == true
        @test [p["kind"] for p in idx["points"]] == ["merger", "control", "merger", "control"]
        @test idx["points"][4]["control_of"] == pts[3].id
        columns, rows = sweep_summary(sdir)
        @test columns[1:5] == ["index", "id", "kind", "control_of", "seed"]
        @test rows[2]["kind"] == "control" && rows[2]["control_of"] == pts[1].id
        @test_throws ErrorException load_sweep_config(
            (
                write(
                    joinpath(work, "bad.toml"),
                    replace(
                        read(joinpath(work, "sweep.toml"), String),
                        "controls = true" => "controls = \"yes\"",
                    ),
                );
                joinpath(work, "bad.toml")
            ),
        )

        # Paired figure on fixture output: mergers solid, controls dashed
        fixtures = joinpath(@__DIR__, "fixtures")
        status = Dict{Int,Dict{String,Any}}()
        for p in pts
            out = joinpath(sdir, p.id, "run", "output")
            mkpath(out)
            cp(joinpath(fixtures, "out1000"), joinpath(out, "out1000"))
            cp(joinpath(fixtures, "lagr.7"), joinpath(out, "lagr.7"))
            status[p.index] =
                Dict{String,Any}("status" => "done", "exit_status" => 0, "elapsed_seconds" => 1.0)
        end
        write_sweep_index(sdir, scfg, pts; status = status)
        _, done_m = Nbody6Dynamics._sweep_done_points(sdir)
        _, done_c = Nbody6Dynamics._sweep_done_points(sdir; kind = "control")
        _, done_all = Nbody6Dynamics._sweep_done_points(sdir; kind = "")
        @test length(done_m) == 2 && length(done_c) == 2 && length(done_all) == 4
        vis_c = VisualizationConfig(;
            format = "png",
            column = "single",
            output_dir = joinpath(work, "plots"),
        )
        @test isfile(plot_control_comparison(sdir, vis_c))
        @test isfile(
            plot_control_comparison(sdir, vis_c; quantity = :energy, filename = "ctrl_energy"),
        )
        figs = sweep_figures(sdir, vis_c)
        @test length(figs) == 3 && all(isfile, figs)      # one seed: comparison + control figures
        @test length(sweep_ensembles(sdir, :lagrangian)) == 2   # mergers only
        status[2]["status"] = "failed"
        status[4]["status"] = "failed"
        write_sweep_index(sdir, scfg, pts; status = status)
        @test_throws ErrorException plot_control_comparison(sdir, vis_c; filename = "none")
        rm(work; recursive = true, force = true)
    end

    # =====================================================================
    @testset "Stellar type labels" begin
        # Standard Hurley et al. (2000) SSE/BSE table used by this fork
        @test startswith(STELLAR_TYPE_LABELS[0], "MS")
        @test startswith(STELLAR_TYPE_LABELS[1], "MS")
        @test startswith(STELLAR_TYPE_LABELS[10], "HeWD")
        @test startswith(STELLAR_TYPE_LABELS[12], "ONeWD")
        @test startswith(STELLAR_TYPE_LABELS[13], "NS")
        @test startswith(STELLAR_TYPE_LABELS[14], "BH")
        @test length(STELLAR_TYPE_LABELS) == 16   # K* = 0..15 complete
    end

    # =====================================================================
    # Regression tests against verbatim Nbody6PPGPU-beijing output.
    # The synthetic fixtures above mirror the readers' assumptions by
    # construction; these fixtures were copied from a real run and would
    # have caught the esc.11 column-map, M*-vs-<M>, and K*-label bugs.
    @testset "Real-output fixtures" begin
        FIXDIR = joinpath(@__DIR__, "fixtures")

        @testset "esc.11" begin
            escs = read_escapers(joinpath(FIXDIR, "esc.11"))
            @test length(escs) == 88          # every data line parses
            e1 = escs[1]
            @test e1.time_myr ≈ 8.11195 rtol = 1e-5     # T[Myr], not TTOT
            @test e1.mass_solar ≈ 18.1687 rtol = 1e-5   # M[M*], not NB mass
            @test e1.escape_energy ≈ 1017.22 rtol = 1e-5
            @test e1.velocity_kms ≈ 209.841 rtol = 1e-5
            @test e1.stellar_type == 14                 # BH escaper
            @test e1.name == 2248
            @test e1.phi_deg ≈ 22.2302 rtol = 1e-5      # ANGLE PHI [deg]
            @test e1.theta_deg ≈ 28.7010 rtol = 1e-5    # ANGLE THETA [deg]
            # escape.F angle conventions: φ ∈ [0, 360], θ ∈ [-90, 90]
            @test all(0.0 ≤ e.phi_deg ≤ 360.0 for e in escs)
            @test all(-90.0 ≤ e.theta_deg ≤ 90.0 for e in escs)
        end

        @testset "out1000 scaling + ADJUST" begin
            diag = read_diagnostics(joinpath(FIXDIR, "out1000"))
            u = extract_scaling(diag)
            @test u.rbar ≈ 5.503639 rtol = 1e-6
            @test u.zmbar ≈ 27695.934588 rtol = 1e-6   # M* (scale), not <M>
            @test u.tscale ≈ 1.15885054 rtol = 1e-6
            @test u.vstar ≈ 4.65224425 rtol = 1e-6
            @test length(diag.adjust) ≥ 6
            a0 = diag.adjust[1]
            @test a0.time_nb == 0.0
            @test a0.qvir ≈ 0.321 atol = 1e-3           # Q = T/|W|
            @test a0.e_tot ≈ -0.5408 atol = 1e-4
        end

        @testset "lagr.7" begin
            lagr = read_lagr(joinpath(FIXDIR, "lagr.7"))
            @test length(lagr.mass_fractions) == 18
            @test lagr.time[1] == 0.0
            i50 = argmin(abs.(lagr.mass_fractions .- 0.5))
            @test lagr.radii[i50, 1] ≈ 1.2442064 rtol = 1e-6
        end

        @testset "sev.83" begin
            sev = read_stellar_evolution(joinpath(FIXDIR, "sev.83_0"))
            @test sev.n_stars == 40
            @test sev.time_myr == 0.0
            @test length(sev.records) == 40
            r1 = sev.records[1]
            @test r1.time_nb == 0.0
            @test r1.name == Int32(1)
            @test r1.stellar_type == Int32(1)            # massive MS star, Hurley K*=1
            @test r1.ri ≈ 0.694596 rtol = 1e-5           # RI [pc]
            @test r1.mass_solar ≈ 52.5373 rtol = 1e-5
            @test r1.log_teff ≈ 4.71026 rtol = 1e-5
            @test r1.ms_lifetime_myr ≈ 4.46907 rtol = 1e-5   # TM
            @test r1.mass_core == 0.0                    # MC (MS star)
            @test r1.radius_envelope ≈ 1e-10 rtol = 1e-3 # RE placeholder on MS
        end

        @testset "bev.82" begin
            # Truncated excerpt (40 of 224 pairs) of the t = 0 record of a
            # two-cluster run with 20 % primordial binaries per cluster.
            bev = read_binary_evolution(joinpath(FIXDIR, "bev.82_0"))
            @test bev.n_pairs == 40
            @test bev.time_myr == 0.0
            @test length(bev.records) == 40
            r1 = bev.records[1]
            @test (r1.name1, r1.name2) == (1, 2)          # primordial pairs lead the body list
            @test r1.eccentricity ≈ 0.92188 rtol = 1e-5
            @test r1.log_period_days ≈ 1.20647 rtol = 1e-5
            @test r1.mass1 ≈ 0.221134 rtol = 1e-5
            @test r1.ms_lifetime1_myr ≈ 4.98658e5 rtol = 1e-5
            # Every record is a bound orbit whose own columns obey Kepler III.
            for r in bev.records
                P_yr = 10^r.log_period_days / 365.25
                a_au = semi_major_axis_pc(r) / Nbody6Dynamics._AU_IN_PC
                @test a_au^3 / P_yr^2 ≈ r.mass1 + r.mass2 rtol = 2e-2
                @test 0 ≤ r.eccentricity < 1
            end
            @test all(r -> r.name1 < r.name2, bev.records)
        end
    end

    # =====================================================================
    @testset "Types and accessors" begin
        params = zeros(Float32, 20)
        params[1] = 1.5f0   # time
        params[3] = 2.0f0   # rbar
        params[4] = 0.6f0   # zmbar
        params[11] = 8.0f0  # tscale
        params[12] = 4.0f0  # vstar

        hdr = SnapshotHeader(Int32(100), Int32(1), Int32(1), Int32(20), params)

        @test time_nb(hdr) ≈ 1.5
        @test rbar(hdr) ≈ 2.0
        @test tscale(hdr) ≈ 8.0
        @test time_myr(hdr) ≈ 12.0   # 1.5 * 8.0

        u = UnitScaling(2.0, 0.6, 8.0, 4.0)
        @test Nbody6Dynamics.to_pc(u, 1.0) ≈ 2.0
        @test Nbody6Dynamics.to_msun(u, 1.0) ≈ 0.6
        @test Nbody6Dynamics.to_myr(u, 1.0) ≈ 8.0
        @test Nbody6Dynamics.to_kms(u, 1.0) ≈ 4.0
    end

    # =====================================================================
    @testset "Plotting (smoke tests)" begin
        Nbody6Dynamics.set_publication_theme!()
        # The theme's faces must be live in this (fresh) process: a face
        # captured at precompile time has a null FreeType pointer and Makie
        # falls back to its sans default without warning.
        theme_fonts = publication_theme().fonts
        for key in (:regular, :bold, :italic)
            @test getfield(theme_fonts[key][], :ft_ptr) != C_NULL
        end
        @test theme_fonts[:regular][] === Nbody6Dynamics.texfont(:text)

        vis = VisualizationConfig(;
            output_dir = joinpath(TESTDIR, "test_plots"),
            format = "png",
            dpi = 72,
            figsize = (6, 4),
        )

        # --- Snapshot plots ---
        n = 50
        params = zeros(Float32, 20)
        params[1] = 1.0f0
        hdr = SnapshotHeader(Int32(n), Int32(1), Int32(1), Int32(20), params)
        snap = Snapshot(
            hdr,
            Int32.(1:n),
            Float32.(rand(n)),
            Float32.(randn(3, n)),
            Float32.(0.1 .* randn(3, n)),
            Float32[],
            Float32[],
        )

        plot_snapshot(snap, vis; filename = "test_snap")
        @test isfile(joinpath(TESTDIR, "test_plots", "test_snap_xy.png"))
        @test isfile(joinpath(TESTDIR, "test_plots", "test_snap_xz.png"))

        # --- Energy plot ---
        adj = [
            AdjustRecord(0.0, 0.0, 1.0, 0.0, -0.25, 100, 10, 1.0),
            AdjustRecord(0.5, 50.0, 0.98, 1e-6, -0.249, 99, 9, 1.1),
            AdjustRecord(1.0, 100.0, 0.99, 2e-6, -0.248, 98, 8, 1.2),
        ]
        diag = DiagnosticsData(adj, Dict{String,Float64}())
        plot_energy(diag, vis; filename = "test_energy")
        @test isfile(joinpath(TESTDIR, "test_plots", "test_energy.png"))

        # --- Lagrangian plot ---
        lagr = LagrangianData(
            [0.0, 0.5, 1.0],
            [0.1, 0.5, 1.0],
            [0.05 0.06 0.07; 0.2 0.22 0.25; 1.0 1.05 1.1],
        )
        plot_lagrangian(lagr, vis; filename = "test_lagr", selected_fractions = [0.1, 0.5, 1.0])
        @test isfile(joinpath(TESTDIR, "test_plots", "test_lagr.png"))

        # --- HR diagram ---
        sev = StellarEvolutionSnapshot(
            0.5,
            4,
            [
                StellarRecord(0.5, Int32(1), Int32(101), Int32(0), 1.0, 0.8, 0.1, -0.5, 3.75),
                StellarRecord(0.5, Int32(2), Int32(102), Int32(1), 0.8, 1.2, 1.5, 0.2, 4.10),
                StellarRecord(0.5, Int32(3), Int32(103), Int32(2), 1.5, 0.7, 2.0, 0.8, 3.60),
                StellarRecord(0.5, Int32(4), Int32(104), Int32(13), 2.0, 10.0, 5.0, 1.5, 4.50),
            ],
        )
        plot_hr(sev, vis; filename = "test_hr")
        @test isfile(joinpath(TESTDIR, "test_plots", "test_hr.png"))

        # --- HR evolution ---
        sev2 = StellarEvolutionSnapshot(
            1.0,
            2,
            [
                StellarRecord(1.0, Int32(1), Int32(101), Int32(2), 1.5, 0.75, 0.5, -0.2, 3.60),
                StellarRecord(1.0, Int32(2), Int32(102), Int32(4), 0.9, 1.10, 2.0, 0.4, 3.90),
            ],
        )
        plot_hr_evolution([sev, sev2], vis; filename = "test_hr_evo")
        @test isfile(joinpath(TESTDIR, "test_plots", "test_hr_evo.png"))
    end

    # =====================================================================
    @testset "Escaper and SSE plots (smoke tests)" begin
        Nbody6Dynamics.set_publication_theme!()

        plots_dir = joinpath(TESTDIR, "test_plots_f23")
        vis = VisualizationConfig(;
            output_dir = plots_dir,
            format = "png",
            dpi = 72,
            figsize = (6, 4),
        )

        # --- Escaper plots ---
        escs = [
            EscaperRecord(1.0, 0.5, -0.1, 15.0, 0, 101, 10.0, -20.0),
            EscaperRecord(2.0, 0.3, 0.2, 30.0, 1, 202, 120.0, 35.0),
            EscaperRecord(3.5, 1.4, 0.5, 250.0, 14, 303, 300.0, -60.0),
            EscaperRecord(4.0, 0.6, 0.1, 22.0, 11, 404, 200.0, 5.0),
        ]

        plot_escapers(escs, vis; filename = "test_escapers")
        @test isfile(joinpath(plots_dir, "test_escapers.png"))

        plot_escape_anisotropy(escs, vis; filename = "test_esc_aniso")
        @test isfile(joinpath(plots_dir, "test_esc_aniso.png"))

        # Empty-input guards: warn, return nothing, produce no file
        ret = @test_logs (:warn, r"No escaper records") plot_escapers(
            EscaperRecord[],
            vis;
            filename = "test_escapers_empty",
        )
        @test ret === nothing
        @test !isfile(joinpath(plots_dir, "test_escapers_empty.png"))

        ret = @test_logs (:warn, r"No escaper records") plot_escape_anisotropy(
            EscaperRecord[],
            vis;
            filename = "test_esc_aniso_empty",
        )
        @test ret === nothing
        @test !isfile(joinpath(plots_dir, "test_esc_aniso_empty.png"))

        # --- SSE plots ---
        # Mixed snapshot: two MS stars (finite TM, one past turnoff), one
        # giant with a partial core, one NS remnant (MC = M).  Records built
        # with the 9-argument convenience constructor carry NaN SSE fields
        # and must be skipped gracefully.
        sev_mix = StellarEvolutionSnapshot(
            50.0,
            5,
            [
                StellarRecord(
                    1.0,
                    Int32(1),
                    Int32(1),
                    Int32(0),
                    0.5,
                    0.4,
                    -0.6,
                    -0.3,
                    3.65,
                    8.0e4,
                    0.0,
                    0.0,
                    1e-10,
                ),
                StellarRecord(
                    1.0,
                    Int32(2),
                    Int32(2),
                    Int32(1),
                    1.2,
                    1.0,
                    0.0,
                    0.0,
                    3.76,
                    30.0,
                    0.0,
                    0.0,
                    1e-10,
                ),
                StellarRecord(
                    1.0,
                    Int32(3),
                    Int32(3),
                    Int32(3),
                    2.0,
                    1.8,
                    1.8,
                    1.2,
                    3.68,
                    40.0,
                    0.25,
                    0.02,
                    15.0,
                ),
                StellarRecord(
                    1.0,
                    Int32(4),
                    Int32(4),
                    Int32(13),
                    3.0,
                    1.4,
                    -5.0,
                    -5.0,
                    5.0,
                    10.0,
                    1.4,
                    1e-5,
                    1e-10,
                ),
                StellarRecord(1.0, Int32(5), Int32(5), Int32(1), 0.9, 0.8, -0.1, -0.1, 3.70),
            ],
        )

        plot_mass_segregation(sev_mix, vis; filename = "test_mass_seg")
        @test isfile(joinpath(plots_dir, "test_mass_seg.png"))

        plot_evolutionary_clock(sev_mix, vis; filename = "test_evo_clock")
        @test isfile(joinpath(plots_dir, "test_evo_clock.png"))

        plot_core_mass([sev_mix], vis; filename = "test_core_mass")
        @test isfile(joinpath(plots_dir, "test_core_mass.png"))

        # MS-only snapshot: no evolved stars — core-mass plot must warn,
        # return nothing, and create no file.
        sev_ms = StellarEvolutionSnapshot(
            0.5,
            2,
            [
                StellarRecord(
                    0.5,
                    Int32(1),
                    Int32(1),
                    Int32(0),
                    1.0,
                    0.8,
                    0.1,
                    -0.5,
                    3.75,
                    90.0,
                    0.0,
                    0.0,
                    1e-10,
                ),
                StellarRecord(
                    0.5,
                    Int32(2),
                    Int32(2),
                    Int32(1),
                    0.8,
                    1.2,
                    1.5,
                    0.2,
                    4.10,
                    5.2,
                    0.0,
                    0.0,
                    1e-10,
                ),
            ],
        )
        ret = @test_logs (:warn, r"No evolved stars") plot_core_mass(
            [sev_ms],
            vis;
            filename = "test_core_mass_ms",
        )
        @test ret === nothing
        @test !isfile(joinpath(plots_dir, "test_core_mass_ms.png"))

        # NaN-TM main-sequence snapshot (pre-v2026.07 data): the clock has
        # no valid TM and must warn + return nothing without a file.
        sev_nan = StellarEvolutionSnapshot(
            0.5,
            1,
            [StellarRecord(0.5, Int32(1), Int32(1), Int32(1), 1.0, 0.8, 0.1, -0.5, 3.75)],
        )
        ret = @test_logs (:warn, r"No main-sequence records") plot_evolutionary_clock(
            sev_nan,
            vis;
            filename = "test_evo_clock_nan",
        )
        @test ret === nothing
        @test !isfile(joinpath(plots_dir, "test_evo_clock_nan.png"))
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
    @testset "Diagnostics — key-value ADJUST format" begin
        diag_path = joinpath(TESTDIR, "out1000_kv")
        write(
            diag_path,
            """
 PHYSICAL SCALING:  R* = 2.500  M* = 1000.0  V* = 4.200  T* = 12.00
                    <M> = 0.500  SU = 1.0  AU = 1.0

 ADJUST:  TIME   0.0000  T[Myr]      0.00  Q  1.000  DE  0.00E+00  ETOT  -0.2500
 RMIN =    0.001 RSCALE =    1.234
 TIME[NB]    0.0000 N    10000 <NB>      0 NPAIRS    500

 ADJUST:  TIME   0.5000  T[Myr]     50.00  Q  0.987  DE  1.23E-06  ETOT  -0.2499
 RMIN =    0.002 RSCALE =    1.300
 TIME[NB]    0.5000 N     9990 <NB>      0 NPAIRS    495

 END RUN
""",
        )

        diag = read_diagnostics(diag_path)

        @test length(diag.adjust) == 2
        # First epoch
        @test diag.adjust[1].time_nb ≈ 0.0
        @test diag.adjust[1].time_myr ≈ 0.0
        @test diag.adjust[1].qvir ≈ 1.0
        @test diag.adjust[1].n == 10000
        @test diag.adjust[1].npairs == 500
        @test diag.adjust[1].rscale ≈ 1.234

        # Second epoch — TIME[NB] values override ADJUST defaults
        @test diag.adjust[2].time_nb ≈ 0.5
        @test diag.adjust[2].n == 9990
        @test diag.adjust[2].npairs == 495
        @test diag.adjust[2].rscale ≈ 1.300

        # Physical scaling
        @test diag.physical_scaling["R*"] ≈ 2.5
        @test diag.physical_scaling["T*"] ≈ 12.0
    end

    # =====================================================================
    @testset "Diagnostics — mixed format with partial epochs" begin
        diag_path = joinpath(TESTDIR, "out1000_mixed")
        write(
            diag_path,
            """
 ADJUST:    0.0000      0.00  1.000  0.00E+00 -0.2500  10000    500   1.234
 ADJUST:  TIME   1.0000  T[Myr]    100.00  Q  0.990  DE  5.00E-06  ETOT  -0.2480
 TIME[NB]    1.0000 N     9500 <NB>      0 NPAIRS    450
""",
        )
        diag = read_diagnostics(diag_path)
        @test length(diag.adjust) == 2

        # First: positional format
        @test diag.adjust[1].time_nb ≈ 0.0
        @test diag.adjust[1].n == 10000
        @test diag.adjust[1].rscale ≈ 1.234

        # Second: key-value + TIME[NB] merge
        @test diag.adjust[2].time_nb ≈ 1.0
        @test diag.adjust[2].n == 9500
        @test diag.adjust[2].npairs == 450
    end

    # =====================================================================
    @testset "Auto FPS calculation" begin
        _auto_fps = Nbody6Dynamics._auto_fps

        # Few frames → clamped to min
        @test _auto_fps(5; target_duration = 12.0, min_fps = 1, max_fps = 10) == 1
        # Many frames → clamped to max
        @test _auto_fps(500; target_duration = 12.0, min_fps = 2, max_fps = 30) == 30
        # Medium range → computed value
        fps = _auto_fps(120; target_duration = 12.0, min_fps = 2, max_fps = 30)
        @test fps == 10
        # Exact target
        @test _auto_fps(24; target_duration = 12.0, min_fps = 1, max_fps = 30) == 2
    end

    # =====================================================================
    @testset "Animations (smoke tests)" begin
        Nbody6Dynamics.set_publication_theme!()

        vis = VisualizationConfig(;
            output_dir = joinpath(TESTDIR, "test_anims"),
            format = "png",
            dpi = 72,
            figsize = (4, 3),
        )

        # Build two simple test snapshots
        n = 30
        params1 = zeros(Float32, 20)
        params1[1] = 0.0f0
        params2 = zeros(Float32, 20)
        params2[1] = 1.0f0
        hdr1 = SnapshotHeader(Int32(n), Int32(1), Int32(1), Int32(20), params1)
        hdr2 = SnapshotHeader(Int32(n), Int32(2), Int32(1), Int32(20), params2)
        snap1 = Snapshot(
            hdr1,
            Int32.(1:n),
            Float32.(rand(n)),
            Float32.(randn(3, n)),
            Float32.(0.1 .* randn(3, n)),
            Float32[],
            Float32[],
        )
        snap2 = Snapshot(
            hdr2,
            Int32.(1:n),
            Float32.(rand(n)),
            Float32.(randn(3, n) .+ 0.5),
            Float32.(0.1 .* randn(3, n)),
            Float32[],
            Float32[],
        )

        # --- Cluster animation (explicit fps) ---
        outpaths = animate_cluster([snap1, snap2], vis; filename = "test_cluster_anim", fps = 2)
        @test all(isfile, outpaths)
        @test all(p -> endswith(p, ".gif"), outpaths)

        # --- Cluster animation (auto fps) ---
        outpaths = animate_cluster([snap1, snap2], vis; filename = "test_cluster_anim_auto")
        @test all(isfile, outpaths)
        @test all(p -> endswith(p, ".gif"), outpaths)

        # --- Lagrangian animation ---
        lagr = LagrangianData(
            [0.0, 0.5, 1.0],
            [0.1, 0.5, 1.0],
            [0.05 0.06 0.07; 0.2 0.22 0.25; 1.0 1.05 1.1],
        )
        outpath = animate_lagrangian(
            lagr,
            vis;
            filename = "test_lagr_anim",
            fps = 2,
            selected_fractions = [0.1, 0.5, 1.0],
        )
        @test isfile(outpath)
        @test endswith(outpath, ".gif")

        # --- Lagrangian animation (auto fps) ---
        outpath = animate_lagrangian(
            lagr,
            vis;
            filename = "test_lagr_anim_auto",
            selected_fractions = [0.1, 0.5, 1.0],
        )
        @test isfile(outpath)
        @test endswith(outpath, ".gif")

        # --- HR animation ---
        sev1 = StellarEvolutionSnapshot(
            0.5,
            2,
            [
                StellarRecord(0.5, Int32(1), Int32(101), Int32(0), 1.0, 0.8, 0.1, -0.5, 3.75),
                StellarRecord(0.5, Int32(2), Int32(102), Int32(1), 0.8, 1.2, 1.5, 0.2, 4.10),
            ],
        )
        sev2 = StellarEvolutionSnapshot(
            1.0,
            2,
            [
                StellarRecord(1.0, Int32(1), Int32(101), Int32(2), 1.5, 0.75, 0.5, -0.2, 3.60),
                StellarRecord(1.0, Int32(2), Int32(102), Int32(4), 0.9, 1.10, 2.0, 0.4, 3.90),
            ],
        )
        outpath = animate_hr([sev1, sev2], vis; filename = "test_hr_anim", fps = 2)
        @test isfile(outpath)
        @test endswith(outpath, ".gif")

        # --- HR animation (auto fps) ---
        outpath = animate_hr([sev1, sev2], vis; filename = "test_hr_anim_auto")
        @test isfile(outpath)
        @test endswith(outpath, ".gif")
    end

    # =====================================================================
    @testset "External post-processing" begin
        # Create a fake output directory with some recognisable files
        ext_dir = joinpath(TESTDIR, "external_output")
        mkpath(ext_dir)

        # --- scan_output on empty directory ---
        scan_empty = scan_output(ext_dir)
        @test scan_empty.dir == abspath(ext_dir)
        @test scan_empty.available[:snapshots_conf3] == false
        @test scan_empty.available[:diagnostics] == false
        @test scan_empty.available[:lagr] == false
        @test scan_empty.available[:escapers] == false
        @test scan_empty.available[:stellar_evo] == false
        @test scan_empty.available[:binary_evo] == false
        @test scan_empty.available[:snapshots_hdf5] == false

        # --- Populate with dummy files ---
        # Write a minimal valid out1000
        open(joinpath(ext_dir, "out1000"), "w") do io
            println(io, " ADJUST: TIME  1.0  T[MYR]  0.5  Q  0.5  DE  1e-8  E  -0.25")
        end
        # Write a dummy lagr.7
        open(joinpath(ext_dir, "lagr.7"), "w") do io
            println(io, "# TIME  0.01  0.10  0.50  1.00")
            println(io, "0.0  0.1  0.3  1.0  5.0")
            println(io, "1.0  0.1  0.3  1.1  5.2")
        end
        # Write a dummy esc.11
        open(joinpath(ext_dir, "esc.11"), "w") do io
            println(io, "# escaper data")
            println(io, "1  1.0  0.5  1.0 2.0 3.0  0.1 0.2 0.3")
        end

        # The engine's regularised-binary record (real fixture)
        cp(joinpath(@__DIR__, "fixtures", "bev.82_0"), joinpath(ext_dir, "bev.82_0"); force = true)

        # --- scan_output with partial data ---
        scan = scan_output(ext_dir)
        @test scan.available[:diagnostics] == true
        @test scan.available[:lagr] == true
        @test scan.available[:escapers] == true
        @test scan.available[:snapshots_conf3] == false
        @test scan.available[:stellar_evo] == false
        @test scan.available[:binary_evo] == true && length(scan.binary_evo_files) == 1

        # --- OutputScan display ---
        buf = IOBuffer()
        show(buf, MIME("text/plain"), scan)
        output_str = String(take!(buf))
        @test occursin("Diagnostics", output_str)
        @test occursin("out1000", output_str)
        @test occursin("not found", output_str)  # for missing categories
        @test occursin("Regularised binaries", output_str)
        @test occursin("binary_population", output_str)

        # --- scan_output error on non-existent directory ---
        @test_throws ErrorException scan_output("/nonexistent/path")

        # --- postprocess_external (data only, no plots) ---
        results = postprocess_external(ext_dir; make_plots = false)
        @test haskey(results, :scan)
        @test results[:scan] isa OutputScan
        @test haskey(results, :diagnostics)
        @test haskey(results, :lagr)
        @test haskey(results, :binary_evo)
        bevs_ext = results[:binary_evo]
        @test bevs_ext isa Vector{BinaryEvolutionSnapshot} && length(bevs_ext) == 1
        @test bevs_ext[1].n_pairs == length(bevs_ext[1].records) > 0
        # A bev.82 record without diagnostics is still read, and an absent one is absent
        rm(joinpath(ext_dir, "bev.82_0"))
        @test !haskey(postprocess_external(ext_dir; make_plots = false), :binary_evo)
    end

    # =====================================================================
    # Adversarial external post-processing tests (included from separate file)
    # =====================================================================
    include("test_external_adversarial_inner.jl")

    # =====================================================================
    @testset "Initial Conditions — Merger IC Generator" begin
        # StableRNGs guarantees an identical stream across Julia versions, so
        # reference values in these tests survive upgrades.
        using StableRNGs
        rng = StableRNG(12345)

        # --- Plummer sampler ---
        @testset "Plummer sampler" begin
            N = 1000
            a = 1.0
            pos, vel = sample_plummer(N, a; rng = rng)
            @test size(pos) == (3, N)
            @test size(vel) == (3, N)
            # Particles should be centred near origin
            cm = vec(sum(pos, dims = 2)) ./ N
            @test all(abs.(cm) .< 0.5)
            # Radii should be finite and positive
            radii = [sqrt(sum(pos[:, i] .^ 2)) for i in 1:N]
            @test all(radii .> 0)
            @test all(isfinite.(radii))
            # Half-mass radius should be roughly 1.305 × a for Plummer
            r_sorted = sort(radii)
            r_hm = r_sorted[N ÷ 2]
            @test 0.5 * a < r_hm < 3.0 * a
        end

        # --- King sampler ---
        @testset "King sampler" begin
            N = 1000
            W0 = 6.0
            rt = 10.0
            pos, vel = sample_king(N, W0, rt; rng = rng)
            @test size(pos) == (3, N)
            @test size(vel) == (3, N)
            # All particles should be within the tidal radius (with some tolerance
            # from the sampling/scaling)
            radii = [sqrt(sum(pos[:, i] .^ 2)) for i in 1:N]
            @test all(isfinite.(radii))
            # Invalid W0
            @test_throws ArgumentError sample_king(10, -1.0, 5.0)
            @test_throws ArgumentError sample_king(10, 6.0, -1.0)
        end

        # --- King ODE solver ---
        @testset "King ODE solver" begin
            rhat, What, rho = Nbody6Dynamics._solve_king(6.0)
            @test length(rhat) == length(What) == length(rho)
            @test What[1] ≈ 6.0
            @test What[end] ≈ 0.0 atol = 0.01
            @test rhat[1] ≈ 0.0
            @test rhat[end] > 0.0
            # Density should be monotonically decreasing
            @test all(diff(rho[1:(end - 1)]) .≤ 0.01)
        end

        # --- King concentration vs published values (physics validation) ---
        @testset "King concentration c(W0)" begin
            # c = log10(r_t/r_0) for the standard dimensionless King equation;
            # reference values from the King (1966) model tables.
            for (W0, c_ref) in [
                (3.0, 0.672),
                (5.0, 1.029),
                (6.0, 1.255),
                (7.0, 1.528),
                (9.0, 2.119),
                (12.0, 2.739),
            ]
                rhat, _, _ = Nbody6Dynamics._solve_king(W0)
                c = log10(rhat[end])
                @test isapprox(c, c_ref; rtol = 0.01)
            end
        end

        # --- Kroupa IMF ---
        @testset "Kroupa IMF" begin
            masses = sample_kroupa(5000; m_low = 0.08, m_up = 100.0, rng = rng)
            @test length(masses) == 5000
            @test all(masses .≥ 0.08)
            @test all(masses .≤ 100.0)
            # Most stars should be low-mass (IMF is bottom-heavy)
            @test count(m -> m < 1.0, masses) > 3000
            # Custom mass limits
            masses_narrow = sample_kroupa(500; m_low = 1.0, m_up = 10.0, rng = rng)
            @test all(1.0 .≤ masses_narrow .≤ 10.0)
        end

        # --- Virialisation ---
        @testset "Virialise" begin
            N = 200
            mass = ones(N) ./ N
            pos = randn(rng, 3, N)
            vel = randn(rng, 3, N) .* 0.1
            virialise!(mass, pos, vel)
            # CM should be at origin
            cm_pos = vec(sum(mass' .* pos, dims = 2))
            cm_vel = vec(sum(mass' .* vel, dims = 2))
            @test all(abs.(cm_pos) .< 1e-10)
            @test all(abs.(cm_vel) .< 1e-10)
            # Virial ratio should be ~0.5
            T = 0.5 * sum(mass[i] * sum(vel[:, i] .^ 2) for i in 1:N)
            W = 0.0
            for i in 1:N, j in (i + 1):N
                dr = pos[:, i] .- pos[:, j]
                W -= mass[i] * mass[j] / sqrt(sum(dr .^ 2))
            end
            Q = T / abs(W)
            @test Q ≈ 0.5 atol = 0.01
        end

        # --- Kepler velocity ---
        @testset "Kepler velocity" begin
            # Equal mass, circular orbit
            v1, v2 = kepler_velocity(1.0, 1.0, 10.0, 0.0)
            @test v1 ≈ v2  # symmetric
            @test v1 > 0
            # Eccentric orbit has lower apocentre velocity
            v1e, v2e = kepler_velocity(1.0, 1.0, 10.0, 0.9)
            @test v1e < v1
            # Invalid inputs
            @test_throws ArgumentError kepler_velocity(1.0, 1.0, 10.0, 1.0)
            @test_throws ArgumentError kepler_velocity(1.0, 1.0, -1.0, 0.5)
        end

        # --- Jacobi radius ---
        @testset "Jacobi radius" begin
            rj = jacobi_radius(10.0, 1.0, 3.0)
            @test rj ≈ 10.0 * (1.0 / 9.0)^(1/3)
            # Equal mass: rJ = d × (1/3)^{1/3}
            rj_eq = jacobi_radius(10.0, 1.0, 1.0)
            @test rj_eq ≈ 10.0 * (1/3)^(1/3) atol = 1e-10
        end

        # --- Two-cluster orbit setup ---
        @testset "Two-cluster orbit" begin
            N1, N2 = 100, 100
            pos1 = randn(rng, 3, N1) .* 0.5
            vel1 = randn(rng, 3, N1) .* 0.01
            mass1 = ones(N1) ./ N1
            pos2 = randn(rng, 3, N2) .* 0.5
            vel2 = randn(rng, 3, N2) .* 0.01
            mass2 = ones(N2) ./ N2
            # Centre each cluster exactly (as virialise! does in the real
            # pipeline) so the COM/momentum/energy identities below are exact.
            for (p, v, m) in ((pos1, vel1, mass1), (pos2, vel2, mass2))
                p .-= sum(m' .* p, dims = 2) ./ sum(m)
                v .-= sum(m' .* v, dims = 2) ./ sum(m)
            end

            pos, vel, mass, ranges, _ = Nbody6Dynamics.setup_two_cluster_orbit(
                pos1,
                vel1,
                mass1,
                pos2,
                vel2,
                mass2,
                20.0,
                0.5;
                truncate_jacobi_flag = false,
            )
            @test length(mass) == N1 + N2
            @test size(pos, 2) == N1 + N2
            @test ranges == [1:N1, (N1 + 1):(N1 + N2)]
            # Combined CM should be near origin
            M = sum(mass)
            cm = vec(sum(mass' .* pos, dims = 2)) ./ M
            @test all(abs.(cm) .< 0.1)

            # Total momentum must vanish and the two-COM orbit must have the
            # Keplerian energy of the requested (d_apo, e) orbit.
            p = vec(sum(mass' .* vel, dims = 2))
            @test all(abs.(p) .< 1e-10)
            r1 = 1:N1
            r2 = (N1 + 1):(N1 + N2)
            M1 = sum(mass[r1])
            M2 = sum(mass[r2])
            com1 = vec(sum(mass[r1]' .* pos[:, r1], dims = 2)) ./ M1
            com2 = vec(sum(mass[r2]' .* pos[:, r2], dims = 2)) ./ M2
            vcom1 = vec(sum(mass[r1]' .* vel[:, r1], dims = 2)) ./ M1
            vcom2 = vec(sum(mass[r2]' .* vel[:, r2], dims = 2)) ./ M2
            d = sqrt(sum((com1 .- com2) .^ 2))
            @test d ≈ 20.0 rtol = 1e-10
            v_rel2 = sum((vcom1 .- vcom2) .^ 2)
            a_orb = 20.0 / (1.0 + 0.5)
            ε = 0.5 * v_rel2 - (M1 + M2) / d          # specific orbital energy, G = 1
            @test ε ≈ -(M1 + M2) / (2a_orb) rtol = 1e-10
        end

        # --- Plummer half-mass relation (physics validation) ---
        @testset "Plummer r_hm = 1.305 a" begin
            N = 20_000
            a = 1.0
            pos, _ = sample_plummer(N, a; rng = rng)
            mass = fill(1.0 / N, N)
            r_hm = Nbody6Dynamics.half_mass_radius(mass, pos; centre = zeros(3))
            @test r_hm ≈ 1.3048 rtol = 0.05    # statistical tolerance
        end

        # --- Kroupa mean mass vs analytic (physics validation) ---
        @testset "Kroupa mean mass" begin
            @test kroupa_mean_mass(0.08, 100.0) ≈ 0.58 rtol = 0.05
            N = 50_000
            m = sample_kroupa(N; m_low = 0.08, m_up = 100.0, rng = rng)
            μ = kroupa_mean_mass(0.08, 100.0)
            # Sample mean within a generous statistical band of the analytic mean
            @test abs(sum(m) / N - μ) / μ < 0.1
        end

        # --- verif_triorbit.toml Lagrange equilibrium (config regression) ---
        @testset "verif_triorbit config equilibrium" begin
            path = joinpath(@__DIR__, "..", "input_files", "verif_triorbit.toml")
            cfg = load_merger_config(path)
            @test length(cfg.clusters) == 3
            m = expected_mass(cfg.clusters[1].imf, cfg.clusters[1].N)   # 900 M☉, equal bodies
            @test m ≈ 900.0
            rc = 6.0
            ω = sqrt(Nbody6Dynamics._G_PC_KMS2_MSUN * m / (sqrt(3) * rc^3))   # km s⁻¹ pc⁻¹
            for spec in cfg.clusters
                r⃗ = spec.position
                v⃗ = spec.velocity
                @test sqrt(sum(r⃗ .^ 2)) ≈ rc rtol = 1e-3
                # Speed = ω r and velocity ⟂ radius (circular Lagrange orbit)
                @test sqrt(sum(v⃗ .^ 2)) ≈ ω * rc rtol = 1e-3
                @test abs(sum(r⃗ .* v⃗)) / (rc * ω * rc) < 1e-3
            end
            # Zero net momentum for equal-mass clusters
            vsum = sum(spec.velocity for spec in cfg.clusters)
            @test all(abs.(vsum) .< 0.05 * ω * rc)
        end

        # --- dat.10 writer ---
        @testset "dat.10 writer" begin
            N = 50
            mass = ones(N) ./ N
            pos = randn(rng, 3, N)
            vel = randn(rng, 3, N)
            dat10_path = joinpath(TESTDIR, "test_dat10.dat")
            write_dat10(dat10_path, mass, pos, vel)
            lines = readlines(dat10_path)
            @test length(lines) == N
            # Each line should have 7 columns
            cols = split(lines[1])
            @test length(cols) == 7
            # First column should be mass
            @test parse(Float64, cols[1]) ≈ 1.0 / N
        end

        # --- N-body unit conversion ---
        @testset "N-body unit conversion" begin
            mass = [0.5, 0.5]
            pos = [1.0 -1.0; 0.0 0.0; 0.0 0.0]
            vel = [0.0 0.0; 1.0 -1.0; 0.0 0.0]
            M_tot = 1e5  # M☉
            rbar = 2.0   # pc
            to_nbody_units!(mass, pos, vel, M_tot, rbar)
            @test sum(mass) ≈ 1e-5  # each mass = 0.5/1e5
            @test pos[1, 1] ≈ 0.5   # 1.0 / 2.0
        end

        # --- TOML config loading (kepler mode) ---
        @testset "Merger config TOML — kepler" begin
            cfg_path = joinpath(TESTDIR, "test_merger.toml")
            write(
                cfg_path,
                """
[merger]
n_clusters = 2
orbit_mode = "kepler"

[merger.cluster1]
model = "plummer"
N = 100
mass_total = 60.0
rbar = 1.0

[merger.cluster2]
model = "king"
N = 100
W0 = 5.0
mass_total = 60.0
rbar = 1.0

[merger.orbit]
apocentre = 10.0
eccentricity = 0.5

[merger.output]
format = "nbody"
truncate_jacobi = false
""",
            )
            cfg = load_merger_config(cfg_path)
            @test length(cfg.clusters) == 2
            @test cfg.orbit_mode == "kepler"
            @test cfg.clusters[1].profile isa PlummerProfile
            @test cfg.clusters[2].profile isa KingProfile
            @test cfg.clusters[2].profile.W0 == 5.0
            @test cfg.orbit.apocentre == 10.0
            @test cfg.orbit.eccentricity == 0.5
            @test cfg.output.format == "nbody"
            @test cfg.output.truncate_jacobi == false
        end

        # --- TOML config loading (explicit mode) ---
        @testset "Merger config TOML — explicit" begin
            cfg_path = joinpath(TESTDIR, "test_merger_explicit.toml")
            write(
                cfg_path,
                """
[merger]
n_clusters = 3
orbit_mode = "explicit"

[merger.cluster1]
model = "plummer"
N = 50
mass_total = 30.0
rbar = 1.0
position = [-5.0, 0.0, 0.0]
velocity = [1.0, 0.0, 0.0]

[merger.cluster2]
model = "king"
N = 50
W0 = 4.0
mass_total = 30.0
rbar = 1.0
position = [2.5, 4.33, 0.0]
velocity = [-0.5, -0.87, 0.0]

[merger.cluster3]
model = "plummer"
N = 50
mass_total = 30.0
rbar = 0.8
position = [2.5, -4.33, 0.0]
velocity = [-0.5, 0.87, 0.0]

[merger.output]
format = "nbody"
truncate_jacobi = false
""",
            )
            cfg = load_merger_config(cfg_path)
            @test length(cfg.clusters) == 3
            @test cfg.orbit_mode == "explicit"
            @test cfg.clusters[1].position == [-5.0, 0.0, 0.0]
            @test cfg.clusters[3].velocity == [-0.5, 0.87, 0.0]
        end

        # --- Fail-fast validation of merger TOML values ---
        @testset "Merger config validation" begin
            # Template with one overridable block per constrained section
            function _merger_toml(;
                cluster1::String = "model = \"king\"\nN = 100\nW0 = 5.0\nrbar = 1.0",
                orbit::String = "apocentre = 10.0\neccentricity = 0.5",
                output::String = "tcrit = 10.0\ndtadj = 1.0\ndeltat = 1.0",
            )
                """
                [merger]
                n_clusters = 2
                orbit_mode = "kepler"

                [merger.cluster1]
                $cluster1

                [merger.cluster2]
                model = "plummer"
                N = 100
                rbar = 1.0

                [merger.orbit]
                $orbit

                [merger.output]
                $output
                """
            end

            # Valid template loads cleanly — smoke test
            ok_path = joinpath(TESTDIR, "merger_val_ok.toml")
            write(ok_path, _merger_toml())
            @test load_merger_config(ok_path) isa MergerConfig

            cases = [
                (
                    "N",
                    _merger_toml(cluster1 = "model = \"king\"\nN = 1\nrbar = 1.0"),
                    "merger.cluster1.N",
                ),
                (
                    "rbar",
                    _merger_toml(cluster1 = "model = \"king\"\nN = 100\nrbar = -1.0"),
                    "merger.cluster1.rbar",
                ),
                (
                    "W0",
                    _merger_toml(cluster1 = "model = \"king\"\nN = 100\nW0 = -2.0\nrbar = 1.0"),
                    "merger.cluster1.W0",
                ),
                (
                    "kroupa_bounds",
                    _merger_toml(
                        cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n" *
                                   "imf = \"kroupa\"\nbodyn = 50.0\nbody1 = 0.1",
                    ),
                    "0 < bodyn < body1",
                ),
                (
                    "target_mass",
                    _merger_toml(
                        cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n" *
                                   "imf = \"kroupa\"\nmass_total = -500.0",
                    ),
                    "target_mass",
                ),
                (
                    "particle_mass",
                    _merger_toml(
                        cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n" *
                                   "imf = \"equal\"\nmass_total = 0.0",
                    ),
                    "particle_mass",
                ),
                (
                    "apocentre",
                    _merger_toml(orbit = "apocentre = -5.0\neccentricity = 0.5"),
                    "merger.orbit.apocentre",
                ),
                (
                    "eccentricity",
                    _merger_toml(orbit = "apocentre = 10.0\neccentricity = 1.0"),
                    "merger.orbit.eccentricity",
                ),
                (
                    "tcrit",
                    _merger_toml(output = "tcrit = 0.0\ndtadj = 1.0\ndeltat = 1.0"),
                    "merger.output.tcrit",
                ),
                (
                    "dtadj",
                    _merger_toml(output = "tcrit = 10.0\ndtadj = -0.5\ndeltat = 1.0"),
                    "merger.output.dtadj",
                ),
                (
                    "deltat",
                    _merger_toml(output = "tcrit = 10.0\ndtadj = 1.0\ndeltat = 0.0"),
                    "merger.output.deltat",
                ),
                ("nbody6_qe", _merger_toml() * "\n[merger.nbody6]\nqe = 0.0\n", "merger.nbody6.qe"),
                (
                    "nbody6_kz16",
                    _merger_toml() * "\n[merger.nbody6]\nkz16 = 5\n",
                    "merger.nbody6.kz16",
                ),
                (
                    "nbody6_rs0",
                    _merger_toml() * "\n[merger.nbody6]\nrs0 = -0.1\n",
                    "merger.nbody6.rs0",
                ),
                (
                    "nbody6_tcomp",
                    _merger_toml() * "\n[merger.nbody6]\ntcomp = 0.0\n",
                    "merger.nbody6.tcomp",
                ),
                (
                    "nbody6_gmin_gmax",
                    _merger_toml() * "\n[merger.nbody6]\ngmin = 0.1\ngmax = 0.01\n",
                    "gmin < gmax",
                ),
                (
                    "nbody6_ncrit",
                    _merger_toml() * "\n[merger.nbody6]\nncrit = 0\n",
                    "merger.nbody6.ncrit",
                ),
                (
                    "nbody6_kz_index",
                    _merger_toml() * "\n[merger.nbody6.kz]\n\"99\" = 1\n",
                    "merger.nbody6.kz index",
                ),
                (
                    "stellar_zmet",
                    _merger_toml() * "\n[merger.stellar]\nzmet = 0.1\n",
                    "merger.stellar.zmet",
                ),
                (
                    "stellar_level",
                    _merger_toml() * "\n[merger.stellar]\nlevel = \"X\"\n",
                    "merger.stellar.level",
                ),
                (
                    "stellar_dtplot",
                    _merger_toml() * "\n[merger.stellar]\ndtplot = 0.5\n",
                    "merger.stellar.dtplot",
                ),
                (
                    "stellar_epoch0",
                    _merger_toml() * "\n[merger.stellar]\nepoch0 = 5.0\n",
                    "merger.stellar.epoch0",
                ),
                ("tidal_kz14_3", _merger_toml() * "\n[merger.tidal]\nkz14 = 3\n", "not supported"),
                (
                    "tidal_kz14_bad",
                    _merger_toml() * "\n[merger.tidal]\nkz14 = 7\n",
                    "merger.tidal.kz14",
                ),
                (
                    "tidal_qe",
                    _merger_toml() * "\n[merger.tidal]\nkz14 = 2\ngmg = 1.0e11\nrg0 = 8.5\n",
                    "tidal work",
                ),
                (
                    "tidal_gmg",
                    _merger_toml() * "\n[merger.tidal]\nkz14 = 2\nrg0 = 8.5\n",
                    "merger.tidal.gmg",
                ),
                (
                    "tidal_rg",
                    _merger_toml() * "\n[merger.tidal]\nkz14 = 5\nvg = [0.0, 220.0, 0.0]\n",
                    "merger.tidal.rg",
                ),
                (
                    "binaries_fraction",
                    _merger_toml(
                        cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n[merger.cluster1.binaries]\nfraction = 1.5",
                    ),
                    "binaries.fraction",
                ),
                (
                    "imf_rescale_factor",
                    _merger_toml(
                        cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n" *
                                   "imf = \"kroupa\"\nmass_total = 5000.0",
                    ),
                    "rescale factor",
                ),
            ]
            for (label, body, expected) in cases
                path = joinpath(TESTDIR, "merger_val_bad_$label.toml")
                write(path, body)
                @test_throws expected load_merger_config(path)
            end

            # Moderate rescale factors load with a warning (×1.6 here)
            warn_path = joinpath(TESTDIR, "merger_val_warn_rescale.toml")
            write(
                warn_path,
                _merger_toml(
                    cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\nimf = \"kroupa\"\n" *
                               "mass_total = $(round(1.6 * 100 * kroupa_mean_mass(0.08, 100.0); digits = 1))",
                ),
            )
            @test_logs (:warn, r"rescales the Kroupa IMF") match_mode = :any load_merger_config(
                warn_path,
            )
        end

        # --- [merger.nbody6] integration parameters and regime guards ---
        @testset "Nbody6 integration parameters" begin
            p0 = Nbody6ParameterSpec()
            @test p0.qe == 2.0e-4 && p0.nnbopt == 0 && p0.rs0 == 0.0 && p0.kz16 == 0
            ρ̂_plummer = Nbody6Dynamics._central_density_contrast(PlummerProfile())
            @test ρ̂_plummer ≈ 2 * 1.305^3 rtol = 0.01
            c6 = Nbody6Dynamics._central_density_contrast(KingProfile(W0 = 6.0))
            c3 = Nbody6Dynamics._central_density_contrast(KingProfile(W0 = 3.0))
            @test isfinite(c6) && c6 > c3 > 1        # monotonic in W0 (3.3 and 9.3)
            @test c3 < ρ̂_plummer < c6                # Plummer's extended halo sits between W0 = 3 and 6

            clusters = [
                ClusterSpec(profile = PlummerProfile(), N = 100, rbar = 1.0, imf = KroupaIMF()),
                ClusterSpec(
                    profile = KingProfile(W0 = 6.0),
                    N = 120,
                    rbar = 2.0,
                    imf = KroupaIMF(),
                ),
            ]
            ranges = [1:100, 101:220]
            r = resolve_nbody6_parameters(p0, clusters, ranges, 220, 4.0)
            @test r.nnbopt == 20                          # clamp(round(√220) = 15, 20, 300)
            @test r.rs0 ≈ min(2 * 0.25 * cbrt(2 * 20 / 100), 0.25)   # r_h = 1/4, N_min = 100, doubled rule
            @test r.rmin ≈ 4 * 0.25 / (100 * cbrt(ρ̂_plummer))
            @test r.dtmin ≈ 0.04 * sqrt(r.rmin^3 * 220)
            @test r.qe == 2.0e-4 && r.kz16 == 0
            user = Nbody6ParameterSpec(;
                qe = 1e-3,
                nnbopt = 50,
                rs0 = 0.1,
                rmin = 1e-4,
                dtmin = 1e-6,
                kz16 = 2,
            )
            ru = resolve_nbody6_parameters(user, clusters, ranges, 220, 4.0)
            @test ru.nnbopt == 50 && ru.rs0 == 0.1 && ru.rmin == 1e-4 && ru.dtmin == 1e-6
            @test ru.kz16 == 2 && ru.qe == 1e-3
            @test_throws ArgumentError Nbody6Dynamics._assert_resolved(p0)

            # Guards: oversized RS0 refused, cold collapse and unresolved members warned
            wide = Nbody6ParameterSpec(; nnbopt = 20, rs0 = 0.5, rmin = 1e-4, dtmin = 1e-6)
            @test_throws ErrorException Nbody6Dynamics._check_multicluster_regime(
                0.5,
                2.0,
                wide,
                0.25,
            )
            @test_logs (:warn, r"cold-collapse") Nbody6Dynamics._check_multicluster_regime(
                0.1,
                2.0,
                r,
                0.25,
            )
            @test_logs (:warn, r"unresolved") Nbody6Dynamics._check_multicluster_regime(
                0.5,
                8.0,
                r,
                0.25,
            )
            @test_logs Nbody6Dynamics._check_multicluster_regime(0.5, 2.0, r, 0.25)
            @test_logs (:warn, r"undersampled") Nbody6Dynamics._check_multicluster_regime(
                0.5,
                2.0,
                r,
                0.25;
                deltat = 1.0,
                t_cr_member_min_nb = 0.2,
            )
            @test_logs Nbody6Dynamics._check_multicluster_regime(
                0.5,
                2.0,
                r,
                0.25;
                deltat = 0.1,
                t_cr_member_min_nb = 0.2,
            )

            # Writer: unresolved specs refused, resolved values written
            inp = joinpath(mktempdir(), "t.inp")
            @test_throws ArgumentError generate_merger_inp(inp, 220, 4.0, 0.6; nbody6 = p0)
            generate_merger_inp(inp, 220, 4.0, 0.6; nbody6 = ru, tcrit = 5.0)
            txt = read(inp, String)
            @test occursin("NNBOPT=50,", txt) && occursin("QE=1.000E-03", txt)
            @test occursin("RS0=0.1,", txt) && occursin("DTMIN=1.000E-06,RMIN=1.000E-04,", txt)
            @test occursin("KZ(11:20)=0 1 0 0 0 2 0 0 3 0", txt)

            # [merger.nbody6] round trip
            nb_path = joinpath(TESTDIR, "merger_nbody6.toml")
            write(
                nb_path,
                """
[merger]
n_clusters = 2
orbit_mode = "kepler"

[merger.cluster1]
model = "king"
N = 100
rbar = 1.0

[merger.cluster2]
model = "plummer"
N = 100
rbar = 1.0

[merger.nbody6]
qe = 1.0e-3
nnbopt = 30
kz16 = 1
tcomp = 7200.0
tcrtp0 = 500.0
isernb = 20
ncrit = 5
smax = 0.5

[merger.nbody6.kz]
"40" = 2
"19" = 4

[merger.stellar]
kz19 = 0
level = "0"
zmet = 0.02
epoch0 = -1.0
dtplot = 2.0
""",
            )
            cfg_nb = load_merger_config(nb_path)
            @test cfg_nb.nbody6.qe == 1.0e-3 && cfg_nb.nbody6.nnbopt == 30
            @test cfg_nb.nbody6.kz16 == 1 && cfg_nb.nbody6.rs0 == 0.0
            @test cfg_nb.nbody6.tcomp == 7200.0 && cfg_nb.nbody6.tcrtp0 == 500.0
            @test cfg_nb.nbody6.isernb == 20 &&
                  cfg_nb.nbody6.ncrit == 5 &&
                  cfg_nb.nbody6.smax == 0.5
            @test cfg_nb.nbody6.kz == Dict(40 => 2, 19 => 4)
            @test cfg_nb.stellar.kz19 == 0 && cfg_nb.stellar.level == "0"
            @test cfg_nb.stellar.zmet == 0.02 &&
                  cfg_nb.stellar.epoch0 == -1.0 &&
                  cfg_nb.stellar.dtplot == 2.0

            # Writer: run control, stellar settings, KZ overrides (an override of a named
            # index warns), and mass bounds reach the file
            rnb = resolve_nbody6_parameters(cfg_nb.nbody6, clusters, ranges, 220, 4.0)
            inp2 = joinpath(mktempdir(), "t2.inp")
            @test_logs (:warn, r"overrides KZ\(19\)") match_mode = :any generate_merger_inp(
                inp2,
                220,
                4.0,
                0.6;
                nbody6 = rnb,
                stellar = cfg_nb.stellar,
                mass_bounds = (0.5, 20.0),
            )
            txt2 = read(inp2, String)
            @test occursin("KSTART=1,TCOMP=7200,TCRTP0=500,isernb=20,iserreg=40,iserks=0 /", txt2)
            @test occursin("N=220,NFIX=1,NCRIT=5,", txt2) &&
                  occursin(",NNBOPT=30,NRUN=1,NCOMM=10,", txt2)
            @test occursin("SMAX=0.5,", txt2) && occursin("Level='0' /", txt2)
            @test occursin("KZ(11:20)=0 0 0 0 0 1 0 0 4 0", txt2)   # KZ(12) off with kz19 = 0; override 19 → 4
            @test occursin("KZ(31:40)=0 0 0 0 0 0 0 0 0 2", txt2)
            @test occursin("BODY1=20,BODYN=0.5,NBIN0=0,NHI0=0,ZMET=0.02,EPOCH0=-1,DTPLOT=2 /", txt2)

            # Tidal field: point-mass and MWPotential2014 namelists, isolated by default
            @test cfg_nb.tidal.kz14 == 0
            td2 = TidalSpec(; kz14 = 2, gmg = 1.0e11, rg0 = 8.5)
            inp3 = joinpath(mktempdir(), "t3.inp")
            generate_merger_inp(inp3, 220, 4.0, 0.6; nbody6 = ru, tidal = td2)
            txt3 = read(inp3, String)
            @test occursin("KZ(11:20)=0 1 0 2 0 2 0 0 3 0", txt3)
            @test occursin("&INXTRNL0\nGMG=1E+11,RG0=8.5 /", txt3)
            td5 = TidalSpec(; kz14 = 5, rg = [8.0, 0.0, 0.0], vg = [0.0, 220.0, 0.0])
            generate_merger_inp(inp3, 220, 4.0, 0.6; nbody6 = ru, tidal = td5)
            txt5 = read(inp3, String)
            @test occursin("KZ(11:20)=0 1 0 5 0 2 0 0 3 0", txt5)
            @test occursin("&INXTRNL0\nRG=8,0,0,VG=0,220,0 /", txt5)
            generate_merger_inp(inp3, 220, 4.0, 0.6; nbody6 = ru, tidal = TidalSpec(; kz14 = 1))
            @test !occursin("&INXTRNL0", read(inp3, String))
            @test_throws ErrorException generate_merger_inp(
                inp3,
                220,
                4.0,
                0.6;
                nbody6 = ru,
                tidal = TidalSpec(; kz14 = 4),
            )
            @test_throws ErrorException Nbody6Dynamics._validate_tidal(
                TidalSpec(; kz14 = 2, gmg = 1e11),
            )
            @test_throws ErrorException Nbody6Dynamics._validate_tidal(
                TidalSpec(; kz14 = 5, rg = [8.0, 0, 0]),
            )
            @test_throws ErrorException Nbody6Dynamics._validate_tidal_tolerance(
                td2,
                Nbody6ParameterSpec(),
            )
            @test Nbody6Dynamics._validate_tidal_tolerance(
                td2,
                Nbody6ParameterSpec(; qe = 0.05),
            ) === nothing
            @test Nbody6Dynamics._validate_tidal_tolerance(TidalSpec(), Nbody6ParameterSpec()) ===
                  nothing

            # Crossing time: G = 1 virial system with E = -M²/(4 r_v) has t_cr = (2 r_v)^{3/2}/√M
            @test crossing_time(1.0, -0.25) ≈ 2.0^1.5
            @test isnan(crossing_time(1.0, 0.1))
            @test Nbody6Dynamics._nbody_time_myr(1.0, 1.0) ≈ 14.91 rtol = 0.01
        end

        # --- Full pipeline — kepler mode (small N) ---
        @testset "Full merger IC pipeline — kepler" begin
            out_dir = mktempdir()
            cfg = MergerConfig(
                [
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 100,
                        rbar = 1.0,
                        imf = RescaledKroupaIMF(
                            bodyn = 0.1,
                            body1 = 50.0,
                            target_mass = 100 * kroupa_mean_mass(0.1, 50.0),
                        ),
                    ),
                    ClusterSpec(
                        profile = KingProfile(W0 = 5.0),
                        N = 100,
                        rbar = 1.0,
                        imf = RescaledKroupaIMF(
                            bodyn = 0.1,
                            body1 = 50.0,
                            target_mass = 100 * kroupa_mean_mass(0.1, 50.0),
                        ),
                    ),
                ],
                "kepler",
                OrbitSpec(apocentre = 10.0, eccentricity = 0.5),
                MergerOutputSpec(format = "nbody", truncate_jacobi = true, output_dir = out_dir),
            )
            result = generate_merger_ic(cfg; rng = rng)
            @test result isa MergerICResult
            @test isdir(result.output_dir)
            @test isfile(joinpath(result.output_dir, "dat.10"))
            @test isfile(joinpath(result.output_dir, "merger.inp"))
            @test isfile(joinpath(result.output_dir, "merger_summary.txt"))
            @test result.N_total > 0
            @test result.M_total > 0
            @test result.orbit_mode == "kepler"
            @test length(result.cluster_ranges) == 2

            # Verify dat.10 mass normalisation
            lines = readlines(joinpath(result.output_dir, "dat.10"))
            masses = [parse(Float64, split(l)[1]) for l in lines]
            @test sum(masses) ≈ 1.0 atol = 1e-10

            # Verify .inp has KZ(22)=2 and the resolved integration parameters
            inp_text = read(joinpath(result.output_dir, "merger.inp"), String)
            @test occursin("KZ(21:30)=0 2 2 0 0 1 0 0 0 1", inp_text)
            @test occursin("QE=2.000E-04", inp_text) && !occursin("RS0=0.5,", inp_text)
            ic_meta = Nbody6Dynamics.TOML.parsefile(joinpath(result.output_dir, "merger_ic.toml"))
            @test 0 < ic_meta["nbody6"]["rs0"] ≤ minimum(c.rbar for c in cfg.clusters) / result.rbar
            @test ic_meta["nbody6"]["rmin"] > 0 && ic_meta["nbody6"]["dtmin"] > 0
            @test isfinite(ic_meta["meta"]["q_virial"]) && ic_meta["meta"]["q_virial"] > 0
            @test ic_meta["meta"]["rbar_over_rhm_min"] > 1
            @test ic_meta["meta"]["t_cr_config_nb"] > ic_meta["meta"]["t_cr_member_min_nb"] > 0
            @test ic_meta["meta"]["t_star_myr"] > 0
            @test ic_meta["stellar"]["level"] == "C" && ic_meta["nbody6"]["tcrtp0"] == 3600.0
            @test ic_meta["nbody6"]["kz"] == Dict{String,Any}()
            summary_text = read(joinpath(result.output_dir, "merger_summary.txt"), String)
            @test occursin("virial ratio Q = T/|W|", summary_text)
            @test occursin("crossing time: configuration", summary_text)
            @test occursin("Integration (merger.inp", summary_text)
            @test occursin("Stellar evolution: KZ(19)=3 Level=C", summary_text)
            @test occursin("ZMET=0.001,EPOCH0=0,DTPLOT=1 /", inp_text)
            @test occursin("BODY1=50,BODYN=0.1,", inp_text)
        end

        # --- Full pipeline — explicit mode (3 clusters) ---
        @testset "Full merger IC pipeline — explicit 3-cluster" begin
            out_dir = mktempdir()
            cfg = MergerConfig(
                [
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 80,
                        rbar = 1.0,
                        imf = RescaledKroupaIMF(
                            bodyn = 0.1,
                            body1 = 50.0,
                            target_mass = 80 * kroupa_mean_mass(0.1, 50.0),
                        ),
                        position = [-5.0, 0.0, 0.0],
                        velocity = [1.0, 0.0, 0.0],
                    ),
                    ClusterSpec(
                        profile = KingProfile(W0 = 5.0),
                        N = 80,
                        rbar = 1.0,
                        imf = RescaledKroupaIMF(
                            bodyn = 0.1,
                            body1 = 50.0,
                            target_mass = 80 * kroupa_mean_mass(0.1, 50.0),
                        ),
                        position = [2.5, 4.33, 0.0],
                        velocity = [-0.5, -0.87, 0.0],
                    ),
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 60,
                        rbar = 0.8,
                        imf = EqualMassIMF(particle_mass = 400.0 / 60),
                        position = [2.5, -4.33, 0.0],
                        velocity = [-0.5, 0.87, 0.0],
                    ),
                ],
                "explicit",
                OrbitSpec(),  # ignored in explicit mode
                MergerOutputSpec(format = "nbody", truncate_jacobi = false, output_dir = out_dir),
            )
            result = generate_merger_ic(cfg; rng = rng)
            @test result isa MergerICResult
            @test result.orbit_mode == "explicit"
            @test length(result.cluster_ranges) == 3
            @test result.N_total == 80 + 80 + 60
            # Total mass should sum to 1 in N-body units
            lines = readlines(joinpath(result.output_dir, "dat.10"))
            @test sum(parse(Float64, split(l)[1]) for l in lines) ≈ 1.0 atol = 1e-10
        end

        # --- Equal-mass IMF ---
        @testset "Equal mass IMF" begin
            out_dir = mktempdir()
            cfg = MergerConfig(
                [
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 50,
                        rbar = 1.0,
                        imf = EqualMassIMF(particle_mass = 500.0 / 50),
                    ),
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 50,
                        rbar = 1.0,
                        imf = EqualMassIMF(particle_mass = 500.0 / 50),
                    ),
                ],
                "kepler",
                OrbitSpec(apocentre = 8.0, eccentricity = 0.3),
                MergerOutputSpec(format = "nbody", truncate_jacobi = false, output_dir = out_dir),
            )
            generate_merger_ic(cfg; rng = rng)
            lines = readlines(joinpath(out_dir, "dat.10"))
            @test length(lines) == 100
        end

        # --- Primordial binaries ---
        @testset "Primordial binaries" begin
            rng_b = StableRNG(777)
            m = sample_kroupa(1000; rng = rng_b)
            b = sample_binaries(BinarySpec(; fraction = 0.3), m, rng_b)
            n_b = round(Int, 0.3 * 1000 / 1.3)
            @test length(b.primary) == n_b == length(b.a_pc) == length(b.e)
            @test isempty(intersect(b.primary, b.secondary)) &&
                  allunique(vcat(b.primary, b.secondary))
            @test all(b.m1 .≥ b.m2)
            @test all(0 .≤ b.e .< 1) && 0.55 < sum(b.e) / n_b < 0.78            # thermal mean 2/3
            @test all(5e-8 .< b.a_pc .< 0.2)                                      # 0.01 AU … 4×10⁴ AU
            bc = sample_binaries(
                BinarySpec(;
                    fraction = 0.5,
                    period = "loguniform",
                    a_min = 1.0,
                    a_max = 10.0,
                    eccentricity = "circular",
                ),
                m,
                rng_b,
            )
            au = Nbody6Dynamics._AU_PC
            @test all(bc.e .== 0) && all(au .≤ bc.a_pc .≤ 10au)
            @test length(bc.primary) == round(Int, 0.5 * 1000 / 1.5)
            bq = sample_binaries(
                BinarySpec(; fraction = 0.2, pairing = "uniform_q", q_min = 0.5),
                m,
                rng_b,
            )
            @test all(0.5 .≤ bq.m2 ./ bq.m1 .≤ 1.0)
            @test isempty(sample_binaries(BinarySpec(), m, rng_b).primary)
            @test_throws ErrorException Nbody6Dynamics._validate_binaries(
                BinarySpec(; fraction = 1.0),
                "c",
            )
            @test_throws ErrorException Nbody6Dynamics._validate_binaries(
                BinarySpec(; pairing = "x"),
                "c",
            )
            @test_throws ErrorException Nbody6Dynamics._validate_binaries(
                BinarySpec(; a_min = 5.0, a_max = 1.0),
                "c",
            )
            # Kroupa (1995) periods stay within the distribution's support
            lp = [Nbody6Dynamics._sample_log_period_kroupa1995(rng_b) for _ in 1:2000]
            @test all(1.0 .≤ lp .≤ 8.43) && 3.0 < sum(lp) / length(lp) < 6.0

            # Kepler relative orbit: energy −GM/(2a) and r within [a(1−e), a(1+e)]
            for e_test in (0.0, 0.5, 0.9)
                r, v = Nbody6Dynamics._kepler_relative_orbit(2.0, 1e-3, e_test, rng_b)
                rn = sqrt(sum(r .^ 2))
                @test isapprox(0.5 * sum(v .^ 2) - 2.0 / rn, -2.0 / (2 * 1e-3); rtol = 1e-8)
                @test 1e-3 * (1 - e_test) - 1e-12 ≤ rn ≤ 1e-3 * (1 + e_test) + 1e-12
            end
            R = Nbody6Dynamics._random_rotation(rng_b)
            @test isapprox(R * R', [1 0 0; 0 1 0; 0 0 1]; atol = 1e-12)

            # Expansion of a tiny system set: pairs first, then singles; centres of mass kept
            pos_s = [0.0 1.0 2.0 10.0 11.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]
            vel_s = zeros(3, 5)
            mass_s = [2.0, 1.0, 1.5, 1.0, 1.0]
            cb1 = (
                system_binary = [1, 0, 2],
                m1 = [1.2, 0.9],
                m2 = [0.8, 0.6],
                a_pc = [1e-3, 2e-3],
                e = [0.0, 0.5],
            )
            cb2 = (
                system_binary = [0, 0],
                m1 = Float64[],
                m2 = Float64[],
                a_pc = Float64[],
                e = Float64[],
            )
            ex = expand_binaries(
                pos_s,
                vel_s,
                mass_s,
                [1:3, 4:5],
                [cb1, cb2],
                [[1, 2, 3], [1, 2]];
                rng = rng_b,
            )
            @test ex.n_pairs == [2, 0] && length(ex.mass) == 7
            @test ex.mass[1:4] ≈ [1.2, 0.8, 0.9, 0.6] && ex.mass[5:7] ≈ [1.0, 1.0, 1.0]
            @test ex.cluster_blocks[1] == [1:4, 5:5] && ex.cluster_blocks[2] == [5:4, 6:7]
            for (k, sys) in ((1, 1), (2, 3))
                i1, i2 = 2k - 1, 2k
                M = ex.mass[i1] + ex.mass[i2]
                com = (ex.mass[i1] .* ex.pos[:, i1] .+ ex.mass[i2] .* ex.pos[:, i2]) ./ M
                @test isapprox(com, pos_s[:, sys]; atol = 1e-12)
                r_rel = ex.pos[:, i1] .- ex.pos[:, i2]
                v_rel = ex.vel[:, i1] .- ex.vel[:, i2]
                @test 0.5 * sum(v_rel .^ 2) - M / sqrt(sum(r_rel .^ 2)) < 0
            end
            @test sum(ex.mass) ≈ sum(mass_s)
            @test isnan(ex.hard_fraction[2]) && 0 ≤ ex.hard_fraction[1] ≤ 1
            @test Nbody6Dynamics._members_from_blocks([1:4, 5:5]) == [1, 2, 3, 4, 5]

            # Summary parsing with and without pair counts
            sp = joinpath(TESTDIR, "summary_bin.txt")
            write(
                sp,
                "  Cluster 1: king, imf=kroupa, N=10 (after trunc: 9, binaries: 2), M=1 M☉\n" *
                "  Cluster 2: plummer, imf=kroupa, N=6 (after trunc: 5, binaries: 1), M=1 M☉\n",
            )
            mem = parse_merger_summary(sp)
            @test mem[1] == vcat(1:4, 7:11) && mem[2] == vcat(5:6, 12:14)
            write(
                sp,
                "  Cluster 1: king, imf=kroupa, N=10 (after trunc: 9), M=1 M☉\n" *
                "  Cluster 2: plummer, imf=kroupa, N=6 (after trunc: 5), M=1 M☉\n",
            )
            mem0 = parse_merger_summary(sp)
            @test mem0[1] == collect(1:9) && mem0[2] == collect(10:14)

            # Full generation with a binary-rich cluster: ordering, input, metadata round trip
            cfg_bin = MergerConfig(
                [
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 200,
                        rbar = 1.0,
                        imf = KroupaIMF(),
                        binaries = BinarySpec(; fraction = 0.3),
                    ),
                    ClusterSpec(
                        profile = KingProfile(W0 = 5.0),
                        N = 200,
                        rbar = 1.0,
                        imf = KroupaIMF(),
                    ),
                ],
                "kepler",
                OrbitSpec(apocentre = 10.0, eccentricity = 0.5),
                MergerOutputSpec(; output_dir = mktempdir());
                seed = 5,
            )
            res_b = generate_merger_ic(cfg_bin)
            @test res_b.n_pairs[2] == 0 && res_b.n_pairs[1] ≥ 30
            nb = sum(res_b.n_pairs)
            inp_b = read(joinpath(res_b.output_dir, "merger.inp"), String)
            @test occursin("NBIN0=$(nb),", inp_b) &&
                  occursin("KZ(1:10)=1 -1 2 0 0 0 3 2 0 0", inp_b)
            dat =
                [parse.(Float64, split(l)) for l in eachline(joinpath(res_b.output_dir, "dat.10"))]
            @test length(dat) == res_b.N_total
            seps = [sqrt(sum((dat[2k - 1][2:4] .- dat[2k][2:4]) .^ 2)) for k in 1:nb]
            @test maximum(seps) < 0.05                          # pairs are far smaller than a cluster
            @test length(res_b.cluster_ranges[1]) + length(res_b.cluster_ranges[2]) == res_b.N_total
            @test res_b.cluster_ranges[1][1:(2nb)] == collect(1:(2nb))
            ic_b = load_merger_ic_result(res_b.output_dir)
            @test ic_b.cluster_ranges == res_b.cluster_ranges && ic_b.n_pairs == res_b.n_pairs
            @test ic_b.cluster_specs[1].binaries.fraction == 0.3
            @test parse_merger_summary(joinpath(res_b.output_dir, "merger_summary.txt")) ==
                  res_b.cluster_ranges
            @test occursin(
                "Primordial binaries: NBIN0 = $(nb)",
                read(joinpath(res_b.output_dir, "merger_summary.txt"), String),
            )
            ic_meta_b = Nbody6Dynamics.TOML.parsefile(joinpath(res_b.output_dir, "merger_ic.toml"))
            @test ic_meta_b["meta"]["nbin0"] == nb && length(ic_meta_b["cluster_blocks"][1]) == 2
            # TOML round trip of the binaries table
            bt_path = joinpath(TESTDIR, "merger_bin.toml")
            write(
                bt_path,
                """
[merger]
n_clusters = 2
orbit_mode = "kepler"

[merger.cluster1]
model = "plummer"
N = 100
rbar = 1.0

[merger.cluster1.binaries]
fraction = 0.25
pairing = "uniform_q"
period = "loguniform"
a_min = 0.5
a_max = 50.0
eccentricity = "circular"

[merger.cluster2]
model = "king"
N = 100
rbar = 1.0
""",
            )
            cfg_bt = load_merger_config(bt_path)
            @test cfg_bt.clusters[1].binaries.fraction == 0.25 &&
                  cfg_bt.clusters[1].binaries.pairing == "uniform_q"
            @test cfg_bt.clusters[1].binaries.a_max == 50.0 &&
                  cfg_bt.clusters[2].binaries.fraction == 0.0
        end

        # --- Summary write→parse round-trip (guards the format/regex coupling) ---
        @testset "Merger summary round-trip" begin
            out_dir = mktempdir()
            cfg = MergerConfig(
                [
                    ClusterSpec(
                        profile = KingProfile(W0 = 5.0),
                        N = 150,
                        rbar = 1.0,
                        imf = RescaledKroupaIMF(
                            bodyn = 0.1,
                            body1 = 50.0,
                            target_mass = 150 * kroupa_mean_mass(0.1, 50.0),
                        ),
                    ),
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 100,
                        rbar = 1.0,
                        imf = EqualMassIMF(particle_mass = 8e2 / 100),
                    ),
                ],
                "kepler",
                OrbitSpec(apocentre = 10.0, eccentricity = 0.4),
                MergerOutputSpec(format = "nbody", truncate_jacobi = true, output_dir = out_dir),
            )
            result = generate_merger_ic(cfg; rng = rng)
            ranges = parse_merger_summary(joinpath(out_dir, "merger_summary.txt"))
            @test ranges == result.cluster_ranges
            @test sum(length, ranges) == result.N_total
        end

        # --- MergerPipelineConfig in Nbody6Config ---
        @testset "MergerPipelineConfig in load_config" begin
            cfg = load_config(joinpath(TESTDIR, "test_config.toml"))
            @test cfg.merger.enabled == false
            @test cfg.merger.config_file == ""
        end
    end

    # =====================================================================
    @testset "Per-cluster structure" begin
        rng_s = StableRNG(2024)
        function _synthetic_cluster(N, a, centre, vcm, rng)
            pos, vel = sample_plummer(N, a; rng = rng)
            mass = fill(1.0 / N, N)
            virialise!(mass, pos, vel)
            pos .+= centre
            vel .+= vcm
            return pos, vel, mass
        end
        N1, N2 = 400, 300
        p1, v1, m1 = _synthetic_cluster(N1, 0.2, [-3.0, 0.0, 0.0], [0.0, 0.3, 0.0], rng_s)
        p2, v2, m2 = _synthetic_cluster(N2, 0.15, [3.0, 0.0, 0.0], [0.0, -0.3, 0.0], rng_s)
        # Unbind the last 5 % of cluster 1: radial speed well above escape (√(2/a) ≈ 3.2)
        n_unb = 20
        for j in (N1 - n_unb + 1):N1
            dir = p1[:, j] .- [-3.0, 0.0, 0.0]
            dir ./= max(sqrt(sum(dir .^ 2)), 1e-6)
            v1[:, j] .= [0.0, 0.3, 0.0] .+ 5.0 .* dir
        end
        pos = hcat(p1, p2)
        vel = hcat(v1, v2)
        mass = vcat(m1, m2)
        names = Int32.(1:(N1 + N2))
        function _snap(t)
            params = zeros(Float32, 20)
            params[1] = t
            params[3] = 1.0f0
            params[4] = 1.0f0
            params[11] = 1.0f0
            params[12] = 1.0f0
            Snapshot(
                SnapshotHeader(Int32(N1 + N2), Int32(1), Int32(1), Int32(20), params),
                names,
                Float32.(mass),
                Float32.(pos),
                Float32.(vel),
                Float32[],
                Float32[],
            )
        end
        snaps = [_snap(0.0), _snap(1.0)]
        ranges = [1:N1, (N1 + 1):(N1 + N2)]

        st = cluster_structure(snaps, ranges)
        @test st isa ClusterStructure
        @test st.time == [0.0, 1.0]
        @test st.n_members[:, 1] == [N1, N2]
        @test N1 - n_unb - 20 ≤ st.n_bound[1, 1] ≤ N1 - n_unb
        @test st.n_bound[2, 1] ≥ 0.9 * N2
        @test 0.88 ≤ st.bound_mass_fraction[1, 1] ≤ 0.95
        @test isapprox(st.centre[:, 1, 1], [-3.0, 0.0, 0.0]; atol = 0.05)
        @test isapprox(st.centre[:, 2, 1], [3.0, 0.0, 0.0]; atol = 0.05)
        @test isapprox(st.r_lagr[2, 1, 1], 1.305 * 0.2; rtol = 0.15)
        @test isapprox(st.r_lagr[2, 2, 1], 1.305 * 0.15; rtol = 0.15)
        @test st.r_lagr[1, 1, 1] < st.r_lagr[2, 1, 1] < st.r_lagr[3, 1, 1]
        @test 0.4 ≤ st.q_virial[2, 1] ≤ 0.6
        @test isfinite(st.sigma_1d[1, 1]) && st.sigma_1d[1, 1] > 0
        @test st.sigma_1d[1, 1] > st.sigma_1d[2, 1] * 0.5   # both of order √(M/r)

        # Unbound fast members inflate the all-member ratio; the bound selection removes them
        Q_all, _ = per_cluster_virial(snaps, ranges; bound_only = false)
        Q_bound, nm = per_cluster_virial(snaps, ranges)
        @test Q_all[1, 1] > Q_bound[1, 1]
        @test 0.4 ≤ Q_bound[1, 1] ≤ 0.65
        @test nm[1, 1] == N1
        st_all = cluster_structure(snaps, ranges; bound_only = false)
        @test st_all.n_bound[1, 1] == N1 && st_all.bound_mass_fraction[1, 1] == 1.0

        # Whole-system bound fraction: the 20 kicked stars are unbound to the pair as well
        fb = bound_fraction(snaps[1])
        @test 0.9 ≤ fb ≤ 0.98

        vis_p = VisualizationConfig(;
            format = "png",
            column = "single",
            units = "nbody",
            output_dir = mktempdir(),
        )
        # Radial profiles against the generating models
        rng_p = StableRNG(9090)
        Np = 6000
        pp, vp = sample_plummer(Np, 0.2; rng = rng_p)
        mp = fill(1.0 / Np, Np)
        virialise!(mp, pp, vp)
        prof = radial_profile(pp, vp, mp, zeros(3))
        @test prof isa RadialProfile && length(prof.r) == 12 && sum(prof.n) ≥ 0.98 * Np
        @test isapprox(prof.r_h, 1.305 * 0.2; rtol = 0.08)
        ρ_pl = model_density(PlummerProfile(), 1.0, prof.r_h)
        good = prof.n .≥ 80
        ratios = prof.rho[good] ./ [ρ_pl(x) for x in prof.r[good]]
        @test all(0.7 .≤ ratios .≤ 1.3)
        @test all(abs.(prof.beta[prof.n .≥ 150]) .< 0.25)           # isotropic Plummer
        σ_iso = prof.sigma_r[prof.n .≥ 150]
        @test all(0.3 .< σ_iso .< 3.0)
        @test_throws ArgumentError radial_profile(pp[:, 1:10], vp[:, 1:10], mp[1:10], zeros(3))
        # King W0 = 6, scaled to r_h = 0.3: the model integrates to the mass and half-mass radius
        pk, vk = sample_king(Np, 6.0, 1.0; rng = rng_p)
        mk = fill(1.0 / Np, Np)
        r_hk = Nbody6Dynamics.half_mass_radius(mk, pk; centre = zeros(3))
        pk .*= 0.3 / r_hk
        virialise!(mk, pk, vk)
        profk = radial_profile(pk, vk, mk, zeros(3))
        ρ_k = model_density(KingProfile(W0 = 6.0), 1.0, 0.3)
        goodk = profk.n .≥ 80
        ratk = profk.rho[goodk] ./ [ρ_k(x) for x in profk.r[goodk]]
        @test all(0.65 .≤ ratk .≤ 1.35)
        @test ρ_k(100.0) == 0.0                                       # beyond the tidal radius
        r_grid = exp10.(range(-3, 1; length = 400))
        m_int = sum(
            4π * r_grid[i]^2 * ρ_k(r_grid[i]) * (r_grid[i + 1] - r_grid[i]) for
            i in 1:(length(r_grid) - 1)
        )
        @test isapprox(m_int, 1.0; rtol = 0.05)
        # Per-cluster profiles on the synthetic pair and the system profile
        profs = cluster_profiles(snaps[1], ranges)
        @test length(profs) == 2 && all(!isnothing, profs)
        @test isapprox(profs[2].r_h, 1.305 * 0.15; rtol = 0.15)
        sp = system_profile(snaps[1])
        @test sp.M ≈ sum(mass) && sp.r_h > profs[1].r_h                 # the pair is wider than a member
        # Figures: with and without the generating models, physical off
        specs_s = [
            ClusterSpec(profile = PlummerProfile(), N = N1, rbar = 1.305 * 0.2, imf = KroupaIMF()),
            ClusterSpec(profile = PlummerProfile(), N = N2, rbar = 1.305 * 0.15, imf = KroupaIMF()),
        ]
        pd1 = plot_density_profiles(snaps[1], ranges, vis_p; specs = specs_s)
        @test isfile(pd1)
        pd2 = plot_density_profiles(snaps[1], ranges, vis_p; filename = "density_nomodel")
        @test isfile(pd2)
        pv = plot_velocity_dispersion(snaps[1], ranges, vis_p)
        @test isfile(pv)
        @test_throws ArgumentError plot_density_profiles(
            snaps[1],
            ranges,
            vis_p;
            specs = specs_s[1:1],
        )

        # Figure: with and without the engine overlay (NB units, raster draft)
        vis_s = VisualizationConfig(;
            format = "png",
            column = "single",
            units = "nbody",
            output_dir = mktempdir(),
        )
        path = plot_cluster_structure(snaps, ranges, vis_s)
        @test isfile(path)
        lagr_s = LagrangianData([0.0, 1.0], [0.1, 0.5, 0.9], [0.5 0.5; 3.0 3.0; 6.0 6.0])
        path2 = plot_cluster_structure(
            snaps,
            ranges,
            vis_s;
            lagr = lagr_s,
            filename = "structure_overlay",
        )
        @test isfile(path2)
    end

    # =====================================================================
    @testset "Merger plot suite smoke" begin
        # plot_merger_ic on a small generated IC (all five figure families)
        ic_dir = mktempdir()
        rng_plot = StableRNG(4242)
        cfg_ic = MergerConfig(
            [
                ClusterSpec(
                    profile = PlummerProfile(),
                    N = 60,
                    rbar = 1.0,
                    imf = RescaledKroupaIMF(
                        bodyn = 0.1,
                        body1 = 50.0,
                        target_mass = 60 * kroupa_mean_mass(0.1, 50.0),
                    ),
                ),
                ClusterSpec(
                    profile = KingProfile(W0 = 5.0),
                    N = 60,
                    rbar = 1.0,
                    imf = EqualMassIMF(particle_mass = 5.0),
                ),
            ],
            "kepler",
            OrbitSpec(apocentre = 10.0, eccentricity = 0.5),
            MergerOutputSpec(format = "nbody", truncate_jacobi = false, output_dir = ic_dir),
        )
        ic_result = generate_merger_ic(cfg_ic; rng = rng_plot)
        plots_dir = mktempdir()
        vis_smoke = VisualizationConfig(; format = "png", output_dir = plots_dir)
        plot_merger_ic(ic_result, vis_smoke)
        for stem in (
            "merger_ic_xy",
            "merger_ic_xz",
            "merger_ic_yz",
            "merger_ic_overview",
            "merger_ic_velocity",
            "merger_ic_imf",
            "merger_ic_density",
        )
            @test isfile(joinpath(plots_dir, stem * ".png"))
        end

        # Synthetic two-cluster snapshots: converging COMs, NB units
        function _smoke_snapshot(t, d)
            n_half = 40
            params = zeros(Float32, 20)
            params[1] = Float32(t)
            pos = zeros(Float32, 3, 2 * n_half)
            vel = Float32.(0.05 .* randn(rng_plot, 3, 2 * n_half))
            pos[:, 1:n_half] .= Float32.(0.3 .* randn(rng_plot, 3, n_half))
            pos[1, 1:n_half] .-= Float32(d)
            pos[:, (n_half + 1):end] .= Float32.(0.3 .* randn(rng_plot, 3, n_half))
            pos[1, (n_half + 1):end] .+= Float32(d)
            Snapshot(
                SnapshotHeader(Int32(2 * n_half), Int32(1), Int32(1), Int32(20), params),
                Int32.(1:(2 * n_half)),
                fill(Float32(1 / (2 * n_half)), 2 * n_half),
                pos,
                vel,
                Float32[],
                Float32[],
            )
        end
        snaps = [_smoke_snapshot(t, d) for (t, d) in ((0.0, 4.0), (1.0, 2.0), (2.0, 0.5))]
        ranges = [1:40, 41:80]

        Q, n_mem = per_cluster_virial(snaps, ranges)
        @test size(Q) == (2, 3)
        @test all(n_mem .== 40)
        @test all(isfinite, Q)

        plot_cluster_separation(snaps, ranges, vis_smoke)
        @test isfile(joinpath(plots_dir, "merger_cluster_separation.png"))
        plot_cluster_virial(snaps, ranges, vis_smoke)
        @test isfile(joinpath(plots_dir, "merger_cluster_virial.png"))
        # Membership as body-index blocks (what parse_merger_summary returns)
        Q_blocks, n_blocks = per_cluster_virial(snaps, collect.(ranges))
        @test Q_blocks == Q && n_blocks == n_mem
        plot_cluster_virial(
            snaps,
            collect.(ranges),
            vis_smoke;
            filename = "merger_cluster_virial_blocks",
        )
        @test isfile(joinpath(plots_dir, "merger_cluster_virial_blocks.png"))

        # Envelope statistics helper: NaNs are skipped, all-NaN columns stay NaN
        env = [1.0 NaN 3.0; 5.0 NaN 1.0]
        lo, hi, mean_vals = Nbody6Dynamics._envelope_stats(env)
        @test lo == [1.0, NaN, 1.0] || (lo[1] == 1.0 && isnan(lo[2]) && lo[3] == 1.0)
        @test hi[1] == 5.0 && isnan(hi[2]) && hi[3] == 3.0
        @test mean_vals[1] == 3.0 && isnan(mean_vals[2]) && mean_vals[3] == 2.0
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
        )
            @test haskey(hw, key)
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
    @testset "Multi-panel canvas stays one column wide" begin
        for col in ("single", "double")
            cfg = VisualizationConfig(; column = col)
            pw, ph = Nbody6Dynamics._figsize_px(cfg)
            # Grids keep the column width whatever the column count; stacks are unchanged
            @test Nbody6Dynamics._fig_multipanel(cfg, 2, 3)[1] == pw
            @test Nbody6Dynamics._fig_multipanel(cfg, 1, 2; inner_ticks = false)[1] == pw
            @test Nbody6Dynamics._fig_multipanel(cfg, 2, 1) ==
                  (pw, 2 * ph + Nbody6Dynamics._MULTIPANEL_VGAP)
            # Three boxes of the preset aspect, compact gaps, one decoration strip each way
            w3, h3 = Nbody6Dynamics._fig_multipanel(cfg, 2, 3; inner_ticks = false)
            gap = Nbody6Dynamics._MULTIPANEL_GAP_COMPACT
            strip = Nbody6Dynamics._AXIS_PROTRUSION
            box_w = (pw - strip - 2 * gap) / 3
            @test box_w == Nbody6Dynamics._multipanel_box_width(cfg, 3; inner_ticks = false)
            @test h3 == round(Int, 2 * box_w * ph / pw + gap + strip)
            @test Nbody6Dynamics._fig_multipanel(
                cfg,
                2,
                3;
                inner_ticks = false,
                extra_height = 30,
            )[2] == round(Int, 2 * box_w * ph / pw + gap + strip + 30)
            # Square panels are taller; a reserved colorbar column narrows them
            @test Nbody6Dynamics._fig_multipanel(
                cfg,
                2,
                3;
                inner_ticks = false,
                panel_aspect = 1.0,
            )[2] > h3
            @test Nbody6Dynamics._fig_multipanel(
                cfg,
                2,
                3;
                inner_ticks = false,
                panel_aspect = 1.0,
                extra_width = 110,
            )[2] < Nbody6Dynamics._fig_multipanel(
                cfg,
                2,
                3;
                inner_ticks = false,
                panel_aspect = 1.0,
            )[2]
        end
        @test Nbody6Dynamics._multipanel_gap(3; inner_ticks = false) ==
              Nbody6Dynamics._MULTIPANEL_GAP_COMPACT
        @test Nbody6Dynamics._multipanel_gap(3; inner_ticks = true) ==
              Nbody6Dynamics._MULTIPANEL_HGAP
        @test Nbody6Dynamics._multipanel_gap(1; inner_ticks = false) ==
              Nbody6Dynamics._MULTIPANEL_HGAP
        # Marker scale follows the panel width; the annotation band is a data-free strip
        cfg_s = VisualizationConfig(; column = "single")
        @test Nbody6Dynamics._multipanel_scale(cfg_s, 1) == 1.0
        s3 = Nbody6Dynamics._multipanel_scale(cfg_s, 3; inner_ticks = false)
        @test 0.25 < s3 < 1 / 3
        @test Nbody6Dynamics._multipanel_scale(cfg_s, 3; inner_ticks = false, extra_width = 110) <
              s3
        @test Nbody6Dynamics._multipanel_scale(cfg_s, 3; inner_ticks = true) < s3
        @test 0 < Nbody6Dynamics._MONTAGE_BAND_FRAC < 0.5
    end

    # =====================================================================
    # Edge cases for pure helpers and degenerate reader inputs
    # =====================================================================
    @testset "Degenerate axis ranges" begin
        # Identical stars (equal-mass, unevolved) give zero-span HR data
        @test Nbody6Dynamics._padded_range(3.678, 3.678) == (3.578, 3.778)
        lo, hi = Nbody6Dynamics._padded_range(3.0, 4.0)
        @test lo ≈ 2.94 && hi ≈ 4.06
        @test !isempty(Nbody6Dynamics._logval_ticks(3.62, 3.74))
        @test Nbody6Dynamics._logval_ticks(3.678, 3.678) == [3.678]
        @test !isempty(Nbody6Dynamics._nice_ticks(1.0, 1.0001))
        @test Nbody6Dynamics._nice_ticks(0.0, 10.0) == collect(0.0:1.0:10.0) ||
              !isempty(Nbody6Dynamics._nice_ticks(0.0, 10.0))
    end

    @testset "Log-tick generator edge cases" begin
        # >2 in-range decades → decades only
        vals, _ = Nbody6Dynamics._log_ticks(0.05, 50.0)
        @test vals == [0.1, 1.0, 10.0]
        # ≤2 in-range decades → 2×/5× intermediates appear
        vals2, labels2 = Nbody6Dynamics._log_ticks(0.5, 30.0)
        @test all(v -> v in vals2, (0.5, 1.0, 2.0, 5.0, 10.0, 20.0))
        # 10^0 renders as plain "1"
        lab1 = String(labels2[findfirst(==(1.0), vals2)])
        @test occursin("1", lab1) && !occursin("10", lab1)
        # Degenerate equal endpoints still yield ≥ 2 ticks
        vals3, _ = Nbody6Dynamics._log_ticks(2.0, 2.0)
        @test length(vals3) ≥ 2
        # Plain decimals throughout on a short axis within 10⁻³–10⁴ (the
        # virial-ratio panel: 0.5 … 10, never "5 × 10⁻¹ … 10¹")
        vals4, labels4 = Nbody6Dynamics._log_ticks(0.45, 13.0)
        @test vals4 == [0.5, 1.0, 2.0, 5.0, 10.0]
        @test [Nbody6Dynamics._log_tick_label(v, true) for v in vals4] == ["0.5", "1", "2", "5", "10"]
        @test all(l -> !occursin("times", String(l)) && !occursin("^", String(l)), labels4)
        @test [
            Nbody6Dynamics._log_tick_label(v, true) for v in (0.001, 0.01, 0.1, 100.0, 50000.0)
        ] == ["0.001", "0.01", "0.1", "100", "50000"]
        # Exponent form beyond the plain window, with the mandatory collapses
        _, labels5 = Nbody6Dynamics._log_ticks(1e-8, 1e-3)
        @test occursin("10^{-8}", String(labels5[1])) && occursin("10^{-3}", String(labels5[end]))
        _, labels6 = Nbody6Dynamics._log_ticks(0.5, 1e6)
        @test occursin("1", String(labels6[1])) && !occursin("10", String(labels6[1]))
        @test String(labels6[2]) == "\$10\$" && occursin("10^{2}", String(labels6[3]))
        @test Nbody6Dynamics._log_tick_label(2e-7, false) == "2\\times 10^{-7}"
        @test Nbody6Dynamics._log_tick_label(0.2, false) == "0.2"
        @test Nbody6Dynamics._log_tick_label(20.0, false) == "20"
        @test Nbody6Dynamics._log_tick_label(200.0, false) == "2\\times 10^{2}"
    end

    @testset "Degenerate reader inputs" begin
        # sev.83 with zero stars: header-only file parses to an empty snapshot
        p_empty = joinpath(TESTDIR, "sev.83_empty")
        write(p_empty, "  0  0.0\n")
        sev0 = read_stellar_evolution(p_empty)
        @test sev0.n_stars == 0 && isempty(sev0.records)

        # esc.11 line with only 12 tokens (angle column truncated) is skipped
        p_trunc = joinpath(TESTDIR, "esc_trunc.11")
        write(p_trunc, "  1.0 5.0e-4 20.0 40.0 2.0e-3 1.234 0.5 -0.1 15.6 0 101 22.2\n")
        @test isempty(read_escapers(p_trunc))
    end

    @testset "half_mass_radius boundaries" begin
        # Single particle: the half-mass radius is that particle's radius
        @test Nbody6Dynamics.half_mass_radius([2.0], zeros(3, 1); centre = zeros(3)) == 0.0
        # Two equal masses: cumulative ≥ M/2 at the inner particle
        pos2 = [1.0 2.0; 0.0 0.0; 0.0 0.0]
        @test Nbody6Dynamics.half_mass_radius([0.5, 0.5], pos2; centre = zeros(3)) ≈ 1.0
    end

    @testset "Kepler near-parabolic limit" begin
        # Apocentre speed → 0 as e → 1 (vis-viva); must stay finite/nonneg
        v1, v2 = kepler_velocity(1.0, 1.0, 10.0, 1.0 - 1e-12)
        @test 0.0 ≤ v1 < 1e-5 && 0.0 ≤ v2 < 1e-5
    end

    # =====================================================================
    # Static QA ships with the tests.
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
        script0 =
            read(Nbody6Dynamics._write_launch_script(launch_dir, args..., make_cfg(0)), String)
        @test !occursin("OMP_NUM_THREADS", script0)
        @test occursin("OMP_STACKSIZE", script0)
        script4 =
            read(Nbody6Dynamics._write_launch_script(launch_dir, args..., make_cfg(4)), String)
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
        @test isnan(Nbody6Dynamics._dump_time("comm.1")) &&
              isnan(Nbody6Dynamics._dump_time("conf.3_1"))
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
        s_app = read(
            Nbody6Dynamics._write_launch_script(run_dir, args..., cfg_r; append = true),
            String,
        )
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
            read(wd_out, String) *
            " ADJUST:  TIME    5.00000E-01  T[Myr]   0.338E+01  Q   0.457E+00\n",
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
        @test all(
            l -> count(==(','), l) == fieldcount(Nbody6Dynamics.TelemetrySample) - 1,
            csv_lines,
        )
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
    @testset "Threaded snapshot reading" begin
        dir = mktempdir()
        write_snap(path, t, n) = open(path, "w") do io
            _write_fortran_record(io, Int32[n, 1, 1, 20])
            params = zeros(Float32, 20)
            params[1] = Float32(t)
            params[3] = 1.0f0
            params[4] = 0.5f0
            params[11] = 10.0f0
            params[12] = 5.0f0
            params[18] = 2.0f0
            _write_fortran_record(io, params)
            for i in 1:n
                buf = IOBuffer()
                write(
                    buf,
                    Float32(1 / n),
                    Float32(i),
                    Float32(-i),
                    0.0f0,
                    0.1f0,
                    0.0f0,
                    0.0f0,
                    Int32(i),
                )
                _write_fortran_record(io, take!(buf))
            end
        end
        for (k, t) in enumerate((0.0, 0.5, 1.0, 1.5, 2.0))
            write_snap(joinpath(dir, "conf.3_$(t)"), t, 4 + k)
        end
        write(joinpath(dir, "conf.3_2.5"), "corrupt")
        serial =
            @test_logs (:warn, r"Skipping corrupt snapshot") (:warn, r"1 / 6 snapshot") read_all_conf3(
                dir;
                threaded = false,
            )
        threaded =
            @test_logs (:warn, r"Skipping corrupt snapshot") (:warn, r"1 / 6 snapshot") read_all_conf3(
                dir;
                threaded = true,
            )
        @test length(serial) == length(threaded) == 5
        @test [time_nb(s.header) for s in threaded] == [0.0, 0.5, 1.0, 1.5, 2.0]
        @test all(
            a.pos == b.pos && a.mass == b.mass && a.name == b.name for
            (a, b) in zip(serial, threaded)
        )
        @test [nparticles(s) for s in threaded] == [5, 6, 7, 8, 9]
        # The ordered reader is generic over the element type and the reader
        names = Nbody6Dynamics._read_ordered(
            String,
            p -> uppercase(basename(p)),
            dir,
            ["conf.3_0.0", "conf.3_1.0"];
            desc = "",
            what = "name",
            threaded = true,
        )
        @test names == ["CONF.3_0.0", "CONF.3_1.0"]
        @test occursin(
            "--threads=$(Threads.nthreads())",
            string(Nbody6Dynamics._sweep_worker_command(dir)),
        )
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
        write(
            one,
            " ADJUST:  TIME    0.00000E+00  T[Myr]   0.000E+00  Q   0.450E+00  DE   0.000E+00\n",
        )
        @test Nbody6Dynamics._live_diagnostics_panel(one) === nothing
        @test Nbody6Dynamics._live_diagnostics_panel(joinpath(dir, "absent")) === nothing
        # Config keys
        @test SimulationConfig().live_diagnostics == false &&
              SimulationConfig().live_interval == 30.0
        cfg_path = joinpath(dir, "c.toml")
        write(cfg_path, "[simulation]\nlive_diagnostics = true\nlive_interval = 5\n")
        c = load_config(cfg_path)
        @test c.simulation.live_diagnostics && c.simulation.live_interval == 5.0
        write(cfg_path, "[simulation]\nlive_interval = 0.5\n")
        @test_throws ErrorException load_config(cfg_path)
    end

    # =====================================================================
    @testset "Telemetry readers and figure" begin
        dir = mktempdir()
        header = join(string.(fieldnames(TelemetrySample)), ",")
        row(
            i;
            gpu = "NaN",
        ) = "$(5.0 * i),3,$(200 + i),$(210 + i),$(10.0 * i),$(2.0 + 0.1 * i),1.5,$gpu,$gpu,$gpu,$gpu,$gpu"
        write(
            joinpath(dir, "telemetry.csv"),
            header * "\n" * join([row(i) for i in 1:6], "\n") * "\n",
        )
        s = read_telemetry(joinpath(dir, "telemetry.csv"))
        @test length(s) == 6 &&
              s[1].elapsed_s == 5.0 &&
              s[6].rss_mib == 206.0 &&
              s[1].n_processes == 3
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
        @test read_telemetry(joinpath(dir, "empty.csv") |> p -> (write(p, ""); p)) ==
              TelemetrySample[]
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
        @test Nbody6Dynamics._telemetry_time_axis([0.0, 60.0])[1] == [0.0, 60.0]
        @test Nbody6Dynamics._telemetry_time_axis([0.0, 600.0])[1] == [0.0, 10.0]
        @test Nbody6Dynamics._telemetry_time_axis([0.0, 7200.0])[1] == [0.0, 2.0]
    end

    # =====================================================================
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
            base_toml(
                runs_dir,
                extra_sim = "",
                merger = "enabled = false\nconfig_file = \"\"",
            ) = """
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
                strip.(
                    split(read(`nvidia-smi --query-gpu=index --format=csv,noheader`, String), '\n'),
                ),
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
    @testset "Static QA — Aqua" begin
        using Aqua
        # Method ambiguities are checked for this package only — recursing
        # into the Makie/SciML dependency tree reports upstream noise.
        Aqua.test_all(Nbody6Dynamics; ambiguities = false, persistent_tasks = false)
        Aqua.test_ambiguities(Nbody6Dynamics)
    end

    @testset "Static QA — ExplicitImports" begin
        using ExplicitImports
        # Pragmatic subset: the package uses plain `using` for its small,
        # stable dependency surface (full explicit-import migration is
        # tracked in the roadmap). These checks catch the real hazards:
        # stale explicit imports, self-qualified names, and accesses of
        # non-owning modules.
        @test check_no_stale_explicit_imports(Nbody6Dynamics) === nothing
        @test check_no_self_qualified_accesses(Nbody6Dynamics) === nothing
    end

    @testset "Static QA — JET" begin
        using JET
        # Reports scoped to this package's own frames — Base/dependency
        # internals (e.g. @sync's sync_end, tuple broadcasting) produce
        # known false positives outside our control.
        JET.test_package(Nbody6Dynamics; target_modules = (Nbody6Dynamics,))
    end
end  # top-level testset
