# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

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

    # Git state that cannot restore the files → set aside, never
    # deleted, and cloned again
    rm(joinpath(src, "configure"))
    rm(joinpath(src, ".git", "HEAD"))
    @test_logs (:warn, r"cannot be restored") match_mode = :any Nbody6Dynamics._ensure_source_tree(
        src,
        inst,
    )
    @test Nbody6Dynamics._source_tree_ready(src)
    aside = filter(startswith("engine.incomplete-"), readdir(dirname(src)))
    @test length(aside) == 1
    @test read(joinpath(dirname(src), only(aside), "untracked.txt"), String) == "keep me"

    # A non-empty directory that is no checkout at all (a mistyped
    # install_dir) is refused and left untouched
    foreign = joinpath(mktempdir(), "results")
    mkpath(foreign)
    write(joinpath(foreign, "data.csv"), "1,2,3")
    @test_throws ArgumentError Nbody6Dynamics._ensure_source_tree(foreign, inst)
    @test read(joinpath(foreign, "data.csv"), String) == "1,2,3"

    # An empty directory is simply cloned into
    empty_dir = joinpath(mktempdir(), "engine")
    mkpath(empty_dir)
    Nbody6Dynamics._ensure_source_tree(empty_dir, inst)
    @test Nbody6Dynamics._source_tree_ready(empty_dir)
end

@testset "GPU build target" begin
    # configure arguments: the engine's configure finds nvcc only on PATH
    # (its --with-cuda fallback reuses the cached PATH check), so the
    # build exports the toolkit's bin and names the prefix explicitly.
    build_of(body) = (p = joinpath(mktempdir(), "cfg.toml"); write(p, body); load_config(p).build)
    cpu_args =
        Nbody6Dynamics._configure_args(build_of("[build]\nenable_gpu = false\n"), "/usr/local/cuda")
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
    @test Nbody6Dynamics._parse_nvcc_release("Cuda compilation tools, release 12.8, V12.8.93") ==
          "12.8"
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
    @test Nbody6Dynamics._check_cuda_arch_support(["sm_90"], ["sm_80", "sm_90"], "12.4") === nothing
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
    @test resolve_cuda_arch(BuildConfig(; enable_gpu = true, cuda_arch = ["sm_90", "sm_120"])) ==
          ["sm_90", "sm_120"]
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
    @test_throws ErrorException load_config(write_cfg("[build]\ncuda_arch = [\"compute_90\"]\n"))
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
    @test occursin("cudaDevAttrClockRate", helper_src) && occursin("checkCudaErrors", helper_src)
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
    write(err, "[R.0 AVX Reg.F ] x Perf.(Gflops) 100.0\n[R.0 GPU Reg.F ] x Perf.(Gflops) 300.0\n")
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
