# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

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
export_width = 5.0
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
    @test cfg.visualization.export_width == 5.0
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
export_width = 7.0
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
    @test cfg2.visualization.export_width == cfg.visualization.export_width
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
        ("startup_timeout", "[simulation]\nstartup_timeout = -1.0\n", "simulation.startup_timeout"),
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
        ("export_width", "[visualization]\nexport_width = 0.0\n", "visualization.export_width"),
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
        ("zoom_frac", "[visualization.style]\nzoom_frac = 0.0\n", "visualization.style.zoom_frac"),
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

    # Unknown keys and wrong types are refused, naming the key.
    unknown = joinpath(TESTDIR, "unknown_key.toml")
    write(unknown, "[simulation]\nomp_thread = 4\n")
    err = try
        load_config(unknown)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("simulation.omp_thread", err.msg)
    write(unknown, "[simulaton]\nomp_threads = 4\n")
    err = try
        load_config(unknown)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("simulaton", err.msg)
    write(unknown, "[simulation]\nomp_threads = \"four\"\n")
    err = try
        load_config(unknown)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError &&
          occursin("simulation.omp_threads", err.msg) &&
          occursin("integer", err.msg)
    write(unknown, "[build]\ncuda_arch = [\"sm_90\", 90]\n")
    err = try
        load_config(unknown)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("build.cuda_arch", err.msg)
    write(unknown, "[visualization.style]\nmarker_budgett = 1.0\n")
    err = try
        load_config(unknown)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError && occursin("visualization.style.marker_budgett", err.msg)
    # Defaults come from the structs: an empty file is the default configuration.
    write(unknown, "")
    @test load_config(unknown).simulation.omp_threads == SimulationConfig().omp_threads
end

# =====================================================================

# =====================================================================
@testset "Project directory and path resolution" begin
    proj = mktempdir()
    cfg_path = joinpath(proj, "config.toml")
    write(
        cfg_path,
        """
[install]
enabled = false
install_dir = "engine/tree"

[simulation]
run_test = false
input_file = "inputs/case.inp"
runs_dir = "results"

[postprocess]
data_dir = "external/output"

[merger]
enabled = false
config_file = "merger.toml"
""",
    )
    cfg = load_config(cfg_path)
    # The directory of the file is the resolution base and the default base_dir.
    @test cfg.config_dir == proj
    # A configuration built in memory takes the working directory.
    @test Nbody6Config(
        InstallConfig(),
        BuildConfig(),
        SimulationConfig(),
        PostprocessConfig(),
        VisualizationConfig(),
        MergerPipelineConfig(),
    ).config_dir == pwd()
    # config_dir is not a configuration key: it survives save/load through the path only.
    frozen = joinpath(proj, "sub", "frozen.toml")
    mkpath(dirname(frozen))
    save_config(cfg, frozen)
    @test !haskey(Nbody6Dynamics.TOML.parsefile(frozen), "config_dir")
    @test load_config(frozen).config_dir == dirname(frozen)

    # One resolution base for every relative path; absolute paths are honoured.
    @test Nbody6Dynamics._resolve_path(proj, "a/b") == joinpath(proj, "a", "b")
    @test Nbody6Dynamics._resolve_path(proj, "/abs/c") == "/abs/c"
    @test Nbody6Dynamics._resolve_path(proj, "") == ""
    @test Nbody6Dynamics._resolve_path(joinpath(proj, "x"), "../y") == joinpath(proj, "y")
    abs_cfg = Nbody6Dynamics._with_absolute_paths(cfg, proj)
    @test abs_cfg.install.install_dir == joinpath(proj, "engine", "tree")
    @test abs_cfg.simulation.input_file == joinpath(proj, "inputs", "case.inp")
    @test abs_cfg.simulation.runs_dir == joinpath(proj, "results")
    @test abs_cfg.postprocess.data_dir == joinpath(proj, "external", "output")
    @test abs_cfg.merger.config_file == joinpath(proj, "merger.toml")
    @test abs_cfg.config_dir == proj
    @test abs_cfg.build == cfg.build && abs_cfg.visualization.style == cfg.visualization.style
    # An already absolute form is a fixed point.
    again = Nbody6Dynamics._with_absolute_paths(abs_cfg, "/elsewhere")
    @test again.install.install_dir == abs_cfg.install.install_dir
    # An empty data_dir stays empty (it means "no external directory").
    @test Nbody6Dynamics._with_absolute_paths(load_config(frozen), proj).postprocess.data_dir ==
          joinpath(proj, "external", "output")
    plain = joinpath(proj, "plain.toml")
    write(plain, "[postprocess]\ndata_dir = \"\"\n")
    @test Nbody6Dynamics._with_absolute_paths(load_config(plain), proj).postprocess.data_dir == ""

    # The latest-run lookup and the input file follow base_dir, not the package tree.
    mkpath(joinpath(proj, "results", "run_20260101_000000_abcd"))
    @test Nbody6Dynamics._find_latest_run(cfg, proj) ==
          joinpath(proj, "results", "run_20260101_000000_abcd")
    @test isempty(Nbody6Dynamics._find_latest_run(cfg, mktempdir()))
    err = try
        run_simulation(cfg)
    catch e
        e
    end
    @test err isa ErrorException && occursin(joinpath(proj, "inputs", "case.inp"), err.msg)

    # Shipped inputs are reachable without a checkout.
    @test isfile(example_input("N1k_quick.inp"))
    @test example_input("showcase/equal_pipeline.toml") == normpath(
        joinpath(Nbody6Dynamics._PACKAGE_ROOT, "input_files", "showcase", "equal_pipeline.toml"),
    )
    @test_throws ArgumentError example_input("does_not_ship.inp")
end
