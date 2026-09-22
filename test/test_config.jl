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
end

# =====================================================================
