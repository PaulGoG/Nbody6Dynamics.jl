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
pair_sum_max_n = 5000

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
    @test cfg.postprocess.pair_sum_max_n == 5000
    @test PostprocessConfig().pair_sum_max_n == 2_000_000
    @test PostprocessConfig().snapshot_integrals
    @test PostprocessConfig().bound_max_passes == 100
    @test PostprocessConfig().bound_min_members == 10
    @test PostprocessConfig().energy_events && PostprocessConfig().mass_budget
    @test PostprocessConfig().energy_jump_threshold == 1.0e-4
    @test PostprocessConfig().energy_jump_sigma == 5.0
    @test PostprocessConfig().events_file == "event.35"
    pp = PostprocessConfig()
    @test (pp.centre_n_stop, pp.centre_stop_fraction, pp.centre_n_floor, pp.centre_shrink) ==
          (1000, 0.01, 500, 0.025)
    @test pp.lagr_n_min == 100 && pp.profile_nbins == 12
    @test (pp.profile_n_min_sigma, pp.profile_n_min_beta, pp.profile_n_min_density) ==
          (50, 300, 100)
    @test pp.density_neighbours == 6 && pp.msr_n_tracer == 50
    @test pp.msr_scan == [10, 20, 50, 100, 200] && pp.msr_n_random == 200
    @test (pp.msr_n_random_small, pp.msr_n_small, pp.msr_seed, pp.msr_threshold) ==
          (1000, 20, 0, 2.0)
    @test (pp.rotation_n_min, pp.rotation_nbins_max, pp.rotation_n_null) == (1000, 16, 20)
    @test (pp.shape_n_min_shell, pp.shape_n_min_enclosed) == (3000, 7000)
    @test pp.pair_k == 10 && pp.pair_tidal_factor ≈ 1 / 3
    @test pp.relaxation_gamma == [0.02, 0.11]
    @test cfg.visualization.dpi == 150
    @test cfg.visualization.export_width == 5.0
    @test cfg.install.source_url == "https://github.com/PaulGoG/Nbody6PPGPU-beijing.git"
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
    @test InstallConfig().ref == "v2026.07.10"
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
        ("pair_sum_max_n", "[postprocess]\npair_sum_max_n = -1\n", "postprocess.pair_sum_max_n"),
        ("span_n_periods", "[postprocess]\nspan_n_periods = 0.0\n", "postprocess.span_n_periods"),
        (
            "bound_gamma_max",
            "[postprocess]\nbound_gamma_max = 0.0\n",
            "postprocess.bound_gamma_max",
        ),
        ("bound_n_miss", "[postprocess]\nbound_n_miss = -1\n", "postprocess.bound_n_miss"),
        (
            "population_frame",
            "[postprocess]\npopulation_frame = \"shells\"\n",
            "postprocess.population_frame",
        ),
        (
            "population_m1_max",
            "[postprocess]\npopulation_m1_min = 0.9\npopulation_m1_max = 0.5\n",
            "postprocess.population_m1_max",
        ),
        (
            "population_edges",
            "[postprocess]\npopulation_edges = [0.0, 0.6, 0.5, 1.0]\n",
            "postprocess.population_edges",
        ),
        (
            "bound_max_passes",
            "[postprocess]\nbound_max_passes = 0\n",
            "postprocess.bound_max_passes",
        ),
        (
            "bound_min_members",
            "[postprocess]\nbound_min_members = 1\n",
            "postprocess.bound_min_members",
        ),
        (
            "energy_jump_threshold",
            "[postprocess]\nenergy_jump_threshold = -1.0e-4\n",
            "postprocess.energy_jump_threshold",
        ),
        (
            "energy_jump_sigma",
            "[postprocess]\nenergy_jump_sigma = -5.0\n",
            "postprocess.energy_jump_sigma",
        ),
        (
            "events_file",
            "[postprocess]\nmass_budget = true\nevents_file = \"\"\n",
            "postprocess.events_file",
        ),
        ("centre_shrink", "[postprocess]\ncentre_shrink = 1.0\n", "postprocess.centre_shrink"),
        ("centre_n_stop", "[postprocess]\ncentre_n_stop = 0\n", "postprocess.centre_n_stop"),
        ("msr_scan", "[postprocess]\nmsr_scan = [20, 10]\n", "postprocess.msr_scan"),
        ("msr_scan_empty", "[postprocess]\nmsr_scan = []\n", "postprocess.msr_scan"),
        (
            "profile_n_min_beta",
            "[postprocess]\nprofile_n_min_sigma = 50\nprofile_n_min_beta = 10\n",
            "postprocess.profile_n_min_beta",
        ),
        (
            "relaxation_gamma",
            "[postprocess]\nrelaxation_gamma = [0.02, -0.1]\n",
            "postprocess.relaxation_gamma",
        ),
        ("pair_k", "[postprocess]\npair_k = 0\n", "postprocess.pair_k"),
        ("rotation_n_null", "[postprocess]\nrotation_n_null = -1\n", "postprocess.rotation_n_null"),
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
            "snapshot_render",
            "[visualization.style]\nsnapshot_render = \"hexbin\"\n",
            "visualization.style.snapshot_render",
        ),
        (
            "density_min_n",
            "[visualization.style]\ndensity_min_n = 0\n",
            "visualization.style.density_min_n",
        ),
        (
            "density_bins",
            "[visualization.style]\ndensity_bins = 8\n",
            "visualization.style.density_bins",
        ),
        (
            "density_decades",
            "[visualization.style]\ndensity_decades = 0.0\n",
            "visualization.style.density_decades",
        ),
        (
            "density_mass_frac",
            "[visualization.style]\ndensity_mass_frac = 1.5\n",
            "visualization.style.density_mass_frac",
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
    # A bare name resolves across the folders by purpose (the flat names of 0.3.0 stay valid).
    @test example_input("N1k_quick.inp") == example_input("engine/N1k_quick.inp")
    @test example_input("merger_demo_small.toml") == normpath(
        joinpath(Nbody6Dynamics._PACKAGE_ROOT, "input_files", "mergers", "merger_demo_small.toml"),
    )
end

# =====================================================================

@testset "Every shipped configuration loads" begin
    root = joinpath(Nbody6Dynamics._PACKAGE_ROOT, "input_files")
    # The layout by purpose: nothing loose at the top, every folder present.
    @test isempty(
        filter(isfile, readdir(root; join = true)) ∩
        filter(p -> endswith(p, ".toml") || endswith(p, ".inp"), readdir(root; join = true)),
    )
    for sub in ("engine", "clusters", "mergers", "verification", "sweeps", "showcase", "gpu")
        @test isdir(joinpath(root, sub))
    end
    # Every TOML parses with the loader its top-level table selects, and the
    # relative references inside the sweep files point at files that ship.
    n = 0
    for (dir, _, files) in walkdir(root), f in sort(files)
        endswith(f, ".toml") || continue
        path = joinpath(dir, f)
        raw = Nbody6Dynamics.TOML.parsefile(path)
        if haskey(raw, "sweep")
            @test load_sweep_config(path) isa Nbody6Dynamics.SweepConfig
        elseif haskey(raw, "install") || haskey(raw, "simulation")
            @test load_config(path) isa Nbody6Dynamics.Nbody6Config
        else
            @test load_merger_config(path) isa Nbody6Dynamics.MergerConfig
        end
        n += 1
    end
    @test n ≥ 30
    # The engine inputs: the size ladder with engine-generated initial
    # conditions, every one with a dump at every adjustment and, where a
    # point-mass field is requested, its namelist.
    inps = filter(f -> endswith(f, ".inp"), readdir(joinpath(root, "engine")))
    @test length(inps) == 5
    for f in inps
        text = read(joinpath(root, "engine", f), String)
        @test occursin(",NCOMM=1,", text)
        kz = Dict{Int,Int}()
        for m in eachmatch(r"KZ\((\d+):(\d+)\)=\s*([-\d ]+)", text)
            vals = parse.(Int, split(strip(m.captures[3])))
            for (k, v) in zip(parse(Int, m.captures[1]):parse(Int, m.captures[2]), vals)
                kz[k] = v
            end
        end
        @test kz[2] == 1 && kz[22] in (0, 1)   # engine-generated initial conditions
        @test kz[8] != 2                       # pairs from dat.10 need KZ(22) ≥ 2
        @test (kz[14] == 0) == !occursin("&INXTRNL0", text)
        nbin0 = parse(Int, match(r"NBIN0=(\d+)", text).captures[1])
        @test (nbin0 > 0) == (kz[8] == 1)
    end
    # The single-cluster science cases: single mode, baseline numerics.
    for f in readdir(joinpath(root, "clusters"))
        m = load_merger_config(joinpath(root, "clusters", f))
        @test m.orbit_mode == "single" && length(m.clusters) == 1
        @test m.nbody6.preset == "baseline" && m.nbody6.checkpoint
        @test m.output.tcrit_myr > 0 && m.output.dtadj_myr > 0 && m.output.deltat_myr > 0
    end
end

# =====================================================================

@testset "Cluster and cost sections" begin
    # Field-wise equality: `==` on these structs compares their vectors by identity.
    same(a, b) = all(getfield(a, k) == getfield(b, k) for k in fieldnames(typeof(a)))
    # Defaults and round trip
    cfg0 = load_config(joinpath(TESTDIR, "empty_sections.toml") |> p -> (write(p, ""); p))
    @test same(cfg0.cluster, ClusterConfig()) && same(cfg0.cost, CostConfig())
    @test !cfg0.cluster.enabled && cfg0.cost.enabled
    rt = joinpath(TESTDIR, "sections_rt.toml")
    write(
        rt,
        """
        [simulation]
        wall_budget = 7200.0
        stop_margin = 300.0

        [cluster]
        enabled = true
        partition = "GPU"
        gres = "gpu:1"
        mem = "64G"
        time_limit = "1-12:00:00"
        partition_limit = "1-12:00:00"
        cpus_per_task = 16
        sbatch_options = ["--hint=nomultithread"]
        prologue = ["module load cuda", "export JULIA_DEPOT_PATH=/data/depot:"]
        chain = "dependency"
        max_segments = 3
        overhead = 600.0

        [cost]
        table = [[1.0e5, 64.0], [5.0e5, 281.0], [1.0e6, 1087.0]]
        budget_hours = 100.0
        on_exceed = "refuse"
        """,
    )
    cfg = load_config(rt)
    @test cfg.cluster.partition == "GPU" && cfg.cluster.gres == "gpu:1"
    @test cfg.cluster.prologue == ["module load cuda", "export JULIA_DEPOT_PATH=/data/depot:"]
    @test cfg.cost.table == [[1.0e5, 64.0], [5.0e5, 281.0], [1.0e6, 1087.0]]
    rt2 = joinpath(TESTDIR, "sections_rt2.toml")
    save_config(cfg, rt2)
    cfg2 = load_config(rt2)
    @test same(cfg2.cluster, cfg.cluster) && same(cfg2.cost, cfg.cost)
    # Absolute-path form keeps both sections
    abs_cfg = Nbody6Dynamics._with_absolute_paths(cfg, TESTDIR)
    @test same(abs_cfg.cluster, cfg.cluster) && same(abs_cfg.cost, cfg.cost)
    # Positional constructors
    six = Nbody6Config(
        cfg.install,
        cfg.build,
        cfg.simulation,
        cfg.postprocess,
        cfg.visualization,
        cfg.merger,
    )
    @test same(six.cluster, ClusterConfig()) && six.config_dir == pwd()
    seven = Nbody6Config(
        cfg.install,
        cfg.build,
        cfg.simulation,
        cfg.postprocess,
        cfg.visualization,
        cfg.merger,
        "/x",
    )
    @test same(seven.cost, CostConfig()) && seven.config_dir == "/x"

    # Slurm time formats
    sec = Nbody6Dynamics._slurm_seconds
    @test sec("90") == 5400
    @test sec("1:30") == 90
    @test sec("1:30:00") == 5400
    @test sec("2-12") == 2 * 86400 + 12 * 3600
    @test sec("2-12:30") == 2 * 86400 + 12 * 3600 + 1800
    @test sec("2-12:30:15") == 2 * 86400 + 12 * 3600 + 1815
    @test sec("36:00:00") == 129600
    for bad in ("", "1:2:3:4", "12h", "1-", "a-1", "1:xx")
        @test_throws ArgumentError sec(bad)
    end
    @test Nbody6Dynamics._slurm_time(129600) == "1-12:00:00"
    @test Nbody6Dynamics._slurm_time(5400) == "01:30:00"
    @test Nbody6Dynamics._slurm_time(sec("2-12:30:15")) == "2-12:30:15"
    @test_throws ArgumentError Nbody6Dynamics._slurm_time(-1)

    # Refusals name the key; a disabled section is not checked
    refused(body, key) = begin
        p = joinpath(TESTDIR, "sections_bad.toml")
        write(p, body)
        err = try
            load_config(p)
            nothing
        catch e
            e
        end
        err isa ArgumentError && occursin(key, err.msg)
    end
    base = "[simulation]\nwall_budget = 7200.0\nstop_margin = 300.0\n\n[cluster]\nenabled = true\n"
    @test refused(base * "scheduler = \"pbs\"\n", "cluster.scheduler")
    @test refused(base * "time_limit = \"12h\"\n", "cluster.time_limit")
    @test refused(
        base * "time_limit = \"2-00:00:00\"\npartition_limit = \"1-12:00:00\"\n",
        "cluster.time_limit",
    )
    @test refused(base * "partition_limit = \"soon\"\n", "cluster.partition_limit")
    @test refused(base * "chain = \"array\"\n", "cluster.chain")
    @test refused(base * "chain = \"dependency\"\n", "cluster.max_segments")
    @test refused(base * "mem = \"64 GB\"\n", "cluster.mem")
    @test refused(base * "sbatch_options = [\"#SBATCH --exclusive\"]\n", "cluster.sbatch_options")
    @test refused(base * "cpus_per_task = 0\n", "cluster.cpus_per_task")
    @test refused(base * "sbatch = \"\"\n", "cluster.sbatch")
    @test refused("[cluster]\nenabled = true\n", "cluster.enabled")               # no wall budget
    @test refused(base * "time_limit = \"2:00:00\"\n", "cluster.time_limit")     # 7200 + 300 + 1800 > 7200
    @test !refused(base * "time_limit = \"2:35:00\"\n", "cluster")              # 9300 ≤ 9300
    @test !refused("[cluster]\nscheduler = \"pbs\"\nchain = \"array\"\nmem = \"x\"\n", "cluster")
    @test refused("[cost]\ntable = [[1.0e5, 64.0]]\n", "cost.table")
    @test refused("[cost]\ntable = [[1.0e5, 64.0, 1.0], [5.0e5, 281.0]]\n", "cost.table")
    @test refused("[cost]\ntable = [[1.0e5, -64.0], [5.0e5, 281.0]]\n", "cost.table")
    @test refused("[cost]\ntable = [[5.0e5, 281.0], [1.0e5, 64.0]]\n", "cost.table")
    @test refused("[cost]\ntable = [1.0e5, 64.0]\n", "cost.table")
    @test refused("[cost]\nmerger_factor = 0.0\n", "cost.merger_factor")
    @test refused("[cost]\nbudget_hours = -1.0\n", "cost.budget_hours")
    @test refused("[cost]\non_exceed = \"stop\"\n", "cost.on_exceed")
    @test refused("[cluster]\nenable = true\n", "cluster.enable")
end
