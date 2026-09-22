module Nbody6Dynamics

using TOML
using Dates
using Printf
using LaTeXStrings
using ProgressMeter
using Random
using SpecialFunctions
using OrdinaryDiffEqTsit5
using Logging
using LoggingExtras: FormatLogger, MinLevelLogger, TeeLogger
using LinearAlgebra: BLAS
using InteractiveUtils: InteractiveUtils
using UnicodePlots: UnicodePlots
using PrecompileTools: @setup_workload, @compile_workload

# Package root directory — all relative config paths resolve against this.
# Computed at precompile time: @__DIR__ = src/, dirname = Nbody6Dynamics/.
"""
Root of the package source tree. Used for read-only purposes only — the
shipped `input_files/`, the CUDA helper headers under `deps/cuda/` and the
package's own provenance stamp. Nothing is ever written under it: run and
build directories resolve against the configuration's `config_dir` or an
explicit `base_dir`, so an installation in a read-only depot works as a
checkout does.
"""
const _PACKAGE_ROOT = dirname(@__DIR__)

"""
    example_input(name) -> String

Absolute path of a file shipped under the package's `input_files/` directory
(`example_input("N1k_quick.inp")`, `example_input("showcase/equal_pipeline.toml")`),
for a user who installed the package by URL and has no checkout to point a
configuration at. Raises an `ArgumentError` when no such file ships.
"""
function example_input(name::AbstractString)::String
    path = normpath(joinpath(_PACKAGE_ROOT, "input_files", name))
    isfile(path) || throw(ArgumentError("no shipped input file named \"$name\" under input_files/"))
    return path
end

# ---------------------------------------------------------------------------
# Core types (must come first)
# ---------------------------------------------------------------------------
include("types.jl")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
include("config.jl")

# ---------------------------------------------------------------------------
# Shared utilities (safesave-style backups, …)
# ---------------------------------------------------------------------------
include("util.jl")

# ---------------------------------------------------------------------------
# Platform detection & dependency checking
# ---------------------------------------------------------------------------
include("platform.jl")

# ---------------------------------------------------------------------------
# Install / build pipeline
# ---------------------------------------------------------------------------
include("install.jl")

# ---------------------------------------------------------------------------
# Runtime hardware telemetry (used by the simulation executor)
# ---------------------------------------------------------------------------
include("telemetry.jl")

# ---------------------------------------------------------------------------
# Simulation execution
# ---------------------------------------------------------------------------
include("run.jl")

# ---------------------------------------------------------------------------
# GPU validation of a CUDA host (logged stages under runs/)
# ---------------------------------------------------------------------------
include("gpu_validation.jl")

# ---------------------------------------------------------------------------
# I/O readers
# ---------------------------------------------------------------------------
include("io/io.jl")

# ---------------------------------------------------------------------------
# Per-cluster structure from snapshots (used by the merger plots)
# ---------------------------------------------------------------------------
include("cluster_structure.jl")
include("binary_population.jl")
include("stellar_population.jl")
include("remnant.jl")

# ---------------------------------------------------------------------------
# Parameter sweeps
# ---------------------------------------------------------------------------
include("sweep.jl")
include("ensemble.jl")

# ---------------------------------------------------------------------------
# Figure interface (implemented by the Makie extension, see ext/)
# ---------------------------------------------------------------------------
include("plotting_api.jl")

# ---------------------------------------------------------------------------
# External post-processing (standalone, config-free)
# ---------------------------------------------------------------------------
include("external.jl")

# ---------------------------------------------------------------------------
# Initial conditions generator (merger ICs)
# ---------------------------------------------------------------------------
include("ic/ic.jl")

# ---------------------------------------------------------------------------
# High-level orchestration
# ---------------------------------------------------------------------------

"""
    postprocess(cfg::Nbody6Config; run_dir = "", base_dir = cfg.config_dir) -> Dict{Symbol,Any}

Run all enabled post-processing steps and return collected results.

Data directory resolution (in priority order):
1. `run_dir` keyword — `run_dir/output/` (from a simulation run)
2. `cfg.postprocess.data_dir` — explicit external directory from config.toml,
   resolved against `base_dir` when relative
3. Falls back to the most recent run under `base_dir/<runs_dir>/`
   (via `_find_latest_run`); errors if none exists.

`base_dir` defaults to the configuration's `config_dir` (the directory of the
file it was loaded from).

Besides the readers, the derived data products are computed here: the class
census (`:stellar_census`) and, for a merger run — one whose output holds a
`merger_summary.txt` and at least two snapshots — the remnant diagnostics
(`:remnant`, a [`RemnantDiagnostics`](@ref)).
"""
function postprocess(
    cfg::Nbody6Config;
    run_dir::AbstractString = "",
    base_dir::AbstractString = cfg.config_dir,
)::Dict{Symbol,Any}
    pp = cfg.postprocess

    # Resolve the directory containing simulation output files.
    # run_dir (explicit argument) takes priority over pp.data_dir (config).
    sim_dir = if !isempty(run_dir)
        # run_dir points to runs/<run_id>/ — the data lives in output/
        joinpath(run_dir, "output")
    elseif !isempty(pp.data_dir)
        # Explicit external directory from config
        _resolve_path(base_dir, pp.data_dir)
    else
        latest = _find_latest_run(cfg, base_dir)
        isempty(latest) && error(
            "No run_dir/data_dir given and no matching run found under " *
            joinpath(base_dir, cfg.simulation.runs_dir),
        )
        @info "Post-processing most recent run: $(basename(latest))"
        joinpath(latest, "output")
    end

    isdir(sim_dir) || error("Data directory does not exist: $sim_dir")

    results = Dict{Symbol,Any}()

    # Snapshots (snapshot_format is validated to "conf3" at config load)
    if pp.snapshot_format == "conf3"
        snaps = read_all_conf3(sim_dir, pp.snapshot_pattern)
        isempty(snaps) || (results[:snapshots] = snaps)
    end

    # Diagnostics from stdout
    if pp.parse_stdout
        diag_path = joinpath(sim_dir, pp.stdout_file)
        if isfile(diag_path)
            results[:diagnostics] = read_diagnostics(diag_path)
        else
            @warn "Stdout file not found: $diag_path"
        end
    end

    # Lagrangian radii
    if pp.read_lagr
        lagr_path = joinpath(sim_dir, pp.lagr_file)
        if isfile(lagr_path)
            results[:lagr] = read_lagr(lagr_path)
        else
            @warn "Lagrangian radii file not found: $lagr_path"
        end
    end

    # Escapers
    if pp.read_escapers
        esc_path = joinpath(sim_dir, pp.escapers_file)
        if isfile(esc_path)
            results[:escapers] = read_escapers(esc_path)
        else
            @warn "Escaper file not found: $esc_path"
        end
    end

    # Stellar evolution
    if pp.read_stellar_evo
        sevs = read_all_stellar_evolution(sim_dir, pp.stellar_evo_pattern)
        isempty(sevs) || (results[:stellar_evo] = sevs)
    end

    # Regularised binaries
    if pp.read_binary_evo
        bevs = read_all_binary_evolution(sim_dir, pp.binary_evo_pattern)
        isempty(bevs) || (results[:binary_evo] = bevs)
    end

    # Class census over single stars and binary members
    if haskey(results, :stellar_evo)
        results[:stellar_census] = stellar_census(
            results[:stellar_evo],
            get(results, :binary_evo, BinaryEvolutionSnapshot[]),
        )
    end

    # Remnant diagnostics of a merger run: a data product, so it is computed
    # here and not in the figure layer, and a headless host gets it too.
    summary_path = joinpath(sim_dir, "merger_summary.txt")
    if haskey(results, :snapshots) && length(results[:snapshots]) ≥ 2 && isfile(summary_path)
        ranges = parse_merger_summary(summary_path)
        if !isempty(ranges)
            @info "Remnant diagnostics (bound set, core radius, rotation, segregation)..."
            results[:remnant] = remnant_diagnostics(results[:snapshots], ranges)
        end
    end

    return results
end

"""
    run_pipeline(cfg::Nbody6Config; base_dir = cfg.config_dir, run_id = "") -> Dict{Symbol,Any}

Top-level orchestrator that runs the full pipeline (or any subset) based on
the config flags.  This is the **single entry point** for config-driven workflows.

# Pipeline phases (controlled by config flags)

| Phase          | Controlled by                           | Notes                                   |
|:---------------|:----------------------------------------|:----------------------------------------|
| Install/Build  | `install.enabled`                       | Clone + compile Nbody6++                |
| Simulation     | `simulation.run_test`                   | Run the N-body integration              |
| Post-process   | `postprocess.enabled`                   | Read output data + sanity checks        |
| Plots          | `visualization.enabled`                 | Generate static plots + GIF animations  |

# Post-process–only mode

To run **only** post-processing and plotting on existing output (from this
project or an external source), set in `config.toml`:

```toml
[simulation]
run_test = false          # skip simulation

[postprocess]
enabled   = true
data_dir  = "/path/to/output"   # point at any dir with Nbody6++ files
```

Then call:
```julia
cfg = load_config("config.toml")
results = run_pipeline(cfg)
```

If `data_dir` is empty but `run_test = false`, the most recent run in
`runs/` is post-processed.

`run_id` fixes the run directory name (`<runs_dir>/<run_id>`); when empty,
`generate_run_id` derives it from `simulation.run_id_prefix` and the time.

`base_dir` is the project directory: every relative path of the configuration
(`install.install_dir`, `simulation.input_file`, `simulation.runs_dir`,
`postprocess.data_dir`, `merger.config_file`) resolves against it. It defaults
to `cfg.config_dir`, the directory of the file `load_config` read, so
`julia scripts/run_setup.jl path/to/config.toml` keeps its engine and its
runs next to that file, never inside the package.

# Returns
A `Dict{Symbol,Any}` with keys `:snapshots`, `:diagnostics`, `:lagr`,
`:escapers`, `:stellar_evo`, `:binary_evo`, `:stellar_census`, `:remnant`
(present only when corresponding data exists).
Returns an empty dict if post-processing is disabled.
"""
function run_pipeline(
    cfg::Nbody6Config;
    base_dir::AbstractString = cfg.config_dir,
    run_id::AbstractString = "",
)::Dict{Symbol,Any}
    # A run that must end in figures is refused before the build and the
    # integration, not after them: the backend is a property of the session,
    # and discovering it missing at phase 4 would cost the whole run.
    cfg.visualization.enabled && _require_plotting(
        :run_pipeline,
        "Set `[visualization] enabled = false` to run the pipeline without figures.",
    )

    # Phases actually performed, stamped into RUN_INFO.toml at the end so an
    # interrupted pipeline is distinguishable from a finished one.
    t_pipeline = time()
    phases = String[]

    # ── Phase 1: Install / Build ──
    if cfg.install.enabled
        @info "Phase 1: Installing Nbody6++..."
        setup_nbody6(cfg; base_dir = base_dir)
        push!(phases, "install")
    end

    # ── Phase 1.5: Merger IC Generation ──
    merger_result = nothing
    run_dir = ""
    if cfg.merger.enabled
        @info "Phase 1.5: Generating merger initial conditions..."
        merger_cfg_path = _resolve_path(base_dir, cfg.merger.config_file)
        isempty(merger_cfg_path) && error("merger.enabled=true but merger.config_file is empty")
        isfile(merger_cfg_path) || error("Merger config not found: $merger_cfg_path")

        merger_cfg = load_merger_config(merger_cfg_path)

        # Create a run directory for this merger. generate_run_id adds the
        # timestamp + a 4-hex uniqueness suffix, so two merger pipelines
        # started in the same second cannot collide; the configured prefix
        # is preserved so _find_latest_run can locate merger runs.
        isempty(run_id) && (run_id = generate_run_id("merger_" * cfg.simulation.run_id_prefix))
        run_dir = joinpath(base_dir, cfg.simulation.runs_dir, run_id)
        ic_dir = joinpath(run_dir, "output")
        mkpath(ic_dir)

        merger_result = generate_merger_ic(merger_cfg; output_dir = ic_dir)
        @info "  Merger ICs written to: $ic_dir"
        push!(phases, "merger_ic")
    end

    # ── Phase 2: Simulation ──
    if cfg.simulation.run_test
        if merger_result !== nothing
            # Merger + simulation: run inside the merger's output directory
            # so the binary finds dat.10 in its working directory
            @info "Phase 2: Running simulation with merger ICs..."
            run_dir = _run_merger_simulation(cfg, merger_result; base_dir = base_dir)
        else
            @info "Phase 2: Running simulation..."
            run_dir = run_simulation(cfg; base_dir = base_dir, run_id = run_id)
        end
        push!(phases, "simulation")
    else
        @info "Phase 2: Simulation skipped (run_test = false)"
        if isempty(run_dir) && isempty(cfg.postprocess.data_dir)
            run_dir = _find_latest_run(cfg, base_dir)
        end
    end

    # ── Phase 3: Post-processing ──
    results = Dict{Symbol,Any}()
    if cfg.postprocess.enabled
        if !isempty(cfg.postprocess.data_dir)
            @info "Phase 3: Post-processing external data: $(cfg.postprocess.data_dir)"
            results = postprocess(cfg; base_dir = base_dir)
            push!(phases, "postprocess")
        elseif !isempty(run_dir)
            @info "Phase 3: Post-processing run: $(basename(run_dir))"
            results = postprocess(cfg; run_dir = run_dir, base_dir = base_dir)
            push!(phases, "postprocess")
        else
            @warn "Phase 3: No data to post-process (no run_dir and no data_dir)"
        end
    end

    # Thread merger IC result into results for plotting
    if merger_result !== nothing
        results[:merger_ic] = merger_result
    end

    # Directory that receives the derived products of this data set
    products_dir = if !isempty(cfg.postprocess.data_dir)
        # normpath strips any trailing slash so dirname yields the parent
        dirname(_resolve_path(base_dir, cfg.postprocess.data_dir))
    else
        run_dir
    end

    # Data products are written here, not by the figure layer: a host
    # without a plotting backend must get them too.
    if haskey(results, :stellar_census) && !isempty(products_dir)
        write_stellar_census(joinpath(products_dir, "stellar_census.csv"), results[:stellar_census])
    end
    if haskey(results, :remnant) && !isempty(products_dir)
        write_remnant_diagnostics(
            joinpath(products_dir, "remnant_diagnostics.csv"),
            results[:remnant],
        )
    end

    # ── Phase 4: Plots ──
    if cfg.visualization.enabled && !isempty(results)
        @info "Phase 4: Generating plots..."
        generate_plots(results, cfg; run_dir = products_dir)
        push!(phases, "plots")
    end

    _stamp_pipeline_completion(run_dir, phases, time() - t_pipeline)
    return results
end

"""
Run Nbody6++ inside the merger IC output directory so `dat.10` is found in
the working directory: thin wrapper around [`_execute_simulation`](@ref)
pointing at `merger.inp` instead of the config's input file.
"""
function _run_merger_simulation(
    cfg::Nbody6Config,
    merger_result::MergerICResult;
    base_dir::AbstractString = cfg.config_dir,
)::String
    # The run directory is the merger's parent (merger writes to run_dir/output/).
    # Absolute paths required — the launch script runs from out_dir via cd().
    run_dir = abspath(dirname(merger_result.output_dir))
    out_dir = abspath(merger_result.output_dir)  # contains dat.10 and merger.inp
    input_path = joinpath(out_dir, "merger.inp")

    @info "Run dir: $run_dir"
    @info "Input:   $input_path"

    return _execute_simulation(
        cfg,
        run_dir,
        out_dir,
        input_path;
        base_dir = base_dir,
        label = "merger simulation",
    )
end

"""Find the most recent run directory under `base_dir/runs_dir/` (`runs_dir` resolved against `base_dir`).

Matches both plain runs (`<prefix>_*`) and merger runs
(`merger_<prefix>_*`), which share the timestamp-based naming from
`generate_run_id` so a lexicographic sort on the timestamp part yields
the most recent run.
"""
function _find_latest_run(cfg::Nbody6Config, base_dir::AbstractString)::String
    runs_base = _resolve_path(base_dir, cfg.simulation.runs_dir)
    isdir(runs_base) || return ""
    prefix = cfg.simulation.run_id_prefix * "_"
    merger_prefix = "merger_" * prefix
    dirs = filter(readdir(runs_base; join = true)) do p
        isdir(p) && (startswith(basename(p), prefix) || startswith(basename(p), merger_prefix))
    end
    isempty(dirs) && return ""
    # Sort by the timestamp part of the run ID and return the latest
    timestamp_key(p) = begin
        b = basename(p)
        startswith(b, merger_prefix) ? b[(length("merger_") + 1):end] : b
    end
    return sort(dirs; by = timestamp_key)[end]
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
export Nbody6Config,
    InstallConfig,
    BuildConfig,
    SimulationConfig,
    PostprocessConfig,
    VisualizationConfig,
    MergerPipelineConfig
export Snapshot, SnapshotHeader, DiagnosticsData, AdjustRecord, LagrangianData, UnitScaling
export EscaperRecord, StellarRecord, StellarEvolutionSnapshot, STELLAR_TYPE_LABELS
export load_config,
    save_config, setup_nbody6, run_simulation, postprocess, generate_plots, run_pipeline
export scan_output, postprocess_external, OutputScan
export read_conf3, read_all_conf3
export read_diagnostics, extract_scaling
export read_lagr
export read_escapers
export read_stellar_evolution, read_all_stellar_evolution
export BinaryRecord, BinaryEvolutionSnapshot, read_binary_evolution, read_all_binary_evolution
export BinaryPopulation, binary_population, binary_hardness, hardness_scale, binary_scales
export StellarClass, STELLAR_CLASSES, stellar_class, stellar_class_index
export HRPopulation, hr_population, hr_populations
export StellarCensus, stellar_census, class_counts, classes_present, write_stellar_census
export semi_major_axis_pc, binding_energy
export plot_binary_population, plot_binary_orbital_elements, plot_binary_period_distribution
export TelemetrySample, read_telemetry, read_run_telemetry, plot_telemetry
export run_gpu_validation
export SweepConfig, SweepPoint, load_sweep_config, sweep_points, prepare_sweep, run_sweep
export run_sweep_point, read_sweep_index, write_sweep_index, sweep_summary, write_sweep_summary
export sweep_visualization, plot_sweep_lagrangian, plot_sweep_energy, sweep_figures
export EnsembleStatistics, ensemble_statistics, sweep_ensembles, plot_sweep_ensemble
export RemnantDiagnostics, RotationProfile, remnant_diagnostics, write_remnant_diagnostics
export coalescence_time, core_radius, rotation_analysis, mass_segregation
export plot_remnant_rotation, plot_rotation_profile, plot_remnant_structure
export plot_mass_segregation_evolution, remnant_figures
export control_merger_dict, write_control_merger_config, plot_control_comparison
export plot_snapshot, plot_snapshot_evolution
export plot_lagrangian, plot_energy, plot_particle_count
export plot_hr, plot_hr_evolution
export plot_escapers, plot_escape_anisotropy
export plot_mass_segregation, plot_evolutionary_clock, plot_core_mass
export plot_cluster_separation, plot_cluster_virial, per_cluster_virial, parse_merger_summary
export ClusterStructure, cluster_structure, plot_cluster_structure, bound_fraction
export RadialProfile, radial_profile, cluster_profiles, system_profile, model_density
export plot_density_profiles, plot_velocity_dispersion
export animate_cluster, animate_hr, animate_lagrangian
export set_publication_theme!, publication_theme, plotting_available
export generate_run_id, restart_simulation, export_for_paper, example_input
export nparticles, time_nb, time_myr, rbar, zmbar, tscale, vstar, rscale, rc
export detect_platform, check_dependencies, detect_cuda_path
export engine_interval
export detect_compute_capabilities, cuda_arch_from_compute_cap, cuda_gencode_flags
export resolve_cuda_arch, nvcc_release, nvcc_supported_archs
export ClusterSpec,
    BinarySpec,
    OrbitSpec,
    MergerOutputSpec,
    Nbody6ParameterSpec,
    StellarSpec,
    TidalSpec,
    MergerConfig,
    MergerICResult
export DensityProfile,
    KingProfile,
    PlummerProfile,
    IMFSpec,
    KroupaIMF,
    RescaledKroupaIMF,
    EqualMassIMF,
    profile_name,
    imf_name,
    expected_mass,
    sample_masses,
    kroupa_mean_mass
export load_merger_config, generate_merger_ic, run_merger_pipeline, load_merger_ic_result
export plot_merger_ic
export sample_plummer, sample_king, sample_kroupa
export virialise!, kepler_velocity, jacobi_radius
export write_dat10, generate_merger_inp, resolve_nbody6_parameters, crossing_time, to_nbody_units!
export sample_binaries, expand_binaries

# ---------------------------------------------------------------------------
# Precompile workload: the configuration and I/O paths every session
# hits first, plus a small merger initial-condition generation. Plotting is
# left out; the CairoMakie precompile is its own. Runs at package
# precompilation only, silently.
# ---------------------------------------------------------------------------
@setup_workload begin
    _pc_dir = mktempdir()
    _pc_cfg = joinpath(_pc_dir, "config.toml")
    write(_pc_cfg, "[simulation]\nomp_threads = 2\n\n[merger]\nenabled = false\n")
    _pc_out = joinpath(_pc_dir, "out1000")
    write(
        _pc_out,
        " ADJUST: TIME  0.0  T[MYR]  0.0  Q  0.5  DE  0.0  E  -0.25\n" *
        " ADJUST: TIME  1.0  T[MYR]  0.5  Q  0.5  DE  1e-8  E  -0.25\n",
    )
    _pc_merger = joinpath(_pc_dir, "merger.toml")
    write(
        _pc_merger,
        """
        [merger]
        n_clusters = 2
        orbit_mode = "kepler"
        seed = 1

        [merger.cluster1]
        model = "plummer"
        N = 40
        rbar = 1.0
        imf = "kroupa"

        [merger.cluster2]
        model = "plummer"
        N = 40
        rbar = 1.0
        imf = "kroupa"

        [merger.orbit]
        apocentre = 5.0
        eccentricity = 0.5

        [merger.output]
        format = "nbody"
        truncate_jacobi = true
        output_dir = "$(_pc_dir)"
        tcrit = 1.0
        dtadj = 0.25
        deltat = 0.5
        """,
    )
    @compile_workload begin
        Logging.with_logger(Logging.NullLogger()) do
            cfg = load_config(_pc_cfg)
            save_config(cfg, joinpath(_pc_dir, "frozen.toml"))
            read_diagnostics(_pc_out)
            engine_interval(0.63)
            generate_run_id("pc")
            _format_elapsed(125.0)
            generate_merger_ic(load_merger_config(_pc_merger))
        end
    end
    rm(_pc_dir; recursive = true, force = true)
end

end # module
