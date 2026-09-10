module Nbody6Dynamics

using TOML
using Dates
using Printf
using CairoMakie
using LaTeXStrings
using ProgressMeter
using Random
using SpecialFunctions
using OrdinaryDiffEqTsit5
using Logging
using LoggingExtras: FormatLogger, MinLevelLogger, TeeLogger
using LinearAlgebra: BLAS
using MathTeXEngine: texfont

# Package root directory — all relative config paths resolve against this.
# Computed at precompile time: @__DIR__ = src/, dirname = Nbody6Dynamics/.
const _PROJECT_ROOT = dirname(@__DIR__)

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
# I/O readers
# ---------------------------------------------------------------------------
include("io/io.jl")

# ---------------------------------------------------------------------------
# Per-cluster structure from snapshots (used by the merger plots)
# ---------------------------------------------------------------------------
include("cluster_structure.jl")
include("binary_population.jl")
include("remnant.jl")

# ---------------------------------------------------------------------------
# Parameter sweeps
# ---------------------------------------------------------------------------
include("sweep.jl")
include("ensemble.jl")

# ---------------------------------------------------------------------------
# Plotting (sets publication theme on load)
# ---------------------------------------------------------------------------
include("plotting/plotting.jl")

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
    postprocess(cfg::Nbody6Config; run_dir = "", base_dir = _PROJECT_ROOT) -> Dict{Symbol,Any}

Run all enabled post-processing steps and return collected results.

Data directory resolution (in priority order):
1. `run_dir` keyword — `run_dir/output/` (from a simulation run)
2. `cfg.postprocess.data_dir` — explicit external directory from config.toml
3. Falls back to the most recent run under `base_dir/runs/`
   (via `_find_latest_run`); errors if none exists.

`base_dir` defaults to the package root directory, making the pipeline path-agnostic.
"""
function postprocess(
    cfg::Nbody6Config;
    run_dir::AbstractString = "",
    base_dir::AbstractString = _PROJECT_ROOT,
)::Dict{Symbol,Any}
    pp = cfg.postprocess

    # Resolve the directory containing simulation output files.
    # run_dir (explicit argument) takes priority over pp.data_dir (config).
    sim_dir = if !isempty(run_dir)
        # run_dir points to runs/<run_id>/ — the data lives in output/
        joinpath(run_dir, "output")
    elseif !isempty(pp.data_dir)
        # Explicit external directory from config
        abspath(pp.data_dir)
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

    return results
end

"""
    _plot_binary_diagnostics(bevs, results, vis)

Binary-population figures for `bevs`. The hard/soft energy scale and the
stellar count per epoch come from the conf.3 snapshots when present
([`binary_scales`](@ref)); without snapshots the stellar count falls back
to the ADJUST diagnostics and the pairs are left unclassified.
"""
function _plot_binary_diagnostics(
    bevs::Vector{BinaryEvolutionSnapshot},
    results::Dict{Symbol,Any},
    vis::VisualizationConfig,
)
    snaps = get(results, :snapshots, Snapshot[])::Vector{Snapshot}
    scales = binary_scales(bevs, snaps)
    if scales === nothing
        n_stars = nothing
        if haskey(results, :diagnostics)
            adj = (results[:diagnostics]::DiagnosticsData).adjust
            if !isempty(adj)
                t_adj = [a.time_myr for a in adj]
                n_stars = [adj[argmin(abs.(t_adj .- b.time_myr))].n for b in bevs]
            end
        end
        pop = binary_population(bevs; n_stars = n_stars)
        m_mean = fill(NaN, length(bevs))
        sigma = fill(NaN, length(bevs))
    else
        pop = binary_population(
            bevs;
            n_stars = scales.n_stars,
            m_mean = scales.m_mean,
            sigma_kms = scales.sigma_kms,
        )
        m_mean = scales.m_mean
        sigma = scales.sigma_kms
    end
    plot_binary_population(pop, vis; filename = "binary_population")
    plot_binary_period_distribution(bevs, vis; filename = "binary_period_distribution")
    epochs =
        length(bevs) > 1 ?
        ((1, "binary_orbital_elements_initial"), (length(bevs), "binary_orbital_elements_final")) :
        ((1, "binary_orbital_elements_initial"),)
    for (idx, fname) in epochs
        plot_binary_orbital_elements(
            bevs[idx],
            vis;
            m_mean = m_mean[idx],
            sigma_kms = sigma[idx],
            filename = fname,
        )
    end
    return nothing
end

"""
    generate_plots(results::Dict{Symbol,Any}, vis::VisualizationConfig;
                   sim_dir::AbstractString = "", animations::Bool = true)

Generate all available plots (and, when `animations = true`, GIF
animations) from post-processing results, saving into `vis.output_dir`.
`sim_dir` names the directory holding the raw simulation output; when given,
it is searched for `merger_summary.txt` to produce the merger-specific
figures (inter-cluster separation, per-cluster virial ratio).

This is the single plot dispatcher — both the config-driven pipeline
(via the `Nbody6Config` method) and [`postprocess_external`](@ref) route
through it.
"""
function generate_plots(
    results::Dict{Symbol,Any},
    vis::VisualizationConfig;
    sim_dir::AbstractString = "",
    animations::Bool = true,
)
    set_publication_theme!()

    if haskey(results, :snapshots)
        snaps = results[:snapshots]::Vector{Snapshot}
        if !isempty(snaps)
            @info "Plotting final snapshot..."
            plot_snapshot(snaps[end], vis; filename = "snapshot_final")
            if length(snaps) > 1
                @info "Plotting snapshot evolution..."
                plot_snapshot_evolution(snaps, vis; filename = "snapshot_evolution")
            end

            # Merger-specific: inter-cluster separation and per-cluster virial.
            # merger_summary.txt lives in the simulation output directory.
            if !isempty(sim_dir) && isdir(sim_dir)
                summary_path = joinpath(sim_dir, "merger_summary.txt")
                if isfile(summary_path) && length(snaps) ≥ 2
                    ranges = parse_merger_summary(summary_path)
                    if !isempty(ranges)
                        if length(ranges) ≥ 2
                            @info "Plotting inter-cluster separation..."
                            plot_cluster_separation(snaps, ranges, vis)
                        end
                        @info "Plotting per-cluster virial ratio..."
                        plot_cluster_virial(snaps, ranges, vis)
                        @info "Plotting per-cluster structure..."
                        plot_cluster_structure(
                            snaps,
                            ranges,
                            vis;
                            lagr = get(results, :lagr, nothing),
                        )
                        # Radial profiles against the generating models (merger_ic.toml)
                        ic_meta = joinpath(sim_dir, "merger_ic.toml")
                        specs =
                            isfile(ic_meta) ? load_merger_ic_result(sim_dir).cluster_specs : nothing
                        specs === nothing || length(specs) == length(ranges) || (specs = nothing)
                        @info "Plotting density profiles (initial and final snapshots)..."
                        plot_density_profiles(
                            snaps[1],
                            ranges,
                            vis;
                            specs = specs,
                            filename = "merger_density_profiles_initial",
                        )
                        plot_density_profiles(
                            snaps[end],
                            ranges,
                            vis;
                            specs = specs,
                            filename = "merger_density_profiles_final",
                        )
                        @info "Plotting velocity dispersion profiles (final snapshot)..."
                        plot_velocity_dispersion(snaps[end], ranges, vis)
                        @info "Remnant diagnostics (bound set, core radius, rotation, segregation)..."
                        diag = remnant_diagnostics(snaps, ranges)
                        write_remnant_diagnostics(
                            joinpath(dirname(abspath(sim_dir)), "remnant_diagnostics.csv"),
                            diag,
                        )
                        remnant_figures(diag, vis)
                    end
                end
            end
        end
    end

    if haskey(results, :diagnostics)
        diag = results[:diagnostics]::DiagnosticsData
        if !isempty(diag.adjust)
            @info "Plotting energy diagnostics..."
            plot_energy(diag, vis; filename = "energy")
            plot_particle_count(diag, vis; filename = "particle_count")
        end
    end

    # Unit scaling for readers whose files carry no header (lagr.7):
    # derived from the diagnostics when available, else NB units.
    scaling =
        haskey(results, :diagnostics) ? extract_scaling(results[:diagnostics]::DiagnosticsData) :
        nothing

    if haskey(results, :lagr)
        lagr = results[:lagr]::LagrangianData
        if !isempty(lagr.time)
            @info "Plotting Lagrangian radii..."
            plot_lagrangian(lagr, vis; filename = "lagrangian_radii", units = scaling)
        end
    end

    if haskey(results, :escapers)
        escs = results[:escapers]::Vector{EscaperRecord}
        if !isempty(escs)
            @info "Plotting escaper analysis..."
            plot_escapers(escs, vis; filename = "escapers")
            plot_escape_anisotropy(escs, vis; filename = "escape_anisotropy")
        end
    end

    if haskey(results, :stellar_evo)
        sevs = results[:stellar_evo]::Vector{StellarEvolutionSnapshot}
        if !isempty(sevs)
            # Three HR diagrams: beginning, middle, end
            mid = max(1, length(sevs) ÷ 2)
            hr_epochs = [
                (1, "hr_diagram_early"),
                (mid, "hr_diagram_mid"),
                (length(sevs), "hr_diagram_final"),
            ]
            for (idx, fname) in hr_epochs
                @info "Plotting HR diagram (epoch $idx/$(length(sevs)))..."
                plot_hr(sevs[idx], vis; filename = fname)
            end
            if length(sevs) > 1
                @info "Plotting HR evolution..."
                plot_hr_evolution(sevs, vis; filename = "hr_evolution")
            end
            @info "Plotting SSE quantities..."
            plot_mass_segregation(sevs[end], vis; filename = "mass_segregation")
            plot_evolutionary_clock(sevs[end], vis; filename = "evolutionary_clock")
            plot_core_mass(sevs, vis; filename = "core_mass_growth")
        end
    end

    if haskey(results, :binary_evo)
        bevs = results[:binary_evo]::Vector{BinaryEvolutionSnapshot}
        if !isempty(bevs)
            @info "Plotting binary population..."
            _plot_binary_diagnostics(bevs, results, vis)
        end
    end

    # --- Animations (GIF) ---
    if animations
        if haskey(results, :snapshots)
            snaps = results[:snapshots]::Vector{Snapshot}
            if length(snaps) > 1
                @info "Animating cluster evolution..."
                animate_cluster(snaps, vis; filename = "cluster_evolution")
            end
        end

        if haskey(results, :lagr)
            lagr = results[:lagr]::LagrangianData
            if length(lagr.time) > 1
                @info "Animating Lagrangian radii..."
                animate_lagrangian(lagr, vis; filename = "lagrangian_anim", units = scaling)
            end
        end

        if haskey(results, :stellar_evo)
            sevs = results[:stellar_evo]::Vector{StellarEvolutionSnapshot}
            if length(sevs) > 1
                @info "Animating HR diagram evolution..."
                animate_hr(sevs, vis; filename = "hr_evolution_anim")
            end
        end
    end

    # --- Merger IC diagnostic plots ---
    if haskey(results, :merger_ic)
        merger_res = results[:merger_ic]::MergerICResult
        @info "Plotting merger IC diagnostics..."
        plot_merger_ic(merger_res, vis)
    end

    @info "All plots and animations saved to: $(vis.output_dir)"
    return nothing
end

"""
    generate_plots(results::Dict{Symbol,Any}, cfg::Nbody6Config;
                   run_dir::AbstractString = "")

Config-driven wrapper around the `VisualizationConfig` method.

When `run_dir` is provided (e.g. `runs/run_XXXX/`), plots are saved to
`run_dir/<visualization.output_dir>/` (typically `runs/run_XXXX/plots/`) and
`run_dir/output/` is searched for merger metadata. Otherwise, falls back to
`visualization.output_dir` relative to the package root.
"""
function generate_plots(results::Dict{Symbol,Any}, cfg::Nbody6Config; run_dir::AbstractString = "")
    # Build a VisualizationConfig with the output_dir resolved to the run
    vis = if !isempty(run_dir)
        plots_dir = joinpath(run_dir, cfg.visualization.output_dir)
        # Forward EVERY field except output_dir — a missed field here means
        # the user's [visualization] settings are silently dropped on the
        # main pipeline path (this has happened twice; see git history).
        VisualizationConfig(;
            enabled = cfg.visualization.enabled,
            format = cfg.visualization.format,
            dpi = cfg.visualization.dpi,
            column = cfg.visualization.column,
            figsize = cfg.visualization.figsize,
            units = cfg.visualization.units,
            output_dir = plots_dir,
            style = cfg.visualization.style,
        )
    else
        cfg.visualization
    end
    sim_dir = !isempty(run_dir) ? joinpath(run_dir, "output") : ""
    return generate_plots(results, vis; sim_dir = sim_dir)
end

"""
    run_pipeline(cfg::Nbody6Config; base_dir = _PROJECT_ROOT, run_id = "") -> Dict{Symbol,Any}

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

# Returns
A `Dict{Symbol,Any}` with keys `:snapshots`, `:diagnostics`, `:lagr`,
`:escapers`, `:stellar_evo`, `:binary_evo` (present only when corresponding data exists).
Returns an empty dict if post-processing is disabled.
"""
function run_pipeline(
    cfg::Nbody6Config;
    base_dir::AbstractString = _PROJECT_ROOT,
    run_id::AbstractString = "",
)::Dict{Symbol,Any}
    # ── Phase 1: Install / Build ──
    if cfg.install.enabled
        @info "Phase 1: Installing Nbody6++..."
        setup_nbody6(cfg; base_dir = base_dir)
    end

    # ── Phase 1.5: Merger IC Generation ──
    merger_result = nothing
    run_dir = ""
    if cfg.merger.enabled
        @info "Phase 1.5: Generating merger initial conditions..."
        merger_cfg_path = cfg.merger.config_file
        isempty(merger_cfg_path) && error("merger.enabled=true but merger.config_file is empty")
        # Resolve relative paths against project root
        if !isabspath(merger_cfg_path)
            merger_cfg_path = joinpath(base_dir, merger_cfg_path)
        end
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
            results = postprocess(cfg)
        elseif !isempty(run_dir)
            @info "Phase 3: Post-processing run: $(basename(run_dir))"
            results = postprocess(cfg; run_dir = run_dir)
        else
            @warn "Phase 3: No data to post-process (no run_dir and no data_dir)"
        end
    end

    # Thread merger IC result into results for plotting
    if merger_result !== nothing
        results[:merger_ic] = merger_result
    end

    # ── Phase 4: Plots ──
    if cfg.visualization.enabled && !isempty(results)
        @info "Phase 4: Generating plots..."
        plot_run_dir = if !isempty(cfg.postprocess.data_dir)
            # normpath strips any trailing slash so dirname yields the parent
            dirname(abspath(normpath(cfg.postprocess.data_dir)))
        elseif !isempty(run_dir)
            run_dir
        else
            ""
        end
        generate_plots(results, cfg; run_dir = plot_run_dir)
    end

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
    base_dir::AbstractString = _PROJECT_ROOT,
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

"""Find the most recent run directory under `base_dir/runs_dir/`.

Matches both plain runs (`<prefix>_*`) and merger runs
(`merger_<prefix>_*`), which share the timestamp-based naming from
`generate_run_id` so a lexicographic sort on the timestamp part yields
the most recent run.
"""
function _find_latest_run(cfg::Nbody6Config, base_dir::AbstractString)::String
    runs_base = joinpath(base_dir, cfg.simulation.runs_dir)
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
export semi_major_axis_pc, binding_energy
export plot_binary_population, plot_binary_orbital_elements, plot_binary_period_distribution
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
export set_publication_theme!
export generate_run_id, restart_simulation, export_for_paper
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

end # module
