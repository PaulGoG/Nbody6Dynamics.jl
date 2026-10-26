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
using KernelAbstractions: KernelAbstractions, @kernel, @index, @Const, CPU
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

Absolute path of a file shipped under the package's `input_files/` directory,
for a user who installed the package by URL and has no checkout to point a
configuration at. `name` is a path below `input_files/`
(`example_input("engine/N1k_quick.inp")`, `example_input("showcase/equal_pipeline.toml")`);
a bare file name is looked up across the sub-directories
(`example_input("N1k_quick.inp")`), which keeps the names of the flat layout
before 0.4.0 valid. Raises an `ArgumentError` when no such file ships, or
when a bare name matches in more than one place.
"""
function example_input(name::AbstractString)::String
    root = joinpath(_PACKAGE_ROOT, "input_files")
    path = normpath(joinpath(root, name))
    isfile(path) && return path
    if !occursin('/', name) && !occursin('\\', name)
        hits = String[]
        for (dir, _, files) in walkdir(root)
            name in files && push!(hits, normpath(joinpath(dir, name)))
        end
        length(hits) == 1 && return first(hits)
        length(hits) > 1 && throw(
            ArgumentError(
                "\"$name\" ships in more than one place under input_files/ (" *
                join(relpath.(hits, root), ", ") *
                "); give the sub-directory",
            ),
        )
    end
    throw(ArgumentError("no shipped input file named \"$name\" under input_files/"))
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
include("pipeline/platform.jl")

# ---------------------------------------------------------------------------
# Install / build pipeline
# ---------------------------------------------------------------------------
include("pipeline/install.jl")

# ---------------------------------------------------------------------------
# Runtime hardware telemetry (used by the simulation executor)
# ---------------------------------------------------------------------------
include("pipeline/telemetry.jl")

# ---------------------------------------------------------------------------
# Simulation execution
# ---------------------------------------------------------------------------
include("pipeline/run.jl")
include("pipeline/resume.jl")
include("pipeline/cost.jl")
include("pipeline/cluster.jl")

# ---------------------------------------------------------------------------
# GPU validation of a CUDA host (logged stages under runs/)
# ---------------------------------------------------------------------------
include("pipeline/gpu_validation.jl")

# ---------------------------------------------------------------------------
# I/O readers
# ---------------------------------------------------------------------------
include("io/io.jl")

# ---------------------------------------------------------------------------
# Potential, bound sets and integrals of snapshots
# ---------------------------------------------------------------------------
include("analysis/potential.jl")
include("analysis/integrals.jl")
include("analysis/energy_events.jl")
include("analysis/budgets.jl")

# ---------------------------------------------------------------------------
# Per-cluster structure from snapshots (used by the merger plots)
# ---------------------------------------------------------------------------
include("analysis/cluster_structure.jl")
include("analysis/binary_population.jl")
include("analysis/stellar_population.jl")
include("analysis/remnant.jl")
include("analysis/shape_mixing.jl")
include("analysis/pair_finder.jl")
include("analysis/equipartition.jl")
include("analysis/angular_momentum.jl")
include("analysis/merger_tree.jl")
include("analysis/validity.jl")
include("analysis/binary_spans.jl")
include("analysis/bound_spans.jl")
include("analysis/binary_census.jl")
include("analysis/gw_observables.jl")
include("analysis/binary_population_series.jl")
include("analysis/gw_catalogues.jl")
include("analysis/collision_catalogue.jl")
include("analysis/survival.jl")

# ---------------------------------------------------------------------------
# Parameter sweeps
# ---------------------------------------------------------------------------
include("pipeline/sweep.jl")
include("analysis/ensemble.jl")

# ---------------------------------------------------------------------------
# Figure interface (implemented by the Makie extension, see ext/)
# ---------------------------------------------------------------------------
include("plotting_api.jl")

# ---------------------------------------------------------------------------
# External post-processing (standalone, config-free)
# ---------------------------------------------------------------------------
include("pipeline/external.jl")

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
census (`:stellar_census`), the integrals of motion of every snapshot
(`:integrals`, [`SnapshotIntegrals`](@ref), with `postprocess.snapshot_integrals`,
including the Jacobi integral when the stdout capture echoes a linearised tidal
field) and, for a merger run — one whose output holds a `merger_summary.txt` and
at least two snapshots — the remnant diagnostics (`:remnant`, a
[`RemnantDiagnostics`](@ref)), which reuse the integrals' bound sets. Both
need an O(N²) pair sum per snapshot and run only while no snapshot holds more
than `postprocess.pair_sum_max_n` stars. From the stdout capture the stamped
event lines (`:engine_events`, [`EngineEvent`](@ref)) and the energy-jump
localisation (`:energy_events`, [`EnergyEventAnalysis`](@ref), with
`postprocess.energy_events`) follow; with the integrals, the escapers and the
counters of `postprocess.events_file` (`:event_counters`, [`EventCounters`](@ref))
the mass budget (`:mass_budget`, [`MassBudgetRow`](@ref), with
`postprocess.mass_budget`). Under a linearised tidal field ([`tidal_mode`](@ref)
`:linearised`) the engine's `DE` is not an error measure and no interval is
flagged.
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
    # Wall time of each stage, in order (`results[:timings]`, RUN_INFO [pipeline]).
    timings = Pair{String,Float64}[]
    t_lap = Ref(time())
    lap!(label) = (t = time(); push!(timings, label => t - t_lap[]); t_lap[] = t; nothing)

    # Snapshots (snapshot_format is validated to "conf3" at config load)
    if pp.snapshot_format == "conf3"
        snaps = read_all_conf3(sim_dir, pp.snapshot_pattern)
        isempty(snaps) || (results[:snapshots] = snaps)
    end
    lap!("snapshots")

    # Diagnostics from stdout, and its stamped event lines
    if pp.parse_stdout
        diag_path = joinpath(sim_dir, pp.stdout_file)
        if isfile(diag_path)
            results[:diagnostics] = read_diagnostics(diag_path)
            if pp.energy_events || pp.mass_budget
                results[:engine_events] = read_engine_events(diag_path)
            end
        else
            @warn "Stdout file not found: $diag_path"
        end
    end
    lap!("stdout")

    # Engine counters per output time (written with stellar evolution on)
    if pp.mass_budget
        events_path = joinpath(sim_dir, pp.events_file)
        if isfile(events_path)
            results[:event_counters] = read_event_counters(events_path)
        else
            @info "No counters file $(pp.events_file): stellar mass loss taken as zero in the mass budget"
        end
    end
    lap!("event_counters")

    # Lagrangian radii
    if pp.read_lagr
        lagr_path = joinpath(sim_dir, pp.lagr_file)
        if isfile(lagr_path)
            results[:lagr] = read_lagr(lagr_path)
        else
            @warn "Lagrangian radii file not found: $lagr_path"
        end
    end
    lap!("lagr")

    # Escapers
    if pp.read_escapers
        esc_path = joinpath(sim_dir, pp.escapers_file)
        if isfile(esc_path)
            results[:escapers] = read_escapers(esc_path)
        else
            @warn "Escaper file not found: $esc_path"
        end
    end
    lap!("escapers")

    # Stellar evolution
    if pp.read_stellar_evo
        sevs = read_all_stellar_evolution(sim_dir, pp.stellar_evo_pattern)
        isempty(sevs) || (results[:stellar_evo] = sevs)
    end
    lap!("stellar_evolution")

    # Regularised binaries
    if pp.read_binary_evo
        bevs = read_all_binary_evolution(sim_dir, pp.binary_evo_pattern)
        isempty(bevs) || (results[:binary_evo] = bevs)
    end
    lap!("binary_evolution")

    # Class census over single stars and binary members
    if haskey(results, :stellar_evo)
        results[:stellar_census] = stellar_census(
            results[:stellar_evo],
            get(results, :binary_evo, BinaryEvolutionSnapshot[]),
        )
    end
    lap!("census")

    # Integrals of motion and remnant diagnostics: data products, so they are
    # computed here and not in the figure layer, and a headless host gets them
    # too. Each needs an O(N²) pair sum per snapshot (the bound set), so both
    # are subject to the configured limit; the remnant diagnostics reuse the
    # integrals' bound sets.
    bound = (; max_passes = pp.bound_max_passes, min_members = pp.bound_min_members)
    # External field: the engine never adds the potential energy of a
    # linearised tide to its total, so the Jacobi integral of the snapshots
    # is the conserved quantity there and the engine's DE is not an error.
    mode = haskey(results, :diagnostics) ? tidal_mode(results[:diagnostics]) : :none
    tidal = haskey(results, :diagnostics) ? results[:diagnostics].tidal : (0.0, 0.0, 0.0, 0.0)
    tidal_energy_amended = _tidal_energy_amended(sim_dir)
    de_valid = engine_energy_valid(mode; tidal_energy_amended)
    if mode == :linearised && !de_valid
        @info "Linearised tidal field (KZ(14) = $(results[:diagnostics].kz[14])): the engine's DE " *
              "excludes the tidal potential; e_jacobi of snapshot_integrals.csv is the conserved energy"
    elseif mode == :linearised
        @info "Linearised tidal field (KZ(14) = $(results[:diagnostics].kz[14])) with the " *
              "tidal-energy-total amendment: the engine's DE measures the error of the Jacobi integral"
    elseif mode != :none
        @info "External field mode $mode (KZ(14) = $(results[:diagnostics].kz[14])): snapshot " *
              "energies are not comparable with the engine's total"
    end
    if haskey(results, :snapshots) && pp.snapshot_integrals
        if _within_pair_sum_limit(results[:snapshots], pp, "Snapshot integrals")
            @info "Snapshot integrals (momenta, energies, bound set)..."
            results[:integrals] = snapshot_integrals(results[:snapshots]; bound..., tidal)
        end
    end
    lap!("integrals")
    if pp.energy_events && haskey(results, :diagnostics) && haskey(results, :engine_events)
        results[:energy_events] = energy_intervals(
            results[:diagnostics],
            results[:engine_events];
            threshold = pp.energy_jump_threshold,
            sigma = pp.energy_jump_sigma,
            engine_energy_valid = de_valid,
        )
        @info "Energy-jump localisation: " * energy_events_summary(results[:energy_events])
    end
    lap!("energy_events")
    if pp.mass_budget && haskey(results, :integrals)
        if haskey(results, :diagnostics)
            results[:mass_budget] = mass_budget(
                results[:integrals],
                get(results, :escapers, EscaperRecord[]),
                get(results, :event_counters, EventCounters[]),
                get(results, :engine_events, EngineEvent[]),
                extract_scaling(results[:diagnostics]),
            )
            @info "Mass budget: " * mass_budget_summary(results[:mass_budget])
            results[:energy_budget] = energy_budget(
                results[:integrals],
                get(results, :escapers, EscaperRecord[]),
                get(results, :event_counters, EventCounters[]),
                results[:diagnostics],
                extract_scaling(results[:diagnostics]);
                bevs = get(results, :binary_evo, BinaryEvolutionSnapshot[]),
                # The integrals are those of the snapshots, in their order
                snapshots = get(results, :snapshots, Snapshot[]),
            )
            @info "Energy budget: " * energy_budget_summary(results[:energy_budget])
        else
            @warn "Mass and energy budgets skipped: the unit scaling needs the stdout capture"
        end
    end
    lap!("budgets")
    summary_path = joinpath(sim_dir, "merger_summary.txt")
    if haskey(results, :snapshots) && length(results[:snapshots]) ≥ 2 && isfile(summary_path)
        ranges = parse_merger_summary(summary_path)
        if !isempty(ranges) &&
           _within_pair_sum_limit(results[:snapshots], pp, "Remnant diagnostics")
            @info "Remnant diagnostics (bound set, core radius, rotation, segregation)..."
            bound_sets =
                haskey(results, :integrals) ? [s.bound for s in results[:integrals]] : nothing
            results[:remnant] = remnant_diagnostics(
                results[:snapshots],
                ranges;
                settings = pp,
                bound_sets,
                bound...,
            )
            lap!("remnant")
            if pp.merger_tree && bound_sets !== nothing
                @info "Progenitor budget" * (length(ranges) ≥ 2 ? " and merger tree..." : "...")
                results[:progenitor_budget] =
                    progenitor_budget(results[:snapshots], ranges, bound_sets; settings = pp)
                length(ranges) ≥ 2 && (
                    results[:merger_tree] =
                        merger_tree(results[:snapshots], ranges; bound_sets, bound...)
                )
            end
        end
    end
    lap!("progenitor_budget_merger_tree")
    # Angular-momentum budget and equipartition: the bound sets of the integrals
    # with the unit scaling of the stdout capture.
    if haskey(results, :integrals) && haskey(results, :diagnostics) && haskey(results, :snapshots)
        bsets = [s.bound for s in results[:integrals]]
        units = extract_scaling(results[:diagnostics])
        if pp.angular_momentum
            @info "Angular-momentum budget and rotation by radius and mass class..."
            results[:angular_momentum] = angular_momentum_budget(
                results[:snapshots],
                bsets,
                units;
                events = get(results, :engine_events, EngineEvent[]),
                class_edges_solar = Tuple(pp.mass_class_edges),
                settings = pp,
                n_min_class = pp.mass_class_n_min,
            )
        end
        lap!("angular_momentum")
        if pp.equipartition
            @info "Equipartition mass of the luminous stars within r_h..."
            results[:equipartition] = equipartition_series(
                results[:snapshots],
                bsets,
                units;
                sevs = get(results, :stellar_evo, StellarEvolutionSnapshot[]),
                settings = pp,
            )
        end
    end
    lap!("equipartition")
    # NBIN0 of the primordial name rule: the &INDATA echo, else the IC record
    # of a generated run (output/merger_ic.toml), else none (labels off)
    nbin0 = haskey(results, :diagnostics) ? results[:diagnostics].nbin0 : -1
    if nbin0 < 0 && (pp.binary_spans || pp.bound_spans)
        ic_path = joinpath(sim_dir, "merger_ic.toml")
        if isfile(ic_path)
            ic = try
                TOML.parsefile(ic_path)
            catch e
                e isa Base.TOML.ParserError || rethrow()
                Dict{String,Any}()
            end
            meta = get(ic, "meta", Dict{String,Any}())
            v = meta isa AbstractDict ? get(meta, "nbin0", get(meta, "n_pairs", nothing)) : nothing
            v isa Integer && v ≥ 0 && (nbin0 = Int(v))
        end
        nbin0 < 0 &&
            @warn "No NBIN0 in the stdout capture or in output/merger_ic.toml: primordial pair " *
                  "labels of the binary spans are off (NBIN0 taken as 0)"
    end
    nbin0 = max(nbin0, 0)
    # Engine-regularised binary spans from the binary lines of the stdout capture
    if pp.binary_spans && !pp.parse_stdout
        @warn "postprocess.binary_spans needs the stdout capture, but parse_stdout = false: " *
              "binary spans skipped"
    elseif pp.binary_spans && haskey(results, :diagnostics)
        diag = results[:diagnostics]
        kz10 = length(diag.kz) ≥ 10 ? diag.kz[10] : -1
        tstar = extract_scaling(diag).tscale
        if kz10 < 2
            @warn "KZ(10) = $(kz10 < 0 ? "unknown (no OPTIONS echo)" : kz10): the engine writes " *
                  "END KSREG lines only for KZ(10) > 1, so no KS episode can be closed; binary " *
                  "spans skipped"
        elseif diag.nzero <= 0
            @warn "No N in the input echo of $(pp.stdout_file): binary spans skipped"
        elseif !(isfinite(tstar) && tstar > 0)
            @warn "No T* in the stdout capture: KS periods cannot be converted, binary spans skipped"
        else
            @info "Binary spans from the KS, few-body and GW records..."
            coll_path = joinpath(sim_dir, "coll.13")
            coal_path = joinpath(sim_dir, "coal.24")
            collisions = isfile(coll_path) ? read_collision_log(coll_path) : CollisionRecord[]
            coalescences = isfile(coal_path) ? read_coalescence_log(coal_path) : CoalescenceRecord[]
            # A restarted segment appends its re-run of the interval after its dump
            restarts = _span_restart_starts(
                _read_toml_or_empty(joinpath(dirname(sim_dir), "RUN_INFO.toml")),
            )
            collisions = _span_restart_dedup(collisions, restarts)
            coalescences = _span_restart_dedup(coalescences, restarts)
            results[:binary_event_log] = read_binary_events(joinpath(sim_dir, pp.stdout_file))
            results[:coalescences] = coalescences
            results[:binary_spans] = binary_spans(
                results[:binary_event_log];
                listings = read_all_ks_pair_listings(sim_dir),
                collisions = collisions,
                coalescences = coalescences,
                nzero = diag.nzero,
                nbin0 = nbin0,
                tstar_myr = tstar,
                n_p = pp.span_n_periods,
            )
        end
    end
    lap!("binary_spans")
    # Physically bound pairs and their spans over the snapshots (level P)
    if pp.bound_spans && haskey(results, :snapshots)
        @info "Bound pairs and their spans over $(length(results[:snapshots])) snapshots..."
        results[:bound_spans] = bound_spans(
            results[:snapshots];
            nbin0 = nbin0,
            gamma_max = pp.bound_gamma_max,
            n_miss = pp.bound_n_miss,
        )
    end
    lap!("bound_spans")
    if haskey(results, :binary_spans) && haskey(results, :bound_spans)
        results[:binary_census] =
            binary_census(results[:binary_spans], results[:bound_spans]; nbin0 = nbin0)
    end
    lap!("binary_census")
    # Binary fractions per shell, class and census level
    if pp.binary_population && haskey(results, :snapshots)
        @info "Binary population over $(length(results[:snapshots])) snapshots..."
        kstar =
            haskey(results, :stellar_evo) ?
            _bpop_pipeline_kstar(
                results[:snapshots],
                results[:stellar_evo],
                get(results, :binary_evo, BinaryEvolutionSnapshot[]),
            ) : nothing
        results[:binary_population] = binary_population_series(
            results[:snapshots];
            radius_frame = Symbol(pp.population_frame),
            edges = isempty(pp.population_edges) ? nothing : pp.population_edges,
            centre_kwargs = _centre_keywords(pp),
            kstar = kstar,
            q_min = pp.population_q_min,
            primary_window_msun = (pp.population_m1_min, pp.population_m1_max),
            gamma_max = pp.bound_gamma_max,
        )
    end
    lap!("binary_population")
    # Pair-type composition, hardness and the orbital-element distributions
    if (pp.binary_pair_types || pp.binary_distributions) && haskey(results, :snapshots)
        kstar_pt =
            haskey(results, :stellar_evo) ?
            _bpop_pipeline_kstar(
                results[:snapshots],
                results[:stellar_evo],
                get(results, :binary_evo, BinaryEvolutionSnapshot[]),
            ) : nothing
        r_tide = pp.population_r_tide_nb > 0 ? pp.population_r_tide_nb : NaN
        if pp.binary_pair_types
            @info "Binary pair types and hardness over $(length(results[:snapshots])) snapshots..."
            results[:binary_pair_types] = binary_pair_type_series(
                results[:snapshots];
                radius_frame = Symbol(pp.population_frame),
                edges = isempty(pp.population_edges) ? nothing : pp.population_edges,
                centre_kwargs = _centre_keywords(pp),
                kstar = kstar_pt,
                gamma_max = pp.bound_gamma_max,
                r_tide_nb = r_tide,
                projection_axis = pp.population_projection_axis,
            )
        end
        if pp.binary_distributions
            @info "Binary orbital-element distributions over $(length(results[:snapshots])) snapshots..."
            results[:binary_distributions] = binary_distribution_series(
                results[:snapshots];
                centre_kwargs = _centre_keywords(pp),
                gamma_max = pp.bound_gamma_max,
                r_tide_nb = r_tide,
            )
        end
    end
    lap!("binary_pair_types")
    # Gravitational-wave catalogues of the compact pairs (event lines, spans, coal.24, esc.11)
    if pp.gw_catalogues && !pp.parse_stdout
        @warn "postprocess.gw_catalogues needs the stdout capture, but parse_stdout = false: " *
              "GW catalogues skipped"
    elseif pp.gw_catalogues && haskey(results, :diagnostics)
        scaling = extract_scaling(results[:diagnostics])
        if !(isfinite(scaling.tscale) && scaling.tscale > 0)
            @warn "No T* in the stdout capture: GW catalogues skipped"
        else
            @info "Gravitational-wave catalogues of the compact pairs..."
            log = get!(results, :binary_event_log) do
                read_binary_events(joinpath(sim_dir, pp.stdout_file))
            end
            coalescences = get!(results, :coalescences) do
                coal_path = joinpath(sim_dir, "coal.24")
                isfile(coal_path) ? read_coalescence_log(coal_path) : CoalescenceRecord[]
            end
            spans = get(results, :binary_spans, BinarySpan[])
            band = (pp.gw_band_hz[1], pp.gw_band_hz[2])
            t_now = pp.gw_t_now_myr > 0 ? pp.gw_t_now_myr : NaN
            results[:gw_catalogues] = (
                compact = compact_binary_catalogue(
                    log,
                    spans;
                    tstar_myr = scaling.tscale,
                    rbar_pc = scaling.rbar,
                    binary_evo = get(results, :binary_evo, BinaryEvolutionSnapshot[]),
                    band_hz = band,
                ),
                escaped = escaped_binary_catalogue(
                    log;
                    tstar_myr = scaling.tscale,
                    rbar_pc = scaling.rbar,
                    vstar_kms = scaling.vstar,
                    escapers = get(results, :escapers, EscaperRecord[]),
                    t_now_myr = t_now,
                    band_hz = band,
                ),
                mergers = gw_merger_catalogue(
                    log,
                    coalescences,
                    spans;
                    tstar_myr = scaling.tscale,
                    rbar_pc = scaling.rbar,
                    vstar_kms = scaling.vstar,
                ),
            )
        end
    end
    lap!("gw_catalogues")
    # Collisions and coalescences with the progenitor cluster of each star
    if pp.collision_catalogue && haskey(results, :diagnostics)
        coll_path = joinpath(sim_dir, "coll.13")
        coal_path = joinpath(sim_dir, "coal.24")
        scaling = extract_scaling(results[:diagnostics])
        if (isfile(coll_path) || isfile(coal_path)) &&
           isfinite(scaling.tscale) &&
           scaling.tscale > 0
            @info "Collision catalogue with the progenitor cluster of each star..."
            collisions = isfile(coll_path) ? read_collision_log(coll_path) : CollisionRecord[]
            coalescences = get!(results, :coalescences) do
                isfile(coal_path) ? read_coalescence_log(coal_path) : CoalescenceRecord[]
            end
            ranges = parse_merger_summary(joinpath(sim_dir, "merger_summary.txt"))
            log = get(results, :binary_event_log, nothing)
            results[:collision_catalogue] = collision_catalogue(
                collisions,
                coalescences,
                ranges;
                tstar_myr = scaling.tscale,
                rbar_pc = scaling.rbar,
                nbin0 = nbin0,
                mass_loss = log === nothing ? nothing : log.mass_loss,
            )
        end
    end
    lap!("collision_catalogue")
    # Lifetimes of the regularised binaries: Kaplan–Meier on the age axis and the
    # cumulative incidence of the end classes (escape and the open ends censor)
    if pp.binary_survival && haskey(results, :binary_spans) && haskey(results, :diagnostics)
        spans = results[:binary_spans]
        t_last = maximum((s.t_end_myr for s in spans); init = 0.0)
        if !isempty(spans) && isfinite(t_last)
            @info "Binary survival statistics over $(length(spans)) spans..."
            table = event_table(spans; time_axis = :age, t_end = t_last)
            results[:binary_survival] = (
                curve = kaplan_meier(table),
                incidence = aalen_johansen(table),
                n_spans = length(spans),
            )
        end
    end
    lap!("binary_survival")
    # Validity record: built from the products above, with the build record of
    # the output directory and the run record of its parent.
    if pp.validity
        results[:validity] = validity_record(
            results,
            pp;
            build_info = _read_toml_or_empty(joinpath(sim_dir, "BUILD_INFO.toml")),
            run_info = _read_toml_or_empty(joinpath(dirname(sim_dir), "RUN_INFO.toml")),
        )
        energy = results[:validity].sections["energy"]
        @info "Validity record: energy $(energy["verdict"])" *
              (isempty(energy["reasons"]) ? "" : " ($(join(energy["reasons"], "; ")))") *
              ", $(length(results[:validity].epochs)) epochs"
    end
    lap!("validity")
    results[:timings] = timings
    @info "Post-processing stage times [s]: " *
          join(("$k $(round(v; digits = 1))" for (k, v) in timings if v ≥ 0.05), ", ") *
          "; total $(round(sum(last, timings; init = 0.0); digits = 1))"

    return results
end

"""
    _tidal_energy_amended(sim_dir) -> Bool

Whether the build record copied into the run's output directory
(`BUILD_INFO.toml`) lists the `tidal-energy-total` source amendment among
its `source_amendments`; false without a record.
"""
function _tidal_energy_amended(sim_dir::AbstractString)::Bool
    path = joinpath(sim_dir, "BUILD_INFO.toml")
    isfile(path) || return false
    record = try
        TOML.parsefile(path)
    catch e
        e isa Base.TOML.ParserError || rethrow()
        @warn "Build record $path is not valid TOML; the engine is taken as unamended"
        return false
    end
    return "tidal-energy-total" in get(record, "source_amendments", String[])
end

"""
    _within_pair_sum_limit(snaps, pp, what) -> Bool

Whether the largest snapshot of `snaps` holds at most `pp.pair_sum_max_n`
stars; warns, naming `what`, the key and the size of the pair sum otherwise.
"""
function _within_pair_sum_limit(
    snaps::Vector{Snapshot},
    pp::PostprocessConfig,
    what::AbstractString,
)::Bool
    n_max = maximum(nparticles, snaps)
    n_max ≤ pp.pair_sum_max_n && return true
    @warn "$what skipped: snapshots of up to $n_max stars exceed postprocess.pair_sum_max_n = " *
          "$(pp.pair_sum_max_n) (each snapshot costs an exact pair sum of " *
          "$(@sprintf("%.1e", float(n_max)^2)) pair evaluations)"
    return false
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
An explicit `run_id` whose run directory already holds a `RUN_INFO.toml` (a
run that was started or finished) is refused with an `ArgumentError` before
anything is built or written; such a run is continued with
[`resume_pipeline`](@ref).

`base_dir` is the project directory: every relative path of the configuration
(`install.install_dir`, `simulation.input_file`, `simulation.runs_dir`,
`postprocess.data_dir`, `merger.config_file`) resolves against it. It defaults
to `cfg.config_dir`, the directory of the file `load_config` read, so
`julia scripts/run_setup.jl path/to/config.toml` keeps its engine and its
runs next to that file, never inside the package.

A simulation phase that ends on a stop request (`simulation.wall_budget`, or
a signal) ends the pipeline there, with an empty result and no completion
stamp; [`resume_pipeline`](@ref) continues it.

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
    # An explicit run ID must not reuse the directory of a started or finished
    # run: its generator and engine phases would overwrite that run's files.
    if !isempty(run_id)
        existing = joinpath(_resolve_path(base_dir, cfg.simulation.runs_dir), run_id)
        isfile(joinpath(existing, "RUN_INFO.toml")) && throw(
            ArgumentError(
                "run_pipeline: the run directory $existing already holds a RUN_INFO.toml " *
                "(a started or finished run); choose another run_id, or continue that run " *
                "with resume_pipeline(\"$existing\")",
            ),
        )
    end

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

    # A run stopped on request is continued by `resume_pipeline`; the phases
    # below belong to the completed run.
    if cfg.simulation.run_test && _run_status(run_dir) == "stopped"
        @info "The simulation stopped on request before its end time; post-processing and " *
              "figures are left to the resumed run (resume_pipeline, scripts/run_resume.jl)"
        return Dict{Symbol,Any}()
    end
    return _finish_pipeline(cfg, run_dir, base_dir, merger_result, phases, t_pipeline)
end

"""
    _finish_pipeline(cfg, run_dir, base_dir, merger_result, phases, t_pipeline)
        -> Dict{Symbol,Any}

Phases 3 and 4 of the pipeline for the data of `run_dir`: post-processing,
the data products, the figures, and the completion stamp. `phases` lists
the phases already performed and is extended; `t_pipeline` is the start
time the stamp's elapsed time is counted from. Shared by
[`run_pipeline`](@ref) and [`resume_pipeline`](@ref).
"""
function _finish_pipeline(
    cfg::Nbody6Config,
    run_dir::AbstractString,
    base_dir::AbstractString,
    merger_result,
    phases::Vector{String},
    t_pipeline::Real,
)::Dict{Symbol,Any}
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
        # An absolute data_dir is not normalised by _resolve_path: normpath and
        # strip a trailing separator so that dirname yields the parent
        data_path = normpath(_resolve_path(base_dir, cfg.postprocess.data_dir))
        length(data_path) > 1 && (data_path = rstrip(data_path, ['/', '\\']))
        dirname(data_path)
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
    if haskey(results, :integrals) && !isempty(products_dir)
        write_snapshot_integrals(
            joinpath(products_dir, "snapshot_integrals.csv"),
            results[:integrals],
        )
    end
    if haskey(results, :energy_events) && !isempty(products_dir)
        write_energy_intervals(
            joinpath(products_dir, "energy_intervals.csv"),
            results[:energy_events],
        )
        write_engine_events(
            joinpath(products_dir, "engine_events.csv"),
            results[:energy_events].events,
        )
    end
    if haskey(results, :mass_budget) && !isempty(products_dir)
        write_mass_budget(joinpath(products_dir, "mass_budget.csv"), results[:mass_budget])
    end
    if haskey(results, :energy_budget) && !isempty(products_dir)
        write_energy_budget(joinpath(products_dir, "energy_budget.csv"), results[:energy_budget])
    end
    if haskey(results, :angular_momentum) && !isempty(products_dir)
        write_angular_momentum(
            joinpath(products_dir, "angular_momentum.csv"),
            results[:angular_momentum],
        )
    end
    if haskey(results, :equipartition) && !isempty(products_dir)
        write_equipartition(joinpath(products_dir, "equipartition.csv"), results[:equipartition])
    end
    if haskey(results, :merger_tree) && !isempty(products_dir)
        write_merger_tree(joinpath(products_dir, "merger_tree.csv"), results[:merger_tree])
    end
    if haskey(results, :binary_spans) && !isempty(products_dir)
        write_binary_spans(joinpath(products_dir, "binary_spans.csv"), results[:binary_spans])
    end
    if haskey(results, :bound_spans) && !isempty(products_dir)
        write_bound_spans(joinpath(products_dir, "bound_spans.csv"), results[:bound_spans])
    end
    if haskey(results, :binary_census) && !isempty(products_dir)
        write_binary_census(joinpath(products_dir, "binary_census.csv"), results[:binary_census])
    end
    if haskey(results, :binary_population) && !isempty(products_dir)
        write_binary_population(
            joinpath(products_dir, "binary_population.csv"),
            results[:binary_population],
        )
    end
    if haskey(results, :binary_pair_types) && !isempty(products_dir)
        write_binary_pair_types(
            joinpath(products_dir, "binary_pair_types.csv"),
            results[:binary_pair_types],
        )
    end
    if haskey(results, :binary_distributions) && !isempty(products_dir)
        write_binary_distributions(
            joinpath(products_dir, "binary_distributions.csv"),
            results[:binary_distributions],
        )
    end
    if haskey(results, :binary_survival) && !isempty(products_dir)
        write_survival_curve(
            joinpath(products_dir, "binary_survival.csv"),
            results[:binary_survival].curve,
        )
        write_cumulative_incidence(
            joinpath(products_dir, "binary_incidence.csv"),
            results[:binary_survival].incidence,
        )
    end
    if haskey(results, :collision_catalogue) && !isempty(products_dir)
        write_collision_catalogue(
            joinpath(products_dir, "collision_catalogue.csv"),
            results[:collision_catalogue],
        )
    end
    if haskey(results, :gw_catalogues) && !isempty(products_dir)
        cat = results[:gw_catalogues]
        write_compact_binary_catalogue(
            joinpath(products_dir, "compact_binary_catalogue.csv"),
            cat.compact,
        )
        write_escaped_binary_catalogue(
            joinpath(products_dir, "escaped_binary_catalogue.csv"),
            cat.escaped,
        )
        write_gw_merger_catalogue(joinpath(products_dir, "gw_merger_catalogue.csv"), cat.mergers)
    end
    if haskey(results, :progenitor_budget) && !isempty(products_dir)
        write_progenitor_budget(
            joinpath(products_dir, "progenitor_budget.csv"),
            results[:progenitor_budget],
        )
    end
    if haskey(results, :validity) && !isempty(products_dir)
        write_validity(joinpath(products_dir, "validity.toml"), results[:validity])
        write_validity_snapshots(
            joinpath(products_dir, "validity_snapshots.csv"),
            results[:validity],
        )
    end

    # ── Phase 4: Plots ──
    if cfg.visualization.enabled && !isempty(results)
        @info "Phase 4: Generating plots..."
        generate_plots(results, cfg; run_dir = products_dir)
        push!(phases, "plots")
    end

    _stamp_pipeline_completion(
        run_dir,
        phases,
        time() - t_pipeline;
        postprocess_seconds = get(results, :timings, nothing),
    )
    return results
end

"""
    resume_pipeline(run_dir; base_dir = nothing) -> Dict{Symbol,Any}

Continue the pipeline of an interrupted run from the configuration frozen
in `run_dir`: one more engine segment through [`resume_run`](@ref), and,
once the run has reached its end time, the post-processing and figure
phases with the completion stamp. The call can be repeated: a run that is
still short of its end time returns an empty result and waits for the next
call, and a run whose pipeline has finished is left alone, so one queue job
per segment needs no knowledge of how many segments the run takes.
"""
function resume_pipeline(
    run_dir::AbstractString;
    base_dir::Union{Nothing,AbstractString} = nothing,
)::Dict{Symbol,Any}
    run_dir = abspath(run_dir)
    cfg = load_config(joinpath(run_dir, "config.toml"))
    base = base_dir === nothing ? cfg.config_dir : String(base_dir)
    cfg.visualization.enabled && _require_plotting(
        :resume_pipeline,
        "Set `[visualization] enabled = false` in the run's config.toml to finish without figures.",
    )
    if _pipeline_completed(run_dir) && _run_status(run_dir) == "completed"
        @info "Run $(basename(run_dir)) is complete and its pipeline has finished; nothing to do"
        return Dict{Symbol,Any}()
    end

    t_pipeline = time()
    resume_run(run_dir; base_dir = base)
    status = _run_status(run_dir)
    if status != "completed"
        @info "Run $(basename(run_dir)) is not at its end time yet (status $(status)); the " *
              "remaining phases are left to a later call"
        return Dict{Symbol,Any}()
    end

    out_dir = joinpath(run_dir, "output")
    merger_result =
        isfile(joinpath(out_dir, "merger_ic.toml")) && isfile(joinpath(out_dir, "dat.10")) ?
        load_merger_ic_result(out_dir) : nothing
    return _finish_pipeline(cfg, run_dir, base, merger_result, String["simulation"], t_pipeline)
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

function __init__()
    Base.Experimental.register_error_hint(_plotting_error_hint, MethodError)
    atexit(_stop_engines_at_exit)
    return nothing
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
# The exported surface is the workflow: configuration, orchestration, the
# readers and their record types, the diagnostics, the merger generator's
# specification types and entry points, and the figure routines. Lower-level
# helpers that a session would rarely call unqualified — accessors with short
# names (`rbar`, `rc`), build internals, IC samplers, sweep-index plumbing —
# are `public`: documented and reachable as `Nbody6Dynamics.name`, not
# brought into scope by `using`.
export Nbody6Config,
    InstallConfig,
    BuildConfig,
    SimulationConfig,
    PostprocessConfig,
    VisualizationConfig,
    MergerPipelineConfig,
    ClusterConfig,
    CostConfig,
    PlotStyle
export load_config, save_config, example_input
export setup_nbody6,
    run_simulation, restart_simulation, resume_run, postprocess, run_pipeline, resume_pipeline
export generate_plots, generate_run_id, export_for_paper, run_gpu_validation
export CostEstimate, estimate_cost, cost_rate
export submit_run, cluster_segment, job_script, read_job_ledger
export scan_output, postprocess_external, OutputScan
export Snapshot, SnapshotHeader, DiagnosticsData, AdjustRecord, LagrangianData, UnitScaling
export EscaperRecord, StellarRecord, StellarEvolutionSnapshot, STELLAR_TYPE_LABELS
export BinaryRecord, BinaryEvolutionSnapshot, TelemetrySample
export read_conf3, read_all_conf3, read_diagnostics, extract_scaling, read_lagr, read_escapers
export cumulative_energy_error, read_energy_totals
export read_stellar_evolution, read_all_stellar_evolution
export read_binary_evolution, read_all_binary_evolution
export read_telemetry, read_run_telemetry
export BinaryPopulation, binary_population
export StellarClass, STELLAR_CLASSES, stellar_class
export HRPopulation, hr_population, hr_populations
export StellarCensus, stellar_census, write_stellar_census
export pair_potential, SystemSet, reduce_pairs, BoundSet, bound_set
export SnapshotIntegrals, snapshot_integrals, write_snapshot_integrals
export tidal_mode, engine_energy_valid
export EngineEvent, EventCounters, ENGINE_EVENT_KINDS, KSTAR_ABSENT
export read_engine_events, read_event_counters
export BinaryEventLog, read_binary_events, DYNAMICS_EVENT_KINDS, EVOLUTION_EVENT_KINDS
export KSEvent, FewBodyEvent, HierarchyEvent, ExchangeEvent, TerminationEvent
export GWEvent, BinaryEscapeEvent, MassLossEvent, TidalEvent, RocheEvent
export KSPairListing, KSPairRecord, read_ks_pair_listing, read_all_ks_pair_listings
export WidePairListing, WidePairRecord, read_wide_pair_listing, read_all_wide_pair_listings
export HierarchyListingRecord, read_hierarchy_listing, HierarchyLogRecord, read_hierarchy_log
export CollisionRecord, read_collision_log, CoalescenceRecord, read_coalescence_log
export CMBodyRecord, read_cmbody_log, BinaryStageRecord, read_binary_stage_log
export DegenerateRecord, read_degenerate_log, EscapedBinaryRecord, read_escaped_binary_log
export HardestBinaryRecord, read_hardest_binary_log
export EnergyInterval, EnergyEventAnalysis, energy_intervals, energy_events_summary
export write_energy_intervals, write_engine_events
export MassBudgetRow, mass_budget, mass_budget_summary, write_mass_budget
export EnergyBudgetRow, energy_budget, energy_budget_summary, write_energy_budget
export ClusterStructure, cluster_structure, bound_fraction, per_cluster_virial, parse_merger_summary
export RadialProfile, radial_profile, cluster_profiles, system_profile, model_density
export RemnantDiagnostics, RotationProfile, remnant_diagnostics, write_remnant_diagnostics
export coalescence_time, core_radius, rotation_analysis, mass_segregation
export shrinking_sphere_centre, lagrangian_grid, relaxation_times, local_relaxation_time
export core_collapse_markers
export ShapeShell, shape_tensor, projected_ellipticity, aplus, two_sample_d, mixing
export SnapshotPair, local_density, snapshot_pairs, reconcile_pairs
export EquipartitionFit, EquipartitionRow, equipartition_fit, equipartition_series
export write_equipartition
export AngularMomentumBudget, angular_momentum_budget, write_angular_momentum
export MergerEvent, MergerTree, merger_tree, write_merger_tree
export ProgenitorBudgetRow, progenitor_budget, write_progenitor_budget
export ValidityRecord, ValidityEpoch, ValiditySnapshot, validity_record
export write_validity, write_validity_snapshots
export BinarySpan, binary_spans, write_binary_spans
export BoundPair, BoundSpan, bound_pairs, bound_spans, write_bound_spans
export CensusSpan, binary_census, write_binary_census
export gw_enhancement,
    gw_time_factor,
    gw_inspiral_time,
    gw_inspiral_time_mandel,
    gw_zwick_factor,
    chirp_mass,
    orbital_frequency,
    gw_frequency,
    eccentricity_at_frequency
export BinaryPopulationRow, binary_population_series, write_binary_population, stellar_type_map
export BinaryPairTypeRow, binary_pair_type_series, write_binary_pair_types
export BinaryDistributionRow, binary_distribution_series, write_binary_distributions
export CompactBinaryRow, EscapedBinaryRow, GWMergerRow, pair_class
export compact_binary_catalogue, escaped_binary_catalogue, gw_merger_catalogue
export write_compact_binary_catalogue, write_escaped_binary_catalogue, write_gw_merger_catalogue
export peters_semi_major_axis, band_residence, eccentricity_after
export RotatingKingModel, solve_rotating_king, rotating_king_moments, sample_rotating_king
export RotatingKingProfile, rotation_moments
export CollisionRow,
    collision_catalogue, collision_summary, progenitor_of, write_collision_catalogue
export EventTable, SurvivalCurve, CumulativeIncidence, event_table, kaplan_meier, nelson_aalen
export aalen_johansen, turnbull, logrank_test, cox_snell_residuals, survival_at, confidence_interval
export write_survival_curve, write_cumulative_incidence
export CalibrationRow, CalibrationTable, CalibratedEstimate, load_calibration, calibrated
export sample_plummer_om, lynden_bell_flip!, stretch_positions!, rank_masses_by_energy!
export SweepConfig, SweepPoint, load_sweep_config, sweep_points, prepare_sweep, run_sweep
export sweep_summary, write_sweep_summary, sweep_visualization
export EnsembleStatistics, ensemble_statistics, sweep_ensembles
export control_merger_dict, write_control_merger_config
export ClusterSpec,
    SegregationSpec,
    BinarySpec,
    MassiveBinarySpec,
    OrbitSpec,
    ComplexSpec,
    MergerOutputSpec,
    Nbody6ParameterSpec,
    StellarSpec,
    TidalSpec,
    MergerConfig,
    MergerICResult
export DensityProfile,
    KingProfile, PlummerProfile, IMFSpec, KroupaIMF, RescaledKroupaIMF, EqualMassIMF, PowerLawIMF
export load_merger_config, generate_merger_ic, run_merger_pipeline, load_merger_ic_result
export engine_interval
export publication_theme, plotting_available, PlottingUnavailable
export plot_snapshot, plot_snapshot_evolution
export plot_lagrangian, plot_energy, plot_particle_count
export plot_hr, plot_hr_evolution
export plot_escapers, plot_escape_anisotropy
export plot_mass_segregation, plot_evolutionary_clock, plot_core_mass
export plot_binary_population, plot_binary_orbital_elements, plot_binary_period_distribution
export plot_cluster_separation, plot_cluster_virial, plot_cluster_structure
export plot_density_profiles, plot_velocity_dispersion
export plot_remnant_rotation, plot_rotation_profile, plot_remnant_structure
export plot_mass_segregation_evolution, remnant_figures
export plot_sweep_lagrangian, plot_sweep_energy, plot_sweep_ensemble, plot_control_comparison
export sweep_figures, plot_merger_ic, plot_telemetry
export animate_cluster, animate_hr, animate_lagrangian

public nparticles, pair_indices, time_nb, time_myr, rbar, zmbar, tscale, vstar, rscale, rc
public to_pc, to_msun, to_myr, to_kms
public detect_platform, check_dependencies, detect_cuda_path
public detect_compute_capabilities, cuda_arch_from_compute_cap, cuda_gencode_flags
public resolve_cuda_arch, nvcc_release, nvcc_supported_archs
public binary_hardness, hardness_scale, binary_scales, semi_major_axis_pc, binding_energy
public stellar_class_index, class_counts, classes_present
public run_sweep_point, read_sweep_index, write_sweep_index
public profile_name, imf_name, expected_mass, sample_masses, kroupa_mean_mass
public sample_plummer, sample_king, sample_kroupa, sample_binaries, expand_binaries
public virialise!,
    kepler_velocity, jacobi_radius, roche_truncation_radius, crossing_time, to_nbody_units!
public write_dat10, generate_merger_inp, resolve_nbody6_parameters
public set_publication_theme!

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
