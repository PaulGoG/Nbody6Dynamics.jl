# =============================================================================
# Figure dispatcher
# =============================================================================
# The single entry point every caller reaches: the config-driven pipeline
# (`run_pipeline`), `postprocess_external`, and users driving the readers by
# hand all route through `generate_plots`.

"""
    _plot_binary_diagnostics(bevs, results, vis; pair_sum_max_n)

Binary-population figures for `bevs`. The hard/soft energy scale and the
stellar count per epoch come from the conf.3 snapshots when present
([`binary_scales`](@ref)); without snapshots the stellar count falls back
to the ADJUST diagnostics and the pairs are left unclassified. The scale
is taken over the bound systems only while no snapshot holds more than
`pair_sum_max_n` particles (an O(N²) selection), over every system above
that, with a warning.
"""
function _plot_binary_diagnostics(
    bevs::Vector{BinaryEvolutionSnapshot},
    results::Dict{Symbol,Any},
    vis::VisualizationConfig;
    pair_sum_max_n::Integer,
)
    snaps = get(results, :snapshots, Snapshot[])::Vector{Snapshot}
    n_max = isempty(snaps) ? 0 : maximum(nparticles, snaps)
    n_max ≤ pair_sum_max_n ||
        @warn "Binary hard/soft scale taken over all systems, not the bound set: " *
              "snapshots of up to $n_max particles exceed postprocess.pair_sum_max_n = " *
              "$pair_sum_max_n (O(N²) pair sums)"
    scales = binary_scales(bevs, snaps; pair_sum_max_n)
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
                   sim_dir::AbstractString = "", animations::Bool = true,
                   pair_sum_max_n::Integer = PostprocessConfig().pair_sum_max_n)

Generate all available plots (and, when `animations = true`, GIF
animations) from post-processing results, saving into `vis.output_dir`.
`sim_dir` names the directory holding the raw simulation output; when given,
it is searched for `merger_summary.txt` to produce the merger-specific
figures (inter-cluster separation, per-cluster virial ratio).

`pair_sum_max_n` bounds the merger figures built on O(N²) pair sums, the
per-cluster virial ratio and structure, the bound-member density and
velocity-dispersion profiles and the remnant figures: they are drawn only
when no snapshot holds more particles than this, with a warning otherwise;
the inter-cluster separation is drawn regardless. The bound-set selection
behind the binary figures' hard/soft scale obeys the same limit, above
which the scale is taken over every system. The default is that of
`PostprocessConfig`, and the `Nbody6Config` method passes
`cfg.postprocess.pair_sum_max_n`.

This is the single plot dispatcher — both the config-driven pipeline
(via the `Nbody6Config` method) and [`postprocess_external`](@ref) route
through it.
"""
function Nbody6Dynamics.generate_plots(
    results::Dict{Symbol,Any},
    vis::VisualizationConfig;
    sim_dir::AbstractString = "",
    animations::Bool = true,
    pair_sum_max_n::Integer = PostprocessConfig().pair_sum_max_n,
)
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
                        n_max = maximum(nparticles, snaps)
                        if n_max ≤ pair_sum_max_n
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
                                isfile(ic_meta) ? load_merger_ic_result(sim_dir).cluster_specs :
                                nothing
                            specs === nothing ||
                                length(specs) == length(ranges) ||
                                (specs = nothing)
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
                            # Computed by `postprocess` for pipeline runs (and written
                            # there as remnant_diagnostics.csv); derived here only for
                            # results that did not pass through it.
                            diag = get(results, :remnant, nothing)
                            if diag === nothing
                                @info "Remnant diagnostics (bound set, core radius, rotation, segregation)..."
                                diag = remnant_diagnostics(snaps, ranges)
                            end
                            remnant_figures(diag::RemnantDiagnostics, vis)
                        else
                            @warn "Per-cluster virial ratio, structure, profiles and " *
                                  "remnant figures skipped: snapshots of up to $n_max " *
                                  "particles exceed postprocess.pair_sum_max_n = " *
                                  "$pair_sum_max_n (O(N²) pair sums)"
                        end
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
            hr_bevs = get(results, :binary_evo, BinaryEvolutionSnapshot[])
            # Three HR diagrams: beginning, middle, end
            mid = max(1, length(sevs) ÷ 2)
            hr_epochs = [
                (1, "hr_diagram_early"),
                (mid, "hr_diagram_mid"),
                (length(sevs), "hr_diagram_final"),
            ]
            for (idx, fname) in hr_epochs
                @info "Plotting HR diagram (epoch $idx/$(length(sevs)))..."
                plot_hr(sevs, vis; epoch = idx, bevs = hr_bevs, filename = fname)
            end
            if length(sevs) > 1
                @info "Plotting HR evolution..."
                plot_hr_evolution(sevs, vis; bevs = hr_bevs, filename = "hr_evolution")
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
            _plot_binary_diagnostics(bevs, results, vis; pair_sum_max_n)
        end
    end

    # Run telemetry: the sampler's CSVs live in the run directory above the
    # output directory (absent for output produced outside the pipeline).
    if !isempty(sim_dir) && isdir(sim_dir)
        samples = read_run_telemetry(dirname(abspath(sim_dir)))
        if length(samples) ≥ 2
            @info "Plotting run telemetry..."
            plot_telemetry(samples, vis; filename = "telemetry")
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
                animate_hr(
                    sevs,
                    vis;
                    bevs = get(results, :binary_evo, BinaryEvolutionSnapshot[]),
                    filename = "hr_evolution_anim",
                )
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
function Nbody6Dynamics.generate_plots(
    results::Dict{Symbol,Any},
    cfg::Nbody6Config;
    run_dir::AbstractString = "",
)
    # Build a VisualizationConfig with the output_dir resolved to the run
    vis = if !isempty(run_dir)
        plots_dir = joinpath(run_dir, cfg.visualization.output_dir)
        # Every field is carried over by construction, so a setting added to
        # the struct cannot be dropped here.
        fields = (f => getfield(cfg.visualization, f) for f in fieldnames(VisualizationConfig))
        VisualizationConfig(; fields..., output_dir = plots_dir)
    else
        cfg.visualization
    end
    sim_dir = !isempty(run_dir) ? joinpath(run_dir, "output") : ""
    return generate_plots(
        results,
        vis;
        sim_dir = sim_dir,
        pair_sum_max_n = cfg.postprocess.pair_sum_max_n,
    )
end
