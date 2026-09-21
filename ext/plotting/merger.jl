# =============================================================================
# Merger-specific diagnostics
# =============================================================================
# Track each initial cluster's centre-of-mass over the course of the run.
# Gives meaningful pre-coalescence diagnostics for multi-cluster setups where
# global Lagrangian radii measure cluster separation rather than internal
# structure.

"""Coalescence heuristic: max/mean pairwise-separation ratio below which the ensemble is considered merged."""
const _COALESCENCE_RATIO_MAX = 2.0

"""Coalescence heuristic: mean separation must also drop below this multiple of the initial minimum separation."""
const _COALESCENCE_MEAN_SEP_FACTOR = 1.5

"""
    plot_cluster_separation(snaps, cluster_ranges, cfg;
                            filename = "merger_cluster_separation")

Plot the evolution of pairwise centre-of-mass separations between the initial
clusters. For small n_clusters (≤ 5) every pair is shown individually; for
larger N, min/max/mean envelopes are plotted instead.  Separations are shown
in pc and times in Myr when `cfg.units == "physical"` (header AS scaling).

The time at which the ratio `max_sep / mean_sep` drops below a heuristic
threshold is marked as an estimated "coalescence time". After that epoch the
initial cluster decomposition is no longer physically meaningful.
"""
@publication function Nbody6Dynamics.plot_cluster_separation(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    filename::AbstractString = "merger_cluster_separation",
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 2 || (@info "Skipping separation plot: need ≥ 2 clusters (got $n_cl)"; return nothing)
    n_t = length(snaps)
    n_t ≥ 2 || (@info "Skipping separation plot: need ≥ 2 snapshots"; return nothing)

    physical = cfg.units == "physical" && all(_has_physical_scaling(s.header) for s in snaps)
    t = physical ? [time_myr(s.header) for s in snaps] : [time_nb(s.header) for s in snaps]
    coms, present, r_rms = _cluster_com_trajectories(snaps, cluster_ranges)

    # Pairwise separations over time: pairs × n_times
    pairs = [(i, j) for i in 1:(n_cl - 1) for j in (i + 1):n_cl]
    n_pairs = length(pairs)
    seps = fill(NaN, n_pairs, n_t)
    for (p, (i, j)) in enumerate(pairs)
        for k in 1:n_t
            (present[i, k] && present[j, k]) || continue
            dx = coms[1, i, k] - coms[1, j, k]
            dy = coms[2, i, k] - coms[2, j, k]
            dz = coms[3, i, k] - coms[3, j, k]
            seps[p, k] = sqrt(dx*dx + dy*dy + dz*dz)
        end
    end
    # NB length → pc (uniform scaling; ratio-based coalescence heuristics
    # below are unaffected)
    if physical
        for k in 1:n_t
            seps[:, k] .*= rbar(snaps[k].header)
        end
    end

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = physical ? L"d_{ij} \; [\mathrm{pc}]" : L"d_{ij} \; [\mathrm{NB}]",
        xticks = _time_ticks(first(t), last(t)),
    )

    detail = n_cl ≤ 5
    if detail
        # One line per pair, labelled
        n_lines = 0
        for (p, (i, j)) in enumerate(pairs)
            valid = .!isnan.(seps[p, :])
            any(valid) || continue
            lines!(
                ax,
                t[valid],
                seps[p, valid];
                color = _OKABE_ITO[mod1(p, length(_OKABE_ITO))],
                linewidth = _STYLE.data,
                label = "$(i)–$(j)",
            )
            n_lines += 1
        end
        n_lines ≥ 2 && _top_legend!(fig, ax; nbanks = min(3, cld(n_pairs, 4)))
        # First-class coalescence time: all initial clusters spatially merged
        coal = coalescence_time(snaps, cluster_ranges)
        if coal.index !== nothing
            t_c = physical ? coal.time_myr : coal.time_nb
            vlines!(ax, [t_c]; color = :black, linestyle = :dash, linewidth = _STYLE.guide)
            text!(
                ax,
                t_c,
                maximum(filter(!isnan, vec(seps)));
                text = latexstring(
                    "t_\\mathrm{coalesce} = $(_fmt_latex_sig3(t_c))\\;",
                    physical ? "\\mathrm{Myr}" : "[\\mathrm{NB}]",
                ),
                align = (:left, :top),
                offset = (4, 0),
                fontsize = _ANNOTATION_FONTSIZE,
            )
        end
    else
        # Envelope: show min/max band plus mean line
        d_min, d_max, d_mean = _envelope_stats(seps)

        valid = .!isnan.(d_mean)
        sep_color = _SEMANTIC_COLORS[:separation]
        band_plot = band!(ax, t[valid], d_min[valid], d_max[valid]; color = (sep_color, 0.25))
        # Darker same-hue edges on the band fill
        sep_edge = _band_edge(sep_color)
        lines!(ax, t[valid], d_min[valid]; color = sep_edge, linewidth = _STYLE.envelope)
        lines!(ax, t[valid], d_max[valid]; color = sep_edge, linewidth = _STYLE.envelope)
        mean_plot = lines!(ax, t[valid], d_mean[valid]; color = sep_color, linewidth = _STYLE.data)

        # Mark estimated coalescence: max/mean ratio and mean separation both
        # drop below their heuristic thresholds (consts at top of file)
        d_init_min = minimum(filter(!isnan, d_min))
        ratio = [isnan(d_mean[k]) ? NaN : d_max[k] / d_mean[k] for k in 1:n_t]
        idx_merge = findfirst(
            k ->
                !isnan(ratio[k]) &&
                ratio[k] < _COALESCENCE_RATIO_MAX &&
                d_mean[k] < _COALESCENCE_MEAN_SEP_FACTOR * d_init_min,
            1:n_t,
        )

        merge_plot = nothing
        merge_label = ""
        if idx_merge !== nothing && t[idx_merge] > first(t)
            t_merge = t[idx_merge]
            merge_plot =
                vlines!(ax, [t_merge]; color = :black, linestyle = :dash, linewidth = _STYLE.guide)
            merge_label = latexstring("t_\\mathrm{merge} \\approx $(round(t_merge; digits=2))")
        end

        # Twin axis: count of spatially-distinct clusters. Two clusters are
        # merged when |COM_i − COM_j| < (r_rms_i + r_rms_j), i.e. their
        # member spheres overlap. Connected components via union-find.
        # Tick/label colour matches the count series (semantic N colour).
        survivor_count = _count_spatial_clusters(coms, r_rms, present; overlap_factor = 1.0)
        count_color = _SEMANTIC_COLORS[:n_particles]
        ax2 = Axis(
            fig[1, 1];
            yaxisposition = :right,
            ylabel = L"N_\mathrm{clusters}\;\mathrm{(\geq 3\;members)}",
            ylabelcolor = count_color,
            yticklabelcolor = count_color,
            ytickcolor = count_color,
            # Overlay axis: grid off so it doesn't double-draw over ax
            ygridvisible = false,
            xgridvisible = false,
            xticklabelsvisible = false,
            xticksvisible = false,
            topspinevisible = false,
            bottomspinevisible = false,
            leftspinevisible = false,
            rightspinevisible = true,
            yticks = _integer_ticks(0, n_cl),
            limits = ((nothing, nothing), (0, n_cl + max(1, ceil(Int, 0.1 * n_cl)))),
        )
        hidespines!(ax2, :t, :b, :l)
        surv_plot = stairs!(
            ax2,
            t,
            Float64.(survivor_count);
            color = count_color,
            linewidth = _STYLE.data,
            step = :post,
        )

        # Manual legend combining primary and twin axes: horizontal, above
        legend_elems = Any[band_plot, mean_plot]
        legend_labels = [L"\min\;-\;\max\;\mathrm{range}", L"\mathrm{mean}"]
        if merge_plot !== nothing
            push!(legend_elems, merge_plot)
            push!(legend_labels, merge_label)
        end
        push!(legend_elems, surv_plot)
        push!(legend_labels, L"\mathrm{surviving\;clusters}")
        _top_legend!(fig, legend_elems, legend_labels)
    end

    return _save_fig(cfg, filename, fig)
end

"""
    plot_cluster_virial(snaps, cluster_ranges, cfg;
                        filename = "merger_cluster_virial")

Plot the internal virial ratio Q_i(t) of each initial cluster. For
n_clusters ≤ 5 every cluster is shown individually; for larger N, a min/max
envelope with the mean is drawn. A reference line at Q = 0.5 marks virial
equilibrium.  The time axis is in Myr when `cfg.units == "physical"`
(Q is dimensionless).

**Interpretation note:** the metric tracks particles by their original
cluster ID and, by default, only those still bound to that group
(`per_cluster_virial(...; bound_only = true)`), so stripped stars and
kicked stellar remnants do not enter. It is physically meaningful *before*
coalescence (which clusters are still in virial equilibrium, which undergo
violent relaxation). After merger the bound subset of an ID-group is
whatever remains self-bound inside the remnant, so the curve loses its
meaning as a cluster diagnostic; the departure from `Q ≈ 0.5` is a rough
proxy for the time of coalescence.
"""
@publication function Nbody6Dynamics.plot_cluster_virial(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    filename::AbstractString = "merger_cluster_virial",
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping virial plot: no cluster ranges"; return nothing)
    n_t = length(snaps)
    n_t ≥ 2 || (@info "Skipping virial plot: need ≥ 2 snapshots"; return nothing)

    physical = cfg.units == "physical" && all(_has_physical_scaling(s.header) for s in snaps)
    t = physical ? [time_myr(s.header) for s in snaps] : [time_nb(s.header) for s in snaps]
    @info "Computing per-cluster virial ratios ($(n_cl) clusters × $(n_t) snapshots)..."
    Q, _ = per_cluster_virial(snaps, cluster_ranges)

    fig = Figure(; size = _figsize_px(cfg))

    # Pick a log scale when any cluster Q spikes high (mergers routinely do)
    q_all = filter(!isnan, vec(Q))
    use_log = !isempty(q_all) && maximum(q_all) > cfg.style.q_log_threshold

    # Clamp tiny/zero values on the log scale only
    Q_plot = use_log ? max.(Q, cfg.style.q_floor) : copy(Q)
    q_plotted = filter(!isnan, vec(Q_plot))

    ax = Axis(
        fig[1, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = L"Q_i = T_i / |W_i|",
        xticks = _time_ticks(first(t), last(t)),
        yscale = use_log ? log10 : identity,
        yticks = use_log && !isempty(q_plotted) ? _log_ticks(extrema(q_plotted)...) :
                 Makie.automatic,
    )

    detail = n_cl ≤ 5
    n_series = 0
    if detail
        for i in 1:n_cl
            valid = .!isnan.(@view Q[i, :])
            any(valid) || continue
            lines!(
                ax,
                t[valid],
                Q_plot[i, valid];
                color = _OKABE_ITO[mod1(i, length(_OKABE_ITO))],
                linewidth = _STYLE.data,
                label = latexstring("\\mathrm{Cluster}\\;$(i)"),
            )
            n_series += 1
        end
    else
        q_min, q_max, q_mean = _envelope_stats(Q_plot)
        valid = .!isnan.(q_mean)
        vir_color = _SEMANTIC_COLORS[:virial]
        band!(
            ax,
            t[valid],
            q_min[valid],
            q_max[valid];
            color = (vir_color, 0.25),
            label = L"\min\;-\;\max\;\mathrm{range}",
        )
        # Darker same-hue edges on the band fill
        vir_edge = _band_edge(vir_color)
        lines!(ax, t[valid], q_min[valid]; color = vir_edge, linewidth = _STYLE.envelope)
        lines!(ax, t[valid], q_max[valid]; color = vir_edge, linewidth = _STYLE.envelope)
        lines!(
            ax,
            t[valid],
            q_mean[valid];
            color = vir_color,
            linewidth = _STYLE.data,
            label = L"\mathrm{mean}",
        )
        n_series = 2
    end

    # Virial equilibrium reference
    hlines!(
        ax,
        [0.5];
        color = :gray50,
        linestyle = :dash,
        linewidth = _STYLE.guide,
        label = L"Q = 0.5\;\mathrm{(virial\;equilibrium)}",
    )

    # ≥ 2 legend entries whenever at least one data series was drawn
    n_series ≥ 1 && _top_legend!(fig, ax; nbanks = detail ? min(3, cld(n_cl + 1, 4)) : 1)

    return _save_fig(cfg, filename, fig)
end

# -----------------------------------------------------------------------------
# Per-cluster structure: bound half-mass radius and bound mass fraction
# -----------------------------------------------------------------------------

"""
    plot_cluster_structure(snaps, cluster_ranges, cfg; lagr = nothing,
                           filename = "merger_cluster_structure")

Two stacked panels sharing the time axis: the half-mass radius of every
initial cluster measured about its own centre from its bound members
([`cluster_structure`](@ref)), and the bound mass fraction. With `lagr`
given, the engine's global 50 % Lagrangian radius about its single density
centre is overlaid dashed grey on the radius panel, which makes the
difference between the configuration and its members visible. For more
than five clusters the min–max envelope and mean are drawn. Radii in pc and
times in Myr when `cfg.units == "physical"`.
"""
@publication function Nbody6Dynamics.plot_cluster_structure(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    lagr::Union{Nothing,LagrangianData} = nothing,
    filename::AbstractString = "merger_cluster_structure",
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping cluster structure plot: no cluster ranges"; return nothing)
    n_t = length(snaps)
    n_t ≥ 2 || (@info "Skipping cluster structure plot: need ≥ 2 snapshots"; return nothing)

    physical = cfg.units == "physical" && all(_has_physical_scaling(s.header) for s in snaps)
    t = physical ? [time_myr(s.header) for s in snaps] : [time_nb(s.header) for s in snaps]
    r_unit = physical ? rbar(snaps[1].header) : 1.0
    @info "Computing per-cluster structure ($(n_cl) clusters × $(n_t) snapshots)..."
    st = cluster_structure(snaps, cluster_ranges)
    r_h = st.r_lagr[2, :, :] .* r_unit
    f_bound = st.bound_mass_fraction

    fig = Figure(; size = _fig_two_panel(cfg))
    ttk = _time_ticks(first(t), last(t))
    ax1 = Axis(
        fig[1, 1];
        ylabel = physical ? L"r_{h,i} \; [\mathrm{pc}]" : L"r_{h,i} \; [\mathrm{NB}]",
        xticklabelsvisible = false,
        xticks = ttk,
    )
    ax2 = Axis(
        fig[2, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = L"M_{\mathrm{bound},i} \, / \, M_i",
        xticks = ttk,
    )

    detail = n_cl ≤ 5
    n_series = 0
    if detail
        for i in 1:n_cl
            valid = .!isnan.(@view r_h[i, :])
            any(valid) || continue
            color = _OKABE_ITO[mod1(i, length(_OKABE_ITO))]
            lines!(
                ax1,
                t[valid],
                r_h[i, valid];
                color = color,
                linewidth = _STYLE.data,
                label = latexstring("\\mathrm{Cluster}\\;$(i)"),
            )
            lines!(ax2, t[valid], f_bound[i, valid]; color = color, linewidth = _STYLE.data)
            n_series += 1
        end
    else
        sep_color = _SEMANTIC_COLORS[:separation]
        sep_edge = _band_edge(sep_color)
        for (ax, data, labelled) in ((ax1, r_h, true), (ax2, f_bound, false))
            lo, hi, mean = _envelope_stats(data)
            valid = .!isnan.(mean)
            any(valid) || continue
            if labelled
                band!(
                    ax,
                    t[valid],
                    lo[valid],
                    hi[valid];
                    color = (sep_color, 0.25),
                    label = L"\min\;-\;\max\;\mathrm{range}",
                )
                lines!(
                    ax,
                    t[valid],
                    mean[valid];
                    color = sep_color,
                    linewidth = _STYLE.data,
                    label = L"\mathrm{mean}",
                )
                n_series = 2
            else
                band!(ax, t[valid], lo[valid], hi[valid]; color = (sep_color, 0.25))
                lines!(ax, t[valid], mean[valid]; color = sep_color, linewidth = _STYLE.data)
            end
            lines!(ax, t[valid], lo[valid]; color = sep_edge, linewidth = _STYLE.envelope)
            lines!(ax, t[valid], hi[valid]; color = sep_edge, linewidth = _STYLE.envelope)
        end
    end

    # Engine's global 50 % radius about its single density centre
    if lagr !== nothing && !isempty(lagr.time)
        i50 = findfirst(==(0.5), lagr.mass_fractions)
        if i50 !== nothing
            t_l = physical ? lagr.time .* tscale(snaps[1].header) : lagr.time
            lines!(
                ax1,
                t_l,
                lagr.radii[i50, :] .* r_unit;
                color = :gray40,
                linestyle = :dash,
                linewidth = _STYLE.fit,
                label = L"\mathrm{Engine}\;r_{50}\;\mathrm{(global\;centre)}",
            )
            n_series += 1
        end
    end

    hlines!(ax2, [1.0]; color = :gray50, linestyle = :dash, linewidth = _STYLE.guide)
    _annotate!(ax2, L"M_{\mathrm{bound},i} = M_i"; corner = :tl, dy = 0.10, color = :gray40)
    ylims!(ax2, 0.0, 1.08)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    n_series ≥ 2 && _top_legend!(fig, ax1; nbanks = min(3, cld(n_series, 4)))

    return _save_fig(cfg, filename, fig)
end

# -----------------------------------------------------------------------------
# Radial profiles: density against the generating model, velocity dispersions
# -----------------------------------------------------------------------------

"""
    plot_density_profiles(snap, cluster_ranges, cfg; specs = nothing,
                          filename = "merger_density_profiles")

Density profiles `ρ(r)` of every initial cluster about its own centre from
its bound members ([`cluster_profiles`](@ref)), log–log, with the
generating model overlaid dashed when `specs` (the clusters'
[`ClusterSpec`](@ref)s) are given: each model is evaluated at the cluster's
initial half-mass radius and its bound mass at this snapshot. A ratio strip
`ρ / ρ_model` with a labelled unity guide sits beneath the main panel. For
more than five clusters the min–max envelope and mean of the ratio are
drawn. Radii in pc and densities in M☉ pc⁻³ when `cfg.units == "physical"`.
"""
@publication function Nbody6Dynamics.plot_density_profiles(
    snap::Snapshot,
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    specs::Union{Nothing,AbstractVector} = nothing,   # Vector{ClusterSpec}; the type is defined later in the module
    filename::AbstractString = "merger_density_profiles",
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping density profiles: no cluster ranges"; return nothing)
    specs === nothing ||
        length(specs) == n_cl ||
        throw(ArgumentError("specs must have one entry per cluster range"))
    h = snap.header
    physical = cfg.units == "physical" && _has_physical_scaling(h)
    r_unit = physical ? rbar(h) : 1.0
    ρ_unit = physical ? zmbar(h) / rbar(h)^3 : 1.0
    profiles = cluster_profiles(snap, cluster_ranges)
    any(!isnothing, profiles) ||
        (@info "Skipping density profiles: no cluster with enough members"; return nothing)
    with_model = specs !== nothing

    fig = Figure(; size = with_model ? _fig_two_panel(cfg) : _figsize_px(cfg))
    ax1 = Axis(
        fig[1, 1];
        xlabel = with_model ? "" : (physical ? L"r \; [\mathrm{pc}]" : L"r \; [\mathrm{NB}]"),
        ylabel = physical ? L"\rho \; [\mathrm{M}_\odot\,\mathrm{pc}^{-3}]" :
                 L"\rho \; [\mathrm{NB}]",
        xscale = log10,
        yscale = log10,
        xticklabelsvisible = !with_model,
    )
    ax2 =
        with_model ?
        Axis(
            fig[2, 1];
            xlabel = physical ? L"r \; [\mathrm{pc}]" : L"r \; [\mathrm{NB}]",
            ylabel = L"\rho \, / \, \rho_{\mathrm{model}}",
            xscale = log10,
            yscale = log10,
            yticks = _log_ticks(0.1, 10.0),
        ) : nothing

    detail = n_cl ≤ 5
    r_all = Float64[]
    ρ_all = Float64[]
    ratios = Vector{Vector{Float64}}()
    r_ratio = Vector{Vector{Float64}}()
    n_series = 0
    for i in 1:n_cl
        prof = profiles[i]
        prof === nothing && continue
        valid = .!isnan.(prof.rho) .& (prof.rho .> 0)
        any(valid) || continue
        r = prof.r[valid] .* r_unit
        ρ = prof.rho[valid] .* ρ_unit
        append!(r_all, r)
        append!(ρ_all, ρ)
        color = _OKABE_ITO[mod1(i, length(_OKABE_ITO))]
        if detail
            lines!(
                ax1,
                r,
                ρ;
                color = color,
                linewidth = _STYLE.data,
                label = latexstring("\\mathrm{Cluster}\\;$(i)"),
            )
            n_series += 1
        else
            lines!(ax1, r, ρ; color = (:gray50, 0.5), linewidth = _STYLE.ghost)
        end
        if with_model
            spec = specs[i]
            r_h_model = spec.rbar / rbar(h)              # initial half-mass radius [NB]
            ρ_model = model_density(spec.profile, prof.M, r_h_model)
            ρm = [ρ_model(x) for x in prof.r[valid]] .* ρ_unit
            ok = ρm .> 0
            if detail && any(ok)
                lines!(ax1, r[ok], ρm[ok]; color = color, linestyle = :dash, linewidth = _STYLE.fit)
                lines!(ax2, r[ok], ρ[ok] ./ ρm[ok]; color = color, linewidth = _STYLE.data)
            end
            push!(ratios, ρ[ok] ./ ρm[ok])
            push!(r_ratio, r[ok])
        end
    end
    if !detail && with_model && !isempty(ratios)
        # Envelope on a common radial grid: interpolate each ratio onto the union grid
        grid = exp10.(
            range(log10(minimum(minimum, r_ratio)), log10(maximum(maximum, r_ratio)); length = 24),
        )
        mat = fill(NaN, length(ratios), length(grid))
        for (i, (rr, q)) in enumerate(zip(r_ratio, ratios))
            for (k, g) in enumerate(grid)
                j = searchsortedlast(rr, g)
                (j < 1 || j ≥ length(rr)) && continue
                f = (g - rr[j]) / (rr[j + 1] - rr[j])
                mat[i, k] = (1 - f) * q[j] + f * q[j + 1]
            end
        end
        lo, hi, mean = _envelope_stats(mat)
        valid = .!isnan.(mean)
        c = _SEMANTIC_COLORS[:separation]
        band!(
            ax2,
            grid[valid],
            lo[valid],
            hi[valid];
            color = (c, 0.25),
            label = L"\min\;-\;\max\;\mathrm{range}",
        )
        lines!(
            ax2,
            grid[valid],
            mean[valid];
            color = c,
            linewidth = _STYLE.data,
            label = L"\mathrm{mean}",
        )
        lines!(ax1, r_all[1:1], ρ_all[1:1]; color = :gray50, label = L"\mathrm{Clusters}")
        n_series = 2
    end
    if with_model && detail
        lines!(
            ax1,
            r_all[1:1],
            ρ_all[1:1];
            color = :black,
            linestyle = :dash,
            linewidth = _STYLE.fit,
            label = L"\mathrm{Generating\;model}",
        )
        n_series += 1
    end
    if with_model
        hlines!(ax2, [1.0]; color = :gray50, linestyle = :dash, linewidth = _STYLE.guide)
        _annotate!(ax2, L"\rho = \rho_{\mathrm{model}}"; corner = :tr, color = :gray40)
        # Within a factor of ten of the model; the outermost shells at a King
        # model's tidal edge run away and are clipped.
        ylims!(ax2, 0.1, 10.0)
        linkxaxes!(ax1, ax2)
        rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    end
    if !isempty(ρ_all)
        ax1.yticks = _log_ticks(minimum(ρ_all), maximum(ρ_all))
        xt = _log_ticks(minimum(r_all), maximum(r_all))
        ax1.xticks = xt
        ax2 === nothing || (ax2.xticks = xt)
    end
    n_series ≥ 2 && _top_legend!(fig, ax1; nbanks = min(3, cld(n_series, 4)))
    return _save_fig(cfg, filename, fig)
end

"""
    plot_velocity_dispersion(snap, cluster_ranges, cfg;
                             filename = "merger_velocity_dispersion")

Radial (solid) and one-dimensional tangential (dashed) velocity dispersion
profiles of every initial cluster from its bound members, with the
anisotropy `β(r) = 1 − σ_t²/σ_r²` beneath and the isotropic `β = 0` guide
labelled. Velocities in km s⁻¹ and radii in pc when `cfg.units ==
"physical"`. For more than five clusters the min–max envelope and mean of
`σ_r` and `β` are drawn.
"""
@publication function Nbody6Dynamics.plot_velocity_dispersion(
    snap::Snapshot,
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    filename::AbstractString = "merger_velocity_dispersion",
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping velocity dispersion: no cluster ranges"; return nothing)
    h = snap.header
    physical = cfg.units == "physical" && _has_physical_scaling(h)
    r_unit = physical ? rbar(h) : 1.0
    v_unit = physical ? vstar(h) : 1.0
    profiles = cluster_profiles(snap, cluster_ranges)
    any(!isnothing, profiles) ||
        (@info "Skipping velocity dispersion: no cluster with enough members"; return nothing)

    fig = Figure(; size = _fig_two_panel(cfg))
    ax1 = Axis(
        fig[1, 1];
        ylabel = physical ? L"\sigma \; [\mathrm{km\,s^{-1}}]" : L"\sigma \; [\mathrm{NB}]",
        xscale = log10,
        xticklabelsvisible = false,
    )
    ax2 = Axis(
        fig[2, 1];
        xlabel = physical ? L"r \; [\mathrm{pc}]" : L"r \; [\mathrm{NB}]",
        ylabel = L"\beta = 1 - \sigma_t^2 / \sigma_r^2",
        xscale = log10,
    )
    detail = n_cl ≤ 5
    r_all = Float64[]
    n_series = 0
    rows_sr = Vector{Vector{Float64}}()
    rows_beta = Vector{Vector{Float64}}()
    rows_r = Vector{Vector{Float64}}()
    for i in 1:n_cl
        prof = profiles[i]
        prof === nothing && continue
        valid = .!isnan.(prof.sigma_r)
        any(valid) || continue
        r = prof.r[valid] .* r_unit
        append!(r_all, r)
        color = _OKABE_ITO[mod1(i, length(_OKABE_ITO))]
        if detail
            lines!(
                ax1,
                r,
                prof.sigma_r[valid] .* v_unit;
                color = color,
                linewidth = _STYLE.data,
                label = latexstring("\\mathrm{Cluster}\\;$(i)"),
            )
            lines!(
                ax1,
                r,
                prof.sigma_t[valid] .* v_unit;
                color = color,
                linestyle = :dash,
                linewidth = _STYLE.fit,
            )
            lines!(ax2, r, prof.beta[valid]; color = color, linewidth = _STYLE.data)
            n_series += 1
        else
            push!(rows_sr, prof.sigma_r[valid] .* v_unit)
            push!(rows_beta, prof.beta[valid])
            push!(rows_r, r)
        end
    end
    if !detail && !isempty(rows_r)
        grid = exp10.(
            range(log10(minimum(minimum, rows_r)), log10(maximum(maximum, rows_r)); length = 24),
        )
        function onto(rows_y)
            mat = fill(NaN, length(rows_y), length(grid))
            for (i, (rr, q)) in enumerate(zip(rows_r, rows_y))
                for (k, g) in enumerate(grid)
                    j = searchsortedlast(rr, g)
                    (j < 1 || j ≥ length(rr)) && continue
                    f = (g - rr[j]) / (rr[j + 1] - rr[j])
                    mat[i, k] = (1 - f) * q[j] + f * q[j + 1]
                end
            end
            return mat
        end
        c = _SEMANTIC_COLORS[:virial]
        for (ax, mat, labelled) in ((ax1, onto(rows_sr), true), (ax2, onto(rows_beta), false))
            lo, hi, mean = _envelope_stats(mat)
            valid = .!isnan.(mean)
            any(valid) || continue
            if labelled
                band!(
                    ax,
                    grid[valid],
                    lo[valid],
                    hi[valid];
                    color = (c, 0.25),
                    label = L"\sigma_r\;\min\;-\;\max",
                )
                lines!(
                    ax,
                    grid[valid],
                    mean[valid];
                    color = c,
                    linewidth = _STYLE.data,
                    label = L"\sigma_r\;\mathrm{mean}",
                )
            else
                band!(ax, grid[valid], lo[valid], hi[valid]; color = (c, 0.25))
                lines!(ax, grid[valid], mean[valid]; color = c, linewidth = _STYLE.data)
            end
        end
        n_series = 2
    end
    if detail && n_series ≥ 1
        lines!(
            ax1,
            r_all[1:1],
            [NaN];
            color = :black,
            linestyle = :dash,
            linewidth = _STYLE.fit,
            label = L"\sigma_t\;\mathrm{(dashed)}",
        )
        n_series += 1
    end
    hlines!(ax2, [0.0]; color = :gray50, linestyle = :dash, linewidth = _STYLE.guide)
    _annotate!(ax2, L"\beta = 0\;\mathrm{(isotropic)}"; corner = :tr, color = :gray40)
    if !isempty(r_all)
        xt = _log_ticks(minimum(r_all), maximum(r_all))
        ax1.xticks = xt
        ax2.xticks = xt
    end
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    n_series ≥ 2 && _top_legend!(fig, ax1; nbanks = min(3, cld(n_series, 4)))
    return _save_fig(cfg, filename, fig)
end
