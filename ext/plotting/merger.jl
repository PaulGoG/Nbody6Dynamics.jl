# =============================================================================
# Merger-specific diagnostics
# =============================================================================
# Track each initial cluster through its own bound set over the course of the
# run. Gives meaningful pre-coalescence diagnostics for multi-cluster setups
# where global Lagrangian radii measure cluster separation rather than
# internal structure.

"""
    plot_cluster_separation(snaps, cluster_ranges, cfg;
                            filename = "merger_cluster_separation")

Plot the separation of the density centres of every pair of initial clusters
against time, with the half-mass threshold `min(r_h,i, r_h,j)` below which
the pair counts as merged, both from [`coalescence_time`](@ref). For
n_clusters ≤ 5 every pair is shown individually, its threshold dashed in the
same colour; for more clusters the min/max envelope and mean of the
separations are drawn with the mean threshold dashed and the number of
distinct clusters (`n_distinct`) on a twin axis. Separations are shown in pc
and times in Myr when `cfg.units == "physical"` (header AS scaling).

The coalescence time of [`coalescence_time`](@ref) is marked when the
clusters coalesce within the snapshots. After that epoch the initial cluster
decomposition is no longer physically meaningful.
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
    # One bound set per cluster and snapshot, shared by every series below
    coal = coalescence_time(snaps, cluster_ranges)

    # Density-centre separations and half-mass thresholds: pairs × n_times
    pairs = [(i, j) for i in 1:(n_cl - 1) for j in (i + 1):n_cl]
    seps = copy(coal.separation)
    thr = copy(coal.threshold)
    # NB length → pc
    if physical
        for k in 1:n_t
            s = rbar(snaps[k].header)
            seps[:, k] .*= s
            thr[:, k] .*= s
        end
    end
    t_c = coal.index === nothing ? NaN : (physical ? coal.time_myr : coal.time_nb)
    t_c_label =
        isfinite(t_c) ?
        latexstring(
            "t_\\mathrm{coalesce} = $(_fmt_latex_sig3(t_c))\\;",
            physical ? "\\mathrm{Myr}" : "[\\mathrm{NB}]",
        ) : L""
    d_all = filter(isfinite, vcat(vec(seps), vec(thr)))

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = physical ? L"d_{ij} \; [\mathrm{pc}]" : L"d_{ij} \; [\mathrm{NB}]",
        xticks = _time_ticks(first(t), last(t)),
    )
    # Upper limit of the separation axis, from which the event label hangs
    y_top = NaN
    if isempty(d_all)
        _no_data_note!(ax, "No cluster with enough bound members")
    else
        y_lo, y_top = _pad_limits(0.0, maximum(d_all); floor = 0)
        ylims!(ax, y_lo, y_top)
        ax.yticks = _nice_ticks(y_lo, y_top)
    end
    # Two legend families; headers above their entries keep the legend within
    # the figure width
    family_titles = ["Separation of density centres:", "Half-mass threshold:"]
    legend_kw = (;
        orientation = :horizontal,
        framevisible = false,
        titleposition = :top,
        tellheight = true,
        padding = (0, 0, 0, 0),
    )

    detail = n_cl ≤ 5
    if detail
        # One colour per pair: separation solid, its threshold dashed
        pair_entries = _LegendElement[]
        pair_labels = AbstractString[]
        for (p, (i, j)) in enumerate(pairs)
            valid = isfinite.(@view seps[p, :])
            any(valid) || continue
            color = _OKABE_ITO[mod1(p, length(_OKABE_ITO))]
            lines!(ax, t[valid], seps[p, valid]; color = color, linewidth = _STYLE.data)
            lines!(
                ax,
                t[valid],
                thr[p, valid];
                color = color,
                linestyle = :dash,
                linewidth = _STYLE.guide,
            )
            push!(pair_entries, LineElement(; color = color, linewidth = _STYLE.data))
            push!(pair_labels, "$(i)–$(j)")
        end
        isfinite(y_top) && _event_marker!(ax, t_c, t_c_label, y_top, (first(t), last(t)))
        if !isempty(pair_entries)
            Legend(
                fig[0, :],
                [
                    pair_entries,
                    _LegendElement[LineElement(;
                        color = :black,
                        linestyle = :dash,
                        linewidth = _STYLE.guide,
                    )],
                ],
                [pair_labels, AbstractString[L"\min(r_{h,i},\, r_{h,j})"]],
                family_titles;
                nbanks = min(3, cld(length(pair_entries), 4)),
                legend_kw...,
            )
        end
    else
        # Envelope: min/max band plus mean line, mean threshold dashed
        d_min, d_max, d_mean = _envelope_stats(seps)
        _, _, thr_mean = _envelope_stats(thr)

        valid = .!isnan.(d_mean)
        sep_color = _SEMANTIC_COLORS[:separation]
        # Darker same-hue edges on the band fill
        sep_edge = _band_edge(sep_color)
        if any(valid)
            band!(ax, t[valid], d_min[valid], d_max[valid]; color = (sep_color, 0.25))
            lines!(ax, t[valid], d_min[valid]; color = sep_edge, linewidth = _STYLE.envelope)
            lines!(ax, t[valid], d_max[valid]; color = sep_edge, linewidth = _STYLE.envelope)
            lines!(ax, t[valid], d_mean[valid]; color = sep_color, linewidth = _STYLE.data)
            lines!(
                ax,
                t[valid],
                thr_mean[valid];
                color = sep_color,
                linestyle = :dash,
                linewidth = _STYLE.guide,
            )
        end
        isfinite(y_top) && _event_marker!(ax, t_c, t_c_label, y_top, (first(t), last(t)))

        # Twin axis: number of distinct clusters, the connected components of
        # the merged-pair graph of `coalescence_time`. Tick/label colour
        # matches the count series (semantic N colour).
        survivor_count = coal.n_distinct
        count_color = _SEMANTIC_COLORS[:n_particles]
        ax2 = Axis(
            fig[1, 1];
            yaxisposition = :right,
            ylabel = L"N_\mathrm{distinct}",
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
        stairs!(
            ax2,
            t,
            Float64.(survivor_count);
            color = count_color,
            linewidth = _STYLE.data,
            step = :post,
        )

        # The count is named by its twin axis in its own colour; a third
        # family would carry the legend past the figure width
        if any(valid)
            Legend(
                fig[0, :],
                [
                    _LegendElement[
                        PolyElement(;
                            color = (sep_color, 0.25),
                            strokecolor = sep_edge,
                            strokewidth = _STYLE.envelope,
                        ),
                        LineElement(; color = sep_color, linewidth = _STYLE.data),
                    ],
                    _LegendElement[LineElement(;
                        color = sep_color,
                        linestyle = :dash,
                        linewidth = _STYLE.guide,
                    )],
                ],
                [AbstractString["Min–max range", "Mean"], AbstractString["Mean"]],
                family_titles;
                nbanks = 1,
                legend_kw...,
            )
        end
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
proxy for the time of coalescence. `bound` passes `max_passes` and
`min_members` to [`bound_set`](@ref).
"""
@publication function Nbody6Dynamics.plot_cluster_virial(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    filename::AbstractString = "merger_cluster_virial",
    bound::NamedTuple = (;),
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping virial plot: no cluster ranges"; return nothing)
    n_t = length(snaps)
    n_t ≥ 2 || (@info "Skipping virial plot: need ≥ 2 snapshots"; return nothing)

    physical = cfg.units == "physical" && all(_has_physical_scaling(s.header) for s in snaps)
    t = physical ? [time_myr(s.header) for s in snaps] : [time_nb(s.header) for s in snaps]
    @info "Computing per-cluster virial ratios ($(n_cl) clusters × $(n_t) snapshots)..."
    Q, _ = per_cluster_virial(snaps, cluster_ranges; bound...)

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
times in Myr when `cfg.units == "physical"`. `bound` passes `max_passes` and
`min_members` to [`bound_set`](@ref).
"""
@publication function Nbody6Dynamics.plot_cluster_structure(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    lagr::Union{Nothing,LagrangianData} = nothing,
    filename::AbstractString = "merger_cluster_structure",
    bound::NamedTuple = (;),
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping cluster structure plot: no cluster ranges"; return nothing)
    n_t = length(snaps)
    n_t ≥ 2 || (@info "Skipping cluster structure plot: need ≥ 2 snapshots"; return nothing)

    physical = cfg.units == "physical" && all(_has_physical_scaling(s.header) for s in snaps)
    t = physical ? [time_myr(s.header) for s in snaps] : [time_nb(s.header) for s in snaps]
    r_unit = physical ? rbar(snaps[1].header) : 1.0
    @info "Computing per-cluster structure ($(n_cl) clusters × $(n_t) snapshots)..."
    st = cluster_structure(snaps, cluster_ranges; bound...)
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

    # Radii of every series of the upper panel, for its limits
    r_panel = filter(isfinite, vec(r_h))

    # Engine's global 50 % radius about its single density centre
    if lagr !== nothing && !isempty(lagr.time)
        i50 = findfirst(==(0.5), lagr.mass_fractions)
        if i50 !== nothing
            t_l = physical ? lagr.time .* tscale(snaps[1].header) : lagr.time
            r_50 = lagr.radii[i50, :] .* r_unit
            append!(r_panel, filter(isfinite, r_50))
            lines!(
                ax1,
                t_l,
                r_50;
                color = :gray40,
                linestyle = :dash,
                linewidth = _STYLE.fit,
                label = L"\mathrm{Engine}\;r_{50}\;\mathrm{(global\;centre)}",
            )
            n_series += 1
        end
    end

    hlines!(ax2, [1.0]; color = :gray50, linestyle = :dash, linewidth = _STYLE.guide)
    # Occupancy of the bound-fraction panel, for the guide label and the limits
    f_x = Float64[]
    f_y = Float64[]
    for i in axes(f_bound, 1), k in axes(f_bound, 2)
        isfinite(f_bound[i, k]) || continue
        push!(f_x, t[k])
        push!(f_y, f_bound[i, k])
    end
    isempty(f_y) ||
        _guide_label!(ax2, L"M_{\mathrm{bound},i} = M_i", 1.0; xs = f_x, ys = f_y, color = :gray40)
    f_lo, f_hi = _pad_limits(isempty(f_y) ? 0.0 : minimum(f_y), 1.0; floor = 0)
    ylims!(ax2, f_lo, f_hi)
    ax2.yticks = _nice_ticks(f_lo, f_hi)
    if !isempty(r_panel)
        r_lo, r_hi = _pad_limits(0.0, maximum(r_panel); floor = 0)
        ylims!(ax1, r_lo, r_hi)
        ax1.yticks = _nice_ticks(r_lo, r_hi)
    end
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
`bound` passes `max_passes` and `min_members` to [`bound_set`](@ref).
"""
@publication function Nbody6Dynamics.plot_density_profiles(
    snap::Snapshot,
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    specs::Union{Nothing,AbstractVector} = nothing,   # Vector{ClusterSpec}; the type is defined later in the module
    filename::AbstractString = "merger_density_profiles",
    bound::NamedTuple = (;),
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
    profiles = cluster_profiles(snap, cluster_ranges; bound...)
    any(!isnothing, profiles) ||
        (@info "Skipping density profiles: no cluster with enough members"; return nothing)
    with_model = specs !== nothing

    # The ratio panel is an auxiliary strip, not a second main panel.
    fig = Figure(; size = with_model ? _fig_with_strip(cfg) : _figsize_px(cfg))
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
        ) : nothing
    with_model && _strip_rows!(fig, 1, 2)

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
        q_x = reduce(vcat, r_ratio; init = Float64[])
        q_y = reduce(vcat, ratios; init = Float64[])
        isempty(q_y) || _guide_label!(
            ax2,
            L"\rho = \rho_{\mathrm{model}}",
            1.0;
            xs = q_x,
            ys = q_y,
            color = :gray40,
            log = true,
            xlog = true,
        )
        # Symmetric about the model in log, holding the data: a factor of two
        # at least, so that agreement is not magnified into disagreement.
        q_fin = [q for q in q_y if isfinite(q) && q > 0]
        m = (isempty(q_fin) ? log10(2) : max(maximum(abs.(log10.(q_fin))), log10(2))) + 0.1
        ylims!(ax2, 10.0^-m, 10.0^m)
        ax2.yticks = _log_ticks(10.0^-m, 10.0^m)
        linkxaxes!(ax1, ax2)
        rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    end
    if !isempty(ρ_all)
        # One tick style for every epoch of a run, whatever range it spans.
        ax1.yticks = _log_ticks(minimum(ρ_all), maximum(ρ_all); plain = false)
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
`σ_r` and `β` are drawn. `bound` passes `max_passes` and `min_members` to
[`bound_set`](@ref).
"""
@publication function Nbody6Dynamics.plot_velocity_dispersion(
    snap::Snapshot,
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    filename::AbstractString = "merger_velocity_dispersion",
    bound::NamedTuple = (;),
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping velocity dispersion: no cluster ranges"; return nothing)
    h = snap.header
    physical = cfg.units == "physical" && _has_physical_scaling(h)
    r_unit = physical ? rbar(h) : 1.0
    v_unit = physical ? vstar(h) : 1.0
    profiles = cluster_profiles(snap, cluster_ranges; bound...)
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
    # Occupancy of both panels, for their limits and the anisotropy guide label
    sig_all = Float64[]
    beta_x = Float64[]
    beta_y = Float64[]
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
            append!(sig_all, filter(isfinite, prof.sigma_r[valid] .* v_unit))
            append!(sig_all, filter(isfinite, prof.sigma_t[valid] .* v_unit))
            for (x, b) in zip(r, prof.beta[valid])
                isfinite(b) || continue
                push!(beta_x, x)
                push!(beta_y, b)
            end
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
            if ax === ax1
                append!(sig_all, filter(isfinite, vcat(lo[valid], hi[valid])))
            else
                for (x, b) in zip(vcat(grid[valid], grid[valid]), vcat(lo[valid], hi[valid]))
                    isfinite(b) || continue
                    push!(beta_x, x)
                    push!(beta_y, b)
                end
            end
        end
        n_series = 2
    end
    if detail && n_series ≥ 1
        # Colour names the cluster, line style the component of the dispersion.
        lines!(
            ax1,
            r_all[1:1],
            [NaN];
            color = :black,
            linewidth = _STYLE.data,
            label = L"\sigma_r\;\mathrm{(solid)}",
        )
        lines!(
            ax1,
            r_all[1:1],
            [NaN];
            color = :black,
            linestyle = :dash,
            linewidth = _STYLE.fit,
            label = L"\sigma_t\;\mathrm{(dashed)}",
        )
        n_series += 2
    end
    hlines!(ax2, [0.0]; color = :gray50, linestyle = :dash, linewidth = _STYLE.guide)
    isempty(beta_y) || _guide_label!(
        ax2,
        L"\beta = 0\;\mathrm{(isotropic)}",
        0.0;
        xs = beta_x,
        ys = beta_y,
        color = :gray40,
        xlog = true,
    )
    isempty(sig_all) || ylims!(ax1, _pad_limits(0.0, maximum(sig_all); floor = 0)...)
    isempty(beta_y) ||
        ylims!(ax2, _pad_limits(min(0.0, minimum(beta_y)), max(0.0, maximum(beta_y)))...)
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
