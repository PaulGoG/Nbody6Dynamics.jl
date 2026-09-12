# =============================================================================
# Binary-population figures (bev.82 diagnostics)
# =============================================================================

"""Bin width in dex for the log-period histograms."""
const _PERIOD_BIN_DEX = 0.5

"""Upper limit of the eccentricity axis; the band above e = 1 holds the annotations."""
const _E_AXIS_TOP = 1.14

"""Legend entry types accepted by the grouped legends (population and ensemble figures)."""
const _LegendElement = Union{LineElement,MarkerElement,PolyElement}

"""Marker of a hardness class: circles for hard pairs, triangles for soft ones."""
_binary_marker(hard::Bool) = hard ? :circle : :utriangle

"""
    plot_binary_population(pop::BinaryPopulation, cfg::VisualizationConfig;
                           filename = "binary_population") -> String

Two stacked panels sharing the time axis: regularised pair counts (total
and, when classified, hard and soft) above the binary fraction
``f_\\mathrm{b} = N_\\mathrm{b} / (N_\\mathrm{s} + N_\\mathrm{b})`` and the hard
fraction ``N_\\mathrm{hard} / N_\\mathrm{pairs}`` in per cent. The net change
of the pair count is annotated in the upper panel.
"""
function Nbody6Dynamics.plot_binary_population(
    pop::BinaryPopulation,
    cfg::VisualizationConfig;
    filename::AbstractString = "binary_population",
)::String
    n = length(pop.time_myr)
    n > 0 || error("No binary population data to plot")

    t = pop.time_myr
    single = n == 1
    tmin, tmax = single ? (t[1] - 0.5, t[1] + 0.5) : extrema(t)
    tticks = single ? _nice_ticks(tmin, tmax; target_n = 5) : _time_ticks(tmin, tmax)
    c_tot = _SEMANTIC_COLORS[:n_pairs]
    c_hard = _SEMANTIC_COLORS[:binary_hard]
    c_soft = _SEMANTIC_COLORS[:binary_soft]

    fig = Figure(; size = _fig_multipanel(cfg, 2, 1))
    ax1 = Axis(fig[1, 1]; ylabel = "Regularised pairs", xticks = tticks, xticklabelsvisible = false)
    ax2 = Axis(fig[2, 1]; xlabel = L"t\;[\mathrm{Myr}]", ylabel = "Fraction [%]", xticks = tticks)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)

    # Counts (markers alone for a single epoch)
    if single
        scatter!(ax1, t, pop.n_pairs; color = c_tot, markersize = 14, label = "All pairs")
        if pop.classified
            scatter!(
                ax1,
                t,
                pop.n_hard;
                color = c_hard,
                marker = :circle,
                markersize = 12,
                label = "Hard",
            )
            scatter!(
                ax1,
                t,
                pop.n_soft;
                color = c_soft,
                marker = :utriangle,
                markersize = 12,
                label = "Soft",
            )
        end
    else
        lines!(ax1, t, pop.n_pairs; color = c_tot, linewidth = 2.5, label = "All pairs")
        if pop.classified
            lines!(ax1, t, pop.n_hard; color = c_hard, linewidth = 2.0, label = "Hard")
            lines!(
                ax1,
                t,
                pop.n_soft;
                color = c_soft,
                linewidth = 2.0,
                linestyle = :dash,
                label = "Soft",
            )
        end
    end
    n0, n1 = pop.n_pairs[1], pop.n_pairs[end]
    if !single && n0 > 0
        Δ = n1 - n0
        pct = @sprintf("%.2g", 100 * abs(Δ) / n0)
        sign = Δ < 0 ? "-" : "+"
        _annotate!(
            ax1,
            latexstring("N_\\mathrm{pairs}: $(n0) \\rightarrow $(n1)\\;($(sign)$(pct)\\,\\%)");
            corner = :br,
            color = c_tot,
        )
    end

    # Fractions
    fb = 100 .* pop.binary_fraction
    any(isfinite, fb) && (
        single ? scatter!(ax2, t, fb; color = c_tot, markersize = 14, label = L"f_\mathrm{b}") :
        lines!(ax2, t, fb; color = c_tot, linewidth = 2.5, label = L"f_\mathrm{b}")
    )
    if pop.classified
        fh = [np > 0 ? 100 * nh / np : NaN for (nh, np) in zip(pop.n_hard, pop.n_pairs)]
        single ?
        scatter!(
            ax2,
            t,
            fh;
            color = c_hard,
            marker = :circle,
            markersize = 12,
            label = "Hard fraction",
        ) : lines!(ax2, t, fh; color = c_hard, linewidth = 2.0, label = "Hard fraction")
    end
    if !any(isfinite, fb) && !pop.classified
        _no_data_note!(ax2, "No stellar count or energy scale available")
    end
    ylims!(ax2, 0.0, nothing)

    # One legend, two families: the pair counts of the upper panel and the
    # fractions of the lower one (colour encodes the quantity in both).
    # Concretely typed vectors: Makie reads `Vector{Any}` labels as multi-line text.
    _el(c; marker = :circle, linestyle = :solid, lw = 2.0) =
        single ? MarkerElement(; color = c, marker = marker, markersize = 12) :
        LineElement(; color = c, linewidth = lw, linestyle = linestyle)
    count_entries = _LegendElement[_el(c_tot; lw = 2.5)]
    count_labels = AbstractString["All pairs"]
    fraction_entries = _LegendElement[]
    fraction_labels = AbstractString[]
    any(isfinite, fb) &&
        (push!(fraction_entries, _el(c_tot; lw = 2.5)); push!(fraction_labels, L"f_\mathrm{b}"))
    if pop.classified
        push!(count_entries, _el(c_hard))
        push!(count_labels, L"Hard ($E_\mathrm{b} > \langle m \rangle \sigma^2$)")
        push!(count_entries, _el(c_soft; marker = :utriangle, linestyle = :dash))
        push!(count_labels, "Soft")
        push!(fraction_entries, _el(c_hard))
        push!(fraction_labels, L"N_\mathrm{hard} / N_\mathrm{pairs}")
    end
    groups = [count_entries]
    group_labels = [count_labels]
    titles = ["Pairs:"]
    if !isempty(fraction_entries)
        push!(groups, fraction_entries)
        push!(group_labels, fraction_labels)
        push!(titles, "Fractions:")
    end
    if sum(length, groups) ≥ 2
        Legend(
            fig[0, :],
            groups,
            group_labels,
            titles;
            orientation = :horizontal,
            nbanks = sum(length, groups) > 3 ? 2 : 1,   # keep the legend within the figure width
            framevisible = false,
            titleposition = :left,
            tellheight = true,
            padding = (0, 0, 0, 0),
        )
    end

    return _save_fig(cfg, filename, fig)
end

"""
    plot_binary_orbital_elements(bev::BinaryEvolutionSnapshot, cfg::VisualizationConfig;
                                 m_mean = NaN, sigma_kms = NaN,
                                 filename = "binary_orbital_elements") -> String

Semi-major axis [au] against eccentricity of every regularised pair at one
epoch. With a hard/soft energy scale (`m_mean` [M☉], `sigma_kms` [km s⁻¹])
the pairs are classed by the Heggie criterion, the class fractions are
annotated, and the boundary semi-major axis of a pair of the mean
component-mass product is marked.
"""
function Nbody6Dynamics.plot_binary_orbital_elements(
    bev::BinaryEvolutionSnapshot,
    cfg::VisualizationConfig;
    m_mean::Real = NaN,
    sigma_kms::Real = NaN,
    filename::AbstractString = "binary_orbital_elements",
)::String
    recs = bev.records
    classified = isfinite(m_mean) && isfinite(sigma_kms) && m_mean > 0 && sigma_kms > 0
    t_val = @sprintf("%.3g", bev.time_myr)

    fig = Figure(; size = _figsize_px(cfg))
    a_au = [semi_major_axis_pc(r) / _AU_IN_PC for r in recs]
    e = [r.eccentricity for r in recs]
    finite = [isfinite(a) && a > 0 && isfinite(ee) for (a, ee) in zip(a_au, e)]
    a_au, e = a_au[finite], e[finite]
    lo, hi = isempty(a_au) ? (1e-2, 1e4) : (minimum(a_au) / 1.5, maximum(a_au) * 1.5)

    ax = Axis(
        fig[1, 1];
        xlabel = L"a\;[\mathrm{au}]",
        ylabel = L"e",
        xscale = log10,
        xticks = _log_ticks(lo, hi),
        yticks = 0.0:0.2:1.0,
    )
    xlims!(ax, lo, hi)
    # e ≤ 1 by construction: the band above unity is a data-free annotation strip.
    ylims!(ax, -0.02, _E_AXIS_TOP)
    _annotate!(ax, latexstring("t = $(t_val)\\;\\mathrm{Myr}"); corner = :tl)

    if isempty(a_au)
        _no_data_note!(ax, "No regularised pairs")
        return _save_fig(cfg, filename, fig)
    end

    c_tot = _SEMANTIC_COLORS[:n_pairs]
    c_hard = _SEMANTIC_COLORS[:binary_hard]
    c_soft = _SEMANTIC_COLORS[:binary_soft]
    n = length(a_au)
    ms = clamp(cfg.style.marker_budget / max(n, 1), cfg.style.marker_min, cfg.style.marker_max)

    if classified
        x = binary_hardness(bev, m_mean, sigma_kms)[finite]
        hard = x .> 1.0
        n_hard = count(hard)
        n_soft = n - n_hard
        p_hard = @sprintf("%.3g", 100 * n_hard / n)
        p_soft = @sprintf("%.3g", 100 * n_soft / n)
        classes = (
            (hard, c_hard, true, "Hard: $(n_hard) ($(p_hard) %)"),
            (.!hard, c_soft, false, "Soft: $(n_soft) ($(p_soft) %)"),
        )
        for (mask, c, is_hard, label) in classes
            any(mask) || continue
            scatter!(
                ax,
                a_au[mask],
                e[mask];
                color = (c, 0.75),
                marker = _binary_marker(is_hard),
                markersize = ms,
                label = label,
            )
        end
        # Boundary for a pair of the mean component-mass product
        m1m2 = sum(r.mass1 * r.mass2 for r in recs[finite]) / n
        a_hs_au = _G_PC_KMS2_MSUN * m1m2 / (2 * m_mean * sigma_kms^2) / _AU_IN_PC
        if lo < a_hs_au < hi
            vlines!(
                ax,
                [a_hs_au];
                color = :black,
                linestyle = :dash,
                linewidth = 1.5,
                label = L"Hard/soft boundary, $\langle m_1 m_2 \rangle$",
            )
        end
        n_series = (n_hard > 0) + (n_soft > 0) + (lo < a_hs_au < hi)
        n_series ≥ 2 && _top_legend!(fig, ax)
    else
        scatter!(ax, a_au, e; color = (c_tot, 0.75), markersize = ms)
        _annotate!(ax, latexstring("N_\\mathrm{pairs} = $(n)"); corner = :tl, color = c_tot)
    end

    return _save_fig(cfg, filename, fig)
end

"""
    plot_binary_period_distribution(bevs::Vector{BinaryEvolutionSnapshot},
                                    cfg::VisualizationConfig;
                                    filename = "binary_period_distribution") -> String

Histogram of ``\\log_{10}(P/\\mathrm{d})`` of the regularised pairs at the
first epoch (filled) and, when more than one snapshot is given, at the last
epoch (dashed outline), with the pair counts annotated.
"""
function Nbody6Dynamics.plot_binary_period_distribution(
    bevs::AbstractVector{BinaryEvolutionSnapshot},
    cfg::VisualizationConfig;
    filename::AbstractString = "binary_period_distribution",
)::String
    isempty(bevs) && error("No binary snapshots to plot")
    first_bev, last_bev = bevs[1], bevs[end]
    logp0 = filter(isfinite, [r.log_period_days for r in first_bev.records])
    logp1 = filter(isfinite, [r.log_period_days for r in last_bev.records])
    two = length(bevs) > 1

    fig = Figure(; size = _figsize_px(cfg))
    all_p = two ? vcat(logp0, logp1) : logp0
    lo = isempty(all_p) ? 0.0 : floor(minimum(all_p) / _PERIOD_BIN_DEX) * _PERIOD_BIN_DEX
    hi = isempty(all_p) ? 1.0 : ceil(maximum(all_p) / _PERIOD_BIN_DEX) * _PERIOD_BIN_DEX
    hi ≤ lo && (hi = lo + _PERIOD_BIN_DEX)
    edges = collect(lo:_PERIOD_BIN_DEX:hi)

    ax = Axis(
        fig[1, 1];
        xlabel = L"\log_{10}(P \, / \, \mathrm{d})",
        ylabel = "Number of pairs",
        xticks = _nice_ticks(lo, hi),
    )
    xlims!(ax, lo, hi)

    if isempty(all_p)
        _no_data_note!(ax, "No regularised pairs")
        return _save_fig(cfg, filename, fig)
    end

    c0 = _SEMANTIC_COLORS[:n_pairs]
    t0 = @sprintf("%.3g", first_bev.time_myr)
    t1 = @sprintf("%.3g", last_bev.time_myr)
    label0 = latexstring("t = $(t0)\\;\\mathrm{Myr}")
    hist!(
        ax,
        logp0;
        bins = edges,
        color = (c0, 0.35),
        strokecolor = _band_edge(c0),
        strokewidth = 1.5,
        label = label0,
    )
    if two
        label1 = latexstring("t = $(t1)\\;\\mathrm{Myr}")
        stephist!(
            ax,
            logp1;
            bins = edges,
            color = :black,
            linestyle = :dash,
            linewidth = 2.0,
            label = label1,
        )
        _top_legend!(fig, ax)
        _annotate!(
            ax,
            latexstring("N_\\mathrm{pairs}: $(length(logp0)) \\rightarrow $(length(logp1))");
            corner = :tl,
            color = c0,
        )
    else
        _annotate!(
            ax,
            latexstring("t = $(t0)\\;\\mathrm{Myr},\\; N_\\mathrm{pairs} = $(length(logp0))");
            corner = :tl,
            color = c0,
        )
    end
    ylims!(ax, 0.0, nothing)

    return _save_fig(cfg, filename, fig)
end
