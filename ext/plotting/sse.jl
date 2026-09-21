# =============================================================================
# SSE-quantity plots from sev.83_* snapshots (v2026.07+ hrplot.F columns)
# =============================================================================

"""
    plot_mass_segregation(sev::StellarEvolutionSnapshot, cfg::VisualizationConfig;
                          filename = "mass_segregation")

Mass-segregation diagnostic from one stellar evolution snapshot: distance
from the density centre RI [pc] (log x) against stellar mass [M☉] (log y),
with the snapshot epoch annotated.  Records with non-finite or non-positive
RI or mass are skipped (both axes are logarithmic).

Warns and returns `nothing` when no record passes the validity filter;
otherwise saves the figure and returns `nothing`.
"""
function Nbody6Dynamics.plot_mass_segregation(
    sev::StellarEvolutionSnapshot,
    cfg::VisualizationConfig;
    filename::AbstractString = "mass_segregation",
)
    valid = [
        r for r in sev.records if
        isfinite(r.ri) && r.ri > 0 && isfinite(r.mass_solar) && r.mass_solar > 0
    ]
    isempty(valid) && (@warn "No stellar records with positive RI and mass"; return nothing)

    ri = [r.ri for r in valid]
    m = [r.mass_solar for r in valid]

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = L"r\;[\mathrm{pc}]",
        ylabel = L"m\;[\mathrm{M}_\odot]",
        xscale = log10,
        yscale = log10,
        xticks = _log_ticks(extrema(ri)...),
        yticks = _log_ticks(extrema(m)...),
        xgridvisible = false,
        ygridvisible = false,
    )
    scatter!(
        ax,
        ri,
        m;
        color = _OKABE_ITO[1],
        markersize = _marker_size(cfg, length(valid)),
        strokewidth = 0,
    )

    t_val = _fmt_latex_sig3(sev.time_myr)
    _annotate!(ax, latexstring("t = $(t_val)\\;\\mathrm{Myr}"); corner = :tr)

    _save_fig(cfg, filename, fig)
    return nothing
end

"""
    plot_evolutionary_clock(sev::StellarEvolutionSnapshot, cfg::VisualizationConfig;
                            filename = "evolutionary_clock")

Evolutionary-clock histogram for the main-sequence population (K* ≤ 1) of
one snapshot: the age fraction t/T_MS, where t is the snapshot epoch [Myr]
and T_MS the per-star main-sequence lifetime TM [Myr] from SSE, clamped to
[0, 2].  Stars past turnoff (t/T_MS > 1) are counted in an annotation and
the T_MS boundary is marked.  Records with non-finite or non-positive TM
are skipped.

Warns and returns `nothing` when no main-sequence record carries a valid
TM; otherwise saves the figure and returns `nothing`.
"""
function Nbody6Dynamics.plot_evolutionary_clock(
    sev::StellarEvolutionSnapshot,
    cfg::VisualizationConfig;
    filename::AbstractString = "evolutionary_clock",
)
    ms_records = [
        r for r in sev.records if
        r.stellar_type ≤ 1 && isfinite(r.ms_lifetime_myr) && r.ms_lifetime_myr > 0
    ]
    isempty(ms_records) &&
        (@warn "No main-sequence records with valid MS lifetime TM"; return nothing)

    frac = [clamp(sev.time_myr / r.ms_lifetime_myr, 0.0, 2.0) for r in ms_records]
    n_ms = length(frac)

    # Histogram over [0, hi]: hi sits just above the largest fraction so the
    # last bin is right-open like the others (floor 0.05 guards the t = 0
    # snapshot where every fraction is zero).
    hi = max(maximum(frac) * 1.05, 0.05)
    nbins = 24
    edges = collect(range(0.0, hi; length = nbins + 1))
    counts = [count(f -> edges[j] ≤ f < edges[j + 1], frac) for j in 1:nbins]
    c_max = maximum(counts)

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = L"t / T_\mathrm{MS}",
        ylabel = L"N",
        xticks = _nice_ticks(0.0, hi),
        yticks = _integer_ticks(0, c_max),
        limits = (0.0, hi, 0.0, max(c_max, 1) * 1.18),
    )

    # Step-filled band with a darker same-hue edge (IMF histogram style).
    hist_color = _OKABE_ITO[1]
    step_x = Float64[]
    step_top = Float64[]
    for j in 1:nbins
        push!(step_x, edges[j])
        push!(step_top, counts[j])
        push!(step_x, edges[j + 1])
        push!(step_top, counts[j])
    end
    band!(ax, step_x, zeros(length(step_x)), step_top; color = (hist_color, 0.55))
    lines!(ax, step_x, step_top; color = _band_edge(hist_color), linewidth = 1.5)

    # Main-sequence turnoff boundary — labelled reference line.
    if hi > 1.0
        vlines!(ax, [1.0]; color = :gray50, linestyle = :dash, linewidth = 1.0)
        # The label goes on whichever side of the line has room; anchoring it
        # always to the right clips it when the axis ends just past t = T_MS.
        x_rel = 1.0 / hi
        to_right = x_rel < 0.8
        text!(
            ax,
            to_right ? x_rel + 0.01 : x_rel - 0.01,
            0.80;
            text = L"t = T_\mathrm{MS}",
            space = :relative,
            align = (to_right ? :left : :right, :top),
            fontsize = 14,
            color = :gray30,
        )
    end

    t_val = _fmt_latex_sig3(sev.time_myr)
    _annotate!(ax, latexstring("t = $(t_val)\\;\\mathrm{Myr}"))

    n_over = count(>(1.0), frac)
    if n_over > 0
        pct_str = @sprintf("%.2g", 100 * n_over / n_ms)
        _annotate!(
            ax,
            latexstring("$(n_over)\\;($(pct_str)\\,\\%)\\;\\mathrm{past}\\;T_\\mathrm{MS}");
            corner = :tr,
            color = _band_edge(hist_color),
        )
    end

    _save_fig(cfg, filename, fig)
    return nothing
end

"""
    plot_core_mass(sevs::Vector{StellarEvolutionSnapshot}, cfg::VisualizationConfig;
                   filename = "core_mass_growth")

Core mass MC [M☉] against total stellar mass [M☉] for the evolved stars
(K* ≥ 2) of the last snapshot, log–log, colour/marker-coded by stellar
class as in the HR figures, with the epoch annotated and the MC = M identity guide marked.
Records with non-finite or non-positive core or total mass are skipped
(SSE fields are NaN in pre-v2026.07 data).

Warns and returns `nothing` when there are no evolved stars with a valid
core mass — early snapshots are MS-only; otherwise saves the figure and
returns `nothing`.
"""
function Nbody6Dynamics.plot_core_mass(
    sevs::Vector{StellarEvolutionSnapshot},
    cfg::VisualizationConfig;
    filename::AbstractString = "core_mass_growth",
)
    isempty(sevs) && (@warn "No stellar evolution snapshots to plot"; return nothing)

    sev = sevs[end]
    evolved = [
        r for r in sev.records if r.stellar_type ≥ 2 &&
            isfinite(r.mass_core) &&
            r.mass_core > 0 &&
            isfinite(r.mass_solar) &&
            r.mass_solar > 0
    ]
    isempty(evolved) && (
        @warn "No evolved stars (K* ≥ 2) with positive core mass — skipping core-mass plot";
        return nothing
    )

    m = [r.mass_solar for r in evolved]
    mc = [r.mass_core for r in evolved]
    m_lo, m_hi = extrema(m)
    mc_lo, mc_hi = extrema(mc)

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = L"M\;[\mathrm{M}_\odot]",
        ylabel = L"M_\mathrm{c}\;[\mathrm{M}_\odot]",
        xscale = log10,
        yscale = log10,
        xticks = _log_ticks(m_lo, m_hi),
        yticks = _log_ticks(mc_lo, mc_hi),
    )

    # MC = M identity guide (core mass cannot exceed the total mass;
    # stripped remnants sit on the line).
    g_lo = min(m_lo, mc_lo) * 0.8
    g_hi = max(m_hi, mc_hi) * 1.2
    lines!(ax, [g_lo, g_hi], [g_lo, g_hi]; color = :gray50, linestyle = :dash, linewidth = 1.0)
    _annotate!(ax, L"M_\mathrm{c} = M"; corner = :br, fontsize = 14, color = :gray30)

    # Classes, colours and markers of the HR figures
    classes_present = sort(unique(stellar_class_index(r.stellar_type) for r in evolved))
    ms = _marker_size(cfg, length(evolved))
    for k in classes_present
        sub = [r for r in evolved if stellar_class_index(r.stellar_type) == k]
        style = _hr_style(k)
        scatter!(
            ax,
            [r.mass_solar for r in sub],
            [r.mass_core for r in sub];
            color = style.color,
            marker = style.marker,
            markersize = ms,
            strokecolor = _band_edge(style.color),
            strokewidth = 1,
            label = STELLAR_CLASSES[k].label,
        )
    end
    length(classes_present) ≥ 2 && _top_legend!(fig, ax; nbanks = cld(length(classes_present), 4))

    # Points accumulate along/below the identity diagonal — the upper-left
    # triangle (MC > M) is empty by construction, so the epoch goes there.
    t_val = _fmt_latex_sig3(sev.time_myr)
    _annotate!(ax, latexstring("t = $(t_val)\\;\\mathrm{Myr}"))

    _save_fig(cfg, filename, fig)
    return nothing
end
