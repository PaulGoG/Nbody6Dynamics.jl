# =============================================================================
# Tidal disruptions at the static black hole (KZ(44))
# =============================================================================

"""Colours of the event kinds: full disruptions vermillion, partial blue, plunges black (Okabe–Ito)."""
const _TDE_KIND_STYLE = (
    full = (label = "Full", colour = _OKABE_ITO[6], marker = :circle),
    partial = (label = "Partial", colour = _OKABE_ITO[1], marker = :utriangle),
    plunge = (label = "Plunge", colour = _OKABE_ITO[8], marker = :rect),
)

"""
    plot_tde_events(events::AbstractVector{TDEEvent}, cfg::VisualizationConfig;
                    filename = "tde_events") -> String

Events at the static black hole in two stacked panels sharing the time axis.
Upper panel: the cumulative number of full disruptions, partial disruptions
and plunges against time (steps, one colour per kind, the counts in the
legend, or stated in the panel for a single kind). Lower panel, logarithmic: the penetration factor β = r_t/q of
every event (marker shape and colour by kind) with a guide line at β = 1.
The mass removed from the cluster is stated in the upper panel.

Throws an error when no event has a finite `time_myr` (read the events with
`tstar_myr`).
"""
@publication function Nbody6Dynamics.plot_tde_events(
    events::AbstractVector{TDEEvent},
    cfg::VisualizationConfig;
    filename::AbstractString = "tde_events",
)::String
    evs = sort!([e for e in events if isfinite(e.time_myr)]; by = e -> e.time_myr)
    isempty(evs) && error("No TDE event with a physical time to plot")
    t_all = [e.time_myr for e in evs]
    tmin = 0.0
    tmax = t_all[end] > 0 ? t_all[end] * 1.02 : 1.0
    tticks = _time_ticks(tmin, tmax)

    kinds = [k for k in (:full, :partial, :plunge) if any(e -> e.kind === k, evs)]
    beta = [e.beta for e in evs if e.beta > 0]
    b_lo, b_hi = isempty(beta) ? (0.5, 2.0) : _imbh_log_limits(vcat(beta, 1.0))
    n_max = maximum(k -> count(e -> e.kind === k, evs), kinds)

    fig = Figure(; size = _fig_multipanel(cfg, 2, 1))
    ax1 = Axis(fig[1, 1]; ylabel = "Cumulative events", xticks = tticks, xticklabelsvisible = false)
    ax2 = Axis(
        fig[2, 1];
        xlabel = L"t\;[\mathrm{Myr}]",
        ylabel = L"\beta = r_\mathrm{t}/q",
        yscale = log10,
        yticks = _log_ticks(b_lo, b_hi),
        xticks = tticks,
    )
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    xlims!(ax2, tmin, tmax)
    ylims!(ax1, _pad_limits(0.0, Float64(n_max); floor = 0.0)...)
    ylims!(ax2, b_lo, b_hi)

    hlines!(ax2, [1.0]; color = :gray45, linestyle = :dash, linewidth = _STYLE.guide)
    elements = _LegendElement[]
    labels = AbstractString[]
    for k in kinds
        st = _TDE_KIND_STYLE[k]
        tk = [e.time_myr for e in evs if e.kind === k]
        n = length(tk)
        stairs!(
            ax1,
            vcat(tmin, tk, tmax),
            vcat(0, 1:n, n);
            step = :post,
            color = st.colour,
            linewidth = _STYLE.data,
        )
        bk = [(e.time_myr, e.beta) for e in evs if e.kind === k && e.beta > 0]
        isempty(bk) || scatter!(
            ax2,
            first.(bk),
            last.(bk);
            color = (st.colour, 0.75),
            marker = st.marker,
            markersize = _STYLE.marker,
            strokecolor = _band_edge(st.colour),
            strokewidth = _STYLE.marker_stroke,
        )
        push!(elements, LineElement(; color = st.colour, linewidth = _STYLE.data))
        push!(labels, "$(st.label) ($n)")
    end

    # Mass removed in the upper left, where the cumulative curves start at zero
    s = tde_summary(evs)
    single = length(elements) < 2
    _annotate!(
        ax1,
        latexstring("\$\\Delta M = $(_fmt_latex_sig(s.dm_msun, 3))\\;\\mathrm{M}_\\odot\$ removed");
        corner = :tl,
        dy = single ? 0.1 : 0.0,
    )

    # One kind: no legend, its count in the upper panel
    if length(elements) ≥ 2
        _top_legend!(fig, elements, labels)
    else
        _annotate!(ax1, labels[1]; corner = :tl, color = _TDE_KIND_STYLE[kinds[1]].colour)
    end

    return _save_fig(cfg, filename, fig)
end
