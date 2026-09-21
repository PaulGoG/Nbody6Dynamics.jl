# =============================================================================
# Escaper analysis plots from esc.11 data
# =============================================================================

# Compact-remnant threshold on the Hurley K* scale: white dwarfs and heavier
# (K* ≥ 10) versus luminous stars (K* < 10).
const _ESC_REMNANT_KSTAR_MIN = 10

# Two-class encoding shared by all escaper figures: colour paired with
# marker so the distinction survives grayscale.
const _ESC_LUMINOUS_COLOR = _OKABE_ITO[1]   # blue
const _ESC_REMNANT_COLOR = _OKABE_ITO[6]    # vermillion
const _ESC_LUMINOUS_LABEL = L"\mathrm{Luminous}\;(K^* < 10)"
const _ESC_REMNANT_LABEL = L"\mathrm{Compact\;remnant}\;(K^* \geq 10)"

"""Whether an escaper is a compact remnant (K* ≥ 10) on the Hurley scale."""
_esc_is_remnant(r::EscaperRecord) = r.stellar_type ≥ _ESC_REMNANT_KSTAR_MIN

"""
    _scatter_escaper_classes!(ax, records, xf, yf, ms) -> (elements, labels)

Scatter `records` split into the two stellar classes (luminous vs compact
remnant) with the shared colour + marker encoding, using `xf`/`yf` to
extract the plotted quantities. Returns the legend elements and labels for
the classes actually present.
"""
function _scatter_escaper_classes!(ax, records, xf, yf, ms)
    legend_elems = Any[]
    legend_labels = LaTeXString[]
    for (sel, color, marker, label) in (
        (
            [r for r in records if !_esc_is_remnant(r)],
            _ESC_LUMINOUS_COLOR,
            :circle,
            _ESC_LUMINOUS_LABEL,
        ),
        (
            [r for r in records if _esc_is_remnant(r)],
            _ESC_REMNANT_COLOR,
            :utriangle,
            _ESC_REMNANT_LABEL,
        ),
    )
        isempty(sel) && continue
        p = scatter!(
            ax,
            [xf(r) for r in sel],
            [yf(r) for r in sel];
            color = color,
            marker = marker,
            markersize = ms,
            strokewidth = 0,
        )
        push!(legend_elems, p)
        push!(legend_labels, label)
    end
    return legend_elems, legend_labels
end

"""
    plot_escapers(escapers::Vector{EscaperRecord}, cfg::VisualizationConfig;
                  filename = "escapers")

Two-panel escaper summary sharing the escape-time axis [Myr]: (top)
cumulative escaped mass [M☉] as a step curve with the total escaper count
and mass annotated; (bottom) escape velocity [km/s] on a log axis,
distinguishing luminous stars (K* < 10) from compact remnants (K* ≥ 10).

Escaper quantities are intrinsically physical (esc.11 tokens 6–13), so the
axes are Myr / M☉ / km s⁻¹ regardless of `cfg.units`.  Warns and returns
`nothing` for empty input; otherwise saves the figure and returns `nothing`.
"""
@publication function Nbody6Dynamics.plot_escapers(
    escapers::Vector{EscaperRecord},
    cfg::VisualizationConfig;
    filename::AbstractString = "escapers",
)
    isempty(escapers) && (@warn "No escaper records to plot"; return nothing)

    esc = sort(escapers; by = r -> r.time_myr)
    t = [r.time_myr for r in esc]
    m_cum = cumsum([r.mass_solar for r in esc])
    n_esc = length(esc)

    fig = Figure(; size = _fig_two_panel(cfg))
    # Fewer ticks than the _time_ticks default: escape times often span a
    # narrow window and dense decimal labels collide at column width.
    ttk = _nice_ticks(first(t), last(t); target_n = 6)
    esc_color = _SEMANTIC_COLORS[:escapers]

    # --- Panel 1: cumulative escaped mass (step curve) ---
    ax1 = Axis(
        fig[1, 1];
        ylabel = L"M_\mathrm{esc}\;[\mathrm{M}_\odot]",
        xticklabelsvisible = false,
        xticks = ttk,
    )
    stairs!(ax1, t, m_cum; color = esc_color, step = :post)
    # The curve rises towards the lower right corner's diagonal — the
    # top-left corner stays free for the totals annotation.
    m_str = _fmt_latex_sig3(m_cum[end])
    _annotate!(
        ax1,
        latexstring(
            "N_\\mathrm{esc} = $(n_esc),\\;M_\\mathrm{esc} = $(m_str)\\;\\mathrm{M}_\\odot",
        );
        color = esc_color,
    )

    # --- Panel 2: escape velocity vs time (log y), two stellar classes ---
    # Non-positive or non-finite velocities are invalid on the log axis.
    vok = [r for r in esc if isfinite(r.velocity_kms) && r.velocity_kms > 0]
    v_all = [r.velocity_kms for r in vok]
    ax2 = Axis(
        fig[2, 1];
        xlabel = L"t\;[\mathrm{Myr}]",
        ylabel = L"v_\mathrm{esc}\;[\mathrm{km\,s^{-1}}]",
        xticks = ttk,
        yscale = log10,
        yticks = isempty(v_all) ? Makie.automatic : _log_ticks(extrema(v_all)...),
    )

    ms = _marker_size(cfg, length(vok))
    legend_elems, legend_labels =
        _scatter_escaper_classes!(ax2, vok, r -> r.time_myr, r -> r.velocity_kms, ms)
    isempty(vok) && _no_data_note!(ax2, L"\mathrm{no\;valid\;escape\;velocities}")
    # Single-entry legends are suppressed per the house standard.
    length(legend_labels) ≥ 2 && _top_legend!(fig, legend_elems, legend_labels)

    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)

    _save_fig(cfg, filename, fig)
    return nothing
end

"""
    plot_escape_anisotropy(escapers::Vector{EscaperRecord}, cfg::VisualizationConfig;
                           filename = "escape_anisotropy")

Sky projection of the escape directions from esc.11: azimuth φ [deg]
(from the x-axis, native range [0°, 360°]) against elevation θ [deg]
(from the xy-plane, [-90°, 90°]), on plain Cartesian axes with equal
data margins.  Luminous stars (K* < 10) and compact remnants (K* ≥ 10)
are distinguished by colour and marker; the escaper count is annotated.

Warns and returns `nothing` when no record carries finite direction
angles; otherwise saves the figure and returns `nothing`.
"""
@publication function Nbody6Dynamics.plot_escape_anisotropy(
    escapers::Vector{EscaperRecord},
    cfg::VisualizationConfig;
    filename::AbstractString = "escape_anisotropy",
)
    isempty(escapers) && (@warn "No escaper records to plot"; return nothing)

    valid = [r for r in escapers if isfinite(r.phi_deg) && isfinite(r.theta_deg)]
    isempty(valid) && (@warn "No escaper records with finite direction angles"; return nothing)

    # escape.F writes φ ∈ [0°, 360°] and θ ∈ [-90°, 90°]; keep the native
    # ranges and frame the data with equal fractional margins on both axes.
    phi = [r.phi_deg for r in valid]
    theta = [r.theta_deg for r in valid]
    phi_lo, phi_hi = extrema(phi)
    th_lo, th_hi = extrema(theta)
    dphi = max(phi_hi - phi_lo, 1.0) * 0.05
    dth = max(th_hi - th_lo, 1.0) * 0.05

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = L"\phi\;[\mathrm{deg}]",
        ylabel = L"\theta\;[\mathrm{deg}]",
        limits = (phi_lo - dphi, phi_hi + dphi, th_lo - dth, th_hi + dth),
        xticks = _nice_ticks(phi_lo - dphi, phi_hi + dphi),
        yticks = _nice_ticks(th_lo - dth, th_hi + dth; target_n = 6),
    )

    ms = _marker_size(cfg, length(valid))
    legend_elems, legend_labels =
        _scatter_escaper_classes!(ax, valid, r -> r.phi_deg, r -> r.theta_deg, ms)
    length(legend_labels) ≥ 2 && _top_legend!(fig, legend_elems, legend_labels)

    _annotate!(ax, latexstring("N_\\mathrm{esc} = $(length(valid))"))

    _save_fig(cfg, filename, fig)
    return nothing
end
