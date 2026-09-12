# =============================================================================
# HR diagram (Hertzsprung-Russell) from stellar evolution data
# =============================================================================

# Unphysical values used by SSE/BSE as placeholders for undefined luminosity/
# temperature (e.g. massless remnants, stars past the end of the grid). Points
# at these values would otherwise force the axes to extend into empty regions
# and visually compress the true stellar distribution.
const _HR_MIN_LOG_L = -5.0
const _HR_MIN_LOG_TEFF = 3.0

"""
    _hr_valid_records(records)

Return the subset of stellar records whose `log_luminosity` and `log_teff` are
within physically meaningful ranges for an HR diagram. Placeholders such as
`log_L = -10` are discarded.
"""
@inline function _hr_valid_records(records)
    return [r for r in records if r.log_luminosity > _HR_MIN_LOG_L && r.log_teff > _HR_MIN_LOG_TEFF]
end

# Okabe–Ito colour + marker encoding for stellar types K* ∈ 0:15: the eight
# palette colours cover one cycle (K* 0–7); the second cycle (K* ≥ 8, white
# dwarfs and compact objects) repeats the colours with a distinct marker so
# the two cycles remain separable (and survive grayscale).
const _HR_COLORS =
    Dict{Int,Makie.RGBAf}(k => _OKABE_ITO[mod1(k + 1, length(_OKABE_ITO))] for k in 0:15)

"""Colour for stellar type `kt` from the Okabe–Ito cycle (grey fallback)."""
_hr_color(kt::Integer) = get(_HR_COLORS, Int(kt), Makie.RGBAf(0.5, 0.5, 0.5, 1))

"""Marker class for stellar type `kt`: circles for K* < 8, triangles for K* ≥ 8."""
_hr_marker(kt::Integer)::Symbol = Int(kt) < 8 ? :circle : :utriangle

"""
    _hr_limits(valid_sets) -> ((t_lo, t_hi), (l_lo, l_hi)) or nothing

Axis limits for HR diagrams: extrema of the valid records' (log Teff, log L)
across one or more record subsets, with a 6 % data margin on each side and
a minimum span of `_HR_MIN_SPAN_DEX` per axis, so a population of identical
stars (equal-mass bodies, no evolution) still yields a finite axis.
Returns `nothing` when no record survives the validity filter.
"""
function _hr_limits(valid_sets::Vector{<:Vector})
    all_teff = reduce(vcat, [[r.log_teff for r in v] for v in valid_sets]; init = Float64[])
    all_lum = reduce(vcat, [[r.log_luminosity for r in v] for v in valid_sets]; init = Float64[])
    isempty(all_teff) && return nothing
    return (_padded_range(extrema(all_teff)...), _padded_range(extrema(all_lum)...))
end

const _HR_MIN_SPAN_DEX = 0.2

"""`(lo, hi)` widened by 6 % on each side, or to `_HR_MIN_SPAN_DEX` about the midpoint when narrower."""
function _padded_range(lo::Real, hi::Real)
    span = hi - lo
    if span < _HR_MIN_SPAN_DEX
        mid = 0.5 * (lo + hi)
        return (mid - 0.5 * _HR_MIN_SPAN_DEX, mid + 0.5 * _HR_MIN_SPAN_DEX)
    end
    return (lo - 0.06 * span, hi + 0.06 * span)
end

"""
    plot_hr(sev::StellarEvolutionSnapshot, cfg::VisualizationConfig;
            filename = "hr_diagram") -> String

Plot a Hertzsprung-Russell diagram (log Teff vs log L) from a single stellar
evolution snapshot, coloured by stellar type K*.

Returns the output file path.
"""
function Nbody6Dynamics.plot_hr(
    sev::StellarEvolutionSnapshot,
    cfg::VisualizationConfig;
    filename::AbstractString = "hr_diagram",
)::String
    isempty(sev.records) && error("No stellar records to plot")

    valid = _hr_valid_records(sev.records)
    isempty(valid) && error("No stellar records pass HR validity filter")

    t_val = @sprintf("%.3g", sev.time_myr)

    # Data ranges (with margin) for tick placement
    tlims, llims = _hr_limits([valid])

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = L"\log_{10}(T_\mathrm{eff} \, / \, \mathrm{K})",
        ylabel = L"\log_{10}(L \, / \, L_\odot)",
        xreversed = true,   # hot → cool from left to right
        xticks = _logval_ticks(tlims...),
        yticks = _logval_ticks(llims...),
        xgridvisible = false,
        ygridvisible = false,
    )
    # sev.time_myr is in Myr, not NB units.  Top-right corner: the sequence
    # enters at top-left, so the upper-right above the ridge line is empty.
    _annotate!(ax, latexstring("t = $(t_val)\\;\\mathrm{Myr}"); corner = :tr)

    # Group by stellar type for legend
    types_present = sort(unique(r.stellar_type for r in valid))

    for kt in types_present
        mask = [r for r in valid if r.stellar_type == kt]
        teff = [r.log_teff for r in mask]
        lum = [r.log_luminosity for r in mask]
        label = get(STELLAR_TYPE_LABELS, Int(kt), "K*=$kt")
        scatter!(
            ax,
            teff,
            lum;
            color = _hr_color(kt),
            marker = _hr_marker(kt),
            markersize = 14,
            label = label,
        )
    end

    if 2 ≤ length(types_present) ≤ 12
        _top_legend!(fig, ax; nbanks = 2)
    end

    return _save_fig(cfg, filename, fig)
end

"""
    plot_hr_evolution(sevs::Vector{StellarEvolutionSnapshot},
                      cfg::VisualizationConfig;
                      filename = "hr_evolution",
                      max_panels = 6) -> String

Plot an HR diagram panel grid showing evolution over multiple epochs.
Selects up to `max_panels` snapshots spaced evenly in time.

Returns the output file path.
"""
function Nbody6Dynamics.plot_hr_evolution(
    sevs::Vector{StellarEvolutionSnapshot},
    cfg::VisualizationConfig;
    filename::AbstractString = "hr_evolution",
    max_panels::Int = 6,
)::String
    isempty(sevs) && error("No stellar evolution snapshots to plot")

    # Select snapshots evenly spaced in time
    n = length(sevs)
    indices = if n ≤ max_panels
        collect(1:n)
    else
        unique(round.(Int, range(1, n; length = max_panels)))
    end

    ncols = min(length(indices), 3)
    nrows = cld(length(indices), ncols)

    # One column wide; inner tick labels are hidden, so the compact gap applies.
    gap = _multipanel_gap(ncols; inner_ticks = false)
    fig = Figure(; size = _fig_multipanel(cfg, nrows, ncols; inner_ticks = false))
    target_ticks = ncols > 1 ? 3 : 5

    # Precompute the valid-record subset for each panel once
    valid_per_panel = [_hr_valid_records(sevs[i].records) for i in indices]

    # Consistent axis limits across panels (using valid records only);
    # normal order for limits — xreversed handles the flip
    lims = _hr_limits(valid_per_panel)
    lims === nothing && error("No valid HR records across any panel")
    tlims, llims = lims
    xtk_hr = _logval_ticks(tlims[1], tlims[2]; target_n = target_ticks)
    ytk_hr = _logval_ticks(llims[1], llims[2]; target_n = target_ticks)
    # Data-free band above the data, where the time annotation sits
    llims_plot = (llims[1], llims[2] + _MONTAGE_BAND_FRAC * (llims[2] - llims[1]))
    marker_scale = max(_multipanel_scale(cfg, ncols; inner_ticks = false), 0.6)

    for (panel_idx, si) in enumerate(indices)
        row = cld(panel_idx, ncols)
        col = mod1(panel_idx, ncols)
        sev = sevs[si]

        # Only show axis labels on border panels
        show_xlab = row == nrows
        show_ylab = col == 1

        t_val = @sprintf("%.3g", sev.time_myr)
        ax = Axis(
            fig[row, col];
            xlabel = show_xlab ? L"\log_{10}(T_\mathrm{eff} \, / \, \mathrm{K})" : "",
            ylabel = show_ylab ? L"\log_{10}(L \, / \, L_\odot)" : "",
            xlabelsize = 22,
            ylabelsize = 22,
            xticklabelsize = 18,
            yticklabelsize = 18,
            xreversed = true,
            limits = (tlims..., llims_plot...),
            xticks = xtk_hr,
            yticks = ytk_hr,
            xticklabelsvisible = show_xlab,
            yticklabelsvisible = show_ylab,
            xgridvisible = false,
            ygridvisible = false,
        )
        # sev.time_myr is in Myr, not NB units.  Top-right in-axis corner is
        # empty on an HR diagram (the sequence enters at top-left).
        _annotate!(ax, latexstring("t = $(t_val)\\;\\mathrm{Myr}"); corner = :tr)

        v = valid_per_panel[panel_idx]
        isempty(v) || scatter!(
            ax,
            [r.log_teff for r in v],
            [r.log_luminosity for r in v];
            color = [_hr_color(r.stellar_type) for r in v],
            marker = [_hr_marker(r.stellar_type) for r in v],
            markersize = 14 * marker_scale,
            strokewidth = 0,
        )
    end

    colgap!(fig.layout, gap)
    rowgap!(fig.layout, gap)

    return _save_fig(cfg, filename, fig)
end
