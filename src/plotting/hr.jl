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
    plot_hr(sev::StellarEvolutionSnapshot, cfg::VisualizationConfig;
            filename = "hr_diagram") -> String

Plot a Hertzsprung-Russell diagram (log Teff vs log L) from a single stellar
evolution snapshot, coloured by stellar type K*.

Returns the output file path.
"""
function plot_hr(
    sev::StellarEvolutionSnapshot,
    cfg::VisualizationConfig;
    filename::AbstractString = "hr_diagram",
)::String
    isempty(sev.records) && error("No stellar records to plot")

    valid = _hr_valid_records(sev.records)
    isempty(valid) && error("No stellar records pass HR validity filter")

    t_val = @sprintf("%.3g", sev.time_myr)

    # Compute data ranges for tick placement
    all_teff = [r.log_teff for r in valid]
    all_lum = [r.log_luminosity for r in valid]
    teff_min, teff_max = extrema(all_teff)
    lum_min, lum_max = extrema(all_lum)
    dt = (teff_max - teff_min) * 0.06
    dl = (lum_max - lum_min) * 0.06

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = L"\log_{10}(T_\mathrm{eff} \, / \, \mathrm{K})",
        ylabel = L"\log_{10}(L \, / \, L_\odot)",
        xreversed = true,   # hot → cool from left to right
        xticks = _logval_ticks(teff_min - dt, teff_max + dt),
        yticks = _logval_ticks(lum_min - dl, lum_max + dl),
        xgridvisible = false,
        ygridvisible = false,
    )
    # sev.time_myr is in Myr, not NB units.  Top-right corner: the sequence
    # enters at top-left, so the upper-right above the ridge line is empty.
    text!(
        ax,
        0.96,
        0.96;
        text = latexstring("t = $(t_val)\\;\\mathrm{Myr}"),
        space = :relative,
        align = (:right, :top),
        fontsize = 16,
    )

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
function plot_hr_evolution(
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

    fig = Figure(; size = _fig_multipanel(cfg, nrows, ncols))

    # Precompute the valid-record subset for each panel once
    valid_per_panel = [_hr_valid_records(sevs[i].records) for i in indices]

    # Consistent axis limits across panels (using valid records only)
    all_teff = reduce(vcat, [[r.log_teff for r in v] for v in valid_per_panel]; init = Float64[])
    all_lum =
        reduce(vcat, [[r.log_luminosity for r in v] for v in valid_per_panel]; init = Float64[])
    isempty(all_teff) && error("No valid HR records across any panel")
    tmin, tmax = extrema(all_teff)
    lmin, lmax = extrema(all_lum)
    dt = (tmax - tmin) * 0.06
    dl = (lmax - lmin) * 0.06
    # Normal order for limits; xreversed handles the flip
    tlims = (tmin - dt, tmax + dt)
    llims = (lmin - dl, lmax + dl)
    xtk_hr = _logval_ticks(tlims[1], tlims[2]; target_n = 5)
    ytk_hr = _logval_ticks(llims[1], llims[2]; target_n = 5)

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
            limits = (tlims..., llims...),
            xticks = xtk_hr,
            yticks = ytk_hr,
            xticklabelsvisible = show_xlab,
            yticklabelsvisible = show_ylab,
            xgridvisible = false,
            ygridvisible = false,
        )
        # sev.time_myr is in Myr, not NB units.  Top-right in-axis corner is
        # empty on an HR diagram (the sequence enters at top-left).
        text!(
            ax,
            0.96,
            0.96;
            text = latexstring("t = $(t_val)\\;\\mathrm{Myr}"),
            space = :relative,
            align = (:right, :top),
            fontsize = 16,
        )

        v = valid_per_panel[panel_idx]
        isempty(v) || scatter!(
            ax,
            [r.log_teff for r in v],
            [r.log_luminosity for r in v];
            color = [_hr_color(r.stellar_type) for r in v],
            marker = [_hr_marker(r.stellar_type) for r in v],
            markersize = 14,
            strokewidth = 0,
        )
    end

    colgap!(fig.layout, _MULTIPANEL_HGAP)
    rowgap!(fig.layout, _MULTIPANEL_VGAP)

    return _save_fig(cfg, filename, fig)
end
