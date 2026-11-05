# =============================================================================
# Calibration report: time-resolved paired differences and the forest plot of the
# window-mean shifts per arm, against the equivalence margin
# =============================================================================

"""Arm colours of a calibration report: Okabe–Ito hues other than the verdict colours
(blue, reddish purple, orange, sky blue, black, yellow), a viridis sampling beyond six."""
const _ARM_PALETTE_INDICES = (1, 4, 2, 5, 8, 7)

"""
    _arm_colors(rows) -> Dict{String,RGBAf}

One colour per arm of the table `rows` (arms in name order), shared by every figure
drawn from it.
"""
function _arm_colors(rows::AbstractVector)
    arms = sort!(unique(String(r.arm) for r in rows))
    n = length(arms)
    colors = if n ≤ length(_ARM_PALETTE_INDICES)
        [_OKABE_ITO[_ARM_PALETTE_INDICES[i]] for i in 1:n]
    else
        cmap = Makie.to_colormap(:viridis)
        [cmap[round(Int, 1 + (i - 1) / max(n - 1, 1) * (length(cmap) - 1))] for i in 1:n]
    end
    return Dict{String,Makie.RGBAf}(a => Makie.to_color(c) for (a, c) in zip(arms, colors))
end

"""Marker and colour of each verdict: shape and hue both carry it, so it survives grey scale."""
_verdict_style(verdict::AbstractString) =
    verdict == "equivalent" ? (marker = :circle, color = _OKABE_ITO[3]) :
    verdict == "different" ? (marker = :diamond, color = _OKABE_ITO[6]) :
    (marker = :rect, color = Makie.to_color(:gray55))

"""Verdicts in legend order with their labels."""
const _VERDICTS =
    (("equivalent", "Equivalent"), ("different", "Different"), ("inconclusive", "Inconclusive"))

"""LaTeX symbol of an observable of the paired analysis and whether it is a radius [pc]."""
function _observable_symbol(observable::AbstractString)
    observable == "r10" && return ("r_{10\\,\\%}", true)
    observable == "r50" && return ("r_{50\\,\\%}", true)
    observable == "r90" && return ("r_{90\\,\\%}", true)
    observable == "energy" && return ("|\\Delta E / E|", false)
    observable == "n_stars" && return ("N", false)
    observable == "n_pairs" && return ("N_\\mathrm{pairs}", false)
    throw(ArgumentError("unknown observable \"$observable\""))
end

"""
    _power_of_ten(values) -> Int

Power of ten factored out of an axis whose largest finite magnitude among `values` lies
outside 10⁻²–10⁴ (`0` inside, for an empty or all-zero set).
"""
function _power_of_ten(values)::Int
    v = maximum((abs(x) for x in values if isfinite(x)); init = 0.0)
    v == 0 && return 0
    e = floor(Int, log10(v))
    return -2 ≤ e ≤ 3 ? 0 : e
end

"""Axis label of the difference of an observable (`window`: its window mean ⟨·⟩_W), with the unit and a factored power of ten `e`."""
function _difference_label(observable::AbstractString, e::Int; window::Bool = false)
    sym, radius = _observable_symbol(observable)
    parts = String[]
    e == 0 || push!(parts, "10^{$(e)}")
    radius && push!(parts, "\\mathrm{pc}")
    unit = isempty(parts) ? "" : " \\; [" * join(parts, "\\,") * "]"
    body = window ? "\\langle \\Delta $(sym) \\rangle_W" : "\\Delta $(sym)"
    return latexstring(body * unit)
end

"""Whether an observable is a count (`n_stars`, `n_pairs`), whose axes take integer ticks."""
_count_observable(observable::AbstractString) = observable in ("n_stars", "n_pairs")

"""The grid point of a figure: `grid` when given, the only grid point of `rows` otherwise."""
function _calibration_grid(rows::AbstractVector, grid)
    grids = sort!(unique(String(r.grid) for r in rows))
    if grid === nothing
        length(grids) == 1 && return only(grids)
        throw(
            ArgumentError(
                "the table spans several grid points; pass grid = one of " *
                join(repr.(grids), ", "),
            ),
        )
    end
    String(grid) in grids || throw(
        ArgumentError("grid point \"$grid\" not in the table; one of " * join(repr.(grids), ", ")),
    )
    return String(grid)
end

"""
    _margin_text(row) -> Union{String,LaTeXString}

The margin of a table row for an annotation: `ε = κ σ_ch = value` for a chaos margin,
`ε = value` (marked relative for a relative one), `No margin` when it is not finite.
"""
function _margin_text(row)
    isfinite(row.margin) || return "No margin"
    _, radius = _observable_symbol(row.observable)
    unit = radius ? "\\,\\mathrm{pc}" : ""
    value = _fmt_latex_sig(row.margin, 3) * unit
    if row.margin_kind == "chaos" && isfinite(row.chaos_floor) && row.chaos_floor > 0
        κ = _fmt_latex_sig(row.margin / row.chaos_floor, 2)
        return latexstring("\\varepsilon = $(κ)\\,\\sigma_\\mathrm{ch} = $(value)")
    end
    suffix = row.margin_kind == "relative" ? "\\; (\\mathrm{relative})" : ""
    return latexstring("\\varepsilon = $(value)$(suffix)")
end

"""Line spacing of stacked annotations, in relative axis units."""
const _ANNOTATION_STEP = 0.075

"""Fraction of the axis height kept free above the data for `n` stacked annotation lines."""
_strip_fraction(n::Integer) = 0.05 + _ANNOTATION_STEP * n

"""
    _annotate_lines!(ax, lines; x = 0.04)

Stacked single-line annotations from the top of `ax`, left-aligned at the relative
abscissa `x`, one text per line so that plain text and LaTeX strings mix; the first line
is the uppermost. The strip they take is [`_strip_fraction`](@ref) of the axis height.
"""
function _annotate_lines!(ax, lines::AbstractVector; x::Real = 0.04)
    for (i, line) in enumerate(lines)
        text!(
            ax,
            x,
            0.96 - _ANNOTATION_STEP * (i - 1);
            text = line,
            space = :relative,
            align = (:left, :top),
            fontsize = _ANNOTATION_FONTSIZE,
            color = :black,
        )
    end
    return nothing
end

"""
    _limits_with_strip(a, b, frac; pad = 0.06) -> (lo, hi)

Axis limits of the data range `[a, b]` with a margin of `pad` of its width on both
sides and a strip of `frac` of the axis height free above the data, for the annotation.
"""
function _limits_with_strip(a::Real, b::Real, frac::Real; pad::Real = 0.06)
    w = b > a ? b - a : max(abs(a), 1.0)
    lo = a - pad * w
    return (lo, lo + (b + pad * w - lo) / (1 - frac))
end

"""
    plot_paired_differences(series, rows, cfg::VisualizationConfig; observable = "r50",
                            grid = nothing, reference = "base", t_settle_myr = 0.0,
                            filename = "paired_difference_<observable>") -> String

Time-resolved mean paired difference, arm minus `reference`, of `observable` for every
arm ([`paired_difference_series`](@ref) gives `series`), the mean as a line and its
standard error as a band in the arm's colour, over the equivalence margin ±ε of the
table `rows` ([`paired_comparison`](@ref), [`block_comparison`](@ref)) as a grey band
and a dashed zero line. The interval before `t_settle_myr` (outside the window of the
shifts) is shaded. Arms whose difference vanishes at every node are named in the
annotation instead of drawn. `grid` selects the grid point of a table that spans
several. Colours per arm as in [`plot_paired_shifts`](@ref).
"""
@publication function Nbody6Dynamics.plot_paired_differences(
    series::AbstractVector,
    rows::AbstractVector,
    cfg::VisualizationConfig;
    observable::AbstractString = "r50",
    grid::Union{Nothing,AbstractString} = nothing,
    reference::AbstractString = "base",
    t_settle_myr::Real = 0.0,
    filename::AbstractString = "paired_difference_" * String(observable),
)::String
    g = _calibration_grid(rows, grid)
    trows = [r for r in rows if r.observable == observable && r.grid == g]
    isempty(trows) && throw(ArgumentError("no table row of the observable \"$observable\""))
    sel = [s for s in series if s.observable == observable && s.grid == g && !isempty(s.time)]
    colors = _arm_colors(rows)
    margin = first(trows).margin
    identical =
        sort!([String(s.arm) for s in sel if all(iszero, s.mean) && all(x -> !(x > 0), s.se)])
    drawn = [s for s in sel if !(String(s.arm) in identical)]
    hi =
        [abs(x) for s in drawn for x in vcat(s.mean .+ s.se, s.mean .- s.se, s.mean) if isfinite(x)]
    e = _power_of_ten(vcat(hi, isfinite(margin) ? [margin] : Float64[]))
    f = 10.0^(-e)

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(fig[1, 1]; xlabel = L"t \; [\mathrm{Myr}]", ylabel = _difference_label(observable, e))
    t_lo = minimum((first(s.time) for s in sel); init = Inf)
    t_hi = maximum((last(s.time) for s in sel); init = -Inf)
    grey = _OKABE_ITO[8]
    if isfinite(margin) && margin > 0
        hspan!(ax, -margin * f, margin * f; color = (grey, 0.10))
        hlines!(
            ax,
            [-margin * f, margin * f];
            color = (grey, 0.55),
            linestyle = :dash,
            linewidth = _STYLE.guide,
        )
    end
    hlines!(ax, [0.0]; color = grey, linestyle = :dash, linewidth = _STYLE.guide)
    if isfinite(t_lo) && t_settle_myr > t_lo
        vspan!(ax, t_lo, t_settle_myr; color = (:gray50, 0.12))
        vlines!(ax, [t_settle_myr]; color = :gray45, linestyle = :dash, linewidth = _STYLE.guide)
    end
    occ_y = Float64[]
    for s in drawn
        c = colors[String(s.arm)]
        m = s.mean .* f
        if all(isfinite, s.se)
            lo = (s.mean .- s.se) .* f
            up = (s.mean .+ s.se) .* f
            band!(ax, s.time, lo, up; color = (c, _STYLE.band_alpha))
            append!(occ_y, lo, up)
        end
        lines!(ax, s.time, m; color = c, linewidth = _STYLE.data)
        append!(occ_y, m)
    end
    isempty(sel) && _no_data_note!(ax, "No paired $(observable) series")
    if isfinite(t_lo) && t_hi > t_lo
        xlims!(ax, t_lo, t_hi)
        ax.xticks = _time_ticks(t_lo, t_hi)
    end
    ns = unique(s.n for s in sel)
    lines_text = Any[]
    isempty(ns) || push!(
        lines_text,
        (length(ns) == 1 ? "$(only(ns))" : "$(minimum(ns))–$(maximum(ns))") *
        " seeds per arm" *
        (t_settle_myr > 0 ? "; window from $(_fmt_plain_sig(t_settle_myr, 3)) Myr" : ""),
    )
    push!(lines_text, _margin_text(first(trows)))
    isempty(identical) || push!(lines_text, "Identical to $(reference): " * join(identical, ", "))
    # The data and the margin band below, the annotation in a free strip above, starting
    # right of the settling interval so that its boundary never crosses the text.
    ε = isfinite(margin) ? margin * f : 0.0
    y_a = min(minimum(occ_y; init = 0.0), -ε)
    y_b = max(maximum(occ_y; init = 0.0), ε)
    y_a == y_b && ((y_a, y_b) = (-1.0, 1.0))
    y_lims = _limits_with_strip(y_a, y_b, _strip_fraction(length(lines_text)))
    ylims!(ax, y_lims...)
    _count_observable(observable) && e == 0 && (ax.yticks = _integer_ticks(y_lims...; target_n = 5))
    x_text = 0.04
    if isfinite(t_lo) && t_hi > t_lo && t_settle_myr > t_lo
        x_text = max(x_text, (t_settle_myr - t_lo) / (t_hi - t_lo) + 0.02)
    end
    _annotate_lines!(ax, lines_text; x = x_text)

    entries = Vector{Vector{_LegendElement}}()
    labels = Vector{Vector{AbstractString}}()
    titles = String[]
    arms_drawn = sort!([String(s.arm) for s in drawn])
    if !isempty(arms_drawn)
        push!(
            entries,
            _LegendElement[
                LineElement(; color = colors[a], linewidth = _STYLE.data) for a in arms_drawn
            ],
        )
        push!(labels, AbstractString[a for a in arms_drawn])
        push!(titles, "Arm − $(reference):")
    end
    push!(
        entries,
        _LegendElement[
            PolyElement(;
                color = (grey, _STYLE.band_alpha),
                strokecolor = _band_edge(grey),
                strokewidth = _STYLE.envelope,
            ),
            PolyElement(;
                color = (grey, 0.10),
                strokecolor = (grey, 0.55),
                strokewidth = _STYLE.guide,
                linestyle = :dash,
            ),
        ],
    )
    push!(labels, AbstractString["±SE", "±ε"])
    push!(titles, "Bands:")
    Legend(
        fig[0, :],
        entries,
        labels,
        titles;
        orientation = :horizontal,
        nbanks = length(arms_drawn) > 3 ? 2 : 1,
        framevisible = false,
        titleposition = :left,
        tellheight = true,
        padding = (0, 0, 0, 0),
    )
    return _save_fig(cfg, filename, fig)
end

"""
    plot_paired_shifts(rows, cfg::VisualizationConfig; observable = "r50", grid = nothing,
                       reference = "base", α = 0.05, holm = true,
                       filename = "paired_shift_<observable>") -> String

Forest plot of the table `rows` ([`paired_comparison`](@ref),
[`block_comparison`](@ref)) for `observable`: per arm (one row each, in name order) the
window-mean shift against `reference` with its (1 − 2α) interval as a bar in the arm's
colour, the verdict as the marker (circle and green for `equivalent`, diamond and
vermillion for `different`, square and grey for `inconclusive`), the margin ±ε as a grey
band and a dashed zero line. With `holm` the Holm-adjusted p-values of the equivalence
test stand on the right-hand axis. `α` must be the level of the table.
"""
@publication function Nbody6Dynamics.plot_paired_shifts(
    rows::AbstractVector,
    cfg::VisualizationConfig;
    observable::AbstractString = "r50",
    grid::Union{Nothing,AbstractString} = nothing,
    reference::AbstractString = "base",
    α::Real = 0.05,
    holm::Bool = true,
    filename::AbstractString = "paired_shift_" * String(observable),
)::String
    g = _calibration_grid(rows, grid)
    trows = sort!([r for r in rows if r.observable == observable && r.grid == g]; by = r -> r.arm)
    isempty(trows) && throw(ArgumentError("no table row of the observable \"$observable\""))
    colors = _arm_colors(rows)
    k = length(trows)
    y = Float64[k - i + 1 for i in 1:k]               # first arm on top
    margin = first(trows).margin
    vals = vcat(
        [r.shift for r in trows],
        [r.ci_low for r in trows],
        [r.ci_high for r in trows],
        isfinite(margin) ? [margin] : Float64[],
    )
    e = _power_of_ten(vals)
    f = 10.0^(-e)
    span = maximum((abs(v) * f for v in vals if isfinite(v)); init = 0.0)
    span = span > 0 ? span : 1.0

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = _difference_label(observable, e; window = true),
        yticks = (y, [String(r.arm) for r in trows]),
        yticksvisible = false,
        ygridvisible = false,
    )
    grey = _OKABE_ITO[8]
    # Rows from 1 (bottom) to k and a free strip above them for the annotation; the
    # margin band and the guides end at the rows, below the strip.
    y_lims = (0.5, 0.5 + (k - 0.15) / (1 - _strip_fraction(2)))
    y_rows = (y_lims[1], k + 0.35)
    if isfinite(margin) && margin > 0
        poly!(
            ax,
            Rect2f(-margin * f, y_rows[1], 2 * margin * f, y_rows[2] - y_rows[1]);
            color = (grey, 0.10),
        )
        linesegments!(
            ax,
            [
                Point2f(-margin * f, y_rows[1]),
                Point2f(-margin * f, y_rows[2]),
                Point2f(margin * f, y_rows[1]),
                Point2f(margin * f, y_rows[2]),
            ];
            color = (grey, 0.55),
            linestyle = :dash,
            linewidth = _STYLE.guide,
        )
    end
    linesegments!(
        ax,
        [Point2f(0, y_rows[1]), Point2f(0, y_rows[2])];
        color = grey,
        linestyle = :dash,
        linewidth = _STYLE.guide,
    )
    for (yi, r) in zip(y, trows)
        c = colors[String(r.arm)]
        if isfinite(r.ci_low) && isfinite(r.ci_high)
            rangebars!(
                ax,
                [yi],
                [r.ci_low * f],
                [r.ci_high * f];
                direction = :x,
                color = c,
                linewidth = _STYLE.data,
                whiskerwidth = 16,
            )
        end
        if isfinite(r.shift)
            v = _verdict_style(r.verdict)
            scatter!(
                ax,
                [r.shift * f],
                [yi];
                marker = v.marker,
                color = v.color,
                markersize = 1.3 * _STYLE.marker,
                strokecolor = _band_edge(v.color),
                strokewidth = _STYLE.marker_stroke,
            )
        end
    end
    xlims!(ax, -1.15 * span, 1.15 * span)
    _count_observable(observable) &&
        e == 0 &&
        (ax.xticks = _integer_ticks(-1.15 * span, 1.15 * span; target_n = 5))
    ylims!(ax, y_lims...)
    ns = unique(r.n for r in trows)
    level = 100 * (1 - 2α)
    _annotate_lines!(
        ax,
        Any[
            (length(ns) == 1 ? "$(only(ns))" :
             "$(minimum(ns))–$(maximum(ns))") * " seeds; $(_fmt_plain_sig(level, 3)) % intervals",
            _margin_text(first(trows)),
        ],
    )
    if holm
        ax_p = Axis(
            fig[1, 1];
            yaxisposition = :right,
            backgroundcolor = :transparent,
            yticks = (y, [_fmt_plain_sig(r.p_equivalence_holm, 2) for r in trows]),
            ylabel = rich(rich("p"; font = :italic), subscript("TOST"), " (Holm)"),
            yticksvisible = false,
            ygridvisible = false,
            xgridvisible = false,
        )
        hidexdecorations!(ax_p)
        hidespines!(ax_p)
        linkyaxes!(ax, ax_p)
        xlims!(ax_p, -1.15 * span, 1.15 * span)
        # The link takes the limits of the empty right-hand axis: set both again.
        ylims!(ax, y_lims...)
        ylims!(ax_p, y_lims...)
    end

    markers = _LegendElement[
        MarkerElement(;
            marker = _verdict_style(v).marker,
            color = _verdict_style(v).color,
            strokecolor = _band_edge(_verdict_style(v).color),
            strokewidth = _STYLE.marker_stroke,
            markersize = 1.3 * _STYLE.marker,
        ) for (v, _) in _VERDICTS
    ]
    band = _LegendElement[PolyElement(;
        color = (grey, 0.10),
        strokecolor = (grey, 0.55),
        strokewidth = _STYLE.guide,
        linestyle = :dash,
    ),]
    Legend(
        fig[0, :],
        [markers, band],
        [AbstractString[l for (_, l) in _VERDICTS], AbstractString["±ε"]],
        ["Verdict:", ""];   # ε is defined in the annotation; "" keeps the group gap
        orientation = :horizontal,
        nbanks = 1,
        framevisible = false,
        titleposition = :left,
        tellheight = true,
        padding = (0, 0, 0, 0),
    )
    return _save_fig(cfg, filename, fig)
end

"""
    calibration_figures(rows, series, cfg::VisualizationConfig; observables = <those of rows>,
                        grid = nothing, reference = "base", t_settle_myr = 0.0, α = 0.05,
                        holm = true) -> Vector{String}

The figures of the calibration report: per observable the time-resolved differences
([`plot_paired_differences`](@ref)) and the forest plot of the shifts
([`plot_paired_shifts`](@ref)) from one table `rows` and its `series`
([`paired_difference_series`](@ref)), so that every figure shares the arm colours.
Returns the paths in that order.
"""
function Nbody6Dynamics.calibration_figures(
    rows::AbstractVector,
    series::AbstractVector,
    cfg::VisualizationConfig;
    observables::AbstractVector{<:AbstractString} = unique(String(r.observable) for r in rows),
    grid::Union{Nothing,AbstractString} = nothing,
    reference::AbstractString = "base",
    t_settle_myr::Real = 0.0,
    α::Real = 0.05,
    holm::Bool = true,
)::Vector{String}
    paths = String[]
    for o in observables
        push!(
            paths,
            plot_paired_differences(
                series,
                rows,
                cfg;
                observable = o,
                grid,
                reference,
                t_settle_myr,
            ),
        )
        push!(paths, plot_paired_shifts(rows, cfg; observable = o, grid, reference, α, holm))
    end
    return paths
end
