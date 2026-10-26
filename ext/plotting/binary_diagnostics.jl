# =============================================================================
# Binary diagnostics: span timeline, population series, survival curves and
# the a–e / f–e planes of the compact pairs
# =============================================================================

# ---------------------------------------------------------------------------
# Binary timeline (span bars)
# ---------------------------------------------------------------------------

"""Bar colour of a span origin: `(origin, legend label, index into _OKABE_ITO)`; any other origin is drawn grey as `Other`."""
const _SPAN_ORIGIN_STYLE = (
    (:primordial, "Primordial", 1),             # blue
    (:primordial_reformed, "Primordial, re-formed", 5),   # sky blue
    (:exchange, "Exchange", 2),                 # orange
    (:dynamical, "Dynamical", 3),               # bluish green
)

"""
End marker of a span by its `end_class` (see [`binary_spans`](@ref)): a star
for `:gw_merger`; a diamond for the collision- and disruption-type ends
`:stellar_merger`, `:collision`, `:disruption`, `:sn_disruption`; a triangle
for `:exchange`; a square for `:escape_bound`; none for `:flyby` and for
`:censored_run_end`, whose span carries the arrow head of a right-censored
span instead; a circle for `:unresolved` and any class not listed.
"""
const _SPAN_END_MARKER = Dict{Symbol,Union{Nothing,Symbol}}(
    :gw_merger => :star5,
    :stellar_merger => :diamond,
    :collision => :diamond,
    :disruption => :diamond,
    :sn_disruption => :diamond,
    :exchange => :utriangle,
    :escape_bound => :rect,
    :flyby => nothing,
    :censored_run_end => nothing,
    :unresolved => :circle,
)

"""Legend labels of the end markers, in drawing order; `:rtriangle` is the arrow head of a right-censored span."""
const _SPAN_END_LEGEND = (
    (:star5, "GW merger"),
    (:diamond, "Collision/disruption"),
    (:utriangle, "Exchange"),
    (:rect, "Escape"),
    (:circle, "Unresolved"),
    (:rtriangle, "Open at end"),
)

"""Height of a span bar in units of the row spacing."""
const _SPAN_BAR_HEIGHT = 0.6

"""Layout units per bar of the timeline canvas."""
const _SPAN_ROW_UNITS = 12

"""Most bars whose pair names are written as y tick labels."""
const _SPAN_MAX_LABELLED = 30

"""Both component types of the opening record are white dwarfs, neutron stars or black holes (K* 10–14)."""
_timeline_compact(s::BinarySpan)::Bool = s.kstar1_form in 10:14 && s.kstar2_form in 10:14

"""
    _timeline_selection(spans, select) -> Vector{BinarySpan}

The spans of `spans` with finite formation and end times that `select` keeps:
`:compact` those of two compact components at formation
(`_timeline_compact`), `:primordial` those of primordial pairs
(`primordial`, origins `:primordial` and `:primordial_reformed`), `:all`
every span. Throws an `ArgumentError` for any other `select`.
"""
function _timeline_selection(spans::AbstractVector{BinarySpan}, select::Symbol)::Vector{BinarySpan}
    select in (:compact, :primordial, :all) ||
        throw(ArgumentError("select must be :compact, :primordial or :all; got :$select"))
    keep(s) =
        isfinite(s.t_form_myr) &&
        isfinite(s.t_end_myr) &&
        (select === :all || (select === :compact ? _timeline_compact(s) : s.primordial))
    return BinarySpan[s for s in spans if keep(s)]
end

"""Colour and legend label of a span origin."""
function _span_origin_style(origin::Symbol)
    for (key, label, idx) in _SPAN_ORIGIN_STYLE
        key === origin && return (_OKABE_ITO[idx], label)
    end
    return (Makie.to_color(:gray55), "Other")
end

"""
    plot_binary_timeline(spans::AbstractVector{BinarySpan}, cfg::VisualizationConfig;
                         filename = "binary_timeline", max_bars = 60,
                         select = :compact, t_end_myr = NaN) -> String

Swimmer plot of the engine-regularised binaries: one horizontal bar per
[`BinarySpan`](@ref) from its formation to its end, rows in order of
formation time, the earliest at the top, labelled by the pair names
`name1–name2` (the labels are dropped above $(_SPAN_MAX_LABELLED) bars, the ticks
kept). The bar colour encodes the origin (primordial blue, re-formed
primordial sky blue, exchange orange, dynamical bluish green), the marker at
the end of a bar the end class: a star for `:gw_merger`; a diamond for
`:stellar_merger`, `:collision`, `:disruption` and `:sn_disruption`; a
triangle for `:exchange`; a square for `:escape_bound`; none for `:flyby`; a
circle for `:unresolved` and any other class. A right-pointing
arrow head ends a right-censored span, an open circle starts a left-censored
one.

# Arguments
- `spans`: the spans of [`binary_spans`](@ref)
- `max_bars`: most bars drawn; beyond it the `max_bars` longest-lived spans are
  kept and their number stated in the panel
- `select`: `:compact` (both components white dwarfs, neutron stars or black
  holes at formation), `:primordial` (primordial pairs, re-formed spans
  included) or `:all`
- `t_end_myr`: right end of the time axis [Myr]; `NaN` takes the latest span end

The canvas grows by $(_SPAN_ROW_UNITS) layout units per bar beyond the
single-panel height. Throws an `ArgumentError` for an unknown `select` or
`max_bars < 1`, and an error when no span is selected.
"""
@publication function Nbody6Dynamics.plot_binary_timeline(
    spans::AbstractVector{BinarySpan},
    cfg::VisualizationConfig;
    filename::AbstractString = "binary_timeline",
    max_bars::Integer = 60,
    select::Symbol = :compact,
    t_end_myr::Real = NaN,
)::String
    max_bars ≥ 1 || throw(ArgumentError("max_bars must be ≥ 1; got $max_bars"))
    selected = _timeline_selection(spans, select)
    isempty(selected) && error("No spans to plot")
    n_selected = length(selected)
    lifetime = [s.t_end_myr - s.t_form_myr for s in selected]
    pool = n_selected > max_bars ? selected[sortperm(lifetime; rev = true)[1:max_bars]] : selected
    sort!(pool; by = s -> (s.t_form_myr, s.name1, s.name2))
    n = length(pool)

    ends = [s.t_end_myr for s in spans if isfinite(s.t_end_myr)]
    t_hi =
        isfinite(t_end_myr) && t_end_myr > 0 ? Float64(t_end_myr) :
        (isempty(ends) ? 0.0 : maximum(ends))
    t_hi > 0 || (t_hi = 1.0)

    width, height = _figsize_px(cfg)
    height = max(height, _SPAN_ROW_UNITS * n)
    # Row spacing in layout units, less the legend rows and the x decorations
    row_units = (height - 220) / n
    ms = clamp(1.1 * row_units, 6.0, _STYLE.marker)
    labelled = n ≤ _SPAN_MAX_LABELLED
    names = ["$(s.name1)–$(s.name2)" for s in pool]

    fig = Figure(; size = (width, height))
    ax = Axis(
        fig[1, 1];
        xlabel = L"t\;[\mathrm{Myr}]",
        ylabel = labelled ? "Pair" : "Binary, by formation time",
        xticks = _time_ticks(0, t_hi),
        yticks = (collect(1:n), names),
        yticklabelsvisible = labelled,
        yticklabelsize = clamp(0.8 * row_units, 10.0, Float64(_STYLE.tick)),
        yreversed = true,
        ygridvisible = false,
    )
    xlims!(ax, -0.02 * t_hi, 1.02 * t_hi)
    # Explicit limits set the axis direction: the earliest span at the top
    ylims!(ax, n + 0.7, 0.3)

    # Bars, one poly per origin
    h = _SPAN_BAR_HEIGHT
    origin_elements = _LegendElement[]
    origin_labels = AbstractString[]
    origins = unique([s.origin for s in pool])
    order = [key for (key, _, _) in _SPAN_ORIGIN_STYLE]
    sort!(origins; by = o -> something(findfirst(==(o), order), length(order) + 1))
    other_drawn = false
    for origin in origins
        c, label = _span_origin_style(origin)
        rects = [
            Rect2f(s.t_form_myr, i - h / 2, max(s.t_end_myr - s.t_form_myr, 0.0), h) for
            (i, s) in enumerate(pool) if s.origin === origin
        ]
        poly!(ax, rects; color = c, strokecolor = _band_edge(c), strokewidth = 1.0)
        label == "Other" && other_drawn && continue
        other_drawn |= label == "Other"
        push!(origin_elements, LineElement(; color = c, linewidth = 2 * _STYLE.data))
        push!(origin_labels, label)
    end

    # End markers, the arrow heads of right-censored spans, the open circles
    # at the start of left-censored ones
    points = Dict{Symbol,Vector{Point2f}}()
    starts = Point2f[]
    for (i, s) in enumerate(pool)
        marker = s.right_censored ? :rtriangle : get(_SPAN_END_MARKER, s.end_class, :circle)
        marker === nothing || push!(get!(Vector{Point2f}, points, marker), Point2f(s.t_end_myr, i))
        s.left_censored && push!(starts, Point2f(s.t_form_myr, i))
    end
    end_elements = _LegendElement[]
    end_labels = AbstractString[]
    for (marker, label) in _SPAN_END_LEGEND
        haskey(points, marker) || continue
        scatter!(
            ax,
            points[marker];
            marker = marker,
            color = :black,
            markersize = ms,
            strokewidth = 0,
        )
        push!(
            end_elements,
            MarkerElement(;
                marker = marker,
                color = :black,
                markersize = _STYLE.marker,
                strokewidth = 0,
            ),
        )
        push!(end_labels, label)
    end
    if !isempty(starts)
        scatter!(
            ax,
            starts;
            marker = :circle,
            color = :white,
            strokecolor = :black,
            strokewidth = _STYLE.marker_stroke,
            markersize = ms,
        )
        push!(
            end_elements,
            MarkerElement(;
                marker = :circle,
                color = :white,
                strokecolor = :black,
                strokewidth = _STYLE.marker_stroke,
                markersize = _STYLE.marker,
            ),
        )
        push!(end_labels, "Open at start")
    end

    if n_selected > n
        xs = Float64[]
        ys = Float64[]
        for (i, s) in enumerate(pool)
            append!(xs, (s.t_form_myr, (s.t_form_myr + s.t_end_myr) / 2, s.t_end_myr))
            append!(ys, (-i, -i, -i))
        end
        _annotate!(
            ax,
            "$(n) of $(n_selected) shown";
            corner = _emptiest_corner(
                xs,
                ys;
                xlims = (-0.02 * t_hi, 1.02 * t_hi),
                ylims = (-(n + 0.7), -0.3),
            ),
        )
    end

    # Two legend rows on top: the origins over the end markers
    if length(origin_elements) + length(end_elements) ≥ 2
        if isempty(end_elements)
            _top_legend!(
                fig,
                origin_elements,
                origin_labels;
                title = "Origin:",
                nbanks = cld(length(origin_elements), 3),
            )
        else
            _top_legend!(
                fig,
                end_elements,
                end_labels;
                title = "End:",
                nbanks = cld(length(end_elements), 2),
            )
            Legend(
                fig[-1, :],
                origin_elements,
                origin_labels,
                "Origin:";
                orientation = :horizontal,
                nbanks = cld(length(origin_elements), 3),
                framevisible = false,
                titleposition = :left,
                tellheight = true,
                padding = (0, 0, 0, 0),
            )
        end
    end

    return _save_fig(cfg, filename, fig)
end

# ---------------------------------------------------------------------------
# Binary population series
# ---------------------------------------------------------------------------

"""Census level letter of the annotation: E for `:ks`, P for `:bound`."""
_census_level_letter(level::Symbol)::String =
    level === :ks ? "E" : level === :bound ? "P" : String(level)

"""Legend label of the shell of a [`BinaryPopulationRow`](@ref), in the frame of the row."""
function _population_shell_label(r::BinaryPopulationRow)::AbstractString
    r.bin == 0 && return "All systems"
    lo, hi = r.edge_lo, r.edge_hi
    if r.radius_frame === :lagrangian
        return "$(_fmt_plain_sig(100 * lo, 3))–$(_fmt_plain_sig(100 * hi, 3)) % of mass"
    elseif r.radius_frame === :milone && r.bin == 1
        return L"R < R_\mathrm{c}"
    elseif r.radius_frame === :milone && r.bin == 2
        return L"R_\mathrm{c} \le R < R_\mathrm{h}"
    end
    x = r.radius_frame === :rh ? "r/r_\\mathrm{h}" : "R/R_\\mathrm{h}"
    lo == 0 && isinf(hi) && return "All radii"
    lo == 0 && return latexstring("$x < $(_fmt_latex_sig(hi, 3))")
    isinf(hi) && return latexstring("$x \\ge $(_fmt_latex_sig(lo, 3))")
    return latexstring("$(_fmt_latex_sig(lo, 3)) \\le $x < $(_fmt_latex_sig(hi, 3))")
end

"""
    plot_binary_population_series(rows::AbstractVector{BinaryPopulationRow},
                                  cfg::VisualizationConfig;
                                  filename = "binary_population_series",
                                  level = :bound, class = :all, bins = (0, 1, 2)) -> String

Two stacked panels sharing the time axis: the number fraction of binaries
``f_\\mathrm{b} = N_\\mathrm{b} / (N_\\mathrm{s} + N_\\mathrm{b})`` in per cent
with its binomial band ``f_\\mathrm{b} \\pm \\sigma_{f_\\mathrm{b}}`` above, the
mass fraction ``f_M = M_\\mathrm{b} / (M_\\mathrm{s} + M_\\mathrm{b})`` below.
One series per radial shell of `bins` (bin 0 holds every system, the others
are labelled by their edges in the frame of the rows), coloured in
Okabe–Ito order, with markers at the snapshots when there are fewer than 40.
The census level (E for `:ks`, P for `:bound`) and the stellar class are
stated in the upper panel.

# Arguments
- `rows`: the rows of [`binary_population_series`](@ref)
- `level`: census level, `:ks` or `:bound`
- `class`: `:all` or the `key` of a [`STELLAR_CLASSES`](@ref) entry
- `bins`: the shells drawn, in legend order

Throws an error when no row matches `level`, `class` and `bins`.
"""
@publication function Nbody6Dynamics.plot_binary_population_series(
    rows::AbstractVector{BinaryPopulationRow},
    cfg::VisualizationConfig;
    filename::AbstractString = "binary_population_series",
    level::Symbol = :bound,
    class::Symbol = :all,
    bins = (0, 1, 2),
)::String
    sel = [r for r in rows if r.level === level && r.class === class && r.bin in bins]
    isempty(sel) && error(
        "No binary population rows at level :$level, class :$class and bins $(Tuple(bins)) to plot",
    )
    shells = unique(Int[b for b in bins if any(r -> r.bin == b, sel)])
    times = sort!(unique([r.t_myr for r in sel]))
    single = length(times) == 1
    show_markers = length(times) < 40
    tmin, tmax = single ? (times[1] - 0.5, times[1] + 0.5) : (times[1], times[end])
    tticks = _time_ticks(tmin, tmax)

    fig = Figure(; size = _fig_multipanel(cfg, 2, 1))
    ax1 =
        Axis(fig[1, 1]; ylabel = L"f_\mathrm{b}\;[\%]", xticks = tticks, xticklabelsvisible = false)
    ax2 = Axis(fig[2, 1]; xlabel = L"t\;[\mathrm{Myr}]", ylabel = L"f_M\;[\%]", xticks = tticks)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    xlims!(ax2, tmin, tmax)

    elements = Vector{Vector{_LegendElement}}()
    labels = AbstractString[]
    upper_x = Float64[]
    upper_y = Float64[]
    lower_band = Float64[]
    lower_y = Float64[]
    for (k, b) in enumerate(shells)
        c = _OKABE_ITO[mod1(k, length(_OKABE_ITO))]
        rs = sort!([r for r in sel if r.bin == b]; by = r -> r.t_myr)
        t = [r.t_myr for r in rs]
        fb = [100 * r.f_b_number for r in rs]
        se = [100 * r.f_b_number_se for r in rs]
        fm = [100 * r.f_b_mass for r in rs]
        ok = isfinite.(fb) .& isfinite.(se)
        if count(ok) ≥ 2
            lo = clamp.(fb[ok] .- se[ok], 0.0, 100.0)
            hi = clamp.(fb[ok] .+ se[ok], 0.0, 100.0)
            band!(ax1, t[ok], lo, hi; color = (c, _STYLE.band_alpha))
            # Dense series (no markers) take thin band edges, or the edges hide the data lines
            edge_width = show_markers ? _STYLE.band_edge : 1.0
            lines!(ax1, t[ok], lo; color = _band_edge(c), linewidth = edge_width)
            lines!(ax1, t[ok], hi; color = _band_edge(c), linewidth = edge_width)
            append!(upper_x, t[ok])
            append!(upper_y, hi)
            append!(lower_band, lo)
        end
        for (ax, f) in ((ax1, fb), (ax2, fm))
            single || lines!(ax, t, f; color = c, linewidth = _STYLE.data)
            show_markers && scatter!(
                ax,
                t,
                f;
                color = c,
                markersize = _STYLE.marker,
                strokecolor = _band_edge(c),
                strokewidth = _STYLE.marker_stroke,
            )
        end
        fin = isfinite.(fb)
        append!(upper_x, t[fin])
        append!(upper_y, fb[fin])
        append!(lower_band, fb[fin])
        append!(lower_y, filter(isfinite, fm))

        entry = _LegendElement[]
        single || push!(entry, LineElement(; color = c, linewidth = _STYLE.data))
        show_markers && push!(
            entry,
            MarkerElement(;
                color = c,
                marker = :circle,
                markersize = _STYLE.marker,
                strokecolor = _band_edge(c),
                strokewidth = _STYLE.marker_stroke,
            ),
        )
        push!(elements, entry)
        push!(labels, _population_shell_label(rs[1]))
    end

    # Tight limits about the data (the lower bound at zero only when the data reach it)
    f_lo, f_hi = _pad_limits(
        isempty(lower_band) ? 0.0 : minimum(lower_band),
        isempty(upper_y) ? 100.0 : maximum(upper_y);
        floor = 0,
    )
    ylims!(ax1, f_lo, f_hi)
    ax1.yticks = _nice_ticks(f_lo, f_hi)
    m_lo, m_hi = _pad_limits(
        isempty(lower_y) ? 0.0 : minimum(lower_y),
        isempty(lower_y) ? 100.0 : maximum(lower_y);
        floor = 0,
    )
    ylims!(ax2, m_lo, m_hi)
    ax2.yticks = _nice_ticks(m_lo, m_hi)

    class_index = findfirst(c -> c.key === class, STELLAR_CLASSES)
    class_text =
        class === :all ? "all stars" :
        "class: " * (class_index === nothing ? String(class) : STELLAR_CLASSES[class_index].label)
    corner = _emptiest_corner(upper_x, upper_y; xlims = (tmin, tmax), ylims = (f_lo, f_hi))
    _annotate!(ax1, "Level $(_census_level_letter(level)), $(class_text)"; corner = corner)
    if length(shells) ≥ 2
        _top_legend!(fig, elements, labels; title = "Shell:", nbanks = length(shells) > 3 ? 2 : 1)
    else
        # A single shell needs no legend: it is named under the level
        _annotate!(ax1, labels[1]; corner = corner, dy = 0.08)
    end

    return _save_fig(cfg, filename, fig)
end

# ---------------------------------------------------------------------------
# Survival curves with a number-at-risk strip
# ---------------------------------------------------------------------------

"""Opacity of the confidence bands of the survival curves, kept low so that overlapping bands stay legible."""
const _SURVIVAL_BAND_ALPHA = 0.18

"""Fields every survival curve must carry."""
const _SURVIVAL_FIELDS = (:label, :time, :survival, :lower, :upper, :n_risk, :censor_times)
"""Largest number of censor ticks drawn per curve."""
const _SURVIVAL_MAX_TICKS = 60

"""Weight of the number-at-risk strip relative to the main panel."""
const _SURVIVAL_STRIP_WEIGHT = 0.4

"""Check the fields and lengths of survival curve `k`; throws `ArgumentError` or `DimensionMismatch`."""
function _check_survival_curve(c::NamedTuple, k::Int)
    for f in _SURVIVAL_FIELDS
        haskey(c, f) || throw(ArgumentError("survival curve $k has no field `$f`"))
    end
    n = length(c.time)
    all(length(getfield(c, f)) == n for f in (:survival, :lower, :upper, :n_risk)) || throw(
        DimensionMismatch(
            "survival curve $k: time, survival, lower, upper and n_risk must have equal lengths",
        ),
    )
    issorted(c.time) || throw(ArgumentError("survival curve $k: time must be sorted"))
    return nothing
end

"""Latest time of a survival curve: its last step or its last censoring, 0 when it has neither."""
_survival_extent(c::NamedTuple)::Float64 = max(
    isempty(c.time) ? 0.0 : Float64(maximum(c.time)),
    isempty(c.censor_times) ? 0.0 : Float64(maximum(c.censor_times)),
)

"""
    _survival_steps(c, x_end) -> (t, s, lower, upper)

Step nodes of a Kaplan–Meier curve for `step = :post`: the curve starts at
(0, 1) and is held at its last value up to `x_end`; non-finite band limits
collapse onto the curve.
"""
function _survival_steps(c::NamedTuple, x_end::Float64)
    t = Float64.(c.time)
    s = Float64.(c.survival)
    lo = [isfinite(l) ? Float64(l) : v for (l, v) in zip(c.lower, s)]
    hi = [isfinite(u) ? Float64(u) : v for (u, v) in zip(c.upper, s)]
    if isempty(t) || t[1] > 0
        pushfirst!(t, 0.0)
        pushfirst!(s, 1.0)
        pushfirst!(lo, 1.0)
        pushfirst!(hi, 1.0)
    end
    if x_end > t[end]
        push!(t, x_end)
        push!(s, s[end])
        push!(lo, lo[end])
        push!(hi, hi[end])
    end
    return t, s, lo, hi
end

"""Vertices of a post-step function through `(t, y)`, for a stepped `band!`."""
function _post_step_vertices(t::Vector{Float64}, y::Vector{Float64})
    xs = Float64[]
    ys = Float64[]
    for i in eachindex(t)
        push!(xs, t[i])
        push!(ys, y[i])
        if i < lastindex(t)
            push!(xs, t[i + 1])
            push!(ys, y[i])
        end
    end
    return xs, ys
end

"""
    _number_at_risk(c, x) -> Int

Subjects of survival curve `c` still at risk at time `x`: the number at risk
of the first step at or after `x`, plus the censorings in between (none after
the last step).
"""
function _number_at_risk(c::NamedTuple, x::Real)::Int
    i = searchsortedfirst(c.time, x)
    inside = i ≤ length(c.time)
    n_next = inside ? Int(c.n_risk[i]) : 0
    t_next = inside ? Float64(c.time[i]) : Inf
    return n_next + count(τ -> x ≤ τ < t_next, c.censor_times)
end

"""
    plot_survival_curves(curves::AbstractVector{<:NamedTuple}, cfg::VisualizationConfig;
                         filename = "binary_survival",
                         xlabel = L"\\tau\\;[\\mathrm{Myr}]") -> String

Kaplan–Meier survival curves of binary lifetimes. Main panel: one step curve
per element of `curves`, starting at (0, 1), with its 95 % confidence band as
a stepped low-opacity band and its censorings as short vertical ticks on the
curve. Strip below, sharing the time axis: the number at risk of every curve
at the tick positions, one row per curve in its colour. Colours follow the
Okabe–Ito order; above four curves the line style alternates between solid
and dashed.

# Arguments
- `curves`: NamedTuples `(label, time, survival, lower, upper, n_risk,
  censor_times)`; `time` sorted, with `survival`, the band limits `lower`,
  `upper` and the number at risk `n_risk` (just before each time) aligned to
  it, and `censor_times` the censored durations
- `xlabel`: label of the duration axis

Throws an error for an empty `curves`, an `ArgumentError` for a curve
without one of the fields or with unsorted times, and a `DimensionMismatch`
for misaligned vectors.
"""
@publication function Nbody6Dynamics.plot_survival_curves(
    curves::AbstractVector{<:NamedTuple},
    cfg::VisualizationConfig;
    filename::AbstractString = "binary_survival",
    xlabel = L"\tau\;[\mathrm{Myr}]",
)::String
    isempty(curves) && error("No survival curves to plot")
    for (k, c) in enumerate(curves)
        _check_survival_curve(c, k)
    end
    n_curves = length(curves)
    extent = maximum(_survival_extent, curves)
    x_max = extent > 0 ? extent : 1.0
    ticks = filter(x -> 0 ≤ x ≤ x_max, _nice_ticks(0, x_max))

    fig = Figure(; size = _fig_with_strip(cfg))
    ax1 = Axis(
        fig[1, 1];
        ylabel = L"Survival probability $S$",
        xticks = ticks,
        xticklabelsvisible = false,
    )
    ax2 = Axis(fig[2, 1]; xlabel = xlabel, xticks = ticks, yreversed = true, ygridvisible = false)
    hideydecorations!(ax2)
    linkxaxes!(ax1, ax2)
    rowsize!(fig.layout, 1, Auto(1.0))
    rowsize!(fig.layout, 2, Auto(_SURVIVAL_STRIP_WEIGHT))
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    xlims!(ax1, -0.04 * x_max, 1.04 * x_max)
    ylims!(ax2, 0.5, n_curves + 0.5)

    colours = [_OKABE_ITO[mod1(k, length(_OKABE_ITO))] for k in 1:n_curves]
    styles = [n_curves > 4 && iseven(k) ? :dash : :solid for k in 1:n_curves]
    steps = [_survival_steps(c, x_max) for c in curves]
    # Tight about the curves and their bands, the top always at unity
    s_min = minimum(
        min(minimum(st[2]), minimum(filter(isfinite, st[3]); init = 1.0)) for st in steps;
        init = 1.0,
    )
    y_floor = s_min > 0.3 ? max(0.0, floor(10 * (s_min - 0.04)) / 10) : 0.0
    ylims!(ax1, y_floor - 0.02 * (1 - y_floor), 1.0 + 0.02 * (1 - y_floor))
    y_step = 1 - y_floor ≤ 0.35 ? 0.05 : 1 - y_floor ≤ 0.6 ? 0.1 : 0.2
    ax1.yticks = y_floor:y_step:1.0

    # Bands beneath every curve
    for (k, (t, _, lo, hi)) in enumerate(steps)
        xs, lo_v = _post_step_vertices(t, lo)
        _, hi_v = _post_step_vertices(t, hi)
        band!(ax1, xs, lo_v, hi_v; color = (colours[k], _SURVIVAL_BAND_ALPHA))
    end
    strip_fontsize = min(_ANNOTATION_FONTSIZE, 0.7 * _CANVAS.strip / n_curves)
    for (k, c) in enumerate(curves)
        t, s, _, _ = steps[k]
        stairs!(
            ax1,
            t,
            s;
            step = :post,
            color = colours[k],
            linewidth = _STYLE.data,
            linestyle = styles[k],
        )
        if !isempty(c.censor_times)
            # At most 60 ticks per curve (evenly spaced in rank), or they hide the curve
            ct = sort!(Float64.(c.censor_times))
            if length(ct) > _SURVIVAL_MAX_TICKS
                ct = ct[round.(Int, range(1, length(ct); length = _SURVIVAL_MAX_TICKS))]
            end
            sc = [s[max(searchsortedlast(t, τ), 1)] for τ in ct]
            scatter!(ax1, ct, sc; marker = :vline, markersize = 12, color = colours[k])
        end
        for x in ticks
            text!(
                ax2,
                x,
                k;
                text = string(_number_at_risk(c, x)),
                align = (:center, :center),
                color = colours[k],
                fontsize = strip_fontsize,
            )
        end
    end

    if n_curves ≥ 2
        elements = _LegendElement[
            LineElement(; color = colours[k], linestyle = styles[k], linewidth = _STYLE.data)
            for k in 1:n_curves
        ]
        labels = AbstractString[String(c.label) for c in curves]
        _top_legend!(fig, elements, labels; nbanks = cld(n_curves, 3))
    end

    return _save_fig(cfg, filename, fig)
end

# ---------------------------------------------------------------------------
# a–e and f–e planes of the compact pairs
# ---------------------------------------------------------------------------

"""
Colour of a compact pair class: `(pair_class, legend label, index into
_OKABE_ITO)` — BH–BH vermillion, BH–NS reddish purple, NS–NS blue, BH–WD
bluish green, NS–WD sky blue, WD–WD yellow; any other class is grey.
"""
const _GW_CLASS_STYLE = (
    (:bhbh, "BH–BH", 6),
    (:bhns, "BH–NS", 4),
    (:nsns, "NS–NS", 1),
    (:bhwd, "BH–WD", 3),
    (:nswd, "NS–WD", 5),
    (:wdwd, "WD–WD", 7),
)

"""Marker of a formation channel: `(channel, legend label, marker)`; any other channel is a diamond."""
const _GW_CHANNEL_MARKER = (
    (:primordial, "Primordial", :circle),
    (:dynamical, "Dynamical", :utriangle),
    (:ejected, "Ejected", :rect),
)

"""Peters coalescence times of the guide lines of the a–e plane [Myr]."""
const _GW_CONTOUR_TIMES_MYR = (10.0, 100.0, 1000.0, 13_800.0)

"""Grid points per axis of the log₁₀ t_GW field behind the guide-line contours."""
const _GW_CONTOUR_GRID = 160

"""Margin around the data of the logarithmic axes [dex]."""
const _GW_MARGIN_DEX = 0.15

"""Lower edge of the ground-based band [Hz] and the width shown above it [dex]."""
const _GW_LVK_HZ = 10.0
const _GW_LVK_SPAN_DEX = 1.0

"""A record with a bound orbit (a > 0, 0 ≤ e < 1), positive masses and a positive f₂₂."""
_gw_plane_record(r::CompactBinaryRow)::Bool =
    isfinite(r.a_au) &&
    r.a_au > 0 &&
    isfinite(r.e) &&
    0 ≤ r.e < 1 &&
    isfinite(r.f_gw22_hz) &&
    r.f_gw22_hz > 0 &&
    r.m1_msun > 0 &&
    r.m2_msun > 0

function _gw_class_style(cls::Symbol)
    for (key, label, idx) in _GW_CLASS_STYLE
        key === cls && return (_OKABE_ITO[idx], label)
    end
    return (Makie.to_color(:gray55), "Other")
end

function _gw_channel_style(channel::Symbol)
    for (key, label, marker) in _GW_CHANNEL_MARKER
        key === channel && return (marker, label)
    end
    return (:diamond, "Other")
end

"""Logarithmic limits of `values` with a margin of `_GW_MARGIN_DEX`, opened to at least one decade."""
function _gw_log_limits(values::AbstractVector{<:Real})
    lo = minimum(values) / 10^_GW_MARGIN_DEX
    hi = maximum(values) * 10^_GW_MARGIN_DEX
    span = log10(hi / lo)
    if span < 1
        widen = 10^((1 - span) / 2)
        lo /= widen
        hi *= widen
    end
    return (lo, hi)
end

"""Median of a non-empty vector."""
function _gw_median(v::AbstractVector{<:Real})::Float64
    s = sort(Float64.(v))
    n = length(s)
    return isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2
end

"""Inline label of a Peters guide line: the time in Myr below 10⁴ Myr, in Gyr above."""
_gw_contour_label(t_myr::Float64) =
    t_myr < 1.0e4 ? latexstring("t_\\mathrm{GW} = $(_fmt_latex_sig(t_myr, 3))\\;\\mathrm{Myr}") :
    latexstring("t_\\mathrm{GW} = $(_fmt_latex_sig(t_myr / 1000, 3))\\;\\mathrm{Gyr}")

"""
    plot_gw_planes(rows::AbstractVector{CompactBinaryRow}, cfg::VisualizationConfig;
                   filename = "gw_planes", last_only = true,
                   band_hz = (1e-4, 1e-1)) -> String

The compact pairs in two planes side by side: semi-major axis against
``1 - e`` (both logarithmic) on the left, GW frequency ``f_{22}`` against
eccentricity on the right. One marker per pair at its last record (by
`time_nb`); with `last_only = false` every record is drawn faint and the last
one bold. Colour encodes the pair class (`pair_class`: BH–BH
vermillion, BH–NS reddish purple, NS–NS blue, BH–WD bluish green, NS–WD sky
blue, WD–WD yellow, others grey), marker shape the channel (primordial
circle, dynamical triangle, ejected square, other diamond).

The left panel carries the curves of constant Peters (1964,
DOI 10.1103/PhysRev.136.B1224) coalescence time
([`gw_inspiral_time`](@ref)) for the median component masses of the pairs,
at 10, 100 and 1000 Myr and 13.8 Gyr, labelled along the lines. The
right panel shades the space-borne band `band_hz` (LISA) and, when inside the
axis range, the ground-based band ``f \\ge 10`` Hz (LVK). The frequency axis
spans the data with a 0.15-dex margin and is extended to the band edges when
the data lie within two decades of them.

# Arguments
- `rows`: the rows of [`compact_binary_catalogue`](@ref); records without a
  bound orbit or a GW frequency are left out
- `last_only`: draw only the last record of every pair
- `band_hz`: (lo, hi) of the space-borne band [Hz]

Throws an error for an empty `rows` or when no record has a bound orbit, and
an `ArgumentError` for a band without 0 < lo < hi.
"""
@publication function Nbody6Dynamics.plot_gw_planes(
    rows::AbstractVector{CompactBinaryRow},
    cfg::VisualizationConfig;
    filename::AbstractString = "gw_planes",
    last_only::Bool = true,
    band_hz = (1e-4, 1e-1),
)::String
    isempty(rows) && error("No compact-binary records to plot")
    band_lo, band_hi = Float64(band_hz[1]), Float64(band_hz[2])
    (0 < band_lo < band_hi && isfinite(band_hi)) ||
        throw(ArgumentError("band_hz must satisfy 0 < lo < hi; got $(band_hz)"))
    recs = filter(_gw_plane_record, rows)
    isempty(recs) &&
        error("No compact-binary records with a bound orbit and a GW frequency to plot")

    # The last record of every pair
    last_of = Dict{Tuple{Int,Int},Int}()
    for (i, r) in enumerate(recs)
        key = minmax(r.name1, r.name2)
        j = get(last_of, key, 0)
        (j == 0 || r.time_nb ≥ recs[j].time_nb) && (last_of[key] = i)
    end
    is_last = falses(length(recs))
    foreach(j -> is_last[j] = true, values(last_of))
    bold = recs[is_last]
    faint = last_only ? CompactBinaryRow[] : recs[.!is_last]
    drawn = vcat(bold, faint)

    a = [r.a_au for r in drawn]
    y = [1 - r.e for r in drawn]
    f = [r.f_gw22_hz for r in drawn]
    a_lo, a_hi = _gw_log_limits(a)
    y_lo, y_hi = _gw_log_limits(y)
    # Room above 1 − e = 1 for the markers of the circular pairs and the guide-line labels
    maximum(y) ≥ 0.5 && (y_hi = max(y_hi, 10^0.3))
    # The data, extended to the band edges within two decades of them, and to
    # one decade of the ground-based band once the data reach 0.1 Hz
    f_span = [minimum(f), maximum(f)]
    maximum(f) ≥ band_lo / 100 && minimum(f) ≤ band_hi * 100 && append!(f_span, (band_lo, band_hi))
    maximum(f) ≥ _GW_LVK_HZ / 100 &&
        push!(f_span, _GW_LVK_HZ * 10^(_GW_LVK_SPAN_DEX - _GW_MARGIN_DEX))
    f_lo, f_hi = _gw_log_limits(f_span)

    # Legend content first: its rows are added to the canvas height
    classes = unique([r.pair_class for r in drawn])
    class_order = [key for (key, _, _) in _GW_CLASS_STYLE]
    sort!(classes; by = c -> something(findfirst(==(c), class_order), length(class_order) + 1))
    channels = unique([r.channel for r in drawn])
    channel_order = [key for (key, _, _) in _GW_CHANNEL_MARKER]
    sort!(channels; by = c -> something(findfirst(==(c), channel_order), length(channel_order) + 1))
    class_elements = _LegendElement[]
    class_labels = AbstractString[]
    for cls in classes
        c, label = _gw_class_style(cls)
        label in class_labels && continue
        push!(
            class_elements,
            MarkerElement(;
                color = c,
                marker = :circle,
                markersize = _STYLE.marker,
                strokecolor = _band_edge(c),
                strokewidth = _STYLE.marker_stroke,
            ),
        )
        push!(class_labels, label)
    end
    channel_elements = _LegendElement[]
    channel_labels = AbstractString[]
    for ch in channels
        marker, label = _gw_channel_style(ch)
        label in channel_labels && continue
        push!(
            channel_elements,
            MarkerElement(;
                color = :gray60,
                marker = marker,
                markersize = _STYLE.marker,
                strokecolor = _band_edge(:gray60),
                strokewidth = _STYLE.marker_stroke,
            ),
        )
        push!(channel_labels, label)
    end
    n_entries = length(class_elements) + length(channel_elements)
    nbanks = n_entries > 6 ? 2 : 1
    legend_height = n_entries ≥ 2 ? 34 * nbanks + 10 : 0

    fig = Figure(; size = _fig_multipanel(cfg, 1, 2; extra_height = legend_height))
    ax1 = Axis(
        fig[1, 1];
        xlabel = L"a\;[\mathrm{au}]",
        ylabel = L"1 - e",
        xscale = log10,
        yscale = log10,
        xticks = _log_ticks(a_lo, a_hi),
        yticks = _log_ticks(y_lo, y_hi),
    )
    ax2 = Axis(
        fig[1, 2];
        xlabel = L"f_{22}\;[\mathrm{Hz}]",
        ylabel = L"e",
        xscale = log10,
        xticks = _log_ticks(f_lo, f_hi),
        yticks = 0.0:0.2:1.0,
    )
    colgap!(fig.layout, _multipanel_gap(2))
    xlims!(ax1, a_lo, a_hi)
    ylims!(ax1, y_lo, y_hi)
    xlims!(ax2, f_lo, f_hi)
    # e < 1 by construction: the band above unity holds the band labels.
    ylims!(ax2, -0.03, _E_AXIS_TOP)

    # Peters guide lines beneath the data, for the median masses: contours of
    # log₁₀ t_GW over the plane, labelled inline along the lines
    m1 = _gw_median([r.m1_msun for r in bold])
    m2 = _gw_median([r.m2_msun for r in bold])
    a_grid = 10 .^ range(log10(a_lo), log10(a_hi); length = _GW_CONTOUR_GRID)
    y_grid = 10 .^ range(log10(y_lo), log10(min(y_hi, 1.0)); length = _GW_CONTOUR_GRID)
    log_t = [log10(gw_inspiral_time(m1, m2, a, 1 - yy) / 1.0e6) for a in a_grid, yy in y_grid]
    contour!(
        ax1,
        a_grid,
        y_grid,
        log_t;
        levels = log10.(collect(_GW_CONTOUR_TIMES_MYR)),
        color = :gray45,
        linestyle = :dash,
        linewidth = _STYLE.guide,
        labels = true,
        labelsize = _ANNOTATION_FONTSIZE,
        labelcolor = :gray35,
        labelformatter = level -> _gw_contour_label(10.0^level),
    )

    # Detector bands beneath the data
    y_label = (1 + _E_AXIS_TOP) / 2
    for (lo, hi, name) in ((band_lo, band_hi, "LISA"), (_GW_LVK_HZ, Inf, "LVK"))
        v_lo, v_hi = max(lo, f_lo), min(hi, f_hi)
        v_lo < v_hi || continue
        vspan!(ax2, v_lo, v_hi; color = (:gray50, 0.12))
        edges = filter(x -> f_lo < x < f_hi, [lo, hi])
        isempty(edges) ||
            vlines!(ax2, edges; color = :gray40, linestyle = :dash, linewidth = _STYLE.guide)
        text!(
            ax2,
            sqrt(v_lo * v_hi),
            y_label;
            text = name,
            align = (:center, :center),
            color = :gray30,
            fontsize = _ANNOTATION_FONTSIZE,
        )
    end

    # Records: the faint history first, the last records bold on top
    scale = max(_multipanel_scale(cfg, 2), 0.8)
    for (set, alpha, stroke) in ((faint, 0.25, 0.0), (bold, 1.0, _STYLE.marker_stroke))
        isempty(set) && continue
        for cls in unique([r.pair_class for r in set]), ch in unique([r.channel for r in set])
            sub = [r for r in set if r.pair_class === cls && r.channel === ch]
            isempty(sub) && continue
            c, _ = _gw_class_style(cls)
            marker, _ = _gw_channel_style(ch)
            kw = (
                color = (c, alpha),
                marker = marker,
                markersize = _STYLE.marker * scale,
                strokecolor = _band_edge(c),
                strokewidth = stroke,
            )
            scatter!(ax1, [r.a_au for r in sub], [1 - r.e for r in sub]; kw...)
            scatter!(ax2, [r.f_gw22_hz for r in sub], [r.e for r in sub]; kw...)
        end
    end

    if n_entries ≥ 2
        groups = Vector{Vector{_LegendElement}}()
        group_labels = Vector{Vector{AbstractString}}()
        titles = String[]
        if !isempty(class_elements)
            push!(groups, class_elements)
            push!(group_labels, class_labels)
            push!(titles, "Class:")
        end
        if !isempty(channel_elements)
            push!(groups, channel_elements)
            push!(group_labels, channel_labels)
            push!(titles, "Channel:")
        end
        Legend(
            fig[0, :],
            groups,
            group_labels,
            titles;
            orientation = :horizontal,
            nbanks = nbanks,
            framevisible = false,
            titleposition = :left,
            tellheight = true,
            padding = (0, 0, 0, 0),
        )
    end

    return _save_fig(cfg, filename, fig)
end
