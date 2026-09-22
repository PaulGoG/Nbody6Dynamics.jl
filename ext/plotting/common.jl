# =============================================================================
# Shared theme, layout, export and helper routines of the figure layer
# =============================================================================
#
# Every figure is composed in layout units on a canvas fixed by its type, under
# one theme, and exported so that the canvas width becomes the printed width:
# text, line weights and markers keep their proportions whatever the target.
# A journal column is an override of the printed width, never of the layout.

# ---------------------------------------------------------------------------
# Theme
# ---------------------------------------------------------------------------

# Computer Modern via MathTeXEngine's texfont API (direct dependency). The
# faces are resolved when the theme is built, never in a constant: a
# FreeType face created while the package precompiles is serialised with a
# null pointer, and Makie then falls back to its default sans font for every
# plain-text label without a word (LaTeX strings, which MathTeXEngine
# renders with its live faces, were the only text in Computer Modern).

"""
    _cm_fonts() -> NamedTuple

MathTeXEngine's Computer Modern faces (`texfont(:text)`, `:bold`,
`:italic`) for the `fonts` attribute of the theme, taken from the live
registry at call time.
"""
_cm_fonts() = (regular = texfont(:text), bold = texfont(:bold), italic = texfont(:italic))

"""
    _STYLE

Sizes and weights by role, in layout units. Figure routines name the role a
mark plays instead of a number, so one series family looks the same in every
figure: `data` for the series plotted, `fit` for a model or secondary curve
beside it, `envelope` for the extremes about a mean, `ghost` for greyed
background curves, `guide` for reference lines (dashed, labelled), `band_edge`
for the outline of an area fill of opacity `band_alpha`; `marker` with a
`marker_stroke` edge of the same, darker hue; `label`, `tick` and `annotation`
font sizes.
"""
const _STYLE = (
    data = 3.0,
    fit = 2.0,
    envelope = 1.5,
    ghost = 1.5,
    guide = 1.5,
    band_edge = 2.5,
    band_alpha = 0.35,
    marker = 14.0,
    marker_stroke = 1.5,
    spine = 1.5,
    label = 26,
    tick = 22,
    annotation = 21,
)

"""
    publication_theme() -> Theme

The publication theme: Computer Modern fonts from MathTeXEngine, 26-unit
labels and legends over 22-unit tick labels, boxed axes with inward ticks and
no minor ticks, a faint dashed grid, frameless horizontal legends, 3-unit
data lines and 14-unit markers with a 1.5-unit edge. Built on every call so
that the font faces are live ones, not faces captured at precompile time.
Every figure routine of the package draws inside
`with_theme(publication_theme())`, so the session's own theme is neither
needed nor changed; [`set_publication_theme!`](@ref) activates it globally
for figures composed by hand.
"""
function Nbody6Dynamics.publication_theme()
    return Theme(
        fonts = _cm_fonts(),
        fontsize = _STYLE.label,
        figure_padding = (10, 30, 10, 10),   # right: the overhang of the last x tick label
        linewidth = _STYLE.data,
        markersize = _STYLE.marker,
        Axis = (
            spinewidth = _STYLE.spine,
            xticklabelsize = _STYLE.tick,
            yticklabelsize = _STYLE.tick,
            xlabelpadding = 8.0,
            ylabelpadding = 8.0,
            xgridstyle = :dash,
            ygridstyle = :dash,
            xgridcolor = (:grey, 0.12),
            ygridcolor = (:grey, 0.12),
            xminorticksvisible = false,
            yminorticksvisible = false,
            xtickalign = 1,       # ticks face inward
            ytickalign = 1,
            xtickwidth = _STYLE.spine,
            ytickwidth = _STYLE.spine,
            xticksize = 10,
            yticksize = 10,
            topspinevisible = true,
            rightspinevisible = true,
        ),
        Lines = (linewidth = _STYLE.data,),
        Scatter = (strokewidth = _STYLE.marker_stroke,),
        Legend = (
            framevisible = false,
            orientation = :horizontal,
            titlefont = :bold,
            patchsize = (34, 20),
            rowgap = 4,
            padding = (0, 0, 0, 0),
        ),
        Colorbar = (
            labelsize = _STYLE.label,
            ticklabelsize = _STYLE.tick,
            tickalign = 1,
            width = 18,
        ),
    )
end

"""
    set_publication_theme!()

Activate [`publication_theme`](@ref) globally, for figures composed by hand
from Makie calls. The figure routines of the package do not depend on it.
"""
function Nbody6Dynamics.set_publication_theme!()
    set_theme!(publication_theme())
end

"""
    @publication function Nbody6Dynamics.plot_something(args...; kwargs...) … end

Define a figure routine whose body runs inside
`with_theme(publication_theme())`. The theme then holds however the routine
is reached — from the dispatcher, a sweep driver or a user's session — and
the session's global theme is left alone.
"""
macro publication(definition)
    Meta.isexpr(definition, :function) && length(definition.args) == 2 ||
        throw(ArgumentError("@publication expects a long-form function definition"))
    body = definition.args[2]
    definition.args[2] = quote
        with_theme(publication_theme()) do
            $body
        end
    end
    return esc(definition)
end

# ---------------------------------------------------------------------------
# Canvases and export
# ---------------------------------------------------------------------------

"""
Canvas sizes in layout units: a single panel, the height each further stacked
main panel and each auxiliary strip (ratio, residual, census) adds, and the
least width of a grid of panels together with the width it takes per column.
"""
const _CANVAS = (
    width = 900,
    height = 600,
    stacked_panel = 350,
    strip = 180,
    grid_width = 1200,
    grid_column = 500,
)

"""Printed widths of the journal column presets [in]."""
const _COLUMN_WIDTH_IN = Dict("single" => 3.4, "double" => 7.05)

"""Canvas of a single-panel figure in layout units. The configuration sets the printed width of the export, not the canvas."""
_figsize_px(::VisualizationConfig) = (_CANVAS.width, _CANVAS.height)

"""Printed width of an export [in]: the journal column when `cfg.column` names one, `cfg.export_width` otherwise."""
_export_width_in(cfg::VisualizationConfig)::Float64 =
    get(_COLUMN_WIDTH_IN, cfg.column, cfg.export_width)

# ---------------------------------------------------------------------------
# Standardised figure sizes for publication-consistent box dimensions
# ---------------------------------------------------------------------------
# All single-panel plots share the same axis-box area.  Plots with extra
# elements (colorbar, second panel) get a wider or taller canvas so the
# primary axis box remains the same physical size.

"""Extra canvas width reserved for a right-side colorbar column."""
const _COLORBAR_WIDTH = 140

"""Single-panel with a right-side colorbar — extra width keeps axis box size."""
_fig_with_colorbar(cfg::VisualizationConfig) =
    (_figsize_px(cfg)[1] + _COLORBAR_WIDTH, _figsize_px(cfg)[2])

"""Two vertically stacked main panels (e.g. energy + virial)."""
_fig_two_panel(cfg::VisualizationConfig) =
    (_figsize_px(cfg)[1], _figsize_px(cfg)[2] + _CANVAS.stacked_panel)

"""
Grids of panels get a wider canvas than a single panel —
`max(_CANVAS.grid_width, ncols × _CANVAS.grid_column)` layout units — shared by
the `ncols` axis boxes after one axis-decoration strip (`_AXIS_PROTRUSION`,
the outer y label and tick labels), the reserved `extra_width` (a shared
colorbar column) and the gaps, at `panel_aspect` (height over width; that of
the single panel by default); the rows stack and one decoration strip is
added below for the x label and tick labels, plus `extra_height` for whatever
sits above or below the grid (a legend row, an auxiliary strip). The panel
gap is [`_multipanel_gap`](@ref): compact when the inner tick labels are
hidden, the stack value otherwise; use it with `colgap!`/`rowgap!`.
Single-column stacks (`ncols = 1`) are a full single panel plus
`_CANVAS.stacked_panel` for every further row. Returns
`(total_width, total_height)` in layout units.
"""
const _MULTIPANEL_HGAP = 84
const _MULTIPANEL_VGAP = 84
const _MULTIPANEL_GAP_COMPACT = 28

"""Layout units one axis label and its tick labels take along one side (label, tick labels, pads)."""
const _AXIS_PROTRUSION = 84

"""Canvas width of a figure with `ncols` columns of panels."""
_grid_canvas_width(ncols::Int)::Int =
    ncols == 1 ? _CANVAS.width : max(_CANVAS.grid_width, ncols * _CANVAS.grid_column)

"""
    _fit_canvas_to_boxes!(fig, nrows, ncols, box_width, aspect)

Give the `nrows × ncols` panel cells of `fig` a fixed size — `box_width` wide,
`aspect` times as high — and resize the canvas to the layout. For panels whose
shape is prescribed (equal-aspect projections) the canvas must follow the
boxes: sized the other way round, the boxes shrink inside their cells and
leave blank gutters, and a colourbar spanning the rows overshoots them.
"""
function _fit_canvas_to_boxes!(fig::Figure, nrows::Int, ncols::Int, box_width::Real, aspect::Real)
    foreach(c -> colsize!(fig.layout, c, Fixed(box_width)), 1:ncols)
    foreach(r -> rowsize!(fig.layout, r, Fixed(box_width * aspect)), 1:nrows)
    resize_to_layout!(fig)
    return nothing
end

"""Gap between the panels of a grid: compact when `ncols > 1` and the inner tick labels are hidden, `_MULTIPANEL_HGAP` otherwise."""
_multipanel_gap(ncols::Int; inner_ticks::Bool = true) =
    (ncols > 1 && !inner_ticks) ? _MULTIPANEL_GAP_COMPACT : _MULTIPANEL_HGAP

"""Fraction of the data range kept free above the data in montage panels, so the in-axis time annotation never meets a marker."""
const _MONTAGE_BAND_FRAC = 0.18

"""Width of one axis box of a grid (`_fig_multipanel` conventions) in canvas units."""
function _multipanel_box_width(
    cfg::VisualizationConfig,
    ncols::Int;
    inner_ticks::Bool = true,
    extra_width::Real = 0,
)
    gap = _multipanel_gap(ncols; inner_ticks)
    return (_grid_canvas_width(ncols) - extra_width - _AXIS_PROTRUSION - (ncols - 1) * gap) / ncols
end

"""Width of one axis box of a grid relative to the single-panel width (1 for stacks); scales markers with the panel."""
function _multipanel_scale(
    cfg::VisualizationConfig,
    ncols::Int;
    inner_ticks::Bool = true,
    extra_width::Real = 0,
)
    ncols == 1 && return 1.0
    return _multipanel_box_width(cfg, ncols; inner_ticks, extra_width) / _figsize_px(cfg)[1]
end

function _fig_multipanel(
    cfg::VisualizationConfig,
    nrows::Int,
    ncols::Int;
    inner_ticks::Bool = true,
    panel_aspect::Union{Nothing,Real} = nothing,
    extra_width::Real = 0,
    extra_height::Real = 0,
)
    pw, ph = _figsize_px(cfg)
    ncols == 1 && return (pw, ph + (nrows - 1) * _CANVAS.stacked_panel + round(Int, extra_height))
    gap = _multipanel_gap(ncols; inner_ticks)
    box_w = _multipanel_box_width(cfg, ncols; inner_ticks, extra_width)
    aspect = panel_aspect === nothing ? ph / pw : Float64(panel_aspect)
    height = nrows * box_w * aspect + (nrows - 1) * gap + _AXIS_PROTRUSION + extra_height
    return (_grid_canvas_width(ncols), round(Int, height))
end

# ---------------------------------------------------------------------------
# Tick and limit helpers for consistent axis formatting
# ---------------------------------------------------------------------------

"""
Compute nice round tick values spanning `[lo, hi]`.
Steps are chosen from multiples of 1, 2, 2.5, 5, 10 (×10^n), targeting
approximately `target_n` ticks.  Returns a `Vector{Float64}`.
"""
function _nice_ticks(lo::Real, hi::Real; target_n::Int = 8)
    span = Float64(hi - lo)
    span ≤ 0 && return [Float64(lo)]
    raw_step = span / target_n
    # Round to a nice step: 1, 2, 2.5, 5, 10, 20, 25, 50, ...
    pow = 10.0^floor(log10(raw_step))
    candidates = [1.0, 2.0, 2.5, 5.0, 10.0]
    step = pow * candidates[argmin(abs.(candidates .* pow .- raw_step))]
    t0 = ceil(lo / step) * step
    t1 = floor(hi / step) * step
    ticks = collect(t0:step:t1)
    # A range narrower than one step holds no rounded tick: use its ends.
    return isempty(ticks) ? [Float64(lo), Float64(hi)] : ticks
end

"""
Compute square (equal-span) axis limits centered on the origin (0,0):
data extent plus a small margin (tight limits per the economy-of-space
standard). Returns `(lo, hi, lo, hi)` — symmetric on both axes.
"""
function _square_limits(xs, ys; pad_frac::Float64 = 0.03)
    r = max(maximum(abs, xs), maximum(abs, ys))
    hs = r * (1 + pad_frac)
    return (-hs, hs, -hs, hs)
end

"""
Compute nice colorbar tick values for a `[cmin, cmax]` range.
"""
function _nice_colorbar_ticks(cmin::Real, cmax::Real; target_n::Int = 6)
    _nice_ticks(cmin, cmax; target_n = target_n)
end

"""
Compute nice time-axis ticks with higher density than the default.
"""
function _time_ticks(tmin::Real, tmax::Real)
    _nice_ticks(tmin, tmax; target_n = 10)
end

"""
Compute ticks for axes whose values are already logarithmic (e.g. log Teff,
log L).  Steps are anchored to multiples of 0.5 or 1.0 for narrow ranges,
or 2/5 for wide ranges.
"""
function _logval_ticks(lo::Real, hi::Real; target_n::Int = 8)
    span = Float64(hi - lo)
    span ≤ 0 && return [Float64(lo)]
    raw_step = span / target_n
    # Choose from 0.5, 1, 2, 5 (suited for log-value axes)
    candidates = [0.5, 1.0, 2.0, 5.0]
    step = candidates[argmin(abs.(candidates .- raw_step))]
    t0 = ceil(lo / step) * step
    t1 = floor(hi / step) * step
    ticks = collect(t0:step:t1)
    # Narrow log-value ranges (identical stars) fall between multiples of the
    # coarsest step: fall back to a finer step, then to the range ends.
    if isempty(ticks)
        ticks = collect((ceil(lo / 0.1) * 0.1):0.1:(floor(hi / 0.1) * 0.1))
    end
    return isempty(ticks) ? [Float64(lo), Float64(hi)] : ticks
end

"""
Compute integer-only ticks for a `[lo, hi]` range.  Useful for count data
(N, N_pairs) where fractional ticks are physically meaningless.
"""
function _integer_ticks(lo::Real, hi::Real; target_n::Int = 6)
    ilo = ceil(Int, lo)
    ihi = floor(Int, hi)
    ihi < ilo && return [Float64(round(Int, (lo + hi) / 2))]
    span = ihi - ilo
    span == 0 && return [Float64(ilo)]
    # Choose a nice integer step
    raw_step = max(1, span ÷ target_n)
    candidates = [1, 2, 5, 10, 20, 25, 50, 100, 200, 500, 1000]
    step = candidates[argmin(abs.(candidates .- raw_step))]
    t0 = ceil(Int, ilo / step) * step
    t1 = floor(Int, ihi / step) * step
    return Float64.(collect(t0:step:t1))
end

"""Build output file path with the configured format extension."""
function _output_path(cfg::VisualizationConfig, basename::AbstractString)::String
    mkpath(cfg.output_dir)
    return joinpath(cfg.output_dir, basename * "." * cfg.format)
end

"""Build output path for animations (always .gif)."""
function _anim_output_path(cfg::VisualizationConfig, basename::AbstractString)::String
    mkpath(cfg.output_dir)
    return joinpath(cfg.output_dir, basename * ".gif")
end

"""Least raster density of an export, in pixels per layout unit."""
const _MIN_PX_PER_UNIT = 4.0

"""
    _export_scale(cfg, canvas_width) -> (pt_per_unit, px_per_unit)

Scales that map a canvas `canvas_width` layout units wide onto the printed
width ([`_export_width_in`](@ref)): points per unit for vector output, and
pixels per unit for raster output at `cfg.dpi` over that width, never
coarser than `_MIN_PX_PER_UNIT`.
"""
function _export_scale(cfg::VisualizationConfig, canvas_width::Real)
    width_in = _export_width_in(cfg)
    return (
        pt_per_unit = 72 * width_in / canvas_width,
        px_per_unit = max(_MIN_PX_PER_UNIT, cfg.dpi * width_in / canvas_width),
    )
end

"""
    _save_fig(cfg, basename, fig) -> String

Resolve the output path, back up any existing file (never-overwrite policy),
save `fig` and log the location. The canvas width of `fig` becomes the
printed width ([`_export_scale`](@ref)), so the export enters a document at
native size. Returns the path.
"""
function _save_fig(cfg::VisualizationConfig, basename::AbstractString, fig)::String
    path = _output_path(cfg, basename)
    _backup_existing(path)
    scale = _export_scale(cfg, size(fig.scene)[1])
    if cfg.format in ("pdf", "svg")
        save(path, fig; pt_per_unit = scale.pt_per_unit)
    else
        save(path, fig; px_per_unit = scale.px_per_unit)
    end
    @info "Saved: $path"
    return path
end

"""Legend entry types accepted by the grouped legends; Makie needs the entry vectors concretely typed."""
const _LegendElement = Union{LineElement,MarkerElement,PolyElement}

"""
    _top_legend!(fig, ax; title = nothing, nbanks = 1, kwargs...)
    _top_legend!(fig, elements, labels; title = nothing, nbanks = 1, kwargs...)

Standard legend placement: horizontal, above the axes, outside the plot
area (`fig[0, :]`), optional bold family header via `title` shown to the
left. Callers must still guard against single-entry legends.
"""
function _top_legend!(fig::Figure, ax::Axis; title = nothing, nbanks::Int = 1, kwargs...)
    # The family title must be passed positionally — Legend's convenience
    # constructors take (layout, ax, title); a `title` kwarg is ignored.
    args = title === nothing ? (ax,) : (ax, title)
    Legend(
        fig[0, :],
        args...;
        orientation = :horizontal,
        nbanks = nbanks,
        framevisible = false,
        titleposition = :left,
        tellheight = true,
        padding = (0, 0, 0, 0),
        kwargs...,
    )
    return nothing
end

function _top_legend!(
    fig::Figure,
    elements::AbstractVector,
    labels::AbstractVector;
    title = nothing,
    nbanks::Int = 1,
    kwargs...,
)
    args = title === nothing ? (elements, labels) : (elements, labels, title)
    Legend(
        fig[0, :],
        args...;
        orientation = :horizontal,
        nbanks = nbanks,
        framevisible = false,
        titleposition = :left,
        tellheight = true,
        padding = (0, 0, 0, 0),
        kwargs...,
    )
    return nothing
end

"""
    _log_ticks(lo, hi) -> (values, labels)

Decade-anchored ticks for `log10`-scaled axes: `10^n` at every decade in
range, with 2× and 5× intermediates when the range spans ≤ 2 decades.
Labels follow the axis-typography standard: plain decimals throughout
(`0.01, 0.1, 1, 10, 100`) when every tick lies within 10⁻³–10⁴ and the
ticks span at most four decades; otherwise the exponent form, in which
`10^0`, `10^1` and the 2×/5× multiples of `10^{-1}`–`10^{1}` still
collapse to `1`, `10`, `0.2`, `5`, `20`.
"""
function _log_ticks(lo::Real, hi::Real)
    lo, hi = min(lo, hi), max(lo, hi)
    lo > 0 || (lo = hi / 1e3)          # guard: log axes need positive range
    e_lo = floor(Int, log10(lo) + 1e-12)
    e_hi = ceil(Int, log10(hi) - 1e-12)
    # Sparse-decade test counts decades actually inside [lo, hi] — the
    # exponent-bin span overcounts when the range endpoints sit mid-decade.
    n_dec = count(e -> lo * (1 - 1e-9) ≤ 10.0^e ≤ hi * (1 + 1e-9), e_lo:e_hi)
    mults = n_dec ≤ 2 ? (1.0, 2.0, 5.0) : (1.0,)
    vals = Float64[]
    for e in e_lo:e_hi, m in mults
        v = m * 10.0^e
        lo * (1 - 1e-9) ≤ v ≤ hi * (1 + 1e-9) && push!(vals, v)
    end
    length(vals) < 2 && (vals = [10.0^e_lo, 10.0^e_hi])
    exponents = [floor(Int, log10(v) + 1e-9) for v in vals]
    plain = all(e -> -3 ≤ e ≤ 4, exponents) && maximum(exponents) - minimum(exponents) ≤ 4
    labels = [latexstring(_log_tick_label(v, plain)) for v in vals]
    return (vals, labels)
end

"""
    _log_tick_label(v, plain) -> String

LaTeX body of a log-axis tick label for `v = m × 10^e` with `m` ∈ {1, 2, 5}:
the plain decimal when `plain` is set or `-1 ≤ e ≤ 1` (the mandatory
collapses `10^0 → 1`, `10^1 → 10`, `2 × 10^{-1} → 0.2`), otherwise
`10^{e}` or `m \\times 10^{e}`.
"""
function _log_tick_label(v::Real, plain::Bool)::String
    e = floor(Int, log10(v) + 1e-9)
    m = round(Int, v / 10.0^e)
    if plain || -1 ≤ e ≤ 1
        return e < 0 ? string(round(v; digits = -e)) : string(round(Int, v))
    end
    return m == 1 ? "10^{$(e)}" : "$(m)\\times 10^{$(e)}"
end

"""
    _fmt_latex_sig(x, n = 3) -> String

Format a non-negative quantity to `n` significant digits for labels and
annotations, following the power-of-ten typography standard: plain decimal
within 10⁻²–10⁴, mantissa `\\times 10^{e}` outside (`10^{e}` alone when the
mantissa is 1) — never computer notation, which `@sprintf("%.3g", 1200)`
and `@sprintf("%.2g", 100)` both produce.
"""
function _fmt_latex_sig(x::Real, n::Int = 3)::String
    x == 0 && return "0"
    v = round(Float64(x); sigdigits = n)
    e = floor(Int, log10(abs(v)))
    if -2 ≤ e ≤ 3
        return isinteger(v) ? string(Int(v)) : string(v)
    end
    m = round(v / 10.0^e; sigdigits = n)
    m == 1 && return "10^{$(e)}"
    m_str = isinteger(m) ? string(Int(m)) : string(m)
    return "$(m_str) \\times 10^{$(e)}"
end

"""`_fmt_latex_sig(x, 3)`."""
_fmt_latex_sig3(x::Real)::String = _fmt_latex_sig(x, 3)

"""
    _superscript(e) -> String

The integer `e` in Unicode superscript digits, with U+207B for a negative
sign (`-5` → `⁻⁵`).
"""
function _superscript(e::Integer)::String
    digits_sup = ('⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹')
    s = string(abs(e))
    body = join(digits_sup[Int(c - '0') + 1] for c in s)
    return e < 0 ? "⁻" * body : body
end

"""
    _fmt_plain_sig(x, n = 3) -> String

`_fmt_latex_sig` for plain-text contexts (legend entries, multi-line
annotations): the same rules, with the power of ten written in Unicode
superscripts (`2.5×10⁻⁵`, `10⁴`) instead of TeX.
"""
function _fmt_plain_sig(x::Real, n::Int = 3)::String
    x == 0 && return "0"
    sign = x < 0 ? "-" : ""
    v = round(abs(Float64(x)); sigdigits = n)
    e = floor(Int, log10(v))
    if -2 ≤ e ≤ 3
        return sign * (isinteger(v) ? string(Int(v)) : string(v))
    end
    m = round(v / 10.0^e; sigdigits = n)
    m == 1 && return sign * "10" * _superscript(e)
    m_str = isinteger(m) ? string(Int(m)) : string(m)
    return sign * m_str * "×10" * _superscript(e)
end

"""Legend or annotation text of a grid-axis value: booleans and integers verbatim, reals to three significant digits (`_fmt_plain_sig`), strings as they are."""
_legend_value(v::Bool) = string(v)
_legend_value(v::Integer) = string(v)
_legend_value(v::Real) = _fmt_plain_sig(v, 3)
_legend_value(v) = String(string(v))

"""
    _marker_size(cfg, n) -> Float64

Scatter marker size for `n` particles: `marker_budget / n` clamped to
`[marker_min, marker_max]` (all from `cfg.style`).
"""
_marker_size(cfg::VisualizationConfig, n::Integer) =
    clamp(cfg.style.marker_budget / max(n, 1), cfg.style.marker_min, cfg.style.marker_max)

# ---------------------------------------------------------------------------
# Shared annotation, layout, and data-preparation helpers
# ---------------------------------------------------------------------------

"""Standard fontsize for in-axis annotations (quantitative takeaways)."""
const _ANNOTATION_FONTSIZE = _STYLE.annotation

"""Row gap between the stacked panels of two-panel figures."""
const _TWO_PANEL_ROWGAP = 12

"""Column gap between an axis and its colorbar."""
const _COLORBAR_COLGAP = 10

"""
    _annotate!(ax, text; corner = :tl, color = :black,
               fontsize = _ANNOTATION_FONTSIZE, dy = 0.0)

In-axis annotation at a standard corner in relative coordinates, coloured to
the relevant series per the annotation standard. `corner` is one of `:tl`,
`:tr`, `:br`; `dy` shifts the anchor downward for stacked annotations.
`text` may be an `Observable` (animations).
"""
function _annotate!(
    ax,
    text;
    corner::Symbol = :tl,
    color = :black,
    fontsize::Real = _ANNOTATION_FONTSIZE,
    dy::Real = 0.0,
)
    x, y, align = if corner === :tl
        (0.04, 0.96 - dy, (:left, :top))
    elseif corner === :tr
        (0.96, 0.96 - dy, (:right, :top))
    elseif corner === :br
        (0.96, 0.04 + dy, (:right, :bottom))
    else
        throw(ArgumentError("corner must be :tl, :tr, or :br; got $corner"))
    end
    text!(
        ax,
        x,
        y;
        text = text,
        space = :relative,
        align = align,
        fontsize = fontsize,
        color = color,
    )
    return nothing
end

"""
    _emptiest_corner(xs, ys; corners = (:tl, :tr, :br), width = 0.35, height = 0.2)

The corner of the data's bounding box (fractions `width` × `height` of the
ranges) holding the fewest points among `corners`, for placing an
annotation clear of series whose course is not known in advance. Ties
resolve in the order given.
"""
function _emptiest_corner(
    xs::AbstractVector{<:Real},
    ys::AbstractVector{<:Real};
    corners = (:tl, :tr, :br),
    width::Real = 0.35,
    height::Real = 0.2,
    avoid_x::Real = NaN,
)
    isempty(xs) && return first(corners)
    x0, x1 = extrema(xs)
    y0, y1 = extrema(ys)
    dx = max(x1 - x0, eps(Float64))
    dy = max(y1 - y0, eps(Float64))
    u_avoid = isfinite(avoid_x) ? (avoid_x - x0) / dx : NaN
    counts = map(corners) do c
        right = c in (:tr, :br)
        top = c in (:tl, :tr)
        # A vertical event marker spans the whole panel height, so it rules out
        # both corners of its own side however empty the data leave them; one
        # data point at the marker's position cannot express that.
        if isfinite(u_avoid) && (right ? u_avoid ≥ 1 - width : u_avoid ≤ width)
            return typemax(Int)
        end
        count(zip(xs, ys)) do (x, y)
            u = (x - x0) / dx
            v = (y - y0) / dy
            (right ? u ≥ 1 - width : u ≤ width) && (top ? v ≥ 1 - height : v ≤ height)
        end
    end
    return corners[argmin(counts)]
end

"""Centred grey note for panels with no plottable data."""
function _no_data_note!(ax, text)
    text!(
        ax,
        0.5,
        0.55;
        text = text,
        space = :relative,
        align = (:center, :center),
        color = :gray30,
        fontsize = _STYLE.annotation,
    )
    return nothing
end

"""
    _log_color_range(values) -> (log_vals, cmin, cmax)

`log10` colour scale for mass colouring: floors values at `1e-30` and widens
a degenerate (single-value) range by ±0.5 so the colormap stays defined.
"""
function _log_color_range(values)
    log_vals = log10.(max.(values, 1e-30))
    cmin, cmax = extrema(log_vals)
    if cmin ≈ cmax
        cmin -= 0.5
        cmax += 0.5
    end
    return log_vals, cmin, cmax
end

"""
    _envelope_stats(values::AbstractMatrix) -> (lo, hi, mean)

Per-column (per-timestep) minimum / maximum / mean over the non-NaN entries
of a series × times matrix; all-NaN columns stay NaN in every output.
"""
function _envelope_stats(values::AbstractMatrix{<:Real})
    n_series, n_t = size(values)
    lo = fill(NaN, n_t)
    hi = fill(NaN, n_t)
    mean_vals = fill(NaN, n_t)
    for k in 1:n_t
        vals = Float64[]
        for i in 1:n_series
            x = values[i, k]
            isnan(x) || push!(vals, x)
        end
        isempty(vals) && continue
        lo[k] = minimum(vals)
        hi[k] = maximum(vals)
        mean_vals[k] = sum(vals) / length(vals)
    end
    return lo, hi, mean_vals
end

# ---------------------------------------------------------------------------
# Semantic colour table
# ---------------------------------------------------------------------------
# NOTE: must precede the includes below — hr.jl builds its stellar-type
# colour map from _OKABE_ITO at include time.

# Okabe–Ito colourblind-safe palette: Makie's Wong colours (7 entries) plus
# black, completing the 8-colour Okabe–Ito set.  Order:
# 1 blue, 2 orange, 3 bluish green, 4 reddish purple, 5 sky blue,
# 6 vermillion, 7 yellow, 8 black.
const _OKABE_ITO = vcat(Makie.wong_colors(), Makie.RGBAf(0, 0, 0, 1))

"""One consistent colour per physical quantity across every figure of the
project (series family encoded by colour, role by line style)."""
const _SEMANTIC_COLORS = Dict{Symbol,Makie.RGBAf}(
    :energy_error => _OKABE_ITO[1],  # blue          |ΔE/E|
    :virial => _OKABE_ITO[6],  # vermillion    Q = T/|W|
    :n_particles => _OKABE_ITO[3],  # bluish green  N (counts)
    :n_pairs => _OKABE_ITO[2],  # orange        N_pairs
    :separation => _OKABE_ITO[5],  # sky blue      pairwise separations
    :escapers => _OKABE_ITO[4],  # reddish purple escaper counts & cumulative mass
    :binary_hard => _OKABE_ITO[6],  # vermillion    hard pairs (E_b > ⟨m⟩σ²)
    :binary_soft => _OKABE_ITO[5],  # sky blue      soft pairs
    :rotation => _OKABE_ITO[1],  # blue          λ_R, v_rot/σ
    :spin => _OKABE_ITO[4],  # reddish purple Peebles spin λ_P
    :core_radius => _OKABE_ITO[6],  # vermillion    r_c
    :half_mass_radius => _OKABE_ITO[3],  # bluish green  r_h
    :segregation => _OKABE_ITO[2],  # orange        Λ_MSR, r_h ratio
)

"""
Darkened same-hue edge colour, at full opacity, for area fills and marker
edges. Accepts any colour Makie does (`:black`, `(colour, alpha)`, `RGBAf`).
"""
function _band_edge(colour)::Makie.RGBAf
    c = Makie.to_color(colour)
    return Makie.RGBAf(0.7 * c.r, 0.7 * c.g, 0.7 * c.b, 1.0)
end
