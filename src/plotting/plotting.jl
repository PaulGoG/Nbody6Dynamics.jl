# =============================================================================
# Plotting submodule — includes all visualisation routines
# =============================================================================

# ---------------------------------------------------------------------------
# Publication-quality theme with Computer Modern (LaTeX) fonts
# ---------------------------------------------------------------------------

# Makie's MathTeXEngine renders L"..." strings in Computer Modern automatically.

# Computer Modern via MathTeXEngine's texfont API (direct dependency) —
# no depot scanning, no silent fallback: if the fonts are missing this
# fails loudly at load time rather than degrading to serif.
const _CM_FONT = (regular = texfont(:text), bold = texfont(:bold), italic = texfont(:italic))

const PUBLICATION_THEME = Theme(
    fontsize = 22,
    fonts = _CM_FONT,
    figure_padding = 16,
    Axis = (
        xlabelsize = 20,
        ylabelsize = 20,
        xticklabelsize = 16,
        yticklabelsize = 16,
        xlabelpadding = 10.0,
        ylabelpadding = 10.0,
        spinewidth = 1.5,
        xtickwidth = 1.2,
        ytickwidth = 1.2,
        xtickalign = 1.0,     # ticks face inward
        ytickalign = 1.0,
        xticksize = 8,
        yticksize = 8,
        # No minor ticks; grey dashed major grid at very low opacity.
        # Dense scatter plots (cluster projections, HR) disable the grid
        # locally via x/ygridvisible = false.
        xminorticksvisible = false,
        yminorticksvisible = false,
        xgridvisible = true,
        ygridvisible = true,
        xgridstyle = :dash,
        ygridstyle = :dash,
        xgridcolor = (:grey, 0.12),
        ygridcolor = (:grey, 0.12),
        topspinevisible = true,
        rightspinevisible = true,
    ),
    Legend = (
        framevisible = true,
        framewidth = 1.0,
        labelsize = 15,
        patchsize = (25, 14),
        padding = (8, 8, 6, 6),
        rowgap = 4,
    ),
    Lines = (linewidth = 2.2,),
    Colorbar = (labelsize = 18, ticklabelsize = 14, tickalign = 1.0, width = 14),
)

"""
    set_publication_theme!()

Activate the publication-quality Makie theme with Computer Modern fonts globally.
"""
function set_publication_theme!()
    set_theme!(PUBLICATION_THEME)
end

# ---------------------------------------------------------------------------
# Common helpers
# ---------------------------------------------------------------------------

# Journal column-width presets (final printed size, inches) and the internal
# render scale: figures are designed at final×_PRINT_SCALE canvas units so
# the 22 pt theme text lands at 22/_PRINT_SCALE ≈ 8.8 pt and 2.2-unit lines
# at ≈ 0.9 pt when the export is reduced to the true column width.
const _PRINT_SCALE = 2.5
const _COLUMN_PRESETS = Dict(
    "single" => (3.4, 2.6),    # ≈ 86–90 mm single column
    "double" => (7.05, 4.35),  # ≈ 178–183 mm double column
)

"""Render scale of the canvas relative to the final printed size (1.0 for
free-form `figsize` canvases, `_PRINT_SCALE` for column presets)."""
_render_scale(cfg::VisualizationConfig) = haskey(_COLUMN_PRESETS, cfg.column) ? _PRINT_SCALE : 1.0

"""Canvas size in Makie units (1 unit = 1 pt): column preset × render scale,
or the free-form `figsize` inches when `column` is empty/unknown."""
function _figsize_px(cfg::VisualizationConfig)
    if haskey(_COLUMN_PRESETS, cfg.column)
        w, h = _COLUMN_PRESETS[cfg.column]
        return (w * _PRINT_SCALE * 72, h * _PRINT_SCALE * 72)
    end
    return (cfg.figsize[1] * 72, cfg.figsize[2] * 72)
end

# ---------------------------------------------------------------------------
# Standardised figure sizes for publication-consistent box dimensions
# ---------------------------------------------------------------------------
# All single-panel plots share the same axis-box area.  Plots with extra
# elements (colorbar, second panel) get a wider or taller canvas so the
# primary axis box remains the same physical size.

"""Extra canvas width reserved for a right-side colorbar column."""
const _COLORBAR_WIDTH = 110

"""Single-panel with a right-side colorbar — extra width keeps axis box size."""
_fig_with_colorbar(cfg::VisualizationConfig) =
    (_figsize_px(cfg)[1] + _COLORBAR_WIDTH, _figsize_px(cfg)[2])

"""Two vertically stacked panels (e.g. energy + virial)."""
_fig_two_panel(cfg::VisualizationConfig) =
    (_figsize_px(cfg)[1], round(Int, _figsize_px(cfg)[2] * 1.45))

"""
Multi-panel grids stay one column wide: the canvas keeps the preset width,
the `ncols` panels share it minus the gaps at `panel_aspect` (height over
width; the preset's ratio by default), and the rows stack. A montage
therefore enters a document at native size like every other figure; at the
`single` preset a three-column grid has 25 mm panels, so montage callers
use three ticks per axis. `extra_width` is reserved on the right (a shared
colorbar) and taken from the panel area. The panel gap is
[`_multipanel_gap`](@ref): compact when the inner tick labels are hidden,
the stack value otherwise; use it with `colgap!`/`rowgap!`. Single-column
stacks (`ncols = 1`) keep panels of the full single-panel size. Returns
`(total_width, total_height)` in Makie screen units.
"""
const _MULTIPANEL_HGAP = 70
const _MULTIPANEL_VGAP = 70
const _MULTIPANEL_GAP_COMPACT = 24

"""Gap between the panels of a grid: compact when `ncols > 1` and the inner tick labels are hidden, `_MULTIPANEL_HGAP` otherwise."""
_multipanel_gap(ncols::Int; inner_ticks::Bool = true) =
    (ncols > 1 && !inner_ticks) ? _MULTIPANEL_GAP_COMPACT : _MULTIPANEL_HGAP

"""Fraction of the data range kept free above the data in montage panels, so the in-axis time annotation never meets a marker."""
const _MONTAGE_BAND_FRAC = 0.18

"""Width of one panel of a grid relative to the single-panel width (1 for stacks); scales markers with the panel."""
function _multipanel_scale(
    cfg::VisualizationConfig,
    ncols::Int;
    inner_ticks::Bool = true,
    extra_width::Real = 0,
)
    ncols == 1 && return 1.0
    pw = _figsize_px(cfg)[1]
    gap = _multipanel_gap(ncols; inner_ticks)
    return (pw - extra_width - (ncols - 1) * gap) / ncols / pw
end

function _fig_multipanel(
    cfg::VisualizationConfig,
    nrows::Int,
    ncols::Int;
    inner_ticks::Bool = true,
    panel_aspect::Union{Nothing,Real} = nothing,
    extra_width::Real = 0,
)
    pw, ph = _figsize_px(cfg)
    ncols == 1 && return (pw, nrows * ph + (nrows - 1) * _MULTIPANEL_VGAP)
    gap = _multipanel_gap(ncols; inner_ticks)
    panel_w = (pw - extra_width - (ncols - 1) * gap) / ncols
    aspect = panel_aspect === nothing ? ph / pw : Float64(panel_aspect)
    return (pw, round(Int, nrows * panel_w * aspect + (nrows - 1) * gap))
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

"""
    _save_fig(cfg, basename, fig) -> String

Resolve the output path, back up any existing file (never-overwrite policy),
save `fig` sized for print, and log the location. Vector formats (pdf/svg)
are scaled so the document width equals the true column width; raster
output gets `cfg.dpi` at that final size. Returns the path.
"""
function _save_fig(cfg::VisualizationConfig, basename::AbstractString, fig)::String
    path = _output_path(cfg, basename)
    _backup_existing(path)
    s = _render_scale(cfg)
    if cfg.format in ("pdf", "svg")
        save(path, fig; pt_per_unit = 1 / s)
    else
        save(path, fig; px_per_unit = cfg.dpi / 72 / s)
    end
    @info "Saved: $path"
    return path
end

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
`10^0` renders as `1`, per the axis-typography standard.
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
    labels = map(vals) do v
        e = floor(Int, log10(v) + 1e-9)
        m = round(Int, v / 10.0^e)
        if m == 1
            e == 0 ? L"1" : latexstring("10^{$(e)}")
        else
            e == 0 ? latexstring("$(m)") : latexstring("$(m)\\times 10^{$(e)}")
        end
    end
    return (vals, labels)
end

"""
    _fmt_latex_sig3(x) -> String

Format a non-negative quantity to 3 significant digits for LaTeX
annotations, following the power-of-ten typography standard: plain
decimal within 10⁻²–10⁴, mantissa `\\times 10^{e}` outside — never
computer notation.
"""
function _fmt_latex_sig3(x::Real)::String
    x == 0 && return "0"
    m_str, e_str = split(@sprintf("%.2e", x), 'e')
    e = parse(Int, e_str)
    if -2 ≤ e ≤ 3
        v = round(x; sigdigits = 3)
        return isinteger(v) ? string(Int(v)) : string(v)
    end
    return "$(m_str) \\times 10^{$(e)}"
end

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

"""Standard fontsize for in-axis annotations (quantitative takeaways, §10)."""
const _ANNOTATION_FONTSIZE = 16

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
)
    isempty(xs) && return first(corners)
    x0, x1 = extrema(xs)
    y0, y1 = extrema(ys)
    dx = max(x1 - x0, eps(Float64))
    dy = max(y1 - y0, eps(Float64))
    counts = map(corners) do c
        count(zip(xs, ys)) do (x, y)
            u = (x - x0) / dx
            v = (y - y0) / dy
            right = c in (:tr, :br)
            top = c in (:tl, :tr)
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
        fontsize = 18,
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

"""Darkened same-hue edge colour for `band!` fills (edge at full opacity)."""
_band_edge(c::Makie.RGBAf) = Makie.RGBAf(0.7 * c.r, 0.7 * c.g, 0.7 * c.b, 1.0)

# ---------------------------------------------------------------------------
# Include plot source files
# ---------------------------------------------------------------------------

include("snapshots.jl")
include("lagrangian.jl")
include("energy.jl")
include("hr.jl")
include("escapers.jl")
include("sse.jl")
include("animation.jl")
include("merger.jl")
include("binaries.jl")
include("sweep.jl")
include("ensemble.jl")
include("remnant.jl")
include("control.jl")
