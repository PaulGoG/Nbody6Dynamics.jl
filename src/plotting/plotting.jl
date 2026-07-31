# =============================================================================
# Plotting submodule — includes all visualisation routines
# =============================================================================

# ---------------------------------------------------------------------------
# Publication-quality theme with Computer Modern (LaTeX) fonts
# ---------------------------------------------------------------------------

# Makie's MathTeXEngine renders L"..." strings in Computer Modern automatically.
# For regular (non-LaTeX) text, we use the "NewComputerModern" font family which
# ships with MathTeXEngine (a CairoMakie dependency), so no system font install
# is needed.

const _CM_FONT = let
    # MathTeXEngine (a Makie dependency) bundles NewComputerModern OTF fonts.
    # Search the Julia package depot for them.
    _found = ""
    for depot in Base.DEPOT_PATH
        candidate = joinpath(depot, "packages", "MathTeXEngine")
        isdir(candidate) || continue
        for entry in readdir(candidate; join = true)
            font_dir = joinpath(entry, "assets", "fonts", "NewComputerModern")
            if isdir(font_dir) && isfile(joinpath(font_dir, "NewCM10-Regular.otf"))
                _found = font_dir
                break
            end
        end
        isempty(_found) || break
    end

    if !isempty(_found)
        (
            regular     = joinpath(_found, "NewCM10-Regular.otf"),
            bold        = joinpath(_found, "NewCM10-Bold.otf"),
            italic      = joinpath(_found, "NewCM10-Italic.otf"),
            bold_italic = joinpath(_found, "NewCM10-BoldItalic.otf"),
        )
    else
        # Fallback: let Makie use its defaults if fonts are not found
        (regular = "serif",)
    end
end

const PUBLICATION_THEME = Theme(
    fontsize = 22,
    fonts    = _CM_FONT,
    figure_padding = 16,
    Axis = (
        xlabelsize         = 20,
        ylabelsize         = 20,
        titlesize          = 20,
        xticklabelsize     = 16,
        yticklabelsize     = 16,
        xlabelpadding      = 10.0,
        ylabelpadding      = 10.0,
        spinewidth         = 1.5,
        xtickwidth         = 1.2,
        ytickwidth         = 1.2,
        xminortickwidth    = 0.8,
        yminortickwidth    = 0.8,
        xtickalign         = 1.0,     # ticks face inward
        ytickalign         = 1.0,
        xticksize          = 8,
        yticksize          = 8,
        # No minor ticks; grey dashed major grid at very low opacity.
        # Dense scatter plots (cluster projections, HR) disable the grid
        # locally via x/ygridvisible = false.
        xminorticksvisible = false,
        yminorticksvisible = false,
        xgridvisible       = true,
        ygridvisible       = true,
        xgridstyle         = :dash,
        ygridstyle         = :dash,
        xgridcolor         = (:grey, 0.12),
        ygridcolor         = (:grey, 0.12),
        topspinevisible    = true,
        rightspinevisible  = true,
    ),
    Legend = (
        framevisible = true,
        framewidth   = 1.0,
        labelsize    = 15,
        patchsize    = (25, 14),
        padding      = (8, 8, 6, 6),
        rowgap       = 4,
    ),
    Lines = (
        linewidth = 2.2,
    ),
    Colorbar = (
        labelsize     = 18,
        ticklabelsize = 14,
        tickalign     = 1.0,
        width         = 14,
    ),
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

"""Convert (width_in, height_in) → Makie screen units (1 unit = 1/72 inch)."""
_figsize_px(cfg::VisualizationConfig) = (cfg.figsize[1] * 72, cfg.figsize[2] * 72)

# ---------------------------------------------------------------------------
# Standardised figure sizes for publication-consistent box dimensions
# ---------------------------------------------------------------------------
# All single-panel plots share the same axis-box area.  Plots with extra
# elements (colorbar, second panel) get a wider or taller canvas so the
# primary axis box remains the same physical size.

"""Single-panel with a right-side colorbar — extra width keeps axis box size."""
_fig_with_colorbar(cfg::VisualizationConfig) =
    (_figsize_px(cfg)[1] + 110, _figsize_px(cfg)[2])

"""Two vertically stacked panels (e.g. energy + virial)."""
_fig_two_panel(cfg::VisualizationConfig) =
    (_figsize_px(cfg)[1], round(Int, _figsize_px(cfg)[2] * 1.45))

"""
Multi-panel grid — each panel matches the single-panel size.
Returns `(total_width, total_height)` in Makie screen units.
Use the same `hgap`/`vgap` values with `colgap!`/`rowgap!` on the layout.
"""
const _MULTIPANEL_HGAP = 70
const _MULTIPANEL_VGAP = 70

function _fig_multipanel(cfg::VisualizationConfig, nrows::Int, ncols::Int)
    pw, ph = _figsize_px(cfg)
    return (ncols * pw + (ncols - 1) * _MULTIPANEL_HGAP,
            nrows * ph + (nrows - 1) * _MULTIPANEL_VGAP)
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
    return collect(t0:step:t1)
end

"""
Compute square (equal-span) axis limits centered on the origin (0,0),
snapped outward to a nice round value with a visual buffer.
Returns `(lo, hi, lo, hi)` — symmetric on both axes so ticks land on
clean multiples and the outermost tick never overlaps the spine.
"""
function _square_limits(xs, ys; pad_frac::Float64 = 0.08)
    # Maximum absolute extent across both coordinates
    r = max(maximum(abs, xs), maximum(abs, ys))
    # Snap outward to a nice round value: 10, 20, 25, 50, 100, 200, ...
    pow = 10.0^floor(log10(max(r, 1.0)))
    candidates = [1.0, 2.0, 2.5, 5.0, 10.0]
    # Find the smallest nice value ≥ r
    hs = pow * candidates[findfirst(c -> c * pow >= r, candidates)]
    # Add visual buffer so the outermost tick sits inside the axis
    hs += hs * pad_frac
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
    return collect(t0:step:t1)
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
save `fig` at the configured DPI, and log the location. Returns the path.
"""
function _save_fig(cfg::VisualizationConfig, basename::AbstractString, fig)::String
    path = _output_path(cfg, basename)
    _backup_existing(path)
    save(path, fig; px_per_unit = cfg.dpi / 72)
    @info "Saved: $path"
    return path
end

"""
    _marker_size(cfg, n) -> Float64

Scatter marker size for `n` particles: `marker_budget / n` clamped to
`[marker_min, marker_max]` (all from `cfg.style`).
"""
_marker_size(cfg::VisualizationConfig, n::Integer) =
    clamp(cfg.style.marker_budget / max(n, 1), cfg.style.marker_min, cfg.style.marker_max)


# ---------------------------------------------------------------------------
# Include plot source files
# ---------------------------------------------------------------------------

include("snapshots.jl")
include("lagrangian.jl")
include("energy.jl")
include("hr.jl")
include("animation.jl")
include("merger.jl")
