# =============================================================================
# Snapshot visualisation — 2D projections and density maps
# =============================================================================

"""
    plot_snapshot(snap::Snapshot, cfg::VisualizationConfig;
                  filename::AbstractString = "snapshot",
                  projections::Vector{Symbol} = [:xy, :xz])

Generate publication-quality plots of particle positions in the requested
projections (:xy, :xz, :yz): a scatter coloured by mass, or the projected
surface mass density when `cfg.style.snapshot_render` selects it for this
particle count (see [`PlotStyle`](@ref)).  Positions are shown in pc, the
density in M☉ pc⁻² and the time annotation in Myr when
`cfg.units == "physical"` (header AS scaling); N-body units otherwise.
"""
@publication function Nbody6Dynamics.plot_snapshot(
    snap::Snapshot,
    cfg::VisualizationConfig;
    filename::AbstractString = "snapshot",
    projections::Vector{Symbol} = [:xy, :xz, :yz],
)
    n = nparticles(snap)
    density = _density_render(cfg, n)
    physical = _physical_units(cfg, [snap], density)
    unit_str = physical ? "pc" : "NB"
    r_scale = physical ? rbar(snap.header) : 1.0
    t_val = physical ? time_myr(snap.header) : time_nb(snap.header)

    ms = _marker_size(cfg, n)

    # Colour by mass (log scale for visual contrast)
    log_m, cmin, cmax = _log_color_range(Float64.(snap.mass))

    for proj in projections
        fig = Figure(; size = _fig_with_colorbar(cfg))

        ix, iy, xsym, ysym = _proj_indices(proj)
        px = snap.pos[ix, :] .* r_scale
        py = snap.pos[iy, :] .* r_scale

        # Square limits and nice ticks for consistent box sizes
        m = density ? _body_masses(snap, physical) : Float64[]
        xlo, xhi, ylo, yhi =
            density ? _centred_limits(_mass_extent(px, py, m, cfg.style.density_mass_frac)) :
            _square_limits(px, py)
        xtk = _nice_ticks(xlo, xhi; target_n = 5)
        ytk = _nice_ticks(ylo, yhi; target_n = 5)

        ax = Axis(
            fig[1, 1];
            xlabel = _coord_label(xsym, unit_str),
            ylabel = _coord_label(ysym, unit_str),
            aspect = DataAspect(),
            limits = (xlo, xhi, ylo, yhi),
            xticks = xtk,
            yticks = ytk,
            xgridvisible = false,
            ygridvisible = false,
            _density_axis_attributes(density)...,
        )

        if density
            Σ = _surface_density(px, py, m, xhi, cfg.style.density_bins)
            hi = _log_peak(Σ)
            lo = hi - cfg.style.density_decades
            _density_heatmap!(ax, Σ, xhi, lo, hi)
            _annotate!(ax, _time_annotation(t_val, physical); _DENSITY_ANNOTATION...)
            _density_colorbar!(fig[1, 2], lo, hi, physical)
        else
            _annotate!(ax, _time_annotation(t_val, physical))
            sc = scatter!(
                ax,
                px,
                py;
                color = log_m,
                colormap = :viridis,
                colorrange = (cmin, cmax),
                markersize = ms,
                strokewidth = 0,
            )
            Colorbar(
                fig[1, 2],
                sc;
                label = L"\log_{10}(m \, / \, M_\mathrm{tot})",
                ticks = _nice_colorbar_ticks(cmin, cmax),
            )
        end

        colgap!(fig.layout, _COLORBAR_COLGAP)
        # A square axis in the 3:2 canvas would leave blank margins either
        # side: fix the box and let the canvas follow it.
        _fit_canvas_to_boxes!(fig, 1, 1, _figsize_px(cfg)[2] - _AXIS_PROTRUSION, 1.0)

        _save_fig(cfg, "$(filename)_$(proj)", fig)
    end

    return nothing
end

"""
    plot_snapshot_evolution(snaps::Vector{Snapshot}, cfg::VisualizationConfig;
                           filename::AbstractString = "snapshot_evolution",
                           projections::Vector{Symbol} = [:xy, :xz, :yz],
                           max_panels::Int = 6)

Multi-panel figure showing cluster evolution across snapshots, as mass-coloured
scatter panels or as projected surface mass density on one colour scale
(`cfg.style.snapshot_render`, decided once for the figure by its largest
snapshot).  Positions are shown in pc and time annotations in Myr when
`cfg.units == "physical"`.

When the spatial extent varies dramatically across panels (e.g. merger runs),
each panel is zoomed on its own data and carries its own tick values.  Under
global limits (the extent is comparable throughout) only the border panels
are labelled.  The extent of a density panel is that of the mass fraction
`cfg.style.density_mass_frac`, so escapers do not set the scale.
"""
@publication function Nbody6Dynamics.plot_snapshot_evolution(
    snaps::Vector{Snapshot},
    cfg::VisualizationConfig;
    filename::AbstractString = "snapshot_evolution",
    projections::Vector{Symbol} = [:xy, :xz, :yz],
    max_panels::Int = 6,
)
    ns = length(snaps)
    ns == 0 && return nothing

    # Select evenly-spaced snapshots
    indices = ns <= max_panels ? (1:ns) : round.(Int, range(1, ns; length = max_panels))

    density = _density_render(cfg, maximum(nparticles(snaps[i]) for i in indices))
    physical = _physical_units(cfg, snaps[indices], density)
    unit_str = physical ? "pc" : "NB"
    r_scales = [physical ? rbar(s.header) : 1.0 for s in snaps]
    masses = Dict(i => (density ? _body_masses(snaps[i], physical) : Float64[]) for i in indices)

    ncols = min(length(indices), 3)
    nrows = ceil(Int, length(indices) / ncols)

    # Global mass colour scale across all selected snapshots
    all_m = reduce(vcat, [Float64.(snaps[i].mass) for i in indices])
    _, cmin, cmax = _log_color_range(all_m)

    for projection in projections
        ix, iy, xsym, ysym = _proj_indices(projection)
        xlab = _coord_label(xsym, unit_str)
        ylab = _coord_label(ysym, unit_str)

        # Check if adaptive zoom is needed
        extents =
            density ?
            [
                _mass_extent(
                    snaps[i].pos[ix, :] .* r_scales[i],
                    snaps[i].pos[iy, :] .* r_scales[i],
                    masses[i],
                    cfg.style.density_mass_frac,
                ) for i in indices
            ] : Float64[]
        use_adaptive =
            density ? minimum(extents) / maximum(extents) < cfg.style.zoom_frac :
            _needs_adaptive_zoom(snaps, indices, ix, iy, cfg.style.zoom_frac)

        # Panels zoomed on their own data each need their tick values. A
        # density map fills its panel, so the time annotation sits on the map
        # instead of in a data-free band above it.
        inner_ticks = use_adaptive
        band = density ? 0.0 : _MONTAGE_BAND_FRAC

        # Square panels, the shared colorbar column taken from the panel area;
        # inner tick labels show only under adaptive zoom.
        gap = _multipanel_gap(ncols; inner_ticks)
        fig = Figure(;
            size = _fig_multipanel(
                cfg,
                nrows,
                ncols;
                inner_ticks,
                panel_aspect = 1 + band,
                extra_width = _COLORBAR_WIDTH,
            ),
        )
        target_ticks = ncols > 1 ? 3 : 5
        marker_scale = _multipanel_scale(cfg, ncols; inner_ticks, extra_width = _COLORBAR_WIDTH)

        # Global limits (used when adaptive is off)
        gxlo, gxhi, gylo, gyhi = if density
            _centred_limits(maximum(extents))
        else
            all_x = reduce(vcat, [snaps[i].pos[ix, :] .* r_scales[i] for i in indices])
            all_y = reduce(vcat, [snaps[i].pos[iy, :] .* r_scales[i] for i in indices])
            _square_limits(all_x, all_y)
        end

        # Surface-density maps on one colour scale across the panels
        maps = Dict{Int,Matrix{Float64}}()
        if density
            for (k, i) in enumerate(indices)
                hw = use_adaptive ? _centred_limits(extents[k])[2] : gxhi
                maps[i] = _surface_density(
                    snaps[i].pos[ix, :] .* r_scales[i],
                    snaps[i].pos[iy, :] .* r_scales[i],
                    masses[i],
                    hw,
                    cfg.style.density_bins,
                )
            end
        end
        dhi = density ? maximum(_log_peak, values(maps)) : 0.0
        dlo = dhi - cfg.style.density_decades

        for (panel, idx) in enumerate(indices)
            snap = snaps[idx]
            px = snap.pos[ix, :] .* r_scales[idx]
            py = snap.pos[iy, :] .* r_scales[idx]
            row = div(panel - 1, ncols) + 1
            col = mod(panel - 1, ncols) + 1

            # Axis labels: only on border panels
            show_xlab = row == nrows
            show_ylab = col == 1

            # Tick labels: on every panel under adaptive zoom, only on the
            # border panels under global limits.
            show_xtick = inner_ticks || show_xlab
            show_ytick = inner_ticks || show_ylab

            # Per-panel or global limits
            if use_adaptive
                xlo, xhi, ylo, yhi =
                    density ? _centred_limits(extents[panel]) : _square_limits(px, py)
            else
                xlo, xhi, ylo, yhi = gxlo, gxhi, gylo, gyhi
            end
            xtk = _nice_ticks(xlo, xhi; target_n = target_ticks)
            ytk = _nice_ticks(ylo, yhi; target_n = target_ticks)
            # Data-free band above the data, where the time annotation sits
            yhi += band * (yhi - ylo)

            t_val = physical ? time_myr(snap.header) : time_nb(snap.header)

            ax = Axis(
                fig[row, col];
                xlabel = show_xlab ? xlab : "",
                ylabel = show_ylab ? ylab : "",
                aspect = DataAspect(),
                limits = (xlo, xhi, ylo, yhi),
                xticks = xtk,
                yticks = ytk,
                xticklabelsvisible = show_xtick,
                yticklabelsvisible = show_ytick,
                xgridvisible = false,
                ygridvisible = false,
                _density_axis_attributes(density)...,
            )

            if density
                _density_heatmap!(ax, maps[idx], xhi, dlo, dhi)
                _annotate!(ax, _time_annotation(t_val, physical); _DENSITY_ANNOTATION...)
            else
                _annotate!(ax, _time_annotation(t_val, physical))
                log_m = log10.(max.(Float64.(snap.mass), 1e-30))
                ms = max(cfg.style.marker_min, _marker_size(cfg, nparticles(snap)) * marker_scale)
                scatter!(
                    ax,
                    px,
                    py;
                    color = log_m,
                    colormap = :viridis,
                    colorrange = (cmin, cmax),
                    markersize = ms,
                    strokewidth = 0,
                )
            end
        end

        # Shared colorbar spanning all rows, right of the panel grid
        if density
            _density_colorbar!(fig[1:nrows, ncols + 1], dlo, dhi, physical)
        else
            Colorbar(
                fig[1:nrows, ncols + 1];
                colormap = :viridis,
                colorrange = (cmin, cmax),
                label = L"\log_{10}(m \, / \, M_\mathrm{tot})",
                ticks = _nice_colorbar_ticks(cmin, cmax),
            )
        end

        # Makie adds the protrusions of the inner tick labels to the gap, so
        # the compact gap serves either way; fixed boxes keep the equal-aspect
        # panels flush with their cells and the colourbar level with the grid.
        colgap!(fig.layout, _MULTIPANEL_GAP_COMPACT)
        rowgap!(fig.layout, _MULTIPANEL_GAP_COMPACT)
        colgap!(fig.layout, ncols, _COLORBAR_COLGAP)
        box_w =
            _multipanel_box_width(cfg, ncols; inner_ticks = false, extra_width = _COLORBAR_WIDTH)
        _fit_canvas_to_boxes!(fig, nrows, ncols, box_w, 1 + band)

        _save_fig(cfg, "$(filename)_$(projection)", fig)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _proj_indices(proj::Symbol)
    proj == :xy && return (1, 2, "x", "y")
    proj == :xz && return (1, 3, "x", "z")
    proj == :yz && return (2, 3, "y", "z")
    error("Unknown projection: $proj.  Use :xy, :xz, or :yz.")
end

"""Axis label for a spatial coordinate `sym` with unit string (e.g. `x [pc]`):
italic variable, upright bracketed unit."""
_coord_label(sym::AbstractString, unit::AbstractString) =
    latexstring("$(sym) \\; [\\mathrm{$(unit)}]")

"""Whether a snapshot header carries a valid physical scaling (rbar and
tscale positive).  Guards the `cfg.units == "physical"` branches against
headers without AS scaling data (fall back to N-body units)."""
_has_physical_scaling(h::SnapshotHeader) = rbar(h) > 0 && tscale(h) > 0

"""In-axis time annotation: `t = … Myr` (physical) or `t = … [NB]`."""
function _time_annotation(t_val::Real, physical::Bool)
    t_str = _fmt_latex_sig(t_val, 3)
    return physical ? latexstring("t = $(t_str)\\;\\mathrm{Myr}") :
           latexstring("t = $(t_str)\\;[\\mathrm{NB}]")
end

# ---------------------------------------------------------------------------
# Surface-density rendering
# ---------------------------------------------------------------------------

"""Perceptually uniform colormap of the surface-density maps; cells below the
colour range take its darkest colour, so empty sky reads as background."""
const _DENSITY_COLORMAP = :inferno

"""Gaussian smoothing length of the surface-density maps [grid cells]: enough
to suppress the shot noise of sparsely populated cells without widening the
core at the default resolution."""
const _DENSITY_SMOOTH_CELLS = 1.0

"""Placement of the time annotation on a density map, which has no data-free
band: white, in the top-right corner, where no tick mark can sit beside it and
read as a minus sign (ticks are drawn on the left and bottom spines only)."""
const _DENSITY_ANNOTATION = (; color = :white, corner = :tr)

"""Whether projections of `n` particles render as surface density under
`cfg.style.snapshot_render` (`"auto"`: from `density_min_n` particles on)."""
function _density_render(cfg::VisualizationConfig, n::Integer)::Bool
    mode = cfg.style.snapshot_render
    return mode == "density" || (mode == "auto" && n ≥ cfg.style.density_min_n)
end

"""Whether a figure of `snaps` is drawn in physical units: `cfg.units`, the
header scaling of every snapshot and, for a density map, the mass scale."""
function _physical_units(cfg::VisualizationConfig, snaps, density::Bool)::Bool
    cfg.units == "physical" || return false
    all(_has_physical_scaling(s.header) for s in snaps) || return false
    return !density || all(zmbar(s.header) > 0 for s in snaps)
end

"""Particle masses in M☉ (`physical`) or N-body units."""
_body_masses(snap::Snapshot, physical::Bool) =
    Float64.(snap.mass) .* (physical ? zmbar(snap.header) : 1.0)

"""
    _mass_extent(px, py, m, frac) -> Float64

Half-width of the smallest origin-centred square that holds the mass
fraction `frac`: the mass-weighted `frac` quantile of max(|x|, |y|).
"""
function _mass_extent(px, py, m, frac::Real)::Float64
    a = max.(abs.(px), abs.(py))
    isempty(a) && return 1.0
    order = sortperm(a)
    target = frac * sum(m)
    acc = 0.0
    r = Float64(a[order[end]])
    for k in order
        acc += m[k]
        if acc ≥ target
            r = Float64(a[k])
            break
        end
    end
    return r > 0 ? r : 1.0
end

"""Square limits `±half_width`, padded like [`_square_limits`](@ref)."""
function _centred_limits(half_width::Real; pad_frac::Float64 = 0.03)
    hs = half_width * (1 + pad_frac)
    return (-hs, hs, -hs, hs)
end

"""
    _surface_density(px, py, m, half_width, nbins) -> Matrix{Float64}

Projected surface mass density on an `nbins × nbins` grid over
[-half_width, half_width]², `Σ[i, j]` at the `i`-th x and `j`-th y cell: each
mass deposited in its cell and divided by the cell area, then smoothed with a
separable Gaussian of `_DENSITY_SMOOTH_CELLS` cells. Particles outside the
square (and non-finite positions) are dropped; mass smoothed across the edge
is lost, so the total is conserved for mass well inside the square.
"""
function _surface_density(px, py, m, half_width::Real, nbins::Int)::Matrix{Float64}
    h = 2 * half_width / nbins
    Σ = zeros(nbins, nbins)
    @inbounds for k in eachindex(px, py, m)
        x = (px[k] + half_width) / h
        y = (py[k] + half_width) / h
        (0 ≤ x < nbins && 0 ≤ y < nbins) || continue
        Σ[floor(Int, x) + 1, floor(Int, y) + 1] += m[k]
    end
    r = ceil(Int, 3 * _DENSITY_SMOOTH_CELLS)
    w = [exp(-0.5 * (d / _DENSITY_SMOOTH_CELLS)^2) for d in (-r):r]
    w ./= sum(w)
    smoothed_x = zeros(nbins, nbins)
    @inbounds for j in 1:nbins, i in 1:nbins, d in (-r):r
        1 ≤ i + d ≤ nbins && (smoothed_x[i, j] += w[d + r + 1] * Σ[i + d, j])
    end
    fill!(Σ, 0.0)
    @inbounds for j in 1:nbins, d in (-r):r, i in 1:nbins
        1 ≤ j + d ≤ nbins && (Σ[i, j] += w[d + r + 1] * smoothed_x[i, j + d])
    end
    return Σ ./= h^2
end

"""Top of a density colour range: log₁₀ of the largest cell value."""
_log_peak(Σ::AbstractMatrix) = log10(max(maximum(Σ), floatmin(Float64)))

"""Heatmap of `log₁₀ Σ` on the square `±half_width`, cells below `lo` in the
colormap's darkest colour."""
function _density_heatmap!(ax, Σ::AbstractMatrix, half_width::Real, lo::Real, hi::Real)
    return heatmap!(
        ax,
        _cell_centres(half_width, size(Σ, 1)),
        _cell_centres(half_width, size(Σ, 2)),
        _log_density(Σ, lo);
        colormap = _DENSITY_COLORMAP,
        colorrange = (lo, hi),
        lowclip = first(Makie.to_colormap(_DENSITY_COLORMAP)),
    )
end

_cell_centres(half_width::Real, n::Int) =
    range(-half_width + half_width / n, half_width - half_width / n; length = n)

"""`log₁₀ Σ`, empty cells one decade below `lo` so that they take the low clip."""
_log_density(Σ::AbstractMatrix, lo::Real) = log10.(max.(Σ, 10.0^(lo - 1)))

"""Axis attributes of a density panel: tick marks drawn inward stay visible on
the dark map."""
_density_axis_attributes(density::Bool) =
    density ? (; xtickcolor = :white, ytickcolor = :white) : (;)

function _density_colorbar!(position, lo::Real, hi::Real, physical::Bool)
    label =
        physical ? L"\log_{10}\,\Sigma \; [\mathrm{M_\odot \, pc^{-2}}]" :
        L"\log_{10}\,\Sigma \; [\mathrm{NB}]"
    return Colorbar(
        position;
        colormap = _DENSITY_COLORMAP,
        colorrange = (lo, hi),
        label = label,
        ticks = _nice_colorbar_ticks(lo, hi),
    )
end

"""
    _needs_adaptive_zoom(snaps, indices, ix, iy, zoom_frac) -> Bool

Return `true` if the spatial extent varies by more than `1/zoom_frac` across
the selected snapshots, meaning per-panel adaptive zoom should be used.
"""
function _needs_adaptive_zoom(
    snaps::Vector{Snapshot},
    indices,
    ix::Int,
    iy::Int,
    zoom_frac::Real,
)::Bool
    extents = Float64[]
    for idx in indices
        snap = snaps[idx]
        sx = snap.pos[ix, :]
        sy = snap.pos[iy, :]
        push!(extents, max(maximum(abs, sx), maximum(abs, sy), 0.1))
    end
    return (minimum(extents) / maximum(extents)) < zoom_frac
end
