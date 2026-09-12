# =============================================================================
# Snapshot visualisation — 2D projections and density maps
# =============================================================================

"""
    plot_snapshot(snap::Snapshot, cfg::VisualizationConfig;
                  filename::AbstractString = "snapshot",
                  projections::Vector{Symbol} = [:xy, :xz])

Generate publication-quality scatter plots of particle positions in the
requested projections (:xy, :xz, :yz).  Positions are shown in pc and the
time annotation in Myr when `cfg.units == "physical"` (header AS scaling);
N-body units otherwise.
"""
function Nbody6Dynamics.plot_snapshot(
    snap::Snapshot,
    cfg::VisualizationConfig;
    filename::AbstractString = "snapshot",
    projections::Vector{Symbol} = [:xy, :xz, :yz],
)
    physical = cfg.units == "physical" && _has_physical_scaling(snap.header)
    unit_str = physical ? "pc" : "NB"
    r_scale = physical ? rbar(snap.header) : 1.0
    t_val = physical ? time_myr(snap.header) : time_nb(snap.header)
    n = nparticles(snap)

    ms = _marker_size(cfg, n)

    # Colour by mass (log scale for visual contrast)
    log_m, cmin, cmax = _log_color_range(Float64.(snap.mass))

    for proj in projections
        fig = Figure(; size = _fig_with_colorbar(cfg))

        ix, iy, xsym, ysym = _proj_indices(proj)
        px = snap.pos[ix, :] .* r_scale
        py = snap.pos[iy, :] .* r_scale

        # Square limits and nice ticks for consistent box sizes
        xlo, xhi, ylo, yhi = _square_limits(px, py)
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
        )
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

        colgap!(fig.layout, _COLORBAR_COLGAP)

        _save_fig(cfg, "$(filename)_$(proj)", fig)
    end

    return nothing
end

"""
    plot_snapshot_evolution(snaps::Vector{Snapshot}, cfg::VisualizationConfig;
                           filename::AbstractString = "snapshot_evolution",
                           projections::Vector{Symbol} = [:xy, :xz, :yz],
                           max_panels::Int = 6)

Multi-panel figure showing cluster evolution across snapshots.  Positions
are shown in pc and time annotations in Myr when `cfg.units == "physical"`.

When the spatial extent varies dramatically across panels (e.g. merger runs),
each panel is zoomed on its own data.  Panels then need their own tick values,
which fit only in a grid of at most `_MAX_TICKED_MONTAGE_COLS` columns; in a
wider grid the panels carry no tick values and each annotation states that
panel's half-width instead, so the scale is still readable and the labels of
adjacent panels cannot run into each other.  Under global limits (the extent is
comparable throughout) only the border panels are labelled.
"""
function Nbody6Dynamics.plot_snapshot_evolution(
    snaps::Vector{Snapshot},
    cfg::VisualizationConfig;
    filename::AbstractString = "snapshot_evolution",
    projections::Vector{Symbol} = [:xy, :xz, :yz],
    max_panels::Int = 6,
)
    ns = length(snaps)
    ns == 0 && return nothing

    physical = cfg.units == "physical" && all(_has_physical_scaling(s.header) for s in snaps)
    unit_str = physical ? "pc" : "NB"
    r_scales = [physical ? rbar(s.header) : 1.0 for s in snaps]

    # Select evenly-spaced snapshots
    indices = ns <= max_panels ? (1:ns) : round.(Int, range(1, ns; length = max_panels))

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
        use_adaptive = _needs_adaptive_zoom(snaps, indices, ix, iy, cfg.style.zoom_frac)

        # Per-panel zoom means every panel needs its own tick values, and three
        # of those do not fit across one column width: the labels of adjacent
        # panels run into each other and the time annotation is clipped. Beyond
        # _MAX_TICKED_MONTAGE_COLS the panels therefore carry no tick values and
        # the annotation states each panel's own half-width instead.
        inner_ticks = use_adaptive && ncols <= _MAX_TICKED_MONTAGE_COLS
        scale_in_annotation = use_adaptive && !inner_ticks
        # The half-width line needs a taller data-free band than the time alone.
        band = scale_in_annotation ? _MONTAGE_BAND_FRAC_TWO_LINE : _MONTAGE_BAND_FRAC

        # One column wide, square panels, the shared colorbar column taken
        # from the panel area; inner tick labels show only under adaptive zoom.
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
        all_x = reduce(vcat, [snaps[i].pos[ix, :] .* r_scales[i] for i in indices])
        all_y = reduce(vcat, [snaps[i].pos[iy, :] .* r_scales[i] for i in indices])
        gxlo, gxhi, gylo, gyhi = _square_limits(all_x, all_y)

        for (panel, idx) in enumerate(indices)
            snap = snaps[idx]
            px = snap.pos[ix, :] .* r_scales[idx]
            py = snap.pos[iy, :] .* r_scales[idx]
            row = div(panel - 1, ncols) + 1
            col = mod(panel - 1, ncols) + 1

            # Axis labels: only on border panels
            show_xlab = row == nrows
            show_ylab = col == 1

            # Tick labels: per panel when adaptive zoom can afford them, only on
            # border panels under global limits, and nowhere when the panels are
            # zoomed individually but too narrow to label (the border panel's
            # values would then be read as everyone's).
            show_xtick = scale_in_annotation ? false : (inner_ticks || show_xlab)
            show_ytick = scale_in_annotation ? false : (inner_ticks || show_ylab)

            # Per-panel or global limits
            if use_adaptive
                xlo, xhi, ylo, yhi = _square_limits(px, py)
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
                xlabelsize = 22,
                ylabelsize = 22,
                xticklabelsize = 18,
                yticklabelsize = 18,
                aspect = DataAspect(),
                limits = (xlo, xhi, ylo, yhi),
                xticks = xtk,
                yticks = ytk,
                xticklabelsvisible = show_xtick,
                yticklabelsvisible = show_ytick,
                xgridvisible = false,
                ygridvisible = false,
            )
            _annotate!(ax, _time_annotation(t_val, physical))
            # xhi is the half-width: _square_limits is symmetric about zero.
            scale_in_annotation &&
                _annotate!(ax, _extent_annotation(xhi, unit_str); dy = _MONTAGE_ANNOTATION_DY)

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

        # Shared colorbar spanning all rows, right of the panel grid
        Colorbar(
            fig[1:nrows, ncols + 1];
            colormap = :viridis,
            colorrange = (cmin, cmax),
            label = L"\log_{10}(m \, / \, M_\mathrm{tot})",
            ticks = _nice_colorbar_ticks(cmin, cmax),
        )

        colgap!(fig.layout, gap)
        rowgap!(fig.layout, gap)

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
    t_str = @sprintf("%.3g", t_val)
    return physical ? latexstring("t = $(t_str)\\;\\mathrm{Myr}") :
           latexstring("t = $(t_str)\\;[\\mathrm{NB}]")
end

"""
Half-extent annotation (`± 42 pc`) of a montage panel zoomed on its own data.
It sits under the time annotation and replaces the tick values, which such a
panel is too narrow to carry; two significant digits keep it inside the panel.
"""
_extent_annotation(half_width::Real, unit::AbstractString) =
    latexstring("\\pm $(@sprintf("%.2g", half_width))\\;\\mathrm{$(unit)}")

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
