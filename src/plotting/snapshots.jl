# =============================================================================
# Snapshot visualisation — 2D projections and density maps
# =============================================================================

"""
    plot_snapshot(snap::Snapshot, cfg::VisualizationConfig;
                  filename::AbstractString = "snapshot",
                  projections::Vector{Symbol} = [:xy, :xz])

Generate publication-quality scatter plots of particle positions in the
requested projections (:xy, :xz, :yz).
"""
function plot_snapshot(
    snap::Snapshot, cfg::VisualizationConfig;
    filename::AbstractString = "snapshot",
    projections::Vector{Symbol} = [:xy, :xz, :yz],
)
    t_nb = time_nb(snap.header)
    n = nparticles(snap)

    # Particle size: visible but not overlapping
    ms = clamp(18000 / n, 4.0, 20.0)

    # Colour by mass (log scale for visual contrast)
    m = Float64.(snap.mass)
    log_m = log10.(max.(m, 1e-30))
    cmin, cmax = extrema(log_m)
    if cmin ≈ cmax; cmin -= 0.5; cmax += 0.5; end

    for proj in projections
        fig = Figure(; size = _fig_with_colorbar(cfg))

        ix, iy, xlabel, ylabel = _proj_indices(proj)

        # Square limits and nice ticks for consistent box sizes
        xlo, xhi, ylo, yhi = _square_limits(snap.pos[ix, :], snap.pos[iy, :])
        xtk = _nice_ticks(xlo, xhi)
        ytk = _nice_ticks(ylo, yhi)

        t_str = @sprintf("%.4f", t_nb)
        ax = Axis(fig[1, 1];
            xlabel = xlabel,
            ylabel = ylabel,
            title  = latexstring("\\mathrm{N} = $(n), \\;\\; \\mathrm{t} = $(t_str) \\; \\mathrm{[NB]}"),
            aspect = DataAspect(),
            limits = (xlo, xhi, ylo, yhi),
            xticks = xtk,
            yticks = ytk,
        )

        sc = scatter!(ax, snap.pos[ix, :], snap.pos[iy, :];
            color      = log_m,
            colormap   = :viridis,
            colorrange = (cmin, cmax),
            markersize = ms,
            strokewidth = 0,
        )

        Colorbar(fig[1, 2], sc;
            label = L"\log_{10}(\mathrm{m} \, / \, \mathrm{M}_\mathrm{tot})",
            ticks = _nice_colorbar_ticks(cmin, cmax),
        )

        colgap!(fig.layout, 10)

        save(_output_path(cfg, "$(filename)_$(proj)"), fig; px_per_unit = cfg.dpi / 72)
    end

    return nothing
end

"""
    plot_snapshot_evolution(snaps::Vector{Snapshot}, cfg::VisualizationConfig;
                           filename::AbstractString = "snapshot_evolution",
                           projections::Vector{Symbol} = [:xy, :xz, :yz],
                           max_panels::Int = 6)

Multi-panel figure showing cluster evolution across snapshots.
When the spatial extent varies dramatically across panels (e.g. merger runs),
per-panel adaptive zoom is enabled and tick labels are shown on every panel
so the reader can infer the scale from the axis values.
"""
function plot_snapshot_evolution(
    snaps::Vector{Snapshot}, cfg::VisualizationConfig;
    filename::AbstractString = "snapshot_evolution",
    projections::Vector{Symbol} = [:xy, :xz, :yz],
    max_panels::Int = 6,
)
    ns = length(snaps)
    ns == 0 && return nothing

    # Select evenly-spaced snapshots
    indices = ns <= max_panels ? (1:ns) : round.(Int, range(1, ns; length = max_panels))

    ncols = min(length(indices), 3)
    nrows = ceil(Int, length(indices) / ncols)

    # Global mass colour scale across all selected snapshots
    all_m = reduce(vcat, [Float64.(snaps[i].mass) for i in indices])
    log_m_global = log10.(max.(all_m, 1e-30))
    cmin, cmax = extrema(log_m_global)
    if cmin ≈ cmax; cmin -= 0.5; cmax += 0.5; end

    for projection in projections
        # Extra width for the shared colorbar column
        pw, ph = _fig_multipanel(cfg, nrows, ncols)
        fig = Figure(; size = (pw + 110, ph))

        ix, iy, xlab, ylab = _proj_indices(projection)

        # Check if adaptive zoom is needed
        use_adaptive = _needs_adaptive_zoom(snaps, indices, ix, iy)

        # Global limits (used when adaptive is off)
        all_x = reduce(vcat, [snaps[i].pos[ix, :] for i in indices])
        all_y = reduce(vcat, [snaps[i].pos[iy, :] for i in indices])
        gxlo, gxhi, gylo, gyhi = _square_limits(all_x, all_y)

        for (panel, idx) in enumerate(indices)
            snap = snaps[idx]
            row = div(panel - 1, ncols) + 1
            col = mod(panel - 1, ncols) + 1

            # Axis labels: only on border panels
            show_xlab = row == nrows
            show_ylab = col == 1

            # Tick labels: always visible when adaptive (scales differ),
            # only on border panels otherwise
            show_xtick = use_adaptive || show_xlab
            show_ytick = use_adaptive || show_ylab

            # Per-panel or global limits
            if use_adaptive
                xlo, xhi, ylo, yhi = _square_limits(snap.pos[ix, :], snap.pos[iy, :])
            else
                xlo, xhi, ylo, yhi = gxlo, gxhi, gylo, gyhi
            end
            xtk = _nice_ticks(xlo, xhi)
            ytk = _nice_ticks(ylo, yhi)

            t_str = @sprintf("%.3f", time_nb(snap.header))
            panel_title = latexstring("\\mathrm{t} = $(t_str) \\; \\mathrm{[NB]}")

            ax = Axis(fig[row, col];
                xlabel = show_xlab ? xlab : "",
                ylabel = show_ylab ? ylab : "",
                title  = panel_title,
                titlesize = 22,
                xlabelsize = 22,
                ylabelsize = 22,
                xticklabelsize = use_adaptive ? 15 : 18,
                yticklabelsize = use_adaptive ? 15 : 18,
                aspect = DataAspect(),
                limits = (xlo, xhi, ylo, yhi),
                xticks = xtk,
                yticks = ytk,
                xticklabelsvisible = show_xtick,
                yticklabelsvisible = show_ytick,
            )

            log_m = log10.(max.(Float64.(snap.mass), 1e-30))
            ms = clamp(18000 / nparticles(snap), 4.0, 18.0)
            scatter!(ax, snap.pos[ix, :], snap.pos[iy, :];
                color      = log_m,
                colormap   = :viridis,
                colorrange = (cmin, cmax),
                markersize = ms,
                strokewidth = 0,
            )
        end

        # Shared colorbar spanning all rows, right of the panel grid
        Colorbar(fig[1:nrows, ncols + 1];
            colormap   = :viridis,
            colorrange = (cmin, cmax),
            label      = L"\log_{10}(\mathrm{m} \, / \, \mathrm{M}_\mathrm{tot})",
            ticks      = _nice_colorbar_ticks(cmin, cmax),
        )

        colgap!(fig.layout, _MULTIPANEL_HGAP)
        rowgap!(fig.layout, _MULTIPANEL_VGAP)

        save(_output_path(cfg, "$(filename)_$(projection)"), fig; px_per_unit = cfg.dpi / 72)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _proj_indices(proj::Symbol)
    proj == :xy && return (1, 2, L"\mathrm{x} \; \mathrm{[NB]}", L"\mathrm{y} \; \mathrm{[NB]}")
    proj == :xz && return (1, 3, L"\mathrm{x} \; \mathrm{[NB]}", L"\mathrm{z} \; \mathrm{[NB]}")
    proj == :yz && return (2, 3, L"\mathrm{y} \; \mathrm{[NB]}", L"\mathrm{z} \; \mathrm{[NB]}")
    error("Unknown projection: $proj.  Use :xy, :xz, or :yz.")
end

"""
    _needs_adaptive_zoom(snaps, indices, ix, iy) -> Bool

Return `true` if the spatial extent varies by more than ~7× across the
selected snapshots, meaning per-panel adaptive zoom should be used.
"""
function _needs_adaptive_zoom(snaps::Vector{Snapshot}, indices, ix::Int, iy::Int)::Bool
    extents = Float64[]
    for idx in indices
        snap = snaps[idx]
        sx = snap.pos[ix, :]; sy = snap.pos[iy, :]
        push!(extents, max(maximum(abs, sx), maximum(abs, sy), 0.1))
    end
    return (minimum(extents) / maximum(extents)) < 0.15
end
