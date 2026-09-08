# =============================================================================
# Sweep comparison figures: one series per point on common axes, coloured
# by the value of one grid axis (seeds share the colour)
# =============================================================================

"""Completed points of a sweep as `(entry, run_dir)` pairs."""
function _sweep_done_points(sweep_dir::AbstractString)
    idx = read_sweep_index(sweep_dir)
    return idx,
    [(p, joinpath(p["dir"], _SWEEP_RUN_ID)) for p in idx["points"] if p["status"] == "done"]
end

"""Colour per axis value: the Okabe–Ito cycle up to seven values, a viridis
sampling beyond."""
function _sweep_axis_colors(values::AbstractVector)
    n = length(values)
    colors = if n ≤ 7
        [_OKABE_ITO[i] for i in 1:n]
    else
        cmap = Makie.to_colormap(:viridis)
        [cmap[round(Int, 1 + (i - 1) / max(n - 1, 1) * (length(cmap) - 1))] for i in 1:n]
    end
    return Dict{Any,Any}(v => c for (v, c) in zip(values, colors))
end

"""Physical unit scaling of a run from its stdout; identity when absent."""
function _run_scaling(out_dir::AbstractString)
    path = joinpath(out_dir, "out1000")
    isfile(path) || return UnitScaling(1.0, 1.0, 1.0, 1.0)
    return extract_scaling(read_diagnostics(path))
end

"""Resolve the comparison axis: the given key, or the first grid axis."""
function _sweep_axis(idx, axis::AbstractString)
    axes = String[idx["sweep"]["axes"]...]
    isempty(axis) && return first(axes)
    axis in axes || throw(
        ArgumentError("axis \"$axis\" is not a grid axis of this sweep; axes: $(join(axes, ", "))"),
    )
    return axis
end

"""Draw one series per point, coloured by its value of `axis`, with one
legend entry per distinct value; `series(run_dir) -> (x, y)` or `nothing`."""
function _sweep_overlay!(ax, pts, axis::AbstractString, series)
    values = sort(unique(p["values"][axis] for (p, _) in pts))
    colors = _sweep_axis_colors(values)
    labelled = Set{Any}()
    n_drawn = 0
    all_x = Float64[]
    all_y = Float64[]
    for (p, run_dir) in pts
        xy = series(run_dir)
        xy === nothing && continue
        x, y = xy
        append!(all_x, x)
        append!(all_y, y)
        v = p["values"][axis]
        if v in labelled
            lines!(ax, x, y; color = (colors[v], 0.85), linewidth = 1.8)
        else
            lines!(
                ax,
                x,
                y;
                color = (colors[v], 0.85),
                linewidth = 1.8,
                label = _format_axis_value(v),
            )
            push!(labelled, v)
        end
        n_drawn += 1
    end
    return n_drawn, length(labelled), _emptiest_corner(all_x, all_y)
end

"""
    plot_sweep_lagrangian(sweep_dir, cfg::VisualizationConfig;
                          axis = "", fraction = 0.5, filename = "sweep_lagrangian") -> String

Lagrangian radius of mass fraction `fraction` against physical time for
every completed point of the sweep, coloured by the value of `axis` (the
first grid axis by default); seeds share the colour.
"""
function plot_sweep_lagrangian(
    sweep_dir::AbstractString,
    cfg::VisualizationConfig;
    axis::AbstractString = "",
    fraction::Real = 0.5,
    filename::AbstractString = "sweep_lagrangian",
)::String
    idx, pts = _sweep_done_points(sweep_dir)
    isempty(pts) && error("No completed points in sweep $sweep_dir")
    axis = _sweep_axis(idx, axis)
    pct = round(Int, 100 * fraction)

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = L"t \; [\mathrm{Myr}]",
        ylabel = latexstring("r_{$(pct)\\,\\%} \\; [\\mathrm{pc}]"),
    )
    n_drawn, n_values, corner = _sweep_overlay!(
        ax,
        pts,
        axis,
        run_dir -> begin
            out = joinpath(run_dir, "output")
            path = joinpath(out, "lagr.7")
            isfile(path) || return nothing
            lagr = read_lagr(path)
            k = argmin(abs.(lagr.mass_fractions .- fraction))
            u = _run_scaling(out)
            (to_myr(u, lagr.time), to_pc(u, lagr.radii[k, :]))
        end,
    )
    n_drawn == 0 && _no_data_note!(ax, "No Lagrangian radii in the completed runs")
    n_seeds = length(idx["sweep"]["seeds"])
    _annotate!(
        ax,
        "$(n_drawn) runs, $(n_seeds) seed$(n_seeds == 1 ? "" : "s") per point";
        corner = corner,
    )
    n_values ≥ 2 && _top_legend!(fig, ax; title = _axis_short(axis) * ":")
    return _save_fig(cfg, filename, fig)
end

"""
    plot_sweep_energy(sweep_dir, cfg::VisualizationConfig;
                      axis = "", filename = "sweep_energy") -> String

Relative energy error |ΔE/E| per adjustment against physical time for
every completed point, on a logarithmic axis, coloured by the value of
`axis`.
"""
function plot_sweep_energy(
    sweep_dir::AbstractString,
    cfg::VisualizationConfig;
    axis::AbstractString = "",
    filename::AbstractString = "sweep_energy",
)::String
    idx, pts = _sweep_done_points(sweep_dir)
    isempty(pts) && error("No completed points in sweep $sweep_dir")
    axis = _sweep_axis(idx, axis)

    fig = Figure(; size = _figsize_px(cfg))
    ax =
        Axis(fig[1, 1]; xlabel = L"t \; [\mathrm{Myr}]", ylabel = L"|\Delta E / E|", yscale = log10)
    lo, hi = Inf, 0.0
    n_drawn, n_values, _ = _sweep_overlay!(
        ax,
        pts,
        axis,
        run_dir -> begin
            path = joinpath(run_dir, "output", "out1000")
            isfile(path) || return nothing
            adj = read_diagnostics(path).adjust
            keep = [a for a in adj if isfinite(a.de_rel) && a.de_rel != 0]
            isempty(keep) && return nothing
            y = abs.([a.de_rel for a in keep])
            lo = min(lo, minimum(y))
            hi = max(hi, maximum(y))
            ([a.time_myr for a in keep], y)
        end,
    )
    if n_drawn == 0
        _no_data_note!(ax, "No ADJUST records in the completed runs")
    else
        ax.yticks = _log_ticks(lo, hi)
        ylims!(ax, lo / 2, hi * 2)
    end
    n_values ≥ 2 && _top_legend!(fig, ax; title = _axis_short(axis) * ":")
    return _save_fig(cfg, filename, fig)
end

"""
    sweep_figures(sweep_dir, cfg::VisualizationConfig; axis = "") -> Vector{String}

The comparison figures of a sweep ([`plot_sweep_lagrangian`](@ref),
[`plot_sweep_energy`](@ref)) for one grid axis.
"""
function sweep_figures(
    sweep_dir::AbstractString,
    cfg::VisualizationConfig;
    axis::AbstractString = "",
)
    return [
        plot_sweep_lagrangian(sweep_dir, cfg; axis = axis),
        plot_sweep_energy(sweep_dir, cfg; axis = axis),
    ]
end
