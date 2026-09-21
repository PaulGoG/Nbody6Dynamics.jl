# =============================================================================
# Sweep comparison figures: one series per point on common axes, coloured
# by the value of one grid axis (seeds share the colour)
# =============================================================================

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

"""Resolve the comparison axis: the given key, the first grid axis, or `""`
for a sweep without axes (seeds only)."""
function _sweep_axis(idx, axis::AbstractString)
    axes = String[idx["sweep"]["axes"]...]
    isempty(axes) && isempty(axis) && return ""
    isempty(axes) && throw(ArgumentError("axis \"$axis\" given, but this sweep has no grid axes"))
    isempty(axis) && return first(axes)
    axis in axes || throw(
        ArgumentError("axis \"$axis\" is not a grid axis of this sweep; axes: $(join(axes, ", "))"),
    )
    return axis
end

"""Value of `axis` at a point; `nothing` for a sweep without axes."""
_point_axis_value(p, axis::AbstractString) = isempty(axis) ? nothing : p["values"][axis]

"""Draw one series per point, coloured by its value of `axis`, with one
legend entry per distinct value; `series(run_dir) -> (x, y)` or `nothing`."""
function _sweep_overlay!(ax, pts, axis::AbstractString, series)
    values =
        sort(unique(_point_axis_value(p, axis) for (p, _) in pts); by = v -> v === nothing ? 0 : v)
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
        v = _point_axis_value(p, axis)
        if v === nothing || v in labelled
            lines!(ax, x, y; color = (colors[v], 0.85), linewidth = _STYLE.data)
        else
            lines!(
                ax,
                x,
                y;
                color = (colors[v], 0.85),
                linewidth = _STYLE.data,
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
@publication function Nbody6Dynamics.plot_sweep_lagrangian(
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
        run_dir -> _run_series(run_dir, :lagrangian; fraction = fraction),
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
@publication function Nbody6Dynamics.plot_sweep_energy(
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
    n_drawn, n_values, _ = _sweep_overlay!(ax, pts, axis, run_dir -> begin
        s = _run_series(run_dir, :energy)
        s === nothing && return nothing
        lo = min(lo, minimum(s[2]))
        hi = max(hi, maximum(s[2]))
        s
    end)
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

The comparison figures of a sweep for one grid axis
([`plot_sweep_lagrangian`](@ref), [`plot_sweep_energy`](@ref)) and, when
the sweep has more than one seed, the ensemble figures of the same
quantities ([`plot_sweep_ensemble`](@ref)), and the merger–control pairs
when the sweep carries controls ([`plot_control_comparison`](@ref)).
"""
function Nbody6Dynamics.sweep_figures(
    sweep_dir::AbstractString,
    cfg::VisualizationConfig;
    axis::AbstractString = "",
)
    paths = [
        plot_sweep_lagrangian(sweep_dir, cfg; axis = axis),
        plot_sweep_energy(sweep_dir, cfg; axis = axis),
    ]
    idx = read_sweep_index(sweep_dir)
    if length(idx["sweep"]["seeds"]) ≥ 2
        push!(paths, plot_sweep_ensemble(sweep_dir, cfg; quantity = :lagrangian, axis = axis))
        push!(paths, plot_sweep_ensemble(sweep_dir, cfg; quantity = :energy, axis = axis))
    end
    if get(idx["sweep"], "controls", false)
        push!(paths, plot_control_comparison(sweep_dir, cfg; quantity = :lagrangian, axis = axis))
    end
    return paths
end
