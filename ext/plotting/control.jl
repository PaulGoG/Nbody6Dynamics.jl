# =============================================================================
# Merger against control: paired series of every sweep point and its
# isolated single-cluster companion
# =============================================================================

"""
    plot_control_comparison(sweep_dir, cfg::VisualizationConfig;
                            quantity = :lagrangian, axis = "", fraction = 0.5,
                            filename = "control_<quantity>") -> String

For every completed merger point with a completed control, the series of
`quantity` (`:lagrangian`, `:energy`, `:n_stars`, `:n_pairs`) of the merger
(solid) and of its control (dashed) in the colour of the point's value of
`axis`; seeds share the colour. Legend families: the axis values and the
line roles.
"""
@publication function Nbody6Dynamics.plot_control_comparison(
    sweep_dir::AbstractString,
    cfg::VisualizationConfig;
    quantity::Symbol = :lagrangian,
    axis::AbstractString = "",
    fraction::Real = 0.5,
    filename::AbstractString = "control_" * String(quantity),
)::String
    idx, mergers = _sweep_done_points(sweep_dir; kind = "merger")
    _, controls = _sweep_done_points(sweep_dir; kind = "control")
    control_of = Dict(p["control_of"] => run_dir for (p, run_dir) in controls)
    pairs = [
        (p, run_dir, control_of[p["id"]]) for (p, run_dir) in mergers if haskey(control_of, p["id"])
    ]
    isempty(pairs) && error("No merger point with a completed control in sweep $sweep_dir")
    axis = _sweep_axis(idx, axis)
    values = sort(
        unique(_point_axis_value(p, axis) for (p, _, _) in pairs);
        by = v -> v === nothing ? 0 : v,
    )
    colors = _sweep_axis_colors(values)

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = L"t \; [\mathrm{Myr}]",
        ylabel = _series_label(quantity, fraction),
        yscale = quantity === :energy ? log10 : identity,
    )
    lo, hi = Inf, 0.0
    n_pairs = 0
    for (p, run_m, run_c) in pairs
        sm = _run_series(run_m, quantity; fraction = fraction)
        sc = _run_series(run_c, quantity; fraction = fraction)
        (sm === nothing || sc === nothing) && continue
        c = colors[_point_axis_value(p, axis)]
        lines!(ax, sm[1], sm[2]; color = (c, 0.9), linewidth = _STYLE.data)
        lines!(ax, sc[1], sc[2]; color = (c, 0.9), linewidth = _STYLE.data, linestyle = :dash)
        lo = min(lo, minimum(sm[2]), minimum(sc[2]))
        hi = max(hi, maximum(sm[2]), maximum(sc[2]))
        n_pairs += 1
    end
    n_pairs == 0 && _no_data_note!(ax, "No $(quantity) series in the paired runs")
    if quantity === :energy && n_pairs > 0 && lo > 0
        ax.yticks = _log_ticks(lo, hi)
        ylims!(ax, lo / 2, hi * 2)
    end
    _annotate!(ax, "$(n_pairs) merger–control pairs"; corner = :tl)

    entries = Vector{Vector{_LegendElement}}()
    labels = Vector{Vector{AbstractString}}()
    titles = String[]
    if length(values) ≥ 2 && !isempty(axis)
        push!(
            entries,
            _LegendElement[
                LineElement(; color = colors[v], linewidth = _STYLE.data) for v in values
            ],
        )
        push!(labels, AbstractString[_format_axis_value(v) for v in values])
        push!(titles, _axis_short(axis) * ":")
    end
    grey = _OKABE_ITO[8]
    push!(
        entries,
        _LegendElement[
            LineElement(; color = grey, linewidth = _STYLE.data),
            LineElement(; color = grey, linewidth = _STYLE.data, linestyle = :dash),
        ],
    )
    # "Control" rather than "Isolated control": with a grid-axis family beside
    # it the longer label pushes the horizontal legend past the figure width and
    # is clipped. The caption carries what the control is.
    push!(labels, AbstractString["Merger", "Control"])
    push!(titles, "Run:")
    Legend(
        fig[0, :],
        entries,
        labels,
        titles;
        orientation = :horizontal,
        # Two families on one row overrun the figure width; banking each of them
        # keeps the legend inside it.
        nbanks = (length(titles) ≥ 2 || length(values) > 4) ? 2 : 1,
        framevisible = false,
        titleposition = :left,
        tellheight = true,
        padding = (0, 0, 0, 0),
    )
    return _save_fig(cfg, filename, fig)
end
