# =============================================================================
# Lagrangian radii evolution plot
# =============================================================================

"""
    plot_lagrangian(lagr::LagrangianData, cfg::VisualizationConfig;
                    filename::AbstractString = "lagrangian_radii",
                    selected_fractions::Vector{Float64} = Float64[],
                    units::Union{UnitScaling,Nothing} = nothing)

Plot Lagrangian radii as a function of time.
If `selected_fractions` is empty, a sensible default subset is plotted.
When `cfg.units == "physical"` and `units` is provided, times are shown in
Myr and radii in pc; N-body units otherwise (lagr.7 carries no header, so
the scaling must be supplied by the caller).
"""
function plot_lagrangian(
    lagr::LagrangianData, cfg::VisualizationConfig;
    filename::AbstractString = "lagrangian_radii",
    selected_fractions::Vector{Float64} = Float64[],
    units::Union{UnitScaling,Nothing} = nothing,
)
    isempty(lagr.time) && (@warn "No Lagrangian data to plot"; return nothing)

    # Default selection: 1%, 10%, 50%, 90%, 100%
    if isempty(selected_fractions)
        selected_fractions = [0.01, 0.1, 0.5, 0.9, 1.0]
    end

    physical = cfg.units == "physical" && units !== nothing &&
               units.rbar > 0 && units.tscale > 0
    ts       = physical ? lagr.time .* units.tscale : lagr.time
    r_scale  = physical ? units.rbar : 1.0

    # Closest available mass fractions and the plotted (positive, scaled)
    # radii — the log-axis tick range comes from the actual data extents.
    sel_idx = [argmin(abs.(lagr.mass_fractions .- f)) for f in selected_fractions]
    r_pos = [r * r_scale for idx in sel_idx
             for r in @view(lagr.radii[idx, :]) if r > 0]

    fig = Figure(; size = _figsize_px(cfg))
    ttk = _time_ticks(first(ts), last(ts))
    ax = Axis(fig[1, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = physical ? L"r_\mathrm{L} \; [\mathrm{pc}]" : L"r_\mathrm{L} \; [\mathrm{NB}]",
        yscale = log10,
        xticks = ttk,
        yticks = isempty(r_pos) ? Makie.automatic : _log_ticks(extrema(r_pos)...),
    )

    for (ci, idx) in enumerate(sel_idx)
        actual_frac = lagr.mass_fractions[idx]

        # Format as clean integer percentage where possible
        pct_val = actual_frac * 100
        pct = isinteger(pct_val) ? @sprintf("%d", Int(pct_val)) : @sprintf("%.1f", pct_val)

        # Mask non-positive radii (empty shells at early times) — they are
        # invalid on the log axis; NaN points are skipped by Makie.
        ys = [r > 0 ? r * r_scale : NaN for r in @view lagr.radii[idx, :]]
        lines!(ax, ts, ys;
            label = latexstring("$(pct)\\%"),
            color = _OKABE_ITO[mod1(ci, length(_OKABE_ITO))],
        )
    end

    if length(selected_fractions) ≥ 2
        _top_legend!(fig, ax; title = L"M(r)/M_\mathrm{tot}:")
    end

    _save_fig(cfg, filename, fig)
    return nothing
end
