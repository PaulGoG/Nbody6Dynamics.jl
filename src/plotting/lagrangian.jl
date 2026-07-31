# =============================================================================
# Lagrangian radii evolution plot
# =============================================================================

"""
    plot_lagrangian(lagr::LagrangianData, cfg::VisualizationConfig;
                    filename::AbstractString = "lagrangian_radii",
                    selected_fractions::Vector{Float64} = Float64[])

Plot Lagrangian radii as a function of time.
If `selected_fractions` is empty, a sensible default subset is plotted.
"""
function plot_lagrangian(
    lagr::LagrangianData, cfg::VisualizationConfig;
    filename::AbstractString = "lagrangian_radii",
    selected_fractions::Vector{Float64} = Float64[],
)
    isempty(lagr.time) && (@warn "No Lagrangian data to plot"; return nothing)

    # Default selection: 1%, 10%, 50%, 90%, 100%
    if isempty(selected_fractions)
        selected_fractions = [0.01, 0.1, 0.5, 0.9, 1.0]
    end

    fig = Figure(; size = _figsize_px(cfg))
    ttk = _time_ticks(first(lagr.time), last(lagr.time))
    ax = Axis(fig[1, 1];
        xlabel = L"\mathrm{t} \; \mathrm{[NB]}",
        ylabel = L"\mathrm{r}_\mathrm{L} \; \mathrm{[NB]}",
        yscale = log10,
        xticks = ttk,
    )

    colors = Makie.wong_colors()

    ci = 1
    for frac in selected_fractions
        # Find closest available mass fraction
        idx = argmin(abs.(lagr.mass_fractions .- frac))
        actual_frac = lagr.mass_fractions[idx]

        # Format as clean integer percentage where possible
        pct_val = actual_frac * 100
        pct = isinteger(pct_val) ? @sprintf("%d", Int(pct_val)) : @sprintf("%.1f", pct_val)
        label = latexstring("\\mathrm{M}(\\mathrm{r})/\\mathrm{M}_\\mathrm{tot} = $(pct)\\%")

        # Mask non-positive radii (empty shells at early times) — they are
        # invalid on the log axis; NaN points are skipped by Makie.
        ys = [r > 0 ? r : NaN for r in @view lagr.radii[idx, :]]
        lines!(ax, lagr.time, ys;
            label = label,
            color = colors[mod1(ci, length(colors))],
        )
        ci += 1
    end

    if length(selected_fractions) ≥ 2
        axislegend(ax; position = :rt, framevisible = true,
                   backgroundcolor = (:white, 0.7))
    end

    _save_fig(cfg, filename, fig)
    return nothing
end
