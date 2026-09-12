# =============================================================================
# Lagrangian radii evolution plot
# =============================================================================

"""Default mass-fraction subset for Lagrangian figures: 1%, 10%, 50%, 90%, 100%."""
const _LAGR_DEFAULT_FRACTIONS = [0.01, 0.1, 0.5, 0.9, 1.0]

"""
    _lagrangian_axis(fig, lagr, cfg, selected_fractions, units)
        -> (ax, ts, sel_idx, r_scale)

Shared prologue of the Lagrangian plot and animation: resolves units
(Myr/pc when `cfg.units == "physical"` and a `UnitScaling` is supplied,
N-body otherwise), selects the closest available mass fractions, and
constructs the log-radius axis with decade ticks from the positive scaled
radii.
"""
function _lagrangian_axis(
    fig,
    lagr::LagrangianData,
    cfg::VisualizationConfig,
    selected_fractions::Vector{Float64},
    units::Union{UnitScaling,Nothing},
)
    physical = cfg.units == "physical" && units !== nothing && units.rbar > 0 && units.tscale > 0
    ts = physical ? lagr.time .* units.tscale : lagr.time
    r_scale = physical ? units.rbar : 1.0

    # Closest available mass fractions and the plotted (positive, scaled)
    # radii — the log-axis tick range comes from the actual data extents.
    sel_idx = [argmin(abs.(lagr.mass_fractions .- f)) for f in selected_fractions]
    r_pos = [r * r_scale for idx in sel_idx for r in @view(lagr.radii[idx, :]) if r > 0]

    ax = Axis(
        fig[1, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = physical ? L"r_\mathrm{L} \; [\mathrm{pc}]" : L"r_\mathrm{L} \; [\mathrm{NB}]",
        yscale = log10,
        xticks = _time_ticks(first(ts), last(ts)),
        yticks = isempty(r_pos) ? Makie.automatic : _log_ticks(extrema(r_pos)...),
    )
    return ax, ts, sel_idx, r_scale
end

"""Percentage legend label for a mass fraction (clean integer when exact)."""
function _fraction_pct_label(frac::Real)
    pct_val = frac * 100
    pct = isinteger(pct_val) ? @sprintf("%d", Int(pct_val)) : @sprintf("%.1f", pct_val)
    return latexstring("$(pct)\\%")
end

"""Radii row `idx` scaled by `r_scale`, with non-positive values (empty
shells at early times) masked to NaN — invalid on the log axis."""
_masked_radii(lagr::LagrangianData, idx::Int, r_scale::Real) =
    [r > 0 ? r * r_scale : NaN for r in @view lagr.radii[idx, :]]

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
function Nbody6Dynamics.plot_lagrangian(
    lagr::LagrangianData,
    cfg::VisualizationConfig;
    filename::AbstractString = "lagrangian_radii",
    selected_fractions::Vector{Float64} = Float64[],
    units::Union{UnitScaling,Nothing} = nothing,
)
    isempty(lagr.time) && (@warn "No Lagrangian data to plot"; return nothing)

    if isempty(selected_fractions)
        selected_fractions = _LAGR_DEFAULT_FRACTIONS
    end

    fig = Figure(; size = _figsize_px(cfg))
    ax, ts, sel_idx, r_scale = _lagrangian_axis(fig, lagr, cfg, selected_fractions, units)

    for (ci, idx) in enumerate(sel_idx)
        lines!(
            ax,
            ts,
            _masked_radii(lagr, idx, r_scale);
            label = _fraction_pct_label(lagr.mass_fractions[idx]),
            color = _OKABE_ITO[mod1(ci, length(_OKABE_ITO))],
        )
    end

    if length(selected_fractions) ≥ 2
        _top_legend!(fig, ax; title = L"M(r)/M_\mathrm{tot}:")
    end

    _save_fig(cfg, filename, fig)
    return nothing
end
