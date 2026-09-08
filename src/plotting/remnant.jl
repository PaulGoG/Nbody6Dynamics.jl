# =============================================================================
# Remnant figures: rotation, structure, mass segregation
# =============================================================================

"""Time axis of a `RemnantDiagnostics` in the configured units: values and label."""
function _remnant_time(diag::RemnantDiagnostics, cfg::VisualizationConfig)
    physical = cfg.units == "physical" && all(isfinite, diag.time_myr) && any(!=(1.0), diag.rbar)
    t = physical ? diag.time_myr : diag.time
    label = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]"
    t_coal = physical ? diag.coalescence_time_myr : diag.coalescence_time
    t_seg = physical ? diag.segregation_time_myr : diag.segregation_time
    return (; physical, t, label, t_coal, t_seg)
end

"""Occupancy points of a series extended by one point standing for an event
marker's label (`x` the marker time, `y` the top of the data range), so
[`_emptiest_corner`](@ref) keeps annotations away from it. Returns the
series unchanged when the marker is absent."""
_with_marker(v::AbstractVector{<:Real}, extra::Real) =
    isfinite(extra) ? vcat(Float64.(v), Float64(extra)) : Float64.(v)

"""Fraction of the time range within which a marker label flips to the left of its line."""
const _MARKER_EDGE_FRACTION = 0.2

"""Dashed vertical marker at `x` with `label` hanging from `y_top` (data
coordinates, the caller's upper axis limit), on the right of the line
unless `x` lies in the last `_MARKER_EDGE_FRACTION` of `t_range`; no-op
for `NaN`."""
function _event_marker!(
    ax,
    x::Real,
    label::AbstractString,
    y_top::Real,
    t_range::Tuple{<:Real,<:Real};
    color = :black,
)
    isfinite(x) || return nothing
    vlines!(ax, [x]; color = color, linestyle = :dash, linewidth = 1.5)
    span = max(t_range[2] - t_range[1], eps(Float64))
    right_side = (x - t_range[1]) / span < 1 - _MARKER_EDGE_FRACTION
    text!(
        ax,
        x,
        y_top;
        text = label,
        align = (right_side ? :left : :right, :top),
        offset = (right_side ? 4 : -4, -2),
        fontsize = _ANNOTATION_FONTSIZE,
        color = color,
    )
    return nothing
end

"""
    plot_remnant_rotation(diag::RemnantDiagnostics, cfg::VisualizationConfig;
                          filename = "remnant_rotation") -> String

Two stacked panels sharing the time axis: the ordered-motion parameter
``λ_R`` with the Peebles spin ``λ_P`` on a twin axis, and the alignment
``\\cos θ`` of the spin axis with the initial orbital angular momentum,
with the coalescence time marked.
"""
function plot_remnant_rotation(
    diag::RemnantDiagnostics,
    cfg::VisualizationConfig;
    filename::AbstractString = "remnant_rotation",
)::String
    isempty(diag.time) && error("No remnant diagnostics to plot")
    ax_t = _remnant_time(diag, cfg)
    t = ax_t.t
    ttk = _time_ticks(first(t), last(t))
    c_r = _SEMANTIC_COLORS[:rotation]
    c_p = _SEMANTIC_COLORS[:spin]

    fig = Figure(; size = _fig_multipanel(cfg, 2, 1))
    ax1 = Axis(fig[1, 1]; ylabel = L"\lambda_R", xticks = ttk, xticklabelsvisible = false)
    valid = isfinite.(diag.lambda_r)
    lines!(ax1, t[valid], diag.lambda_r[valid]; color = c_r, linewidth = 2.2, label = L"\lambda_R")
    ylims!(ax1, 0.0, 1.05)
    ax1b = Axis(
        fig[1, 1];
        yaxisposition = :right,
        ylabel = L"\lambda_P",
        yticklabelcolor = c_p,
        ylabelcolor = c_p,
        xticksvisible = false,
        xticklabelsvisible = false,
        xgridvisible = false,
        ygridvisible = false,
    )
    hidespines!(ax1b)
    validp = isfinite.(diag.lambda_peebles)
    lines!(
        ax1b,
        t[validp],
        diag.lambda_peebles[validp];
        color = c_p,
        linewidth = 2.0,
        linestyle = :dash,
        label = L"\lambda_P",
    )
    linkxaxes!(ax1, ax1b)
    any(validp) && ylims!(ax1b, 0.0, 1.15 * maximum(diag.lambda_peebles[validp]))

    ax2 = Axis(
        fig[2, 1];
        xlabel = ax_t.label,
        ylabel = L"\cos\theta_{\mathrm{spin},\,L_\mathrm{orb}}",
        xticks = ttk,
    )
    valida = isfinite.(diag.spin_alignment)
    lines!(ax2, t[valida], diag.spin_alignment[valida]; color = :black, linewidth = 2.0)
    hlines!(ax2, [1.0]; color = (:grey, 0.6), linestyle = :dot, linewidth = 1.2)
    # The guide label sits at the end of the axis furthest from the
    # coalescence marker, so the two never cross.
    label_left =
        !isfinite(ax_t.t_coal) ||
        (ax_t.t_coal - first(t)) > 0.5 * max(last(t) - first(t), eps(Float64))
    text!(
        ax2,
        label_left ? first(t) : last(t),
        1.0;
        text = "Aligned with orbital L",
        align = (label_left ? :left : :right, :top),
        offset = (label_left ? 4 : -4, -2),
        fontsize = _ANNOTATION_FONTSIZE,
        color = :grey,
    )
    ylims!(ax2, -1.1, 1.15)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    tr = (first(t), last(t))
    _event_marker!(ax1, ax_t.t_coal, L"t_\mathrm{coalesce}", 1.05, tr)
    _event_marker!(ax2, ax_t.t_coal, L"t_\mathrm{coalesce}", 1.15, tr)
    if any(valid)
        _annotate!(
            ax1,
            latexstring(
                "\\lambda_R = $(_fmt_latex_sig3(diag.lambda_r[findlast(valid)]))\\;(\\mathrm{final})",
            );
            corner = _emptiest_corner(
                _with_marker(t[valid], ax_t.t_coal),
                _with_marker(diag.lambda_r[valid], maximum(diag.lambda_r[valid]));
                corners = (:tl, :tr),
            ),
            color = c_r,
        )
    end
    _top_legend!(
        fig,
        _LegendElement[
            LineElement(; color = c_r, linewidth = 2.2),
            LineElement(; color = c_p, linewidth = 2.0, linestyle = :dash),
        ],
        AbstractString[L"\lambda_R", L"\lambda_P"],
    )
    return _save_fig(cfg, filename, fig)
end

"""
    plot_rotation_profile(profile::RotationProfile, cfg::VisualizationConfig;
                          rbar = 1.0, lambda_r = NaN, filename = "remnant_rotation_profile") -> String

``v_\\mathrm{rot}/σ`` against cylindrical radius for the shells of a
[`RotationProfile`](@ref); `rbar` converts NB lengths to pc.
"""
function plot_rotation_profile(
    profile::RotationProfile,
    cfg::VisualizationConfig;
    rbar::Real = 1.0,
    lambda_r::Real = NaN,
    filename::AbstractString = "remnant_rotation_profile",
)::String
    isempty(profile.radius) && error("Empty rotation profile")
    physical = cfg.units == "physical" && rbar != 1.0
    R = physical ? profile.radius .* rbar : profile.radius
    c_r = _SEMANTIC_COLORS[:rotation]
    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = physical ? L"R \; [\mathrm{pc}]" : L"R \; [\mathrm{NB}]",
        ylabel = L"v_\mathrm{rot} / \sigma",
    )
    valid = isfinite.(profile.v_rot_over_sigma)
    scatterlines!(
        ax,
        R[valid],
        profile.v_rot_over_sigma[valid];
        color = c_r,
        linewidth = 2.0,
        markersize = 12,
    )
    hlines!(ax, [0.0]; color = (:grey, 0.6), linestyle = :dot, linewidth = 1.2)
    if isfinite(lambda_r)
        _annotate!(
            ax,
            latexstring("\\lambda_R = $(_fmt_latex_sig3(lambda_r))");
            corner = _emptiest_corner(R[valid], profile.v_rot_over_sigma[valid]),
            color = c_r,
        )
    end
    return _save_fig(cfg, filename, fig)
end

"""
    plot_remnant_structure(diag::RemnantDiagnostics, cfg::VisualizationConfig;
                           filename = "remnant_structure") -> String

Core and half-mass radii of the bound remnant against time above the ratio
``r_h / r_c``, with the coalescence time marked.
"""
function plot_remnant_structure(
    diag::RemnantDiagnostics,
    cfg::VisualizationConfig;
    filename::AbstractString = "remnant_structure",
)::String
    isempty(diag.time) && error("No remnant diagnostics to plot")
    ax_t = _remnant_time(diag, cfg)
    t = ax_t.t
    ttk = _time_ticks(first(t), last(t))
    scale = ax_t.physical ? diag.rbar : ones(length(t))
    r_c = diag.r_core .* scale
    r_h = diag.r_half .* scale
    c_core = _SEMANTIC_COLORS[:core_radius]
    c_half = _SEMANTIC_COLORS[:half_mass_radius]

    fig = Figure(; size = _fig_multipanel(cfg, 2, 1))
    ax1 = Axis(
        fig[1, 1];
        ylabel = ax_t.physical ? L"r \; [\mathrm{pc}]" : L"r \; [\mathrm{NB}]",
        xticks = ttk,
        xticklabelsvisible = false,
    )
    vc = isfinite.(r_c)
    vh = isfinite.(r_h)
    lines!(ax1, t[vc], r_c[vc]; color = c_core, linewidth = 2.2, label = L"r_c")
    lines!(ax1, t[vh], r_h[vh]; color = c_half, linewidth = 2.2, linestyle = :dash, label = L"r_h")
    r_top = 1.15 * maximum(vcat(r_c[vc], r_h[vh]); init = 1.0)
    ylims!(ax1, 0.0, r_top)
    ax2 = Axis(fig[2, 1]; xlabel = ax_t.label, ylabel = L"r_h / r_c", xticks = ttk)
    ratio = r_h ./ r_c
    vr = isfinite.(ratio)
    lines!(ax2, t[vr], ratio[vr]; color = :black, linewidth = 2.0)
    q_top = 1.15 * maximum(ratio[vr]; init = 1.0)
    ylims!(ax2, 0.0, q_top)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    tr = (first(t), last(t))
    _event_marker!(ax1, ax_t.t_coal, L"t_\mathrm{coalesce}", r_top, tr)
    _event_marker!(ax2, ax_t.t_coal, L"t_\mathrm{coalesce}", q_top, tr)
    _top_legend!(fig, ax1)
    return _save_fig(cfg, filename, fig)
end

"""
    plot_mass_segregation_evolution(diag::RemnantDiagnostics, cfg::VisualizationConfig;
                                    lambda_threshold = 2.0,
                                    filename = "remnant_mass_segregation") -> String

``Λ_\\mathrm{MSR}`` with its error band against time above the half-mass
radius ratio of the massive subset, with the threshold, coalescence and
segregation times marked.
"""
function plot_mass_segregation_evolution(
    diag::RemnantDiagnostics,
    cfg::VisualizationConfig;
    lambda_threshold::Real = _MSR_THRESHOLD,
    filename::AbstractString = "remnant_mass_segregation",
)::String
    isempty(diag.time) && error("No remnant diagnostics to plot")
    ax_t = _remnant_time(diag, cfg)
    t = ax_t.t
    ttk = _time_ticks(first(t), last(t))
    c_msr = _SEMANTIC_COLORS[:segregation]

    fig = Figure(; size = _fig_multipanel(cfg, 2, 1))
    ax1 =
        Axis(fig[1, 1]; ylabel = L"\Lambda_\mathrm{MSR}", xticks = ttk, xticklabelsvisible = false)
    v = isfinite.(diag.lambda_msr) .& isfinite.(diag.lambda_msr_err)
    l_top = 1.3 * lambda_threshold
    if any(v)
        band!(
            ax1,
            t[v],
            diag.lambda_msr[v] .- diag.lambda_msr_err[v],
            diag.lambda_msr[v] .+ diag.lambda_msr_err[v];
            color = (c_msr, 0.25),
        )
        lines!(ax1, t[v], diag.lambda_msr[v]; color = c_msr, linewidth = 2.2)
        hi = maximum(diag.lambda_msr[v] .+ diag.lambda_msr_err[v])
        l_top = 1.15 * max(hi, lambda_threshold)
    else
        _no_data_note!(ax1, "Too few members for the segregation estimate")
    end
    ylims!(ax1, 0.0, l_top)
    hlines!(ax1, [lambda_threshold]; color = (:grey, 0.7), linestyle = :dot, linewidth = 1.4)
    text!(
        ax1,
        t[end],
        lambda_threshold;
        text = latexstring(
            "\\Lambda_\\mathrm{MSR} = $(_fmt_latex_sig3(lambda_threshold))\\;\\mathrm{(threshold)}",
        ),
        align = (:right, :bottom),
        offset = (0, 3),
        fontsize = _ANNOTATION_FONTSIZE,
        color = :grey,
    )
    hlines!(ax1, [1.0]; color = (:grey, 0.5), linestyle = :dash, linewidth = 1.0)

    ax2 =
        Axis(fig[2, 1]; xlabel = ax_t.label, ylabel = L"r_{h,\mathrm{massive}} / r_h", xticks = ttk)
    vs = isfinite.(diag.segregation_ratio)
    any(vs) && lines!(ax2, t[vs], diag.segregation_ratio[vs]; color = c_msr, linewidth = 2.0)
    hlines!(ax2, [1.0]; color = (:grey, 0.5), linestyle = :dash, linewidth = 1.0)
    s_top = any(vs) ? 1.15 * max(maximum(diag.segregation_ratio[vs]), 1.0) : 1.2
    ylims!(ax2, 0.0, s_top)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    tr = (first(t), last(t))
    for (ax, top) in ((ax1, l_top), (ax2, s_top))
        _event_marker!(ax, ax_t.t_coal, L"t_\mathrm{coalesce}", top, tr)
        _event_marker!(ax, ax_t.t_seg, L"t_\mathrm{seg}", 0.85 * top, tr; color = c_msr)
    end
    return _save_fig(cfg, filename, fig)
end

"""
    remnant_figures(diag::RemnantDiagnostics, cfg::VisualizationConfig) -> Vector{String}

The four remnant figures: rotation parameters, rotation profile of the
last analysed snapshot, structure, and mass segregation.
"""
function remnant_figures(diag::RemnantDiagnostics, cfg::VisualizationConfig)
    paths = String[
        plot_remnant_rotation(diag, cfg),
        plot_remnant_structure(diag, cfg),
        plot_mass_segregation_evolution(diag, cfg),
    ]
    if !isempty(diag.profile.radius)
        k = findlast(isfinite, diag.lambda_r)
        push!(
            paths,
            plot_rotation_profile(
                diag.profile,
                cfg;
                rbar = k === nothing ? 1.0 : diag.rbar[k],
                lambda_r = k === nothing ? NaN : diag.lambda_r[k],
            ),
        )
    end
    return paths
end
