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

"""Clearance of an event label from its own line, as a fraction of the time range."""
const _MARKER_LABEL_CLEARANCE = 0.01

"""Dashed vertical marker at `x` with `label` hanging from `y_top` (data
coordinates, the caller's upper axis limit) clear of the line — one per cent
of `t_range` to its right, or to its left when `x` lies in the last
`_MARKER_EDGE_FRACTION` of `t_range`; no-op for `NaN`."""
function _event_marker!(
    ax,
    x::Real,
    label::AbstractString,
    y_top::Real,
    t_range::Tuple{<:Real,<:Real};
    color = :black,
)
    isfinite(x) || return nothing
    vlines!(ax, [x]; color = color, linestyle = :dash, linewidth = _STYLE.guide)
    span = max(t_range[2] - t_range[1], eps(Float64))
    right_side = (x - t_range[1]) / span < 1 - _MARKER_EDGE_FRACTION
    dx = _MARKER_LABEL_CLEARANCE * span
    text!(
        ax,
        right_side ? x + dx : x - dx,
        y_top;
        text = label,
        align = (right_side ? :left : :right, :top),
        offset = (0, -2),
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
@publication function Nbody6Dynamics.plot_remnant_rotation(
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
    lines!(
        ax1,
        t[valid],
        diag.lambda_r[valid];
        color = c_r,
        linewidth = _STYLE.data,
        label = L"\lambda_R",
    )
    λ_hi = any(valid) ? max(1.0, maximum(diag.lambda_r[valid])) : 1.0
    l_top = _pad_limits(0.0, λ_hi; floor = 0)[2]
    ylims!(ax1, 0.0, l_top)
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
        linewidth = _STYLE.data,
        linestyle = :dash,
        label = L"\lambda_P",
    )
    linkxaxes!(ax1, ax1b)
    any(validp) && ylims!(ax1b, 0.0, 1.15 * maximum(diag.lambda_peebles[validp]))

    valida = isfinite.(diag.spin_alignment)
    # The remnant spin stays close to the orbital axis: the panel holds the
    # data down to cos θ = 0.8 at least, instead of the full [-1, 1] range
    # the alignment could in principle take.
    a_lo = any(valida) ? min(minimum(diag.spin_alignment[valida]), 0.8) : 0.8
    # Headroom for the guide label, which clears its line by 3 % of the data range
    a_top = max(1.05, 1.0 + 0.05 * (1.1 - a_lo))
    ax2 = Axis(
        fig[2, 1];
        xlabel = ax_t.label,
        ylabel = L"\cos\theta_{\mathrm{spin,\,orb}}",
        xticks = ttk,
        yticks = _nice_ticks(a_lo - 0.1, a_top),
    )
    lines!(ax2, t[valida], diag.spin_alignment[valida]; color = :black, linewidth = _STYLE.data)
    hlines!(ax2, [1.0]; color = (:grey, 0.6), linestyle = :dot, linewidth = _STYLE.guide)
    any(valida) && _guide_label!(
        ax2,
        "Aligned with orbital L",
        1.0;
        xs = t[valida],
        ys = diag.spin_alignment[valida],
        color = :grey,
    )
    ylims!(ax2, a_lo - 0.1, a_top)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    tr = (first(t), last(t))
    _event_marker!(ax1, ax_t.t_coal, L"t_\mathrm{coalesce}", l_top, tr)
    _event_marker!(ax2, ax_t.t_coal, L"t_\mathrm{coalesce}", a_top, tr)
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
                avoid_x = ax_t.t_coal,
            ),
            color = c_r,
        )
    end
    _top_legend!(
        fig,
        _LegendElement[
            LineElement(; color = c_r, linewidth = _STYLE.data),
            LineElement(; color = c_p, linewidth = _STYLE.data, linestyle = :dash),
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
@publication function Nbody6Dynamics.plot_rotation_profile(
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
        linewidth = _STYLE.data,
        markersize = _STYLE.marker,
        strokecolor = _band_edge(c_r),
        strokewidth = _STYLE.marker_stroke,
    )
    hlines!(ax, [0.0]; color = (:grey, 0.6), linestyle = :dot, linewidth = _STYLE.guide)
    any(valid) && _guide_label!(
        ax,
        L"v_\mathrm{rot} = 0",
        0.0;
        xs = R[valid],
        ys = profile.v_rot_over_sigma[valid],
        color = :grey,
    )
    if isfinite(lambda_r)
        _annotate!(
            ax,
            latexstring("\\lambda_R = $(_fmt_latex_sig3(lambda_r))");
            corner = _emptiest_corner(
                R[valid],
                profile.v_rot_over_sigma[valid];
                corners = (:tl, :tr, :br),
            ),
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
@publication function Nbody6Dynamics.plot_remnant_structure(
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
    lines!(ax1, t[vc], r_c[vc]; color = c_core, linewidth = _STYLE.data, label = L"r_c")
    lines!(
        ax1,
        t[vh],
        r_h[vh];
        color = c_half,
        linewidth = _STYLE.data,
        linestyle = :dash,
        label = L"r_h",
    )
    r_top = 1.15 * maximum(vcat(r_c[vc], r_h[vh]); init = 1.0)
    ylims!(ax1, 0.0, r_top)
    ax2 = Axis(fig[2, 1]; xlabel = ax_t.label, ylabel = L"r_h / r_c", xticks = ttk)
    ratio = r_h ./ r_c
    vr = isfinite.(ratio)
    lines!(ax2, t[vr], ratio[vr]; color = :black, linewidth = _STYLE.data)
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
@publication function Nbody6Dynamics.plot_mass_segregation_evolution(
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
        lo_band = diag.lambda_msr[v] .- diag.lambda_msr_err[v]
        hi_band = diag.lambda_msr[v] .+ diag.lambda_msr_err[v]
        band!(ax1, t[v], lo_band, hi_band; color = (c_msr, _STYLE.band_alpha))
        # Darker same-hue edges on the error band
        msr_edge = _band_edge(c_msr)
        lines!(ax1, t[v], lo_band; color = msr_edge, linewidth = _STYLE.band_edge)
        lines!(ax1, t[v], hi_band; color = msr_edge, linewidth = _STYLE.band_edge)
        lines!(ax1, t[v], diag.lambda_msr[v]; color = c_msr, linewidth = _STYLE.data)
        l_top = 1.15 * max(maximum(hi_band), lambda_threshold)
    else
        _no_data_note!(ax1, "Too few members for the segregation estimate")
    end
    ylims!(ax1, 0.0, l_top)
    hlines!(
        ax1,
        [lambda_threshold];
        color = (:grey, 0.7),
        linestyle = :dot,
        linewidth = _STYLE.guide,
    )
    hlines!(ax1, [1.0]; color = (:grey, 0.5), linestyle = :dash, linewidth = _STYLE.guide)
    if any(v)
        _guide_label!(
            ax1,
            latexstring(
                "\\Lambda_\\mathrm{MSR} = $(_fmt_latex_sig3(lambda_threshold))\\;\\mathrm{(segregated\\;above)}",
            ),
            lambda_threshold;
            xs = t[v],
            ys = diag.lambda_msr[v],
            color = :grey,
        )
        _guide_label!(
            ax1,
            L"\Lambda_\mathrm{MSR} = 1",
            1.0;
            xs = t[v],
            ys = diag.lambda_msr[v],
            color = :grey,
        )
    end

    ax2 =
        Axis(fig[2, 1]; xlabel = ax_t.label, ylabel = L"r_{h,\mathrm{massive}} / r_h", xticks = ttk)
    vs = isfinite.(diag.segregation_ratio)
    any(vs) &&
        lines!(ax2, t[vs], diag.segregation_ratio[vs]; color = c_msr, linewidth = _STYLE.data)
    hlines!(ax2, [1.0]; color = (:grey, 0.5), linestyle = :dash, linewidth = _STYLE.guide)
    any(vs) && _guide_label!(
        ax2,
        L"r_{h,\mathrm{massive}} = r_h",
        1.0;
        xs = t[vs],
        ys = diag.segregation_ratio[vs],
        color = :grey,
    )
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
function Nbody6Dynamics.remnant_figures(diag::RemnantDiagnostics, cfg::VisualizationConfig)
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
