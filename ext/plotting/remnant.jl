# =============================================================================
# Remnant figures: rotation, structure, mass segregation
# =============================================================================

"""Time axis of a `RemnantDiagnostics` in the configured units: values, label and the epochs."""
function _remnant_time(diag::RemnantDiagnostics, cfg::VisualizationConfig)
    physical = cfg.units == "physical" && all(isfinite, diag.time_myr) && any(!=(1.0), diag.rbar)
    t = physical ? diag.time_myr : diag.time
    label = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]"
    t_coal = physical ? diag.coalescence_time_myr : diag.coalescence_time
    t_seg = physical ? diag.segregation_time_myr : diag.segregation_time
    t_contact = physical ? diag.t_contact_myr : diag.t_contact
    t_core = physical ? diag.t_core_myr : diag.t_core
    return (; physical, t, label, t_coal, t_seg, t_contact, t_core)
end

"""Height of each panel of the two-panel remnant figures in layout units, for the placement of the epoch labels."""
const _REMNANT_PANEL_HEIGHT = 400

"""Markers of the merger epochs on a remnant figure: first contact, core coalescence, centre coalescence."""
_remnant_epoch_markers(ax_t) = [
    (ax_t.t_contact, L"t_\mathrm{contact}", _EPOCH_STYLE.contact...),
    (ax_t.t_core, L"t_\mathrm{core}", _EPOCH_STYLE.core...),
    (ax_t.t_coal, L"t_\mathrm{coalesce}", _EPOCH_STYLE.join...),
]

"""
    _corner_clear_of(xs, ys, marker_xs, xlims, ylims; corners = (:tl, :tr), width = 0.35)

The emptiest corner ([`_emptiest_corner`](@ref) within the axis limits) among
those of `corners` whose side holds none of the vertical markers at
`marker_xs`, whose labels hang from the top of the panel; among all of
`corners` when every side holds one.
"""
function _corner_clear_of(xs, ys, marker_xs, xlims, ylims; corners = (:tl, :tr), width = 0.35)
    x0, x1 = Float64.(xlims)
    span = max(x1 - x0, eps(Float64))
    us = [(x - x0) / span for x in marker_xs if isfinite(x)]
    clear =
        Tuple(c for c in corners if !any(u -> (c in (:tr, :br) ? u ≥ 1 - width : u ≤ width), us))
    return _emptiest_corner(
        xs,
        ys;
        corners = isempty(clear) ? corners : clear,
        xlims = (x0, x1),
        ylims = Float64.(ylims),
    )
end

"""
    plot_remnant_rotation(diag::RemnantDiagnostics, cfg::VisualizationConfig;
                          filename = "remnant_rotation") -> String

Two stacked panels sharing the time axis: the ordered-motion parameter
``λ_R`` with the Peebles spin ``λ_P`` on a twin axis, and the alignment
``\\cos θ`` of the spin axis with the orbital angular momentum of the
most-bound cores of the progenitors at the last snapshot before the core coalescence (`MergerEpochs`;
the initial orbital angular momentum for a configuration that never
coalesces). The first contact, the core coalescence and the centre
coalescence are marked and the `interacting` and `coalescing` regimes
shaded; the series start at the core coalescence, before which the bound set
holds several bodies.
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
    # Headroom above cos θ = 1 for the two lines of the guide label
    a_top = max(1.07, 1.0 + 0.16 * (1.1 - a_lo))
    ax2 = Axis(
        fig[2, 1];
        xlabel = ax_t.label,
        ylabel = L"\cos\theta_{\mathrm{spin,\,orb}}",
        xticks = ttk,
        # No tick above cos θ = 1: the headroom holds the guide label only
        yticks = filter(≤(1.0 + 1e-9), _nice_ticks(a_lo - 0.1, a_top)),
    )
    lines!(ax2, t[valida], diag.spin_alignment[valida]; color = :black, linewidth = _STYLE.data)
    hlines!(ax2, [1.0]; color = (:grey, 0.6), linestyle = :dot, linewidth = _STYLE.guide)
    # The reference axis: the orbital angular momentum of the most-bound cores at
    # the core coalescence (the initial orbit when the cores never coalesce)
    orbit_lines =
        isfinite(ax_t.t_core) ?
        [
            L"Aligned with $\mathbf{L}_\mathrm{orb}$ of the",
            L"most-bound cores before $t_\mathrm{core}$",
        ] : [L"Aligned with the initial $\mathbf{L}_\mathrm{orb}$"]
    # Box of the guide label in data coordinates, which the epoch labels keep clear of
    label_box = Tuple{Vector{Float64},Vector{Float64}}[]
    if any(valida)
        # Right-aligned at the end of the series, below the guide when the series
        # leaves that stretch free, above it otherwise
        xs_a = t[valida]
        ys_a = diag.spin_alignment[valida]
        x_end = maximum(xs_a)
        width = maximum(_label_length, orbit_lines) / _MARKER_AXIS_WIDTH * (last(t) - first(t))
        a_span = a_top - (a_lo - 0.1)
        line_height = 1.25 * _ANNOTATION_FONTSIZE / _REMNANT_PANEL_HEIGHT * a_span
        height = length(orbit_lines) * line_height
        below =
            !any(
                i -> xs_a[i] ≥ x_end - width && 1.0 - 0.03 * a_span - height ≤ ys_a[i] ≤ 1.0,
                eachindex(xs_a),
            )
        y_first = below ? 1.0 - 0.02 * a_span : 1.0 + 0.02 * a_span + height
        for (j, line) in enumerate(orbit_lines)
            text!(
                ax2,
                x_end,
                y_first - (j - 1) * line_height;
                text = line,
                align = (:right, :top),
                color = :grey,
                fontsize = _ANNOTATION_FONTSIZE,
            )
        end
        # Rows across the box, closer than the padding of the occupied intervals
        x_box = [x_end - width, x_end]
        for y in range(y_first - height, y_first; length = 6)
            push!(label_box, (x_box, [y, y]))
        end
    end
    ylims!(ax2, a_lo - 0.1, a_top)
    tr = (first(t), last(t))
    bands = _regime_bands!((ax1, ax2), t, diag.regime; label = false)
    for ax in (ax1, ax1b, ax2)
        _time_limits!(ax, t)
    end
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    markers = _remnant_epoch_markers(ax_t)
    # λ_P on the scale of the λ_R axis, for the placement of the epoch labels
    p_scale = any(validp) ? l_top / (1.15 * maximum(diag.lambda_peebles[validp])) : 1.0
    _epoch_markers!(
        ax1,
        markers,
        (0.0, l_top),
        tr;
        series = [
            (t[valid], diag.lambda_r[valid]),
            (t[validp], diag.lambda_peebles[validp] .* p_scale),
        ],
        axis_height = _REMNANT_PANEL_HEIGHT,
    )
    _epoch_markers!(
        ax2,
        markers,
        (a_lo - 0.1, a_top),
        tr;
        series = vcat(
            [(t[valida], diag.spin_alignment[valida]), ([first(t), last(t)], [1.0, 1.0])],
            label_box,
        ),
        axis_height = _REMNANT_PANEL_HEIGHT,
    )
    if any(valid)
        _annotate!(
            ax1,
            latexstring(
                "\\lambda_R = $(_fmt_latex_sig3(diag.lambda_r[findlast(valid)]))\\;(\\mathrm{final})",
            );
            corner = _corner_clear_of(
                t[valid],
                diag.lambda_r[valid],
                (m[1] for m in markers),
                tr,
                (0.0, l_top),
            ),
            color = c_r,
        )
    end
    band_elements, band_labels = _regime_legend(bands)
    _top_legend!(
        fig,
        _LegendElement[
            LineElement(; color = c_r, linewidth = _STYLE.data),
            LineElement(; color = c_p, linewidth = _STYLE.data, linestyle = :dash),
            band_elements...,
        ],
        AbstractString[L"\lambda_R", L"\lambda_P", band_labels...],
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
                corners = (:tl, :tr),
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
``r_h / r_c``, from the core coalescence on, with the first contact, the core
coalescence and the centre coalescence marked and the `interacting` and
`coalescing` regimes shaded.
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
    _regime_bands!((ax1, ax2), t, diag.regime)
    _time_limits!(ax1, t)
    _time_limits!(ax2, t)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    tr = (first(t), last(t))
    markers = _remnant_epoch_markers(ax_t)
    _epoch_markers!(
        ax1,
        markers,
        (0.0, r_top),
        tr;
        series = [(t[vc], r_c[vc]), (t[vh], r_h[vh])],
        axis_height = _REMNANT_PANEL_HEIGHT,
    )
    _epoch_markers!(
        ax2,
        markers,
        (0.0, q_top),
        tr;
        series = [(t[vr], ratio[vr])],
        axis_height = _REMNANT_PANEL_HEIGHT,
    )
    _top_legend!(fig, ax1)
    return _save_fig(cfg, filename, fig)
end

"""
    plot_mass_segregation_evolution(diag::RemnantDiagnostics, cfg::VisualizationConfig;
                                    lambda_threshold = 2.0,
                                    filename = "remnant_mass_segregation") -> String

``Λ_\\mathrm{MSR}`` with its error band against time above the half-mass
radius ratio of the massive subset, from the core coalescence on, with the
threshold, the first contact, the core and centre coalescences and the
segregation time marked and the `interacting` and `coalescing` regimes
shaded.
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
        # The guide labels keep clear of the error band, not only of its centre line,
        # and of each other: the unity label sees the threshold guide and the strip
        # below it, where the threshold label may sit, as occupied
        band_x = vcat(t[v], t[v], t[v])
        band_y = vcat(diag.lambda_msr[v], lo_band, hi_band)
        guide_x = collect(range(first(t), last(t); length = 50))
        _guide_label!(
            ax1,
            latexstring(
                "\\Lambda_\\mathrm{MSR} = $(_fmt_latex_sig3(lambda_threshold))\\;\\mathrm{(segregated\\;above)}",
            ),
            lambda_threshold;
            xs = band_x,
            ys = band_y,
            color = :grey,
        )
        _guide_label!(
            ax1,
            L"\Lambda_\mathrm{MSR} = 1",
            1.0;
            xs = vcat(band_x, guide_x, guide_x),
            ys = vcat(
                band_y,
                fill(Float64(lambda_threshold), 50),
                fill(0.5 * (1.0 + lambda_threshold), 50),
            ),
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
    bands = _regime_bands!((ax1, ax2), t, diag.regime; label = false)
    _time_limits!(ax1, t)
    _time_limits!(ax2, t)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    tr = (first(t), last(t))
    markers = vcat(_remnant_epoch_markers(ax_t), [(ax_t.t_seg, L"t_\mathrm{seg}", c_msr, :dash)])
    msr_series =
        any(v) ?
        [
            (t[v], diag.lambda_msr[v] .- diag.lambda_msr_err[v]),
            (t[v], diag.lambda_msr[v] .+ diag.lambda_msr_err[v]),
        ] : Tuple{Vector{Float64},Vector{Float64}}[]
    _epoch_markers!(
        ax1,
        markers,
        (0.0, l_top),
        tr;
        series = msr_series,
        axis_height = _REMNANT_PANEL_HEIGHT,
    )
    _epoch_markers!(
        ax2,
        markers,
        (0.0, s_top),
        tr;
        series = [(t[vs], diag.segregation_ratio[vs])],
        axis_height = _REMNANT_PANEL_HEIGHT,
    )
    # The regimes are the only regions to name: one entry per regime shaded
    isempty(bands) || _top_legend!(fig, _regime_legend(bands)...; title = "Regime:")
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
