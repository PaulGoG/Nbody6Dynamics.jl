# =============================================================================
# Merger IC diagnostic plots
# =============================================================================
# Uses the same publication theme, figure sizing, axis helpers, and output
# conventions as the main simulation plotting pipeline (ext/plotting/).

"""Fraction of a frame's half-width kept free beyond the mass quantile."""
const _IC_PAD_FRAC = 0.05
"""Lowest height-over-width of a projection box (a flat configuration is kept readable)."""
const _IC_MIN_ASPECT = 0.4
"""Largest number of stars drawn as velocity arrows."""
const _IC_ARROWS_MAX = 400
"""Length of the longest velocity arrow as a fraction of the frame width."""
const _IC_ARROW_FRAC = 0.12

"""
    plot_merger_ic(result::MergerICResult, vis::VisualizationConfig)

Generate all diagnostic plots for a merger IC:
- 3× spatial projections (XY, XZ, YZ) with mass colourmap
- 3-panel overview
- Velocity field quiver plot coloured by cluster membership
- Kroupa IMF histogram with reference power-law slopes
- Radial mass-density profile per cluster (ρ in M☉ pc⁻³)

All plots use the publication theme and are saved to `vis.output_dir`.
"""
@publication function Nbody6Dynamics.plot_merger_ic(
    result::MergerICResult,
    vis::VisualizationConfig,
)
    mkpath(vis.output_dir)

    N = result.N_total
    n_cl = length(result.cluster_ranges)

    # Physical positions and masses
    x_pc = result.pos_physical[1, :]
    y_pc = result.pos_physical[2, :]
    z_pc = result.pos_physical[3, :]
    mass_solar = result.mass_physical

    # Mass colour scale (log of physical mass in solar masses)
    log_m, cmin, cmax = _log_color_range(mass_solar)

    # Marker size — same formula as snapshot plots
    # Particles are not data points: the budget's size for a sparse sample
    # is capped at the data-marker size, so a few hundred stars stay dots.
    ms = min(_marker_size(vis, N), _STYLE.marker)
    density = _density_render(vis, N)
    # Half of the mass fraction outside the snapshot frames stays outside
    # these: a sample of a few hundred stars shows every star but a sampler
    # tail, a large one keeps its halo scale.
    frac = 1 - (1 - vis.style.density_mass_frac) / 2

    # ── Frames ───────────────────────────────────────────────────
    # The frame of a projection follows the mass, not the outermost star: a
    # half-width per axis from the mass quantile `density_mass_frac` of the
    # coordinate's magnitude (the quantile the snapshot figures take on the
    # square), so a flat configuration gets a wide box instead of empty rows
    # and an outlier of the sampler does not set the scale. A box never gets
    # flatter than `_IC_MIN_ASPECT`. The overview shares one frame across its
    # three panels (common scale); each projection figure takes its own.
    half = Dict(
        k => _mass_extent(c, c, mass_solar, frac) for
        (k, c) in ((:x, x_pc), (:y, y_pc), (:z, z_pc))
    )
    ic_projections = [
        (:xy, x_pc, y_pc, :x, :y, L"x \; [\mathrm{pc}]", L"y \; [\mathrm{pc}]"),
        (:xz, x_pc, z_pc, :x, :z, L"x \; [\mathrm{pc}]", L"z \; [\mathrm{pc}]"),
        (:yz, y_pc, z_pc, :y, :z, L"y \; [\mathrm{pc}]", L"z \; [\mathrm{pc}]"),
    ]
    frame_half(hx, hy) =
        (hx * (1 + _IC_PAD_FRAC), max(hy, _IC_MIN_ASPECT * hx) * (1 + _IC_PAD_FRAC))

    # One panel: the surface-density map above `density_min_n` particles (the
    # colour range is the caller's, shared by the panels of a montage), the
    # mass-coloured scatter otherwise. Returns the scatter plot or `nothing`.
    function draw_projection!(ax, px, py, hx, hy, lo, hi, marker)
        if density
            Σ = _surface_density(px, py, mass_solar, max(hx, hy), vis.style.density_bins)
            _density_heatmap!(ax, Σ, max(hx, hy), lo, hi)
            return nothing
        end
        return scatter!(
            ax,
            px,
            py;
            color = log_m,
            colormap = :viridis,
            colorrange = (cmin, cmax),
            markersize = marker,
            strokewidth = 0,
            rasterize = true,
        )
    end
    # Colour range of the density maps of a set of projections: `density_decades`
    # below the highest peak, so the panels of one figure share a scale.
    function density_range(projs, hx, hy)
        density || return (0.0, 0.0)
        peak = maximum(
            _log_peak(_surface_density(px, py, mass_solar, max(hx, hy), vis.style.density_bins)) for (_, px, py) in projs
        )
        return (peak - vis.style.density_decades, peak)
    end
    projection_axis(position, xlab, ylab, hx, hy, band; kwargs...) = Axis(
        position;
        xlabel = xlab,
        ylabel = ylab,
        aspect = DataAspect(),
        limits = (-hx, hx, -hy, hy + band),
        xticks = _nice_ticks(-hx, hx; target_n = 5),
        yticks = _nice_ticks(-hy, hy; target_n = 5),
        xgridvisible = false,
        ygridvisible = false,
        _density_axis_attributes(density)...,
        kwargs...,
    )
    annotate_projection!(ax, text) =
        density ? _annotate!(ax, text; _DENSITY_ANNOTATION...) : _annotate!(ax, text)
    mass_colorbar!(position, sc) = Colorbar(
        position,
        sc;
        label = L"\log_{10}(m \; [\mathrm{M}_\odot])",
        ticks = _nice_colorbar_ticks(cmin, cmax),
    )

    # ── Single-projection plots ──────────────────────────────────
    # A data-free band above the data holds the projection label here too.
    box_w = _figsize_px(vis)[1] - _AXIS_PROTRUSION
    for (proj, px, py, kx, ky, xlab, ylab) in ic_projections
        hx, hy = frame_half(half[kx], half[ky])
        band = _MONTAGE_BAND_FRAC * 2hy
        lo, hi = density_range(((proj, px, py),), hx, hy)
        fig = Figure(; size = _fig_with_colorbar(vis))
        ax = projection_axis(fig[1, 1], xlab, ylab, hx, hy, band)
        sc = draw_projection!(ax, px, py, hx, hy, lo, hi, ms)
        annotate_projection!(ax, uppercase(string(proj)))
        density ? _density_colorbar!(fig[1, 2], lo, hi, true) : mass_colorbar!(fig[1, 2], sc)
        colgap!(fig.layout, _COLORBAR_COLGAP)
        # The box has the frame's shape: fix it and let the canvas follow.
        _fit_canvas_to_boxes!(fig, 1, 1, box_w, (2hy + band) / (2hx))
        _save_fig(vis, "merger_ic_$(proj)", fig)
    end

    # ── 3-panel overview ─────────────────────────────────────────
    # One row of three boxes of the same shape on one frame (the largest
    # horizontal and vertical extents of the three projections) with a
    # data-free band above the data for the projection label; the inner y
    # tick labels are hidden, so the compact gap applies, the markers scale
    # with the box, and one colourbar serves the three panels.
    hx3, hy3 = frame_half(max(half[:x], half[:y]), max(half[:y], half[:z]))
    band3 = _MONTAGE_BAND_FRAC * 2hy3
    panel_aspect = (2hy3 + band3) / (2hx3)
    gap3 = _multipanel_gap(3; inner_ticks = false)
    fig3 = Figure(;
        size = _fig_multipanel(
            vis,
            1,
            3;
            inner_ticks = false,
            panel_aspect = panel_aspect,
            extra_width = _COLORBAR_WIDTH,
            extra_height = 2 * _ANNOTATION_FONTSIZE,   # the orbit annotation row
        ),
    )
    box3 = _multipanel_box_width(vis, 3; inner_ticks = false, extra_width = _COLORBAR_WIDTH)
    ms3 = max(vis.style.marker_min, ms * box3 / _figsize_px(vis)[1])
    lo3, hi3 = density_range(((proj, px, py) for (proj, px, py, _...) in ic_projections), hx3, hy3)
    # At the single-column preset a y label per panel would leave stamp-sized
    # boxes: the projection labels name both axes, so the first panel carries
    # one shared label and the others none. The double preset labels each.
    own_ylabels = vis.column != "single"

    orbit_annot = if result.orbit_mode == "kepler"
        d_apo = Nbody6Dynamics.orbit_apocentre(result.orbit_spec, result.cluster_specs[1].rbar)
        ecc = result.orbit_spec.eccentricity
        latexstring(
            "d_\\mathrm{apo} = $(round(d_apo; digits=1))\\;\\mathrm{pc},\\; e = $(round(ecc; digits=2))",
        )
    else
        nothing
    end
    # The orbit parameters once, in an annotation row above the panels (row 1,
    # so that the panel row keeps its index for the box sizing below)
    orbit_annot === nothing || Label(
        fig3[1, 1:3],
        orbit_annot;
        tellwidth = false,
        fontsize = _ANNOTATION_FONTSIZE,
        halign = :right,
    )
    sc3 = nothing
    for (col, (proj, px, py, _, _, xlab, ylab)) in enumerate(ic_projections)
        ax = projection_axis(
            fig3[2, col],
            xlab,
            own_ylabels ? ylab : (col == 1 ? L"y,\; z \; [\mathrm{pc}]" : ""),
            hx3,
            hy3,
            band3;
            xticks = _nice_ticks(-hx3, hx3; target_n = 3),
            yticks = _nice_ticks(-hy3, hy3; target_n = 3),
            yticklabelsvisible = col == 1,
        )
        sc = draw_projection!(ax, px, py, hx3, hy3, lo3, hi3, ms3)
        sc === nothing || (sc3 = sc)
        annotate_projection!(ax, uppercase(string(proj)))
    end
    if density
        _density_colorbar!(fig3[2, 4], lo3, hi3, true)
    else
        # A colorbar as tall as one small box holds few tick labels
        Colorbar(
            fig3[2, 4],
            sc3;
            label = L"\log_{10}(m \; [\mathrm{M}_\odot])",
            ticks = _nice_colorbar_ticks(cmin, cmax; target_n = vis.column == "single" ? 3 : 6),
        )
    end
    colgap!(fig3.layout, gap3)
    foreach(c -> colsize!(fig3.layout, c, Fixed(box3)), 1:3)
    rowsize!(fig3.layout, 2, Fixed(box3 * panel_aspect))
    resize_to_layout!(fig3)
    _save_fig(vis, "merger_ic_overview", fig3)

    # ── Velocity field ───────────────────────────────────────────
    # The xy frame of the overview; a subsample of the stars as arrows whose
    # longest one spans `_IC_ARROW_FRAC` of the frame, over the stars as faint
    # dots. The legend sits above the axes, so the single-panel canvas holds.
    hxv, hyv = frame_half(half[:x], half[:y])
    bandv = _MONTAGE_BAND_FRAC * 2hyv
    fig_v = Figure(; size = _figsize_px(vis))
    ax_v = Axis(
        fig_v[1, 1];
        xlabel = L"x \; [\mathrm{pc}]",
        ylabel = L"y \; [\mathrm{pc}]",
        aspect = DataAspect(),
        limits = (-hxv, hxv, -hyv, hyv + bandv),
        xticks = _nice_ticks(-hxv, hxv; target_n = 5),
        yticks = _nice_ticks(-hyv, hyv; target_n = 5),
        xgridvisible = false,
        ygridvisible = false,
    )
    _annotate!(ax_v, "XY")

    vx_kms = result.vel_physical[1, :]
    vy_kms = result.vel_physical[2, :]

    # Faint background particles
    scatter!(
        ax_v,
        x_pc,
        y_pc;
        color = (:grey60, 0.5),
        markersize = max(2.0, 0.25 * ms),
        strokewidth = 0,
        rasterize = true,
    )

    stride = max(1, N ÷ _IC_ARROWS_MAX)
    sample = 1:stride:N
    v_max = maximum(sqrt.(vx_kms[sample] .^ 2 .+ vy_kms[sample] .^ 2))
    arrow_scale = _IC_ARROW_FRAC * 2hxv / max(v_max, 1e-10)

    # Cluster colours cycle the Okabe–Ito palette; when n_cl is large and a
    # per-cluster legend would overflow the figure, colour by density model
    # instead (first two palette colours: king → blue, plummer → orange).
    model_color(model) = model == "king" ? _OKABE_ITO[1] : _OKABE_ITO[2]
    per_cluster_legend = n_cl ≤ 8

    legend_elems = []
    legend_labels = String[]
    seen_models = Set{String}()

    for (ci, rng) in enumerate(result.cluster_ranges)
        model = profile_name(result.cluster_specs[ci].profile)
        col = per_cluster_legend ? _OKABE_ITO[mod1(ci, length(_OKABE_ITO))] : model_color(model)
        in_cluster = falses(N)
        in_cluster[rng] .= true
        sub = filter(i -> in_cluster[i], collect(sample))
        arrows!(
            ax_v,
            x_pc[sub],
            y_pc[sub],
            vx_kms[sub] .* arrow_scale,
            vy_kms[sub] .* arrow_scale;
            color = col,
            linewidth = _STYLE.envelope,
            arrowsize = 7,
        )
        if per_cluster_legend
            push!(legend_elems, [PolyElement(color = col)])
            push!(legend_labels, "Cluster $ci")
        elseif !(model in seen_models)
            push!(seen_models, model)
            push!(legend_elems, [PolyElement(color = col)])
            n_same = count(s -> profile_name(s.profile) == model, result.cluster_specs)
            push!(legend_labels, "$(titlecase(model)) ($n_same)")
        end
    end
    length(legend_labels) ≥ 2 && _top_legend!(fig_v, legend_elems, legend_labels)
    _fit_canvas_to_boxes!(fig_v, 1, 1, box_w, (2hyv + bandv) / (2hxv))

    _save_fig(vis, "merger_ic_velocity", fig_v)

    # ── IMF histogram ────────────────────────────────────────────
    # Histogram data first: the log-log axis takes explicit decade-anchored
    # ticks from the actual data extents.
    m_min = max(minimum(mass_solar), 1e-3)
    m_max = maximum(mass_solar)
    # A single mass (equal-mass IMF without a seeded IMBH) gets one bar on a
    # factor-of-two window instead of zero-width bins, and no Kroupa slopes.
    single_mass = m_max ≤ m_min * (1 + 1e-9)
    if single_mass
        m_min /= 2
        m_max *= 2
    end
    m_edges = 10.0 .^ range(log10(m_min), log10(m_max); length = 40)
    m_centres = sqrt.(m_edges[1:(end - 1)] .* m_edges[2:end])
    dm_log = diff(log10.(m_edges))
    counts = [count(m -> m_edges[j] ≤ m < m_edges[j + 1], mass_solar) for j in eachindex(dm_log)]
    dn = counts ./ dm_log
    mask = dn .> 0

    # Compute a sensible y-floor: one-third-decade below the smallest populated
    # bin, but never more than ~4 decades below the tallest bar.
    dn_min = any(mask) ? minimum(dn[mask]) : 1.0
    dn_max = any(mask) ? maximum(dn[mask]) : 10.0
    y_floor = max(dn_min / 3, dn_max / 1e4)

    fig_h = Figure(; size = _figsize_px(vis))
    ax_h = Axis(
        fig_h[1, 1];
        xlabel = L"m \; [\mathrm{M}_\odot]",
        ylabel = L"\mathrm{d}N / \mathrm{d}\log m",
        xscale = log10,
        yscale = log10,
        xticks = _log_ticks(m_min, m_max),
        yticks = _log_ticks(y_floor, dn_max),
    )

    # Render the histogram as a step-filled band: empty bins drop to `y_floor`.
    # `band!` fills between two step-expanded y-curves and is robust on log-log
    # axes, unlike per-bin `poly!` or barplot-with-log-x.
    dn_plot = [m ? d : y_floor for (m, d) in zip(mask, dn)]
    step_x = Float64[]
    step_top = Float64[]
    for j in eachindex(dn_plot)
        push!(step_x, m_edges[j])
        push!(step_top, dn_plot[j])
        push!(step_x, m_edges[j + 1])
        push!(step_top, dn_plot[j])
    end
    step_bot = fill(y_floor, length(step_x))
    hist_color = _OKABE_ITO[1]
    band!(ax_h, step_x, step_bot, step_top; color = (hist_color, 0.55))
    # Darker same-hue edge on the band fill
    lines!(ax_h, step_x, step_top; color = _band_edge(hist_color), linewidth = _STYLE.band_edge)

    # Reference slopes (the canonical Kroupa segment exponents and break
    # masses from imf.jl), normalised to the most populated bin (more robust
    # than picking m≈0.3 when the sample doesn't span the full Kroupa range).
    if any(mask) && !single_mass
        i_ref = argmax(dn)
        norm_val = dn[i_ref]
        m_ref = m_centres[i_ref]
        x_ref = 10.0 .^ range(log10(m_min / 1.5), log10(m_max * 1.5); length = 200)
        n_slopes = 0
        # Segments are drawn in order of increasing mass, each continuing the
        # previous one at the break mass, so the reference curve is unbroken
        # where the Kroupa exponent changes.
        anchor_m, anchor_y = m_ref, norm_val
        for (α, lbl, lo, hi, c) in [
            (
                _KROUPA_ALPHAS[2],
                L"\alpha = 1.3",
                _KROUPA_BREAKS[2],
                _KROUPA_BREAKS[3],
                _OKABE_ITO[2],
            ),
            (_KROUPA_ALPHAS[3], L"\alpha = 2.3", _KROUPA_BREAKS[3], 150.0, _OKABE_ITO[6]),
        ]
            seg = filter(x -> lo ≤ x ≤ hi, x_ref)
            isempty(seg) && continue
            lines!(
                ax_h,
                seg,
                anchor_y .* (seg ./ anchor_m) .^ (1.0 - α);
                color = c,
                linewidth = _STYLE.fit,
                linestyle = :dash,
                label = lbl,
            )
            anchor_y *= (hi / anchor_m)^(1.0 - α)
            anchor_m = hi
            n_slopes += 1
        end
        ylims!(ax_h, y_floor, dn_max * 2)
        n_slopes ≥ 2 && _top_legend!(fig_h, ax_h)
    end
    single_mass && ylims!(ax_h, y_floor, dn_max * 2)

    _save_fig(vis, "merger_ic_imf", fig_h)

    # ── Radial mass-density profile per cluster ──────────────────
    # True mass density: shell mass / shell volume (positions in pc, masses
    # in M☉ → ρ in M☉ pc⁻³).
    fig_r = Figure(; size = _figsize_px(vis))
    ax_r = Axis(
        fig_r[1, 1];
        xlabel = L"r \; [\mathrm{pc}]",
        ylabel = L"\rho(r) \; [\mathrm{M}_\odot\,\mathrm{pc}^{-3}]",
        xscale = log10,
        yscale = log10,
    )

    # When many clusters are present, label one line per model rather than per
    # cluster; keep all individual curves but make them semi-transparent so the
    # family (King vs Plummer) is visually dominant.
    per_cluster_density = n_cl ≤ 8
    labelled_models = Set{String}()
    line_alpha = per_cluster_density ? 1.0 : 0.55
    n_labels = 0

    # Plotted data extents for decade-anchored log-log ticks
    r_ext = (Inf, -Inf)
    ρ_ext = (Inf, -Inf)

    for (ci, rng) in enumerate(result.cluster_ranges)
        model = profile_name(result.cluster_specs[ci].profile)
        col = per_cluster_density ? _OKABE_ITO[mod1(ci, length(_OKABE_ITO))] : model_color(model)
        m_cl = mass_solar[rng]
        px_cl = x_pc[rng]
        py_cl = y_pc[rng]
        pz_cl = z_pc[rng]
        M_cl = sum(m_cl)
        cx = sum(m_cl .* px_cl) / M_cl
        cy = sum(m_cl .* py_cl) / M_cl
        cz = sum(m_cl .* pz_cl) / M_cl
        r = sqrt.((px_cl .- cx) .^ 2 .+ (py_cl .- cy) .^ 2 .+ (pz_cl .- cz) .^ 2)
        # Shells follow the sample: √N bins between 8 and 30 keep the profile
        # readable for a few hundred stars as for 10⁵
        n_bins = clamp(round(Int, sqrt(length(r))), 8, 30)
        r_lo = max(minimum(r), 0.01)
        r_hi = maximum(r)
        r_hi > r_lo || continue
        r_edges = 10.0 .^ range(log10(r_lo), log10(r_hi); length = n_bins + 1)
        r_mid = sqrt.(r_edges[1:(end - 1)] .* r_edges[2:end])
        shell_vol = (4π / 3) .* (r_edges[2:end] .^ 3 .- r_edges[1:(end - 1)] .^ 3)
        shell_mass = zeros(n_bins)
        for (ri, mi) in zip(r, m_cl)
            j = min(searchsortedlast(r_edges, ri), n_bins)
            j ≥ 1 && (shell_mass[j] += mi)
        end
        ρ = shell_mass ./ shell_vol
        m_r = ρ .> 0

        label = if per_cluster_density
            "Cluster $ci"
        elseif !(model in labelled_models)
            push!(labelled_models, model)
            n_same = count(s -> profile_name(s.profile) == model, result.cluster_specs)
            "$(titlecase(model)) ($n_same)"
        else
            nothing  # no legend entry for subsequent lines of same model
        end

        if any(m_r)
            r_ext = (min(r_ext[1], minimum(r_mid[m_r])), max(r_ext[2], maximum(r_mid[m_r])))
            ρ_ext = (min(ρ_ext[1], minimum(ρ[m_r])), max(ρ_ext[2], maximum(ρ[m_r])))
        end
        lines!(
            ax_r,
            r_mid[m_r],
            ρ[m_r];
            color = (col, line_alpha),
            linewidth = _STYLE.data,
            label = isnothing(label) ? nothing : label,
        )
        isnothing(label) || (n_labels += 1)
    end
    if isfinite(r_ext[1]) && r_ext[1] < r_ext[2]
        ax_r.xticks = _log_ticks(r_ext...)
    end
    if isfinite(ρ_ext[1]) && ρ_ext[1] < ρ_ext[2]
        ax_r.yticks = _log_ticks(ρ_ext...)
    end
    n_labels ≥ 2 && _top_legend!(fig_r, ax_r)

    _save_fig(vis, "merger_ic_density", fig_r)

    @info "All merger IC plots saved to: $(vis.output_dir)"
    return nothing
end
