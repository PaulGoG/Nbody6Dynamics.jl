# =============================================================================
# Merger IC diagnostic plots
# =============================================================================
# Uses the same publication theme, figure sizing, axis helpers, and output
# conventions as the main simulation plotting pipeline (src/plotting/).

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
function plot_merger_ic(result::MergerICResult, vis::VisualizationConfig)
    set_publication_theme!()
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
    ms = _marker_size(vis, N)

    # ── Single-projection plots ──────────────────────────────────
    ic_projections = [
        (:xy, x_pc, y_pc, L"x \; [\mathrm{pc}]", L"y \; [\mathrm{pc}]"),
        (:xz, x_pc, z_pc, L"x \; [\mathrm{pc}]", L"z \; [\mathrm{pc}]"),
        (:yz, y_pc, z_pc, L"y \; [\mathrm{pc}]", L"z \; [\mathrm{pc}]"),
    ]

    for (proj, px, py, xlab, ylab) in ic_projections
        label_str = uppercase(string(proj))
        fig = Figure(; size = _fig_with_colorbar(vis))

        xlo, xhi, ylo, yhi = _square_limits(px, py)
        xtk = _nice_ticks(xlo, xhi; target_n = 5)
        ytk = _nice_ticks(ylo, yhi; target_n = 5)

        ax = Axis(
            fig[1, 1];
            xlabel = xlab,
            ylabel = ylab,
            aspect = DataAspect(),
            limits = (xlo, xhi, ylo, yhi),
            xticks = xtk,
            yticks = ytk,
            xgridvisible = false,
            ygridvisible = false,
        )
        _annotate!(ax, label_str)
        sc = scatter!(
            ax,
            px,
            py;
            color = log_m,
            colormap = :viridis,
            colorrange = (cmin, cmax),
            markersize = ms,
            strokewidth = 0,
            rasterize = true,
        )
        Colorbar(
            fig[1, 2],
            sc;
            label = L"\log_{10}(m \; [\mathrm{M}_\odot])",
            ticks = _nice_colorbar_ticks(cmin, cmax),
        )
        colgap!(fig.layout, _COLORBAR_COLGAP)

        _save_fig(vis, "merger_ic_$(proj)", fig)
    end

    # ── 3-panel overview ─────────────────────────────────────────
    pw, ph = _fig_multipanel(vis, 1, 3)
    fig3 = Figure(; size = (pw + _COLORBAR_WIDTH, ph))   # extra width for shared colorbar

    # Global limits across all 3 projections
    all_coords = vcat(x_pc, y_pc, z_pc)
    xlo, xhi, _, _ = _square_limits(all_coords, all_coords)
    tk = _nice_ticks(xlo, xhi; target_n = 5)

    orbit_annot = if result.orbit_mode == "kepler"
        d_apo = result.orbit_spec.apocentre
        ecc = result.orbit_spec.eccentricity
        latexstring(
            "d_\\mathrm{apo} = $(round(d_apo; digits=1))\\;\\mathrm{pc},\\; e = $(round(ecc; digits=2))",
        )
    else
        nothing
    end

    for (col, (proj, px, py, xlab, ylab)) in enumerate(ic_projections)
        label_str = uppercase(string(proj))
        # Every panel carries its own ylabel — the y-quantities differ
        # (Y, Z, Z); only the tick labels are shared/suppressed.
        ax = Axis(
            fig3[1, col];
            xlabel = xlab,
            ylabel = ylab,
            aspect = DataAspect(),
            limits = (xlo, xhi, xlo, xhi),
            xticks = tk,
            yticks = tk,
            yticklabelsvisible = col == 1,
            xgridvisible = false,
            ygridvisible = false,
        )
        _annotate!(ax, label_str)
        if col == 1 && orbit_annot !== nothing
            _annotate!(ax, orbit_annot; dy = 0.07)
        end
        sc = scatter!(
            ax,
            px,
            py;
            color = log_m,
            colormap = :viridis,
            colorrange = (cmin, cmax),
            markersize = ms * 0.7,
            strokewidth = 0,
            rasterize = true,
        )
        col == 3 && Colorbar(
            fig3[1, 4],
            sc;
            label = L"\log_{10}(m \; [\mathrm{M}_\odot])",
            ticks = _nice_colorbar_ticks(cmin, cmax),
        )
    end

    colgap!(fig3.layout, _MULTIPANEL_HGAP)

    _save_fig(vis, "merger_ic_overview", fig3)

    # ── Velocity field ───────────────────────────────────────────
    # No side panel any more (legend sits above the axes) → single-panel size
    fig_v = Figure(; size = _figsize_px(vis))
    ax_v = Axis(
        fig_v[1, 1];
        xlabel = L"x \; [\mathrm{pc}]",
        ylabel = L"y \; [\mathrm{pc}]",
        aspect = DataAspect(),
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
        color = (:grey70, 0.3),
        markersize = 0.8,
        strokewidth = 0,
        rasterize = true,
    )

    stride = max(1, N ÷ 500)
    v_max = maximum(sqrt.(vx_kms[1:stride:N] .^ 2 .+ vy_kms[1:stride:N] .^ 2))
    arrow_scale = 2.0 / max(v_max, 1e-10)

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
        sub = filter(i -> i in rng, collect(1:stride:N))
        arrows!(
            ax_v,
            x_pc[sub],
            y_pc[sub],
            vx_kms[sub] .* arrow_scale,
            vy_kms[sub] .* arrow_scale;
            color = col,
            linewidth = 1.0,
            arrowsize = 4,
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

    _save_fig(vis, "merger_ic_velocity", fig_v)

    # ── IMF histogram ────────────────────────────────────────────
    # Histogram data first: the log-log axis takes explicit decade-anchored
    # ticks from the actual data extents.
    m_min = max(minimum(mass_solar), 1e-3)
    m_max = maximum(mass_solar)
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
    lines!(ax_h, step_x, step_top; color = _band_edge(hist_color), linewidth = 1.5)

    # Reference slopes (the canonical Kroupa segment exponents and break
    # masses from imf.jl), normalised to the most populated bin (more robust
    # than picking m≈0.3 when the sample doesn't span the full Kroupa range).
    if any(mask)
        i_ref = argmax(dn)
        norm_val = dn[i_ref]
        m_ref = m_centres[i_ref]
        x_ref = 10.0 .^ range(log10(m_min / 1.5), log10(m_max * 1.5); length = 200)
        n_slopes = 0
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
                norm_val .* (seg ./ m_ref) .^ (1.0 - α);
                color = c,
                linewidth = 2,
                linestyle = :dash,
                label = lbl,
            )
            n_slopes += 1
        end
        ylims!(ax_h, y_floor, dn_max * 2)
        n_slopes ≥ 2 && _top_legend!(fig_h, ax_h)
    end

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
        n_bins = 30
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
            linewidth = 2,
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
