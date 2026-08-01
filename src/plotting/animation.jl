# =============================================================================
# Animation of cluster evolution — GIF output
# =============================================================================

"""
    _auto_fps(nframes; target_duration = 12.0, min_fps = 1, max_fps = 30) -> Int

Compute a frame rate so the animation lasts approximately `target_duration`
seconds, clamped to `[min_fps, max_fps]`.
"""
function _auto_fps(nframes::Int;
                   target_duration::Float64 = 12.0,
                   min_fps::Int = 1, max_fps::Int = 30)::Int
    return clamp(round(Int, nframes / target_duration), min_fps, max_fps)
end

"""
    animate_cluster(snaps::Vector{Snapshot}, cfg::VisualizationConfig;
                    filename::AbstractString = "cluster_evolution",
                    projections::Vector{Symbol} = [:xy, :xz, :yz],
                    fps::Union{Int,Nothing} = nothing) -> Vector{String}

Create animated GIFs of the cluster's spatial evolution across snapshots,
one per projection.  Positions are shown in pc and the time annotation in
Myr when `cfg.units == "physical"` (header AS scaling); N-body otherwise.

# Arguments
- `snaps`: ordered vector of `Snapshot`s
- `cfg`: visualisation configuration (figsize, dpi, output_dir)
- `filename`: output filename stem (`.gif` appended automatically)
- `projections`: spatial projections to animate (any of :xy, :xz, :yz)
- `fps`: frames per second (`nothing` = use `cfg.style.anim_fps`; `0` there
  auto-calculates for ~`cfg.style.anim_target_seconds` s, clamped to 1–10 fps)

Returns a vector of output file paths.
"""
function animate_cluster(
    snaps::Vector{Snapshot}, cfg::VisualizationConfig;
    filename::AbstractString = "cluster_evolution",
    projections::Vector{Symbol} = [:xy, :xz, :yz],
    fps::Union{Int,Nothing} = nothing,
)::Vector{String}
    isempty(snaps) && error("No snapshots to animate")

    nframes = length(snaps)
    # Cluster snapshots are discrete, each worth studying individually:
    # fewer snapshots → slower pace, many → speed up.
    fps = something(fps, cfg.style.anim_fps)
    fps > 0 || (fps = _auto_fps(nframes;
        target_duration = cfg.style.anim_target_seconds, min_fps = 1, max_fps = 10))

    physical = cfg.units == "physical" &&
               all(_has_physical_scaling(s.header) for s in snaps)
    unit_str = physical ? "pc" : "NB"
    scales   = [physical ? rbar(snap.header) : 1.0 for snap in snaps]

    ms = _marker_size(cfg, nparticles(snaps[1]))

    # Mass colour scale — global across all frames for consistency
    all_m = reduce(vcat, [Float64.(snap.mass) for snap in snaps])
    log_m_global = log10.(max.(all_m, 1e-30))
    cmin, cmax = extrema(log_m_global)
    if cmin ≈ cmax
        cmin -= 0.5
        cmax += 0.5
    end

    outpaths = String[]

    for projection in projections
        ix, iy, xsym, ysym = _proj_indices(projection)

        # Check if adaptive zoom is needed
        all_indices = collect(1:nframes)
        use_adaptive = _needs_adaptive_zoom(snaps, all_indices, ix, iy, cfg.style.zoom_frac)

        # Pre-compute per-frame limits for adaptive mode
        frame_limits = if use_adaptive
            [_square_limits(snaps[i].pos[ix, :] .* scales[i],
                            snaps[i].pos[iy, :] .* scales[i]) for i in 1:nframes]
        else
            # Global limits — same for every frame
            all_x = reduce(vcat, [snaps[i].pos[ix, :] .* scales[i] for i in 1:nframes])
            all_y = reduce(vcat, [snaps[i].pos[iy, :] .* scales[i] for i in 1:nframes])
            lims = _square_limits(all_x, all_y)
            fill(lims, nframes)
        end

        outpath = _anim_output_path(cfg, "$(filename)_$(projection)")

        mode_str = use_adaptive ? "adaptive zoom" : "global limits"
        @info "Animating $nframes frames ($projection, $mode_str) → $outpath  ($(fps) fps, ~$(round(Int, nframes/fps)) s)"

        # Observable for the frame index
        frame_idx = Observable(1)

        time_text = @lift begin
            h = snaps[$frame_idx].header
            _time_annotation(physical ? time_myr(h) : time_nb(h), physical)
        end

        # Initial limits
        xlo0, xhi0, ylo0, yhi0 = frame_limits[1]

        fig = Figure(; size = _fig_with_colorbar(cfg))
        ax = Axis(fig[1, 1];
            xlabel = _coord_label(xsym, unit_str),
            ylabel = _coord_label(ysym, unit_str),
            aspect = DataAspect(),
            limits = (xlo0, xhi0, ylo0, yhi0),
            xticks = _nice_ticks(xlo0, xhi0; target_n = 5),
            yticks = _nice_ticks(ylo0, yhi0; target_n = 5),
            xgridvisible = false,
            ygridvisible = false,
        )
        text!(ax, 0.04, 0.96; text = time_text,
            space = :relative, align = (:left, :top), fontsize = 16)

        # Single source observable so positions and colours update atomically
        # even when the particle count changes between frames (Point2f handles
        # the varying length).
        frame_data = @lift begin
            snap = snaps[$frame_idx]
            sc = scales[$frame_idx]
            (points = Point2f.(snap.pos[ix, :] .* sc, snap.pos[iy, :] .* sc),
             colors = log10.(max.(Float64.(snap.mass), 1e-30)))
        end
        pts    = @lift($frame_data.points)
        colors = @lift($frame_data.colors)

        scatter!(ax, pts;
            color      = colors,
            colormap   = :viridis,
            colorrange = (cmin, cmax),
            markersize = ms,
            strokewidth = 0,
        )

        # Add colorbar
        Colorbar(fig[1, 2]; colormap = :viridis, colorrange = (cmin, cmax),
                 label = L"\log_{10}(m \, / \, M_\mathrm{tot})",
                 ticks = _nice_colorbar_ticks(cmin, cmax))
        colgap!(fig.layout, 10)

        _backup_existing(outpath)
        record(fig, outpath, 1:nframes; framerate = fps) do i
            frame_idx[] = i
            if use_adaptive
                xlo, xhi, ylo, yhi = frame_limits[i]
                xlims!(ax, xlo, xhi)
                ylims!(ax, ylo, yhi)
                ax.xticks = _nice_ticks(xlo, xhi; target_n = 5)
                ax.yticks = _nice_ticks(ylo, yhi; target_n = 5)
            end
        end

        @info "Animation saved: $outpath  ($nframes frames, $(fps) fps)"
        push!(outpaths, outpath)
    end

    return outpaths
end

"""
    animate_hr(sevs::Vector{StellarEvolutionSnapshot}, cfg::VisualizationConfig;
               filename::AbstractString = "hr_evolution_anim",
               fps::Union{Int,Nothing} = nothing) -> String

Animate HR diagram evolution across stellar evolution snapshots.

# Arguments
- `fps`: frames per second (`nothing` = use `cfg.style.anim_fps`; `0` there
  auto-calculates for ~`cfg.style.anim_target_seconds` s, clamped to 1–8 fps).
  HR frames are information-dense, so the auto rate favours a slower pace
  than cluster animations.

Returns the output file path.
"""
function animate_hr(
    sevs::Vector{StellarEvolutionSnapshot}, cfg::VisualizationConfig;
    filename::AbstractString = "hr_evolution_anim",
    fps::Union{Int,Nothing} = nothing,
)::String
    isempty(sevs) && error("No stellar evolution snapshots to animate")

    # Filter placeholder values per frame; use for both limits and plotted data
    valid_per_frame = [_hr_valid_records(sev.records) for sev in sevs]

    # Compute global axis limits from valid records only
    all_teff = reduce(vcat, [[r.log_teff for r in v] for v in valid_per_frame]; init = Float64[])
    all_lum  = reduce(vcat, [[r.log_luminosity for r in v] for v in valid_per_frame]; init = Float64[])
    isempty(all_teff) && error("No valid HR records in any snapshot")

    tmin, tmax = extrema(all_teff)
    lmin, lmax = extrema(all_lum)
    dt = (tmax - tmin) * 0.06
    dl = (lmax - lmin) * 0.06
    xlims = (tmin - dt, tmax + dt)
    ylims = (lmin - dl, lmax + dl)

    nframes = length(sevs)
    # HR frames are information-dense (stellar types, population structure);
    # keep pace slow so each epoch is readable.
    fps = something(fps, cfg.style.anim_fps)
    fps > 0 || (fps = _auto_fps(nframes;
        target_duration = cfg.style.anim_target_seconds, min_fps = 1, max_fps = 8))

    outpath = _anim_output_path(cfg, filename)

    @info "Animating HR diagram: $nframes frames → $outpath  ($(fps) fps, ~$(round(Int, nframes/fps)) s)"

    frame_idx = Observable(1)

    # sev.time_myr is in Myr, not NB units
    time_text = @lift begin
        t_str = @sprintf("%.3g", sevs[$frame_idx].time_myr)
        latexstring("t = $(t_str)\\;\\mathrm{Myr}")
    end

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(fig[1, 1];
        xlabel = L"\log_{10}(T_\mathrm{eff} \, / \, \mathrm{K})",
        ylabel = L"\log_{10}(L \, / \, L_\odot)",
        limits = (xlims..., ylims...),
        xreversed = true,
        xticks = _logval_ticks(xlims[1], xlims[2]),
        yticks = _logval_ticks(ylims[1], ylims[2]),
        xgridvisible = false,
        ygridvisible = false,
    )
    # Top-right in-axis corner is empty on an HR diagram (the sequence
    # enters at top-left).
    text!(ax, 0.96, 0.96; text = time_text,
        space = :relative, align = (:right, :top), fontsize = 16)

    # Single source observable so positions, colours, and markers update
    # atomically even when the star count changes between epochs (Point2f
    # handles the varying length).
    frame_data = @lift begin
        v = valid_per_frame[$frame_idx]
        (points = Point2f.([r.log_teff for r in v],
                           [r.log_luminosity for r in v]),
         colors  = [_hr_color(r.stellar_type) for r in v],
         markers = [_hr_marker(r.stellar_type) for r in v])
    end
    pts      = @lift($frame_data.points)
    col_data = @lift($frame_data.colors)
    mk_data  = @lift($frame_data.markers)

    scatter!(ax, pts;
        color = col_data, marker = mk_data, markersize = 14, strokewidth = 0)

    _backup_existing(outpath)
    record(fig, outpath, 1:nframes; framerate = fps) do i
        frame_idx[] = i
    end

    @info "HR animation saved: $outpath  ($nframes frames, $(fps) fps)"
    return outpath
end

"""
    animate_lagrangian(lagr::LagrangianData, cfg::VisualizationConfig;
                       filename::AbstractString = "lagrangian_anim",
                       fps::Union{Int,Nothing} = nothing,
                       selected_fractions::Vector{Float64} = Float64[],
                       units::Union{UnitScaling,Nothing} = nothing) -> String

Animate Lagrangian radii evolution with a sweeping time cursor.

# Arguments
- `fps`: frames per second (`nothing` = use `cfg.style.anim_fps`; `0` there
  auto-calculates for ~`cfg.style.anim_target_seconds` s, clamped to 2–30 fps).
  Lagrangian data is a dense time series so the auto rate allows smooth,
  fast playback.
- `units`: when `cfg.units == "physical"` and a `UnitScaling` is provided,
  times are shown in Myr and radii in pc; N-body units otherwise.

Returns the output file path.
"""
function animate_lagrangian(
    lagr::LagrangianData, cfg::VisualizationConfig;
    filename::AbstractString = "lagrangian_anim",
    fps::Union{Int,Nothing} = nothing,
    selected_fractions::Vector{Float64} = Float64[],
    units::Union{UnitScaling,Nothing} = nothing,
)::String
    isempty(lagr.time) && error("No Lagrangian data to animate")

    if isempty(selected_fractions)
        selected_fractions = [0.01, 0.1, 0.5, 0.9, 1.0]
    end

    nt = length(lagr.time)
    # Lagrangian data is a dense time series (often hundreds of steps)
    # with smooth playback; allow up to 30 fps.
    fps = something(fps, cfg.style.anim_fps)
    fps > 0 || (fps = _auto_fps(nt;
        target_duration = cfg.style.anim_target_seconds, min_fps = 2, max_fps = 30))

    physical = cfg.units == "physical" && units !== nothing &&
               units.rbar > 0 && units.tscale > 0
    ts       = physical ? lagr.time .* units.tscale : lagr.time
    r_scale  = physical ? units.rbar : 1.0

    # Closest available mass fractions and the plotted (positive, scaled)
    # radii — the log-axis tick range comes from the actual data extents.
    frac_indices = [argmin(abs.(lagr.mass_fractions .- f)) for f in selected_fractions]
    r_pos = [r * r_scale for idx in frac_indices
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

    # Pre-plot all lines (full data) in light gray as ghost background.
    # Non-positive radii (empty shells) are invalid on the log axis → NaN.
    for idx in frac_indices
        lines!(ax, ts, [r > 0 ? r * r_scale : NaN for r in @view lagr.radii[idx, :]];
            color = :gray82, linewidth = 1.0)
    end

    # Animated lines — use Point2f Observables to avoid x/y length mismatch
    frame_idx = Observable(1)

    for (ci, fidx) in enumerate(frac_indices)
        actual_frac = lagr.mass_fractions[fidx]
        pct_val = actual_frac * 100
        pct = isinteger(pct_val) ? @sprintf("%d", Int(pct_val)) : @sprintf("%.1f", pct_val)

        ys = [r > 0 ? r * r_scale : NaN for r in @view lagr.radii[fidx, :]]
        pts = @lift(Point2f.(ts[1:$frame_idx], ys[1:$frame_idx]))

        lines!(ax, pts;
            color = _OKABE_ITO[mod1(ci, length(_OKABE_ITO))],
            label = latexstring("$(pct)\\%"),
        )
    end

    # Vertical cursor line
    vlines!(ax, @lift(ts[$frame_idx]);
        color = :gray40, linestyle = :dash, linewidth = 1.0)

    # Static layout: the top legend row is added before record() starts
    if length(selected_fractions) ≥ 2
        _top_legend!(fig, ax; title = L"M(r)/M_\mathrm{tot}:")
    end

    outpath = _anim_output_path(cfg, filename)

    @info "Animating Lagrangian radii: $nt frames → $outpath  ($(fps) fps, ~$(round(Int, nt/fps)) s)"

    _backup_existing(outpath)
    record(fig, outpath, 1:nt; framerate = fps) do i
        frame_idx[] = i
    end

    @info "Lagrangian animation saved: $outpath  ($nt frames, $(fps) fps)"
    return outpath
end
