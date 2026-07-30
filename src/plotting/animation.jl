# =============================================================================
# Animation of cluster evolution — GIF output
# =============================================================================

"""
    _auto_fps(nframes; target_duration, min_fps, max_fps) -> Int

Compute a frame rate so the animation lasts approximately `target_duration`
seconds, clamped to `[min_fps, max_fps]`.
"""
function _auto_fps(nframes::Int;
                   target_duration::Float64, min_fps::Int, max_fps::Int)::Int
    return clamp(round(Int, nframes / target_duration), min_fps, max_fps)
end

"""
    animate_cluster(snaps::Vector{Snapshot}, cfg::VisualizationConfig;
                    filename::AbstractString = "cluster_evolution",
                    projections::Vector{Symbol} = [:xy, :xz, :yz],
                    fps::Union{Int,Nothing} = nothing,
                    trail_frac::Float64 = 0.0) -> Vector{String}

Create animated GIFs of the cluster's spatial evolution across snapshots,
one per projection.

# Arguments
- `snaps`: ordered vector of `Snapshot`s
- `cfg`: visualisation configuration (figsize, dpi, output_dir)
- `filename`: output filename stem (`.gif` appended automatically)
- `projections`: spatial projections to animate (any of :xy, :xz, :yz)
- `fps`: frames per second (`nothing` = auto-calculate for ~12 s total duration,
  clamped to 1–10 fps; manual `Int` overrides)
- `trail_frac`: fraction of previous positions to overlay as a fading trail
  (0.0 = no trail, 0.3 = overlay last 30% of elapsed frames)

Returns a vector of output file paths.
"""
function animate_cluster(
    snaps::Vector{Snapshot}, cfg::VisualizationConfig;
    filename::AbstractString = "cluster_evolution",
    projections::Vector{Symbol} = [:xy, :xz, :yz],
    fps::Union{Int,Nothing} = nothing,
    trail_frac::Float64 = 0.0,
)::Vector{String}
    isempty(snaps) && error("No snapshots to animate")

    nframes = length(snaps)
    # Cluster snapshots are discrete, each worth studying individually.
    # Target ~12 s total; fewer snapshots → slower pace, many → speed up.
    fps = something(fps,
        _auto_fps(nframes; target_duration = 12.0, min_fps = 1, max_fps = 10))

    # Marker size: visible but not overlapping
    n_typical = nparticles(snaps[1])
    ms = clamp(18000 / n_typical, 4.0, 20.0)

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
        ix, iy, xlab, ylab = _proj_indices(projection)

        # Check if adaptive zoom is needed
        all_indices = collect(1:nframes)
        use_adaptive = _needs_adaptive_zoom(snaps, all_indices, ix, iy)

        # Pre-compute per-frame limits for adaptive mode
        frame_limits = if use_adaptive
            [_square_limits(snaps[i].pos[ix, :], snaps[i].pos[iy, :]) for i in 1:nframes]
        else
            # Global limits — same for every frame
            all_x = reduce(vcat, [snap.pos[ix, :] for snap in snaps])
            all_y = reduce(vcat, [snap.pos[iy, :] for snap in snaps])
            lims = _square_limits(all_x, all_y)
            fill(lims, nframes)
        end

        outpath = _anim_output_path(cfg, "$(filename)_$(projection)")

        mode_str = use_adaptive ? "adaptive zoom" : "global limits"
        @info "Animating $nframes frames ($projection, $mode_str) → $outpath  ($(fps) fps, ~$(round(Int, nframes/fps)) s)"

        # Observable for the frame index
        frame_idx = Observable(1)

        # Title as observable
        title_text = @lift begin
            snap = snaps[$frame_idx]
            n = nparticles(snap)
            t_str = @sprintf("%.4f", time_nb(snap.header))
            latexstring("\\mathrm{N} = $(n), \\;\\; \\mathrm{t} = $(t_str) \\; \\mathrm{[NB]}")
        end

        # Initial limits
        xlo0, xhi0, ylo0, yhi0 = frame_limits[1]

        fig = Figure(; size = _fig_with_colorbar(cfg))
        ax = Axis(fig[1, 1];
            xlabel = xlab,
            ylabel = ylab,
            title  = title_text,
            aspect = DataAspect(),
            limits = (xlo0, xhi0, ylo0, yhi0),
            xticks = _nice_ticks(xlo0, xhi0),
            yticks = _nice_ticks(ylo0, yhi0),
        )

        # Use Point2f observable to handle varying particle counts across frames
        pts = @lift(Point2f.(snaps[$frame_idx].pos[ix, :], snaps[$frame_idx].pos[iy, :]))
        colors = @lift(log10.(max.(Float64.(snaps[$frame_idx].mass), 1e-30)))

        scatter!(ax, pts;
            color      = colors,
            colormap   = :viridis,
            colorrange = (cmin, cmax),
            markersize = ms,
            strokewidth = 0,
        )

        # Add colorbar
        Colorbar(fig[1, 2]; colormap = :viridis, colorrange = (cmin, cmax),
                 label = L"\log_{10}(\mathrm{m} \, / \, \mathrm{M}_\mathrm{tot})",
                 ticks = _nice_colorbar_ticks(cmin, cmax))
        colgap!(fig.layout, 10)

        record(fig, outpath, 1:nframes; framerate = fps) do i
            frame_idx[] = i
            if use_adaptive
                xlo, xhi, ylo, yhi = frame_limits[i]
                xlims!(ax, xlo, xhi)
                ylims!(ax, ylo, yhi)
                ax.xticks = _nice_ticks(xlo, xhi)
                ax.yticks = _nice_ticks(ylo, yhi)
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
- `fps`: frames per second (`nothing` = auto-calculate for ~15 s total duration,
  clamped to 1–8 fps; manual `Int` overrides). HR frames are information-dense,
  so the auto rate favours a slower pace than cluster animations.

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
    # HR frames are information-dense (stellar types, population structure).
    # Target ~15 s; keep pace slow so each epoch is readable.
    fps = something(fps,
        _auto_fps(nframes; target_duration = 15.0, min_fps = 1, max_fps = 8))

    outpath = _anim_output_path(cfg, filename)

    @info "Animating HR diagram: $nframes frames → $outpath  ($(fps) fps, ~$(round(Int, nframes/fps)) s)"

    frame_idx = Observable(1)

    title_text = @lift begin
        sev = sevs[$frame_idx]
        t_str = @sprintf("%.4f", sev.time_myr)
        latexstring("\\mathrm{t}_\\mathrm{NB} = $(t_str), \\;\\; \\mathrm{N}_\\star = $(sev.n_stars)")
    end

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(fig[1, 1];
        xlabel = L"\log_{10}(\mathrm{T}_\mathrm{eff} \, / \, \mathrm{K})",
        ylabel = L"\log_{10}(\mathrm{L} \, / \, \mathrm{L}_\odot)",
        title  = title_text,
        limits = (xlims..., ylims...),
        xreversed = true,
        xticks = _logval_ticks(xlims[1], xlims[2]),
        yticks = _logval_ticks(ylims[1], ylims[2]),
    )

    # Use Point2f observable to handle varying star counts across epochs
    pts = @lift(Point2f.(
        [r.log_teff for r in valid_per_frame[$frame_idx]],
        [r.log_luminosity for r in valid_per_frame[$frame_idx]]))
    col_data = @lift([get(_HR_COLORS, Int(r.stellar_type), :gray50)
                      for r in valid_per_frame[$frame_idx]])

    scatter!(ax, pts;
        color = col_data, markersize = 14, strokewidth = 0)

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
                       selected_fractions::Vector{Float64} = Float64[]) -> String

Animate Lagrangian radii evolution with a sweeping time cursor.

# Arguments
- `fps`: frames per second (`nothing` = auto-calculate for ~12 s total duration,
  clamped to 2–30 fps; manual `Int` overrides). Lagrangian data is a dense
  time series so the auto rate allows smooth, fast playback.

Returns the output file path.
"""
function animate_lagrangian(
    lagr::LagrangianData, cfg::VisualizationConfig;
    filename::AbstractString = "lagrangian_anim",
    fps::Union{Int,Nothing} = nothing,
    selected_fractions::Vector{Float64} = Float64[],
)::String
    isempty(lagr.time) && error("No Lagrangian data to animate")

    if isempty(selected_fractions)
        selected_fractions = [0.01, 0.1, 0.5, 0.9, 1.0]
    end

    nt = length(lagr.time)
    # Lagrangian data is a dense time series (often hundreds of steps).
    # Target ~12 s with smooth playback; allow up to 30 fps.
    fps = something(fps,
        _auto_fps(nt; target_duration = 12.0, min_fps = 2, max_fps = 30))

    fig = Figure(; size = _figsize_px(cfg))
    ttk = _time_ticks(first(lagr.time), last(lagr.time))
    ax = Axis(fig[1, 1];
        xlabel = L"\mathrm{t} \; \mathrm{[NB]}",
        ylabel = L"\mathrm{r}_\mathrm{L} \; \mathrm{[NB]}",
        title  = L"\textbf{Lagrangian Radii Evolution}",
        yscale = log10,
        yminorticksvisible = false,
        xticks = ttk,
    )

    colors = Makie.wong_colors()

    # Pre-plot all lines (full data) in light gray as ghost background
    frac_indices = Int[]
    for frac in selected_fractions
        idx = argmin(abs.(lagr.mass_fractions .- frac))
        push!(frac_indices, idx)
        lines!(ax, lagr.time, lagr.radii[idx, :];
            color = :gray82, linewidth = 1.0)
    end

    # Animated lines — use Point2f Observables to avoid x/y length mismatch
    frame_idx = Observable(1)

    for (ci, (frac, fidx)) in enumerate(zip(selected_fractions, frac_indices))
        actual_frac = lagr.mass_fractions[fidx]
        pct_val = actual_frac * 100
        pct = isinteger(pct_val) ? @sprintf("%d", Int(pct_val)) : @sprintf("%.1f", pct_val)

        pts = @lift(Point2f.(lagr.time[1:$frame_idx], lagr.radii[fidx, 1:$frame_idx]))

        lines!(ax, pts;
            color = colors[mod1(ci, length(colors))],
            label = latexstring("\\mathrm{M}(\\mathrm{r})/\\mathrm{M}_\\mathrm{tot} = $(pct)\\%"),
        )
    end

    # Vertical cursor line
    vlines!(ax, @lift(lagr.time[$frame_idx]);
        color = :gray40, linestyle = :dash, linewidth = 1.0)

    axislegend(ax; position = :rt, framevisible = true,
               backgroundcolor = (:white, 0.7))

    outpath = _anim_output_path(cfg, filename)

    @info "Animating Lagrangian radii: $nt frames → $outpath  ($(fps) fps, ~$(round(Int, nt/fps)) s)"

    record(fig, outpath, 1:nt; framerate = fps) do i
        frame_idx[] = i
    end

    @info "Lagrangian animation saved: $outpath  ($nt frames, $(fps) fps)"
    return outpath
end
