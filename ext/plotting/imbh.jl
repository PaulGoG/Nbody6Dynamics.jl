# =============================================================================
# IMBH figures: wandering about the stellar centres against the Brownian
# predictions, and the most bound companion per snapshot
# =============================================================================

"""Colour of the IMBH offsets about the stellar centre of mass (Okabe–Ito blue)."""
const _IMBH_CM_COLOUR = _OKABE_ITO[1]

"""Colour of the IMBH offsets about the stellar density centre (Okabe–Ito orange)."""
const _IMBH_DC_COLOUR = _OKABE_ITO[2]

"""Guide lines of the Brownian predictions: `(legend label, index into _OKABE_ITO)`."""
const _IMBH_PREDICTION_STYLE = (("Plummer (CHL02)", 3), ("Harmonic core", 4))

"""Colours of the most frequent companions are taken cyclically from the first seven Okabe–Ito entries; black is kept for the IMRI stars."""
const _IMBH_COMPANION_COLOURS = 7

"""Layout units one further legend row adds to the canvas."""
const _IMBH_LEGEND_ROW = 36

"""The samples `(x, y)` with both finite and `y > 0`, as `Float64` vectors (for logarithmic axes)."""
function _imbh_positive(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    keep = [isfinite(xi) && isfinite(yi) && yi > 0 for (xi, yi) in zip(x, y)]
    return (Float64.(x[keep]), Float64.(y[keep]))
end

"""
    _imbh_log_limits(values; bottom_dex = 0.12, top_dex = 0.25) -> (lo, hi)

Limits of a logarithmic axis about the positive `values`: margins of
`bottom_dex` below and `top_dex` above (room for the guide labels), opened
symmetrically to at least one decade.
"""
function _imbh_log_limits(
    values::AbstractVector{<:Real};
    bottom_dex::Real = 0.12,
    top_dex::Real = 0.25,
)
    lo = minimum(values) / 10^bottom_dex
    hi = maximum(values) * 10^top_dex
    span = log10(hi / lo)
    if span < 1
        widen = 10^((1 - span) / 2)
        lo /= widen
        hi *= widen
    end
    return (Float64(lo), Float64(hi))
end

"""
    _imbh_label_above(y, xs, ys, x_lo, x_hi) -> Bool

Whether the inline label of a horizontal guide line at `y` on a logarithmic
axis goes above the line at its right end: true unless the box above it (the
right 30 % of `[x_lo, x_hi]`, a factor 1.5 in `y`) holds more of the points
`(xs, ys)` than the box below.
"""
function _imbh_label_above(
    y::Real,
    xs::AbstractVector{<:Real},
    ys::AbstractVector{<:Real},
    x_lo::Real,
    x_hi::Real,
)::Bool
    x_cut = x_hi - 0.3 * (x_hi - x_lo)
    n_above = count(i -> xs[i] ≥ x_cut && y ≤ ys[i] ≤ 1.5 * y, eachindex(xs))
    n_below = count(i -> xs[i] ≥ x_cut && y / 1.5 ≤ ys[i] ≤ y, eachindex(xs))
    return n_above ≤ n_below
end

# ---------------------------------------------------------------------------
# Wandering
# ---------------------------------------------------------------------------

"""
    plot_imbh_wandering(traj::AbstractVector{IMBHState}, wandering::AbstractVector{<:NamedTuple},
                        cfg::VisualizationConfig; prediction_plummer_pc = NaN,
                        prediction_core_pc = NaN, filename = "imbh_wandering") -> String

Brownian wandering of the IMBH in two stacked panels sharing the time axis.
Upper panel, logarithmic: the distance of the IMBH from the stellar centre
of mass (`r_cm_pc`, blue) and, thinner, from the stellar density centre
(`r_pc`, orange), over the core radius ``r_\\mathrm{c}`` (grey, dashed) and
the half-mass radius ``r_\\mathrm{h}`` (grey, dotted) of the stars, each
labelled at its end; non-positive values are left out. Lower panel,
logarithmic: the rms offset per coordinate ``x_\\mathrm{rms} = \\sqrt{\\langle x^2 \\rangle}``
of every window of `wandering` as a marker at the window centre with a
horizontal bar over the window, in the centre-of-mass frame when the rows
hold it (the frame of Chatterjee, Hernquist & Loeb 2002,
DOI 10.1086/339174), in the density-centre frame otherwise, coloured as the
distance from that centre above. The equilibrium predictions are dashed
horizontal guide lines labelled at their right end, and the rms offset
pooled over all windows (weighted by the states per window) is stated in
the panel.

# Arguments
- `traj`: the states of [`imbh_trajectory`](@ref)
- `wandering`: the rows of [`imbh_wandering`](@ref), of either frame or both
- `prediction_plummer_pc`: [`brownian_rms_plummer`](@ref) [pc]; drawn when
  finite and positive
- `prediction_core_pc`: [`brownian_rms_core`](@ref) (Lin & Tremaine 1980,
  DOI 10.1086/157803) [pc]; drawn when finite and positive

An empty `wandering` leaves the lower panel with the predictions and the
note `No window`. Throws an error when `traj` holds no state with a finite
time.
"""
@publication function Nbody6Dynamics.plot_imbh_wandering(
    traj::AbstractVector{IMBHState},
    wandering::AbstractVector{<:NamedTuple},
    cfg::VisualizationConfig;
    prediction_plummer_pc::Real = NaN,
    prediction_core_pc::Real = NaN,
    filename::AbstractString = "imbh_wandering",
)::String
    states = sort!([s for s in traj if isfinite(s.time_myr)]; by = s -> s.time_myr)
    isempty(states) && error("No IMBH states to plot")
    t = [s.time_myr for s in states]
    tmin, tmax = t[1] == t[end] ? (t[1] - 0.5, t[end] + 0.5) : (t[1], t[end])
    tticks = _time_ticks(tmin, tmax)

    # Upper panel: distances and the radii of the stars
    t_cm, r_cm = _imbh_positive(t, [s.r_cm_pc for s in states])
    t_dc, r_dc = _imbh_positive(t, [s.r_pc for s in states])
    t_c, r_c = _imbh_positive(t, [s.r_core_pc for s in states])
    t_h, r_h = _imbh_positive(t, [s.r_half_pc for s in states])
    xs_upper = vcat(t_cm, t_dc, t_c, t_h)
    ys_upper = vcat(r_cm, r_dc, r_c, r_h)

    # Lower panel: the windows of one frame, the predictions
    frame = any(r -> r.frame === :centre_of_mass, wandering) ? :centre_of_mass : :density_centre
    windows = [
        r for r in wandering if r.frame === frame &&
            r.n > 0 &&
            isfinite(r.x2_mean_pc2) &&
            r.x2_mean_pc2 > 0 &&
            isfinite(r.t_start_myr) &&
            isfinite(r.t_end_myr)
    ]
    x_rms = Float64[sqrt(r.x2_mean_pc2) for r in windows]
    t_start = Float64[r.t_start_myr for r in windows]
    t_end = Float64[r.t_end_myr for r in windows]
    t_mid = (t_start .+ t_end) ./ 2
    c_frame = frame === :centre_of_mass ? _IMBH_CM_COLOUR : _IMBH_DC_COLOUR
    predictions = Tuple{Float64,String,Makie.RGBAf}[]
    for (v, (label, idx)) in
        zip((prediction_plummer_pc, prediction_core_pc), _IMBH_PREDICTION_STYLE)
        isfinite(v) && v > 0 || continue
        push!(predictions, (Float64(v), label, Makie.to_color(_OKABE_ITO[idx])))
    end

    # Legend content first: a second row is added to the canvas height
    elements = _LegendElement[]
    labels = AbstractString[]
    if !isempty(t_cm)
        push!(elements, LineElement(; color = _IMBH_CM_COLOUR, linewidth = _STYLE.data))
        push!(labels, "About c.m.")
    end
    if !isempty(t_dc)
        push!(elements, LineElement(; color = _IMBH_DC_COLOUR, linewidth = _STYLE.fit))
        push!(labels, "About density centre")
    end
    if !isempty(t_c)
        push!(elements, LineElement(; color = :gray45, linestyle = :dash, linewidth = _STYLE.guide))
        push!(labels, L"r_\mathrm{c}")
    end
    if !isempty(t_h)
        push!(elements, LineElement(; color = :gray45, linestyle = :dot, linewidth = _STYLE.guide))
        push!(labels, L"r_\mathrm{h}")
    end
    for (_, label, c) in predictions
        push!(elements, LineElement(; color = c, linestyle = :dash, linewidth = _STYLE.guide))
        push!(labels, label)
    end
    nbanks = length(elements) > 4 ? 2 : 1

    upper_values = isempty(ys_upper) ? vcat(x_rms, [p[1] for p in predictions]) : ys_upper
    # The first states sit at the centre (r → 0): the lower limit follows the first percentile
    y1_lo, y1_hi = if isempty(upper_values)
        (0.1, 1.0)
    else
        sorted_upper = sort(upper_values)
        lo_robust = sorted_upper[max(1, ceil(Int, 0.01 * length(sorted_upper)))]
        _imbh_log_limits([lo_robust, sorted_upper[end]])
    end
    lower_values = vcat(x_rms, [p[1] for p in predictions])
    y2_lo, y2_hi = isempty(lower_values) ? (y1_lo, y1_hi) : _imbh_log_limits(lower_values)

    frame_text = frame === :centre_of_mass ? "c.m. frame" : "density-centre frame"
    fig =
        Figure(; size = _fig_multipanel(cfg, 2, 1; extra_height = (nbanks - 1) * _IMBH_LEGEND_ROW))
    ax1 = Axis(
        fig[1, 1];
        ylabel = L"r\;[\mathrm{pc}]",
        yscale = log10,
        yticks = _log_ticks(y1_lo, y1_hi),
        xticks = tticks,
        xticklabelsvisible = false,
    )
    ax2 = Axis(
        fig[2, 1];
        xlabel = L"t\;[\mathrm{Myr}]",
        ylabel = latexstring("\$x_\\mathrm{rms}\\;[\\mathrm{pc}]\$, $(frame_text)"),
        yscale = log10,
        yticks = _log_ticks(y2_lo, y2_hi),
        xticks = tticks,
    )
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    # A right margin keeps the inline labels of the predictions clear of the windows
    x_hi = isempty(predictions) ? tmax : tmax + 0.16 * (tmax - tmin)
    xlims!(ax2, tmin, x_hi)
    ylims!(ax1, y1_lo, y1_hi)
    ylims!(ax2, y2_lo, y2_hi)

    # Radii of the stars beneath the distances (named in the legend)
    isempty(t_h) ||
        lines!(ax1, t_h, r_h; color = :gray45, linestyle = :dot, linewidth = _STYLE.guide)
    isempty(t_c) ||
        lines!(ax1, t_c, r_c; color = :gray45, linestyle = :dash, linewidth = _STYLE.guide)
    isempty(t_dc) || lines!(ax1, t_dc, r_dc; color = _IMBH_DC_COLOUR, linewidth = _STYLE.fit)
    isempty(t_cm) || lines!(ax1, t_cm, r_cm; color = _IMBH_CM_COLOUR, linewidth = _STYLE.data)

    # Predictions beneath the windows
    for (v, _, c) in predictions
        hlines!(ax2, [v]; color = c, linestyle = :dash, linewidth = _STYLE.guide)
    end
    if !isempty(windows)
        rangebars!(
            ax2,
            x_rms,
            t_start,
            t_end;
            direction = :x,
            color = c_frame,
            linewidth = _STYLE.fit,
            whiskerwidth = 10,
        )
        scatter!(
            ax2,
            t_mid,
            x_rms;
            color = c_frame,
            markersize = _STYLE.marker,
            strokecolor = _band_edge(c_frame),
            strokewidth = _STYLE.marker_stroke,
        )
    end
    win_x = vcat(t_start, t_mid, t_end)
    win_y = vcat(x_rms, x_rms, x_rms)
    # Inline labels at the right end; of two close predictions the lower is labelled beneath its line
    sorted = sort(predictions; by = first)
    close_pair = length(sorted) == 2 && sorted[2][1] / sorted[1][1] < 1.6
    for (k, (v, label, c)) in enumerate(sorted)
        up = close_pair ? k == 2 : _imbh_label_above(v, win_x, win_y, tmin, tmax)
        text!(
            ax2,
            x_hi,
            v;
            text = label,
            align = (:right, up ? :bottom : :top),
            offset = (-6, up ? 3 : -3),
            color = c,
            fontsize = _ANNOTATION_FONTSIZE,
        )
    end

    # Pooled rms offset in the emptiest corner, clear of the windows and the guide lines
    corner_x = copy(win_x)
    corner_y = log10.(win_y)
    for (v, _, _) in predictions, u in range(0, 1; length = 11)
        push!(corner_x, tmin + u * (tmax - tmin))
        push!(corner_y, log10(v))
    end
    corner = _emptiest_corner(
        corner_x,
        corner_y;
        width = 0.45,
        xlims = (tmin, tmax),
        ylims = (log10(y2_lo), log10(y2_hi)),
    )
    if isempty(windows)
        _annotate!(ax2, "No window"; corner = corner)
    else
        n_total = sum(r -> r.n, windows)
        pooled = sqrt(sum(r -> r.n * r.x2_mean_pc2, windows) / n_total)
        _annotate!(
            ax2,
            latexstring("\$x_\\mathrm{rms} = $(_fmt_latex_sig(pooled, 3))\$ pc (all windows)");
            corner = corner,
            color = c_frame,
        )
    end

    length(elements) ≥ 2 && _top_legend!(fig, elements, labels; nbanks = nbanks)

    return _save_fig(cfg, filename, fig)
end

# ---------------------------------------------------------------------------
# Bound companions
# ---------------------------------------------------------------------------

"""
    plot_imbh_companions(rows::AbstractVector{IMBHCompanion}, cfg::VisualizationConfig;
                         max_named = 6, filename = "imbh_companions") -> String

The most bound two-body companion of the IMBH per snapshot in two stacked
panels sharing the time axis: the semi-major axis `a` [au] above
(logarithmic), the eccentricity `e` below. One marker per snapshot with a
bound companion (`companion_name > 0`); the colour encodes the companion's
identity — the `max_named` most frequent companions in Okabe–Ito order of
frequency, labelled `Body <name>`, every other one grey as `Other` — and the
fill its regularisation: filled for the KS partner of the IMBH (`in_pair`),
open-faced otherwise. A black star on top marks an IMRI candidate, a pair
whose Peters coalescence time `t_gw_myr` is shorter than the time span of
the rows. The upper panel states the bound fraction, the number of
snapshots and the exchanges of [`imbh_companion_summary`](@ref).

# Arguments
- `rows`: the rows of [`imbh_companions`](@ref)
- `max_named`: most companions given their own colour and legend entry

Without a bound row the axes are drawn empty with the note
`No bound companion`. Throws an error for an empty `rows` and an
`ArgumentError` for `max_named < 0`.
"""
@publication function Nbody6Dynamics.plot_imbh_companions(
    rows::AbstractVector{IMBHCompanion},
    cfg::VisualizationConfig;
    max_named::Integer = 6,
    filename::AbstractString = "imbh_companions",
)::String
    isempty(rows) && error("No companion rows to plot")
    max_named ≥ 0 || throw(ArgumentError("max_named must be ≥ 0; got $max_named"))
    times = Float64[r.time_myr for r in rows if isfinite(r.time_myr)]
    t_lo, t_hi = isempty(times) ? (0.0, 1.0) : extrema(times)
    tmin, tmax = t_lo == t_hi ? (t_lo - 0.5, t_hi + 0.5) : (t_lo, t_hi)
    tticks = _time_ticks(tmin, tmax)
    t_span = t_hi - t_lo
    bound = IMBHCompanion[r for r in rows if r.companion_name > 0 && isfinite(r.time_myr)]

    # Companions ranked by frequency, ties by name
    counts = Dict{Int,Int}()
    for r in bound
        counts[r.companion_name] = get(counts, r.companion_name, 0) + 1
    end
    ranked = sort!(collect(keys(counts)); by = nm -> (-counts[nm], nm))
    named = ranked[1:min(Int(max_named), length(ranked))]
    rank_of = Dict(nm => k for (k, nm) in enumerate(named))
    grey = Makie.to_color(:gray55)
    colour(k::Int) = k == 0 ? grey : Makie.to_color(_OKABE_ITO[mod1(k, _IMBH_COMPANION_COLOURS)])
    imri = IMBHCompanion[r for r in bound if isfinite(r.t_gw_myr) && r.t_gw_myr < t_span]

    # Legend content first: its rows beyond the first are added to the canvas height
    companion_elements = _LegendElement[]
    companion_labels = AbstractString[]
    filled(c) = MarkerElement(;
        color = c,
        marker = :circle,
        markersize = _STYLE.marker,
        strokecolor = _band_edge(c),
        strokewidth = _STYLE.marker_stroke,
    )
    for (k, nm) in enumerate(named)
        push!(companion_elements, filled(colour(k)))
        push!(companion_labels, "Body $(nm)")
    end
    if length(ranked) > length(named)
        push!(companion_elements, filled(grey))
        push!(companion_labels, "Other")
    end
    state_elements = _LegendElement[]
    state_labels = AbstractString[]
    if any(r -> r.in_pair, bound)
        push!(state_elements, filled(grey))
        push!(state_labels, "KS pair")
    end
    if any(r -> !r.in_pair, bound)
        push!(
            state_elements,
            MarkerElement(;
                color = :white,
                marker = :circle,
                markersize = _STYLE.marker,
                strokecolor = grey,
                strokewidth = _STYLE.marker_stroke,
            ),
        )
        push!(state_labels, "Unregularised")
    end
    if !isempty(imri)
        push!(
            state_elements,
            MarkerElement(;
                color = :black,
                marker = :star5,
                markersize = _STYLE.marker,
                strokewidth = 0,
            ),
        )
        push!(state_labels, "IMRI candidate")
    end
    show_legend = length(companion_elements) + length(state_elements) ≥ 2
    companion_banks = max(cld(length(companion_elements), 4), 1)
    legend_rows =
        show_legend ?
        (isempty(state_elements) ? 0 : 1) + (isempty(companion_elements) ? 0 : companion_banks) : 0

    a_points = [
        (r.time_myr, r.semi_major_axis_au) for
        r in bound if isfinite(r.semi_major_axis_au) && r.semi_major_axis_au > 0
    ]
    a_lo, a_hi =
        isempty(a_points) ? (1.0, 1.0e3) :
        _imbh_log_limits([p[2] for p in a_points]; bottom_dex = 0.2, top_dex = 0.3)

    fig = Figure(;
        size = _fig_multipanel(
            cfg,
            2,
            1;
            extra_height = max(legend_rows - 1, 0) * _IMBH_LEGEND_ROW,
        ),
    )
    ax1 = Axis(
        fig[1, 1];
        ylabel = L"a\;[\mathrm{au}]",
        yscale = log10,
        yticks = _log_ticks(a_lo, a_hi),
        xticks = tticks,
        xticklabelsvisible = false,
    )
    ax2 = Axis(
        fig[2, 1];
        xlabel = L"t\;[\mathrm{Myr}]",
        ylabel = L"e",
        yticks = 0.0:0.2:1.0,
        xticks = tticks,
    )
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    xlims!(ax2, tmin, tmax)
    ylims!(ax1, a_lo, a_hi)
    ylims!(ax2, -0.04, 1.04)

    if isempty(bound)
        _no_data_note!(ax1, "No bound companion")
        return _save_fig(cfg, filename, fig)
    end

    # Markers: the grey companions beneath, the named ones on top
    for k in 0:length(named), paired in (true, false)
        sub = IMBHCompanion[
            r for r in bound if get(rank_of, r.companion_name, 0) == k && r.in_pair == paired
        ]
        isempty(sub) && continue
        c = colour(k)
        kw = (
            marker = :circle,
            markersize = _STYLE.marker,
            color = paired ? c : Makie.to_color(:white),
            strokecolor = paired ? _band_edge(c) : c,
            strokewidth = _STYLE.marker_stroke,
        )
        ta, a = _imbh_positive([r.time_myr for r in sub], [r.semi_major_axis_au for r in sub])
        isempty(ta) || scatter!(ax1, ta, a; kw...)
        te = Float64[r.time_myr for r in sub if isfinite(r.eccentricity)]
        e = Float64[r.eccentricity for r in sub if isfinite(r.eccentricity)]
        isempty(te) || scatter!(ax2, te, e; kw...)
    end
    if !isempty(imri)
        star = (marker = :star5, color = :black, markersize = 0.75 * _STYLE.marker, strokewidth = 0)
        ta, a = _imbh_positive([r.time_myr for r in imri], [r.semi_major_axis_au for r in imri])
        isempty(ta) || scatter!(ax1, ta, a; star...)
        te = Float64[r.time_myr for r in imri if isfinite(r.eccentricity)]
        e = Float64[r.eccentricity for r in imri if isfinite(r.eccentricity)]
        isempty(te) || scatter!(ax2, te, e; star...)
    end

    bh_summary = imbh_companion_summary(rows)
    n_snap = bh_summary.n_snapshots
    n_exch = bh_summary.n_exchanges
    summary_text =
        "Bound in $(_fmt_plain_sig(100 * bh_summary.fraction_bound, 3)) % of $(n_snap) " *
        (n_snap == 1 ? "snapshot" : "snapshots") *
        ", $(n_exch) " *
        (n_exch == 1 ? "exchange" : "exchanges")
    corner = _emptiest_corner(
        Float64[p[1] for p in a_points],
        Float64[log10(p[2]) for p in a_points];
        width = 0.55,
        xlims = (tmin, tmax),
        ylims = (log10(a_lo), log10(a_hi)),
    )
    _annotate!(ax1, summary_text; corner = corner)

    # Two legend rows on top: the states over the companions
    if show_legend
        if isempty(state_elements) || isempty(companion_elements)
            only_state = isempty(companion_elements)
            _top_legend!(
                fig,
                only_state ? state_elements : companion_elements,
                only_state ? state_labels : companion_labels;
                title = only_state ? "State:" : "Companion:",
                nbanks = only_state ? 1 : companion_banks,
            )
        else
            _top_legend!(fig, state_elements, state_labels; title = "State:")
            Legend(
                fig[-1, :],
                companion_elements,
                companion_labels,
                "Companion:";
                orientation = :horizontal,
                nbanks = companion_banks,
                framevisible = false,
                titleposition = :left,
                tellheight = true,
                padding = (0, 0, 0, 0),
            )
        end
    end

    return _save_fig(cfg, filename, fig)
end
