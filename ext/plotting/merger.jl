# =============================================================================
# Merger-specific diagnostics
# =============================================================================
# Track each initial cluster through its own bound set over the course of the
# run. Gives meaningful pre-coalescence diagnostics for multi-cluster setups
# where global Lagrangian radii measure cluster separation rather than
# internal structure. The merger epochs decide where these per-cluster
# quantities hold: the series end where a cluster is retired and the regimes
# of the system are shaded.

# -----------------------------------------------------------------------------
# Merger epochs on a time axis: regime bands and epoch markers
# -----------------------------------------------------------------------------

"""
Shading of the system regimes of the merger epochs on a time axis, in their
order: regime, legend label and fill. Low-alpha greys, the `coalescing` band
the darker, so that the bands survive greyscale and never compete with the
colours of the series; `detached` and `settled` stay unshaded.
"""
const _REGIME_BANDS = (
    (:interacting, "Interacting", Makie.RGBAf(0.55, 0.55, 0.55, 0.16)),
    (:coalescing, "Coalescing", Makie.RGBAf(0.2, 0.2, 0.2, 0.18)),
)

"""
    _regime_bands!(axes, t, regime; label = true) -> Vector{Tuple{String,Makie.RGBAf}}

Full-height bands behind the series over the `interacting` and `coalescing`
snapshots of `regime` (one regime per time of `t`) on every axis of `axes`.
A run of equal regimes spans from its first snapshot to the first snapshot
of the next run, or to the last snapshot for the final run. With `label` the
first band of each regime on the first axis carries the legend label.
Returns the labels and fills of the regimes drawn, in order of appearance,
for legends assembled by hand ([`_regime_legend`](@ref)).
"""
function _regime_bands!(
    axes,
    t::AbstractVector{<:Real},
    regime::AbstractVector{Symbol};
    label::Bool = true,
)
    n = length(t)
    n == length(regime) ||
        throw(DimensionMismatch("_regime_bands!: $(length(regime)) regimes for $n times"))
    drawn = Tuple{String,Makie.RGBAf}[]
    k = 1
    while k ≤ n
        r = regime[k]
        j = k
        while j < n && regime[j + 1] === r
            j += 1
        end
        x_lo = Float64(t[k])
        x_hi = Float64(t[j < n ? j + 1 : n])
        b = findfirst(e -> e[1] === r, _REGIME_BANDS)
        if b !== nothing && x_hi > x_lo
            _, name, fill = _REGIME_BANDS[b]
            new = all(d -> d[1] != name, drawn)
            for (a, ax) in enumerate(axes)
                p = vspan!(
                    ax,
                    [x_lo],
                    [x_hi];
                    color = fill,
                    label = (label && new && a == 1) ? name : nothing,
                )
                # Behind every series of the axis
                translate!(p, 0, 0, -10)
            end
            new && push!(drawn, (name, fill))
        end
        k = j + 1
    end
    return drawn
end

"""Legend entries (elements, labels) of the regime bands returned by [`_regime_bands!`](@ref)."""
function _regime_legend(drawn)
    elements = _LegendElement[
        PolyElement(; color = f, strokecolor = _band_edge(:gray50), strokewidth = 1) for
        (_, f) in drawn
    ]
    return elements, AbstractString[name for (name, _) in drawn]
end

"""Width a vertical epoch label takes beside its line, in layout units: one annotation text line with clearance."""
const _MARKER_LABEL_WIDTH = 1.2 * _ANNOTATION_FONTSIZE

"""Width of the time axis of a single-panel figure in layout units (the canvas less the axis protrusions), for spacing the epoch labels."""
const _MARKER_AXIS_WIDTH = 780

"""Height of the axis of a single panel in layout units, for the length of the epoch labels in data units."""
const _MARKER_AXIS_HEIGHT = 340

"""The mathematics of a `LaTeXString` without its `\$` delimiters; plain text wrapped in `\\mathrm`."""
function _math_body(s::AbstractString)::String
    str = String(s)
    return startswith(str, "\$") && endswith(str, "\$") ? str[2:(end - 1)] :
           "\\mathrm{" * replace(str, " " => "\\;") * "}"
end

"""Rendered length of a label in layout units, estimated from the characters left after its LaTeX markup at the annotation size (an upper bound: subscripts are smaller)."""
function _label_length(s::AbstractString)::Float64
    str = replace(String(s), r"\\(mathrm|mathbf|text)" => "", r"\\[;,]" => " ", r"[{}$_^\\]" => "")
    return 0.45 * _ANNOTATION_FONTSIZE * length(str)
end

"""
    _occupied_intervals(series, x_lo, x_hi, pad) -> Vector{Tuple{Float64,Float64}}

The y-intervals the polylines `series` (pairs of x and y vectors, non-finite
points breaking nothing but skipped) cover within `x_lo ≤ x ≤ x_hi`, each
widened by `pad`, merged and sorted.
"""
function _occupied_intervals(series, x_lo::Real, x_hi::Real, pad::Real)
    raw = Tuple{Float64,Float64}[]
    for (xs, ys) in series
        idx = [k for k in eachindex(xs, ys) if isfinite(xs[k]) && isfinite(ys[k])]
        for (a, b) in zip(idx, Iterators.drop(idx, 1))
            xa, xb, ya, yb = xs[a], xs[b], ys[a], ys[b]
            (max(xa, xb) < x_lo || min(xa, xb) > x_hi) && continue
            at(x) = xb == xa ? ya : ya + (yb - ya) * (x - xa) / (xb - xa)
            u, v =
                at(clamp(x_lo, min(xa, xb), max(xa, xb))), at(clamp(x_hi, min(xa, xb), max(xa, xb)))
            push!(raw, (min(u, v) - pad, max(u, v) + pad))
        end
        for k in idx
            x_lo ≤ xs[k] ≤ x_hi && push!(raw, (ys[k] - pad, ys[k] + pad))
        end
    end
    sort!(raw)
    merged = Tuple{Float64,Float64}[]
    for (a, b) in raw
        if !isempty(merged) && a ≤ merged[end][2]
            merged[end] = (merged[end][1], max(merged[end][2], b))
        else
            push!(merged, (a, b))
        end
    end
    return merged
end

"""
    _epoch_markers!(ax, markers, y_range, t_range; series = (), axis_height = _MARKER_AXIS_HEIGHT)

Vertical lines at the epochs `markers`, tuples `(x, label, color, linestyle)`
(entries with a non-finite `x` skipped), each labelled vertically beside its
line, reading upwards, so that epochs a label width apart never collide. A
label sits left of its line, or right of it when the line lies within one
label width of the start of `t_range`; lines closer than one label width
share one label, their labels joined by commas. The label hangs from the top
of the axis (`y_range`, data coordinates) unless the polylines `series` (x
and y vectors as drawn) cross that stretch beside the line: it then hangs from
the top of the highest gap between the series long enough to hold it
(estimated for an axis of `axis_height` layout units). `label` is a
`LaTeXString`.
"""
function _epoch_markers!(
    ax,
    markers,
    y_range::Tuple{<:Real,<:Real},
    t_range::Tuple{<:Real,<:Real};
    series = (),
    axis_height::Real = _MARKER_AXIS_HEIGHT,
)
    ms = sort!([m for m in markers if isfinite(m[1])]; by = first)
    isempty(ms) && return nothing
    span = max(Float64(t_range[2] - t_range[1]), eps(Float64))
    gap = _MARKER_LABEL_WIDTH / _MARKER_AXIS_WIDTH * span
    y_lo, y_hi = Float64.(y_range)
    y_span = max(y_hi - y_lo, eps(Float64))
    for (x, _, color, style) in ms
        vlines!(ax, [Float64(x)]; color = color, linestyle = style, linewidth = _STYLE.guide)
    end
    T = eltype(ms)
    groups = [T[ms[1]]]
    for m in ms[2:end]
        m[1] - last(groups[end])[1] < gap ? push!(groups[end], m) : push!(groups, T[m])
    end
    for g in groups
        left = first(g)[1] - t_range[1] ≥ gap
        x = Float64(left ? first(g)[1] : last(g)[1])
        text =
            length(g) == 1 ? first(g)[2] : latexstring(join((_math_body(m[2]) for m in g), ",\\;"))
        # Highest gap of the series beside the line that holds the label
        len = (_label_length(text) + 8) / axis_height * y_span
        # The text column beside the line; a series that only touches the line is clear of it
        occupied = _occupied_intervals(
            series,
            left ? x - gap : x + 0.1 * gap,
            left ? x - 0.1 * gap : x + gap,
            0.02 * y_span,
        )
        top = y_hi
        for (a, b) in Iterators.reverse(occupied)
            if b ≥ top
                top = min(top, a)
            elseif top - max(b, y_lo) ≥ len
                break
            else
                top = a
            end
        end
        y_anchor = top - y_lo ≥ len ? top : y_hi
        text!(
            ax,
            x,
            y_anchor;
            text = text,
            rotation = π / 2,
            align = (:right, left ? :bottom : :top),
            offset = (left ? -3 : 3, -4),
            fontsize = _ANNOTATION_FONTSIZE,
            color = first(g)[3],
        )
    end
    return nothing
end

"""Height a second row of a top legend (bold family titles over one row of entries) adds to the canvas, in layout units."""
const _LEGEND_ROW_HEIGHT = 70

"""Time limits of an axis over the whole run (Makie's default 5 % margin), whatever span its series cover."""
function _time_limits!(ax, t::AbstractVector{<:Real})
    lo, hi = Float64(first(t)), Float64(last(t))
    pad = 0.05 * max(hi - lo, eps(Float64))
    xlims!(ax, lo - pad, hi + pad)
    return nothing
end

"""Line colours and styles of the epoch markers: contact, core coalescence, merger-tree join."""
const _EPOCH_STYLE =
    (contact = (:gray40, :dash), core = (:black, :dash), join = (:gray45, :dashdot))

"""Merger-tree groups of event `e` as `a{+}b`, tags of a group comma-separated."""
_event_groups(ep::MergerEpochs, e::Int)::String =
    join(ep.group_a[e], ",") * "{+}" * join(ep.group_b[e], ",")

"""Times of the epochs: N-body, or Myr when `physical`; `NaN` for index 0."""
_epoch_time(ep::MergerEpochs, k::Integer, physical::Bool)::Float64 =
    k > 0 ? (physical ? ep.time_myr[k] : ep.time[k]) : NaN

"""
    _separation_markers(ep, physical) -> Vector

Epoch markers of the separation figure: the first contact and the core
coalescence of the system with their times, and every merger-tree event: its
join (the centre coalescence, labelled `t_coalesce` with its time for two
clusters) and, for three or more clusters, its core coalescence unless it is
the system's.
"""
function _separation_markers(ep::MergerEpochs, physical::Bool)
    unit = physical ? "\\;\\mathrm{Myr}" : "\\;[\\mathrm{NB}]"
    valued(name, k) =
        latexstring("t_\\mathrm{$name} = ", _fmt_latex_sig3(_epoch_time(ep, k, physical)), unit)
    markers = Tuple{Float64,LaTeXString,Symbol,Symbol}[]
    ep.contact_index > 0 && push!(
        markers,
        (
            _epoch_time(ep, ep.contact_index, physical),
            valued("contact", ep.contact_index),
            _EPOCH_STYLE.contact...,
        ),
    )
    ep.core_index > 0 && push!(
        markers,
        (
            _epoch_time(ep, ep.core_index, physical),
            valued("core", ep.core_index),
            _EPOCH_STYLE.core...,
        ),
    )
    two = size(ep.tag_regime, 1) == 2
    for e in eachindex(ep.event_join)
        k = ep.event_join[e]
        if k > 0
            label =
                two ? valued("coalesce", k) :
                latexstring("t_\\mathrm{join}\\,(", _event_groups(ep, e), ")")
            push!(markers, (_epoch_time(ep, k, physical), label, _EPOCH_STYLE.join...))
        end
        c = ep.event_core[e]
        if !two && c > 0 && c != ep.core_index
            push!(
                markers,
                (
                    _epoch_time(ep, c, physical),
                    latexstring("t_\\mathrm{core}\\,(", _event_groups(ep, e), ")"),
                    _EPOCH_STYLE.core...,
                ),
            )
        end
    end
    return markers
end

"""
    _distinct_clusters(ep::MergerEpochs) -> Vector{Int}

Number of distinct clusters at every snapshot under the group rule of
[`coalescence_time`](@ref), rebuilt from the separation and threshold
matrices of the merger epochs: the connected components of the pairs closer
than their threshold, over the clusters with a finite separation to another
cluster (a lone present cluster is not counted).
"""
function _distinct_clusters(ep::MergerEpochs)::Vector{Int}
    n_cl = size(ep.tag_regime, 1)
    n_t = length(ep.time)
    parent = collect(1:n_cl)
    function root(a::Int)
        while parent[a] != a
            parent[a] = parent[parent[a]]
            a = parent[a]
        end
        return a
    end
    present = falses(n_cl)
    out = zeros(Int, n_t)
    for k in 1:n_t
        parent .= 1:n_cl
        fill!(present, false)
        for (q, (i, j)) in enumerate(ep.pairs)
            d = ep.separation[q, k]
            isfinite(d) || continue
            present[i] = true
            present[j] = true
            if d < ep.threshold[q, k]
                ri, rj = root(i), root(j)
                ri == rj || (parent[ri] = rj)
            end
        end
        out[k] = count(i -> present[i] && root(i) == i, 1:n_cl)
    end
    return out
end

"""
    _pair_noise_floor(snaps, cluster_ranges, ep, bound_sets; bound = (;)) -> Matrix{Float64}

Centre-noise floor of the distance of every pair of clusters from its core
coalescence (`pair_core` of the epochs) on, `n_pairs × n_t`, `NaN` elsewhere:
the RMS distance `√(3 (σ_i² + σ_j²))` of two density centres that coincide,
with `σ` the per-axis noise of the shrinking-sphere centre of a cluster's
members in the bound set of the whole snapshot (the centres of
[`coalescence_time`](@ref), default centre keys; the `sigma` of
[`shrinking_sphere_centre`](@ref), a lower bound of the scatter). The bound
sets come from `bound_sets`, or are evaluated with the keywords `bound` at
the snapshots needed. N-body units.
"""
function _pair_noise_floor(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    ep::MergerEpochs,
    bound_sets::Union{Nothing,AbstractVector{BoundSet}};
    bound::NamedTuple = (;),
)::Matrix{Float64}
    n_pairs = length(ep.pairs)
    n_t = length(snaps)
    out = fill(NaN, n_pairs, n_t)
    k0 = minimum(filter(>(0), ep.pair_core); init = n_t + 1)
    for k in k0:n_t
        merged = [p for p in 1:n_pairs if 0 < ep.pair_core[p] ≤ k && isfinite(ep.separation[p, k])]
        isempty(merged) && continue
        snap = snaps[k]
        stars = bound_sets === nothing ? bound_set(snap; bound...).stars : bound_sets[k].stars
        in_bound = falses(nparticles(snap))
        in_bound[stars] .= true
        σ = fill(NaN, length(cluster_ranges))
        for p in merged, i in ep.pairs[p]
            isnan(σ[i]) || continue
            sel =
                [s for s in Nbody6Dynamics._member_indices(snap, cluster_ranges[i]) if in_bound[s]]
            length(sel) ≥ 3 || continue
            σ[i] =
                shrinking_sphere_centre(Float64.(snap.pos[:, sel]), Float64.(snap.mass[sel])).sigma
        end
        for p in merged
            i, j = ep.pairs[p]
            out[p, k] = sqrt(3 * (σ[i]^2 + σ[j]^2))
        end
    end
    return out
end

"""Largest finite entry of every column of `m`, `NaN` for a column without one."""
function _column_max(m::AbstractMatrix{<:Real})::Vector{Float64}
    out = fill(NaN, size(m, 2))
    for k in axes(m, 2), x in view(m, :, k)
        isfinite(x) && (isnan(out[k]) || x > out[k]) && (out[k] = x)
    end
    return out
end

"""
    plot_cluster_separation(snaps, cluster_ranges, cfg; coalescence = nothing, epochs = nothing,
                            settings = PostprocessConfig(), bound = (;), bound_sets = nothing,
                            filename = "merger_cluster_separation")

Plot the separation of the density centres of every pair of initial clusters
against time, with the half-mass threshold `min(r_h,i, r_h,j)` below which
the pair counts as merged ([`coalescence_time`](@ref)), over the merger
epochs ([`merger_epochs`](@ref)). `coalescence` and `epochs` take the
results of the pipeline (`results[:coalescence]`, `results[:merger_epochs]`);
when `epochs` is not given both are evaluated here on the bound sets
(`bound_sets`, or [`bound_set`](@ref) of every snapshot with the keywords
`bound`) with the `merger_*` keys of `settings`. Without `coalescence` the
number of distinct clusters is rebuilt from the separations of the epochs.

A pair's separation and threshold end at its core coalescence (`pair_core`):
from there the distance is the offset of two sub-centres inside one remnant
and continues as a thin line of the pair's colour, muted, against the
centre-noise floor of the two centres (dotted, the largest over the merged
pairs). The first contact `t_contact` and the core coalescence `t_core` are
marked with their times, and so is every join of the merger tree (the centre
coalescence, `t_coalesce` for two clusters) and, for three or more clusters,
the core coalescence of every join; the `interacting` and `coalescing`
regimes are shaded. For n_clusters ≤ 5 every pair is shown individually, its
threshold dashed in the same colour; for more clusters the min/max envelope
and mean of the separations are drawn with the mean threshold dashed, the
mean offset muted and the number of distinct clusters (`n_distinct`) on a
twin axis. Separations are shown in pc and times in Myr when
`cfg.units == "physical"` (header AS scaling).
"""
@publication function Nbody6Dynamics.plot_cluster_separation(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    filename::AbstractString = "merger_cluster_separation",
    coalescence::Union{Nothing,NamedTuple} = nothing,
    epochs::Union{Nothing,MergerEpochs} = nothing,
    settings::PostprocessConfig = PostprocessConfig(),
    bound::NamedTuple = (;),
    bound_sets::Union{Nothing,AbstractVector{BoundSet}} = nothing,
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 2 || (@info "Skipping separation plot: need ≥ 2 clusters (got $n_cl)"; return nothing)
    n_t = length(snaps)
    n_t ≥ 2 || (@info "Skipping separation plot: need ≥ 2 snapshots"; return nothing)
    bound_sets === nothing ||
        length(bound_sets) == n_t ||
        throw(DimensionMismatch("$(length(bound_sets)) bound sets for $n_t snapshots"))

    physical = cfg.units == "physical" && all(_has_physical_scaling(s.header) for s in snaps)
    t = physical ? [time_myr(s.header) for s in snaps] : [time_nb(s.header) for s in snaps]
    # The pipeline's centre coalescence and merger epochs, evaluated here with the
    # configured settings when the caller has none
    bsets = bound_sets
    co = coalescence
    ep = epochs
    if ep === nothing
        bsets === nothing && (bsets = BoundSet[bound_set(s; bound...) for s in snaps])
        co === nothing &&
            (co = coalescence_time(snaps, cluster_ranges; bound_sets = bsets, bound...))
        ep = merger_epochs(snaps, cluster_ranges, bsets, co; settings, bound...)
    end
    (length(ep.time) == n_t && size(ep.tag_regime, 1) == n_cl) || throw(
        DimensionMismatch(
            "merger epochs of $(size(ep.tag_regime, 1)) clusters and $(length(ep.time)) " *
            "snapshots for $n_cl clusters and $n_t snapshots",
        ),
    )
    n_distinct = co === nothing ? _distinct_clusters(ep) : co.n_distinct

    # Separations, thresholds, offsets after the core coalescence and the noise
    # floor, pairs × n_times; NB length → pc
    scale = physical ? [rbar(s.header) for s in snaps] : ones(n_t)
    n_pairs = length(ep.pairs)
    seps = ep.separation .* scale'
    thr = ep.threshold .* scale'
    retired = BitMatrix([0 < ep.pair_core[p] ≤ k for p in 1:n_pairs, k in 1:n_t])
    noise =
        any(retired) ? _pair_noise_floor(snaps, cluster_ranges, ep, bsets; bound) .* scale' :
        fill(NaN, n_pairs, n_t)
    floor_max = _column_max(noise)
    d_all = filter(isfinite, vcat(vec(seps), vec(thr), floor_max))

    # A second legend row (offsets and regimes) widens the canvas by its height,
    # so that the axis keeps the size of a single panel
    second_row = any(retired) || any(r -> r in (:interacting, :coalescing), ep.regime)
    w, h = _figsize_px(cfg)
    fig = Figure(; size = (w, second_row ? h + _LEGEND_ROW_HEIGHT : h))
    ax = Axis(
        fig[1, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = physical ? L"d_{ij} \; [\mathrm{pc}]" : L"d_{ij} \; [\mathrm{NB}]",
        xticks = _time_ticks(first(t), last(t)),
    )
    _time_limits!(ax, t)
    # Limits of the separation axis, between which the epoch labels are placed
    y_lo = 0.0
    y_top = NaN
    if isempty(d_all)
        _no_data_note!(ax, "No cluster with enough bound members")
    else
        y_lo, y_top = _pad_limits(0.0, maximum(d_all); floor = 0)
        ylims!(ax, y_lo, y_top)
        ax.yticks = _nice_ticks(y_lo, y_top)
    end
    t_epochs = physical ? ep.time_myr : ep.time
    bands = _regime_bands!((ax,), t_epochs, ep.regime; label = false)
    # Families of the legend, in two rows: the separations and their thresholds,
    # then what follows the core coalescence and the regimes
    family_titles = ["Separation of density centres:", "Half-mass threshold:"]
    legend_kw = (;
        orientation = :horizontal,
        framevisible = false,
        titleposition = :top,
        tellheight = true,
        padding = (0, 0, 0, 0),
    )
    offset_color = (:gray30, 0.45)
    floor_style = (; color = :gray40, linestyle = :dot, linewidth = _STYLE.guide)

    detail = n_cl ≤ 5
    # The polylines drawn, which the epoch labels keep clear of
    drawn = Tuple{Vector{Float64},Vector{Float64}}[]
    top_entries = Vector{_LegendElement}[]
    top_labels = Vector{AbstractString}[]
    if detail
        # One colour per pair: separation solid and its threshold dashed up to the
        # pair's core coalescence, the centre offset muted after it
        pair_entries = _LegendElement[]
        pair_labels = AbstractString[]
        for (p, (i, j)) in enumerate(ep.pairs)
            any(isfinite, @view seps[p, :]) || continue
            color = _OKABE_ITO[mod1(p, length(_OKABE_ITO))]
            k_core = ep.pair_core[p]
            pre = k_core > 0 ? (1:k_core) : (1:n_t)
            vs = [k for k in pre if isfinite(seps[p, k])]
            vt = [k for k in pre if isfinite(thr[p, k])]
            lines!(ax, t[vs], seps[p, vs]; color = color, linewidth = _STYLE.data)
            lines!(
                ax,
                t[vt],
                thr[p, vt];
                color = color,
                linestyle = :dash,
                linewidth = _STYLE.guide,
            )
            push!(drawn, (t[vs], seps[p, vs]), (t[vt], thr[p, vt]))
            if k_core > 0
                vo = [k for k in k_core:n_t if isfinite(seps[p, k])]
                lines!(ax, t[vo], seps[p, vo]; color = (color, 0.45), linewidth = _STYLE.fit)
                push!(drawn, (t[vo], seps[p, vo]))
            end
            push!(pair_entries, LineElement(; color = color, linewidth = _STYLE.data))
            push!(pair_labels, "$(i)–$(j)")
            # A single pair names its offset in its own colour, several in grey
            offset_color = length(pair_entries) == 1 ? (color, 0.45) : (:gray30, 0.45)
        end
        if !isempty(pair_entries)
            push!(top_entries, pair_entries)
            push!(top_labels, pair_labels)
        end
    else
        # Envelope over the pairs before their core coalescence: min/max band
        # plus mean line, mean threshold dashed; mean offset after it, muted
        pre = copy(seps)
        pre[retired] .= NaN
        thr_pre = copy(thr)
        thr_pre[retired] .= NaN
        post = fill(NaN, n_pairs, n_t)
        post[retired] .= seps[retired]
        d_min, d_max, d_mean = _envelope_stats(pre)
        _, _, thr_mean = _envelope_stats(thr_pre)
        _, _, off_mean = _envelope_stats(post)

        valid = .!isnan.(d_mean)
        sep_color = _SEMANTIC_COLORS[:separation]
        # Darker same-hue edges on the band fill
        sep_edge = _band_edge(sep_color)
        if any(valid)
            band!(ax, t[valid], d_min[valid], d_max[valid]; color = (sep_color, 0.25))
            lines!(ax, t[valid], d_min[valid]; color = sep_edge, linewidth = _STYLE.envelope)
            lines!(ax, t[valid], d_max[valid]; color = sep_edge, linewidth = _STYLE.envelope)
            lines!(ax, t[valid], d_mean[valid]; color = sep_color, linewidth = _STYLE.data)
            push!(drawn, (t[valid], d_min[valid]), (t[valid], d_max[valid]))
            vt = .!isnan.(thr_mean)
            push!(drawn, (t[vt], thr_mean[vt]))
            lines!(
                ax,
                t[vt],
                thr_mean[vt];
                color = sep_color,
                linestyle = :dash,
                linewidth = _STYLE.guide,
            )
            push!(
                top_entries,
                _LegendElement[
                    PolyElement(;
                        color = (sep_color, 0.25),
                        strokecolor = sep_edge,
                        strokewidth = _STYLE.envelope,
                    ),
                    LineElement(; color = sep_color, linewidth = _STYLE.data),
                ],
            )
            push!(top_labels, AbstractString["Min–max range", "Mean"])
        end
        vo = .!isnan.(off_mean)
        any(vo) &&
            lines!(ax, t[vo], off_mean[vo]; color = (sep_color, 0.45), linewidth = _STYLE.fit)
        push!(drawn, (t[vo], off_mean[vo]))
        offset_color = (sep_color, 0.45)

        # Twin axis: number of distinct clusters, the connected components of
        # the merged-pair graph of `coalescence_time`. Tick/label colour
        # matches the count series (semantic N colour).
        count_color = _SEMANTIC_COLORS[:n_particles]
        ax2 = Axis(
            fig[1, 1];
            yaxisposition = :right,
            ylabel = L"N_\mathrm{distinct}",
            ylabelcolor = count_color,
            yticklabelcolor = count_color,
            ytickcolor = count_color,
            # Overlay axis: grid off so it doesn't double-draw over ax
            ygridvisible = false,
            xgridvisible = false,
            xticklabelsvisible = false,
            xticksvisible = false,
            topspinevisible = false,
            bottomspinevisible = false,
            leftspinevisible = false,
            rightspinevisible = true,
            yticks = _integer_ticks(0, n_cl),
            limits = ((nothing, nothing), (0, n_cl + max(1, ceil(Int, 0.1 * n_cl)))),
        )
        hidespines!(ax2, :t, :b, :l)
        linkxaxes!(ax, ax2)
        stairs!(
            ax2,
            t,
            Float64.(n_distinct);
            color = count_color,
            linewidth = _STYLE.data,
            step = :post,
        )
    end
    vf = isfinite.(floor_max)
    any(vf) && lines!(ax, t[vf], floor_max[vf]; floor_style...)
    isfinite(y_top) && _epoch_markers!(
        ax,
        _separation_markers(ep, physical),
        (y_lo, y_top),
        (first(t), last(t));
        series = drawn,
    )

    if !isempty(top_entries)
        threshold_entry = _LegendElement[LineElement(;
            color = detail ? :black : _SEMANTIC_COLORS[:separation],
            linestyle = :dash,
            linewidth = _STYLE.guide,
        )]
        push!(top_entries, threshold_entry)
        push!(top_labels, AbstractString[detail ? L"\min(r_{h,i},\, r_{h,j})" : "Mean"])
        # The second row: the offsets after the core coalescence and the regimes
        low_entries = Vector{_LegendElement}[]
        low_labels = Vector{AbstractString}[]
        low_titles = String[]
        if any(retired)
            push!(
                low_entries,
                _LegendElement[
                    LineElement(; color = offset_color, linewidth = _STYLE.fit),
                    LineElement(; floor_style...),
                ],
            )
            push!(low_labels, AbstractString["Centre offset", "Noise floor"])
            push!(low_titles, "After core coalescence:")
        end
        if !isempty(bands)
            elements, names = _regime_legend(bands)
            push!(low_entries, elements)
            push!(low_labels, names)
            push!(low_titles, "Regime:")
        end
        legends = GridLayout(fig[0, 1])
        Legend(
            legends[1, 1],
            top_entries,
            top_labels,
            family_titles;
            nbanks = detail ? min(3, cld(length(top_labels[1]), 4)) : 1,
            legend_kw...,
        )
        if !isempty(low_entries)
            Legend(legends[2, 1], low_entries, low_labels, low_titles; nbanks = 1, legend_kw...)
            rowgap!(legends, 1, 8)
        end
    end

    return _save_fig(cfg, filename, fig)
end

"""
    plot_cluster_virial(snaps, cluster_ranges, cfg; epochs = nothing, bound = (;),
                        filename = "merger_cluster_virial")

Plot the internal virial ratio Q_i(t) of each initial cluster. For
n_clusters ≤ 5 every cluster is shown individually; for larger N, a min/max
envelope with the mean is drawn. A reference line at Q = 0.5 marks virial
equilibrium.  The time axis is in Myr when `cfg.units == "physical"`
(Q is dimensionless).

The metric tracks particles by their original cluster ID and, by default,
only those still bound to that group (`per_cluster_virial(...; bound_only =
true)`), so stripped stars and kicked stellar remnants do not enter. It holds
before the core coalescence (which clusters are still in virial equilibrium,
which undergo violent relaxation): after it the stars of a cluster are bound
to the remnant, not to each other. With `epochs` (the `MergerEpochs` of the
run) a cluster's curve therefore ends where [`per_cluster_virial`](@ref)
retires it (`coalescing`, `settled` or `dissolved`) and the `interacting` and
`coalescing` regimes of the system are shaded; without them the curves run
over the whole run. `bound` passes `max_passes` and `min_members` to
[`bound_set`](@ref).
"""
@publication function Nbody6Dynamics.plot_cluster_virial(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    filename::AbstractString = "merger_cluster_virial",
    bound::NamedTuple = (;),
    epochs::Union{Nothing,MergerEpochs} = nothing,
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping virial plot: no cluster ranges"; return nothing)
    n_t = length(snaps)
    n_t ≥ 2 || (@info "Skipping virial plot: need ≥ 2 snapshots"; return nothing)

    physical = cfg.units == "physical" && all(_has_physical_scaling(s.header) for s in snaps)
    t = physical ? [time_myr(s.header) for s in snaps] : [time_nb(s.header) for s in snaps]
    @info "Computing per-cluster virial ratios ($(n_cl) clusters × $(n_t) snapshots)..."
    Q, _ = per_cluster_virial(snaps, cluster_ranges; bound..., epochs)

    fig = Figure(; size = _figsize_px(cfg))

    # Pick a log scale when any cluster Q spikes high (mergers routinely do)
    q_all = filter(!isnan, vec(Q))
    use_log = !isempty(q_all) && maximum(q_all) > cfg.style.q_log_threshold

    # Clamp tiny/zero values on the log scale only
    Q_plot = use_log ? max.(Q, cfg.style.q_floor) : copy(Q)
    q_plotted = filter(!isnan, vec(Q_plot))

    ax = Axis(
        fig[1, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = L"Q_i = T_i / |W_i|",
        xticks = _time_ticks(first(t), last(t)),
        yscale = use_log ? log10 : identity,
        yticks = use_log && !isempty(q_plotted) ? _log_ticks(extrema(q_plotted)...) :
                 Makie.automatic,
    )
    _time_limits!(ax, t)

    detail = n_cl ≤ 5
    n_series = 0
    if detail
        for i in 1:n_cl
            valid = .!isnan.(@view Q[i, :])
            any(valid) || continue
            lines!(
                ax,
                t[valid],
                Q_plot[i, valid];
                color = _OKABE_ITO[mod1(i, length(_OKABE_ITO))],
                linewidth = _STYLE.data,
                label = latexstring("\\mathrm{Cluster}\\;$(i)"),
            )
            n_series += 1
        end
    else
        q_min, q_max, q_mean = _envelope_stats(Q_plot)
        valid = .!isnan.(q_mean)
        vir_color = _SEMANTIC_COLORS[:virial]
        band!(
            ax,
            t[valid],
            q_min[valid],
            q_max[valid];
            color = (vir_color, 0.25),
            label = L"\min\;-\;\max\;\mathrm{range}",
        )
        # Darker same-hue edges on the band fill
        vir_edge = _band_edge(vir_color)
        lines!(ax, t[valid], q_min[valid]; color = vir_edge, linewidth = _STYLE.envelope)
        lines!(ax, t[valid], q_max[valid]; color = vir_edge, linewidth = _STYLE.envelope)
        lines!(
            ax,
            t[valid],
            q_mean[valid];
            color = vir_color,
            linewidth = _STYLE.data,
            label = L"\mathrm{mean}",
        )
        n_series = 2
    end

    # Virial equilibrium reference
    hlines!(
        ax,
        [0.5];
        color = :gray50,
        linestyle = :dash,
        linewidth = _STYLE.guide,
        label = L"Q = 0.5\;\mathrm{(virial\;equilibrium)}",
    )
    # Regimes of the system, behind the series
    n_bands =
        epochs === nothing ? 0 :
        length(_regime_bands!((ax,), physical ? epochs.time_myr : epochs.time, epochs.regime))

    # ≥ 2 legend entries whenever at least one data series was drawn
    n_entries = n_series + 1 + n_bands
    n_series ≥ 1 && _top_legend!(fig, ax; nbanks = min(3, cld(n_entries, 4)))

    return _save_fig(cfg, filename, fig)
end

# -----------------------------------------------------------------------------
# Per-cluster structure: bound half-mass radius and bound mass fraction
# -----------------------------------------------------------------------------

"""
    plot_cluster_structure(snaps, cluster_ranges, cfg; lagr = nothing, epochs = nothing,
                           settings = PostprocessConfig(), bound = (;),
                           filename = "merger_cluster_structure")

Two stacked panels sharing the time axis: the half-mass radius of every
initial cluster measured about its own centre from its bound members
([`cluster_structure`](@ref) with the centre keys and `merger_min_members`
of `settings`), and the bound mass fraction. With `lagr` given, the engine's
global 50 % Lagrangian radius about its single density centre is overlaid
dashed grey on the radius panel, which makes the difference between the
configuration and its members visible. With `epochs` (the `MergerEpochs` of
the run) a cluster's curves end where `cluster_structure` retires it
(`coalescing`, `settled` or `dissolved`: its stars are bound to the remnant,
not to each other) and the `interacting` and `coalescing` regimes of the
system are shaded on both panels; without them the curves run over the whole
run. For more than five clusters the min–max envelope and mean are drawn.
Radii in pc and times in Myr when `cfg.units == "physical"`. `bound` passes
`max_passes` and `min_members` to [`bound_set`](@ref).
"""
@publication function Nbody6Dynamics.plot_cluster_structure(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    lagr::Union{Nothing,LagrangianData} = nothing,
    filename::AbstractString = "merger_cluster_structure",
    bound::NamedTuple = (;),
    epochs::Union{Nothing,MergerEpochs} = nothing,
    settings::PostprocessConfig = PostprocessConfig(),
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping cluster structure plot: no cluster ranges"; return nothing)
    n_t = length(snaps)
    n_t ≥ 2 || (@info "Skipping cluster structure plot: need ≥ 2 snapshots"; return nothing)

    physical = cfg.units == "physical" && all(_has_physical_scaling(s.header) for s in snaps)
    t = physical ? [time_myr(s.header) for s in snaps] : [time_nb(s.header) for s in snaps]
    r_unit = physical ? rbar(snaps[1].header) : 1.0
    @info "Computing per-cluster structure ($(n_cl) clusters × $(n_t) snapshots)..."
    st = cluster_structure(snaps, cluster_ranges; settings, bound..., epochs)
    r_h = st.r_lagr[2, :, :] .* r_unit
    f_bound = st.bound_mass_fraction

    fig = Figure(; size = _fig_two_panel(cfg))
    ttk = _time_ticks(first(t), last(t))
    ax1 = Axis(
        fig[1, 1];
        ylabel = physical ? L"r_{h,i} \; [\mathrm{pc}]" : L"r_{h,i} \; [\mathrm{NB}]",
        xticklabelsvisible = false,
        xticks = ttk,
    )
    ax2 = Axis(
        fig[2, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = L"M_{\mathrm{bound},i} \, / \, M_i",
        xticks = ttk,
    )

    detail = n_cl ≤ 5
    n_series = 0
    if detail
        for i in 1:n_cl
            valid = .!isnan.(@view r_h[i, :])
            any(valid) || continue
            color = _OKABE_ITO[mod1(i, length(_OKABE_ITO))]
            lines!(
                ax1,
                t[valid],
                r_h[i, valid];
                color = color,
                linewidth = _STYLE.data,
                label = latexstring("\\mathrm{Cluster}\\;$(i)"),
            )
            lines!(ax2, t[valid], f_bound[i, valid]; color = color, linewidth = _STYLE.data)
            n_series += 1
        end
    else
        sep_color = _SEMANTIC_COLORS[:separation]
        sep_edge = _band_edge(sep_color)
        for (ax, data, labelled) in ((ax1, r_h, true), (ax2, f_bound, false))
            lo, hi, mean = _envelope_stats(data)
            valid = .!isnan.(mean)
            any(valid) || continue
            if labelled
                band!(
                    ax,
                    t[valid],
                    lo[valid],
                    hi[valid];
                    color = (sep_color, 0.25),
                    label = L"\min\;-\;\max\;\mathrm{range}",
                )
                lines!(
                    ax,
                    t[valid],
                    mean[valid];
                    color = sep_color,
                    linewidth = _STYLE.data,
                    label = L"\mathrm{mean}",
                )
                n_series = 2
            else
                band!(ax, t[valid], lo[valid], hi[valid]; color = (sep_color, 0.25))
                lines!(ax, t[valid], mean[valid]; color = sep_color, linewidth = _STYLE.data)
            end
            lines!(ax, t[valid], lo[valid]; color = sep_edge, linewidth = _STYLE.envelope)
            lines!(ax, t[valid], hi[valid]; color = sep_edge, linewidth = _STYLE.envelope)
        end
    end

    # Radii of every series of the upper panel, for its limits
    r_panel = filter(isfinite, vec(r_h))

    # Engine's global 50 % radius about its single density centre
    if lagr !== nothing && !isempty(lagr.time)
        i50 = findfirst(==(0.5), lagr.mass_fractions)
        if i50 !== nothing
            t_l = physical ? lagr.time .* tscale(snaps[1].header) : lagr.time
            r_50 = lagr.radii[i50, :] .* r_unit
            append!(r_panel, filter(isfinite, r_50))
            lines!(
                ax1,
                t_l,
                r_50;
                color = :gray40,
                linestyle = :dash,
                linewidth = _STYLE.fit,
                label = L"\mathrm{Engine}\;r_{50}\;\mathrm{(global\;centre)}",
            )
            n_series += 1
        end
    end

    hlines!(ax2, [1.0]; color = :gray50, linestyle = :dash, linewidth = _STYLE.guide)
    # Occupancy of the bound-fraction panel, for the guide label and the limits
    f_x = Float64[]
    f_y = Float64[]
    for i in axes(f_bound, 1), k in axes(f_bound, 2)
        isfinite(f_bound[i, k]) || continue
        push!(f_x, t[k])
        push!(f_y, f_bound[i, k])
    end
    isempty(f_y) ||
        _guide_label!(ax2, L"M_{\mathrm{bound},i} = M_i", 1.0; xs = f_x, ys = f_y, color = :gray40)
    f_lo, f_hi = _pad_limits(isempty(f_y) ? 0.0 : minimum(f_y), 1.0; floor = 0)
    # Headroom above the unity guide for its label
    f_hi = 1.0 + 1.6 * (f_hi - 1.0)
    ylims!(ax2, f_lo, f_hi)
    # No tick above unity, which the fraction cannot exceed: the headroom holds
    # the guide label, and a tick there would meet the panel above at the joint
    ax2.yticks = filter(≤(1.0 + 1e-9), _nice_ticks(f_lo, f_hi))
    if !isempty(r_panel)
        r_lo, r_hi = _pad_limits(0.0, maximum(r_panel); floor = 0)
        ylims!(ax1, r_lo, r_hi)
        ax1.yticks = _nice_ticks(r_lo, r_hi)
    end
    # Regimes of the system on both panels, behind the series
    if epochs !== nothing
        t_ep = physical ? epochs.time_myr : epochs.time
        n_series += length(_regime_bands!((ax1, ax2), t_ep, epochs.regime))
    end
    _time_limits!(ax1, t)
    _time_limits!(ax2, t)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    n_series ≥ 2 && _top_legend!(fig, ax1; nbanks = min(3, cld(n_series, 4)))

    return _save_fig(cfg, filename, fig)
end

# -----------------------------------------------------------------------------
# Radial profiles: density against the generating model, velocity dispersions
# -----------------------------------------------------------------------------

"""
    _remnant_frame_tags(epochs, snap, n_cl) -> BitVector

Clusters retired by the merger epochs at `snap` (`coalescing`, `settled` or
`dissolved`), whose profile [`cluster_profiles`](@ref) with `retired =
:remnant` takes about the remnant centre; none without epochs.
"""
function _remnant_frame_tags(epochs, snap::Snapshot, n_cl::Int)::BitVector
    epochs === nothing && return falses(n_cl)
    k = Nbody6Dynamics._epoch_snapshot_index(epochs, snap)
    return Nbody6Dynamics._retired_tags(epochs, n_cl, length(epochs.time))[:, k]
end

"""Annotation naming the clusters drawn about the remnant centre; `nothing` when there are none."""
function _remnant_frame_note(tags::AbstractVector{Int})::Union{Nothing,String}
    isempty(tags) && return nothing
    one = length(tags) == 1
    return (one ? "Cluster $(only(tags))" : "Clusters " * join(tags, ", ")) *
           " inside the remnant:\n" *
           (one ? "profile" : "profiles") *
           " about its centre"
end

"""
    plot_density_profiles(snap, cluster_ranges, cfg; specs = nothing, epochs = nothing,
                          settings = PostprocessConfig(), bound = (;),
                          filename = "merger_density_profiles")

Density profiles `ρ(r)` of every initial cluster about its own centre from
its bound members ([`cluster_profiles`](@ref) with the centre, profile and
`merger_min_members` keys of `settings`), log–log, with the generating model
overlaid dashed when `specs` (the clusters' [`ClusterSpec`](@ref)s) are
given: each model is evaluated at the cluster's initial half-mass radius and
its bound mass at this snapshot. A ratio strip `ρ / ρ_model` with a labelled
unity guide sits beneath the main panel. For more than five clusters the
min–max envelope and mean of the ratio are drawn. Radii in pc and densities
in M☉ pc⁻³ when `cfg.units == "physical"`. `bound` passes `max_passes` and
`min_members` to [`bound_set`](@ref).

With `epochs` (the `MergerEpochs` of the run, `snap` among its snapshots) only
the clusters still `detached`, `interacting` or a `satellite` at `snap` are
drawn about their own centre against their model; a cluster retired there
(`coalescing`, `settled` or `dissolved`: its stars are bound to the remnant,
not to each other) is drawn from its members in the bound set of the whole
snapshot about the remnant centre, without a model, and an annotation names
it.
"""
@publication function Nbody6Dynamics.plot_density_profiles(
    snap::Snapshot,
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    specs::Union{Nothing,AbstractVector} = nothing,   # Vector{ClusterSpec}; the type is defined later in the module
    filename::AbstractString = "merger_density_profiles",
    bound::NamedTuple = (;),
    epochs::Union{Nothing,MergerEpochs} = nothing,
    settings::PostprocessConfig = PostprocessConfig(),
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping density profiles: no cluster ranges"; return nothing)
    specs === nothing ||
        length(specs) == n_cl ||
        throw(ArgumentError("specs must have one entry per cluster range"))
    h = snap.header
    physical = cfg.units == "physical" && _has_physical_scaling(h)
    r_unit = physical ? rbar(h) : 1.0
    ρ_unit = physical ? zmbar(h) / rbar(h)^3 : 1.0
    profiles =
        cluster_profiles(snap, cluster_ranges; settings, bound..., epochs, retired = :remnant)
    any(!isnothing, profiles) ||
        (@info "Skipping density profiles: no cluster with enough members"; return nothing)
    in_remnant = _remnant_frame_tags(epochs, snap, n_cl)
    # The generating model describes a cluster about its own centre only
    with_model = specs !== nothing && any(i -> profiles[i] !== nothing && !in_remnant[i], 1:n_cl)

    # The ratio panel is an auxiliary strip, not a second main panel.
    fig = Figure(; size = with_model ? _fig_with_strip(cfg) : _figsize_px(cfg))
    ax1 = Axis(
        fig[1, 1];
        xlabel = with_model ? "" : (physical ? L"r \; [\mathrm{pc}]" : L"r \; [\mathrm{NB}]"),
        ylabel = physical ? L"\rho \; [\mathrm{M}_\odot\,\mathrm{pc}^{-3}]" :
                 L"\rho \; [\mathrm{NB}]",
        xscale = log10,
        yscale = log10,
        xticklabelsvisible = !with_model,
    )
    ax2 =
        with_model ?
        Axis(
            fig[2, 1];
            xlabel = physical ? L"r \; [\mathrm{pc}]" : L"r \; [\mathrm{NB}]",
            ylabel = L"\rho \, / \, \rho_{\mathrm{model}}",
            xscale = log10,
            yscale = log10,
        ) : nothing
    with_model && _strip_rows!(fig, 1, 2)

    detail = n_cl ≤ 5
    r_all = Float64[]
    ρ_all = Float64[]
    ratios = Vector{Vector{Float64}}()
    r_ratio = Vector{Vector{Float64}}()
    n_series = 0
    for i in 1:n_cl
        prof = profiles[i]
        prof === nothing && continue
        valid = .!isnan.(prof.rho) .& (prof.rho .> 0)
        any(valid) || continue
        r = prof.r[valid] .* r_unit
        ρ = prof.rho[valid] .* ρ_unit
        append!(r_all, r)
        append!(ρ_all, ρ)
        color = _OKABE_ITO[mod1(i, length(_OKABE_ITO))]
        if detail
            lines!(
                ax1,
                r,
                ρ;
                color = color,
                linewidth = _STYLE.data,
                label = latexstring("\\mathrm{Cluster}\\;$(i)"),
            )
            n_series += 1
        else
            lines!(ax1, r, ρ; color = (:gray50, 0.5), linewidth = _STYLE.ghost)
        end
        if with_model && !in_remnant[i]
            spec = specs[i]
            r_h_model = spec.rbar / rbar(h)              # initial half-mass radius [NB]
            ρ_model = model_density(spec.profile, prof.M, r_h_model)
            ρm = [ρ_model(x) for x in prof.r[valid]] .* ρ_unit
            ok = ρm .> 0
            if detail && any(ok)
                lines!(ax1, r[ok], ρm[ok]; color = color, linestyle = :dash, linewidth = _STYLE.fit)
                lines!(ax2, r[ok], ρ[ok] ./ ρm[ok]; color = color, linewidth = _STYLE.data)
            end
            push!(ratios, ρ[ok] ./ ρm[ok])
            push!(r_ratio, r[ok])
        end
    end
    if !detail && with_model && !isempty(ratios)
        # Envelope on a common radial grid: interpolate each ratio onto the union grid
        grid = exp10.(
            range(log10(minimum(minimum, r_ratio)), log10(maximum(maximum, r_ratio)); length = 24),
        )
        mat = fill(NaN, length(ratios), length(grid))
        for (i, (rr, q)) in enumerate(zip(r_ratio, ratios))
            for (k, g) in enumerate(grid)
                j = searchsortedlast(rr, g)
                (j < 1 || j ≥ length(rr)) && continue
                f = (g - rr[j]) / (rr[j + 1] - rr[j])
                mat[i, k] = (1 - f) * q[j] + f * q[j + 1]
            end
        end
        lo, hi, mean = _envelope_stats(mat)
        valid = .!isnan.(mean)
        c = _SEMANTIC_COLORS[:separation]
        band!(
            ax2,
            grid[valid],
            lo[valid],
            hi[valid];
            color = (c, 0.25),
            label = L"\min\;-\;\max\;\mathrm{range}",
        )
        lines!(
            ax2,
            grid[valid],
            mean[valid];
            color = c,
            linewidth = _STYLE.data,
            label = L"\mathrm{mean}",
        )
        lines!(ax1, r_all[1:1], ρ_all[1:1]; color = :gray50, label = L"\mathrm{Clusters}")
        n_series = 2
    end
    if with_model && detail
        lines!(
            ax1,
            r_all[1:1],
            ρ_all[1:1];
            color = :black,
            linestyle = :dash,
            linewidth = _STYLE.fit,
            label = L"\mathrm{Generating\;model}",
        )
        n_series += 1
    end
    if with_model
        hlines!(ax2, [1.0]; color = :gray50, linestyle = :dash, linewidth = _STYLE.guide)
        q_x = reduce(vcat, r_ratio; init = Float64[])
        q_y = reduce(vcat, ratios; init = Float64[])
        isempty(q_y) || _guide_label!(
            ax2,
            L"\rho = \rho_{\mathrm{model}}",
            1.0;
            xs = q_x,
            ys = q_y,
            color = :gray40,
            log = true,
            xlog = true,
        )
        # Symmetric about the model in log, holding the data: a factor of two
        # at least, so that agreement is not magnified into disagreement.
        q_fin = [q for q in q_y if isfinite(q) && q > 0]
        m = (isempty(q_fin) ? log10(2) : max(maximum(abs.(log10.(q_fin))), log10(2))) + 0.1
        ylims!(ax2, 10.0^-m, 10.0^m)
        ax2.yticks = _log_ticks(10.0^-m, 10.0^m)
        linkxaxes!(ax1, ax2)
        rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    end
    if !isempty(ρ_all)
        # One tick style for every epoch of a run, whatever range it spans.
        ax1.yticks = _log_ticks(minimum(ρ_all), maximum(ρ_all); plain = false)
        xt = _log_ticks(minimum(r_all), maximum(r_all))
        ax1.xticks = xt
        ax2 === nothing || (ax2.xticks = xt)
    end
    note = _remnant_frame_note([i for i in 1:n_cl if in_remnant[i] && profiles[i] !== nothing])
    if note !== nothing && !isempty(ρ_all)
        _annotate!(
            ax1,
            note;
            corner = _emptiest_corner(log10.(r_all), log10.(ρ_all); corners = (:tr, :br, :tl)),
            color = :gray25,
        )
    end
    n_series ≥ 2 && _top_legend!(fig, ax1; nbanks = min(3, cld(n_series, 4)))
    return _save_fig(cfg, filename, fig)
end

"""
    plot_velocity_dispersion(snap, cluster_ranges, cfg; epochs = nothing,
                             settings = PostprocessConfig(), bound = (;),
                             filename = "merger_velocity_dispersion")

Radial (solid) and one-dimensional tangential (dashed) velocity dispersion
profiles of every initial cluster from its bound members
([`cluster_profiles`](@ref) with the keys of `settings`), with the anisotropy
`β(r) = 1 − σ_t²/σ_r²` beneath and the isotropic `β = 0` guide labelled.
Velocities in km s⁻¹ and radii in pc when `cfg.units == "physical"`. For
more than five clusters the min–max envelope and mean of `σ_r` and `β` are
drawn. `bound` passes `max_passes` and `min_members` to [`bound_set`](@ref).
With `epochs` (the `MergerEpochs` of the run) a cluster retired at `snap`
(`coalescing`, `settled` or `dissolved`) is drawn from its members in the
bound set of the whole snapshot about the remnant centre, as in
[`plot_density_profiles`](@ref), and an annotation names it.
"""
@publication function Nbody6Dynamics.plot_velocity_dispersion(
    snap::Snapshot,
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    cfg::VisualizationConfig;
    filename::AbstractString = "merger_velocity_dispersion",
    bound::NamedTuple = (;),
    epochs::Union{Nothing,MergerEpochs} = nothing,
    settings::PostprocessConfig = PostprocessConfig(),
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping velocity dispersion: no cluster ranges"; return nothing)
    h = snap.header
    physical = cfg.units == "physical" && _has_physical_scaling(h)
    r_unit = physical ? rbar(h) : 1.0
    v_unit = physical ? vstar(h) : 1.0
    profiles =
        cluster_profiles(snap, cluster_ranges; settings, bound..., epochs, retired = :remnant)
    any(!isnothing, profiles) ||
        (@info "Skipping velocity dispersion: no cluster with enough members"; return nothing)
    in_remnant = _remnant_frame_tags(epochs, snap, n_cl)

    fig = Figure(; size = _fig_two_panel(cfg))
    ax1 = Axis(
        fig[1, 1];
        ylabel = physical ? L"\sigma \; [\mathrm{km\,s^{-1}}]" : L"\sigma \; [\mathrm{NB}]",
        xscale = log10,
        xticklabelsvisible = false,
    )
    ax2 = Axis(
        fig[2, 1];
        xlabel = physical ? L"r \; [\mathrm{pc}]" : L"r \; [\mathrm{NB}]",
        ylabel = L"\beta = 1 - \sigma_t^2 / \sigma_r^2",
        xscale = log10,
    )
    detail = n_cl ≤ 5
    r_all = Float64[]
    n_series = 0
    rows_sr = Vector{Vector{Float64}}()
    rows_beta = Vector{Vector{Float64}}()
    rows_r = Vector{Vector{Float64}}()
    # Occupancy of both panels, for their limits and the anisotropy guide label
    sig_all = Float64[]
    beta_x = Float64[]
    beta_y = Float64[]
    for i in 1:n_cl
        prof = profiles[i]
        prof === nothing && continue
        valid = .!isnan.(prof.sigma_r)
        any(valid) || continue
        r = prof.r[valid] .* r_unit
        append!(r_all, r)
        color = _OKABE_ITO[mod1(i, length(_OKABE_ITO))]
        if detail
            lines!(
                ax1,
                r,
                prof.sigma_r[valid] .* v_unit;
                color = color,
                linewidth = _STYLE.data,
                label = latexstring("\\mathrm{Cluster}\\;$(i)"),
            )
            lines!(
                ax1,
                r,
                prof.sigma_t[valid] .* v_unit;
                color = color,
                linestyle = :dash,
                linewidth = _STYLE.fit,
            )
            lines!(ax2, r, prof.beta[valid]; color = color, linewidth = _STYLE.data)
            append!(sig_all, filter(isfinite, prof.sigma_r[valid] .* v_unit))
            append!(sig_all, filter(isfinite, prof.sigma_t[valid] .* v_unit))
            for (x, b) in zip(r, prof.beta[valid])
                isfinite(b) || continue
                push!(beta_x, x)
                push!(beta_y, b)
            end
            n_series += 1
        else
            push!(rows_sr, prof.sigma_r[valid] .* v_unit)
            push!(rows_beta, prof.beta[valid])
            push!(rows_r, r)
        end
    end
    if !detail && !isempty(rows_r)
        grid = exp10.(
            range(log10(minimum(minimum, rows_r)), log10(maximum(maximum, rows_r)); length = 24),
        )
        function onto(rows_y)
            mat = fill(NaN, length(rows_y), length(grid))
            for (i, (rr, q)) in enumerate(zip(rows_r, rows_y))
                for (k, g) in enumerate(grid)
                    j = searchsortedlast(rr, g)
                    (j < 1 || j ≥ length(rr)) && continue
                    f = (g - rr[j]) / (rr[j + 1] - rr[j])
                    mat[i, k] = (1 - f) * q[j] + f * q[j + 1]
                end
            end
            return mat
        end
        c = _SEMANTIC_COLORS[:virial]
        for (ax, mat, labelled) in ((ax1, onto(rows_sr), true), (ax2, onto(rows_beta), false))
            lo, hi, mean = _envelope_stats(mat)
            valid = .!isnan.(mean)
            any(valid) || continue
            if labelled
                band!(
                    ax,
                    grid[valid],
                    lo[valid],
                    hi[valid];
                    color = (c, 0.25),
                    label = L"\sigma_r\;\min\;-\;\max",
                )
                lines!(
                    ax,
                    grid[valid],
                    mean[valid];
                    color = c,
                    linewidth = _STYLE.data,
                    label = L"\sigma_r\;\mathrm{mean}",
                )
            else
                band!(ax, grid[valid], lo[valid], hi[valid]; color = (c, 0.25))
                lines!(ax, grid[valid], mean[valid]; color = c, linewidth = _STYLE.data)
            end
            if ax === ax1
                append!(sig_all, filter(isfinite, vcat(lo[valid], hi[valid])))
            else
                for (x, b) in zip(vcat(grid[valid], grid[valid]), vcat(lo[valid], hi[valid]))
                    isfinite(b) || continue
                    push!(beta_x, x)
                    push!(beta_y, b)
                end
            end
        end
        n_series = 2
    end
    if detail && n_series ≥ 1
        # Colour names the cluster, line style the component of the dispersion.
        lines!(
            ax1,
            r_all[1:1],
            [NaN];
            color = :black,
            linewidth = _STYLE.data,
            label = L"\sigma_r\;\mathrm{(solid)}",
        )
        lines!(
            ax1,
            r_all[1:1],
            [NaN];
            color = :black,
            linestyle = :dash,
            linewidth = _STYLE.fit,
            label = L"\sigma_t\;\mathrm{(dashed)}",
        )
        n_series += 2
    end
    hlines!(ax2, [0.0]; color = :gray50, linestyle = :dash, linewidth = _STYLE.guide)
    isempty(beta_y) || _guide_label!(
        ax2,
        L"\beta = 0\;\mathrm{(isotropic)}",
        0.0;
        xs = beta_x,
        ys = beta_y,
        color = :gray40,
        xlog = true,
    )
    isempty(sig_all) || ylims!(ax1, _pad_limits(0.0, maximum(sig_all); floor = 0)...)
    if isempty(beta_y)
        # No shell holds the members the anisotropy needs: say so, and keep the
        # tick labels of the empty panel clear of the joint
        _no_data_note!(
            ax2,
            latexstring(
                "Fewer than $(settings.profile_n_min_beta) members per shell: no \$\\beta\$",
            ),
        )
        ylims!(ax2, _pad_limits(-1.0, 1.0)...)
    else
        ylims!(ax2, _pad_limits(min(0.0, minimum(beta_y)), max(0.0, maximum(beta_y)))...)
    end
    if !isempty(r_all)
        xt = _log_ticks(minimum(r_all), maximum(r_all))
        ax1.xticks = xt
        ax2.xticks = xt
    end
    note = _remnant_frame_note([i for i in 1:n_cl if in_remnant[i] && profiles[i] !== nothing])
    if note !== nothing
        # Occupancy of the dispersion panel, in log r as the axis draws it
        occ_x = Float64[]
        occ_y = Float64[]
        for prof in profiles
            prof === nothing && continue
            for (x, a, b) in zip(prof.r, prof.sigma_r, prof.sigma_t), s in (a, b)
                (isfinite(s) && x > 0) || continue
                push!(occ_x, log10(x * r_unit))
                push!(occ_y, s * v_unit)
            end
        end
        isempty(occ_x) || _annotate!(
            ax1,
            note;
            corner = _emptiest_corner(
                occ_x,
                occ_y;
                corners = (:tr, :br, :tl),
                ylims = isempty(sig_all) ? nothing : _pad_limits(0.0, maximum(sig_all); floor = 0),
            ),
            color = :gray25,
        )
    end
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    n_series ≥ 2 && _top_legend!(fig, ax1; nbanks = min(3, cld(n_series, 4)))
    return _save_fig(cfg, filename, fig)
end
