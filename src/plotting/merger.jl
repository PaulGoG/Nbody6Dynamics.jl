# =============================================================================
# Merger-specific diagnostics
# =============================================================================
# Track each initial cluster's centre-of-mass over the course of the run.
# Gives meaningful pre-coalescence diagnostics for multi-cluster setups where
# global Lagrangian radii measure cluster separation rather than internal
# structure.

"""Coalescence heuristic: max/mean pairwise-separation ratio below which the ensemble is considered merged."""
const _COALESCENCE_RATIO_MAX = 2.0

"""Coalescence heuristic: mean separation must also drop below this multiple of the initial minimum separation."""
const _COALESCENCE_MEAN_SEP_FACTOR = 1.5

"""
    parse_merger_summary(path::AbstractString) -> Vector{UnitRange{Int}}

Read `merger_summary.txt` and return a vector of index ranges, one per initial
cluster, using the post-truncation particle counts. Ranges refer to the
particle IDs in `dat.10` (1-based, contiguous across clusters).

Returns an empty vector if the file cannot be parsed.
"""
function parse_merger_summary(path::AbstractString)::Vector{UnitRange{Int}}
    isfile(path) || return UnitRange{Int}[]

    trunc_counts = Int[]
    # Match e.g. "Cluster 5: plummer, imf=kroupa, N=1000 (after trunc: 966),"
    # Tolerant of extra comma-separated fields between the profile name and
    # the N=… count (the summary format has grown fields before; the
    # write→parse round-trip test in runtests.jl guards this coupling).
    pattern = r"Cluster\s+\d+:\s+.*?\bN=\d+\s+\(after trunc:\s+(\d+)\)"
    for line in eachline(path)
        m = match(pattern, line)
        m === nothing && continue
        push!(trunc_counts, parse(Int, m.captures[1]))
    end

    isempty(trunc_counts) && return UnitRange{Int}[]

    ranges = Vector{UnitRange{Int}}(undef, length(trunc_counts))
    offset = 0
    for (i, n) in enumerate(trunc_counts)
        ranges[i] = (offset + 1):(offset + n)
        offset += n
    end
    return ranges
end

"""
    _cluster_com_trajectories(snaps, cluster_ranges)
        -> (coms, present, r_rms)

Compute the mass-weighted centre of mass of each initial cluster at every
snapshot. `coms` is `3 × n_clusters × n_times`. `present[i, k]` is false when
fewer than 3 members of cluster `i` remain in snapshot `k`. `r_rms[i, k]` is
the mass-weighted RMS radius of cluster `i`'s members from their COM (a
proxy for cluster extent); `NaN` when `present` is false.
"""
function _cluster_com_trajectories(
    snaps::Vector{Snapshot},
    cluster_ranges::Vector{UnitRange{Int}},
)
    n_cl = length(cluster_ranges)
    n_t  = length(snaps)
    coms    = fill(NaN, 3, n_cl, n_t)
    r_rms   = fill(NaN, n_cl, n_t)
    present = falses(n_cl, n_t)

    for (k, snap) in enumerate(snaps)
        names_k = Int.(snap.name)
        for (i, rng) in enumerate(cluster_ranges)
            mask = [n in rng for n in names_k]
            n_mem = count(mask)
            n_mem < 3 && continue
            m = Float64.(snap.mass[mask])
            M = sum(m)
            M > 0 || continue
            cx = sum(m .* Float64.(snap.pos[1, mask])) / M
            cy = sum(m .* Float64.(snap.pos[2, mask])) / M
            cz = sum(m .* Float64.(snap.pos[3, mask])) / M
            coms[1, i, k] = cx
            coms[2, i, k] = cy
            coms[3, i, k] = cz

            # Mass-weighted RMS radius: sqrt(Σ m_k r_k² / Σ m_k)
            x = Float64.(@view snap.pos[1, mask])
            y = Float64.(@view snap.pos[2, mask])
            z = Float64.(@view snap.pos[3, mask])
            ssum = 0.0
            @inbounds for j in eachindex(m)
                dx = x[j] - cx
                dy = y[j] - cy
                dz = z[j] - cz
                ssum += m[j] * (dx*dx + dy*dy + dz*dz)
            end
            r_rms[i, k] = sqrt(ssum / M)
            present[i, k] = true
        end
    end
    return coms, present, r_rms
end

"""
    _count_spatial_clusters(coms, r_rms, present; overlap_factor = 1.0) -> Vector{Int}

At each snapshot, count spatially distinct initial clusters via union-find.
Two clusters are considered merged when their COMs are closer than
`overlap_factor × (r_ref_i + r_ref_j)`, where `r_ref` is the *initial* RMS
radius of each cluster (held fixed throughout). Using the initial radius
(rather than current) prevents spurious "mergers" caused by infall-driven
puffing — the criterion then measures true spatial convergence of COMs.
"""
function _count_spatial_clusters(
    coms::Array{Float64,3},
    r_rms::AbstractMatrix{<:Real},
    present::AbstractMatrix{Bool};
    overlap_factor::Real = 1.0,
)
    n_cl = size(coms, 2)
    n_t  = size(coms, 3)
    counts = zeros(Int, n_t)

    # Reference radii: each cluster's RMS radius in its first present snapshot
    r_ref = fill(NaN, n_cl)
    for i in 1:n_cl
        for k in 1:n_t
            if present[i, k] && !isnan(r_rms[i, k])
                r_ref[i] = r_rms[i, k]
                break
            end
        end
    end

    # Simple union-find with path compression
    parent = zeros(Int, n_cl)
    function find(p::Vector{Int}, a::Int)
        while p[a] != a
            p[a] = p[p[a]]
            a = p[a]
        end
        return a
    end

    for k in 1:n_t
        active = findall(@view present[:, k])
        isempty(active) && (counts[k] = 0; continue)

        for i in 1:n_cl; parent[i] = i; end
        for ia in eachindex(active)
            i = active[ia]
            isnan(r_ref[i]) && continue
            for jb in ia+1:length(active)
                j = active[jb]
                isnan(r_ref[j]) && continue
                dx = coms[1, i, k] - coms[1, j, k]
                dy = coms[2, i, k] - coms[2, j, k]
                dz = coms[3, i, k] - coms[3, j, k]
                d  = sqrt(dx*dx + dy*dy + dz*dz)
                if d < overlap_factor * (r_ref[i] + r_ref[j])
                    ri = find(parent, i)
                    rj = find(parent, j)
                    ri == rj || (parent[ri] = rj)
                end
            end
        end
        roots = Set{Int}()
        for i in active; push!(roots, find(parent, i)); end
        counts[k] = length(roots)
    end
    return counts
end

"""
    plot_cluster_separation(snaps, cluster_ranges, cfg;
                            filename = "merger_cluster_separation")

Plot the evolution of pairwise centre-of-mass separations between the initial
clusters. For small n_clusters (≤ 5) every pair is shown individually; for
larger N, min/max/mean envelopes are plotted instead.  Separations are shown
in pc and times in Myr when `cfg.units == "physical"` (header AS scaling).

The time at which the ratio `max_sep / mean_sep` drops below a heuristic
threshold is marked as an estimated "coalescence time". After that epoch the
initial cluster decomposition is no longer physically meaningful.
"""
function plot_cluster_separation(
    snaps::Vector{Snapshot},
    cluster_ranges::Vector{UnitRange{Int}},
    cfg::VisualizationConfig;
    filename::AbstractString = "merger_cluster_separation",
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 2 || (@info "Skipping separation plot: need ≥ 2 clusters (got $n_cl)"; return nothing)
    n_t = length(snaps)
    n_t ≥ 2 || (@info "Skipping separation plot: need ≥ 2 snapshots"; return nothing)

    physical = cfg.units == "physical" &&
               all(_has_physical_scaling(s.header) for s in snaps)
    t = physical ? [time_myr(s.header) for s in snaps] : [time_nb(s.header) for s in snaps]
    coms, present, r_rms = _cluster_com_trajectories(snaps, cluster_ranges)

    # Pairwise separations over time: pairs × n_times
    pairs = [(i, j) for i in 1:n_cl-1 for j in i+1:n_cl]
    n_pairs = length(pairs)
    seps = fill(NaN, n_pairs, n_t)
    for (p, (i, j)) in enumerate(pairs)
        for k in 1:n_t
            (present[i, k] && present[j, k]) || continue
            dx = coms[1, i, k] - coms[1, j, k]
            dy = coms[2, i, k] - coms[2, j, k]
            dz = coms[3, i, k] - coms[3, j, k]
            seps[p, k] = sqrt(dx*dx + dy*dy + dz*dz)
        end
    end
    # NB length → pc (uniform scaling; ratio-based coalescence heuristics
    # below are unaffected)
    if physical
        for k in 1:n_t
            seps[:, k] .*= rbar(snaps[k].header)
        end
    end

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(fig[1, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = physical ? L"d_{ij} \; [\mathrm{pc}]" : L"d_{ij} \; [\mathrm{NB}]",
        xticks = _time_ticks(first(t), last(t)),
    )

    detail = n_cl ≤ 5
    if detail
        # One line per pair, labelled
        n_lines = 0
        for (p, (i, j)) in enumerate(pairs)
            valid = .!isnan.(seps[p, :])
            any(valid) || continue
            lines!(ax, t[valid], seps[p, valid];
                color = _OKABE_ITO[mod1(p, length(_OKABE_ITO))],
                linewidth = 1.8,
                label = "$(i)–$(j)")
            n_lines += 1
        end
        n_lines ≥ 2 && _top_legend!(fig, ax; nbanks = min(3, cld(n_pairs, 4)))
    else
        # Envelope: show min/max band plus mean line
        d_min  = fill(NaN, n_t)
        d_max  = fill(NaN, n_t)
        d_mean = fill(NaN, n_t)
        for k in 1:n_t
            vals = Float64[]
            for p in 1:n_pairs
                x = seps[p, k]
                isnan(x) || push!(vals, x)
            end
            isempty(vals) && continue
            d_min[k]  = minimum(vals)
            d_max[k]  = maximum(vals)
            d_mean[k] = sum(vals) / length(vals)
        end

        valid = .!isnan.(d_mean)
        sep_color = _SEMANTIC_COLORS[:separation]
        band_plot = band!(ax, t[valid], d_min[valid], d_max[valid];
                          color = (sep_color, 0.25))
        # Darker same-hue edges on the band fill
        sep_edge = _band_edge(sep_color)
        lines!(ax, t[valid], d_min[valid]; color = sep_edge, linewidth = 1.0)
        lines!(ax, t[valid], d_max[valid]; color = sep_edge, linewidth = 1.0)
        mean_plot = lines!(ax, t[valid], d_mean[valid];
                           color = sep_color, linewidth = 2.4)

        # Mark estimated coalescence: max/mean ratio and mean separation both
        # drop below their heuristic thresholds (consts at top of file)
        d_init_min = minimum(filter(!isnan, d_min))
        ratio = [isnan(d_mean[k]) ? NaN : d_max[k] / d_mean[k] for k in 1:n_t]
        idx_merge = findfirst(k -> !isnan(ratio[k]) && ratio[k] < _COALESCENCE_RATIO_MAX &&
                                    d_mean[k] < _COALESCENCE_MEAN_SEP_FACTOR * d_init_min, 1:n_t)

        merge_plot = nothing
        merge_label = ""
        if idx_merge !== nothing && t[idx_merge] > first(t)
            t_merge = t[idx_merge]
            merge_plot = vlines!(ax, [t_merge]; color = :black,
                                 linestyle = :dash, linewidth = 1.5)
            merge_label = latexstring("t_\\mathrm{merge} \\approx $(round(t_merge; digits=2))")
        end

        # Twin axis: count of spatially-distinct clusters. Two clusters are
        # merged when |COM_i − COM_j| < (r_rms_i + r_rms_j), i.e. their
        # member spheres overlap. Connected components via union-find.
        # Tick/label colour matches the count series (semantic N colour).
        survivor_count = _count_spatial_clusters(coms, r_rms, present;
                                                  overlap_factor = 1.0)
        count_color = _SEMANTIC_COLORS[:n_particles]
        ax2 = Axis(fig[1, 1];
            yaxisposition = :right,
            ylabel = L"N_\mathrm{clusters}\;\mathrm{(\geq 3\;members)}",
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
        surv_plot = stairs!(ax2, t, Float64.(survivor_count);
                            color = count_color, linewidth = 2.2,
                            step = :post)

        # Manual legend combining primary and twin axes: horizontal, above
        legend_elems = Any[band_plot, mean_plot]
        legend_labels = [L"\min\;-\;\max\;\mathrm{range}", L"\mathrm{mean}"]
        if merge_plot !== nothing
            push!(legend_elems, merge_plot)
            push!(legend_labels, merge_label)
        end
        push!(legend_elems, surv_plot)
        push!(legend_labels, L"\mathrm{surviving\;clusters}")
        _top_legend!(fig, legend_elems, legend_labels)
    end

    return _save_fig(cfg, filename, fig)
end

# -----------------------------------------------------------------------------
# Per-cluster virial ratio (COM-subtracted, self-gravity only)
# -----------------------------------------------------------------------------

"""
    _cluster_virial_snapshot(snap, rng) -> (Q, n_mem)

Compute the virial ratio Q = T/|W| for the members of an initial cluster in a
single snapshot, after subtracting the cluster centre-of-mass velocity. Only
self-gravity between cluster members contributes to W. Returns `(NaN, n_mem)`
if the cluster has < 3 members or |W| is not positive.

N-body units are assumed (G = 1).
"""
function _cluster_virial_snapshot(snap::Snapshot, rng::UnitRange{Int})
    names_k = Int.(snap.name)
    mask = [n in rng for n in names_k]
    n_mem = count(mask)
    n_mem < 3 && return (NaN, n_mem)

    m = Float64.(snap.mass[mask])
    M = sum(m)
    M > 0 || return (NaN, n_mem)

    # Positions and velocities for cluster members
    x = Float64.(@view snap.pos[1, mask])
    y = Float64.(@view snap.pos[2, mask])
    z = Float64.(@view snap.pos[3, mask])
    vx = Float64.(@view snap.vel[1, mask])
    vy = Float64.(@view snap.vel[2, mask])
    vz = Float64.(@view snap.vel[3, mask])

    # COM velocity (subtract bulk motion so Q measures internal kinetic energy)
    vcx = sum(m .* vx) / M
    vcy = sum(m .* vy) / M
    vcz = sum(m .* vz) / M

    # Internal kinetic energy
    T = 0.0
    @inbounds for i in eachindex(m)
        dvx = vx[i] - vcx
        dvy = vy[i] - vcy
        dvz = vz[i] - vcz
        T += 0.5 * m[i] * (dvx*dvx + dvy*dvy + dvz*dvz)
    end

    # Self-gravitational potential energy of the cluster (G = 1)
    W = 0.0
    @inbounds for i in 1:n_mem - 1
        xi, yi, zi, mi = x[i], y[i], z[i], m[i]
        for j in i+1:n_mem
            dx = xi - x[j]
            dy = yi - y[j]
            dz = zi - z[j]
            r = sqrt(dx*dx + dy*dy + dz*dz)
            r > 0 || continue
            W -= mi * m[j] / r
        end
    end

    absW = abs(W)
    absW > 0 || return (NaN, n_mem)
    return (T / absW, n_mem)
end

"""
    per_cluster_virial(snaps, cluster_ranges) -> (Q::Matrix{Float64}, n_mem::Matrix{Int})

Compute the internal virial ratio of every initial cluster at every snapshot.
`Q[i, k]` is the ratio for cluster `i` at snapshot `k`; `NaN` when the cluster
has < 3 members left.

Complexity: O(Σ_i N_i²) per snapshot. Manageable for N_i ≲ few × 10³.
"""
function per_cluster_virial(
    snaps::Vector{Snapshot},
    cluster_ranges::Vector{UnitRange{Int}},
)
    n_cl = length(cluster_ranges)
    n_t  = length(snaps)
    Q     = fill(NaN, n_cl, n_t)
    n_mem = zeros(Int, n_cl, n_t)

    for (k, snap) in enumerate(snaps)
        for (i, rng) in enumerate(cluster_ranges)
            q, nm = _cluster_virial_snapshot(snap, rng)
            Q[i, k] = q
            n_mem[i, k] = nm
        end
    end
    return Q, n_mem
end

"""
    plot_cluster_virial(snaps, cluster_ranges, cfg;
                        filename = "merger_cluster_virial")

Plot the internal virial ratio Q_i(t) of each initial cluster. For
n_clusters ≤ 5 every cluster is shown individually; for larger N, a min/max
envelope with the mean is drawn. A reference line at Q = 0.5 marks virial
equilibrium.  The time axis is in Myr when `cfg.units == "physical"`
(Q is dimensionless).

**Interpretation note:** the metric tracks particles by their original
cluster ID, not spatial membership. It is physically meaningful *before*
coalescence (diagnosing which clusters are still in virial equilibrium vs
undergoing violent relaxation). After merger, the same ID-group is spread
throughout the merged system, so T/|W| for its self-gravity only diverges —
this is expected, not a bug, and provides a rough proxy for "time of
coalescence" as the point where Q_i saturates at a large constant.
"""
function plot_cluster_virial(
    snaps::Vector{Snapshot},
    cluster_ranges::Vector{UnitRange{Int}},
    cfg::VisualizationConfig;
    filename::AbstractString = "merger_cluster_virial",
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || (@info "Skipping virial plot: no cluster ranges"; return nothing)
    n_t = length(snaps)
    n_t ≥ 2 || (@info "Skipping virial plot: need ≥ 2 snapshots"; return nothing)

    physical = cfg.units == "physical" &&
               all(_has_physical_scaling(s.header) for s in snaps)
    t = physical ? [time_myr(s.header) for s in snaps] : [time_nb(s.header) for s in snaps]
    @info "Computing per-cluster virial ratios ($(n_cl) clusters × $(n_t) snapshots)..."
    Q, _ = per_cluster_virial(snaps, cluster_ranges)

    fig = Figure(; size = _figsize_px(cfg))

    # Pick a log scale when any cluster Q spikes high (mergers routinely do)
    q_all = filter(!isnan, vec(Q))
    use_log = !isempty(q_all) && maximum(q_all) > cfg.style.q_log_threshold

    # Clamp tiny/zero values on the log scale only
    Q_plot = use_log ? max.(Q, cfg.style.q_floor) : copy(Q)
    q_plotted = filter(!isnan, vec(Q_plot))

    ax = Axis(fig[1, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = L"Q_i = T_i / |W_i|",
        xticks = _time_ticks(first(t), last(t)),
        yscale = use_log ? log10 : identity,
        yticks = use_log && !isempty(q_plotted) ?
                 _log_ticks(extrema(q_plotted)...) : Makie.automatic,
    )

    detail = n_cl ≤ 5
    n_series = 0
    if detail
        for i in 1:n_cl
            valid = .!isnan.(@view Q[i, :])
            any(valid) || continue
            lines!(ax, t[valid], Q_plot[i, valid];
                color = _OKABE_ITO[mod1(i, length(_OKABE_ITO))],
                linewidth = 1.8,
                label = latexstring("\\mathrm{cluster}\\;$(i)"))
            n_series += 1
        end
    else
        q_min = fill(NaN, n_t)
        q_max = fill(NaN, n_t)
        q_mean = fill(NaN, n_t)
        for k in 1:n_t
            vals = Float64[]
            for i in 1:n_cl
                x = Q_plot[i, k]
                isnan(x) || push!(vals, x)
            end
            isempty(vals) && continue
            q_min[k]  = minimum(vals)
            q_max[k]  = maximum(vals)
            q_mean[k] = sum(vals) / length(vals)
        end
        valid = .!isnan.(q_mean)
        vir_color = _SEMANTIC_COLORS[:virial]
        band!(ax, t[valid], q_min[valid], q_max[valid];
              color = (vir_color, 0.25),
              label = L"\min\;-\;\max\;\mathrm{range}")
        # Darker same-hue edges on the band fill
        vir_edge = _band_edge(vir_color)
        lines!(ax, t[valid], q_min[valid]; color = vir_edge, linewidth = 1.0)
        lines!(ax, t[valid], q_max[valid]; color = vir_edge, linewidth = 1.0)
        lines!(ax, t[valid], q_mean[valid];
               color = vir_color, linewidth = 2.4,
               label = L"\mathrm{mean}")
        n_series = 2
    end

    # Virial equilibrium reference
    hlines!(ax, [0.5];
            color = :gray50, linestyle = :dash, linewidth = 1.0,
            label = L"Q = 0.5\;\mathrm{(virial\;equilibrium)}")

    # ≥ 2 legend entries whenever at least one data series was drawn
    n_series ≥ 1 &&
        _top_legend!(fig, ax; nbanks = detail ? min(3, cld(n_cl + 1, 4)) : 1)

    return _save_fig(cfg, filename, fig)
end
