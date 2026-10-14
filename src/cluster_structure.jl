# =============================================================================
# Per-cluster structure from snapshots by original membership
# =============================================================================
# The engine keeps one density centre and one scale radius, so before
# coalescence its lagr.7, RC, and RTIDE describe the configuration, not the
# member clusters. These routines measure every initial cluster about its own
# centre, from the members it still binds.

"""
    ClusterStructure

Per-cluster structural quantities at every snapshot, all by original
membership (`cluster_ranges` of `merger_ic.toml`) and in N-body units.
Arrays are `n_clusters × n_times`; `NaN` marks clusters with fewer than
`_MIN_MEMBERS` bound members.

# Fields
- `time`: snapshot times [NB]
- `n_members`: members present in the snapshot
- `n_bound`: members bound to the cluster (negative energy in its own frame)
- `bound_mass_fraction`: bound mass over present member mass
- `centre`: `3 × n_clusters × n_times` shrinking-sphere centre of the bound members
- `r_lagr`: `3 × n_clusters × n_times` radii enclosing 10 %, 50 %, 90 % of the bound mass
- `sigma_1d`: one-dimensional velocity dispersion of the bound members [NB]
- `q_virial`: `T/|W|` of the bound members (COM velocity subtracted)
"""
struct ClusterStructure
    time::Vector{Float64}
    n_members::Matrix{Int}
    n_bound::Matrix{Int}
    bound_mass_fraction::Matrix{Float64}
    centre::Array{Float64,3}
    r_lagr::Array{Float64,3}
    sigma_1d::Matrix{Float64}
    q_virial::Matrix{Float64}
end

const _LAGR_FRACTIONS_MEMBER = (0.1, 0.5, 0.9)
const _MIN_MEMBERS = 10

"""
    _member_indices(snap, members) -> Vector{Int}

Indices in `snap` of the particles whose original name (body index in
`dat.10`) belongs to `members`, a contiguous range or, with primordial
binaries, a vector of body indices.
"""
function _member_indices(snap::Snapshot, members::UnitRange{Int})::Vector{Int}
    return [i for (i, n) in enumerate(snap.name) if first(members) ≤ n ≤ last(members)]
end
function _member_indices(snap::Snapshot, members::AbstractVector{Int})::Vector{Int}
    isempty(members) && return Int[]
    flag = falses(maximum(members))
    flag[members] .= true
    return [i for (i, n) in enumerate(snap.name) if 1 ≤ n ≤ length(flag) && flag[n]]
end

"""
    _shrinking_sphere_centre(pos, mass; keep = 0.9, n_stop = 20) -> Vector{Float64}

Mass-weighted centre by the shrinking-sphere iteration (Power et al. 2003):
starting from the centre of mass, keep the fraction `keep` of the particles
closest to the current centre and recentre, until `n_stop` particles (or
1 % of the input, whichever is larger) remain. O(N log N) per iteration;
insensitive to unbound outliers and to a second cluster's tidal debris.
"""
function _shrinking_sphere_centre(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    keep::Float64 = 0.9,
    n_stop::Int = 20,
)::Vector{Float64}
    idx = collect(1:length(mass))
    n_min = max(n_stop, ceil(Int, 0.01 * length(mass)))
    c = zeros(3)
    while true
        M = sum(@view mass[idx])
        for k in 1:3
            c[k] = sum(mass[i] * pos[k, i] for i in idx) / M
        end
        length(idx) ≤ n_min && break
        d2 = [(pos[1, i] - c[1])^2 + (pos[2, i] - c[2])^2 + (pos[3, i] - c[3])^2 for i in idx]
        n_keep = max(n_min, round(Int, keep * length(idx)))
        idx = idx[partialsortperm(d2, 1:n_keep)]
    end
    return c
end

"""
    _virial_energies(snap, idx) -> (T, W, M)

Kinetic energy about the mass-weighted mean velocity, potential energy and
mass of the systems formed by the snapshot stars `idx` (regularised pairs at
their centres of mass, components of pairs split by `idx` left out).
"""
function _virial_energies(snap::Snapshot, idx::AbstractVector{Int})
    sys, _ = _snapshot_systems(snap, collect(idx))
    isempty(sys.mass) && return (NaN, NaN, 0.0)
    vref = _mean_velocity(sys, eachindex(sys.mass))
    Φ = pair_potential(sys.pos, sys.mass)
    T = 0.0
    W = 0.0
    for s in eachindex(sys.mass)
        T +=
            0.5 *
            sys.mass[s] *
            (
                (sys.vel[1, s] - vref[1])^2 +
                (sys.vel[2, s] - vref[2])^2 +
                (sys.vel[3, s] - vref[3])^2
            )
        W += 0.5 * sys.mass[s] * Φ[s]
    end
    return (T, W, sum(sys.mass))
end

"""
    _lagrangian_radii(pos, mass, centre, fractions) -> Vector{Float64}

Radii about `centre` enclosing the given fractions of the total mass.
"""
function _lagrangian_radii(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    centre::AbstractVector{Float64},
    fractions,
)::Vector{Float64}
    r = [
        sqrt((pos[1, i] - centre[1])^2 + (pos[2, i] - centre[2])^2 + (pos[3, i] - centre[3])^2)
        for i in eachindex(mass)
    ]
    order = sortperm(r)
    cum = cumsum(mass[order])
    M = cum[end]
    return [r[order[something(findfirst(≥(f * M), cum), length(cum))]] for f in fractions]
end

"""
    RadialProfile

Radial structure of one particle set about a centre, in log-spaced shells
(the innermost shell is the sphere enclosing `n_min_inner` particles):

- `r`: mass-weighted mean radius of each shell
- `r_edges`: shell boundaries (`length(r) + 1`)
- `rho`: shell density `ΔM / (4π/3 (r_{k+1}³ − r_k³))`
- `sigma_r`, `sigma_t`: radial and one-dimensional tangential velocity
  dispersions of the shell about the set's centre-of-mass velocity
  (`σ_t² = ½ ⟨|v_t − ⟨v_t⟩|²⟩`); `NaN` in shells with fewer than 5 members
- `beta`: anisotropy `1 − σ_t²/σ_r²` (0 isotropic, 1 radial)
- `n`: members per shell
- `M`, `r_h`: total mass and half-mass radius of the set
"""
struct RadialProfile
    r::Vector{Float64}
    r_edges::Vector{Float64}
    rho::Vector{Float64}
    sigma_r::Vector{Float64}
    sigma_t::Vector{Float64}
    beta::Vector{Float64}
    n::Vector{Int}
    M::Float64
    r_h::Float64
end

"""
    radial_profile(pos, vel, mass, centre; nbins = 12, n_min_inner = 20,
                   mass_fraction_max = 0.99) -> RadialProfile

Shell profile of the particles about `centre` (N-body or any consistent
units): `nbins` shells whose outer edges are log-spaced between the radius
enclosing `n_min_inner` particles and the radius enclosing
`mass_fraction_max` of the mass.
"""
function radial_profile(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    centre::AbstractVector{Float64};
    nbins::Int = 12,
    n_min_inner::Int = 20,
    mass_fraction_max::Float64 = 0.99,
)::RadialProfile
    N = length(mass)
    N ≥ 2 * n_min_inner ||
        throw(ArgumentError("radial_profile: need at least $(2 * n_min_inner) particles; got $N"))
    nbins ≥ 2 || throw(ArgumentError("radial_profile: nbins must be ≥ 2"))
    d = [
        sqrt((pos[1, i] - centre[1])^2 + (pos[2, i] - centre[2])^2 + (pos[3, i] - centre[3])^2)
        for i in 1:N
    ]
    order = sortperm(d)
    cum = cumsum(mass[order])
    M = cum[end]
    r_h = d[order[findfirst(≥(0.5 * M), cum)]]
    r_min = max(d[order[n_min_inner]], eps())
    r_max = d[order[findfirst(≥(mass_fraction_max * M), cum)]]
    r_max > r_min || (r_max = r_min * 10)
    edges = vcat(0.0, exp10.(range(log10(r_min), log10(r_max); length = nbins)))
    vc = vec(sum(vel .* mass'; dims = 2)) ./ M
    r_mean = fill(NaN, nbins)
    rho = fill(NaN, nbins)
    σr = fill(NaN, nbins)
    σt = fill(NaN, nbins)
    β = fill(NaN, nbins)
    n = zeros(Int, nbins)
    for k in 1:nbins
        idx = [i for i in 1:N if edges[k] ≤ d[i] < edges[k + 1]]
        n[k] = length(idx)
        isempty(idx) && continue
        m = mass[idx]
        ΔM = sum(m)
        rho[k] = ΔM / (4π / 3 * (edges[k + 1]^3 - edges[k]^3))
        r_mean[k] = sum(m .* d[idx]) / ΔM
        n[k] < 5 && continue
        vr = zeros(n[k])
        vt = zeros(3, n[k])
        for (j, i) in enumerate(idx)
            rx = pos[1, i] - centre[1]
            ry = pos[2, i] - centre[2]
            rz = pos[3, i] - centre[3]
            rn = max(d[i], eps())
            ux, uy, uz = rx / rn, ry / rn, rz / rn
            dvx, dvy, dvz = vel[1, i] - vc[1], vel[2, i] - vc[2], vel[3, i] - vc[3]
            vr[j] = dvx * ux + dvy * uy + dvz * uz
            vt[1, j] = dvx - vr[j] * ux
            vt[2, j] = dvy - vr[j] * uy
            vt[3, j] = dvz - vr[j] * uz
        end
        w = m ./ ΔM
        mean_vr = sum(w .* vr)
        σr2 = sum(w .* (vr .- mean_vr) .^ 2)
        mean_vt = vec(sum(vt .* w'; dims = 2))
        σt2 = 0.5 * sum(w[j] * sum((vt[:, j] .- mean_vt) .^ 2) for j in 1:n[k])
        σr[k] = sqrt(σr2)
        σt[k] = sqrt(σt2)
        β[k] = σr2 > 0 ? 1 - σt2 / σr2 : NaN
    end
    return RadialProfile(r_mean, edges, rho, σr, σt, β, n, M, r_h)
end

"""
    cluster_profiles(snap, cluster_ranges; bound_only = true, nbins = 12)
        -> Vector{Union{Nothing,RadialProfile}}

Radial profile of every initial cluster about its own shrinking-sphere
centre from its bound members ([`cluster_structure`](@ref) conventions;
`max_passes` and `min_members` go to [`bound_set`](@ref)); `nothing` for
clusters with too few members.
"""
function cluster_profiles(
    snap::Snapshot,
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    bound_only::Bool = true,
    nbins::Int = 12,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)
    pos_all = Float64.(snap.pos)
    vel_all = Float64.(snap.vel)
    mass_all = Float64.(snap.mass)
    out = Vector{Union{Nothing,RadialProfile}}(nothing, length(cluster_ranges))
    for (i, rng) in enumerate(cluster_ranges)
        idx = _member_indices(snap, rng)
        length(idx) < 2 * _MIN_MEMBERS && continue
        sel = bound_only ? bound_set(snap, idx; max_passes, min_members).stars : idx
        length(sel) < 2 * _MIN_MEMBERS && continue
        p = pos_all[:, sel]
        c = _shrinking_sphere_centre(p, mass_all[sel])
        out[i] = radial_profile(p, vel_all[:, sel], mass_all[sel], c; nbins = nbins)
    end
    return out
end

"""
    system_profile(snap; nbins = 12) -> RadialProfile

Radial profile of the whole snapshot about its shrinking-sphere centre
(the remnant after coalescence).
"""
function system_profile(snap::Snapshot; nbins::Int = 12)::RadialProfile
    pos = Float64.(snap.pos)
    vel = Float64.(snap.vel)
    mass = Float64.(snap.mass)
    c = _shrinking_sphere_centre(pos, mass)
    return radial_profile(pos, vel, mass, c; nbins = nbins)
end

"""
    bound_fraction(snap::Snapshot; max_passes = 100, min_members = 10) -> Float64

Mass fraction of the stars bound to the whole system in its own frame
(iterative energetic unbinding with regularised pairs at their centres of
mass, [`bound_set`](@ref)): a snapshot-based escaper measure independent of
the engine's escape sphere. One O(N²) pair sum plus O(N·N_removed) per pass.
"""
function bound_fraction(
    snap::Snapshot;
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)::Float64
    mass = Float64.(snap.mass)
    isempty(mass) && return NaN
    sel = bound_set(snap; max_passes, min_members).stars
    return sum(mass[sel]) / sum(mass)
end

"""
    cluster_structure(snaps, cluster_ranges; bound_only = true) -> ClusterStructure

Structure of every initial cluster at every snapshot, measured about the
cluster's own shrinking-sphere centre from the members it binds
(`bound_only = false` uses all present members and reproduces the plain
membership statistics). The bound selection removes tidally stripped
stars and kicked stellar remnants, which otherwise dominate the velocity
dispersion and the virial ratio after the first supernovae.

Regularised pairs enter the bound selection, `sigma_1d` and `q_virial` as
their centres of mass; `max_passes` and `min_members` go to
[`bound_set`](@ref). One O(N_i²) pair sum per cluster and snapshot.
"""
function cluster_structure(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    bound_only::Bool = true,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)::ClusterStructure
    n_cl = length(cluster_ranges)
    n_t = length(snaps)
    time = [time_nb(s.header) for s in snaps]
    n_members = zeros(Int, n_cl, n_t)
    n_bound = zeros(Int, n_cl, n_t)
    fbound = fill(NaN, n_cl, n_t)
    centre = fill(NaN, 3, n_cl, n_t)
    r_lagr = fill(NaN, 3, n_cl, n_t)
    σ = fill(NaN, n_cl, n_t)
    q = fill(NaN, n_cl, n_t)
    for (k, snap) in enumerate(snaps)
        pos_all = Float64.(snap.pos)
        mass_all = Float64.(snap.mass)
        for (i, rng) in enumerate(cluster_ranges)
            idx = _member_indices(snap, rng)
            n_members[i, k] = length(idx)
            length(idx) < _MIN_MEMBERS && continue
            if bound_only
                b = bound_set(snap, idx; max_passes, min_members)
                sel = b.stars
                T, W = b.kinetic, b.potential
            else
                sel = idx
                T, W, _ = _virial_energies(snap, idx)
            end
            n_bound[i, k] = length(sel)
            fbound[i, k] = sum(mass_all[sel]) / sum(mass_all[idx])
            length(sel) < _MIN_MEMBERS && continue
            p = pos_all[:, sel]
            m = mass_all[sel]
            M = sum(m)
            c = _shrinking_sphere_centre(p, m)
            centre[:, i, k] = c
            r_lagr[:, i, k] = _lagrangian_radii(p, m, c, _LAGR_FRACTIONS_MEMBER)
            σ[i, k] = sqrt(2 * T / (3 * M))
            q[i, k] = W < 0 ? T / abs(W) : NaN
        end
    end
    return ClusterStructure(time, n_members, n_bound, fbound, centre, r_lagr, σ, q)
end

# ---------------------------------------------------------------------------
# Original-membership diagnostics of a merger run
# ---------------------------------------------------------------------------

"""
    parse_merger_summary(path::AbstractString) -> Vector{Vector{Int}}

Read `merger_summary.txt` and return the body indices (`dat.10` order,
1-based) of every initial cluster from the post-truncation body counts
and, when present, the pair counts: pairs of all clusters come first,
cluster by cluster, then the singles, so a cluster with binaries owns two
contiguous blocks. Returns an empty vector if the file cannot be parsed.
"""
function parse_merger_summary(path::AbstractString)::Vector{Vector{Int}}
    isfile(path) || return Vector{Int}[]

    counts = Tuple{Int,Int}[]   # (bodies after truncation, pairs)
    # Match e.g. "Cluster 5: plummer, imf=kroupa, N=1000 (after trunc: 966, binaries: 12),"
    # Tolerant of extra comma-separated fields between the profile name and
    # the N=… count (the summary format has grown fields before; the
    # write→parse round-trip test in runtests.jl guards this coupling).
    pattern = r"Cluster\s+\d+:\s+.*?\bN=\d+\s+\(after trunc:\s+(\d+)(?:,\s*binaries:\s+(\d+))?\)"
    for line in eachline(path)
        m = match(pattern, line)
        m === nothing && continue
        n_b = m.captures[2] === nothing ? 0 : parse(Int, m.captures[2])
        push!(counts, (parse(Int, m.captures[1]), n_b))
    end
    isempty(counts) && return Vector{Int}[]

    nbin0 = sum(last, counts)
    members = Vector{Vector{Int}}(undef, length(counts))
    pair_offset = 0
    single_offset = 2nbin0
    for (i, (n, n_b)) in enumerate(counts)
        n_s = n - 2n_b
        blocks = [(pair_offset + 1):(pair_offset + 2n_b), (single_offset + 1):(single_offset + n_s)]
        members[i] = _members_from_blocks(blocks)
        pair_offset += 2n_b
        single_offset += n_s
    end
    return members
end

"""
    _cluster_virial_snapshot(snap, rng; bound_only = true, max_passes, min_members) -> (Q, n_mem)

Virial ratio `Q = T/|W|` of the members of an initial cluster in one
snapshot, about the members' mass-weighted mean velocity and with only the
self-gravity among the selected members in `W`; regularised pairs count as
their centres of mass. With `bound_only` the selection is the bound subset
([`bound_set`](@ref)), which removes tidally stripped stars and kicked stellar
remnants; otherwise every present member counts. `n_mem` is the number of
present members. Returns `(NaN, n_mem)` with fewer than 3 selected members or
a non-negative `W`. N-body units (`G = 1`).
"""
function _cluster_virial_snapshot(
    snap::Snapshot,
    rng::AbstractVector{Int};
    bound_only::Bool = true,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)
    idx = _member_indices(snap, rng)
    n_mem = length(idx)
    n_mem < 3 && return (NaN, n_mem)
    if bound_only
        b = bound_set(snap, idx; max_passes, min_members)
        length(b.stars) < 3 && return (NaN, n_mem)
        T, W = b.kinetic, b.potential
    else
        T, W, _ = _virial_energies(snap, idx)
    end
    W < 0 || return (NaN, n_mem)
    return (T / abs(W), n_mem)
end

"""
    per_cluster_virial(snaps, cluster_ranges; bound_only = true)
        -> (Q::Matrix{Float64}, n_mem::Matrix{Int})

Internal virial ratio of every initial cluster at every snapshot, `Q[i, k]`
for cluster `i` at snapshot `k` (`NaN` with fewer than 3 selected members).
By default only the members bound to the cluster enter `T` and `W`
([`bound_set`](@ref), which takes `max_passes` and `min_members`);
`bound_only = false` uses every present member and is dominated by kicked
remnants after the first supernovae. One O(N_i²) pair sum per cluster and
snapshot.
"""
function per_cluster_virial(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    bound_only::Bool = true,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)
    n_cl = length(cluster_ranges)
    n_t = length(snaps)
    Q = fill(NaN, n_cl, n_t)
    n_mem = zeros(Int, n_cl, n_t)

    for (k, snap) in enumerate(snaps)
        for (i, rng) in enumerate(cluster_ranges)
            q, nm = _cluster_virial_snapshot(snap, rng; bound_only, max_passes, min_members)
            Q[i, k] = q
            n_mem[i, k] = nm
        end
    end
    return Q, n_mem
end
