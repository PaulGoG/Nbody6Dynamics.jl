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
const _BOUND_MAX_ITER = 8

"""
    _member_indices(snap, rng) -> Vector{Int}

Indices in `snap` of the particles whose original name lies in `rng`.
"""
function _member_indices(snap::Snapshot, rng::UnitRange{Int})::Vector{Int}
    return [i for (i, n) in enumerate(snap.name) if first(rng) ≤ n ≤ last(rng)]
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
    _self_potential(pos, mass) -> Vector{Float64}

Potential of every particle from all the others, `Φ_i = −Σ_{j≠i} m_j / r_ij`
(`G = 1`), by the O(N²) pair sum.
"""
function _self_potential(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
)::Vector{Float64}
    n = length(mass)
    Φ = zeros(Float64, n)
    @inbounds for i in 1:(n - 1)
        xi, yi, zi, mi = pos[1, i], pos[2, i], pos[3, i], mass[i]
        for j in (i + 1):n
            dx = xi - pos[1, j]
            dy = yi - pos[2, j]
            dz = zi - pos[3, j]
            r = sqrt(dx * dx + dy * dy + dz * dz)
            r > 0 || continue
            Φ[i] -= mass[j] / r
            Φ[j] -= mi / r
        end
    end
    return Φ
end

"""
    _bound_members(pos, vel, mass) -> Vector{Int}

Indices of the particles bound to the group they form: iterate the
self-consistent selection `½|v − v_com|² + Φ < 0`, with the centre-of-mass
velocity and the potential recomputed from the current bound set, until the
set no longer changes (at most `_BOUND_MAX_ITER` passes). O(N²) per pass.
"""
function _bound_members(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
)::Vector{Int}
    bound = collect(1:length(mass))
    for _ in 1:_BOUND_MAX_ITER
        length(bound) < _MIN_MEMBERS && return bound
        p = pos[:, bound]
        v = vel[:, bound]
        m = mass[bound]
        M = sum(m)
        vc = vec(sum(v .* m'; dims = 2)) ./ M
        Φ = _self_potential(p, m)
        keep = [
            0.5 * ((v[1, i] - vc[1])^2 + (v[2, i] - vc[2])^2 + (v[3, i] - vc[3])^2) + Φ[i] < 0
            for i in eachindex(m)
        ]
        all(keep) && return bound
        bound = bound[keep]
    end
    return bound
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
    bound_fraction(snap::Snapshot) -> Float64

Mass fraction of the particles bound to the whole system in its own frame
(self-consistent negative-energy selection, [`_bound_members`](@ref)):
a snapshot-based escaper measure independent of the engine's escape
sphere. O(N²) per pass; intended for N ≲ 3×10⁴.
"""
function bound_fraction(snap::Snapshot)::Float64
    pos = Float64.(snap.pos)
    vel = Float64.(snap.vel)
    mass = Float64.(snap.mass)
    isempty(mass) && return NaN
    sel = _bound_members(pos, vel, mass)
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

Complexity O(Σ_i N_i²) per snapshot for the bound selection and the virial
ratio; intended for N_i ≲ 10⁴.
"""
function cluster_structure(
    snaps::Vector{Snapshot},
    cluster_ranges::Vector{UnitRange{Int}};
    bound_only::Bool = true,
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
        vel_all = Float64.(snap.vel)
        mass_all = Float64.(snap.mass)
        for (i, rng) in enumerate(cluster_ranges)
            idx = _member_indices(snap, rng)
            n_members[i, k] = length(idx)
            length(idx) < _MIN_MEMBERS && continue
            pos = pos_all[:, idx]
            vel = vel_all[:, idx]
            mass = mass_all[idx]
            sel = bound_only ? _bound_members(pos, vel, mass) : collect(1:length(idx))
            n_bound[i, k] = length(sel)
            fbound[i, k] = sum(mass[sel]) / sum(mass)
            length(sel) < _MIN_MEMBERS && continue
            p = pos[:, sel]
            v = vel[:, sel]
            m = mass[sel]
            M = sum(m)
            c = _shrinking_sphere_centre(p, m)
            centre[:, i, k] = c
            r_lagr[:, i, k] = _lagrangian_radii(p, m, c, _LAGR_FRACTIONS_MEMBER)
            vc = vec(sum(v .* m'; dims = 2)) ./ M
            T = 0.0
            for j in eachindex(m)
                T += 0.5 * m[j] * ((v[1, j] - vc[1])^2 + (v[2, j] - vc[2])^2 + (v[3, j] - vc[3])^2)
            end
            σ[i, k] = sqrt(2 * T / (3 * M))
            W = 0.5 * sum(m .* _self_potential(p, m))
            q[i, k] = W < 0 ? T / abs(W) : NaN
        end
    end
    return ClusterStructure(time, n_members, n_bound, fbound, centre, r_lagr, σ, q)
end
