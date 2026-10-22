# =============================================================================
# k-d tree neighbour utilities: local density and snapshot binary finder
# =============================================================================
#
# Both utilities query a k-d tree built on the stellar positions, so the
# neighbour searches cost O(n log n) instead of the O(n²) of a direct scan.
# The local density is the Casertano & Hut (1985) estimator
# (DOI 10.1086/163485). The binary finder joins bound neighbours greedily in
# order of their specific orbital energy, one level of node joining after Bate
# (2009a, DOI 10.1111/j.1365-2966.2008.14106.x) ordered by energy instead of
# separation, and classifies the pairs by the hardness convention of Heggie
# (1975) (DOI 10.1093/mnras/173.3.729).

using NearestNeighbors: KDTree, knn

"""
    _self_slot(idx, self) -> Int

Position in the neighbour list `idx` of a tree query at star `self` that is
dropped to leave the star's neighbours: the entry holding `self` itself, or
the last entry when coincident points pushed `self` out of the list.
"""
function _self_slot(idx::AbstractVector{<:Integer}, self::Int)
    @inbounds for p in eachindex(idx)
        idx[p] == self && return p
    end
    return lastindex(idx)
end

# ---------------------------------------------------------------------------
# Local density (Casertano & Hut 1985)
# ---------------------------------------------------------------------------

"""
    local_density(pos, mass; k = 6, weighting = :mass) -> Vector{Float64}

Casertano & Hut (1985, DOI 10.1086/163485) density estimate at every
particle of `pos` (3 × n) with masses `mass`, from its `k` nearest
neighbours found on a k-d tree in O(n log n). With `r_k` the distance to the
`k`-th neighbour,

- `weighting = :mass`: ``ρ_i = \\sum_{j=1}^{k-1} m_j / (4π r_k^3 / 3)``, the
  mass of the `k − 1` nearest neighbours over the sphere reaching the `k`-th;
  the mass density behind the density centre and core radius
  ([`core_radius`](@ref)) of a cluster with a mass spectrum;
- `weighting = :number`: ``n_i = (k - 1) / (4π r_k^3 / 3)``, the number
  density, for equal-weight tracers or when the mass spectrum must not bias
  the centre.

The particle itself is excluded by index, so coincident points count as
neighbours of one another; a particle whose `k`-th neighbour coincides with
it gets `0`. `k = 6` is the rank recommended by Casertano & Hut and used by
the engine.

Throws `ArgumentError` unless `k ≥ 2`, `n > k` and `weighting` is `:mass` or
`:number`, and `DimensionMismatch` when `pos` is not 3 × `length(mass)`.

# Example
```julia
ρ = local_density(Float64.(snap.pos), Float64.(snap.mass))
core_radius(Float64.(snap.pos), ρ)
```
"""
function local_density(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    k::Int = 6,
    weighting::Symbol = :mass,
)
    n = length(mass)
    size(pos) == (3, n) ||
        throw(DimensionMismatch("local_density: pos is $(size(pos)), expected (3, $n)"))
    k ≥ 2 || throw(ArgumentError("local_density: k must be at least 2, got $k"))
    n > k || throw(ArgumentError("local_density needs more than k = $k particles, got $n"))
    weighting in (:mass, :number) ||
        throw(ArgumentError("local_density: weighting must be :mass or :number, got :$weighting"))
    tree = KDTree(pos)
    idxs, dists = knn(tree, pos, k + 1, true)
    return _density_from_neighbours(idxs, dists, mass, k, weighting === :mass)
end

"""
    _density_from_neighbours(idxs, dists, mass, k, by_mass) -> Vector{Float64}

Densities from sorted `k + 1`-neighbour query results `idxs`, `dists` at
every particle; see [`local_density`](@ref).
"""
function _density_from_neighbours(
    idxs::AbstractVector{<:AbstractVector{<:Integer}},
    dists::AbstractVector{<:AbstractVector{<:Real}},
    mass::AbstractVector{Float64},
    k::Int,
    by_mass::Bool,
)
    n = length(mass)
    ρ = zeros(n)
    @inbounds for i in 1:n
        idx = idxs[i]
        d = dists[i]
        skip = _self_slot(idx, i)
        m_inner = 0.0
        r_k = 0.0
        rank = 0
        for p in eachindex(idx)
            p == skip && continue
            rank += 1
            if rank < k
                m_inner += mass[idx[p]]
            else
                r_k = Float64(d[p])
                break
            end
        end
        numerator = by_mass ? m_inner : Float64(k - 1)
        ρ[i] = r_k > 0 ? numerator / (4π / 3 * r_k^3) : 0.0
    end
    return ρ
end

# ---------------------------------------------------------------------------
# Snapshot binary finder (node joining after Bate 2009a, one level)
# ---------------------------------------------------------------------------

"""
    SnapshotPair

A bound pair found among the stars of a snapshot by
[`snapshot_pairs`](@ref), with the elements of its relative Kepler orbit
(G = 1, the units of the input arrays).

# Fields
- `i`, `j`: star indices, `i` the more massive component
- `a`: semi-major axis of the relative orbit
- `e`: eccentricity
- `binding_energy`: ``m_i m_j / (2a)``
- `hardness`: ``x = E_\\mathrm{b} / kT``; `NaN` when no `kT` was supplied
- `r_third`: distance from the pair's centre of mass to the nearest other
  star; `Inf` when there is none
- `tidal_ok`: whether the apocentre clears the third body,
  ``a(1 + e) < f_\\mathrm{tid} r_\\mathrm{third}``
"""
struct SnapshotPair
    i::Int
    j::Int
    a::Float64
    e::Float64
    binding_energy::Float64
    hardness::Float64
    r_third::Float64
    tidal_ok::Bool
end

"""
    snapshot_pairs(pos, vel, mass; k = 10, tidal_factor = 1/3, kT = NaN)
        -> Vector{SnapshotPair}

Bound pairs among the stars `pos`, `vel` (3 × n) of masses `mass`, G = 1 in
any consistent units, sorted by decreasing binding energy.

Candidates are the unordered pairs `(i, j)` in which `j` is one of the `k`
nearest neighbours of `i` or vice versa, found on a k-d tree. For each
candidate the relative orbit gives the specific energy
``ε = |\\mathbf{v}|^2 / 2 - μ / |\\mathbf{r}|`` with ``μ = m_i + m_j``; a
bound candidate (`ε < 0`) has ``a = -μ / (2ε)`` and
``e = \\sqrt{1 - h^2 / (a μ)}``, `h` the specific angular momentum. Bound
candidates are joined greedily in ascending ``ε = -(m_i + m_j)/2a``, most
bound first, each star used once: a pair is accepted when neither star
belongs to an accepted pair, so the first pairs accepted are mutually most
bound, and a pair is accepted as well when the most-bound choice of one of
its stars was already taken by another pair. This is one level of node
joining after Bate (2009a, DOI 10.1111/j.1365-2966.2008.14106.x), who joins
the closest bound pair of nodes, ordered by ε instead of separation.
Hierarchical nodes (an accepted pair bound to a third body) are not built,
so a triple appears as its inner pair. Coincident stars, whose relative
orbit is undefined, are not paired.

The candidate restriction means that a pair is found only when it is closer
than the `k`-th neighbour of at least one of its stars; wide soft pairs
beyond that distance are not searched.

Two-body boundedness ignores the environment: a pair that is bound in
isolation but sits inside a dense core is tidally unbound. `tidal_ok` flags
the pairs whose apocentre `a(1 + e)` stays within `tidal_factor` times the
distance `r_third` from the pair's centre of mass to the nearest other star.

`kT` is the hard/soft energy scale in the units of the input, typically
``⟨m⟩ σ_\\mathrm{1D}^2`` of the bound set (in N-body units, the
[`hardness_scale`](@ref) product `m_mean * sigma_kms^2` divided by
`zmbar * vstar^2`); `hardness = E_b / kT` and a pair is hard when
`hardness > 1` (Heggie 1975, DOI 10.1093/mnras/173.3.729).

The components of the engine's KS pairs are stored in a snapshot at their
predicted positions and velocities, so the finder's elements for those pairs
agree with the engine's up to the single-precision rounding of the file.

Throws `ArgumentError` unless `k ≥ 1`, `tidal_factor > 0`, `kT` is positive
or `NaN` and there are at least two stars, and `DimensionMismatch` unless
`pos` and `vel` are both 3 × `length(mass)`.

# Example
```julia
pairs = snapshot_pairs(Float64.(snap.pos), Float64.(snap.vel), Float64.(snap.mass); kT)
count(p -> p.tidal_ok && p.hardness > 1, pairs)
```
"""
function snapshot_pairs(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    k::Int = 10,
    tidal_factor::Real = 1 / 3,
    kT::Real = NaN,
)
    n = length(mass)
    size(pos) == size(vel) ||
        throw(DimensionMismatch("snapshot_pairs: pos is $(size(pos)), vel is $(size(vel))"))
    size(pos) == (3, n) ||
        throw(DimensionMismatch("snapshot_pairs: pos is $(size(pos)), expected (3, $n)"))
    n ≥ 2 || throw(ArgumentError("snapshot_pairs needs at least 2 stars, got $n"))
    k ≥ 1 || throw(ArgumentError("snapshot_pairs: k must be at least 1, got $k"))
    tidal_factor > 0 ||
        throw(ArgumentError("snapshot_pairs: tidal_factor must be positive, got $tidal_factor"))
    (isnan(kT) || kT > 0) ||
        throw(ArgumentError("snapshot_pairs: kT must be positive or NaN, got $kT"))
    tree = KDTree(pos)
    kk = min(k, n - 1)
    idxs, _ = knn(tree, pos, kk + 1, true)
    return _snapshot_pairs(tree, idxs, pos, vel, mass, Float64(tidal_factor), Float64(kT))
end

"""
    _relative_orbit(pos, vel, mass, i, j) -> (; r, ε)

Separation and specific orbital energy of the relative orbit of stars `i`
and `j` (G = 1).
"""
@inline function _relative_orbit(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    i::Int,
    j::Int,
)
    @inbounds begin
        rx = pos[1, i] - pos[1, j]
        ry = pos[2, i] - pos[2, j]
        rz = pos[3, i] - pos[3, j]
        vx = vel[1, i] - vel[1, j]
        vy = vel[2, i] - vel[2, j]
        vz = vel[3, i] - vel[3, j]
        μ = mass[i] + mass[j]
    end
    r = sqrt(rx * rx + ry * ry + rz * rz)
    ε = 0.5 * (vx * vx + vy * vy + vz * vz) - μ / r
    return (; r, ε)
end

"""
    _snapshot_pairs(tree, idxs, pos, vel, mass, tidal_factor, kT) -> Vector{SnapshotPair}

Candidate collection, node joining and pair elements behind
[`snapshot_pairs`](@ref), from the sorted neighbour lists `idxs` of the
tree query at every star.
"""
function _snapshot_pairs(
    tree::KDTree,
    idxs::AbstractVector{<:AbstractVector{<:Integer}},
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    tidal_factor::Float64,
    kT::Float64,
)
    n = length(mass)
    candidates = NTuple{2,Int}[]
    sizehint!(candidates, sum(length, idxs))
    @inbounds for i in 1:n
        idx = idxs[i]
        skip = _self_slot(idx, i)
        for p in eachindex(idx)
            p == skip && continue
            j = Int(idx[p])
            push!(candidates, i < j ? (i, j) : (j, i))
        end
    end
    unique!(sort!(candidates))

    bound = Tuple{Float64,Int,Int}[]
    for (i, j) in candidates
        orbit = _relative_orbit(pos, vel, mass, i, j)
        orbit.r > 0 && orbit.ε < 0 && push!(bound, (orbit.ε, i, j))
    end
    sort!(bound)

    paired = falses(n)
    pairs = SnapshotPair[]
    centre = zeros(3)
    n_query = min(3, n)
    for (ε, i, j) in bound
        (paired[i] || paired[j]) && continue
        paired[i] = paired[j] = true
        mi, mj = mass[i], mass[j]
        μ = mi + mj
        a = -μ / (2ε)
        @inbounds begin
            rx = pos[1, i] - pos[1, j]
            ry = pos[2, i] - pos[2, j]
            rz = pos[3, i] - pos[3, j]
            vx = vel[1, i] - vel[1, j]
            vy = vel[2, i] - vel[2, j]
            vz = vel[3, i] - vel[3, j]
            for c in 1:3
                centre[c] = (mi * pos[c, i] + mj * pos[c, j]) / μ
            end
        end
        hx = ry * vz - rz * vy
        hy = rz * vx - rx * vz
        hz = rx * vy - ry * vx
        h² = hx * hx + hy * hy + hz * hz
        e = sqrt(max(1 - h² / (a * μ), 0.0))

        nbr, dist = knn(tree, centre, n_query, true)
        r_third = Inf
        for p in eachindex(nbr)
            (nbr[p] == i || nbr[p] == j) && continue
            r_third = Float64(dist[p])
            break
        end

        E_b = mi * mj / (2a)
        primary, secondary = mj > mi ? (j, i) : (i, j)
        push!(
            pairs,
            SnapshotPair(
                primary,
                secondary,
                a,
                e,
                E_b,
                E_b / kT,
                r_third,
                a * (1 + e) < tidal_factor * r_third,
            ),
        )
    end
    sort!(pairs; by = p -> p.binding_energy, rev = true)
    return pairs
end

# ---------------------------------------------------------------------------
# Reconciliation with the engine's regularised pairs
# ---------------------------------------------------------------------------

"""
    reconcile_pairs(pairs, snap::Snapshot)
        -> (; n_snapshot, n_ks, n_both, n_snapshot_only, n_ks_only, ks_found::BitVector)

Compare the pairs found by [`snapshot_pairs`](@ref) on `snap` with the
engine's KS pairs ([`pair_indices`](@ref)). A KS pair counts as found when a
`SnapshotPair` holds the same two star indices in either order;
`ks_found[k]` records this for KS pair k. `n_snapshot` and `n_ks` are the
sizes of the two lists, `n_both` the found KS pairs, `n_ks_only` the KS pairs
not found and `n_snapshot_only` the finder pairs that are not KS pairs.

The KS list is the hard-and-close subset of the binary population, since
regularisation requires a separation below RMIN. The finder adds the pairs
that are bound but wider than RMIN; for massive primaries these are
themselves hard.
"""
function reconcile_pairs(pairs::AbstractVector{SnapshotPair}, snap::Snapshot)
    found = Set{NTuple{2,Int}}()
    sizehint!(found, length(pairs))
    for p in pairs
        push!(found, p.i < p.j ? (p.i, p.j) : (p.j, p.i))
    end
    ks = pair_indices(snap)
    ks_set = Set{NTuple{2,Int}}()
    sizehint!(ks_set, length(ks))
    ks_found = falses(length(ks))
    for (k, (a, b)) in enumerate(ks)
        key = a < b ? (a, b) : (b, a)
        push!(ks_set, key)
        ks_found[k] = key in found
    end
    n_both = count(ks_found)
    n_snapshot_only = count(key -> !(key in ks_set), found)
    return (;
        n_snapshot = length(pairs),
        n_ks = length(ks),
        n_both,
        n_snapshot_only,
        n_ks_only = length(ks) - n_both,
        ks_found,
    )
end
