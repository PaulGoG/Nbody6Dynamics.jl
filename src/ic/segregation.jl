# =============================================================================
# Primordial mass segregation of the sampled systems
# =============================================================================
#
# Baumgardt, De Marchi & Kroupa (2008), Appendix A, with the degree of
# segregation S of McLuster (Küpper et al. 2011, Eq. A1). Included before
# config.jl because `SegregationSpec` is a field type of `ClusterSpec`; the
# profile tags and samplers it calls are resolved at run time.

"""Values of `[merger.clusterN.segregation] method`."""
const _SEGREGATION_METHODS = ("none", "bdk08")

"""
    SegregationSpec(; method = "none", degree = 0.0)

Primordial mass segregation of a cluster (`[merger.clusterN.segregation]`).

# Fields
- `method::String`: `"none"` (unsegregated sampling) or `"bdk08"`.
- `degree::Float64`: degree of segregation `S`, `0 ≤ S ≤ 1`; `0` is a random
  assignment, `1` complete ordering.

`"bdk08"` is the scheme of Baumgardt, De Marchi & Kroupa (2008, Appendix A;
DOI 10.1086/590488) with the partial ordering of McLuster (Küpper et al. 2011,
Eq. A1; DOI 10.1111/j.1365-2966.2011.19412.x). `N′ ≥ N ⟨m⟩/m_min` phase-space
points of the unsegregated model are ranked by specific energy, most bound
first; the systems, sorted by mass in descending order and permuted by the
McLuster slot order of degree `S`, take one point each from consecutive
energy slices of `N′ m/M` points. The mass-density profile of the model is
preserved, so the cluster stays in virial equilibrium. The scheme acts on
system masses: a primordial pair is placed as its centre of mass, before the
pairs are expanded.
"""
Base.@kwdef struct SegregationSpec
    method::String = "none"
    degree::Float64 = 0.0
end

"""`SegregationSpec` as a TOML table (metadata round trip)."""
_segregation_table(s::SegregationSpec) =
    Dict{String,Any}("method" => s.method, "degree" => s.degree)

"""`[merger.cluster<idx>.segregation]` as a [`SegregationSpec`](@ref)."""
_parse_segregation_table(d::AbstractDict, idx::Int)::SegregationSpec =
    _parse_section(SegregationSpec, d, "merger.cluster$idx.segregation")

"""
    _validate_segregation(s::SegregationSpec, key::String)

Refuse an unknown `method` and a `degree` outside `[0, 1]`, naming
`<key>.segregation.<field>`.
"""
function _validate_segregation(s::SegregationSpec, key::String)
    s.method in _SEGREGATION_METHODS || _config_error(
        "$key.segregation.method",
        "must be one of \"none\", \"bdk08\"; got \"$(s.method)\"",
    )
    (0 ≤ s.degree ≤ 1) ||
        _config_error("$key.segregation.degree", "must satisfy 0 ≤ S ≤ 1; got $(s.degree)")
    return nothing
end

function _check_segregation_arguments(n::Integer, S::Real)
    n ≥ 1 || throw(ArgumentError("segregation: the number of systems must be ≥ 1; got $n"))
    (0 ≤ S ≤ 1) || throw(ArgumentError("segregation: the degree S must satisfy 0 ≤ S ≤ 1; got $S"))
    return nothing
end

"""
    _segregation_relative_slots(n, S, rng) -> Vector{Float64}

The `n` relative slots `yᵢ = 1 − Xᵢ^(1−S)`, `Xᵢ` uniform in `[0, 1)`, of the
McLuster slot order, drawn with the same random-number use as
[`segregation_order`](@ref) for `S < 1`. They follow
`P(y ≤ Y) = 1 − (1 − Y)^{1/(1−S)}`.
"""
function _segregation_relative_slots(n::Integer, S::Real, rng::AbstractRNG)::Vector{Float64}
    _check_segregation_arguments(n, S)
    e = 1.0 - Float64(S)
    return Float64[1.0 - rand(rng)^e for _ in 1:n]
end

# Fenwick tree over the slot occupancy: `tree[i]` counts the free slots in
# (i − lowbit(i), i]; all n slots start free.
_fenwick_all_free(n::Int)::Vector{Int} = Int[i & -i for i in 1:n]

# Slot (1-based) of the k-th free slot in increasing order, which is then
# marked occupied; O(log n).
function _fenwick_take!(tree::Vector{Int}, k::Int)::Int
    n = length(tree)
    pos = 0
    rem = k
    step = prevpow(2, n)
    while step > 0
        nxt = pos + step
        if nxt ≤ n && tree[nxt] < rem
            pos = nxt
            rem -= tree[nxt]
        end
        step >>= 1
    end
    slot = pos + 1
    i = slot
    while i ≤ n
        tree[i] -= 1
        i += i & -i
    end
    return slot
end

"""
    segregation_order(n::Integer, S::Real, rng::AbstractRNG) -> Vector{Int}

Slot order of degree `S` after McLuster (Küpper et al. 2011, Eq. A1): with the
systems sorted by mass in descending order, system `i` takes the free slot of
rank `j = min(⌊yᵢ (n − i + 1)⌋, n − i)` (0-based, increasing slot order) among
the `n − i + 1` still free, `yᵢ = 1 − Xᵢ^(1−S)`. Returns `slot`, `slot[i]` the
slot (1-based) of the `i`-th most massive system. `S = 1` is the identity
(nothing is drawn), `S = 0` a uniformly random permutation. O(n log n).
Throws an `ArgumentError` unless `n ≥ 1` and `0 ≤ S ≤ 1`.

# Example
```julia
segregation_order(5, 1.0, Xoshiro(1))   # [1, 2, 3, 4, 5]
```
"""
function segregation_order(n::Integer, S::Real, rng::AbstractRNG)::Vector{Int}
    _check_segregation_arguments(n, S)
    m = Int(n)
    S == 1 && return collect(1:m)
    e = 1.0 - Float64(S)
    tree = _fenwick_all_free(m)
    slot = Vector{Int}(undef, m)
    for i in 1:m
        y = 1.0 - rand(rng)^e
        free = m - i + 1
        j = min(floor(Int, y * free), free - 1)
        slot[i] = _fenwick_take!(tree, j + 1)
    end
    return slot
end

# Phase-space points of the unsegregated model with their specific energies in
# the model's units: Plummer with G = M = 1, E = v²/2 − 1/√(r² + a²); King in
# units of σ², E = v²/2 − Ŵ(r̂).
function _segregation_model_points(profile, n_pts::Int, rh::Float64, rng::AbstractRNG)
    if profile isa PlummerProfile
        a = rh / _PLUMMER_RHM_OVER_A
        pos, vel = sample_plummer(n_pts, a; rng = rng)
        energy = Vector{Float64}(undef, n_pts)
        for k in 1:n_pts
            r2 = pos[1, k]^2 + pos[2, k]^2 + pos[3, k]^2
            v2 = vel[1, k]^2 + vel[2, k]^2 + vel[3, k]^2
            energy[k] = 0.5 * v2 - 1.0 / sqrt(r2 + a^2)
        end
        return pos, vel, energy
    elseif profile isa KingProfile
        return _sample_king_with_energy(n_pts, profile.W0, 1.0; rng = rng)
    else
        throw(ArgumentError("segregate_systems: unsupported profile $(typeof(profile))"))
    end
end

"""
    segregate_systems(profile, mass, rh; degree, rng) -> (pos, vel, energy)

Phase-space coordinates of the systems of masses `mass` (any order) drawn from
`profile` with primordial mass segregation of degree `degree` after
Baumgardt, De Marchi & Kroupa (2008, Appendix A):

1. `N′ = max(n, ⌈n ⟨m⟩/m_min⌉)` points of the unsegregated model — Plummer with
   scale radius `rh / 1.305`, King with unit tidal radius (`rh` unused), as
   `_sample_cluster` samples them — ranked by specific energy, most bound first.
2. The systems, sorted by mass in descending order (stable), are arranged in
   the slot order [`segregation_order`](@ref)`(n, degree, rng)`; the system in
   slot `k` takes a point drawn uniformly among the energy ranks
   `⌊N′ M_cum[k−1]⌋ + 1 … max(⌊N′ M_cum[k−1]⌋ + 1, ⌊N′ M_cum[k]⌋)`, with `M_cum`
   the cumulative mass fraction in slot order.

Returns positions and velocities (3 × n, the sampler's units) and the specific
energies of the chosen points, all in the order of `mass`. `degree = 1` gives
energies non-decreasing in descending mass; `degree = 0` a random assignment.
Throws an `ArgumentError` for an empty `mass`, a non-positive mass or
`degree ∉ [0, 1]`.
"""
function segregate_systems(
    profile,
    mass::AbstractVector{Float64},
    rh::Real;
    degree::Real,
    rng::AbstractRNG,
)::Tuple{Matrix{Float64},Matrix{Float64},Vector{Float64}}
    n = length(mass)
    _check_segregation_arguments(n, degree)
    all(m -> isfinite(m) && m > 0, mass) ||
        throw(ArgumentError("segregate_systems: every mass must be finite and > 0"))
    M = sum(mass)
    n_pts = max(n, ceil(Int, M / minimum(mass)))

    pts_pos, pts_vel, pts_energy = _segregation_model_points(profile, n_pts, Float64(rh), rng)
    by_energy = sortperm(pts_energy)

    by_mass = sortperm(mass; rev = true)
    slot = segregation_order(n, degree, rng)
    in_slot = Vector{Int}(undef, n)
    for i in 1:n
        in_slot[slot[i]] = by_mass[i]
    end

    pos = Matrix{Float64}(undef, 3, n)
    vel = Matrix{Float64}(undef, 3, n)
    energy = Vector{Float64}(undef, n)
    m_cum = 0.0
    for k in 1:n
        s = in_slot[k]
        lo = min(floor(Int, n_pts * (m_cum / M)) + 1, n_pts)
        m_cum += mass[s]
        hi = min(max(lo, floor(Int, n_pts * (m_cum / M))), n_pts)
        p = by_energy[rand(rng, lo:hi)]
        for d in 1:3
            pos[d, s] = pts_pos[d, p]
            vel[d, s] = pts_vel[d, p]
        end
        energy[s] = pts_energy[p]
    end
    return pos, vel, energy
end
