# =============================================================================
# Merger epochs: contact, core coalescence, dissolution and the regimes
# =============================================================================
#
# A merger is a sequence of epochs rather than one time: contact of the
# self-bound cores (Roche-lobe overflow; Priyatikanto et al. 2016, DOI
# 10.1093/mnras/stw060), coalescence of the most-bound cores in phase space
# (Han et al. 2018, DOI 10.1093/mnras/stx2792; Karam & Sills 2022, DOI
# 10.1093/mnras/stac1298), then spatial and kinematic mixing. The epochs set
# the regime of every snapshot, progenitor and pair, and the regime decides in
# which frame a tag-resolved quantity is physical.

"""Regimes of a snapshot, a progenitor or a pair, in the order of a merger."""
const _MERGER_REGIMES = (:detached, :interacting, :coalescing, :settled, :dissolved, :satellite)

"""Regimes in which the progenitors form one body (the remnant frame applies)."""
const _REMNANT_REGIMES = (:coalescing, :settled)

"""
    MergerEpochs

Epochs of the merger of the initial clusters (the tags of `cluster_ranges`) and
the regime of every snapshot, progenitor and pair, from [`merger_epochs`](@ref).
Indices are snapshot indices, `0` for "never"; times are N-body times with a
`_myr` companion, `NaN` for "never". The regimes are `:detached`, `:interacting`,
`:coalescing`, `:settled`, `:dissolved` and `:satellite` ([`regime_label`](@ref)).

# Fields
System:
- `time`, `time_myr`: snapshot times
- `contact_index`, `t_contact`, `t_contact_myr`: first contact of any pair
- `core_index`, `t_core`, `t_core_myr`: first snapshot from which the
  progenitors form one remnant for the rest of the run: every progenitor
  present is `coalescing`, `settled`, a `satellite` or `dissolved`, and at
  least one is `coalescing` or `settled`
- `settled_index`: first `settled` snapshot from `core_index` on
- `coalescence_index`, `t_coal`, `t_coal_myr`: the centre coalescence of
  [`coalescence_time`](@ref), kept for continuity
- `t_cr_core`: crossing time of the bound set at `core_index`
- `l_orb`: orbital angular momentum of the most-bound cores about their centre
  of mass at the last snapshot before `core_index`, where the cores are still
  distinct (NB units; `NaN` without a core coalescence)
- `regime`: regime of every snapshot (`n_t`)

Pairs, in the order of [`coalescence_time`](@ref):
- `pairs`, `separation`, `threshold`: as returned by `coalescence_time`
- `pair_contact`: first snapshot with a Roche-lobe overflow of either core
- `pair_join`: snapshot from which the two clusters share a group of the
  coalescence test at every later snapshot ([`merger_tree`](@ref))
- `pair_core`: core coalescence of the merger-tree event that joins the groups
  of the two clusters
- `roche_fill`: `max(r_f,i / R_L,i, r_f,j / R_L,j)`, `n_pairs × n_t`
- `delta_s`: phase-space offset δs of the two most-bound cores, `n_pairs × n_t`
- `pair_regime`: `n_pairs × n_t`

Merger-tree events, in the order of [`merger_tree`](@ref):
- `group_a`, `group_b`: the tags of the two groups joined
- `event_join`: join index of the merger tree
- `event_contact`, `event_core`, `event_settled`: first contact of any pair
  across the two groups, core coalescence of their main progenitors, first
  `settled` snapshot after it
- `main_progenitor`: the most massive tag of the more massive group
- `mass_ratio`: bound mass of the lighter group over that of the heavier one,
  at the last snapshot before `event_contact`

Progenitors:
- `tag_contact`, `tag_core`, `tag_dissolved`: first contact with any other
  progenitor, first core coalescence of an event that includes the tag, and the
  dissolution index ([`dissolution_index`](@ref))
- `bound_mass`: mass of the tag's members in the global bound set, `n_cl × n_t`
- `self_bound_mass`: mass of the tag's self-bound subset, `n_cl × n_t`; `0`
  below `min_members`
- `tag_regime`: `n_cl × n_t`

Thresholds (the `merger_*` keys of [`PostprocessConfig`](@ref)):
`min_members`, `contact_lagrangian`, `core_n`, `core_delta`,
`dissolved_fraction`, `coalescing_tcr`.
"""
struct MergerEpochs
    time::Vector{Float64}
    time_myr::Vector{Float64}
    contact_index::Int
    t_contact::Float64
    t_contact_myr::Float64
    core_index::Int
    t_core::Float64
    t_core_myr::Float64
    settled_index::Int
    coalescence_index::Int
    t_coal::Float64
    t_coal_myr::Float64
    t_cr_core::Float64
    l_orb::Vector{Float64}
    regime::Vector{Symbol}
    pairs::Vector{Tuple{Int,Int}}
    separation::Matrix{Float64}
    threshold::Matrix{Float64}
    pair_contact::Vector{Int}
    pair_join::Vector{Int}
    pair_core::Vector{Int}
    roche_fill::Matrix{Float64}
    delta_s::Matrix{Float64}
    pair_regime::Matrix{Symbol}
    group_a::Vector{Vector{Int}}
    group_b::Vector{Vector{Int}}
    event_join::Vector{Int}
    event_contact::Vector{Int}
    event_core::Vector{Int}
    event_settled::Vector{Int}
    main_progenitor::Vector{Int}
    mass_ratio::Vector{Float64}
    tag_contact::Vector{Int}
    tag_core::Vector{Int}
    tag_dissolved::Vector{Int}
    bound_mass::Matrix{Float64}
    self_bound_mass::Matrix{Float64}
    tag_regime::Matrix{Symbol}
    min_members::Int
    contact_lagrangian::Float64
    core_n::Int
    core_delta::Float64
    dissolved_fraction::Float64
    coalescing_tcr::Float64
end

# ---------------------------------------------------------------------------
# Detectors
# ---------------------------------------------------------------------------

"""
    _roche_lobe_radius(d, q) -> Float64

Roche-lobe radius of a body of mass ratio `q = M_1 / M_2` at separation `d`
(Eggleton 1983, DOI 10.1086/160960):
`R_L = d · 0.49 q^{2/3} / (0.6 q^{2/3} + ln(1 + q^{1/3}))`, accurate to 1 %
for every `0 < q < ∞`.
"""
function _roche_lobe_radius(d::Real, q::Real)::Float64
    q23 = cbrt(q)^2
    return d * 0.49 * q23 / (0.6 * q23 + log1p(cbrt(q)))
end

"""
    _roche_fill(r_i, r_j, m_i, m_j, d) -> Float64

Larger Roche-lobe filling factor `r / R_L` of two bodies with Lagrangian radii
`r_i`, `r_j`, masses `m_i`, `m_j` and separation `d`; `NaN` when an input is
not finite or a mass not positive.
"""
function _roche_fill(r_i::Real, r_j::Real, m_i::Real, m_j::Real, d::Real)::Float64
    all(isfinite, (r_i, r_j, m_i, m_j, d)) || return NaN
    (m_i > 0 && m_j > 0 && d ≥ 0) || return NaN
    d == 0 && return Inf
    return max(r_i / _roche_lobe_radius(d, m_i / m_j), r_j / _roche_lobe_radius(d, m_j / m_i))
end

"""
    contact_index(radius_i, radius_j, mass_i, mass_j, separation) -> Int

First snapshot at which either of two self-bound cores overflows its Roche
lobe: `radius_i[k] > R_L(separation[k], mass_i[k] / mass_j[k])` or the same
with `i` and `j` exchanged, `R_L` the Roche-lobe radius of Eggleton (1983, DOI
10.1086/160960). The radii are the Lagrangian radii of the self-bound cores
at the mass fraction `merger_contact_lagrangian`, the masses their self-bound
masses and `separation` the distance of their density centres; the contact
epoch follows Priyatikanto et al. (2016, DOI 10.1093/mnras/stw060). Snapshots
with a non-finite entry or a non-positive mass are skipped. No persistence is
required: a fly-by that separates again stays past contact. Returns `0` when
the cores never touch.

Throws `DimensionMismatch` when the series differ in length.

# Example
```julia
contact_index([0.5, 0.9], [0.5, 0.9], [1.0, 1.0], [1.0, 1.0], [4.0, 1.5])   # 2
```
"""
function contact_index(
    radius_i::AbstractVector{<:Real},
    radius_j::AbstractVector{<:Real},
    mass_i::AbstractVector{<:Real},
    mass_j::AbstractVector{<:Real},
    separation::AbstractVector{<:Real},
)::Int
    n = length(separation)
    all(v -> length(v) == n, (radius_i, radius_j, mass_i, mass_j)) ||
        throw(DimensionMismatch("contact_index: the five series must have equal length"))
    for k in 1:n
        f = _roche_fill(radius_i[k], radius_j[k], mass_i[k], mass_j[k], separation[k])
        isfinite(f) || f == Inf || continue
        f > 1 && return k
    end
    return 0
end

"""
    core_index(delta_s; threshold = 2.0) -> Int

First snapshot from which the phase-space offset of two most-bound cores,
`δs = Δx / σ_x + Δv / σ_v` (Han et al. 2018, DOI 10.1093/mnras/stx2792, Eq. 3),
stays below `threshold` for the rest of the run. `Δx`, `Δv` are the offsets of
the two cores in position and velocity and `σ_x`, `σ_v` the position and
velocity dispersions of the host core. Non-finite entries (a core absent) are
skipped. The scan runs back from the last snapshot over the finite entries
and stops at the first `δs ≥ threshold`; the earliest snapshot of that run is
the core coalescence, provided it is followed by a later finite entry or is
the last snapshot, so that a pericentre passage that separates again does not
count. Returns `0` when the cores never coalesce.

Throws `ArgumentError` for `threshold ≤ 0`.

# Example
```julia
core_index([5.0, 3.0, 1.5, 1.2, 0.8])   # 3
```
"""
function core_index(delta_s::AbstractVector{<:Real}; threshold::Real = 2.0)::Int
    threshold > 0 || throw(ArgumentError("core_index: threshold must be > 0; got $threshold"))
    n = length(delta_s)
    idx = 0
    later = false
    for k in n:-1:1
        isfinite(delta_s[k]) || continue
        delta_s[k] < threshold || break
        (later || k == n) && (idx = k)
        later = true
    end
    return idx
end

"""
    dissolution_index(bound_mass, initial_mass; fraction = 0.1) -> Int

First snapshot from which the self-bound mass of a progenitor stays below
`fraction × initial_mass` for the rest of the run: the secondary has been
shredded (de la Fuente Marcos & de la Fuente Marcos 2010, DOI
10.1088/0004-637X/719/1/104; Han et al. 2018, DOI 10.1093/mnras/stx2792).
A non-finite entry counts as no self-bound core (mass 0). Returns `0` when the
progenitor keeps a self-bound core above the limit at the last snapshot.

Throws `ArgumentError` unless `0 < fraction < 1` and `initial_mass > 0`.

# Example
```julia
dissolution_index([1.0, 0.6, 0.05, 0.02], 1.0)   # 3
```
"""
function dissolution_index(
    bound_mass::AbstractVector{<:Real},
    initial_mass::Real;
    fraction::Real = 0.1,
)::Int
    0 < fraction < 1 ||
        throw(ArgumentError("dissolution_index: fraction must lie in (0, 1); got $fraction"))
    initial_mass > 0 ||
        throw(ArgumentError("dissolution_index: initial_mass must be > 0; got $initial_mass"))
    limit = fraction * initial_mass
    idx = 0
    for k in length(bound_mass):-1:1
        m = bound_mass[k]
        (isfinite(m) ? Float64(m) : 0.0) < limit || break
        idx = k
    end
    return idx
end

"""
    regime_label(k, contact, core, settled; dissolved = 0, satellite = false) -> Symbol

Regime of snapshot `k` from the epoch indices of a progenitor, a pair or the
whole system (`0` for an epoch never reached), in this order of precedence:

1. `:coalescing` from `core` on, `:settled` from `settled` on;
2. before `core`: `:dissolved` from `dissolved` on (a shredded progenitor);
3. `:satellite` when `satellite` is set (a progenitor in contact with an
   already coalesced group that has not coalesced with it);
4. `:interacting` from `contact` on;
5. `:detached` otherwise.

# Example
```julia
regime_label(5, 3, 6, 9)    # :interacting
regime_label(7, 3, 6, 9)    # :coalescing
```
"""
function regime_label(
    k::Integer,
    contact::Integer,
    core::Integer,
    settled::Integer;
    dissolved::Integer = 0,
    satellite::Bool = false,
)::Symbol
    if core > 0 && k ≥ core
        return (settled > 0 && k ≥ settled) ? :settled : :coalescing
    end
    dissolved > 0 && k ≥ dissolved && return :dissolved
    satellite && return :satellite
    contact > 0 && k ≥ contact && return :interacting
    return :detached
end

# ---------------------------------------------------------------------------
# Groups of the coalescence test and the events of the merger tree
# ---------------------------------------------------------------------------

"""
    _coalescence_groups(co) -> Matrix{Int}

Group of every cluster at every snapshot under the test of
[`coalescence_time`](@ref): the root tag of its connected component in the
graph of pairs closer than their threshold, over the clusters present;
`0` where a cluster is absent. `co` is the result of `coalescence_time`.
"""
function _coalescence_groups(co)::Matrix{Int}
    present = co.present
    n_cl, n_t = size(present)
    groups = zeros(Int, n_cl, n_t)
    parent = collect(1:n_cl)
    for k in 1:n_t
        for i in 1:n_cl
            parent[i] = i
        end
        for (q, (i, j)) in enumerate(co.pairs)
            (present[i, k] && present[j, k]) || continue
            co.separation[q, k] < co.threshold[q, k] || continue
            ri = _union_root!(parent, i)
            rj = _union_root!(parent, j)
            ri == rj || (parent[ri] = rj)
        end
        for i in 1:n_cl
            present[i, k] && (groups[i, k] = _union_root!(parent, i))
        end
    end
    return groups
end

"""
    _tree_events(co) -> (; groups, joined, events)

Joins of the merger tree from the result `co` of [`coalescence_time`](@ref)
under its group rule. The pair `q = (i, j)` is joined from snapshot
`joined[q]`, the earliest `k` from which the two clusters belong to the same
group of `_coalescence_groups` at every later snapshot at which both are
present, provided `k` is followed by a later co-present snapshot or is the
last; `0` when never. The joined pairs, in increasing `joined[q]` and then pair
order, unite the groups by union–find; every pair that unites two distinct
groups is an event `(; pair, index, group_a, group_b)`, `group_a` holding the
smallest tag.
"""
function _tree_events(co)
    groups = _coalescence_groups(co)
    present = co.present
    pairs = co.pairs
    n_cl, n_t = size(present)
    joined = zeros(Int, length(pairs))
    for (q, (i, j)) in enumerate(pairs)
        later = false
        @inbounds for k in n_t:-1:1
            (present[i, k] && present[j, k]) || continue
            groups[i, k] == groups[j, k] || break
            (later || k == n_t) && (joined[q] = k)
            later = true
        end
    end
    order = sort!(findall(>(0), joined); by = q -> (joined[q], q))
    parent = collect(1:n_cl)
    events =
        NamedTuple{(:pair, :index, :group_a, :group_b),Tuple{Int,Int,Vector{Int},Vector{Int}}}[]
    for q in order
        i, j = pairs[q]
        ri = _union_root!(parent, i)
        rj = _union_root!(parent, j)
        ri == rj && continue
        group_a = [t for t in 1:n_cl if _union_root!(parent, t) == ri]
        group_b = [t for t in 1:n_cl if _union_root!(parent, t) == rj]
        if first(group_b) < first(group_a)
            group_a, group_b = group_b, group_a
        end
        push!(events, (; pair = q, index = joined[q], group_a, group_b))
        parent[rj] = ri
    end
    return (; groups, joined, events)
end

# ---------------------------------------------------------------------------
# Per-progenitor series
# ---------------------------------------------------------------------------

"""
    _most_bound_core(snap, sel, n_core) -> (; pos, vel, sigma_x, sigma_v)

Core of the snapshot stars `sel`: the `n_core` systems (regularised pairs at
their centres of mass) of largest specific binding energy
`½|v − v̄|² + Φ`, `Φ` the potential of the systems of `sel` alone and `v̄` their
mass-weighted mean velocity. Returns the number-weighted centroid `pos` and
mean velocity `vel` of the core and its three-dimensional RMS dispersions
`sigma_x`, `sigma_v` about them; `NaN` for fewer than two systems. Number
weighting follows the shrinking-sphere centre: with a broad mass spectrum a
mass-weighted core locks onto the few massive remnants. O(n²) in `length(sel)`.
"""
function _most_bound_core(snap::Snapshot, sel::Vector{Int}, n_core::Integer)
    sys, _ = _snapshot_systems(snap, sel)
    n_sys = length(sys.mass)
    if n_sys < 2
        return (; pos = fill(NaN, 3), vel = fill(NaN, 3), sigma_x = NaN, sigma_v = NaN)
    end
    Φ = pair_potential(sys.pos, sys.mass)
    v̄ = _mean_velocity(sys, 1:n_sys)
    e = Vector{Float64}(undef, n_sys)
    @inbounds for s in 1:n_sys
        e[s] =
            0.5 * ((sys.vel[1, s] - v̄[1])^2 + (sys.vel[2, s] - v̄[2])^2 + (sys.vel[3, s] - v̄[3])^2) +
            Φ[s]
    end
    n = min(Int(n_core), n_sys)
    top = partialsortperm(e, 1:n)
    x_c = zeros(3)
    v_c = zeros(3)
    @inbounds for s in top, a in 1:3
        x_c[a] += sys.pos[a, s]
        v_c[a] += sys.vel[a, s]
    end
    x_c ./= n
    v_c ./= n
    sx = 0.0
    sv = 0.0
    @inbounds for s in top, a in 1:3
        sx += (sys.pos[a, s] - x_c[a])^2
        sv += (sys.vel[a, s] - v_c[a])^2
    end
    return (; pos = x_c, vel = v_c, sigma_x = sqrt(sx / n), sigma_v = sqrt(sv / n))
end

"""
    _epoch_tag_series(snaps, cluster_ranges, bound_sets, settings; max_passes, min_members)

Per progenitor and snapshot: the members in the global bound set (count
`n_bound`, mass `bound_mass`), the self-bound subset of the tag's members
([`bound_set`](@ref) with `max_passes`, `min_members`; mass
`self_bound_mass`, shrinking-sphere density centre `self_centre` with the
centre keys of `settings`, Lagrangian radius `self_radius` at
`merger_contact_lagrangian`), and the most-bound core of the members in the
global bound set (`_most_bound_core` with `merger_core_n`). A subset
with fewer than `merger_min_members` members counts as absent (mass 0, other
entries `NaN`). Also the initial mass of every tag, `m_initial`.
"""
function _epoch_tag_series(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    bound_sets::AbstractVector{BoundSet},
    settings::PostprocessConfig;
    max_passes::Integer,
    min_members::Integer,
)
    n_cl = length(cluster_ranges)
    n_t = length(snaps)
    n_min = settings.merger_min_members
    keywords = _centre_keywords(settings)
    fraction = (settings.merger_contact_lagrangian,)
    n_bound = zeros(Int, n_cl, n_t)
    bound_mass = zeros(n_cl, n_t)
    self_bound_mass = zeros(n_cl, n_t)
    self_centre = fill(NaN, 3, n_cl, n_t)
    self_radius = fill(NaN, n_cl, n_t)
    core_pos = fill(NaN, 3, n_cl, n_t)
    core_vel = fill(NaN, 3, n_cl, n_t)
    core_sigma_x = fill(NaN, n_cl, n_t)
    core_sigma_v = fill(NaN, n_cl, n_t)
    m_initial = [
        sum(Float64(snaps[1].mass[s]) for s in _member_indices(snaps[1], rng); init = 0.0) for
        rng in cluster_ranges
    ]
    for (k, snap) in enumerate(snaps)
        in_bound = falses(nparticles(snap))
        in_bound[bound_sets[k].stars] .= true
        for (i, rng) in enumerate(cluster_ranges)
            members = _member_indices(snap, rng)
            length(members) < n_min && continue
            # Self-bound subset of the tag (set B)
            b = bound_set(snap, members; max_passes, min_members)
            if b.status !== :dissolved && length(b.stars) ≥ n_min
                p = Float64.(snap.pos[:, b.stars])
                m = Float64.(snap.mass[b.stars])
                c = _shrinking_sphere_centre(p, m; keywords...)
                self_bound_mass[i, k] = sum(m)
                self_centre[:, i, k] .= c
                self_radius[i, k] = _lagrangian_radii(p, m, c, fraction)[1]
            end
            # Members in the global bound set (set A) and their most-bound core
            sel = [s for s in members if in_bound[s]]
            n_bound[i, k] = length(sel)
            bound_mass[i, k] = sum(Float64(snap.mass[s]) for s in sel; init = 0.0)
            length(sel) < n_min && continue
            core = _most_bound_core(snap, sel, settings.merger_core_n)
            core_pos[:, i, k] .= core.pos
            core_vel[:, i, k] .= core.vel
            core_sigma_x[i, k] = core.sigma_x
            core_sigma_v[i, k] = core.sigma_v
        end
    end
    return (;
        m_initial,
        n_bound,
        bound_mass,
        self_bound_mass,
        self_centre,
        self_radius,
        core_pos,
        core_vel,
        core_sigma_x,
        core_sigma_v,
    )
end

"""Distance of the columns `i` and `j` of the slice `x[:, :, k]`."""
@inline function _column_distance(x::Array{Float64,3}, i::Int, j::Int, k::Int)::Float64
    return sqrt(
        (x[1, i, k] - x[1, j, k])^2 + (x[2, i, k] - x[2, j, k])^2 + (x[3, i, k] - x[3, j, k])^2,
    )
end

"""Crossing time `M^{5/2} / (−2E)^{3/2}` of the bound set `b` of `snap`; `NaN` when unbound or empty."""
function _bound_crossing_time(snap::Snapshot, b::BoundSet)::Float64
    isempty(b.stars) && return NaN
    M = sum(Float64(snap.mass[i]) for i in b.stars)
    return crossing_time(M, b.kinetic + b.potential)
end

# ---------------------------------------------------------------------------
# Epochs
# ---------------------------------------------------------------------------

"""
    merger_epochs(snaps, cluster_ranges, bound_sets, coalescence = nothing;
                  settings = PostprocessConfig(), max_passes = settings.bound_max_passes,
                  min_members = settings.bound_min_members) -> MergerEpochs

Epochs of the merger of the initial clusters and the regime of every
snapshot, progenitor and pair.

`bound_sets` holds the global bound set of every snapshot (the `bound` of
[`snapshot_integrals`](@ref)); `coalescence` the result of
[`coalescence_time`](@ref) on them, evaluated here when not given. The groups
and joins are those of [`merger_tree`](@ref) (the group rule of
`coalescence_time`). The detectors:

- Contact ([`contact_index`](@ref)): the self-bound subset of every tag
  ([`bound_set`](@ref) of its members with `max_passes`, `min_members`) gives
  the self-bound mass, the shrinking-sphere density centre (centre keys of
  `settings`) and the Lagrangian radius at the mass fraction
  `merger_contact_lagrangian`; a pair is in contact from the first snapshot at
  which either radius exceeds the Roche-lobe radius of Eggleton (1983, DOI
  10.1086/160960) at the separation of the centres and the ratio of the
  self-bound masses (Priyatikanto et al. 2016, DOI 10.1093/mnras/stw060).
- Core coalescence ([`core_index`](@ref)): the core of a tag is the
  `merger_core_n` systems of largest binding energy among its members in the
  global bound set, in the potential of those members alone; `δs = Δx / σ_x +
  Δv / σ_v` with the offsets of the two cores and the dispersions of the core
  of the tag with the larger bound mass (Han et al. 2018, DOI
  10.1093/mnras/stx2792); coalescence when `δs < merger_core_delta` from then
  to the end. For a merger-tree event the core coalescence is that of the
  most massive tag of each group (bound mass at the last snapshot before the
  event's contact, which also sets `main_progenitor` and `mass_ratio`).
- Dissolution ([`dissolution_index`](@ref)): the self-bound mass stays below
  `merger_dissolved_fraction` of the tag's initial mass.

Every subset needs `merger_min_members` members to count as present. The
`coalescing` regime lasts `merger_coalescing_tcr` crossing times of the bound
set (`M^{5/2} / (−2E)^{3/2}`, Heggie & Mathieu 1986) after a core coalescence,
open-ended when that crossing time is undefined; `settled` follows.

The regime of a progenitor at a snapshot follows [`regime_label`](@ref) with
its first contact, the latest core coalescence of an event that includes it
and its dissolution; it is a `satellite` when it has not coalesced but is in
contact with a progenitor that has. A pair takes the core coalescence of the
event that joins its two groups. The system is `detached` before the first
contact of any pair, `interacting` until `core_index` (defined in
[`MergerEpochs`](@ref)), then `coalescing` and `settled` after the latest core
coalescence. A single progenitor is one body throughout: every regime is
`settled` and no epoch is evaluated.

Cost: one self-bound subset and one potential of the bound members per tag
and snapshot, O(Σᵢ Nᵢ²) per snapshot.

# Arguments
- `snaps`: snapshots in time order
- `cluster_ranges`: names of every initial cluster
- `bound_sets`: global bound set of every snapshot
- `coalescence`: result of `coalescence_time` on `snaps` and `bound_sets`
- `settings`: the `merger_*` and centre keys
- `max_passes`, `min_members`: limits of the per-tag self-bound subsets

Throws `ArgumentError` for no snapshots or no clusters, `DimensionMismatch`
when the bound sets or the coalescence result do not match the snapshots.

# Example
```julia
co = coalescence_time(snaps, ranges; bound_sets)
ep = merger_epochs(snaps, ranges, bound_sets, co; settings = cfg.postprocess)
ep.t_core, ep.regime
```
"""
function merger_epochs(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    bound_sets::AbstractVector{BoundSet},
    coalescence::Union{Nothing,NamedTuple} = nothing;
    settings::PostprocessConfig = PostprocessConfig(),
    max_passes::Integer = settings.bound_max_passes,
    min_members::Integer = settings.bound_min_members,
)::MergerEpochs
    isempty(snaps) && throw(ArgumentError("merger_epochs needs at least one snapshot"))
    n_cl = length(cluster_ranges)
    n_cl ≥ 1 || throw(ArgumentError("merger_epochs needs at least one cluster"))
    n_t = length(snaps)
    length(bound_sets) == n_t || throw(
        DimensionMismatch("merger_epochs: $(length(bound_sets)) bound sets for $n_t snapshots"),
    )
    t = [time_nb(s.header) for s in snaps]
    t_myr = [time_myr(s.header) for s in snaps]
    at_nb(k) = k > 0 ? t[k] : NaN
    at_myr(k) = k > 0 ? t_myr[k] : NaN
    n_cl == 1 && return _single_body_epochs(snaps, cluster_ranges, bound_sets, settings, t, t_myr)

    co =
        coalescence === nothing ?
        coalescence_time(snaps, cluster_ranges; bound_sets, max_passes, min_members) : coalescence
    size(co.present) == (n_cl, n_t) || throw(
        DimensionMismatch(
            "merger_epochs: coalescence result of size $(size(co.present)) for $n_cl clusters and $n_t snapshots",
        ),
    )
    tree = _tree_events(co)
    s = _epoch_tag_series(snaps, cluster_ranges, bound_sets, settings; max_passes, min_members)
    pairs = co.pairs
    n_pairs = length(pairs)
    slot = Dict{Tuple{Int,Int},Int}(p => q for (q, p) in enumerate(pairs))
    pair_of(i, j) = slot[(min(i, j), max(i, j))]

    # Pair series and detectors
    roche = fill(NaN, n_pairs, n_t)
    δs = fill(NaN, n_pairs, n_t)
    pair_contact = zeros(Int, n_pairs)
    d_self = fill(NaN, n_t)
    for (q, (i, j)) in enumerate(pairs)
        for k in 1:n_t
            d_self[k] = _column_distance(s.self_centre, i, j, k)
            m_i, m_j = s.self_bound_mass[i, k], s.self_bound_mass[j, k]
            if m_i > 0 && m_j > 0
                roche[q, k] =
                    _roche_fill(s.self_radius[i, k], s.self_radius[j, k], m_i, m_j, d_self[k])
            end
            (isfinite(s.core_sigma_x[i, k]) && isfinite(s.core_sigma_x[j, k])) || continue
            h = s.bound_mass[i, k] ≥ s.bound_mass[j, k] ? i : j
            Δx = _column_distance(s.core_pos, i, j, k)
            Δv = _column_distance(s.core_vel, i, j, k)
            δs[q, k] = Δx / s.core_sigma_x[h, k] + Δv / s.core_sigma_v[h, k]
        end
        pair_contact[q] = contact_index(
            view(s.self_radius, i, :),
            view(s.self_radius, j, :),
            view(s.self_bound_mass, i, :),
            view(s.self_bound_mass, j, :),
            d_self,
        )
    end

    # Merger-tree events
    n_ev = length(tree.events)
    group_a = [e.group_a for e in tree.events]
    group_b = [e.group_b for e in tree.events]
    event_join = [e.index for e in tree.events]
    event_contact = zeros(Int, n_ev)
    event_core = zeros(Int, n_ev)
    event_settled = zeros(Int, n_ev)
    event_tcr = fill(NaN, n_ev)
    main_progenitor = zeros(Int, n_ev)
    mass_ratio = fill(NaN, n_ev)
    pair_event = zeros(Int, n_pairs)
    for (e, ev) in enumerate(tree.events)
        first_contact = typemax(Int)
        for a in ev.group_a, b in ev.group_b
            q = pair_of(a, b)
            pair_event[q] = e
            pair_contact[q] > 0 && (first_contact = min(first_contact, pair_contact[q]))
        end
        event_contact[e] = first_contact == typemax(Int) ? 0 : first_contact
        k_ref = max(1, (event_contact[e] > 0 ? event_contact[e] : ev.index) - 1)
        main_a = ev.group_a[argmax([s.bound_mass[a, k_ref] for a in ev.group_a])]
        main_b = ev.group_b[argmax([s.bound_mass[b, k_ref] for b in ev.group_b])]
        M_a = sum(s.bound_mass[a, k_ref] for a in ev.group_a)
        M_b = sum(s.bound_mass[b, k_ref] for b in ev.group_b)
        main_progenitor[e] = M_a ≥ M_b ? main_a : main_b
        mass_ratio[e] = max(M_a, M_b) > 0 ? min(M_a, M_b) / max(M_a, M_b) : NaN
        c = core_index(view(δs, pair_of(main_a, main_b), :); threshold = settings.merger_core_delta)
        event_core[e] = c
        c > 0 || continue
        event_tcr[e] = _bound_crossing_time(snaps[c], bound_sets[c])
        window =
            settings.merger_coalescing_tcr == 0 ? 0.0 :
            (isfinite(event_tcr[e]) ? settings.merger_coalescing_tcr * event_tcr[e] : Inf)
        k_s = findfirst(k -> k ≥ c && t[k] ≥ t[c] + window, 1:n_t)
        event_settled[e] = k_s === nothing ? 0 : k_s
    end
    pair_core = [e == 0 ? 0 : event_core[e] for e in pair_event]

    # Progenitor epochs
    events_of = [[e for e in 1:n_ev if i in group_a[e] || i in group_b[e]] for i in 1:n_cl]
    tag_contact = zeros(Int, n_cl)
    tag_core = zeros(Int, n_cl)
    for i in 1:n_cl
        contacts = [pair_contact[pair_of(i, j)] for j in 1:n_cl if j != i]
        tag_contact[i] = minimum(filter(>(0), contacts); init = typemax(Int))
        tag_contact[i] == typemax(Int) && (tag_contact[i] = 0)
        cores = [event_core[e] for e in events_of[i] if event_core[e] > 0]
        tag_core[i] = isempty(cores) ? 0 : minimum(cores)
    end
    tag_dissolved = [
        s.m_initial[i] > 0 ?
        dissolution_index(
            view(s.self_bound_mass, i, :),
            s.m_initial[i];
            fraction = settings.merger_dissolved_fraction,
        ) : 0 for i in 1:n_cl
    ]
    # Latest event core at or before k among the events `es`: (core, settled), (0, 0) if none
    function latest(es, k)
        best = 0
        for e in es
            c = event_core[e]
            (c > 0 && c ≤ k) || continue
            (best == 0 || c ≥ event_core[best]) && (best = e)
        end
        return best == 0 ? (0, 0) : (event_core[best], event_settled[best])
    end
    coalesced(j, k) = tag_core[j] > 0 && k ≥ tag_core[j]
    tag_regime = Matrix{Symbol}(undef, n_cl, n_t)
    for i in 1:n_cl, k in 1:n_t
        core_k, settled_k = latest(events_of[i], k)
        satellite =
            core_k == 0 &&
            any(j -> j != i && 0 < pair_contact[pair_of(i, j)] ≤ k && coalesced(j, k), 1:n_cl)
        tag_regime[i, k] = regime_label(
            k,
            tag_contact[i],
            core_k,
            settled_k;
            dissolved = tag_dissolved[i],
            satellite,
        )
    end

    # Pair regimes
    pair_regime = Matrix{Symbol}(undef, n_pairs, n_t)
    for (q, (i, j)) in enumerate(pairs), k in 1:n_t
        e = pair_event[q]
        if e > 0 && event_core[e] > 0 && k ≥ event_core[e]
            pair_regime[q, k] = regime_label(k, pair_contact[q], event_core[e], event_settled[e])
        else
            ri, rj = tag_regime[i, k], tag_regime[j, k]
            in_contact = 0 < pair_contact[q] ≤ k
            pair_regime[q, k] = regime_label(
                k,
                pair_contact[q],
                0,
                0;
                dissolved = (ri === :dissolved || rj === :dissolved) ? k : 0,
                satellite = in_contact && (ri === :satellite || rj === :satellite),
            )
        end
    end

    # System epochs
    sys_contact = minimum(filter(>(0), pair_contact); init = typemax(Int))
    sys_contact == typemax(Int) && (sys_contact = 0)
    one_body(k) =
        any(i -> tag_regime[i, k] in _REMNANT_REGIMES, 1:n_cl) && all(
            i ->
                s.n_bound[i, k] < settings.merger_min_members ||
                tag_regime[i, k] in (:coalescing, :settled, :satellite, :dissolved),
            1:n_cl,
        )
    sys_core = 0
    for k in n_t:-1:1
        one_body(k) || break
        sys_core = k
    end
    all_events = collect(1:n_ev)
    regime = Vector{Symbol}(undef, n_t)
    for k in 1:n_t
        core_k, settled_k = (sys_core > 0 && k ≥ sys_core) ? latest(all_events, k) : (0, 0)
        regime[k] = regime_label(k, sys_contact, core_k, settled_k)
    end
    settled_index = 0
    if sys_core > 0
        k_s = findfirst(k -> k ≥ sys_core && regime[k] === :settled, 1:n_t)
        settled_index = k_s === nothing ? 0 : k_s
    end
    t_cr_core = sys_core > 0 ? _bound_crossing_time(snaps[sys_core], bound_sets[sys_core]) : NaN
    # Orbital angular momentum of the cores at the last snapshot before they coalesce: at
    # `core_index` itself the cores coincide and only a noise residual is left.
    l_orb = sys_core > 0 ? _core_orbital_angular_momentum(s, max(sys_core - 1, 1)) : fill(NaN, 3)
    k_coal = co.index === nothing ? 0 : co.index

    return MergerEpochs(
        t,
        t_myr,
        sys_contact,
        at_nb(sys_contact),
        at_myr(sys_contact),
        sys_core,
        at_nb(sys_core),
        at_myr(sys_core),
        settled_index,
        k_coal,
        at_nb(k_coal),
        at_myr(k_coal),
        t_cr_core,
        l_orb,
        regime,
        pairs,
        co.separation,
        co.threshold,
        pair_contact,
        tree.joined,
        pair_core,
        roche,
        δs,
        pair_regime,
        group_a,
        group_b,
        event_join,
        event_contact,
        event_core,
        event_settled,
        main_progenitor,
        mass_ratio,
        tag_contact,
        tag_core,
        tag_dissolved,
        s.bound_mass,
        s.self_bound_mass,
        tag_regime,
        settings.merger_min_members,
        settings.merger_contact_lagrangian,
        settings.merger_core_n,
        settings.merger_core_delta,
        settings.merger_dissolved_fraction,
        settings.merger_coalescing_tcr,
    )
end

"""
    _core_orbital_angular_momentum(series, k) -> Vector{Float64}

Orbital angular momentum `Σ Mᵢ (xᵢ − X) × (vᵢ − V)` of the most-bound cores
of `_epoch_tag_series` at snapshot `k`, `xᵢ`, `vᵢ` the core centroid
and mean velocity, `Mᵢ` the tag's mass in the global bound set and `X`, `V`
their mass-weighted means; `NaN` with fewer than two cores.
"""
function _core_orbital_angular_momentum(series, k::Int)::Vector{Float64}
    n_cl = size(series.bound_mass, 1)
    sel = [i for i in 1:n_cl if isfinite(series.core_sigma_x[i, k]) && series.bound_mass[i, k] > 0]
    length(sel) ≥ 2 || return fill(NaN, 3)
    M = [series.bound_mass[i, k] for i in sel]
    Mt = sum(M)
    X = zeros(3)
    V = zeros(3)
    for (w, i) in zip(M, sel), a in 1:3
        X[a] += w * series.core_pos[a, i, k] / Mt
        V[a] += w * series.core_vel[a, i, k] / Mt
    end
    L = zeros(3)
    for (w, i) in zip(M, sel)
        r = series.core_pos[:, i, k] .- X
        v = series.core_vel[:, i, k] .- V
        L .+= w .* (r × v)
    end
    return L
end

"""
    _single_body_epochs(snaps, cluster_ranges, bound_sets, settings, t, t_myr) -> MergerEpochs

Epochs of a run with one progenitor: one body throughout, every regime
`settled`, no epoch reached, the bound mass of the tag at every snapshot.
"""
function _single_body_epochs(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    bound_sets::AbstractVector{BoundSet},
    settings::PostprocessConfig,
    t::Vector{Float64},
    t_myr::Vector{Float64},
)::MergerEpochs
    n_t = length(snaps)
    bound_mass = zeros(1, n_t)
    for (k, snap) in enumerate(snaps)
        in_bound = falses(nparticles(snap))
        in_bound[bound_sets[k].stars] .= true
        bound_mass[1, k] = sum(
            Float64(snap.mass[s]) for s in _member_indices(snap, cluster_ranges[1]) if in_bound[s];
            init = 0.0,
        )
    end
    return MergerEpochs(
        t,
        t_myr,
        0,
        NaN,
        NaN,
        0,
        NaN,
        NaN,
        0,
        0,
        NaN,
        NaN,
        NaN,
        fill(NaN, 3),
        fill(:settled, n_t),
        Tuple{Int,Int}[],
        zeros(0, n_t),
        zeros(0, n_t),
        Int[],
        Int[],
        Int[],
        zeros(0, n_t),
        zeros(0, n_t),
        Matrix{Symbol}(undef, 0, n_t),
        Vector{Int}[],
        Vector{Int}[],
        Int[],
        Int[],
        Int[],
        Int[],
        Int[],
        Float64[],
        [0],
        [0],
        [0],
        bound_mass,
        fill(NaN, 1, n_t),
        fill(:settled, 1, n_t),
        settings.merger_min_members,
        settings.merger_contact_lagrangian,
        settings.merger_core_n,
        settings.merger_core_delta,
        settings.merger_dissolved_fraction,
        settings.merger_coalescing_tcr,
    )
end

# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

"""
    write_merger_epochs(path, epochs::MergerEpochs) -> String

Write the epochs as CSV. A leading comment line records the thresholds and
the crossing time of the `coalescing` window
(`# n_clusters=… n_snapshots=… merger_min_members=… … t_cr_core_nb=…`); one
row per entity follows, with the columns

    kind,group_a,group_b,contact_index,t_contact_nb,t_contact_myr,core_index,t_core_nb,
    t_core_myr,settled_index,t_settled_nb,t_settled_myr,join_index,t_join_nb,t_join_myr,
    dissolved_index,t_dissolved_nb,t_dissolved_myr,main_progenitor,mass_ratio

`kind` is `system` (`group_a` all tags; `join` the centre coalescence of
[`coalescence_time`](@ref)), `event` (one per merger-tree join), `pair`
(`group_a`, `group_b` the two tags; `core` and `settled` those of the event
joining their groups) or `tag` (`group_a` the tag; `core` its first core
coalescence; `dissolved` its dissolution). Groups are `;`-joined tags, indices
`0` for "never", non-finite values and fields that do not apply empty cells.
An existing file is backed up, never overwritten. Returns the path.
"""
function write_merger_epochs(path::AbstractString, epochs::MergerEpochs)::String
    _backup_existing(path)
    t, tm = epochs.time, epochs.time_myr
    nb(k) = k > 0 ? _csv_precise(t[k]) : ""
    my(k) = k > 0 ? _csv_precise(tm[k]) : ""
    epoch(k) = (string(k), nb(k), my(k))
    blank = ("", "", "")
    n_cl = size(epochs.tag_regime, 1)
    open(path, "w") do io
        println(
            io,
            "# n_clusters=$(n_cl) n_snapshots=$(length(t)) merger_min_members=$(epochs.min_members)" *
            " merger_contact_lagrangian=$(epochs.contact_lagrangian) merger_core_n=$(epochs.core_n)" *
            " merger_core_delta=$(epochs.core_delta) merger_dissolved_fraction=$(epochs.dissolved_fraction)" *
            " merger_coalescing_tcr=$(epochs.coalescing_tcr) t_cr_core_nb=$(_csv_precise(epochs.t_cr_core))",
        )
        println(
            io,
            "kind,group_a,group_b,contact_index,t_contact_nb,t_contact_myr,core_index,t_core_nb," *
            "t_core_myr,settled_index,t_settled_nb,t_settled_myr,join_index,t_join_nb,t_join_myr," *
            "dissolved_index,t_dissolved_nb,t_dissolved_myr,main_progenitor,mass_ratio",
        )
        row(cells...) = println(io, join(cells, ","))
        row(
            "system",
            join(1:n_cl, ";"),
            "",
            epoch(epochs.contact_index)...,
            epoch(epochs.core_index)...,
            epoch(epochs.settled_index)...,
            epoch(epochs.coalescence_index)...,
            blank...,
            "",
            "",
        )
        for e in eachindex(epochs.event_join)
            row(
                "event",
                join(epochs.group_a[e], ";"),
                join(epochs.group_b[e], ";"),
                epoch(epochs.event_contact[e])...,
                epoch(epochs.event_core[e])...,
                epoch(epochs.event_settled[e])...,
                epoch(epochs.event_join[e])...,
                blank...,
                string(epochs.main_progenitor[e]),
                _csv_precise(epochs.mass_ratio[e]),
            )
        end
        for (q, (i, j)) in enumerate(epochs.pairs)
            e = findfirst(
                x ->
                    (i in epochs.group_a[x] && j in epochs.group_b[x]) ||
                    (j in epochs.group_a[x] && i in epochs.group_b[x]),
                eachindex(epochs.event_join),
            )
            settled = e === nothing ? 0 : epochs.event_settled[e]
            row(
                "pair",
                string(i),
                string(j),
                epoch(epochs.pair_contact[q])...,
                epoch(epochs.pair_core[q])...,
                epoch(settled)...,
                epoch(epochs.pair_join[q])...,
                blank...,
                "",
                "",
            )
        end
        if !isempty(epochs.pairs)
            for i in 1:n_cl
                row(
                    "tag",
                    string(i),
                    "",
                    epoch(epochs.tag_contact[i])...,
                    epoch(epochs.tag_core[i])...,
                    blank...,
                    blank...,
                    epoch(epochs.tag_dissolved[i])...,
                    "",
                    "",
                )
            end
        end
    end
    return String(path)
end

"""
    write_pair_separation(path, epochs::MergerEpochs) -> String

Write the separation series of every cluster pair as CSV, one row per
(snapshot, pair), snapshots outer: `time_nb,time_myr,pair,separation,threshold,centre_offset,regime`.
`pair` is `i;j`; `separation` and `threshold` are the density-centre distance
and the merge threshold `min(r_h,i, r_h,j)` of [`coalescence_time`](@ref)
while the pair is `detached`, `interacting`, `dissolved` or a `satellite`, and
empty from its core coalescence on, where the same distance is written as
`centre_offset`, the offset of two sub-centres inside one remnant; `regime`
the pair's regime. Non-finite values are empty cells. An existing file is
backed up, never overwritten. Returns the path.
"""
function write_pair_separation(path::AbstractString, epochs::MergerEpochs)::String
    _backup_existing(path)
    open(path, "w") do io
        println(io, "time_nb,time_myr,pair,separation,threshold,centre_offset,regime")
        for k in eachindex(epochs.time), (q, (i, j)) in enumerate(epochs.pairs)
            merged = epochs.pair_regime[q, k] in _REMNANT_REGIMES
            d = epochs.separation[q, k]
            cells = (
                _csv_precise(epochs.time[k]),
                _csv_precise(epochs.time_myr[k]),
                "$i;$j",
                merged ? "" : _csv_precise(d),
                merged ? "" : _csv_precise(epochs.threshold[q, k]),
                merged ? _csv_precise(d) : "",
                String(epochs.pair_regime[q, k]),
            )
            println(io, join(cells, ","))
        end
    end
    return String(path)
end
