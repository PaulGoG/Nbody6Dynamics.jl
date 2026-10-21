# =============================================================================
# Bound spans: physically bound pairs of the snapshots joined across snapshots
# =============================================================================
#
# Level P of the binary census. At every snapshot the mutually most-bound pairs
# (one level of the node joining of Bate 2009, DOI 10.1111/j.1365-2966.2009.14165.x)
# are kept when the tidal perturbation of their nearest neighbours at apocentre
# is weak, or when the dominant perturber is the outer body of a triple that is
# stable by the criterion of Mardling & Aarseth (2001, DOI
# 10.1046/j.1365-8711.2001.03974.x). The engine's KS pairs are members
# regardless, since single-precision snapshots cannot resolve hard pairs. The
# memberships of a pair of names are then joined over the snapshot sequence.

# ---------------------------------------------------------------------------
# Record types
# ---------------------------------------------------------------------------

"""
    BoundPair

A physically bound pair of one snapshot from [`bound_pairs`](@ref), in the
N-body units of the snapshot (G = 1).

| field | meaning | unit |
|---|---|---|
| `name1`, `name2` | component names, `name1 < name2` | — |
| `a`, `e` | semi-major axis and eccentricity of the relative orbit; `NaN` for an unbound KS pair | NB, — |
| `eps` | specific two-body energy ``ε = \\frac{1}{2}\\lvert Δ\\mathbf{v}\\rvert^2 - (m_1 + m_2)/r`` | NB |
| `gamma_apo` | apocentre perturbation ``γ_\\mathrm{apo}``; `NaN` for a KS pair kept only as such | — |
| `hierarchy` | kept as the inner pair of a stable triple despite ``γ_\\mathrm{apo} ≥ γ_\\mathrm{max}`` | — |
| `outer_name` | name of the outer body of that triple, 0 when none | — |
| `ks` | the pair is a KS pair of the snapshot | — |
| `m1`, `m2` | masses of `name1` and `name2` | NB |
"""
struct BoundPair
    name1::Int
    name2::Int
    a::Float64
    e::Float64
    eps::Float64
    gamma_apo::Float64
    hierarchy::Bool
    outer_name::Int
    ks::Bool
    m1::Float64
    m2::Float64
end

"""
    BoundSpan

One physically bound pair of names over a snapshot sequence (level P) from
[`bound_spans`](@ref): the member snapshots `j_first … j_last`, with up to
`n_miss` consecutive excused misses bridged. Times and elements are in N-body
units.

| field | meaning | unit |
|---|---|---|
| `name1`, `name2` | component names, `name1 < name2` | — |
| `primordial` | the pair is (2k − 1, 2k) with 2k ≤ 2·NBIN0 | — |
| `j_first`, `j_last` | indices of the first and last member snapshot | — |
| `t_form_lo`, `t_form_hi` | formation interval ``(t_{j_\\mathrm{first} - 1}, t_{j_\\mathrm{first}}]``; `t_form_lo` is `NaN` when left-censored | NB |
| `t_end_lo`, `t_end_hi` | end interval ``[t_{j_\\mathrm{last}}, t_{j_\\mathrm{last} + 1})``; `t_end_hi` is `NaN` when right-censored | NB |
| `left_censored` | `j_first` is the first snapshot | — |
| `right_censored` | `j_last` is the last snapshot | — |
| `n_snapshots` | member snapshots | — |
| `n_missed` | bridged (excused) misses | — |
| `ks_fraction` | member snapshots in which the pair is a KS pair, over `n_snapshots` | — |
| `hierarchy_flag` | the pair is the inner pair of a stable triple in some member snapshot | — |
| `outer_name` | last nonzero outer body of such a triple, 0 when none | — |
| `a_first`, `e_first` | elements at `j_first` | NB, — |
| `a_last`, `e_last` | elements at `j_last` | NB, — |
| `m1`, `m2` | masses of `name1` and `name2` at `j_first` | NB |
| `end_unbound` | not right-censored, and at `j_last + 1` both components are present, unbound and separated by more than `k_sep · a_last` | — |
"""
Base.@kwdef struct BoundSpan
    name1::Int
    name2::Int
    primordial::Bool
    j_first::Int
    j_last::Int
    t_form_lo::Float64
    t_form_hi::Float64
    t_end_lo::Float64
    t_end_hi::Float64
    left_censored::Bool
    right_censored::Bool
    n_snapshots::Int
    n_missed::Int
    ks_fraction::Float64
    hierarchy_flag::Bool
    outer_name::Int
    a_first::Float64
    e_first::Float64
    a_last::Float64
    e_last::Float64
    m1::Float64
    m2::Float64
    end_unbound::Bool
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Separation, specific energy, semi-major axis and eccentricity of the relative
# orbit of stars i and j (G = 1); ε is NaN for coincident stars, a and e NaN
# unless the orbit is bound
function _bound_orbit(
    pos::AbstractMatrix{<:Real},
    vel::AbstractMatrix{<:Real},
    mass::AbstractVector{<:Real},
    i::Int,
    j::Int,
)::NTuple{4,Float64}
    @inbounds begin
        rx = Float64(pos[1, i]) - Float64(pos[1, j])
        ry = Float64(pos[2, i]) - Float64(pos[2, j])
        rz = Float64(pos[3, i]) - Float64(pos[3, j])
        vx = Float64(vel[1, i]) - Float64(vel[1, j])
        vy = Float64(vel[2, i]) - Float64(vel[2, j])
        vz = Float64(vel[3, i]) - Float64(vel[3, j])
        μ = Float64(mass[i]) + Float64(mass[j])
    end
    r = sqrt(rx * rx + ry * ry + rz * rz)
    r > 0 || return (r, NaN, NaN, NaN)
    ε = 0.5 * (vx * vx + vy * vy + vz * vz) - μ / r
    ε < 0 || return (r, ε, NaN, NaN)
    a = -μ / (2ε)
    hx = ry * vz - rz * vy
    hy = rz * vx - rx * vz
    hz = rx * vy - ry * vx
    e = sqrt(max(1 - (hx * hx + hy * hy + hz * hz) / (a * μ), 0.0))
    return (r, ε, a, e)
end

# γ_apo of the pair (i, j) with apocentre distance r_apo and c.m. `centre`, and
# the index of the body of its largest term (0 when there is none)
function _bound_perturbation(
    tree::KDTree,
    mass::Vector{Float64},
    i::Int,
    j::Int,
    centre::Vector{Float64},
    m_b::Float64,
    r_apo::Float64,
    n_perturbers::Int,
)::Tuple{Float64,Int}
    n_perturbers == 0 && return (0.0, 0)
    nbr, dist = knn(tree, centre, min(n_perturbers + 2, length(mass)), true)
    γ = 0.0
    dominant = 0
    term_max = -Inf
    taken = 0
    for p in eachindex(nbr)
        taken == n_perturbers && break
        l = Int(nbr[p])
        (l == i || l == j) && continue
        taken += 1
        d = Float64(dist[p])
        term = d > 0 ? 2 * (mass[l] / m_b) * (r_apo / d)^3 : Inf
        γ += term
        if term > term_max
            term_max = term
            dominant = l
        end
    end
    return (γ, dominant)
end

# Whether body l is bound to the c.m. of the pair (i, j) and the triple
# satisfies the Mardling & Aarseth (2001) stability limit
function _bound_stable_triple(
    pos::Matrix{Float64},
    vel::Matrix{Float64},
    mass::Vector{Float64},
    i::Int,
    j::Int,
    l::Int,
    a_in::Float64,
    centre::Vector{Float64},
    v_b::Vector{Float64},
    m_b::Float64,
)::Bool
    m_l = mass[l]
    μ = m_b + m_l
    @inbounds begin
        rx = pos[1, l] - centre[1]
        ry = pos[2, l] - centre[2]
        rz = pos[3, l] - centre[3]
        vx = vel[1, l] - v_b[1]
        vy = vel[2, l] - v_b[2]
        vz = vel[3, l] - v_b[3]
        sx = pos[1, i] - pos[1, j]
        sy = pos[2, i] - pos[2, j]
        sz = pos[3, i] - pos[3, j]
        ux = vel[1, i] - vel[1, j]
        uy = vel[2, i] - vel[2, j]
        uz = vel[3, i] - vel[3, j]
    end
    d = sqrt(rx * rx + ry * ry + rz * rz)
    d > 0 || return false
    ε₃ = 0.5 * (vx * vx + vy * vy + vz * vz) - μ / d
    ε₃ < 0 || return false
    a_out = -μ / (2ε₃)
    # Outer and inner specific angular momenta
    hx = ry * vz - rz * vy
    hy = rz * vx - rx * vz
    hz = rx * vy - ry * vx
    kx = sy * uz - sz * uy
    ky = sz * ux - sx * uz
    kz = sx * uy - sy * ux
    h_out = sqrt(hx * hx + hy * hy + hz * hz)
    h_in = sqrt(kx * kx + ky * ky + kz * kz)
    e_out = sqrt(max(1 - h_out^2 / (a_out * μ), 0.0))
    e_out < 1 || return false
    # Inclination [deg]; an undefined angle (radial orbit) counts as coplanar
    h_prod = h_out * h_in
    incl = h_prod > 0 ? acosd(clamp((hx * kx + hy * ky + hz * kz) / h_prod, -1.0, 1.0)) : 0.0
    q = m_l / m_b
    limit = 2.8 * (1 + q)^0.4 * (1 + e_out)^0.4 * (1 - e_out)^(-0.2) * (1 - 0.3 * incl / 180)
    return a_out * (1 - e_out) / a_in > limit
end

# A BoundPair of stars i, j of `snap`, names ordered
function _bound_pair(
    snap::Snapshot,
    mass::Vector{Float64},
    i::Int,
    j::Int,
    a::Float64,
    e::Float64,
    ε::Float64,
    γ::Float64,
    hierarchy::Bool,
    outer_name::Int,
    ks::Bool,
)::BoundPair
    n_i = Int(snap.name[i])
    n_j = Int(snap.name[j])
    first_i = n_i < n_j
    return BoundPair(
        first_i ? n_i : n_j,
        first_i ? n_j : n_i,
        a,
        e,
        ε,
        γ,
        hierarchy,
        outer_name,
        ks,
        first_i ? mass[i] : mass[j],
        first_i ? mass[j] : mass[i],
    )
end

# A missed snapshot is excused when both components are present and bound to
# each other and neither is in another pair of the snapshot
function _bound_excused(
    key::Tuple{Int,Int},
    snap::Snapshot,
    index::Dict{Int,Int},
    paired::Set{Int},
)::Bool
    (key[1] in paired || key[2] in paired) && return false
    i = get(index, key[1], 0)
    j = get(index, key[2], 0)
    (i == 0 || j == 0) && return false
    ε = _bound_orbit(snap.pos, snap.vel, snap.mass, i, j)[2]
    return ε < 0
end

# The span of `key` over the member snapshots js
function _bound_span(
    key::Tuple{Int,Int},
    js::AbstractVector{Int},
    missed::Int,
    members::Vector{Dict{Tuple{Int,Int},BoundPair}},
    snaps::Vector{Snapshot},
    index::Vector{Dict{Int,Int}},
    t::Vector{Float64},
    k_sep::Float64,
    nbin0::Int,
)::BoundSpan
    ns = length(snaps)
    j1 = first(js)
    j2 = last(js)
    q1 = members[j1][key]
    q2 = members[j2][key]
    n_ks = 0
    hierarchy = false
    outer = 0
    for j in js
        q = members[j][key]
        n_ks += q.ks
        hierarchy |= q.hierarchy
        q.outer_name != 0 && (outer = q.outer_name)
    end
    left = j1 == 1
    right = j2 == ns
    end_unbound = false
    if !right && isfinite(q2.a)
        s = snaps[j2 + 1]
        i = get(index[j2 + 1], key[1], 0)
        l = get(index[j2 + 1], key[2], 0)
        if i != 0 && l != 0
            r, ε, _, _ = _bound_orbit(s.pos, s.vel, s.mass, i, l)
            end_unbound = ε > 0 && r > k_sep * q2.a
        end
    end
    n = length(js)
    return BoundSpan(;
        name1 = key[1],
        name2 = key[2],
        primordial = _span_primordial(key, nbin0),
        j_first = j1,
        j_last = j2,
        t_form_lo = left ? NaN : t[j1 - 1],
        t_form_hi = t[j1],
        t_end_lo = t[j2],
        t_end_hi = right ? NaN : t[j2 + 1],
        left_censored = left,
        right_censored = right,
        n_snapshots = n,
        n_missed = missed,
        ks_fraction = n_ks / n,
        hierarchy_flag = hierarchy,
        outer_name = outer,
        a_first = q1.a,
        e_first = q1.e,
        a_last = q2.a,
        e_last = q2.e,
        m1 = q1.m1,
        m2 = q1.m2,
        end_unbound = end_unbound,
    )
end

# Spans of one pair of names from its ascending member snapshots js
function _bound_walk!(
    spans::Vector{BoundSpan},
    key::Tuple{Int,Int},
    js::Vector{Int},
    members::Vector{Dict{Tuple{Int,Int},BoundPair}},
    paired::Vector{Set{Int}},
    index::Vector{Dict{Int,Int}},
    snaps::Vector{Snapshot},
    t::Vector{Float64},
    n_miss::Int,
    k_sep::Float64,
    nbin0::Int,
)::Vector{BoundSpan}
    p = 1
    while p <= length(js)
        p_first = p
        j_last = js[p]
        missed = 0
        p += 1
        while p <= length(js)
            j_next = js[p]
            gap = j_next - j_last - 1
            gap <= n_miss || break
            all(
                m -> _bound_excused(key, snaps[m], index[m], paired[m]),
                (j_last + 1):(j_next - 1),
            ) || break
            missed += gap
            j_last = j_next
            p += 1
        end
        push!(
            spans,
            _bound_span(
                key,
                view(js, p_first:(p - 1)),
                missed,
                members,
                snaps,
                index,
                t,
                k_sep,
                nbin0,
            ),
        )
    end
    return spans
end

# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

"""
    bound_pairs(snap::Snapshot; k::Integer = 10, n_perturbers::Integer = 4,
                gamma_max::Real = 0.1) -> Vector{BoundPair}

The physically bound pairs of a snapshot, in its N-body units (G = 1, arrays
converted to `Float64`), sorted by `(name1, name2)`.

# Arguments
- `snap`: the snapshot; its first `2·npairs` stars are the KS components
- `k`: neighbour rank of the candidate search of [`snapshot_pairs`](@ref)
- `n_perturbers`: nearest other bodies to a pair's centre of mass entering ``γ_\\mathrm{apo}``
- `gamma_max`: largest accepted ``γ_\\mathrm{apo}``

# Method
1. Candidates are the bound, mutually most-bound pairs of
   [`snapshot_pairs`](@ref), the first level of the node joining of Bate
   (2009, DOI 10.1111/j.1365-2966.2009.14165.x).
2. With ``m_b = m_1 + m_2`` and ``d_l`` the distances of the `n_perturbers`
   nearest other bodies from the pair's centre of mass, the apocentre
   perturbation is
   ``γ_\\mathrm{apo} = \\sum_l 2 (m_l / m_b) [a(1 + e) / d_l]^3``. A candidate
   is accepted when ``γ_\\mathrm{apo} <`` `gamma_max`. Otherwise it is kept
   with `hierarchy = true` when the body ``l^*`` of the largest term is bound
   to the pair's centre of mass,
   ``ε_3 = \\frac{1}{2}\\lvert\\mathbf{v}_l - \\mathbf{v}_b\\rvert^2 - (m_b + m_l)/d < 0``,
   and the triple satisfies the stability limit of Mardling & Aarseth (2001,
   DOI 10.1046/j.1365-8711.2001.03974.x, eq. 90) in the form
   ``R_{p,\\mathrm{out}} / a > 2.8 (1 + q)^{2/5} (1 + e_\\mathrm{out})^{2/5} (1 - e_\\mathrm{out})^{-1/5} (1 - 0.3\\, i / 180°)``,
   with ``q = m_l / m_b``, ``e_\\mathrm{out}`` and
   ``R_{p,\\mathrm{out}} = a_\\mathrm{out}(1 - e_\\mathrm{out})`` of the orbit of
   ``l^*`` about the pair's centre of mass and ``i`` [deg] the angle between
   the inner and outer orbital angular momenta (0 when one of them vanishes);
   it is rejected otherwise.
3. Every KS pair of the snapshot is a member whatever steps 1 and 2 give, with
   `ks = true`; a KS pair not accepted in step 2 carries the elements of its
   relative state when bound (`NaN` when not), `gamma_apo = NaN` and
   `hierarchy = false`. A candidate that joins a KS component to another star
   is dropped.

Throws `ArgumentError` for `n_perturbers < 0` or `gamma_max ≤ 0`, and as
[`snapshot_pairs`](@ref) for `k < 1`.

# Example
```julia
pairs = bound_pairs(read_conf3(path))
count(p -> p.hierarchy, pairs)
```
"""
function bound_pairs(
    snap::Snapshot;
    k::Integer = 10,
    n_perturbers::Integer = 4,
    gamma_max::Real = 0.1,
)::Vector{BoundPair}
    n_perturbers >= 0 ||
        throw(ArgumentError("bound_pairs: n_perturbers must be >= 0; got $n_perturbers"))
    gamma_max > 0 || throw(ArgumentError("bound_pairs: gamma_max must be > 0; got $gamma_max"))
    pos = Float64.(snap.pos)
    vel = Float64.(snap.vel)
    mass = Float64.(snap.mass)
    n = length(mass)
    n_p = Int(n_perturbers)
    γ_max = Float64(gamma_max)
    partner = zeros(Int, n)
    for p in 1:(snap.npairs)
        partner[2p - 1] = 2p
        partner[2p] = 2p - 1
    end
    ks_found = falses(snap.npairs)
    pairs = BoundPair[]
    if n >= 2
        candidates = snapshot_pairs(pos, vel, mass; k = Int(k))
        tree = KDTree(pos)
        centre = zeros(3)
        v_b = zeros(3)
        for c in candidates
            i, j = c.i, c.j
            ks = partner[i] == j
            # A KS component joined to a third star: the KS pair wins
            !ks && (partner[i] != 0 || partner[j] != 0) && continue
            _, ε, a, e = _bound_orbit(pos, vel, mass, i, j)
            isfinite(a) || continue
            m_b = mass[i] + mass[j]
            @inbounds for d in 1:3
                centre[d] = (mass[i] * pos[d, i] + mass[j] * pos[d, j]) / m_b
                v_b[d] = (mass[i] * vel[d, i] + mass[j] * vel[d, j]) / m_b
            end
            γ, l = _bound_perturbation(tree, mass, i, j, centre, m_b, a * (1 + e), n_p)
            hierarchy = false
            outer_name = 0
            if !(γ < γ_max)
                (l != 0 && _bound_stable_triple(pos, vel, mass, i, j, l, a, centre, v_b, m_b)) ||
                    continue
                hierarchy = true
                outer_name = Int(snap.name[l])
            end
            ks && (ks_found[max(i, j) ÷ 2] = true)
            push!(pairs, _bound_pair(snap, mass, i, j, a, e, ε, γ, hierarchy, outer_name, ks))
        end
    end
    for p in 1:(snap.npairs)
        ks_found[p] && continue
        i, j = 2p - 1, 2p
        _, ε, a, e = _bound_orbit(pos, vel, mass, i, j)
        push!(pairs, _bound_pair(snap, mass, i, j, a, e, ε, NaN, false, 0, true))
    end
    sort!(pairs; by = q -> (q.name1, q.name2))
    return pairs
end

"""
    bound_spans(snapshots::AbstractVector{Snapshot}; nbin0::Integer, k::Integer = 10,
                n_perturbers::Integer = 4, gamma_max::Real = 0.1,
                n_miss::Integer = 1, k_sep::Real = 20.0) -> Vector{BoundSpan}

The physically bound pairs (level P) of a snapshot sequence joined to
[`BoundSpan`](@ref)s, in N-body units, sorted by `(t_form_hi, name1, name2)`.

# Arguments
- `snapshots`: the snapshots in increasing time
- `nbin0`: initial number of primordial binaries NBIN0
- `k`, `n_perturbers`, `gamma_max`: membership at a snapshot, see [`bound_pairs`](@ref)
- `n_miss`: largest number of consecutive excused misses bridged within a span
- `k_sep`: separation, in units of the last semi-major axis, beyond which an
  unbound pair counts as separated (`end_unbound`)

# Method
The members of snapshot j are the pairs of [`bound_pairs`](@ref), at time
``t_j`` = [`time_nb`](@ref) of its header. For each pair of names, a span runs
from a member snapshot over the following ones; a snapshot in which the pair is
not a member is an excused miss when both components are present, their
two-body energy ε is negative and neither is a component of another member
pair. Up to `n_miss` consecutive excused misses are bridged when a member
snapshot follows them; misses at the end of the sequence are not. Otherwise the
span ends at its last member snapshot ``j_2``, and the next membership of the
pair opens a new span.

Formation and end are interval-censored: ``t_\\mathrm{form} ∈ (t_{j_1 - 1}, t_{j_1}]``,
``t_\\mathrm{end} ∈ [t_{j_2}, t_{j_2 + 1})``. A span is left-censored when
``j_1`` is the first snapshot (`t_form_lo = NaN`) and right-censored when
``j_2`` is the last (`t_end_hi = NaN`). `end_unbound` records that at
``j_2 + 1`` both components are present with ε > 0 and a separation above
`k_sep · a_last`, false for a right-censored span or a `NaN` `a_last`.
A pair is primordial when its names are (2k − 1, 2k) with 2k ≤ 2·NBIN0.

Throws `ArgumentError` for `nbin0 < 0`, `n_miss < 0`, `k_sep ≤ 0`, an empty
`snapshots` or snapshot times that do not increase strictly.

# Example
```julia
spans = bound_spans(read_all_conf3(dir); nbin0 = 2_000)
write_bound_spans(joinpath(dir, "bound_spans.csv"), spans)
```
"""
function bound_spans(
    snapshots::AbstractVector{Snapshot};
    nbin0::Integer,
    k::Integer = 10,
    n_perturbers::Integer = 4,
    gamma_max::Real = 0.1,
    n_miss::Integer = 1,
    k_sep::Real = 20.0,
)::Vector{BoundSpan}
    nbin0 >= 0 || throw(ArgumentError("bound_spans: nbin0 must be >= 0; got $nbin0"))
    n_miss >= 0 || throw(ArgumentError("bound_spans: n_miss must be >= 0; got $n_miss"))
    k_sep > 0 || throw(ArgumentError("bound_spans: k_sep must be > 0; got $k_sep"))
    isempty(snapshots) && throw(ArgumentError("bound_spans: no snapshots"))
    snaps = Vector{Snapshot}(snapshots)
    ns = length(snaps)
    t = Float64[time_nb(s.header) for s in snaps]
    for j in 2:ns
        t[j] > t[j - 1] || throw(
            ArgumentError(
                "bound_spans: snapshot times must increase; t[$j] = $(t[j]) after $(t[j - 1])",
            ),
        )
    end

    members = Vector{Dict{Tuple{Int,Int},BoundPair}}(undef, ns)
    paired = Vector{Set{Int}}(undef, ns)
    index = Vector{Dict{Int,Int}}(undef, ns)
    memberships = Dict{Tuple{Int,Int},Vector{Int}}()
    for j in 1:ns
        s = snaps[j]
        d = Dict{Tuple{Int,Int},BoundPair}()
        names = Set{Int}()
        for q in bound_pairs(s; k, n_perturbers, gamma_max)
            key = (q.name1, q.name2)
            d[key] = q
            push!(names, q.name1, q.name2)
            push!(get!(Vector{Int}, memberships, key), j)
        end
        members[j] = d
        paired[j] = names
        index[j] = Dict{Int,Int}(Int(s.name[i]) => i for i in eachindex(s.name))
    end

    spans = BoundSpan[]
    for key in sort!(collect(keys(memberships)))
        _bound_walk!(
            spans,
            key,
            memberships[key],
            members,
            paired,
            index,
            snaps,
            t,
            Int(n_miss),
            Float64(k_sep),
            Int(nbin0),
        )
    end
    sort!(spans; by = s -> (s.t_form_hi, s.name1, s.name2))
    return spans
end

"""
    write_bound_spans(path, spans::AbstractVector{BoundSpan}) -> String

Write the spans of [`bound_spans`](@ref) as CSV, one column per field of
[`BoundSpan`](@ref) in declaration order (N-body units, non-finite values as
empty cells). An existing file is backed up, never overwritten. Returns the
path.
"""
function write_bound_spans(path::AbstractString, spans::AbstractVector{BoundSpan})::String
    _backup_existing(path)
    n_fields = fieldcount(BoundSpan)
    open(path, "w") do io
        println(io, join(fieldnames(BoundSpan), ","))
        for s in spans
            println(io, join(ntuple(f -> _csv_precise(getfield(s, f)), n_fields), ","))
        end
    end
    return String(path)
end
