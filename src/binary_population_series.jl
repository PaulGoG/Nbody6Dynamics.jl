# =============================================================================
# Binary population series (design W8): number, star and mass binary fractions
# per radial shell, stellar class and census level of a snapshot sequence
# =============================================================================

"""
    BinaryPopulationRow

Binary fractions of one snapshot, census level, radial bin and stellar class
from [`binary_population_series`](@ref). Radii are in the N-body units of the
snapshot, masses in M☉.

| field | meaning | unit |
|---|---|---|
| `snapshot` | index of the snapshot in the input vector | — |
| `t_nb`, `t_myr` | snapshot time | NB, Myr |
| `level` | census level: `:ks` (E, KS pairs) or `:bound` (P, bound pairs) | — |
| `radius_frame` | `:lagrangian` (enclosed-mass bins) or `:rh` (bins in r/r_h) | — |
| `bin` | 0 for all systems, else the shell index 1, 2, … | — |
| `edge_lo`, `edge_hi` | shell edges: mass fractions (`:lagrangian`) or multiples of r_h (`:rh`) | — |
| `r_lo_nb`, `r_hi_nb` | shell radii; `NaN` for an empty Lagrangian shell | NB |
| `class` | `:all` or the `key` of a [`STELLAR_CLASSES`](@ref) entry | — |
| `n_singles` | single stars (of the class) | — |
| `n_binaries` | binaries (with a primary of the class) | — |
| `n_in_binaries` | binary components (of the class) | — |
| `n_hierarchies` | binaries of `n_binaries` in a hierarchy: at `:bound` the inner pair of a stable triple, at `:ks` a KS pair with a merged inner binary (negative `pair_name`) | — |
| `m_singles_msun` | mass of the singles of `n_singles` | M☉ |
| `m_binaries_msun` | mass of the binaries (`:all`) or of the components of the class in binaries | M☉ |
| `f_b_number` | ``f_b = N_b / (N_s + N_b)`` | — |
| `f_b_stars` | ``f_⋆ = N_\\mathrm{in} / (N_s + N_\\mathrm{in})`` | — |
| `f_b_mass` | ``f_M = M_b / (M_s + M_b)`` | — |
| `f_b_number_se` | binomial standard error ``\\sqrt{f_b (1 - f_b) / (N_s + N_b)}`` | — |
| `f_b_q` | fraction of binaries with ``q > q_\\mathrm{min}`` among the systems with single or primary mass in the window | — |
| `r_h_binaries_nb`, `r_h_singles_nb` | half-mass radii of the binaries and singles of the row; `NaN` on shells | NB |

Fractions with a vanishing denominator are `NaN`.
"""
struct BinaryPopulationRow
    snapshot::Int
    t_nb::Float64
    t_myr::Float64
    level::Symbol
    radius_frame::Symbol
    bin::Int
    edge_lo::Float64
    edge_hi::Float64
    r_lo_nb::Float64
    r_hi_nb::Float64
    class::Symbol
    n_singles::Int
    n_binaries::Int
    n_in_binaries::Int
    n_hierarchies::Int
    m_singles_msun::Float64
    m_binaries_msun::Float64
    f_b_number::Float64
    f_b_stars::Float64
    f_b_mass::Float64
    f_b_number_se::Float64
    f_b_q::Float64
    r_h_binaries_nb::Float64
    r_h_singles_nb::Float64
end

const _BPOP_LEVELS = (:ks, :bound)
const _BPOP_FRAMES = (:lagrangian, :rh)
const _BPOP_EDGES_LAGRANGIAN = (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
const _BPOP_EDGES_RH = (0.0, 0.5, 1.0, 2.0, 4.0, Inf)
# Relative tolerance of an epoch match in N-body time: the sev.83/bev.82 data
# lines carry six significant digits, conf.3 a Float32.
const _BPOP_EPOCH_RTOL = 1.0e-5

# Systems of a snapshot at one census level: the pairs first, then the singles.
# `secondary` is 0 for a single, whose `primary` is the star itself.
struct _BpopSystems
    mass::Vector{Float64}
    radius::Vector{Float64}
    primary::Vector{Int}
    secondary::Vector{Int}
    hierarchy::Vector{Bool}
end

# Counts and masses per (bin + 1, class column); column 1 is `:all`.
struct _BpopTally
    n_singles::Matrix{Int}
    n_binaries::Matrix{Int}
    n_in_binaries::Matrix{Int}
    n_hierarchies::Matrix{Int}
    m_singles::Matrix{Float64}
    m_binaries::Matrix{Float64}
    n_window::Matrix{Int}
    n_window_q::Matrix{Int}
end

_BpopTally(nr::Int, nc::Int) = _BpopTally(
    zeros(Int, nr, nc),
    zeros(Int, nr, nc),
    zeros(Int, nr, nc),
    zeros(Int, nr, nc),
    zeros(nr, nc),
    zeros(nr, nc),
    zeros(Int, nr, nc),
    zeros(Int, nr, nc),
)

"""
    binary_population_series(snapshots::AbstractVector{Snapshot}; levels = (:ks, :bound),
                             radius_frame::Symbol = :lagrangian, edges = nothing,
                             centre = :shrinking_sphere, centre_kwargs::NamedTuple = (;),
                             kstar = nothing, q_min::Real = 0.5,
                             primary_window_msun::Tuple{<:Real,<:Real} = (0.5, 0.9),
                             k::Integer = 10, n_perturbers::Integer = 4,
                             gamma_max::Real = 0.1, max_untyped_fraction::Real = 0.01)
        -> Vector{BinaryPopulationRow}

Binary fractions of a snapshot sequence per census level, radial bin and
stellar class, one [`BinaryPopulationRow`](@ref) per snapshot × level × bin ×
class, in that order.

# Arguments
- `snapshots`: the snapshots
- `levels`: census levels, any of `:ks` (level E) and `:bound` (level P)
- `radius_frame`: `:lagrangian` (shells in enclosed mass fraction) or `:rh`
  (shells in units of the half-mass radius)
- `edges`: shell edges, strictly increasing from 0; mass fractions ending at 1
  for `:lagrangian` (default `[0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]`), multiples
  of r_h for `:rh` (default `[0, 0.5, 1, 2, 4, Inf]`)
- `centre`: `:shrinking_sphere` ([`shrinking_sphere_centre`](@ref) of the
  stars of each snapshot, with `centre_kwargs`), `:origin`, or a fixed point
  (length-3 vector, NB)
- `kstar`: `nothing` (class `:all` only) or a vector aligned with `snapshots`
  of name → K* maps, e.g. from [`stellar_type_map`](@ref). A star without a
  type counts under `:all` only: the engine lists neither the outer body of a
  merged hierarchy (a KS pair with a negative `pair_name`) nor the body that
  represents a chain subsystem in sev.83 or bev.82
- `q_min`, `primary_window_msun`: mass-ratio threshold and inclusive primary
  mass window [M☉] of `f_b_q`
- `k`, `n_perturbers`, `gamma_max`: pair membership at level P, see
  [`bound_pairs`](@ref)
- `max_untyped_fraction`: largest fraction of the stars of a snapshot without a
  stellar type; more indicates mismatched epochs and raises an `ArgumentError`

Massless bodies (the ghosts of merged and chain subsystems) are left out of
every count.

# Method
At a level the stars of a snapshot form systems. Level `:ks` takes the KS
pairs (stars 2k − 1 and 2k), level `:bound` the pairs of [`bound_pairs`](@ref);
the inner pair of a stable triple is a binary and its outer body a single.
Every other star is a single. A binary has the mass m₁ + m₂ and the position
of its centre of mass; its primary is the heavier component and
q = m₂/m₁ ≤ 1. Masses are converted to M☉ with [`zmbar`](@ref).

Radii are measured from the centre of each snapshot, computed once for all
levels. The half-mass radius r_h of a set of systems is the radius of the
first system, in radius order, at which the enclosed mass reaches half the
total. With `:rh` shell b holds the systems with
``e_b r_h ≤ r < e_{b+1} r_h``; with `:lagrangian` those whose enclosed mass
fraction lies in ``(e_b, e_{b+1}]``. Bin 0 holds every system.

With ``N_s`` singles and ``N_b`` binaries of total masses ``M_s`` and ``M_b``,
the number fraction is ``f_b = N_b/(N_s + N_b)``, with binomial standard error
``\\sqrt{f_b(1 - f_b)/(N_s + N_b)}``, the fraction of stars in binaries
``f_⋆ = 2N_b/(N_s + 2N_b) = 2f_b/(1 + f_b)`` and the mass fraction
``f_M = M_b/(M_s + M_b)``. A class row counts the singles of the class and the
binaries whose primary belongs to it; its `n_in_binaries` and
`m_binaries_msun` count the components of the class in any binary, so that its
``f_⋆`` is the fraction of the members of the class bound in pairs, whatever
their partner. `f_b_q` is the fraction of binaries with q > `q_min` among the
systems whose single or primary mass lies in `primary_window_msun`; on
`:main_sequence` rows it is the analogue of the q-restricted fractions
measured photometrically in globular clusters (Milone et al. 2012,
DOI 10.1051/0004-6361/201016384).

A fraction is meaningful only together with the period range of its level:
level E is restricted to the pairs the integrator regularised, level P extends
to every pair bound against its nearest perturbers, and observed fractions
refer to their own period windows (Hurley et al. 2005,
DOI 10.1111/j.1365-2966.2005.09448.x). Photometric fractions of unresolved
systems, which include chance superpositions, are not computed.

Throws `ArgumentError` for an empty `snapshots`, empty or unknown `levels`, an
unknown `radius_frame`, invalid `edges`, an unknown `centre`, a `kstar` not
aligned with `snapshots` or lacking a star, `q_min` outside (0, 1), a window
with lo ≥ hi or lo < 0, and for a star in two pairs.

# Example
```julia
snaps = read_all_conf3(dir)
rows = binary_population_series(snaps; levels = (:ks,), radius_frame = :rh)
write_binary_population(joinpath(dir, "binary_population.csv"), rows)
```
"""
function binary_population_series(
    snapshots::AbstractVector{Snapshot};
    levels = (:ks, :bound),
    radius_frame::Symbol = :lagrangian,
    edges = nothing,
    centre = :shrinking_sphere,
    centre_kwargs::NamedTuple = (;),
    kstar = nothing,
    q_min::Real = 0.5,
    primary_window_msun::Tuple{<:Real,<:Real} = (0.5, 0.9),
    k::Integer = 10,
    n_perturbers::Integer = 4,
    gamma_max::Real = 0.1,
    max_untyped_fraction::Real = 0.01,
)::Vector{BinaryPopulationRow}
    isempty(snapshots) && throw(ArgumentError("binary_population_series: no snapshots"))
    levs = _bpop_levels(levels)
    radius_frame in _BPOP_FRAMES || throw(
        ArgumentError(
            "binary_population_series: radius_frame must be :lagrangian or :rh; got :$radius_frame",
        ),
    )
    shell_edges = _bpop_edges(edges, radius_frame)
    c_fixed = _bpop_fixed_centre(centre)
    kmaps = _bpop_kstar_maps(kstar, length(snapshots))
    0 < q_min < 1 ||
        throw(ArgumentError("binary_population_series: q_min must lie in (0, 1); got $q_min"))
    window = (Float64(primary_window_msun[1]), Float64(primary_window_msun[2]))
    (window[1] ≥ 0 && window[1] < window[2]) || throw(
        ArgumentError(
            "binary_population_series: primary_window_msun must satisfy 0 ≤ lo < hi; got $primary_window_msun",
        ),
    )
    q_thr = Float64(q_min)
    0 ≤ max_untyped_fraction ≤ 1 || throw(
        ArgumentError(
            "binary_population_series: max_untyped_fraction must lie in [0, 1]; got $max_untyped_fraction",
        ),
    )
    untyped = Float64(max_untyped_fraction)

    rows = BinaryPopulationRow[]
    for (j, snap) in enumerate(snapshots)
        mass = Float64.(snap.mass)
        pos = Float64.(snap.pos)
        c =
            c_fixed === nothing ? shrinking_sphere_centre(pos, mass; centre_kwargs...).centre :
            c_fixed
        kcls = kmaps === nothing ? Int[] : _bpop_classes(snap, kmaps[j], j, untyped)
        t_nb = time_nb(snap.header)
        t_myr = time_myr(snap.header)
        m_unit = zmbar(snap.header)
        for level in levs
            comps, hierarchy = _bpop_components(snap, level, k, n_perturbers, gamma_max)
            sys = _bpop_systems(snap, mass, pos, c, comps, hierarchy, j)
            shells = _bpop_shells(sys.radius, sys.mass, shell_edges, radius_frame)
            _bpop_rows!(
                rows,
                sys,
                shells,
                mass,
                kcls,
                j,
                t_nb,
                t_myr,
                level,
                radius_frame,
                m_unit,
                window,
                q_thr,
            )
        end
    end
    return rows
end

function _bpop_levels(levels)::Vector{Symbol}
    levs = Symbol[]
    for l in levels
        (l isa Symbol && l in _BPOP_LEVELS) || throw(
            ArgumentError(
                "binary_population_series: levels must be drawn from (:ks, :bound); got $(repr(l))",
            ),
        )
        push!(levs, l)
    end
    isempty(levs) && throw(ArgumentError("binary_population_series: levels is empty"))
    return levs
end

function _bpop_edges(edges, frame::Symbol)::Vector{Float64}
    e = if edges === nothing
        collect(Float64, frame === :rh ? _BPOP_EDGES_RH : _BPOP_EDGES_LAGRANGIAN)
    else
        collect(Float64, edges)
    end
    length(e) ≥ 2 ||
        throw(ArgumentError("binary_population_series: edges needs at least two entries; got $e"))
    e[1] == 0 || throw(ArgumentError("binary_population_series: edges must start at 0; got $e"))
    all(i -> e[i + 1] > e[i], 1:(length(e) - 1)) ||
        throw(ArgumentError("binary_population_series: edges must increase strictly; got $e"))
    frame === :lagrangian &&
        e[end] != 1 &&
        throw(
            ArgumentError(
                "binary_population_series: Lagrangian edges are mass fractions ending at 1; got $e",
            ),
        )
    return e
end

# The centre shared by all snapshots, or `nothing` for a shrinking-sphere centre
# computed per snapshot.
function _bpop_fixed_centre(centre)::Union{Nothing,Vector{Float64}}
    centre === :origin && return zeros(3)
    centre === :shrinking_sphere && return nothing
    (centre isa AbstractVector{<:Real} && length(centre) == 3) && return Vector{Float64}(centre)
    throw(
        ArgumentError(
            "binary_population_series: centre must be :origin, :shrinking_sphere or a length-3 vector; got $(repr(centre))",
        ),
    )
end

function _bpop_kstar_maps(kstar, n::Int)
    kstar === nothing && return nothing
    (kstar isa AbstractVector && length(kstar) == n) || throw(
        ArgumentError(
            "binary_population_series: kstar must be a vector of $n name → K* maps, one per snapshot",
        ),
    )
    return collect(kstar)
end

# Class index of every star, 0 for a star without a stellar type. The engine
# types neither the outer body of a merged hierarchy nor the representative of
# a chain subsystem in sev.83/bev.82; such stars count under `:all` only. More
# untyped stars than `max_untyped` (a fraction of the stars) means mismatched
# epochs and is an error.
function _bpop_classes(snap::Snapshot, kmap, j::Int, max_untyped::Float64)::Vector{Int}
    kmap isa AbstractDict ||
        throw(ArgumentError("binary_population_series: kstar[$j] must be a name → K* dictionary"))
    cls = Vector{Int}(undef, length(snap.name))
    n_untyped = 0
    first_untyped = 0
    for i in eachindex(snap.name)
        nm = Int(snap.name[i])
        kt = get(kmap, nm, nothing)
        if kt === nothing
            cls[i] = 0
            n_untyped += 1
            first_untyped == 0 && (first_untyped = nm)
            continue
        end
        cls[i] = stellar_class_index(kt)
    end
    n_untyped > max_untyped * length(cls) && throw(
        ArgumentError(
            "binary_population_series: $n_untyped of the $(length(cls)) stars of snapshot $j " *
            "(e.g. $first_untyped) have no stellar type in kstar, more than the fraction " *
            "max_untyped_fraction = $max_untyped",
        ),
    )
    return cls
end

function _bpop_components(
    snap::Snapshot,
    level::Symbol,
    k::Integer,
    n_perturbers::Integer,
    gamma_max::Real,
)::Tuple{Vector{NTuple{2,Int}},Vector{Bool}}
    if level === :ks
        comps = pair_indices(snap)
        # Merged hierarchies: KS pairs with a negative centre-of-mass name
        merged = [!isempty(snap.pair_name) && snap.pair_name[k] < 0 for k in eachindex(comps)]
        return comps, merged
    end
    pairs = bound_pairs(snap; k, n_perturbers, gamma_max)
    index = Dict{Int,Int}(Int(snap.name[i]) => i for i in eachindex(snap.name))
    comps = NTuple{2,Int}[(index[p.name1], index[p.name2]) for p in pairs]
    return comps, Bool[p.hierarchy for p in pairs]
end

function _bpop_systems(
    snap::Snapshot,
    mass::Vector{Float64},
    pos::Matrix{Float64},
    c::Vector{Float64},
    comps::Vector{NTuple{2,Int}},
    hierarchy::Vector{Bool},
    j::Int,
)::_BpopSystems
    n = length(mass)
    paired = falses(n)
    for (a, b) in comps
        for i in (a, b)
            paired[i] && throw(
                ArgumentError(
                    "binary_population_series: star $(snap.name[i]) of snapshot $j is a component of two pairs",
                ),
            )
            paired[i] = true
        end
    end
    n_pairs = length(comps)
    # Massless singles are ghosts of merged or chain subsystems, not systems
    n_ghosts = count(i -> !paired[i] && mass[i] == 0, 1:n)
    n_sys = n - n_pairs - n_ghosts
    m = Vector{Float64}(undef, n_sys)
    r = Vector{Float64}(undef, n_sys)
    primary = Vector{Int}(undef, n_sys)
    secondary = zeros(Int, n_sys)
    hier = fill(false, n_sys)
    @inbounds for (s, (a, b)) in enumerate(comps)
        m_s = mass[a] + mass[b]
        d2 = 0.0
        for d in 1:3
            x = (mass[a] * pos[d, a] + mass[b] * pos[d, b]) / m_s
            d2 += (x - c[d])^2
        end
        m[s] = m_s
        r[s] = sqrt(d2)
        primary[s], secondary[s] = mass[b] > mass[a] ? (b, a) : (a, b)
        hier[s] = hierarchy[s]
    end
    s = n_pairs
    @inbounds for i in 1:n
        (paired[i] || mass[i] == 0) && continue
        s += 1
        d2 = 0.0
        for d in 1:3
            d2 += (pos[d, i] - c[d])^2
        end
        m[s] = mass[i]
        r[s] = sqrt(d2)
        primary[s] = i
    end
    return _BpopSystems(m, r, primary, secondary, hier)
end

function _bpop_half_mass_radius(radius::Vector{Float64}, mass::Vector{Float64})::Float64
    isempty(radius) && return NaN
    order = sortperm(radius)
    cum = cumsum(mass[order]) ./ sum(mass)
    i = findfirst(≥(0.5), cum)
    return i === nothing ? NaN : radius[order[i]]
end

_bpop_scaled(e::Float64, r_h::Float64)::Float64 = isinf(e) ? e : e * r_h

# Shell of every system (0 when outside all shells) and the edges and radii of
# bin 0 and the shells, indexed bin + 1.
function _bpop_shells(
    radius::Vector{Float64},
    mass::Vector{Float64},
    edges::Vector{Float64},
    frame::Symbol,
)
    n = length(radius)
    nb = length(edges) - 1
    shell = zeros(Int, n)
    edge_lo = [0.0; edges[1:nb]]
    edge_hi = [frame === :lagrangian ? 1.0 : Inf; edges[2:end]]
    r_lo = zeros(nb + 1)
    r_hi = fill(NaN, nb + 1)
    r_hi[1] = n == 0 ? NaN : maximum(radius)
    if frame === :rh
        r_h = _bpop_half_mass_radius(radius, mass)
        for b in 1:nb
            r_lo[b + 1] = _bpop_scaled(edges[b], r_h)
            r_hi[b + 1] = _bpop_scaled(edges[b + 1], r_h)
        end
        for s in 1:n
            r = radius[s]
            for b in 1:nb
                if r_lo[b + 1] ≤ r < r_hi[b + 1]
                    shell[s] = b
                    break
                end
            end
        end
    else
        order = sortperm(radius)
        cum = cumsum(mass[order]) ./ sum(mass)
        # e_b < cum ≤ e_{b+1}; round-off beyond the last edge joins the last shell
        for (p, s) in enumerate(order)
            b = clamp(searchsortedfirst(edges, cum[p]) - 1, 1, nb)
            shell[s] = b
            r_hi[b + 1] = radius[s]
        end
        prev = 0.0
        for b in 1:nb
            r_lo[b + 1] = prev
            isnan(r_hi[b + 1]) || (prev = r_hi[b + 1])
        end
    end
    return (; shell, edge_lo, edge_hi, r_lo, r_hi)
end

function _bpop_add_single!(t::_BpopTally, row::Int, col::Int, m::Float64, in_window::Bool)
    t.n_singles[row, col] += 1
    t.m_singles[row, col] += m
    in_window && (t.n_window[row, col] += 1)
    return nothing
end

function _bpop_add_binary!(
    t::_BpopTally,
    row::Int,
    col::Int,
    hierarchy::Bool,
    in_window::Bool,
    q_above::Bool,
)
    t.n_binaries[row, col] += 1
    hierarchy && (t.n_hierarchies[row, col] += 1)
    if in_window
        t.n_window[row, col] += 1
        q_above && (t.n_window_q[row, col] += 1)
    end
    return nothing
end

function _bpop_tally!(
    t::_BpopTally,
    row::Int,
    sys::_BpopSystems,
    s::Int,
    mass::Vector{Float64},
    kcls::Vector{Int},
    colof::Vector{Int},
    m_unit::Float64,
    window::NTuple{2,Float64},
    q_min::Float64,
)
    p = sys.primary[s]
    sec = sys.secondary[s]
    in_window = window[1] ≤ m_unit * mass[p] ≤ window[2]
    if sec == 0
        _bpop_add_single!(t, row, 1, mass[p], in_window)
        col = _bpop_column(colof, kcls, p)
        col > 0 && _bpop_add_single!(t, row, col, mass[p], in_window)
        return nothing
    end
    q_above = mass[sec] / mass[p] > q_min
    _bpop_add_binary!(t, row, 1, sys.hierarchy[s], in_window, q_above)
    t.n_in_binaries[row, 1] += 2
    t.m_binaries[row, 1] += sys.mass[s]
    col = _bpop_column(colof, kcls, p)
    col > 0 && _bpop_add_binary!(t, row, col, sys.hierarchy[s], in_window, q_above)
    for i in (p, sec)
        col = _bpop_column(colof, kcls, i)
        col > 0 || continue
        t.n_in_binaries[row, col] += 1
        t.m_binaries[row, col] += mass[i]
    end
    return nothing
end

# Class column of star `i`; 0 without classes or for an unclassified star
_bpop_column(colof::Vector{Int}, kcls::Vector{Int}, i::Int)::Int =
    (isempty(kcls) || kcls[i] == 0) ? 0 : colof[kcls[i]]

_bpop_ratio(a::Real, d::Real)::Float64 = d == 0 ? NaN : Float64(a) / Float64(d)

function _bpop_rows!(
    rows::Vector{BinaryPopulationRow},
    sys::_BpopSystems,
    shells,
    mass::Vector{Float64},
    kcls::Vector{Int},
    j::Int,
    t_nb::Float64,
    t_myr::Float64,
    level::Symbol,
    frame::Symbol,
    m_unit::Float64,
    window::NTuple{2,Float64},
    q_min::Float64,
)
    n_cls = length(STELLAR_CLASSES)
    cols = [0]
    if !isempty(kcls)
        present = falses(n_cls)
        for c in kcls
            c > 0 && (present[c] = true)
        end
        append!(cols, findall(present))
    end
    colof = zeros(Int, n_cls)
    for (q, c) in enumerate(cols)
        c > 0 && (colof[c] = q)
    end
    n_rows = length(shells.edge_lo)
    n_cols = length(cols)
    t = _BpopTally(n_rows, n_cols)
    singles = [Int[] for _ in 1:n_cols]
    binaries = [Int[] for _ in 1:n_cols]
    for s in eachindex(sys.mass)
        b = shells.shell[s]
        _bpop_tally!(t, 1, sys, s, mass, kcls, colof, m_unit, window, q_min)
        b > 0 && _bpop_tally!(t, b + 1, sys, s, mass, kcls, colof, m_unit, window, q_min)
        members = sys.secondary[s] == 0 ? singles : binaries
        push!(members[1], s)
        col = _bpop_column(colof, kcls, sys.primary[s])
        col > 0 && push!(members[col], s)
    end
    r_h_bin = [_bpop_half_mass_radius(sys.radius[l], sys.mass[l]) for l in binaries]
    r_h_sin = [_bpop_half_mass_radius(sys.radius[l], sys.mass[l]) for l in singles]
    for row in 1:n_rows, q in 1:n_cols
        key = cols[q] == 0 ? :all : STELLAR_CLASSES[cols[q]].key
        n_s = t.n_singles[row, q]
        n_b = t.n_binaries[row, q]
        n_in = t.n_in_binaries[row, q]
        m_s = t.m_singles[row, q] * m_unit
        m_b = t.m_binaries[row, q] * m_unit
        f = _bpop_ratio(n_b, n_s + n_b)
        push!(
            rows,
            BinaryPopulationRow(
                j,
                t_nb,
                t_myr,
                level,
                frame,
                row - 1,
                shells.edge_lo[row],
                shells.edge_hi[row],
                shells.r_lo[row],
                shells.r_hi[row],
                key,
                n_s,
                n_b,
                n_in,
                t.n_hierarchies[row, q],
                m_s,
                m_b,
                f,
                _bpop_ratio(n_in, n_s + n_in),
                _bpop_ratio(m_b, m_s + m_b),
                n_s + n_b == 0 ? NaN : sqrt(f * (1 - f) / (n_s + n_b)),
                _bpop_ratio(t.n_window_q[row, q], t.n_window[row, q]),
                row == 1 ? r_h_bin[q] : NaN,
                row == 1 ? r_h_sin[q] : NaN,
            ),
        )
    end
    return rows
end

"""
    write_binary_population(path, rows::AbstractVector{BinaryPopulationRow}) -> String

Write the rows of [`binary_population_series`](@ref) as CSV, one column per
field of [`BinaryPopulationRow`](@ref) in declaration order (symbols as their
names, non-finite values as empty cells). An existing file is backed up, never
overwritten. Returns the path.
"""
function write_binary_population(
    path::AbstractString,
    rows::AbstractVector{BinaryPopulationRow},
)::String
    _backup_existing(path)
    n_fields = fieldcount(BinaryPopulationRow)
    open(path, "w") do io
        println(io, join(fieldnames(BinaryPopulationRow), ","))
        for r in rows
            println(io, join(ntuple(f -> _bpop_cell(getfield(r, f)), n_fields), ","))
        end
    end
    return String(path)
end

_bpop_cell(v::Symbol) = String(v)
_bpop_cell(v) = _csv_precise(v)

"""
    stellar_type_map(sev::StellarEvolutionSnapshot,
                     bev::Union{Nothing,BinaryEvolutionSnapshot}) -> Dict{Int,Int}

Stellar type K* of every star of one epoch, keyed by name: the single stars of
the sev.83 records, then both components of every bev.82 record, which take
precedence (sev.83 omits the members of KS pairs). The `kstar` input of
[`binary_population_series`](@ref).

# Example
```julia
kmap = stellar_type_map(read_stellar_evolution(sev_path), read_binary_evolution(bev_path))
```
"""
function stellar_type_map(
    sev::StellarEvolutionSnapshot,
    bev::Union{Nothing,BinaryEvolutionSnapshot},
)::Dict{Int,Int}
    kmap = Dict{Int,Int}()
    for r in sev.records
        kmap[Int(r.name)] = Int(r.stellar_type)
    end
    if bev !== nothing
        for r in bev.records
            kmap[Int(r.name1)] = Int(r.stellar_type1)
            kmap[Int(r.name2)] = Int(r.stellar_type2)
        end
    end
    return kmap
end

# Index of the epoch nearest to `t` [NB] within the matching tolerance, or
# `nothing`. Epochs are timed by their records' N-body time: the file headers
# carry the time in Myr with one decimal only. Epochs without records are skipped.
function _bpop_nearest_epoch(
    epochs::AbstractVector{<:Union{StellarEvolutionSnapshot,BinaryEvolutionSnapshot}},
    t::Float64,
)::Union{Nothing,Int}
    best = 0
    d_best = Inf
    for (i, e) in enumerate(epochs)
        isempty(e.records) && continue
        d = abs(e.records[1].time_nb - t)
        if d < d_best
            best = i
            d_best = d
        end
    end
    return (best > 0 && d_best ≤ _BPOP_EPOCH_RTOL * max(1.0, abs(t))) ? best : nothing
end

# Name → K* maps aligned with `snapshots` from the sev.83 and bev.82 epochs of
# the same time, or `nothing` (warned) when a snapshot has no sev.83 epoch.
function _bpop_pipeline_kstar(
    snapshots::AbstractVector{Snapshot},
    sevs::AbstractVector{StellarEvolutionSnapshot},
    bevs::AbstractVector{BinaryEvolutionSnapshot},
)::Union{Nothing,Vector{Dict{Int,Int}}}
    maps = Dict{Int,Int}[]
    for (j, snap) in enumerate(snapshots)
        t = time_nb(snap.header)
        i = _bpop_nearest_epoch(sevs, t)
        if i === nothing
            @warn "Binary population: snapshot $j (t = $t NB) has no sev.83 epoch at its time; " *
                  "stellar classes omitted"
            return nothing
        end
        ib = _bpop_nearest_epoch(bevs, t)
        push!(maps, stellar_type_map(sevs[i], ib === nothing ? nothing : bevs[ib]))
    end
    return maps
end
