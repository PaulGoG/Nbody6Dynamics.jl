# =============================================================================
# Collision catalogue: stellar collisions (coll.13) and coalescences (coal.24)
# of a multi-cluster run, tagged with the progenitor cluster of each component
# =============================================================================
#
# Product names follow the engine. MIX (mix.f) orders the colliding pair so
# that I1 is the more evolved star, writes coll.13 in that order and returns
# the order to CMBODY through its dummy arguments; CMBODY (cmbody.f) then gives
# the product, which occupies I1, the name of the heavier body (BODY, the
# masses of the record). COAL (coal.f) places the product in I1 and writes
# NAME(I1), NAME(I2) to coal.24 before it exchanges the two names when
# BODY0(I1) < BODY0(I2); BODY0 is not part of the record.

"""R☉ per pc, for the coal.24 separation in N-body units."""
const _COLLCAT_RSUN_PER_PC = 4.435e7

"""
    CollisionRow

One stellar collision (`coll.13`) or coalescence (`coal.24`) of the collision
catalogue, from [`collision_catalogue`](@ref). The component listed first is
the one whose name the product keeps (engine convention, see
[`collision_catalogue`](@ref)); `NaN` and `0` mark unavailable values.

| field | meaning | unit |
|---|---|---|
| `kind` | `:collision` (coll.13) or `:coalescence` (coal.24) | — |
| `channel` | `channel` of the [`CoalescenceRecord`](@ref); `:none` for coll.13 | — |
| `time_nb`, `time_myr` | time of the record | NB, Myr |
| `name1`, `name2` | component names; the product keeps `name1` | — |
| `kstar1`, `kstar2`, `kstar_product` | component types before the event, product type | — |
| `m1_msun`, `m2_msun` | component masses before the event (BODY·M*) | M☉ |
| `m_product_msun`, `m_lost_msun` | product mass and mass lost, as recorded | M☉ |
| `progenitor1`, `progenitor2` | index of the progenitor cluster of each component in `cluster_ranges`; `0` when in no range | — |
| `mixed` | both progenitors known and different | — |
| `primordial_pair` | the two names are a primordial `dat.10` pair (2k − 1, 2k), k ≤ `nbin0` | — |
| `generation` | 1 + the number of earlier rows whose product may carry `name1` or `name2` | — |
| `separation_rsun` | separation at the event; coal.24: R12 · R̄ in R☉ (`NaN` without R̄) | R☉ |
| `eccentricity`, `period_d` | orbit at the event, as recorded | —, d |
| `r_centre` | distance from the density centre: coll.13 in core radii (RI/RC), coal.24 in NB (RI) | RC or NB |
| `kicked_before` | a `BINARY KICK` line named `name1` or `name2` before `time_nb` | — |
"""
Base.@kwdef struct CollisionRow
    kind::Symbol
    channel::Symbol
    time_nb::Float64
    time_myr::Float64
    name1::Int
    name2::Int
    kstar1::Int
    kstar2::Int
    kstar_product::Int
    m1_msun::Float64
    m2_msun::Float64
    m_product_msun::Float64
    m_lost_msun::Float64
    progenitor1::Int
    progenitor2::Int
    mixed::Bool
    primordial_pair::Bool
    generation::Int
    separation_rsun::Float64
    eccentricity::Float64
    period_d::Float64
    r_centre::Float64
    kicked_before::Bool
end

"""
    progenitor_of(name::Integer, cluster_ranges) -> Int

Index of the progenitor cluster of the body `name`: the position of the first
member list of `cluster_ranges` (the `dat.10` names of each initial cluster,
as returned by [`parse_merger_summary`](@ref)) that contains `name`, `0` when
none does. A linear search per call; [`collision_catalogue`](@ref) builds a
name table once instead.

# Example
```julia
cluster_ranges = parse_merger_summary(joinpath(dir, "merger_summary.txt"))
progenitor_of(150, cluster_ranges)
```
"""
function progenitor_of(
    name::Integer,
    cluster_ranges::AbstractVector{<:AbstractVector{<:Integer}},
)::Int
    for (i, members) in enumerate(cluster_ranges)
        name in members && return i
    end
    return 0
end

"""Name → progenitor-cluster table of `cluster_ranges`; `ArgumentError` for a name in two clusters."""
function _collcat_progenitor_table(
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
)::Dict{Int,Int}
    table = Dict{Int,Int}()
    for (i, members) in enumerate(cluster_ranges), name in members
        j = get(table, name, 0)
        if j != 0 && j != i
            throw(
                ArgumentError(
                    "collision_catalogue: cluster_ranges overlap; name $name lies in clusters $j and $i",
                ),
            )
        end
        table[name] = i
    end
    return table
end

"""The record lists the product's name second (coll.13: the heavier body is the second one)."""
_collcat_product_second(c::CollisionRecord)::Bool = c.mass2_msun > c.mass1_msun
_collcat_product_second(::CoalescenceRecord)::Bool = false

"""Components `(name, kstar, mass)` of a record, the one whose name the product keeps first."""
function _collcat_components(c::Union{CollisionRecord,CoalescenceRecord})
    a = (c.name1, c.kstar1, c.mass1_msun)
    b = (c.name2, c.kstar2, c.mass2_msun)
    return _collcat_product_second(c) ? (b, a) : (a, b)
end

_collcat_kind(::CollisionRecord)::Symbol = :collision
_collcat_kind(::CoalescenceRecord)::Symbol = :coalescence

_collcat_channel(::CollisionRecord)::Symbol = :none
_collcat_channel(c::CoalescenceRecord)::Symbol = c.channel

_collcat_separation_rsun(c::CollisionRecord, ::Float64)::Float64 = c.separation_rsun
_collcat_separation_rsun(c::CoalescenceRecord, rbar_pc::Float64)::Float64 =
    c.separation_nb * rbar_pc * _COLLCAT_RSUN_PER_PC

_collcat_r_centre(c::CollisionRecord)::Float64 = c.distance_over_core_radius
_collcat_r_centre(c::CoalescenceRecord)::Float64 = c.distance_density_centre_nb

"""Earliest `BINARY KICK` time [NB] of every name in `mass_loss` (time from TTOT·T* when TIME[NB] is absent)."""
function _collcat_first_kicks(
    mass_loss::Union{Nothing,AbstractVector{MassLossEvent}},
    tstar_myr::Float64,
)::Dict{Int,Float64}
    first_kick = Dict{Int,Float64}()
    mass_loss === nothing && return first_kick
    for ev in mass_loss
        ev.kind === :binary_kick || continue
        t = isfinite(ev.time_nb) ? ev.time_nb : ev.time_myr / tstar_myr
        isfinite(t) || continue
        for name in (ev.name1, ev.name2)
            name == 0 && continue
            first_kick[name] = min(get(first_kick, name, Inf), t)
        end
    end
    return first_kick
end

"""One [`CollisionRow`](@ref) of a coll.13 or coal.24 record."""
function _collcat_row(
    c::Union{CollisionRecord,CoalescenceRecord},
    table::Dict{Int,Int},
    first_kick::Dict{Int,Float64},
    generation::Int;
    tstar_myr::Float64,
    rbar_pc::Float64,
    nbin0::Int,
)::CollisionRow
    (n1, k1, m1), (n2, k2, m2) = _collcat_components(c)
    p1 = get(table, n1, 0)
    p2 = get(table, n2, 0)
    t = c.time_nb
    return CollisionRow(;
        kind = _collcat_kind(c),
        channel = _collcat_channel(c),
        time_nb = t,
        time_myr = t * tstar_myr,
        name1 = n1,
        name2 = n2,
        kstar1 = k1,
        kstar2 = k2,
        kstar_product = c.kstar_product,
        m1_msun = m1,
        m2_msun = m2,
        m_product_msun = c.mass_product_msun,
        m_lost_msun = c.mass_lost_msun,
        progenitor1 = p1,
        progenitor2 = p2,
        mixed = p1 != 0 && p2 != 0 && p1 != p2,
        primordial_pair = _span_primordial(_span_key(n1, n2), nbin0),
        generation = generation,
        separation_rsun = _collcat_separation_rsun(c, rbar_pc),
        eccentricity = c.eccentricity,
        period_d = c.period_d,
        r_centre = _collcat_r_centre(c),
        kicked_before = get(first_kick, n1, Inf) < t || get(first_kick, n2, Inf) < t,
    )
end

"""
    collision_catalogue(collisions::AbstractVector{CollisionRecord},
                        coalescences::AbstractVector{CoalescenceRecord},
                        cluster_ranges; tstar_myr, rbar_pc = NaN, nbin0 = 0,
                        mass_loss = nothing) -> Vector{CollisionRow}

Collision catalogue of a run: one [`CollisionRow`](@ref) per stellar
collision of `coll.13` and per coalescence of `coal.24`, each component tagged
with its progenitor cluster, sorted by `(time_nb, name1, name2)`.

# Arguments
- `collisions`: coll.13 records, from [`read_collision_log`](@ref)
- `coalescences`: coal.24 records, from [`read_coalescence_log`](@ref)
- `cluster_ranges`: the `dat.10` names of each initial cluster, from [`parse_merger_summary`](@ref); empty for a single-cluster run (every progenitor `0`)
- `tstar_myr`: time unit T* [Myr]
- `rbar_pc`: length unit R̄ [pc], for the coal.24 separation; `NaN` leaves it `NaN`
- `nbin0`: number of primordial pairs, the first `2nbin0` names of `dat.10`
- `mass_loss`: the `mass_loss` records of [`read_binary_events`](@ref), for `kicked_before`; `nothing` sets it `false`

# Conventions
The product keeps the name of the component listed first in the row. For
coll.13 this is the engine rule of `mix.f` and `cmbody.f`: the product takes
the name of the heavier body (record masses, the record order on a tie), so
the components are swapped when `mass2_msun > mass1_msun`. For coal.24 the
record order is kept; `coal.f` gives the product the name of the progenitor
with the larger initial mass (BODY0), which the record does not carry, so
either name may continue. The absorbed body becomes a ghost and its name does
not recur, so `generation` counts a coal.24 row for both of its names; it
counts earlier rows in the sorted order. `kicked_before` uses the `:binary_kick`
lines ([`MassLossEvent`](@ref)), which name both members of the kicked KS
pair; a kick of the companion therefore flags the row as well. `r_centre` is
in core radii for coll.13 and in N-body units for coal.24, as recorded.

Throws an `ArgumentError` for `tstar_myr ≤ 0`, `nbin0 < 0` or a name that
lies in two clusters of `cluster_ranges`.

# Example
```julia
rows = collision_catalogue(read_collision_log(joinpath(dir, "coll.13")),
                           read_coalescence_log(joinpath(dir, "coal.24")),
                           parse_merger_summary(joinpath(dir, "merger_summary.txt"));
                           tstar_myr = 0.87, rbar_pc = 2.5, nbin0 = 500,
                           mass_loss = read_binary_events(log_path).mass_loss)
collision_summary(rows).n_mixed
```
"""
function collision_catalogue(
    collisions::AbstractVector{CollisionRecord},
    coalescences::AbstractVector{CoalescenceRecord},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    tstar_myr::Real,
    rbar_pc::Real = NaN,
    nbin0::Integer = 0,
    mass_loss::Union{Nothing,AbstractVector{MassLossEvent}} = nothing,
)::Vector{CollisionRow}
    tstar_myr > 0 ||
        throw(ArgumentError("collision_catalogue: tstar_myr must be > 0; got $tstar_myr"))
    nbin0 ≥ 0 || throw(ArgumentError("collision_catalogue: nbin0 must be ≥ 0; got $nbin0"))
    tstar = Float64(tstar_myr)
    rbar = Float64(rbar_pc)
    table = _collcat_progenitor_table(cluster_ranges)
    first_kick = _collcat_first_kicks(mass_loss, tstar)

    # (time, name1, name2, source: 1 coll.13 / 2 coal.24, index), names product first
    order = Tuple{Float64,Int,Int,Int,Int}[]
    for (i, c) in enumerate(collisions)
        (n1, _, _), (n2, _, _) = _collcat_components(c)
        push!(order, (c.time_nb, n1, n2, 1, i))
    end
    for (i, c) in enumerate(coalescences)
        (n1, _, _), (n2, _, _) = _collcat_components(c)
        push!(order, (c.time_nb, n1, n2, 2, i))
    end
    sort!(order)

    units = (; tstar_myr = tstar, rbar_pc = rbar, nbin0 = Int(nbin0))
    products = Dict{Int,Int}()   # product name → number of rows so far
    rows = CollisionRow[]
    sizehint!(rows, length(order))
    for (_, n1, n2, source, i) in order
        generation = 1 + get(products, n1, 0) + (n2 == n1 ? 0 : get(products, n2, 0))
        products[n1] = get(products, n1, 0) + 1
        # coal.24: either name may be the product's (coal.f, BODY0); the other is a ghost's
        (source == 2 && n2 != n1) && (products[n2] = get(products, n2, 0) + 1)
        row = if source == 1
            _collcat_row(collisions[i], table, first_kick, generation; units...)
        else
            _collcat_row(coalescences[i], table, first_kick, generation; units...)
        end
        push!(rows, row)
    end
    return rows
end

"""
    collision_summary(rows::AbstractVector{CollisionRow}) -> NamedTuple

Summary of a collision catalogue from [`collision_catalogue`](@ref):

- `n_collisions`, `n_coalescences`: rows of coll.13 and of coal.24
- `n_mixed`, `n_primordial_pairs`, `n_kicked_before`: rows with the flag set
- `m_lost_total_msun`: sum of the finite `m_lost_msun` [M☉]
- `t_first_myr`, `t_last_myr`: earliest and latest finite `time_myr`; `NaN` without one
- `by_progenitor_pair`: `Dict{Tuple{Int,Int},Int}`, rows per unordered progenitor pair `(min, max)` (`0` for an unknown progenitor)
- `by_pair_class`: `Dict{Symbol,Int}`, rows per [`pair_class`](@ref) of the two component types
- `max_generation`: highest `generation`; `0` for an empty catalogue

# Example
```julia
s = collision_summary(rows)
s.by_progenitor_pair[(1, 2)]   # collisions between members of clusters 1 and 2
```
"""
function collision_summary(rows::AbstractVector{CollisionRow})::NamedTuple
    by_progenitor_pair = Dict{Tuple{Int,Int},Int}()
    by_pair_class = Dict{Symbol,Int}()
    m_lost = 0.0
    t_first = Inf
    t_last = -Inf
    for r in rows
        key = (min(r.progenitor1, r.progenitor2), max(r.progenitor1, r.progenitor2))
        by_progenitor_pair[key] = get(by_progenitor_pair, key, 0) + 1
        class = pair_class(r.kstar1, r.kstar2)
        by_pair_class[class] = get(by_pair_class, class, 0) + 1
        isfinite(r.m_lost_msun) && (m_lost += r.m_lost_msun)
        if isfinite(r.time_myr)
            t_first = min(t_first, r.time_myr)
            t_last = max(t_last, r.time_myr)
        end
    end
    return (;
        n_collisions = count(r -> r.kind === :collision, rows),
        n_coalescences = count(r -> r.kind === :coalescence, rows),
        n_mixed = count(r -> r.mixed, rows),
        n_primordial_pairs = count(r -> r.primordial_pair, rows),
        n_kicked_before = count(r -> r.kicked_before, rows),
        m_lost_total_msun = m_lost,
        t_first_myr = isfinite(t_first) ? t_first : NaN,
        t_last_myr = isfinite(t_last) ? t_last : NaN,
        by_progenitor_pair = by_progenitor_pair,
        by_pair_class = by_pair_class,
        max_generation = maximum(r -> r.generation, rows; init = 0),
    )
end

"""CSV cell: floats with 12 significant digits and `NaN` as `NaN`, Bools as `true`/`false`, Symbols by name."""
_collcat_cell(v::AbstractFloat)::String = isnan(v) ? "NaN" : @sprintf("%.12g", v)
_collcat_cell(v)::String = string(v)

"""
    write_collision_catalogue(path, rows::AbstractVector{CollisionRow}) -> String

Write the collision catalogue as CSV: a header of the field names of
[`CollisionRow`](@ref) in declaration order, one line per row. Floats carry
12 significant digits (`NaN` written as `NaN`), Bools `true`/`false`, Symbols
their names. An existing file is backed up, never overwritten. Returns the
path.

# Example
```julia
write_collision_catalogue(joinpath(dir, "collision_catalogue.csv"), rows)
```
"""
function write_collision_catalogue(path::AbstractString, rows::AbstractVector{CollisionRow})::String
    _backup_existing(path)
    n_fields = fieldcount(CollisionRow)
    open(path, "w") do io
        println(io, join(fieldnames(CollisionRow), ","))
        for r in rows
            println(io, join(ntuple(f -> _collcat_cell(getfield(r, f)), n_fields), ","))
        end
    end
    return String(path)
end
