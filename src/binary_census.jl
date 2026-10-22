# =============================================================================
# Binary census: union of the engine-regularised (level E) and the physically
# bound (level P) spans of each pair
# =============================================================================
#
# A pair of component names can be a KS pair without being bound at a snapshot
# (a hyperbolic regularisation, a pair resolved only by KS) and bound without
# being regularised (a wide or soft binary). The census joins, per pair, the
# overlapping or adjacent spans of both levels into one record that keeps the
# limits and counts of each level (definition D5 of the binary-census note).

"""
    CensusSpan

One binary of the census from [`binary_census`](@ref): the union of the
overlapping or adjacent level-E ([`BinarySpan`](@ref)) and level-P
([`BoundSpan`](@ref)) spans of a pair of component names.

| field | meaning | unit |
|---|---|---|
| `name1`, `name2` | component names, `name1 < name2` | — |
| `primordial` | the primordial name rule of [`binary_spans`](@ref) | — |
| `levels` | `:E`, `:P` or `:EP`: the levels that contribute | — |
| `n_level_e`, `n_level_p` | contributing spans of each level | — |
| `t_form_nb`, `t_end_nb` | earliest formation (level-E opening, level-P first member snapshot) and latest end (level-E close, level-P last member snapshot) | NB |
| `t_form_e`, `t_end_e` | limits of the level-E spans; `NaN` without one | NB |
| `t_form_p`, `t_end_p` | first and last member snapshot of the level-P spans; `NaN` without one | NB |
| `left_censored`, `right_censored` | the earliest contribution is left-censored; a contribution is open at the end | — |
| `origin` | origin of the earliest level-E span when it opens no later than the level-P span, else `:primordial` (primordial pair bound at the first snapshot) or `:unresolved` | — |
| `end_class` | end class of the last level-E span when it ends no earlier than the level-P span; else `:censored_run_end`, `:separated` (unbound and apart, level P) or `:unresolved` | — |
| `ks_time_fraction` | KS time over the level-E extent (level-E fraction weighted by duration); `NaN` without level E | — |
"""
Base.@kwdef struct CensusSpan
    name1::Int
    name2::Int
    primordial::Bool
    levels::Symbol
    n_level_e::Int
    n_level_p::Int
    t_form_nb::Float64
    t_end_nb::Float64
    t_form_e::Float64
    t_end_e::Float64
    t_form_p::Float64
    t_end_p::Float64
    left_censored::Bool
    right_censored::Bool
    origin::Symbol
    end_class::Symbol
    ks_time_fraction::Float64
end

# Interval of a span for the union: level E from its opening to its close,
# level P widened to the midpoints between its boundary snapshots and their
# neighbours, so that adjacent spans of the two levels join while two level-P
# spans separated by one unexcused miss stay apart.
_census_interval(s::BinarySpan) = (s.t_form_nb, s.t_end_nb)
_census_interval(s::BoundSpan) = (
    isnan(s.t_form_lo) ? s.t_form_hi : (s.t_form_lo + s.t_form_hi) / 2,
    isnan(s.t_end_hi) ? s.t_end_lo : (s.t_end_lo + s.t_end_hi) / 2,
)

"""
    binary_census(level_e::AbstractVector{BinarySpan}, level_p::AbstractVector{BoundSpan};
                  nbin0::Integer) -> Vector{CensusSpan}

Union of the engine-regularised spans ([`binary_spans`](@ref)) and the
physically bound spans ([`bound_spans`](@ref)) of each pair: spans of the same
pair whose intervals overlap — a level-P span widened to the midpoints between
its first member snapshot and the one before, and between its last member
snapshot and the one after — are joined, transitively, into one
[`CensusSpan`](@ref). A pair can have several census spans (it separates and
re-forms). Flyby spans of level E (`:flyby`, never bound) are left out.
Returns the spans sorted by `(t_form_nb, name1, name2)`; `ArgumentError` for
`nbin0 < 0`.

# Example
```julia
e = binary_spans(read_binary_events(log); nzero, nbin0, tstar_myr)
p = bound_spans(read_all_conf3(dir); nbin0)
census = binary_census(e, p; nbin0)
```
"""
function binary_census(
    level_e::AbstractVector{BinarySpan},
    level_p::AbstractVector{BoundSpan};
    nbin0::Integer,
)::Vector{CensusSpan}
    nbin0 ≥ 0 || throw(ArgumentError("binary_census: nbin0 must be ≥ 0; got $nbin0"))
    groups = Dict{Tuple{Int,Int},Tuple{Vector{BinarySpan},Vector{BoundSpan}}}()
    for s in level_e
        s.end_class === :flyby && continue
        push!(get!(() -> (BinarySpan[], BoundSpan[]), groups, (s.name1, s.name2))[1], s)
    end
    for s in level_p
        push!(get!(() -> (BinarySpan[], BoundSpan[]), groups, (s.name1, s.name2))[2], s)
    end
    out = CensusSpan[]
    for (key, (es, ps)) in groups
        # Interval items: (start, end, level, index)
        items = Tuple{Float64,Float64,Int,Int}[]
        for (i, s) in enumerate(es)
            push!(items, (_census_interval(s)..., 1, i))
        end
        for (i, s) in enumerate(ps)
            push!(items, (_census_interval(s)..., 2, i))
        end
        sort!(items)
        group = Tuple{Float64,Float64,Int,Int}[]
        reach = -Inf
        for it in items
            if !isempty(group) && it[1] > reach
                push!(out, _census_span(key, group, es, ps, Int(nbin0)))
                empty!(group)
                reach = -Inf
            end
            push!(group, it)
            reach = max(reach, it[2])
        end
        isempty(group) || push!(out, _census_span(key, group, es, ps, Int(nbin0)))
    end
    sort!(out; by = s -> (s.t_form_nb, s.name1, s.name2))
    return out
end

function _census_span(
    key::Tuple{Int,Int},
    group::Vector{Tuple{Float64,Float64,Int,Int}},
    es::Vector{BinarySpan},
    ps::Vector{BoundSpan},
    nbin0::Int,
)::CensusSpan
    e = BinarySpan[es[it[4]] for it in group if it[3] == 1]
    p = BoundSpan[ps[it[4]] for it in group if it[3] == 2]
    t_form_e = isempty(e) ? NaN : minimum(s.t_form_nb for s in e)
    t_end_e = isempty(e) ? NaN : maximum(s.t_end_nb for s in e)
    t_form_p = isempty(p) ? NaN : minimum(s.t_form_hi for s in p)
    t_end_p = isempty(p) ? NaN : maximum(s.t_end_lo for s in p)
    first_e = isempty(e) ? nothing : e[argmin([s.t_form_nb for s in e])]
    last_e = isempty(e) ? nothing : e[argmax([s.t_end_nb for s in e])]
    first_p = isempty(p) ? nothing : p[argmin([s.t_form_hi for s in p])]
    last_p = isempty(p) ? nothing : p[argmax([s.t_end_lo for s in p])]
    e_first = first_e !== nothing && (first_p === nothing || t_form_e ≤ t_form_p)
    e_last = last_e !== nothing && (last_p === nothing || t_end_e ≥ t_end_p)
    primordial = _span_primordial(key, nbin0)
    origin = if e_first
        first_e.origin
    elseif primordial && first_p.left_censored
        :primordial
    else
        :unresolved
    end
    end_class = if e_last
        last_e.end_class
    elseif last_p.right_censored
        :censored_run_end
    elseif last_p.end_unbound
        :separated
    else
        :unresolved
    end
    left = e_first ? first_e.left_censored : first_p.left_censored
    right = any(s -> s.right_censored, e) || any(s -> s.right_censored, p)
    ks = if isempty(e)
        NaN
    else
        dur = sum(s.t_end_nb - s.t_form_nb for s in e)
        dur > 0 ? sum(s.ks_time_fraction * (s.t_end_nb - s.t_form_nb) for s in e) / dur : 1.0
    end
    return CensusSpan(;
        name1 = key[1],
        name2 = key[2],
        primordial = primordial,
        levels = isempty(e) ? :P : isempty(p) ? :E : :EP,
        n_level_e = length(e),
        n_level_p = length(p),
        t_form_nb = min(isnan(t_form_e) ? Inf : t_form_e, isnan(t_form_p) ? Inf : t_form_p),
        t_end_nb = max(isnan(t_end_e) ? -Inf : t_end_e, isnan(t_end_p) ? -Inf : t_end_p),
        t_form_e = t_form_e,
        t_end_e = t_end_e,
        t_form_p = t_form_p,
        t_end_p = t_end_p,
        left_censored = left,
        right_censored = right,
        origin = origin,
        end_class = end_class,
        ks_time_fraction = ks,
    )
end

"""
    write_binary_census(path, spans::AbstractVector{CensusSpan}) -> String

Write the census as CSV, one column per field of [`CensusSpan`](@ref) in
declaration order (non-finite values as empty cells). An existing file is
backed up, never overwritten. Returns the path.
"""
function write_binary_census(path::AbstractString, spans::AbstractVector{CensusSpan})::String
    _backup_existing(path)
    n_fields = fieldcount(CensusSpan)
    open(path, "w") do io
        println(io, join(fieldnames(CensusSpan), ","))
        for s in spans
            println(io, join(ntuple(f -> _csv_precise(getfield(s, f)), n_fields), ","))
        end
    end
    return String(path)
end
