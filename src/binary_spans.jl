# =============================================================================
# Binary spans: engine-regularised pairs joined across their KS episodes
# =============================================================================
#
# A pair that leaves KS regularisation and is regularised again within a few
# orbital periods (a soft binary terminated near apocentre, a pair carried
# through a chain or a hierarchical merger) is one binary, not several. The
# records of the engine log, of coll.13, coal.24 and of the bdat.9 listings
# are reduced to KS episodes per pair of names, the episodes are joined to
# spans, and every span receives an origin and an end class.

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------

# Two records are simultaneous when abs(Δt) ≤ rtol·max(1, abs(t)) (E17.10 stamps).
const _SPAN_DUP_RTOL = 1.0e-9
# abs(Δt) [NB] within which records of different files or families are joined (F9.2 stamps).
const _SPAN_JOIN_TOL = 1.0e-2
# Log lines before a closing record searched for a GW merger trigger that carries no time.
const _SPAN_TRIGGER_WINDOW = 200
# Days per Myr (Julian year of 365.25 d).
const _DAYS_PER_MYR = 3.6525e8

# KS kinds that open an episode
const _SPAN_KS_OPENING = (:ks_new, :ks_new_primordial, :ks_new_compact)
# Opening records whose elements give way to a simultaneous NEW KSREG line
const _SPAN_WEAK_OPENING = (:ks_new_compact, :ks_new_mass_loss)
# Few-body kinds that open and close a chain episode
const _SPAN_CHAIN_OPENING = (:chain_new, :triple_new, :quad_new, :chain_init)
const _SPAN_CHAIN_CLOSING = (:chain_end, :chterm)
# Few-body kinds that mark a three-body formation
const _SPAN_THREE_BODY_KINDS = (:chain_end, :chterm, :chterm_names, :chain_binary)

# ---------------------------------------------------------------------------
# Record type
# ---------------------------------------------------------------------------

"""
    BinarySpan

One engine-regularised binary (level E) from [`binary_spans`](@ref): a pair of
component names from the regularisation that opens it to the termination that
ends it, with the KS episodes in between joined.

| field | meaning | unit |
|---|---|---|
| `name1`, `name2` | component names, `name1 < name2` | — |
| `primordial` | the pair is (2k − 1, 2k) with 2k ≤ 2·NBIN0 | — |
| `hierarchy_outer` | a name is ≤ 0 or > NZERO: outer orbit of a merged hierarchy | — |
| `origin` | `:primordial`, `:primordial_reformed`, `:exchange` or `:dynamical` | — |
| `origin_sublabel` | `:gw_capture`, `:capture` or `:three_body` for a dynamical origin, `:none` otherwise | — |
| `mixed` | origin not `:primordial`, and a component was earlier in a span of another primordial pair | — |
| `t_form_nb`, `t_end_nb` | start and end of the span | NB |
| `t_form_myr`, `t_end_myr` | the same | Myr |
| `left_censored` | the first episode was already open at the first record of the pair | — |
| `right_censored` | the span is open at the end of the log | — |
| `end_class` | how the span ends, see [`binary_spans`](@ref) | — |
| `n_ks_episodes` | KS episodes joined into the span | — |
| `n_gaps_merged` | gaps between them, `n_ks_episodes − 1` | — |
| `max_gap_nb` | largest joined gap, 0 when none | NB |
| `ks_time_fraction` | Σ episode durations over the span duration; 1 for a span of zero duration | — |
| `in_chain_count` | joined gaps covered by a chain holding both components | — |
| `hierarchy_flag` | `hierarchy_outer`, or a joined gap covered by a hierarchy with the pair as inner binary | — |
| `outer_partner_name` | outer body of the last such hierarchy, 0 otherwise | — |
| `ecc_form`, `semi_form_nb`, `period_form_d`, `eb_form_nb` | e, a, period, binding energy of the opening record | —, NB, d, NB |
| `m1_form_msun`, `m2_form_msun` | component masses of the opening record | M☉ |
| `kstar1_form`, `kstar2_form`, `kstar_cm_form` | component and c.m. types of the opening record | — |
| `r_cm_form_nb`, `gamma_form` | c.m. distance and relative perturbation of the opening record | NB, — |
| `ecc_end`, `semi_end_nb`, `period_end_d`, `eb_end_nb` | the same elements of the closing record of the last episode | —, NB, d, NB |
| `m1_end_msun`, `m2_end_msun` | masses of that record | M☉ |
| `kstar1_end`, `kstar2_end`, `kstar_cm_end` | types of that record | — |
| `r_cm_end_nb`, `gamma_end` | c.m. distance and perturbation of that record | NB, — |
| `line_form`, `line_end` | log lines of the opening and the closing record | — |

The opening record is the NEW KSREG, NS/BH BINARY or MDNEW KSREG line of the
first episode, the closing record the END KSREG or MDEND KSREG line of the
last; index 1 of masses and types belongs to `name1` whatever the order of the
names in the record. MDEND/MDNEW lines carry no c.m. distance. Absent values are
`NaN` for reals, [`KSTAR_ABSENT`](@ref) for types and 0 for lines: the opening
elements of a left-censored span, the closing elements of a span whose last
episode ended without a closing record (a collision, a coalescence, an escape,
a new regularisation of the pair) or is still open.
"""
Base.@kwdef struct BinarySpan
    name1::Int
    name2::Int
    primordial::Bool
    hierarchy_outer::Bool
    origin::Symbol
    origin_sublabel::Symbol
    mixed::Bool
    t_form_nb::Float64
    t_end_nb::Float64
    t_form_myr::Float64
    t_end_myr::Float64
    left_censored::Bool
    right_censored::Bool
    end_class::Symbol
    n_ks_episodes::Int
    n_gaps_merged::Int
    max_gap_nb::Float64
    ks_time_fraction::Float64
    in_chain_count::Int
    hierarchy_flag::Bool
    outer_partner_name::Int
    ecc_form::Float64
    semi_form_nb::Float64
    period_form_d::Float64
    eb_form_nb::Float64
    m1_form_msun::Float64
    m2_form_msun::Float64
    kstar1_form::Int
    kstar2_form::Int
    kstar_cm_form::Int
    r_cm_form_nb::Float64
    gamma_form::Float64
    ecc_end::Float64
    semi_end_nb::Float64
    period_end_d::Float64
    eb_end_nb::Float64
    m1_end_msun::Float64
    m2_end_msun::Float64
    kstar1_end::Int
    kstar2_end::Int
    kstar_cm_end::Int
    r_cm_end_nb::Float64
    gamma_end::Float64
    line_form::Int
    line_end::Int
end

# ---------------------------------------------------------------------------
# Internal records
# ---------------------------------------------------------------------------

# Orbital elements of one KS record, masses and types ordered by name
struct _SpanElements
    ecc::Float64
    semi_nb::Float64
    period_d::Float64
    eb_nb::Float64
    m1_msun::Float64
    m2_msun::Float64
    kstar1::Int
    kstar2::Int
    kstar_cm::Int
    r_cm_nb::Float64
    gamma::Float64
end

const _SPAN_NO_ELEMENTS =
    _SpanElements(NaN, NaN, NaN, NaN, NaN, NaN, KSTAR_ABSENT, KSTAR_ABSENT, KSTAR_ABSENT, NaN, NaN)

# One marker of a pair: priority 1 closing, 2 opening, 3 terminal
struct _SpanMarker
    time::Float64
    line::Int
    priority::Int
    kind::Symbol
    elements::_SpanElements
end

# One KS episode of a pair; t_close = +Inf while it is open at the end of the log
mutable struct _SpanEpisode
    t_open::Float64
    t_close::Float64
    t_latest_open::Float64
    open_at_end::Bool
    left_censored::Bool
    kind_form::Symbol
    line_form::Int
    form::_SpanElements
    kind_end::Symbol
    line_end::Int
    final::_SpanElements
end

# One chain episode: its members (names 1…NZERO) and its time interval
mutable struct _SpanChain
    t_start::Float64
    t_end::Float64
    members::Set{Int}
end

# One hierarchy episode of an inner pair, with the name of the outer body
struct _SpanHierarchy
    t_start::Float64
    t_end::Float64
    outer::Int
end

# Joined episodes first_episode…last_episode of a pair, before origin and end class
struct _SpanRaw
    key::Tuple{Int,Int}
    first_episode::Int
    last_episode::Int
    in_chain_count::Int
    hierarchy_flag::Bool
    outer_partner_name::Int
    max_gap_nb::Float64
    t_form_nb::Float64
    t_end_nb::Float64
end

# Time stamps and record indices of the origin and end-class tests
struct _SpanLookup
    collisions::Dict{Tuple{Int,Int},Vector{Float64}}
    mergers::Dict{Tuple{Int,Int},Vector{Float64}}
    escapes::Dict{Tuple{Int,Int},Vector{Float64}}
    captures::Dict{Tuple{Int,Int},Vector{Float64}}
    gw_captures::Dict{Tuple{Int,Int},Vector{Float64}}
    listed::Dict{Tuple{Int,Int},Vector{Float64}}
    exchanges_of_name::Dict{Int,Vector{Int}}
    chain_binaries_of_name::Dict{Int,Vector{Int}}
    three_body_of_name::Dict{Int,Vector{Int}}
    gw_triggers::Vector{Tuple{Float64,Int}}
end

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

_span_valid(name1::Int, name2::Int)::Bool = name1 != 0 && name2 != 0 && name1 != name2
_span_key(name1::Int, name2::Int)::Tuple{Int,Int} = (min(name1, name2), max(name1, name2))

_span_primordial(key::Tuple{Int,Int}, nbin0::Int)::Bool =
    key[1] >= 1 && isodd(key[1]) && key[2] == key[1] + 1 && key[2] <= 2 * nbin0

_span_simultaneous(a::Float64, b::Float64)::Bool = abs(a - b) <= _SPAN_DUP_RTOL * max(1.0, abs(a))

_span_near(times::Vector{Float64}, t::Float64)::Bool = any(x -> abs(x - t) <= _SPAN_JOIN_TOL, times)

# A time stamp within [t_open, t_close] of an episode, widened by the join tolerance
_span_within(times::Vector{Float64}, t_open::Float64, t_close::Float64)::Bool =
    any(x -> t_open - _SPAN_JOIN_TOL <= x <= t_close + _SPAN_JOIN_TOL, times)

# Unbound elements: binding energy > 0 (μH of a hyperbolic pair) or a < 0
_span_unbound(e::_SpanElements)::Bool = e.eb_nb > 0 || e.semi_nb < 0

_span_covers(s::Float64, f::Float64, t_c::Float64, t_o::Float64)::Bool =
    s <= t_c + _SPAN_JOIN_TOL && f >= t_o - _SPAN_JOIN_TOL

# Exactly one of the two names of `key` is among {a, b}
_span_shares_one(key::Tuple{Int,Int}, a::Int, b::Int)::Bool =
    ((key[1] == a || key[1] == b) + (key[2] == a || key[2] == b)) == 1

_span_is_member(r::FewBodyEvent, n::Int, nzero::Int)::Bool = 1 <= n <= nzero && n in r.names

_span_members(r::FewBodyEvent, nzero::Int)::Set{Int} =
    Set{Int}(n for n in r.names if 1 <= n <= nzero)

_span_intersects(r::FewBodyEvent, members::Set{Int}, nzero::Int)::Bool =
    any(n -> 1 <= n <= nzero && n in members, r.names)

function _span_push_time!(
    d::Dict{Tuple{Int,Int},Vector{Float64}},
    name1::Int,
    name2::Int,
    t::Float64,
)::Nothing
    (_span_valid(name1, name2) && isfinite(t)) || return nothing
    push!(get!(Vector{Float64}, d, _span_key(name1, name2)), t)
    return nothing
end

_span_push_index!(d::Dict{Int,Vector{Int}}, n::Int, i::Int) = push!(get!(Vector{Int}, d, n), i)

function _span_elements(r::KSEvent)::_SpanElements
    swap = r.name1 > r.name2
    return _SpanElements(
        r.ecc,
        r.semi_nb,
        r.period_d,
        r.eb_nb,
        swap ? r.m2_msun : r.m1_msun,
        swap ? r.m1_msun : r.m2_msun,
        swap ? r.kstar2 : r.kstar1,
        swap ? r.kstar1 : r.kstar2,
        r.kstar_cm,
        r.r_cm_nb,
        r.gamma,
    )
end

# MDEND/MDNEW lines carry no c.m. distance relative to the density centre
function _span_elements(r::MassLossEvent)::_SpanElements
    swap = r.name1 > r.name2
    return _SpanElements(
        r.ecc,
        r.semi_nb,
        r.period_d,
        r.eb_nb,
        swap ? r.m2_msun : r.m1_msun,
        swap ? r.m1_msun : r.m2_msun,
        swap ? r.kstar2 : r.kstar1,
        swap ? r.kstar1 : r.kstar2,
        r.kstar_cm,
        NaN,
        r.gamma,
    )
end

# Latest finite time stamp of `records`, at least `t`
_span_latest(t::Float64, records::AbstractVector)::Float64 = foldl(
    (a, r) -> isfinite(r.time_nb) && r.time_nb > a ? Float64(r.time_nb) : a,
    records;
    init = t,
)

function _span_latest_time(
    log::BinaryEventLog,
    collisions::AbstractVector{CollisionRecord},
    coalescences::AbstractVector{CoalescenceRecord},
)::Float64
    t = -Inf
    t = _span_latest(t, log.ks)
    t = _span_latest(t, log.fewbody)
    t = _span_latest(t, log.hierarchy)
    t = _span_latest(t, log.exchange)
    t = _span_latest(t, log.termination)
    t = _span_latest(t, log.gw)
    t = _span_latest(t, log.escape)
    t = _span_latest(t, log.mass_loss)
    t = _span_latest(t, log.tidal)
    t = _span_latest(t, log.roche)
    t = _span_latest(t, collisions)
    t = _span_latest(t, coalescences)
    return t
end

# ---------------------------------------------------------------------------
# Step 1: KS episodes per pair
# ---------------------------------------------------------------------------

function _span_add_marker!(
    markers::Dict{Tuple{Int,Int},Vector{_SpanMarker}},
    name1::Int,
    name2::Int,
    m::_SpanMarker,
)::Nothing
    (_span_valid(name1, name2) && isfinite(m.time)) || return nothing
    push!(get!(Vector{_SpanMarker}, markers, _span_key(name1, name2)), m)
    return nothing
end

function _span_markers(
    log::BinaryEventLog,
    collisions::AbstractVector{CollisionRecord},
    coalescences::AbstractVector{CoalescenceRecord},
)::Dict{Tuple{Int,Int},Vector{_SpanMarker}}
    markers = Dict{Tuple{Int,Int},Vector{_SpanMarker}}()
    for r in log.ks
        priority = r.kind === :ks_end ? 1 : r.kind in _SPAN_KS_OPENING ? 2 : 0
        priority == 0 && continue
        m = _SpanMarker(r.time_nb, r.line, priority, r.kind, _span_elements(r))
        _span_add_marker!(markers, r.name1, r.name2, m)
    end
    for r in log.mass_loss
        priority = r.kind === :ks_end_mass_loss ? 1 : r.kind === :ks_new_mass_loss ? 2 : 0
        priority == 0 && continue
        m = _SpanMarker(r.time_nb, r.line, priority, r.kind, _span_elements(r))
        _span_add_marker!(markers, r.name1, r.name2, m)
    end
    for r in collisions
        m = _SpanMarker(r.time_nb, 0, 3, :collision, _SPAN_NO_ELEMENTS)
        _span_add_marker!(markers, r.name1, r.name2, m)
    end
    for r in coalescences
        m = _SpanMarker(r.time_nb, 0, 3, :coalescence, _SPAN_NO_ELEMENTS)
        _span_add_marker!(markers, r.name1, r.name2, m)
    end
    for r in log.roche
        r.kind === :roche_coal || continue
        m = _SpanMarker(r.time_nb, r.line, 3, r.kind, _SPAN_NO_ELEMENTS)
        _span_add_marker!(markers, r.name1, r.name2, m)
    end
    for r in log.escape
        r.kind === :binary_escape || continue
        m = _SpanMarker(r.time_nb, r.line, 3, r.kind, _SPAN_NO_ELEMENTS)
        _span_add_marker!(markers, r.name1, r.name2, m)
    end
    return markers
end

# Sorted times of the bdat.9 listings that hold each pair
function _span_listed(listings::AbstractVector{KSPairListing})::Dict{Tuple{Int,Int},Vector{Float64}}
    listed = Dict{Tuple{Int,Int},Vector{Float64}}()
    for l in listings
        for r in l.records
            _span_push_time!(listed, r.name1, r.name2, l.time_nb)
        end
    end
    foreach(sort!, values(listed))
    return listed
end

# Opening time of a left-censored episode closed at t_close
function _span_censored_open(
    listed::Vector{Float64},
    t_prev::Float64,
    t_close::Float64,
    t_start::Float64,
)::Float64
    i = searchsortedlast(listed, t_prev) + 1
    i <= length(listed) && listed[i] <= t_close && return listed[i]
    return min(max(t_start, t_prev), t_close)
end

# Episodes of one pair from its markers, sorted by (time, priority, line)
function _span_walk(
    markers::Vector{_SpanMarker},
    listed::Vector{Float64},
    t_start::Float64,
)::Vector{_SpanEpisode}
    episodes = _SpanEpisode[]
    current = 0     # index of the open episode, 0 while closed
    t_prev = -Inf   # close time of the last episode
    for m in markers
        if m.priority == 2
            if current != 0
                ep = episodes[current]
                if _span_simultaneous(m.time, ep.t_latest_open)
                    # The same regularisation written twice: keep the NEW KSREG elements
                    if ep.kind_form in _SPAN_WEAK_OPENING && m.kind in (:ks_new, :ks_new_primordial)
                        ep.kind_form = m.kind
                        ep.line_form = m.line
                        ep.form = m.elements
                    end
                    ep.t_latest_open = m.time
                    continue
                end
                # A new regularisation without termination: close without a record
                ep.t_close = m.time
                t_prev = m.time
            end
            push!(
                episodes,
                _SpanEpisode(
                    m.time,
                    Inf,
                    m.time,
                    false,
                    false,
                    m.kind,
                    m.line,
                    m.elements,
                    :none,
                    0,
                    _SPAN_NO_ELEMENTS,
                ),
            )
            current = length(episodes)
        elseif m.priority == 1
            if current != 0
                ep = episodes[current]
                ep.t_close = m.time
                ep.kind_end = m.kind
                ep.line_end = m.line
                ep.final = m.elements
                current = 0
                t_prev = m.time
            elseif m.time - t_prev > _SPAN_JOIN_TOL
                # A termination without its opening: the episode was open before the log
                t_open = _span_censored_open(listed, t_prev, m.time, t_start)
                push!(
                    episodes,
                    _SpanEpisode(
                        t_open,
                        m.time,
                        t_open,
                        false,
                        true,
                        :none,
                        0,
                        _SPAN_NO_ELEMENTS,
                        m.kind,
                        m.line,
                        m.elements,
                    ),
                )
                t_prev = m.time
            end
        elseif current != 0
            episodes[current].t_close = m.time
            current = 0
            t_prev = m.time
        end
    end
    current != 0 && (episodes[current].open_at_end = true)
    return episodes
end

# Period of an episode [NB]: closing record, else opening record, else NaN
function _span_period_nb(ep::_SpanEpisode, tstar::Float64)::Float64
    scale = tstar * _DAYS_PER_MYR
    p = ep.final.period_d
    isfinite(p) && p > 0 && return p / scale
    p = ep.form.period_d
    isfinite(p) && p > 0 && return p / scale
    return NaN
end

# ---------------------------------------------------------------------------
# Step 2: chain and hierarchy episodes
# ---------------------------------------------------------------------------

function _span_chains(
    log::BinaryEventLog,
    nzero::Int,
)::Tuple{Vector{_SpanChain},Dict{Int,Vector{Int}}}
    records = sort!(filter(r -> isfinite(r.time_nb), log.fewbody); by = r -> r.line)
    chains = _SpanChain[]
    current = 0
    for r in records
        if current == 0
            r.kind in _SPAN_CHAIN_OPENING || continue
            push!(chains, _SpanChain(r.time_nb, Inf, _span_members(r, nzero)))
            current = length(chains)
            continue
        end
        chain = chains[current]
        _span_intersects(r, chain.members, nzero) || continue
        union!(chain.members, _span_members(r, nzero))
        if r.kind in _SPAN_CHAIN_CLOSING
            chain.t_end = r.time_nb
            current = 0
        end
    end
    of_name = Dict{Int,Vector{Int}}()
    for (i, c) in enumerate(chains), n in c.members
        _span_push_index!(of_name, n, i)
    end
    return chains, of_name
end

function _span_hierarchies(log::BinaryEventLog)::Dict{Tuple{Int,Int},Vector{_SpanHierarchy}}
    records = sort!(filter(r -> isfinite(r.time_nb), log.hierarchy); by = r -> r.line)
    episodes = Dict{Tuple{Int,Int},Vector{_SpanHierarchy}}()
    pending = Dict{Tuple{Int,Int},Vector{Int}}()
    for r in records
        n1, n2 = r.names[1], r.names[2]
        _span_valid(n1, n2) || continue
        key = _span_key(n1, n2)
        if r.kind === :hierarchy_new || r.kind === :hierarchy_new_merge2
            v = get!(Vector{_SpanHierarchy}, episodes, key)
            push!(v, _SpanHierarchy(r.time_nb, Inf, r.names[3]))
            push!(get!(Vector{Int}, pending, key), length(v))
        elseif r.kind === :hierarchy_end && haskey(pending, key)
            v = episodes[key]
            for i in pending[key]
                v[i] = _SpanHierarchy(v[i].t_start, r.time_nb, v[i].outer)
            end
            empty!(pending[key])
        end
    end
    return episodes
end

function _span_chain_covers(
    chains::Vector{_SpanChain},
    chains_of_name::Dict{Int,Vector{Int}},
    key::Tuple{Int,Int},
    t_c::Float64,
    t_o::Float64,
)::Bool
    for i in get(chains_of_name, key[1], Int[])
        c = chains[i]
        key[2] in c.members && _span_covers(c.t_start, c.t_end, t_c, t_o) && return true
    end
    return false
end

# Index of the last hierarchy episode covering the gap, 0 when none
function _span_hierarchy_cover(hierarchy::Vector{_SpanHierarchy}, t_c::Float64, t_o::Float64)::Int
    found = 0
    for (i, h) in enumerate(hierarchy)
        _span_covers(h.t_start, h.t_end, t_c, t_o) && (found = i)
    end
    return found
end

# ---------------------------------------------------------------------------
# Step 3: spans
# ---------------------------------------------------------------------------

function _span_episodes_by_name(
    episodes::Dict{Tuple{Int,Int},Vector{_SpanEpisode}},
)::Dict{Int,Vector{Tuple{Tuple{Int,Int},Int}}}
    of_name = Dict{Int,Vector{Tuple{Tuple{Int,Int},Int}}}()
    for (key, eps) in episodes, j in eachindex(eps), n in key
        push!(get!(Vector{Tuple{Tuple{Int,Int},Int}}, of_name, n), (key, j))
    end
    return of_name
end

# Whether an episode of another pair holding a component overlaps the gap (t_c, t_o)
function _span_blocked(
    key::Tuple{Int,Int},
    t_c::Float64,
    t_o::Float64,
    episodes::Dict{Tuple{Int,Int},Vector{_SpanEpisode}},
    episodes_of_name::Dict{Int,Vector{Tuple{Tuple{Int,Int},Int}}},
)::Bool
    lo = t_o - _SPAN_DUP_RTOL * max(1.0, abs(t_o))
    hi = t_c + _SPAN_DUP_RTOL * max(1.0, abs(t_c))
    for n in key, (other, j) in get(episodes_of_name, n, Tuple{Tuple{Int,Int},Int}[])
        other == key && continue
        ep = episodes[other][j]
        ep.t_open < lo && ep.t_close > hi && return true
    end
    return false
end

# Spans of one pair, appended to `raw` in time order
function _span_join!(
    raw::Vector{_SpanRaw},
    key::Tuple{Int,Int},
    episodes::Dict{Tuple{Int,Int},Vector{_SpanEpisode}},
    episodes_of_name::Dict{Int,Vector{Tuple{Tuple{Int,Int},Int}}},
    chains::Vector{_SpanChain},
    chains_of_name::Dict{Int,Vector{Int}},
    hierarchies::Dict{Tuple{Int,Int},Vector{_SpanHierarchy}},
    nzero::Int,
    tstar::Float64,
    n_p::Float64,
    t_open_end::Float64,
)::Vector{_SpanRaw}
    eps = episodes[key]
    outer_orbit = key[1] <= 0 || key[2] > nzero
    hierarchy = get(hierarchies, key, _SpanHierarchy[])
    start = 1
    n_chain = 0
    in_hierarchy = false
    outer = 0
    max_gap = 0.0
    for e in eachindex(eps)
        if e < length(eps)
            a = eps[e]
            b = eps[e + 1]
            gap = b.t_open - a.t_close
            period = _span_period_nb(a, tstar)
            g1 = gap <= 0 || (isfinite(period) && gap <= n_p * period)
            chain = _span_chain_covers(chains, chains_of_name, key, a.t_close, b.t_open)
            h = _span_hierarchy_cover(hierarchy, a.t_close, b.t_open)
            if (g1 || chain || h > 0) &&
               !_span_blocked(key, a.t_close, b.t_open, episodes, episodes_of_name)
                n_chain += chain
                if h > 0
                    in_hierarchy = true
                    outer = hierarchy[h].outer
                end
                max_gap = max(max_gap, gap)
                continue
            end
        end
        last_ep = eps[e]
        t_end = last_ep.open_at_end ? t_open_end : last_ep.t_close
        push!(
            raw,
            _SpanRaw(
                key,
                start,
                e,
                n_chain,
                outer_orbit || in_hierarchy,
                outer,
                max_gap,
                eps[start].t_open,
                t_end,
            ),
        )
        start = e + 1
        n_chain = 0
        in_hierarchy = false
        outer = 0
        max_gap = 0.0
    end
    return raw
end

# ---------------------------------------------------------------------------
# Steps 4 and 5: origin and end class
# ---------------------------------------------------------------------------

function _span_lookup(
    log::BinaryEventLog,
    collisions::AbstractVector{CollisionRecord},
    coalescences::AbstractVector{CoalescenceRecord},
    listed::Dict{Tuple{Int,Int},Vector{Float64}},
    nzero::Int,
)::_SpanLookup
    coll = Dict{Tuple{Int,Int},Vector{Float64}}()
    for r in collisions
        _span_push_time!(coll, r.name1, r.name2, r.time_nb)
    end
    mergers = Dict{Tuple{Int,Int},Vector{Float64}}()
    for r in coalescences
        _span_push_time!(mergers, r.name1, r.name2, r.time_nb)
    end
    for r in log.roche
        r.kind === :roche_coal && _span_push_time!(mergers, r.name1, r.name2, r.time_nb)
    end
    escapes = Dict{Tuple{Int,Int},Vector{Float64}}()
    for r in log.escape
        r.kind === :binary_escape && _span_push_time!(escapes, r.name1, r.name2, r.time_nb)
    end
    captures = Dict{Tuple{Int,Int},Vector{Float64}}()
    for r in log.tidal
        r.kind === :capture && _span_push_time!(captures, r.name1, r.name2, r.time_nb)
    end
    gw_captures = Dict{Tuple{Int,Int},Vector{Float64}}()
    for r in log.gw
        r.kind === :gw_capture && _span_push_time!(gw_captures, r.name1, r.name2, r.time_nb)
    end
    exchanges_of_name = Dict{Int,Vector{Int}}()
    for (i, r) in enumerate(log.exchange)
        isfinite(r.time_nb) || continue
        for n in (r.name1_initial, r.name2_initial, r.name1_final, r.name2_final)
            n == 0 || _span_push_index!(exchanges_of_name, n, i)
        end
    end
    chain_binaries_of_name = Dict{Int,Vector{Int}}()
    three_body_of_name = Dict{Int,Vector{Int}}()
    for (i, r) in enumerate(log.fewbody)
        isfinite(r.time_nb) || continue
        if r.kind === :chain_binary
            for n in (r.names[1], r.names[2])
                n == 0 || _span_push_index!(chain_binaries_of_name, n, i)
            end
        end
        if r.kind in _SPAN_THREE_BODY_KINDS
            for n in r.names
                1 <= n <= nzero && _span_push_index!(three_body_of_name, n, i)
            end
        end
    end
    gw_triggers = Tuple{Float64,Int}[
        (r.time_nb, r.line) for
        r in log.gw if r.kind === :gw_ksint_coll || r.kind === :gw_unpert_coll
    ]
    return _SpanLookup(
        coll,
        mergers,
        escapes,
        captures,
        gw_captures,
        listed,
        exchanges_of_name,
        chain_binaries_of_name,
        three_body_of_name,
        gw_triggers,
    )
end

# An exchange or chain binary that leaves the pair `key` bound at time t
function _span_exchange_formed(
    key::Tuple{Int,Int},
    t::Float64,
    lookup::_SpanLookup,
    log::BinaryEventLog,
)::Bool
    for j in get(lookup.exchanges_of_name, key[1], Int[])
        r = log.exchange[j]
        _span_key(r.name1_final, r.name2_final) == key &&
            abs(r.time_nb - t) <= _SPAN_JOIN_TOL &&
            return true
    end
    for j in get(lookup.chain_binaries_of_name, key[1], Int[])
        r = log.fewbody[j]
        _span_key(r.names[1], r.names[2]) == key &&
            abs(r.time_nb - t) <= _SPAN_JOIN_TOL &&
            return true
    end
    return false
end

function _span_origin(
    s::_SpanRaw,
    t_first_close::Float64,
    primordial::Bool,
    first_of_key::Bool,
    before::Bool,
    lookup::_SpanLookup,
    log::BinaryEventLog,
    t_start::Float64,
    nzero::Int,
)::Tuple{Symbol,Symbol}
    key = s.key
    t = s.t_form_nb
    if primordial
        listed = get(lookup.listed, key, Float64[])
        at_start =
            _span_simultaneous(t, t_start) ||
            (!isempty(listed) && listed[1] <= t_start + _SPAN_JOIN_TOL)
        return (first_of_key && at_start ? :primordial : :primordial_reformed, :none)
    end
    before && _span_exchange_formed(key, t, lookup, log) && return (:exchange, :none)
    # A capture at a pericentre of the first episode (opened by the hyperbolic NEW KSREG)
    _span_within(get(lookup.gw_captures, key, Float64[]), t, t_first_close) &&
        return (:dynamical, :gw_capture)
    _span_within(get(lookup.captures, key, Float64[]), t, t_first_close) &&
        return (:dynamical, :capture)
    for j in get(lookup.three_body_of_name, key[1], Int[])
        r = log.fewbody[j]
        abs(r.time_nb - t) <= _SPAN_JOIN_TOL &&
            _span_is_member(r, key[2], nzero) &&
            return (:dynamical, :three_body)
    end
    return (:dynamical, :none)
end

# A GW merger trigger at t_c, or without time within the lines before line_end
function _span_gw_trigger(triggers::Vector{Tuple{Float64,Int}}, t_c::Float64, line_end::Int)::Bool
    for (t, l) in triggers
        if isfinite(t)
            abs(t - t_c) <= _SPAN_JOIN_TOL && return true
        elseif line_end > 0 && line_end - _SPAN_TRIGGER_WINDOW <= l < line_end
            return true
        end
    end
    return false
end

# An exchange or chain binary in [t_c − tol, t_next] that takes the pair apart
function _span_exchange_ended(
    key::Tuple{Int,Int},
    t_c::Float64,
    t_next::Float64,
    lookup::_SpanLookup,
    log::BinaryEventLog,
)::Bool
    lo = t_c - _SPAN_JOIN_TOL
    for n in key, j in get(lookup.exchanges_of_name, n, Int[])
        r = log.exchange[j]
        lo <= r.time_nb <= t_next || continue
        _span_key(r.name1_initial, r.name2_initial) == key && return true
        _span_shares_one(key, r.name1_final, r.name2_final) && return true
    end
    for n in key, j in get(lookup.chain_binaries_of_name, n, Int[])
        r = log.fewbody[j]
        lo <= r.time_nb <= t_next && _span_shares_one(key, r.names[1], r.names[2]) && return true
    end
    return false
end

function _span_end_class(
    ep::_SpanEpisode,
    form::_SpanElements,
    key::Tuple{Int,Int},
    t_next::Float64,
    lookup::_SpanLookup,
    log::BinaryEventLog,
)::Symbol
    ep.open_at_end && return :censored_run_end
    t_c = ep.t_close
    collision = _span_near(get(lookup.collisions, key, Float64[]), t_c)
    collision && _span_gw_trigger(lookup.gw_triggers, t_c, ep.line_end) && return :gw_merger
    _span_near(get(lookup.mergers, key, Float64[]), t_c) && return :stellar_merger
    collision && return :collision
    _span_near(get(lookup.escapes, key, Float64[]), t_c) && return :escape_bound
    _span_exchange_ended(key, t_c, t_next, lookup, log) && return :exchange
    ep.kind_end === :ks_end_mass_loss && return :sn_disruption
    if ep.kind_end === :ks_end && _span_unbound(ep.final)
        # Unbound at both ends: a hyperbolic encounter, not a disrupted binary
        return _span_unbound(form) ? :flyby : :disruption
    end
    return :unresolved
end

# The span raw[i]; the spans of one pair are consecutive in `raw`, in time order
function _span_record(
    i::Int,
    raw::Vector{_SpanRaw},
    episodes::Dict{Tuple{Int,Int},Vector{_SpanEpisode}},
    spans_of_name::Dict{Int,Vector{Int}},
    lookup::_SpanLookup,
    log::BinaryEventLog,
    nbin0::Int,
    nzero::Int,
    tstar::Float64,
    t_start::Float64,
)::BinarySpan
    s = raw[i]
    key = s.key
    eps = episodes[key]
    first_ep = eps[s.first_episode]
    last_ep = eps[s.last_episode]
    primordial = _span_primordial(key, nbin0)
    first_of_key = i == 1 || raw[i - 1].key != key
    t_next = i < length(raw) && raw[i + 1].key == key ? raw[i + 1].t_form_nb : Inf

    # Components that were in a span of another pair formed earlier
    before = false
    before_primordial = false
    for n in key, j in get(spans_of_name, n, Int[])
        o = raw[j]
        (o.key != key && o.t_form_nb < s.t_form_nb) || continue
        before = true
        before_primordial |= _span_primordial(o.key, nbin0)
    end
    t_first_close = first_ep.open_at_end ? s.t_end_nb : first_ep.t_close
    origin, sublabel = _span_origin(
        s,
        t_first_close,
        primordial,
        first_of_key,
        before,
        lookup,
        log,
        t_start,
        nzero,
    )
    end_class = _span_end_class(last_ep, first_ep.form, key, t_next, lookup, log)

    t_end = s.t_end_nb
    ks_time = 0.0
    for e in (s.first_episode):(s.last_episode)
        ep = eps[e]
        ks_time += (ep.open_at_end ? t_end : ep.t_close) - ep.t_open
    end
    duration = t_end - s.t_form_nb
    n_episodes = s.last_episode - s.first_episode + 1
    f = first_ep.form
    z = last_ep.final
    return BinarySpan(;
        name1 = key[1],
        name2 = key[2],
        primordial = primordial,
        hierarchy_outer = key[1] <= 0 || key[2] > nzero,
        origin = origin,
        origin_sublabel = sublabel,
        mixed = origin !== :primordial && before_primordial,
        t_form_nb = s.t_form_nb,
        t_end_nb = t_end,
        t_form_myr = s.t_form_nb * tstar,
        t_end_myr = t_end * tstar,
        left_censored = first_ep.left_censored,
        right_censored = last_ep.open_at_end,
        end_class = end_class,
        n_ks_episodes = n_episodes,
        n_gaps_merged = n_episodes - 1,
        max_gap_nb = s.max_gap_nb,
        ks_time_fraction = duration == 0 ? 1.0 : ks_time / duration,
        in_chain_count = s.in_chain_count,
        hierarchy_flag = s.hierarchy_flag,
        outer_partner_name = s.outer_partner_name,
        ecc_form = f.ecc,
        semi_form_nb = f.semi_nb,
        period_form_d = f.period_d,
        eb_form_nb = f.eb_nb,
        m1_form_msun = f.m1_msun,
        m2_form_msun = f.m2_msun,
        kstar1_form = f.kstar1,
        kstar2_form = f.kstar2,
        kstar_cm_form = f.kstar_cm,
        r_cm_form_nb = f.r_cm_nb,
        gamma_form = f.gamma,
        ecc_end = z.ecc,
        semi_end_nb = z.semi_nb,
        period_end_d = z.period_d,
        eb_end_nb = z.eb_nb,
        m1_end_msun = z.m1_msun,
        m2_end_msun = z.m2_msun,
        kstar1_end = z.kstar1,
        kstar2_end = z.kstar2,
        kstar_cm_end = z.kstar_cm,
        r_cm_end_nb = z.r_cm_nb,
        gamma_end = z.gamma,
        line_form = first_ep.line_form,
        line_end = last_ep.line_end,
    )
end

# ---------------------------------------------------------------------------
# Public interface
# ---------------------------------------------------------------------------

"""
    binary_spans(log::BinaryEventLog;
                 listings::AbstractVector{KSPairListing} = KSPairListing[],
                 collisions::AbstractVector{CollisionRecord} = CollisionRecord[],
                 coalescences::AbstractVector{CoalescenceRecord} = CoalescenceRecord[],
                 nzero::Integer, nbin0::Integer, tstar_myr::Real,
                 t_start_nb::Real = 0.0, t_end_nb::Real = NaN, n_p::Real = 2.0)
        -> Vector{BinarySpan}

The engine-regularised binaries of a run (level E): the KS episodes of every
pair of component names in the engine log, joined to [`BinarySpan`](@ref)s,
each with its origin and its end class.

# Arguments
- `log`: binary and few-body records of the engine log, from [`read_binary_events`](@ref)
- `listings`: the `bdat.9_<t>` listings ([`read_all_ks_pair_listings`](@ref)); they date the
  opening of an episode whose NEW KSREG line is not in the log
- `collisions`, `coalescences`: `coll.13` and `coal.24` ([`read_collision_log`](@ref),
  [`read_coalescence_log`](@ref))
- `nzero`: initial number of bodies NZERO; a name ≤ 0 or above it is a c.m. body
- `nbin0`: initial number of primordial binaries NBIN0
- `tstar_myr`: time unit T* [Myr]
- `t_start_nb`: time at which the log starts [NB]
- `t_end_nb`: end given to spans still open at the end of the log [NB]; with `NaN` the
  latest time stamp of all records
- `n_p`: largest gap, in orbital periods of the earlier episode, that is joined

# Method
Names are compared as unordered pairs; records with a name 0 or two equal names
are ignored.

1. Episodes. A NEW KSREG, NS/BH BINARY or MDNEW KSREG line opens an episode of
   its pair, an END KSREG or MDEND KSREG line closes it, and a collision
   (`coll.13`), a coalescence (`coal.24`, ROCHE COAL) or a BINARY ESCAPE of the
   pair ends an open episode without a closing record. An opening simultaneous
   with the one that opened the episode is the same regularisation written twice
   (its NEW KSREG elements replace those of an NS/BH BINARY or MDNEW line); a
   later opening of an open episode ends it. A closing record of a closed pair
   is ignored within $(_SPAN_JOIN_TOL) NB of the last close, and otherwise closes a
   left-censored episode opened at the first listing that holds the pair after
   the last close, else at `t_start_nb` or the last close.
2. Coverage. Chain episodes run from NEW CHAIN/TRIPLE/QUAD to the first END
   CHAIN or CHTERM record whose members (names 1…NZERO) intersect theirs, and
   collect the members of every record in between that intersects them;
   hierarchy episodes run from NEW MERGER or NEW MERGE2 of an inner pair to the
   next END MERGER of the same pair. Unclosed episodes run to the end.
3. Spans. Consecutive episodes of a pair are joined when the gap g between them
   is ≤ 0 or at most `n_p` periods of the earlier episode (its closing record,
   else its opening record), or when a chain holding both components or a
   hierarchy with the pair as inner binary covers the gap, and in either case
   no episode of another pair sharing a component overlaps the gap.
4. Origin. `:primordial` for the first span of a primordial pair, (2k − 1, 2k)
   with 2k ≤ 2·NBIN0, formed at `t_start_nb` or listed at the start;
   `:primordial_reformed` for its later spans; `:exchange` when an EXCHANGE or
   CHAIN BINARY line leaves the pair bound at the formation time and a
   component was in a span of another pair before; `:dynamical` otherwise, with
   the sublabel `:gw_capture` (GR CAPTURE of the pair) or `:capture` (NEW
   CAPTURE, tidal) at a pericentre of the first episode, which the hyperbolic
   NEW KSREG opened, else `:three_body` (END CHAIN, CHTERM or CHAIN BINARY
   with both components as members) at the formation time. `mixed` marks a non-primordial origin with a component that
   was earlier in a span of another primordial pair.

The end class is that of the last episode, the first that applies of:

1. `:censored_run_end`: open at the end of the log;
2. `:gw_merger`: a collision of the pair at the close with a GR KSINT COLL or
   COLL DUE TO GW UNPERT trigger, at the same time or, the trigger carrying no
   time, within $(_SPAN_TRIGGER_WINDOW) log lines before the closing record;
3. `:stellar_merger`: a coalescence (`coal.24`, ROCHE COAL) of the pair;
4. `:collision`: a collision of the pair without a GW trigger;
5. `:escape_bound`: a BINARY ESCAPE of the pair;
6. `:exchange`: an EXCHANGE line between the close and the next span of the pair
   whose initial pair is this one or whose final pair holds one of its
   components, or a CHAIN BINARY line holding one of its components;
7. `:sn_disruption`: closed by MDEND KSREG;
8. `:flyby`: closed by END KSREG unbound (positive binding energy or a < 0)
   and opened unbound: a hyperbolic encounter that never was a binary;
9. `:disruption`: closed by END KSREG unbound after a bound opening;
10. `:unresolved`: none of these.

Tolerances: records of one family are simultaneous within a relative 10⁻⁹
(E17.10 stamps); records of different files or families
(F9.2 stamps) are joined within $(_SPAN_JOIN_TOL) NB.

Returns the spans sorted by `(t_form_nb, name1, name2)`, an empty vector for an
empty log. Throws `ArgumentError` for `nzero ≤ 0`, `nbin0 < 0`, a `tstar_myr` or
`n_p` that is not finite and positive, or a non-finite `t_start_nb`.

# Examples
```julia
dir = "runs/<run_id>/output"
log = read_binary_events(joinpath(dir, "out1000"))
spans = binary_spans(
    log;
    listings = read_all_ks_pair_listings(dir),
    collisions = read_collision_log(joinpath(dir, "coll.13")),
    coalescences = read_coalescence_log(joinpath(dir, "coal.24")),
    nzero = 20_000,
    nbin0 = 2_000,
    tstar_myr = 0.87,
)
write_binary_spans(joinpath(dir, "binary_spans.csv"), spans)
```
"""
function binary_spans(
    log::BinaryEventLog;
    listings::AbstractVector{KSPairListing} = KSPairListing[],
    collisions::AbstractVector{CollisionRecord} = CollisionRecord[],
    coalescences::AbstractVector{CoalescenceRecord} = CoalescenceRecord[],
    nzero::Integer,
    nbin0::Integer,
    tstar_myr::Real,
    t_start_nb::Real = 0.0,
    t_end_nb::Real = NaN,
    n_p::Real = 2.0,
)::Vector{BinarySpan}
    nzero > 0 || throw(ArgumentError("binary_spans: nzero must be > 0; got $nzero"))
    nbin0 >= 0 || throw(ArgumentError("binary_spans: nbin0 must be >= 0; got $nbin0"))
    (isfinite(tstar_myr) && tstar_myr > 0) ||
        throw(ArgumentError("binary_spans: tstar_myr must be finite and > 0; got $tstar_myr"))
    (isfinite(n_p) && n_p > 0) ||
        throw(ArgumentError("binary_spans: n_p must be finite and > 0; got $n_p"))
    isfinite(t_start_nb) ||
        throw(ArgumentError("binary_spans: t_start_nb must be finite; got $t_start_nb"))
    nz = Int(nzero)
    nb = Int(nbin0)
    tstar = Float64(tstar_myr)
    t_start = Float64(t_start_nb)
    periods = Float64(n_p)

    markers = _span_markers(log, collisions, coalescences)
    isempty(markers) && return BinarySpan[]
    listed = _span_listed(listings)

    # Step 1
    episodes = Dict{Tuple{Int,Int},Vector{_SpanEpisode}}()
    for (key, ms) in markers
        sort!(ms; by = m -> (m.time, m.priority, m.line))
        eps = _span_walk(ms, get(listed, key, Float64[]), t_start)
        isempty(eps) || (episodes[key] = eps)
    end

    # Steps 2 and 3
    t_open_end =
        isfinite(t_end_nb) ? Float64(t_end_nb) : _span_latest_time(log, collisions, coalescences)
    chains, chains_of_name = _span_chains(log, nz)
    hierarchies = _span_hierarchies(log)
    episodes_of_name = _span_episodes_by_name(episodes)
    raw = _SpanRaw[]
    for key in sort!(collect(keys(episodes)))
        _span_join!(
            raw,
            key,
            episodes,
            episodes_of_name,
            chains,
            chains_of_name,
            hierarchies,
            nz,
            tstar,
            periods,
            t_open_end,
        )
    end

    # Steps 4 and 5
    lookup = _span_lookup(log, collisions, coalescences, listed, nz)
    spans_of_name = Dict{Int,Vector{Int}}()
    for (i, s) in enumerate(raw), n in s.key
        _span_push_index!(spans_of_name, n, i)
    end
    spans = BinarySpan[
        _span_record(i, raw, episodes, spans_of_name, lookup, log, nb, nz, tstar, t_start) for
        i in eachindex(raw)
    ]
    sort!(spans; by = s -> (s.t_form_nb, s.name1, s.name2))
    return spans
end

"""
    write_binary_spans(path, spans::AbstractVector{BinarySpan}) -> String

Write the spans of [`binary_spans`](@ref) as CSV, one column per field of
[`BinarySpan`](@ref) in declaration order (times in N-body units and Myr,
symbols as text, non-finite values as empty cells, absent stellar types as
$(KSTAR_ABSENT)). An existing file is backed up, never overwritten. Returns the
path.
"""
function write_binary_spans(path::AbstractString, spans::AbstractVector{BinarySpan})::String
    _backup_existing(path)
    n_fields = fieldcount(BinarySpan)
    open(path, "w") do io
        println(io, join(fieldnames(BinarySpan), ","))
        for s in spans
            println(io, join(ntuple(f -> _csv_precise(getfield(s, f)), n_fields), ","))
        end
    end
    return String(path)
end
