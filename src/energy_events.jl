# =============================================================================
# Localisation of energy jumps in the adjustment intervals
# =============================================================================
# The engine's ADJUST line prints the energy change of the adjustment interval
# (DELTA, its relative form DE) and the running values of its correction
# channels (EBIN, EMERGE, ESUB, ECOLL, EMDOT, ECDOT). The routine evaluates the
# energy, then removes escapers, then writes output, so every event line
# between two ADJUST lines belongs to the interval that ends at the second.
# One record per interval holds the increments of the channels and the counts
# of the events the log places in it; the intervals whose |DE| stands out are
# flagged.

"""
    EnergyInterval

One adjustment interval `(t_start, t_end]` of the engine log, from
[`energy_intervals`](@ref): `index` is the ADJUST record that ends it (`k ≥ 2`
in `diag.adjust`); `de_rel` and `de_abs` its `DE` and `DELTA`; `d_ebin`,
`d_emerge`, `d_esub`, `d_ecoll`, `d_emdot`, `d_ecdot` the increments of the
engine's correction channels over the interval (N-body units; NaN when either
end lacks the term); `counts` the number of events of each kind of
`ENGINE_EVENT_KINDS`, in that order, and `n_events` their sum; `score` = |DE|
over the robust scale of the run's DE series (NaN when the scale is zero);
`flagged` whether |DE| reaches the flagging level.
"""
struct EnergyInterval
    index::Int
    t_start::Float64
    t_end::Float64
    de_rel::Float64
    de_abs::Float64
    d_ebin::Float64
    d_emerge::Float64
    d_esub::Float64
    d_ecoll::Float64
    d_emdot::Float64
    d_ecdot::Float64
    counts::NTuple{16,Int}
    n_events::Int
    score::Float64
    flagged::Bool
end

"""
    EnergyEventAnalysis

Result of [`energy_intervals`](@ref): the `intervals`, the `events` they were
built from, the flagging parameters `threshold` and `sigma`, the robust
`scale` of the DE series (1.4826 × MAD over the intervals with finite `DE`;
0 when fewer than two) and `engine_energy_valid`, false when the engine's DE
is not an error measure (a linearised tidal field, see [`tidal_mode`](@ref)),
in which case no interval is flagged.
"""
struct EnergyEventAnalysis
    intervals::Vector{EnergyInterval}
    events::Vector{EngineEvent}
    threshold::Float64
    sigma::Float64
    scale::Float64
    engine_energy_valid::Bool
end

"""
    energy_intervals(diag::DiagnosticsData, events::AbstractVector{EngineEvent};
                     threshold = 1.0e-4, sigma = 5.0, engine_energy_valid = true)
        -> EnergyEventAnalysis

Split the run into its adjustment intervals, one [`EnergyInterval`](@ref) per
ADJUST record `k ≥ 2` of `diag.adjust`, with the increments of the engine's
correction channels and the counts of the `events` whose `interval` is `k`
(events of other intervals or of an unknown kind are ignored).

An interval is flagged when

    |DE| ≥ max(threshold, sigma × 1.4826 × MAD),

the MAD taken over the finite `DE` of all intervals of the run. The floor
`threshold` keeps a run whose errors sit at round-off level unflagged; the
robust scale keeps the rule relative to the run's own noise, so a few large
jumps neither inflate it nor hide each other. With `engine_energy_valid =
false` (see [`engine_energy_valid`](@ref)) no interval is flagged. Throws an
`ArgumentError` for a negative `threshold` or `sigma`.

# Example
```julia
a = energy_intervals(diag, events; engine_energy_valid = engine_energy_valid(tidal_mode(diag)))
energy_events_summary(a)
```
"""
function energy_intervals(
    diag::DiagnosticsData,
    events::AbstractVector{EngineEvent};
    threshold::Real = 1.0e-4,
    sigma::Real = 5.0,
    engine_energy_valid::Bool = true,
)::EnergyEventAnalysis
    threshold >= 0 || throw(ArgumentError("threshold must be ≥ 0, got $threshold"))
    sigma >= 0 || throw(ArgumentError("sigma must be ≥ 0, got $sigma"))
    adj = diag.adjust
    nint = max(length(adj) - 1, 0)

    # Event counts per interval in one pass; column k − 1 holds interval k.
    nkinds = length(ENGINE_EVENT_KINDS)
    C = zeros(Int, nkinds, nint)
    for ev in events
        k = ev.interval
        2 <= k <= length(adj) || continue
        j = findfirst(==(ev.kind), ENGINE_EVENT_KINDS)
        j === nothing && continue
        C[j, k - 1] += 1
    end

    x = Float64[adj[k].de_rel for k in 2:length(adj) if isfinite(adj[k].de_rel)]
    scale = length(x) < 2 ? 0.0 : 1.4826 * _median(abs.(x .- _median(x)))
    level = max(Float64(threshold), Float64(sigma) * scale)

    intervals = Vector{EnergyInterval}(undef, nint)
    for i in 1:nint
        a, b = adj[i], adj[i + 1]
        counts = ntuple(j -> C[j, i], Val(16))
        de = b.de_rel
        intervals[i] = EnergyInterval(
            i + 1,
            a.time_nb,
            b.time_nb,
            de,
            b.de_abs,
            b.ebin - a.ebin,
            b.emerge - a.emerge,
            b.esub - a.esub,
            b.ecoll - a.ecoll,
            b.emdot - a.emdot,
            b.ecdot - a.ecdot,
            counts,
            sum(counts),
            scale > 0 ? abs(de) / scale : NaN,
            engine_energy_valid && isfinite(de) && abs(de) >= level,
        )
    end
    return EnergyEventAnalysis(
        intervals,
        collect(EngineEvent, events),
        Float64(threshold),
        Float64(sigma),
        scale,
        engine_energy_valid,
    )
end

"""Median of a non-empty real vector (mean of the two middle values for an even length)."""
function _median(v::AbstractVector{<:Real})::Float64
    isempty(v) && throw(ArgumentError("median of an empty vector"))
    s = sort(v)
    n = length(s)
    h = n ÷ 2
    return isodd(n) ? Float64(s[h + 1]) : 0.5 * (Float64(s[h]) + Float64(s[h + 1]))
end

"""Whether no correction channel of `iv` changed: every finite increment is zero."""
function _no_channel_increment(iv::EnergyInterval)::Bool
    return all(
        d -> !isfinite(d) || d == 0,
        (iv.d_ebin, iv.d_emerge, iv.d_esub, iv.d_ecoll, iv.d_emdot, iv.d_ecdot),
    )
end

"""
    energy_events_summary(a::EnergyEventAnalysis) -> String

One-line summary of [`energy_intervals`](@ref): the number of intervals and of
flagged ones with the flagging level and the robust scale, how many flagged
intervals contain events and how many contain neither an event nor a channel
increment (a jump the log does not explain), and the event counts of the
flagged intervals by kind, in the order of `ENGINE_EVENT_KINDS`. When the
engine's DE is not an error measure, only the number of intervals is given.
"""
function energy_events_summary(a::EnergyEventAnalysis)::String
    n = length(a.intervals)
    if !a.engine_energy_valid
        return "engine DE is not an error measure (linearised tidal field): no interval flagged; " *
               "$n adjustment intervals"
    end
    level = max(a.threshold, a.sigma * a.scale)
    m = p = q = 0
    totals = zeros(Int, length(ENGINE_EVENT_KINDS))
    for iv in a.intervals
        iv.flagged || continue
        m += 1
        if iv.n_events > 0
            p += 1
            totals .+= iv.counts
        elseif _no_channel_increment(iv)
            q += 1
        end
    end
    kinds = ["$(ENGINE_EVENT_KINDS[j])=$(totals[j])" for j in eachindex(totals) if totals[j] > 0]
    listed = isempty(kinds) ? "none" : join(kinds, ", ")
    return "$n adjustment intervals, $m flagged (|DE| ≥ $(@sprintf("%.3g", level)); " *
           "robust scale $(@sprintf("%.3g", a.scale))); $p flagged with events, " *
           "$q flagged with neither an event nor a channel increment; " *
           "events in flagged intervals: $listed"
end

"""
    write_energy_intervals(path, a::EnergyEventAnalysis) -> String

Write the intervals of [`energy_intervals`](@ref) as CSV, one row per
interval: `index`, `t_start`, `t_end`, `de_rel`, `de_abs`, the channel
increments `d_*`, `score`, `flagged`, `n_events` and one `n_<kind>` column per
kind of `ENGINE_EVENT_KINDS` (N-body units; non-finite values as empty cells).
An existing file is backed up, never overwritten. Returns the path.
"""
function write_energy_intervals(path::AbstractString, a::EnergyEventAnalysis)::String
    _backup_existing(path)
    open(path, "w") do io
        header = (
            "index",
            "t_start",
            "t_end",
            "de_rel",
            "de_abs",
            "d_ebin",
            "d_emerge",
            "d_esub",
            "d_ecoll",
            "d_emdot",
            "d_ecdot",
            "score",
            "flagged",
            "n_events",
            map(k -> "n_" * String(k), ENGINE_EVENT_KINDS)...,
        )
        println(io, join(header, ","))
        for iv in a.intervals
            cells = (
                _csv_cell(iv.index),
                _csv_precise(iv.t_start),
                _csv_precise(iv.t_end),
                _csv_precise(iv.de_rel),
                _csv_precise(iv.de_abs),
                _csv_precise(iv.d_ebin),
                _csv_precise(iv.d_emerge),
                _csv_precise(iv.d_esub),
                _csv_precise(iv.d_ecoll),
                _csv_precise(iv.d_emdot),
                _csv_precise(iv.d_ecdot),
                _csv_precise(iv.score),
                iv.flagged ? "true" : "false",
                _csv_cell(iv.n_events),
                map(_csv_cell, iv.counts)...,
            )
            println(io, join(cells, ","))
        end
    end
    return String(path)
end

"""
    write_engine_events(path, events::AbstractVector{EngineEvent}) -> String

Write the stamped event lines (see [`EngineEvent`](@ref)) as CSV, one row per
event, in the field order of the struct; a stellar type the line does not
carry and non-finite values are written as empty cells. An existing file is
backed up, never overwritten. Returns the path.
"""
function write_engine_events(path::AbstractString, events::AbstractVector{EngineEvent})::String
    _backup_existing(path)
    kstar(k::Int) = k == KSTAR_ABSENT ? "" : _csv_cell(k)
    open(path, "w") do io
        println(
            io,
            "line,interval,time_nb,kind,name1,name2,name3,kstar1,kstar2,kstar3,count,n_total," *
            "energy_nb,mass_solar,dmass_solar,velocity_kms",
        )
        for ev in events
            cells = (
                _csv_cell(ev.line),
                _csv_cell(ev.interval),
                _csv_precise(ev.time_nb),
                _csv_cell(String(ev.kind)),
                _csv_cell(ev.name1),
                _csv_cell(ev.name2),
                _csv_cell(ev.name3),
                kstar(ev.kstar1),
                kstar(ev.kstar2),
                kstar(ev.kstar3),
                _csv_cell(ev.count),
                _csv_cell(ev.n_total),
                _csv_precise(ev.energy_nb),
                _csv_precise(ev.mass_solar),
                _csv_precise(ev.dmass_solar),
                _csv_precise(ev.velocity_kms),
            )
            println(io, join(cells, ","))
        end
    end
    return String(path)
end
