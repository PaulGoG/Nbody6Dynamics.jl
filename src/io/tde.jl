# =============================================================================
# Tidal disruptions and plunges at the static black hole (engine KZ(44) = 1)
# =============================================================================
#
# The engine fork (v2026.07.13 on, mbh.F) writes one stamped line per event:
#   ' TDE        T[NB] NAME K* M[M*] R*[R*] Q[PC] RT[PC] BETA E[NB] MBH[M*] t name k m r q rt beta e mbh'
#   ' TDE PLUNGE T[NB] …  (same fields)'
#   ' TDE PARTIAL T[NB] NAME K* M[M*] R*[R*] Q[PC] RT[PC] BETA E[NB] DM[M*] t name k m r q rt beta e dm'

"""
    TDEEvent

One event at the static black hole, from [`read_tde_events`](@ref).

| field | meaning | unit |
|---|---|---|
| `kind` | `:full` (the star disrupted, removed), `:plunge` (pericentre inside 4 G M/c² above r_t, removed), `:partial` (the star loses `dm_msun` and survives) | — |
| `time_nb`, `time_myr` | event time (the pericentre passage within one step) | NB, Myr |
| `name`, `kstar` | body name and stellar type | — |
| `mass_msun` | stellar mass before the event | M☉ |
| `radius_rsun` | stellar radius (0 for point masses) | R☉ |
| `q_pc` | pericentre distance of the two-body orbit about the black hole | pc |
| `r_t_pc` | disruption radius of the star | pc |
| `beta` | penetration factor r_t/q | — |
| `e_nb` | energy booked by the engine (of the removed or stripped mass) | NB |
| `dm_msun` | mass removed from the cluster: the star's mass for `:full` and `:plunge`, the stripped mass for `:partial` | M☉ |
| `mbh_msun` | black-hole mass at the event; `NaN` on partial lines | M☉ |
"""
Base.@kwdef struct TDEEvent
    kind::Symbol
    time_nb::Float64
    time_myr::Float64 = NaN
    name::Int
    kstar::Int
    mass_msun::Float64
    radius_rsun::Float64
    q_pc::Float64
    r_t_pc::Float64
    beta::Float64
    e_nb::Float64
    dm_msun::Float64
    mbh_msun::Float64 = NaN
end

"""
    read_tde_events(path; tstar_myr = NaN) -> Vector{TDEEvent}

Events at the static black hole from the engine stdout `path`: every line
beginning `TDE` (full disruptions), `TDE PLUNGE` and `TDE PARTIAL`, in file
order, with `time_myr = time_nb · tstar_myr` when `tstar_myr` is given.
Lines whose numbers do not parse are skipped with a count in a warning.
Throws an `ErrorException` when the file does not exist.

# Example
```julia
events = read_tde_events(joinpath(run, "output", "out1000"); tstar_myr = 0.21)
count(e -> e.kind === :full, events)
```
"""
function read_tde_events(path::AbstractString; tstar_myr::Real = NaN)::Vector{TDEEvent}
    isfile(path) || error("read_tde_events: no file $path")
    out = TDEEvent[]
    bad = 0
    for ln in eachline(path)
        s = lstrip(ln)
        startswith(s, "TDE") || continue
        kind =
            startswith(s, "TDE PARTIAL") ? :partial : startswith(s, "TDE PLUNGE") ? :plunge : :full
        tok = split(s)
        k0 = findfirst(t -> tryparse(Float64, t) !== nothing, tok)
        vals =
            k0 === nothing ? Float64[] : [something(tryparse(Float64, t), NaN) for t in tok[k0:end]]
        if length(vals) != 10 || any(isnan, vals)
            bad += 1
            continue
        end
        t, name, kw, m, r, q, rt, β, e, last = vals
        push!(
            out,
            TDEEvent(;
                kind = kind,
                time_nb = t,
                time_myr = t * Float64(tstar_myr),
                name = round(Int, name),
                kstar = round(Int, kw),
                mass_msun = m,
                radius_rsun = r,
                q_pc = q,
                r_t_pc = rt,
                beta = β,
                e_nb = e,
                dm_msun = kind === :partial ? last : m,
                mbh_msun = kind === :partial ? NaN : last,
            ),
        )
    end
    bad > 0 && @warn "read_tde_events: $bad TDE lines of $path did not parse"
    return out
end

"""
    tde_summary(events::AbstractVector{TDEEvent}) -> NamedTuple

Counts of the events by kind (`n_full`, `n_partial`, `n_plunge`), the mass
removed from the cluster `dm_msun` (all kinds), the time span of the events
`t_first_myr`, `t_last_myr` and the event rate `rate_per_myr` (events over the
span; `NaN` for fewer than two events).
"""
function tde_summary(events::AbstractVector{TDEEvent})
    n = length(events)
    ts = Float64[e.time_myr for e in events]
    t0 = n > 0 ? minimum(ts) : NaN
    t1 = n > 0 ? maximum(ts) : NaN
    return (
        n_full = count(e -> e.kind === :full, events),
        n_partial = count(e -> e.kind === :partial, events),
        n_plunge = count(e -> e.kind === :plunge, events),
        dm_msun = sum(e -> e.dm_msun, events; init = 0.0),
        t_first_myr = t0,
        t_last_myr = t1,
        rate_per_myr = (n ≥ 2 && t1 > t0) ? (n - 1) / (t1 - t0) : NaN,
    )
end

"""
    write_tde_events(path, events::AbstractVector{TDEEvent}) -> String

Write the events as CSV, one column per field of [`TDEEvent`](@ref) in
declaration order (`NaN` as empty cells, Symbols by name), after a `#` line
with the counts of [`tde_summary`](@ref). An existing file is backed up, never
overwritten. Returns the path.
"""
function write_tde_events(path::AbstractString, events::AbstractVector{TDEEvent})::String
    s = tde_summary(events)
    _backup_existing(path)
    open(path, "w") do io
        println(
            io,
            "# n_full=$(s.n_full) n_partial=$(s.n_partial) n_plunge=$(s.n_plunge) " *
            "dm_msun=$(s.dm_msun)",
        )
        println(io, join(fieldnames(TDEEvent), ","))
        for e in events
            cells = (
                v isa AbstractFloat ? (isfinite(v) ? repr(v) : "") : string(v) for
                v in (getfield(e, f) for f in fieldnames(TDEEvent))
            )
            println(io, join(cells, ","))
        end
    end
    return String(path)
end

"""
    tde_rate_series(events::AbstractVector{TDEEvent}; bin_myr, t_end_myr = NaN)
        -> Vector{NamedTuple}

Events per time bin of width `bin_myr` [Myr] from 0 to `t_end_myr` (default:
the last event): one row per bin with `t_start_myr`, `t_end_myr`, the counts
`n_full`, `n_partial`, `n_plunge`, the mass removed `dm_msun` and the rate of
all kinds `rate_per_myr` (events per Myr). Throws an `ArgumentError` for a
non-positive `bin_myr` or events without a physical time (`time_myr` NaN).

# Example
```julia
rates = tde_rate_series(events; bin_myr = 1.0)
```
"""
function tde_rate_series(
    events::AbstractVector{TDEEvent};
    bin_myr::Real,
    t_end_myr::Real = NaN,
)::Vector{NamedTuple}
    (isfinite(bin_myr) && bin_myr > 0) ||
        throw(ArgumentError("tde_rate_series: bin_myr must be finite and > 0; got $bin_myr"))
    any(e -> !isfinite(e.time_myr), events) &&
        throw(ArgumentError("tde_rate_series: events without time_myr (read with tstar_myr)"))
    t_end =
        isfinite(t_end_myr) ? Float64(t_end_myr) : maximum((e.time_myr for e in events); init = 0.0)
    nbins = max(1, ceil(Int, t_end / bin_myr))
    rows = NamedTuple[]
    for k in 1:nbins
        t0 = (k - 1) * bin_myr
        t1 = k * bin_myr
        inbin = [e for e in events if (t0 ≤ e.time_myr < t1) || (k == nbins && e.time_myr == t1)]
        n = length(inbin)
        push!(
            rows,
            (
                t_start_myr = t0,
                t_end_myr = t1,
                n_full = count(e -> e.kind === :full, inbin),
                n_partial = count(e -> e.kind === :partial, inbin),
                n_plunge = count(e -> e.kind === :plunge, inbin),
                dm_msun = sum(e -> e.dm_msun, inbin; init = 0.0),
                rate_per_myr = n / bin_myr,
            ),
        )
    end
    return rows
end
