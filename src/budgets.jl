# =============================================================================
# Mass and energy budgets of the snapshots
# =============================================================================
# The engine removes mass from the cluster by escaper removal, by stellar-
# evolution mass loss (accumulated in ZMDOT, M☉, written to event.35) and by
# the mass lost in collisions and coalescences, which leaves ZMASS but not
# ZMDOT (cmbody.f, coal.f). Escapers are removed before the output of the same
# adjustment and event.35 is written before the snapshot of the same output,
# so at a snapshot time t the stars of the snapshot, the escapers with
# t_esc ≤ t, the ZMDOT of the last output ≤ t and the collision losses up to t
# make up the initial mass.

"""
    MassBudgetRow

The mass budget of one snapshot, N-body mass units, from [`mass_budget`](@ref):
`m_snapshot` the stars' mass; `m_escaped` the mass of the escapers removed at or
before the snapshot time (`esc.11`, binary escapers as their centre of mass);
`m_lost_stellar` the stellar-evolution mass loss ZMDOT of the last `event.35`
row at or before it (0 when there is none); `m_lost_collision` the mass lost in
the collisions and coalescences up to it (the `DM` of the `NEW STAR` and `COAL`
lines); `m_engine` the engine's own total ZMASS after the last escape at or
before it (NaN when none); `residual = m0 − (m_snapshot + m_escaped +
m_lost_stellar + m_lost_collision)` with `m0` the first snapshot's mass, and
`residual_rel = residual / m0`.
"""
struct MassBudgetRow
    time_nb::Float64
    time_myr::Float64
    n_stars::Int
    m_snapshot::Float64
    m_escaped::Float64
    m_lost_stellar::Float64
    m_lost_collision::Float64
    m_engine::Float64
    residual::Float64
    residual_rel::Float64
end

"""
    mass_budget(ints::AbstractVector{SnapshotIntegrals}, escapers::AbstractVector{EscaperRecord},
                counters::AbstractVector{EventCounters}, events::AbstractVector{EngineEvent},
                units::UnitScaling; time_tolerance = 1.0e-5) -> Vector{MassBudgetRow}

Mass budget of every snapshot (see [`MassBudgetRow`](@ref)), with the mass of
the first snapshot as the initial mass `m0`. Escaper and counter times are
converted from Myr with `units.tscale`, masses from M☉ with `units.zmbar`; a
record stamped in Myr enters the budget of a snapshot at time `t` when its time
is `≤ t + time_tolerance × max(t, 1)` (N-body units), which absorbs the
six-digit rounding of the Myr stamps; the event lines carry exact N-body
stamps and enter when their time is `≤ t` to round-off (`1e-9 × max(t, 1)`).
Of `events` only the `:merged_star` and `:coalescence` lines (their mass loss
`dmass_solar`) and the `:escape` lines (the engine's total mass) are used.
Throws an `ArgumentError` without snapshots, for non-positive `units.zmbar`
or `units.tscale`, or a negative `time_tolerance`.
"""
function mass_budget(
    ints::AbstractVector{SnapshotIntegrals},
    escapers::AbstractVector{EscaperRecord},
    counters::AbstractVector{EventCounters},
    events::AbstractVector{EngineEvent},
    units::UnitScaling;
    time_tolerance::Real = 1.0e-5,
)::Vector{MassBudgetRow}
    isempty(ints) && throw(ArgumentError("mass_budget needs at least one snapshot"))
    units.zmbar > 0 || throw(ArgumentError("units.zmbar must be positive, got $(units.zmbar)"))
    units.tscale > 0 || throw(ArgumentError("units.tscale must be positive, got $(units.tscale)"))
    time_tolerance >= 0 || throw(ArgumentError("time_tolerance must be ≥ 0, got $time_tolerance"))
    zmbar, tscale = units.zmbar, units.tscale

    # Escapers sorted by time with the cumulative escaped mass.
    t_esc_raw = Float64[e.time_myr / tscale for e in escapers]
    p_esc = sortperm(t_esc_raw)
    t_esc = t_esc_raw[p_esc]
    m_esc_cum = cumsum(Float64[escapers[j].mass_solar / zmbar for j in p_esc])

    # event.35 rows sorted by time; ZMDOT is already cumulative.
    t_cnt_raw = Float64[c.time_myr / tscale for c in counters]
    p_cnt = sortperm(t_cnt_raw)
    t_cnt = t_cnt_raw[p_cnt]
    zmdot = Float64[counters[j].zmdot / zmbar for j in p_cnt]

    # Collision and coalescence losses with finite stamp and DM, cumulative.
    coll = EngineEvent[
        ev for ev in events if (ev.kind === :merged_star || ev.kind === :coalescence) &&
            isfinite(ev.time_nb) &&
            isfinite(ev.dmass_solar)
    ]
    p_coll = sortperm(Float64[ev.time_nb for ev in coll])
    t_coll = Float64[coll[j].time_nb for j in p_coll]
    m_coll_cum = cumsum(Float64[coll[j].dmass_solar / zmbar for j in p_coll])

    # Escape lines carry the engine's total mass after the removal.
    esc_ev = EngineEvent[ev for ev in events if ev.kind === :escape && isfinite(ev.time_nb)]
    p_ev = sortperm(Float64[ev.time_nb for ev in esc_ev])
    t_ev = Float64[esc_ev[j].time_nb for j in p_ev]
    m_ev = Float64[esc_ev[j].mass_solar / zmbar for j in p_ev]

    m0 = Float64(ints[1].mass)
    rows = Vector{MassBudgetRow}(undef, length(ints))
    for (i, s) in enumerate(ints)
        bound_myr, bound_nb = _budget_bounds(s.time_nb, time_tolerance)
        j = searchsortedlast(t_esc, bound_myr)
        m_escaped = j == 0 ? 0.0 : m_esc_cum[j]
        j = searchsortedlast(t_cnt, bound_myr)
        m_stellar = j == 0 ? 0.0 : zmdot[j]
        j = searchsortedlast(t_coll, bound_nb)
        m_collision = j == 0 ? 0.0 : m_coll_cum[j]
        j = searchsortedlast(t_ev, bound_nb)
        m_engine = j == 0 ? NaN : m_ev[j]
        residual = m0 - (s.mass + m_escaped + m_stellar + m_collision)
        rows[i] = MassBudgetRow(
            s.time_nb,
            s.time_myr,
            s.n_stars,
            s.mass,
            m_escaped,
            m_stellar,
            m_collision,
            m_engine,
            residual,
            m0 == 0 ? NaN : residual / m0,
        )
    end
    return rows
end

"""
    mass_budget_summary(rows::AbstractVector{MassBudgetRow}) -> String

One-line summary of [`mass_budget`](@ref): the largest |residual| relative to
the first snapshot's mass and its time, and the escaped, stellar-evolution and
collision losses of the last snapshot as percentages of that mass.
"""
function mass_budget_summary(rows::AbstractVector{MassBudgetRow})::String
    isempty(rows) && return "mass budget: no snapshots"
    m0 = rows[1].m_snapshot
    worst, t_worst = NaN, NaN
    for r in rows
        x = abs(r.residual) / m0
        if isfinite(x) && !(x <= worst)
            worst, t_worst = x, r.time_nb
        end
    end
    last_row = rows[end]
    f(v) = @sprintf("%.3g", v)
    return "mass budget over $(length(rows)) snapshots: max |residual|/m0 = $(f(worst)) " *
           "at t = $(f(t_worst)) NB; by the last snapshot escaped " *
           "$(f(100 * last_row.m_escaped / m0)) %, stellar mass loss " *
           "$(f(100 * last_row.m_lost_stellar / m0)) %, collisions " *
           "$(f(100 * last_row.m_lost_collision / m0)) % of m0"
end

"""
    write_mass_budget(path, rows::AbstractVector{MassBudgetRow}) -> String

Write the rows of [`mass_budget`](@ref) as CSV, one per snapshot (N-body units
plus `time_myr`; non-finite values as empty cells). An existing file is backed
up, never overwritten. Returns the path.
"""
function write_mass_budget(path::AbstractString, rows::AbstractVector{MassBudgetRow})::String
    _backup_existing(path)
    open(path, "w") do io
        println(
            io,
            "time_nb,time_myr,n_stars,m_snapshot,m_escaped,m_lost_stellar,m_lost_collision," *
            "m_engine,residual,residual_rel",
        )
        for r in rows
            cells = (
                _csv_precise(r.time_nb),
                _csv_precise(r.time_myr),
                _csv_cell(r.n_stars),
                _csv_precise(r.m_snapshot),
                _csv_precise(r.m_escaped),
                _csv_precise(r.m_lost_stellar),
                _csv_precise(r.m_lost_collision),
                _csv_precise(r.m_engine),
                _csv_precise(r.residual),
                _csv_precise(r.residual_rel),
            )
            println(io, join(cells, ","))
        end
    end
    return String(path)
end

"""Upper time bounds of a snapshot at `t` for Myr-stamped records (`tol` relative) and for exact N-body stamps."""
function _budget_bounds(t::Float64, tol::Real)::NTuple{2,Float64}
    return (t + Float64(tol) * max(t, 1.0), t + 1.0e-9 * max(t, 1.0))
end

# ---------------------------------------------------------------------------
# Energy budget
# ---------------------------------------------------------------------------
# The engine conserves ETOT = ZKIN − POT + ETIDE + EBIN + ESUB + EMERGE + ECOLL
# + EMDOT + ECDOT up to its integration error and subtracts the energy of every
# removed escaper from it (escape.F; the sums ESESC, EBESC, EMESC in event.35).
# A snapshot holds ZKIN − POT (with the tidal term as e_jacobi) and EBIN as the
# pairs' internal energy; the members of chains and hierarchies appear with
# their physical coordinates, so ESUB and EMERGE are inside the snapshot sums.
# Hence e_jacobi(t) + E_escaped(t) + ECOLL + EMDOT + ECDOT is the conserved
# total the engine checks, reconstructed from the snapshot.

"""
    EnergyBudgetRow

The energy budget of one snapshot, N-body units, from [`energy_budget`](@ref):
`e_snapshot` the Jacobi integral of the snapshot, `T_cm + W + e_pairs +
W_tide`; `e_pairs` the regularised pairs' internal energy inside it, from
`bev.82` at the snapshot times it covers (`pair_source = :bev`, the engine's
double-precision elements) and from the single-precision snapshot otherwise
(`:snapshot`, where the tightest pairs are quantised); `e_escaped` the energy
the engine subtracted for the escapers removed at or before the snapshot time
— `ESESC + EBESC + EMESC` of the last `event.35` row (`source = :counters`),
or ``Σ m\\,E_\\mathrm{esc}`` over the `esc.11` records when there is no
counters file (`:escapers`; the internal energy of escaping pairs is then
missing), or 0 (`:none`); `e_channels` the engine's `ECOLL + EMDOT + ECDOT`,
from the same `event.35` row with `:counters` and from the last ADJUST record
at or before the snapshot time otherwise (four printed digits; NaN when there
is none); `e_total = e_snapshot + e_escaped + e_channels`; `residual` its
change since the first snapshot and `residual_rel` that change over
`|e_total|` of the first snapshot; `detot` the engine's own cumulative energy
change `DETOT` at the last ADJUST record at or before the snapshot, which the
residual tracks when the reconstruction is exact.
"""
struct EnergyBudgetRow
    time_nb::Float64
    time_myr::Float64
    e_snapshot::Float64
    e_pairs::Float64
    e_escaped::Float64
    e_channels::Float64
    e_total::Float64
    residual::Float64
    residual_rel::Float64
    detot::Float64
    source::Symbol
    pair_source::Symbol
end

"""
    energy_budget(ints::AbstractVector{SnapshotIntegrals}, escapers::AbstractVector{EscaperRecord},
                  counters::AbstractVector{EventCounters}, diag::DiagnosticsData,
                  units::UnitScaling; bevs = BinaryEvolutionSnapshot[], time_tolerance = 1.0e-5)
        -> Vector{EnergyBudgetRow}

Energy budget of every snapshot (see [`EnergyBudgetRow`](@ref)): the Jacobi
integral of the snapshot plus the energy of the removed escapers plus the
engine's collision, mass-loss and kick corrections, which together are the
total the engine conserves. The escaper energies and the corrections come
from `counters` (`event.35`, full precision) when given, else from
`escapers` (`esc.11`) and the ADJUST records of `diag` (four digits); the
pairs' internal energy from `bevs` (`bev.82`, matched by the N-body stamp
of its records, the header time being rounded) at the snapshot times it
covers. Times in Myr are converted with `units.tscale`, masses with
`units.zmbar`, the `bev.82` energies with `units.zmbar × units.vstar²`;
time bounds as in [`mass_budget`](@ref). Throws an `ArgumentError` without
snapshots, for non-positive `units.zmbar`, `units.tscale` or `units.vstar`,
or a negative `time_tolerance`.
"""
function energy_budget(
    ints::AbstractVector{SnapshotIntegrals},
    escapers::AbstractVector{EscaperRecord},
    counters::AbstractVector{EventCounters},
    diag::DiagnosticsData,
    units::UnitScaling;
    bevs::AbstractVector{BinaryEvolutionSnapshot} = BinaryEvolutionSnapshot[],
    time_tolerance::Real = 1.0e-5,
)::Vector{EnergyBudgetRow}
    isempty(ints) && throw(ArgumentError("energy_budget needs at least one snapshot"))
    units.zmbar > 0 || throw(ArgumentError("units.zmbar must be positive, got $(units.zmbar)"))
    units.tscale > 0 || throw(ArgumentError("units.tscale must be positive, got $(units.tscale)"))
    units.vstar > 0 || throw(ArgumentError("units.vstar must be positive, got $(units.vstar)"))
    time_tolerance >= 0 || throw(ArgumentError("time_tolerance must be ≥ 0, got $time_tolerance"))
    zmbar, tscale = units.zmbar, units.tscale

    # Escaper energies and correction channels: event.35 at full precision
    # when present (the ADJUST line prints four digits, which bound the
    # closure once the kick energies exceed the total by orders of magnitude).
    source = !isempty(counters) ? :counters : (!isempty(escapers) ? :escapers : :none)
    if source === :counters
        p = sortperm(Float64[c.time_myr for c in counters])
        t_esc = Float64[counters[j].time_myr / tscale for j in p]
        e_esc = Float64[counters[j].esesc + counters[j].ebesc + counters[j].emesc for j in p]
        t_chan = t_esc
        e_chan = Float64[counters[j].ecoll + counters[j].emdot + counters[j].ecdot for j in p]
    else
        p = sortperm(Float64[e.time_myr for e in escapers])
        t_esc = Float64[escapers[j].time_myr / tscale for j in p]
        e_esc =
            cumsum(Float64[escapers[j].mass_solar / zmbar * escapers[j].escape_energy for j in p])
        q = sortperm(Float64[a.time_nb for a in diag.adjust])
        t_chan = Float64[diag.adjust[j].time_nb for j in q]
        e_chan =
            Float64[diag.adjust[j].ecoll + diag.adjust[j].emdot + diag.adjust[j].ecdot for j in q]
    end
    q = sortperm(Float64[a.time_nb for a in diag.adjust])
    t_adj = Float64[diag.adjust[j].time_nb for j in q]
    detot_adj = Float64[diag.adjust[j].detot for j in q]

    # Pair energies from bev.82 (double precision in the engine) at the
    # snapshot times it covers; the single-precision snapshot otherwise.
    # The bev.82 header time is TPHYS rounded to few digits; the records carry
    # the exact N-body stamp of the line.
    bev_t(b) = isempty(b.records) ? b.time_myr / tscale : b.records[1].time_nb
    pb = sortperm(Float64[bev_t(b) for b in bevs])
    t_bev = Float64[bev_t(bevs[j]) for j in pb]
    e_bev = Float64[
        -sum(binding_energy(r) for r in bevs[j].records; init = 0.0) / (zmbar * units.vstar^2)
        for j in pb
    ]

    rows = Vector{EnergyBudgetRow}(undef, length(ints))
    e_ref = NaN
    for (i, s) in enumerate(ints)
        bound_myr, _ = _budget_bounds(s.time_nb, time_tolerance)
        lower_myr = s.time_nb - (bound_myr - s.time_nb)
        j = searchsortedlast(t_esc, bound_myr)
        e_escaped = j == 0 ? 0.0 : e_esc[j]
        j = searchsortedlast(t_chan, bound_myr)
        e_channels = j == 0 ? NaN : e_chan[j]
        j = searchsortedlast(t_adj, bound_myr)
        detot = j == 0 ? NaN : detot_adj[j]
        j = searchsortedlast(t_bev, bound_myr)
        pair_source = j > 0 && t_bev[j] >= lower_myr ? :bev : :snapshot
        e_pairs = pair_source === :bev ? e_bev[j] : s.pair_internal
        e_snapshot = s.kinetic + s.potential + e_pairs + s.tidal_potential
        e_total = e_snapshot + e_escaped + e_channels
        i == 1 && (e_ref = e_total)
        residual = e_total - e_ref
        rows[i] = EnergyBudgetRow(
            s.time_nb,
            s.time_myr,
            e_snapshot,
            e_pairs,
            e_escaped,
            e_channels,
            e_total,
            residual,
            e_ref == 0 ? NaN : residual / abs(e_ref),
            detot,
            source,
            pair_source,
        )
    end
    return rows
end

"""
    energy_budget_summary(rows::AbstractVector{EnergyBudgetRow}) -> String

One-line summary of [`energy_budget`](@ref): the largest |residual| relative
to the first snapshot's total and its time, and the source of the escaper
energy.
"""
function energy_budget_summary(rows::AbstractVector{EnergyBudgetRow})::String
    isempty(rows) && return "energy budget: no snapshots"
    worst, t_worst = NaN, NaN
    for r in rows
        x = abs(r.residual_rel)
        if isfinite(x) && !(x <= worst)
            worst, t_worst = x, r.time_nb
        end
    end
    f(v) = @sprintf("%.3g", v)
    n_bev = count(r -> r.pair_source === :bev, rows)
    return "energy budget over $(length(rows)) snapshots: max |residual|/|E0| = $(f(worst)) " *
           "at t = $(f(t_worst)) NB; escaper energy and corrections from $(rows[1].source); " *
           "pair energies from bev.82 at $n_bev snapshots"
end

"""
    write_energy_budget(path, rows::AbstractVector{EnergyBudgetRow}) -> String

Write the rows of [`energy_budget`](@ref) as CSV, one per snapshot (N-body
units plus `time_myr`; non-finite values as empty cells). An existing file is
backed up, never overwritten. Returns the path.
"""
function write_energy_budget(path::AbstractString, rows::AbstractVector{EnergyBudgetRow})::String
    _backup_existing(path)
    open(path, "w") do io
        println(
            io,
            "time_nb,time_myr,e_snapshot,e_pairs,e_escaped,e_channels,e_total,residual," *
            "residual_rel,detot,source,pair_source",
        )
        for r in rows
            cells = (
                _csv_precise(r.time_nb),
                _csv_precise(r.time_myr),
                _csv_precise(r.e_snapshot),
                _csv_precise(r.e_pairs),
                _csv_precise(r.e_escaped),
                _csv_precise(r.e_channels),
                _csv_precise(r.e_total),
                _csv_precise(r.residual),
                _csv_precise(r.residual_rel),
                _csv_precise(r.detot),
                String(r.source),
                String(r.pair_source),
            )
            println(io, join(cells, ","))
        end
    end
    return String(path)
end
