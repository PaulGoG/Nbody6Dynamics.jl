# =============================================================================
# One pass over the engine log for the binary, few-body and GW event records
# =============================================================================
#
# The line families are parsed in events_dynamics.jl (KS episodes, chains and
# other few-body regularisations, hierarchies, exchanges, termination
# diagnostics) and events_evolution.jl (gravitational radiation, escapes of
# binaries, mass-loss terminations and kicks, tidal capture and
# circularisation, Roche-lobe overflow). The layouts are those of the fork's
# FORMAT statements; the hazards each parser corrects are in the docstrings of
# the record types.

"""
    BinaryEventLog

Records of the binary, few-body and gravitational-wave lines of one engine
log, by family, each vector in file order:

| Field | Records | Lines |
|---|---|---|
| `ks` | [`KSEvent`](@ref) | NEW/END KSREG, NS/BH BINARY |
| `fewbody` | [`FewBodyEvent`](@ref) | NEW TRIPLE/QUAD/CHAIN, END CHAIN, CHTERM, CHAIN B-B/ESCAPE/BINARY |
| `hierarchy` | [`HierarchyEvent`](@ref) | NEW/END MERGER and variants |
| `exchange` | [`ExchangeEvent`](@ref) | EXCHANGE |
| `termination` | [`TerminationEvent`](@ref) | END RUN DEFERRED, ENFORCED CHAIN TINY STEPS |
| `gw` | [`GWEvent`](@ref) | GR KSTIDE/KSINT/TIDES/BRAKE, GW mergers, GW recoil, GR CAPTURE/FLYBY |
| `escape` | [`BinaryEscapeEvent`](@ref) | BINARY/HIARCH/QUAD ESCAPE |
| `mass_loss` | [`MassLossEvent`](@ref) | MDEND/MDNEW KSREG, BINARY KICK |
| `tidal` | [`TidalEvent`](@ref) | NEW/END CIRC, CIRCULARIZED, NEW CAPTURE |
| `roche` | [`RocheEvent`](@ref) | NEW/END ROCHE, ROCHE COAL |

`skipped` counts, by kind, the lines of a known kind that did not parse.
"""
struct BinaryEventLog
    ks::Vector{KSEvent}
    fewbody::Vector{FewBodyEvent}
    hierarchy::Vector{HierarchyEvent}
    exchange::Vector{ExchangeEvent}
    termination::Vector{TerminationEvent}
    gw::Vector{GWEvent}
    escape::Vector{BinaryEscapeEvent}
    mass_loss::Vector{MassLossEvent}
    tidal::Vector{TidalEvent}
    roche::Vector{RocheEvent}
    skipped::Dict{Symbol,Int}
end

BinaryEventLog() = BinaryEventLog(
    KSEvent[],
    FewBodyEvent[],
    HierarchyEvent[],
    ExchangeEvent[],
    TerminationEvent[],
    GWEvent[],
    BinaryEscapeEvent[],
    MassLossEvent[],
    TidalEvent[],
    RocheEvent[],
    Dict{Symbol,Int}(),
)

_store!(log::BinaryEventLog, r::KSEvent) = push!(log.ks, r)
_store!(log::BinaryEventLog, r::FewBodyEvent) = push!(log.fewbody, r)
_store!(log::BinaryEventLog, r::HierarchyEvent) = push!(log.hierarchy, r)
_store!(log::BinaryEventLog, r::ExchangeEvent) = push!(log.exchange, r)
_store!(log::BinaryEventLog, r::TerminationEvent) = push!(log.termination, r)
_store!(log::BinaryEventLog, r::GWEvent) = push!(log.gw, r)
_store!(log::BinaryEventLog, r::BinaryEscapeEvent) = push!(log.escape, r)
_store!(log::BinaryEventLog, r::MassLossEvent) = push!(log.mass_loss, r)
_store!(log::BinaryEventLog, r::TidalEvent) = push!(log.tidal, r)
_store!(log::BinaryEventLog, r::RocheEvent) = push!(log.roche, r)

"""
    read_binary_events(path::AbstractString) -> BinaryEventLog

Read the binary, few-body and gravitational-wave event lines of an engine log
(`out1000`) in one pass and return them by family ([`BinaryEventLog`](@ref)).
Every record carries its line number, its adjustment interval (one more than
the number of `ADJUST:` lines above it, the convention of
[`read_engine_events`](@ref)), its time in N-body units and its kind. A line
of a known kind that does not parse is counted in `skipped` and reported in
one warning. The capture is held in memory while it is read (about 3 MiB per
N-body time unit at 2×10⁴ bodies with `KZ(10) = 2`).
"""
function read_binary_events(path::AbstractString)::BinaryEventLog
    isfile(path) || error("Engine log not found: $path")
    lines = readlines(path)
    log = BinaryEventLog()
    n_adjust = 0
    i = 1
    while i <= length(lines)
        s = lstrip(lines[i])
        if startswith(s, "ADJUST:")
            n_adjust += 1
            i += 1
            continue
        end
        family = :dynamics
        kind = _dynamics_line_kind(s)
        if kind === :none
            family = :evolution
            kind = _evolution_line_kind(s)
        end
        if kind === :none
            i += 1
            continue
        end
        consumed = 1
        try
            record, consumed = if family === :dynamics
                _parse_dynamics_line(kind, lines, i, n_adjust + 1)
            else
                _parse_evolution_line(kind, lines, i, n_adjust + 1)
            end
            _store!(log, record)
        catch e
            e isa Union{ArgumentError,BoundsError,InexactError} || rethrow()
            log.skipped[kind] = get(log.skipped, kind, 0) + 1
            consumed = 1
        end
        i += max(consumed, 1)
    end
    n_skipped = sum(values(log.skipped); init = 0)
    n_skipped > 0 &&
        @warn "$(basename(path)): $n_skipped binary-event line(s) did not parse" skipped =
            log.skipped
    return log
end
