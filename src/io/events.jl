# =============================================================================
# Readers for the stamped event lines of the engine log (out1000) and for the
# event counters of event.35
# =============================================================================
#
# Each event kind and the engine routine that writes its line:
#
#   kind                    line prefix                  writer
#   :ks_rectification       BEG KSRECT TIME[NB]          ksrect.f
#   :ks_new                 NEW KSREG                    ksinit.F
#   :ks_new_compact         NS/BH BINARY T[NB]           ksinit.F
#   :ks_new_mdot            MDNEW KSREG TIME[NB]         mdot.F
#   :ks_end_mdot            MDEND KSREG TIME[NB]         mdot.F
#   :remnant_formation      WD/NS/BH FORMATION T=        mdot.F
#   :kick_full              FULL VELOCITY KICK:          kick.F
#   :kick_zero              ZERO VELOCITY KICK:          kick.F
#   :kick_low               LOW VELOCITY KICK:           kick.F
#   :collision              <word> COLLISION NSYS        cmbody.f
#   :merged_star            NEW STAR TIME[NB]            mix.f
#   :coalescence            <word> COAL: IQCOLL          coal.f
#   :escape                 ESCAPE T,NESC =              escape.F
#   :supernova_no_remnant   MDOT WARNING! SN KW=15       mdot.F
#   :chain_start            NEW CHAIN                    chain.f
#   :chain_end              CHTERM:                      chterm.f
#
# Lines without a stamp of their own (ENFORCED MASS LOSS, MIX TYPE CHANGE,
# WD/NS BINARY, BINARY COAL, BINARY/QUAD/HIARCH ESCAPE, RECTIFY K*, WARNING!,
# NEW TZ) are not read; the stamped line of the same event is.
#
# Interval rule: the ADJUST routine evaluates the energy, then removes the
# escapers, then writes its output. Every line between two ADJUST lines
# therefore belongs to the interval that ends at the second, and the ESCAPE
# line, printed after the first, carries the first's time stamp.
#
# event.35 (events.f) holds one row of cumulative counters per output time,
# preceded by a TIME[Myr] header row when the run started fresh.

"""
    read_engine_events(path::AbstractString) -> Vector{EngineEvent}

Parse the stamped event lines of the engine log (`out1000`) and return one
[`EngineEvent`](@ref) per line, in file order. The lines are written by
`ksrect.f`, `ksinit.F`, `mdot.F`, `kick.F`, `cmbody.f`, `mix.f`, `coal.f`,
`escape.F`, `chain.f` and `chterm.f`; the kinds are listed in
[`ENGINE_EVENT_KINDS`](@ref) and the fields each kind carries in the
docstring of `EngineEvent`. `interval` is one more than the number of
`ADJUST:` lines above the event, so a line printed between the (k−1)-th and
the k-th ADJUST line has `interval == k`; `line` is the 1-based line number.
Lines that do not parse (overflowed Fortran fields, glued exponents such as
`1.08726-157`) are skipped.
"""
function read_engine_events(path::AbstractString)::Vector{EngineEvent}
    isfile(path) || error("Diagnostics file not found: $path")

    events = EngineEvent[]
    n_adjust = 0
    line = 0

    for text in eachline(path)
        line += 1
        s = lstrip(text)

        if startswith(s, "ADJUST:")
            n_adjust += 1
            continue
        end

        kind = _engine_event_kind(s)
        kind === :none && continue

        try
            push!(events, _parse_engine_event(kind, split(s), n_adjust + 1, line))
        catch e
            @debug "Skipping unparseable event line" line exception = e
        end
    end

    return events
end

"""
    read_event_counters(path::AbstractString) -> Vector{EventCounters}

Parse `event.35`, the engine's event counters (`events.f`), and return one
[`EventCounters`](@ref) per output time. Each data row holds 58 tokens: the
time in Myr, 13 interaction counters, 12 energy terms and the maximum mass,
9 stellar-evolution counters, 7 cumulative mass-loss terms and the 16
stellar-type counts. The `TIME[Myr]` header row (absent after a restart),
blank lines, and rows that are short or do not parse are skipped.
"""
function read_event_counters(path::AbstractString)::Vector{EventCounters}
    isfile(path) || error("Event counters file not found: $path")

    records = EventCounters[]

    for text in eachline(path)
        tokens = split(text)
        isempty(tokens) && continue
        tokens[1] == "TIME[Myr]" && continue
        if length(tokens) < 58
            @debug "Skipping short event.35 row" text ntokens = length(tokens)
            continue
        end
        try
            push!(records, _parse_event_counters(tokens))
        catch e
            @debug "Skipping unparseable event.35 row" text exception = e
        end
    end

    return records
end

# ---------------------------------------------------------------------------
# Internal parsers
# ---------------------------------------------------------------------------

"""Index of the first token equal to `label`; `ArgumentError` when absent."""
function _label_index(tokens::AbstractVector{<:AbstractString}, label::AbstractString)::Int
    i = findfirst(==(label), tokens)
    isnothing(i) && throw(ArgumentError("Label $label absent from event line"))
    return i
end

"""The `n` tokens following the first token equal to `label`."""
function _after(tokens::AbstractVector{<:AbstractString}, label::AbstractString, n::Integer)
    i = _label_index(tokens, label)
    i + n <= length(tokens) || throw(ArgumentError("Fewer than $n tokens follow $label"))
    return view(tokens, (i + 1):(i + n))
end

"""The `n` tokens following the `=` that follows the first token equal to `label`."""
function _after_equals(tokens::AbstractVector{<:AbstractString}, label::AbstractString, n::Integer)
    i = _label_index(tokens, label) + 1
    i <= length(tokens) && tokens[i] == "=" ||
        throw(ArgumentError("No `=` follows $label in event line"))
    i + n <= length(tokens) || throw(ArgumentError("Fewer than $n tokens follow $label ="))
    return view(tokens, (i + 1):(i + n))
end

"""Event kind of a stripped log line from its prefix; `:none` for other lines."""
function _engine_event_kind(s::AbstractString)::Symbol
    startswith(s, "BEG KSRECT TIME[NB]") && return :ks_rectification
    startswith(s, "NEW KSREG") && return :ks_new
    startswith(s, "NS/BH BINARY T[NB]") && return :ks_new_compact
    startswith(s, "MDNEW KSREG TIME[NB]") && return :ks_new_mdot
    startswith(s, "MDEND KSREG TIME[NB]") && return :ks_end_mdot
    startswith(s, "WD/NS/BH FORMATION T=") && return :remnant_formation
    startswith(s, "FULL VELOCITY KICK:") && return :kick_full
    startswith(s, "ZERO VELOCITY KICK:") && return :kick_zero
    startswith(s, "LOW VELOCITY KICK:") && return :kick_low
    occursin(r"^\S+\s+COLLISION\s+NSYS", s) && return :collision
    startswith(s, "NEW STAR TIME[NB]") && return :merged_star
    occursin(r"^\S+\s+COAL:\s+IQCOLL", s) && return :coalescence
    startswith(s, "ESCAPE T,NESC =") && return :escape
    startswith(s, "MDOT WARNING! SN KW=15") && return :supernova_no_remnant
    startswith(s, "NEW CHAIN") && return :chain_start
    startswith(s, "CHTERM:") && return :chain_end
    return :none
end

"""Build the `EngineEvent` of a line of known `kind` from its tokens."""
function _parse_engine_event(
    kind::Symbol,
    tokens::AbstractVector{<:AbstractString},
    interval::Int,
    line::Int,
)::EngineEvent
    if kind === :ks_rectification ||
       kind === :ks_new ||
       kind === :ks_new_compact ||
       kind === :ks_new_mdot ||
       kind === :ks_end_mdot
        time_label = kind === :ks_new_compact ? "T[NB]" : "TIME[NB]"
        nm = _after(tokens, "NM1,2,S=", 3)
        kw = _after(tokens, "KW1,2,S=", 3)
        m = _after(tokens, "M1,2[*]", 2)
        return EngineEvent(;
            time_nb = _parse_fortran_float(_after(tokens, time_label, 1)[1]),
            kind = kind,
            interval = interval,
            line = line,
            name1 = parse(Int, nm[1]),
            name2 = parse(Int, nm[2]),
            name3 = parse(Int, nm[3]),
            kstar1 = parse(Int, kw[1]),
            kstar2 = parse(Int, kw[2]),
            kstar3 = parse(Int, kw[3]),
            energy_nb = _parse_fortran_float(_after(tokens, "e,a,eb[NB]=", 3)[3]),
            mass_solar = _parse_fortran_float(m[1]) + _parse_fortran_float(m[2]),
        )
    elseif kind === :remnant_formation
        ik = _after(tokens, "I,K,N,KW=", 5)
        m = _after(tokens, "M0/1[*]", 2)
        m0 = _parse_fortran_float(m[1])
        m1 = _parse_fortran_float(m[2])
        return EngineEvent(;
            time_nb = _parse_fortran_float(_after(tokens, "T=", 1)[1]),
            kind = kind,
            interval = interval,
            line = line,
            name1 = parse(Int, ik[3]),
            kstar1 = parse(Int, ik[4]),
            kstar2 = parse(Int, ik[5]),
            mass_solar = m1,
            dmass_solar = m0 - m1,
        )
    elseif kind === :kick_full || kind === :kick_zero || kind === :kick_low
        ids = _after(tokens, "I,N=", 2)
        ks = _after(tokens, "ICASE", 5)
        mn = _parse_fortran_float(_after(tokens, "MN[M*]", 1)[1])
        m0 = _parse_fortran_float(_after(tokens, "M0[M*]", 1)[1])
        return EngineEvent(;
            time_nb = _parse_fortran_float(_after(tokens, "TIME[NB]", 1)[1]),
            kind = kind,
            interval = interval,
            line = line,
            name1 = parse(Int, ids[2]),
            kstar1 = parse(Int, ks[1]),
            kstar2 = parse(Int, ks[2]),
            kstar3 = parse(Int, ks[3]),
            mass_solar = mn,
            dmass_solar = m0 - mn,
            velocity_kms = _parse_fortran_float(_after(tokens, "VK[km/s]", 1)[1]),
        )
    elseif kind === :collision
        nm = _after(tokens, "NM1,2,S=", 3)
        kw = _after(tokens, "KW1,2,S=", 3)
        return EngineEvent(;
            time_nb = _parse_fortran_float(_after(tokens, "TIME[NB]", 1)[1]),
            kind = kind,
            interval = interval,
            line = line,
            name1 = parse(Int, nm[1]),
            name2 = parse(Int, nm[2]),
            name3 = parse(Int, nm[3]),
            kstar1 = parse(Int, kw[1]),
            kstar2 = parse(Int, kw[2]),
            kstar3 = parse(Int, kw[3]),
            count = parse(Int, _after_equals(tokens, "NSYS", 1)[1]),
            energy_nb = _parse_fortran_float(_after(tokens, "e,a,eb[NB],P[d]=", 4)[3]),
            mass_solar = _parse_fortran_float(_after(tokens, "M1,2,3[*]", 3)[3]),
        )
    elseif kind === :merged_star
        nm = _after(tokens, "NM1,2=", 2)
        return EngineEvent(;
            time_nb = _parse_fortran_float(_after(tokens, "TIME[NB]", 1)[1]),
            kind = kind,
            interval = interval,
            line = line,
            name1 = parse(Int, nm[1]),
            name2 = parse(Int, nm[2]),
            kstar1 = parse(Int, _after(tokens, "K*(I1)", 1)[1]),
            kstar2 = parse(Int, _after(tokens, "K*(I2)", 1)[1]),
            kstar3 = parse(Int, _after(tokens, "K*(INEW)", 1)[1]),
            mass_solar = _parse_fortran_float(_after(tokens, "M(INEW)[M*]", 1)[1]),
            dmass_solar = _parse_fortran_float(_after(tokens, "DM[M*]", 1)[1]),
        )
    elseif kind === :coalescence
        nn = _after(tokens, "N1,2", 2)
        kw = _after(tokens, "KW1,2,S", 3)
        mm = _after(tokens, "M12S,DM[*]", 4)
        return EngineEvent(;
            time_nb = _parse_fortran_float(_after(tokens, "TIME[NB]", 1)[1]),
            kind = kind,
            interval = interval,
            line = line,
            name1 = parse(Int, nn[1]),
            name2 = parse(Int, nn[2]),
            kstar1 = parse(Int, kw[1]),
            kstar2 = parse(Int, kw[2]),
            kstar3 = parse(Int, kw[3]),
            count = parse(Int, _after(tokens, "IQCOLL", 1)[1]),
            energy_nb = _parse_fortran_float(_after(tokens, "e,a,eb,dp[NB]=", 4)[3]),
            mass_solar = _parse_fortran_float(mm[3]),
            dmass_solar = _parse_fortran_float(mm[4]),
        )
    elseif kind === :escape
        tn = _after_equals(tokens, "T,NESC", 2)
        return EngineEvent(;
            time_nb = _parse_fortran_float(tn[1]),
            kind = kind,
            interval = interval,
            line = line,
            count = parse(Int, tn[2]),
            n_total = parse(Int, _after_equals(tokens, "N,NSESC,NBESC", 3)[1]),
            energy_nb = _parse_fortran_float(_after_equals(tokens, "NCRIT1", 1)[1]),
            mass_solar = _parse_fortran_float(_after(tokens, "M*", 1)[1]),
        )
    elseif kind === :supernova_no_remnant
        # Ten label tokens, then TTOT, I, NM, K* before, KW, M0, M1, MC, AGE, TIME.
        length(tokens) >= 17 ||
            throw(ArgumentError("Fewer than 17 tokens in MDOT WARNING! SN line"))
        m0 = _parse_fortran_float(tokens[16])
        m1 = _parse_fortran_float(tokens[17])
        return EngineEvent(;
            time_nb = _parse_fortran_float(tokens[11]),
            kind = kind,
            interval = interval,
            line = line,
            name1 = parse(Int, tokens[13]),
            kstar1 = parse(Int, tokens[14]),
            kstar2 = parse(Int, tokens[15]),
            mass_solar = m1,
            dmass_solar = m0 - m1,
        )
    elseif kind === :chain_start
        # After `UNIT]`: NCH, NP, ENERGY, RSUM, RGRAV, TCR, RMAXS, NAMEC(1:NM).
        v = _after(tokens, "UNIT]", 7)
        first_name = _label_index(tokens, "UNIT]") + 8
        names = ntuple(Val(3)) do k
            j = first_name + k - 1
            j <= length(tokens) ? parse(Int, tokens[j]) : 0
        end
        return EngineEvent(;
            time_nb = NaN,
            kind = kind,
            interval = interval,
            line = line,
            name1 = names[1],
            name2 = names[2],
            name3 = names[3],
            count = parse(Int, v[1]),
            energy_nb = _parse_fortran_float(v[3]),
        )
    elseif kind === :chain_end
        # After `STEP(I3,I4,ICM1)[NB]`: time, I3, I4, ICM1, NI3, NI4, NCM.
        v = _after(tokens, "STEP(I3,I4,ICM1)[NB]", 7)
        return EngineEvent(;
            time_nb = _parse_fortran_float(v[1]),
            kind = kind,
            interval = interval,
            line = line,
            name1 = parse(Int, v[5]),
            name2 = parse(Int, v[6]),
            name3 = parse(Int, v[7]),
        )
    end
    throw(ArgumentError("Unknown engine event kind $kind"))
end

"""Build an `EventCounters` from the 58 tokens of an event.35 data row."""
function _parse_event_counters(tokens::AbstractVector{<:AbstractString})::EventCounters
    f(k) = _parse_fortran_float(tokens[k])
    i(k) = parse(Int, tokens[k])
    return EventCounters(
        f(1),
        # ndiss … nkick
        i(2),
        i(3),
        i(4),
        i(5),
        i(6),
        i(7),
        i(8),
        i(9),
        i(10),
        i(11),
        i(12),
        i(13),
        i(14),
        # ebin … mmax
        f(15),
        f(16),
        f(17),
        f(18),
        f(19),
        f(20),
        f(21),
        f(22),
        f(23),
        f(24),
        f(25),
        f(26),
        # nmdot … nbs
        i(27),
        i(28),
        i(29),
        i(30),
        i(31),
        i(32),
        i(33),
        i(34),
        i(35),
        # zmrg … zmdot
        f(36),
        f(37),
        f(38),
        f(39),
        f(40),
        f(41),
        f(42),
        # ntype(1:16)
        ntuple(k -> i(42 + k), Val(16)),
    )
end
