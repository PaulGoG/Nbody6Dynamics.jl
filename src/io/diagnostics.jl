# =============================================================================
# Parser for Nbody6++ simulation stdout (out1000)
# =============================================================================
#
# Key patterns extracted (three related lines per epoch):
#   ADJUST:  TIME val  T[Myr] val  Q val  DE val  DELTA val  DETOT val  … ETOT val …
#   RMIN = …  RSCALE = …
#   TIME[NB]  val  N  val  NPAIRS  val  …
#
# The ADJUST line also carries the energy terms EKIN, POT, ETIDE, EBIN,
# EMERGE, ESUB, ECOLL, EMDOT, ECDOT, which are kept in the record.
#
# Also:
#   PHYSICAL SCALING:  R* = …  M* = …  V* = …  T* = …
#   OPTIONS, then the row of indices 1 … 50 and the row of the 50 KZ values
#   (written once at startup).
#   TIDAL PARAMETERS:  TIDAL(1:4)  TSCALE = …  RTIDE = …  (written once when KZ(14) ≤ 2)

"""
    read_diagnostics(path::AbstractString) -> DiagnosticsData

Parse the simulation stdout file for ADJUST lines, supplementary RSCALE
and TIME[NB] lines, and physical scaling info; the TIME[NB] lines give N,
NPAIRS and NMERGE.

Each ADJUST epoch in Nbody6++ output consists of up to three lines:
1. `ADJUST:` — energy, virial ratio, time
2. `RMIN = … RSCALE = …` — scale radius (may appear between ADJUST and TIME[NB])
3. `TIME[NB] … N … NPAIRS …` — particle counts

These are merged into a single `AdjustRecord` per epoch. An ADJUST line that
does not parse (too few fields, an overflowed `****` field) is skipped
together with its supplementary lines; the count is reported in one warning
per file.

The energy terms the key-value ADJUST line prints after `DETOT` (`EKIN`,
`POT`, `ETIDE`, `EBIN`, `EMERGE`, `ESUB`, `ECOLL`, `EMDOT`, `ECDOT`) are kept
in the record, NaN when the line lacks them.
The first `OPTIONS` echo gives `kz`, the 50 option values, left empty when
the capture lacks the echo or its value row does not hold 50 integers.
The `TIDAL PARAMETERS:` echo gives `tidal`, the four coefficients
`TIDAL(1:4)` in N-body units, zeros when the capture lacks it.
`nbin0` from the `&INDATA` echo, `nzero` (N) from the value row under the
`N  NFIX  NCRIT` header of the `&ININPUT` echo (`-1` when absent).
"""
function read_diagnostics(path::AbstractString)::DiagnosticsData
    isfile(path) || error("Diagnostics file not found: $path")

    adjust_records = AdjustRecord[]
    scaling = Dict{String,Float64}()

    # Accumulator for the current epoch being assembled
    pending_time_nb = NaN
    pending_time_myr = NaN
    pending_qvir = NaN
    pending_de_rel = NaN
    pending_e_tot = NaN
    pending_n = 0
    pending_npairs = 0
    pending_nmerge = 0
    pending_rscale = 0.0
    pending_de_abs = NaN
    pending_detot = NaN
    pending_ekin = NaN
    pending_pot = NaN
    pending_etide = NaN
    pending_ebin = NaN
    pending_emerge = NaN
    pending_esub = NaN
    pending_ecoll = NaN
    pending_emdot = NaN
    pending_ecdot = NaN
    have_adjust = false   # true once we've seen ADJUST: for this epoch
    n_unparsed = 0        # ADJUST lines skipped as unparseable

    kz = Int[]
    tidal = (0.0, 0.0, 0.0, 0.0)
    options_seen = false
    options_pending = 0   # non-empty lines after OPTIONS still to read; the second holds the values
    nbin0 = -1
    nzero = -1
    nzero_pending = false   # the N NFIX NCRIT header was seen, its value row is next

    for line in eachline(path)
        stripped = lstrip(line)

        # ── N of the &ININPUT echo: first non-empty row after its header ──
        if nzero_pending && !isempty(stripped)
            m = match(r"^(\d+)\s", stripped)
            m === nothing || (nzero = parse(Int, something(m.captures[1])))
            nzero_pending = false
        elseif nzero < 0 && startswith(stripped, "N  NFIX  NCRIT")
            nzero_pending = true
            continue
        end

        # ── NBIN0 of the &INDATA echo; no `continue`, the line may serve other branches ──
        if nbin0 < 0
            m = match(r"NBIN0\s*=\s*(\d+)", stripped)
            m === nothing || (nbin0 = parse(Int, something(m.captures[1])))
        end

        # ── OPTIONS echo: index row, then value row ──
        if options_pending > 0
            isempty(stripped) && continue
            options_pending -= 1
            if options_pending == 0
                _parse_options_row!(kz, stripped)
            end
            continue
        end
        if !options_seen && strip(stripped) == "OPTIONS"
            options_seen = true
            options_pending = 2
            continue
        end

        # ── TIDAL PARAMETERS echo ──
        if occursin("TIDAL PARAMETERS:", stripped)
            tidal = _parse_tidal_line(stripped)
            continue
        end

        # ── ADJUST line: start a new epoch ──
        if startswith(stripped, "ADJUST:")
            # Flush any previous pending epoch
            if have_adjust
                push!(
                    adjust_records,
                    AdjustRecord(
                        pending_time_nb,
                        pending_time_myr,
                        pending_qvir,
                        pending_de_rel,
                        pending_e_tot,
                        pending_n,
                        pending_npairs,
                        pending_rscale,
                        pending_de_abs,
                        pending_detot,
                        pending_ekin,
                        pending_pot,
                        pending_etide,
                        pending_ebin,
                        pending_emerge,
                        pending_esub,
                        pending_ecoll,
                        pending_emdot,
                        pending_ecdot,
                        pending_nmerge,
                    ),
                )
            end
            # Parse the new ADJUST line; an unparseable one is skipped with
            # its RSCALE and TIME[NB] lines, never kept as a zero record
            vals = _parse_adjust_fields(stripped)
            if vals === nothing
                n_unparsed += 1
                have_adjust = false
                continue
            end
            pending_time_nb = vals.time_nb
            pending_time_myr = vals.time_myr
            pending_qvir = vals.qvir
            pending_de_rel = vals.de_rel
            pending_e_tot = vals.e_tot
            pending_n = vals.n
            pending_npairs = vals.npairs
            pending_nmerge = 0
            pending_rscale = vals.rscale
            pending_de_abs = vals.de_abs
            pending_detot = vals.detot
            pending_ekin = vals.ekin
            pending_pot = vals.pot
            pending_etide = vals.etide
            pending_ebin = vals.ebin
            pending_emerge = vals.emerge
            pending_esub = vals.esub
            pending_ecoll = vals.ecoll
            pending_emdot = vals.emdot
            pending_ecdot = vals.ecdot
            have_adjust = true
            continue
        end

        # ── RMIN / RSCALE line: supplements current epoch ──
        if have_adjust && occursin("RSCALE", stripped) && occursin("RMIN", stripped)
            m = match(r"RSCALE\s*=\s*([\d.Ee\+\-]+)", stripped)
            if !isnothing(m)
                val = tryparse(Float64, something(m.captures[1]))
                !isnothing(val) && pending_rscale == 0.0 && (pending_rscale = val)
            end
            continue
        end

        # ── TIME[NB] line: N, NPAIRS for current epoch ──
        if have_adjust && startswith(stripped, "TIME[NB]")
            info = _parse_time_nb_line(stripped)
            pending_n = info.n
            pending_npairs = info.npairs
            pending_nmerge = info.nmerge
            continue
        end

        # ── PHYSICAL SCALING block ──
        if occursin("PHYSICAL SCALING", stripped) || occursin("R* =", stripped)
            _parse_scaling_line!(scaling, stripped)
            continue
        end
    end

    # Flush the last pending epoch
    if have_adjust
        push!(
            adjust_records,
            AdjustRecord(
                pending_time_nb,
                pending_time_myr,
                pending_qvir,
                pending_de_rel,
                pending_e_tot,
                pending_n,
                pending_npairs,
                pending_rscale,
                pending_de_abs,
                pending_detot,
                pending_ekin,
                pending_pot,
                pending_etide,
                pending_ebin,
                pending_emerge,
                pending_esub,
                pending_ecoll,
                pending_emdot,
                pending_ecdot,
                pending_nmerge,
            ),
        )
    end

    # Forward-fill N and NPAIRS: TIME[NB] lines appear less frequently than
    # ADJUST lines, so many records have n=0.  Carry forward the last known
    # particle count so plots don't show sawtooth artefacts.
    _forward_fill_particle_counts!(adjust_records)
    n_unparsed > 0 &&
        @warn "$path: $n_unparsed ADJUST lines skipped as unparseable (too few fields or a " *
              "Fortran field overflow)"

    return DiagnosticsData(adjust_records, scaling, kz, tidal, nbin0, nzero)
end

# ---------------------------------------------------------------------------
# Internal parsers
# ---------------------------------------------------------------------------

"""Parsed fields from an ADJUST line (before merging with supplementary lines)."""
struct _AdjustFields
    time_nb::Float64
    time_myr::Float64
    qvir::Float64
    de_rel::Float64
    e_tot::Float64
    n::Int
    npairs::Int
    rscale::Float64
    de_abs::Float64
    detot::Float64
    ekin::Float64
    pot::Float64
    etide::Float64
    ebin::Float64
    emerge::Float64
    esub::Float64
    ecoll::Float64
    emdot::Float64
    ecdot::Float64
end

function _AdjustFields(
    time_nb::Real,
    time_myr::Real,
    qvir::Real,
    de_rel::Real,
    e_tot::Real,
    n::Integer,
    npairs::Integer,
    rscale::Real,
    de_abs::Real,
    detot::Real,
)
    return _AdjustFields(
        time_nb,
        time_myr,
        qvir,
        de_rel,
        e_tot,
        n,
        npairs,
        rscale,
        de_abs,
        detot,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
    )
end

"""
Parse a single ADJUST line.  Supports two formats:
  - Positional:  T_NB  T_Myr  QV  DE  E_TOT  N  NPAIRS  RSCALE …
  - Key-value:   TIME val  T[Myr] val  Q val  DE val … ETOT val …
Reals are read as Fortran reals (`E` or `D` exponent). Returns `nothing` for
a line that does not parse (too few tokens, an overflowed `****` field), which
the caller skips.
"""
function _parse_adjust_fields(line::AbstractString)::Union{Nothing,_AdjustFields}
    body = replace(line, r"^ADJUST:\s*" => "")
    tokens = split(body)
    length(tokens) < 4 && return nothing

    try
        if uppercase(tokens[1]) == "TIME"
            return _parse_adjust_kv(tokens)
        else
            return _parse_adjust_positional(tokens)
        end
    catch e
        e isa Union{ArgumentError,ErrorException,InexactError} || rethrow()
        @debug "Skipping unparseable ADJUST line" line exception = e
        return nothing
    end
end

"""Parse positional ADJUST: T_NB T_Myr QV DE E_TOT N NPAIRS RSCALE …"""
function _parse_adjust_positional(tokens)::_AdjustFields
    length(tokens) >= 8 || error("Too few tokens")
    _AdjustFields(
        _parse_fortran_float(tokens[1]),
        _parse_fortran_float(tokens[2]),
        _parse_fortran_float(tokens[3]),
        _parse_fortran_float(tokens[4]),
        _parse_fortran_float(tokens[5]),
        parse(Int, tokens[6]),
        parse(Int, tokens[7]),
        _parse_fortran_float(tokens[8]),
        NaN,
        NaN,
    )
end

"""Parse key-value ADJUST: TIME val T[Myr] val Q val DE val … ETOT val …"""
function _parse_adjust_kv(tokens)::_AdjustFields
    kv = Dict{String,String}()
    i = 1
    while i < length(tokens)
        key = uppercase(tokens[i])
        if i + 1 <= length(tokens)
            kv[key] = tokens[i + 1]
            i += 2
        else
            break
        end
    end

    kv_real(key, default) = _parse_fortran_float(get(kv, key, default))
    time_nb = kv_real("TIME", "0")
    time_myr = kv_real("T[MYR]", "0")
    qvir = kv_real("Q", "0")
    de_rel = kv_real("DE", "0")
    e_tot = _parse_fortran_float(get(kv, "ETOT", get(kv, "E", "0")))
    n = round(Int, kv_real("N", "0"))
    npairs = round(Int, kv_real("NPAIRS", "0"))
    rscale = _parse_fortran_float(get(kv, "RSCALE", get(kv, "RSCL", "0")))
    de_abs = kv_real("DELTA", "NaN")
    detot = kv_real("DETOT", "NaN")
    ekin = kv_real("EKIN", "NaN")
    pot = kv_real("POT", "NaN")
    etide = kv_real("ETIDE", "NaN")
    ebin = kv_real("EBIN", "NaN")
    emerge = kv_real("EMERGE", "NaN")
    esub = kv_real("ESUB", "NaN")
    ecoll = kv_real("ECOLL", "NaN")
    emdot = kv_real("EMDOT", "NaN")
    ecdot = kv_real("ECDOT", "NaN")

    return _AdjustFields(
        time_nb,
        time_myr,
        qvir,
        de_rel,
        e_tot,
        n,
        npairs,
        rscale,
        de_abs,
        detot,
        ekin,
        pot,
        etide,
        ebin,
        emerge,
        esub,
        ecoll,
        emdot,
        ecdot,
    )
end

"""
Fill `kz` from the value row of the `OPTIONS` echo when it holds exactly 50
integers; leave it empty otherwise.
"""
function _parse_options_row!(kz::Vector{Int}, line::AbstractString)
    tokens = split(line)
    parsed = Int[]
    ok = length(tokens) == 50
    if ok
        for t in tokens
            v = tryparse(Int, t)
            if isnothing(v)
                ok = false
                break
            end
            push!(parsed, v)
        end
    end
    if ok
        append!(kz, parsed)
    else
        @debug "Unparseable OPTIONS echo" line
    end
    return kz
end

"""
Parse `TIDAL(1:4)` from the `TIDAL PARAMETERS:` echo; zeros when the line
does not match or a value does not parse.
"""
function _parse_tidal_line(line::AbstractString)::NTuple{4,Float64}
    m = match(r"TIDAL PARAMETERS:\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line)
    if isnothing(m)
        @debug "Unparseable TIDAL PARAMETERS echo" line
        return (0.0, 0.0, 0.0, 0.0)
    end
    try
        return (
            _parse_fortran_float(something(m.captures[1])),
            _parse_fortran_float(something(m.captures[2])),
            _parse_fortran_float(something(m.captures[3])),
            _parse_fortran_float(something(m.captures[4])),
        )
    catch e
        @debug "Unparseable TIDAL PARAMETERS echo" line exception = e
        return (0.0, 0.0, 0.0, 0.0)
    end
end

"""
Parse a TIME[NB] line for N, NPAIRS and NMERGE.

Format: `TIME[NB]  val  N  val  <NB>  val  NPAIRS  val  …`
"""
function _parse_time_nb_line(line::AbstractString)
    n = 0
    npairs = 0
    nmerge = 0

    m_n = match(r"\bN\s+(\d+)", line)
    if !isnothing(m_n)
        n = parse(Int, something(m_n.captures[1]))
    end

    m_np = match(r"NPAIRS\s+(\d+)", line)
    if !isnothing(m_np)
        npairs = parse(Int, something(m_np.captures[1]))
    end

    m_nm = match(r"NMERGE\s+(\d+)", line)
    if !isnothing(m_nm)
        nmerge = parse(Int, something(m_nm.captures[1]))
    end

    return (; n, npairs, nmerge)
end

"""
Extract key=value pairs from PHYSICAL SCALING output lines.
Patterns:  `R* = 1.234`  `T* = 5.678`  etc.
"""
function _parse_scaling_line!(scaling::Dict{String,Float64}, line::AbstractString)
    for m in eachmatch(r"([A-Z<>\*]+)\s*=\s*([\d.Ee\+\-]+)", line)
        # Both capture groups are non-optional — `something` narrows the
        # Union{Nothing,SubString} eltype for type stability (JET-clean).
        key = String(something(m.captures[1]))
        val = tryparse(Float64, String(something(m.captures[2])))
        !isnothing(val) && (scaling[key] = val)
    end
end

"""
Forward-fill N, NPAIRS and NMERGE across ADJUST records.  TIME[NB] lines (which
carry particle counts) appear only at snapshot intervals, not at every
ADJUST step.  Records without a TIME[NB] line have n=0; we replace those
with the most recent known value.
"""
function _forward_fill_particle_counts!(records::Vector{AdjustRecord})
    last_n = 0
    last_np = 0
    last_nm = 0
    for i in eachindex(records)
        r = records[i]
        if r.n > 0
            last_n = r.n
            last_np = r.npairs
            last_nm = r.nmerge
        elseif last_n > 0
            records[i] = AdjustRecord(
                r.time_nb,
                r.time_myr,
                r.qvir,
                r.de_rel,
                r.e_tot,
                last_n,
                last_np,
                r.rscale,
                r.de_abs,
                r.detot,
                r.ekin,
                r.pot,
                r.etide,
                r.ebin,
                r.emerge,
                r.esub,
                r.ecoll,
                r.emdot,
                r.ecdot,
                last_nm,
            )
        end
    end
end

"""
    extract_scaling(diag::DiagnosticsData) -> UnitScaling

Build a `UnitScaling` from parsed physical scaling data. A missing `M*` or
`T*` gives `NaN` for that factor, with one warning, so that every quantity
converted with it is `NaN` rather than a number in the wrong units; a missing
`R*` or `V*` falls back to 1.

The mass scale is `M*` (ZMBAR, the NB→M☉ conversion factor for the *total*
mass) — NOT `<M>`, which is the mean stellar mass. Nbody6++ prints both in
the PHYSICAL SCALING line; `start.F` redefines ZMBAR as the mass scale
factor at startup, and `units.f` converts masses as `ZMBAR*M`.
"""
function extract_scaling(diag::DiagnosticsData)::UnitScaling
    s = diag.physical_scaling
    missing_keys = filter(k -> !haskey(s, k), ["M*", "T*"])
    isempty(missing_keys) ||
        @warn "PHYSICAL SCALING without $(join(missing_keys, " and ")): the converted " *
              "quantities (masses in M☉, times in Myr) are NaN" maxlog = 1
    UnitScaling(get(s, "R*", 1.0), get(s, "M*", NaN), get(s, "T*", NaN), get(s, "V*", 1.0))
end

"""
    cumulative_energy_error(diag::DiagnosticsData)
        -> (; time_nb::Vector{Float64}, errtot::Vector{Float64}, detot::Vector{Float64})

Cumulative energy-conservation record of a run, one entry per adjustment.
`errtot` is the running sum of the relative interval errors `de_rel` (the
engine's `ERRTOT`, dimensionless); `detot` is the running sum of the
interval energy changes in N-body energy units (the engine's `DETOT`),
taken from the ADJUST records when every record carries it, otherwise
accumulated from `de_abs`, and `NaN` throughout when neither is available.
The ADJUST line prints four significant digits, so `errtot` agrees with the
END RUN total of [`read_energy_totals`](@ref) only to that precision.
"""
function cumulative_energy_error(diag::DiagnosticsData)
    isempty(diag.adjust) && return (; time_nb = Float64[], errtot = Float64[], detot = Float64[])

    n = length(diag.adjust)
    time_nb = [r.time_nb for r in diag.adjust]
    de_rel_vals = [r.de_rel for r in diag.adjust]
    errtot = cumsum(de_rel_vals)

    detot_vals = [r.detot for r in diag.adjust]
    de_abs_vals = [r.de_abs for r in diag.adjust]

    detot = if all(isfinite, detot_vals)
        detot_vals
    elseif all(isfinite, de_abs_vals)
        cumsum(de_abs_vals)
    else
        fill(NaN, n)
    end

    return (; time_nb, errtot, detot)
end

"""Parse a Fortran real with an `E` or `D` exponent."""
_parse_fortran_float(s::AbstractString)::Float64 = parse(Float64, replace(s, r"[dD]" => "E"))

"""
    read_energy_totals(path::AbstractString)
        -> Union{Nothing,@NamedTuple{time_nb::Float64,errtot::Float64,detot::Float64}}

Totals the engine prints on its `END RUN` line: the final time in N-body
units, `ERRTOT` (sum of the relative interval errors) and `DETOT` (sum of
the interval energy changes, N-body energy units). The last `END RUN` line
of the stdout capture is used, so a restarted run reports its final
segment; `nothing` when the capture has no such line.
"""
function read_energy_totals(path::AbstractString)
    isfile(path) || error("Diagnostics file not found: $path")

    re = r"TOFF/TIME/TTOT=\s*(\S+)\s+(\S+)\s+(\S+).*?ERRTOT\s*=\s*(\S+)\s+DETOT\s*=\s*(\S+)"
    last_totals = nothing

    for line in eachline(path)
        occursin("END RUN", line) || continue
        m = match(re, line)
        isnothing(m) && continue
        ttot = _parse_fortran_float(something(m.captures[3]))
        errtot = _parse_fortran_float(something(m.captures[4]))
        detot = _parse_fortran_float(something(m.captures[5]))
        last_totals = (; time_nb = ttot, errtot = errtot, detot = detot)
    end

    return last_totals
end
