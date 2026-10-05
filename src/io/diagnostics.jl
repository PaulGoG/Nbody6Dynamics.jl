# =============================================================================
# Parser for Nbody6++ simulation stdout (out1000)
# =============================================================================
#
# Key patterns extracted (three related lines per epoch):
#   ADJUST:  TIME val  T[Myr] val  Q val  DE val  DELTA val  DETOT val  … ETOT val …
#   RMIN = …  RSCALE = …
#   TIME[NB]  val  N  val  NPAIRS  val  …
#
# Also:
#   PHYSICAL SCALING:  R* = …  M* = …  V* = …  T* = …

"""
    read_diagnostics(path::AbstractString) -> DiagnosticsData

Parse the simulation stdout file for ADJUST lines, supplementary RSCALE
and TIME[NB] lines (N, NPAIRS), and physical scaling info.

Each ADJUST epoch in Nbody6++ output consists of up to three lines:
1. `ADJUST:` — energy, virial ratio, time
2. `RMIN = … RSCALE = …` — scale radius (may appear between ADJUST and TIME[NB])
3. `TIME[NB] … N … NPAIRS …` — particle counts

These are merged into a single `AdjustRecord` per epoch.
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
    pending_rscale = 0.0
    pending_de_abs = NaN
    pending_detot = NaN
    have_adjust = false   # true once we've seen ADJUST: for this epoch

    for line in eachline(path)
        stripped = lstrip(line)

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
                    ),
                )
            end
            # Parse the new ADJUST line
            vals = _parse_adjust_fields(stripped)
            pending_time_nb = vals.time_nb
            pending_time_myr = vals.time_myr
            pending_qvir = vals.qvir
            pending_de_rel = vals.de_rel
            pending_e_tot = vals.e_tot
            pending_n = vals.n
            pending_npairs = vals.npairs
            pending_rscale = vals.rscale
            pending_de_abs = vals.de_abs
            pending_detot = vals.detot
            have_adjust = true
            continue
        end

        # ── RMIN / RSCALE line: supplements current epoch ──
        if have_adjust && occursin("RSCALE", stripped) && occursin("RMIN", stripped)
            m = match(r"RSCALE\s*=\s*([\d.Ee\+\-]+)", stripped)
            if !isnothing(m)
                val = tryparse(Float64, m.captures[1])
                !isnothing(val) && pending_rscale == 0.0 && (pending_rscale = val)
            end
            continue
        end

        # ── TIME[NB] line: N, NPAIRS for current epoch ──
        if have_adjust && startswith(stripped, "TIME[NB]")
            info = _parse_time_nb_line(stripped)
            pending_n = info.n
            pending_npairs = info.npairs
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
            ),
        )
    end

    # Forward-fill N and NPAIRS: TIME[NB] lines appear less frequently than
    # ADJUST lines, so many records have n=0.  Carry forward the last known
    # particle count so plots don't show sawtooth artefacts.
    _forward_fill_particle_counts!(adjust_records)

    return DiagnosticsData(adjust_records, scaling)
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
end

"""
Parse a single ADJUST line.  Supports two formats:
  - Positional:  T_NB  T_Myr  QV  DE  E_TOT  N  NPAIRS  RSCALE …
  - Key-value:   TIME val  T[Myr] val  Q val  DE val … ETOT val …
"""
function _parse_adjust_fields(line::AbstractString)::_AdjustFields
    body = replace(line, r"^ADJUST:\s*" => "")
    tokens = split(body)
    length(tokens) < 4 && return _AdjustFields(0, 0, 0, 0, 0, 0, 0, 0, NaN, NaN)

    try
        if uppercase(tokens[1]) == "TIME"
            return _parse_adjust_kv(tokens)
        else
            return _parse_adjust_positional(tokens)
        end
    catch e
        @debug "Skipping unparseable ADJUST line" line exception = e
        return _AdjustFields(0, 0, 0, 0, 0, 0, 0, 0, NaN, NaN)
    end
end

"""Parse positional ADJUST: T_NB T_Myr QV DE E_TOT N NPAIRS RSCALE …"""
function _parse_adjust_positional(tokens)::_AdjustFields
    length(tokens) >= 8 || error("Too few tokens")
    _AdjustFields(
        parse(Float64, tokens[1]),
        parse(Float64, tokens[2]),
        parse(Float64, tokens[3]),
        parse(Float64, tokens[4]),
        parse(Float64, tokens[5]),
        parse(Int, tokens[6]),
        parse(Int, tokens[7]),
        parse(Float64, tokens[8]),
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

    time_nb = parse(Float64, get(kv, "TIME", "0"))
    time_myr = parse(Float64, get(kv, "T[MYR]", "0"))
    qvir = parse(Float64, get(kv, "Q", "0"))
    de_rel = parse(Float64, get(kv, "DE", "0"))
    e_tot = parse(Float64, get(kv, "ETOT", get(kv, "E", "0")))
    n = round(Int, parse(Float64, get(kv, "N", "0")))
    npairs = round(Int, parse(Float64, get(kv, "NPAIRS", "0")))
    rscale = parse(Float64, get(kv, "RSCALE", get(kv, "RSCL", "0")))
    de_abs = parse(Float64, get(kv, "DELTA", "NaN"))
    detot = parse(Float64, get(kv, "DETOT", "NaN"))

    return _AdjustFields(time_nb, time_myr, qvir, de_rel, e_tot, n, npairs, rscale, de_abs, detot)
end

"""
Parse a TIME[NB] line for N and NPAIRS.

Format: `TIME[NB]  val  N  val  <NB>  val  NPAIRS  val  …`
"""
function _parse_time_nb_line(line::AbstractString)
    n = 0
    npairs = 0

    m_n = match(r"\bN\s+(\d+)", line)
    if !isnothing(m_n)
        n = parse(Int, m_n.captures[1])
    end

    m_np = match(r"NPAIRS\s+(\d+)", line)
    if !isnothing(m_np)
        npairs = parse(Int, m_np.captures[1])
    end

    return (; n, npairs)
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
Forward-fill N and NPAIRS across ADJUST records.  TIME[NB] lines (which
carry particle counts) appear only at snapshot intervals, not at every
ADJUST step.  Records without a TIME[NB] line have n=0; we replace those
with the most recent known value.
"""
function _forward_fill_particle_counts!(records::Vector{AdjustRecord})
    last_n = 0
    last_np = 0
    for i in eachindex(records)
        r = records[i]
        if r.n > 0
            last_n = r.n
            last_np = r.npairs
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
            )
        end
    end
end

"""
    extract_scaling(diag::DiagnosticsData) -> UnitScaling

Build a `UnitScaling` from parsed physical scaling data.
Falls back to unit values if keys are missing.

The mass scale is `M*` (ZMBAR, the NB→M☉ conversion factor for the *total*
mass) — NOT `<M>`, which is the mean stellar mass. Nbody6++ prints both in
the PHYSICAL SCALING line; `start.F` redefines ZMBAR as the mass scale
factor at startup, and `units.f` converts masses as `ZMBAR*M`.
"""
function extract_scaling(diag::DiagnosticsData)::UnitScaling
    s = diag.physical_scaling
    UnitScaling(get(s, "R*", 1.0), get(s, "M*", 1.0), get(s, "T*", 1.0), get(s, "V*", 1.0))
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
        ttot = _parse_fortran_float(m.captures[3])
        errtot = _parse_fortran_float(m.captures[4])
        detot = _parse_fortran_float(m.captures[5])
        last_totals = (; time_nb = ttot, errtot = errtot, detot = detot)
    end

    return last_totals
end
