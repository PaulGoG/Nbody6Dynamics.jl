# =============================================================================
# Reader for Nbody6++ sev.83_* (single-star evolution snapshots)
# =============================================================================
#
# Each sev.83_<t> file is an ASCII snapshot of stellar properties at one
# epoch, written by `hrplot.F`. Two layouts are read:
#   Header line:  "NS  TPHYS"        — star count and physical time [Myr]
#   Data lines (15 tokens, `hrplot.F` of upstream v2026.07+):
#                 TTOT  I  NAME  K*  RI[pc]  M[M☉]  LOG10(L)  LOG10(R)
#                 LOG10(Teff)  AGE  EPOCH  TM[Myr]  MC[M☉]  RCC[R☉]  RE[R☉]
#   Data lines (11 tokens, the layout of the engine manual): the same through
#                 AGE  EPOCH, with no SSE tail; TM, MC, RCC, RE become NaN.
#
# Header time is TPHYS in **Myr**, while token 1 of every data line is TTOT
# in **NB units** — both clocks are kept.

"""
    read_stellar_evolution(path::AbstractString) -> StellarEvolutionSnapshot

Parse a single sev.83_* file into a `StellarEvolutionSnapshot`. The snapshot
time is the header TPHYS [Myr]; each record additionally stores its own
per-line NB time (TTOT).

Both `hrplot.F` layouts are accepted: 15 tokens per data line (upstream
v2026.07+, ending at TM MC RCC RE) and the 11 tokens documented by the engine
manual (ending at AGE EPOCH), for which `ms_lifetime_myr`, `mass_core`,
`radius_core` and `radius_envelope` are filled with `NaN`. Lines with any
other token count, or with unparseable fields (Fortran field overflow), are
skipped and reported once per file as a warning.
"""
function read_stellar_evolution(path::AbstractString)::StellarEvolutionSnapshot
    isfile(path) || error("Stellar evolution file not found: $path")

    lines = readlines(path)
    isempty(lines) && return StellarEvolutionSnapshot(0.0, 0, StellarRecord[])

    # Parse header — first non-empty line: "NS  TPHYS[Myr]"
    header_tokens = split(strip(lines[1]))
    n_stars = parse(Int, header_tokens[1])
    time_myr = parse(Float64, header_tokens[2])

    records = StellarRecord[]
    sizehint!(records, n_stars)

    n_data = 0
    n_skipped = 0
    for i in 2:length(lines)
        stripped = strip(lines[i])
        isempty(stripped) && continue
        n_data += 1

        tokens = split(stripped)
        if length(tokens) != 15 && length(tokens) != 11
            n_skipped += 1
            @debug "Skipping stellar line of unexpected layout" line = stripped n_tokens =
                length(tokens)
            continue
        end

        try
            t_nb = parse(Float64, tokens[1])
            idx = parse(Int32, tokens[2])
            name = parse(Int32, tokens[3])
            kstar = parse(Int32, tokens[4])
            ri = parse(Float64, tokens[5])
            mass = parse(Float64, tokens[6])
            logl = parse(Float64, tokens[7])
            logr = parse(Float64, tokens[8])
            logt = parse(Float64, tokens[9])
            # The 11-token layout stops at AGE EPOCH: no SSE tail to read.
            tm, mc, rcc, re = if length(tokens) == 15
                (
                    parse(Float64, tokens[12]),
                    parse(Float64, tokens[13]),
                    parse(Float64, tokens[14]),
                    parse(Float64, tokens[15]),
                )
            else
                (NaN, NaN, NaN, NaN)
            end
            push!(
                records,
                StellarRecord(t_nb, idx, name, kstar, ri, mass, logl, logr, logt, tm, mc, rcc, re),
            )
        catch e
            n_skipped += 1
            @debug "Skipping unparseable stellar line" line = stripped exception = e
        end
    end

    n_skipped > 0 &&
        @warn "sev file: $n_skipped of $n_data data lines skipped (unexpected layout or unparseable)" path

    return StellarEvolutionSnapshot(time_myr, n_stars, records)
end

"""
    read_all_stellar_evolution(dir::AbstractString, pattern::AbstractString = "sev.83_*")
        -> Vector{StellarEvolutionSnapshot}

Read all sev.83_* files matching `pattern` in `dir`, sorted by snapshot time.
The default pattern matches what this fork's `hrplot.F` actually writes
(`sev.83_<time>`, e.g. `sev.83_0`, `sev.83_0.5`).
"""
function read_all_stellar_evolution(
    dir::AbstractString,
    pattern::AbstractString = "sev.83_*",
)::Vector{StellarEvolutionSnapshot}
    isdir(dir) || return StellarEvolutionSnapshot[]

    # Glob the pattern; order files by numeric time suffix where possible
    # (final ordering is by parsed snapshot time anyway).
    files = filter(f -> _glob_match(f, pattern), readdir(dir))
    prefix = split(pattern, '*'; limit = 2)[1]
    sort!(files; by = f -> begin
        suffix = replace(f, prefix => ""; count = 1)
        something(tryparse(Float64, suffix), Inf)
    end)

    snapshots = StellarEvolutionSnapshot[]
    @showprogress desc = "Reading stellar evolution..." for f in files
        path = joinpath(dir, f)
        try
            snap = read_stellar_evolution(path)
            push!(snapshots, snap)
        catch e
            @debug "Skipping unreadable sev file" path exception = e
        end
    end

    sort!(snapshots; by = s -> s.time_myr)
    return snapshots
end

"""
Simple glob matching for patterns like `sev.83_*` (only supports a single `*`).
"""
function _glob_match(filename::AbstractString, pattern::AbstractString)::Bool
    if !occursin('*', pattern)
        return filename == pattern
    end
    parts = split(pattern, '*'; limit = 2)
    return startswith(filename, parts[1]) && endswith(filename, parts[2])
end
