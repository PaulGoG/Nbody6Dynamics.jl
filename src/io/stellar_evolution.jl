# =============================================================================
# Reader for Nbody6++ sev.83_* (single-star evolution snapshots)
# =============================================================================
#
# Each sev.83_<t> file is an ASCII snapshot of stellar properties at one
# epoch, written by `hrplot.F` (upstream v2026.07+, 15 tokens per line):
#   Header line:  "NS  TPHYS"        — star count and physical time [Myr]
#   Data lines:   TTOT  I  NAME  K*  RI[pc]  M[M☉]  LOG10(L)  LOG10(R)
#                 LOG10(Teff)  AGE  EPOCH  TM[Myr]  MC[M☉]  RCC[R☉]  RE[R☉]
#
# Header time is TPHYS in **Myr**, while token 1 of every data line is TTOT
# in **NB units** — both clocks are kept.

"""
    read_stellar_evolution(path::AbstractString) -> StellarEvolutionSnapshot

Parse a single sev.83_* file into a `StellarEvolutionSnapshot`. The snapshot
time is the header TPHYS [Myr]; each record additionally stores its own
per-line NB time (TTOT).
"""
function read_stellar_evolution(path::AbstractString)::StellarEvolutionSnapshot
    isfile(path) || error("Stellar evolution file not found: $path")

    lines = readlines(path)
    isempty(lines) && return StellarEvolutionSnapshot(0.0, 0, StellarRecord[])

    # Parse header — first non-empty line: "NS  TPHYS[Myr]"
    header_tokens = split(strip(lines[1]))
    n_stars  = parse(Int, header_tokens[1])
    time_myr = parse(Float64, header_tokens[2])

    records = StellarRecord[]
    sizehint!(records, n_stars)

    for i in 2:length(lines)
        stripped = strip(lines[i])
        isempty(stripped) && continue

        tokens = split(stripped)
        length(tokens) < 15 && continue

        try
            t_nb  = parse(Float64, tokens[1])
            idx   = parse(Int32, tokens[2])
            name  = parse(Int32, tokens[3])
            kstar = parse(Int32, tokens[4])
            ri    = parse(Float64, tokens[5])
            mass  = parse(Float64, tokens[6])
            logl  = parse(Float64, tokens[7])
            logr  = parse(Float64, tokens[8])
            logt  = parse(Float64, tokens[9])
            tm    = parse(Float64, tokens[12])
            mc    = parse(Float64, tokens[13])
            rcc   = parse(Float64, tokens[14])
            re    = parse(Float64, tokens[15])
            push!(records, StellarRecord(t_nb, idx, name, kstar, ri, mass,
                                         logl, logr, logt, tm, mc, rcc, re))
        catch e
            @debug "Skipping unparseable stellar line" line = stripped exception = e
        end
    end

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
