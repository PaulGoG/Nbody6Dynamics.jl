# =============================================================================
# Reader for Nbody6++ sev*.83 (single-star evolution snapshots)
# =============================================================================
#
# Each sev*.83 file is an ASCII snapshot of stellar properties at one epoch.
# Format:
#   Header line:   "NS  Time[NB]"  followed by  NS (number of stars) and time
#   Data lines:    I  NAME  K*  RI/RC  M[M☉]  LOG10(L)  LOG10(R)  LOG10(Teff)
#
# The first line typically has the count and NB time; subsequent lines are
# per-star records.

"""
    read_stellar_evolution(path::AbstractString) -> StellarEvolutionSnapshot

Parse a single sev*.83 file into a `StellarEvolutionSnapshot`.
"""
function read_stellar_evolution(path::AbstractString)::StellarEvolutionSnapshot
    isfile(path) || error("Stellar evolution file not found: $path")

    lines = readlines(path)
    isempty(lines) && return StellarEvolutionSnapshot(0.0, 0, StellarRecord[])

    # Parse header — first non-empty line: "NS  Time"
    header_tokens = split(strip(lines[1]))
    n_stars  = parse(Int, header_tokens[1])
    time_nb  = parse(Float64, header_tokens[2])

    records = StellarRecord[]
    sizehint!(records, n_stars)

    for i in 2:length(lines)
        stripped = strip(lines[i])
        isempty(stripped) && continue

        tokens = split(stripped)
        length(tokens) < 9 && continue

        try
            # Format: TIME_NB  INDEX  NAME  K*  RI/RC  MASS  LOGL  LOGR  LOGT  [extras...]
            _t_nb = parse(Float64, tokens[1])  # per-line time (use header time)
            idx   = parse(Int32, tokens[2])
            name  = parse(Int32, tokens[3])
            kstar = parse(Int32, tokens[4])
            ri_rc = parse(Float64, tokens[5])
            mass  = parse(Float64, tokens[6])
            logl  = parse(Float64, tokens[7])
            logr  = parse(Float64, tokens[8])
            logt  = parse(Float64, tokens[9])
            push!(records, StellarRecord(time_nb, idx, name, kstar, ri_rc, mass, logl, logr, logt))
        catch e
            @debug "Skipping unparseable stellar line" line = stripped exception = e
        end
    end

    return StellarEvolutionSnapshot(time_nb, n_stars, records)
end

"""
    read_all_stellar_evolution(dir::AbstractString, pattern::AbstractString = "sev*.83")
        -> Vector{StellarEvolutionSnapshot}

Read all sev*.83 files matching `pattern` in `dir`, sorted by time.
"""
function read_all_stellar_evolution(
    dir::AbstractString,
    pattern::AbstractString = "sev*.83",
)::Vector{StellarEvolutionSnapshot}
    isdir(dir) || return StellarEvolutionSnapshot[]

    # Glob the pattern
    files = filter(f -> _glob_match(f, pattern), readdir(dir))
    sort!(files)

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
Simple glob matching for patterns like `sev*.83` (only supports a single `*`).
"""
function _glob_match(filename::AbstractString, pattern::AbstractString)::Bool
    if !occursin('*', pattern)
        return filename == pattern
    end
    parts = split(pattern, '*'; limit = 2)
    return startswith(filename, parts[1]) && endswith(filename, parts[2])
end
