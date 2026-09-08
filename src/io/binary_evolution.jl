# =============================================================================
# Reader for Nbody6++ bev.82_* (regularised-binary evolution snapshots)
# =============================================================================
#
# Each bev.82_<t> file is an ASCII snapshot of the KS-regularised binaries at
# one epoch, written by `hrplot.F` (upstream v2026.07, FORMAT 5, 32 tokens
# per line):
#   Header line:  "NPAIRS  TPHYS"      — pair count and physical time [Myr]
#   Data lines:   TTOT  J1  J2  NAME1  NAME2  K*1  K*2  K*cm  RI[pc]  ECC
#                 LOG10(P/d)  LOG10(a/R☉)  M1[M☉]  M2[M☉]  LOG10(L1)  LOG10(L2)
#                 LOG10(R1)  LOG10(R2)  LOG10(Teff1)  LOG10(Teff2)  AGE1  AGE2
#                 EPOCH1  EPOCH2  TM1  TM2  MC1  MC2  RCC1  RCC2  RE1  RE2
#
# As for sev.83, the header time is TPHYS in Myr while token 1 of every data
# line is TTOT in NB units. Only KS pairs (component index below IFIRST) are
# written: wide binaries beyond the regularisation distance never appear.

"""Tokens per data line of a bev.82 record (`hrplot.F` FORMAT 5)."""
const _BEV_TOKENS = 32

"""
    read_binary_evolution(path::AbstractString) -> BinaryEvolutionSnapshot

Parse a single bev.82_* file into a [`BinaryEvolutionSnapshot`](@ref). The
snapshot time is the header TPHYS [Myr]; each record additionally stores its
own per-line NB time (TTOT). Lines with fewer than 32 tokens, or with
unparseable fields (Fortran field overflow), are skipped.
"""
function read_binary_evolution(path::AbstractString)::BinaryEvolutionSnapshot
    isfile(path) || error("Binary evolution file not found: $path")

    lines = readlines(path)
    isempty(lines) && return BinaryEvolutionSnapshot(0.0, 0, BinaryRecord[])

    header_tokens = split(strip(lines[1]))
    length(header_tokens) ≥ 2 || error("Malformed bev.82 header in $path: \"$(strip(lines[1]))\"")
    n_pairs = parse(Int, header_tokens[1])
    time_myr = parse(Float64, header_tokens[2])

    records = BinaryRecord[]
    sizehint!(records, max(n_pairs, 0))

    for i in 2:length(lines)
        stripped = strip(lines[i])
        isempty(stripped) && continue

        tokens = split(stripped)
        length(tokens) < _BEV_TOKENS && continue

        try
            t_nb = parse(Float64, tokens[1])
            ints = ntuple(k -> parse(Int32, tokens[1 + k]), 7)
            reals = ntuple(k -> parse(Float64, tokens[8 + k]), 24)
            push!(records, BinaryRecord(t_nb, ints..., reals...))
        catch e
            @debug "Skipping unparseable binary line" line = stripped exception = e
        end
    end

    return BinaryEvolutionSnapshot(time_myr, n_pairs, records)
end

"""
    read_all_binary_evolution(dir::AbstractString, pattern::AbstractString = "bev.82_*")
        -> Vector{BinaryEvolutionSnapshot}

Read all bev.82_* files matching `pattern` in `dir`, sorted by snapshot time.
The default pattern matches what `hrplot.F` writes (`bev.82_<time>`, e.g.
`bev.82_0`, `bev.82_0.5`).
"""
function read_all_binary_evolution(
    dir::AbstractString,
    pattern::AbstractString = "bev.82_*",
)::Vector{BinaryEvolutionSnapshot}
    isdir(dir) || return BinaryEvolutionSnapshot[]

    files = filter(f -> _glob_match(f, pattern), readdir(dir))
    prefix = split(pattern, '*'; limit = 2)[1]
    sort!(files; by = f -> begin
        suffix = replace(f, prefix => ""; count = 1)
        something(tryparse(Float64, suffix), Inf)
    end)

    snapshots = BinaryEvolutionSnapshot[]
    @showprogress desc = "Reading binary evolution..." for f in files
        path = joinpath(dir, f)
        try
            push!(snapshots, read_binary_evolution(path))
        catch e
            @debug "Skipping unreadable bev file" path exception = e
        end
    end

    sort!(snapshots; by = s -> s.time_myr)
    return snapshots
end
