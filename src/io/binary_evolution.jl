# =============================================================================
# Reader for Nbody6++ bev.82_* (regularised-binary evolution snapshots)
# =============================================================================
#
# Each bev.82_<t> file is an ASCII snapshot of the KS-regularised binaries at
# one epoch, written by `hrplot.F`. Two layouts are read:
#   Header line:  "NPAIRS  TPHYS"      — pair count and physical time [Myr]
#   Data lines (32 tokens, upstream v2026.07, FORMAT 5):
#                 TTOT  J1  J2  NAME1  NAME2  K*1  K*2  K*cm  RI[pc]  ECC
#                 LOG10(P/d)  LOG10(a/R☉)  M1[M☉]  M2[M☉]  LOG10(L1)  LOG10(L2)
#                 LOG10(R1)  LOG10(R2)  LOG10(Teff1)  LOG10(Teff2)  AGE1  AGE2
#                 EPOCH1  EPOCH2  TM1  TM2  MC1  MC2  RCC1  RCC2  RE1  RE2
#   Data lines (24 tokens, the layout of the engine manual): the same through
#                 EPOCH1 EPOCH2, with no SSE tail; TM, MC, RCC, RE of both
#                 components become NaN.
#
# As for sev.83, the header time is TPHYS in Myr while token 1 of every data
# line is TTOT in NB units. Only KS pairs (component index below IFIRST) are
# written: wide binaries beyond the regularisation distance never appear.

"""Tokens per data line of a bev.82 record (`hrplot.F` FORMAT 5)."""
const _BEV_TOKENS = 32

"""Tokens per data line of the older bev.82 layout of the engine manual (no SSE tail)."""
const _BEV_TOKENS_SHORT = 24

"""
    read_binary_evolution(path::AbstractString) -> BinaryEvolutionSnapshot

Parse a single bev.82_* file into a [`BinaryEvolutionSnapshot`](@ref). The
snapshot time is the header TPHYS [Myr]; each record additionally stores its
own per-line NB time (TTOT).

Both `hrplot.F` layouts are accepted: 32 tokens per data line (upstream
v2026.07, ending at the SSE tail TM MC RCC RE of both components) and the 24
tokens documented by the engine manual (ending at EPOCH(I1) EPOCH(I2)), for
which `ms_lifetime1_myr`, `ms_lifetime2_myr`, `mass_core1`, `mass_core2`,
`radius_core1`, `radius_core2`, `radius_envelope1` and `radius_envelope2` are
filled with `NaN`. Lines with any other token count, or with unparseable
fields (Fortran field overflow), are skipped and reported once per file as a
warning.
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

    n_data = 0
    n_skipped = 0
    for i in 2:length(lines)
        stripped = strip(lines[i])
        isempty(stripped) && continue
        n_data += 1

        tokens = split(stripped)
        if length(tokens) != _BEV_TOKENS && length(tokens) != _BEV_TOKENS_SHORT
            n_skipped += 1
            @debug "Skipping binary line of unexpected layout" line = stripped n_tokens =
                length(tokens)
            continue
        end

        try
            t_nb = parse(Float64, tokens[1])
            ints = ntuple(k -> parse(Int32, tokens[1 + k]), 7)
            # The short layout stops at EPOCH(I1) EPOCH(I2): the eight trailing
            # SSE reals (TM, MC, RCC, RE per component) are unknown.
            reals = if length(tokens) == _BEV_TOKENS
                ntuple(k -> parse(Float64, tokens[8 + k]), 24)
            else
                (ntuple(k -> parse(Float64, tokens[8 + k]), 16)..., ntuple(_ -> NaN, 8)...)
            end
            push!(records, BinaryRecord(t_nb, ints..., reals...))
        catch e
            n_skipped += 1
            @debug "Skipping unparseable binary line" line = stripped exception = e
        end
    end

    n_skipped > 0 &&
        @warn "bev file: $n_skipped of $n_data data lines skipped (unexpected layout or unparseable)" path

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
