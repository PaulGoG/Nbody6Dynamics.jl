# =============================================================================
# Reader for Nbody6++ lagr.7 (Lagrangian radii, formatted ASCII)
# =============================================================================
#
# Supports two formats:
#
# 1. Legacy block format: each epoch produces `rows_per_block` rows.
#    The first row per block is TIME followed by radii values.
#
# 2. Modern single-row format: header lines starting with ## or TIME,
#    then one row per epoch with TIME + all sub-block columns concatenated.
#    The first 18 values after TIME are the total Lagrangian radii.
#
# Standard mass fractions (18 bins):
#   0.001, 0.003, 0.005, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3,
#   0.4,   0.5,   0.6,   0.7,  0.8,  0.9,  0.95, 0.99, 1.0

const LAGR_MASS_FRACTIONS = [
    0.001, 0.003, 0.005, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3,
    0.4,   0.5,   0.6,   0.7,  0.8,  0.9,  0.95, 0.99, 1.0,
]

# Number of radii columns in the total-particle block (18 fractions + RC)
const _LAGR_NCOLS = 19

"""
    read_lagr(path::AbstractString; rows_per_block::Int = 15) -> LagrangianData

Parse the lagr.7 file and return Lagrangian radii evolution.

Auto-detects the file format (legacy block vs modern single-row).
For the legacy format, `rows_per_block` controls how many rows per epoch.
"""
function read_lagr(path::AbstractString; rows_per_block::Int = 15)::LagrangianData
    isfile(path) || error("Lagrangian radii file not found: $path")

    lines = readlines(path)
    isempty(lines) && return LagrangianData(Float64[], LAGR_MASS_FRACTIONS, Matrix{Float64}(undef, 0, 0))

    # Detect format: modern format has header lines starting with ## or non-numeric text
    if _is_modern_format(lines)
        return _read_lagr_modern(lines)
    else
        return _read_lagr_legacy(lines, rows_per_block)
    end
end

"""Detect modern format by checking if first line starts with ## or 'TIME'."""
function _is_modern_format(lines::Vector{String})::Bool
    for line in lines
        s = strip(line)
        isempty(s) && continue
        return startswith(s, '#') || startswith(uppercase(s), "TIME")
    end
    return false
end

"""
Read modern single-row format:
  ## header comment
  TIME  0.001  0.003  ...  (column labels)
  0.0   val1   val2   ...  (data rows, one per epoch)
"""
function _read_lagr_modern(lines::Vector{String})::LagrangianData
    times = Float64[]
    radii_rows = Vector{Vector{Float64}}()

    for line in lines
        s = strip(line)
        isempty(s) && continue
        startswith(s, '#') && continue

        tokens = split(s)
        isempty(tokens) && continue

        # Skip label rows (first token is not a number)
        t = tryparse(Float64, tokens[1])
        isnothing(t) && continue

        push!(times, t)

        # Extract the first 18 values after TIME as total Lagrangian radii
        # (the 19th is typically the core radius <RC, which we skip)
        nlagr = length(LAGR_MASS_FRACTIONS)
        ncols = min(nlagr, length(tokens) - 1)
        push!(radii_rows, [parse(Float64, tokens[1 + k]) for k in 1:ncols])
    end

    isempty(times) && return LagrangianData(Float64[], LAGR_MASS_FRACTIONS, Matrix{Float64}(undef, 0, 0))

    nf = length(radii_rows[1])
    nt = length(times)
    mat = Matrix{Float64}(undef, nf, nt)
    for j in 1:nt
        nr = length(radii_rows[j])
        mat[1:min(nr, nf), j] .= radii_rows[j][1:min(nr, nf)]
        nr < nf && (mat[nr+1:nf, j] .= NaN)
    end

    fracs = nf == length(LAGR_MASS_FRACTIONS) ? LAGR_MASS_FRACTIONS :
            collect(range(0, 1; length = nf))

    return LagrangianData(times, fracs, mat)
end

"""
Read legacy block format: rows_per_block rows per epoch,
first row of each block is TIME followed by radii values.
"""
function _read_lagr_legacy(lines::Vector{String}, rows_per_block::Int)::LagrangianData
    times = Float64[]
    radii_rows = Vector{Vector{Float64}}()

    if rows_per_block <= 0
        rows_per_block = _detect_block_size(lines)
    end

    for block_start in 1:rows_per_block:length(lines)
        line = strip(lines[block_start])
        isempty(line) && continue

        tokens = split(line)
        length(tokens) < 2 && continue

        t = tryparse(Float64, tokens[1])
        isnothing(t) && continue

        push!(times, t)
        push!(radii_rows, [parse(Float64, tok) for tok in tokens[2:end]])
    end

    isempty(times) && return LagrangianData(Float64[], LAGR_MASS_FRACTIONS, Matrix{Float64}(undef, 0, 0))

    nf = length(radii_rows[1])
    nt = length(times)
    mat = Matrix{Float64}(undef, nf, nt)
    for j in 1:nt
        nr = length(radii_rows[j])
        mat[1:min(nr, nf), j] .= radii_rows[j][1:min(nr, nf)]
        nr < nf && (mat[nr+1:nf, j] .= NaN)
    end

    fracs = length(LAGR_MASS_FRACTIONS) == nf ? LAGR_MASS_FRACTIONS : collect(range(0, 1; length = nf))

    return LagrangianData(times, fracs, mat)
end

"""
Auto-detect block size by finding consecutive lines starting with the same time value.
"""
function _detect_block_size(lines::Vector{String})::Int
    isempty(lines) && return 1
    first_tok = split(strip(lines[1]))[1]
    count = 1
    for i in 2:length(lines)
        tokens = split(strip(lines[i]))
        isempty(tokens) && break
        tokens[1] == first_tok ? (count += 1) : break
    end
    return max(count, 1)
end
