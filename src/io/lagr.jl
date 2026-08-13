# =============================================================================
# Reader for Nbody6++ lagr.7 (Lagrangian radii, formatted ASCII)
# =============================================================================
#
# File format (upstream lagr.f): header lines starting with ## or TIME,
# then one row per epoch with TIME + all sub-block columns concatenated.
# The first 18 values after TIME are the total Lagrangian radii; the 19th
# is the core radius <RC (skipped here).
#
# Standard mass fractions (18 bins):
#   0.001, 0.003, 0.005, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3,
#   0.4,   0.5,   0.6,   0.7,  0.8,  0.9,  0.95, 0.99, 1.0

const LAGR_MASS_FRACTIONS = [
    0.001,
    0.003,
    0.005,
    0.01,
    0.03,
    0.05,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    0.95,
    0.99,
    1.0,
]

"""
    read_lagr(path::AbstractString) -> LagrangianData

Parse a lagr.7 file (upstream `lagr.f` single-row format: `##`/`TIME`
header lines, then one row per epoch) and return the Lagrangian radii
evolution. Raises an error on files without the expected header, rather
than mis-parsing an unrecognised layout.
"""
function read_lagr(path::AbstractString)::LagrangianData
    isfile(path) || error("Lagrangian radii file not found: $path")

    lines = readlines(path)
    isempty(lines) &&
        return LagrangianData(Float64[], LAGR_MASS_FRACTIONS, Matrix{Float64}(undef, 0, 0))

    # Fail fast on unrecognised layouts: the first non-empty line must be a
    # header (## comment or TIME column-label row, per upstream lagr.f).
    first_content = findfirst(l -> !isempty(strip(l)), lines)
    if first_content !== nothing
        s = strip(lines[first_content])
        startswith(s, '#') ||
            startswith(uppercase(s), "TIME") ||
            error(
                "Unrecognised lagr.7 format in $path: expected a ##/TIME header " *
                "line (upstream lagr.f layout), found: $(first(s, 40))",
            )
    end

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

    isempty(times) &&
        return LagrangianData(Float64[], LAGR_MASS_FRACTIONS, Matrix{Float64}(undef, 0, 0))

    nf = length(radii_rows[1])
    nt = length(times)
    mat = Matrix{Float64}(undef, nf, nt)
    for j in 1:nt
        nr = length(radii_rows[j])
        mat[1:min(nr, nf), j] .= radii_rows[j][1:min(nr, nf)]
        nr < nf && (mat[(nr + 1):nf, j] .= NaN)
    end

    fracs =
        nf == length(LAGR_MASS_FRACTIONS) ? LAGR_MASS_FRACTIONS : collect(range(0, 1; length = nf))

    return LagrangianData(times, fracs, mat)
end
