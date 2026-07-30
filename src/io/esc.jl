# =============================================================================
# Reader for Nbody6++ esc.11 (escaper data, formatted ASCII)
# =============================================================================
#
# esc.11 records particles that escape the cluster.  Each line contains:
#   TESC  BESC  EESC  VKM  KSTARI  NAMEI
#
# Where:
#   TESC   — escape time [Myr]
#   BESC   — mass at escape [M☉]
#   EESC   — dimensionless escape energy
#   VKM    — escape velocity [km/s]
#   KSTARI — stellar type (K*)
#   NAMEI  — particle identifier (NAME)

"""
    read_escapers(path::AbstractString) -> Vector{EscaperRecord}

Parse the esc.11 file and return a vector of `EscaperRecord`.
Skips blank lines and lines with fewer than 6 tokens.
"""
function read_escapers(path::AbstractString)::Vector{EscaperRecord}
    isfile(path) || error("Escaper file not found: $path")

    records = EscaperRecord[]

    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, '#') && continue

        tokens = split(stripped)
        length(tokens) < 6 && continue

        try
            tesc   = parse(Float64, tokens[1])
            besc   = parse(Float64, tokens[2])
            eesc   = parse(Float64, tokens[3])
            vkm    = parse(Float64, tokens[4])
            kstar  = parse(Int, tokens[5])
            namei  = parse(Int, tokens[6])
            push!(records, EscaperRecord(tesc, besc, eesc, vkm, kstar, namei))
        catch e
            @debug "Skipping unparseable escaper line" line exception = e
        end
    end

    return records
end
