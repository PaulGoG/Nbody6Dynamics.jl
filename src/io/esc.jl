# =============================================================================
# Reader for Nbody6++ esc.11 (escaper data, formatted ASCII)
# =============================================================================
#
# esc.11 records particles that escape the cluster. The Fortran writer
# (escape.F, format `1X,1P,9E13.5,I4,I10,16E13.5`) emits per line:
#
#   token  column        unit
#   1      TTOT          NB time
#   2      BODY          NB mass
#   3      RI            NB distance at escape
#   4      VI            NB speed at escape
#   5      STEP          NB timestep
#   6      TESC = T[Myr] escape time [Myr]
#   7      BESC = M[M*]  mass at escape [M☉]
#   8      EESC          escape energy
#   9      VKM  = VI[km/s] escape velocity [km/s]
#   10     KSTARI = K*   stellar type
#   11     NAMEI         particle identifier (NAME)
#   12+    angles, SSE quantities, energy budget, tidal terms
#
# The physical-unit quantities live in tokens 6–11 — NOT 1–6 (the first five
# columns are NB-unit diagnostics). Binary escapers go to unit 31 (escbin.31),
# so every esc.11 data line has this single shape.

"""
    read_escapers(path::AbstractString) -> Vector{EscaperRecord}

Parse the esc.11 file and return a vector of `EscaperRecord` (physical
units: Myr, M☉, km/s). Skips the header, blank lines, and lines with fewer
than 11 tokens.
"""
function read_escapers(path::AbstractString)::Vector{EscaperRecord}
    isfile(path) || error("Escaper file not found: $path")

    records = EscaperRecord[]

    for line in eachline(path)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, '#') && continue

        tokens = split(stripped)
        length(tokens) < 11 && continue

        try
            tesc = parse(Float64, tokens[6])
            besc = parse(Float64, tokens[7])
            eesc = parse(Float64, tokens[8])
            vkm = parse(Float64, tokens[9])
            kstar = parse(Int, tokens[10])
            namei = parse(Int, tokens[11])
            push!(records, EscaperRecord(tesc, besc, eesc, vkm, kstar, namei))
        catch e
            @debug "Skipping unparseable escaper line" line exception = e
        end
    end

    return records
end
