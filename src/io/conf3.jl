# =============================================================================
# Reader for Nbody6++ conf.3 snapshot files (Fortran unformatted binary)
# =============================================================================
#
# Layout:
#   Record 1:  NTOT, MODEL, NRUN, NK                     (4 × Int32)
#   Record 2:  Bulk array containing params + all particle data
#
# Record 2 layout (Nbody6PPGPU-beijing output.F):
#   AS(1:NK)            — header parameters
#   BODYS(1:NTOT)       — masses
#   RHOS(1:NTOT)        — local densities
#   XNS(1:NTOT)         — neighbour numbers (unused)
#   XS(1:3, 1:NTOT)     — positions (Fortran column-major)
#   VS(1:3, 1:NTOT)     — velocities (Fortran column-major)
#   PHI(1:NTOT)         — potentials
#   NAME(1:NTOT)        — particle identifiers (Int32)

"""
    read_conf3(path::AbstractString) -> Snapshot

Read a single conf.3 snapshot from a Fortran unformatted binary file.
Supports both the bulk-array format (one record for all data) and the
legacy per-particle record format.
"""
function read_conf3(path::AbstractString)::Snapshot
    open(path, "r") do io
        return _read_conf3_snapshot(io)
    end
end

"""
    read_all_conf3(dir::AbstractString, pattern::AbstractString = "conf.3_*") -> Vector{Snapshot}

Read all conf.3 snapshot files matching `pattern` in `dir`, sorted by filename.
"""
function read_all_conf3(dir::AbstractString, pattern::AbstractString = "conf.3_*")::Vector{Snapshot}
    prefix = replace(pattern, "*" => "")

    files = filter(readdir(dir; join = false)) do f
        startswith(f, prefix) || f == replace(pattern, "_*" => "")
    end
    # Sort numerically by the suffix after the prefix. Suffixes may be integer
    # (conf.3_10) or decimal (conf.3_0.5, conf.3_10.0) depending on DELTAT; both
    # must parse so that time-order is preserved for downstream evolution plots.
    sort!(files; by = f -> begin
        suffix = replace(f, prefix => ""; count = 1)
        something(tryparse(Float64, suffix), Inf)
    end)

    isempty(files) && (@warn "No snapshot files matching '$pattern' found in $dir"; return Snapshot[])

    snapshots = Snapshot[]
    sizehint!(snapshots, length(files))
    n_skipped = 0
    p = Progress(length(files); desc = "Reading snapshots: ", showspeed = true)
    for fname in files
        try
            push!(snapshots, read_conf3(joinpath(dir, fname)))
        catch e
            n_skipped += 1
            @warn "Skipping corrupt snapshot: $fname" exception=(e, catch_backtrace())
        end
        next!(p)
    end
    n_skipped > 0 && @warn "$n_skipped / $(length(files)) snapshot files were corrupt and skipped"
    return snapshots
end

# ---------------------------------------------------------------------------
# Internal
# ---------------------------------------------------------------------------

function _read_conf3_snapshot(io::IO)::Snapshot
    # Record 1: header integers
    hdr_data = read_fortran_record(io)
    hdr_ints = reinterpret(Int32, hdr_data) |> collect
    ntot  = hdr_ints[1]
    model = hdr_ints[2]
    nrun  = hdr_ints[3]
    nk    = hdr_ints[4]
    n     = Int(ntot)

    # Record 2: peek at size to determine format
    rec2_size = peek_record_size(io)
    expected_bulk = (Int(nk) + 11 * n) * 4  # NK params + 11 arrays × NTOT

    if rec2_size == expected_bulk
        return _read_bulk_format(io, ntot, model, nrun, nk, n)
    else
        # Legacy per-particle records: record 2 = params, then N particle records
        return _read_per_particle_format(io, ntot, model, nrun, nk, n)
    end
end

"""
Read the bulk-array format where record 2 contains params + all particle data.
"""
function _read_bulk_format(io::IO, ntot, model, nrun, nk, n::Int)::Snapshot
    data = read_fortran_record(io)
    buf  = IOBuffer(data)

    # 1. AS parameters
    params = Vector{Float32}(undef, nk)
    read!(buf, params)

    header = SnapshotHeader(ntot, model, nrun, nk, params)

    # 2. BODYS — masses
    mass = Vector{Float32}(undef, n)
    read!(buf, mass)

    # 3. RHOS — densities
    rho = Vector{Float32}(undef, n)
    read!(buf, rho)

    # 4. XNS — neighbour counts (skip)
    skip(buf, n * 4)

    # 5. XS(1:3, 1:NTOT) — positions (Fortran column-major = x1,y1,z1,x2,y2,z2,...)
    pos = Matrix{Float32}(undef, 3, n)
    read!(buf, pos)

    # 6. VS(1:3, 1:NTOT) — velocities
    vel = Matrix{Float32}(undef, 3, n)
    read!(buf, vel)

    # 7. PHI — potentials
    phi = Vector{Float32}(undef, n)
    read!(buf, phi)

    # 8. NAME — particle identifiers
    names = Vector{Int32}(undef, n)
    read!(buf, names)

    return Snapshot(header, names, mass, pos, vel, rho, phi)
end

"""
Read the legacy per-particle record format where each particle has its own record.
"""
function _read_per_particle_format(io::IO, ntot, model, nrun, nk, n::Int)::Snapshot
    params = read_fortran_record(io, Float32, Int(nk))
    header = SnapshotHeader(ntot, model, nrun, nk, params)

    # Detect format from first particle record
    rec_size = peek_record_size(io)
    extended = (rec_size == 44)

    names = Vector{Int32}(undef, n)
    mass  = Vector{Float32}(undef, n)
    pos   = Matrix{Float32}(undef, 3, n)
    vel   = Matrix{Float32}(undef, 3, n)
    rho   = extended ? Vector{Float32}(undef, n) : Float32[]
    phi   = extended ? Vector{Float32}(undef, n) : Float32[]

    for i in 1:n
        data = read_fortran_record(io)
        buf = IOBuffer(data)
        if extended
            mass[i]   = read(buf, Float32)
            rho[i]    = read(buf, Float32)
            _         = read(buf, Float32)  # XNS
            pos[1, i] = read(buf, Float32)
            pos[2, i] = read(buf, Float32)
            pos[3, i] = read(buf, Float32)
            vel[1, i] = read(buf, Float32)
            vel[2, i] = read(buf, Float32)
            vel[3, i] = read(buf, Float32)
            phi[i]    = read(buf, Float32)
            names[i]  = read(buf, Int32)
        else
            mass[i]   = read(buf, Float32)
            pos[1, i] = read(buf, Float32)
            pos[2, i] = read(buf, Float32)
            pos[3, i] = read(buf, Float32)
            vel[1, i] = read(buf, Float32)
            vel[2, i] = read(buf, Float32)
            vel[3, i] = read(buf, Float32)
            names[i]  = read(buf, Int32)
        end
    end

    return Snapshot(header, names, mass, pos, vel, rho, phi)
end
