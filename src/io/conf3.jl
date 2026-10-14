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
#
# NTOT = N + NPAIRS: entries 1:N are the stars, the components of regularised
# pair k being entries 2k − 1 and 2k; entry N + k is the centre of mass of
# pair k. NPAIRS is AS(2). While a chain regularisation is active, its
# centre-of-mass body (NAME = 0) stays among the stars beside its members.
# The reader keeps the stars only.

"""
    read_conf3(path::AbstractString) -> Snapshot

Read a single conf.3 snapshot from a Fortran unformatted binary file.
Supports both the bulk-array format (one record for all data) and the
legacy per-particle record format. The centre-of-mass entries of regularised
pairs, which follow the stars in the file, are dropped (their count and names
are kept as `npairs` and `pair_name` of the returned [`Snapshot`](@ref)), and so
is the centre-of-mass body of an active chain regularisation (NAME = 0), whose
members the file carries individually. Throws `ArgumentError` when the header
pair count is inconsistent with the number of entries.
"""
function read_conf3(path::AbstractString)::Snapshot
    open(path, "r") do io
        return _read_conf3_snapshot(io)
    end
end

"""
    read_all_conf3(dir, pattern = "conf.3_*"; threaded = Threads.nthreads() > 1) -> Vector{Snapshot}

Read every conf.3 snapshot matching `pattern` in `dir`, ordered by the
numeric time suffix. With `threaded = true` (the default in a Julia session
started with more than one thread) the files are read in chunks on
`Threads.@spawn` tasks and assembled in file order, which matters once a
run holds hundreds of snapshots of 10⁵ particles; the result is identical
to the serial read. Corrupt files are skipped with a warning each.
"""
function read_all_conf3(
    dir::AbstractString,
    pattern::AbstractString = "conf.3_*";
    threaded::Bool = Threads.nthreads() > 1,
)::Vector{Snapshot}
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

    isempty(files) &&
        (@warn "No snapshot files matching '$pattern' found in $dir"; return Snapshot[])

    return _read_ordered(
        Snapshot,
        read_conf3,
        dir,
        files;
        desc = "Reading snapshots: ",
        what = "snapshot",
        threaded = threaded,
    )
end

"""
    _read_ordered(T, reader, dir, files; desc, what, threaded) -> Vector{T}

Apply `reader(path)` to every file of `files` (already in the wanted order)
and return the results in that order as a `Vector{T}`, skipping the files
whose read throws (a warning naming each, then a count). Serial, or chunked
over one `Threads.@spawn` task per thread when `threaded`; the progress bar
is shared, `next!` being thread-safe. Reading is pure, so the threaded
result equals the serial one.
"""
function _read_ordered(
    ::Type{T},
    reader,
    dir::AbstractString,
    files::AbstractVector{<:AbstractString};
    desc::AbstractString,
    what::AbstractString,
    threaded::Bool,
)::Vector{T} where {T}
    n = length(files)
    results = Vector{Union{Nothing,T}}(nothing, n)
    failures = Vector{Any}(nothing, n)
    p = Progress(n; desc = desc, showspeed = true)
    function read_one!(i)
        try
            results[i] = reader(joinpath(dir, files[i]))
        catch e
            failures[i] = (e, catch_backtrace())
        end
        next!(p)
        return nothing
    end
    if threaded && n > 1
        chunk = cld(n, min(Threads.nthreads(), n))
        tasks = [Threads.@spawn foreach(read_one!, lo:min(lo + chunk - 1, n)) for lo in 1:chunk:n]
        foreach(wait, tasks)
    else
        foreach(read_one!, 1:n)
    end
    out = T[]
    sizehint!(out, n)
    n_skipped = 0
    for i in 1:n
        r = results[i]
        if r === nothing
            n_skipped += 1
            @warn "Skipping corrupt $what: $(files[i])" exception = failures[i]
        else
            push!(out, r)
        end
    end
    n_skipped > 0 && @warn "$n_skipped / $n $what files were corrupt and skipped"
    return out
end

# ---------------------------------------------------------------------------
# Internal
# ---------------------------------------------------------------------------

function _read_conf3_snapshot(io::IO)::Snapshot
    # Record 1: header integers
    hdr_data = read_fortran_record(io)
    hdr_ints = reinterpret(Int32, hdr_data) |> collect
    ntot = hdr_ints[1]
    model = hdr_ints[2]
    nrun = hdr_ints[3]
    nk = hdr_ints[4]
    n = Int(ntot)

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
    buf = IOBuffer(data)

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

    return _star_snapshot(header, names, mass, pos, vel, rho, phi)
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
    mass = Vector{Float32}(undef, n)
    pos = Matrix{Float32}(undef, 3, n)
    vel = Matrix{Float32}(undef, 3, n)
    rho = extended ? Vector{Float32}(undef, n) : Float32[]
    phi = extended ? Vector{Float32}(undef, n) : Float32[]

    for i in 1:n
        data = read_fortran_record(io)
        buf = IOBuffer(data)
        if extended
            mass[i] = read(buf, Float32)
            rho[i] = read(buf, Float32)
            _ = read(buf, Float32)  # XNS
            pos[1, i] = read(buf, Float32)
            pos[2, i] = read(buf, Float32)
            pos[3, i] = read(buf, Float32)
            vel[1, i] = read(buf, Float32)
            vel[2, i] = read(buf, Float32)
            vel[3, i] = read(buf, Float32)
            phi[i] = read(buf, Float32)
            names[i] = read(buf, Int32)
        else
            mass[i] = read(buf, Float32)
            pos[1, i] = read(buf, Float32)
            pos[2, i] = read(buf, Float32)
            pos[3, i] = read(buf, Float32)
            vel[1, i] = read(buf, Float32)
            vel[2, i] = read(buf, Float32)
            vel[3, i] = read(buf, Float32)
            names[i] = read(buf, Int32)
        end
    end

    return _star_snapshot(header, names, mass, pos, vel, rho, phi)
end

"""
    _star_snapshot(header, name, mass, pos, vel, rho, phi) -> Snapshot

Keep the stars of the engine's NTOT entries. Dropped are the centre-of-mass
entries of regularised pairs behind the stars (`output.F` writes `J = 1, NTOT`
with NTOT = N + NPAIRS, NPAIRS in AS(2)), whose names are kept as `pair_name`,
and the centre-of-mass body of an active chain regularisation (NAME = 0,
`reduce.f`), whose members `output.F` restores from `CHDATA` beside it.
"""
function _star_snapshot(header::SnapshotHeader, name, mass, pos, vel, rho, phi)::Snapshot
    ntot = length(name)
    np_raw = header.nk ≥ 2 ? header.params[2] : 0.0f0
    (isinteger(np_raw) && np_raw ≥ 0) ||
        throw(ArgumentError("conf.3 header NPAIRS = $np_raw is not a non-negative integer"))
    npairs = Int(np_raw)
    n = ntot - npairs
    2npairs ≤ n ||
        throw(ArgumentError("conf.3 header NPAIRS = $npairs is inconsistent with NTOT = $ntot"))
    pair_name = name[(n + 1):ntot]
    chain = [i for i in (2npairs + 1):n if name[i] == 0]
    npairs == 0 &&
        isempty(chain) &&
        return Snapshot(header, name, mass, pos, vel, rho, phi, 0, pair_name)
    keep = isempty(chain) ? (1:n) : setdiff(1:n, chain)
    pick(a) = isempty(a) ? a : a[keep]
    return Snapshot(
        header,
        name[keep],
        mass[keep],
        pos[:, keep],
        vel[:, keep],
        pick(rho),
        pick(phi),
        npairs,
        pair_name,
    )
end
