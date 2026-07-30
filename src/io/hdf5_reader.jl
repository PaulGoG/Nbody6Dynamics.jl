# =============================================================================
# Reader for Nbody6++ HDF5/H5Part snapshot files (data.40.h5part)
# =============================================================================
#
# H5Part layout:
#   /Step#0/  NAME[N], TIME[N], X[N], Y[N], Z[N], VX[N], VY[N], VZ[N],
#             M[N], K*[N], RHO[N], PHI[N], L[N], RS[N], Teff[N], ...
#   /Step#1/  ...

"""
    read_hdf5_snapshot(path::AbstractString, step::Int) -> Snapshot

Read a single time step from an H5Part file.
"""
function read_hdf5_snapshot(path::AbstractString, step::Int)::Snapshot
    h5open(path, "r") do fid
        return _read_h5part_step(fid, step)
    end
end

"""
    read_hdf5_snapshots(path::AbstractString) -> Vector{Snapshot}

Read all time steps from an H5Part file.
"""
function read_hdf5_snapshots(path::AbstractString)::Vector{Snapshot}
    h5open(path, "r") do fid
        step_names = sort(filter(k -> startswith(k, "Step#"), keys(fid));
                          by = k -> parse(Int, replace(k, "Step#" => "")))
        isempty(step_names) && error("No Step# groups found in $path")

        snapshots = Vector{Snapshot}(undef, length(step_names))
        p = Progress(length(step_names); desc = "Reading HDF5 steps: ", showspeed = true)
        for (i, sname) in enumerate(step_names)
            step_idx = parse(Int, replace(sname, "Step#" => ""))
            snapshots[i] = _read_h5part_step(fid, step_idx)
            next!(p)
        end
        return snapshots
    end
end

"""
    list_hdf5_steps(path::AbstractString) -> Vector{Int}

Return sorted list of step indices available in the file.
"""
function list_hdf5_steps(path::AbstractString)::Vector{Int}
    h5open(path, "r") do fid
        step_names = filter(k -> startswith(k, "Step#"), keys(fid))
        return sort([parse(Int, replace(k, "Step#" => "")) for k in step_names])
    end
end

# ---------------------------------------------------------------------------
# Internal
# ---------------------------------------------------------------------------

function _read_h5part_step(fid, step::Int)::Snapshot
    group_name = "Step#$step"
    haskey(fid, group_name) || error("Step group '$group_name' not found in HDF5 file")
    g = fid[group_name]

    gkeys = keys(g)

    # Required datasets
    names_raw = _h5read_dataset(g, gkeys, ["NAME", "Name", "name"])
    x   = _h5read_dataset(g, gkeys, ["X", "x"])
    y   = _h5read_dataset(g, gkeys, ["Y", "y"])
    z   = _h5read_dataset(g, gkeys, ["Z", "z"])
    vx  = _h5read_dataset(g, gkeys, ["VX", "Vx", "vx"])
    vy  = _h5read_dataset(g, gkeys, ["VY", "Vy", "vy"])
    vz  = _h5read_dataset(g, gkeys, ["VZ", "Vz", "vz"])
    m   = _h5read_dataset(g, gkeys, ["M", "m", "Mass", "mass"])

    n = length(x)

    # Construct position/velocity matrices
    pos = Matrix{Float32}(undef, 3, n)
    vel = Matrix{Float32}(undef, 3, n)
    pos[1, :] .= Float32.(x)
    pos[2, :] .= Float32.(y)
    pos[3, :] .= Float32.(z)
    vel[1, :] .= Float32.(vx)
    vel[2, :] .= Float32.(vy)
    vel[3, :] .= Float32.(vz)

    # Optional datasets
    rho = _h5try_dataset(g, gkeys, ["RHO", "Rho", "rho"], n)
    phi = _h5try_dataset(g, gkeys, ["PHI", "Phi", "phi"], n)

    # Extract time from TIME dataset or attribute
    t_nb = _h5read_time(g, gkeys)

    # Build a minimal header from HDF5 metadata
    params = zeros(Float32, 20)
    params[1] = Float32(t_nb)

    header = SnapshotHeader(Int32(n), Int32(step), Int32(0), Int32(20), params)

    return Snapshot(header, Int32.(names_raw), Float32.(m), pos, vel, rho, phi)
end

"""
Try multiple possible dataset names (HDF5 field naming varies across versions).
"""
function _h5read_dataset(g, gkeys, candidates::Vector{String})
    for name in candidates
        name in gkeys && return read(g[name])
    end
    error("None of $(candidates) found in HDF5 group. Available: $(gkeys)")
end

"""
Try to read an optional dataset; return empty Float32[] if not found.
"""
function _h5try_dataset(g, gkeys, candidates::Vector{String}, n::Int)::Vector{Float32}
    for name in candidates
        if name in gkeys
            return Float32.(read(g[name]))
        end
    end
    return Float32[]
end

"""
Extract simulation time from the step group.
"""
function _h5read_time(g, gkeys)::Float64
    # TIME may be a dataset (H5Part convention) or an attribute
    for name in ["TIME", "Time", "time"]
        if name in gkeys
            t = read(g[name])
            return Float64(t isa AbstractArray ? t[1] : t)
        end
    end
    # Try attributes
    attrs = HDF5.attributes(g)
    for name in ["TIME", "Time", "time"]
        if haskey(attrs, name)
            return Float64(read(attrs[name]))
        end
    end
    return 0.0
end
