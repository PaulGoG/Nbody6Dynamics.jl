# =============================================================================
# Integrals of motion of snapshots
# =============================================================================
# Mass, centre of mass, linear and angular momentum, kinetic and potential
# energy of every snapshot, of the whole system and of its bound set, in the
# engine frame and N-body units (G = 1). The kinetic and potential energies of
# the reduced system (regularised pairs as centres of mass) are the engine's
# EKIN and −POT of the ADJUST line, the pairs' internal energy its EBIN, so a
# snapshot taken at an adjustment time checks the snapshot data against the
# engine's double-precision bookkeeping.

"""
    SnapshotIntegrals

Integrals of one snapshot in N-body units (`G = 1`), engine frame, from
[`snapshot_integrals`](@ref):

- `time_nb`, `time_myr`; `n_stars`, `n_pairs` (regularised), `n_hierarchies`
  (pairs whose centre-of-mass name is negative, i.e. hierarchical mergers; their
  components enter as single stars)
- `mass`; `com`, `vcom`: centre of mass and its velocity
- `momentum`: ``Σ m v``
- `angular_momentum`: ``Σ m x × v`` of the stars about the origin;
  `angular_momentum_com` about the centre of mass in its frame
- `kinetic_stars`: ``½ Σ m v²`` of the stars
- `kinetic`, `potential`: ``½ Σ M V²`` and ``W = ½ Σ M Φ`` of the reduced
  system (compare `EKIN` and `−POT` of the ADJUST line)
- `pair_internal`: ``Σ (½ μ v_\\mathrm{rel}² − m_1 m_2 / r)`` of the regularised
  pairs (compare `EBIN`)
- `bound`: the [`BoundSet`](@ref) of the whole snapshot;
  `bound_mass`, `bound_angular_momentum` (its stars, about the origin)
"""
struct SnapshotIntegrals
    time_nb::Float64
    time_myr::Float64
    n_stars::Int
    n_pairs::Int
    n_hierarchies::Int
    mass::Float64
    com::NTuple{3,Float64}
    vcom::NTuple{3,Float64}
    momentum::NTuple{3,Float64}
    angular_momentum::NTuple{3,Float64}
    angular_momentum_com::NTuple{3,Float64}
    kinetic_stars::Float64
    kinetic::Float64
    potential::Float64
    pair_internal::Float64
    bound::BoundSet
    bound_mass::Float64
    bound_angular_momentum::NTuple{3,Float64}
end

"""
    snapshot_integrals(snap::Snapshot; max_passes = 100, min_members = 10,
                       backend = KernelAbstractions.CPU()) -> SnapshotIntegrals
    snapshot_integrals(snaps::Vector{Snapshot}; kwargs...) -> Vector{SnapshotIntegrals}

Integrals of motion of a snapshot (see [`SnapshotIntegrals`](@ref)). One exact
pair sum over the reduced system gives `potential` and seeds the bound-set
iteration ([`bound_set`](@ref), whose keywords are passed on). O(n²) in the
number of systems.
"""
function snapshot_integrals(
    snap::Snapshot;
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
    backend = CPU(),
)::SnapshotIntegrals
    n = nparticles(snap)
    m = Float64.(snap.mass)
    x = Float64.(snap.pos)
    v = Float64.(snap.vel)
    M = sum(m)
    com = ntuple(k -> sum(m[i] * x[k, i] for i in 1:n; init = 0.0) / M, 3)
    vcom = ntuple(k -> sum(m[i] * v[k, i] for i in 1:n; init = 0.0) / M, 3)
    P = ntuple(k -> sum(m[i] * v[k, i] for i in 1:n; init = 0.0), 3)
    L = _angular_momentum(m, x, v, 1:n, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    Lcm = _angular_momentum(m, x, v, 1:n, com, vcom)
    Tstars = 0.5 * sum(m[i] * (v[1, i]^2 + v[2, i]^2 + v[3, i]^2) for i in 1:n; init = 0.0)

    sys, _ = _snapshot_systems(snap, collect(1:n))
    hierarchies = count(<(0), snap.pair_name)
    Φ = pair_potential(sys.pos, sys.mass; backend = backend)
    T = 0.0
    W = 0.0
    @inbounds for s in eachindex(sys.mass)
        T += 0.5 * sys.mass[s] * (sys.vel[1, s]^2 + sys.vel[2, s]^2 + sys.vel[3, s]^2)
        W += 0.5 * sys.mass[s] * Φ[s]
    end
    Ebin = 0.0
    for s in eachindex(sys.mass)
        j = sys.star2[s]
        j == 0 && continue
        i = sys.star1[s]
        μ = m[i] * m[j] / (m[i] + m[j])
        r = sqrt((x[1, i] - x[1, j])^2 + (x[2, i] - x[2, j])^2 + (x[3, i] - x[3, j])^2)
        w2 = (v[1, i] - v[1, j])^2 + (v[2, i] - v[2, j])^2 + (v[3, i] - v[3, j])^2
        Ebin += 0.5 * μ * w2 - m[i] * m[j] / r
    end

    b = bound_set(
        sys;
        potential = Φ,
        max_passes = max_passes,
        min_members = min_members,
        backend = backend,
    )
    Mb = sum(m[i] for i in b.stars; init = 0.0)
    Lb = _angular_momentum(m, x, v, b.stars, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    return SnapshotIntegrals(
        time_nb(snap.header),
        time_myr(snap.header),
        n,
        snap.npairs,
        hierarchies,
        M,
        com,
        vcom,
        P,
        L,
        Lcm,
        Tstars,
        T,
        W,
        Ebin,
        b,
        Mb,
        Lb,
    )
end

function snapshot_integrals(snaps::AbstractVector{Snapshot}; kwargs...)::Vector{SnapshotIntegrals}
    return [snapshot_integrals(s; kwargs...) for s in snaps]
end

"""``Σ m (x − x₀) × (v − v₀)`` over the stars `idx`."""
function _angular_momentum(
    m,
    x,
    v,
    idx,
    x0::NTuple{3,Float64},
    v0::NTuple{3,Float64},
)::NTuple{3,Float64}
    lx = ly = lz = 0.0
    @inbounds for i in idx
        rx, ry, rz = x[1, i] - x0[1], x[2, i] - x0[2], x[3, i] - x0[3]
        ux, uy, uz = v[1, i] - v0[1], v[2, i] - v0[2], v[3, i] - v0[3]
        lx += m[i] * (ry * uz - rz * uy)
        ly += m[i] * (rz * ux - rx * uz)
        lz += m[i] * (rx * uy - ry * ux)
    end
    return (lx, ly, lz)
end

"""
    write_snapshot_integrals(path, ints::Vector{SnapshotIntegrals}) -> String

Write the integrals as CSV, one row per snapshot (N-body units plus
`time_myr`). `l*_removed` is the angular momentum of the stars lost since the
first row, `L(t₀) − L(t)` (escapers and stellar mass loss together). An
existing file is backed up, never overwritten.
"""
function write_snapshot_integrals(
    path::AbstractString,
    ints::AbstractVector{SnapshotIntegrals},
)::String
    _backup_existing(path)
    L0 = isempty(ints) ? (0.0, 0.0, 0.0) : ints[1].angular_momentum
    open(path, "w") do io
        println(
            io,
            "time_nb,time_myr,n_stars,n_pairs,n_hierarchies,mass,x_cm,y_cm,z_cm,vx_cm,vy_cm,vz_cm," *
            "px,py,pz,lx,ly,lz,lx_cm,ly_cm,lz_cm,t_stars,t_systems,w_systems,e_pairs,e_total,q_virial," *
            "n_bound,mass_bound,t_bound,w_bound,lx_bound,ly_bound,lz_bound,bound_passes,bound_status," *
            "lx_removed,ly_removed,lz_removed",
        )
        for s in ints
            b = s.bound
            cells = Any[
                s.time_nb,
                s.time_myr,
                s.n_stars,
                s.n_pairs,
                s.n_hierarchies,
                s.mass,
                s.com...,
                s.vcom...,
                s.momentum...,
                s.angular_momentum...,
                s.angular_momentum_com...,
                s.kinetic_stars,
                s.kinetic,
                s.potential,
                s.pair_internal,
                s.kinetic + s.potential + s.pair_internal,
                s.potential < 0 ? s.kinetic / abs(s.potential) : NaN,
                length(b.stars),
                s.bound_mass,
                b.kinetic,
                b.potential,
                s.bound_angular_momentum...,
                b.passes,
                String(b.status),
                (L0 .- s.angular_momentum)...,
            ]
            println(io, join(_csv_precise.(cells), ","))
        end
    end
    return String(path)
end

"""CSV cell with 12 significant digits: conservation checks need more than the 6 of `_csv_cell`."""
_csv_precise(v::AbstractFloat) = isfinite(v) ? @sprintf("%.12g", v) : ""
_csv_precise(v) = _csv_cell(v)
