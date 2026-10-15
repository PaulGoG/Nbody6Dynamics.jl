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
- `tidal_potential`: ``−½ Σ m (T₁ x² + T₃ z²)`` over the stars, the potential
  energy of the linearised external field of `KZ(14) ≤ 2` (`xtrnlv.F`), 0
  without a field
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
    tidal_potential::Float64
    bound::BoundSet
    bound_mass::Float64
    bound_angular_momentum::NTuple{3,Float64}
end

"""
    snapshot_integrals(snap::Snapshot; max_passes = 100, min_members = 10,
                       backend = KernelAbstractions.CPU(),
                       tidal = (0.0, 0.0, 0.0, 0.0)) -> SnapshotIntegrals
    snapshot_integrals(snaps::Vector{Snapshot}; kwargs...) -> Vector{SnapshotIntegrals}

Integrals of motion of a snapshot (see [`SnapshotIntegrals`](@ref)). One exact
pair sum over the reduced system gives `potential` and seeds the bound-set
iteration ([`bound_set`](@ref), whose keywords are passed on). O(n²) in the
number of systems. `tidal` is the linearised field `TIDAL(1:4)` of
[`DiagnosticsData`](@ref); its `T₁` and `T₃` give `tidal_potential` (zeros:
no external field).
"""
function snapshot_integrals(
    snap::Snapshot;
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
    backend = CPU(),
    tidal::NTuple{4,Float64} = (0.0, 0.0, 0.0, 0.0),
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
    # Potential of the linearised field of KZ(14) ≤ 2; T₂ and T₄ do not enter it.
    Wt = 0.0
    @inbounds for i in 1:n
        Wt += m[i] * (tidal[1] * x[1, i]^2 + tidal[3] * x[3, i]^2)
    end
    Wt *= -0.5

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
        Wt,
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
first row, `L(t₀) − L(t)` (escapers and stellar mass loss together).
`e_jacobi` is the Jacobi integral `e_total + w_tidal`, equal to `e_total`
without an external field. An existing file is backed up, never overwritten.
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
            "px,py,pz,lx,ly,lz,lx_cm,ly_cm,lz_cm,t_stars,t_systems,w_systems,e_pairs,w_tidal,e_total," *
            "e_jacobi,q_virial," *
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
                s.tidal_potential,
                s.kinetic + s.potential + s.pair_internal,
                s.kinetic + s.potential + s.pair_internal + s.tidal_potential,
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

"""
    tidal_mode(kz14::Integer) -> Symbol
    tidal_mode(diag::DiagnosticsData) -> Symbol

The class of the engine's external field option `KZ(14)`: `:none` (0),
`:linearised` (1, 2: a static field in the rotating frame, whose potential
energy the engine never adds to its total — the ADJUST `DE` drifts and the
Jacobi integral `e_jacobi` of [`snapshot_integrals`](@ref) is the conserved
quantity), `:orbit` (3, 5: the engine accumulates the work of the field in
`ETIDE`, so its `DE` is an error measure, while snapshots cannot reconstruct
it), `:external` (4: neither). Any other value throws an `ArgumentError`. For a
`DiagnosticsData` without the options echo (`kz` empty) the mode is `:none`.
"""
function tidal_mode(kz14::Integer)::Symbol
    kz14 == 0 && return :none
    (kz14 == 1 || kz14 == 2) && return :linearised
    (kz14 == 3 || kz14 == 5) && return :orbit
    kz14 == 4 && return :external
    throw(ArgumentError("KZ(14) must be one of 0, 1, 2, 3, 4, 5, got $kz14"))
end

function tidal_mode(diag::DiagnosticsData)::Symbol
    isempty(diag.kz) && return :none
    length(diag.kz) >= 14 || throw(
        ArgumentError("options echo holds $(length(diag.kz)) KZ values, KZ(14) needs at least 14"),
    )
    return tidal_mode(diag.kz[14])
end

"""
    engine_energy_valid(mode::Symbol; tidal_energy_amended = false) -> Bool

Whether the engine's ADJUST `DE` measures an integration error under the
tidal mode of [`tidal_mode`](@ref): true for `:none` and `:orbit`, and for
`:linearised` when the engine was built with the `tidal-energy-total`
source amendment (`tidal_energy_amended`, read from the `source_amendments`
of the build record), which restores the tidal potential energy in the
engine's total; false for `:linearised` otherwise and for `:external`.
"""
function engine_energy_valid(mode::Symbol; tidal_energy_amended::Bool = false)::Bool
    (mode === :none || mode === :orbit) && return true
    mode === :linearised && return tidal_energy_amended
    mode === :external && return false
    throw(
        ArgumentError(
            "tidal mode must be one of :none, :linearised, :orbit, :external, got :$mode",
        ),
    )
end
