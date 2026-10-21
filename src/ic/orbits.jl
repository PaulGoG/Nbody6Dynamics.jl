# =============================================================================
# Orbital mechanics and multi-cluster assembly
# =============================================================================

# ---------------------------------------------------------------------------
# Code units of the IC generator: G = 1 with masses in M☉ and lengths in pc.
# ---------------------------------------------------------------------------

"""Gravitational constant in pc (km s⁻¹)² M☉⁻¹ (CODATA 2018 G with the IAU 2015 nominal solar mass parameter)."""
const _G_PC_KMS2_MSUN = 4.30091e-3

"""Code velocity unit √(G M☉ / pc) in km s⁻¹ (≈ 0.06558)."""
const _CODE_VSTAR_KMS = sqrt(_G_PC_KMS2_MSUN)

# Code time unit pc / (km s⁻¹) / _CODE_VSTAR_KMS in Myr (≈ 14.91 Myr), the
# crossing-time unit of a 1 M☉, 1 pc system with G = 1.
const _PC_KM = 3.0857e13
const _MYR_S = 3.15576e13
const _CODE_TIME_MYR = _PC_KM / _MYR_S / _CODE_VSTAR_KMS

"""
    kepler_velocity(M1, M2, d_apo, ecc) -> (v1, v2)

Compute centre-of-mass velocities at apocentre for a two-body Keplerian orbit.

# Arguments
- `M1`, `M2`: cluster masses (any consistent unit system, G=1 assumed)
- `d_apo`: apocentre separation
- `ecc`: orbital eccentricity (0 = circular, <1 = bound)

# Returns
Tuple `(v1, v2)` — scalar speeds along the tangential direction at apocentre.
The direction is perpendicular to the line connecting the two centres.

At apocentre, the velocity is purely tangential:
```math
v_{\\rm apo} = \\sqrt{\\frac{G(M_1+M_2)(1-e)}{a(1+e)}}
```
where ``a = d_{\\rm apo}/(1+e)`` is the semi-major axis.
"""
function kepler_velocity(M1::Float64, M2::Float64, d_apo::Float64, ecc::Float64)
    0.0 ≤ ecc < 1.0 || throw(ArgumentError("Eccentricity must be in [0, 1), got $ecc"))
    d_apo > 0.0 || throw(ArgumentError("Apocentre must be positive, got $d_apo"))

    a = d_apo / (1.0 + ecc)  # semi-major axis
    M_total = M1 + M2

    # Orbital speed at apocentre (vis-viva at r = d_apo)
    v_apo = sqrt(M_total * (1.0 - ecc) / (a * (1.0 + ecc)))

    # CM decomposition: v1 = (M2/M_tot) v_apo, v2 = (M1/M_tot) v_apo
    v1 = (M2 / M_total) * v_apo
    v2 = (M1 / M_total) * v_apo

    return v1, v2
end

"""
    jacobi_radius(d, M_self, M_other) -> Float64

Instantaneous Jacobi (tidal) radius at separation `d`:

```math
r_J = d \\left(\\frac{M_{\\rm self}}{3 M_{\\rm other}}\\right)^{1/3}
```
"""
function jacobi_radius(d::Float64, M_self::Float64, M_other::Float64)
    return d * (M_self / (3.0 * M_other))^(1.0 / 3.0)
end

"""
    truncate_jacobi(pos, vel, mass, r_trunc) -> (pos_t, vel_t, mass_t)

Return copies of `(pos, vel, mass)` with particles beyond `r_trunc` from the
coordinate origin removed, and the indices of the kept particles. Non-mutating (the inputs are left untouched), so
no `!` — the cluster is assumed centred at the origin (pre-offset).
"""
function truncate_jacobi(
    pos::Matrix{Float64},
    vel::Matrix{Float64},
    mass::Vector{Float64},
    r_trunc::Float64,
)
    N = length(mass)
    keep = Bool[]
    sizehint!(keep, N)
    for i in 1:N
        r = sqrt(pos[1, i]^2 + pos[2, i]^2 + pos[3, i]^2)
        push!(keep, r ≤ r_trunc)
    end
    idx = findall(keep)
    n_removed = N - length(idx)
    n_removed > 0 &&
        @info "Jacobi truncation: removed $n_removed / $N particles (r_trunc = $(round(r_trunc; digits=4)))"
    return pos[:, idx], vel[:, idx], mass[idx], idx
end

"""
    setup_two_cluster_orbit(pos1, vel1, mass1, pos2, vel2, mass2,
                            d_apo, ecc; truncate_jacobi_flag = true,
                            truncation_radii = (0.0, 0.0),
                            orbit_normal = [0.0, 0.0, -1.0], apocentre_angle = 0.0)
        -> (pos, vel, mass, cluster_ranges, kept)

Place two clusters on a Keplerian orbit at apocentre (code units, `G = 1`).

With `n = orbit_normal / |orbit_normal|`, the reference direction
`u ∝ x̂ − (x̂·n) n` (`ŷ` in place of `x̂` when `n ∥ x̂`) and
`φ = apocentre_angle` [degrees], the apocentre direction is
`e_apo = u cos φ + (n × u) sin φ` and the velocity direction
`e_vel = n × e_apo`. Cluster 1 is shifted by `−d1 e_apo` in position and by
`−v1 e_vel` in velocity, cluster 2 by `+d2 e_apo` and `+v2 e_vel`, where
`d1, d2` are the CM distances and `v1, v2` the apocentre speeds
([`kepler_velocity`](@ref)). The separation of the two centres at the start
is `d_apo e_apo` and the orbital angular momentum of the pair is along `n`.
The defaults give `e_apo = x̂`, `e_vel = −ŷ`: cluster 1 at `(−d1, 0, 0)`
moving along `+ŷ`, cluster 2 at `(+d2, 0, 0)` moving along `−ŷ`, bit for bit
the layout of earlier versions.

If `truncate_jacobi_flag = true`, cluster `i` is truncated before combining
at `truncation_radii[i]` when that is positive, else at its instantaneous
Jacobi radius ([`jacobi_radius`](@ref)).

Returns combined `(pos, vel, mass)` arrays in the system CM frame plus a
vector of `UnitRange{Int}` giving each cluster's (post-truncation) particle
index range in the output and the indices of the kept particles in each
input set — the same contract as [`combine_clusters_explicit`](@ref).
"""
function setup_two_cluster_orbit(
    pos1::Matrix{Float64},
    vel1::Matrix{Float64},
    mass1::Vector{Float64},
    pos2::Matrix{Float64},
    vel2::Matrix{Float64},
    mass2::Vector{Float64},
    d_apo::Float64,
    ecc::Float64;
    truncate_jacobi_flag::Bool = true,
    truncation_radii::Tuple{Float64,Float64} = (0.0, 0.0),
    orbit_normal::AbstractVector{<:Real} = [0.0, 0.0, -1.0],
    apocentre_angle::Real = 0.0,
)
    e_apo, e_vel = _orbit_directions(orbit_normal, apocentre_angle)
    M1 = sum(mass1)
    M2 = sum(mass2)
    M_total = M1 + M2

    # Optional truncation: fixed radii where given, else the Jacobi radii
    kept1 = collect(1:length(mass1))
    kept2 = collect(1:length(mass2))
    if truncate_jacobi_flag
        r_trunc_1 = truncation_radii[1] > 0 ? truncation_radii[1] : jacobi_radius(d_apo, M1, M2)
        r_trunc_2 = truncation_radii[2] > 0 ? truncation_radii[2] : jacobi_radius(d_apo, M2, M1)
        pos1, vel1, mass1, kept1 = truncate_jacobi(pos1, vel1, mass1, r_trunc_1)
        pos2, vel2, mass2, kept2 = truncate_jacobi(pos2, vel2, mass2, r_trunc_2)
        # Update totals after truncation
        M1 = sum(mass1)
        M2 = sum(mass2)
        M_total = M1 + M2
    end

    # CM distances
    d1 = (M2 / M_total) * d_apo  # cluster 1 distance from CM
    d2 = (M1 / M_total) * d_apo  # cluster 2 distance from CM

    # Apocentre velocities
    v1, v2 = kepler_velocity(M1, M2, d_apo, ecc)

    N1 = length(mass1)
    N2 = length(mass2)
    N_total = N1 + N2

    pos_out = zeros(Float64, 3, N_total)
    vel_out = zeros(Float64, 3, N_total)
    mass_out = zeros(Float64, N_total)

    # Cluster 1: shift by −d1 e_apo with velocity −v1 e_vel; with the default
    # orientation the products are by ±1 and ±0, so the result is exact.
    for i in 1:N1
        for k in 1:3
            pos_out[k, i] = pos1[k, i] - d1 * e_apo[k]
            vel_out[k, i] = vel1[k, i] - v1 * e_vel[k]
        end
        mass_out[i] = mass1[i]
    end

    # Cluster 2: shift by +d2 e_apo with velocity +v2 e_vel
    for i in 1:N2
        j = N1 + i
        for k in 1:3
            pos_out[k, j] = pos2[k, i] + d2 * e_apo[k]
            vel_out[k, j] = vel2[k, i] + v2 * e_vel[k]
        end
        mass_out[j] = mass2[i]
    end

    return pos_out, vel_out, mass_out, [1:N1, (N1 + 1):N_total], [kept1, kept2]
end

"""
    _orbit_directions(orbit_normal, apocentre_angle) -> (e_apo, e_vel)

Unit apocentre direction `e_apo = u cos φ + (n × u) sin φ` and velocity
direction `e_vel = n × e_apo` of a two-cluster orbit with unit normal
`n = orbit_normal / |orbit_normal|`, reference direction `u ∝ x̂ − (x̂·n) n`
(`ŷ − (ŷ·n) n` when that is shorter than 10⁻⁸) and `φ = apocentre_angle`
[degrees], as 3-tuples. The defaults `[0, 0, −1]`, `0` give exactly
`(1, 0, 0)` and `(0, −1, 0)`. Throws an `ArgumentError` unless the normal has
three finite components of non-zero norm and the angle is finite.
"""
function _orbit_directions(orbit_normal::AbstractVector{<:Real}, apocentre_angle::Real)
    length(orbit_normal) == 3 ||
        throw(ArgumentError("orbit_normal must have 3 components; got $(length(orbit_normal))"))
    nx, ny, nz = Float64(orbit_normal[1]), Float64(orbit_normal[2]), Float64(orbit_normal[3])
    n_norm = sqrt(nx^2 + ny^2 + nz^2)
    (isfinite(n_norm) && n_norm > 0) ||
        throw(ArgumentError("orbit_normal must be finite with non-zero norm; got $orbit_normal"))
    isfinite(apocentre_angle) ||
        throw(ArgumentError("apocentre_angle must be finite; got $apocentre_angle"))
    nx, ny, nz = nx / n_norm, ny / n_norm, nz / n_norm
    # u = x̂ − (x̂·n) n, or ŷ − (ŷ·n) n for a normal along x̂
    ux, uy, uz = 1.0 - nx * nx, 0.0 - nx * ny, 0.0 - nx * nz
    if sqrt(ux^2 + uy^2 + uz^2) < 1e-8
        ux, uy, uz = 0.0 - ny * nx, 1.0 - ny * ny, 0.0 - ny * nz
    end
    u_norm = sqrt(ux^2 + uy^2 + uz^2)
    ux, uy, uz = ux / u_norm, uy / u_norm, uz / u_norm
    # w = n × u
    wx, wy, wz = ny * uz - nz * uy, nz * ux - nx * uz, nx * uy - ny * ux
    φ = deg2rad(Float64(apocentre_angle))
    c, s = cos(φ), sin(φ)
    e_apo = (ux * c + wx * s, uy * c + wy * s, uz * c + wz * s)
    e_vel = (
        ny * e_apo[3] - nz * e_apo[2],
        nz * e_apo[1] - nx * e_apo[3],
        nx * e_apo[2] - ny * e_apo[1],
    )
    return e_apo, e_vel
end

"""
    two_cluster_orbit_quantities(M1, M2, d_apo, ecc, rh1; sigma1 = NaN, e_total = NaN,
                                 r_out = (NaN, NaN)) -> NamedTuple

Derived quantities of a two-cluster Kepler orbit started at apocentre, in pc,
km s⁻¹, M☉ and Myr with `G = _G_PC_KMS2_MSUN` [pc km² s⁻² M☉⁻¹].

# Arguments
- `M1`, `M2`: cluster masses [M☉]; > 0
- `d_apo`: apocentre separation [pc]; > 0
- `ecc`: eccentricity, `0 ≤ e < 1`
- `rh1`: half-mass radius of cluster 1 [pc]; > 0
- `sigma1`: velocity dispersion of cluster 1 [km s⁻¹] (3-D, mass-weighted rms
  speed about the mean)
- `e_total`: total energy of the configuration, internal plus orbital
  [M☉ km² s⁻²]
- `r_out`: outer radii of the two clusters [pc]

# Returns
With `M = M1 + M2`, `μ = M1 M2 / M`, `a = d_apo / (1 + e)`, `r_p = a (1 − e)`,
`p = a (1 − e²)`:
- `a_pc = a`, `apocentre_pc = d_apo`, `pericentre_pc = r_p` [pc]
- `apocentre_rh = d_apo / rh1`, `pericentre_rh = r_p / rh1`
- `period_myr = 2π √(a³ / (G M))`, converted from pc / (km s⁻¹) to Myr
- `e_orb = −G M1 M2 / (2a)`, orbital energy [M☉ km² s⁻²]
- `l_orb = μ √(G M p)`, orbital angular momentum [M☉ pc km s⁻¹]
- `j_orb = √(p / rh1)`, specific orbital angular momentum in units of
  `√(G M rh1)`
- `v_peri_kms = √(G M (1 + e) / r_p)`, relative speed at pericentre [km s⁻¹]
- `v_peri_over_sigma = v_peri_kms / sigma1`
- `lambda_peebles = l_orb √|e_total| / (G M^{5/2})`, the spin parameter of
  Peebles (1969, DOI 10.1086/149876) of the orbital angular momentum; with
  `l` in M☉ pc km s⁻¹, `E` in M☉ km² s⁻² and `G` in pc km² s⁻² M☉⁻¹ the
  numerator and the denominator both carry M☉^{3/2} pc km² s⁻², so it is
  dimensionless; `NaN` when `e_total` is `NaN`
- `roche_lobe_pc = (R_L,1, R_L,2)`, Roche-lobe radii at apocentre after
  Eggleton (1983, DOI 10.1086/160960),
  `R_L,i = d_apo 0.49 q^{2/3} / (0.6 q^{2/3} + ln(1 + q^{1/3}))` with
  `q = M_i / M_other` [pc]
- `roche_filling = (r_out[1] / R_L,1, r_out[2] / R_L,2)`

Throws an `ArgumentError` for non-positive masses, `d_apo` or `rh1`, and for
`e` outside `[0, 1)`.

# Example
```julia
q = two_cluster_orbit_quantities(1000.0, 1000.0, 12.0, 0.5, 2.0; sigma1 = 2.0)
q.pericentre_pc   # 4.0
```
"""
function two_cluster_orbit_quantities(
    M1::Real,
    M2::Real,
    d_apo::Real,
    ecc::Real,
    rh1::Real;
    sigma1::Real = NaN,
    e_total::Real = NaN,
    r_out::Tuple{Real,Real} = (NaN, NaN),
)::NamedTuple
    M1 > 0 || throw(ArgumentError("two_cluster_orbit_quantities: M1 must be > 0; got $M1"))
    M2 > 0 || throw(ArgumentError("two_cluster_orbit_quantities: M2 must be > 0; got $M2"))
    d_apo > 0 || throw(ArgumentError("two_cluster_orbit_quantities: d_apo must be > 0; got $d_apo"))
    rh1 > 0 || throw(ArgumentError("two_cluster_orbit_quantities: rh1 must be > 0; got $rh1"))
    0 ≤ ecc < 1 || throw(
        ArgumentError("two_cluster_orbit_quantities: eccentricity must be in [0, 1); got $ecc"),
    )
    G = _G_PC_KMS2_MSUN
    m1, m2, d, e, r_h = Float64(M1), Float64(M2), Float64(d_apo), Float64(ecc), Float64(rh1)
    M = m1 + m2
    μ = m1 * m2 / M
    a = d / (1 + e)
    r_p = a * (1 - e)
    p = a * (1 - e^2)
    l_orb = μ * sqrt(G * M * p)
    v_peri = sqrt(G * M * (1 + e) / r_p)
    E = Float64(e_total)
    λ = isnan(E) ? NaN : l_orb * sqrt(abs(E)) / (G * M^2.5)
    roche(q) = d * 0.49 * cbrt(q)^2 / (0.6 * cbrt(q)^2 + log(1 + cbrt(q)))
    R_L1, R_L2 = roche(m1 / m2), roche(m2 / m1)
    return (
        a_pc = a,
        apocentre_pc = d,
        pericentre_pc = r_p,
        apocentre_rh = d / r_h,
        pericentre_rh = r_p / r_h,
        period_myr = 2π * sqrt(a^3 / (G * M)) * (_PC_KM / _MYR_S),
        e_orb = -G * m1 * m2 / (2a),
        l_orb = l_orb,
        j_orb = sqrt(p / r_h),
        v_peri_kms = v_peri,
        v_peri_over_sigma = v_peri / Float64(sigma1),
        lambda_peebles = λ,
        roche_lobe_pc = (R_L1, R_L2),
        roche_filling = (Float64(r_out[1]) / R_L1, Float64(r_out[2]) / R_L2),
    )
end

"""
    combine_clusters_explicit(cluster_data, specs; truncate_jacobi_flag = true,
                              truncation_radii = Float64[])
        -> (pos, vel, mass, cluster_ranges, kept)

Combine N ≥ 2 clusters using explicit per-cluster position/velocity offsets.

Each `ClusterSpec` in `specs` must have `position` and `velocity` vectors set.
Particles in each cluster are shifted by the cluster's CM offset.

If `truncate_jacobi_flag=true`, each cluster is truncated at the Jacobi radius
with respect to its nearest neighbour before combining; when
`truncation_radii` is non-empty (one entry per cluster) and its `i`-th entry
is positive, cluster `i` is truncated at that radius instead.

Returns combined `(pos, vel, mass)` arrays and a vector of `UnitRange{Int}`
giving the particle index range for each cluster in the output.
"""
function combine_clusters_explicit(
    cluster_data::Vector{<:NamedTuple},
    specs::Vector{ClusterSpec};
    truncate_jacobi_flag::Bool = true,
    truncation_radii::AbstractVector{Float64} = Float64[],
)
    n = length(cluster_data)
    n == length(specs) || error("cluster_data and specs must have same length")
    (isempty(truncation_radii) || length(truncation_radii) == n) || throw(
        DimensionMismatch(
            "truncation_radii has $(length(truncation_radii)) entries for $n clusters",
        ),
    )

    # Optional Jacobi truncation against nearest neighbour
    truncated = Vector{
        NamedTuple{(:pos, :vel, :mass),Tuple{Matrix{Float64},Matrix{Float64},Vector{Float64}}},
    }()
    kept = Vector{Vector{Int}}()

    for i in 1:n
        # Keep only the particle arrays; sampled cluster data may carry
        # extra fields (energies) that the truncated vector does not hold.
        cd = (pos = cluster_data[i].pos, vel = cluster_data[i].vel, mass = cluster_data[i].mass)
        if truncate_jacobi_flag && !isempty(truncation_radii) && truncation_radii[i] > 0
            pos_t, vel_t, mass_t, idx_t =
                truncate_jacobi(cd.pos, cd.vel, cd.mass, Float64(truncation_radii[i]))
            push!(truncated, (pos = pos_t, vel = vel_t, mass = mass_t))
            push!(kept, idx_t)
        elseif truncate_jacobi_flag && n ≥ 2
            M_self = sum(cd.mass)
            pos_i = specs[i].position

            # Find nearest-neighbour distance
            d_min = Inf
            M_nearest = 0.0
            for j in 1:n
                j == i && continue
                d = sqrt(sum((pos_i[k] - specs[j].position[k])^2 for k in 1:3))
                if d < d_min
                    d_min = d
                    M_nearest = sum(cluster_data[j].mass)
                end
            end

            r_jacobi = jacobi_radius(d_min, M_self, M_nearest)
            pos_t, vel_t, mass_t, idx_t = truncate_jacobi(cd.pos, cd.vel, cd.mass, r_jacobi)
            push!(truncated, (pos = pos_t, vel = vel_t, mass = mass_t))
            push!(kept, idx_t)
        else
            push!(truncated, cd)
            push!(kept, collect(1:length(cd.mass)))
        end
    end

    # Count total particles and build ranges
    N_total = sum(length(cd.mass) for cd in truncated)
    pos_out = zeros(Float64, 3, N_total)
    vel_out = zeros(Float64, 3, N_total)
    mass_out = zeros(Float64, N_total)
    cluster_ranges = UnitRange{Int}[]

    offset = 0
    for i in 1:n
        cd = truncated[i]
        Ni = length(cd.mass)
        rng = (offset + 1):(offset + Ni)
        push!(cluster_ranges, rng)

        cm_pos = specs[i].position
        cm_vel = specs[i].velocity ./ _CODE_VSTAR_KMS   # km s⁻¹ → code units (G = 1)

        for j in 1:Ni
            idx = offset + j
            for k in 1:3
                pos_out[k, idx] = cd.pos[k, j] + cm_pos[k]
                vel_out[k, idx] = cd.vel[k, j] + cm_vel[k]
            end
            mass_out[idx] = cd.mass[j]
        end

        offset += Ni
    end

    # Shift to system centre of mass
    M_total = sum(mass_out)
    cm = zeros(3)
    cmv = zeros(3)
    for i in 1:N_total
        for k in 1:3
            cm[k] += mass_out[i] * pos_out[k, i]
            cmv[k] += mass_out[i] * vel_out[k, i]
        end
    end
    cm ./= M_total
    cmv ./= M_total
    for i in 1:N_total
        for k in 1:3
            pos_out[k, i] -= cm[k]
            vel_out[k, i] -= cmv[k]
        end
    end

    return pos_out, vel_out, mass_out, cluster_ranges, kept
end

"""
    combine_single_cluster(cluster_data::NamedTuple)
        -> (pos, vel, mass, cluster_ranges, kept)

The single-cluster counterpart of [`combine_clusters_explicit`](@ref): the
sampled cluster is shifted to its centre of mass (position and velocity,
mass-weighted) with no truncation. Returns copies of the arrays, the index
range `[1:N]` and `kept = [collect(1:N)]`.
"""
function combine_single_cluster(cluster_data::NamedTuple)
    pos = copy(cluster_data.pos)
    vel = copy(cluster_data.vel)
    mass = copy(cluster_data.mass)
    N = length(mass)

    # Shift to the centre of mass
    M_total = sum(mass)
    cm = zeros(3)
    cmv = zeros(3)
    for i in 1:N
        for k in 1:3
            cm[k] += mass[i] * pos[k, i]
            cmv[k] += mass[i] * vel[k, i]
        end
    end
    cm ./= M_total
    cmv ./= M_total
    for i in 1:N
        for k in 1:3
            pos[k, i] -= cm[k]
            vel[k, i] -= cmv[k]
        end
    end

    return pos, vel, mass, UnitRange{Int}[1:N], [collect(1:N)]
end

# Oort constants of the engine's solar-neighbourhood tide (xtrnl0.F, KZ(14) = 1) [km s⁻¹ pc⁻¹]
const _OORT_A_KMS_PC = 14.4e-3
const _OORT_B_KMS_PC = -12.0e-3

"""
    tidal_frame_rotation(tidal::TidalSpec) -> Float64

Angular velocity Ω [km s⁻¹ pc⁻¹] of the frame in which the engine integrates
under `tidal.kz14`, along +z: `A − B` with the Oort constants of `xtrnl0.F`
(A = 14.4, B = −12.0 km s⁻¹ kpc⁻¹) for `kz14 = 1`, `√(G M_G / R_G³)` of the
circular orbit in the point-mass galaxy for `kz14 = 2` (`gmg` [M☉], `rg0`
[kpc]); 0 for the inertial integrations (`kz14 = 0, 5`). The Coriolis terms
of `xtrnld.f` (F_x += 2Ω ẏ, F_y −= 2Ω ẋ) fix the sense: a pair orbit with its
angular momentum along +z is prograde with respect to the Galactic rotation.
"""
function tidal_frame_rotation(tidal::TidalSpec)::Float64
    tidal.kz14 == 1 && return _OORT_A_KMS_PC - _OORT_B_KMS_PC
    tidal.kz14 == 2 && return sqrt(_G_PC_KMS2_MSUN * tidal.gmg / (1000 * tidal.rg0)^3)
    return 0.0
end

"""
    pair_jacobi_radius(M, tidal::TidalSpec) -> Float64

Jacobi radius [pc] of a bound system of total mass `M` [M☉] in the
linearised tide of the engine (Hill problem): `(G M / (4A(A − B)))^{1/3}` for
`kz14 = 1` (≈ 1.41 (M/M☉)^{1/3} pc) and `(G M / (3Ω²))^{1/3}` for the
point-mass galaxy of `kz14 = 2` (A = 3Ω/4, B = −Ω/4); `Inf` without a
rotating-frame tide. A pair whose apocentre reaches it is not bound to itself.
"""
function pair_jacobi_radius(M::Real, tidal::TidalSpec)::Float64
    M > 0 || throw(ArgumentError("pair_jacobi_radius: M must be > 0; got $M"))
    if tidal.kz14 == 1
        A, B = _OORT_A_KMS_PC, _OORT_B_KMS_PC
        return cbrt(_G_PC_KMS2_MSUN * M / (4A * (A - B)))
    elseif tidal.kz14 == 2
        Ω = tidal_frame_rotation(tidal)
        return cbrt(_G_PC_KMS2_MSUN * M / (3Ω^2))
    end
    return Inf
end

"""
    to_rotating_frame!(pos, vel, Ω)

Replace inertial velocities `vel` [km s⁻¹] of bodies at `pos` [pc] (3 × N,
about the rotation axis through the origin) by the velocities in the frame
rotating with Ω ẑ [km s⁻¹ pc⁻¹]: v − Ω × r, i.e. vₓ + Ω y, v_y − Ω x.
"""
function to_rotating_frame!(pos::AbstractMatrix{Float64}, vel::AbstractMatrix{Float64}, Ω::Real)
    size(pos) == size(vel) || throw(DimensionMismatch("pos is $(size(pos)), vel is $(size(vel))"))
    for i in axes(pos, 2)
        vel[1, i] += Ω * pos[2, i]
        vel[2, i] -= Ω * pos[1, i]
    end
    return vel
end
