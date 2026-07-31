# =============================================================================
# Orbital mechanics and multi-cluster assembly
# =============================================================================

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
coordinate origin removed. Non-mutating (the inputs are left untouched), so
no `!` — the cluster is assumed centred at the origin (pre-offset).
"""
function truncate_jacobi(pos::Matrix{Float64}, vel::Matrix{Float64},
                         mass::Vector{Float64}, r_trunc::Float64)
    N = length(mass)
    keep = Bool[]
    sizehint!(keep, N)
    for i in 1:N
        r = sqrt(pos[1, i]^2 + pos[2, i]^2 + pos[3, i]^2)
        push!(keep, r ≤ r_trunc)
    end
    idx = findall(keep)
    n_removed = N - length(idx)
    n_removed > 0 && @info "Jacobi truncation: removed $n_removed / $N particles (r_trunc = $(round(r_trunc; digits=4)))"
    return pos[:, idx], vel[:, idx], mass[idx]
end

"""
    setup_two_cluster_orbit(pos1, vel1, mass1, pos2, vel2, mass2,
                            d_apo, ecc; truncate_jacobi_flag=true)
        -> (pos, vel, mass)

Place two clusters on a Keplerian orbit at apocentre.

Cluster 1 is placed at `(-d1, 0, 0)` with velocity `(0, +v1, 0)`,
cluster 2 at `(+d2, 0, 0)` with velocity `(0, -v2, 0)`, where
`d1, d2` are the CM distances and `v1, v2` are the apocentre speeds.

If `truncate_jacobi_flag=true`, each cluster is truncated at its
instantaneous Jacobi radius before combining.

Returns combined `(pos, vel, mass)` arrays in the system CM frame.
"""
function setup_two_cluster_orbit(
    pos1::Matrix{Float64}, vel1::Matrix{Float64}, mass1::Vector{Float64},
    pos2::Matrix{Float64}, vel2::Matrix{Float64}, mass2::Vector{Float64},
    d_apo::Float64, ecc::Float64;
    truncate_jacobi_flag::Bool = true
)
    M1 = sum(mass1)
    M2 = sum(mass2)
    M_total = M1 + M2

    # Optional Jacobi truncation
    if truncate_jacobi_flag
        rJ1 = jacobi_radius(d_apo, M1, M2)
        rJ2 = jacobi_radius(d_apo, M2, M1)
        pos1, vel1, mass1 = truncate_jacobi(pos1, vel1, mass1, rJ1)
        pos2, vel2, mass2 = truncate_jacobi(pos2, vel2, mass2, rJ2)
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

    # Cluster 1: shift to (-d1, 0, 0) with velocity (0, +v1, 0)
    for i in 1:N1
        pos_out[1, i] = pos1[1, i] - d1
        pos_out[2, i] = pos1[2, i]
        pos_out[3, i] = pos1[3, i]
        vel_out[1, i] = vel1[1, i]
        vel_out[2, i] = vel1[2, i] + v1
        vel_out[3, i] = vel1[3, i]
        mass_out[i] = mass1[i]
    end

    # Cluster 2: shift to (+d2, 0, 0) with velocity (0, -v2, 0)
    for i in 1:N2
        j = N1 + i
        pos_out[1, j] = pos2[1, i] + d2
        pos_out[2, j] = pos2[2, i]
        pos_out[3, j] = pos2[3, i]
        vel_out[1, j] = vel2[1, i]
        vel_out[2, j] = vel2[2, i] - v2
        vel_out[3, j] = vel2[3, i]
        mass_out[j] = mass2[i]
    end

    return pos_out, vel_out, mass_out
end

"""
    combine_clusters_explicit(cluster_data, specs; truncate_jacobi_flag=true)
        -> (pos, vel, mass, cluster_ranges)

Combine N ≥ 2 clusters using explicit per-cluster position/velocity offsets.

Each `ClusterSpec` in `specs` must have `position` and `velocity` vectors set.
Particles in each cluster are shifted by the cluster's CM offset.

If `truncate_jacobi_flag=true`, each cluster is truncated at the Jacobi radius
with respect to its nearest neighbour before combining.

Returns combined `(pos, vel, mass)` arrays and a vector of `UnitRange{Int}`
giving the particle index range for each cluster in the output.
"""
function combine_clusters_explicit(
    cluster_data::Vector{<:NamedTuple},
    specs::Vector{ClusterSpec};
    truncate_jacobi_flag::Bool = true
)
    n = length(cluster_data)
    n == length(specs) || error("cluster_data and specs must have same length")

    # Optional Jacobi truncation against nearest neighbour
    truncated = Vector{NamedTuple{(:pos, :vel, :mass), Tuple{Matrix{Float64}, Matrix{Float64}, Vector{Float64}}}}()

    for i in 1:n
        cd = cluster_data[i]
        if truncate_jacobi_flag && n ≥ 2
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

            rJ = jacobi_radius(d_min, M_self, M_nearest)
            pos_t, vel_t, mass_t = truncate_jacobi(cd.pos, cd.vel, cd.mass, rJ)
            push!(truncated, (pos = pos_t, vel = vel_t, mass = mass_t))
        else
            push!(truncated, cd)
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
        cm_vel = specs[i].velocity

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
            cm[k]  += mass_out[i] * pos_out[k, i]
            cmv[k] += mass_out[i] * vel_out[k, i]
        end
    end
    cm  ./= M_total
    cmv ./= M_total
    for i in 1:N_total
        for k in 1:3
            pos_out[k, i] -= cm[k]
            vel_out[k, i] -= cmv[k]
        end
    end

    return pos_out, vel_out, mass_out, cluster_ranges
end
