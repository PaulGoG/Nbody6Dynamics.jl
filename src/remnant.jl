# =============================================================================
# Remnant diagnostics: coalescence, core radius, rotation, mass segregation
# =============================================================================
#
# The "remnant" is the self-consistent bound set of the whole system at a
# snapshot (see `_bound_members`); before coalescence it is simply the bound
# part of the whole configuration. Cluster centre-of-mass trajectories and
# the union-find count of spatially distinct clusters live here because the
# coalescence time is derived from them; the merger figures reuse them.

using LinearAlgebra: dot, ×

"""
    _cluster_com_trajectories(snaps, cluster_ranges)
        -> (coms, present, r_rms)

Compute the mass-weighted centre of mass of each initial cluster at every
snapshot. `coms` is `3 × n_clusters × n_times`. `present[i, k]` is false when
fewer than 3 members of cluster `i` remain in snapshot `k`. `r_rms[i, k]` is
the mass-weighted RMS radius of cluster `i`'s members from their COM (a
proxy for cluster extent); `NaN` when `present` is false.
"""
function _cluster_com_trajectories(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
)
    n_cl = length(cluster_ranges)
    n_t = length(snaps)
    coms = fill(NaN, 3, n_cl, n_t)
    r_rms = fill(NaN, n_cl, n_t)
    present = falses(n_cl, n_t)

    for (k, snap) in enumerate(snaps)
        for (i, rng) in enumerate(cluster_ranges)
            mask = falses(length(snap.name))
            mask[_member_indices(snap, rng)] .= true
            n_mem = count(mask)
            n_mem < 3 && continue
            m = Float64.(snap.mass[mask])
            M = sum(m)
            M > 0 || continue
            cx = sum(m .* Float64.(snap.pos[1, mask])) / M
            cy = sum(m .* Float64.(snap.pos[2, mask])) / M
            cz = sum(m .* Float64.(snap.pos[3, mask])) / M
            coms[1, i, k] = cx
            coms[2, i, k] = cy
            coms[3, i, k] = cz

            # Mass-weighted RMS radius: sqrt(Σ m_k r_k² / Σ m_k)
            x = Float64.(@view snap.pos[1, mask])
            y = Float64.(@view snap.pos[2, mask])
            z = Float64.(@view snap.pos[3, mask])
            ssum = 0.0
            @inbounds for j in eachindex(m)
                dx = x[j] - cx
                dy = y[j] - cy
                dz = z[j] - cz
                ssum += m[j] * (dx*dx + dy*dy + dz*dz)
            end
            r_rms[i, k] = sqrt(ssum / M)
            present[i, k] = true
        end
    end
    return coms, present, r_rms
end

"""
    _count_spatial_clusters(coms, r_rms, present; overlap_factor = 1.0) -> Vector{Int}

At each snapshot, count spatially distinct initial clusters via union-find.
Two clusters are considered merged when their COMs are closer than
`overlap_factor × (r_ref_i + r_ref_j)`, where `r_ref` is the *initial* RMS
radius of each cluster (held fixed throughout). Using the initial radius
(rather than current) prevents spurious "mergers" caused by infall-driven
puffing — the criterion then measures true spatial convergence of COMs.
"""
function _count_spatial_clusters(
    coms::Array{Float64,3},
    r_rms::AbstractMatrix{<:Real},
    present::AbstractMatrix{Bool};
    overlap_factor::Real = 1.0,
)
    n_cl = size(coms, 2)
    n_t = size(coms, 3)
    counts = zeros(Int, n_t)

    # Reference radii: each cluster's RMS radius in its first present snapshot
    r_ref = fill(NaN, n_cl)
    for i in 1:n_cl
        for k in 1:n_t
            if present[i, k] && !isnan(r_rms[i, k])
                r_ref[i] = r_rms[i, k]
                break
            end
        end
    end

    # Simple union-find with path compression
    parent = zeros(Int, n_cl)
    function find(p::Vector{Int}, a::Int)
        while p[a] != a
            p[a] = p[p[a]]
            a = p[a]
        end
        return a
    end

    for k in 1:n_t
        active = findall(@view present[:, k])
        isempty(active) && (counts[k] = 0; continue)

        for i in 1:n_cl
            parent[i] = i
        end
        for ia in eachindex(active)
            i = active[ia]
            isnan(r_ref[i]) && continue
            for jb in (ia + 1):length(active)
                j = active[jb]
                isnan(r_ref[j]) && continue
                dx = coms[1, i, k] - coms[1, j, k]
                dy = coms[2, i, k] - coms[2, j, k]
                dz = coms[3, i, k] - coms[3, j, k]
                d = sqrt(dx*dx + dy*dy + dz*dz)
                if d < overlap_factor * (r_ref[i] + r_ref[j])
                    ri = find(parent, i)
                    rj = find(parent, j)
                    ri == rj || (parent[ri] = rj)
                end
            end
        end
        roots = Set{Int}()
        for i in active
            push!(roots, find(parent, i))
        end
        counts[k] = length(roots)
    end
    return counts
end

# ---------------------------------------------------------------------------
# Coalescence
# ---------------------------------------------------------------------------

"""
    coalescence_time(snaps, cluster_ranges; overlap_factor = 1.0)
        -> (; index, time_nb, time_myr, n_distinct)

First snapshot at which every initial cluster has joined one spatial group:
two clusters count as merged when their centres of mass are closer than
`overlap_factor` times the sum of their initial RMS radii
([`_count_spatial_clusters`](@ref)). `index` is `nothing` and the times are
`NaN` when the clusters never coalesce within the snapshots; `n_distinct`
is the count of spatially distinct clusters at every snapshot.
"""
function coalescence_time(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    overlap_factor::Real = 1.0,
)
    length(cluster_ranges) ≥ 2 ||
        throw(ArgumentError("coalescence_time needs at least two clusters"))
    isempty(snaps) && throw(ArgumentError("coalescence_time needs at least one snapshot"))
    coms, present, r_rms = _cluster_com_trajectories(snaps, cluster_ranges)
    n_distinct = _count_spatial_clusters(coms, r_rms, present; overlap_factor = overlap_factor)
    index = findfirst(k -> n_distinct[k] == 1 && count(@view present[:, k]) ≥ 2, eachindex(snaps))
    if index === nothing
        return (; index = nothing, time_nb = NaN, time_myr = NaN, n_distinct)
    end
    h = snaps[index].header
    return (; index, time_nb = time_nb(h), time_myr = time_myr(h), n_distinct)
end

"""
    _orbital_angular_momentum(snap, cluster_ranges) -> Vector{Float64}

Orbital angular momentum of the initial clusters about the system's centre
of mass, from their member centres of mass and mean velocities (NB units).
"""
function _orbital_angular_momentum(
    snap::Snapshot,
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
)
    n_cl = length(cluster_ranges)
    M = zeros(n_cl)
    R = zeros(3, n_cl)
    V = zeros(3, n_cl)
    for (i, rng) in enumerate(cluster_ranges)
        idx = _member_indices(snap, rng)
        isempty(idx) && continue
        m = Float64.(snap.mass[idx])
        M[i] = sum(m)
        for k in 1:3
            R[k, i] = sum(m[j] * snap.pos[k, idx[j]] for j in eachindex(idx)) / M[i]
            V[k, i] = sum(m[j] * snap.vel[k, idx[j]] for j in eachindex(idx)) / M[i]
        end
    end
    Mtot = sum(M)
    Rc = vec(sum(R .* M'; dims = 2)) ./ Mtot
    Vc = vec(sum(V .* M'; dims = 2)) ./ Mtot
    L = zeros(3)
    for i in 1:n_cl
        r = R[:, i] .- Rc
        v = V[:, i] .- Vc
        L .+= M[i] .* (r × v)
    end
    return L
end

# ---------------------------------------------------------------------------
# Local density and core radius (Casertano & Hut 1985)
# ---------------------------------------------------------------------------

"""Neighbour rank of the local density estimate (Casertano & Hut 1985)."""
const _DENSITY_NEIGHBOURS = 6

"""
    _local_density(pos, mass; k = _DENSITY_NEIGHBOURS) -> Vector{Float64}

Casertano & Hut (1985) density estimate at every particle: the mass of
its `k − 1` nearest neighbours divided by the volume of the sphere reaching
the `k`-th. O(N²).
"""
function _local_density(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    k::Int = _DENSITY_NEIGHBOURS,
)
    n = length(mass)
    n > k || throw(ArgumentError("_local_density needs more than $k particles, got $n"))
    ρ = zeros(n)
    d2 = Vector{Float64}(undef, n)
    order = Vector{Int}(undef, n)
    @inbounds for i in 1:n
        for j in 1:n
            dx = pos[1, i] - pos[1, j]
            dy = pos[2, i] - pos[2, j]
            dz = pos[3, i] - pos[3, j]
            d2[j] = dx * dx + dy * dy + dz * dz
        end
        d2[i] = Inf
        partialsortperm!(order, d2, 1:k)
        m_inner = 0.0
        for j in 1:(k - 1)
            m_inner += mass[order[j]]
        end
        r_k = sqrt(d2[order[k]])
        ρ[i] = r_k > 0 ? m_inner / (4π / 3 * r_k^3) : 0.0
    end
    return ρ
end

"""
    _snapshot_density(snap, idx) -> Vector{Float64}

Local densities of the particles `idx`: the engine's values when the
snapshot carries a density column with positive entries for most of them,
otherwise [`_local_density`](@ref) on the subset.
"""
function _snapshot_density(snap::Snapshot, idx::AbstractVector{Int})
    if length(snap.rho) == nparticles(snap)
        ρ = Float64.(snap.rho[idx])
        count(>(0), ρ) ≥ length(idx) ÷ 2 && return ρ
    end
    return _local_density(Float64.(snap.pos[:, idx]), Float64.(snap.mass[idx]))
end

"""
    core_radius(pos, rho) -> (; r_core, centre)

Density-weighted centre `Σ ρ_i r_i / Σ ρ_i` and core radius
`√(Σ ρ_i² |r_i − r_d|² / Σ ρ_i²)` of Casertano & Hut (1985), the
definition behind the engine's `RC`. Summed over particles the weights are
`ρ³` in the continuum limit, so for a Plummer sphere of scale length `a`
the estimator gives `√0.3 a ≈ 0.548 a` with exact densities; the
sixth-neighbour estimate lies ≈ 4 % below that.
"""
function core_radius(pos::AbstractMatrix{Float64}, rho::AbstractVector{Float64})
    length(rho) == size(pos, 2) || throw(DimensionMismatch("rho and pos disagree in length"))
    w = sum(rho)
    w > 0 || throw(ArgumentError("core_radius: densities must not all vanish"))
    centre = [sum(rho[i] * pos[k, i] for i in eachindex(rho)) / w for k in 1:3]
    w2 = sum(abs2, rho)
    s = 0.0
    for i in eachindex(rho)
        s +=
            rho[i]^2 *
            ((pos[1, i] - centre[1])^2 + (pos[2, i] - centre[2])^2 + (pos[3, i] - centre[3])^2)
    end
    return (; r_core = sqrt(s / w2), centre)
end

# ---------------------------------------------------------------------------
# Rotation
# ---------------------------------------------------------------------------

"""
    RotationProfile

Rotation of a stellar system about its spin axis in cylindrical shells of
equal membership: `radius` (mean cylindrical radius per shell), `v_rot`
(mass-weighted mean azimuthal velocity), `sigma` (one-dimensional
dispersion of the shell about its mean velocity), all in NB units, and the
ratio `v_rot_over_sigma`.
"""
struct RotationProfile
    radius::Vector{Float64}
    v_rot::Vector{Float64}
    sigma::Vector{Float64}
    v_rot_over_sigma::Vector{Float64}
end

"""Cylindrical shells of the rotation profile."""
const _ROTATION_SHELLS = 8

"""
    rotation_analysis(pos, vel, mass; nbins = _ROTATION_SHELLS)
        -> (; L, axis, lambda_r, lambda_peebles, profile)

Rotation of a bound set about its shrinking-sphere centre and mean velocity:
angular momentum `L`, unit spin `axis`, the intrinsic ordered-to-total
motion parameter
``λ_R = Σ m R |v_rot| / Σ m R √(v_rot² + σ²)`` (the λ_R of Emsellem et al.
2007 evaluated with the shell values of the [`RotationProfile`](@ref)),
the Peebles (1969) spin parameter ``λ_P = |L| √|E| / (G M^{5/2})`` with
`E` the total energy of the set (`G = 1`), and the profile itself.
"""
function rotation_analysis(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    nbins::Int = _ROTATION_SHELLS,
)
    n = length(mass)
    n ≥ _MIN_MEMBERS ||
        throw(ArgumentError("rotation_analysis needs at least $(_MIN_MEMBERS) particles, got $n"))
    nbins ≥ 1 || throw(ArgumentError("nbins must be positive"))
    M = sum(mass)
    c = _shrinking_sphere_centre(pos, mass)
    vc = vec(sum(vel .* mass'; dims = 2)) ./ M
    L = zeros(3)
    for i in 1:n
        r = (pos[1, i] - c[1], pos[2, i] - c[2], pos[3, i] - c[3])
        v = (vel[1, i] - vc[1], vel[2, i] - vc[2], vel[3, i] - vc[3])
        L[1] += mass[i] * (r[2] * v[3] - r[3] * v[2])
        L[2] += mass[i] * (r[3] * v[1] - r[1] * v[3])
        L[3] += mass[i] * (r[1] * v[2] - r[2] * v[1])
    end
    Lnorm = sqrt(sum(abs2, L))
    axis = Lnorm > 0 ? L ./ Lnorm : [0.0, 0.0, 1.0]

    # Cylindrical radius and azimuthal velocity of every particle
    R = zeros(n)
    vφ = zeros(n)
    for i in 1:n
        r = [pos[k, i] - c[k] for k in 1:3]
        v = [vel[k, i] - vc[k] for k in 1:3]
        r_par = dot(r, axis)
        r_perp = r .- r_par .* axis
        Ri = sqrt(sum(abs2, r_perp))
        R[i] = Ri
        if Ri > 0
            φ̂ = axis × (r_perp ./ Ri)
            vφ[i] = dot(v, φ̂)
        end
    end

    # Shells of equal membership in cylindrical radius
    order = sortperm(R)
    nb = min(nbins, n ÷ max(_MIN_MEMBERS, 1))
    nb = max(nb, 1)
    edges = round.(Int, range(0, n; length = nb + 1))
    radius = zeros(nb)
    v_rot = zeros(nb)
    sigma = zeros(nb)
    shell_of = zeros(Int, n)
    for b in 1:nb
        idx = order[(edges[b] + 1):edges[b + 1]]
        shell_of[idx] .= b
        m = mass[idx]
        Mb = sum(m)
        radius[b] = sum(m .* R[idx]) / Mb
        v_rot[b] = sum(m .* vφ[idx]) / Mb
        vmean = [sum(m[j] * (vel[k, idx[j]] - vc[k]) for j in eachindex(idx)) / Mb for k in 1:3]
        s2 = 0.0
        for j in eachindex(idx)
            i = idx[j]
            s2 += m[j] * sum((vel[k, i] - vc[k] - vmean[k])^2 for k in 1:3)
        end
        sigma[b] = sqrt(s2 / (3 * Mb))
    end
    ratio = [sigma[b] > 0 ? v_rot[b] / sigma[b] : NaN for b in 1:nb]
    profile = RotationProfile(radius, v_rot, sigma, ratio)

    num = 0.0
    den = 0.0
    for i in 1:n
        b = shell_of[i]
        num += mass[i] * R[i] * abs(v_rot[b])
        den += mass[i] * R[i] * sqrt(v_rot[b]^2 + sigma[b]^2)
    end
    lambda_r = den > 0 ? num / den : NaN

    T = 0.0
    for i in 1:n
        T += 0.5 * mass[i] * sum((vel[k, i] - vc[k])^2 for k in 1:3)
    end
    W = 0.5 * sum(mass .* _self_potential(pos, mass))
    E = T + W
    lambda_peebles = Lnorm * sqrt(abs(E)) / M^2.5

    return (; L, axis, lambda_r, lambda_peebles, profile)
end

# ---------------------------------------------------------------------------
# Mass segregation (Allison et al. 2009)
# ---------------------------------------------------------------------------

"""Number of massive stars and of random comparison sets of the Λ_MSR estimate."""
const _MSR_N_MASSIVE = 50
const _MSR_N_RANDOM = 50

"""Length of the Euclidean minimum spanning tree of the points `pos[:, idx]` (Prim, O(n²))."""
function _mst_length(pos::AbstractMatrix{Float64}, idx::AbstractVector{Int})
    n = length(idx)
    n ≤ 1 && return 0.0
    in_tree = falses(n)
    best = fill(Inf, n)
    in_tree[1] = true
    for j in 2:n
        best[j] = sqrt(sum((pos[k, idx[1]] - pos[k, idx[j]])^2 for k in 1:3))
    end
    total = 0.0
    for _ in 2:n
        j_min = 0
        d_min = Inf
        for j in 1:n
            if !in_tree[j] && best[j] < d_min
                d_min = best[j]
                j_min = j
            end
        end
        in_tree[j_min] = true
        total += d_min
        for j in 1:n
            in_tree[j] && continue
            d = sqrt(sum((pos[k, idx[j_min]] - pos[k, idx[j]])^2 for k in 1:3))
            d < best[j] && (best[j] = d)
        end
    end
    return total
end

"""
    mass_segregation(pos, mass; n_massive = 50, n_random = 50, seed = 0)
        -> (; lambda_msr, lambda_err, r_half_ratio)

Mass-segregation ratio of Allison et al. (2009): the mean minimum-spanning-
tree length of `n_random` random sets of `n_massive` stars over that of the
`n_massive` most massive stars, `Λ_MSR = ⟨l_random⟩ / l_massive`, with the
error `σ_random / l_massive`; `Λ_MSR ≈ 1` without segregation, `> 1` with
the massive stars concentrated. `r_half_ratio` is the half-mass radius of
the massive subset over that of the whole set, both about the
shrinking-sphere centre. The random sets come from a `MersenneTwister`
seeded with `seed`, so the estimate is reproducible.
"""
function mass_segregation(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    n_massive::Int = _MSR_N_MASSIVE,
    n_random::Int = _MSR_N_RANDOM,
    seed::Integer = 0,
)
    n = length(mass)
    n_massive ≥ 2 || throw(ArgumentError("n_massive must be ≥ 2"))
    n_random ≥ 2 || throw(ArgumentError("n_random must be ≥ 2"))
    n > n_massive ||
        throw(ArgumentError("mass_segregation needs more than $n_massive particles, got $n"))
    massive = partialsortperm(mass, 1:n_massive; rev = true)
    l_massive = _mst_length(pos, massive)
    rng = Random.MersenneTwister(seed)
    lengths = [_mst_length(pos, Random.randperm(rng, n)[1:n_massive]) for _ in 1:n_random]
    mean_l = sum(lengths) / n_random
    std_l = sqrt(sum((lengths .- mean_l) .^ 2) / (n_random - 1))
    lambda_msr = l_massive > 0 ? mean_l / l_massive : NaN
    lambda_err = l_massive > 0 ? std_l / l_massive : NaN
    c = _shrinking_sphere_centre(pos, mass)
    r_half_all = _lagrangian_radii(pos, mass, c, (0.5,))[1]
    r_half_massive = _lagrangian_radii(pos[:, massive], mass[massive], c, (0.5,))[1]
    return (; lambda_msr, lambda_err, r_half_ratio = r_half_massive / r_half_all)
end

# ---------------------------------------------------------------------------
# Time series
# ---------------------------------------------------------------------------

"""
    RemnantDiagnostics

Diagnostics of the bound remnant of the whole system at every snapshot
(NB units unless noted), from [`remnant_diagnostics`](@ref).

# Fields
- `time`, `time_myr`: snapshot times
- `rbar`: NB length unit in pc at every snapshot (1 when unscaled)
- `n_bound`, `bound_mass_fraction`: bound members and their mass fraction
- `r_core`, `r_half`: Casertano–Hut core radius and half-mass radius
- `lambda_r`, `lambda_peebles`: rotation parameters ([`rotation_analysis`](@ref))
- `spin_axis`: unit spin vectors, `3 × n_t`
- `spin_alignment`: cosine of the angle between the spin axis and the
  orbital angular momentum of the initial clusters at the first snapshot
- `lambda_msr`, `lambda_msr_err`, `segregation_ratio`: mass segregation
  ([`mass_segregation`](@ref))
- `coalescence_time`, `coalescence_time_myr`: first snapshot at which all
  initial clusters overlap ([`coalescence_time`](@ref)); `NaN` if never
- `segregation_time`, `segregation_time_myr`: first snapshot with
  `Λ_MSR − σ > lambda_threshold`; `NaN` if never
- `profile`: rotation profile of the last snapshot with enough members
"""
struct RemnantDiagnostics
    time::Vector{Float64}
    time_myr::Vector{Float64}
    rbar::Vector{Float64}
    n_bound::Vector{Int}
    bound_mass_fraction::Vector{Float64}
    r_core::Vector{Float64}
    r_half::Vector{Float64}
    lambda_r::Vector{Float64}
    lambda_peebles::Vector{Float64}
    spin_axis::Matrix{Float64}
    spin_alignment::Vector{Float64}
    lambda_msr::Vector{Float64}
    lambda_msr_err::Vector{Float64}
    segregation_ratio::Vector{Float64}
    coalescence_time::Float64
    coalescence_time_myr::Float64
    segregation_time::Float64
    segregation_time_myr::Float64
    profile::RotationProfile
end

"""Λ_MSR level (minus its error) above which the remnant counts as mass-segregated."""
const _MSR_THRESHOLD = 2.0

"""
    remnant_diagnostics(snaps, cluster_ranges; n_massive = 50, n_random = 50,
                        seed = 0, lambda_threshold = 2.0, nbins = 8,
                        overlap_factor = 1.0) -> RemnantDiagnostics

Bound remnant of the whole system at every snapshot (the self-consistent
bound set of all particles), its core and half-mass radii, rotation and
mass segregation, plus the coalescence and segregation times. O(N²) per
snapshot; intended for N ≲ 2×10⁴.
"""
function remnant_diagnostics(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    n_massive::Int = _MSR_N_MASSIVE,
    n_random::Int = _MSR_N_RANDOM,
    seed::Integer = 0,
    lambda_threshold::Real = _MSR_THRESHOLD,
    nbins::Int = _ROTATION_SHELLS,
    overlap_factor::Real = 1.0,
)::RemnantDiagnostics
    isempty(snaps) && throw(ArgumentError("remnant_diagnostics needs at least one snapshot"))
    n_t = length(snaps)
    t = [time_nb(s.header) for s in snaps]
    t_myr = [time_myr(s.header) for s in snaps]
    rbars = [rbar(s.header) for s in snaps]
    n_bound = zeros(Int, n_t)
    fbound = fill(NaN, n_t)
    r_core = fill(NaN, n_t)
    r_half = fill(NaN, n_t)
    λr = fill(NaN, n_t)
    λp = fill(NaN, n_t)
    axis = fill(NaN, 3, n_t)
    align = fill(NaN, n_t)
    λmsr = fill(NaN, n_t)
    λerr = fill(NaN, n_t)
    seg = fill(NaN, n_t)
    profile = RotationProfile(Float64[], Float64[], Float64[], Float64[])

    L_orb = _orbital_angular_momentum(snaps[1], cluster_ranges)
    L_norm = sqrt(sum(abs2, L_orb))
    L_hat = L_norm > 0 ? L_orb ./ L_norm : fill(NaN, 3)

    for (k, snap) in enumerate(snaps)
        pos_all = Float64.(snap.pos)
        vel_all = Float64.(snap.vel)
        mass_all = Float64.(snap.mass)
        length(mass_all) < _MIN_MEMBERS && continue
        sel = _bound_members(pos_all, vel_all, mass_all)
        n_bound[k] = length(sel)
        fbound[k] = sum(mass_all[sel]) / sum(mass_all)
        length(sel) < _MIN_MEMBERS && continue
        pos = pos_all[:, sel]
        vel = vel_all[:, sel]
        mass = mass_all[sel]

        if length(sel) > _DENSITY_NEIGHBOURS
            ρ = _snapshot_density(snap, sel)
            cr = core_radius(pos, ρ)
            r_core[k] = cr.r_core
        end
        c = _shrinking_sphere_centre(pos, mass)
        r_half[k] = _lagrangian_radii(pos, mass, c, (0.5,))[1]

        rot = rotation_analysis(pos, vel, mass; nbins = nbins)
        λr[k] = rot.lambda_r
        λp[k] = rot.lambda_peebles
        axis[:, k] = rot.axis
        align[k] = all(isfinite, L_hat) ? dot(rot.axis, L_hat) : NaN
        profile = rot.profile

        if length(sel) > n_massive
            ms = mass_segregation(pos, mass; n_massive, n_random, seed)
            λmsr[k] = ms.lambda_msr
            λerr[k] = ms.lambda_err
            seg[k] = ms.r_half_ratio
        end
    end

    coal =
        length(cluster_ranges) ≥ 2 ?
        coalescence_time(snaps, cluster_ranges; overlap_factor = overlap_factor) :
        (; index = nothing, time_nb = NaN, time_myr = NaN)
    k_seg = findfirst(k -> isfinite(λmsr[k]) && λmsr[k] - λerr[k] > lambda_threshold, 1:n_t)
    t_seg = k_seg === nothing ? NaN : t[k_seg]
    t_seg_myr = k_seg === nothing ? NaN : t_myr[k_seg]

    return RemnantDiagnostics(
        t,
        t_myr,
        rbars,
        n_bound,
        fbound,
        r_core,
        r_half,
        λr,
        λp,
        axis,
        align,
        λmsr,
        λerr,
        seg,
        coal.time_nb,
        coal.time_myr,
        t_seg,
        t_seg_myr,
        profile,
    )
end

"""
    write_remnant_diagnostics(path, diag::RemnantDiagnostics) -> String

Write the time series of `diag` as CSV (one row per snapshot, NB units
plus `time_myr`) with the coalescence and segregation times in a leading
comment line. An existing file is backed up, never overwritten.
"""
function write_remnant_diagnostics(path::AbstractString, diag::RemnantDiagnostics)
    _backup_existing(path)
    open(path, "w") do io
        println(
            io,
            "# coalescence_time_nb=$(diag.coalescence_time) coalescence_time_myr=$(diag.coalescence_time_myr) segregation_time_nb=$(diag.segregation_time) segregation_time_myr=$(diag.segregation_time_myr)",
        )
        println(
            io,
            "time_nb,time_myr,rbar_pc,n_bound,bound_mass_fraction,r_core,r_half,lambda_r,lambda_peebles,spin_x,spin_y,spin_z,spin_alignment,lambda_msr,lambda_msr_err,segregation_ratio",
        )
        for k in eachindex(diag.time)
            cells = [
                diag.time[k],
                diag.time_myr[k],
                diag.rbar[k],
                diag.n_bound[k],
                diag.bound_mass_fraction[k],
                diag.r_core[k],
                diag.r_half[k],
                diag.lambda_r[k],
                diag.lambda_peebles[k],
                diag.spin_axis[1, k],
                diag.spin_axis[2, k],
                diag.spin_axis[3, k],
                diag.spin_alignment[k],
                diag.lambda_msr[k],
                diag.lambda_msr_err[k],
                diag.segregation_ratio[k],
            ]
            println(io, join(_csv_cell.(cells), ","))
        end
    end
    return String(path)
end
