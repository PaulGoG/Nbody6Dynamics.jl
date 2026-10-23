# =============================================================================
# Remnant diagnostics: coalescence, core radius, rotation, mass segregation
# =============================================================================
#
# The "remnant" is the self-consistent bound set of the whole system at a
# snapshot (see `bound_set`); before coalescence it is simply the bound
# part of the whole configuration. The coalescence time follows from the
# density-centre separations of the surviving per-cluster memberships of the
# global bound set.

using LinearAlgebra: dot, ×

# ---------------------------------------------------------------------------
# Coalescence
# ---------------------------------------------------------------------------

"""Root of `a` in the union–find forest `parent`, with path halving."""
function _union_root!(parent::Vector{Int}, a::Int)::Int
    @inbounds while parent[a] != a
        parent[a] = parent[parent[a]]
        a = parent[a]
    end
    return a
end

"""
    coalescence_time(snaps, cluster_ranges; bound_sets = nothing, max_passes = 100, min_members = 10)
        -> (; index, time_nb, time_myr, n_distinct::Vector{Int},
             separation::Matrix{Float64}, threshold::Matrix{Float64},
             pairs::Vector{Tuple{Int,Int}}, present::BitMatrix)

First snapshot from which the initial clusters form a single spatial group
for the rest of the run.

At every snapshot each cluster (tag of `cluster_ranges`) is represented by
its members that belong to the bound set of the whole system
([`bound_set`](@ref) of the snapshot with `max_passes` and `min_members`, or
the entries of `bound_sets`, one per snapshot, when the caller has them); a
cluster with fewer than `2 × _MIN_MEMBERS` bound members is absent at that
snapshot. Its centre is the shrinking-sphere density centre of those members
and `r_h` their half-mass radius about it. The members are taken from the
global bound set, not from a bound set of the tag alone: after the merger the
stars of the smaller cluster are bound to the remnant, not to each other, and
a per-tag unbinding would strip them away. Two present
clusters are merged when their centres are closer than
`min(r_h,i, r_h,j)`; `n_distinct` counts the connected components of the
merged-pair graph over the present clusters (union–find). The clusters have
coalesced at snapshot `k` when `n_distinct == 1` with at least two clusters
present there and `n_distinct == 1` at every later snapshot, so a first
passage through the pericentre of an orbit that separates again does not
count (Fellhauer et al. 2002, doi:10.1046/j.1365-8711.2002.05143.x).

The criterion replaces one based on member centres of mass and RMS radii:
the RMS radius of a Plummer sphere diverges logarithmically with the
truncation radius, and the member centre of mass follows the stripped
debris rather than the cluster.

`index` is `nothing` and the times are `NaN` when the clusters never
coalesce within the snapshots. `separation` and `threshold` are
`n_pairs × n_t` in the pair order `(1, 2), (1, 3), …, (2, 3), …`, `NaN`
where a cluster of the pair is absent. One bound-set evaluation per snapshot
unless `bound_sets` is given.
"""
function coalescence_time(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    bound_sets::Union{Nothing,AbstractVector{BoundSet}} = nothing,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)
    n_cl = length(cluster_ranges)
    n_cl ≥ 2 || throw(ArgumentError("coalescence_time needs at least two clusters"))
    isempty(snaps) && throw(ArgumentError("coalescence_time needs at least one snapshot"))
    bound_sets === nothing ||
        length(bound_sets) == length(snaps) ||
        throw(ArgumentError("coalescence_time: one bound set per snapshot is needed"))
    n_t = length(snaps)
    pairs = [(i, j) for i in 1:(n_cl - 1) for j in (i + 1):n_cl]
    separation = fill(NaN, length(pairs), n_t)
    threshold = fill(NaN, length(pairs), n_t)
    n_distinct = zeros(Int, n_t)
    n_present = zeros(Int, n_t)
    present_all = falses(n_cl, n_t)
    centres = zeros(3, n_cl)
    r_h = fill(NaN, n_cl)
    present = falses(n_cl)
    parent = collect(1:n_cl)

    for (k, snap) in enumerate(snaps)
        fill!(present, false)
        bound =
            bound_sets === nothing ? bound_set(snap; max_passes, min_members).stars :
            bound_sets[k].stars
        in_bound = falses(nparticles(snap))
        in_bound[bound] .= true
        for (i, rng) in enumerate(cluster_ranges)
            sel = [s for s in _member_indices(snap, rng) if in_bound[s]]
            length(sel) < 2 * _MIN_MEMBERS && continue
            p = Float64.(snap.pos[:, sel])
            m = Float64.(snap.mass[sel])
            c = _shrinking_sphere_centre(p, m)
            centres[:, i] .= c
            r_h[i] = _lagrangian_radii(p, m, c, (0.5,))[1]
            present[i] = true
        end
        for i in 1:n_cl
            parent[i] = i
        end
        for (q, (i, j)) in enumerate(pairs)
            (present[i] && present[j]) || continue
            dx = centres[1, i] - centres[1, j]
            dy = centres[2, i] - centres[2, j]
            dz = centres[3, i] - centres[3, j]
            d = sqrt(dx * dx + dy * dy + dz * dz)
            thr = min(r_h[i], r_h[j])
            separation[q, k] = d
            threshold[q, k] = thr
            if d < thr
                ri = _union_root!(parent, i)
                rj = _union_root!(parent, j)
                ri == rj || (parent[ri] = rj)
            end
        end
        n_present[k] = count(present)
        present_all[:, k] .= present
        n_distinct[k] = count(i -> present[i] && _union_root!(parent, i) == i, 1:n_cl)
    end

    # Persistence: scan back from the last snapshot while a single group remains
    index = nothing
    for k in n_t:-1:1
        n_distinct[k] == 1 || break
        n_present[k] ≥ 2 && (index = k)
    end
    if index === nothing
        return (;
            index = nothing,
            time_nb = NaN,
            time_myr = NaN,
            n_distinct,
            separation,
            threshold,
            pairs,
            present = present_all,
        )
    end
    h = snaps[index].header
    return (;
        index,
        time_nb = time_nb(h),
        time_myr = time_myr(h),
        n_distinct,
        separation,
        threshold,
        pairs,
        present = present_all,
    )
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
    _local_density(pos, mass; k = _DENSITY_NEIGHBOURS, weighting = :mass) -> Vector{Float64}

Casertano & Hut (1985, doi:10.1086/163485) density estimate at every
particle: the mass of its `k − 1` nearest neighbours divided by the volume
of the sphere reaching the `k`-th (`weighting = :mass`), or the number
density `(k − 1) / (4π/3 r_k³)` (`weighting = :number`, the proxy behind
the engine's `XNDBL`). O(N²).
"""
function _local_density(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    k::Int = _DENSITY_NEIGHBOURS,
    weighting::Symbol = :mass,
)
    weighting in (:mass, :number) ||
        throw(ArgumentError("weighting must be one of :mass | :number, got :$weighting"))
    number = weighting === :number
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
        if number
            m_inner = Float64(k - 1)
        else
            for j in 1:(k - 1)
                m_inner += mass[order[j]]
            end
        end
        r_k = sqrt(d2[order[k]])
        ρ[i] = r_k > 0 ? m_inner / (4π / 3 * r_k^3) : 0.0
    end
    return ρ
end

"""
    _snapshot_density(snap, idx; weighting = :mass, k = _DENSITY_NEIGHBOURS) -> Vector{Float64}

Local densities of the particles `idx`. Mass densities are the engine's
values when the snapshot carries a density column with positive entries for
most of them, otherwise [`_local_density`](@ref) on the subset with
neighbour rank `k`; number densities (`weighting = :number`) are always
computed locally, the engine column being a mass density.
"""
function _snapshot_density(
    snap::Snapshot,
    idx::AbstractVector{Int};
    weighting::Symbol = :mass,
    k::Int = _DENSITY_NEIGHBOURS,
)
    if weighting === :mass && length(snap.rho) == nparticles(snap)
        ρ = Float64.(snap.rho[idx])
        count(>(0), ρ) ≥ length(idx) ÷ 2 && return ρ
    end
    return local_density(Float64.(snap.pos[:, idx]), Float64.(snap.mass[idx]); k, weighting)
end

"""
    core_radius(pos, rho) -> (; r_core, centre)

Density-weighted centre `Σ ρ_i r_i / Σ ρ_i` and core radius
`√(Σ ρ_i² |r_i − r_d|² / Σ ρ_i²)` of Casertano & Hut (1985), the
definition behind the engine's `RC`. Summed over particles the weights are
`ρ³` in the continuum limit, so for a Plummer sphere of scale length `a`
the estimator gives `√0.3 a ≈ 0.548 a` with exact densities. With
sixth-neighbour number densities the calibration table
(`assets/calibration/estimators.csv`) finds it unbiased from N ≈ 3×10⁴ on
and 1–2 % high below (replica scatter 9 % at 10³, 0.3 % at 10⁶); with
sixth-neighbour mass densities and a Kroupa spectrum it is 1.5–3 % high
with ten times the scatter.
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
equal membership, all in NB units: `radius` (mass-weighted mean cylindrical
radius per shell), `v_rot` (mass-weighted mean azimuthal velocity), `sigma`
(one-dimensional dispersion of the shell about its mean velocity), the
ratio `v_rot_over_sigma`, `n` (members per shell), `sigma_phi`
(one-dimensional azimuthal dispersion of the shell about `v_rot`),
`v_rot_debiased` (``\\mathrm{sign}(\\bar v_φ) \\sqrt{\\max(\\bar v_φ^2 − σ_φ^2/n, 0)}``)
and `v_rot_floor` (rms of `v_rot` over the sign-flip null draws of
[`rotation_analysis`](@ref)).
"""
struct RotationProfile
    radius::Vector{Float64}
    v_rot::Vector{Float64}
    sigma::Vector{Float64}
    v_rot_over_sigma::Vector{Float64}
    n::Vector{Int}
    sigma_phi::Vector{Float64}
    v_rot_debiased::Vector{Float64}
    v_rot_floor::Vector{Float64}
end

"""Maximum number of cylindrical shells of the rotation profile."""
const _ROTATION_SHELLS = 16
"""Minimum effective membership per shell of the rotation profile."""
const _ROTATION_SHELL_MEMBERS = 1000
"""Sign-flip null draws of the rotation analysis."""
const _ROTATION_NULL_DRAWS = 20

"""Mass-weighted λ_R from the shell values, `w` holding `Σ m R` per shell."""
function _shell_lambda_r(
    w::AbstractVector{Float64},
    v_rot::AbstractVector{Float64},
    sigma::AbstractVector{Float64},
)::Float64
    num = 0.0
    den = 0.0
    @inbounds for b in eachindex(w)
        num += w[b] * abs(v_rot[b])
        den += w[b] * sqrt(v_rot[b]^2 + sigma[b]^2)
    end
    return den > 0 ? num / den : NaN
end

"""
    _rotation_curve_fit(radius, v_rot, n, r_max) -> (v_peak, r_peak)

Weighted least-squares fit of ``V(R) = 2 V_\\mathrm{pk} (R/R_\\mathrm{pk}) /
(1 + (R/R_\\mathrm{pk})^2)`` to the shells with `0 < radius ≤ r_max`, weights
`n`: `R_pk` on a logarithmic grid of 64 values spanning the fitted shell
radii, `V_pk` the linear amplitude at each. `(NaN, NaN)` with fewer than two
shells.
"""
function _rotation_curve_fit(
    radius::AbstractVector{Float64},
    v_rot::AbstractVector{Float64},
    n::AbstractVector{Int},
    r_max::Float64,
)
    fit = [b for b in eachindex(radius) if 0 < radius[b] ≤ r_max]
    length(fit) < 2 && return (NaN, NaN)
    r_lo = minimum(radius[fit])
    r_hi = maximum(radius[fit])
    v_peak = NaN
    r_peak = NaN
    best = Inf
    for r_pk in exp.(range(log(r_lo), log(r_hi); length = 64))
        sfv = 0.0
        sff = 0.0
        for b in fit
            x = radius[b] / r_pk
            f = 2x / (1 + x^2)
            sfv += n[b] * f * v_rot[b]
            sff += n[b] * f * f
        end
        V = sff > 0 ? sfv / sff : 0.0
        res = 0.0
        for b in fit
            x = radius[b] / r_pk
            res += n[b] * (v_rot[b] - V * 2x / (1 + x^2))^2
        end
        if res < best
            best = res
            v_peak = V
            r_peak = r_pk
        end
    end
    return (v_peak, r_peak)
end

"""
    rotation_analysis(pos, vel, mass; n_min = 1000, nbins_max = 16, n_null = 20, seed = 0,
                      energy = nothing, peebles = true, r_h = NaN, shrink = 0.025, n_stop = 1000,
                      stop_fraction = 0.01, n_floor = 500)
        -> (; L, axis, lambda_r, lambda_r_floor, lambda_peebles, lambda_peebles_floor,
             lambda_peebles_floor_analytic, v_peak, r_peak, profile)

Rotation of a bound set about its shrinking-sphere centre (`shrink`,
`n_stop`, `stop_fraction` and `n_floor` go to `shrinking_sphere_centre`).

`L` is the angular momentum of the whole set about its centre and
mass-weighted mean velocity. The spin `axis` is the unit vector of the
angular momentum of the particles within `2 r_h` of the centre, about their
own mean velocity (`[0, 0, 1]` when it vanishes), so that loosely bound
outer stars and tidal debris do not tilt it; `r_h` is the half-mass radius
about the centre unless passed.

The [`RotationProfile`](@ref) uses cylindrical shells of equal membership,
`max(1, min(nbins_max, ⌊N / n_min⌋))` of them; the mass-weighting noise
of a shell is measured by the sign-flip null rather than assumed through
`N_eff = (Σm)² / Σm²`, which enters the analytic λ_P floor only. Each
carries the mass-weighted mean azimuthal
velocity `v_rot`, the dispersions `sigma` and `sigma_phi`, and the debiased
amplitude ``\\sqrt{\\bar v_φ^2 − σ_φ^2/n}``, which removes the noise bias of
a mean of `n` azimuthal velocities. From the shells follow the
ordered-to-total motion parameter
``λ_R = Σ m R |v_\\mathrm{rot}| / Σ m R \\sqrt{v_\\mathrm{rot}^2 + σ^2}``
(Emsellem et al. 2007, doi:10.1111/j.1365-2966.2007.11752.x) and the
Peebles (1969, doi:10.1086/149876) spin parameter
``λ_P = |L| \\sqrt{|E|} / (G M^{5/2})`` with `E` the total energy of the
set (`G = 1`).

Finite-N noise sets a floor to all three. The sign-flip null reverses the
azimuthal velocity of every particle with probability ½ (`n_null` draws
from a `MersenneTwister` seeded with `seed`), leaving the shells, `σ` and
`E` unchanged, and recomputes `v_rot`, `λ_R` and `|L|`: a non-rotating
realisation of the same particles, the N-body analogue of the non-rotating
Monte-Carlo clusters of Sollima, Baumgardt & Hilker (2019,
doi:10.1093/mnras/stz505). `v_rot_floor`, `lambda_r_floor` and
`lambda_peebles_floor` are the rms values over the draws (`NaN` for
`n_null = 0`). `lambda_peebles_floor_analytic = 3π / (32 √N_eff) ≈
0.295 / √N_eff` is the corresponding analytic floor of λ_P for
uncorrelated azimuthal signs.

`v_peak` and `r_peak` fit the rotation curve
``V(R) = 2 V_\\mathrm{pk} (R/R_\\mathrm{pk}) / (1 + (R/R_\\mathrm{pk})^2)``
of Bianchini et al. (2018, doi:10.1093/mnras/sty2225) to the shells within
`3 r_h`, weighted by membership: `R_pk` on a logarithmic grid of 64 values
across the shell radii, `V_pk` the linear least-squares amplitude; `NaN`
with fewer than two shells.

`spin_axis` fixes the axis instead of deriving it from the particles within
`2 r_h`, so that subsets (a radial range, a mass class) are measured about
the axis of the whole set.

`energy` passes `E` when it is known, as for a [`BoundSet`](@ref)
(`kinetic + potential`, regularised pairs at their centres of mass);
otherwise it is the star-level `T + W` of the given particles, which
includes the internal energy of any pair among them and costs an O(n²) pair
sum. `peebles = false` skips λ_P altogether: `lambda_peebles` and
`lambda_peebles_floor` are then `NaN`, for callers that use the shells and
λ_R only.
"""
function rotation_analysis(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    n_min::Int = _ROTATION_SHELL_MEMBERS,
    nbins_max::Int = _ROTATION_SHELLS,
    n_null::Int = _ROTATION_NULL_DRAWS,
    seed::Integer = 0,
    energy::Union{Nothing,Real} = nothing,
    peebles::Bool = true,
    r_h::Real = NaN,
    shrink::Real = 0.025,
    n_stop::Integer = 1000,
    stop_fraction::Real = 0.01,
    n_floor::Integer = 500,
    spin_axis::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    n = length(mass)
    n ≥ _MIN_MEMBERS ||
        throw(ArgumentError("rotation_analysis needs at least $(_MIN_MEMBERS) particles, got $n"))
    if spin_axis !== nothing
        length(spin_axis) == 3 && all(isfinite, spin_axis) && sqrt(sum(abs2, spin_axis)) > 0 ||
            throw(ArgumentError("spin_axis must be a finite non-zero 3-vector"))
    end
    size(pos, 2) == n && size(vel, 2) == n ||
        throw(DimensionMismatch("pos, vel and mass disagree in length"))
    n_min ≥ 2 || throw(ArgumentError("n_min must be ≥ 2, got $n_min"))
    nbins_max ≥ 1 || throw(ArgumentError("nbins_max must be ≥ 1, got $nbins_max"))
    n_null ≥ 0 || throw(ArgumentError("n_null must be ≥ 0, got $n_null"))
    isnan(r_h) || r_h > 0 || throw(ArgumentError("r_h must be positive or NaN, got $r_h"))
    M = sum(mass)
    c = _shrinking_sphere_centre(pos, mass; shrink, n_stop, stop_fraction, n_floor)
    rh = isnan(r_h) ? _lagrangian_radii(pos, mass, c, (0.5,))[1] : Float64(r_h)
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

    # Spin axis from the particles within 2 r_h, about their mean velocity
    r_in² = (2 * rh)^2
    M_in = 0.0
    px = py = pz = 0.0
    @inbounds for i in 1:n
        rx, ry, rz = pos[1, i] - c[1], pos[2, i] - c[2], pos[3, i] - c[3]
        rx * rx + ry * ry + rz * rz ≤ r_in² || continue
        M_in += mass[i]
        px += mass[i] * vel[1, i]
        py += mass[i] * vel[2, i]
        pz += mass[i] * vel[3, i]
    end
    L_in = zeros(3)
    if M_in > 0
        vix, viy, viz = px / M_in, py / M_in, pz / M_in
        @inbounds for i in 1:n
            rx, ry, rz = pos[1, i] - c[1], pos[2, i] - c[2], pos[3, i] - c[3]
            rx * rx + ry * ry + rz * rz ≤ r_in² || continue
            vx, vy, vz = vel[1, i] - vix, vel[2, i] - viy, vel[3, i] - viz
            L_in[1] += mass[i] * (ry * vz - rz * vy)
            L_in[2] += mass[i] * (rz * vx - rx * vz)
            L_in[3] += mass[i] * (rx * vy - ry * vx)
        end
    end
    L_in_norm = sqrt(sum(abs2, L_in))
    axis = if spin_axis !== nothing
        collect(Float64, spin_axis) ./ sqrt(sum(abs2, spin_axis))
    else
        L_in_norm > 0 ? L_in ./ L_in_norm : [0.0, 0.0, 1.0]
    end
    a1, a2, a3 = axis[1], axis[2], axis[3]

    # Cylindrical radius and azimuthal velocity of every particle, and the
    # part m v_φ (R ẑ − z R̂) of its angular momentum that the sign flip reverses
    R = zeros(n)
    vφ = zeros(n)
    Lφ = zeros(3, n)
    @inbounds for i in 1:n
        rx, ry, rz = pos[1, i] - c[1], pos[2, i] - c[2], pos[3, i] - c[3]
        vx, vy, vz = vel[1, i] - vc[1], vel[2, i] - vc[2], vel[3, i] - vc[3]
        z = rx * a1 + ry * a2 + rz * a3
        qx, qy, qz = rx - z * a1, ry - z * a2, rz - z * a3
        Ri = sqrt(qx * qx + qy * qy + qz * qz)
        R[i] = Ri
        Ri > 0 || continue
        ux, uy, uz = qx / Ri, qy / Ri, qz / Ri
        fx = a2 * uz - a3 * uy
        fy = a3 * ux - a1 * uz
        fz = a1 * uy - a2 * ux
        vφi = vx * fx + vy * fy + vz * fz
        vφ[i] = vφi
        w = mass[i] * vφi
        Lφ[1, i] = w * (Ri * a1 - z * ux)
        Lφ[2, i] = w * (Ri * a2 - z * uy)
        Lφ[3, i] = w * (Ri * a3 - z * uz)
    end

    # Shells of equal membership in cylindrical radius
    N_eff = M^2 / sum(abs2, mass)
    nb = max(1, min(nbins_max, length(mass) ÷ n_min))
    order = sortperm(R)
    edges = round.(Int, range(0, n; length = nb + 1))
    radius = zeros(nb)
    v_rot = zeros(nb)
    sigma = zeros(nb)
    sigma_phi = zeros(nb)
    v_rot_debiased = zeros(nb)
    counts = zeros(Int, nb)
    M_shell = zeros(nb)
    w_R = zeros(nb)
    shell_of = zeros(Int, n)
    @inbounds for b in 1:nb
        lo, hi = edges[b] + 1, edges[b + 1]
        Mb = 0.0
        sR = 0.0
        sφ = 0.0
        sx = sy = sz = 0.0
        for q in lo:hi
            i = order[q]
            m = mass[i]
            shell_of[i] = b
            Mb += m
            sR += m * R[i]
            sφ += m * vφ[i]
            sx += m * (vel[1, i] - vc[1])
            sy += m * (vel[2, i] - vc[2])
            sz += m * (vel[3, i] - vc[3])
        end
        counts[b] = hi - lo + 1
        M_shell[b] = Mb
        w_R[b] = sR
        radius[b] = sR / Mb
        v_rot[b] = sφ / Mb
        mx, my, mz = sx / Mb, sy / Mb, sz / Mb
        s2 = 0.0
        sp2 = 0.0
        for q in lo:hi
            i = order[q]
            m = mass[i]
            s2 +=
                m * (
                    (vel[1, i] - vc[1] - mx)^2 +
                    (vel[2, i] - vc[2] - my)^2 +
                    (vel[3, i] - vc[3] - mz)^2
                )
            sp2 += m * (vφ[i] - v_rot[b])^2
        end
        sigma[b] = sqrt(s2 / (3 * Mb))
        sigma_phi[b] = sqrt(sp2 / Mb)
        v_rot_debiased[b] = sign(v_rot[b]) * sqrt(max(v_rot[b]^2 - sigma_phi[b]^2 / counts[b], 0.0))
    end
    ratio = [sigma[b] > 0 ? v_rot[b] / sigma[b] : NaN for b in 1:nb]
    lambda_r = _shell_lambda_r(w_R, v_rot, sigma)

    E = if !peebles
        NaN
    elseif energy === nothing
        T = 0.0
        for i in 1:n
            T += 0.5 * mass[i] * sum((vel[k, i] - vc[k])^2 for k in 1:3)
        end
        T + 0.5 * sum(mass .* pair_potential(pos, mass))
    else
        Float64(energy)
    end
    peebles = sqrt(abs(E)) / M^2.5
    lambda_peebles = Lnorm * peebles
    lambda_peebles_floor_analytic = 0.295 / sqrt(N_eff)

    # Sign-flip null: shells, σ and E fixed, only the signs of v_φ change
    v_rot_floor = fill(NaN, nb)
    lambda_r_floor = NaN
    lambda_peebles_floor = NaN
    if n_null > 0
        rng = Random.MersenneTwister(seed)
        Lfix = L .- vec(sum(Lφ; dims = 2))
        acc = zeros(nb)
        v_null = zeros(nb)
        sv2 = zeros(nb)
        sλr2 = 0.0
        sλp2 = 0.0
        for _ in 1:n_null
            fill!(acc, 0.0)
            lx, ly, lz = Lfix[1], Lfix[2], Lfix[3]
            @inbounds for i in 1:n
                s = rand(rng, Bool) ? -1.0 : 1.0
                acc[shell_of[i]] += s * mass[i] * vφ[i]
                lx += s * Lφ[1, i]
                ly += s * Lφ[2, i]
                lz += s * Lφ[3, i]
            end
            for b in 1:nb
                v_null[b] = acc[b] / M_shell[b]
                sv2[b] += v_null[b]^2
            end
            sλr2 += _shell_lambda_r(w_R, v_null, sigma)^2
            sλp2 += (sqrt(lx * lx + ly * ly + lz * lz) * peebles)^2
        end
        v_rot_floor .= sqrt.(sv2 ./ n_null)
        lambda_r_floor = sqrt(sλr2 / n_null)
        lambda_peebles_floor = sqrt(sλp2 / n_null)
    end

    v_peak, r_peak = _rotation_curve_fit(radius, v_rot, counts, 3 * rh)
    profile =
        RotationProfile(radius, v_rot, sigma, ratio, counts, sigma_phi, v_rot_debiased, v_rot_floor)

    return (;
        L,
        axis,
        lambda_r,
        lambda_r_floor,
        lambda_peebles,
        lambda_peebles_floor,
        lambda_peebles_floor_analytic,
        v_peak,
        r_peak,
        profile,
    )
end

# ---------------------------------------------------------------------------
# Mass segregation (Allison et al. 2009; Maschberger & Clarke 2011;
# Olczak, Spurzem & Henning 2011)
# ---------------------------------------------------------------------------

"""Random comparison sets per N_MST above `n_small`."""
const _MSR_N_RANDOM = 200
"""Random comparison sets per N_MST up to `n_small`, and that limit."""
const _MSR_N_RANDOM_SMALL = 1000
const _MSR_N_SMALL = 20
"""Default scan of minimum-spanning-tree sizes N_MST."""
const _MSR_N_MST = (10, 20, 50, 100, 200)
"""N_MST of the headline segregation values when it is part of the scan."""
const _MSR_N_HEADLINE = 20

"""Mean and sample standard deviation of `v` (`NaN` deviation for fewer than two entries)."""
function _mean_std(v::AbstractVector{Float64})
    n = length(v)
    μ = sum(v) / n
    n < 2 && return (μ, NaN)
    s = 0.0
    @inbounds for x in v
        s += (x - μ)^2
    end
    return (μ, sqrt(s / (n - 1)))
end

"""
    _mst_edges(pos, idx) -> Vector{Float64}

Edge lengths of the Euclidean minimum spanning tree of the points
`pos[:, idx]` in the order Prim's algorithm adds them (`n − 1` edges,
O(n²)).
"""
function _mst_edges(pos::AbstractMatrix{Float64}, idx::AbstractVector{Int})::Vector{Float64}
    n = length(idx)
    edges = Vector{Float64}(undef, max(n - 1, 0))
    n ≤ 1 && return edges
    in_tree = falses(n)
    best = fill(Inf, n)
    in_tree[1] = true
    p = idx[1]
    @inbounds for j in 2:n
        q = idx[j]
        best[j] =
            sqrt((pos[1, p] - pos[1, q])^2 + (pos[2, p] - pos[2, q])^2 + (pos[3, p] - pos[3, q])^2)
    end
    @inbounds for e in 1:(n - 1)
        j_min = 0
        d_min = Inf
        for j in 1:n
            if !in_tree[j] && best[j] < d_min
                d_min = best[j]
                j_min = j
            end
        end
        in_tree[j_min] = true
        edges[e] = d_min
        p = idx[j_min]
        for j in 1:n
            in_tree[j] && continue
            q = idx[j]
            d = sqrt(
                (pos[1, p] - pos[1, q])^2 + (pos[2, p] - pos[2, q])^2 + (pos[3, p] - pos[3, q])^2,
            )
            d < best[j] && (best[j] = d)
        end
    end
    return edges
end

"""Length of the Euclidean minimum spanning tree of the points `pos[:, idx]` (Prim, O(n²))."""
_mst_length(pos::AbstractMatrix{Float64}, idx::AbstractVector{Int}) = sum(_mst_edges(pos, idx))

"""
    mass_segregation(pos, mass; tracer = nothing, n_mst = (10, 20, 50, 100, 200),
                     n_random = 200, n_random_small = 1000, n_small = 20, seed = 0)
        -> (; lambda_msr, lambda_err, lambda_median, lambda_median_err, gamma_msr, gamma_err,
             r_half_ratio, n_mst_used::Vector{Int}, lambda_scan::Vector{Float64},
             lambda_err_scan::Vector{Float64}, lambda_median_scan::Vector{Float64},
             lambda_median_err_scan::Vector{Float64}, gamma_scan::Vector{Float64},
             gamma_err_scan::Vector{Float64})

Mass-segregation ratios of a tracer class against random sets of the same
size, from the edges of their Euclidean minimum spanning trees (MST).

- `Λ_MSR = ⟨Σ l_random⟩ / Σ l_tracer` with the error `σ(Σ l_random) / Σ l_tracer`
  (Allison et al. 2009, doi:10.1111/j.1365-2966.2009.14508.x);
- `Λ̃_MSR` (`lambda_median`), the same ratio with the median edge length of
  each tree, which is insensitive to single long edges to outlying tracer
  stars (Maschberger & Clarke 2011, doi:10.1111/j.1365-2966.2010.17602.x);
- `Γ_MSR = exp⟨⟨ln l⟩_random⟩ / exp⟨ln l_tracer⟩`, the ratio of geometric
  mean edge lengths, with the multiplicative interval factor
  `exp σ(⟨ln l⟩_random)` as `gamma_err` (Olczak, Spurzem & Henning 2011,
  doi:10.1051/0004-6361/201015264).

All three are ≈ 1 without segregation and > 1 with the tracers concentrated.

`tracer` holds the indices (columns of `pos`) of the tracer class in the
order of decreasing initial mass, chosen once by the caller from a fixed set
of identities. The `N` most massive stars by present mass change identity
through stellar evolution: the massive stars become remnants and drop out of
the set, so a present-mass selection follows whichever stars currently lead
the mass function instead of one population. `nothing` selects the
`maximum(n_mst)` most massive stars by present mass.

The ratios are evaluated for every `N_MST` in `n_mst` with
`2 ≤ N_MST ≤ length(tracer)` (`n_mst_used`, ascending), the tracer subset
being the first `N_MST` entries of `tracer`. Each uses `n_random` random sets
drawn from all particles, `n_random_small` when `N_MST ≤ n_small`: small
trees fluctuate strongly and need ≥ 1000 sets, larger ones ≥ 200. The
headline values are those at `N_MST = 20` when it is part of the scan,
otherwise at the largest `N_MST`; Olczak et al. find the clearest signature
for 10–20 tracers. Without segregation `Λ_MSR` exceeds unity by more than
its error in about one case in six, so a single 1σ excess is no detection.

`r_half_ratio` is the half-mass radius of the tracer set over that of the
whole set, both about the shrinking-sphere centre. The random sets come
from a `MersenneTwister` seeded with `seed`, so the estimates are
reproducible. Throws `ArgumentError` when no `N_MST` fits the tracer class
or the particles do not exceed the largest `N_MST` used.
"""
function mass_segregation(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    tracer::Union{Nothing,AbstractVector{<:Integer}} = nothing,
    n_mst::Union{Tuple{Vararg{Integer}},AbstractVector{<:Integer}} = _MSR_N_MST,
    n_random::Int = _MSR_N_RANDOM,
    n_random_small::Int = _MSR_N_RANDOM_SMALL,
    n_small::Int = _MSR_N_SMALL,
    seed::Integer = 0,
)
    n = size(pos, 2)
    length(mass) == n || throw(DimensionMismatch("mass and pos disagree in length"))
    isempty(n_mst) && throw(ArgumentError("n_mst must not be empty"))
    n_random ≥ 2 || throw(ArgumentError("n_random must be ≥ 2, got $n_random"))
    n_random_small ≥ 2 || throw(ArgumentError("n_random_small must be ≥ 2, got $n_random_small"))
    tr = if tracer === nothing
        collect(Int, partialsortperm(mass, 1:min(Int(maximum(n_mst)), n); rev = true))
    else
        collect(Int, tracer)
    end
    all(i -> 1 ≤ i ≤ n, tr) || throw(ArgumentError("tracer indices must lie in 1:$n"))
    allunique(tr) || throw(ArgumentError("tracer indices must be unique"))
    n_used = sort!(unique!(Int[m for m in n_mst if 2 ≤ m ≤ length(tr)]))
    isempty(n_used) && throw(
        ArgumentError("no N_MST of $(Tuple(n_mst)) lies in 2:$(length(tr)) (tracer class size)"),
    )
    n > last(n_used) ||
        throw(ArgumentError("mass_segregation needs more than $(last(n_used)) particles, got $n"))

    n_s = length(n_used)
    Λ = fill(NaN, n_s)
    Λ_err = fill(NaN, n_s)
    Λ̃ = fill(NaN, n_s)
    Λ̃_err = fill(NaN, n_s)
    Γ = fill(NaN, n_s)
    Γ_err = fill(NaN, n_s)
    rng = Random.MersenneTwister(seed)
    perm = Vector{Int}(undef, n)
    for (s, m) in enumerate(n_used)
        e_t = _mst_edges(pos, view(tr, 1:m))
        l_t = sum(e_t)
        med_t = _median(e_t)
        g_t = sum(log, e_t) / length(e_t)
        k = m ≤ n_small ? n_random_small : n_random
        l_r = Vector{Float64}(undef, k)
        med_r = Vector{Float64}(undef, k)
        g_r = Vector{Float64}(undef, k)
        for d in 1:k
            Random.randperm!(rng, perm)
            e = _mst_edges(pos, view(perm, 1:m))
            l_r[d] = sum(e)
            med_r[d] = _median(e)
            g_r[d] = sum(log, e) / length(e)
        end
        μ_l, σ_l = _mean_std(l_r)
        μ_med, σ_med = _mean_std(med_r)
        μ_g, σ_g = _mean_std(g_r)
        if l_t > 0
            Λ[s] = μ_l / l_t
            Λ_err[s] = σ_l / l_t
        end
        if med_t > 0
            Λ̃[s] = μ_med / med_t
            Λ̃_err[s] = σ_med / med_t
        end
        isfinite(g_t) && (Γ[s] = exp(μ_g - g_t))
        Γ_err[s] = exp(σ_g)
    end
    h = something(findfirst(==(_MSR_N_HEADLINE), n_used), n_s)

    c = _shrinking_sphere_centre(pos, mass)
    r_half_all = _lagrangian_radii(pos, mass, c, (0.5,))[1]
    r_half_tracer = _lagrangian_radii(pos[:, tr], mass[tr], c, (0.5,))[1]
    return (;
        lambda_msr = Λ[h],
        lambda_err = Λ_err[h],
        lambda_median = Λ̃[h],
        lambda_median_err = Λ̃_err[h],
        gamma_msr = Γ[h],
        gamma_err = Γ_err[h],
        r_half_ratio = r_half_tracer / r_half_all,
        n_mst_used = n_used,
        lambda_scan = Λ,
        lambda_err_scan = Λ_err,
        lambda_median_scan = Λ̃,
        lambda_median_err_scan = Λ̃_err,
        gamma_scan = Γ,
        gamma_err_scan = Γ_err,
    )
end

# ---------------------------------------------------------------------------
# Core collapse (Tanikawa, Hut & Makino 2012)
# ---------------------------------------------------------------------------

"""
Running mean of `v` over `window` entries centred on each one (truncated at
the ends), averaging the non-`NaN` entries; `NaN` where `v` itself is `NaN`.
"""
function _running_mean(v::AbstractVector{Float64}, window::Int)::Vector{Float64}
    n = length(v)
    h = window ÷ 2
    out = fill(NaN, n)
    @inbounds for i in 1:n
        isnan(v[i]) && continue
        s = 0.0
        c = 0
        for j in max(1, i - h):min(n, i + h)
            isnan(v[j]) && continue
            s += v[j]
            c += 1
        end
        out[i] = s / c
    end
    return out
end

"""
    core_collapse_markers(time, r_core, rho_core; window = 5)
        -> (; i_rc_min, i_rho_max, t_rc_min, t_rho_max)

Structural markers of core collapse in a snapshot series: the index and
time of the minimum of the core radius `r_core` (`i_rc_min`, `t_rc_min`)
and of the maximum of the core density `rho_core` (`i_rho_max`,
`t_rho_max`), both after a running mean over `window` snapshots (odd,
`≥ 1`, truncated at the ends). `NaN` entries are ignored; an index is
`nothing` and its time `NaN` when every entry is `NaN`.

These are two of the set of collapse markers compared by Tanikawa, Hut &
Makino (2012, doi:10.1016/j.newast.2011.09.001) and Fujii & Portegies Zwart
(2014); the markers based on the first hard binaries come from the snapshot
pair finder. The markers generally fall at different snapshots, so core
collapse is reported as this set, never as one moment.
"""
function core_collapse_markers(
    time::AbstractVector{<:Real},
    r_core::AbstractVector{<:Real},
    rho_core::AbstractVector{<:Real};
    window::Integer = 5,
)
    n = length(time)
    length(r_core) == n && length(rho_core) == n ||
        throw(DimensionMismatch("time, r_core and rho_core disagree in length"))
    window ≥ 1 && isodd(window) ||
        throw(ArgumentError("window must be an odd integer ≥ 1, got $window"))
    t = collect(Float64, time)
    rc = _running_mean(collect(Float64, r_core), Int(window))
    ρc = _running_mean(collect(Float64, rho_core), Int(window))
    i_rc_min = nothing
    i_rho_max = nothing
    for i in 1:n
        if !isnan(rc[i]) && (i_rc_min === nothing || rc[i] < rc[i_rc_min])
            i_rc_min = i
        end
        if !isnan(ρc[i]) && (i_rho_max === nothing || ρc[i] > ρc[i_rho_max])
            i_rho_max = i
        end
    end
    t_rc_min = i_rc_min === nothing ? NaN : t[i_rc_min]
    t_rho_max = i_rho_max === nothing ? NaN : t[i_rho_max]
    return (; i_rc_min, i_rho_max, t_rc_min, t_rho_max)
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
- `lambda_msr`, `lambda_msr_err`, `segregation_ratio`: mass segregation of
  the tracer class, headline values ([`mass_segregation`](@ref))
- `coalescence_time`, `coalescence_time_myr`: first snapshot from which all
  initial clusters form one group ([`coalescence_time`](@ref)); `NaN` if never
- `segregation_time`, `segregation_time_myr`: first snapshot with
  `Λ_MSR − σ > msr_threshold`; `NaN` if never
- `profile`: rotation profile of the last snapshot with enough members
- `r_core_number`: Casertano–Hut core radius from the number density
- `rho_core`: maximum local mass density among the bound members
- `lambda_r_floor`, `lambda_peebles_floor`: sign-flip null floors of λ_R and λ_P
- `lambda_median`, `gamma_msr`, `gamma_msr_err`: median-edge ratio `Λ̃_MSR`
  and geometric-mean ratio `Γ_MSR` with its interval factor
- `v_peak`, `r_peak`: rotation-curve fit ([`rotation_analysis`](@ref))
- `tracer_names`: names of the tracer class of the mass segregation, in
  decreasing mass at the first snapshot
- `n_pairs_snapshot`: pairs of the snapshot pair finder among the bound
  stars that pass its tidal criterion (`tidal_ok`)
- `n_hard_pairs`: those of them with hardness `E_b / kT > 1`
- `hardness_max`: largest hardness among them; `NaN` when there is none
- `kt`: hard/soft energy scale `kT = 2K / (3n)` of the bound set, `K` its
  kinetic energy about its mean velocity and `n` its star count
- `f_max`: largest stellar mass of the bound set over its mean stellar mass
- `t_first_10kt`, `t_first_10kt_myr`: first snapshot with
  `hardness_max > 10`; `NaN` if never
- `t_first_10fmax_kt`, `t_first_10fmax_kt_myr`: first snapshot with
  `hardness_max > 10 f_max`; `NaN` if never
- `t_rc_min`, `t_rc_min_myr`: minimum of the smoothed core radius `r_core`
- `t_rho_max`, `t_rho_max_myr`: maximum of the smoothed core density
  `rho_core`
- `coalescence_pairs`, `coalescence_separation`, `coalescence_threshold`: the
  cluster pairs and, per pair (row) and snapshot (column), the separation of
  the per-cluster density centres and the coalescence threshold of
  [`coalescence_time`](@ref) (`NaN` where a cluster is absent); empty with
  fewer than two clusters

The last four pairs of times are the core-collapse markers: the first hard
binaries (Tanikawa, Hut & Makino 2012, doi:10.1016/j.newast.2011.09.001;
Fujii & Portegies Zwart 2014) and the structural
extrema of `core_collapse_markers` (running mean over five snapshots).
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
    r_core_number::Vector{Float64}
    rho_core::Vector{Float64}
    lambda_r_floor::Vector{Float64}
    lambda_peebles_floor::Vector{Float64}
    lambda_median::Vector{Float64}
    gamma_msr::Vector{Float64}
    gamma_msr_err::Vector{Float64}
    v_peak::Vector{Float64}
    r_peak::Vector{Float64}
    tracer_names::Vector{Int32}
    n_pairs_snapshot::Vector{Int}
    n_hard_pairs::Vector{Int}
    hardness_max::Vector{Float64}
    kt::Vector{Float64}
    f_max::Vector{Float64}
    t_first_10kt::Float64
    t_first_10kt_myr::Float64
    t_first_10fmax_kt::Float64
    t_first_10fmax_kt_myr::Float64
    t_rc_min::Float64
    t_rc_min_myr::Float64
    t_rho_max::Float64
    t_rho_max_myr::Float64
    coalescence_pairs::Vector{Tuple{Int,Int}}
    coalescence_separation::Matrix{Float64}
    coalescence_threshold::Matrix{Float64}
end

"""Λ_MSR level (minus its error) above which the remnant counts as mass-segregated."""
const _MSR_THRESHOLD = 2.0

"""
    remnant_diagnostics(snaps, cluster_ranges; settings = PostprocessConfig(),
                        bound_sets = nothing, max_passes = 100,
                        min_members = 10) -> RemnantDiagnostics

Bound remnant of the whole system at every snapshot (the bound set of all
stars, regularised pairs at their centres of mass, [`bound_set`](@ref) with
`max_passes` and `min_members`), its core radii from the mass and the number
density, core density, half-mass radius, rotation, mass segregation and
snapshot pairs, plus the coalescence, segregation and core-collapse times.
`bound_sets` passes the bound sets when they are already known, one per
snapshot (the `bound` of [`snapshot_integrals`](@ref)); otherwise each costs
one O(N²) pair sum.

Every estimator parameter comes from `settings`, a
[`PostprocessConfig`](@ref). The local densities use the neighbour rank
`density_neighbours`; the half-mass radius and the rotation are measured
about the shrinking-sphere centre with the centre keys. The mass segregation
follows a fixed tracer class: the names of the `msr_n_tracer` most massive
stars of the first snapshot's bound set, in decreasing mass. At every
snapshot the tracers still in the bound set enter [`mass_segregation`](@ref)
in that order, with the N_MST scan `msr_scan`, `msr_n_random` random sets per
N_MST above `msr_n_small` (`msr_n_random_small` up to it) and the seed
`msr_seed`; the segregation time is the first snapshot with
`Λ_MSR − σ > msr_threshold`. [`rotation_analysis`](@ref) takes
`rotation_n_min`, `rotation_nbins_max`, `rotation_n_null` and `msr_seed`.

The snapshot pair finder (`snapshot_pairs` with `pair_k` and
`pair_tidal_factor`) runs on the bound stars of every snapshot with the
energy scale `kT = 2K / (3n)`, `K` the kinetic energy of the bound set about
its mean velocity and `n` its star count; the hardness series count the
pairs that pass the tidal criterion. The first snapshots with a pair harder
than `10 kT` and than `10 f_max kT` (`f_max` the largest over the mean
stellar mass) and the extrema of the core radius and core density
(`core_collapse_markers`, `window = 5`) are the core-collapse markers of
[`RemnantDiagnostics`](@ref). The pair finder costs O(N log N) per snapshot.
"""
function remnant_diagnostics(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    settings::PostprocessConfig = PostprocessConfig(),
    bound_sets::Union{Nothing,AbstractVector{BoundSet}} = nothing,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)::RemnantDiagnostics
    isempty(snaps) && throw(ArgumentError("remnant_diagnostics needs at least one snapshot"))
    settings.msr_n_tracer ≥ 2 ||
        throw(ArgumentError("settings.msr_n_tracer must be ≥ 2, got $(settings.msr_n_tracer)"))
    bound_sets === nothing ||
        length(bound_sets) == length(snaps) ||
        throw(DimensionMismatch("$(length(bound_sets)) bound sets for $(length(snaps)) snapshots"))
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
    r_core_n = fill(NaN, n_t)
    ρ_core = fill(NaN, n_t)
    λr_floor = fill(NaN, n_t)
    λp_floor = fill(NaN, n_t)
    λmed = fill(NaN, n_t)
    γ = fill(NaN, n_t)
    γ_err = fill(NaN, n_t)
    v_pk = fill(NaN, n_t)
    r_pk = fill(NaN, n_t)
    n_pairs = zeros(Int, n_t)
    n_hard = zeros(Int, n_t)
    x_max = fill(NaN, n_t)
    kt = fill(NaN, n_t)
    f_max = fill(NaN, n_t)
    centre_kw = _centre_keywords(settings)
    tracer_names = Int32[]
    profile = RotationProfile(
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Int[],
        Float64[],
        Float64[],
        Float64[],
    )

    L_orb = _orbital_angular_momentum(snaps[1], cluster_ranges)
    L_norm = sqrt(sum(abs2, L_orb))
    L_hat = L_norm > 0 ? L_orb ./ L_norm : fill(NaN, 3)
    bsets = Vector{BoundSet}(undef, n_t)

    for (k, snap) in enumerate(snaps)
        pos_all = Float64.(snap.pos)
        vel_all = Float64.(snap.vel)
        mass_all = Float64.(snap.mass)
        length(mass_all) < _MIN_MEMBERS && continue
        b = bound_sets === nothing ? bound_set(snap; max_passes, min_members) : bound_sets[k]
        bsets[k] = b
        sel = b.stars
        n_bound[k] = length(sel)
        fbound[k] = sum(mass_all[sel]) / sum(mass_all)
        if k == 1
            # Tracer class: fixed identities, in decreasing mass at the first snapshot
            n_tr = min(settings.msr_n_tracer, length(sel))
            top = n_tr > 0 ? partialsortperm(mass_all[sel], 1:n_tr; rev = true) : Int[]
            tracer_names = snap.name[sel[top]]
        end
        length(sel) < _MIN_MEMBERS && continue
        pos = pos_all[:, sel]
        vel = vel_all[:, sel]
        mass = mass_all[sel]

        # Snapshot pairs on the bound stars, hardness against kT = 2K/(3n)
        n_sel = length(sel)
        kt[k] = b.kinetic > 0 ? 2 * b.kinetic / (3 * n_sel) : NaN
        f_max[k] = maximum(mass) / (sum(mass) / n_sel)
        found = snapshot_pairs(
            pos,
            vel,
            mass;
            k = settings.pair_k,
            tidal_factor = settings.pair_tidal_factor,
            kT = kt[k],
        )
        for p in found
            p.tidal_ok || continue
            n_pairs[k] += 1
            p.hardness > 1 && (n_hard[k] += 1)
            if !isnan(p.hardness) && (isnan(x_max[k]) || p.hardness > x_max[k])
                x_max[k] = p.hardness
            end
        end

        k_nb = settings.density_neighbours
        if n_sel > k_nb
            ρ = _snapshot_density(snap, sel; k = k_nb)
            cr = core_radius(pos, ρ)
            r_core[k] = cr.r_core
            ρ_core[k] = maximum(ρ)
            ρn = _snapshot_density(snap, sel; weighting = :number, k = k_nb)
            r_core_n[k] = core_radius(pos, ρn).r_core
        end
        c = _shrinking_sphere_centre(pos, mass; centre_kw...)
        r_half[k] = _lagrangian_radii(pos, mass, c, (0.5,))[1]

        rot = rotation_analysis(
            pos,
            vel,
            mass;
            n_min = settings.rotation_n_min,
            nbins_max = settings.rotation_nbins_max,
            n_null = settings.rotation_n_null,
            seed = settings.msr_seed,
            energy = b.kinetic + b.potential,
            r_h = r_half[k],
            centre_kw...,
        )
        λr[k] = rot.lambda_r
        λp[k] = rot.lambda_peebles
        λr_floor[k] = rot.lambda_r_floor
        λp_floor[k] = rot.lambda_peebles_floor
        v_pk[k] = rot.v_peak
        r_pk[k] = rot.r_peak
        axis[:, k] = rot.axis
        align[k] = all(isfinite, L_hat) ? dot(rot.axis, L_hat) : NaN
        profile = rot.profile

        # Positions of the tracer names among the bound stars, order kept
        slot = Dict{Int32,Int}(snap.name[sel[j]] => j for j in eachindex(sel))
        tracer = Int[]
        for nm in tracer_names
            j = get(slot, nm, 0)
            j > 0 && push!(tracer, j)
        end
        n_scan = maximum((m for m in settings.msr_scan if 2 ≤ m ≤ length(tracer)); init = 0)
        if n_scan ≥ 2 && length(sel) > n_scan
            ms = mass_segregation(
                pos,
                mass;
                tracer,
                n_mst = Tuple(settings.msr_scan),
                n_random = settings.msr_n_random,
                n_random_small = settings.msr_n_random_small,
                n_small = settings.msr_n_small,
                seed = settings.msr_seed,
            )
            λmsr[k] = ms.lambda_msr
            λerr[k] = ms.lambda_err
            seg[k] = ms.r_half_ratio
            λmed[k] = ms.lambda_median
            γ[k] = ms.gamma_msr
            γ_err[k] = ms.gamma_err
        end
    end

    coal =
        length(cluster_ranges) ≥ 2 ?
        coalescence_time(snaps, cluster_ranges; bound_sets = bsets, max_passes, min_members) :
        (;
            index = nothing,
            time_nb = NaN,
            time_myr = NaN,
            pairs = Tuple{Int,Int}[],
            separation = zeros(Float64, 0, n_t),
            threshold = zeros(Float64, 0, n_t),
        )
    k_seg = findfirst(k -> isfinite(λmsr[k]) && λmsr[k] - λerr[k] > settings.msr_threshold, 1:n_t)
    t_seg = k_seg === nothing ? NaN : t[k_seg]
    t_seg_myr = k_seg === nothing ? NaN : t_myr[k_seg]

    # Core-collapse markers: first hard pairs, then the structural extrema
    k_10 = findfirst(k -> x_max[k] > 10, 1:n_t)
    k_10f = findfirst(k -> x_max[k] > 10 * f_max[k], 1:n_t)
    cc = core_collapse_markers(t, r_core, ρ_core; window = 5)
    t_at(i, times) = i === nothing ? NaN : times[i]

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
        r_core_n,
        ρ_core,
        λr_floor,
        λp_floor,
        λmed,
        γ,
        γ_err,
        v_pk,
        r_pk,
        tracer_names,
        n_pairs,
        n_hard,
        x_max,
        kt,
        f_max,
        t_at(k_10, t),
        t_at(k_10, t_myr),
        t_at(k_10f, t),
        t_at(k_10f, t_myr),
        cc.t_rc_min,
        t_at(cc.i_rc_min, t_myr),
        cc.t_rho_max,
        t_at(cc.i_rho_max, t_myr),
        coal.pairs,
        coal.separation,
        coal.threshold,
    )
end

"""
    write_remnant_diagnostics(path, diag::RemnantDiagnostics) -> String

Write the time series of `diag` as CSV (one row per snapshot, NB units
plus `time_myr`) with the coalescence, segregation and core-collapse times
in a leading comment line. An existing file is backed up, never overwritten.
"""
function write_remnant_diagnostics(path::AbstractString, diag::RemnantDiagnostics)
    _backup_existing(path)
    open(path, "w") do io
        println(
            io,
            "# coalescence_time_nb=$(diag.coalescence_time) coalescence_time_myr=$(diag.coalescence_time_myr) segregation_time_nb=$(diag.segregation_time) segregation_time_myr=$(diag.segregation_time_myr)" *
            " t_first_10kt_nb=$(diag.t_first_10kt) t_first_10kt_myr=$(diag.t_first_10kt_myr) t_first_10fmax_kt_nb=$(diag.t_first_10fmax_kt) t_first_10fmax_kt_myr=$(diag.t_first_10fmax_kt_myr)" *
            " t_rc_min_nb=$(diag.t_rc_min) t_rc_min_myr=$(diag.t_rc_min_myr) t_rho_max_nb=$(diag.t_rho_max) t_rho_max_myr=$(diag.t_rho_max_myr)",
        )
        println(
            io,
            "time_nb,time_myr,rbar_pc,n_bound,bound_mass_fraction,r_core,r_half,lambda_r,lambda_peebles,spin_x,spin_y,spin_z,spin_alignment,lambda_msr,lambda_msr_err,segregation_ratio,r_core_number,rho_core,lambda_r_floor,lambda_peebles_floor,lambda_median,gamma_msr,gamma_msr_err,v_peak,r_peak,n_pairs_snapshot,n_hard_pairs,hardness_max,kt,f_max",
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
                diag.r_core_number[k],
                diag.rho_core[k],
                diag.lambda_r_floor[k],
                diag.lambda_peebles_floor[k],
                diag.lambda_median[k],
                diag.gamma_msr[k],
                diag.gamma_msr_err[k],
                diag.v_peak[k],
                diag.r_peak[k],
                diag.n_pairs_snapshot[k],
                diag.n_hard_pairs[k],
                diag.hardness_max[k],
                diag.kt[k],
                diag.f_max[k],
            ]
            println(io, join(_csv_cell.(cells), ","))
        end
    end
    return String(path)
end
