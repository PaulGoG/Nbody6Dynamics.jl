# =============================================================================
# Static analytic test models of the estimator calibration
# =============================================================================
# Osipkov–Merritt anisotropic Plummer sphere (Merritt 1985,
# DOI 10.1086/113810), Lynden-Bell rotation by sign flips (Lynden-Bell 1960,
# DOI 10.1093/mnras/120.3.204), stretched positions for the shape estimator
# (Zemp et al. 2011, DOI 10.1088/0067-0049/197/2/30) and an energy-ranked mass
# segregation of tunable degree (a calibration model of its own). Each model has a known answer for one
# estimator; all of them are also usable as initial conditions. Units as in
# `sample_plummer`: G = 1, total mass 1, lengths in units of the scale radius.

"""
    _om_plummer_df(q, x) -> Float64

Osipkov–Merritt Plummer distribution function, up to a constant, in units
of `GM/a`: `g(q) = (128/7)(1 − x) q^{7/2} + 8x q^{3/2}` with `x = a²/r_a²`.
Requires `q ≥ 0`.
"""
@inline function _om_plummer_df(q::Float64, x::Float64)::Float64
    return (128 / 7) * (1.0 - x) * q^3.5 + 8.0 * x * q^1.5
end

"""
    sample_plummer_om(N::Int, a::Float64, r_a::Float64; rng = Random.default_rng()) -> (pos, vel)

Osipkov–Merritt anisotropic Plummer model of `N` equal-mass particles with
scale radius `a` and anisotropy radius `r_a` (Merritt 1985,
DOI 10.1086/113810). The density is that of the Plummer sphere,
`ρ = 3/(4π a³) (1 + r²/a²)^{-5/2}`; the velocity anisotropy is
`β(r) = r²/(r² + r_a²)`, isotropic inside `r_a` and radially biased outside.
`r_a` is in the same length units as `a`; `r_a = Inf` gives the isotropic
Plummer model.

The distribution function follows from Eddington inversion in the variable
`Q = E + L²/(2r_a²)`: the augmented density `ρ_Q = (1 + r²/r_a²) ρ` is a
polynomial in the relative potential `ψ = (1 + r²/a²)^{-1/2}`, so with
`x = a²/r_a²` the inversion gives the closed form
`f(Q) ∝ (128/7)(1 − x) Q^{7/2} + 8x Q^{3/2}`, which for `x = 0` reduces to
the isotropic Plummer distribution function `∝ (−E)^{7/2}`. It is
non-negative for all bound orbits only when `x ≤ 16/9`, i.e. `r_a ≥ 0.75 a`;
smaller `r_a` is rejected.

Positions are drawn by the inverse CDF as in [`sample_plummer`](@ref);
speed and pitch angle by von Neumann rejection from the distribution
function at the particle's radius, with a uniformly random azimuth of the
tangential velocity.

Returns `(pos::Matrix{Float64}, vel::Matrix{Float64})`, each `3 × N`, in
units where `G = M_total = 1`.

# Example
```julia
pos, vel = sample_plummer_om(10_000, 1.0, 1.5; rng = Xoshiro(1))
```
"""
function sample_plummer_om(
    N::Int,
    a::Float64,
    r_a::Float64;
    rng::AbstractRNG = Random.default_rng(),
)
    N ≥ 1 || throw(ArgumentError("N must be at least 1, got $N"))
    a > 0 || throw(ArgumentError("a must be positive, got $a"))
    r_a ≥ 0.75 * a || throw(
        ArgumentError(
            "r_a must be at least 0.75 a: the Osipkov–Merritt Plummer distribution function " *
            "is negative below it (x = a²/r_a² > 16/9); got r_a = $(r_a / a) a",
        ),
    )

    x = a^2 / r_a^2
    inv_sqrt_a = 1 / sqrt(a)
    pos = zeros(Float64, 3, N)
    vel = zeros(Float64, 3, N)

    for i in 1:N
        # Inverse CDF for radius: r = a / sqrt(u^{-2/3} - 1)
        u = rand(rng)
        r = a / sqrt(u^(-2/3) - 1.0)

        # Uniform direction on sphere
        cosθ = 2.0 * rand(rng) - 1.0
        sinθ = sqrt(1.0 - cosθ^2)
        φ = 2π * rand(rng)
        r̂x = sinθ * cos(φ)
        r̂y = sinθ * sin(φ)
        r̂z = cosθ

        pos[1, i] = r * r̂x
        pos[2, i] = r * r̂y
        pos[3, i] = r * r̂z

        # Relative potential and escape speed in units of GM/a and sqrt(GM/a)
        ψ = 1 / sqrt(1.0 + (r / a)^2)
        ṽ_esc² = 2.0 * ψ
        ṽ_esc = sqrt(ṽ_esc²)
        ra_term = x == 0.0 ? 0.0 : (r / r_a)^2

        # Envelope: g is increasing on [0, ψ] for x ≤ 1, otherwise it may
        # peak at the interior stationary point q*.
        g_max = _om_plummer_df(ψ, x)
        if x > 1.0
            q_star = sqrt(12.0 * x / (64.0 * (x - 1.0)))
            if 0.0 < q_star < ψ
                g_max = max(g_max, _om_plummer_df(q_star, x))
            end
        end
        envelope = ṽ_esc² * g_max

        # Rejection sampling of (ṽ, cos η) with weight ṽ² g(q)
        ṽ = 0.0
        cosη = 0.0
        while true
            ṽ = ṽ_esc * rand(rng)
            cosη = 2.0 * rand(rng) - 1.0
            q = ψ - 0.5 * ṽ^2 * (1.0 + ra_term * (1.0 - cosη^2))
            q > 0.0 && rand(rng) * envelope ≤ ṽ^2 * _om_plummer_df(q, x) && break
        end
        sinη = sqrt(1.0 - cosη^2)

        # Tangential unit vector: orthogonalise the Cartesian axis least
        # aligned with r̂, then rotate by a uniform azimuth about r̂.
        ax, ay, az = abs(r̂x), abs(r̂y), abs(r̂z)
        ex, ey, ez = if ax ≤ ay && ax ≤ az
            (1.0, 0.0, 0.0)
        elseif ay ≤ az
            (0.0, 1.0, 0.0)
        else
            (0.0, 0.0, 1.0)
        end
        er = ex * r̂x + ey * r̂y + ez * r̂z
        t1x = ex - er * r̂x
        t1y = ey - er * r̂y
        t1z = ez - er * r̂z
        t1n = sqrt(t1x^2 + t1y^2 + t1z^2)
        t1x /= t1n
        t1y /= t1n
        t1z /= t1n
        t2x = r̂y * t1z - r̂z * t1y
        t2y = r̂z * t1x - r̂x * t1z
        t2z = r̂x * t1y - r̂y * t1x
        χ = 2π * rand(rng)
        cχ, sχ = cos(χ), sin(χ)
        t̂x = cχ * t1x + sχ * t2x
        t̂y = cχ * t1y + sχ * t2y
        t̂z = cχ * t1z + sχ * t2z

        v = ṽ * inv_sqrt_a
        vel[1, i] = v * (cosη * r̂x + sinη * t̂x)
        vel[2, i] = v * (cosη * r̂y + sinη * t̂y)
        vel[3, i] = v * (cosη * r̂z + sinη * t̂z)
    end

    return pos, vel
end

"""
    lynden_bell_flip!(vel, pos, axis, p; rng = Random.default_rng()) -> vel

Lynden-Bell (1960) rotation by sign flips (DOI 10.1093/mnras/120.3.204).
For every particle off the rotation axis `axis` (any non-zero 3-vector), the
azimuthal velocity `v_φ` about that axis is given the positive sign with
probability `p` and the negative sign otherwise, by reflecting the velocity
in the meridional plane (`v ← v − 2 v_φ φ̂`). Particles on the axis are left
untouched.

The density, the speeds, the radial velocities and every velocity
dispersion are unchanged, so an equilibrium stays an equilibrium; the mean
azimuthal velocity in any shell becomes `(2p − 1)⟨|v_φ|⟩`. `p = ½` is the
non-rotating null, `p = 1` the maximal rotation the given model admits.
`vel` is overwritten and returned.

# Example
```julia
pos, vel = sample_plummer_om(10_000, 1.0, Inf)
lynden_bell_flip!(vel, pos, [0.0, 0.0, 1.0], 0.8)
```
"""
function lynden_bell_flip!(
    vel::AbstractMatrix{Float64},
    pos::AbstractMatrix{Float64},
    axis::AbstractVector{<:Real},
    p::Real;
    rng::AbstractRNG = Random.default_rng(),
)
    0 ≤ p ≤ 1 || throw(ArgumentError("p must lie in [0, 1], got $p"))
    size(pos, 1) == 3 || throw(ArgumentError("pos must have 3 rows, got size $(size(pos))"))
    size(vel) == size(pos) ||
        throw(ArgumentError("vel must have the size of pos $(size(pos)), got $(size(vel))"))
    length(axis) == 3 || throw(ArgumentError("axis must have length 3, got $(length(axis))"))
    a1, a2, a3 = Float64(axis[1]), Float64(axis[2]), Float64(axis[3])
    anorm = sqrt(a1^2 + a2^2 + a3^2)
    anorm > 0 || throw(ArgumentError("axis must be non-zero"))
    âx, ây, âz = a1 / anorm, a2 / anorm, a3 / anorm
    pf = Float64(p)

    for i in axes(pos, 2)
        rx, ry, rz = pos[1, i], pos[2, i], pos[3, i]
        r_par = rx * âx + ry * ây + rz * âz
        px = rx - r_par * âx
        py = ry - r_par * ây
        pz = rz - r_par * âz
        R = sqrt(px^2 + py^2 + pz^2)
        R == 0 && continue
        φx = (ây * pz - âz * py) / R
        φy = (âz * px - âx * pz) / R
        φz = (âx * py - ây * px) / R
        v_φ = vel[1, i] * φx + vel[2, i] * φy + vel[3, i] * φz
        positive = rand(rng) < pf
        if (positive && v_φ < 0) || (!positive && v_φ > 0)
            vel[1, i] -= 2 * v_φ * φx
            vel[2, i] -= 2 * v_φ * φy
            vel[3, i] -= 2 * v_φ * φz
        end
    end

    return vel
end

"""
    stretch_positions!(pos, q, s) -> pos

Scale the `y` coordinates by `q` and the `z` coordinates by `s`, turning a
spherical model into a static one with axis ratios `b/a = q` and `c/a = s`
(for `1 ≥ q ≥ s`) about the `x`, `y`, `z` axes, the known answer for the
shape estimator (Zemp et al. 2011, DOI 10.1088/0067-0049/197/2/30). The
result is not an equilibrium; velocities are left to the caller. `pos`
(`3 × N`) is overwritten and returned.

# Example
```julia
pos, _ = sample_plummer(10_000, 1.0)
stretch_positions!(pos, 0.8, 0.6)
```
"""
function stretch_positions!(pos::AbstractMatrix{Float64}, q::Real, s::Real)
    q > 0 || throw(ArgumentError("q must be positive, got $q"))
    s > 0 || throw(ArgumentError("s must be positive, got $s"))
    size(pos, 1) == 3 || throw(ArgumentError("pos must have 3 rows, got size $(size(pos))"))
    qf, sf = Float64(q), Float64(s)
    for i in axes(pos, 2)
        pos[2, i] *= qf
        pos[3, i] *= sf
    end
    return pos
end

"""
    rank_masses_by_energy!(mass, energy, S; rng = Random.default_rng()) -> mass

Mass segregation of tunable degree `S` for the estimator calibration, a scheme
of its own: it is neither the `S` of Šubr, Kroupa & Baumgardt (2008, DOI
10.1111/j.1365-2966.2008.12993.x) nor that of Baumgardt, De Marchi & Kroupa
(2008, DOI 10.1086/590488), and with positions kept and masses permuted it
changes the mass-density profile, so it is not an initial-condition
segregation method. A random fraction `S` of the particles (`round(Int, S n)` of `n`) receives
the largest masses, ranked so that the most bound of them (lowest `energy`)
carries the largest mass; the remaining masses go to the remaining particles
in random order. `S = 0` is a random permutation of the masses, `S = 1` the
full energy ranking.

The mass set is preserved exactly — the result is a permutation of the
input — and only its assignment to particles changes; `energy` (per-particle
specific or total energy, any consistent units) is not modified. `mass` is
overwritten and returned.

# Example
```julia
rank_masses_by_energy!(mass, energy, 0.5; rng = Xoshiro(1))
```
"""
function rank_masses_by_energy!(
    mass::AbstractVector{Float64},
    energy::AbstractVector{Float64},
    S::Real;
    rng::AbstractRNG = Random.default_rng(),
)
    0 ≤ S ≤ 1 || throw(ArgumentError("S must lie in [0, 1], got $S"))
    length(energy) == length(mass) || throw(
        ArgumentError(
            "energy must have the length of mass ($(length(mass))), got $(length(energy))",
        ),
    )
    n = length(mass)
    n_s = round(Int, S * n)

    order = Random.randperm(rng, n)
    chosen = order[1:n_s]
    rest = order[(n_s + 1):n]
    chosen = chosen[sortperm(view(energy, chosen))]
    sorted = sort(collect(mass); rev = true)
    remaining = Random.shuffle(rng, sorted[(n_s + 1):n])

    assigned = similar(sorted)
    for k in 1:n_s
        assigned[chosen[k]] = sorted[k]
    end
    for k in eachindex(rest)
        assigned[rest[k]] = remaining[k]
    end
    mass .= assigned
    return mass
end
