# =============================================================================
# Cluster density profile samplers: Plummer and King models
# =============================================================================

"""
    sample_plummer(N::Int, a::Float64; rng=Random.default_rng()) -> (pos, vel)

Generate `N` particles from a Plummer model with scale radius `a` using
inverse CDF sampling (Aarseth, Hénon & Wielen 1974).

Returns `(pos::Matrix{Float64}, vel::Matrix{Float64})` each `3 × N`,
in internal units where `G = M_total = 1` and virial radius `r_v = (3π/16) a`.
"""
function sample_plummer(N::Int, a::Float64; rng::AbstractRNG = Random.default_rng())
    pos = zeros(Float64, 3, N)
    vel = zeros(Float64, 3, N)

    for i in 1:N
        # Inverse CDF for radius: r = a / sqrt(X^{-2/3} - 1)
        X = rand(rng)
        r = a / sqrt(X^(-2/3) - 1.0)

        # Uniform direction on sphere
        cosθ = 2.0 * rand(rng) - 1.0
        sinθ = sqrt(1.0 - cosθ^2)
        φ = 2π * rand(rng)

        pos[1, i] = r * sinθ * cos(φ)
        pos[2, i] = r * sinθ * sin(φ)
        pos[3, i] = r * cosθ

        # Escape velocity at radius r
        v_esc = sqrt(2.0) * (r^2 + a^2)^(-0.25)

        # Rejection sampling for velocity magnitude (von Neumann)
        while true
            q = rand(rng)
            g = rand(rng)
            # f(q) ∝ q^2 (1 - q^2)^{7/2}, max at q = sqrt(1/9) → f_max = 0.0920
            if g ≤ q^2 * (1.0 - q^2)^3.5
                v = q * v_esc
                # Isotropic velocity direction
                cosθv = 2.0 * rand(rng) - 1.0
                sinθv = sqrt(1.0 - cosθv^2)
                φv = 2π * rand(rng)
                vel[1, i] = v * sinθv * cos(φv)
                vel[2, i] = v * sinθv * sin(φv)
                vel[3, i] = v * cosθv
                break
            end
        end
    end

    return pos, vel
end

# ---------------------------------------------------------------------------
# King model via ODE integration + rejection sampling
# ---------------------------------------------------------------------------

"""
    _king_ode!(du, u, p, ρ̂)

Right-hand side for the King model ODE in dimensionless form.
`u = [Ŵ, dŴ/dρ̂]` where `Ŵ = (ψ - ψ_t) / σ²` is the lowered potential and
`ρ̂ = r / r_0` is the dimensionless radius.

The King density is:
```math
\\hat{\\rho}(\\hat{W}) = e^{\\hat{W}} \\operatorname{erf}(\\sqrt{\\hat{W}})
    - \\sqrt{\\frac{4\\hat{W}}{\\pi}} \\left(1 + \\frac{2\\hat{W}}{3}\\right)
```
"""
function _king_density(W::Float64)::Float64
    W ≤ 0.0 && return 0.0
    sqW = sqrt(W)
    return exp(W) * erf(sqW) - sqrt(4.0 * W / π) * (1.0 + 2.0 * W / 3.0)
end

"""
    _solve_king(W0; N_radial=10000) -> (r̂, Ŵ, ρ̂_king)

Integrate the King model ODE from the centre (`Ŵ(0) = W0`) outward until
`Ŵ → 0` (the tidal radius). Returns dimensionless radius, potential, and
density arrays.
"""
function _solve_king(W0::Float64; N_radial::Int = 10000)
    # Poisson equation in spherical symmetry:
    #   d²Ŵ/dρ̂² + (2/ρ̂) dŴ/dρ̂ = -9 ρ̂_king(Ŵ)
    # where the factor 9 comes from 4πGρ₀ / (9σ²/4πGρ₀r₀²) normalisation.

    # Central boundary: Ŵ(0) = W0, dŴ/dρ̂(0) = 0
    # Near origin: Ŵ ≈ W0 - (3/2) ρ̂_king(W0) ρ̂²

    ρ0_king = _king_density(W0)

    # Adaptive step: estimate tidal radius from concentration
    # For W0 ~ 1-12, r_t/r_0 ranges from ~2 to ~1000+
    rhat_max = 3.0 * 10.0^(0.6 * W0 / 3.0)  # generous upper bound
    dr = rhat_max / N_radial

    rhat = zeros(Float64, N_radial + 1)
    What = zeros(Float64, N_radial + 1)
    ρ_arr = zeros(Float64, N_radial + 1)

    What[1] = W0
    ρ_arr[1] = ρ0_king

    # Start slightly off centre to avoid 2/r singularity
    rhat[1] = 0.0
    # Use Taylor expansion for the first step
    What[2] = W0 - 1.5 * ρ0_king * dr^2
    rhat[2] = dr

    i_tidal = N_radial + 1  # index of tidal radius

    for i in 2:N_radial
        r = rhat[i]
        W = What[i]
        if W ≤ 0.0
            i_tidal = i
            What[i] = 0.0
            break
        end
        ρ_arr[i] = _king_density(W)

        # Finite difference for 2nd order ODE (Störmer-Verlet like)
        # Ŵ_{i+1} = 2Ŵ_i - Ŵ_{i-1} + dr² [-9 ρ_king(Ŵ_i) - (2/r_i)(Ŵ_i - Ŵ_{i-1})/dr]
        dWdr = (What[i] - What[i-1]) / dr
        d2Wdr2 = -9.0 * ρ_arr[i] - 2.0 / r * dWdr

        What[i+1] = What[i] + dWdr * dr + 0.5 * d2Wdr2 * dr^2
        rhat[i+1] = r + dr

        if What[i+1] ≤ 0.0
            # Linear interpolation for tidal radius
            frac = What[i] / (What[i] - What[i+1])
            rhat[i+1] = rhat[i] + frac * dr
            What[i+1] = 0.0
            ρ_arr[i+1] = 0.0
            i_tidal = i + 1
            break
        end
    end

    return rhat[1:i_tidal], What[1:i_tidal], ρ_arr[1:i_tidal]
end

"""
    sample_king(N::Int, W0::Float64, rt::Float64;
                rng=Random.default_rng()) -> (pos, vel)

Generate `N` particles from a King (1966) model with central dimensionless
potential `W0` and tidal radius `rt`.

Uses the standard recipe:
1. Integrate the King ODE to get the density-potential pair `(ρ̂, Ŵ)(r̂)`
2. Rejection-sample particle positions from the density profile
3. At each radius, rejection-sample velocities from the lowered Maxwellian

Returns `(pos::Matrix{Float64}, vel::Matrix{Float64})` each `3 × N`,
in units where the tidal radius equals `rt`.
"""
function sample_king(N::Int, W0::Float64, rt::Float64;
                     rng::AbstractRNG = Random.default_rng())
    W0 > 0.0 || throw(ArgumentError("W0 must be positive, got $W0"))
    rt > 0.0 || throw(ArgumentError("rt must be positive, got $rt"))

    rhat, What, ρ_arr = _solve_king(W0)

    rhat_t = rhat[end]  # dimensionless tidal radius
    scale = rt / rhat_t  # physical scale: r_phys = scale * r̂

    # Build CDF for radial position sampling
    # P(< r̂) ∝ ∫₀^r̂ ρ(r̂') r̂'² dr̂'
    n_pts = length(rhat)
    cdf = zeros(Float64, n_pts)
    for i in 2:n_pts
        dr = rhat[i] - rhat[i-1]
        # Trapezoidal integration
        cdf[i] = cdf[i-1] + 0.5 * dr * (ρ_arr[i-1] * rhat[i-1]^2 + ρ_arr[i] * rhat[i]^2)
    end
    cdf ./= cdf[end]  # normalise

    # Interpolation helper: find r̂ from CDF value via binary search + linear interp
    function _interp_r(u::Float64)
        # Binary search for the bracket
        lo, hi = 1, n_pts
        while hi - lo > 1
            mid = (lo + hi) >> 1
            cdf[mid] < u ? (lo = mid) : (hi = mid)
        end
        # Linear interpolation within bracket
        frac = (u - cdf[lo]) / (cdf[hi] - cdf[lo] + 1e-30)
        return rhat[lo] + frac * (rhat[hi] - rhat[lo])
    end

    # Interpolation helper: find Ŵ at given r̂
    function _interp_W(r::Float64)
        r ≥ rhat[end] && return 0.0
        lo, hi = 1, n_pts
        while hi - lo > 1
            mid = (lo + hi) >> 1
            rhat[mid] < r ? (lo = mid) : (hi = mid)
        end
        frac = (r - rhat[lo]) / (rhat[hi] - rhat[lo] + 1e-30)
        return What[lo] + frac * (What[hi] - What[lo])
    end

    # σ² in physical units: from the King model, σ² = GM/(9 r_0) × (some factor)
    # We work in units where σ = 1 during sampling, then rescale at the end.
    # The velocity dispersion σ² relates to the King parameter as:
    # v_esc(r) = sqrt(2 Ŵ(r̂)) × σ

    pos = zeros(Float64, 3, N)
    vel = zeros(Float64, 3, N)

    for i in 1:N
        # Sample radius from CDF
        u = rand(rng)
        rhat_i = _interp_r(u)
        r_phys = rhat_i * scale

        # Uniform direction
        cosθ = 2.0 * rand(rng) - 1.0
        sinθ = sqrt(1.0 - cosθ^2)
        φ = 2π * rand(rng)

        pos[1, i] = r_phys * sinθ * cos(φ)
        pos[2, i] = r_phys * sinθ * sin(φ)
        pos[3, i] = r_phys * cosθ

        # Local escape speed (in σ units)
        W_local = _interp_W(rhat_i)
        v_esc = sqrt(2.0 * max(W_local, 0.0))

        # Rejection-sample velocity from lowered Maxwellian:
        # f(v) ∝ v² [exp(-(v²/2 - W)) - 1] for v < v_esc, 0 otherwise
        # = v² [exp(W - v²/2) - 1]
        # Maximum at v* where d/dv [v² (e^{W-v²/2} - 1)] = 0
        if W_local ≤ 1e-10
            # Near tidal radius: essentially zero velocity
            vel[:, i] .= 0.0
            continue
        end

        # Envelope: f_max ≤ v_esc² × exp(W_local)
        f_max = v_esc^2 * exp(W_local)

        while true
            v = v_esc * rand(rng)
            f_v = v^2 * (exp(W_local - 0.5 * v^2) - 1.0)
            if rand(rng) * f_max ≤ f_v
                cosθv = 2.0 * rand(rng) - 1.0
                sinθv = sqrt(1.0 - cosθv^2)
                φv = 2π * rand(rng)
                vel[1, i] = v * sinθv * cos(φv)
                vel[2, i] = v * sinθv * sin(φv)
                vel[3, i] = v * cosθv
                break
            end
        end
    end

    # Rescale velocities to be consistent with positions.
    # In N-body units with G=1, M_total=1:
    # σ² = M / (6 r_0) for King models, where r_0 = scale × 1
    # But we need to set the velocity scale self-consistently.
    # For now, velocities are in units of σ; the virialisation step in
    # `combine_clusters!` will handle the final scaling.

    return pos, vel
end

"""
    virialise!(mass, pos, vel)

Shift to centre-of-mass frame and scale velocities so that the virial
ratio `Q = T/|W| = 0.5` (virial equilibrium). Operates in-place.

Assumes `G = 1`. Uses the exact N-body potential energy (O(N²) — fine for
IC generation with N ≤ 10⁶).
"""
function virialise!(mass::Vector{Float64}, pos::Matrix{Float64}, vel::Matrix{Float64})
    N = length(mass)
    M_total = sum(mass)

    # Centre of mass correction
    cm_pos = zeros(Float64, 3)
    cm_vel = zeros(Float64, 3)
    for i in 1:N
        for k in 1:3
            cm_pos[k] += mass[i] * pos[k, i]
            cm_vel[k] += mass[i] * vel[k, i]
        end
    end
    cm_pos ./= M_total
    cm_vel ./= M_total

    for i in 1:N
        for k in 1:3
            pos[k, i] -= cm_pos[k]
            vel[k, i] -= cm_vel[k]
        end
    end

    # Kinetic energy
    T = 0.0
    for i in 1:N
        v2 = vel[1, i]^2 + vel[2, i]^2 + vel[3, i]^2
        T += 0.5 * mass[i] * v2
    end

    # Potential energy (O(N²))
    W = 0.0
    for i in 1:N
        for j in (i+1):N
            dx = pos[1, i] - pos[1, j]
            dy = pos[2, i] - pos[2, j]
            dz = pos[3, i] - pos[3, j]
            r = sqrt(dx^2 + dy^2 + dz^2)
            W -= mass[i] * mass[j] / r
        end
    end

    # Scale velocities: Q_target = 0.5 → T_new = 0.5 |W|
    # v_new = v_old × sqrt(0.5 |W| / T)
    if T > 0.0 && W < 0.0
        scale_v = sqrt(0.5 * abs(W) / T)
        vel .*= scale_v
    end

    return nothing
end
