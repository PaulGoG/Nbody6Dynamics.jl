# =============================================================================
# Rotating King model (Einsel & Spurzem 1999): self-consistent axisymmetric
# equilibrium by Legendre–Poisson iteration, and exact phase-space sampling
# =============================================================================
#
# Units throughout: G = 1, σ = β^{-1/2} = 1, r_K = 1 with r_K² = 9σ²/(4πGρ_c);
# densities in units of the central density ρ_c, masses in ρ_c r_K³,
# potentials and energies per unit mass in σ².

# Gauss–Legendre nodes of the speed quadrature of the local moments.
const _RK_SPEED_NODES = 64
# Below this |z| the kernels of the first and second v_φ moments are summed
# as power series (the closed forms lose digits to cancellation).
const _RK_SERIES_THRESHOLD = 0.5
# Above this z the e^{−z} parts of the kernels are dropped (no Inf − Inf).
const _RK_ASYMPTOTIC_THRESHOLD = 20.0
# Points of the scan that bounds the King speed law in the velocity sampler,
# and the margin applied to the scanned maximum.
const _RK_SCAN_POINTS = 200
const _RK_SCAN_MARGIN = 1.001
# Iterations after which the truncation surface must lie inside r_max.
const _RK_CLOSURE_ITERATIONS = 50

"""
    RotatingKingModel

Self-consistent rotating King model (Einsel & Spurzem 1999, DOI
10.1046/j.1365-8711.1999.02083.x) with distribution function
`f(E, J_z) = C (e^{−βE} − 1) e^{βΩ₀J_z}` for `E < 0`, tabulated on a polar
grid `(r̂, μ = cos θ)`, as returned by [`solve_rotating_king`](@ref). All
quantities are dimensionless (`G = σ = r_K = 1`).

# Fields
- `W0::Float64`: central dimensionless potential, `ψ(0) = W₀`.
- `ω0::Float64`: dimensionless angular velocity `ω₀ = Ω₀ r_K/σ`.
- `r::Vector{Float64}`: uniform radial grid `r̂ = r/r_K`, `n_r` points from 0 to `r̂_max`.
- `μ::Vector{Float64}`: Gauss–Legendre nodes in `μ = cos θ` on `[−1, 1]`, ascending, `n_μ` points.
- `ψ::Matrix{Float64}`: `ψ = −β(Φ − Φ_t)` at `(r̂_i, μ_j)`, `n_r × n_μ`; zero on the truncation surface,
  negative outside it.
- `ρ::Matrix{Float64}`: density `ρ̂(ψ, ω₀R̂)/ρ̂(W₀, 0)` at `(r̂_i, μ_j)`, unity at the centre and zero
  outside the truncation surface.
- `l_max::Int`: highest (even) Legendre order of the potential expansion.
- `r_t_equator::Float64`: radius of the surface `ψ = 0` in the equatorial plane (`μ = 0`).
- `r_t_pole::Float64`: radius of the surface `ψ = 0` on the rotation axis (`μ = ±1`).
- `mass::Float64`: total mass `M̂ = ∫ρ̂ dV̂` in units of `ρ_c r_K³`.
- `r_h::Float64`: radius of the sphere enclosing `M̂/2`.
- `iterations::Int`: fixed-point iterations performed.
- `residual::Float64`: final `max |ψ^{n+1} − ψ^n|` of the undamped update.
"""
struct RotatingKingModel
    W0::Float64
    ω0::Float64
    r::Vector{Float64}
    μ::Vector{Float64}
    ψ::Matrix{Float64}
    ρ::Matrix{Float64}
    l_max::Int
    r_t_equator::Float64
    r_t_pole::Float64
    mass::Float64
    r_h::Float64
    iterations::Int
    residual::Float64
end

# ---------------------------------------------------------------------------
# Quadrature and Legendre polynomials
# ---------------------------------------------------------------------------

"""
    _rk_legendre(l, x) -> Float64

Legendre polynomial `P_l(x)` by the three-term (Bonnet) recurrence.
"""
function _rk_legendre(l::Integer, x::Float64)::Float64
    l == 0 && return 1.0
    p_prev, p = 1.0, x
    for n in 2:l
        p_prev, p = p, ((2n - 1) * x * p - (n - 1) * p_prev) / n
    end
    return p
end

"""
    _rk_legendre_derivative(n, z) -> (P_n(z), P_n′(z))

`P_n` and its derivative from `P_n′ = n (z P_n − P_{n−1})/(z² − 1)`, for `|z| < 1`.
"""
function _rk_legendre_derivative(n::Integer, z::Float64)
    p_prev, p = 1.0, z
    for k in 2:n
        p_prev, p = p, ((2k - 1) * z * p - (k - 1) * p_prev) / k
    end
    return p, n * (z * p - p_prev) / (z^2 - 1)
end

"""
    _rk_gauss_legendre(n) -> (x, w)

Nodes `x` (ascending) and weights `w` of the `n`-point Gauss–Legendre rule on
`[−1, 1]`, by Newton iteration on `P_n` from the Tricomi initial guesses
`cos(π(i − 1/4)/(n + 1/2))`.
"""
function _rk_gauss_legendre(n::Integer)
    n ≥ 1 || throw(ArgumentError("_rk_gauss_legendre: n must be ≥ 1; got $n"))
    x = zeros(Float64, n)
    w = zeros(Float64, n)
    for i in 1:((n + 1) ÷ 2)
        z = cos(π * (i - 0.25) / (n + 0.5))
        for _ in 1:100
            p, dp = _rk_legendre_derivative(n, z)
            δ = p / dp
            z -= δ
            abs(δ) ≤ 4 * eps(Float64) && break
        end
        _, dp = _rk_legendre_derivative(n, z)
        x[i] = -z
        x[n + 1 - i] = z
        w[i] = 2 / ((1 - z^2) * dp^2)
        w[n + 1 - i] = w[i]
    end
    return x, w
end

# ---------------------------------------------------------------------------
# Local velocity moments at (ψ, k = ω₀R̂)
# ---------------------------------------------------------------------------

"""
    _rk_sinhc(z) -> Float64

`½∫_{−1}^{1} e^{zx} dx = sinh(z)/z`, by the series `1 + z²/6` for `|z| < 10⁻⁴`.
"""
function _rk_sinhc(z::Float64)::Float64
    a = abs(z)
    a < 1e-4 && return 1 + a^2 / 6
    return sinh(a) / a
end

"""
    _rk_vphi_kernel(z) -> Float64

`½∫_{−1}^{1} x e^{zx} dx = cosh(z)/z − sinh(z)/z²` (odd in `z`, `→ z/3` as
`z → 0`); power series `Σ z^{2n+1} (2n + 2)/(2n + 3)!` for `|z| < 0.5`.
"""
function _rk_vphi_kernel(z::Float64)::Float64
    a = abs(z)
    val = if a < _RK_SERIES_THRESHOLD
        a2 = a^2
        a * (
            1 / 3 +
            a2 * (
                1 / 30 +
                a2 * (
                    1 / 840 +
                    a2 * (1 / 45360 + a2 * (1 / 3991680 + a2 * (1 / 518918400 + a2 / 93405312000)))
                )
            )
        )
    elseif a > _RK_ASYMPTOTIC_THRESHOLD
        exp(a) * (a - 1) / (2 * a^2)
    else
        (a * cosh(a) - sinh(a)) / a^2
    end
    return copysign(val, z)
end

"""
    _rk_vphi2_kernel(z) -> Float64

`½∫_{−1}^{1} x² e^{zx} dx = [(z² + 2) sinh z − 2z cosh z]/z³` (even in `z`,
`→ 1/3` as `z → 0`); power series `Σ z^{2n}/((2n)! (2n + 3))` for `|z| < 0.5`.
"""
function _rk_vphi2_kernel(z::Float64)::Float64
    a = abs(z)
    if a < _RK_SERIES_THRESHOLD
        a2 = a^2
        return 1 / 3 +
               a2 * (
            1 / 10 +
            a2 *
            (1 / 168 + a2 * (1 / 6480 + a2 * (1 / 443520 + a2 * (1 / 47174400 + a2 / 7185024000))))
        )
    elseif a > _RK_ASYMPTOTIC_THRESHOLD
        return exp(a) * (a^2 - 2a + 2) / (2 * a^3)
    end
    return ((a^2 + 2) * sinh(a) - 2a * cosh(a)) / a^3
end

"""
    _rk_density(ψ, k, x, w) -> Float64

Unnormalised local density of the rotating King distribution function,
`ρ̂(ψ, k) = 4π ∫₀^{√(2ψ)} u² (e^{ψ − u²/2} − 1) sinh(ku)/(ku) du`, by the
Gauss–Legendre rule `(x, w)` mapped onto `[0, √(2ψ)]`; zero for `ψ ≤ 0`. At
`k = 0` it equals `(2π)^{3/2}` times [`_king_density`](@ref).
"""
function _rk_density(
    ψ::Real,
    k::Real,
    x::AbstractVector{Float64},
    w::AbstractVector{Float64},
)::Float64
    ψf = Float64(ψ)
    kf = Float64(k)
    ψf > 0 || return 0.0
    a = sqrt(2ψf)
    s = 0.0
    @inbounds for q in eachindex(x, w)
        u = 0.5 * a * (1 + x[q])
        s += w[q] * u^2 * expm1(ψf - 0.5 * u^2) * _rk_sinhc(kf * u)
    end
    return 2π * a * s
end

"""
    _rk_moments(ψ, k, x, w) -> NamedTuple{(:ρ, :ρv_φ, :ρv_φ², :ρv_R², :ρu²)}

Unnormalised local moments of the rotating King distribution function at
`(ψ, k = ω₀R̂)`, in units of `σ` for the velocities, by the Gauss–Legendre
rule `(x, w)` on `[0, √(2ψ)]` with `g(u) = e^{ψ − u²/2} − 1`:

- `ρ = 4π ∫ u² g sinh(ku)/(ku) du`;
- `ρv_φ = ρ⟨v_φ⟩ = 4π ∫ u³ g [cosh(ku)/(ku) − sinh(ku)/(ku)²] du`;
- `ρv_φ² = ρ⟨v_φ²⟩ = 4π ∫ u⁴ g ½∫x² e^{kux} dx du`;
- `ρu² = ρ⟨u²⟩ = 4π ∫ u⁴ g sinh(ku)/(ku) du`;
- `ρv_R² = ρ⟨v_R²⟩ = ρ⟨v_z²⟩ = ½(ρu² − ρv_φ²)`.

All vanish for `ψ ≤ 0`.
"""
function _rk_moments(ψ::Real, k::Real, x::AbstractVector{Float64}, w::AbstractVector{Float64})
    ψf = Float64(ψ)
    kf = Float64(k)
    ψf > 0 || return (ρ = 0.0, ρv_φ = 0.0, ρv_φ² = 0.0, ρv_R² = 0.0, ρu² = 0.0)
    a = sqrt(2ψf)
    s0 = 0.0
    s1 = 0.0
    s2 = 0.0
    s3 = 0.0
    @inbounds for q in eachindex(x, w)
        u = 0.5 * a * (1 + x[q])
        g = w[q] * u^2 * expm1(ψf - 0.5 * u^2)
        z = kf * u
        sc = _rk_sinhc(z)
        s0 += g * sc
        s1 += g * u * _rk_vphi_kernel(z)
        s2 += g * u^2 * _rk_vphi2_kernel(z)
        s3 += g * u^2 * sc
    end
    c = 2π * a
    ρu² = c * s3
    ρv_φ² = c * s2
    return (ρ = c * s0, ρv_φ = c * s1, ρv_φ² = ρv_φ², ρv_R² = 0.5 * (ρu² - ρv_φ²), ρu² = ρu²)
end

# ---------------------------------------------------------------------------
# Radial quadrature and the multipole potential
# ---------------------------------------------------------------------------

"""
    _rk_ghost(f, j, parity) -> Float64

`f[j]` on the uniform grid `s_j = (j − 1)h`, continued to `j < 1` by the
parity of the integrand about `s = 0` and by zero beyond the last point.
"""
@inline function _rk_ghost(f::AbstractVector{Float64}, j::Int, parity::Float64)::Float64
    j < 1 && return parity * f[2 - j]
    j > length(f) && return 0.0
    return f[j]
end

"""
    _rk_interval_integrals!(I, f, h, parity) -> I

`I[i] = ∫_{s_i}^{s_{i+1}} f ds` on a uniform grid of spacing `h`, from the
quintic through `s_{i−2} … s_{i+3}` (weights `(11, −93, 802, 802, −93, 11) h/1440`,
sixth order). Points beyond the grid are supplied by [`_rk_ghost`](@ref):
`parity = ±1` for an even or odd integrand at the origin, zero outside `s_max`.
"""
function _rk_interval_integrals!(
    I::AbstractVector{Float64},
    f::AbstractVector{Float64},
    h::Float64,
    parity::Float64,
)
    c = h / 1440
    @inbounds for i in 1:(length(f) - 1)
        I[i] =
            c * (
                11 * (_rk_ghost(f, i - 2, parity) + _rk_ghost(f, i + 3, parity)) -
                93 * (_rk_ghost(f, i - 1, parity) + _rk_ghost(f, i + 2, parity)) +
                802 * (f[i] + f[i + 1])
            )
    end
    return I
end

"""
    _rk_radial_integral(f, h) -> Float64

`∫₀^{s_max} f ds` of an integrand even about the origin, by
[`_rk_interval_integrals!`](@ref).
"""
function _rk_radial_integral(f::AbstractVector{Float64}, h::Float64)::Float64
    I = zeros(Float64, length(f) - 1)
    _rk_interval_integrals!(I, f, h, 1.0)
    return sum(I)
end

"""
    _rk_multipole_potential(r, μ, wμ, ρ, l_max) -> Matrix{Float64}

Gravitational potential, in units of `σ²` and zero at infinity, of the density
`ρ` (in units of `ρ_c`) tabulated on the uniform radial grid `r` (in `r_K`)
and the Gauss–Legendre nodes `μ` with weights `wμ`. The density is projected
on the even Legendre polynomials, `ρ_l(r) = (2l + 1)/2 Σ_j wμ_j ρ(r, μ_j) P_l(μ_j)`,
and each multipole is integrated with the Green's function

`Φ̂_l(r) = −(9/(2l + 1)) [r^{−l−1} ∫₀^r ρ_l s^{l+2} ds + r^l ∫_r^{r_max} ρ_l s^{1−l} ds]`,

the factor `9 = 4πGρ_c r_K²/σ²` converting the Newtonian potential of
`ρ/ρ_c` into units of `σ²`. Both integrals use product-trapezoid rules (`ρ_l`
linear on each grid interval, the power of `s` integrated exactly), the inner
ones accumulated outward and the outer ones inward; on the first interval the
regular behaviour `ρ_l ∝ s^l` is used, and `r = 0` takes the limits
`Φ̂_0(0) = −9 ∫ ρ_0 s ds`, `Φ̂_{l≥2}(0) = 0`.
"""
function _rk_multipole_potential(
    r::AbstractVector{Float64},
    μ::AbstractVector{Float64},
    wμ::AbstractVector{Float64},
    ρ::AbstractMatrix{Float64},
    l_max::Integer,
)::Matrix{Float64}
    n_r, n_μ = size(ρ)
    Φ = zeros(Float64, n_r, n_μ)
    ρ_l = zeros(Float64, n_r)
    A = zeros(Float64, n_r)
    B = zeros(Float64, n_r)
    Φ_l = zeros(Float64, n_r)
    P = zeros(Float64, n_μ)
    for l in 0:2:l_max
        for j in 1:n_μ
            P[j] = _rk_legendre(l, μ[j])
        end
        fill!(ρ_l, 0.0)
        for j in 1:n_μ
            p = (2l + 1) / 2 * wμ[j] * P[j]
            @inbounds for i in 1:n_r
                ρ_l[i] += p * ρ[i, j]
            end
        end
        # Product-trapezoid rules: ρ_l linear on each interval, the powers of s
        # integrated exactly. A polynomial rule on s^{1-l} near the origin
        # (singular for l ≥ 2) distorts the kernel at the innermost radii and
        # feeds a growing mode of the iteration there.
        A[1] = 0.0
        @inbounds for i in 1:(n_r - 1)
            s0, s1 = r[i], r[i + 1]
            q = l + 3
            A[i + 1] = A[i] + 0.5 * (ρ_l[i] + ρ_l[i + 1]) * (s1^q - s0^q) / q
        end
        B[n_r] = 0.0
        @inbounds for i in (n_r - 1):-1:2
            s0, s1 = r[i], r[i + 1]
            q = 2 - l
            w = q == 0 ? log(s1 / s0) : (s1^q - s0^q) / q
            B[i] = B[i + 1] + 0.5 * (ρ_l[i] + ρ_l[i + 1]) * w
        end
        # First interval [0, h]: ρ_l ∝ s^l by regularity, so the integrand is
        # ∝ s and the interval contributes ρ_l(h) h^{2-l}/2 (l = 0: the
        # trapezoid of ρ_l s).
        h1 = r[2]
        B[1] = B[2] + (l == 0 ? 0.25 * (ρ_l[1] + ρ_l[2]) * h1^2 : 0.5 * ρ_l[2] * h1^(2 - l))
        coefficient = 9 / (2l + 1)
        Φ_l[1] = l == 0 ? -coefficient * B[1] : 0.0
        @inbounds for i in 2:n_r
            s = r[i]
            Φ_l[i] = -coefficient * (A[i] / s^(l + 1) + s^l * B[i])
        end
        for j in 1:n_μ
            @inbounds for i in 1:n_r
                Φ[i, j] += Φ_l[i] * P[j]
            end
        end
    end
    return Φ
end

"""
    _rk_coefficients(ψ, μ, wμ, l_max) -> Matrix{Float64}

Even Legendre coefficients `ψ_l(r_i) = (2l + 1)/2 Σ_j wμ_j ψ(r_i, μ_j) P_l(μ_j)`,
`l = 0, 2, …, l_max` in columns `1, 2, …, l_max/2 + 1`.
"""
function _rk_coefficients(
    ψ::AbstractMatrix{Float64},
    μ::AbstractVector{Float64},
    wμ::AbstractVector{Float64},
    l_max::Integer,
)::Matrix{Float64}
    n_r, n_μ = size(ψ)
    n_l = l_max ÷ 2 + 1
    c = zeros(Float64, n_r, n_l)
    for q in 1:n_l
        l = 2 * (q - 1)
        for j in 1:n_μ
            p = (2l + 1) / 2 * wμ[j] * _rk_legendre(l, μ[j])
            @inbounds for i in 1:n_r
                c[i, q] += p * ψ[i, j]
            end
        end
    end
    return c
end

"""
    _rk_legendre_coefficients(m::RotatingKingModel) -> Matrix{Float64}

Even Legendre coefficients of `m.ψ` on the radial grid (see
[`_rk_coefficients`](@ref)), the input of [`_rk_psi_at`](@ref).
"""
function _rk_legendre_coefficients(m::RotatingKingModel)::Matrix{Float64}
    _, wμ = _rk_gauss_legendre(length(m.μ))
    return _rk_coefficients(m.ψ, m.μ, wμ, m.l_max)
end

"""
    _rk_psi_at(m, c, r, μ) -> Float64

`ψ` of the model `m` at radius `r` (in `r_K`) and `μ = cos θ`: the Legendre
series with coefficients `c` (from [`_rk_legendre_coefficients`](@ref))
interpolated linearly in `r` between the bracketing grid rows. Beyond `r_max`
the value at `r_max` is returned.
"""
function _rk_psi_at(m::RotatingKingModel, c::AbstractMatrix{Float64}, r::Real, μ::Real)::Float64
    n_r = length(m.r)
    h = m.r[2] - m.r[1]
    s = clamp(Float64(r) / h, 0.0, Float64(n_r))
    i = clamp(floor(Int, s) + 1, 1, n_r - 1)
    f = clamp(s - (i - 1), 0.0, 1.0)
    x = Float64(μ)
    acc = (1 - f) * c[i, 1] + f * c[i + 1, 1]
    p_prev, p = 1.0, x
    for l in 2:m.l_max
        p_prev, p = p, ((2l - 1) * x * p - (l - 1) * p_prev) / l
        if iseven(l)
            q = l ÷ 2 + 1
            acc += ((1 - f) * c[i, q] + f * c[i + 1, q]) * p
        end
    end
    return acc
end

"""
    _rk_surface_radius(c, r, l_max, μ0) -> Float64

Radius at which `ψ(r, μ0)`, the Legendre series with coefficients `c`, first
falls to zero, by linear interpolation between grid rows; `NaN` when it stays
positive out to the last row.
"""
function _rk_surface_radius(
    c::AbstractMatrix{Float64},
    r::AbstractVector{Float64},
    l_max::Integer,
    μ0::Float64,
)::Float64
    P = [_rk_legendre(l, μ0) for l in 0:2:l_max]
    value(i) = sum(c[i, q] * P[q] for q in eachindex(P))
    v_prev = value(1)
    for i in 2:length(r)
        v = value(i)
        if v ≤ 0
            return r[i - 1] + (r[i] - r[i - 1]) * v_prev / (v_prev - v)
        end
        v_prev = v
    end
    return NaN
end

# ---------------------------------------------------------------------------
# Self-consistent solution
# ---------------------------------------------------------------------------

"""
    _rk_no_equilibrium(ω0, reason) -> ArgumentError

Error raised by [`solve_rotating_king`](@ref) when no equilibrium is found.
"""
_rk_no_equilibrium(ω0::Float64, reason::AbstractString) =
    ArgumentError("solve_rotating_king: no equilibrium at ω0 = $(ω0): $(reason)")

"""
    _rk_initial_potential(r, rk, wk, ρk, n_μ) -> Matrix{Float64}

Starting iterate: the non-rotating King solution `(rk, wk, ρk)` of
[`_solve_king`](@ref) interpolated linearly onto the grid `r`, continued
outside its tidal radius `r_t` by the point-mass tail
`ψ = −9 m_K (1/r_t − 1/r)`, `m_K = ∫₀^{r_t} (ρ/ρ_c) s² ds`; identical for all
`n_μ` angular nodes.
"""
function _rk_initial_potential(
    r::AbstractVector{Float64},
    rk::AbstractVector{Float64},
    wk::AbstractVector{Float64},
    ρk::AbstractVector{Float64},
    n_μ::Integer,
)::Matrix{Float64}
    r_t = rk[end]
    m_K = 0.0
    for i in 2:length(rk)
        m_K += 0.5 * (rk[i] - rk[i - 1]) * (ρk[i - 1] * rk[i - 1]^2 + ρk[i] * rk[i]^2)
    end
    m_K /= ρk[1]
    ψ = zeros(Float64, length(r), n_μ)
    for i in eachindex(r)
        x = r[i]
        value = if x < r_t
            j = clamp(searchsortedlast(rk, x), 1, length(rk) - 1)
            f = (x - rk[j]) / (rk[j + 1] - rk[j])
            (1 - f) * wk[j] + f * wk[j + 1]
        else
            -9 * m_K * (1 / r_t - 1 / x)
        end
        ψ[i, :] .= value
    end
    return ψ
end

"""
    _rk_density_grid!(ρ, ψ, r, μ, ω0, ρ̂₀, x, w) -> ρ

Normalised density `ρ̂(ψ, ω₀ r √(1 − μ²))/ρ̂₀` on the grid, evaluated on the
nodes `μ_j ≤ 0` and mirrored (the solution is symmetric about the equator).
Raises the no-equilibrium `ArgumentError` on a non-finite potential or density.
"""
function _rk_density_grid!(
    ρ::AbstractMatrix{Float64},
    ψ::AbstractMatrix{Float64},
    r::AbstractVector{Float64},
    μ::AbstractVector{Float64},
    ω0::Float64,
    ρ̂₀::Float64,
    x::AbstractVector{Float64},
    w::AbstractVector{Float64},
)
    n_r, n_μ = size(ψ)
    for j in 1:cld(n_μ, 2)
        sinθ = sqrt(max(1 - μ[j]^2, 0.0))
        for i in 1:n_r
            ψ_ij = ψ[i, j]
            isfinite(ψ_ij) || throw(_rk_no_equilibrium(ω0, "non-finite potential"))
            value = ψ_ij > 0 ? _rk_density(ψ_ij, ω0 * r[i] * sinθ, x, w) / ρ̂₀ : 0.0
            isfinite(value) || throw(_rk_no_equilibrium(ω0, "non-finite density"))
            ρ[i, j] = value
            ρ[i, n_μ + 1 - j] = value
        end
    end
    return ρ
end

"""
    _rk_spherical_part!(ψ0, r, μ, wμ, ψh, ω0, ρ̂₀, x, w, W0) -> ψ0

Spherical part `ψ₀(r)` of the potential by outward integration of
`ψ₀'' + (2/r) ψ₀' = −9 ρ̄(r, ψ₀)` with `ψ₀(0) = W0`, `ψ₀'(0) = 0`, where
`ρ̄ = ½ Σ_j wμ_j ρ̂(ψ₀ + ψ_h(r, μ_j), ω₀ r sinθ_j)/ρ̂₀` is the angular average
of the density at the current non-spherical part `ψ_h = Σ_{l≥2} ψ_l P_l`
(interpolated linearly between the grid radii). Classical RK4 with four
sub-steps per grid interval; at the centre `ψ₀'' = −3ρ̄`. The initial-value
problem has no growing mode, whereas the Picard iteration of the `l = 0`
Poisson equation amplifies the core's Jeans mode (`∇²δψ = −9 ρ̂' δψ`, growth
factor ≈ 9/q² > 1 for the modes that fit the core) for `W0 ≳ 4`.
"""
function _rk_spherical_part!(
    ψ0::AbstractVector{Float64},
    r::AbstractVector{Float64},
    μ::AbstractVector{Float64},
    wμ::AbstractVector{Float64},
    ψh::AbstractMatrix{Float64},
    ω0::Float64,
    ρ̂₀::Float64,
    x::AbstractVector{Float64},
    w::AbstractVector{Float64},
    W0::Float64,
)
    nr = length(r)
    nμ = length(μ)
    h = r[2] - r[1]
    sinθ = [sqrt(max(1 - μ[j]^2, 0.0)) for j in 1:nμ]
    # Angular average of the density at radius rr inside grid interval i, for the spherical potential p
    function rhobar(rr::Float64, p::Float64, i::Int)
        f = clamp((rr - r[i]) / h, 0.0, 1.0)
        i2 = min(i + 1, nr)
        acc = 0.0
        @inbounds for j in 1:nμ
            ψt = p + (1 - f) * ψh[i, j] + f * ψh[i2, j]
            ψt > 0 || continue
            acc += wμ[j] * _rk_density(ψt, ω0 * rr * sinθ[j], x, w)
        end
        return 0.5 * acc / ρ̂₀
    end
    rhs(rr, p, u, i) =
        rr < 1e-12 ? (u, -3 * rhobar(rr, p, i)) : (u, -9 * rhobar(rr, p, i) - 2 * u / rr)
    ψ0[1] = W0
    p = W0
    u = 0.0
    hs = h / 4
    for i in 1:(nr - 1)
        rr = r[i]
        for _ in 1:4
            k1p, k1u = rhs(rr, p, u, i)
            k2p, k2u = rhs(rr + hs / 2, p + hs / 2 * k1p, u + hs / 2 * k1u, i)
            k3p, k3u = rhs(rr + hs / 2, p + hs / 2 * k2p, u + hs / 2 * k2u, i)
            k4p, k4u = rhs(rr + hs, p + hs * k3p, u + hs * k3u, i)
            p += hs / 6 * (k1p + 2k2p + 2k3p + k4p)
            u += hs / 6 * (k1u + 2k2u + 2k3u + k4u)
            rr += hs
        end
        ψ0[i + 1] = p
    end
    return ψ0
end

"""
    solve_rotating_king(W0, ω0; l_max = 12, n_r = 600, n_μ = 48, r_max = 0.0,
                        relaxation = 0.4, tol = 1e-9, max_iter = 2000) -> RotatingKingModel

Self-consistent axisymmetric equilibrium of the rotating King distribution
function `f(E, J_z) = C (e^{−βE} − 1) e^{βΩ₀J_z}` (`E < 0`) of Einsel &
Spurzem (1999, DOI 10.1046/j.1365-8711.1999.02083.x), with central potential
`W0 = ψ(0)` and dimensionless angular velocity `ω0 = Ω₀ r_K/σ ≥ 0`
(`ω0 = 0` is the King 1966 model, DOI 10.1086/109930).

Method (Prendergast–Tomer iteration as in Varri & Bertin 2012, DOI
10.1051/0004-6361/201118300). In the units `G = σ = r_K = 1`, the potential
`ψ = −β(Φ − Φ_t)` solves `∇̂²ψ = −9 ρ̂(ψ, R̂)/ρ̂(W₀, 0)`, where the local density
`ρ̂(ψ, k = ω₀R̂)` is a 64-node Gauss–Legendre quadrature over the speed. Starting
from the non-rotating King solution, each iteration evaluates the density on
the grid `(r̂_i, μ_j)` (`n_r` uniform radii from 0 to `r_max`, `n_μ`
Gauss–Legendre nodes in `μ = cos θ`), computes its potential by the Legendre
Green's function for `l = 0, 2, …, l_max` (`_rk_multipole_potential`,
product-trapezoid radial rules), takes the non-spherical part
`ψ_h^{n+1} = −[Φ̂(r̂, μ) − Φ̂₀(r̂)]` (the potential less its angular average)
under-relaxed, `ψ_h ← (1 − α) ψ_h + α ψ_h^{n+1}` with `α = relaxation`, and
then recomputes the spherical part `ψ₀(r̂)` exactly by outward integration of
the `l = 0` equation with that `ψ_h` (`_rk_spherical_part!`; a Picard step on
`l = 0` diverges on the core's Jeans mode for `W0 ≳ 4`), until
`max |ψ^{n+1} − ψ^n| < tol` for the whole `ψ = ψ₀ + ψ_h`. `r_max = 0` selects
1.5 times the tidal radius of the non-rotating model.

Raises `ArgumentError("solve_rotating_king: no equilibrium at ω0 = …")` when
the surface `ψ = 0` does not close inside `r_max` after 50 iterations, when
the density at `r_max` is non-zero at convergence, when the potential or the
density becomes non-finite, or when `max_iter` is reached (critical rotation).

Returns a [`RotatingKingModel`](@ref) in model units (`r_K = 1`, `σ = 1`,
`ρ_c = 1`); [`rotating_king_moments`](@ref) gives its global moments.

# Example
```julia
m = solve_rotating_king(6.0, 0.6)
m.r_t_equator, m.r_t_pole       # equatorial and polar truncation radii [r_K]
rotating_king_moments(m).λ      # T_rot/T_kin; Einsel & Spurzem (1999) Table 1: 0.1981
```
"""
function solve_rotating_king(
    W0::Real,
    ω0::Real;
    l_max::Integer = 12,
    n_r::Integer = 600,
    n_μ::Integer = 48,
    r_max::Real = 0.0,
    relaxation::Real = 0.4,
    tol::Real = 1e-9,
    max_iter::Integer = 2000,
)::RotatingKingModel
    W₀ = Float64(W0)
    ω₀ = Float64(ω0)
    (isfinite(W₀) && W₀ > 0) ||
        throw(ArgumentError("solve_rotating_king: W0 must be positive and finite; got $W0"))
    (isfinite(ω₀) && ω₀ ≥ 0) ||
        throw(ArgumentError("solve_rotating_king: ω0 must be non-negative and finite; got $ω0"))
    (l_max ≥ 0 && iseven(l_max)) ||
        throw(ArgumentError("solve_rotating_king: l_max must be even and ≥ 0; got $l_max"))
    n_r ≥ 16 || throw(ArgumentError("solve_rotating_king: n_r must be ≥ 16; got $n_r"))
    n_μ > l_max || throw(
        ArgumentError("solve_rotating_king: n_μ must exceed l_max; got n_μ = $n_μ, l_max = $l_max"),
    )
    (isfinite(r_max) && r_max ≥ 0) || throw(
        ArgumentError("solve_rotating_king: r_max must be ≥ 0 (0 selects the default); got $r_max"),
    )
    0 < relaxation ≤ 1 ||
        throw(ArgumentError("solve_rotating_king: relaxation must lie in (0, 1]; got $relaxation"))
    tol > 0 || throw(ArgumentError("solve_rotating_king: tol must be > 0; got $tol"))
    max_iter ≥ 1 || throw(ArgumentError("solve_rotating_king: max_iter must be ≥ 1; got $max_iter"))

    α = Float64(relaxation)
    nr = Int(n_r)
    nμ = Int(n_μ)
    rk, wk, ρk = _solve_king(W₀)
    r̂_max = r_max > 0 ? Float64(r_max) : 1.5 * rk[end]
    r = collect(range(0.0, r̂_max; length = nr))
    μ, wμ = _rk_gauss_legendre(nμ)
    xq, wq = _rk_gauss_legendre(_RK_SPEED_NODES)
    ρ̂₀ = _rk_density(W₀, 0.0, xq, wq)

    ψ = _rk_initial_potential(r, rk, wk, ρk, nμ)
    ψh = zeros(Float64, nr, nμ)
    ψ0 = zeros(Float64, nr)
    ρ = zeros(Float64, nr, nμ)
    residual = Inf
    iterations = 0
    converged = false
    for iteration in 1:max_iter
        _rk_density_grid!(ρ, ψ, r, μ, ω₀, ρ̂₀, xq, wq)
        Φ = _rk_multipole_potential(r, μ, wμ, ρ, l_max)
        all(isfinite, Φ) || throw(_rk_no_equilibrium(ω₀, "non-finite potential"))
        # Non-spherical part: the potential less its angular average, under-relaxed
        @inbounds for i in 1:nr
            Φ̄ = 0.0
            for j in 1:nμ
                Φ̄ += wμ[j] * Φ[i, j]
            end
            Φ̄ *= 0.5
            for j in 1:nμ
                ψh[i, j] += α * (-(Φ[i, j] - Φ̄) - ψh[i, j])
            end
        end
        # Spherical part: exact outward integration at the current ψh
        _rk_spherical_part!(ψ0, r, μ, wμ, ψh, ω₀, ρ̂₀, xq, wq, W₀)
        residual = 0.0
        @inbounds for j in 1:nμ, i in 1:nr
            ψ_new = ψ0[i] + ψh[i, j]
            isfinite(ψ_new) || throw(_rk_no_equilibrium(ω₀, "non-finite potential"))
            residual = max(residual, abs(ψ_new - ψ[i, j]))
            ψ[i, j] = ψ_new
        end
        iterations = iteration
        if iteration ≥ _RK_CLOSURE_ITERATIONS && any(≥(0.0), @view ψ[nr, :])
            throw(
                _rk_no_equilibrium(ω₀, "the surface ψ = 0 does not close inside r_max = $(r̂_max)"),
            )
        end
        if residual < tol
            converged = true
            break
        end
    end
    converged ||
        throw(_rk_no_equilibrium(ω₀, "no convergence within max_iter = $(max_iter) iterations"))

    _rk_density_grid!(ρ, ψ, r, μ, ω₀, ρ̂₀, xq, wq)
    all(iszero, @view ρ[nr, :]) ||
        throw(_rk_no_equilibrium(ω₀, "non-zero density at r_max = $(r̂_max)"))

    c = _rk_coefficients(ψ, μ, wμ, l_max)
    r_t_equator = _rk_surface_radius(c, r, l_max, 0.0)
    r_t_pole = _rk_surface_radius(c, r, l_max, 1.0)
    (isfinite(r_t_equator) && isfinite(r_t_pole)) ||
        throw(_rk_no_equilibrium(ω₀, "the surface ψ = 0 does not close inside r_max = $(r̂_max)"))

    h = r[2] - r[1]
    profile = [2π * r[i]^2 * sum(wμ[j] * ρ[i, j] for j in 1:nμ) for i in 1:nr]
    I = zeros(Float64, nr - 1)
    _rk_interval_integrals!(I, profile, h, 1.0)
    cumulative = zeros(Float64, nr)
    for i in 1:(nr - 1)
        cumulative[i + 1] = cumulative[i] + I[i]
    end
    mass = cumulative[end]
    i_h = findfirst(≥(0.5 * mass), cumulative)
    r_h = NaN
    if i_h !== nothing && i_h > 1
        r_h =
            r[i_h - 1] +
            h * (0.5 * mass - cumulative[i_h - 1]) / (cumulative[i_h] - cumulative[i_h - 1])
    end

    return RotatingKingModel(
        W₀,
        ω₀,
        r,
        μ,
        ψ,
        ρ,
        Int(l_max),
        r_t_equator,
        r_t_pole,
        mass,
        r_h,
        iterations,
        residual,
    )
end

"""
    rotating_king_moments(m::RotatingKingModel)
        -> NamedTuple{(:λ, :t_rot, :e_dyn, :σ_c, :r_c, :r_t_over_r_c, :r_h_over_r_c, :r_t_over_r_K,
                       :r_h_over_r_K, :v_rot_rms_over_σ1d, :T_rot, :T_kin, :W)}

Global moments of the rotating King model `m` (Einsel & Spurzem 1999, DOI
10.1046/j.1365-8711.1999.02083.x), as volume integrals over the grid with
`dV̂ = 2π r̂² dr̂ dμ` (sixth-order radial quadrature, Gauss weights in `μ`), in
units of `ρ_c r_K³ σ²` for the energies:

- `T_kin = ½∫ρ̂⟨u²⟩ dV̂` (total kinetic energy, rotation included);
- `T_rot = ½∫ρ̂⟨v_φ⟩² dV̂` (ordered rotational energy);
- `W = ½∫ρ̂ Φ̂ dV̂` with the self-consistent multipole potential `Φ̂` (zero at infinity);
- `λ = T_rot/T_kin`; `t_rot = T_rot/|W|` (Ostriker & Peebles 1973, DOI 10.1086/152513);
- `e_dyn = 1 − r_t_pole/r_t_equator`, the ellipticity of the surface `ψ = 0`;
- `σ_c`, the central one-dimensional velocity dispersion in units of `σ`
  (`⟨u²⟩/3` at `ψ = W₀`, `k = 0`), and `r_c = r_K σ_c/σ`, King's core radius
  formed with the central dispersion (the `r_c` of Einsel & Spurzem 1999 and
  of King 1966, Eq. 9), in units of `r_K`;
- `r_t_over_r_c = r_t_equator/r_c` and `r_h_over_r_c = r_h/r_c` (the table
  convention of Einsel & Spurzem 1999), with `r_t_over_r_K = r_t_equator` and
  `r_h_over_r_K = r_h` in King radii;
- `v_rot_rms_over_σ1d = √(2T_rot/M̂) / √(2T_kin/(3M̂))`.

For a converged model `2T_kin + W = 0` to the discretisation accuracy (scalar
virial theorem).

# Example
```julia
mom = rotating_king_moments(solve_rotating_king(6.0, 0.4))
mom.λ, mom.t_rot, mom.e_dyn
```
"""
function rotating_king_moments(m::RotatingKingModel)
    n_r, n_μ = size(m.ψ)
    h = m.r[2] - m.r[1]
    _, wμ = _rk_gauss_legendre(n_μ)
    xq, wq = _rk_gauss_legendre(_RK_SPEED_NODES)
    ρ̂₀ = _rk_density(m.W0, 0.0, xq, wq)
    Φ = _rk_multipole_potential(m.r, m.μ, wμ, m.ρ, m.l_max)
    f_kin = zeros(Float64, n_r)
    f_rot = zeros(Float64, n_r)
    f_pot = zeros(Float64, n_r)
    for j in 1:n_μ
        sinθ = sqrt(max(1 - m.μ[j]^2, 0.0))
        for i in 1:n_r
            ρ = m.ρ[i, j]
            f_pot[i] += wμ[j] * ρ * Φ[i, j]
            ψ = m.ψ[i, j]
            (ψ > 0 && ρ > 0) || continue
            mo = _rk_moments(ψ, m.ω0 * m.r[i] * sinθ, xq, wq)
            f_kin[i] += wμ[j] * mo.ρu² / ρ̂₀
            f_rot[i] += wμ[j] * (mo.ρv_φ / ρ̂₀)^2 / ρ
        end
    end
    for i in 1:n_r
        dV = 2π * m.r[i]^2
        f_kin[i] *= dV
        f_rot[i] *= dV
        f_pot[i] *= dV
    end
    T_kin = 0.5 * _rk_radial_integral(f_kin, h)
    T_rot = 0.5 * _rk_radial_integral(f_rot, h)
    W = 0.5 * _rk_radial_integral(f_pot, h)
    M = m.mass
    v_rot_rms = sqrt(2 * T_rot / M)
    σ_1d = sqrt(2 * T_kin / (3 * M))
    mo_c = _rk_moments(m.W0, 0.0, xq, wq)
    σ_c = sqrt(mo_c.ρu² / mo_c.ρ / 3)
    return (
        λ = T_rot / T_kin,
        t_rot = T_rot / abs(W),
        e_dyn = 1 - m.r_t_pole / m.r_t_equator,
        σ_c = σ_c,
        r_c = σ_c,
        r_t_over_r_c = m.r_t_equator / σ_c,
        r_h_over_r_c = m.r_h / σ_c,
        r_t_over_r_K = m.r_t_equator,
        r_h_over_r_K = m.r_h,
        v_rot_rms_over_σ1d = v_rot_rms / σ_1d,
        T_rot = T_rot,
        T_kin = T_kin,
        W = W,
    )
end

# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------

"""
    _rotating_king_velocity(ψ, k, rng) -> (v_R, v_φ, v_z)

Exact draw from the local velocity distribution `∝ (e^{ψ − u²/2} − 1) e^{k u x}`
(`u < √(2ψ)`, `x = v_φ/u`, `k = ω₀R̂`), in units of `σ`, cylindrical components.

The speed `u` follows `p(u) ∝ u² (e^{ψ − u²/2} − 1) sinh(ku)/(ku)`: a King
speed is drawn by rejection from the box `[0, √(2ψ)] × [0, h_max]` (`h_max`
from a 200-point scan, raised by 0.1 % so that it bounds the maximum between
scan points), then accepted with probability
`sinh(ku)/(ku) / [sinh(k√(2ψ))/(k√(2ψ))] ≤ 1`. The direction cosine follows
`p(x | u) ∝ e^{kux}` on `[−1, 1]` by inversion,
`x = 1 + ln(ξ + (1 − ξ) e^{−2ku})/(ku)` (`x = 2ξ − 1` for `ku < 10⁻⁶`); the
perpendicular speed `u√(1 − x²)` is split between `v_R` and `v_z` with a
uniform angle.
"""
function _rotating_king_velocity(ψ::Real, k::Real, rng::AbstractRNG)
    ψf = Float64(ψ)
    kf = Float64(k)
    (isfinite(ψf) && ψf > 0) ||
        throw(ArgumentError("_rotating_king_velocity: ψ must be positive and finite; got $ψ"))
    (isfinite(kf) && kf ≥ 0) ||
        throw(ArgumentError("_rotating_king_velocity: k must be non-negative and finite; got $k"))
    a = sqrt(2ψf)
    h_max = 0.0
    for q in 1:_RK_SCAN_POINTS
        u = a * (q - 1) / (_RK_SCAN_POINTS - 1)
        h_max = max(h_max, u^2 * expm1(ψf - 0.5 * u^2))
    end
    h_max *= _RK_SCAN_MARGIN
    s_max = _rk_sinhc(kf * a)
    u = 0.0
    while true
        u = a * rand(rng)
        rand(rng) * h_max < u^2 * expm1(ψf - 0.5 * u^2) || continue
        rand(rng) * s_max < _rk_sinhc(kf * u) && break
    end
    z = kf * u
    ξ = rand(rng)
    x = z < 1e-6 ? 2ξ - 1 : 1 + log(ξ + (1 - ξ) * exp(-2z)) / z
    x = clamp(x, -1.0, 1.0)
    v_perp = u * sqrt(max(1 - x^2, 0.0))
    χ = 2π * rand(rng)
    return (v_perp * cos(χ), u * x, v_perp * sin(χ))
end

"""
    _rk_sample_μ!(cumulative, μ_ext, ρ_ext, i, f, ξ) -> Float64

Inverse-CDF draw of `μ = cos θ` from the density interpolated in radius,
`(1 − f) ρ_ext[i, :] + f ρ_ext[i + 1, :]`, taken as piecewise linear in `μ` on
the extended nodes `μ_ext = (−1, μ_1, …, μ_n, 1)`; `ξ` is the uniform variate
and `cumulative` a buffer of length `length(μ_ext) − 1`. Returns `NaN` when the
interpolated density vanishes.
"""
function _rk_sample_μ!(
    cumulative::AbstractVector{Float64},
    μ_ext::AbstractVector{Float64},
    ρ_ext::AbstractMatrix{Float64},
    i::Int,
    f::Float64,
    ξ::Float64,
)::Float64
    n = length(μ_ext)
    density(s) = (1 - f) * ρ_ext[i, s] + f * ρ_ext[i + 1, s]
    total = 0.0
    d_prev = density(1)
    for s in 1:(n - 1)
        d_next = density(s + 1)
        total += 0.5 * (μ_ext[s + 1] - μ_ext[s]) * (d_prev + d_next)
        cumulative[s] = total
        d_prev = d_next
    end
    total > 0 || return NaN
    target = ξ * total
    s = min(searchsortedlast(cumulative, target) + 1, n - 1)
    base = s == 1 ? 0.0 : cumulative[s - 1]
    q = max(target - base, 0.0)
    d_a = density(s)
    L = μ_ext[s + 1] - μ_ext[s]
    slope = (density(s + 1) - d_a) / L
    denominator = d_a + sqrt(max(d_a^2 + 2 * slope * q, 0.0))
    t = denominator > 0 ? 2q / denominator : 0.0
    return μ_ext[s] + clamp(t, 0.0, L)
end

"""
    sample_rotating_king(N, model::RotatingKingModel; rng = Random.default_rng()) -> (pos, vel)

Draw `N` particles from the rotating King model `model` (Einsel & Spurzem
1999, DOI 10.1046/j.1365-8711.1999.02083.x; see [`solve_rotating_king`](@ref)).

Positions: the radius by inverse CDF of the enclosed mass
`M(r̂) = ∫₀^r̂ 2π s² Σ_j w_j ρ̂(s, μ_j) ds` (trapezoidal, linear interpolation
on the grid); `μ = cos θ` from the conditional density at that radius
(interpolated between the bracketing grid rows, piecewise linear in `μ`,
exact inversion); azimuth uniform. A position at which the interpolated `ψ`
(`_rk_psi_at`) is not positive, possible only in the grid cell that contains
the truncation surface, is redrawn. Velocities: exact draws from the local
distribution function at `(ψ, k = ω₀R̂)` by `_rotating_king_velocity`, rotated
into Cartesian components with the position's azimuth.

Returns `(pos, vel)`, each `3 × N`, in model units (`r_K = 1`, `σ = 1`),
centred on the model centre; the mean velocity is not subtracted. The model
rotates about `+z` (`L_z > 0`, prograde).

# Example
```julia
using Random
m = solve_rotating_king(6.0, 0.6)
pos, vel = sample_rotating_king(10_000, m; rng = Xoshiro(1))
```
"""
function sample_rotating_king(
    N::Integer,
    model::RotatingKingModel;
    rng::AbstractRNG = Random.default_rng(),
)
    N ≥ 1 || throw(ArgumentError("sample_rotating_king: N must be ≥ 1; got $N"))
    r = model.r
    μ = model.μ
    ρ = model.ρ
    n_r, n_μ = size(ρ)
    h = r[2] - r[1]
    _, wμ = _rk_gauss_legendre(n_μ)
    xq, wq = _rk_gauss_legendre(_RK_SPEED_NODES)
    ρ̂₀ = _rk_density(model.W0, 0.0, xq, wq)
    c = _rk_coefficients(model.ψ, μ, wμ, model.l_max)

    μ_ext = vcat(-1.0, μ, 1.0)
    ρ_ext = zeros(Float64, n_r, n_μ + 2)
    for i in 1:n_r
        ψ_pole = _rk_psi_at(model, c, r[i], 1.0)
        ρ_pole = ψ_pole > 0 ? _rk_density(ψ_pole, 0.0, xq, wq) / ρ̂₀ : 0.0
        ρ_ext[i, 1] = ρ_pole
        ρ_ext[i, n_μ + 2] = ρ_pole
        for j in 1:n_μ
            ρ_ext[i, j + 1] = ρ[i, j]
        end
    end

    cdf = zeros(Float64, n_r)
    profile_prev = 0.0
    for i in 2:n_r
        profile = 2π * r[i]^2 * sum(wμ[j] * ρ[i, j] for j in 1:n_μ)
        cdf[i] = cdf[i - 1] + 0.5 * h * (profile_prev + profile)
        profile_prev = profile
    end
    cdf[end] > 0 || throw(ArgumentError("sample_rotating_king: the model carries no mass"))
    cdf ./= cdf[end]

    cumulative = zeros(Float64, n_μ + 1)
    pos = zeros(Float64, 3, N)
    vel = zeros(Float64, 3, N)
    for p in 1:N
        r_p = 0.0
        μ_p = 0.0
        ψ_p = 0.0
        while true
            ξ = rand(rng)
            i = clamp(searchsortedlast(cdf, ξ), 1, n_r - 1)
            Δc = cdf[i + 1] - cdf[i]
            f = Δc > 0 ? clamp((ξ - cdf[i]) / Δc, 0.0, 1.0) : 0.0
            r_p = r[i] + f * h
            μ_p = _rk_sample_μ!(cumulative, μ_ext, ρ_ext, i, f, rand(rng))
            isnan(μ_p) && continue
            ψ_p = _rk_psi_at(model, c, r_p, μ_p)
            ψ_p > 0 && break
        end
        φ = 2π * rand(rng)
        sφ, cφ = sincos(φ)
        R = r_p * sqrt(max(1 - μ_p^2, 0.0))
        pos[1, p] = R * cφ
        pos[2, p] = R * sφ
        pos[3, p] = r_p * μ_p
        v_R, v_φ, v_z = _rotating_king_velocity(ψ_p, model.ω0 * R, rng)
        vel[1, p] = v_R * cφ - v_φ * sφ
        vel[2, p] = v_R * sφ + v_φ * cφ
        vel[3, p] = v_z
    end
    return pos, vel
end

"""
    sample_rotating_king(N, W0, ω0; rng = Random.default_rng(), kwargs...) -> (pos, vel)

Solve the rotating King model with central potential `W0` and dimensionless
angular velocity `ω0` ([`solve_rotating_king`](@ref), keyword arguments
`kwargs` passed through) and draw `N` particles from it
([`sample_rotating_king(N, model)`](@ref sample_rotating_king)). Units as there
(`r_K = 1`, `σ = 1`); reference Einsel & Spurzem (1999), DOI
10.1046/j.1365-8711.1999.02083.x.

# Example
```julia
using Random
pos, vel = sample_rotating_king(10_000, 6.0, 0.6; rng = Xoshiro(1))
```
"""
function sample_rotating_king(
    N::Integer,
    W0::Real,
    ω0::Real;
    rng::AbstractRNG = Random.default_rng(),
    kwargs...,
)
    N ≥ 1 || throw(ArgumentError("sample_rotating_king: N must be ≥ 1; got $N"))
    model = solve_rotating_king(W0, ω0; kwargs...)
    return sample_rotating_king(N, model; rng = rng)
end

# ---------------------------------------------------------------------------
# Model cache, orientation and realised rotation
# ---------------------------------------------------------------------------

# Solved models by (W0, ω0), default solver settings.
const _ROTATING_KING_CACHE = Dict{Tuple{Float64,Float64},RotatingKingModel}()

"""
    _rotating_king_cached(W0, ω0) -> RotatingKingModel

[`solve_rotating_king`](@ref)`(W0, ω0)` with the default solver settings,
stored in `_ROTATING_KING_CACHE` on the first call and returned from it
afterwards (a solve takes several seconds; the generator samples every cluster
of a configuration from the model of its profile). The cache is a plain `Dict`
and not thread-safe: concurrent first calls from several tasks must not occur.
"""
function _rotating_king_cached(W0::Real, ω0::Real)::RotatingKingModel
    key = (Float64(W0), Float64(ω0))
    return get!(() -> solve_rotating_king(key[1], key[2]), _ROTATING_KING_CACHE, key)
end

"""
    _rotate_to_axis!(pos, vel, axis) -> (pos, vel)

Rotate the columns of the `3 × N` matrices `pos` and `vel` in place by the
proper rotation that maps `ẑ` onto `â = axis/‖axis‖`. For `â_z ≥ 0` it is the
rotation about `ẑ × â` (Rodrigues' formula, `R = I + [v]ₓ + [v]ₓ²/(1 + ẑ·â)`
with `v = ẑ × â`; the identity for `â = ẑ`); for `â_z < 0` the rotation by π
about `x̂` followed by the rotation that maps `−ẑ` onto `â` by the same
formula, so that `â = −ẑ` gives the rotation by π about `x̂`. Angular momenta
`r × v` are carried from `ẑ` to `â`.
"""
function _rotate_to_axis!(pos::AbstractMatrix{Float64}, vel::AbstractMatrix{Float64}, axis)
    (size(pos, 1) == 3 && size(vel, 1) == 3 && size(pos, 2) == size(vel, 2)) || throw(
        DimensionMismatch(
            "_rotate_to_axis!: pos and vel must both be 3 × N; got $(size(pos)) and $(size(vel))",
        ),
    )
    length(axis) == 3 ||
        throw(ArgumentError("_rotate_to_axis!: axis must have three components; got $axis"))
    n = sqrt(sum(abs2, axis))
    (isfinite(n) && n > 0) ||
        throw(ArgumentError("_rotate_to_axis!: axis must be finite and non-zero; got $axis"))
    a = (Float64(axis[1]) / n, Float64(axis[2]) / n, Float64(axis[3]) / n)
    flip = a[3] < 0
    # v = u × â with u = ẑ (no flip) or u = −ẑ (after the rotation by π about x̂)
    v = flip ? (a[2], -a[1], 0.0) : (-a[2], a[1], 0.0)
    c = flip ? -a[3] : a[3]
    K = [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]
    R = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0] .+ K .+ (K * K) ./ (1 + c)
    if flip
        R = R * [1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 -1.0]
    end
    for M in (pos, vel)
        @inbounds for p in axes(M, 2)
            x, y, z = M[1, p], M[2, p], M[3, p]
            M[1, p] = R[1, 1] * x + R[1, 2] * y + R[1, 3] * z
            M[2, p] = R[2, 1] * x + R[2, 2] * y + R[2, 3] * z
            M[3, p] = R[3, 1] * x + R[3, 2] * y + R[3, 3] * z
        end
    end
    return pos, vel
end

"""
    rotation_moments(mass, pos, vel; axis = [0, 0, 1], n_bins = 20)
        -> NamedTuple{(:λ, :l_axis, :t_rot, :t_kin)}

Ordered rotation of a particle set about the unit vector `â = axis/‖axis‖`,
measured about the centre of mass `x̄` and in the frame of the mass-weighted
mean velocity `v̄`; `mass` has length `N`, `pos` and `vel` are `3 × N`.

- `l_axis = Σ m [(r − x̄) × (v − v̄)]·â / Σ m`, the specific angular momentum
  along `â` [length × velocity of the input];
- `t_kin = ½ Σ m |v − v̄|²`, the kinetic energy;
- `t_rot = ½ Σ_b M_b ⟨v_φ⟩_b²`, the ordered rotational energy, with the
  particles sorted by the cylindrical radius `R` about `â` through `x̄` and
  divided into `n_bins` bins of equal counts, `M_b` the bin mass and
  `⟨v_φ⟩_b` its mass-weighted mean azimuthal velocity (positive in the sense
  of `â`); the sample estimator of `T_rot = ½∫ρ⟨v_φ⟩² dV` of
  [`rotating_king_moments`](@ref), biased upwards by about `n_bins/(3N)` in `λ`;
- `λ = t_rot/t_kin`.

Particles on the axis (`R = 0`) carry `v_φ = 0`. Throws a
`DimensionMismatch` for inconsistent sizes and an `ArgumentError` for an
empty set, `n_bins ∉ [1, N]`, or an axis that is not a finite non-zero
three-vector.

# Example
```julia
m = solve_rotating_king(6.0, 0.6)
pos, vel = sample_rotating_king(20_000, m)
rotation_moments(ones(20_000), pos, vel).λ   # ≈ rotating_king_moments(m).λ
```
"""
function rotation_moments(
    mass::AbstractVector{<:Real},
    pos::AbstractMatrix{<:Real},
    vel::AbstractMatrix{<:Real};
    axis = [0.0, 0.0, 1.0],
    n_bins::Integer = 20,
)::NamedTuple{(:λ, :l_axis, :t_rot, :t_kin),NTuple{4,Float64}}
    N = length(mass)
    N ≥ 1 || throw(ArgumentError("rotation_moments: the particle set is empty"))
    (size(pos) == (3, N) && size(vel) == (3, N)) || throw(
        DimensionMismatch(
            "rotation_moments: pos and vel must be 3 × $N; got $(size(pos)) and $(size(vel))",
        ),
    )
    1 ≤ n_bins ≤ N ||
        throw(ArgumentError("rotation_moments: n_bins must lie in [1, $N]; got $n_bins"))
    length(axis) == 3 ||
        throw(ArgumentError("rotation_moments: axis must have three components; got $axis"))
    n_a = sqrt(sum(abs2, axis))
    (isfinite(n_a) && n_a > 0) ||
        throw(ArgumentError("rotation_moments: axis must be finite and non-zero; got $axis"))
    â = (Float64(axis[1]) / n_a, Float64(axis[2]) / n_a, Float64(axis[3]) / n_a)

    M = 0.0
    x̄ = zeros(Float64, 3)
    v̄ = zeros(Float64, 3)
    for p in 1:N
        m = Float64(mass[p])
        M += m
        for k in 1:3
            x̄[k] += m * pos[k, p]
            v̄[k] += m * vel[k, p]
        end
    end
    M > 0 || throw(ArgumentError("rotation_moments: the total mass must be > 0; got $M"))
    x̄ ./= M
    v̄ ./= M

    R = zeros(Float64, N)
    v_φ = zeros(Float64, N)
    l_sum = 0.0
    t_kin = 0.0
    for p in 1:N
        m = Float64(mass[p])
        x, y, z = pos[1, p] - x̄[1], pos[2, p] - x̄[2], pos[3, p] - x̄[3]
        u, v, w = vel[1, p] - v̄[1], vel[2, p] - v̄[2], vel[3, p] - v̄[3]
        # Specific angular momentum along â and the cylindrical decomposition
        l_p = â[1] * (y * w - z * v) + â[2] * (z * u - x * w) + â[3] * (x * v - y * u)
        l_sum += m * l_p
        t_kin += 0.5 * m * (u^2 + v^2 + w^2)
        s = â[1] * x + â[2] * y + â[3] * z
        R_p = sqrt(max((x - s * â[1])^2 + (y - s * â[2])^2 + (z - s * â[3])^2, 0.0))
        R[p] = R_p
        v_φ[p] = R_p > 0 ? l_p / R_p : 0.0
    end

    order = sortperm(R)
    t_rot = 0.0
    for b in 1:n_bins
        lo = div((b - 1) * N, n_bins) + 1
        hi = div(b * N, n_bins)
        M_b = 0.0
        p_b = 0.0
        for q in lo:hi
            p = order[q]
            m = Float64(mass[p])
            M_b += m
            p_b += m * v_φ[p]
        end
        M_b > 0 && (t_rot += 0.5 * p_b^2 / M_b)
    end
    λ = t_kin > 0 ? t_rot / t_kin : 0.0
    return (λ = λ, l_axis = l_sum / M, t_rot = t_rot, t_kin = t_kin)
end
