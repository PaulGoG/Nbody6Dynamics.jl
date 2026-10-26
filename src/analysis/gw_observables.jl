# =============================================================================
# Gravitational-wave inspiral time, masses and frequencies of a binary
# =============================================================================
#
# Orbit-averaged quadrupole radiation of a Newtonian Keplerian binary (Peters &
# Mathews 1963, DOI 10.1103/PhysRev.131.435; Peters 1964, DOI
# 10.1103/PhysRev.136.B1224). Masses enter in M☉, semi-major axes in au, times
# leave in yr and frequencies in Hz; the arithmetic is SI throughout.

using QuadGK: quadgk

# Newtonian constant of gravitation [m³ kg⁻¹ s⁻²] (CODATA 2018)
const _GW_G = 6.67430e-11
# Speed of light in vacuum [m s⁻¹]
const _GW_C = 2.99792458e8
# Solar mass [kg]: IAU 2015 nominal GM☉ over the CODATA 2018 G
const _GW_MSUN = 1.98840987e30
# Astronomical unit [m] (IAU 2012)
const _GW_AU = 1.495978707e11
# Julian year [s]
const _GW_YR = 3.15576e7
# Coefficient k = 121/304 of the eccentricity decay and of the Peters track
const _GW_K = 121 / 304

"""
    _gw_check_binary(m1, m2, a, e)

Throw an `ArgumentError` unless both masses and the semi-major axis are
positive and `0 ≤ e < 1`.
"""
function _gw_check_binary(m1::Real, m2::Real, a::Real, e::Real)
    (m1 > 0 && m2 > 0) || throw(ArgumentError("masses must be positive, got m1 = $m1, m2 = $m2"))
    a > 0 || throw(ArgumentError("semi-major axis must be positive, got a = $a"))
    0 ≤ e < 1 || throw(ArgumentError("eccentricity must lie in [0, 1), got e = $e"))
    return nothing
end

"""
    _gw_beta(m1_kg, m2_kg) -> Float64

β = (64/5) G³ m₁ m₂ M / c⁵ [m⁴ s⁻¹] of the circular decay a(t) = (a₀⁴ − 4βt)^{1/4}.
"""
_gw_beta(m1_kg::Float64, m2_kg::Float64)::Float64 =
    64 / 5 * _GW_G^3 * m1_kg * m2_kg * (m1_kg + m2_kg) / _GW_C^5

"""
    _gw_circular_time(m1, m2, a) -> Float64

Circular coalescence time T_c = 5c⁵a⁴/(256 G³ m₁ m₂ M) = a⁴/(4β) [s], masses
in M☉ and `a` in au.
"""
function _gw_circular_time(m1::Real, m2::Real, a::Real)::Float64
    m1_kg = Float64(m1) * _GW_MSUN
    m2_kg = Float64(m2) * _GW_MSUN
    a_m = Float64(a) * _GW_AU
    return 5 * _GW_C^5 * a_m^4 / (256 * _GW_G^3 * m1_kg * m2_kg * (m1_kg + m2_kg))
end

"""
    gw_enhancement(e) -> Float64

Peters–Mathews enhancement of the orbit-averaged GW luminosity of an
eccentric binary over the circular one with the same semi-major axis,

    F(e) = (1 + 73e²/24 + 37e⁴/96) / (1 − e²)^{7/2},

which is also the sum over harmonics Σₙ g(n, e) (Peters & Mathews 1963,
DOI 10.1103/PhysRev.131.435). Throws an `ArgumentError` unless `0 ≤ e < 1`.

# Example
```julia
gw_enhancement(0.5)   # ≈ 4.88
```
"""
function gw_enhancement(e::Real)::Float64
    0 ≤ e < 1 || throw(ArgumentError("eccentricity must lie in [0, 1), got e = $e"))
    ε = Float64(e)
    e² = ε^2
    return (1 + 73 * e² / 24 + 37 * e²^2 / 96) / ((1 - ε) * (1 + ε))^3.5
end

"""
    _gw_time_factor_x(e0) -> Float64

Time factor F(e₀) by the x-form quadrature (e = e₀x), suited to e₀ ≤ 0.5:

    F = (48/19)(1 − e₀²)^{1/2}(1 + ke₀²)^{−3480/2299}
        ∫₀¹ x^{29/19}(1 + ke₀²x²)^{1181/2299}(1 − e₀²x²)^{−3/2} dx,

k = 121/304. Valid for 0 ≤ e₀ < 1; e₀ = 0 returns 1.
"""
function _gw_time_factor_x(e0::Float64)::Float64
    e0 == 0.0 && return 1.0
    ke0² = _GW_K * e0^2
    integrand(x) = x^(29 / 19) * (1 + ke0² * x^2)^(1181 / 2299) / ((1 - e0 * x) * (1 + e0 * x))^1.5
    integral = first(quadgk(integrand, 0.0, 1.0; rtol = 1e-12))
    return 48 / 19 * sqrt((1 - e0) * (1 + e0)) * (1 + ke0²)^(-3480 / 2299) * integral
end

"""
    _gw_time_factor_tau(e0) -> Float64

Time factor F(e₀) by the τ-form quadrature (1 − e² = (1 − e₀²)e^τ), suited to
e₀ > 0.5:

    F = (24/19) e₀^{−48/19}(1 + ke₀²)^{−3480/2299}
        ∫₀^{τ_max} e(τ)^{10/19}(1 + ke(τ)²)^{1181/2299} e^{−τ/2} dτ,

e(τ)² = 1 − (1 − e₀²)e^τ, τ_max = −ln(1 − e₀²), k = 121/304. The integrand is
bounded and tends to F → 768/425 as e₀ → 1. Valid for 0 < e₀ < 1.
"""
function _gw_time_factor_tau(e0::Float64)::Float64
    τ_max = -log((1 - e0) * (1 + e0))
    # e(τ)² = 1 − e^{τ − τ_max}, free of cancellation near τ_max
    function integrand(τ)
        e² = max(-expm1(τ - τ_max), 0.0)
        return e²^(5 / 19) * (1 + _GW_K * e²)^(1181 / 2299) * exp(-τ / 2)
    end
    integral = first(quadgk(integrand, 0.0, τ_max; rtol = 1e-12))
    return 24 / 19 * e0^(-48 / 19) * (1 + _GW_K * e0^2)^(-3480 / 2299) * integral
end

"""
    gw_time_factor(e0) -> Float64

Ratio F(e₀) ≡ T/[T_c(1 − e₀²)^{7/2}] of the exact Peters (1964) coalescence time
T to the circular time T_c of the same semi-major axis, scaled by
(1 − e₀²)^{7/2}. It rises from 1 at e₀ = 0 (minimum ≈ 0.98 near e₀ = 0.5) to
768/425 as e₀ → 1, the high-eccentricity limit of Peters (1964,
DOI 10.1103/PhysRev.136.B1224).

Evaluated by singularity-free quadratures of the Peters integral (adaptive
Gauss–Kronrod, relative tolerance 10⁻¹², k = 121/304): the x-form
(e = e₀x) for e₀ ≤ 0.5 and the τ-form (1 − e² = (1 − e₀²)e^τ) for e₀ > 0.5.
Throws an `ArgumentError` unless `0 ≤ e0 < 1`.

# Example
```julia
gw_time_factor(0.9)   # ≈ 1.139561
```
"""
function gw_time_factor(e0::Real)::Float64
    0 ≤ e0 < 1 || throw(ArgumentError("eccentricity must lie in [0, 1), got e0 = $e0"))
    ε = Float64(e0)
    return ε ≤ 0.5 ? _gw_time_factor_x(ε) : _gw_time_factor_tau(ε)
end

"""
    gw_inspiral_time(m1, m2, a, e) -> Float64

Coalescence time [yr] of a binary of masses `m1`, `m2` [M☉], semi-major axis
`a` [au] and eccentricity `e` under orbit-averaged quadrupole radiation,

    T = T_c (1 − e²)^{7/2} F(e),   T_c = 5c⁵a⁴ / (256 G³ m₁ m₂ M),

with F from [`gw_time_factor`](@ref) (Peters 1964, DOI 10.1103/PhysRev.136.B1224).
No post-Newtonian correction is applied; [`gw_zwick_factor`](@ref) gives the
1PN-recalibrated factor separately. Throws an `ArgumentError` for non-positive
masses or `a`, or `e` outside [0, 1).

# Example
```julia
gw_inspiral_time(30.0, 30.0, 1.0, 0.999)   # ≈ 3.6×10³ yr
```
"""
function gw_inspiral_time(m1::Real, m2::Real, a::Real, e::Real)::Float64
    _gw_check_binary(m1, m2, a, e)
    ε = Float64(e)
    t_c = _gw_circular_time(m1, m2, a)
    return t_c * ((1 - ε) * (1 + ε))^3.5 * gw_time_factor(ε) / _GW_YR
end

"""
    gw_inspiral_time_mandel(m1, m2, a, e) -> Float64

Coalescence time [yr] by the fit of Mandel (2021, DOI 10.3847/2515-5172/ac2d35),

    T ≈ T_c (1 + 0.27e¹⁰ + 0.33e²⁰ + 0.2e¹⁰⁰⁰)(1 − e²)^{7/2},

accurate to better than 3 % for 0 ≤ e ≤ 0.99999 against
[`gw_inspiral_time`](@ref); units and argument checks as there.
"""
function gw_inspiral_time_mandel(m1::Real, m2::Real, a::Real, e::Real)::Float64
    _gw_check_binary(m1, m2, a, e)
    ε = Float64(e)
    t_c = _gw_circular_time(m1, m2, a)
    fit = 1 + 0.27 * ε^10 + 0.33 * ε^20 + 0.2 * ε^1000
    return t_c * fit * ((1 - ε) * (1 + ε))^3.5 / _GW_YR
end

"""
    gw_zwick_factor(m1, m2, a, e) -> Float64

1PN-recalibrated correction factor Q_f = exp(2.5 r_S/p₀) of the Peters
coalescence time (Zwick et al. 2020, DOI 10.1093/mnras/staa1314), with
r_S = 2GM/c² the Schwarzschild radius of the total mass and p₀ = a(1 − e) the
pericentre distance. Masses in M☉, `a` in au. Zwick et al. fitted Q_f to
inspirals evolved with the 3.5PN flux and state it as accurate within 1PN
precision; they give no validity range, the uncertainties of the fitted
coefficients overlap only for p₀ ≲ 35 r_S, and the factor is meaningless near
the ISCO. Identifying their initial (a₀, e₀) with the Newtonian elements of the
engine is ambiguous at O(r_S/p₀). The factor is therefore a sensitivity
measure reported on its own; [`gw_inspiral_time`](@ref) never applies it.
"""
function gw_zwick_factor(m1::Real, m2::Real, a::Real, e::Real)::Float64
    _gw_check_binary(m1, m2, a, e)
    r_s = 2 * _GW_G * (Float64(m1) + Float64(m2)) * _GW_MSUN / _GW_C^2
    p₀ = Float64(a) * _GW_AU * (1 - Float64(e))
    return exp(2.5 * r_s / p₀)
end

"""
    chirp_mass(m1, m2) -> Float64

Chirp mass 𝓜 = (m₁m₂)^{3/5}/M^{1/5} [M☉] of component masses `m1`, `m2`
[M☉], M = m₁ + m₂. Throws an `ArgumentError` for non-positive masses.

# Example
```julia
chirp_mass(30.0, 30.0)   # ≈ 26.12
```
"""
function chirp_mass(m1::Real, m2::Real)::Float64
    (m1 > 0 && m2 > 0) || throw(ArgumentError("masses must be positive, got m1 = $m1, m2 = $m2"))
    μ1 = Float64(m1)
    μ2 = Float64(m2)
    return (μ1 * μ2)^0.6 / (μ1 + μ2)^0.2
end

"""
    orbital_frequency(m1, m2, a) -> Float64

Keplerian orbital frequency f_orb = (2π)⁻¹ √(GM/a³) [Hz] of masses `m1`, `m2`
[M☉] at semi-major axis `a` [au]. Throws an `ArgumentError` for non-positive
masses or `a`.
"""
function orbital_frequency(m1::Real, m2::Real, a::Real)::Float64
    _gw_check_binary(m1, m2, a, 0.0)
    gm = _GW_G * (Float64(m1) + Float64(m2)) * _GW_MSUN
    a_m = Float64(a) * _GW_AU
    return sqrt(gm / a_m^3) / (2π)
end

"""
    _gw_frequency_si(gm, a_m, e, definition) -> Float64

GW frequency [Hz] of the named `definition` for GM [m³ s⁻²] and semi-major
axis `a_m` [m]; see `gw_frequency`. Throws an `ArgumentError` for an
unknown definition.
"""
function _gw_frequency_si(gm::Float64, a_m::Float64, e::Float64, definition::Symbol)::Float64
    if definition === :f22
        return sqrt(gm / a_m^3) / π
    elseif definition === :wen
        return sqrt(gm) / π * (1 + e)^1.1954 / (a_m * (1 - e) * (1 + e))^1.5
    elseif definition === :peri
        return sqrt(gm / (a_m * (1 - e))^3) / π
    end
    throw(ArgumentError("unknown GW frequency definition :$definition; one of :f22, :wen, :peri"))
end

"""
    gw_frequency(m1, m2, a, e; definition = :f22) -> Float64

GW frequency [Hz] of a binary of masses `m1`, `m2` [M☉], semi-major axis `a`
[au] and eccentricity `e`, under one of the definitions in use:

- `:f22`: dominant quadrupole harmonic, f₂₂ = 2f_orb;
- `:wen`: peak harmonic of Wen (2003, DOI 10.1086/378794),
  f = (√(GM)/π)(1 + e)^{1.1954}/[a(1 − e²)]^{1.5};
- `:peri`: circular orbit of the pericentre radius, f = π⁻¹ √(GM/r_p³),
  r_p = a(1 − e).

All three coincide at e = 0. The definitions differ by large factors at high
e, so a product that stores the frequency names the definition in its column.
Throws an `ArgumentError` for an unknown definition or invalid orbit.

# Example
```julia
gw_frequency(30.0, 30.0, 1.0, 0.999; definition = :wen)   # ≈ 12.6×10⁻³ Hz
```
"""
function gw_frequency(m1::Real, m2::Real, a::Real, e::Real; definition::Symbol = :f22)::Float64
    _gw_check_binary(m1, m2, a, e)
    gm = _GW_G * (Float64(m1) + Float64(m2)) * _GW_MSUN
    return _gw_frequency_si(gm, Float64(a) * _GW_AU, Float64(e), definition)
end

"""
    _gw_track_semi_major_axis(c₀, e) -> Float64

Semi-major axis a(e) = c₀ e^{12/19}(1 + ke²)^{870/2299}/(1 − e²) of the Peters
track (Peters 1964, Eq. 5.11), in the unit of `c₀`.
"""
_gw_track_semi_major_axis(c₀::Float64, e::Float64)::Float64 =
    c₀ * e^(12 / 19) * (1 + _GW_K * e^2)^(870 / 2299) / ((1 - e) * (1 + e))

"""
    eccentricity_at_frequency(m1, m2, a0, e0, f_ref; definition = :f22)
        -> NamedTuple{(:e, :a_au, :time_to_ref_yr, :above_at_record)}

Eccentricity of a binary recorded at (`a0` [au], `e0`), masses `m1`, `m2`
[M☉], when its GW frequency of the given `definition` (see
[`gw_frequency`](@ref)) has risen to `f_ref` [Hz] under the Peters (1964)
evolution.

Along the track a(e) = c₀e^{12/19}(1 + ke²)^{870/2299}/(1 − e²), k = 121/304,
with c₀ fixed by (a₀, e₀), each definition's frequency decreases monotonically
with e, so the root of f(e) = f_ref on (0, e₀] is unique; it is found by
bisection in ln e to a relative precision of 10⁻¹² in e, starting from
e₀·10⁻¹² as the lower limit (extended by up to three further factors of
10⁻¹² when needed). For e₀ = 0 the track is circular and e stays 0.

Returns `e`, the semi-major axis `a_au` [au] there, the time `time_to_ref_yr`
[yr] = T(a₀, e₀) − T(a_ref, e_ref) of [`gw_inspiral_time`](@ref), and
`above_at_record`; a binary already at or above `f_ref` at the record returns
`NaN` for the three values and `above_at_record = true`. Throws an
`ArgumentError` for an invalid orbit, a non-positive `f_ref`, an unknown
definition, or a root below the extended lower limit.

# Example
```julia
eccentricity_at_frequency(30.0, 30.0, 1.0, 0.999, 10.0; definition = :wen).e   # ≈ 4.42×10⁻⁴
```
"""
function eccentricity_at_frequency(
    m1::Real,
    m2::Real,
    a0::Real,
    e0::Real,
    f_ref::Real;
    definition::Symbol = :f22,
)::@NamedTuple{e::Float64, a_au::Float64, time_to_ref_yr::Float64, above_at_record::Bool}
    _gw_check_binary(m1, m2, a0, e0)
    f_ref > 0 || throw(ArgumentError("reference frequency must be positive, got f_ref = $f_ref"))
    gm = _GW_G * (Float64(m1) + Float64(m2)) * _GW_MSUN
    a0_m = Float64(a0) * _GW_AU
    ε₀ = Float64(e0)
    f = Float64(f_ref)
    # Validates the definition before any branch
    f_record = _gw_frequency_si(gm, a0_m, ε₀, definition)
    above = (e = NaN, a_au = NaN, time_to_ref_yr = NaN, above_at_record = true)
    f_record ≥ f && return above
    t_record = gw_inspiral_time(m1, m2, a0, ε₀)

    if ε₀ == 0.0
        # Circular track: every definition reduces to f = π⁻¹ √(GM/a³)
        a_ref_au = cbrt(gm / (π * f)^2) / _GW_AU
        t_ref = gw_inspiral_time(m1, m2, a_ref_au, 0.0)
        return (
            e = 0.0,
            a_au = a_ref_au,
            time_to_ref_yr = t_record - t_ref,
            above_at_record = false,
        )
    end

    c₀ = a0_m * (1 - ε₀) * (1 + ε₀) / (ε₀^(12 / 19) * (1 + _GW_K * ε₀^2)^(870 / 2299))
    f_track(ε) = _gw_frequency_si(gm, _gw_track_semi_major_axis(c₀, ε), ε, definition)

    e_lo = ε₀ * 1e-12
    n_extensions = 0
    while f_track(e_lo) < f
        n_extensions == 3 && throw(
            ArgumentError(
                "f_ref = $f_ref Hz is not reached for e ≥ $e_lo along the Peters track of " *
                "a0 = $a0 au, e0 = $e0 (definition :$definition)",
            ),
        )
        e_lo *= 1e-12
        n_extensions += 1
    end

    # f decreases with e: f(e_lo) ≥ f_ref > f(e₀)
    lo = log(e_lo)
    hi = log(ε₀)
    while hi - lo > 1e-12
        mid = (lo + hi) / 2
        if f_track(exp(mid)) ≥ f
            lo = mid
        else
            hi = mid
        end
    end
    e_ref = exp((lo + hi) / 2)
    a_ref_au = _gw_track_semi_major_axis(c₀, e_ref) / _GW_AU
    t_ref = gw_inspiral_time(m1, m2, a_ref_au, e_ref)
    return (e = e_ref, a_au = a_ref_au, time_to_ref_yr = t_record - t_ref, above_at_record = false)
end

"""
    _gw_peters_rhs!(du, u, β, t)

Orbit-averaged Peters equations in u = (ln a, e), SI units:

    d ln a/dt = −β F(e)/a⁴,
    de/dt     = −(19/12) β e (1 + ke²) / [a⁴ (1 − e²)^{5/2}],

with β = (64/5) G³m₁m₂M/c⁵, F the Peters–Mathews enhancement and k = 121/304.
"""
function _gw_peters_rhs!(du, u, β::Float64, t)
    inv_a4 = exp(-4 * u[1])
    ε = u[2]
    e² = ε^2
    one_minus_e² = (1 - ε) * (1 + ε)
    du[1] = -β * inv_a4 * (1 + 73 * e² / 24 + 37 * e²^2 / 96) / one_minus_e²^3.5
    du[2] = -19 / 12 * β * inv_a4 * ε * (1 + _GW_K * e²) / one_minus_e²^2.5
    return nothing
end

"""
    _gw_inspiral_time_ode(m1, m2, a, e0) -> Float64

Coalescence time [yr] by direct integration of the Peters equations
(`_gw_peters_rhs!`) from (`a` [au], `e0`) with Tsit5 at relative and
absolute tolerance 10⁻¹², terminated when a < 10⁻³ a₀; the remaining time of
the end state is added from `gw_inspiral_time`. An independent check
of the quadrature, used by the tests.
"""
function _gw_inspiral_time_ode(m1::Real, m2::Real, a::Real, e0::Real)::Float64
    _gw_check_binary(m1, m2, a, e0)
    β = _gw_beta(Float64(m1) * _GW_MSUN, Float64(m2) * _GW_MSUN)
    a_m = Float64(a) * _GW_AU
    t_c = a_m^4 / (4β)
    ln_a_stop = log(1.0e-3 * a_m)
    cb = ContinuousCallback((u, t, integrator) -> u[1] - ln_a_stop, terminate!)
    # T < (768/425) T_c for every e₀, so 2 T_c brackets the stop
    prob = ODEProblem(_gw_peters_rhs!, [log(a_m), Float64(e0)], (0.0, 2.0 * t_c), β)
    sol = solve(prob, Tsit5(); callback = cb, reltol = 1e-12, abstol = 1e-12, maxiters = 10^7)
    u_end = sol.u[end]
    # The event is located to the root-finding tolerance; the end state's own inspiral time is
    # added below, so only a stop at the end of the time span is an error.
    u_end[1] ≤ ln_a_stop + 1e-3 ||
        error("Peters integration ended at a = $(exp(u_end[1]) / _GW_AU) au before 10⁻³ a₀")
    a_end_au = exp(u_end[1]) / _GW_AU
    return sol.t[end] / _GW_YR + gw_inspiral_time(m1, m2, a_end_au, max(u_end[2], 0.0))
end

# =============================================================================
# Strain amplitudes
# =============================================================================
#
# Orientation-averaged strains of the harmonics of a Keplerian binary (Peters &
# Mathews 1963; Barack & Cutler 2004, DOI 10.1103/PhysRevD.69.082005; Wagg,
# Broekgaarden & Gültekin 2022, DOI 10.3847/1538-4365/ac5c52) and the 22-mode
# amplitude of the LISA Data Challenge galactic binaries. Distances enter in
# kpc; every strain scales as 1/D.

# Kiloparsec [m] (IAU 2015 parsec)
const _GW_KPC = 3.0856775814913673e19

"""Throw an `ArgumentError` unless the distance [kpc] is positive."""
function _gw_check_distance(distance_kpc::Real)
    distance_kpc > 0 ||
        throw(ArgumentError("distance must be positive, got distance_kpc = $distance_kpc"))
    return nothing
end

"""Throw an `ArgumentError` unless the harmonic number is ≥ 1."""
function _gw_check_harmonic(n::Integer)
    n ≥ 1 || throw(ArgumentError("harmonic number must be ≥ 1, got n = $n"))
    return nothing
end

"""Bessel function J_k(x) of any integer order, J₋ₖ = (−1)ᵏ Jₖ."""
function _gw_besselj(k::Int, x::Float64)::Float64
    k ≥ 0 && return besselj(k, x)
    j = besselj(-k, x)
    return isodd(k) ? -j : j
end

"""
    peters_mathews_g(n::Integer, e::Real) -> Float64

Relative GW power g(n, e) radiated in the n-th harmonic of the orbital
frequency by a Keplerian binary of eccentricity `e` (Peters & Mathews 1963,
DOI 10.1103/PhysRev.131.435, Eq. 20). With Jₖ = Jₖ(ne) the Bessel functions of
the first kind,

    g(n, e) = (n⁴/32) { [J_{n−2} − 2e J_{n−1} + (2/n) Jₙ + 2e J_{n+1} − J_{n+2}]²
                        + (1 − e²) [J_{n−2} − 2Jₙ + J_{n+2}]² + (4/(3n²)) Jₙ² }.

A circular orbit radiates only in n = 2, g(2, 0) = 1 and g(n, 0) = 0 for n ≠ 2,
and the harmonics sum to the enhancement, Σ_{n≥1} g(n, e) =
[`gw_enhancement`](@ref)(e). Throws an `ArgumentError` for `n < 1` or `e`
outside [0, 1).

# Example
```julia
peters_mathews_g(2, 0.0)   # 1.0
```
"""
function peters_mathews_g(n::Integer, e::Real)::Float64
    _gw_check_harmonic(n)
    0 ≤ e < 1 || throw(ArgumentError("eccentricity must lie in [0, 1), got e = $e"))
    k = Int(n)
    ε = Float64(e)
    x = k * ε
    j₋₂ = _gw_besselj(k - 2, x)
    j₋₁ = _gw_besselj(k - 1, x)
    j₀ = _gw_besselj(k, x)
    j₊₁ = _gw_besselj(k + 1, x)
    j₊₂ = _gw_besselj(k + 2, x)
    ν = Float64(k)
    b₁ = j₋₂ - 2ε * j₋₁ + 2 / ν * j₀ + 2ε * j₊₁ - j₊₂
    b₂ = j₋₂ - 2j₀ + j₊₂
    return ν^4 / 32 * (b₁^2 + (1 - ε) * (1 + ε) * b₂^2 + 4 / (3ν^2) * j₀^2)
end

"""
    gw_amplitude_ldc(mchirp_msun::Real, f_gw_hz::Real; distance_kpc::Real = 1.0) -> Float64

Strain amplitude of the leading (l = m = 2) mode of a circular binary of chirp
mass `mchirp_msun` [M☉] at GW frequency `f_gw_hz` [Hz] and distance
`distance_kpc` [kpc], in the convention of the LISA Data Challenge galactic
binaries (LDC manual, Eqs. 30–32a and 45–46),

    𝒜 = 2 (G𝓜)^{5/3} (πf)^{2/3} / (c⁴ D).

𝒜 is half the face-on amplitude h₀ = 4(G𝓜)^{5/3}(πf)^{2/3}/(c⁴D) of the
quadrupole formula and carries no orientation average. The amplitude scales as
1/D: the catalogues store the distance-free product 𝒜 × D at D = 1 kpc, and
the value at a distance D is the stored one divided by D in kpc. Throws an
`ArgumentError` for a non-positive chirp mass, frequency or distance.

# Example
```julia
gw_amplitude_ldc(0.5, 1e-3)   # ≈ 1.87×10⁻²² at 1 kpc
```
"""
function gw_amplitude_ldc(mchirp_msun::Real, f_gw_hz::Real; distance_kpc::Real = 1.0)::Float64
    mchirp_msun > 0 ||
        throw(ArgumentError("chirp mass must be positive, got mchirp_msun = $mchirp_msun"))
    f_gw_hz > 0 || throw(ArgumentError("GW frequency must be positive, got f_gw_hz = $f_gw_hz"))
    _gw_check_distance(distance_kpc)
    gm = _GW_G * Float64(mchirp_msun) * _GW_MSUN
    d_m = Float64(distance_kpc) * _GW_KPC
    return 2 * gm^(5 / 3) * (π * Float64(f_gw_hz))^(2 / 3) / (_GW_C^4 * d_m)
end

"""
    gw_harmonic_strain(m1::Real, m2::Real, a_au::Real, e::Real, n::Integer;
                       distance_kpc::Real = 1.0) -> Float64

Orientation- and polarisation-averaged strain amplitude hₙ of the n-th
harmonic (frequency fₙ = n f_orb) of a binary of masses `m1`, `m2` [M☉],
semi-major axis `a_au` [au] and eccentricity `e` at distance `distance_kpc`
[kpc] (Wagg, Broekgaarden & Gültekin 2022, LEGWORK, DOI
10.3847/1538-4365/ac5c52, Eq. 27),

    hₙ² = (2^{28/3}/5) (G𝓜)^{10/3} / (c⁸ D²) · g(n, e)/n² · (π f_orb)^{4/3},

with f_orb = (2π)⁻¹ √(GM/a³), 𝓜 the [`chirp_mass`](@ref) and g from
[`peters_mathews_g`](@ref). A circular binary radiates only in n = 2, where
h₂ = (4/√5) [`gw_amplitude_ldc`](@ref)(𝓜, 2f_orb) = (4/5)^{1/2} h₀, the face-on
amplitude h₀ = 2𝒜 times the root of the orientation average
⟨[(1 + cos²ι)/2]² + cos²ι⟩ = 4/5. The strain scales as 1/D: the
catalogues store hₙ × D at D = 1 kpc, and the value at a distance D is the
stored one divided by D in kpc. Throws an `ArgumentError` for non-positive
masses, `a_au` or distance, `e` outside [0, 1) or `n < 1`.

# Example
```julia
gw_harmonic_strain(0.6, 0.6, 0.01, 0.0, 2; distance_kpc = 8.0)
```
"""
function gw_harmonic_strain(
    m1::Real,
    m2::Real,
    a_au::Real,
    e::Real,
    n::Integer;
    distance_kpc::Real = 1.0,
)::Float64
    _gw_check_binary(m1, m2, a_au, e)
    _gw_check_harmonic(n)
    _gw_check_distance(distance_kpc)
    f_orb = orbital_frequency(m1, m2, a_au)
    gm = _GW_G * chirp_mass(m1, m2) * _GW_MSUN
    ν = Float64(n)
    shape = 2^(28 / 3) / 5 * peters_mathews_g(n, e) / ν^2 * (π * f_orb)^(4 / 3)
    return sqrt(shape) * gm^(5 / 3) / (_GW_C^4 * Float64(distance_kpc) * _GW_KPC)
end

"""
    gw_harmonic_frequency_derivative(m1::Real, m2::Real, a_au::Real, e::Real, n::Integer)
        -> Float64

Orbit-averaged drift ḟₙ = n ḟ_orb [Hz s⁻¹] of the frequency fₙ = n f_orb of the
n-th harmonic under quadrupole radiation, for masses `m1`, `m2` [M☉],
semi-major axis `a_au` [au] and eccentricity `e` (Peters 1964,
DOI 10.1103/PhysRev.136.B1224),

    ḟ_orb = (48/5) π^{8/3} (G𝓜/c³)^{5/3} (2f_orb)^{11/3} F(e),

with F the [`gw_enhancement`](@ref); for n = 2 at e = 0 it is the circular
chirp ḟ = (96/5) π^{8/3} (G𝓜/c³)^{5/3} f^{11/3} of f = 2f_orb. Throws an
`ArgumentError` for non-positive masses or `a_au`, `e` outside [0, 1) or
`n < 1`.
"""
function gw_harmonic_frequency_derivative(
    m1::Real,
    m2::Real,
    a_au::Real,
    e::Real,
    n::Integer,
)::Float64
    _gw_check_binary(m1, m2, a_au, e)
    _gw_check_harmonic(n)
    f_orb = orbital_frequency(m1, m2, a_au)
    t_m = _GW_G * chirp_mass(m1, m2) * _GW_MSUN / _GW_C^3
    fdot_orb = 48 / 5 * π^(8 / 3) * t_m^(5 / 3) * (2 * f_orb)^(11 / 3) * gw_enhancement(e)
    return Float64(n) * fdot_orb
end

"""
    gw_characteristic_strain(m1::Real, m2::Real, a_au::Real, e::Real, n::Integer;
                             distance_kpc::Real = 1.0) -> Float64

Orientation-averaged characteristic strain h_{c,n} of the n-th harmonic of an
inspiralling binary of masses `m1`, `m2` [M☉], semi-major axis `a_au` [au] and
eccentricity `e` at distance `distance_kpc` [kpc] (Wagg, Broekgaarden &
Gültekin 2022, LEGWORK, DOI 10.3847/1538-4365/ac5c52, Eq. 26, after Barack &
Cutler 2004, DOI 10.1103/PhysRevD.69.082005),

    h_{c,n}² = (2^{5/3} / (3π^{4/3})) (G𝓜)^{5/3} / (c³ D²) · f_orb^{−1/3} · g(n, e) / (n F(e)),

with f_orb = (2π)⁻¹ √(GM/a³), g from [`peters_mathews_g`](@ref) and F the
[`gw_enhancement`](@ref). It equals h_{c,n}² = (fₙ²/ḟₙ) hₙ² of the
[`gw_harmonic_strain`](@ref) hₙ and the drift ḟₙ of
[`gw_harmonic_frequency_derivative`](@ref), fₙ = n f_orb. The strain scales as
1/D: the catalogues store h_{c,n} × D at D = 1 kpc, and the value at a
distance D is the stored one divided by D in kpc. Throws an `ArgumentError`
for non-positive masses, `a_au` or distance, `e` outside [0, 1) or `n < 1`.
"""
function gw_characteristic_strain(
    m1::Real,
    m2::Real,
    a_au::Real,
    e::Real,
    n::Integer;
    distance_kpc::Real = 1.0,
)::Float64
    _gw_check_binary(m1, m2, a_au, e)
    _gw_check_harmonic(n)
    _gw_check_distance(distance_kpc)
    f_orb = orbital_frequency(m1, m2, a_au)
    gm = _GW_G * chirp_mass(m1, m2) * _GW_MSUN
    ratio = peters_mathews_g(n, e) / (Float64(n) * gw_enhancement(e))
    hc² = 2^(5 / 3) / (3 * π^(4 / 3)) * gm^(5 / 3) / _GW_C^3 * f_orb^(-1 / 3) * ratio
    return sqrt(hc²) / (Float64(distance_kpc) * _GW_KPC)
end
