# =============================================================================
# Primordial binary populations for the merger initial conditions
# =============================================================================
# Stars are sampled from the IMF, a fraction of them are paired, the
# centres of mass of all systems are placed by the density sampler and
# virialised, and each pair is expanded into two bodies on a Keplerian
# orbit at the end. The engine reads primordial pairs as the first 2·NBIN0
# bodies of dat.10 (bodies 2i−1, 2i), so every cluster's pairs are written
# first, cluster by cluster, followed by every cluster's singles.

"""
    MassiveBinarySpec(; model = "none", m_min = 5.0, fraction = 0.69, pairing = "pool",
                        log_period_min = 0.15, log_period_max = 3.5, period_exponent = -0.55,
                        q_min = 0.1, q_exponent = -0.1, ecc_exponent = -0.45, ecc_max = 0.9,
                        ecc_turnover = true)

Binary population of the massive primaries of one cluster
(`[merger.clusterN.binaries.massive]`). The defaults are the O-star
distributions of Sana et al. (2012, DOI 10.1126/science.1223344):
`f(lP) ∝ lP^π`, `f(q) ∝ q^κ`, `f(e) ∝ e^η` with `lP = log10(P/d)`. π, κ and
`f_bin` are those of Sana et al. (2012); η and the domains of `lP`, `q` and `e`
are as tabulated by Cottaar & Hénault-Brunet (2014, DOI
10.1051/0004-6361/201322445, Table 1) and used by Sana et al. (2013, DOI
10.1051/0004-6361/201219621). `fraction` is the fraction of primaries with a
companion at `log(P/d) ≤ 3.5`: a wider `log_period_max` keeps the density in
`log P` only if `fraction` is raised accordingly. Sana et al. measured O stars
(`M ≳ 15 M☉`); below that mass the branch is an extrapolation (Moe & Di Stefano
2017, DOI 10.3847/1538-4365/aa6fb6).

# Fields
- `model`: `"none"` (no massive branch; every star belongs to the population of
  [`BinarySpec`](@ref)), `"sana2012"` or `"moe2017"` (the joint distribution
  `f(lP, q, e | M₁)` of Moe & Di Stefano 2017 with the inner-binary probability
  of Moe's code; it ignores `fraction`, `period_exponent`, `q_exponent`,
  `ecc_exponent` and `ecc_max` and uses `m_min`, `pairing`, `log_period_min`,
  `log_period_max`, `q_min` and `ecc_turnover`)
- `m_min`: primary-mass threshold [M☉] of the branch; `> 0` (`moe2017` warns
  below 0.8)
- `fraction`: binary fraction of the massive primaries, one Bernoulli trial per
  primary reached unpaired; `[0, 1]` (Sana: `f_bin = 0.69 ± 0.09`)
- `pairing`: `"pool"` (the partner is the unpaired star of mass `≤ m₁` closest to
  `q m₁`; the IMF multiset is kept) or `"replace"` (a random unpaired star below
  `m_min` takes the mass `q m₁`; this alters the IMF)
- `log_period_min`, `log_period_max`: period range [log10 d];
  `0 ≤ log_period_min < log_period_max ≤ 8.43` (Sana: 0.15 and 3.5, secondary);
  `moe2017` requires `log_period_min ≥ 0.15` and `log_period_max ≤ 8.0`
- `period_exponent`: π (Sana: −0.55 ± 0.22)
- `q_min`: lower mass-ratio bound; `0 < q_min < 1` (Sana: 0.1, secondary)
- `q_exponent`: κ (Sana: −0.10 ± 0.58)
- `ecc_exponent`: η; `> −1` (Sana: −0.45 ± 0.17, secondary)
- `ecc_max`: upper eccentricity bound; `0 < ecc_max < 1` (Sana: 0.9, secondary)
- `ecc_turnover`: eccentricity law of `moe2017`. `true` is the variant of Moe's
  code: `f(e) ∝ e^η` below `0.8 e_max`, carrying the probability `8/(9 + η)`,
  and a linear decline to zero at `e_max` above it. `false` is the paper-only law
  `f(e) ∝ e^η` on `[0, e_max]`
"""
Base.@kwdef struct MassiveBinarySpec
    model::String = "none"
    m_min::Float64 = 5.0
    fraction::Float64 = 0.69
    pairing::String = "pool"
    log_period_min::Float64 = 0.15
    log_period_max::Float64 = 3.5
    period_exponent::Float64 = -0.55
    q_min::Float64 = 0.1
    q_exponent::Float64 = -0.1
    ecc_exponent::Float64 = -0.45
    ecc_max::Float64 = 0.9
    ecc_turnover::Bool = true
end

"""
    BinarySpec(; fraction = 0.0, pairing = "random", period = "kroupa1995",
                 a_min = 0.01, a_max = 100.0, q_min = 0.1, eccentricity = "thermal",
                 roche_fill_max = 0.0, massive = MassiveBinarySpec())

Primordial binary population of one cluster (`[merger.clusterN.binaries]`).
With a massive branch (`massive.model ≠ "none"`) the keys `fraction` to
`eccentricity` describe the stars below `massive.m_min` only.

# Fields
- `fraction`: binary fraction by systems, `N_b / (N_s + N_b)`; `0 ≤ fraction < 1`
- `pairing`: `"random"` (both components drawn from the IMF; the more massive
  becomes the primary) or `"uniform_q"` (the secondary's mass is `q m₁` with
  `q` uniform in `[q_min, 1]`, replacing the drawn mass of the partner)
- `period`: `"kroupa1995"` (the birth period distribution of Kroupa 1995b,
  DOI 10.1093/mnras/277.4.1507, Eq. 8,
  `f(log P) ∝ (log P − 1) / (45 + (log P − 1)²)` for `1 ≤ log P/d ≤ 8.43`,
  converted to a semi-major axis with Kepler's third law) or `"loguniform"`
  (semi-major axis log-uniform between `a_min` and `a_max`). Kroupa's birth
  distribution belongs with his pre-main-sequence eigenevolution, which the
  generator does not apply: short eccentric orbits are not circularised and
  the mass ratios are not evolved (Belloni et al. 2017,
  DOI 10.1093/mnras/stx1763)
- `a_min`, `a_max`: semi-major axis bounds [AU] for `"loguniform"`; `0 < a_min < a_max`
- `q_min`: lower mass-ratio bound for `"uniform_q"`; `0 < q_min < 1`
- `eccentricity`: `"thermal"` (`f(e) = 2e`) or `"circular"`
- `roche_fill_max`: Roche-lobe guard at periastron for every pair, the largest
  admitted `φ_i = R_i / (r_L(q_i) a (1 − e))` of either component; `0` (no guard)
  or `0 < roche_fill_max ≤ 1`. `R_i` is the ZAMS radius of `zams_radius`
  (Tout et al. 1996, DOI 10.1093/mnras/281.1.257) and `r_L` the Roche-lobe radius
  of Eggleton (1983, DOI 10.1086/160960); the period and eccentricity laws are
  truncated, not rejected
- `massive`: the branch of the massive primaries, a `MassiveBinarySpec`
"""
Base.@kwdef struct BinarySpec
    fraction::Float64 = 0.0
    pairing::String = "random"
    period::String = "kroupa1995"
    a_min::Float64 = 0.01
    a_max::Float64 = 100.0
    q_min::Float64 = 0.1
    eccentricity::String = "thermal"
    roche_fill_max::Float64 = 0.0
    massive::MassiveBinarySpec = MassiveBinarySpec()
end

const _AU_PC = 1 / 206264.806
const _AU_RSUN = 215.032
const _DAY_YR = 1 / 365.25
# Kroupa (1995b) birth period distribution parameters and range (log10 P in days)
const _KROUPA95_ETA = 2.5
const _KROUPA95_DELTA = 45.0
const _KROUPA95_LOGP_MIN = 1.0
const _KROUPA95_LOGP_MAX = 8.43

"""Fail-fast bounds of a [`BinarySpec`](@ref); `key` names the owning cluster table. Raises an `ArgumentError` naming the offending key."""
function _validate_binaries(b::BinarySpec, key::String)
    (0 ≤ b.fraction < 1) ||
        _config_error("$key.binaries.fraction", "must satisfy 0 ≤ fraction < 1; got $(b.fraction)")
    b.pairing in ("random", "uniform_q") || _config_error(
        "$key.binaries.pairing",
        "must be \"random\" or \"uniform_q\"; got \"$(b.pairing)\"",
    )
    b.period in ("kroupa1995", "loguniform") || _config_error(
        "$key.binaries.period",
        "must be \"kroupa1995\" or \"loguniform\"; got \"$(b.period)\"",
    )
    b.eccentricity in ("thermal", "circular") || _config_error(
        "$key.binaries.eccentricity",
        "must be \"thermal\" or \"circular\"; got \"$(b.eccentricity)\"",
    )
    (0 < b.a_min < b.a_max) || _config_error(
        "$key.binaries",
        "must satisfy 0 < a_min < a_max [AU]; got $(b.a_min), $(b.a_max)",
    )
    (0 < b.q_min < 1) ||
        _config_error("$key.binaries.q_min", "must satisfy 0 < q_min < 1; got $(b.q_min)")
    return nothing
end

"""
    _sample_log_period_kroupa1995(rng) -> Float64

`log10 P` [days] from the Kroupa (1995) birth period distribution by
inverse transform: with `x = log P − 1`, `F(x) = (η/2) ln((δ + x²)/δ)`.
"""
function _sample_log_period_kroupa1995(rng::AbstractRNG)::Float64
    x_max = _KROUPA95_LOGP_MAX - _KROUPA95_LOGP_MIN
    f_max = 0.5 * _KROUPA95_ETA * log((_KROUPA95_DELTA + x_max^2) / _KROUPA95_DELTA)
    u = rand(rng) * f_max
    x = sqrt(_KROUPA95_DELTA * (exp(2u / _KROUPA95_ETA) - 1))
    return _KROUPA95_LOGP_MIN + x
end

"""
    _sample_log_period_kroupa1995(rng, log_period_lo) -> Float64

The Kroupa (1995) law truncated to `log P ≥ log_period_lo` [days]: the uniform
variate of the inverse transform is restricted to `[F(x_lo), F(x_max)]`.
"""
function _sample_log_period_kroupa1995(rng::AbstractRNG, log_period_lo::Float64)::Float64
    x_max = _KROUPA95_LOGP_MAX - _KROUPA95_LOGP_MIN
    f_max = 0.5 * _KROUPA95_ETA * log((_KROUPA95_DELTA + x_max^2) / _KROUPA95_DELTA)
    x_lo = log_period_lo - _KROUPA95_LOGP_MIN
    f_lo = x_lo > 0 ? 0.5 * _KROUPA95_ETA * log((_KROUPA95_DELTA + x_lo^2) / _KROUPA95_DELTA) : 0.0
    u = f_lo + rand(rng) * (f_max - f_lo)
    x = sqrt(_KROUPA95_DELTA * (exp(2u / _KROUPA95_ETA) - 1))
    return _KROUPA95_LOGP_MIN + x
end

"""
    _power_law_quantile(u, x₀, x₁, k) -> Float64

Inverse CDF at `u ∈ [0, 1]` of `p(x) ∝ x^k` on `[x₀, x₁]`:
`[x₀^{k+1} + u (x₁^{k+1} − x₀^{k+1})]^{1/(k+1)}`, and `x₀ (x₁/x₀)^u` for `k = −1`.
"""
function _power_law_quantile(u::Float64, x₀::Float64, x₁::Float64, k::Float64)::Float64
    x = if k == -1
        x₀ * (x₁ / x₀)^u
    else
        s = k + 1
        (x₀^s + u * (x₁^s - x₀^s))^(1 / s)
    end
    return clamp(x, x₀, x₁)
end

"""
    _sana_draw(spec::MassiveBinarySpec, n::Integer, rng) -> (lP, q, e)

`n` independent triples from the truncated power laws of `spec`: `lP` [log10 d]
on `[log_period_min, log_period_max]` with `period_exponent`, `q` on `[q_min, 1]`
with `q_exponent`, `e` on `[0, ecc_max]` with `ecc_exponent`, drawn in the order
`lP, q, e` per binary.
"""
function _sana_draw(spec::MassiveBinarySpec, n::Integer, rng::AbstractRNG)
    lP = Vector{Float64}(undef, n)
    q = Vector{Float64}(undef, n)
    e = Vector{Float64}(undef, n)
    for k in 1:n
        lP[k] = _power_law_quantile(
            rand(rng),
            spec.log_period_min,
            spec.log_period_max,
            spec.period_exponent,
        )
        q[k] = _power_law_quantile(rand(rng), spec.q_min, 1.0, spec.q_exponent)
        e[k] = _power_law_quantile(rand(rng), 0.0, spec.ecc_max, spec.ecc_exponent)
    end
    return lP, q, e
end

"""Roche-lobe radius in units of the separation, Eggleton (1983), `q = m_i/m_other`."""
_roche_lobe_fraction(q::Float64)::Float64 = 0.49 * cbrt(q)^2 / (0.6 * cbrt(q)^2 + log(1 + cbrt(q)))

"""
    _roche_separation_au(m₁, m₂, zmet, φ) -> Float64

Smallest periastron distance [AU] at which neither ZAMS component fills more
than the fraction `φ` of its Roche lobe: `max_i R_i / (φ r_L(q_i))`.
"""
function _roche_separation_au(m₁::Float64, m₂::Float64, zmet::Float64, φ::Float64)::Float64
    a_rsun = max(
        zams_radius(m₁, zmet) / (φ * _roche_lobe_fraction(m₁ / m₂)),
        zams_radius(m₂, zmet) / (φ * _roche_lobe_fraction(m₂ / m₁)),
    )
    return a_rsun / _AU_RSUN
end

# Kepler III in AU, yr, M☉
_semi_major_axis_au(log_period::Float64, M::Float64)::Float64 =
    cbrt(M * (exp10(log_period) * _DAY_YR)^2)
_log_period_days(a_au::Float64, M::Float64)::Float64 = log10(sqrt(a_au^3 / M) / _DAY_YR)

_guard_range_error(m₁::Float64, m₂::Float64, a_rl::Float64, φ::Float64, range::String) = throw(
    ArgumentError(
        "sample_binaries: the pair m₁ = $m₁, m₂ = $m₂ M☉ needs a ≥ $a_rl AU for " *
        "roche_fill_max = $φ, beyond the period range ($range)",
    ),
)

"""
    _guarded_orbit(spec::BinarySpec, m₁, m₂, a_rl, rng) -> (a_au, e, truncated)

Semi-major axis [AU] and eccentricity of a pair of the low-mass population with
the period (or axis) law truncated to `a ≥ a_rl` and the eccentricity law to
`e ≤ 1 − a_rl/a`; `truncated` flags a pair whose period (or axis) range the guard
raised.
"""
function _guarded_orbit(
    spec::BinarySpec,
    m₁::Float64,
    m₂::Float64,
    a_rl::Float64,
    rng::AbstractRNG,
)::Tuple{Float64,Float64,Bool}
    M = m₁ + m₂
    truncated = false
    a_au = if spec.period == "kroupa1995"
        lP_rl = _log_period_days(a_rl, M)
        lP_rl > _KROUPA95_LOGP_MAX && _guard_range_error(
            m₁,
            m₂,
            a_rl,
            spec.roche_fill_max,
            "log(P/d) ≤ $_KROUPA95_LOGP_MAX",
        )
        truncated = lP_rl > _KROUPA95_LOGP_MIN
        _semi_major_axis_au(_sample_log_period_kroupa1995(rng, lP_rl), M)
    else
        a_rl > spec.a_max &&
            _guard_range_error(m₁, m₂, a_rl, spec.roche_fill_max, "a ≤ $(spec.a_max) AU")
        truncated = a_rl > spec.a_min
        a_lo = max(spec.a_min, a_rl)
        a_lo * (spec.a_max / a_lo)^rand(rng)
    end
    a_au = max(a_au, a_rl)
    e = 0.0
    if spec.eccentricity == "thermal"
        e_hi = clamp(1 - a_rl / a_au, 0.0, 1.0)
        e = e_hi * sqrt(rand(rng))
    end
    return a_au, e, truncated
end

"""
    _sana_orbit(spec::MassiveBinarySpec, u_P, u_e, m₁, m₂, a_rl, φ) -> (a_au, e, truncated)

Semi-major axis [AU] and eccentricity of a massive pair from its uniform
variates `u_P` (period) and `u_e` (eccentricity); `φ` is the fill limit that
set `a_rl` (for the error message). For `a_rl > 0` the period law
is truncated to `a ≥ a_rl` and the eccentricity law to `e ≤ min(ecc_max, 1 − a_rl/a)`;
`truncated` flags a pair whose period range the guard raised.
"""
function _sana_orbit(
    spec::MassiveBinarySpec,
    u_P::Float64,
    u_e::Float64,
    m₁::Float64,
    m₂::Float64,
    a_rl::Float64,
    φ::Float64,
)::Tuple{Float64,Float64,Bool}
    M = m₁ + m₂
    truncated = false
    lP_lo = spec.log_period_min
    if a_rl > 0
        lP_rl = _log_period_days(a_rl, M)
        lP_rl > spec.log_period_max &&
            _guard_range_error(m₁, m₂, a_rl, φ, "log(P/d) ≤ $(spec.log_period_max)")
        if lP_rl > lP_lo
            lP_lo = lP_rl
            truncated = true
        end
    end
    lP = _power_law_quantile(u_P, lP_lo, spec.log_period_max, spec.period_exponent)
    a_au = _semi_major_axis_au(lP, M)
    e_hi = spec.ecc_max
    if a_rl > 0
        a_au = max(a_au, a_rl)
        e_rl = max(1 - a_rl / a_au, 0.0)
        if e_rl < e_hi
            e_hi = e_rl
        end
    end
    return a_au, _power_law_quantile(u_e, 0.0, e_hi, spec.ecc_exponent), truncated
end

# -----------------------------------------------------------------------------
# Moe & Di Stefano (2017, DOI 10.3847/1538-4365/aa6fb6): f(lP, q, e | M₁)
# -----------------------------------------------------------------------------
# lP = log10(P/d), L = log10(M₁/M☉). The fits hold for 0.8 < M₁ < 40 M☉; the
# inner-binary probability, the clamp of M₁ and the eccentricity turnover are
# conventions of Moe's code (as ported in COSMIC), not published fits.

"""Primary mass [M☉] clamped to `[0.8, 40]`, the range of the Moe & Di Stefano (2017) statistics (convention of Moe's code)."""
_moe_clamp_m1(m₁::Float64)::Float64 = clamp(m₁, 0.8, 40.0)

"""
    _moe_companion_frequency(lP, m₁) -> Float64

Companions with `q > 0.3` per decade of period, `f_logP;q>0.3(M₁, P)`, Moe &
Di Stefano (2017, DOI 10.3847/1538-4365/aa6fb6, Eqs. 20–23) with `α = 0.018`,
`Δ = 0.7`. The first segment (`f₁`) also covers `lP < 0.2`, as in Moe's code.
"""
function _moe_companion_frequency(lP::Float64, m₁::Float64)::Float64
    L = log10(m₁)
    f₁ = 0.020 + 0.04L + 0.07L^2
    f₂₇ = 0.039 + 0.07L + 0.01L^2
    f₅₅ = 0.078 - 0.05L + 0.04L^2
    α, Δ = 0.018, 0.7
    lP < 1.0 && return f₁
    lP < 2.7 - Δ && return f₁ + (lP - 1) / (1.7 - Δ) * (f₂₇ - f₁ - α * Δ)
    lP < 2.7 + Δ && return f₂₇ + α * (lP - 2.7)
    lP < 5.5 && return f₂₇ + α * Δ + (lP - 2.7 - Δ) / (2.8 - Δ) * (f₅₅ - f₂₇ - α * Δ)
    return f₅₅ * exp(-0.3 * (lP - 5.5))
end

# Linear interpolation in log M₁ between the mass anchors 1.2, 3.5 and 6 M☉ (Moe's code)
function _moe_mass_anchors(m₁::Float64, y₁₂::Float64, y₃₅::Float64, y₆::Float64)::Float64
    m₁ ≤ 1.2 && return y₁₂
    m₁ ≤ 3.5 && return y₁₂ + log(m₁ / 1.2) / log(3.5 / 1.2) * (y₃₅ - y₁₂)
    m₁ ≤ 6.0 && return y₃₅ + log(m₁ / 3.5) / log(6.0 / 3.5) * (y₆ - y₃₅)
    return y₆
end

"""
    _moe_q_parameters(lP, m₁) -> (; γ_small, γ_large, F_twin)

Mass-ratio law of Moe & Di Stefano (2017, DOI 10.3847/1538-4365/aa6fb6, §9.1):
`γ_smallq` on `0.1 < q < 0.3` (Eqs. 13–15), `γ_largeq` on `0.3 < q < 1`
(Eqs. 9–11) and the twin excess `F_twin` above `q = 0.95` (Eqs. 5–7). The
anchors at 0.8–1.2, 3.5 and > 6 M☉ are interpolated linearly in `log M₁`.
"""
function _moe_q_parameters(
    lP::Float64,
    m₁::Float64,
)::NamedTuple{(:γ_small, :γ_large, :F_twin),NTuple{3,Float64}}
    γl₉ = lP < 5 ? -0.5 : -0.5 - 0.3 * (lP - 5)
    γl₁₀ =
        lP < 1 ? -0.5 : lP < 4.5 ? -0.5 - 0.2 * (lP - 1) : lP < 6.5 ? -1.2 - 0.4 * (lP - 4.5) : -2.0
    γl₁₁ = lP < 1 ? -0.5 : lP < 2 ? -0.5 - 0.9 * (lP - 1) : lP < 4 ? -1.4 - 0.3 * (lP - 2) : -2.0
    γs₁₃ = 0.3
    γs₁₄ = lP < 2.5 ? 0.2 : lP < 5.5 ? 0.2 - 0.3 * (lP - 2.5) : -0.7 - 0.2 * (lP - 5.5)
    γs₁₅ = lP < 1 ? 0.1 : lP < 3 ? 0.1 - 0.15 * (lP - 1) : lP < 5.6 ? -0.2 - 0.5 * (lP - 3) : -1.5
    F₁ = 0.30 - 0.15 * log10(m₁)
    lP_twin = m₁ ≤ 6.5 ? 8.0 - m₁ : 1.5
    F = lP < 1 ? F₁ : lP < lP_twin ? F₁ * (1 - (lP - 1) / (lP_twin - 1)) : 0.0
    return (
        γ_small = _moe_mass_anchors(m₁, γs₁₃, γs₁₄, γs₁₅),
        γ_large = _moe_mass_anchors(m₁, γl₉, γl₁₀, γl₁₁),
        F_twin = F,
    )
end

# ∫_{x₀}^{x₁} x^k dx
_power_integral(x₀::Float64, x₁::Float64, k::Float64)::Float64 =
    k == -1 ? log(x₁ / x₀) : (x₁^(k + 1) - x₀^(k + 1)) / (k + 1)

"""
    _moe_q_weights(p, q_lo) -> (w_small, w_twin, w_large)

Integrated weights of the three segments of the Moe & Di Stefano (2017) mass-ratio
law above `q_lo ≥ 0.1`, normalised so that the `q > 0.3` part (power law plus
twin excess) integrates to 1 at `q_lo ≤ 0.3`: the density is
`(1 − F) q^γl / I_l + (F/0.05) 1[q ≥ 0.95]` on `[0.3, 1]` and continues as
`p_l(0.3) (q/0.3)^γs` below 0.3. At `q_lo = 0.1`, `w_small = W_s`.
"""
function _moe_q_weights(
    p::NamedTuple{(:γ_small, :γ_large, :F_twin),NTuple{3,Float64}},
    q_lo::Float64,
)::NTuple{3,Float64}
    γl, γs, F = p.γ_large, p.γ_small, p.F_twin
    I_l = _power_integral(0.3, 1.0, γl)
    p_l03 = (1 - F) * 0.3^γl / I_l
    w_small = q_lo < 0.3 ? p_l03 * 0.3^(-γs) * _power_integral(q_lo, 0.3, γs) : 0.0
    w_twin = F * (1 - max(q_lo, 0.95)) / 0.05
    w_large = (1 - F) * _power_integral(max(q_lo, 0.3), 1.0, γl) / I_l
    return w_small, w_twin, w_large
end

"""
`1 + W_s`, the factor converting `f_logP;q>0.3` into `f_logP;q>0.1` at `(lP, M₁)`
(Moe & Di Stefano 2017, DOI 10.3847/1538-4365/aa6fb6, §9.1; COSMIC's `q_factor`).
"""
_moe_q_factor(lP::Float64, m₁::Float64)::Float64 =
    sum(_moe_q_weights(_moe_q_parameters(lP, m₁), 0.1))

"""
Power-law index η of `f(e) ∝ e^η`, Moe & Di Stefano (2017, DOI
10.3847/1538-4365/aa6fb6, Eqs. 17–18), interpolated in `log M₁` between 3 and
7 M☉ and frozen below `lP = 0.7` at its value there (Moe's code).
"""
function _moe_eta(lP::Float64, m₁::Float64)::Float64
    x = max(lP, 0.7)
    η₃ = 0.6 - 0.7 / (x - 0.5)
    η₇ = 0.9 - 0.2 / (x - 0.5)
    m₁ ≤ 3 && return η₃
    m₁ ≥ 7 && return η₇
    return η₃ + log(m₁ / 3) / log(7 / 3) * (η₇ - η₃)
end

"""
`e_max(P) = 1 − (P/2 d)^{−2/3}`, Moe & Di Stefano (2017, DOI
10.3847/1538-4365/aa6fb6, Eq. 3); `0` (circular) for `P ≤ 2 d`, at most 0.99
(the upper end of the eccentricity grid of Moe's code).
"""
_moe_e_max(lP::Float64)::Float64 = clamp(1 - (exp10(lP) / 2)^(-2 / 3), 0.0, 0.99)

"""
Inner-binary probability of a companion at `lP`, the approximation of Moe's code
(COSMIC port of Moe & Di Stefano 2017, DOI 10.3847/1538-4365/aa6fb6):
`1` for `lP ≤ 1.5`, else `max(0, 1 − 0.11 (lP − 1.5)^{1.43} (M₁/10 M☉)^{0.56})`.
"""
_moe_probbin(lP::Float64, m₁::Float64)::Float64 =
    lP ≤ 1.5 ? 1.0 : max(0.0, 1 - 0.11 * (lP - 1.5)^1.43 * (m₁ / 10)^0.56)

"""
    _moe_period_table(m₁, lP_lo, lP_hi) -> (lP, CumP)

Grid `lP` [log10 d] of step 0.005 from `lP_lo` to `lP_hi` (both included; the
last step may be shorter) and the cumulative trapezoid integral
`CumP(lP) = ∫ f_logP;q>0.3 (1 + W_s) probbin dlP` of Moe & Di Stefano (2017,
DOI 10.3847/1538-4365/aa6fb6) with the inner-binary probability of Moe's code.
`CumP[end]` is the inner-binary fraction `f_bin(M₁)` over `[lP_lo, lP_hi]`.
"""
function _moe_period_table(
    m₁::Float64,
    lP_lo::Float64,
    lP_hi::Float64,
)::Tuple{Vector{Float64},Vector{Float64}}
    n = floor(Int, (lP_hi - lP_lo) / 0.005 + 1e-9)
    lP = [lP_lo + 0.005k for k in 0:n]
    if lP_hi - lP[end] > 1e-9
        push!(lP, lP_hi)
    else
        lP[end] = lP_hi
    end
    cum = zeros(length(lP))
    g_prev =
        _moe_companion_frequency(lP[1], m₁) * _moe_q_factor(lP[1], m₁) * _moe_probbin(lP[1], m₁)
    for k in 2:length(lP)
        g = _moe_companion_frequency(lP[k], m₁) * _moe_q_factor(lP[k], m₁) * _moe_probbin(lP[k], m₁)
        cum[k] = cum[k - 1] + 0.5 * (g_prev + g) * (lP[k] - lP[k - 1])
        g_prev = g
    end
    return lP, cum
end

# Inverse of the tabulated cumulative by linear interpolation (flat stretches skipped)
function _moe_invert_cum(lP::Vector{Float64}, cum::Vector{Float64}, target::Float64)::Float64
    k = searchsortedfirst(cum, target)
    k ≤ 1 && return lP[1]
    k > length(cum) && return lP[end]
    return lP[k - 1] + (target - cum[k - 1]) / (cum[k] - cum[k - 1]) * (lP[k] - lP[k - 1])
end

# Tabulated cumulative at x by linear interpolation
function _moe_cum_at(lP::Vector{Float64}, cum::Vector{Float64}, x::Float64)::Float64
    x ≤ lP[1] && return cum[1]
    x ≥ lP[end] && return cum[end]
    k = searchsortedlast(lP, x)
    return cum[k] + (x - lP[k]) / (lP[k + 1] - lP[k]) * (cum[k + 1] - cum[k])
end

"""
    _moe_draw_q(u₁, u₂, p, q_min) -> Float64

Mass ratio from two uniforms (Moe & Di Stefano 2017, DOI
10.3847/1538-4365/aa6fb6): `u₁` selects the segment (power law `γ_small` on
`0.1–0.3`, twin excess uniform on `0.95–1`, power law `γ_large` on `0.3–1`) by
its weight, `u₂` places `q` in it (logarithmic form for an exponent of −1). The
law is restricted to `q ≥ max(q_min, 0.1)` with the segment weights renormalised
to the admitted range; `q_min ≥ 0.3` drops the small segment.
"""
function _moe_draw_q(
    u₁::Float64,
    u₂::Float64,
    p::NamedTuple{(:γ_small, :γ_large, :F_twin),NTuple{3,Float64}},
    q_min::Float64,
)::Float64
    q_lo = max(q_min, 0.1)
    w_s, w_t, w_l = _moe_q_weights(p, q_lo)
    x = u₁ * (w_s + w_t + w_l)
    x < w_s && return _power_law_quantile(u₂, q_lo, 0.3, p.γ_small)
    if x < w_s + w_t
        q_t = max(q_lo, 0.95)
        return q_t + (1 - q_t) * u₂
    end
    return _power_law_quantile(u₂, max(q_lo, 0.3), 1.0, p.γ_large)
end

"""
    _moe_draw_e(u, η, e_max, turnover) -> Float64

Eccentricity at the quantile `u` of the Moe & Di Stefano (2017, DOI
10.3847/1538-4365/aa6fb6) law. With `turnover` (Moe's code), `f(e) ∝ e^η` below
`0.8 e_max` carries the probability `p = 8/(9 + η)`: `u < p` gives
`e = 0.8 e_max (u/p)^{1/(η+1)}`, otherwise `e = e_max − 0.2 e_max √(1 − u′)`
with `u′ = (u − p)/(1 − p)` (linear decline). Without it, `e = e_max u^{1/(η+1)}`.
Returns 0 for `η ≤ −1` (not normalisable) or `e_max = 0`.
"""
function _moe_draw_e(u::Float64, η::Float64, e_max::Float64, turnover::Bool)::Float64
    (η ≤ -1 || e_max ≤ 0) && return 0.0
    s = 1 / (η + 1)
    turnover || return e_max * u^s
    p = 8 / (9 + η)
    u < p && return 0.8 * e_max * (u / p)^s
    return e_max - 0.2 * e_max * sqrt(max(1 - (u - p) / (1 - p), 0.0))
end

# Cumulative distribution of the eccentricity law of `_moe_draw_e` at e
function _moe_e_cdf(e::Float64, η::Float64, e_max::Float64, turnover::Bool)::Float64
    (η ≤ -1 || e_max ≤ 0) && return 1.0
    e ≤ 0 && return 0.0
    e ≥ e_max && return 1.0
    turnover || return (e / e_max)^(η + 1)
    p = 8 / (9 + η)
    e < 0.8 * e_max && return p * (e / (0.8 * e_max))^(η + 1)
    return p + (1 - p) * (1 - ((e_max - e) / (0.2 * e_max))^2)
end

# The eccentricity law of `_moe_draw_e` truncated to e ≤ e_hi, at the quantile u
function _moe_truncated_e(
    u::Float64,
    η::Float64,
    e_max::Float64,
    e_hi::Float64,
    turnover::Bool,
)::Float64
    e = _moe_draw_e(u * _moe_e_cdf(e_hi, η, e_max, turnover), η, e_max, turnover)
    return min(e, max(e_hi, 0.0))
end

"""
    _moe_draw(spec::MassiveBinarySpec, m₁, rng; table) -> (is_binary, lP, q, e)

One primary of mass `m₁` [M☉] under Moe & Di Stefano (2017, DOI
10.3847/1538-4365/aa6fb6) with the inner-binary sampler of Moe's code. With
`m₁c = _moe_clamp_m1(m₁)` and `table = _moe_period_table(m₁c, log_period_min,
log_period_max)` (passed in by a caller that caches it), one uniform `u` decides
both: `u ≥ f_bin` leaves the star single (`(false, NaN, NaN, NaN)`), otherwise
`lP = CumP⁻¹(u)`. Then `q` follows from `_moe_draw_q` at `(lP, m₁c)`
with `q_min`, and `e` from `_moe_draw_e` with `ecc_turnover`. Draws
one uniform for a single and four for a binary.
"""
function _moe_draw(
    spec::MassiveBinarySpec,
    m₁::Float64,
    rng::AbstractRNG;
    table::Tuple{Vector{Float64},Vector{Float64}} = _moe_period_table(
        _moe_clamp_m1(m₁),
        spec.log_period_min,
        spec.log_period_max,
    ),
)::Tuple{Bool,Float64,Float64,Float64}
    m₁c = _moe_clamp_m1(m₁)
    lPv, cum = table
    u = rand(rng)
    u < cum[end] || return (false, NaN, NaN, NaN)
    lP = _moe_invert_cum(lPv, cum, u)
    u₁ = rand(rng)
    u₂ = rand(rng)
    q = _moe_draw_q(u₁, u₂, _moe_q_parameters(lP, m₁c), spec.q_min)
    e = _moe_draw_e(rand(rng), _moe_eta(lP, m₁c), _moe_e_max(lP), spec.ecc_turnover)
    return (true, lP, q, e)
end

"""
    _moe_orbit(spec::MassiveBinarySpec, table, u_P, u_e, m₁c, m₁, m₂, a_rl, φ) -> (a_au, e, truncated)

Semi-major axis [AU] and eccentricity of a `moe2017` pair from the uniforms
`u_P` (period, over the inner-binary part of `table`) and `u_e` (eccentricity);
`m₁c` is the clamped primary mass of the statistics, `φ` the fill limit that set
`a_rl` (for the error message). For `a_rl > 0` the tabulated `CumP` is restricted
to `[lP(a_rl), log_period_max]` and renormalised, and the eccentricity law of
Moe & Di Stefano (2017, DOI 10.3847/1538-4365/aa6fb6) at the drawn `lP` is
truncated to `e ≤ min(e_max, 1 − a_rl/a)`; `truncated` flags a pair whose period
range the guard raised.
"""
function _moe_orbit(
    spec::MassiveBinarySpec,
    table::Tuple{Vector{Float64},Vector{Float64}},
    u_P::Float64,
    u_e::Float64,
    m₁c::Float64,
    m₁::Float64,
    m₂::Float64,
    a_rl::Float64,
    φ::Float64,
)::Tuple{Float64,Float64,Bool}
    lPv, cum = table
    M = m₁ + m₂
    truncated = false
    c_lo = 0.0
    lP_rl = -Inf
    if a_rl > 0
        lP_rl = _log_period_days(a_rl, M)
        lP_rl > spec.log_period_max &&
            _guard_range_error(m₁, m₂, a_rl, φ, "log(P/d) ≤ $(spec.log_period_max)")
        if lP_rl > lPv[1]
            c_lo = _moe_cum_at(lPv, cum, lP_rl)
            truncated = true
        end
    end
    c_lo < cum[end] || _guard_range_error(
        m₁,
        m₂,
        a_rl,
        φ,
        "log(P/d) ≤ $(spec.log_period_max) with a non-zero inner-binary probability",
    )
    lP = max(_moe_invert_cum(lPv, cum, c_lo + u_P * (cum[end] - c_lo)), lP_rl)
    a_au = _semi_major_axis_au(lP, M)
    e_max = _moe_e_max(lP)
    e_hi = e_max
    if a_rl > 0
        a_au = max(a_au, a_rl)
        e_hi = min(e_max, max(1 - a_rl / a_au, 0.0))
    end
    return a_au, _moe_truncated_e(u_e, _moe_eta(lP, m₁c), e_max, e_hi, spec.ecc_turnover), truncated
end

"""
    _moe_guarded_orbit(spec, table, lP, e, m₁c, m₁, m₂, a_rl, φ, rng) -> (a_au, e, truncated)

Orbit of a `moe2017` pair drawn at `(lP, e)` by `_moe_draw`, at its final
masses. Without the guard (`a_rl = 0`) it is `a(lP)` and `e`. With it, a drawn
`a < a_rl` is redrawn from the truncated laws by `_moe_orbit` (fresh
uniforms), and otherwise an `e > 1 − a_rl/a` is redrawn from the eccentricity law
truncated there (a fresh uniform); both leave the draw distributed as the
truncated law. `truncated` flags `lP(a_rl) > log_period_min`, as for Sana.
"""
function _moe_guarded_orbit(
    spec::MassiveBinarySpec,
    table::Tuple{Vector{Float64},Vector{Float64}},
    lP::Float64,
    e::Float64,
    m₁c::Float64,
    m₁::Float64,
    m₂::Float64,
    a_rl::Float64,
    φ::Float64,
    rng::AbstractRNG,
)::Tuple{Float64,Float64,Bool}
    M = m₁ + m₂
    a_au = _semi_major_axis_au(lP, M)
    a_rl > 0 || return (a_au, e, false)
    a_au < a_rl && return _moe_orbit(spec, table, rand(rng), rand(rng), m₁c, m₁, m₂, a_rl, φ)
    truncated = _log_period_days(a_rl, M) > spec.log_period_min
    e_cap = 1 - a_rl / a_au
    if e > e_cap
        e = _moe_truncated_e(rand(rng), _moe_eta(lP, m₁c), _moe_e_max(lP), e_cap, spec.ecc_turnover)
    end
    return a_au, e, truncated
end

"""
    moe2017_integrals(m₁; lP_lo = 0.2, lP_hi = 8.0, n = 20_001) -> (; f_mult_q03, f_mult_q01, f_bin)

Integrals of the Moe & Di Stefano (2017, DOI 10.3847/1538-4365/aa6fb6) model at
the primary mass `m₁` [M☉], for checking the implementation against the
published values.

# Arguments
- `m₁`: primary mass [M☉]; the fits hold for `0.8 < m₁ < 40` (no clamp applied)
- `lP_lo`, `lP_hi`: period range [log10 d] of the companion frequencies
- `n`: number of trapezoid nodes, `≥ 2`

# Returns
- `f_mult_q03`: `∫ f_logP;q>0.3 dlP` over `[lP_lo, lP_hi]` (their Eq. 1)
- `f_mult_q01`: `∫ f_logP;q>0.3 (1 + W_s) dlP`, the companions with `q > 0.1`
- `f_bin`: the inner-binary fraction `CumP[end]` of Moe's code over `[0.15, 8.0]`
  (cumulative trapezoid on the 0.005-dex grid of the sampler)

# Example
```julia
moe2017_integrals(28.0)   # (f_mult_q03 ≈ 1.03, f_mult_q01 ≈ 2.11, f_bin ≈ 0.96)
```
"""
function moe2017_integrals(
    m₁::Float64;
    lP_lo::Real = 0.2,
    lP_hi::Real = 8.0,
    n::Integer = 20_001,
)::NamedTuple{(:f_mult_q03, :f_mult_q01, :f_bin),NTuple{3,Float64}}
    (isfinite(m₁) && m₁ > 0) ||
        throw(ArgumentError("moe2017_integrals: m₁ must be finite and > 0; got $m₁"))
    lo, hi = Float64(lP_lo), Float64(lP_hi)
    lo < hi || throw(ArgumentError("moe2017_integrals: needs lP_lo < lP_hi; got $lo, $hi"))
    n ≥ 2 || throw(ArgumentError("moe2017_integrals: needs n ≥ 2; got $n"))
    h = (hi - lo) / (n - 1)
    s03 = 0.0
    s01 = 0.0
    for k in 0:(n - 1)
        x = lo + k * h
        w = (k == 0 || k == n - 1) ? 0.5h : h
        f = _moe_companion_frequency(x, m₁)
        s03 += w * f
        s01 += w * f * _moe_q_factor(x, m₁)
    end
    f_bin = last(last(_moe_period_table(m₁, 0.15, 8.0)))
    return (f_mult_q03 = s03, f_mult_q01 = s01, f_bin = f_bin)
end

# Union–find successor: the smallest unremoved position ≥ p (removed p points to p + 1).
function _successor!(parent::Vector{Int}, p::Int)::Int
    root = p
    while parent[root] != root
        root = parent[root]
    end
    while parent[p] != root
        next = parent[p]
        parent[p] = root
        p = next
    end
    return root
end

"""
Unpaired stars ordered by mass (ties by index) with deletion; `up` finds the next
unremoved position upwards, `down` (reversed positions) the next one downwards.
"""
struct _MassPool
    sorted_mass::Vector{Float64}
    order::Vector{Int}
    position::Vector{Int}
    up::Vector{Int}
    down::Vector{Int}
end

function _MassPool(m::Vector{Float64})
    n = length(m)
    order = sortperm(m; alg = MergeSort)
    position = Vector{Int}(undef, n)
    position[order] = 1:n
    return _MassPool(m[order], order, position, collect(1:(n + 1)), collect(1:(n + 1)))
end

function _pool_remove!(pool::_MassPool, star::Int)
    n = length(pool.order)
    p = pool.position[star]
    pool.up[p] = p + 1
    pool.down[n + 1 - p] = n + 2 - p
    return nothing
end

_pool_up!(pool::_MassPool, p::Int)::Int = _successor!(pool.up, p)
function _pool_down!(pool::_MassPool, p::Int)::Int
    n = length(pool.order)
    return n + 1 - _successor!(pool.down, n + 1 - p)
end

"""
    _pool_nearest!(pool, target, m_max) -> Int

Star index of the unremoved star of mass `≤ m_max` closest to `target ≤ m_max`
(ties to the lower index), or 0 when there is none.
"""
function _pool_nearest!(pool::_MassPool, target::Float64, m_max::Float64)::Int
    n = length(pool.order)
    p = searchsortedlast(pool.sorted_mass, target)
    lo = _pool_down!(pool, p)
    if lo > 0
        # The lowest index among the unremoved stars of that mass
        lo = _pool_up!(pool, searchsortedfirst(pool.sorted_mass, pool.sorted_mass[lo]))
    end
    hi = _pool_up!(pool, p + 1)
    (hi ≤ n && pool.sorted_mass[hi] ≤ m_max) || (hi = 0)
    lo == 0 && hi == 0 && return 0
    lo == 0 && return pool.order[hi]
    hi == 0 && return pool.order[lo]
    d_lo = target - pool.sorted_mass[lo]
    d_hi = pool.sorted_mass[hi] - target
    d_lo < d_hi && return pool.order[lo]
    d_hi < d_lo && return pool.order[hi]
    return min(pool.order[lo], pool.order[hi])
end

"""
    _pair_massive!(m, used, s::MassiveBinarySpec, rng, zmet, φ) -> (primary, secondary, a_au, e, n_guarded)

The massive branch: the stars are walked in descending mass (stable); each
unpaired star of mass `≥ m_min` becomes a primary when its Bernoulli trial
`rand(rng) < fraction` succeeds; it then draws inline the period uniform
`u_P`, the mass ratio `q` (power law of exponent `q_exponent` on
`[q_min, 1]`) and the eccentricity uniform `u_e`, and takes a partner by
`pool` or `replace`; without a partner it stays single. `u_P` and `u_e` are
mapped through the period and eccentricity laws at the final masses by
`_sana_orbit`, truncated by the Roche guard when `φ > 0`. Updates `m` (for
`"replace"`) and `used` in place.

With `model = "moe2017"` the Bernoulli trial and the three draws are replaced by
`_moe_draw` (binary decision, `lP`, `q` and `e` at the primary mass clamped to
0.8–40 M☉; a primary that draws single stays unused); the period tables are
cached per clamped mass, and the orbit at the final masses follows
`_moe_guarded_orbit`. Primaries outside 0.8–40 M☉ are counted and reported in
one warning.
"""
function _pair_massive!(
    m::Vector{Float64},
    used::BitVector,
    s::MassiveBinarySpec,
    rng::AbstractRNG,
    zmet::Float64,
    φ::Float64,
)
    primary = Int[]
    secondary = Int[]
    a_au = Float64[]
    e = Float64[]
    n_guarded = 0
    pool_mode = s.pairing == "pool"
    pool = _MassPool(pool_mode ? m : Float64[])
    low = pool_mode ? Int[] : [i for i in eachindex(m) if m[i] < s.m_min]
    moe = s.model == "moe2017"
    tables = Dict{Float64,Tuple{Vector{Float64},Vector{Float64}}}()
    n_clamped = 0
    no_table = (Float64[], Float64[])
    for i in sortperm(m; rev = true, alg = MergeSort)
        m[i] < s.m_min && break
        used[i] && continue
        u_P, u_e, lP, e_draw, m₁c = NaN, NaN, NaN, NaN, NaN
        table = no_table
        if moe
            m₁c = round(_moe_clamp_m1(m[i]); digits = 6)
            n_clamped += !(0.8 ≤ m[i] ≤ 40.0)
            haskey(tables, m₁c) ||
                (tables[m₁c] = _moe_period_table(m₁c, s.log_period_min, s.log_period_max))
            table = tables[m₁c]
            is_binary, lP, q, e_draw = _moe_draw(s, m[i], rng; table = table)
            is_binary || continue
        else
            rand(rng) < s.fraction || continue
            u_P = rand(rng)
            q = _power_law_quantile(rand(rng), s.q_min, 1.0, s.q_exponent)
            u_e = rand(rng)
        end
        j = 0
        if pool_mode
            _pool_remove!(pool, i)
            j = _pool_nearest!(pool, q * m[i], m[i])
            j > 0 && _pool_remove!(pool, j)
        elseif !isempty(low)
            k = rand(rng, eachindex(low))
            j = low[k]
            low[k] = low[end]
            pop!(low)
            m[j] = q * m[i]
        end
        used[i] = true
        j == 0 && continue
        used[j] = true
        a_rl = φ > 0 ? _roche_separation_au(m[i], m[j], zmet, φ) : 0.0
        a_k, e_k, truncated =
            moe ? _moe_guarded_orbit(s, table, lP, e_draw, m₁c, m[i], m[j], a_rl, φ, rng) :
            _sana_orbit(s, u_P, u_e, m[i], m[j], a_rl, φ)
        push!(primary, i)
        push!(secondary, j)
        push!(a_au, a_k)
        push!(e, e_k)
        n_guarded += truncated
    end
    n_clamped > 0 &&
        @warn "moe2017: $n_clamped primaries outside 0.8–40 M☉ use the statistics of the nearest bound"
    return primary, secondary, a_au, e, n_guarded
end

"""
    sample_binaries(spec::BinarySpec, masses, rng; zmet = NaN)
        -> (; masses, primary, secondary, m1, m2, a_pc, e, n_guarded)

Pair the stars of `masses` [M☉]. Without a massive branch
(`spec.massive.model == "none"`), `N_b = round(fraction N / (1 + fraction))` of
the `N` stars are paired by `spec.pairing`, `spec.period` and
`spec.eccentricity`. With `model = "sana2012"` the massive primaries are paired
first (`MassiveBinarySpec`), and the same rule then applies to the
unpaired stars below `massive.m_min`, `N` being their number. With
`model = "moe2017"` each massive primary is a binary with the inner-binary
probability of Moe's code and takes its period, mass ratio and eccentricity from
the joint distribution of Moe & Di Stefano (2017, DOI 10.3847/1538-4365/aa6fb6)
at its mass clamped to 0.8–40 M☉.

With `spec.roche_fill_max = φ > 0` every pair, at its final masses, has its
period law truncated to `a ≥ a_RL = max_i R_i/(φ r_L(q_i))` and its
eccentricity law to `e ≤ 1 − a_RL/a`, so that neither ZAMS component fills more
than `φ` of its Roche lobe at periastron; the radii need the metallicity
`zmet` (finite, `1e-4 ≤ zmet ≤ 0.03`). A pair whose `a_RL` lies beyond the
period range raises an `ArgumentError`.

Returns the (possibly modified, for `"uniform_q"` and `"replace"`) mass vector,
the star indices of the primaries and secondaries (massive pairs first), the
component masses [M☉] with `m1 ≥ m2`, the semi-major axes [pc], the
eccentricities, and `n_guarded`, the number of pairs whose period (or axis)
range the guard raised (its eccentricity cap acts on every pair). Every star
belongs to at most one pair.
Without the massive branch and the guard the draws are those of earlier
versions.
"""
function sample_binaries(
    spec::BinarySpec,
    masses::AbstractVector{Float64},
    rng::AbstractRNG;
    zmet::Real = NaN,
)
    φ = spec.roche_fill_max
    (φ > 0 && !isfinite(zmet)) &&
        throw(ArgumentError("sample_binaries: roche_fill_max = $φ needs a finite metallicity zmet"))
    z = Float64(zmet)
    N = length(masses)
    m = collect(Float64, masses)
    used = falses(N)
    massive_model = spec.massive.model != "none"
    primary_h, secondary_h, a_au_h, e_h, n_guarded =
        massive_model ? _pair_massive!(m, used, spec.massive, rng, z, φ) :
        (Int[], Int[], Float64[], Float64[], 0)

    # Low-mass population: every star without a massive branch, else the
    # unpaired stars below m_min
    subset =
        massive_model ? [i for i in 1:N if !used[i] && m[i] < spec.massive.m_min] : collect(1:N)
    n_sub = length(subset)
    N_b = round(Int, spec.fraction * n_sub / (1 + spec.fraction))
    2N_b ≤ n_sub || (N_b = n_sub ÷ 2)
    perm = subset[randperm(rng, n_sub)]
    primary = perm[1:N_b]
    secondary = perm[(N_b + 1):(2N_b)]
    for k in 1:N_b
        if spec.pairing == "uniform_q"
            q = spec.q_min + (1 - spec.q_min) * rand(rng)
            m[secondary[k]] = q * m[primary[k]]
        elseif m[secondary[k]] > m[primary[k]]
            primary[k], secondary[k] = secondary[k], primary[k]
        end
    end
    m1 = m[primary]
    m2 = m[secondary]
    a_pc = Vector{Float64}(undef, N_b)
    e = Vector{Float64}(undef, N_b)
    for k in 1:N_b
        if φ > 0
            a_rl = _roche_separation_au(m1[k], m2[k], z, φ)
            a_au, e[k], truncated = _guarded_orbit(spec, m1[k], m2[k], a_rl, rng)
            a_pc[k] = a_au * _AU_PC
            n_guarded += truncated
            continue
        end
        a_au = if spec.period == "kroupa1995"
            p_yr = exp10(_sample_log_period_kroupa1995(rng)) * _DAY_YR
            cbrt((m1[k] + m2[k]) * p_yr^2)        # Kepler III in AU, yr, M☉
        else
            spec.a_min * (spec.a_max / spec.a_min)^rand(rng)
        end
        a_pc[k] = a_au * _AU_PC
        e[k] = spec.eccentricity == "thermal" ? sqrt(rand(rng)) : 0.0
    end
    if massive_model
        primary = vcat(primary_h, primary)
        secondary = vcat(secondary_h, secondary)
        m1 = vcat(m[primary_h], m1)
        m2 = vcat(m[secondary_h], m2)
        a_pc = vcat(a_au_h .* _AU_PC, a_pc)
        e = vcat(e_h, e)
    end
    return (
        masses = m,
        primary = primary,
        secondary = secondary,
        m1 = m1,
        m2 = m2,
        a_pc = a_pc,
        e = e,
        n_guarded = n_guarded,
    )
end

"""
    _random_rotation(rng) -> Matrix{Float64}

Uniformly distributed rotation matrix (from a random unit quaternion).
"""
function _random_rotation(rng::AbstractRNG)::Matrix{Float64}
    u1, u2, u3 = rand(rng), rand(rng), rand(rng)
    q0 = sqrt(1 - u1) * sin(2π * u2)
    q1 = sqrt(1 - u1) * cos(2π * u2)
    q2 = sqrt(u1) * sin(2π * u3)
    q3 = sqrt(u1) * cos(2π * u3)
    return [
        1-2(q2^2+q3^2) 2(q1*q2-q0*q3) 2(q1*q3+q0*q2)
        2(q1*q2+q0*q3) 1-2(q1^2+q3^2) 2(q2*q3-q0*q1)
        2(q1*q3-q0*q2) 2(q2*q3+q0*q1) 1-2(q1^2+q2^2)
    ]
end

"""
    _kepler_relative_orbit(M, a, e, rng) -> (r_vec, v_vec)

Relative position and velocity of a two-body orbit of total mass `M`,
semi-major axis `a`, and eccentricity `e` (`G = 1`) at a uniformly random
mean anomaly (Kepler's equation by Newton iteration), in a uniformly random
orientation.
"""
function _kepler_relative_orbit(M::Real, a::Real, e::Real, rng::AbstractRNG)
    mean_anomaly = 2π * rand(rng)
    E = mean_anomaly
    for _ in 1:50
        δ = (E - e * sin(E) - mean_anomaly) / (1 - e * cos(E))
        E -= δ
        abs(δ) < 1e-13 && break
    end
    r = a * (1 - e * cos(E))
    cosν = (cos(E) - e) / (1 - e * cos(E))
    sinν = sqrt(1 - e^2) * sin(E) / (1 - e * cos(E))
    vscale = sqrt(M / (a * (1 - e^2)))
    r_plane = [r * cosν, r * sinν, 0.0]
    v_plane = [-vscale * sinν, vscale * (e + cosν), 0.0]
    R = _random_rotation(rng)
    return R * r_plane, R * v_plane
end

"""
    expand_binaries(pos, vel, mass, system_ranges, cluster_binaries, kept; rng)
        -> (; pos, vel, mass, cluster_blocks, n_pairs, hard_fraction)

Turn the combined *system* arrays (code units, `G = 1`) into *body* arrays:
every pair becomes two bodies on a Keplerian relative orbit about the
system's centre of mass, written first, cluster by cluster, followed by the
singles of every cluster (the engine's primordial-pair convention). For
cluster `i`, `system_ranges[i]` are its systems in the combined arrays,
`kept[i]` their indices in that cluster's sampled system set, and
`cluster_binaries[i]` the sampled pair table (`system_binary[j]` = pair
index of sampled system `j`, 0 for a single, plus `m1`, `m2`, `a_pc`, `e`).
`cluster_blocks[i]` are the contiguous body-index blocks of cluster `i`
(its pairs block and its singles block), `n_pairs[i]` its number of pairs,
and `hard_fraction[i]` the fraction of its pairs with binding energy above
`⟨m⟩ σ²` of the cluster's systems.

`rng` is one `AbstractRNG` for all clusters, drawn from cluster by cluster,
or a vector of one RNG per cluster, cluster `i` drawing from `rng[i]`.
"""
function expand_binaries(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    system_ranges::AbstractVector{<:AbstractVector{Int}},
    cluster_binaries::AbstractVector,
    kept::AbstractVector{<:AbstractVector{Int}};
    rng::Union{AbstractRNG,AbstractVector{<:AbstractRNG}} = Random.default_rng(),
)
    n_cl = length(system_ranges)
    (rng isa AbstractRNG || length(rng) == n_cl) ||
        throw(DimensionMismatch("expand_binaries: $(length(rng)) RNGs for $n_cl clusters"))
    pairs = [Int[] for _ in 1:n_cl]     # combined-array system index of each surviving pair
    pair_ids = [Int[] for _ in 1:n_cl]  # its index in the cluster's pair table
    singles = [Int[] for _ in 1:n_cl]
    for i in 1:n_cl
        cb = cluster_binaries[i]
        for (k, sys) in enumerate(system_ranges[i])
            p = cb.system_binary[kept[i][k]]
            if p > 0
                push!(pairs[i], sys)
                push!(pair_ids[i], p)
            else
                push!(singles[i], sys)
            end
        end
    end
    n_pairs = length.(pairs)
    nbin0 = sum(n_pairs)
    N_bodies = 2nbin0 + sum(length, singles)
    pos_b = zeros(Float64, 3, N_bodies)
    vel_b = zeros(Float64, 3, N_bodies)
    mass_b = zeros(Float64, N_bodies)
    cluster_blocks = Vector{Vector{UnitRange{Int}}}(undef, n_cl)
    hard_fraction = fill(NaN, n_cl)
    b = 0
    for i in 1:n_cl
        cb = cluster_binaries[i]
        rng_i = rng isa AbstractRNG ? rng : rng[i]
        # 1-D velocity dispersion of the cluster's systems for the hard/soft boundary
        idx_sys = system_ranges[i]
        M_sys = sum(mass[idx_sys])
        vc = vec(sum(vel[:, idx_sys] .* mass[idx_sys]'; dims = 2)) ./ M_sys
        σ2 = sum(mass[s] * sum((vel[:, s] .- vc) .^ 2) for s in idx_sys) / (3 * M_sys)
        m_mean = M_sys / length(idx_sys)
        n_hard = 0
        first_pair_body = b + 1
        for (sys, p) in zip(pairs[i], pair_ids[i])
            m1, m2, a, e = cb.m1[p], cb.m2[p], cb.a_pc[p], cb.e[p]
            M = m1 + m2
            r_rel, v_rel = _kepler_relative_orbit(M, a, e, rng_i)
            n_hard += (m1 * m2 / (2a) > m_mean * σ2) ? 1 : 0
            for (body, m, f) in ((b + 1, m1, m2 / M), (b + 2, m2, -m1 / M))
                mass_b[body] = m
                for k in 1:3
                    pos_b[k, body] = pos[k, sys] + f * r_rel[k]
                    vel_b[k, body] = vel[k, sys] + f * v_rel[k]
                end
            end
            b += 2
        end
        hard_fraction[i] = isempty(pairs[i]) ? NaN : n_hard / length(pairs[i])
        cluster_blocks[i] = [first_pair_body:b]
    end
    s = 2nbin0
    for i in 1:n_cl
        first_single = s + 1
        for sys in singles[i]
            s += 1
            mass_b[s] = mass[sys]
            pos_b[:, s] .= pos[:, sys]
            vel_b[:, s] .= vel[:, sys]
        end
        push!(cluster_blocks[i], first_single:s)
    end
    return (
        pos = pos_b,
        vel = vel_b,
        mass = mass_b,
        cluster_blocks = cluster_blocks,
        n_pairs = n_pairs,
        hard_fraction = hard_fraction,
    )
end

"""Body-index vector of a cluster from its contiguous blocks."""
_members_from_blocks(blocks::AbstractVector{<:UnitRange{Int}}) =
    reduce(vcat, (collect(r) for r in blocks); init = Int[])
