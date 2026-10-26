# =============================================================================
# Gravitational-wave catalogues: compact binaries, escaped binaries, mergers
# =============================================================================
#
# Rows of the compact-binary, escaped-binary and GW-merger catalogues, built
# from the event log of the engine, the binary spans, bev.82 and coal.24. Each
# recorded orbit is carried forward along the Peters (1964, DOI
# 10.1103/PhysRev.136.B1224) track of orbit-averaged quadrupole radiation; the
# band residence is the interval between the crossings of the two band edges,
# which does not depend on the output cadence of the run.

# au per pc and R☉ per au
const _GWCAT_AU_PER_PC = 206264.806
const _GWCAT_RSUN_PER_AU = 215.032
# Hubble time [Myr]
const _GWCAT_T_HUBBLE_MYR = 13_800.0
# Default GW band and reference frequencies of f₂₂ [Hz]
const _GWCAT_BAND_HZ = (1.0e-4, 1.0e-1)
const _GWCAT_FREF_HZ = (1.0e-2, 10.0)
# K* of white dwarfs, neutron stars and black holes
const _GWCAT_COMPACT = 10:14
# GWEvent kinds of an applied Peters step (GR KSTIDE, GR KSINT2, GR TIDES)
const _GWCAT_PETERS_KINDS = (:gw_decay_step, :gw_ksint, :gw_tides)
# yr per Myr, s per day
const _GWCAT_YR_PER_MYR = 1.0e6
const _GWCAT_S_PER_DAY = 86_400.0
# Relative perturbation from which the engine applies no GW decay to a perturbed KS pair
const _GWCAT_GAMMA_NO_BRAKE = 0.1
# Lower edge of x = c²p/(GM) in the calibration of Tucker & Will (2021)
const _GWCAT_X_PN = 100.0
# Engine merger trigger: pericentre within 25 Schwarzschild radii
const _GWCAT_TRIGGER_RS = 25.0
# Interval [NB] before an escape searched for a GW recoil or a natal kick
const _GWCAT_KICK_WINDOW_NB = 1.0
# Relative tolerance of the esc.11 time against the escape time
const _GWCAT_ESC_RTOL = 1.0e-6
# Pair classes by the ranks (1 star, 2 WD, 3 NS, 4 BH) of the two components, higher rank first
const _GWCAT_PAIR_CLASSES =
    ((:star_star,), (:wd_star, :wdwd), (:ns_star, :nswd, :nsns), (:bh_star, :bhwd, :bhns, :bhbh))

# ---------------------------------------------------------------------------
# Pair classes
# ---------------------------------------------------------------------------

"""Rank of a component type: 4 BH, 3 NS, 2 WD, 1 non-compact star, 0 absent or out of range."""
function _gwcat_kstar_rank(k::Integer)::Int
    (k == KSTAR_ABSENT || k < 0 || k ≥ 15) && return 0
    k == 14 && return 4
    k == 13 && return 3
    k ≥ 10 && return 2
    return 1
end

"""Both component types are white dwarfs, neutron stars or black holes."""
_gwcat_compact(k1::Integer, k2::Integer)::Bool = k1 in _GWCAT_COMPACT && k2 in _GWCAT_COMPACT

"""
    pair_class(kstar1::Integer, kstar2::Integer) -> Symbol

Class of a pair from the stellar types K* of its components, independent of
their order: `:bhbh`, `:bhns`, `:nsns`, `:bhwd`, `:nswd`, `:wdwd` for two
compact components (BH: K* = 14, NS: 13, WD: 10–12), `:bh_star`, `:ns_star`,
`:wd_star` for a compact component with a star of K* ≤ 9, `:star_star` for
two stars, and `:unknown` when either type is [`KSTAR_ABSENT`](@ref),
negative or ≥ 15.

# Example
```julia
pair_class(13, 14)   # :bhns
pair_class(1, 11)    # :wd_star
```
"""
function pair_class(kstar1::Integer, kstar2::Integer)::Symbol
    r1 = _gwcat_kstar_rank(kstar1)
    r2 = _gwcat_kstar_rank(kstar2)
    (r1 == 0 || r2 == 0) && return :unknown
    return _GWCAT_PAIR_CLASSES[max(r1, r2)][min(r1, r2)]
end

# ---------------------------------------------------------------------------
# Peters track, band residence, state after GW decay
# ---------------------------------------------------------------------------

"""g(e) = e^{12/19}(1 + 121e²/304)^{870/2299}/(1 − e²), the shape of the Peters track a(e) = c₀ g(e)."""
_gwcat_peters_g(e::Float64)::Float64 =
    e^(12 / 19) * (1 + _GW_K * e^2)^(870 / 2299) / ((1 - e) * (1 + e))

"""
    peters_semi_major_axis(a0::Real, e0::Real, e::Real) -> Float64

Semi-major axis [au] at eccentricity `e` on the Peters (1964,
DOI 10.1103/PhysRev.136.B1224) track through the state (`a0` [au], `e0`),

    a(e) = a₀ g(e)/g(e₀),   g(e) = e^{12/19}(1 + 121e²/304)^{870/2299}/(1 − e²).

GW decay lowers e, so the track is defined for 0 ≤ e ≤ e₀; a circular state
(`e0 = 0`) admits only `e = 0` and returns `a0`, and `e = 0` from an eccentric
state returns 0 (coalescence). Throws an `ArgumentError` for `a0 ≤ 0`, `e0`
outside [0, 1), `e < 0` or `e > e0`.

# Example
```julia
peters_semi_major_axis(1.0, 0.9, 0.5)   # ≈ 0.163 au
```
"""
function peters_semi_major_axis(a0::Real, e0::Real, e::Real)::Float64
    a0 > 0 || throw(ArgumentError("semi-major axis must be positive, got a0 = $a0"))
    0 ≤ e0 < 1 || throw(ArgumentError("eccentricity must lie in [0, 1), got e0 = $e0"))
    0 ≤ e ≤ e0 || throw(ArgumentError("eccentricity must lie in [0, e0] = [0, $e0], got e = $e"))
    ε₀ = Float64(e0)
    ε₀ == 0.0 && return Float64(a0)
    return Float64(a0) * (_gwcat_peters_g(Float64(e)) / _gwcat_peters_g(ε₀))
end

"""
    _gwcat_time_to_frequency(m1, m2, a, e, f, definition)
        -> NamedTuple{(:e, :t_yr, :above, :reached)}

Eccentricity and time [yr] at which the GW frequency of `definition` reaches
`f` along the Peters track of a valid state, from
[`eccentricity_at_frequency`](@ref). `above`: the state is already at or above
`f`; `reached = false` when the root lies below the lower limit of that search.
`e` and `t_yr` are `NaN` in both cases.
"""
function _gwcat_time_to_frequency(
    m1::Real,
    m2::Real,
    a::Real,
    e::Real,
    f::Real,
    definition::Symbol,
)::@NamedTuple{e::Float64, t_yr::Float64, above::Bool, reached::Bool}
    r = try
        eccentricity_at_frequency(m1, m2, a, e, f; definition = definition)
    catch err
        # The orbit, f and the definition are validated by the callers: the only
        # remaining ArgumentError is a frequency not reached within the bisection limit.
        err isa ArgumentError || rethrow()
        return (e = NaN, t_yr = NaN, above = false, reached = false)
    end
    r.above_at_record && return (e = NaN, t_yr = NaN, above = true, reached = true)
    return (e = r.e, t_yr = r.time_to_ref_yr, above = false, reached = true)
end

"""
    band_residence(m1, m2, a0, e0, f_lo, f_hi; definition = :f22)
        -> NamedTuple{(:t_in_yr, :t_out_yr, :above_lo, :above_hi)}

Times [yr] from the state (`a0` [au], `e0`) of a binary of masses `m1`, `m2`
[M☉] at which its GW frequency of the given `definition` (see
[`gw_frequency`](@ref)) enters the band [`f_lo`, `f_hi`] [Hz] and leaves it,
under the Peters (1964) evolution.

The frequency rises monotonically along the track. The crossing of a
frequency f is found with [`eccentricity_at_frequency`](@ref), which returns
the eccentricity e_f there; with a_f = a(e_f) from
[`peters_semi_major_axis`](@ref), the crossing time is
T(a₀, e₀) − T(a_f, e_f) of [`gw_inspiral_time`](@ref). A circular state stays
circular, all definitions reduce to f₂₂ and a_f = (GM/(πf)²)^{1/3}.

Returns
- `t_in_yr`: time to `f_lo`, 0 when the state is already at or above it;
- `t_out_yr`: time to `f_hi`, 0 when the state is already above it, the
  coalescence time T(a₀, e₀) when `f_hi` is not reached before coalescence;
- `above_lo`, `above_hi`: the state is already at or above `f_lo`, `f_hi`.

The residence time in band is `t_out_yr − t_in_yr`. Throws an
`ArgumentError` for an invalid orbit, an unknown definition, `f_lo ≤ 0`,
`f_hi ≤ 0` or `f_lo ≥ f_hi`.

# Example
```julia
b = band_residence(30.0, 30.0, 0.1, 0.0, 1e-4, 1e-1)
b.t_out_yr - b.t_in_yr   # ≈ 4.14×10⁶ yr
```
"""
function band_residence(
    m1::Real,
    m2::Real,
    a0::Real,
    e0::Real,
    f_lo::Real,
    f_hi::Real;
    definition::Symbol = :f22,
)::NamedTuple{(:t_in_yr, :t_out_yr, :above_lo, :above_hi),Tuple{Float64,Float64,Bool,Bool}}
    _gw_check_binary(m1, m2, a0, e0)
    (f_lo > 0 && f_hi > 0) ||
        throw(ArgumentError("band edges must be positive, got f_lo = $f_lo, f_hi = $f_hi"))
    f_lo < f_hi || throw(ArgumentError("band requires f_lo < f_hi, got f_lo = $f_lo, f_hi = $f_hi"))
    # Validates the definition
    gw_frequency(m1, m2, a0, e0; definition = definition)
    t_coal = gw_inspiral_time(m1, m2, a0, e0)
    lo = _gwcat_time_to_frequency(m1, m2, a0, e0, f_lo, definition)
    hi = _gwcat_time_to_frequency(m1, m2, a0, e0, f_hi, definition)
    t_in = lo.above ? 0.0 : lo.reached ? lo.t_yr : t_coal
    t_out = hi.above ? 0.0 : hi.reached ? hi.t_yr : t_coal
    return (t_in_yr = t_in, t_out_yr = t_out, above_lo = lo.above, above_hi = hi.above)
end

"""
    eccentricity_after(m1, m2, a0, e0, Δt_yr) -> Tuple{Float64,Float64}

State `(a, e)` (a in au) of a binary of masses `m1`, `m2` [M☉], recorded at
(`a0` [au], `e0`), after `Δt_yr` [yr] of pure GW decay along the Peters
(1964) track; `(NaN, NaN)` when `Δt_yr` reaches the coalescence time
T(a₀, e₀) of [`gw_inspiral_time`](@ref).

The remaining time τ = T(a₀, e₀) − Δt equals T(a(e), e) for one e on the track,
an increasing function of e: it is found by bisection on [0, e₀] to a
relative tolerance of 10⁻¹⁰ in τ (at most 200 iterations), and
a = [`peters_semi_major_axis`](@ref)(a₀, e₀, e). A circular state decays as
a = a₀(τ/T(a₀, 0))^{1/4} with e = 0. Throws an `ArgumentError` for an invalid
orbit or `Δt_yr < 0`.

# Example
```julia
T = gw_inspiral_time(30.0, 30.0, 0.1, 0.9)
a, e = eccentricity_after(30.0, 30.0, 0.1, 0.9, 0.3T)   # gw_inspiral_time(30.0, 30.0, a, e) ≈ 0.7T
```
"""
function eccentricity_after(
    m1::Real,
    m2::Real,
    a0::Real,
    e0::Real,
    Δt_yr::Real,
)::Tuple{Float64,Float64}
    _gw_check_binary(m1, m2, a0, e0)
    Δt_yr ≥ 0 || throw(ArgumentError("elapsed time must be non-negative, got Δt_yr = $Δt_yr"))
    t_coal = gw_inspiral_time(m1, m2, a0, e0)
    Δt_yr ≥ t_coal && return (NaN, NaN)
    α₀ = Float64(a0)
    ε₀ = Float64(e0)
    Δt_yr == 0 && return (α₀, ε₀)
    τ = t_coal - Float64(Δt_yr)
    ε₀ == 0.0 && return (α₀ * (τ / t_coal)^0.25, 0.0)

    remaining(ε) = ε == 0.0 ? 0.0 : gw_inspiral_time(m1, m2, peters_semi_major_axis(α₀, ε₀, ε), ε)
    lo = 0.0
    hi = ε₀
    ε = ε₀
    for _ in 1:200
        ε = (lo + hi) / 2
        t = remaining(ε)
        abs(t - τ) ≤ 1.0e-10 * τ && break
        if t < τ
            lo = ε
        else
            hi = ε
        end
    end
    return (peters_semi_major_axis(α₀, ε₀, ε), ε)
end

# ---------------------------------------------------------------------------
# Derived quantities of one record
# ---------------------------------------------------------------------------

"""
    _gwcat_fdot_circular(mchirp_msun, f_hz) -> Float64

Orbit-averaged GW frequency derivative ḟ = (96/5)π^{8/3}(G𝓜/c³)^{5/3}f^{11/3}
[Hz s⁻¹] of a circular binary of chirp mass `mchirp_msun` [M☉] at GW
frequency `f_hz` (Robson et al. 2019, DOI 10.1088/1361-6382/ab1101).
"""
function _gwcat_fdot_circular(mchirp_msun::Real, f_hz::Real)::Float64
    t_m = _GW_G * Float64(mchirp_msun) * _GW_MSUN / _GW_C^3
    return 96 / 5 * π^(8 / 3) * t_m^(5 / 3) * Float64(f_hz)^(11 / 3)
end

"""Schwarzschild radius 2GM/c² [au] of a total mass `mtot_msun` [M☉]."""
_gwcat_schwarzschild_au(mtot_msun::Float64)::Float64 =
    2 * _GW_G * mtot_msun * _GW_MSUN / _GW_C^2 / _GW_AU

"""Positive finite masses, a positive finite semi-major axis and 0 ≤ e < 1."""
_gwcat_valid_orbit(m1::Float64, m2::Float64, a::Float64, e::Float64)::Bool =
    isfinite(m1) && isfinite(m2) && m1 > 0 && m2 > 0 && isfinite(a) && a > 0 && 0 ≤ e < 1

"""Mass ratio m₂/m₁, symmetric mass ratio, total and chirp mass; `NaN` for invalid masses."""
function _gwcat_mass_state(
    m1::Float64,
    m2::Float64,
)::@NamedTuple{q::Float64, eta::Float64, mtot_msun::Float64, mchirp_msun::Float64}
    (isfinite(m1) && isfinite(m2) && m1 > 0 && m2 > 0) ||
        return (q = NaN, eta = NaN, mtot_msun = NaN, mchirp_msun = NaN)
    mtot = m1 + m2
    return (q = m2 / m1, eta = m1 * m2 / mtot^2, mtot_msun = mtot, mchirp_msun = chirp_mass(m1, m2))
end

const _GWCAT_NO_ORBIT = (
    q_peri_au = NaN,
    period_d = NaN,
    f_orb_hz = NaN,
    f_gw22_hz = NaN,
    f_gw_peak_wen_hz = NaN,
    f_gw_peri_hz = NaN,
    fdot_gw22_hz_s = NaN,
    t_gw_circ_myr = NaN,
    t_gw_peters_myr = NaN,
    t_gw_mandel_myr = NaN,
    q_f_zwick = NaN,
    flag_pn_regime = false,
)

"""Orbit, frequencies and coalescence times of a record (masses M☉, a in au); `NaN` for an invalid orbit."""
function _gwcat_orbit_state(m1::Float64, m2::Float64, a::Float64, e::Float64)
    _gwcat_valid_orbit(m1, m2, a, e) || return _GWCAT_NO_ORBIT
    f_orb = orbital_frequency(m1, m2, a)
    f22 = gw_frequency(m1, m2, a, e; definition = :f22)
    gm = _GW_G * (m1 + m2) * _GW_MSUN
    p_m = a * _GW_AU * (1 - e) * (1 + e)
    return (
        q_peri_au = a * (1 - e),
        period_d = 1 / (f_orb * _GWCAT_S_PER_DAY),
        f_orb_hz = f_orb,
        f_gw22_hz = f22,
        f_gw_peak_wen_hz = gw_frequency(m1, m2, a, e; definition = :wen),
        f_gw_peri_hz = gw_frequency(m1, m2, a, e; definition = :peri),
        fdot_gw22_hz_s = _gwcat_fdot_circular(chirp_mass(m1, m2), f22) * gw_enhancement(e),
        t_gw_circ_myr = _gw_circular_time(m1, m2, a) / _GW_YR / _GWCAT_YR_PER_MYR,
        t_gw_peters_myr = gw_inspiral_time(m1, m2, a, e) / _GWCAT_YR_PER_MYR,
        t_gw_mandel_myr = gw_inspiral_time_mandel(m1, m2, a, e) / _GWCAT_YR_PER_MYR,
        q_f_zwick = gw_zwick_factor(m1, m2, a, e),
        flag_pn_regime = _GW_C^2 * p_m / gm < _GWCAT_X_PN,
    )
end

"""Eccentricities and times [Myr] at which f₂₂ reaches the two reference frequencies; `NaN` when unavailable."""
function _gwcat_reference(
    m1::Float64,
    m2::Float64,
    a::Float64,
    e::Float64,
    f_ref::NTuple{2,Float64},
)::@NamedTuple{e_lo::Float64, e_hi::Float64, t_lo_myr::Float64, t_hi_myr::Float64}
    _gwcat_valid_orbit(m1, m2, a, e) ||
        return (e_lo = NaN, e_hi = NaN, t_lo_myr = NaN, t_hi_myr = NaN)
    lo = _gwcat_time_to_frequency(m1, m2, a, e, f_ref[1], :f22)
    hi = _gwcat_time_to_frequency(m1, m2, a, e, f_ref[2], :f22)
    return (
        e_lo = lo.e,
        e_hi = hi.e,
        t_lo_myr = lo.t_yr / _GWCAT_YR_PER_MYR,
        t_hi_myr = hi.t_yr / _GWCAT_YR_PER_MYR,
    )
end

"""[`band_residence`](@ref) of f₂₂ in Myr; `NaN` and `false` for an invalid orbit."""
function _gwcat_band(
    m1::Float64,
    m2::Float64,
    a::Float64,
    e::Float64,
    band::NTuple{2,Float64},
)::@NamedTuple{t_in_myr::Float64, t_out_myr::Float64, above_lo::Bool, above_hi::Bool}
    _gwcat_valid_orbit(m1, m2, a, e) ||
        return (t_in_myr = NaN, t_out_myr = NaN, above_lo = false, above_hi = false)
    b = band_residence(m1, m2, a, e, band[1], band[2])
    return (
        t_in_myr = b.t_in_yr / _GWCAT_YR_PER_MYR,
        t_out_myr = b.t_out_yr / _GWCAT_YR_PER_MYR,
        above_lo = b.above_lo,
        above_hi = b.above_hi,
    )
end

# ---------------------------------------------------------------------------
# Lookups and argument checks
# ---------------------------------------------------------------------------

"""Unordered pair of names as `(min, max)`."""
_gwcat_key(n1::Integer, n2::Integer)::Tuple{Int,Int} =
    n1 < n2 ? (Int(n1), Int(n2)) : (Int(n2), Int(n1))

"""Spans by unordered pair of names."""
function _gwcat_span_index(
    spans::AbstractVector{BinarySpan},
)::Dict{Tuple{Int,Int},Vector{BinarySpan}}
    index = Dict{Tuple{Int,Int},Vector{BinarySpan}}()
    for s in spans
        push!(get!(() -> BinarySpan[], index, _gwcat_key(s.name1, s.name2)), s)
    end
    return index
end

"""First span of the pair with `t_form_nb ≤ t ≤ t_end_nb + tol`, `nothing` without one."""
function _gwcat_covering_span(
    index::Dict{Tuple{Int,Int},Vector{BinarySpan}},
    n1::Integer,
    n2::Integer,
    t::Float64,
    tol::Float64,
)::Union{Nothing,BinarySpan}
    key = _gwcat_key(n1, n2)
    haskey(index, key) || return nothing
    for s in index[key]
        s.t_form_nb ≤ t ≤ s.t_end_nb + tol && return s
    end
    return nothing
end

"""Times [NB] of the escape lines by unordered pair of names (pairs with a name 0 excluded)."""
function _gwcat_escape_times(log::BinaryEventLog)::Dict{Tuple{Int,Int},Vector{Float64}}
    times = Dict{Tuple{Int,Int},Vector{Float64}}()
    for ev in log.escape
        (ev.name1 != 0 && ev.name2 != 0) || continue
        push!(get!(() -> Float64[], times, _gwcat_key(ev.name1, ev.name2)), ev.time_nb)
    end
    return times
end

"""GW events of the given kinds by unordered pair of names (pairs with a name 0 excluded)."""
function _gwcat_gw_index(log::BinaryEventLog, kinds::Tuple)::Dict{Tuple{Int,Int},Vector{GWEvent}}
    index = Dict{Tuple{Int,Int},Vector{GWEvent}}()
    for ev in log.gw
        (ev.kind in kinds && ev.name1 != 0 && ev.name2 != 0) || continue
        push!(get!(() -> GWEvent[], index, _gwcat_key(ev.name1, ev.name2)), ev)
    end
    return index
end

function _gwcat_check_positive(x::Real, name::String, caller::String)
    (isfinite(x) && x > 0) || throw(ArgumentError("$caller: $name must be finite and > 0; got $x"))
    return nothing
end

function _gwcat_check_band(band, name::String, caller::String)::NTuple{2,Float64}
    length(band) == 2 || throw(ArgumentError("$caller: $name must hold two frequencies; got $band"))
    lo, hi = band
    (isfinite(lo) && isfinite(hi) && 0 < lo < hi) ||
        throw(ArgumentError("$caller: $name must satisfy 0 < lo < hi [Hz]; got $band"))
    return (Float64(lo), Float64(hi))
end

# ---------------------------------------------------------------------------
# Compact-binary catalogue
# ---------------------------------------------------------------------------

"""
    CompactBinaryRow

One record of a pair with two compact components (K* 10–14) in the
compact-binary catalogue, from [`compact_binary_catalogue`](@ref). Times in
Myr unless the suffix says otherwise; `NaN`, `0` and `:none` mark unavailable
values. The GW quantities require a bound orbit with positive masses; all
times to frequency crossings are counted from the record.

| field | meaning | unit |
|---|---|---|
| `record_source` | `:gw_decay` (GR KSTIDE, GR KSINT2, GR TIDES lines), `:ns_bh_binary` (NS/BH BINARY line), `:bev82` (bev.82 record) | — |
| `time_nb`, `time_myr` | time of the record | NB, Myr |
| `name1`, `name2` | component names, `name1` the heavier | — |
| `kstar1`, `kstar2` | component types | — |
| `pair_class` | [`pair_class`](@ref) of the two types | — |
| `m1_msun`, `m2_msun` | component masses, m₁ ≥ m₂ | M☉ |
| `q`, `eta` | m₂/m₁ ≤ 1, m₁m₂/M² | — |
| `mtot_msun`, `mchirp_msun` | M = m₁ + m₂, chirp mass 𝓜 | M☉ |
| `a_au`, `e` | semi-major axis and eccentricity of the record | au, — |
| `q_peri_au` | pericentre a(1 − e) | au |
| `period_d` | Keplerian period 1/f_orb | d |
| `f_orb_hz` | orbital frequency | Hz |
| `f_gw22_hz`, `f_gw_peak_wen_hz`, `f_gw_peri_hz` | GW frequency, definitions `:f22`, `:wen`, `:peri` of [`gw_frequency`](@ref) | Hz |
| `fdot_gw22_hz_s` | (96/5)π^{8/3}(G𝓜/c³)^{5/3}f₂₂^{11/3} times [`gw_enhancement`](@ref)(e) | Hz s⁻¹ |
| `t_gw_circ_myr` | circular coalescence time at the same a | Myr |
| `t_gw_peters_myr`, `t_gw_mandel_myr` | coalescence time: Peters (1964) quadrature, Mandel (2021) fit | Myr |
| `q_f_zwick` | 1PN factor of [`gw_zwick_factor`](@ref) | — |
| `t_merge_abs_myr` | `time_myr + t_gw_peters_myr` | Myr |
| `merges_before_t_hubble` | `t_merge_abs_myr ≤ t_hubble_myr` | — |
| `tgr_engine_nb` | engine TGR = a/abs(ȧ) of a GW decay record | NB |
| `gamma_pert` | relative perturbation γ of the record | — |
| `e_at_fref_lo`, `e_at_fref_hi` | e when f₂₂ reaches the reference frequencies; `NaN` when already above or not reached | — |
| `t_to_fref_lo_myr`, `t_to_fref_hi_myr` | times to these crossings | Myr |
| `t_band_in_myr`, `t_band_out_myr` | times to the entry of f₂₂ into the band and to its exit ([`band_residence`](@ref)) | Myr |
| `above_band_lo`, `above_band_hi` | f₂₂ of the record already above the lower, the upper band edge | — |
| `channel` | `:ejected` (an escape line of the pair precedes the record), else `:primordial` or `:dynamical` by the covering span, `:none` without one | — |
| `span_origin` | `origin` of the covering [`BinarySpan`](@ref), `:none` without one | — |
| `flag_pn_regime` | x = c²a(1 − e²)/(GM) < 100, below the calibration of Tucker & Will (2021) | — |
| `flag_perturbed_no_brake` | γ ≥ 0.1: the engine applies no GW decay to the pair while perturbed | — |
| `record_kind` | the event kind of the record (`:gw_decay_step` for a GR KSTIDE line, `:gw_ksint`/`:gw_tides` for the decision lines, `:ks_new_compact`, `:bev82`) | — |

An application of the decay writes a decision line (GR KSINT2 or GR TIDES) and
the KSTIDE line at the same time with the same state: the decision line is
dropped when the KSTIDE line of the pair at that time is present, so that one
application gives one row (`decision_lines = true` keeps both).
"""
Base.@kwdef struct CompactBinaryRow
    record_source::Symbol
    record_kind::Symbol = :none
    time_nb::Float64
    time_myr::Float64
    name1::Int
    name2::Int
    kstar1::Int
    kstar2::Int
    pair_class::Symbol
    m1_msun::Float64
    m2_msun::Float64
    q::Float64
    eta::Float64
    mtot_msun::Float64
    mchirp_msun::Float64
    a_au::Float64
    e::Float64
    q_peri_au::Float64
    period_d::Float64
    f_orb_hz::Float64
    f_gw22_hz::Float64
    f_gw_peak_wen_hz::Float64
    f_gw_peri_hz::Float64
    fdot_gw22_hz_s::Float64
    t_gw_circ_myr::Float64
    t_gw_peters_myr::Float64
    t_gw_mandel_myr::Float64
    q_f_zwick::Float64
    t_merge_abs_myr::Float64
    merges_before_t_hubble::Bool
    tgr_engine_nb::Float64
    gamma_pert::Float64
    e_at_fref_lo::Float64
    e_at_fref_hi::Float64
    t_to_fref_lo_myr::Float64
    t_to_fref_hi_myr::Float64
    t_band_in_myr::Float64
    t_band_out_myr::Float64
    above_band_lo::Bool
    above_band_hi::Bool
    channel::Symbol
    span_origin::Symbol
    flag_pn_regime::Bool
    flag_perturbed_no_brake::Bool
end

"""Count a record for source `i`: `false` (and counted) for an absent type, `false` for a non-compact pair."""
function _gwcat_admit!(
    kept::Vector{Int},
    absent::Vector{Int},
    i::Int,
    k1::Integer,
    k2::Integer,
)::Bool
    if k1 == KSTAR_ABSENT || k2 == KSTAR_ABSENT
        absent[i] += 1
        return false
    end
    _gwcat_compact(k1, k2) || return false
    kept[i] += 1
    return true
end

"""One [`CompactBinaryRow`](@ref) from the fields of a record, masses in M☉, `a` in au."""
function _gwcat_compact_row(
    ctx,
    source::Symbol,
    kind::Symbol,
    time_nb::Float64,
    n1::Int,
    n2::Int,
    k1::Int,
    k2::Int,
    m1::Float64,
    m2::Float64,
    a::Float64,
    e::Float64,
    gamma::Float64,
    tgr_nb::Float64,
)::CompactBinaryRow
    if m2 > m1
        n1, n2 = n2, n1
        k1, k2 = k2, k1
        m1, m2 = m2, m1
    end
    span = _gwcat_covering_span(ctx.spans, n1, n2, time_nb, 0.0)
    escapes = get(ctx.escapes, _gwcat_key(n1, n2), Float64[])
    ejected = any(t -> t < time_nb, escapes)
    channel =
        ejected ? :ejected : span === nothing ? :none : span.primordial ? :primordial : :dynamical
    ms = _gwcat_mass_state(m1, m2)
    orb = _gwcat_orbit_state(m1, m2, a, e)
    ref = _gwcat_reference(m1, m2, a, e, ctx.f_ref)
    band = _gwcat_band(m1, m2, a, e, ctx.band)
    time_myr = time_nb * ctx.tstar
    t_merge = time_myr + orb.t_gw_peters_myr
    return CompactBinaryRow(;
        record_source = source,
        record_kind = kind,
        time_nb = time_nb,
        time_myr = time_myr,
        name1 = n1,
        name2 = n2,
        kstar1 = k1,
        kstar2 = k2,
        pair_class = pair_class(k1, k2),
        m1_msun = m1,
        m2_msun = m2,
        q = ms.q,
        eta = ms.eta,
        mtot_msun = ms.mtot_msun,
        mchirp_msun = ms.mchirp_msun,
        a_au = a,
        e = e,
        q_peri_au = orb.q_peri_au,
        period_d = orb.period_d,
        f_orb_hz = orb.f_orb_hz,
        f_gw22_hz = orb.f_gw22_hz,
        f_gw_peak_wen_hz = orb.f_gw_peak_wen_hz,
        f_gw_peri_hz = orb.f_gw_peri_hz,
        fdot_gw22_hz_s = orb.fdot_gw22_hz_s,
        t_gw_circ_myr = orb.t_gw_circ_myr,
        t_gw_peters_myr = orb.t_gw_peters_myr,
        t_gw_mandel_myr = orb.t_gw_mandel_myr,
        q_f_zwick = orb.q_f_zwick,
        t_merge_abs_myr = t_merge,
        merges_before_t_hubble = t_merge ≤ ctx.t_hubble,
        tgr_engine_nb = tgr_nb,
        gamma_pert = gamma,
        e_at_fref_lo = ref.e_lo,
        e_at_fref_hi = ref.e_hi,
        t_to_fref_lo_myr = ref.t_lo_myr,
        t_to_fref_hi_myr = ref.t_hi_myr,
        t_band_in_myr = band.t_in_myr,
        t_band_out_myr = band.t_out_myr,
        above_band_lo = band.above_lo,
        above_band_hi = band.above_hi,
        channel = channel,
        span_origin = span === nothing ? :none : span.origin,
        flag_pn_regime = orb.flag_pn_regime,
        flag_perturbed_no_brake = gamma ≥ _GWCAT_GAMMA_NO_BRAKE,
    )
end

"""
    compact_binary_catalogue(log::BinaryEventLog, spans::AbstractVector{BinarySpan};
                             tstar_myr, rbar_pc, binary_evo = BinaryEvolutionSnapshot[],
                             band_hz = (1e-4, 1e-1), f_ref_hz = (1e-2, 10.0),
                             t_hubble_myr = 13_800.0, decision_lines = false) -> Vector{CompactBinaryRow}

Compact-binary catalogue of a run: one [`CompactBinaryRow`](@ref) per record
of a pair whose two components are white dwarfs, neutron stars or black
holes (K* 10–14).

# Arguments
- `log`: event records of the engine log, from [`read_binary_events`](@ref)
- `spans`: binary spans of the run, from [`binary_spans`](@ref)
- `tstar_myr`, `rbar_pc`: time unit T* [Myr] and length unit R̄ [pc]
- `binary_evo`: bev.82 snapshots ([`BinaryEvolutionSnapshot`](@ref)), optional
- `band_hz`: GW band (lo, hi) of the residence times [Hz]
- `f_ref_hz`: the two reference frequencies of f₂₂ [Hz]
- `t_hubble_myr`: Hubble time [Myr] of `merges_before_t_hubble`

# Sources
- every [`GWEvent`](@ref) of an applied Peters step (`:gw_decay_step`,
  `:gw_ksint`, `:gw_tides`): the orbit before the step, its γ and TGR;
- every [`KSEvent`](@ref) of kind `:ks_new_compact` (NS/BH BINARY);
- every bev.82 record, with a = 10^{log₁₀(a/R☉)}/215.032 au.

Lengths convert as a [au] = a [NB]·R̄·206264.806 and times as t [Myr] =
t [NB]·T*. Records with an absent type are dropped and counted, non-compact
pairs skipped; one `@info` reports rows and dropped records per source. The
covering span is the first [`BinarySpan`](@ref) of the same unordered pair
with `t_form_nb ≤ time_nb ≤ t_end_nb`. Rows are sorted by
`(time_nb, name1, name2)`. Throws an `ArgumentError` for a non-positive
`tstar_myr` or `rbar_pc`, or a band or reference pair without 0 < lo < hi.

# Example
```julia
log = read_binary_events(joinpath(dir, "out1000"))
spans = binary_spans(log; nzero = 20_000, nbin0 = 2_000, tstar_myr = 0.87)
rows = compact_binary_catalogue(log, spans; tstar_myr = 0.87, rbar_pc = 2.5,
                                binary_evo = read_all_binary_evolution(dir))
write_compact_binary_catalogue(joinpath(dir, "compact_binary_catalogue.csv"), rows)
```
"""
function compact_binary_catalogue(
    log::BinaryEventLog,
    spans::AbstractVector{BinarySpan};
    tstar_myr::Real,
    rbar_pc::Real,
    binary_evo::AbstractVector = BinaryEvolutionSnapshot[],
    band_hz = _GWCAT_BAND_HZ,
    f_ref_hz = _GWCAT_FREF_HZ,
    t_hubble_myr::Real = _GWCAT_T_HUBBLE_MYR,
    decision_lines::Bool = false,
)::Vector{CompactBinaryRow}
    caller = "compact_binary_catalogue"
    _gwcat_check_positive(tstar_myr, "tstar_myr", caller)
    _gwcat_check_positive(rbar_pc, "rbar_pc", caller)
    band = _gwcat_check_band(band_hz, "band_hz", caller)
    f_ref = _gwcat_check_band(f_ref_hz, "f_ref_hz", caller)
    ctx = (
        tstar = Float64(tstar_myr),
        band = band,
        f_ref = f_ref,
        t_hubble = Float64(t_hubble_myr),
        spans = _gwcat_span_index(spans),
        escapes = _gwcat_escape_times(log),
    )
    au_per_nb = Float64(rbar_pc) * _GWCAT_AU_PER_PC
    kept = zeros(Int, 3)
    absent = zeros(Int, 3)
    rows = CompactBinaryRow[]

    for ev in log.gw
        ev.kind in _GWCAT_PETERS_KINDS || continue
        _gwcat_admit!(kept, absent, 1, ev.kstar1, ev.kstar2) || continue
        push!(
            rows,
            _gwcat_compact_row(
                ctx,
                :gw_decay,
                ev.kind,
                ev.time_nb,
                ev.name1,
                ev.name2,
                ev.kstar1,
                ev.kstar2,
                ev.m1_msun,
                ev.m2_msun,
                ev.semi_nb * au_per_nb,
                ev.ecc,
                ev.gamma,
                ev.t_gr_nb,
            ),
        )
    end
    for ev in log.ks
        ev.kind === :ks_new_compact || continue
        _gwcat_admit!(kept, absent, 2, ev.kstar1, ev.kstar2) || continue
        push!(
            rows,
            _gwcat_compact_row(
                ctx,
                :ns_bh_binary,
                :ks_new_compact,
                ev.time_nb,
                ev.name1,
                ev.name2,
                ev.kstar1,
                ev.kstar2,
                ev.m1_msun,
                ev.m2_msun,
                ev.semi_nb * au_per_nb,
                ev.ecc,
                ev.gamma,
                NaN,
            ),
        )
    end
    for snap in binary_evo, r in snap.records
        _gwcat_admit!(kept, absent, 3, r.stellar_type1, r.stellar_type2) || continue
        push!(
            rows,
            _gwcat_compact_row(
                ctx,
                :bev82,
                :bev82,
                r.time_nb,
                Int(r.name1),
                Int(r.name2),
                Int(r.stellar_type1),
                Int(r.stellar_type2),
                r.mass1,
                r.mass2,
                exp10(r.log_semi_major_axis_rsun) / _GWCAT_RSUN_PER_AU,
                r.eccentricity,
                NaN,
                NaN,
            ),
        )
    end
    if !decision_lines
        applied =
            Set((r.name1, r.name2, r.time_nb) for r in rows if r.record_kind === :gw_decay_step)
        filter!(
            r -> !(
                r.record_kind in (:gw_ksint, :gw_tides) && (r.name1, r.name2, r.time_nb) in applied
            ),
            rows,
        )
    end
    sort!(rows; by = r -> (r.time_nb, r.name1, r.name2))
    @info(
        "Compact-binary catalogue: $(length(rows)) rows",
        gw_decay = (rows = kept[1], kstar_absent = absent[1]),
        ns_bh_binary = (rows = kept[2], kstar_absent = absent[2]),
        bev82 = (rows = kept[3], kstar_absent = absent[3]),
    )
    return rows
end

# ---------------------------------------------------------------------------
# Escaped-binary catalogue
# ---------------------------------------------------------------------------

"""
    EscapedBinaryRow

One escape of a pair with two compact components (K* 10–14) in the
escaped-binary catalogue, from [`escaped_binary_catalogue`](@ref). Times in
Myr unless the suffix says otherwise; `NaN`, `0` and `:none` mark unavailable
values. The orbit is that of the escape line; the GW quantities require a
bound orbit with positive masses, and all times to frequency crossings are
counted from the escape.

| field | meaning | unit |
|---|---|---|
| `t_esc_nb`, `t_esc_myr` | escape time (TTOT of the closing `ESCAPE T,NESC` line) | NB, Myr |
| `name1`, `name2`, `name_cm` | component names (`name1` the heavier), c.m. name | — |
| `kstar1`, `kstar2` | component types | — |
| `pair_class` | [`pair_class`](@ref) of the two types | — |
| `m1_msun`, `m2_msun` | component masses, m₁ ≥ m₂ | M☉ |
| `q`, `mchirp_msun` | m₂/m₁ ≤ 1, chirp mass | —, M☉ |
| `a_au`, `e` | semi-major axis and eccentricity at escape | au, — |
| `period_d` | period as printed | d |
| `eb_nb` | binding energy as printed | NB |
| `r_esc_pc` | c.m. distance from the density centre, r·R̄ | pc |
| `v_esc_kms` | c.m. speed, v·V* | km s⁻¹ |
| `escape_mode` | `:gw_recoil` (a GW recoil of the c.m. or of the same pair within 1 NB before the escape), `:sn_kick` (a BINARY KICK line naming a component in that interval), `:dynamical` | — |
| `t_gw_peters_myr`, `t_gw_mandel_myr` | coalescence time at escape: Peters (1964), Mandel (2021) | Myr |
| `q_f_zwick` | 1PN factor of [`gw_zwick_factor`](@ref) | — |
| `t_merge_abs_myr` | `t_esc_myr + t_gw_peters_myr` | Myr |
| `merges_before_t_hubble` | `t_merge_abs_myr ≤ t_hubble_myr` | — |
| `e_at_fref_lo`, `e_at_fref_hi` | e when f₂₂ reaches the reference frequencies; `NaN` when already above or not reached | — |
| `t_band_in_myr`, `t_band_out_myr` | times to the entry of f₂₂ into the band and to its exit ([`band_residence`](@ref)) | Myr |
| `f_gw22_now_hz`, `e_now`, `a_now_au` | state at `t_now_myr` after GW decay since the escape ([`eccentricity_after`](@ref)) | Hz, —, au |
| `t_to_merge_from_now_myr` | coalescence time of that state | Myr |
| `flag_kept_in_sim` | `false` when esc.11 lists a component within a relative 10⁻⁶ of `t_esc_myr` (the engine removed the pair), `true` otherwise | — |

The present-epoch columns are `NaN` when `t_now_myr` is `NaN` or earlier than
the escape, and when the pair merges before it.
"""
Base.@kwdef struct EscapedBinaryRow
    t_esc_nb::Float64
    t_esc_myr::Float64
    name1::Int
    name2::Int
    name_cm::Int
    kstar1::Int
    kstar2::Int
    pair_class::Symbol
    m1_msun::Float64
    m2_msun::Float64
    q::Float64
    mchirp_msun::Float64
    a_au::Float64
    e::Float64
    period_d::Float64
    eb_nb::Float64
    r_esc_pc::Float64
    v_esc_kms::Float64
    escape_mode::Symbol
    t_gw_peters_myr::Float64
    t_gw_mandel_myr::Float64
    q_f_zwick::Float64
    t_merge_abs_myr::Float64
    merges_before_t_hubble::Bool
    e_at_fref_lo::Float64
    e_at_fref_hi::Float64
    t_band_in_myr::Float64
    t_band_out_myr::Float64
    f_gw22_now_hz::Float64
    e_now::Float64
    a_now_au::Float64
    t_to_merge_from_now_myr::Float64
    flag_kept_in_sim::Bool
end

"""Escape mode of an escape line: a GW recoil, then a natal kick within the preceding interval, else dynamical."""
function _gwcat_escape_mode(
    rec::BinaryEscapeEvent,
    recoils::Vector{GWEvent},
    kicks::Vector{MassLossEvent},
    tstar::Float64,
)::Symbol
    t = rec.time_nb
    isfinite(t) || return :dynamical
    in_window(x) = t - _GWCAT_KICK_WINDOW_NB ≤ x ≤ t
    named = rec.name1 != 0 && rec.name2 != 0
    key = _gwcat_key(rec.name1, rec.name2)
    for ev in recoils
        same_cm = rec.name_cm != 0 && ev.name_cm == rec.name_cm
        same_pair = named && _gwcat_key(ev.name1, ev.name2) == key
        (same_cm || same_pair) && in_window(ev.time_nb) && return :gw_recoil
    end
    for k in kicks
        t_kick = isfinite(k.time_nb) ? k.time_nb : k.time_myr / tstar
        hit =
            (rec.name1 != 0 && (k.name1 == rec.name1 || k.name2 == rec.name1)) ||
            (rec.name2 != 0 && (k.name1 == rec.name2 || k.name2 == rec.name2))
        hit && in_window(t_kick) && return :sn_kick
    end
    return :dynamical
end

"""
    escaped_binary_catalogue(log::BinaryEventLog; tstar_myr, rbar_pc, vstar_kms,
                             escapers = EscaperRecord[], t_now_myr = NaN,
                             band_hz = (1e-4, 1e-1), f_ref_hz = (1e-2, 10.0),
                             t_hubble_myr = 13_800.0) -> Vector{EscapedBinaryRow}

Escaped-binary catalogue of a run: one [`EscapedBinaryRow`](@ref) per
[`BinaryEscapeEvent`](@ref) whose two component types are compact (K* 10–14),
with the GW evolution of the pair after it left the cluster.

# Arguments
- `log`: event records of the engine log, from [`read_binary_events`](@ref)
- `tstar_myr`, `rbar_pc`, `vstar_kms`: time, length and velocity units T* [Myr], R̄ [pc], V* [km s⁻¹]
- `escapers`: esc.11 records ([`EscaperRecord`](@ref)), for `flag_kept_in_sim`
- `t_now_myr`: epoch of the present-epoch columns [Myr]; `NaN` leaves them `NaN`
- `band_hz`, `f_ref_hz`, `t_hubble_myr`: as in [`compact_binary_catalogue`](@ref)

The escape mode looks for a `:gw_recoil` [`GWEvent`](@ref) of the same c.m.
name or the same pair, then for a `:binary_kick` [`MassLossEvent`](@ref)
naming a component, in [t_esc − 1 NB, t_esc]. Rows are sorted by time.
Throws an `ArgumentError` for a non-positive unit or a band or reference
pair without 0 < lo < hi.

# Example
```julia
rows = escaped_binary_catalogue(log; tstar_myr = 0.87, rbar_pc = 2.5, vstar_kms = 2.6,
                                escapers = read_escapers(joinpath(dir, "esc.11")),
                                t_now_myr = 12_000.0)
write_escaped_binary_catalogue(joinpath(dir, "escaped_binary_catalogue.csv"), rows)
```
"""
function escaped_binary_catalogue(
    log::BinaryEventLog;
    tstar_myr::Real,
    rbar_pc::Real,
    vstar_kms::Real,
    escapers::AbstractVector{EscaperRecord} = EscaperRecord[],
    t_now_myr::Real = NaN,
    band_hz = _GWCAT_BAND_HZ,
    f_ref_hz = _GWCAT_FREF_HZ,
    t_hubble_myr::Real = _GWCAT_T_HUBBLE_MYR,
)::Vector{EscapedBinaryRow}
    caller = "escaped_binary_catalogue"
    _gwcat_check_positive(tstar_myr, "tstar_myr", caller)
    _gwcat_check_positive(rbar_pc, "rbar_pc", caller)
    _gwcat_check_positive(vstar_kms, "vstar_kms", caller)
    band = _gwcat_check_band(band_hz, "band_hz", caller)
    f_ref = _gwcat_check_band(f_ref_hz, "f_ref_hz", caller)
    tstar = Float64(tstar_myr)
    rbar = Float64(rbar_pc)
    vstar = Float64(vstar_kms)
    t_now = Float64(t_now_myr)
    t_hubble = Float64(t_hubble_myr)
    recoils = filter(ev -> ev.kind === :gw_recoil, log.gw)
    kicks = filter(ev -> ev.kind === :binary_kick, log.mass_loss)

    rows = EscapedBinaryRow[]
    for rec in log.escape
        _gwcat_compact(rec.kstar1, rec.kstar2) || continue
        n1, n2, k1, k2, m1, m2 =
            rec.name1, rec.name2, rec.kstar1, rec.kstar2, rec.m1_msun, rec.m2_msun
        if m2 > m1
            n1, n2 = n2, n1
            k1, k2 = k2, k1
            m1, m2 = m2, m1
        end
        a = rec.semi_nb * rbar * _GWCAT_AU_PER_PC
        e = rec.ecc
        t_esc_myr = rec.time_nb * tstar
        ms = _gwcat_mass_state(m1, m2)
        orb = _gwcat_orbit_state(m1, m2, a, e)
        ref = _gwcat_reference(m1, m2, a, e, f_ref)
        bnd = _gwcat_band(m1, m2, a, e, band)
        t_merge = t_esc_myr + orb.t_gw_peters_myr

        f_now = e_now = a_now = t_left = NaN
        if _gwcat_valid_orbit(m1, m2, a, e) &&
           isfinite(t_now) &&
           isfinite(t_esc_myr) &&
           t_now ≥ t_esc_myr
            a_now, e_now = eccentricity_after(m1, m2, a, e, (t_now - t_esc_myr) * _GWCAT_YR_PER_MYR)
            if isfinite(a_now)
                f_now = gw_frequency(m1, m2, a_now, e_now; definition = :f22)
                t_left = gw_inspiral_time(m1, m2, a_now, e_now) / _GWCAT_YR_PER_MYR
            end
        end
        listed = any(
            r ->
                (r.name == n1 || r.name == n2) &&
                abs(r.time_myr - t_esc_myr) ≤ _GWCAT_ESC_RTOL * abs(t_esc_myr),
            escapers,
        )
        push!(
            rows,
            EscapedBinaryRow(;
                t_esc_nb = rec.time_nb,
                t_esc_myr = t_esc_myr,
                name1 = n1,
                name2 = n2,
                name_cm = rec.name_cm,
                kstar1 = k1,
                kstar2 = k2,
                pair_class = pair_class(k1, k2),
                m1_msun = m1,
                m2_msun = m2,
                q = ms.q,
                mchirp_msun = ms.mchirp_msun,
                a_au = a,
                e = e,
                period_d = rec.period_d,
                eb_nb = rec.eb_nb,
                r_esc_pc = rec.r_nb * rbar,
                v_esc_kms = rec.v_nb * vstar,
                escape_mode = _gwcat_escape_mode(rec, recoils, kicks, tstar),
                t_gw_peters_myr = orb.t_gw_peters_myr,
                t_gw_mandel_myr = orb.t_gw_mandel_myr,
                q_f_zwick = orb.q_f_zwick,
                t_merge_abs_myr = t_merge,
                merges_before_t_hubble = t_merge ≤ t_hubble,
                e_at_fref_lo = ref.e_lo,
                e_at_fref_hi = ref.e_hi,
                t_band_in_myr = bnd.t_in_myr,
                t_band_out_myr = bnd.t_out_myr,
                f_gw22_now_hz = f_now,
                e_now = e_now,
                a_now_au = a_now,
                t_to_merge_from_now_myr = t_left,
                flag_kept_in_sim = !listed,
            ),
        )
    end
    sort!(rows; by = r -> (r.t_esc_nb, r.name1, r.name2))
    return rows
end

# ---------------------------------------------------------------------------
# GW-merger catalogue
# ---------------------------------------------------------------------------

"""
    GWMergerRow

One coalescence of two compact components (K* 10–14) in the GW-merger
catalogue, from [`gw_merger_catalogue`](@ref). Times in Myr unless the suffix
says otherwise; `NaN`, `0` and `:none` mark unavailable values. Names,
types and masses are in the order of the coalescence record. The trigger
state is the orbit of the last GW decay record (GR KSTIDE, GR KSINT2,
GR TIDES) of the pair at or before the merger.

| field | meaning | unit |
|---|---|---|
| `t_merge_nb`, `t_merge_myr` | time of the coalescence record | NB, Myr |
| `name1`, `name2` | component names | — |
| `kstar1`, `kstar2` | component types before the merger | — |
| `pair_class` | [`pair_class`](@ref) of the two types | — |
| `m1_msun`, `m2_msun` | component masses before the merger | M☉ |
| `m_remnant_msun`, `mass_lost_msun` | product mass and mass lost, as recorded | M☉ |
| `channel` | `channel` of the [`CoalescenceRecord`](@ref) | — |
| `a_trigger_au`, `e_trigger` | orbit of the trigger state | au, — |
| `q_peri_trigger_rs` | its pericentre a(1 − e) over 2GM/c² of the total mass | — |
| `trigger_kind` | `:q25rs` (pericentre within 25 Schwarzschild radii, or below SEPMRG), `:tgr_lt_torb` (TGR < T_orb), `:chain` (chain coalescence), `:unknown` | — |
| `v_kick_kms`, `spin1`, `spin2` | recoil speed and spins of the matching GW recoil line | km s⁻¹, — |
| `v_esc_kms` | escape speed of that line (printed with a KS index that may be stale) | km s⁻¹ |
| `retained` | `v_kick_kms < v_esc_kms`; `false` when either is `NaN` | — |
| `e_at_fref_hi` | e when f₂₂ of the trigger state reaches 10 Hz; `NaN` when already above | — |
| `generation` | 1 + the number of earlier mergers whose remnant is `name1` or `name2` | — |
| `span_primordial` | the covering [`BinarySpan`](@ref) is primordial | — |

The remnant of a merger keeps the `name1` of its coalescence record.
"""
Base.@kwdef struct GWMergerRow
    t_merge_nb::Float64
    t_merge_myr::Float64
    name1::Int
    name2::Int
    kstar1::Int
    kstar2::Int
    pair_class::Symbol
    m1_msun::Float64
    m2_msun::Float64
    m_remnant_msun::Float64
    mass_lost_msun::Float64
    channel::Symbol
    a_trigger_au::Float64
    e_trigger::Float64
    q_peri_trigger_rs::Float64
    trigger_kind::Symbol
    v_kick_kms::Float64
    spin1::Float64
    spin2::Float64
    v_esc_kms::Float64
    retained::Bool
    e_at_fref_hi::Float64
    generation::Int
    span_primordial::Bool
end

"""Trigger class of a merger from the last GW decay record of the pair (`nothing` without one)."""
function _gwcat_trigger_kind(
    ev::Union{Nothing,GWEvent},
    a_au::Float64,
    e::Float64,
    rs_au::Float64,
    channel::Symbol,
)::Symbol
    if ev !== nothing
        ev.qperi_nb < ev.sep_merge_nb && return :q25rs
        a_au * (1 - e) ≤ _GWCAT_TRIGGER_RS * rs_au && return :q25rs
        ev.t_gr_nb < ev.t_orb_nb && return :tgr_lt_torb
    end
    channel === :chain && return :chain
    return :unknown
end

"""
    gw_merger_catalogue(log::BinaryEventLog, coalescences::AbstractVector{CoalescenceRecord},
                        spans::AbstractVector{BinarySpan}; tstar_myr, rbar_pc, vstar_kms)
        -> Vector{GWMergerRow}

GW-merger catalogue of a run: one [`GWMergerRow`](@ref) per coal.24 record
whose two component types are compact (K* 10–14), sorted by time.

# Arguments
- `log`: event records of the engine log, from [`read_binary_events`](@ref)
- `coalescences`: coal.24 records, from [`read_coalescence_log`](@ref)
- `spans`: binary spans of the run, from [`binary_spans`](@ref)
- `tstar_myr`, `rbar_pc`, `vstar_kms`: time, length and velocity units T* [Myr], R̄ [pc], V* [km s⁻¹]

The trigger state is the latest GW decay [`GWEvent`](@ref) (`:gw_decay_step`,
`:gw_ksint`, `:gw_tides`) of the same unordered pair with
`time_nb ≤ t_merge_nb`; the recoil columns come from the `:gw_recoil` event of
the pair nearest in time. The covering span is the first span of the pair
with `t_form_nb ≤ t_merge_nb ≤ t_end_nb + 0.01 NB`, the tolerance of the F9.2
stamps. Generations count earlier rows in time order. Throws an
`ArgumentError` for a non-positive unit.

# Example
```julia
coal = read_coalescence_log(joinpath(dir, "coal.24"))
rows = gw_merger_catalogue(log, coal, spans; tstar_myr = 0.87, rbar_pc = 2.5, vstar_kms = 2.6)
write_gw_merger_catalogue(joinpath(dir, "gw_merger_catalogue.csv"), rows)
```
"""
function gw_merger_catalogue(
    log::BinaryEventLog,
    coalescences::AbstractVector{CoalescenceRecord},
    spans::AbstractVector{BinarySpan};
    tstar_myr::Real,
    rbar_pc::Real,
    vstar_kms::Real,
)::Vector{GWMergerRow}
    caller = "gw_merger_catalogue"
    _gwcat_check_positive(tstar_myr, "tstar_myr", caller)
    _gwcat_check_positive(rbar_pc, "rbar_pc", caller)
    _gwcat_check_positive(vstar_kms, "vstar_kms", caller)
    tstar = Float64(tstar_myr)
    au_per_nb = Float64(rbar_pc) * _GWCAT_AU_PER_PC
    span_index = _gwcat_span_index(spans)
    decays = _gwcat_gw_index(log, _GWCAT_PETERS_KINDS)
    recoils = _gwcat_gw_index(log, (:gw_recoil,))
    records = sort(
        filter(c -> _gwcat_compact(c.kstar1, c.kstar2), collect(coalescences));
        by = c -> c.time_nb,
    )

    remnants = Dict{Int,Int}()
    rows = GWMergerRow[]
    for c in records
        key = _gwcat_key(c.name1, c.name2)
        t = c.time_nb

        trigger = nothing
        for ev in get(decays, key, GWEvent[])
            (isfinite(ev.time_nb) && ev.time_nb ≤ t) || continue
            (trigger === nothing || ev.time_nb ≥ trigger.time_nb) && (trigger = ev)
        end
        a_trig = e_trig = q_rs = e_hi = NaN
        rs = NaN
        if trigger !== nothing
            a_trig = trigger.semi_nb * au_per_nb
            e_trig = trigger.ecc
            mt1, mt2 = trigger.m1_msun, trigger.m2_msun
            if !(isfinite(mt1) && isfinite(mt2))
                mt1, mt2 = c.mass1_msun, c.mass2_msun
            end
            rs = _gwcat_schwarzschild_au(mt1 + mt2)
            q_rs = a_trig * (1 - e_trig) / rs
            if _gwcat_valid_orbit(mt1, mt2, a_trig, e_trig)
                e_hi = _gwcat_time_to_frequency(mt1, mt2, a_trig, e_trig, _GWCAT_FREF_HZ[2], :f22).e
            end
        end

        recoil = nothing
        for ev in get(recoils, key, GWEvent[])
            isfinite(ev.time_nb) || continue
            (recoil === nothing || abs(ev.time_nb - t) < abs(recoil.time_nb - t)) && (recoil = ev)
        end
        v_kick = recoil === nothing ? NaN : recoil.v_kick_kms
        v_esc = recoil === nothing ? NaN : recoil.v_esc_kms

        span = _gwcat_covering_span(span_index, c.name1, c.name2, t, _SPAN_JOIN_TOL)
        generation =
            1 + get(remnants, c.name1, 0) + (c.name2 == c.name1 ? 0 : get(remnants, c.name2, 0))
        remnants[c.name1] = get(remnants, c.name1, 0) + 1

        push!(
            rows,
            GWMergerRow(;
                t_merge_nb = t,
                t_merge_myr = t * tstar,
                name1 = c.name1,
                name2 = c.name2,
                kstar1 = c.kstar1,
                kstar2 = c.kstar2,
                pair_class = pair_class(c.kstar1, c.kstar2),
                m1_msun = c.mass1_msun,
                m2_msun = c.mass2_msun,
                m_remnant_msun = c.mass_product_msun,
                mass_lost_msun = c.mass_lost_msun,
                channel = c.channel,
                a_trigger_au = a_trig,
                e_trigger = e_trig,
                q_peri_trigger_rs = q_rs,
                trigger_kind = _gwcat_trigger_kind(trigger, a_trig, e_trig, rs, c.channel),
                v_kick_kms = v_kick,
                spin1 = recoil === nothing ? NaN : recoil.spin1,
                spin2 = recoil === nothing ? NaN : recoil.spin2,
                v_esc_kms = v_esc,
                retained = v_kick < v_esc,
                e_at_fref_hi = e_hi,
                generation = generation,
                span_primordial = span !== nothing && span.primordial,
            ),
        )
    end
    return rows
end

# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

"""CSV cell: floats with 12 significant digits and `NaN` as `NaN`, Bools as `true`/`false`, Symbols by name."""
_gwcat_cell(v::AbstractFloat)::String = isnan(v) ? "NaN" : @sprintf("%.12g", v)
_gwcat_cell(v)::String = string(v)

"""Write `rows` of type `T` as CSV, one column per field in declaration order; returns the path."""
function _gwcat_write(path::AbstractString, rows::AbstractVector, ::Type{T})::String where {T}
    _backup_existing(path)
    n_fields = fieldcount(T)
    open(path, "w") do io
        println(io, join(fieldnames(T), ","))
        for r in rows
            println(io, join(ntuple(f -> _gwcat_cell(getfield(r, f)), n_fields), ","))
        end
    end
    return String(path)
end

"""
    write_compact_binary_catalogue(path, rows::AbstractVector{CompactBinaryRow}) -> String

Write the compact-binary catalogue as CSV: a header of the field names of
[`CompactBinaryRow`](@ref) in declaration order, one line per row. Floats
carry 12 significant digits (`NaN` written as `NaN`), Bools `true`/`false`,
Symbols their names. An existing file is backed up, never overwritten.
Returns the path.

# Example
```julia
write_compact_binary_catalogue("compact_binary_catalogue.csv", rows)
```
"""
write_compact_binary_catalogue(
    path::AbstractString,
    rows::AbstractVector{CompactBinaryRow},
)::String = _gwcat_write(path, rows, CompactBinaryRow)

"""
    write_escaped_binary_catalogue(path, rows::AbstractVector{EscapedBinaryRow}) -> String

Write the escaped-binary catalogue as CSV, one column per field of
[`EscapedBinaryRow`](@ref) in declaration order, with the cell conventions of
[`write_compact_binary_catalogue`](@ref). An existing file is backed up.
Returns the path.

# Example
```julia
write_escaped_binary_catalogue("escaped_binary_catalogue.csv", rows)
```
"""
write_escaped_binary_catalogue(
    path::AbstractString,
    rows::AbstractVector{EscapedBinaryRow},
)::String = _gwcat_write(path, rows, EscapedBinaryRow)

"""
    write_gw_merger_catalogue(path, rows::AbstractVector{GWMergerRow}) -> String

Write the GW-merger catalogue as CSV, one column per field of
[`GWMergerRow`](@ref) in declaration order, with the cell conventions of
[`write_compact_binary_catalogue`](@ref). An existing file is backed up.
Returns the path.

# Example
```julia
write_gw_merger_catalogue("gw_merger_catalogue.csv", rows)
```
"""
write_gw_merger_catalogue(path::AbstractString, rows::AbstractVector{GWMergerRow})::String =
    _gwcat_write(path, rows, GWMergerRow)
