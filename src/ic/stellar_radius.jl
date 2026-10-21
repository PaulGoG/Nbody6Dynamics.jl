# =============================================================================
# Zero-age main-sequence radii of the engine's stellar evolution
# =============================================================================
# Tout et al. (1996), DOI 10.1093/mnras/281.1.257, with the coefficients of
# the engine: rzamsf (src/Main/zfuncs.f) evaluates a(8)…a(16) of COMMON
# /MSCFF/, which zcnsts (src/Main/zcnsts.f) sets from xz(36:76) of
# include/zdata.h as polynomials in lzs = log10(Z/0.02).

# xz(36:60) and xz(62:76): quartic coefficients of a(8)…a(12) and a(14)…a(16)
const _RZAMS_POLY = (
    (1.715359, 6.2246212e-01, -9.2557761e-01, -1.16996966, -3.0631491e-01),     # a(8)
    (6.597788, -4.2450044e-01, -1.213339427e+01, -1.073509484e+01, -2.51487077), # a(9)
    (1.008855000e+01, -7.11727086, -3.167119479e+01, -2.424848322e+01, -5.33608972), # a(10)
    (1.012495, 3.2699690e-01, -9.23418e-03, -3.876858e-02, -4.12750e-03),        # a(11)
    (7.490166e-02, 2.410413e-02, 7.233664e-02, 3.040467e-02, 1.97741e-03),       # a(12)
    (3.082234, 9.447205e-01, -2.15200882, -2.49219496, -6.3848738e-01),          # a(14)
    (1.784778e+01, -7.4534569, -4.896066856e+01, -4.005386135e+01, -9.09331816), # a(15)
    (2.2582e-04, -1.86899e-03, 3.88783e-03, 1.42402e-03, -7.671e-05),            # a(16)
)
# xz(61): a(13), independent of Z
const _RZAMS_A13 = 1.077422e-02

_horner4(c::NTuple{5,Float64}, x::Float64) = c[1] + x * (c[2] + x * (c[3] + x * (c[4] + x * c[5])))

"""
    _rzams_coefficients(z::Float64) -> NTuple{9,Float64}

The coefficients `a(8)…a(16)` that the engine's `zcnsts` sets for metallicity `z`.
"""
function _rzams_coefficients(z::Float64)::NTuple{9,Float64}
    lzs = log10(z / 0.02)
    p = map(c -> _horner4(c, lzs), _RZAMS_POLY)
    return (p[1], p[2], p[3], p[4], p[5], _RZAMS_A13, p[6], p[7], p[8])
end

"""
    zams_radius(m::Real, z::Real) -> Float64

Zero-age main-sequence radius [R☉] of a star of mass `m` [M☉] at metallicity
`z`: the rational fit of Tout et al. (1996, DOI 10.1093/mnras/281.1.257) with
the coefficients of the engine's SSE package (`rzamsf` after `zcnsts`), so
that the radii agree with those the engine assigns at `t = 0`.

# Arguments
- `m`: stellar mass [M☉]; finite and `> 0`. The fit is calibrated on
  0.1–100 M☉; above 100 M☉ it is an extrapolation.
- `z`: metallicity; `1e-4 ≤ z ≤ 0.03`, the range of the fit and of the engine.

Throws an `ArgumentError` outside these bounds.

# Example
```julia
zams_radius(1.0, 0.02)   # ≈ 0.888 R☉
```
"""
function zams_radius(m::Real, z::Real)::Float64
    (isfinite(m) && m > 0) ||
        throw(ArgumentError("zams_radius: the mass must be finite and > 0 [M☉]; got $m"))
    (1e-4 ≤ z ≤ 0.03) ||
        throw(ArgumentError("zams_radius: the metallicity must satisfy 1e-4 ≤ z ≤ 0.03; got $z"))
    mf = Float64(m)
    a8, a9, a10, a11, a12, a13, a14, a15, a16 = _rzams_coefficients(Float64(z))
    mx = sqrt(mf)
    num = (a8 * mf^2 + a9 * mf^6) * mx + a10 * mf^11 + (a11 + a12 * mx) * mf^19
    den = a13 + a14 * mf^2 + (a15 * mf^8 + mf^18 + a16 * mf^19) * mx
    return num / den
end
