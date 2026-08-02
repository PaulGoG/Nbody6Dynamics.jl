# =============================================================================
# Initial Mass Function (IMF) sampling
# =============================================================================

# Kroupa (2001) continuous, broken power-law representation used throughout.
# The breaks are fixed by the IMF model; the active segments are whatever
# subset overlaps the caller-provided [m_low, m_up] window.
const _KROUPA_BREAKS = (0.01, 0.08, 0.50, Inf)
const _KROUPA_ALPHAS = (0.3, 1.3, 2.3)   # ξ(m) ∝ m^{-α}

# Continuity coefficients c_i such that c_i m^{-α_i} matches at each break.
const _KROUPA_COEFFS = let
    c = ones(Float64, 3)
    c[2] = c[1] * _KROUPA_BREAKS[2]^(_KROUPA_ALPHAS[2] - _KROUPA_ALPHAS[1])
    c[3] = c[2] * _KROUPA_BREAKS[3]^(_KROUPA_ALPHAS[3] - _KROUPA_ALPHAS[2])
    (c[1], c[2], c[3])
end

# Indefinite integral of c m^{-α} from a to b.
@inline function _powlaw_integral(a::Real, b::Real, c::Real, α::Real)
    if abs(α - 1.0) < 1e-10
        return c * log(b / a)
    else
        return c * (b^(1.0 - α) - a^(1.0 - α)) / (1.0 - α)
    end
end

# Indefinite integral of m × c m^{-α} = c m^{1-α} from a to b (used for <m>).
@inline function _powlaw_mass_integral(a::Real, b::Real, c::Real, α::Real)
    if abs(α - 2.0) < 1e-10
        return c * log(b / a)
    else
        return c * (b^(2.0 - α) - a^(2.0 - α)) / (2.0 - α)
    end
end

# Return the list of (a, b, coeff, alpha) tuples that fall inside the
# requested [m_low, m_up] window, in increasing-m order.
function _active_kroupa_segments(m_low::Real, m_up::Real)
    m_low < m_up || error("IMF bounds must satisfy m_low < m_up, got [$m_low, $m_up]")
    segs = Tuple{Float64,Float64,Float64,Float64}[]
    for k in 1:3
        a = max(_KROUPA_BREAKS[k],   m_low)
        b = min(_KROUPA_BREAKS[k+1], m_up)
        a < b || continue
        push!(segs, (a, b, _KROUPA_COEFFS[k], _KROUPA_ALPHAS[k]))
    end
    isempty(segs) && error("No Kroupa IMF segments inside [$m_low, $m_up]")
    return segs
end

"""
    kroupa_mean_mass(m_low = 0.08, m_up = 100.0) -> Float64

Analytical mean stellar mass `⟨m⟩` of the Kroupa (2001) IMF restricted to
`[m_low, m_up]` in M☉. Used for:

- Sanity-checking over-determined configs (`M_total / N` should lie near
  this value for a natural stellar population).
- Estimating the expected cluster mass when only `N` is specified.

For the standard range [0.08, 100] this returns ≈ 0.58 M☉.
"""
function kroupa_mean_mass(m_low::Real = 0.08, m_up::Real = 100.0)::Float64
    segs = _active_kroupa_segments(Float64(m_low), Float64(m_up))
    num   = sum(_powlaw_mass_integral(a, b, c, α) for (a, b, c, α) in segs)
    denom = sum(_powlaw_integral(a, b, c, α)       for (a, b, c, α) in segs)
    return num / denom
end

"""
    sample_kroupa(N::Int; m_low=0.08, m_up=100.0,
                  rng=Random.default_rng()) -> Vector{Float64}

Draw `N` stellar masses from the Kroupa (2001) broken power-law IMF using
inverse-CDF sampling (no rejection, O(N)).

Only the segments of the Kroupa IMF that overlap `[m_low, m_up]` are active.
Masses are returned in M☉, un-normalised — the caller is responsible for
any post-sampling rescale.
"""
function sample_kroupa(N::Int; m_low::Float64 = 0.08, m_up::Float64 = 100.0,
                       rng::AbstractRNG = Random.default_rng())
    segments = _active_kroupa_segments(m_low, m_up)

    # Per-segment probability weights (normalised)
    weights = Float64[_powlaw_integral(a, b, c, α) for (a, b, c, α) in segments]
    total   = sum(weights)
    cum_weights = cumsum(weights) ./ total

    # Inverse CDF within a single power-law segment
    @inline function _sample_segment(a, b, c, α, u)
        if abs(α - 1.0) < 1e-10
            return a * exp(u * log(b / a))
        else
            I_ab = (b^(1.0 - α) - a^(1.0 - α))
            return (a^(1.0 - α) + u * I_ab)^(1.0 / (1.0 - α))
        end
    end

    masses = Vector{Float64}(undef, N)
    for i in 1:N
        r = rand(rng)
        seg_idx = 1
        for j in eachindex(cum_weights)
            if r ≤ cum_weights[j]
                seg_idx = j
                break
            end
        end
        u = rand(rng)
        a, b, c, α = segments[seg_idx]
        masses[i] = _sample_segment(a, b, c, α, u)
    end
    return masses
end

# -----------------------------------------------------------------------------
# Dispatch: sample_masses on IMFSpec
# -----------------------------------------------------------------------------

"""
    sample_masses(imf::IMFSpec, N::Int, rng::AbstractRNG) -> Vector{Float64}

Draw `N` body masses [M☉] according to the IMF specification. Dispatches on
the concrete `IMFSpec` subtype:

- `KroupaIMF`          — natural Kroupa sampling; no rescale.
- `RescaledKroupaIMF`  — Kroupa then uniform rescale to `target_mass`; emits
                         a `@warn` when the rescale factor is outside
                         ×[0.7, 1.4].
- `EqualMassIMF`       — every body gets `particle_mass`.
"""
sample_masses(imf::EqualMassIMF, N::Int, ::AbstractRNG) =
    fill(imf.particle_mass, N)

sample_masses(imf::KroupaIMF, N::Int, rng::AbstractRNG) =
    sample_kroupa(N; m_low = imf.bodyn, m_up = imf.body1, rng = rng)

function sample_masses(imf::RescaledKroupaIMF, N::Int, rng::AbstractRNG)
    raw = sample_kroupa(N; m_low = imf.bodyn, m_up = imf.body1, rng = rng)
    c = imf.target_mass / sum(raw)
    if abs(log10(c)) > log10(1.4)
        eff_lo = c * imf.bodyn
        eff_hi = c * imf.body1
        expected = N * kroupa_mean_mass(imf.bodyn, imf.body1)
        @warn """
        IMF rescale factor is ×$(round(c, digits=2)): the effective mass range shifts
        from [$(imf.bodyn), $(imf.body1)] M☉ to [$(round(eff_lo, digits=3)), $(round(eff_hi, digits=1))] M☉.

        This is "super-particle" mode — individual body masses no longer
        correspond to real stars, and downstream stellar-evolution output
        (SEV/BEV files, HR diagrams) will be non-physical.

        To get a natural Kroupa population, drop `mass_total` from the
        cluster spec (let it emerge as the sum of sampled masses, expected
        ≈ $(round(expected, digits=1)) M☉ for N=$N). If super-particles are
        intended, silence this warning by setting imf="kroupa_rescaled"
        explicitly.
        """ _group=:Nbody6Dynamics
    end
    raw .* c
end
