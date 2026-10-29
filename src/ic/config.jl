# =============================================================================
# Merger IC configuration types and TOML parser
# =============================================================================
#
# Two tag hierarchies drive the IC generation: `DensityProfile` selects the
# spatial sampler (King or Plummer) and `IMFSpec` selects the mass-function
# sampler and its scaling mode. `ClusterSpec` carries the shared parameters
# (N, rbar, position, velocity).
#
# The IMF hierarchy is the load-bearing design decision:
#
#   - `KroupaIMF`           — natural Kroupa sampling in [bodyn, body1].
#                             Total mass is an *output* (= sum of samples).
#   - `RescaledKroupaIMF`   — Kroupa sampling, then a uniform linear rescale
#                             so the sum equals `target_mass`. This is the
#                             legacy "super-particle" mode where each body
#                             represents a mass-weighted lump rather than a
#                             real star. The effective mass range shifts by
#                             the rescale factor; stellar-evolution
#                             prescriptions are not physical in this regime.
#   - `EqualMassIMF`        — all bodies have the same mass.
#
# The TOML loader accepts two equivalent schemas: the flat
# form (`imf = "kroupa"`, `mass_total`, `bodyn`, `body1` at the cluster
# level) is the canonical user-facing one — all shipped configs use it —
# while the structured table form (`profile = {type=...}`, `imf =
# {type=...}`) exists for the lossless metadata round-trip of merger_ic.toml.
# The flat form with `mass_total` set translates to `RescaledKroupaIMF`,
# which warns at sample time when the rescale factor falls outside
# ×[0.7, 1.4].
# =============================================================================

# -----------------------------------------------------------------------------
# Density profile tag hierarchy
# -----------------------------------------------------------------------------

"Tag type for cluster density profile samplers."
abstract type DensityProfile end

"""
    KingProfile(; W0 = 6.0)

King (1966) model with dimensionless central potential `W0`. Larger `W0`
produces a more centrally concentrated cluster (and a larger ratio of tidal
to King radius).
"""
Base.@kwdef struct KingProfile <: DensityProfile
    W0::Float64 = 6.0
end

"""
    PlummerProfile()

Plummer (1911) isotropic model. No free parameters beyond the half-mass
radius carried on `ClusterSpec`.
"""
struct PlummerProfile <: DensityProfile end

"""
    RotatingKingProfile(; W0 = 6.0, omega0 = 0.6, axis = [0.0, 0.0, 1.0])

Rotating King model of Einsel & Spurzem (1999, DOI 10.1046/j.1365-8711.1999.02083.x):
central potential `W0`, dimensionless angular velocity `omega0 = Ω₀ r_K/σ ≥ 0`
(`0` is the King 1966 model), rotation axis `axis` (the angular momentum of the
cluster points along it; `[0, 0, 1]` is prograde with respect to the Galactic
rotation of `KZ(14) = 1, 2`, `[0, 0, -1]` retrograde). Solved and sampled by
[`solve_rotating_king`](@ref) and [`sample_rotating_king`](@ref).
"""
Base.@kwdef struct RotatingKingProfile <: DensityProfile
    W0::Float64 = 6.0
    omega0::Float64 = 0.6
    axis::Vector{Float64} = [0.0, 0.0, 1.0]
end

"Model name used in log messages and metadata TOML."
profile_name(::KingProfile) = "king"
profile_name(::PlummerProfile) = "plummer"
profile_name(::RotatingKingProfile) = "rotating_king"

# -----------------------------------------------------------------------------
# IMF tag hierarchy
# -----------------------------------------------------------------------------

"Tag type for stellar mass-function samplers. Each concrete subtype carries
its own parameters and dispatches `sample_masses`."
abstract type IMFSpec end

"""
    KroupaIMF(; bodyn = 0.08, body1 = 100.0)

Natural Kroupa (2001) sampling in `[bodyn, body1]` M☉. The total cluster
mass is an *output* (the sum of the N samples), not a specified parameter.

Use this mode when the scientific intent is a realistic stellar population.
Stellar-evolution prescriptions in Nbody6++ are well-defined for this case.
"""
Base.@kwdef struct KroupaIMF <: IMFSpec
    bodyn::Float64 = 0.08
    body1::Float64 = 100.0
end

"""
    RescaledKroupaIMF(; bodyn, body1, target_mass)

Super-particle mode: sample N stars from Kroupa in `[bodyn, body1]`, then
apply a uniform linear rescale so the sum equals `target_mass`. The IMF
shape (power-law slope) is preserved but the effective mass range shifts by
the rescale factor.

A warning is emitted at sample time when the rescale factor is outside
×[0.7, 1.4]: in that regime the individual body masses no longer correspond
to real stars and downstream stellar-evolution output is non-physical.

Use this mode deliberately when constructing super-particle simulations.
"""
Base.@kwdef struct RescaledKroupaIMF <: IMFSpec
    bodyn::Float64 = 0.08
    body1::Float64 = 100.0
    target_mass::Float64       # required
end

"""
    EqualMassIMF(; particle_mass)

Every body has identical mass `particle_mass` [M☉].
"""
Base.@kwdef struct EqualMassIMF <: IMFSpec
    particle_mass::Float64     # required
end

"""
    PowerLawIMF(; alpha = 2.35, bodyn = 0.08, body1 = 100.0)

Single power-law mass function `dN/dm ∝ m^{−alpha}` on `[bodyn, body1]` [M☉]
(`alpha = 2.35` is the Salpeter 1955 slope, DOI 10.1086/145971), drawn by the
inverse cumulative distribution; `alpha` is any finite real (the logarithmic
forms apply at 1 and 2).
"""
Base.@kwdef struct PowerLawIMF <: IMFSpec
    alpha::Float64 = 2.35
    bodyn::Float64 = 0.08
    body1::Float64 = 100.0
end

"Short label used in metadata TOML and logs."
imf_name(::KroupaIMF) = "kroupa"
imf_name(::RescaledKroupaIMF) = "kroupa_rescaled"
imf_name(::EqualMassIMF) = "equal"
imf_name(::PowerLawIMF) = "powerlaw"

# -----------------------------------------------------------------------------
# IMBHSpec
# -----------------------------------------------------------------------------

"""Values of `[merger.clusterK.imbh] placement`, `velocity` and `response`."""
const _IMBH_PLACEMENTS = ("centre", "offset")
const _IMBH_VELOCITIES = ("rest", "equipartition", "explicit")
const _IMBH_RESPONSES = ("rigid", "jeans")

"""Values of `[merger.nbody6] imbh_route`."""
const _IMBH_ROUTES = ("auto", "kz24", "flag")

"""
    IMBHSpec(; mass = NaN, mass_fraction = NaN, placement = "centre", offset = zeros(3),
               velocity = "rest", velocity_kms = zeros(3), response = "rigid", spin = 0.0)

Intermediate-mass black hole seeded in a cluster (`[merger.clusterK.imbh]`).

# Fields
- `mass::Float64`: IMBH mass M• [M☉]; exclusive with `mass_fraction`.
- `mass_fraction::Float64`: M•/M_stars of the cluster's systems,
  `0 < mass_fraction ≤ 0.05`; exclusive with `mass`.
- `placement::String`: `"centre"` (the stellar centre of mass) or `"offset"`
  (the centre of mass plus `offset`).
- `offset::Vector{Float64}`: displacement from the stellar centre of mass
  [pc], cluster frame; used by `"offset"`.
- `velocity::String`: `"rest"` (the cluster's mean velocity, nudged by
  10⁻⁶ σ_c in a random direction because the engine's start-up rejects a body
  at exact rest),
  `"equipartition"` (a Gaussian draw of per-component variance
  σ_c² ⟨m⟩/M• about it, Merritt 2005) or `"explicit"` (the mean velocity plus
  `velocity_kms`).
- `velocity_kms::Vector{Float64}`: velocity relative to the cluster [km s⁻¹];
  used by `"explicit"`.
- `response::String`: `"rigid"` (the stellar velocities unchanged) or
  `"jeans"` (rescaled in shells to the isotropic Jeans equilibrium with the
  added point mass).
- `spin::Float64`: dimensionless spin χ, `0 ≤ χ ≤ 1`; recorded only, neither
  route passes it to the engine.

One IMBH per cluster. The IMBHs are the last bodies of `dat.10`, in cluster
order. The `kz24` route (`KZ(24) = 1`) takes one IMBH per configuration and
repeats its state on the engine's `IMBHINIT` line; the `flag` route
(`KZ(24) = 3`, engine fork ≥ v2026.07.11, the package pin) takes 1 to 16 and
passes their `dat.10` names (`[merger.nbody6] imbh_route`).
"""
Base.@kwdef struct IMBHSpec
    mass::Float64 = NaN
    mass_fraction::Float64 = NaN
    placement::String = "centre"
    offset::Vector{Float64} = zeros(3)
    velocity::String = "rest"
    velocity_kms::Vector{Float64} = zeros(3)
    response::String = "rigid"
    spin::Float64 = 0.0
end

# -----------------------------------------------------------------------------
# ClusterSpec
# -----------------------------------------------------------------------------

"""
    ClusterSpec(; N, rbar, profile, imf, position = Float64[], velocity = Float64[])

Specification for a single star cluster.

# Fields
- `N::Int`: number of bodies to sample.
- `rbar::Float64`: target half-mass radius [pc].
- `profile::DensityProfile`: spatial sampler ([`KingProfile`](@ref),
  [`PlummerProfile`](@ref) or [`RotatingKingProfile`](@ref)).
- `imf::IMFSpec`: mass-function sampler ([`KroupaIMF`](@ref),
  [`RescaledKroupaIMF`](@ref), or [`EqualMassIMF`](@ref)).
- `position::Vector{Float64}`: centre-of-mass position [pc] (length 3). For
  `orbit_mode = "kepler"` and `"single"` this is left empty.
- `velocity::Vector{Float64}`: centre-of-mass velocity [km s⁻¹]. For
  `orbit_mode = "kepler"` and `"single"` this is left empty.
- `binaries::BinarySpec`: primordial binary population (none by default).
- `segregation::SegregationSpec`: primordial mass segregation of the systems
  (`SegregationSpec`; none by default).
- `imbh::Union{Nothing,IMBHSpec}`: seeded IMBH ([`IMBHSpec`](@ref)); `nothing` by default.

The keyword constructor takes the structured (tag-typed) fields only; the
flat TOML schema is the concern of [`load_merger_config`](@ref), which
translates it into these types. The seven-argument positional form of
earlier versions (without `segregation`) remains valid.
"""
struct ClusterSpec
    N::Int
    rbar::Float64
    profile::DensityProfile
    imf::IMFSpec
    position::Vector{Float64}
    velocity::Vector{Float64}
    binaries::BinarySpec
    segregation::SegregationSpec
    imbh::Union{Nothing,IMBHSpec}
end

function ClusterSpec(;
    N::Int = 50000,
    rbar::Real = 2.0,
    profile::DensityProfile = KingProfile(),
    imf::IMFSpec = KroupaIMF(),
    position::AbstractVector = Float64[],
    velocity::AbstractVector = Float64[],
    binaries::BinarySpec = BinarySpec(),
    segregation::SegregationSpec = SegregationSpec(),
    imbh::Union{Nothing,IMBHSpec} = nothing,
)
    return ClusterSpec(
        N,
        Float64(rbar),
        profile,
        imf,
        Float64.(position),
        Float64.(velocity),
        binaries,
        segregation,
        imbh,
    )
end

# The eight-field form of earlier versions (before `imbh`).
ClusterSpec(
    N::Integer,
    rbar::Real,
    profile::DensityProfile,
    imf::IMFSpec,
    position::AbstractVector,
    velocity::AbstractVector,
    binaries::BinarySpec,
    segregation::SegregationSpec,
) = ClusterSpec(N, rbar, profile, imf, position, velocity, binaries, segregation, nothing)

# The seven-field form of earlier versions (before `segregation`).
ClusterSpec(
    N::Integer,
    rbar::Real,
    profile::DensityProfile,
    imf::IMFSpec,
    position::AbstractVector,
    velocity::AbstractVector,
    binaries::BinarySpec,
) = ClusterSpec(N, rbar, profile, imf, position, velocity, binaries, SegregationSpec())

"""`BinarySpec` as a TOML table (metadata round trip); `massive` as a subtable."""
function _binaries_table(b::BinarySpec)
    d = Dict{String,Any}(
        String(k) => getfield(b, k) for k in fieldnames(BinarySpec) if k !== :massive
    )
    d["massive"] =
        Dict{String,Any}(String(k) => getfield(b.massive, k) for k in fieldnames(MassiveBinarySpec))
    return d
end

"""
`[<table>.binaries]` as a [`BinarySpec`](@ref), with the subtable
`[<table>.binaries.massive]` as its [`MassiveBinarySpec`](@ref); `table` is the
dotted cluster table (`merger.cluster<i>`, `merger.complex.member`). A
`binaries` value that is not a table is refused.
"""
function _parse_binaries_table(d, table::AbstractString)::BinarySpec
    section = "$table.binaries"
    d isa AbstractDict || _config_error(section, "must be a table; got $d ($(typeof(d)))")
    b = _parse_section(BinarySpec, d, section; skip = ("massive",))
    haskey(d, "massive") || return b
    raw = d["massive"]
    raw isa AbstractDict || _config_error("$section.massive", "must be a table")
    massive = _parse_section(MassiveBinarySpec, raw, "$section.massive")
    return BinarySpec(;
        fraction = b.fraction,
        pairing = b.pairing,
        period = b.period,
        a_min = b.a_min,
        a_max = b.a_max,
        q_min = b.q_min,
        eccentricity = b.eccentricity,
        roche_fill_max = b.roche_fill_max,
        massive = massive,
    )
end

const _MASSIVE_BINARY_MODELS = ("none", "sana2012", "moe2017")
const _MASSIVE_BINARY_PAIRINGS = ("pool", "replace")

"""
    _validate_binary_extensions(b::BinarySpec, key::String)

Fail-fast bounds of the Roche guard (`roche_fill_max`) and of the massive branch
(`massive`) of a [`BinarySpec`](@ref); `key` names the owning cluster table.
Raises an `ArgumentError` naming the offending key. `moe2017` additionally needs
`0.15 ≤ log_period_min` and `log_period_max ≤ 8.0`, and warns for `m_min < 0.8`.
"""
function _validate_binary_extensions(b::BinarySpec, key::String)
    kb = "$key.binaries"
    φ = b.roche_fill_max
    (φ == 0 || 0 < φ ≤ 1) || _config_error(
        "$kb.roche_fill_max",
        "must be 0 (no guard) or satisfy 0 < roche_fill_max ≤ 1; got $φ",
    )
    s = b.massive
    km = "$kb.massive"
    s.model in _MASSIVE_BINARY_MODELS || _config_error(
        "$km.model",
        "must be one of \"none\", \"sana2012\", \"moe2017\"; got \"$(s.model)\"",
    )
    s.pairing in _MASSIVE_BINARY_PAIRINGS ||
        _config_error("$km.pairing", "must be one of \"pool\", \"replace\"; got \"$(s.pairing)\"")
    (isfinite(s.m_min) && s.m_min > 0) ||
        _config_error("$km.m_min", "must be finite and > 0 [M☉]; got $(s.m_min)")
    (0 ≤ s.fraction ≤ 1) ||
        _config_error("$km.fraction", "must satisfy 0 ≤ fraction ≤ 1; got $(s.fraction)")
    (0 ≤ s.log_period_min < s.log_period_max) || _config_error(
        "$km.log_period_min",
        "must satisfy 0 ≤ log_period_min < log_period_max [log10 d]; " *
        "got $(s.log_period_min), $(s.log_period_max)",
    )
    s.log_period_max ≤ _KROUPA95_LOGP_MAX || _config_error(
        "$km.log_period_max",
        "must be ≤ $_KROUPA95_LOGP_MAX [log10 d]; got $(s.log_period_max)",
    )
    isfinite(s.period_exponent) ||
        _config_error("$km.period_exponent", "must be finite; got $(s.period_exponent)")
    (s.period_exponent > -1 || s.log_period_min > 0) || _config_error(
        "$km.period_exponent",
        "≤ −1 needs log_period_min > 0 (the law is not normalisable at 0); " *
        "got $(s.period_exponent) with log_period_min = $(s.log_period_min)",
    )
    (0 < s.q_min < 1) || _config_error("$km.q_min", "must satisfy 0 < q_min < 1; got $(s.q_min)")
    isfinite(s.q_exponent) || _config_error("$km.q_exponent", "must be finite; got $(s.q_exponent)")
    (isfinite(s.ecc_exponent) && s.ecc_exponent > -1) ||
        _config_error("$km.ecc_exponent", "must be finite and > −1; got $(s.ecc_exponent)")
    (0 < s.ecc_max < 1) ||
        _config_error("$km.ecc_max", "must satisfy 0 < ecc_max < 1; got $(s.ecc_max)")
    if s.model == "moe2017"
        s.log_period_min ≥ 0.15 || _config_error(
            "$km.log_period_min",
            "must be ≥ 0.15 [log10 d] for model = \"moe2017\"; got $(s.log_period_min)",
        )
        s.log_period_max ≤ 8.0 || _config_error(
            "$km.log_period_max",
            "must be ≤ 8.0 [log10 d] for model = \"moe2017\" (the domain of Moe & " *
            "Di Stefano 2017); got $(s.log_period_max)",
        )
        s.m_min < 0.8 && @warn "$km.m_min = $(s.m_min) M☉ lies below 0.8 M☉, the lower " *
              "bound of the moe2017 statistics; such primaries use those of 0.8 M☉"
    end
    return nothing
end

# -----------------------------------------------------------------------------
# Tag-type serialisation (metadata TOML round-trip)
# -----------------------------------------------------------------------------

"""
    _profile_table(p::DensityProfile) -> Dict{String,Any}

Serialise a density-profile tag to the structured TOML table form parsed by
[`_parse_profile_table`](@ref), so metadata files round-trip losslessly.
"""
_profile_table(p::KingProfile) = Dict{String,Any}("type" => "king", "W0" => p.W0)
_profile_table(::PlummerProfile) = Dict{String,Any}("type" => "plummer")
_profile_table(p::RotatingKingProfile) = Dict{String,Any}(
    "type" => "rotating_king",
    "W0" => p.W0,
    "omega0" => p.omega0,
    "axis" => copy(p.axis),
)

"""
    _imf_table(imf::IMFSpec) -> Dict{String,Any}

Serialise an IMF tag to the structured TOML table form parsed by
[`_parse_imf_table`](@ref), so metadata files round-trip losslessly.
"""
_imf_table(i::KroupaIMF) =
    Dict{String,Any}("type" => "kroupa", "bodyn" => i.bodyn, "body1" => i.body1)
_imf_table(i::RescaledKroupaIMF) = Dict{String,Any}(
    "type" => "kroupa_rescaled",
    "bodyn" => i.bodyn,
    "body1" => i.body1,
    "target_mass" => i.target_mass,
)
_imf_table(i::EqualMassIMF) =
    Dict{String,Any}("type" => "equal", "particle_mass" => i.particle_mass)

"""
    expected_mass(imf::IMFSpec, N::Int) -> Float64

Expected total cluster mass [M☉] for `N` bodies drawn from `imf` — exact for
`RescaledKroupaIMF`/`EqualMassIMF`, the analytic expectation for `KroupaIMF`.
"""
expected_mass(imf::KroupaIMF, N::Int) = N * kroupa_mean_mass(imf.bodyn, imf.body1)
expected_mass(imf::RescaledKroupaIMF, N::Int) = imf.target_mass
expected_mass(imf::EqualMassIMF, N::Int) = N * imf.particle_mass
expected_mass(imf::PowerLawIMF, N::Int) = N * powerlaw_mean_mass(imf.alpha, imf.bodyn, imf.body1)
_imf_table(i::PowerLawIMF) = Dict{String,Any}(
    "type" => "powerlaw",
    "alpha" => i.alpha,
    "bodyn" => i.bodyn,
    "body1" => i.body1,
)

# -----------------------------------------------------------------------------
# Orbit and output specs
# -----------------------------------------------------------------------------

"""
    OrbitSpec(; apocentre = 15.0, eccentricity = 0.7, apocentre_rh = 0.0, pericentre = 0.0,
                pericentre_rh = 0.0, orbit_normal = [0.0, 0.0, -1.0], apocentre_angle = 0.0)

Orbital parameters for the 2-cluster Kepler convenience mode. Ignored for
`orbit_mode = "explicit"` and `"single"`. The clusters start at apocentre.

# Fields
- `apocentre::Float64`: apocentre separation [pc]; used when none of the
  three keys below is positive.
- `eccentricity::Float64`: orbital eccentricity, `0 ≤ e < 1`.
- `apocentre_rh::Float64`: apocentre separation in units of the input
  half-mass radius `r_h,1` of cluster 1 (`0` = unused).
- `pericentre::Float64`: pericentre separation [pc] (`0` = unused); the
  apocentre is `r_p (1 + e) / (1 − e)`.
- `pericentre_rh::Float64`: pericentre separation in units of `r_h,1`
  (`0` = unused).
- `orbit_normal::Vector{Float64}`: direction of the orbital angular momentum
  of the pair (three components, any non-zero norm).
- `apocentre_angle::Float64`: rotation of the apocentre direction about
  `orbit_normal` [degrees].

At most one of the TOML keys `apocentre`, `apocentre_rh`, `pericentre`,
`pericentre_rh` may be given; none falls back to `apocentre = 15.0`.
[`orbit_apocentre`](@ref) resolves the apocentre in pc.

Orientation: with `n` the unit `orbit_normal`, the reference direction `u`
is the projection of `x̂` onto the plane normal to `n` (of `ŷ` when `n ∥ x̂`);
the separation vector cluster 2 − cluster 1 at the start points along
`e_apo = u cos φ + (n × u) sin φ`, `φ = apocentre_angle`, and the relative
velocity along `n × e_apo`, so the orbital angular momentum is along `n`.
The default `[0, 0, −1]` with `φ = 0` is the orbit of earlier versions: in
the x–y plane, cluster 2 on +x at apocentre and moving towards −y.
"""
Base.@kwdef struct OrbitSpec
    apocentre::Float64 = 15.0
    eccentricity::Float64 = 0.7
    apocentre_rh::Float64 = 0.0
    pericentre::Float64 = 0.0
    pericentre_rh::Float64 = 0.0
    orbit_normal::Vector{Float64} = [0.0, 0.0, -1.0]
    apocentre_angle::Float64 = 0.0
end

"""Keys of `[merger.orbit]` that set the orbit size; at most one may be given."""
const _ORBIT_SIZE_KEYS = ("apocentre", "apocentre_rh", "pericentre", "pericentre_rh")

"""
    orbit_apocentre(orbit::OrbitSpec, rh1::Real) -> Float64

Apocentre separation [pc] of a Kepler-mode orbit, with `rh1` the input
half-mass radius of cluster 1 [pc] and `e = orbit.eccentricity`:
`apocentre_rh > 0` gives `apocentre_rh rh1`; else `pericentre > 0` gives
`pericentre (1 + e) / (1 − e)`; else `pericentre_rh > 0` gives
`pericentre_rh rh1 (1 + e) / (1 − e)`; else `apocentre`. Throws an
`ArgumentError` for `rh1 ≤ 0`.

# Example
```julia
orbit_apocentre(OrbitSpec(; pericentre_rh = 2.0, eccentricity = 0.5), 2.0)   # 12.0
```
"""
function orbit_apocentre(orbit::OrbitSpec, rh1::Real)::Float64
    rh1 > 0 || throw(ArgumentError("orbit_apocentre: rh1 must be > 0 [pc]; got $rh1"))
    r_h = Float64(rh1)
    e = orbit.eccentricity
    orbit.apocentre_rh > 0 && return orbit.apocentre_rh * r_h
    orbit.pericentre > 0 && return orbit.pericentre * (1 + e) / (1 - e)
    orbit.pericentre_rh > 0 && return orbit.pericentre_rh * r_h * (1 + e) / (1 - e)
    return orbit.apocentre
end

"""
    MergerOutputSpec(; format = "nbody", truncate_jacobi = true,
                       output_dir = "", tcrit = 100.0, dtadj = 1.0, deltat = 1.0,
                       truncate_radius_rh = 0.0, truncate_rule = "roche")

Output and integration parameters for merger ICs.

# Fields
- `format::String`: `"nbody"` (KZ(22)=2, N-body units) — the only supported format.
- `truncate_jacobi::Bool`: truncate each cluster at its orbit-dependent
  tidal radius with respect to its nearest neighbour (`truncate_rule`)
  before combining.
- `truncate_radius_rh::Float64`: with `truncate_jacobi = true` and a
  positive value, every cluster is truncated at `truncate_radius_rh` times
  its own input half-mass radius instead of the orbit-dependent radius
  (kepler and explicit modes); `0` keeps the orbit-dependent radius; ≥ 0.
- `truncate_rule::String`: the orbit-dependent radius, `"roche"` (the Roche
  lobe of Eggleton 1983, [`roche_truncation_radius`](@ref), inside the inner
  Lagrange point at every mass ratio) or `"hill"` (the Hill radius
  [`jacobi_radius`](@ref), valid for a cluster much lighter than its
  neighbour; the rule of earlier versions).
- `output_dir::String`: directory for generated files; empty (the default)
  leaves the choice to the caller — `run_merger_pipeline` then creates a
  `runs/merger_<timestamp>/` under the working directory, and the main
  pipeline always writes into its run directory.
- `tcrit::Float64`: simulation end time (NB units when `format = "nbody"`).
- `dtadj::Float64`, `deltat::Float64`: adjustment and snapshot intervals.
- `tcrit_myr`, `dtadj_myr`, `deltat_myr`: the same three in Myr; a positive
  value replaces the NB one, converted at generation with the realised
  N-body time unit of the configuration (`0` = unused). A key may be given
  in one form only.

`dtadj` and `deltat` (and the stellar `dtplot`) are written as the nearest
dyadic rational with an exact decimal expansion ([`engine_interval`](@ref)):
the engine counts the decimal digits of these intervals with a loop that
never terminates on other values. For intervals ≥ 0.25 NB the change is
below 0.4 %; below that the nine-digit cap binds, the spacing stays 2⁻⁹ NB
and the change reaches up to 2⁻¹⁰ NB.
"""
Base.@kwdef struct MergerOutputSpec
    format::String = "nbody"
    truncate_jacobi::Bool = true
    output_dir::String = ""
    tcrit::Float64 = 100.0
    dtadj::Float64 = 1.0
    deltat::Float64 = 1.0
    tcrit_myr::Float64 = 0.0
    dtadj_myr::Float64 = 0.0
    deltat_myr::Float64 = 0.0
    truncate_radius_rh::Float64 = 0.0
    truncate_rule::String = "roche"
end

"""Integration-parameter presets of `[merger.nbody6]` (the `preset` key)."""
const _NBODY6_PRESETS = ("engine", "baseline")

"""Binary-physics option sets of `[merger.nbody6]` (the `binary_physics` key)."""
const _BINARY_PHYSICS = ("inherited", "full")

"""Former values of `preset` and `binary_physics`, read as their new names with a deprecation warning."""
const _NBODY6_RENAMED_VALUES =
    Dict((:preset, "campaign") => "baseline", (:binary_physics, "campaign") => "full")

"""
`(index, value)` pairs of the `KZ` options that `binary_physics = "full"` writes: chain and stable
hierarchies (`KZ(15) = 2`, chain records `KZ(30) = 2`), KS and chain slow-down (`KZ(26) = 2`), Peters
decay of compact KS pairs with sequential tides (`KZ(27) = 3`), GW braking of the inner compact binaries
of hierarchies (`KZ(28) = 1`), BSE Roche transfer with tidal synchronisation (`KZ(34) = 2`), the
hierarchy census (`KZ(18) = 1`), binary records also without primordial pairs (`KZ(8) = 2`), the binary
event files (`KZ(9) = 3`), the KS formation and termination lines (`KZ(10) = 2`) and no per-binary custom
output (`KZ(46) = 0`).
"""
const _FULL_BINARY_KZ = (
    (8, 2),
    (9, 3),
    (10, 2),
    (15, 2),
    (18, 1),
    (26, 2),
    (27, 3),
    (28, 1),
    (30, 2),
    (34, 2),
    (46, 0),
)

"""
    Nbody6ParameterSpec(; qe = 2.0e-4, etai = 0.0, etar = 0.0, etau = 0.0,
                          nnbopt = 0, rs0 = 0.0, rmin = 0.0, dtmin = 0.0,
                          eclose = 1.0, gmin = 1.0e-6, gmax = 0.01, smax = 1.0,
                          kz16 = 0, tcomp = 1.0e8, tcrtp0 = 3600.0,
                          isernb = 40, iserreg = 40, iserks = 0,
                          nfix = 1, ncrit = 10, nrun = 1, ncomm = 1, checkpoint = false,
                          preset = "engine", binary_physics = "inherited",
                          gw_recoil = "off", bh_spin = "zero", escapers = "remove",
                          imbh_route = "auto", imbh_output_dt_myr = 1.0,
                          post_newtonian = "off",
                          kz = Dict{Int,Int}())

Every numerical field of the `&INNBODY6` and `&ININPUT` namelists that the
merger input writer emits (`[merger.nbody6]`), in the N-body units of the
combined system. A zero for `etai`, `etar`, `etau`, `nnbopt`, `rs0`, `rmin`,
or `dtmin` means "the preset's value or rule, resolved at generation time"
([`resolve_nbody6_parameters`](@ref)). Meanings follow the Nbody6++GPU
input manual.

# Fields
- `qe`: energy-error tolerance per adjustment interval; > 0
- `etai`, `etar`, `etau`: irregular, regular, and regularised time-step
  factors; ≥ 0; 0 = preset value
- `nnbopt`: target neighbour number; `0` = `clamp(round(√N_total), 20, 300)`
- `rs0`: initial neighbour-sphere radius; `0` = derived; must not exceed
  the smallest member half-mass radius
- `rmin`, `dtmin`: KS regularisation distance and time step; `0` = derived
- `eclose`: binding energy per unit mass of a hard binary; > 0
- `gmin`, `gmax`: relative perturbation for unperturbed KS motion and the
  soft-binary termination parameter; `0 < gmin < gmax`
- `smax`: maximum time step (a power of two commensurate with 1); > 0
- `kz16`: `KZ(16)`, engine re-derivation of `RMIN`, `DTMIN`, `ECLOSE` every
  `DTADJ`; `0` keeps the written values (recommended for multi-cluster
  systems); 0–3
- `tcomp`: run-time limit [s]; > 0
- `tcrtp0`: termination time [Myr]; > 0
- `isernb`, `iserreg`, `iserks`: MPI block-size thresholds for serial
  execution of irregular, regular, and KS blocks; ≥ 0
- `nfix`: multiplier of `DELTAT` for `conf.3` and binary output; ≥ 1
- `ncrit`: minimum particle number, alternative termination criterion; ≥ 1
- `nrun`: run identification index; ≥ 1
- `ncomm`: the engine writes a restart (`COMMON`) dump on every `ncomm`-th
  dump call, counting the periodic and the terminating calls together; ≥ 1.
  Keep 1: with a larger value the dump of an `END RUN`, of a stop request
  or of an energy halt is written only when the count happens to fit
- `checkpoint`: `true` writes a restart dump at every adjustment
  (`KZ(2) = 1`), from which an interrupted run resumes; requires
  `ncomm = 1`. `false` (`KZ(2) = −1`) leaves the dump at the end of the run
  or at a stop request as the only one
- `preset`: `"engine"` or `"baseline"`; the source of every zero entry.
  `"engine"`: the rules the engine's own example inputs follow (η = 0.02,
  neighbour number √N, `RMIN` from the engine's `adjust.F` rule for the
  smallest member, `DTMIN` from the Aarseth close-encounter step).
  `"baseline"`: the integration baseline of the calibration programme — `ETAI` 0.01,
  `ETAR` 0.02, `ETAU` 0.1, `NNBOPT` 100, `RMIN` three times the 90°
  deflection scale 4 r_h/N of the smallest member, `DTMIN` 0.4 of the
  Aarseth close-encounter step at that `RMIN`; the engine's automatic
  re-derivation of these quantities (`KZ(16)`, `KZ(17)`, `KZ(32)`) is
  refused. Explicit nonzero keys override the preset in both cases.
- `binary_physics`: `"inherited"` or `"full"`; the few-body and binary
  options. `"inherited"` writes the block the generator has always written:
  no chain and no stable hierarchies (`KZ(15) = 0`), no gravitational-wave
  decay of compact pairs (`KZ(27) = 0`), no binary event output.
  `"full"` writes the complete set adopted for the binary programme (second table
  below, `_FULL_BINARY_KZ`) with the Level defaults of `&INSSE`, `&INBSE`
  and `&INCOLL`; it requires stellar evolution (`stellar.kz19 ≥ 3`), since
  the decay, tides and Roche transfer act on the stellar types, and an engine
  from the fork tag `v2026.07.6`: from `v2026.07-campaign.5` on the decay of
  compact KS pairs follows the Peters equations, and earlier tags end the run
  at the first binary–binary chain and can end a run with an unrestartable
  dump inside a chain. The former names `preset = "campaign"` and
  `binary_physics = "campaign"` are read as `"baseline"` and `"full"` with a
  deprecation warning
- `gw_recoil`: recoil of the remnant of a merger of two neutron stars or black
  holes (`&INBSE kicktype`, applied by COAL and MIX under `stellar.kz19 ≥ 3`):
  `"off"` (−2, the Level C default, no recoil), `"fit_q"` (0, the mass-ratio
  fit), `"nr_fits"` (1, the numerical-relativity fits with random spin
  orientations; the recoil direction of the engine is not isotropic) or
  `"test_1e4"` (−1, 10⁴ km s⁻¹, a test value); any value but `"off"` requires
  `stellar.kz19 ≥ 3`
- `bh_spin`: the dimensionless spin given to a black hole at formation
  (`&INBSE bhspinfl`): `"zero"` (0, the Level C default), `"uniform"` (1),
  `"gaussian"` (2, mean 0.5, width 0.2), `"maxwellian"` (3, σ = 0.2),
  `"geneva"` (4) or `"mesa"` (5, tracks in Z and the CO-core mass); the
  draws of 1–3 use an unseeded generator in the engine and are not
  restart-reproducible (warned)
- `escapers`: `"remove"` (`KZ(23) = 2`, the escapers beyond twice the tidal
  radius are removed, the behaviour the generator has always written) or
  `"keep"` (`KZ(23) = 0`, every star stays in the integration with its full
  phase space, for tails and the escaper catalogues; the cost grows with the
  escaped population)
- `imbh_route`: how the seeded IMBHs ([`IMBHSpec`](@ref)) reach the engine:
  `"kz24"` (`KZ(24) = 1`, the IMBH as the last body of `dat.10` repeated on
  the `IMBHINIT` line; one IMBH), `"flag"` (`KZ(24) = 3`, engine fork
  ≥ v2026.07.11: the IMBHs as the last bodies of `dat.10`, one per cluster
  and 1 to 16 in all, flagged by their names on the line after
  `NIMBHS DTBH`) or `"auto"` (`kz24` for one IMBH, `flag` for several);
  see [`resolved_imbh_route`](@ref)
- `imbh_output_dt_myr`: interval of the engine's IMBH record `mbh.45`
  (`DTBH`) [Myr]; > 0
- `post_newtonian`: `"off"` (`KZ(43) = 0`) or `"full"` (`KZ(43) = 7`, engine
  fork ≥ v2026.07.14, whose booking of the PN work closes the energy check;
  v2026.07.12 and .13 halt embedded runs): the 1PN, 2PN and 2.5PN relative
  accelerations of Blanchet (2024, DOI 10.1007/s41114-024-00050-z) in every
  perturbed KS step of a bound pair of two neutron stars or black holes
  (seeded IMBHs included) with a pericentre compactness `m/(c² q) ≥ 10⁻⁶`; such steps run
  without slow-down and replace the orbit-averaged Peters decay, which
  unperturbed pairs keep. `"full"` requires the Peters path, `KZ(27) = 3`
  (from `binary_physics = "full"` or a `kz` override), which also holds the
  merger test
- `kz`: explicit `KZ(i) = v` overrides applied last (`[merger.nbody6.kz]`,
  keys 1–50); an override of an index that also has a named key (8, 14, 16,
  19) warns. `KZ(18) ≥ 2` is refused: it runs the engine's primordial-triple
  generator on the `dat.10` pairs, and the generator writes no triples
  (`NHI0 = 0`). `KZ(46) > 0` (per-binary custom output, several GB per
  N-body time unit at 6×10⁵ bodies) warns

# Presets
| Key | `"engine"` | `"baseline"` |
|---|---|---|
| `etai` | 0.02 | 0.01 |
| `etar` | 0.02 | 0.02 |
| `etau` | 0.1 | 0.1 |
| `nnbopt` | `clamp(round(√N_total), 20, 300)` | 100 |
| `rs0` | `min(2 r_h (2 NNBOPT / N_min)^{1/3}, r_h)` | the same |
| `rmin` | `4 r_h / (N_min ρ̂^{1/3})` | `12 r_h / N_min` |
| `dtmin` | `0.04 √(ETAI/0.02) √(RMIN³ N_total)` | `0.016 √(ETAI/0.02) √(RMIN³ N_total)` |
| `kz16`, `KZ(17)`, `KZ(32)` | free | must be 0 |

`r_h` and `N_min` of the smallest member in N-body units of the combined
system, `ρ̂` its central density contrast, `N_total` the number of bodies.

# Binary physics
| `KZ` | `"inherited"` | `"full"` |
|---|---|---|
| 8 | 2 with primordial pairs, else 0 | 2 |
| 9, 10 | 0, 0 | 3, 2 |
| 15, 18 | 0, 0 | 2, 1 |
| 26, 27, 28 | 1, 0, 0 | 2, 3, 1 |
| 30, 34, 46 | 1, 0, 0 | 2, 2, 0 |
"""
Base.@kwdef struct Nbody6ParameterSpec
    qe::Float64 = 2.0e-4
    etai::Float64 = 0.0
    etar::Float64 = 0.0
    etau::Float64 = 0.0
    nnbopt::Int = 0
    rs0::Float64 = 0.0
    rmin::Float64 = 0.0
    dtmin::Float64 = 0.0
    eclose::Float64 = 1.0
    gmin::Float64 = 1.0e-6
    gmax::Float64 = 0.01
    smax::Float64 = 1.0
    kz16::Int = 0
    tcomp::Float64 = 1.0e8
    tcrtp0::Float64 = 3600.0
    isernb::Int = 40
    iserreg::Int = 40
    iserks::Int = 0
    nfix::Int = 1
    ncrit::Int = 10
    nrun::Int = 1
    ncomm::Int = 1
    checkpoint::Bool = false
    preset::String = "engine"
    binary_physics::String = "inherited"
    gw_recoil::String = "off"
    bh_spin::String = "zero"
    escapers::String = "remove"
    imbh_route::String = "auto"
    imbh_output_dt_myr::Float64 = 1.0
    post_newtonian::String = "off"
    kz::Dict{Int,Int} = Dict{Int,Int}()
end

"""`KZ(43)` per `post_newtonian` value (fork v2026.07.12 kspn.f: bits 1, 2, 4 = 1PN, 2PN, 2.5PN)."""
const _POST_NEWTONIAN_KZ43 = Dict("off" => 0, "full" => 7)

"""`kicktype` of `&INBSE` per `gw_recoil` value (kickGW.f: −2 none, 0 mass-ratio fit, 1 numerical-relativity fits with random spin orientations, −1 the 10⁴ km s⁻¹ test value)."""
const _GW_RECOIL_KICKTYPE = Dict("off" => -2, "fit_q" => 0, "nr_fits" => 1, "test_1e4" => -1)
"""`bhspinfl` of `&INBSE` per `bh_spin` value (kick.F: 0 zero spin, 1 uniform, 2 Gaussian (0.5, 0.2), 3 Maxwellian σ = 0.2, 4 Geneva and 5 MESA tracks of Z and M_CO)."""
const _BH_SPIN_FLAG = Dict(
    "zero" => 0,
    "uniform" => 1,
    "gaussian" => 2,
    "maxwellian" => 3,
    "geneva" => 4,
    "mesa" => 5,
)
"""`KZ(23)` per `escapers` value: `"remove"` removes escapers beyond twice the tidal radius (2), `"keep"` integrates every star (0)."""
const _ESCAPERS_KZ23 = Dict("remove" => 2, "keep" => 0)

"""
    StellarSpec(; kz19 = 3, level = "C", zmet = 0.001, epoch0 = 0.0, dtplot = 1.0,
                  dtplot_myr = 0.0)

Stellar-evolution settings of `merger.inp` (`[merger.stellar]`):

- `kz19`: `KZ(19)`, evolution and mass-loss scheme; `0` = off, `1`–`2`
  supernova schemes, `≥ 3` Eggleton–Tout–Hurley (SSE); ≥ 0
- `level`: SSE/BSE parameter level (Kamlah et al. 2022): `"A"`, `"B"`,
  `"C"`, or `"0"` (no level; the engine's independent defaults)
- `zmet`: metal abundance; `0.0001 ≤ zmet ≤ 0.03` (the engine's own bounds)
- `epoch0`: formation time of the population [Myr]; ≤ 0 (age at start is
  `−epoch0`)
- `dtplot`: interval of the stellar-evolution diagnostics [NB]; > 0 and,
  per the manual, ≥ `deltat`
- `dtplot_myr`: the same in Myr; a positive value replaces `dtplot` at
  generation, converted with the realised time unit
"""
Base.@kwdef struct StellarSpec
    kz19::Int = 3
    level::String = "C"
    zmet::Float64 = 0.001
    epoch0::Float64 = 0.0
    dtplot::Float64 = 1.0
    dtplot_myr::Float64 = 0.0
end

"""
    TidalSpec(; kz14 = 0, gmg = 0.0, rg0 = 0.0, rg = zeros(3), vg = zeros(3),
              velocity_frame = "inertial")

External galactic field of `merger.inp` (`[merger.tidal]`, `KZ(14)` and the
`&INXTRNL0` namelist):

- `kz14`: `0` isolated; `1` standard solar-neighbourhood linearised tide
  (no parameters); `2` point-mass galaxy on a circular orbit (`gmg`,
  `rg0`); `5` MWPotential2014 with the cluster on a galactocentric orbit
  (`rg`, `vg`). Options `3` (point mass + disk + halo + bulge) and `4`
  (Plummer potential) are refused: the engine rescales every velocity to
  the `&INSCALE` virial ratio including the external potential on those
  paths (`xtrnl0.F`), which would destroy the prescribed orbital
  kinematics of a multi-cluster configuration.
- `gmg`: galaxy mass [M☉] (`kz14 = 2`); > 0
- `rg0`: galactocentric distance of the circular orbit [kpc] (`kz14 = 2`); > 0
- `rg`, `vg`: galactocentric position [kpc] and velocity [km s⁻¹] of the
  configuration's centre of mass (`kz14 = 5`); both non-zero
- `velocity_frame`: `"rotating"` or `"inertial"`. Under `kz14 = 1, 2` the
  engine integrates in the frame that rotates with the circular Galactic
  orbit, Ω = Ω ẑ (Ω = A − B = 26.4 km s⁻¹ kpc⁻¹ for `kz14 = 1`, V_G/R_G
  for `kz14 = 2`), and takes the `dat.10` velocities as velocities in that
  frame (`xtrnl0.F`, `xtrnld.f`). The generator builds the clusters and the
  orbit in an inertial frame, so `"inertial"` (the default) converts them,
  writing v − Ω × r for every body. `"rotating"` writes the generated
  velocities as they are (the behaviour of earlier versions): clusters
  non-rotating and a Kepler orbit of the pair in the rotating frame. Without
  a rotating frame (`kz14 = 0, 5`) the value has no effect; the TOML key is
  then refused at load, a programmatic `TidalSpec` may carry either value.

The `&INSCALE` tidal radius stays 0 so the engine derives it from the field
with the generator's `RBAR` (a non-zero value would override `RBAR`).
"""
Base.@kwdef struct TidalSpec
    kz14::Int = 0
    gmg::Float64 = 0.0
    rg0::Float64 = 0.0
    rg::Vector{Float64} = zeros(3)
    vg::Vector{Float64} = zeros(3)
    velocity_frame::String = "inertial"
end

# -----------------------------------------------------------------------------
# MergerConfig
# -----------------------------------------------------------------------------

"""
    MBHSpec(; mass_msun, softening_pc = 0.0, rt_mode = "fixed", r_acc_pc = 0.0,
              partial = "off", beta_full = 1.0, accretion_fraction = 0.0,
              response = "jeans")

Static massive black hole at the coordinate origin (`[merger.mbh]`; engine
fork ≥ v2026.07.13, `KZ(44) = 1` and the `&INMBH` namelist). The black hole
is a fixed point mass acting on every body (irregular force, KS
perturbations as a tidal difference); a single star whose orbit about it
reaches the capture radius is disrupted at its pericentre passage, its energy
booked and its mass removed (engine line `TDE`).

# Fields
- `mass_msun`: black-hole mass [M☉]; > 0
- `softening_pc`: Plummer softening of its force [pc]; ≥ 0 (0: a point mass)
- `rt_mode`: `"fixed"` (every star disrupted at `r_acc_pc`, the practice of
  Cho et al. 2026a) or `"stellar"` (r_t = R★ (M/m★)^{1/3} with the
  stellar-evolution radius; requires `stellar.kz19 ≥ 3`)
- `r_acc_pc`: disruption radius of `"fixed"` [pc]; > 0 there
- `partial`: `"off"` (full disruption when r_t/q ≥ `beta_full`), `"gamma43"`
  or `"gamma53"` (full disruption above β = 1.85 or 0.9 and partial
  disruption down to 0.6 or 0.5, the star losing the fraction of Guillochon &
  Ramirez-Ruiz 2013, DOI 10.1088/0004-637X/767/1/25)
- `beta_full`: threshold of `partial = "off"`; > 0
- `accretion_fraction`: fraction of every disrupted or stripped mass the black
  hole gains (`FACC`; a plunging star is gained whole when it is > 0;
  engine fork ≥ v2026.07.17); 0 ≤ value ≤ 1, 0 keeps the mass fixed
- `response`: `"jeans"` (the stellar velocities rescaled by the isotropic
  Jeans equation with the black hole's potential, as for a seeded IMBH) or
  `"rigid"` (velocities unchanged)

Requires `orbit_mode = "single"` (the cluster at the origin) and no external
field (`tidal.kz14 = 0`).
"""
Base.@kwdef struct MBHSpec
    mass_msun::Float64 = NaN
    softening_pc::Float64 = 0.0
    rt_mode::String = "fixed"
    r_acc_pc::Float64 = 0.0
    partial::String = "off"
    beta_full::Float64 = 1.0
    accretion_fraction::Float64 = 0.0
    response::String = "jeans"
end

"""`IPART` of `&INMBH` per `MBHSpec.partial` value."""
const _MBH_IPART = Dict("off" => 0, "gamma43" => 1, "gamma53" => 2)

"""
    MergerConfig(clusters, orbit_mode, orbit, output[, nbody6, stellar, tidal, seed, virial_max_n,
                 seed_streams, complex])

Top-level configuration for multi-cluster merger initial conditions.

# Orbit modes

- `"kepler"` — exactly 2 clusters placed on a Keplerian orbit at apocentre.
  Requires an [`OrbitSpec`](@ref). Cluster positions/velocities are computed
  automatically.
- `"explicit"` — any number of clusters (≥ 2). Each cluster must specify
  `position = [x, y, z]` and `velocity = [vx, vy, vz]` in its TOML block.
  The `orbit` section is ignored (with a warning if it appears in TOML).
- `"single"` — one cluster at the origin, at rest (`n_clusters = 1`), without
  truncation. `position` and `velocity` are refused; the `orbit` section is
  ignored (with a warning if it appears in TOML).
- `"complex"` — a cluster complex built from a `ComplexSpec`: member
  masses from a mass function, half-mass radii from a size–mass relation,
  centres from a distribution and the bulk velocities scaled to a target
  virial ratio of the members as softened point masses. `clusters` holds the
  member template (`[merger.complex.member]`, `N` and `rbar` placeholders);
  the members are assembled as in explicit mode. `n_clusters` and
  `[merger.clusterN]` tables are refused; the `orbit` section is ignored
  (with a warning if it appears in TOML).

# Fields
- `clusters::Vector{ClusterSpec}`
- `orbit_mode::String`: `"kepler"`, `"explicit"`, `"single"` or `"complex"`
- `orbit::OrbitSpec`
- `output::MergerOutputSpec`
- `nbody6::Nbody6ParameterSpec`: integration and run-control parameters
- `stellar::StellarSpec`: stellar-evolution settings
- `tidal::TidalSpec`: external galactic field
- `seed::Union{Int, Nothing}`: RNG seed. `nothing` = non-deterministic.
- `virial_max_n::Int`: largest N, per cluster and for the combined system, for which the exact
  O(N²) potential is evaluated (virialisation of each cluster, combined virial ratio); default
  `_VIRIAL_NMAX` (200 000)
- `seed_streams::String`: `"shared"` (default) draws every cluster and its
  pairs from one `MersenneTwister(seed)` stream; `"per_cluster"` gives cluster
  `i` its own stream seeded by `(seed, i)`, which samples and expands its
  pairs, so a cluster's realisation does not depend on the other clusters or
  on the orbit
- `complex::ComplexSpec`: parameters of `orbit_mode = "complex"`; default
  `ComplexSpec()`, unused in the other modes
- `mbh::Union{Nothing,MBHSpec}`: static massive black hole at the origin
  (`[merger.mbh]`, [`MBHSpec`](@ref)); `nothing` by default
"""
struct MergerConfig
    clusters::Vector{ClusterSpec}
    orbit_mode::String
    orbit::OrbitSpec
    output::MergerOutputSpec
    nbody6::Nbody6ParameterSpec
    stellar::StellarSpec
    tidal::TidalSpec
    seed::Union{Int,Nothing}
    virial_max_n::Int
    seed_streams::String
    complex::ComplexSpec
    mbh::Union{Nothing,MBHSpec}
end

# The eleven-field form of earlier versions (before `mbh`).
MergerConfig(
    clusters::Vector{ClusterSpec},
    orbit_mode::String,
    orbit::OrbitSpec,
    output::MergerOutputSpec,
    nbody6::Nbody6ParameterSpec,
    stellar::StellarSpec,
    tidal::TidalSpec,
    seed::Union{Int,Nothing},
    virial_max_n::Int,
    seed_streams::String,
    complex::ComplexSpec,
) = MergerConfig(
    clusters,
    orbit_mode,
    orbit,
    output,
    nbody6,
    stellar,
    tidal,
    seed,
    virial_max_n,
    seed_streams,
    complex,
    nothing,
)

"""Values of `[merger] seed_streams`."""
const _SEED_STREAMS = ("shared", "per_cluster")

# Convenience form: derived integration parameters, default stellar
# evolution, and a non-deterministic seed (the drawn seed is still recorded
# in the run metadata).
MergerConfig(
    clusters::Vector{ClusterSpec},
    orbit_mode::AbstractString,
    orbit::OrbitSpec,
    output::MergerOutputSpec;
    nbody6::Nbody6ParameterSpec = Nbody6ParameterSpec(),
    stellar::StellarSpec = StellarSpec(),
    tidal::TidalSpec = TidalSpec(),
    seed::Union{Int,Nothing} = nothing,
    virial_max_n::Int = _VIRIAL_NMAX,
    seed_streams::AbstractString = "shared",
    complex::ComplexSpec = ComplexSpec(),
    mbh::Union{Nothing,MBHSpec} = nothing,
) = MergerConfig(
    clusters,
    String(orbit_mode),
    orbit,
    output,
    nbody6,
    stellar,
    tidal,
    seed,
    virial_max_n,
    String(seed_streams),
    complex,
    mbh,
)

MergerConfig(
    clusters::Vector{ClusterSpec},
    orbit_mode::AbstractString,
    orbit::OrbitSpec,
    output::MergerOutputSpec,
    nbody6::Nbody6ParameterSpec,
    stellar::StellarSpec,
    tidal::TidalSpec,
    seed::Union{Int,Nothing},
    seed_streams::AbstractString = "shared",
) = MergerConfig(
    clusters,
    String(orbit_mode),
    orbit,
    output,
    nbody6,
    stellar,
    tidal,
    seed,
    _VIRIAL_NMAX,
    String(seed_streams),
    ComplexSpec(),
)

# The ten-field form of earlier versions (before `complex`).
MergerConfig(
    clusters::Vector{ClusterSpec},
    orbit_mode::AbstractString,
    orbit::OrbitSpec,
    output::MergerOutputSpec,
    nbody6::Nbody6ParameterSpec,
    stellar::StellarSpec,
    tidal::TidalSpec,
    seed::Union{Int,Nothing},
    virial_max_n::Int,
    seed_streams::AbstractString,
) = MergerConfig(
    clusters,
    String(orbit_mode),
    orbit,
    output,
    nbody6,
    stellar,
    tidal,
    seed,
    virial_max_n,
    String(seed_streams),
    ComplexSpec(),
)

# The nine-field form of earlier versions (before `seed_streams`).
MergerConfig(
    clusters::Vector{ClusterSpec},
    orbit_mode::AbstractString,
    orbit::OrbitSpec,
    output::MergerOutputSpec,
    nbody6::Nbody6ParameterSpec,
    stellar::StellarSpec,
    tidal::TidalSpec,
    seed::Union{Int,Nothing},
    virial_max_n::Int,
) = MergerConfig(
    clusters,
    String(orbit_mode),
    orbit,
    output,
    nbody6,
    stellar,
    tidal,
    seed,
    virial_max_n,
    "shared",
    ComplexSpec(),
)

# -----------------------------------------------------------------------------
# MergerICResult
# -----------------------------------------------------------------------------

"""
    MergerICResult

Structured output from [`generate_merger_ic`](@ref), carrying everything
needed for downstream plotting and inspection without re-reading files.

# Fields
- `output_dir::String`
- `N_total::Int`: total particle count (after optional Jacobi truncation)
- `M_total::Float64`: total mass [M☉]
- `rbar::Float64`: half-mass radius [pc]
- `zmbar::Float64`: mean particle mass [M☉]
- `cluster_ranges::Vector{Vector{Int}}`: body indices of every cluster in
  `dat.10` order (a contiguous range without primordial binaries; with them,
  the cluster's pairs block followed by its singles block)
- `n_pairs::Vector{Int}`: primordial pairs per cluster (`sum` = `NBIN0`)
- `cluster_specs::Vector{ClusterSpec}`
- `orbit_mode::String`
- `orbit_spec::OrbitSpec`
- `mass_physical::Vector{Float64}`: [M☉]
- `pos_physical::Matrix{Float64}`: 3×N [pc]
- `vel_physical::Matrix{Float64}`: 3×N [km/s]
"""
struct MergerICResult
    output_dir::String
    N_total::Int
    M_total::Float64
    rbar::Float64
    zmbar::Float64
    cluster_ranges::Vector{Vector{Int}}
    n_pairs::Vector{Int}
    cluster_specs::Vector{ClusterSpec}
    orbit_mode::String
    orbit_spec::OrbitSpec
    mass_physical::Vector{Float64}
    pos_physical::Matrix{Float64}
    vel_physical::Matrix{Float64}
end

# -----------------------------------------------------------------------------
# TOML loading
# -----------------------------------------------------------------------------

"""
    load_merger_config(path::AbstractString) -> MergerConfig

Load a merger IC configuration from a TOML file. Accepts the canonical flat
schema and the structured metadata schema interchangeably. Keys the schema
does not define are rejected, and parsed values are validated fail-fast
(N ≥ 2, rbar > 0, W0 > 0, IMF mass bounds, Kepler orbit parameters, positive
integration intervals); violations raise an `ArgumentError` naming the
offending `merger.*` key.

# Kepler mode (2 clusters, flat form)

```toml
[merger]
n_clusters = 2
orbit_mode = "kepler"

[merger.cluster1]
model = "king"
N = 50000
rbar = 2.0
W0 = 6.0
imf = "kroupa"
bodyn = 0.08
body1 = 100.0
# mass_total optional: present → super-particle rescale; absent → natural Kroupa

[merger.orbit]
apocentre = 15.0
eccentricity = 0.7
```

# Explicit mode (N clusters)

Each cluster provides `position = [x, y, z]` and `velocity = [vx, vy, vz]`.

# Single mode (1 cluster)

`n_clusters = 1`, `orbit_mode = "single"`: the cluster is placed at the
origin, at rest; `position` and `velocity` must be absent.

# Structured form (equivalent to flat; used by the merger_ic.toml metadata round-trip)

```toml
[merger.cluster1]
N       = 1500
rbar    = 1.0
profile = { type = "king", W0 = 6.0 }
imf     = { type = "kroupa", bounds = [0.08, 100.0] }
```
"""
function load_merger_config(path::AbstractString)::MergerConfig
    raw = TOML.parsefile(path)
    _reject_unknown(raw, ("merger",), basename(path))
    haskey(raw, "merger") || _config_error("merger", "table is missing in $path")
    m = raw["merger"]
    m isa AbstractDict || _config_error("merger", "must be a table; got $m ($(typeof(m)))")

    n_clusters = _typed(get(m, "n_clusters", 2), Int, "merger.n_clusters")
    orbit_mode = _typed(get(m, "orbit_mode", "kepler"), String, "merger.orbit_mode")
    orbit_mode in ("kepler", "explicit", "single", "complex") || _config_error(
        "merger.orbit_mode",
        "must be \"kepler\", \"explicit\", \"single\" or \"complex\", got \"$orbit_mode\"",
    )
    complex_mode = orbit_mode == "complex"
    if complex_mode
        haskey(m, "n_clusters") && _config_error(
            "merger.n_clusters",
            "is not used with orbit_mode = \"complex\" (the members come from [merger.complex]); remove it",
        )
        for k in sort!(collect(String, keys(m)))
            occursin(r"^cluster\d+$", k) && _config_error(
                "merger.$k",
                "is not used with orbit_mode = \"complex\" (the members are built from " *
                "[merger.complex.member]); remove it",
            )
        end
        haskey(m, "complex") ||
            _config_error("merger.complex", "table is missing in $path (orbit_mode = \"complex\")")
        n_clusters = 0
    else
        haskey(m, "complex") && _config_error(
            "merger.complex",
            "is used with orbit_mode = \"complex\" only; got orbit_mode = \"$orbit_mode\"",
        )
    end
    virial_max_n = _typed(get(m, "virial_max_n", _VIRIAL_NMAX), Int, "merger.virial_max_n")
    virial_max_n ≥ 1 || _config_error("merger.virial_max_n", "must be ≥ 1; got $virial_max_n")
    seed_streams = _typed(get(m, "seed_streams", "shared"), String, "merger.seed_streams")
    seed_streams in _SEED_STREAMS || _config_error(
        "merger.seed_streams",
        "must be \"shared\" or \"per_cluster\", got \"$seed_streams\"",
    )

    cluster_keys = ["cluster$k" for k in 1:max(n_clusters, 0)]
    for k in sort!(collect(String, keys(m)))
        (occursin(r"^cluster\d+$", k) && !(k in cluster_keys)) &&
            _config_error("merger.$k", "is given but n_clusters = $n_clusters")
    end
    _reject_unknown(
        m,
        vcat(
            [
                "n_clusters",
                "orbit_mode",
                "seed",
                "seed_streams",
                "virial_max_n",
                "orbit",
                "output",
                "nbody6",
                "stellar",
                "tidal",
                "complex",
                "mbh",
            ],
            cluster_keys,
        ),
        "merger",
    )

    clusters = ClusterSpec[]
    for i in 1:n_clusters
        key = "cluster$i"
        haskey(m, key) || _config_error("merger.$key", "table is missing in $path")
        push!(clusters, _parse_cluster_table(_section_table(m, key, "merger.$key"), i, orbit_mode))
    end
    # The member template joins `clusters` after the per-cluster checks below.
    complex_spec, complex_template = if complex_mode
        _parse_complex_section(m["complex"])
    else
        (ComplexSpec(), nothing)
    end

    if orbit_mode == "kepler"
        n_clusters == 2 || _config_error(
            "merger.n_clusters",
            "must be exactly 2 for orbit_mode = \"kepler\", got $n_clusters",
        )
    end
    if orbit_mode == "single"
        n_clusters == 1 || _config_error(
            "merger.n_clusters",
            "must be exactly 1 for orbit_mode = \"single\", got $n_clusters",
        )
    end
    if orbit_mode == "explicit"
        n_clusters ≥ 2 || _config_error(
            "merger.n_clusters",
            "must be ≥ 2 for orbit_mode = \"explicit\", got $n_clusters",
        )
    end

    orbit_raw = _section_table(m, "orbit", "merger.orbit")
    if orbit_mode in ("explicit", "single", "complex") && !isempty(orbit_raw)
        @warn "[merger.orbit] is defined but orbit_mode=\"$orbit_mode\" — the orbit section will be ignored."
    end
    orbit = _parse_section(OrbitSpec, orbit_raw, "merger.orbit")
    size_keys = String[k for k in _ORBIT_SIZE_KEYS if haskey(orbit_raw, k)]
    length(size_keys) ≤ 1 || _config_error(
        "merger.orbit",
        "give one of apocentre, apocentre_rh, pericentre, pericentre_rh; got " *
        join(size_keys, ", "),
    )

    out_raw = _section_table(m, "output", "merger.output")
    output = _parse_section(MergerOutputSpec, out_raw, "merger.output")
    for key in ("tcrit", "dtadj", "deltat")
        (haskey(out_raw, key) && haskey(out_raw, key * "_myr")) && _config_error(
            "merger.output.$key",
            "and merger.output.$(key)_myr are mutually exclusive",
        )
    end

    nb_raw = _section_table(m, "nbody6", "merger.nbody6")
    nb = _parse_section(Nbody6ParameterSpec, nb_raw, "merger.nbody6"; skip = ("kz",))
    kz_raw = _section_table(nb_raw, "kz", "merger.nbody6.kz")
    kz_over = Dict{Int,Int}()
    for (k, v) in kz_raw
        idx = tryparse(Int, String(k))
        idx === nothing &&
            _config_error("merger.nbody6.kz", "keys must be integers 1–50; got \"$k\"")
        v isa Integer || _config_error("merger.nbody6.kz[\"$k\"]", "must be an integer; got $v")
        kz_over[idx] = Int(v)
    end
    renamed = Dict{Symbol,String}()
    for key in (:preset, :binary_physics)
        value = getfield(nb, key)
        new = get(_NBODY6_RENAMED_VALUES, (key, value), nothing)
        new === nothing && continue
        @warn "merger.nbody6.$key = \"$value\" is deprecated; read as \"$new\""
        renamed[key] = new
    end
    nbody6 = Nbody6ParameterSpec(;
        (
            k => get(renamed, k, getfield(nb, k)) for
            k in fieldnames(Nbody6ParameterSpec) if k != :kz
        )...,
        kz = kz_over,
    )

    st_raw = _section_table(m, "stellar", "merger.stellar")
    stellar = _parse_section(StellarSpec, st_raw, "merger.stellar")

    td_raw = _section_table(m, "tidal", "merger.tidal")
    tidal = _parse_section(TidalSpec, td_raw, "merger.tidal")
    (length(tidal.rg) == 3 && length(tidal.vg) == 3) ||
        _config_error("merger.tidal.rg", "and merger.tidal.vg must have length 3")
    # The frame key acts under a rotating-frame tide only; the default is inert elsewhere.
    (haskey(td_raw, "velocity_frame") && !(tidal.kz14 in (1, 2))) && _config_error(
        "merger.tidal.velocity_frame",
        "needs the rotating frame of kz14 = 1 or 2; got kz14 = $(tidal.kz14); remove it",
    )

    seed_raw = get(m, "seed", nothing)
    seed::Union{Int,Nothing} = if seed_raw === nothing
        nothing
    else
        _typed(seed_raw, Int, "merger.seed")   # any integer is a real seed, including 0
    end

    # density_ratio: rbar of cluster i from the mean density inside the
    # half-mass radius relative to cluster 1, ρ̄ ∝ M / r_h³. Cluster 1 and the
    # cluster itself are validated first (the expected masses need valid IMFs).
    prevalidated = falses(max(n_clusters, 0))
    for i in 1:n_clusters
        c = m["cluster$i"]
        haskey(c, "density_ratio") || continue
        key = "merger.cluster$i.density_ratio"
        i == 1 && _config_error(
            key,
            "is the density relative to cluster 1 and is allowed in [merger.cluster2] and later only",
        )
        ratio = _typed(c["density_ratio"], Float64, key)
        (isfinite(ratio) && ratio > 0) || _config_error(key, "must be finite and > 0; got $ratio")
        haskey(c, "rbar") && _config_error(
            key,
            "and merger.cluster$i.rbar are mutually exclusive (density_ratio derives rbar)",
        )
        for j in (1, i)
            prevalidated[j] && continue
            _validate_cluster_spec(clusters[j], j)
            prevalidated[j] = true
        end
        spec = clusters[i]
        M_ratio = expected_mass(spec.imf, spec.N) / expected_mass(clusters[1].imf, clusters[1].N)
        clusters[i] = ClusterSpec(;
            N = spec.N,
            rbar = clusters[1].rbar * cbrt(M_ratio / ratio),
            profile = spec.profile,
            imf = spec.imf,
            position = spec.position,
            velocity = spec.velocity,
            binaries = spec.binaries,
            segregation = spec.segregation,
            imbh = spec.imbh,
        )
    end

    # Fail-fast validation of the parsed values (bounds documented in the
    # config comments and the manual). kepler_velocity re-checks the orbit
    # parameters at run time for programmatically constructed configs.
    for (i, spec) in enumerate(clusters)
        prevalidated[i] || _validate_cluster_spec(spec, i)
    end
    for (i, spec) in enumerate(clusters)
        spec.N ≤ virial_max_n || _config_error(
            "merger.cluster$i.N",
            "= $(spec.N) exceeds merger.virial_max_n = $virial_max_n: " *
            "virialising a cluster evaluates its exact O(N²) potential " *
            "(≈ $(round(spec.N^2 / 2e9; digits = 1))×10⁹ pair evaluations); " *
            "raise merger.virial_max_n to accept that cost",
        )
    end
    complex_template === nothing || push!(clusters, complex_template)
    if orbit_mode == "kepler"
        size_key = isempty(size_keys) ? "apocentre" : size_keys[1]
        size_value = getfield(orbit, Symbol(size_key))::Float64
        (isfinite(size_value) && size_value > 0) ||
            _config_error("merger.orbit.$size_key", "must be > 0; got $size_value")
        (0 ≤ orbit.eccentricity < 1) || _config_error(
            "merger.orbit.eccentricity",
            "must satisfy 0 ≤ e < 1; got $(orbit.eccentricity)",
        )
        n = orbit.orbit_normal
        (length(n) == 3 && all(isfinite, n) && sqrt(sum(abs2, n)) > 0) || _config_error(
            "merger.orbit.orbit_normal",
            "must be three finite components of non-zero norm; got $n",
        )
        isfinite(orbit.apocentre_angle) || _config_error(
            "merger.orbit.apocentre_angle",
            "must be finite [degrees]; got $(orbit.apocentre_angle)",
        )
    end
    output.format == "nbody" || _config_error(
        "merger.output.format",
        "must be \"nbody\" (KZ(22) = 2, the only supported format); got \"$(output.format)\"",
    )
    output.truncate_rule in _TRUNCATE_RULES || _config_error(
        "merger.output.truncate_rule",
        "must be one of \"roche\", \"hill\"; got \"$(output.truncate_rule)\"",
    )
    output.tcrit > 0 || _config_error("merger.output.tcrit", "must be > 0; got $(output.tcrit)")
    output.dtadj > 0 || _config_error("merger.output.dtadj", "must be > 0; got $(output.dtadj)")
    output.deltat > 0 || _config_error("merger.output.deltat", "must be > 0; got $(output.deltat)")
    output.tcrit_myr ≥ 0 ||
        _config_error("merger.output.tcrit_myr", "must be ≥ 0; got $(output.tcrit_myr)")
    output.dtadj_myr ≥ 0 ||
        _config_error("merger.output.dtadj_myr", "must be ≥ 0; got $(output.dtadj_myr)")
    output.deltat_myr ≥ 0 ||
        _config_error("merger.output.deltat_myr", "must be ≥ 0; got $(output.deltat_myr)")
    (isfinite(output.truncate_radius_rh) && output.truncate_radius_rh ≥ 0) || _config_error(
        "merger.output.truncate_radius_rh",
        "must be finite and ≥ 0 (0 = Jacobi radius); got $(output.truncate_radius_rh)",
    )
    _validate_nbody6(nbody6)
    _validate_stellar(stellar, output)
    _validate_binary_physics(nbody6, stellar)
    _validate_tidal(tidal)
    _validate_tidal_tolerance(tidal, nbody6)
    _validate_imbh_configuration(clusters, nbody6, stellar)
    mbh_raw = _section_table(m, "mbh", "merger.mbh")
    mbh = isempty(mbh_raw) ? nothing : _parse_section(MBHSpec, mbh_raw, "merger.mbh")
    mbh === nothing || _validate_mbh(mbh, orbit_mode, stellar, tidal, nbody6)

    return MergerConfig(
        clusters,
        orbit_mode,
        orbit,
        output,
        nbody6,
        stellar,
        tidal,
        seed,
        virial_max_n,
        seed_streams,
        complex_spec,
        mbh,
    )
end

"""Fail-fast checks of `[merger.mbh]` against the other sections."""
function _validate_mbh(
    b::MBHSpec,
    orbit_mode::AbstractString,
    stellar::StellarSpec,
    tidal::TidalSpec,
    nbody6::Nbody6ParameterSpec,
)
    (isfinite(b.mass_msun) && b.mass_msun > 0) ||
        _config_error("merger.mbh.mass_msun", "must be finite and > 0 [M☉]; got $(b.mass_msun)")
    (isfinite(b.softening_pc) && b.softening_pc ≥ 0) || _config_error(
        "merger.mbh.softening_pc",
        "must be finite and ≥ 0 [pc]; got $(b.softening_pc)",
    )
    b.rt_mode in ("fixed", "stellar") || _config_error(
        "merger.mbh.rt_mode",
        "must be \"fixed\" or \"stellar\"; got \"$(b.rt_mode)\"",
    )
    (b.rt_mode == "fixed" && !(isfinite(b.r_acc_pc) && b.r_acc_pc > 0)) && _config_error(
        "merger.mbh.r_acc_pc",
        "must be finite and > 0 [pc] under rt_mode = \"fixed\"; got $(b.r_acc_pc)",
    )
    (b.rt_mode == "stellar" && stellar.kz19 < 3) && _config_error(
        "merger.mbh.rt_mode",
        "= \"stellar\" needs stellar radii (merger.stellar.kz19 ≥ 3); got kz19 = $(stellar.kz19)",
    )
    haskey(_MBH_IPART, b.partial) || _config_error(
        "merger.mbh.partial",
        "must be one of \"off\", \"gamma43\", \"gamma53\"; got \"$(b.partial)\"",
    )
    (isfinite(b.beta_full) && b.beta_full > 0) ||
        _config_error("merger.mbh.beta_full", "must be finite and > 0; got $(b.beta_full)")
    (isfinite(b.accretion_fraction) && 0 ≤ b.accretion_fraction ≤ 1) || _config_error(
        "merger.mbh.accretion_fraction",
        "must lie in [0, 1]; got $(b.accretion_fraction)",
    )
    b.response in ("jeans", "rigid") || _config_error(
        "merger.mbh.response",
        "must be \"jeans\" or \"rigid\"; got \"$(b.response)\"",
    )
    orbit_mode == "single" || _config_error(
        "merger.mbh",
        "needs orbit_mode = \"single\" (the black hole sits at the origin, at the cluster centre); " *
        "got \"$orbit_mode\"",
    )
    tidal.kz14 == 0 || _config_error(
        "merger.mbh",
        "needs an isolated cluster (merger.tidal.kz14 = 0); got kz14 = $(tidal.kz14)",
    )
    haskey(nbody6.kz, 44) &&
        _config_error("merger.nbody6.kz", "KZ(44) is set by [merger.mbh]; remove the override")
    return nothing
end

"""
    _section_table(parent, name, section) -> AbstractDict

The subtable `parent[name]` of a parsed TOML table, an empty table when it is
absent; a value that is not a table raises `_config_error(section, ...)`.
"""
function _section_table(parent::AbstractDict, name::AbstractString, section::AbstractString)
    haskey(parent, name) || return Dict{String,Any}()
    t = parent[name]
    t isa AbstractDict || _config_error(section, "must be a table; got $t ($(typeof(t)))")
    return t
end

"""Keys of the member template `[merger.complex.member]` that the complex sets itself and that are therefore refused."""
const _COMPLEX_MEMBER_REFUSED = ("N", "rbar", "position", "velocity", "density_ratio", "imbh")

"""
    _parse_complex_section(c::AbstractDict) -> (ComplexSpec, ClusterSpec)

`[merger.complex]` as a validated [`ComplexSpec`](@ref) and its member template
`[merger.complex.member]` as a `ClusterSpec` (parsed as a cluster table whose
errors name `merger.complex.member`) with the placeholders
`N = member_n_min`, `rbar = 1.0`. The template refuses `N`, `rbar`,
`position`, `velocity` and `density_ratio`, which the complex sets, `imbh`, and an
IMF with a target mass (`mass_total`, `target_mass`): the member masses come
from the mass function. A control (`control ≠ "none"`) refuses template
binaries that change the stellar masses (`pairing = "uniform_q"` with
`fraction > 0`, a massive branch with `pairing = "replace"` and
`fraction > 0`): the control needs the star list of the complex.
"""
function _parse_complex_section(c)::Tuple{ComplexSpec,ClusterSpec}
    c isa AbstractDict || _config_error("merger.complex", "must be a table")
    spec = _parse_section(ComplexSpec, c, "merger.complex"; skip = ("member",))
    _validate_complex(spec)
    t = get(c, "member", nothing)
    t isa AbstractDict || _config_error(
        "merger.complex.member",
        "table is missing (the member template of orbit_mode = \"complex\")",
    )
    for k in _COMPLEX_MEMBER_REFUSED
        haskey(t, k) && _config_error(
            "merger.complex.member.$k",
            k == "imbh" ?
            "is not taken by a cluster complex (the kz24 route seeds one IMBH, in a cluster table); remove it" :
            "is set by the complex (mass function, size–mass relation, centres); remove it",
        )
    end
    imf_raw = get(t, "imf", nothing)
    (haskey(t, "mass_total") || (imf_raw isa AbstractDict && haskey(imf_raw, "target_mass"))) &&
        _config_error(
            "merger.complex.member.mass_total",
            "is refused: the member masses come from merger.complex.mass_function",
        )
    parsed = _parse_cluster_table(t, "merger.complex.member", "complex")
    if spec.control != "none"
        b = parsed.binaries
        (b.fraction > 0 && b.pairing == "uniform_q") && _config_error(
            "merger.complex.control",
            "= \"$(spec.control)\" refuses merger.complex.member.binaries.pairing = \"uniform_q\": " *
            "the pairing changes the stellar masses, and the control needs the complex's star list",
        )
        (b.massive.model != "none" && b.massive.fraction > 0 && b.massive.pairing == "replace") &&
            _config_error(
                "merger.complex.control",
                "= \"$(spec.control)\" refuses merger.complex.member.binaries.massive.pairing = " *
                "\"replace\": the pairing changes the stellar masses, and the control needs the " *
                "complex's star list",
            )
    end
    parsed.imf isa RescaledKroupaIMF && _config_error(
        "merger.complex.member.imf",
        "must not be a rescaled Kroupa IMF: the member masses come from merger.complex.mass_function",
    )
    template = ClusterSpec(;
        N = spec.member_n_min,
        rbar = 1.0,
        profile = parsed.profile,
        imf = parsed.imf,
        binaries = parsed.binaries,
        segregation = parsed.segregation,
    )
    _validate_cluster_spec(template, "merger.complex.member")
    return spec, template
end

"""
    _validate_complex(s::ComplexSpec)

Fail-fast bounds of `[merger.complex]` (the table of [`ComplexSpec`](@ref)),
raising an `ArgumentError` that names `merger.complex.<key>`.
"""
function _validate_complex(s::ComplexSpec)
    key(k) = "merger.complex.$k"
    s.mass_function in _COMPLEX_MASS_FUNCTIONS || _config_error(
        key("mass_function"),
        "must be one of \"equal\", \"power_law\"; got \"$(s.mass_function)\"",
    )
    s.size_mass in _COMPLEX_SIZE_MASS || _config_error(
        key("size_mass"),
        "must be one of \"constant\", \"marks_kroupa\", \"frost\", \"dobbs\"; got \"$(s.size_mass)\"",
    )
    s.centres in _COMPLEX_CENTRES || _config_error(
        key("centres"),
        "must be one of \"uniform_sphere\", \"plummer\", \"fractal\", \"filament\"; got \"$(s.centres)\"",
    )
    (s.n_members == 0 || s.n_members ≥ 2) ||
        _config_error(key("n_members"), "must be ≥ 2 (or 0 = unused); got $(s.n_members)")
    (isfinite(s.mass_total) && s.mass_total ≥ 0) || _config_error(
        key("mass_total"),
        "must be finite and > 0 [M☉] (or 0 = unused); got $(s.mass_total)",
    )
    has_n, has_m = s.n_members > 0, s.mass_total > 0
    if s.mass_function == "equal"
        has_n || _config_error(key("n_members"), "must be > 0 for mass_function = \"equal\"")
        has_m || _config_error(key("mass_total"), "must be > 0 [M☉] for mass_function = \"equal\"")
    else
        xor(has_n, has_m) || _config_error(
            key("n_members"),
            "and merger.complex.mass_total: mass_function = \"power_law\" takes exactly one " *
            "of them; got n_members = $(s.n_members), mass_total = $(s.mass_total)",
        )
        isfinite(s.beta) || _config_error(key("beta"), "must be finite; got $(s.beta)")
        (isfinite(s.mass_hi) && 0 < s.mass_lo < s.mass_hi) || _config_error(
            key("mass_lo"),
            "and merger.complex.mass_hi must satisfy 0 < mass_lo < mass_hi [M☉]; " *
            "got $(s.mass_lo), $(s.mass_hi)",
        )
    end
    (isfinite(s.rh) && s.rh > 0) || _config_error(key("rh"), "must be > 0 [pc]; got $(s.rh)")
    (isfinite(s.size_scatter_dex) && s.size_scatter_dex ≥ 0) ||
        _config_error(key("size_scatter_dex"), "must be ≥ 0 [dex]; got $(s.size_scatter_dex)")
    (isfinite(s.radius) && s.radius > 0) ||
        _config_error(key("radius"), "must be > 0 [pc]; got $(s.radius)")
    (isfinite(s.plummer_cutoff) && s.plummer_cutoff > 1) || _config_error(
        key("plummer_cutoff"),
        "must be > 1 (units of the Plummer scale a); got $(s.plummer_cutoff)",
    )
    (1 ≤ s.fractal_dimension ≤ 3) || _config_error(
        key("fractal_dimension"),
        "must satisfy 1 ≤ D ≤ 3; got $(s.fractal_dimension)",
    )
    (isfinite(s.filament_length) && s.filament_length > 0) ||
        _config_error(key("filament_length"), "must be > 0 [pc]; got $(s.filament_length)")
    (isfinite(s.filament_width) && s.filament_width > 0) ||
        _config_error(key("filament_width"), "must be > 0 [pc]; got $(s.filament_width)")
    (isfinite(s.virial_ratio) && s.virial_ratio ≥ 0) ||
        _config_error(key("virial_ratio"), "must be ≥ 0; got $(s.virial_ratio)")
    (0 ≤ s.rotation_fraction < 1) || _config_error(
        key("rotation_fraction"),
        "must satisfy 0 ≤ f < 1; got $(s.rotation_fraction)",
    )
    (s.virial_ratio == 0 && s.rotation_fraction != 0) && _config_error(
        key("rotation_fraction"),
        "must be 0 with merger.complex.virial_ratio = 0 (a cold start); got $(s.rotation_fraction)",
    )
    (isfinite(s.softening) && s.softening ≥ 0) || _config_error(
        key("softening"),
        "must be ≥ 0 [pc] (0 = the mean member r_h); got $(s.softening)",
    )
    s.member_n_min ≥ 2 || _config_error(key("member_n_min"), "must be ≥ 2; got $(s.member_n_min)")
    s.imf_partition in _COMPLEX_IMF_PARTITIONS || _config_error(
        key("imf_partition"),
        "must be one of \"member\", \"global\"; got \"$(s.imf_partition)\"",
    )
    s.control in _COMPLEX_CONTROLS || _config_error(
        key("control"),
        "must be one of \"none\", \"monolithic\", \"cold\"; got \"$(s.control)\"",
    )
    (s.control != "none" && s.imf_partition != "global") && _config_error(
        key("control"),
        "= \"$(s.control)\" needs merger.complex.imf_partition = \"global\" (the control " *
        "takes the complex's star list); got \"$(s.imf_partition)\"",
    )
    return nothing
end

# Energy tolerance below which a tidal-field run halts at its first adjustment:
# this fork never evaluates the tidal potential energy unless it is built with
# the package's tidal-energy-total amendment (the XTRNLV call in energy.F is
# commented out), so an unamended engine's DE measures the tidal work.
const _QE_MIN_TIDAL = 0.01

"""
    _validate_tidal_tolerance(t::TidalSpec, p::Nbody6ParameterSpec)

Refuse an external field with an energy tolerance below `$(_QE_MIN_TIDAL)`:
an engine built without the `tidal-energy-total` source amendment omits the
tidal potential energy from its bookkeeping, so the relative energy change it
checks against `QE` is the work done by the field (10⁻⁴ to 10⁻² per
adjustment interval for realistic fields) and the run would halt at its
first adjustment. The generator does not know which build will run the
input, so the floor applies to every tidal configuration.
"""
function _validate_tidal_tolerance(t::TidalSpec, p::Nbody6ParameterSpec)
    if t.kz14 > 0 && p.qe < _QE_MIN_TIDAL
        _config_error(
            "merger.nbody6.qe",
            "= $(p.qe) with merger.tidal.kz14 = $(t.kz14): an engine built without the " *
            "tidal-energy-total amendment does not evaluate the tidal potential energy (energy.F), " *
            "so its energy check measures the tidal work and halts the run at the first " *
            "adjustment; set qe ≥ $(_QE_MIN_TIDAL) for tidal-field runs",
        )
    end
    return nothing
end

"""Fail-fast checks of `[merger.tidal]`: supported `kz14` and the parameters each option needs."""
function _validate_tidal(t::TidalSpec)
    if t.kz14 in (3, 4)
        _config_error(
            "merger.tidal.kz14",
            "= $(t.kz14) is not supported: on these options the engine " *
            "rescales all velocities to the &INSCALE virial ratio including the external potential, " *
            "which destroys the prescribed orbital kinematics; use 1, 2, or 5",
        )
    end
    t.kz14 in (0, 1, 2, 5) ||
        _config_error("merger.tidal.kz14", "must be one of 0, 1, 2, 5; got $(t.kz14)")
    if t.kz14 == 2
        t.gmg > 0 ||
            _config_error("merger.tidal.gmg", "must be > 0 [M☉] for kz14 = 2; got $(t.gmg)")
        t.rg0 > 0 ||
            _config_error("merger.tidal.rg0", "must be > 0 [kpc] for kz14 = 2; got $(t.rg0)")
    end
    t.velocity_frame in ("rotating", "inertial") || _config_error(
        "merger.tidal.velocity_frame",
        "must be \"rotating\" or \"inertial\"; got \"$(t.velocity_frame)\"",
    )
    if t.kz14 == 5
        any(!=(0.0), t.rg) ||
            _config_error("merger.tidal.rg", "must be non-zero [kpc] for kz14 = 5")
        any(!=(0.0), t.vg) ||
            _config_error("merger.tidal.vg", "must be non-zero [km/s] for kz14 = 5")
    end
    return nothing
end

"""Refuse `binary_physics = "full"` without stellar evolution (`kz19 ≥ 3`): the gravitational-wave decay, tides and Roche transfer of the set act on the stellar types; refuse `post_newtonian` without the Peters path `KZ(27) = 3`, which decays unperturbed pairs and holds the merger test."""
function _validate_binary_physics(p::Nbody6ParameterSpec, s::StellarSpec)
    (p.gw_recoil != "off" && s.kz19 < 3) && _config_error(
        "merger.nbody6.gw_recoil",
        "= \"$(p.gw_recoil)\" requires stellar evolution (merger.stellar.kz19 ≥ 3): the recoil " *
        "is applied by the collision and coalescence routines of the stellar-evolution path; " *
        "got kz19 = $(s.kz19)",
    )
    kz27 = get(p.kz, 27, p.binary_physics == "full" ? 3 : 0)
    (p.post_newtonian != "off" && kz27 < 3) && _config_error(
        "merger.nbody6.post_newtonian",
        "= \"$(p.post_newtonian)\" requires the Peters path KZ(27) = 3 (binary_physics = \"full\" " *
        "or merger.nbody6.kz 27 = 3), which decays unperturbed pairs and holds the merger test; " *
        "got KZ(27) = $kz27",
    )
    (p.binary_physics == "full" && s.kz19 < 3) && _config_error(
        "merger.nbody6.binary_physics",
        "= \"full\" requires stellar evolution (merger.stellar.kz19 ≥ 3): the " *
        "gravitational-wave decay, tides and Roche transfer act on the stellar types; " *
        "got kz19 = $(s.kz19)",
    )
    return nothing
end

"""Largest IMBH count of the flag route: the engine's list holds 16 massive black holes (`imbhinit.F`)."""
const _IMBH_FLAG_MAX = 16

"""
    _resolve_imbh_route(route, n_imbh) -> Symbol

Route of `n_imbh` seeded IMBHs under `imbh_route = route`: `:none` without an
IMBH, `:kz24` or `:flag` as configured, and for `"auto"` `:kz24` for one IMBH
and `:flag` for several. The consistency of the pair is checked by
[`_imbh_route_violation`](@ref).
"""
function _resolve_imbh_route(route::AbstractString, n_imbh::Integer)::Symbol
    n_imbh == 0 && return :none
    route == "kz24" && return :kz24
    route == "flag" && return :flag
    return n_imbh == 1 ? :kz24 : :flag
end

"""
    _imbh_route_violation(route, carriers) -> Union{Nothing,Tuple{String,String}}

`(key, message)` of the IMBH-count rule that `imbh_route = route` and the
indices `carriers` of the clusters with an IMBH violate, `nothing` when they
are consistent: `"kz24"` takes at most one IMBH, `"flag"` at least one and at
most `_IMBH_FLAG_MAX`, `"auto"` at most `_IMBH_FLAG_MAX`.
"""
function _imbh_route_violation(route::AbstractString, carriers::AbstractVector{<:Integer})
    n = length(carriers)
    if route == "kz24" && n > 1
        return (
            "merger.cluster$(carriers[2]).imbh",
            "the kz24 route takes one IMBH; $n clusters carry one",
        )
    elseif route == "flag" && n == 0
        return (
            "merger.nbody6.imbh_route",
            "\"flag\" needs at least one [merger.clusterK.imbh] table",
        )
    elseif n > _IMBH_FLAG_MAX
        return (
            "merger.cluster$(carriers[_IMBH_FLAG_MAX + 1]).imbh",
            "the engine's list holds $(_IMBH_FLAG_MAX) massive black holes; $n clusters carry one",
        )
    end
    return nothing
end

"""
    resolved_imbh_route(cfg::MergerConfig) -> Symbol

Route by which the IMBHs of `cfg` reach the engine: `:none` (no cluster
carries an [`IMBHSpec`](@ref)), `:kz24` (`KZ(24) = 1`, one IMBH repeated on
the `IMBHINIT` line) or `:flag` (`KZ(24) = 3`, engine fork ≥ v2026.07.11, the
IMBH names after `NIMBHS DTBH`). `imbh_route = "auto"` resolves to `:kz24`
for one IMBH and to `:flag` for several. The configuration is assumed
validated ([`load_merger_config`](@ref)).

# Example
```julia
cfg = load_merger_config("merger.toml")
resolved_imbh_route(cfg) === :flag && @info "several IMBHs, KZ(24) = 3"
```
"""
resolved_imbh_route(cfg::MergerConfig)::Symbol =
    _resolve_imbh_route(cfg.nbody6.imbh_route, count(c -> c.imbh !== nothing, cfg.clusters))

"""
    _validate_imbh_configuration(clusters, nbody6, stellar)

Configuration-level checks of IMBH seeding in [`load_merger_config`](@ref):
the IMBH count against `imbh_route` ([`_imbh_route_violation`](@ref): at most
one under `"kz24"`, 1 to 16 under `"flag"`, at most 16 under `"auto"`), and
with an IMBH in some cluster: a `kz` override of index 24 (set by the route)
and, under stellar evolution (`kz19 ≥ 3`), a `mass` below 100 M☉ of any IMBH
are refused; `escapers = "remove"` warns once (an ejected IMBH is removed and
every later restart aborts).
"""
function _validate_imbh_configuration(
    clusters::AbstractVector{ClusterSpec},
    nbody6::Nbody6ParameterSpec,
    stellar::StellarSpec,
)
    carriers = findall(c -> c.imbh !== nothing, clusters)
    violation = _imbh_route_violation(nbody6.imbh_route, carriers)
    violation === nothing || _config_error(violation...)
    isempty(carriers) && return nothing
    haskey(nbody6.kz, 24) &&
        _config_error("merger.nbody6.kz", "KZ(24) is set by the IMBH route; remove the override")
    for k in carriers
        s = clusters[k].imbh::IMBHSpec
        (stellar.kz19 ≥ 3 && isfinite(s.mass) && s.mass < 100) && _config_error(
            "merger.cluster$(k).imbh.mass",
            "an IMBH under stellar evolution must be ≥ 100 M☉ (below, indistinguishable from the " *
            "SSE black holes); got $(s.mass) M☉",
        )
    end
    k = carriers[1]
    nbody6.escapers == "remove" &&
        @warn "merger.cluster$(k).imbh: with escapers = \"remove\" an IMBH ejected beyond twice " *
              "the tidal radius is removed and every later restart aborts (imbhinit.F); use " *
              "escapers = \"keep\" for segmented runs"
    return nothing
end

"""Fail-fast bounds of `[merger.stellar]` (the engine's own limits on `zmet`; `level` one of A, B, C, 0)."""
function _validate_stellar(s::StellarSpec, output::MergerOutputSpec)
    s.kz19 ≥ 0 || _config_error("merger.stellar.kz19", "must be ≥ 0; got $(s.kz19)")
    s.level in ("A", "B", "C", "0") || _config_error(
        "merger.stellar.level",
        "must be one of \"A\", \"B\", \"C\", \"0\"; got \"$(s.level)\"",
    )
    (0.0001 ≤ s.zmet ≤ 0.03) ||
        _config_error("merger.stellar.zmet", "must satisfy 0.0001 ≤ zmet ≤ 0.03; got $(s.zmet)")
    s.epoch0 ≤ 0 || _config_error("merger.stellar.epoch0", "must be ≤ 0 [Myr]; got $(s.epoch0)")
    s.dtplot > 0 || _config_error("merger.stellar.dtplot", "must be > 0 [NB]; got $(s.dtplot)")
    s.dtplot_myr ≥ 0 ||
        _config_error("merger.stellar.dtplot_myr", "must be ≥ 0 [Myr]; got $(s.dtplot_myr)")
    # With NB intervals on both sides the ordering is checked here; with a
    # physical interval on either side it is checked at generation, once the
    # time unit is known.
    if s.dtplot_myr == 0 && output.deltat_myr == 0
        s.dtplot ≥ output.deltat || _config_error(
            "merger.stellar.dtplot",
            "must be ≥ merger.output.deltat; got dtplot = $(s.dtplot), deltat = $(output.deltat)",
        )
    elseif s.dtplot_myr > 0 && output.deltat_myr > 0
        s.dtplot_myr ≥ output.deltat_myr || _config_error(
            "merger.stellar.dtplot_myr",
            "must be ≥ merger.output.deltat_myr; got $(s.dtplot_myr), $(output.deltat_myr)",
        )
    end
    return nothing
end

"""Fail-fast bounds of `[merger.nbody6]`: positive tolerances, non-negative step factors and derivable entries, a known `preset` and its restrictions, `kz16` in 0–3."""
function _validate_nbody6(p::Nbody6ParameterSpec)
    p.qe > 0 || _config_error("merger.nbody6.qe", "must be > 0; got $(p.qe)")
    p.etai ≥ 0 ||
        _config_error("merger.nbody6.etai", "must be ≥ 0 (0 = preset value); got $(p.etai)")
    p.etar ≥ 0 ||
        _config_error("merger.nbody6.etar", "must be ≥ 0 (0 = preset value); got $(p.etar)")
    p.etau ≥ 0 ||
        _config_error("merger.nbody6.etau", "must be ≥ 0 (0 = preset value); got $(p.etau)")
    haskey(_GW_RECOIL_KICKTYPE, p.gw_recoil) || _config_error(
        "merger.nbody6.gw_recoil",
        "must be one of \"off\", \"fit_q\", \"nr_fits\", \"test_1e4\"; got \"$(p.gw_recoil)\"",
    )
    haskey(_BH_SPIN_FLAG, p.bh_spin) || _config_error(
        "merger.nbody6.bh_spin",
        "must be one of \"zero\", \"uniform\", \"gaussian\", \"maxwellian\", \"geneva\", \"mesa\"; " *
        "got \"$(p.bh_spin)\"",
    )
    haskey(_POST_NEWTONIAN_KZ43, p.post_newtonian) || _config_error(
        "merger.nbody6.post_newtonian",
        "must be \"off\" or \"full\"; got \"$(p.post_newtonian)\"",
    )
    haskey(_ESCAPERS_KZ23, p.escapers) || _config_error(
        "merger.nbody6.escapers",
        "must be \"remove\" or \"keep\"; got \"$(p.escapers)\"",
    )
    p.imbh_route in _IMBH_ROUTES || _config_error(
        "merger.nbody6.imbh_route",
        "must be one of \"auto\", \"kz24\", \"flag\"; got \"$(p.imbh_route)\"",
    )
    (isfinite(p.imbh_output_dt_myr) && p.imbh_output_dt_myr > 0) || _config_error(
        "merger.nbody6.imbh_output_dt_myr",
        "must be finite and > 0 [Myr]; got $(p.imbh_output_dt_myr)",
    )
    p.bh_spin in ("uniform", "gaussian", "maxwellian") &&
        @warn "merger.nbody6.bh_spin = \"$(p.bh_spin)\": the engine draws these spins with an unseeded " *
              "generator, so the run is not restart-reproducible"
    p.preset in _NBODY6_PRESETS || _config_error(
        "merger.nbody6.preset",
        "must be one of \"engine\", \"baseline\"; got \"$(p.preset)\"",
    )
    p.nnbopt ≥ 0 ||
        _config_error("merger.nbody6.nnbopt", "must be ≥ 0 (0 = derived); got $(p.nnbopt)")
    p.rs0 ≥ 0 || _config_error("merger.nbody6.rs0", "must be ≥ 0 (0 = derived); got $(p.rs0)")
    p.rmin ≥ 0 || _config_error("merger.nbody6.rmin", "must be ≥ 0 (0 = derived); got $(p.rmin)")
    p.dtmin ≥ 0 || _config_error("merger.nbody6.dtmin", "must be ≥ 0 (0 = derived); got $(p.dtmin)")
    p.kz16 in 0:3 || _config_error("merger.nbody6.kz16", "must be one of 0, 1, 2, 3; got $(p.kz16)")
    p.eclose > 0 || _config_error("merger.nbody6.eclose", "must be > 0; got $(p.eclose)")
    (0 < p.gmin < p.gmax) || _config_error(
        "merger.nbody6",
        "must satisfy 0 < gmin < gmax; got gmin = $(p.gmin), gmax = $(p.gmax)",
    )
    p.smax > 0 || _config_error("merger.nbody6.smax", "must be > 0; got $(p.smax)")
    p.tcomp > 0 || _config_error("merger.nbody6.tcomp", "must be > 0 [s]; got $(p.tcomp)")
    p.tcrtp0 > 0 || _config_error("merger.nbody6.tcrtp0", "must be > 0 [Myr]; got $(p.tcrtp0)")
    for k in (:isernb, :iserreg, :iserks)
        getfield(p, k) ≥ 0 ||
            _config_error("merger.nbody6.$k", "must be ≥ 0; got $(getfield(p, k))")
    end
    for k in (:nfix, :ncrit, :nrun, :ncomm)
        getfield(p, k) ≥ 1 ||
            _config_error("merger.nbody6.$k", "must be ≥ 1; got $(getfield(p, k))")
    end
    p.binary_physics in _BINARY_PHYSICS || _config_error(
        "merger.nbody6.binary_physics",
        "must be one of \"inherited\", \"full\"; got \"$(p.binary_physics)\"",
    )
    for (i, v) in p.kz
        1 ≤ i ≤ 50 || _config_error("merger.nbody6.kz", "index must be within 1–50; got $i")
    end
    get(p.kz, 18, 0) ≥ 2 && _config_error(
        "merger.nbody6.kz",
        "must not set KZ(18) ≥ 2: it runs the engine's primordial-triple generator on the dat.10 " *
        "pairs, and the generator writes no triples (NHI0 = 0); got KZ(18) = $(p.kz[18])",
    )
    if p.preset == "baseline"
        p.kz16 == 0 || _config_error(
            "merger.nbody6.kz16",
            "must be 0 with preset = \"baseline\" (RMIN, DTMIN and ECLOSE stay as written); " *
            "got $(p.kz16)",
        )
        for i in (16, 17, 32)
            (haskey(p.kz, i) && p.kz[i] != 0) && _config_error(
                "merger.nbody6.kz",
                "must not set KZ($i) with preset = \"baseline\" (no automatic adjustment " *
                "of the integration parameters); got KZ($i) = $(p.kz[i])",
            )
        end
    end
    if p.checkpoint
        p.ncomm == 1 || _config_error(
            "merger.nbody6.ncomm",
            "must be 1 with checkpoint = true: the engine counts every dump call and writes " *
            "one in ncomm, the dump of a stop request included; got $(p.ncomm)",
        )
        haskey(p.kz, 2) && _config_error(
            "merger.nbody6.kz",
            "must not set index 2 with checkpoint = true, which sets KZ(2) = 1; got KZ(2) = $(p.kz[2])",
        )
    end
    return nothing
end

"""Keys a `[merger.cluster<k>]` table may carry (flat and structured forms together)."""
const _CLUSTER_TABLE_KEYS = (
    "N",
    "rbar",
    "position",
    "velocity",
    "model",
    "W0",
    "imf",
    "bodyn",
    "body1",
    "mass_total",
    "profile",
    "binaries",
    "segregation",
    "density_ratio",
    "imbh",
)

"""Flat keys that a structured `imf` table (first entry) or `profile` table (second) excludes."""
const _FLAT_IMF_KEYS = ("bodyn", "body1", "mass_total")
const _FLAT_PROFILE_KEYS = ("model", "W0")

"""
    _parse_cluster_table(c, key, orbit_mode) -> ClusterSpec

One cluster table, `[merger.cluster<i>]` or the member template
`[merger.complex.member]` (`orbit_mode = "complex"`), in the flat or the
structured form; `key` is the dotted table name every error carries (an
`Int` index `i` stands for `"merger.cluster<i>"`). Refused: unknown keys (for
the member template, also the keys the complex sets); `profile` or
`binaries` given as a non-table; the flat IMF keys `bodyn`, `body1`,
`mass_total` beside a structured `imf` table and the flat profile keys
`model`, `W0` beside a structured `profile` table; `position` and `velocity`
in kepler and single mode. `orbit_mode = "metadata"` (the cluster tables of
`merger_ic.toml`) applies no mode rule to `position` and `velocity`.
"""
function _parse_cluster_table(
    c::AbstractDict,
    key::AbstractString,
    orbit_mode::AbstractString,
)::ClusterSpec
    # The member template lists only the keys it takes
    allowed = String[
        k for k in _CLUSTER_TABLE_KEYS if
        orbit_mode != "complex" || !(k in _COMPLEX_MEMBER_REFUSED || k == "mass_total")
    ]
    _reject_unknown(c, allowed, key)
    N = _typed(get(c, "N", 50000), Int, "$key.N")
    rbar = _typed(get(c, "rbar", 2.0), Float64, "$key.rbar")

    position = _typed(get(c, "position", Float64[]), Vector{Float64}, "$key.position")
    velocity = _typed(get(c, "velocity", Float64[]), Vector{Float64}, "$key.velocity")
    if orbit_mode == "single" && (haskey(c, "position") || haskey(c, "velocity"))
        _config_error(
            "$key.position",
            "and $key.velocity are not used in single mode (the cluster is placed at the origin, at rest); remove them",
        )
    end
    if orbit_mode == "kepler" && (haskey(c, "position") || haskey(c, "velocity"))
        _config_error(
            "$key.position",
            "and $key.velocity are not used in kepler mode (the orbit places the clusters); remove them",
        )
    end

    if orbit_mode == "explicit"
        length(position) == 3 ||
            _config_error("$key.position", "must be [x, y, z] in explicit orbit mode")
        length(velocity) == 3 ||
            _config_error("$key.velocity", "must be [vx, vy, vz] in explicit orbit mode")
    end

    for (name, flat) in (("imf", _FLAT_IMF_KEYS), ("profile", _FLAT_PROFILE_KEYS))
        (haskey(c, name) && c[name] isa AbstractDict) || continue
        for k in flat
            haskey(c, k) && _config_error(
                "$key.$k",
                "is a key of the flat form and is not mixed with the structured $name table; " *
                "move it into $key.$name or drop the table",
            )
        end
    end
    for name in ("profile", "binaries")
        (haskey(c, name) && !(c[name] isa AbstractDict)) && _config_error(
            "$key.$name",
            "must be a table; got $(c[name]) ($(typeof(c[name])))" *
            (name == "profile" ? "; the flat form names the profile with `model`" : ""),
        )
    end

    profile = if haskey(c, "profile")
        _parse_profile_table(c["profile"], key)
    else
        model = _typed(get(c, "model", "king"), String, "$key.model")
        if model == "king"
            KingProfile(W0 = _typed(get(c, "W0", 6.0), Float64, "$key.W0"))
        elseif model == "plummer"
            haskey(c, "W0") && @warn "[$key] W0 set but model=\"plummer\"; ignored."
            PlummerProfile()
        else
            _config_error("$key.model", "unknown model '$model'; supported: 'king', 'plummer'")
        end
    end

    imf = if haskey(c, "imf") && c["imf"] isa AbstractDict
        _parse_imf_table(c["imf"], N, key)
    else
        imf_str = _typed(get(c, "imf", "kroupa"), String, "$key.imf")
        bodyn = _typed(get(c, "bodyn", 0.08), Float64, "$key.bodyn")
        body1 = _typed(get(c, "body1", 100.0), Float64, "$key.body1")
        mt_raw = get(c, "mass_total", nothing)
        mass_total::Union{Float64,Nothing} =
            mt_raw === nothing ? nothing : _typed(mt_raw, Float64, "$key.mass_total")

        if imf_str == "kroupa"
            mass_total === nothing ? KroupaIMF(bodyn = bodyn, body1 = body1) :
            RescaledKroupaIMF(bodyn = bodyn, body1 = body1, target_mass = mass_total)
        elseif imf_str == "kroupa_rescaled"
            mass_total === nothing &&
                _config_error("$key.mass_total", "is required for imf = \"kroupa_rescaled\"")
            RescaledKroupaIMF(bodyn = bodyn, body1 = body1, target_mass = mass_total)
        elseif imf_str == "equal"
            mass_total === nothing &&
                _config_error("$key.mass_total", "is required for imf = \"equal\"")
            EqualMassIMF(particle_mass = mass_total / N)
        else
            _config_error(
                "$key.imf",
                "unknown imf '$imf_str'; supported: 'kroupa', 'kroupa_rescaled', 'equal'",
            )
        end
    end

    binaries = if haskey(c, "binaries")
        _parse_binaries_table(c["binaries"], key)
    else
        BinarySpec()
    end

    segregation = if haskey(c, "segregation")
        seg_raw = c["segregation"]
        seg_raw isa AbstractDict || _config_error("$key.segregation", "must be a table")
        _parse_segregation_table(seg_raw, key)
    else
        SegregationSpec()
    end

    imbh = if haskey(c, "imbh")
        imbh_raw = c["imbh"]
        imbh_raw isa AbstractDict ||
            _config_error("$key.imbh", "must be a table; got $imbh_raw ($(typeof(imbh_raw)))")
        _parse_imbh_table(imbh_raw, key)
    else
        nothing
    end

    return ClusterSpec(N, rbar, profile, imf, position, velocity, binaries, segregation, imbh)
end

_parse_cluster_table(c::AbstractDict, idx::Int, orbit_mode::AbstractString)::ClusterSpec =
    _parse_cluster_table(c, "merger.cluster$idx", orbit_mode)

"""
    _parse_imbh_table(t::AbstractDict, key::AbstractString) -> IMBHSpec

`[<key>.imbh]` as a validated [`IMBHSpec`](@ref) ([`_validate_imbh`](@ref));
`key` is the dotted cluster table. Keys other than the eight fields are refused.
"""
function _parse_imbh_table(t::AbstractDict, key::AbstractString)::IMBHSpec
    s = _parse_section(IMBHSpec, t, "$key.imbh")
    _validate_imbh(s, key)
    return s
end

"""
    _validate_imbh(s::IMBHSpec, key)

Fail-fast bounds of an [`IMBHSpec`](@ref), naming `<key>.imbh.<field>`:
exactly one of `mass > 0` [M☉] and `0 < mass_fraction ≤ 0.05` (warned above
0.02), the enumerations of `placement`, `velocity` and `response`, finite
three-component `offset` and `velocity_kms`, and `0 ≤ spin ≤ 1`.
"""
function _validate_imbh(s::IMBHSpec, key::AbstractString)
    k = "$key.imbh"
    has_mass, has_fraction = isfinite(s.mass), isfinite(s.mass_fraction)
    xor(has_mass, has_fraction) ||
        _config_error("$k.mass", "exactly one of mass and mass_fraction must be given")
    has_mass && (s.mass > 0 || _config_error("$k.mass", "must be > 0 [M☉]; got $(s.mass)"))
    if has_fraction
        (0 < s.mass_fraction ≤ 0.05) || _config_error(
            "$k.mass_fraction",
            "must satisfy 0 < mass_fraction ≤ 0.05 (BME04a: the cusp extends to M(<r) = M• " *
            "above 5 %); got $(s.mass_fraction)",
        )
        s.mass_fraction > 0.02 &&
            @warn "$k.mass_fraction = $(s.mass_fraction) exceeds 0.02: the IMBH dominates the " *
                  "central potential and the cluster leaves the equilibrium it was sampled in"
    end
    s.placement in _IMBH_PLACEMENTS || _config_error(
        "$k.placement",
        "must be one of \"centre\", \"offset\"; got \"$(s.placement)\"",
    )
    s.velocity in _IMBH_VELOCITIES || _config_error(
        "$k.velocity",
        "must be one of \"rest\", \"equipartition\", \"explicit\"; got \"$(s.velocity)\"",
    )
    s.response in _IMBH_RESPONSES ||
        _config_error("$k.response", "must be one of \"rigid\", \"jeans\"; got \"$(s.response)\"")
    (length(s.offset) == 3 && all(isfinite, s.offset)) ||
        _config_error("$k.offset", "must be three finite components [pc]; got $(s.offset)")
    (length(s.velocity_kms) == 3 && all(isfinite, s.velocity_kms)) || _config_error(
        "$k.velocity_kms",
        "must be three finite components [km s⁻¹]; got $(s.velocity_kms)",
    )
    (0 ≤ s.spin ≤ 1) || _config_error("$k.spin", "must satisfy 0 ≤ spin ≤ 1; got $(s.spin)")
    return nothing
end

# ---------------------------------------------------------------------------
# Fail-fast validation of parsed cluster specifications
# ---------------------------------------------------------------------------

"""
    _validate_cluster_spec(spec::ClusterSpec, key::AbstractString)
    _validate_cluster_spec(spec::ClusterSpec, idx::Int)

Validate one parsed `ClusterSpec` against the documented bounds; raises an
`ArgumentError` naming the offending `<key>.<field>`, `key` the dotted
cluster table (`idx` stands for `"merger.cluster<idx>"`).
"""
_validate_cluster_spec(spec::ClusterSpec, idx::Int) =
    _validate_cluster_spec(spec, "merger.cluster$idx")

function _validate_cluster_spec(spec::ClusterSpec, key_table::AbstractString)
    key = String(key_table)
    spec.N ≥ 2 || _config_error("$key.N", "must be ≥ 2; got $(spec.N)")
    spec.rbar > 0 || _config_error("$key.rbar", "must be > 0; got $(spec.rbar)")
    _validate_profile(spec.profile, key)
    _validate_imf(spec.imf, spec.N, key)
    _validate_binaries(spec.binaries, key)
    _validate_binary_extensions(spec.binaries, key)
    _validate_segregation(spec.segregation, key)
    return nothing
end

"Validate a density-profile tag; `key` names the owning cluster table."
function _validate_profile(p::KingProfile, key::String)
    p.W0 > 0 || _config_error("$key.W0", "must be > 0; got $(p.W0)")
    return nothing
end
_validate_profile(::PlummerProfile, ::String) = nothing
function _validate_profile(p::RotatingKingProfile, key::String)
    p.W0 > 0 || _config_error("$key.W0", "must be > 0; got $(p.W0)")
    (isfinite(p.omega0) && p.omega0 ≥ 0) ||
        _config_error("$key.omega0", "must be finite and ≥ 0; got $(p.omega0)")
    length(p.axis) == 3 ||
        _config_error("$key.axis", "must have three components; got $(length(p.axis))")
    all(isfinite, p.axis) || _config_error("$key.axis", "must be finite; got $(p.axis)")
    sqrt(sum(abs2, p.axis)) > 0 ||
        _config_error("$key.axis", "must have a non-zero norm; got $(p.axis)")
    return nothing
end

"Validate an IMF tag; `key` names the owning cluster table."
function _validate_imf(i::PowerLawIMF, ::Int, key::String)
    (0 < i.bodyn < i.body1) || _config_error(
        "$key.imf",
        "must satisfy 0 < bodyn < body1; got bodyn = $(i.bodyn), body1 = $(i.body1)",
    )
    isfinite(i.alpha) || _config_error("$key.imf.alpha", "must be finite; got $(i.alpha)")
    return nothing
end

function _validate_imf(i::KroupaIMF, ::Int, key::String)
    (0 < i.bodyn < i.body1) || _config_error(
        "$key.imf",
        "must satisfy 0 < bodyn < body1; got bodyn = $(i.bodyn), body1 = $(i.body1)",
    )
    return nothing
end
# Rescale factor c = target_mass / (N ⟨m⟩) of the rescaled Kroupa mode: bodies
# stop corresponding to stars beyond a factor of two either way (refused);
# beyond 1.4 the stellar-evolution output is already unreliable (warned).
const _IMF_RESCALE_WARN = 1.4
const _IMF_RESCALE_MAX = 2.0

function _validate_imf(i::RescaledKroupaIMF, N::Int, key::String)
    (0 < i.bodyn < i.body1) || _config_error(
        "$key.imf",
        "must satisfy 0 < bodyn < body1; got bodyn = $(i.bodyn), body1 = $(i.body1)",
    )
    i.target_mass > 0 || _config_error(
        "$key.imf",
        "target_mass (flat form: mass_total) must be > 0; got $(i.target_mass)",
    )
    mean_m = kroupa_mean_mass(i.bodyn, i.body1)
    c = i.target_mass / (N * mean_m)
    if !(1 / _IMF_RESCALE_MAX ≤ c ≤ _IMF_RESCALE_MAX)
        _config_error(
            key,
            "requests a Kroupa rescale factor of ×$(round(c; digits = 2)) " *
            "(target_mass = $(i.target_mass) M☉ against an expected natural mass of " *
            "$(round(N * mean_m; digits = 1)) M☉ for N = $N); factors outside " *
            "×[$(1 / _IMF_RESCALE_MAX), $(_IMF_RESCALE_MAX)] produce bodies that are not stars. " *
            "Drop mass_total for a natural Kroupa population, choose N and mass_total " *
            "consistent with the IMF mean $(round(mean_m; digits = 2)) M☉, or use " *
            "imf = \"equal\" for a collisionless super-particle model.",
        )
    elseif !(1 / _IMF_RESCALE_WARN ≤ c ≤ _IMF_RESCALE_WARN)
        @warn "config: $key rescales the Kroupa IMF by ×$(round(c; digits = 2)); " *
              "stellar-evolution output will not correspond to real stars"
    end
    return nothing
end
function _validate_imf(i::EqualMassIMF, ::Int, key::String)
    i.particle_mass > 0 || _config_error(
        "$key.imf",
        "particle_mass (flat form: mass_total / N) must be > 0; got $(i.particle_mass)",
    )
    return nothing
end

function _parse_profile_table(p::AbstractDict, table::AbstractString)::DensityProfile
    key = "$table.profile"
    t = _typed(get(p, "type", "king"), String, "$key.type")
    allowed = t == "rotating_king" ? ("type", "W0", "omega0", "axis") : ("type", "W0")
    _reject_unknown(p, allowed, key)
    if t == "king"
        KingProfile(W0 = _typed(get(p, "W0", 6.0), Float64, "$key.W0"))
    elseif t == "plummer"
        haskey(p, "W0") && _config_error("$key.W0", "is not a parameter of type = 'plummer'")
        PlummerProfile()
    elseif t == "rotating_king"
        RotatingKingProfile(
            W0 = _typed(get(p, "W0", 6.0), Float64, "$key.W0"),
            omega0 = _typed(get(p, "omega0", 0.6), Float64, "$key.omega0"),
            axis = _typed(get(p, "axis", [0.0, 0.0, 1.0]), Vector{Float64}, "$key.axis"),
        )
    else
        _config_error(
            "$key.type",
            "unknown type '$t'; supported: 'king', 'plummer', 'rotating_king'",
        )
    end
end

"""Keys a structured `imf` table of each `type` takes; any other key is refused, not dropped."""
const _IMF_TABLE_KEYS = Dict(
    "kroupa" => ("type", "bounds", "bodyn", "body1"),
    "kroupa_rescaled" => ("type", "bounds", "bodyn", "body1", "target_mass"),
    "equal" => ("type", "particle_mass", "target_mass"),
    "powerlaw" => ("type", "alpha", "bounds", "bodyn", "body1"),
)

function _parse_imf_table(imf::AbstractDict, N::Int, table::AbstractString)::IMFSpec
    key = "$table.imf"
    t = _typed(get(imf, "type", "kroupa"), String, "$key.type")
    haskey(_IMF_TABLE_KEYS, t) || _config_error(
        "$key.type",
        "unknown type '$t'; supported: 'kroupa', 'kroupa_rescaled', 'equal', 'powerlaw'",
    )
    _reject_unknown(imf, _IMF_TABLE_KEYS[t], key)
    (haskey(imf, "bounds") && (haskey(imf, "bodyn") || haskey(imf, "body1"))) && _config_error(
        "$key.bounds",
        "and $key.bodyn / $key.body1 give the same bounds; use one form",
    )
    (haskey(imf, "particle_mass") && haskey(imf, "target_mass")) && _config_error(
        "$key.particle_mass",
        "and $key.target_mass are mutually exclusive for type = 'equal'",
    )

    bounds_raw = get(imf, "bounds", nothing)
    bodyn, body1 = if bounds_raw !== nothing
        bounds = _typed(bounds_raw, Vector{Float64}, "$key.bounds")
        length(bounds) == 2 || _config_error("$key.bounds", "must have length 2")
        (bounds[1], bounds[2])
    else
        (
            _typed(get(imf, "bodyn", 0.08), Float64, "$key.bodyn"),
            _typed(get(imf, "body1", 100.0), Float64, "$key.body1"),
        )
    end

    if t == "kroupa"
        return KroupaIMF(bodyn = bodyn, body1 = body1)
    elseif t == "kroupa_rescaled"
        tm = get(imf, "target_mass", nothing)
        tm === nothing &&
            _config_error("$key.target_mass", "is required for type = 'kroupa_rescaled'")
        return RescaledKroupaIMF(
            bodyn = bodyn,
            body1 = body1,
            target_mass = _typed(tm, Float64, "$key.target_mass"),
        )
    elseif t == "equal"
        pm_raw = get(imf, "particle_mass", nothing)
        if pm_raw !== nothing
            return EqualMassIMF(particle_mass = _typed(pm_raw, Float64, "$key.particle_mass"))
        end
        tm = get(imf, "target_mass", nothing)
        tm === nothing && _config_error(
            "$key.particle_mass",
            "or $key.target_mass is required for type = 'equal'",
        )
        return EqualMassIMF(particle_mass = _typed(tm, Float64, "$key.target_mass") / N)
    elseif t == "powerlaw"
        return PowerLawIMF(
            alpha = _typed(get(imf, "alpha", 2.35), Float64, "$key.alpha"),
            bodyn = bodyn,
            body1 = body1,
        )
    else
        _config_error(
            "$key.type",
            "unknown type '$t'; supported: 'kroupa', 'kroupa_rescaled', 'equal', 'powerlaw'",
        )
    end
end
