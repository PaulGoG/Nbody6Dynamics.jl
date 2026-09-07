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
# The TOML loader accepts two equivalent schemas (decision D6): the flat
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

"Model name used in log messages and metadata TOML."
profile_name(::KingProfile) = "king"
profile_name(::PlummerProfile) = "plummer"

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

"Short label used in metadata TOML and logs."
imf_name(::KroupaIMF) = "kroupa"
imf_name(::RescaledKroupaIMF) = "kroupa_rescaled"
imf_name(::EqualMassIMF) = "equal"

# -----------------------------------------------------------------------------
# ClusterSpec
# -----------------------------------------------------------------------------

"""
    ClusterSpec(; N, rbar, profile, imf, position = Float64[], velocity = Float64[])

Specification for a single star cluster.

# Fields
- `N::Int`: number of bodies to sample.
- `rbar::Float64`: target half-mass radius [pc].
- `profile::DensityProfile`: spatial sampler ([`KingProfile`](@ref) or
  [`PlummerProfile`](@ref)).
- `imf::IMFSpec`: mass-function sampler ([`KroupaIMF`](@ref),
  [`RescaledKroupaIMF`](@ref), or [`EqualMassIMF`](@ref)).
- `position::Vector{Float64}`: centre-of-mass position [pc] (length 3). For
  `orbit_mode = "kepler"` this is auto-computed and may be left empty.
- `velocity::Vector{Float64}`: centre-of-mass velocity in internal code units [km s⁻¹]
  (G = 1 with masses in M☉ and lengths in pc, so 1 unit ≈ 0.0656 km/s). For
  `orbit_mode = "kepler"` this is auto-computed and may be left empty.

The keyword constructor takes the structured (tag-typed) fields only; the
flat TOML schema is the concern of [`load_merger_config`](@ref), which
translates it into these types.
"""
struct ClusterSpec
    N::Int
    rbar::Float64
    profile::DensityProfile
    imf::IMFSpec
    position::Vector{Float64}
    velocity::Vector{Float64}
end

function ClusterSpec(;
    N::Int = 50000,
    rbar::Real = 2.0,
    profile::DensityProfile = KingProfile(),
    imf::IMFSpec = KroupaIMF(),
    position::AbstractVector = Float64[],
    velocity::AbstractVector = Float64[],
)
    return ClusterSpec(N, Float64(rbar), profile, imf, Float64.(position), Float64.(velocity))
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

# -----------------------------------------------------------------------------
# Orbit and output specs
# -----------------------------------------------------------------------------

"""
    OrbitSpec(; apocentre = 15.0, eccentricity = 0.7)

Orbital parameters for the 2-cluster Kepler convenience mode. Ignored for
`orbit_mode = "explicit"`.
"""
Base.@kwdef struct OrbitSpec
    apocentre::Float64 = 15.0
    eccentricity::Float64 = 0.7
end

"""
    MergerOutputSpec(; format = "nbody", truncate_jacobi = true,
                       output_dir = ".", tcrit = 100.0, dtadj = 1.0, deltat = 1.0)

Output and integration parameters for merger ICs.

# Fields
- `format::String`: `"nbody"` (KZ(22)=2, N-body units) — the only supported format.
- `truncate_jacobi::Bool`: truncate each cluster at the nearest-neighbour
  Jacobi radius before combining.
- `output_dir::String`: directory for generated files.
- `tcrit::Float64`: simulation end time (NB units when `format = "nbody"`).
- `dtadj::Float64`, `deltat::Float64`: adjustment and snapshot intervals.
"""
Base.@kwdef struct MergerOutputSpec
    format::String = "nbody"
    truncate_jacobi::Bool = true
    output_dir::String = "."
    tcrit::Float64 = 100.0
    dtadj::Float64 = 1.0
    deltat::Float64 = 1.0
end

"""
    Nbody6ParameterSpec(; qe = 2.0e-4, etai = 0.02, etar = 0.02, etau = 0.1,
                          nnbopt = 0, rs0 = 0.0, rmin = 0.0, dtmin = 0.0,
                          eclose = 1.0, gmin = 1.0e-6, gmax = 0.01, smax = 1.0,
                          kz16 = 0, tcomp = 1.0e8, tcrtp0 = 3600.0,
                          isernb = 40, iserreg = 40, iserks = 0,
                          nfix = 1, ncrit = 10, nrun = 1, ncomm = 10,
                          kz = Dict{Int,Int}())

Every numerical field of the `&INNBODY6` and `&ININPUT` namelists that the
merger input writer emits (`[merger.nbody6]`), in the N-body units of the
combined system. A zero for `nnbopt`, `rs0`, `rmin`, or `dtmin` means
"derive from the member clusters at generation time"
([`resolve_nbody6_parameters`](@ref)). Meanings follow the Nbody6++GPU
input manual.

# Fields
- `qe`: energy-error tolerance per adjustment interval; > 0
- `etai`, `etar`, `etau`: irregular, regular, and regularised time-step
  factors; > 0
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
- `ncomm`: multiplier of `DELTAT` for the restart (`COMMON`) dump interval; ≥ 1
- `kz`: explicit `KZ(i) = v` overrides applied last (`[merger.nbody6.kz]`,
  keys 1–50); an override of an index that also has a named key (14, 16,
  19) warns
"""
Base.@kwdef struct Nbody6ParameterSpec
    qe::Float64 = 2.0e-4
    etai::Float64 = 0.02
    etar::Float64 = 0.02
    etau::Float64 = 0.1
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
    ncomm::Int = 10
    kz::Dict{Int,Int} = Dict{Int,Int}()
end

"""
    StellarSpec(; kz19 = 3, level = "C", zmet = 0.001, epoch0 = 0.0, dtplot = 1.0)

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
"""
Base.@kwdef struct StellarSpec
    kz19::Int = 3
    level::String = "C"
    zmet::Float64 = 0.001
    epoch0::Float64 = 0.0
    dtplot::Float64 = 1.0
end

"""
    TidalSpec(; kz14 = 0, gmg = 0.0, rg0 = 0.0, rg = zeros(3), vg = zeros(3))

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

The `&INSCALE` tidal radius stays 0 so the engine derives it from the field
with the generator's `RBAR` (a non-zero value would override `RBAR`).
"""
Base.@kwdef struct TidalSpec
    kz14::Int = 0
    gmg::Float64 = 0.0
    rg0::Float64 = 0.0
    rg::Vector{Float64} = zeros(3)
    vg::Vector{Float64} = zeros(3)
end

# -----------------------------------------------------------------------------
# MergerConfig
# -----------------------------------------------------------------------------

"""
    MergerConfig(clusters, orbit_mode, orbit, output[, nbody6, stellar, tidal, seed])

Top-level configuration for multi-cluster merger initial conditions.

# Orbit modes

- `"kepler"` — exactly 2 clusters placed on a Keplerian orbit at apocentre.
  Requires an [`OrbitSpec`](@ref). Cluster positions/velocities are computed
  automatically.
- `"explicit"` — any number of clusters (≥ 2). Each cluster must specify
  `position = [x, y, z]` and `velocity = [vx, vy, vz]` in its TOML block.
  The `orbit` section is ignored (with a warning if it appears in TOML).

# Fields
- `clusters::Vector{ClusterSpec}`
- `orbit_mode::String`: `"kepler"` or `"explicit"`
- `orbit::OrbitSpec`
- `output::MergerOutputSpec`
- `nbody6::Nbody6ParameterSpec`: integration and run-control parameters
- `stellar::StellarSpec`: stellar-evolution settings
- `tidal::TidalSpec`: external galactic field
- `seed::Union{Int, Nothing}`: RNG seed. `nothing` = non-deterministic.
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
end

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
) = MergerConfig(clusters, String(orbit_mode), orbit, output, nbody6, stellar, tidal, seed)

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
- `cluster_ranges::Vector{UnitRange{Int}}`
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
    cluster_ranges::Vector{UnitRange{Int}}
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
schema and the structured metadata schema interchangeably. Parsed values are
validated fail-fast (N ≥ 2, rbar > 0, W0 > 0, IMF mass bounds, Kepler orbit
parameters, positive integration intervals); violations raise an error
naming the offending `merger.*` key.

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
    haskey(raw, "merger") || error("TOML file missing [merger] section: $path")
    m = raw["merger"]

    n_clusters = get(m, "n_clusters", 2)::Int
    orbit_mode = get(m, "orbit_mode", "kepler")::String
    orbit_mode in ("kepler", "explicit") ||
        error("orbit_mode must be \"kepler\" or \"explicit\", got \"$orbit_mode\"")

    clusters = ClusterSpec[]
    for i in 1:n_clusters
        key = "cluster$i"
        haskey(m, key) || error("Missing [merger.$key] in config")
        push!(clusters, _parse_cluster_table(m[key], i, orbit_mode))
    end

    if orbit_mode == "kepler"
        n_clusters == 2 || error("Kepler orbit mode requires exactly 2 clusters, got $n_clusters")
    end

    orbit_raw = get(m, "orbit", Dict{String,Any}())
    if orbit_mode == "explicit" && !isempty(orbit_raw)
        @warn "[merger.orbit] is defined but orbit_mode=\"explicit\" — the orbit section will be ignored."
    end
    orbit = OrbitSpec(;
        apocentre = Float64(get(orbit_raw, "apocentre", 15.0)),
        eccentricity = Float64(get(orbit_raw, "eccentricity", 0.7)),
    )

    out_raw = get(m, "output", Dict{String,Any}())
    output = MergerOutputSpec(;
        format = get(out_raw, "format", "nbody")::String,
        truncate_jacobi = get(out_raw, "truncate_jacobi", true)::Bool,
        output_dir = get(out_raw, "output_dir", ".")::String,
        tcrit = Float64(get(out_raw, "tcrit", 100.0)),
        dtadj = Float64(get(out_raw, "dtadj", 1.0)),
        deltat = Float64(get(out_raw, "deltat", 1.0)),
    )

    nb_raw = get(m, "nbody6", Dict{String,Any}())
    kz_raw = get(nb_raw, "kz", Dict{String,Any}())
    kz_over = Dict{Int,Int}()
    for (k, v) in kz_raw
        idx = tryparse(Int, String(k))
        idx === nothing && error("config: merger.nbody6.kz keys must be integers 1–50; got \"$k\"")
        v isa Integer || error("config: merger.nbody6.kz[\"$k\"] must be an integer; got $v")
        kz_over[idx] = Int(v)
    end
    nbody6 = Nbody6ParameterSpec(;
        qe = Float64(get(nb_raw, "qe", 2.0e-4)),
        etai = Float64(get(nb_raw, "etai", 0.02)),
        etar = Float64(get(nb_raw, "etar", 0.02)),
        etau = Float64(get(nb_raw, "etau", 0.1)),
        nnbopt = Int(get(nb_raw, "nnbopt", 0)),
        rs0 = Float64(get(nb_raw, "rs0", 0.0)),
        rmin = Float64(get(nb_raw, "rmin", 0.0)),
        dtmin = Float64(get(nb_raw, "dtmin", 0.0)),
        eclose = Float64(get(nb_raw, "eclose", 1.0)),
        gmin = Float64(get(nb_raw, "gmin", 1.0e-6)),
        gmax = Float64(get(nb_raw, "gmax", 0.01)),
        smax = Float64(get(nb_raw, "smax", 1.0)),
        kz16 = Int(get(nb_raw, "kz16", 0)),
        tcomp = Float64(get(nb_raw, "tcomp", 1.0e8)),
        tcrtp0 = Float64(get(nb_raw, "tcrtp0", 3600.0)),
        isernb = Int(get(nb_raw, "isernb", 40)),
        iserreg = Int(get(nb_raw, "iserreg", 40)),
        iserks = Int(get(nb_raw, "iserks", 0)),
        nfix = Int(get(nb_raw, "nfix", 1)),
        ncrit = Int(get(nb_raw, "ncrit", 10)),
        nrun = Int(get(nb_raw, "nrun", 1)),
        ncomm = Int(get(nb_raw, "ncomm", 10)),
        kz = kz_over,
    )
    st_raw = get(m, "stellar", Dict{String,Any}())
    stellar = StellarSpec(;
        kz19 = Int(get(st_raw, "kz19", 3)),
        level = String(get(st_raw, "level", "C")),
        zmet = Float64(get(st_raw, "zmet", 0.001)),
        epoch0 = Float64(get(st_raw, "epoch0", 0.0)),
        dtplot = Float64(get(st_raw, "dtplot", 1.0)),
    )

    td_raw = get(m, "tidal", Dict{String,Any}())
    rg_raw = get(td_raw, "rg", [0.0, 0.0, 0.0])
    vg_raw = get(td_raw, "vg", [0.0, 0.0, 0.0])
    (length(rg_raw) == 3 && length(vg_raw) == 3) ||
        error("config: merger.tidal.rg and merger.tidal.vg must have length 3")
    tidal = TidalSpec(;
        kz14 = Int(get(td_raw, "kz14", 0)),
        gmg = Float64(get(td_raw, "gmg", 0.0)),
        rg0 = Float64(get(td_raw, "rg0", 0.0)),
        rg = Float64.(rg_raw),
        vg = Float64.(vg_raw),
    )

    seed_raw = get(m, "seed", nothing)
    seed::Union{Int,Nothing} = if seed_raw === nothing
        nothing
    else
        Int(seed_raw)   # any integer is a real seed, including 0
    end

    # Fail-fast validation of the parsed values (bounds documented in the
    # config comments and the manual). kepler_velocity re-checks the orbit
    # parameters at run time for programmatically constructed configs.
    for (i, spec) in enumerate(clusters)
        _validate_cluster_spec(spec, i)
    end
    if orbit_mode == "kepler"
        orbit.apocentre > 0 ||
            error("config: merger.orbit.apocentre must be > 0; got $(orbit.apocentre)")
        (0 ≤ orbit.eccentricity < 1) || error(
            "config: merger.orbit.eccentricity must satisfy 0 ≤ e < 1; got $(orbit.eccentricity)",
        )
    end
    output.tcrit > 0 || error("config: merger.output.tcrit must be > 0; got $(output.tcrit)")
    output.dtadj > 0 || error("config: merger.output.dtadj must be > 0; got $(output.dtadj)")
    output.deltat > 0 || error("config: merger.output.deltat must be > 0; got $(output.deltat)")
    _validate_nbody6(nbody6)
    _validate_stellar(stellar, output)
    _validate_tidal(tidal)
    _validate_tidal_tolerance(tidal, nbody6)

    return MergerConfig(clusters, orbit_mode, orbit, output, nbody6, stellar, tidal, seed)
end

# Energy tolerance below which a tidal-field run halts at its first adjustment:
# this fork never evaluates the tidal potential energy (the XTRNLV call in
# energy.F is commented out), so the engine's DE measures the tidal work.
const _QE_MIN_TIDAL = 0.01

"""
    _validate_tidal_tolerance(t::TidalSpec, p::Nbody6ParameterSpec)

Refuse an external field with an energy tolerance below `$(_QE_MIN_TIDAL)`:
the engine's energy bookkeeping omits the tidal potential energy, so the
relative energy change it checks against `QE` is the work done by the field
(10⁻⁴ to 10⁻² per adjustment interval for realistic fields) and the run
would halt at its first adjustment.
"""
function _validate_tidal_tolerance(t::TidalSpec, p::Nbody6ParameterSpec)
    if t.kz14 > 0 && p.qe < _QE_MIN_TIDAL
        error(
            "config: merger.nbody6.qe = $(p.qe) with merger.tidal.kz14 = $(t.kz14): the engine does " *
            "not evaluate the tidal potential energy (energy.F), so its energy check measures the " *
            "tidal work and halts the run at the first adjustment; set qe ≥ $(_QE_MIN_TIDAL) for " *
            "tidal-field runs and assess the integration accuracy against an isolated control run",
        )
    end
    return nothing
end

"""Fail-fast checks of `[merger.tidal]`: supported `kz14` and the parameters each option needs."""
function _validate_tidal(t::TidalSpec)
    if t.kz14 in (3, 4)
        error(
            "config: merger.tidal.kz14 = $(t.kz14) is not supported: on these options the engine " *
            "rescales all velocities to the &INSCALE virial ratio including the external potential, " *
            "which destroys the prescribed orbital kinematics; use 1, 2, or 5",
        )
    end
    t.kz14 in (0, 1, 2, 5) ||
        error("config: merger.tidal.kz14 must be one of 0, 1, 2, 5; got $(t.kz14)")
    if t.kz14 == 2
        t.gmg > 0 || error("config: merger.tidal.gmg must be > 0 [M☉] for kz14 = 2; got $(t.gmg)")
        t.rg0 > 0 || error("config: merger.tidal.rg0 must be > 0 [kpc] for kz14 = 2; got $(t.rg0)")
    elseif t.kz14 == 5
        any(!=(0.0), t.rg) || error("config: merger.tidal.rg must be non-zero [kpc] for kz14 = 5")
        any(!=(0.0), t.vg) || error("config: merger.tidal.vg must be non-zero [km/s] for kz14 = 5")
    end
    return nothing
end

"""Fail-fast bounds of `[merger.stellar]` (the engine's own limits on `zmet`; `level` one of A, B, C, 0)."""
function _validate_stellar(s::StellarSpec, output::MergerOutputSpec)
    s.kz19 ≥ 0 || error("config: merger.stellar.kz19 must be ≥ 0; got $(s.kz19)")
    s.level in ("A", "B", "C", "0") || error(
        "config: merger.stellar.level must be one of \"A\", \"B\", \"C\", \"0\"; got \"$(s.level)\"",
    )
    (0.0001 ≤ s.zmet ≤ 0.03) ||
        error("config: merger.stellar.zmet must satisfy 0.0001 ≤ zmet ≤ 0.03; got $(s.zmet)")
    s.epoch0 ≤ 0 || error("config: merger.stellar.epoch0 must be ≤ 0 [Myr]; got $(s.epoch0)")
    s.dtplot > 0 || error("config: merger.stellar.dtplot must be > 0 [NB]; got $(s.dtplot)")
    s.dtplot ≥ output.deltat || error(
        "config: merger.stellar.dtplot must be ≥ merger.output.deltat; got dtplot = $(s.dtplot), deltat = $(output.deltat)",
    )
    return nothing
end

"""Fail-fast bounds of `[merger.nbody6]`: positive tolerances and step factors, non-negative derivable entries, `kz16` in 0–3."""
function _validate_nbody6(p::Nbody6ParameterSpec)
    p.qe > 0 || error("config: merger.nbody6.qe must be > 0; got $(p.qe)")
    p.etai > 0 || error("config: merger.nbody6.etai must be > 0; got $(p.etai)")
    p.etar > 0 || error("config: merger.nbody6.etar must be > 0; got $(p.etar)")
    p.nnbopt ≥ 0 || error("config: merger.nbody6.nnbopt must be ≥ 0 (0 = derived); got $(p.nnbopt)")
    p.rs0 ≥ 0 || error("config: merger.nbody6.rs0 must be ≥ 0 (0 = derived); got $(p.rs0)")
    p.rmin ≥ 0 || error("config: merger.nbody6.rmin must be ≥ 0 (0 = derived); got $(p.rmin)")
    p.dtmin ≥ 0 || error("config: merger.nbody6.dtmin must be ≥ 0 (0 = derived); got $(p.dtmin)")
    p.kz16 in 0:3 || error("config: merger.nbody6.kz16 must be one of 0, 1, 2, 3; got $(p.kz16)")
    p.etau > 0 || error("config: merger.nbody6.etau must be > 0; got $(p.etau)")
    p.eclose > 0 || error("config: merger.nbody6.eclose must be > 0; got $(p.eclose)")
    (0 < p.gmin < p.gmax) || error(
        "config: merger.nbody6 must satisfy 0 < gmin < gmax; got gmin = $(p.gmin), gmax = $(p.gmax)",
    )
    p.smax > 0 || error("config: merger.nbody6.smax must be > 0; got $(p.smax)")
    p.tcomp > 0 || error("config: merger.nbody6.tcomp must be > 0 [s]; got $(p.tcomp)")
    p.tcrtp0 > 0 || error("config: merger.nbody6.tcrtp0 must be > 0 [Myr]; got $(p.tcrtp0)")
    for k in (:isernb, :iserreg, :iserks)
        getfield(p, k) ≥ 0 || error("config: merger.nbody6.$k must be ≥ 0; got $(getfield(p, k))")
    end
    for k in (:nfix, :ncrit, :nrun, :ncomm)
        getfield(p, k) ≥ 1 || error("config: merger.nbody6.$k must be ≥ 1; got $(getfield(p, k))")
    end
    for (i, v) in p.kz
        1 ≤ i ≤ 50 || error("config: merger.nbody6.kz index must be within 1–50; got $i")
    end
    return nothing
end

# Parse one [merger.clusterN] table handling both flat and structured forms.
function _parse_cluster_table(c::AbstractDict, idx::Int, orbit_mode::AbstractString)::ClusterSpec
    N = Int(get(c, "N", 50000))
    rbar = Float64(get(c, "rbar", 2.0))

    pos_raw = get(c, "position", Float64[])
    vel_raw = get(c, "velocity", Float64[])
    position = isempty(pos_raw) ? Float64[] : Float64.(pos_raw)
    velocity = isempty(vel_raw) ? Float64[] : Float64.(vel_raw)

    if orbit_mode == "explicit"
        length(position) == 3 ||
            error("[merger.cluster$idx] requires position = [x, y, z] in explicit mode")
        length(velocity) == 3 ||
            error("[merger.cluster$idx] requires velocity = [vx, vy, vz] in explicit mode")
    end

    profile = if haskey(c, "profile") && c["profile"] isa AbstractDict
        _parse_profile_table(c["profile"], idx)
    else
        model = get(c, "model", "king")::String
        if model == "king"
            KingProfile(W0 = Float64(get(c, "W0", 6.0)))
        elseif model == "plummer"
            haskey(c, "W0") && @warn "[merger.cluster$idx] W0 set but model=\"plummer\"; ignored."
            PlummerProfile()
        else
            error("[merger.cluster$idx] unknown model '$model'. Supported: 'king', 'plummer'.")
        end
    end

    imf = if haskey(c, "imf") && c["imf"] isa AbstractDict
        _parse_imf_table(c["imf"], N, idx)
    else
        imf_str = get(c, "imf", "kroupa")::String
        bodyn = Float64(get(c, "bodyn", 0.08))
        body1 = Float64(get(c, "body1", 100.0))
        mt_raw = get(c, "mass_total", nothing)
        mass_total::Union{Float64,Nothing} = mt_raw === nothing ? nothing : Float64(mt_raw)

        if imf_str == "kroupa"
            mass_total === nothing ? KroupaIMF(bodyn = bodyn, body1 = body1) :
            RescaledKroupaIMF(bodyn = bodyn, body1 = body1, target_mass = mass_total)
        elseif imf_str == "kroupa_rescaled"
            mass_total === nothing &&
                error("[merger.cluster$idx] imf=\"kroupa_rescaled\" requires mass_total")
            RescaledKroupaIMF(bodyn = bodyn, body1 = body1, target_mass = mass_total)
        elseif imf_str == "equal"
            mass_total === nothing &&
                error("[merger.cluster$idx] imf=\"equal\" requires mass_total")
            EqualMassIMF(particle_mass = mass_total / N)
        else
            error(
                "[merger.cluster$idx] unknown imf '$imf_str'. Supported: 'kroupa', 'kroupa_rescaled', 'equal'.",
            )
        end
    end

    return ClusterSpec(N, rbar, profile, imf, position, velocity)
end

# ---------------------------------------------------------------------------
# Fail-fast validation of parsed cluster specifications
# ---------------------------------------------------------------------------

"""
    _validate_cluster_spec(spec::ClusterSpec, idx::Int)

Validate one parsed `ClusterSpec` against the documented bounds; raises an
`ErrorException` naming the offending `merger.cluster<idx>.<key>`.
"""
function _validate_cluster_spec(spec::ClusterSpec, idx::Int)
    key = "merger.cluster$idx"
    spec.N ≥ 2 || error("config: $key.N must be ≥ 2; got $(spec.N)")
    spec.rbar > 0 || error("config: $key.rbar must be > 0; got $(spec.rbar)")
    _validate_profile(spec.profile, key)
    _validate_imf(spec.imf, spec.N, key)
    return nothing
end

"Validate a density-profile tag; `key` names the owning cluster table."
function _validate_profile(p::KingProfile, key::String)
    p.W0 > 0 || error("config: $key.W0 must be > 0; got $(p.W0)")
    return nothing
end
_validate_profile(::PlummerProfile, ::String) = nothing

"Validate an IMF tag; `key` names the owning cluster table."
function _validate_imf(i::KroupaIMF, ::Int, key::String)
    (0 < i.bodyn < i.body1) || error(
        "config: $key.imf must satisfy 0 < bodyn < body1; " *
        "got bodyn = $(i.bodyn), body1 = $(i.body1)",
    )
    return nothing
end
# Rescale factor c = target_mass / (N ⟨m⟩) of the rescaled Kroupa mode: bodies
# stop corresponding to stars beyond a factor of two either way (refused);
# beyond 1.4 the stellar-evolution output is already unreliable (warned).
const _IMF_RESCALE_WARN = 1.4
const _IMF_RESCALE_MAX = 2.0

function _validate_imf(i::RescaledKroupaIMF, N::Int, key::String)
    (0 < i.bodyn < i.body1) || error(
        "config: $key.imf must satisfy 0 < bodyn < body1; " *
        "got bodyn = $(i.bodyn), body1 = $(i.body1)",
    )
    i.target_mass > 0 || error(
        "config: $key.imf target_mass (flat form: mass_total) must be > 0; " *
        "got $(i.target_mass)",
    )
    mean_m = kroupa_mean_mass(i.bodyn, i.body1)
    c = i.target_mass / (N * mean_m)
    if !(1 / _IMF_RESCALE_MAX ≤ c ≤ _IMF_RESCALE_MAX)
        error(
            "config: $key requests a Kroupa rescale factor of ×$(round(c; digits = 2)) " *
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
    i.particle_mass > 0 || error(
        "config: $key.imf particle_mass (flat form: mass_total / N) must be > 0; " *
        "got $(i.particle_mass)",
    )
    return nothing
end

function _parse_profile_table(p::AbstractDict, idx::Int)::DensityProfile
    t = get(p, "type", "king")::String
    if t == "king"
        KingProfile(W0 = Float64(get(p, "W0", 6.0)))
    elseif t == "plummer"
        PlummerProfile()
    else
        error("[merger.cluster$idx.profile] unknown type '$t'. Supported: 'king', 'plummer'.")
    end
end

function _parse_imf_table(imf::AbstractDict, N::Int, idx::Int)::IMFSpec
    t = get(imf, "type", "kroupa")::String

    bounds_raw = get(imf, "bounds", nothing)
    bodyn, body1 = if bounds_raw !== nothing
        length(bounds_raw) == 2 || error("[merger.cluster$idx.imf] 'bounds' must have length 2")
        Float64(bounds_raw[1]), Float64(bounds_raw[2])
    else
        Float64(get(imf, "bodyn", 0.08)), Float64(get(imf, "body1", 100.0))
    end

    if t == "kroupa"
        return KroupaIMF(bodyn = bodyn, body1 = body1)
    elseif t == "kroupa_rescaled"
        tm = get(imf, "target_mass", nothing)
        tm === nothing &&
            error("[merger.cluster$idx.imf] type='kroupa_rescaled' requires target_mass")
        return RescaledKroupaIMF(bodyn = bodyn, body1 = body1, target_mass = Float64(tm))
    elseif t == "equal"
        pm_raw = get(imf, "particle_mass", nothing)
        if pm_raw !== nothing
            return EqualMassIMF(particle_mass = Float64(pm_raw))
        end
        tm = get(imf, "target_mass", nothing)
        tm === nothing &&
            error("[merger.cluster$idx.imf] type='equal' requires 'particle_mass' or 'target_mass'")
        return EqualMassIMF(particle_mass = Float64(tm) / N)
    else
        error(
            "[merger.cluster$idx.imf] unknown type '$t'. Supported: 'kroupa', 'kroupa_rescaled', 'equal'.",
        )
    end
end
