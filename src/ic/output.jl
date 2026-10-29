# =============================================================================
# dat.10 writer and .inp generator for Nbody6++ merger ICs
# =============================================================================

# Velocity unit of the internal (G = 1, M☉, pc) code-unit system:
# sqrt(G M☉ / pc) in km/s. Sampled and Kepler velocities carry this unit
# until the NB-unit conversion for dat.10.
# Code-unit constants (_G_PC_KMS2_MSUN, _CODE_VSTAR_KMS, _CODE_TIME_MYR)
# are defined in orbits.jl.

"""
    crossing_time(M, E) -> Float64

Crossing time `M^{5/2} / (2|E|)^{3/2}` of a self-gravitating system of mass
`M` and total energy `E < 0` in `G = 1` units (the engine's `TCR`,
`scale.F`/`adjust.F`); `NaN` for `E ≥ 0`.
"""
crossing_time(M::Real, E::Real)::Float64 = E < 0 ? M^2.5 / (2 * abs(E))^1.5 : NaN

"""N-body time unit in Myr for a system of total mass `M_total` [M☉] and length unit `rbar_pc` [pc] (the engine's `T*`)."""
_nbody_time_myr(M_total::Real, rbar_pc::Real)::Float64 = _CODE_TIME_MYR * sqrt(rbar_pc^3 / M_total)

"""Stellar mass bounds `(bodyn, body1)` [M☉] written to `&INDATA` for an IMF specification (inert under `KZ(22) = 2`, kept consistent with the sampled population)."""
_imf_mass_bounds(i::KroupaIMF) = (i.bodyn, i.body1)
_imf_mass_bounds(i::RescaledKroupaIMF) = (i.bodyn, i.body1)
_imf_mass_bounds(i::EqualMassIMF) = (i.particle_mass, i.particle_mass)
_imf_mass_bounds(i::PowerLawIMF) = (i.bodyn, i.body1)
_imf_mass_bounds(i::ExplicitMassesIMF) = extrema(i.masses)

"""
    write_dat10(path, mass, pos, vel)

Write a `dat.10` particle file for Nbody6++ (`KZ(22)=2`).

Each line: `MASS  X  Y  Z  VX  VY  VZ`

All values must be in N-body units (G=1, M_total=1).
"""
function write_dat10(
    path::AbstractString,
    mass::Vector{Float64},
    pos::Matrix{Float64},
    vel::Matrix{Float64},
)
    N = length(mass)
    size(pos) == (3, N) || error("pos must be 3×$N, got $(size(pos))")
    size(vel) == (3, N) || error("vel must be 3×$N, got $(size(vel))")

    _backup_existing(path)
    open(path, "w") do io
        for i in 1:N
            @printf(
                io,
                "%.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e\n",
                mass[i],
                pos[1, i],
                pos[2, i],
                pos[3, i],
                vel[1, i],
                vel[2, i],
                vel[3, i]
            )
        end
    end
    @info "Wrote $N particles to $path"
    return nothing
end

"""
    to_nbody_units!(mass, pos, vel, M_total_solar, rbar_pc)

Convert from physical units (M☉, pc, km/s) to N-body units (G = 1,
M_total = 1, length unit `rbar_pc`) in-place.

The callers pass the half-mass radius of the whole initial configuration as
`rbar_pc`, and the engine reads the result with `KZ(22) = 2`, which skips its
rescaling to E = -1/4. The total energy is therefore not -1/4: an isolated
King W₀ = 6 cluster starts near -0.20 (r_V ≈ 1.25 r_h), and a bound pair lower
still, its orbital binding energy adding to the total.

# N-body scaling
- Mass: `m_nb = m_solar / M_total_solar`
- Length: `r_nb = r_pc / rbar_pc` (rbar = length unit in pc)
- Velocity: `v_nb = v_kms / vstar` where `vstar = sqrt(G M_total / rbar)`
  in km/s, i.e. `_CODE_VSTAR_KMS × sqrt(M_total / rbar)`.
"""
function to_nbody_units!(
    mass::Vector{Float64},
    pos::Matrix{Float64},
    vel::Matrix{Float64},
    M_total_solar::Float64,
    rbar_pc::Float64,
)
    vstar = _CODE_VSTAR_KMS * sqrt(M_total_solar / rbar_pc)
    mass ./= M_total_solar
    pos ./= rbar_pc
    vel ./= vstar
    return nothing
end

"""
    _central_density_contrast(profile::DensityProfile) -> Float64

Central density over the mean density inside the half-mass radius,
`ρ₀ / ρ̄(< r_h)`. Plummer: `2 (r_h/a)³ ≈ 4.45` analytically. King: from the
solved model, with `r_h` located on the cumulative-mass profile. Rotating
King: from the solved axisymmetric model (`ρ₀ = 1` in its units, `M̂` and the
spherical half-mass radius `r_h` of [`RotatingKingModel`](@ref)).
"""
_central_density_contrast(::PlummerProfile) = 2.0 * _PLUMMER_RHM_OVER_A^3
function _central_density_contrast(p::KingProfile)
    rhat, ρ, m_cum = _king_cumulative(p.W0)
    i_h = findfirst(≥(0.5 * m_cum[end]), m_cum)
    r_h = rhat[i_h]
    ρ_mean = 0.5 * m_cum[end] / (4π / 3 * r_h^3)
    return ρ[1] / ρ_mean
end
function _central_density_contrast(p::RotatingKingProfile)
    m = _rotating_king_cached(p.W0, p.omega0)
    return 1.0 / (0.5 * m.mass / (4π / 3 * m.r_h^3))
end

"""
    _engine_digit_counter_terminates(dt) -> Bool

Whether the engine's `string_left.f` terminates on the interval `dt`. That
routine, called at every output (`output.F`), COMMON dump (`mydump.F`) and
stellar-evolution record (`hrplot.F`) with `DELTAT`, `DTADJ` and
`DTPLOT`, counts the decimal digits of `dt` by multiplying by ten until
`dtmp - int(dtmp)` is zero, with a default-kind `int`. A value whose
binary representation never lands on an integer under repeated
multiplication (`0.6302`, `0.1576`, …) grows past 2³¹, where the
conversion overflows, and the loop never ends: the engine sits at 100 %
CPU inside its first output and never prints a second adjustment.
Emulated here in the same arithmetic.
"""
function _engine_digit_counter_terminates(dt::Real)::Bool
    dtmp = Float64(dt)
    ip = 0
    while dtmp - trunc(dtmp) != 0.0
        dtmp ≥ 2.0^31 && return false
        ip += 1
        dtmp *= 10.0
        ip > 60 && return false
    end
    return true
end

"""
    engine_interval(dt; resolution = 1/128, max_digits = 9) -> Float64

`dt` rounded to the nearest dyadic rational `m 2⁻ᵏ` whose spacing `2⁻ᵏ` is
at most `resolution × dt`, with `k ≤ max_digits`. Such a value has an
exact decimal expansion of `k` digits, so the engine's digit counter
([`_engine_digit_counter_terminates`](@ref)) stops after `k` steps, and it
is commensurate with the engine's power-of-two block time steps. The cap
of nine digits is the engine's as well: `string_left.f` formats a count of
at most ten with an `I1` descriptor, so exactly ten digits produce an
invalid format and a runtime error at the next output. Applied to
`DTADJ`, `DELTAT` and `DTPLOT` before they are written to `merger.inp`;
the change is below `resolution / 2` unless `max_digits` binds (intervals
below about `2⁻ᵐᵃˣ⁻ᵈⁱᵍⁱᵗˢ × 128`). `dt ≤ 0` is returned unchanged.

# Example
```julia
engine_interval(0.6302)   # 0.62890625 = 161/256
engine_interval(0.5)      # 0.5
```
"""
function engine_interval(dt::Real; resolution::Real = 1 / 128, max_digits::Integer = 9)::Float64
    dt > 0 || return Float64(dt)
    k = clamp(ceil(Int, -log2(resolution * dt)), 0, Int(max_digits))
    h = 2.0^-k
    return max(round(dt / h), 1.0) * h
end

"""[`engine_interval`](@ref) of `dt`, logging the change when it exceeds rounding noise."""
function _engine_interval_logged(name::AbstractString, dt::Real)::Float64
    v = engine_interval(dt)
    if dt > 0 && !isapprox(v, dt; rtol = 1e-12)
        @info @sprintf(
            "%s = %s NB written as the dyadic %s (requested %.6g; the engine's interval digit counter needs an exact decimal)",
            name,
            _decimal_string(v),
            _decimal_string(v),
            dt
        )
    end
    return v
end

"""
    _decimal_string(x) -> String

Exact decimal text of a dyadic rational with at most ten binary digits
(`0.62890625`, `0.5`, `1`; [`engine_interval`](@ref) stops at nine), for
the interval fields of `merger.inp`: a fixed-precision format would
re-round the value and defeat the rounding.
"""
function _decimal_string(x::Real)::String
    s = @sprintf("%.10f", x)
    s = rstrip(s, '0')
    return String(rstrip(s, '.'))
end

"""
    resolve_nbody6_parameters(spec, clusters, cluster_ranges, N_total, rbar_pc)
        -> Nbody6ParameterSpec

Fill the zero entries of `spec` from its `preset` and the member clusters, in
N-body units of the combined system (length unit `rbar_pc`). With `r_h` the
smallest member half-mass radius in those units, `N_min` the smallest
post-truncation membership, and `ρ̂` the central density contrast of that
cluster's profile ([`_central_density_contrast`](@ref)):

| Key | `"engine"` | `"baseline"` |
|---|---|---|
| `ETAI` | 0.02 | 0.01 |
| `ETAR` | 0.02 | 0.02 |
| `ETAU` | 0.1 | 0.1 |
| `NNBOPT` | `clamp(round(√N_total), 20, 300)` | 100 |
| `RS0` | `min(2 r_h (2 NNBOPT / N_min)^{1/3}, r_h)` | the same |
| `RMIN` | `4 r_h / (N_min ρ̂^{1/3})` | `12 r_h / N_min` |
| `DTMIN` | `0.04 √(ETAI/0.02) √(RMIN³ N_total)` | `0.016 √(ETAI/0.02) √(RMIN³ N_total)` |

- `RS0` is twice the radius enclosing about `NNBOPT` stars at the mean
  density of the half-mass sphere, capped at `r_h`. The factor 2 matches the
  engine's own example inputs (`RS0 ≈ 0.5 r_h` for `NNBOPT ≈ √N`). The
  engine regrows an empty neighbour sphere itself (`fpoly0.F`), so the
  radius is a matter of start-up cost, not of correctness; the start-up hang
  once attributed to it was the interval digit counter
  ([`engine_interval`](@ref)).
- `RMIN` of `"engine"` is the functional form of the engine's own
  re-derivation (`adjust.F`) evaluated for the member cluster instead of the
  whole configuration; that of `"baseline"` is three times the 90°
  deflection scale `4 r_h / N_min`.
- `DTMIN` is the Aarseth close-encounter step with the mean body mass
  `1/N_total` (`"engine"`), or 0.4 of it (`"baseline"`), at the resolved
  `RMIN`; it uses the resolved `ETAI`.

User-supplied (nonzero) values are kept. The result has no zero entries.
"""
function resolve_nbody6_parameters(
    spec::Nbody6ParameterSpec,
    clusters::Vector{ClusterSpec},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    N_total::Int,
    rbar_pc::Float64,
)::Nbody6ParameterSpec
    isempty(clusters) && throw(ArgumentError("resolve_nbody6_parameters: no clusters"))
    rbar_pc > 0 || throw(ArgumentError("resolve_nbody6_parameters: rbar_pc must be > 0"))
    r_h, i_min = findmin([c.rbar for c in clusters] ./ rbar_pc)
    N_min = minimum(length, cluster_ranges)
    N_min ≥ 2 || throw(
        ArgumentError(
            "resolve_nbody6_parameters: a cluster has fewer than 2 members after truncation",
        ),
    )
    spec.preset in _NBODY6_PRESETS ||
        throw(ArgumentError("resolve_nbody6_parameters: unknown preset \"$(spec.preset)\""))
    engine = spec.preset == "engine"
    etai = spec.etai > 0 ? spec.etai : (engine ? 0.02 : 0.01)
    etar = spec.etar > 0 ? spec.etar : 0.02
    etau = spec.etau > 0 ? spec.etau : 0.1
    nnbopt = if spec.nnbopt > 0
        spec.nnbopt
    else
        engine ? clamp(round(Int, sqrt(N_total)), 20, 300) : 100
    end
    rs0 = spec.rs0 > 0 ? spec.rs0 : min(2 * r_h * cbrt(2 * nnbopt / N_min), r_h)
    ρ̂ = _central_density_contrast(clusters[i_min].profile)
    rmin = if spec.rmin > 0
        spec.rmin
    else
        engine ? 4 * r_h / (N_min * cbrt(ρ̂)) : 12 * r_h / N_min
    end
    dtmin = if spec.dtmin > 0
        spec.dtmin
    else
        (engine ? 0.04 : 0.016) * sqrt(etai / 0.02) * sqrt(rmin^3 * N_total)
    end
    kw = Dict{Symbol,Any}(k => getfield(spec, k) for k in fieldnames(Nbody6ParameterSpec))
    kw[:etai] = etai
    kw[:etar] = etar
    kw[:etau] = etau
    kw[:nnbopt] = nnbopt
    kw[:rs0] = rs0
    kw[:rmin] = rmin
    kw[:dtmin] = dtmin
    return Nbody6ParameterSpec(; kw...)
end

"""`Nbody6ParameterSpec` as a TOML-ready table (the `kz` override keys as strings)."""
function _nbody6_table(p::Nbody6ParameterSpec)::Dict{String,Any}
    d = _struct_to_dict(p)
    d["kz"] = Dict{String,Int}(string(k) => v for (k, v) in p.kz)
    return d
end

"""Throw unless every derivable entry of `p` has been resolved (no zeros)."""
function _assert_resolved(p::Nbody6ParameterSpec)
    for k in (:etai, :etar, :etau, :nnbopt, :rs0, :rmin, :dtmin)
        getfield(p, k) > 0 || throw(
            ArgumentError(
                "generate_merger_inp: nbody6.$k is unresolved (0); pass the result of resolve_nbody6_parameters",
            ),
        )
    end
    return nothing
end

# Regime thresholds for the generation-time guards (see the merger docs,
# "Feasibility and limitations").
const _Q_COLD_COLLAPSE = 0.3
const _RBAR_OVER_RHM_UNRESOLVED = 5.0
# Per-pair relative FP32 regular-force error from the coordinate offset above which a warning is given
const _FP32_OFFSET_WARN = 1.0e-4

"""
    _check_multicluster_regime(q_virial, rbar_over_rhm_min, nbody6, r_h_min_nb;
                               deltat = 0.0, t_cr_member_min_nb = NaN)

Generation-time guards against configurations the engine mis-integrates or
mis-diagnoses: refuse a neighbour radius `rs0` wider than the smallest
member half-mass radius (every neighbour list would overflow at start-up);
warn when the combined virial ratio is below `$(_Q_COLD_COLLAPSE)`
(cold-collapse regime, global diagnostics meaningless until the remnant
forms), when the length unit exceeds `$(_RBAR_OVER_RHM_UNRESOLVED)`
member half-mass radii (members unresolved by the single-centre
diagnostics), when the snapshot interval `deltat` exceeds the smallest
member crossing time (member dynamics undersampled), and when the bound
`fp32_offset_bound` = 4u·max|x|/RS0 (u = 2⁻²⁴) on the per-pair relative error
of the single-precision regular force from the coordinate offset exceeds
`$(_FP32_OFFSET_WARN)`.
"""
function _check_multicluster_regime(
    q_virial::Float64,
    rbar_over_rhm_min::Float64,
    nbody6::Nbody6ParameterSpec,
    r_h_min_nb::Float64;
    deltat::Float64 = 0.0,
    t_cr_member_min_nb::Float64 = NaN,
    fp32_offset_bound::Float64 = 0.0,
)
    nbody6.rs0 ≤ r_h_min_nb || error(
        "merger.nbody6.rs0 = $(nbody6.rs0) exceeds the smallest member half-mass radius " *
        "$(round(r_h_min_nb; sigdigits = 3)) (N-body units): a neighbour sphere wider than a " *
        "member cluster overflows every neighbour list at start-up. Set rs0 = 0 to derive it " *
        "or choose a smaller value.",
    )
    if !isnan(q_virial) && q_virial < _Q_COLD_COLLAPSE
        @warn "Combined virial ratio Q = $(round(q_virial; digits = 3)) < $(_Q_COLD_COLLAPSE): " *
              "cold-collapse regime; global infall dominates the evolution and the engine's " *
              "global diagnostics (lagr.7, RC, RTIDE, ADJUST Q) are not meaningful until the " *
              "remnant forms"
    end
    if rbar_over_rhm_min > _RBAR_OVER_RHM_UNRESOLVED
        @warn "RBAR/r_hm,min = $(round(rbar_over_rhm_min; digits = 1)) > $(_RBAR_OVER_RHM_UNRESOLVED): " *
              "member clusters are unresolved by the engine's single-centre diagnostics; use " *
              "the per-cluster snapshot analysis (cluster_ranges) for cluster-level results"
    end
    if fp32_offset_bound > _FP32_OFFSET_WARN
        @warn "FP32 offset bound 4u·max|x|/RS0 = $(round(fp32_offset_bound; sigdigits = 2)) > " *
              "$(_FP32_OFFSET_WARN): the single-precision regular force of a GPU build loses " *
              "accuracy on the outermost bodies; centre the configuration or raise rs0"
    end
    if !isnan(t_cr_member_min_nb) && deltat > t_cr_member_min_nb
        @warn "Snapshot interval deltat = $deltat exceeds the smallest member crossing time " *
              "$(round(t_cr_member_min_nb; sigdigits = 3)) (N-body units): member dynamics are " *
              "undersampled by conf.3; lower merger.output.deltat"
    end
    return nothing
end

"""
    _checkpoint_dump_bytes(n, nnbopt) -> Int

Upper estimate of the size of one restart dump of `n` bodies with target
neighbour number `nnbopt`, in bytes: 576 per body (the per-body records of
the engine's `mydump.F`), 4 per neighbour-list entry with `nnbopt + 1`
entries per body, and 2 MB for the header blocks. The realised mean
neighbour number stays below `nnbopt` (642 MB measured at 6×10⁵ bodies with
`nnbopt = 300`, against 1.07 GB here).
"""
_checkpoint_dump_bytes(n::Integer, nnbopt::Integer)::Int =
    2_000_000 + Int(n) * (576 + 4 * (Int(nnbopt) + 1))

"""
    generate_merger_inp(path, N_total, rbar, zmbar; nbody6, stellar = StellarSpec(),
                        kz14 = 0, tcrit, dtadj, deltat, nrand, mass_bounds = (0.08, 100.0),
                        imbh = nothing)

Generate an Nbody6++ `.inp` file configured for external particle
input via `dat.10` in N-body units (`KZ(22) = 2`, the only supported
input mode).

Follows the Fortran NAMELIST read order expected by `nbody6.F → start.F`:
  `&INNBODY6` → `&ININPUT` → `&INSSE` → `&INBSE` → `&INCOLL` →
  `&INDATA` → `&INSCALE` (→ `&INXTRNL0` if KZ(14) = 2, 5) (→ the
  list-directed `IMBHINIT` line `BIMBH XIMBH(1:3) VIMBH(1:3) DTBH` if
  KZ(24) = 1, or the lines `NIMBHS DTBH` and the `NIMBHS` IMBH names if
  KZ(24) = 3).

`imbh` selects one of two IMBH routes by its `route` field; either sets
`KZ(24)` with the other fixed options (before the `kz` overrides), and an
unknown route raises an `ArgumentError`.

`(route = :kz24, mass_nb, pos_nb, vel_nb, dtbh_nb)`, holding the last
`dat.10` row (the IMBH, N-body units) and the `mbh.45` interval `DTBH` [NB],
sets `KZ(24) = 1` and writes the `IMBHINIT` line in the units the engine
reads it (M☉, pc, km s⁻¹). The engine overwrites particle N with `BODY(N) = BIMBH/ZMBAR`,
`X = XIMBH/RBAR`, `XDOT = VIMBH/VSTAR`, `VSTAR = 1e-5 √(GM/PC) √(ZMASS ZMBAR/RBAR)`,
where under `KZ(22) = 2` its `ZMBAR` is already the mass scale
`ZMBAR N / ZMASS` (`start.F`), not the mean mass written to `&ININPUT`. With
`rbar_w`, `zmbar_w` the values of the `RBAR=%.6f`, `ZMBAR=%.6f` fields and
`zmass_w = zmbar_w N_total` (`ZMASS = 1`), the line carries
`BIMBH = mass_nb zmass_w`, `XIMBH = pos_nb rbar_w` and
`VIMBH = vel_nb VSTAR(zmass_w, rbar_w)`, so the engine's conversion returns
the `dat.10` state to rounding (an identity overwrite).

`(route = :flag, names::Vector{Int}, dtbh_nb)` sets `KZ(24) = 3` (engine
fork ≥ v2026.07.11) and writes, at the same position, the line
`NIMBHS DTBH` (`length(names)` and `dtbh_nb` as `%.15e`) and then the names
separated by single spaces on one line. The engine keeps the `dat.10` mass,
position and velocity of these bodies and flags them as IMBHs (type 14).
`names` must hold 1 to 16 distinct names in `1:N_total`.

An IMBH under stellar evolution (`stellar.kz19 ≥ 3`) needs the Peters decay
of `KZ(27) = 3` on either route; a smaller final `KZ(27)` is refused.

# Key settings
- `nbody6`: resolved [`Nbody6ParameterSpec`](@ref); obtain it from
  [`resolve_nbody6_parameters`](@ref) — unresolved zeros are refused. Its
  `binary_physics = "full"` writes the option set of
  `_FULL_BINARY_KZ` over the inherited block, then its `kz` overrides
  are applied; a resulting `KZ(46) > 0` warns. `checkpoint = true` sets `KZ(2) = 1` (a
  restart dump at every adjustment) and logs the disk space the dumps can
  take.
- `stellar`: [`StellarSpec`](@ref) (`KZ(19)`, `Level`, `ZMET`, `EPOCH0`,
  `DTPLOT`); `KZ(12)` HR diagnostics are switched off with `kz19 = 0`
- `tidal`: [`TidalSpec`](@ref) (`KZ(14)` and the `&INXTRNL0` namelist for
  options 2 and 5)
- `mass_bounds`: `(BODYN, BODY1)` of `&INDATA`, inert under `KZ(22) = 2`
  but kept consistent with the sampled population
"""
function generate_merger_inp(
    path::AbstractString,
    N_total::Int,
    rbar::Float64,
    zmbar::Float64;
    nbody6::Nbody6ParameterSpec,
    stellar::StellarSpec = StellarSpec(),
    tidal::TidalSpec = TidalSpec(),
    tcrit::Float64 = 100.0,
    dtadj::Float64 = 1.0,
    deltat::Float64 = 1.0,
    nrand::Int = 10000,
    mass_bounds::Tuple{Float64,Float64} = (0.08, 100.0),
    nbin0::Int = 0,
    imbh::Union{Nothing,NamedTuple} = nothing,
    mbh::Union{Nothing,NamedTuple} = nothing,
)
    _assert_resolved(nbody6)
    _validate_tidal(tidal)
    nbin0 ≥ 0 || throw(ArgumentError("nbin0 must be ≥ 0; got $nbin0"))
    kz14 = tidal.kz14
    kz = zeros(Int, 50)
    kz[1] = 1
    kz[2] = nbody6.checkpoint ? 1 : -1
    kz[3] = 2
    kz[7] = 3
    kz[8] = nbin0 > 0 ? 2 : 0   # 2: primordial pairs are the first 2·NBIN0 bodies of dat.10
    kz[12] = stellar.kz19 == 0 ? 0 : 1
    kz[14] = kz14
    kz[16] = nbody6.kz16
    kz[19] = stellar.kz19
    kz[22] = 2
    kz[23] = _ESCAPERS_KZ23[nbody6.escapers]
    kz[26] = 1
    kz[30] = 1
    imbh_route = if imbh === nothing
        :none
    else
        hasproperty(imbh, :route) ||
            throw(ArgumentError("imbh: a route field (:kz24 or :flag) is required"))
        imbh.route in (:kz24, :flag) || throw(
            ArgumentError("imbh: unknown route $(repr(imbh.route)); expected :kz24 or :flag"),
        )
        imbh.route
    end
    # kz24: IMBHINIT overwrites particle N; flag: the listed names are the IMBHs
    imbh_route === :kz24 && (kz[24] = 1)
    imbh_route === :flag && (kz[24] = 3)
    mbh === nothing || (kz[44] = 1)
    kz[43] = _POST_NEWTONIAN_KZ43[nbody6.post_newtonian]
    if nbody6.binary_physics == "full"
        for (i, v) in _FULL_BINARY_KZ
            kz[i] = v
        end
    end
    for (i, v) in nbody6.kz
        i in (8, 14, 16, 19) &&
            @warn "merger.nbody6.kz overrides KZ($i) = $v, which also has a named setting (binaries, kz14/tidal, kz16, stellar.kz19)"
        kz[i] = v
    end
    kz[46] > 0 &&
        @warn "KZ(46) = $(kz[46]) writes the per-binary custom output (ASCII, several GB per N-body time unit at 6×10⁵ bodies); meant for test arms"
    # The IMBHINIT line in the units the engine reads it, from the RBAR and
    # ZMBAR it parses and its mass scale ZMBAR N / ZMASS (ZMASS = 1)
    imbh_line = if imbh_route !== :kz24
        nothing
    else
        (length(imbh.pos_nb) == 3 && length(imbh.vel_nb) == 3 && imbh.dtbh_nb > 0) || throw(
            ArgumentError("imbh: pos_nb and vel_nb need three components and dtbh_nb must be > 0"),
        )
        rbar_w = parse(Float64, @sprintf("%.6f", rbar))
        zmbar_w = parse(Float64, @sprintf("%.6f", zmbar))
        zmass_w = zmbar_w * N_total
        vstar_w = _imbh_engine_vstar_kms(zmass_w, rbar_w)
        (
            bimbh = Float64(imbh.mass_nb) * zmass_w,
            ximbh = Float64.(imbh.pos_nb) .* rbar_w,
            vimbh = Float64.(imbh.vel_nb) .* vstar_w,
            dtbh = Float64(imbh.dtbh_nb),
        )
    end
    # The NIMBHS DTBH line and the name line of the flag route
    imbh_flag = if imbh_route !== :flag
        nothing
    else
        names = collect(Int, imbh.names)
        imbh.dtbh_nb > 0 ||
            throw(ArgumentError("imbh: dtbh_nb must be > 0; got $(imbh.dtbh_nb)"))
        1 ≤ length(names) ≤ _IMBH_FLAG_MAX || throw(
            ArgumentError(
                "imbh: the flag route takes 1 to $(_IMBH_FLAG_MAX) names; got $(length(names))",
            ),
        )
        all(n -> 1 ≤ n ≤ N_total, names) ||
            throw(ArgumentError("imbh: names must lie in 1:$N_total; got $names"))
        allunique(names) || throw(ArgumentError("imbh: repeated names in $names"))
        (names = names, dtbh = Float64(imbh.dtbh_nb))
    end
    # The physics refusal after the argument checks
    (imbh !== nothing && stellar.kz19 ≥ 3 && kz[27] < 3) && error(
        "an IMBH under stellar evolution needs the Peters decay of KZ(27) = 3 " *
        "(binary_physics = \"full\" or kz = {27 = 3}); got KZ(27) = $(kz[27])",
    )

    # Build KZ lines: KZ(1:10) = ... etc.
    kz_strs = String[]
    for row in 1:5
        s = (row - 1) * 10 + 1
        push!(kz_strs, "KZ($(s):$(s+9))=$(join(string.(kz[s:s+9]), " "))")
    end

    _backup_existing(path)
    open(path, "w") do io
        # --- 1. &INNBODY6: start/restart, CPU time, checkpointing ---
        println(io, "&INNBODY6")
        @printf(
            io,
            "KSTART=1,TCOMP=%.6G,TCRTP0=%.6G,isernb=%d,iserreg=%d,iserks=%d /\n",
            nbody6.tcomp,
            nbody6.tcrtp0,
            nbody6.isernb,
            nbody6.iserreg,
            nbody6.iserks
        )
        println(io)

        # --- 2. &ININPUT: main simulation parameters + KZ options ---
        println(io, "&ININPUT")
        @printf(
            io,
            "N=%d,NFIX=%d,NCRIT=%d,NRAND=%d,NNBOPT=%d,NRUN=%d,NCOMM=%d,\n",
            N_total,
            nbody6.nfix,
            nbody6.ncrit,
            abs(nrand) % typemax(Int32),
            nbody6.nnbopt,
            nbody6.nrun,
            nbody6.ncomm
        )
        @printf(
            io,
            "ETAI=%.4G,ETAR=%.4G,RS0=%.4G,DTADJ=%s,DELTAT=%s,TCRIT=%.2f,QE=%.3E,RBAR=%.6f,ZMBAR=%.6f,\n",
            nbody6.etai,
            nbody6.etar,
            nbody6.rs0,
            _decimal_string(dtadj),
            _decimal_string(deltat),
            tcrit,
            nbody6.qe,
            rbar,
            zmbar
        )
        println(io, join(kz_strs, "\n"))
        @printf(
            io,
            "DTMIN=%.3E,RMIN=%.3E,ETAU=%.4G,ECLOSE=%.4G,GMIN=%.3E,GMAX=%.4G,SMAX=%.4G,\n",
            nbody6.dtmin,
            nbody6.rmin,
            nbody6.etau,
            nbody6.eclose,
            nbody6.gmin,
            nbody6.gmax,
            nbody6.smax
        )
        println(io, "Level='$(stellar.level)' /")
        println(io)

        # --- 3-5. SSE/BSE/Coll: empty → use Level defaults ---
        println(io, "&INSSE /")
        println(io)
        inbse = String[]
        nbody6.gw_recoil == "off" ||
            push!(inbse, "kicktype=$(_GW_RECOIL_KICKTYPE[nbody6.gw_recoil])")
        nbody6.bh_spin == "zero" || push!(inbse, "bhspinfl=$(_BH_SPIN_FLAG[nbody6.bh_spin])")
        println(io, isempty(inbse) ? "&INBSE /" : "&INBSE " * join(inbse, ",") * " /")
        println(io)
        println(io, "&INCOLL /")
        println(io)

        # --- 6. &INDATA: stellar population ---
        # ALPHAS/BODY1/BODYN are inert under KZ(22)=2 (masses come from
        # dat.10); the bounds are written consistent with the sampled IMF.
        # ZMET, EPOCH0, DTPLOT govern stellar evolution. NBIN0 = number of
        # primordial pairs at the head of dat.10; NHI0 = 0 (no hierarchies).
        println(io, "&INDATA")
        @printf(
            io,
            "ALPHAS=2.35,BODY1=%.4G,BODYN=%.4G,NBIN0=%d,NHI0=0,ZMET=%.4G,EPOCH0=%.4G,DTPLOT=%s /\n",
            mass_bounds[2],
            mass_bounds[1],
            nbin0,
            stellar.zmet,
            stellar.epoch0,
            _decimal_string(stellar.dtplot)
        )
        println(io)

        # --- 7. &INSETUP: placeholder (no Fortran NAMELIST reads this) ---
        println(io, "&INSETUP SEMI=,ECC=,APO=,N2=,SCALE=,ZM1=,ZM2,ZMH,RCUT= /")
        println(io)

        # --- 8. &INSCALE: virial ratio, rotation, tidal ---
        println(io, "&INSCALE")
        println(io, "Q=0.5,VXROT=0.0,VZROT=0.0,RTIDE=0.0 /")
        println(io)

        # --- 9. External potential: read by the engine for KZ(14) = 2 and 5 ---
        if kz14 == 2
            println(io, "&INXTRNL0")
            @printf(io, "GMG=%.6G,RG0=%.6G /\n", tidal.gmg, tidal.rg0)
            println(io)
        elseif kz14 == 5
            println(io, "&INXTRNL0")
            @printf(
                io,
                "RG=%.6G,%.6G,%.6G,VG=%.6G,%.6G,%.6G /\n",
                tidal.rg[1],
                tidal.rg[2],
                tidal.rg[3],
                tidal.vg[1],
                tidal.vg[2],
                tidal.vg[3]
            )
            println(io)
        end

        # --- 10. IMBHINIT (KZ(24) = 1): mass, position, velocity, DTBH ---
        if imbh_line !== nothing
            @printf(
                io,
                "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
                imbh_line.bimbh,
                imbh_line.ximbh[1],
                imbh_line.ximbh[2],
                imbh_line.ximbh[3],
                imbh_line.vimbh[1],
                imbh_line.vimbh[2],
                imbh_line.vimbh[3],
                imbh_line.dtbh
            )
        end
        # --- 10'. Flag route (KZ(24) = 3): NIMBHS DTBH, then the IMBH names ---
        if imbh_flag !== nothing
            @printf(io, "%d %.15e\n", length(imbh_flag.names), imbh_flag.dtbh)
            println(io, join(string.(imbh_flag.names), " "))
        end
        # --- 11. &INMBH (KZ(44) = 1): static black hole at the origin ---
        if mbh !== nothing
            @printf(
                io,
                "&INMBH MBH=%.10E,EPSPC=%.10E,RTMODE=%d,RACCRPC=%.10E,FACC=0.5,BETAF=%.6G,IPART=%d /\n",
                mbh.mass_msun,
                mbh.softening_pc,
                mbh.rtmode,
                mbh.r_acc_pc,
                mbh.beta_full,
                mbh.ipart
            )
        end

        # --- No binaries (NBIN0=0) or hierarchical triples ---
    end

    kz24_note = imbh !== nothing && kz[24] in (1, 3) ? ", KZ(24)=$(kz[24])" : ""
    @info "Wrote .inp file: $path (N=$N_total, NBIN0=$nbin0, KZ(22)=2, KZ(14)=$kz14$kz24_note, binary_physics=$(nbody6.binary_physics), RBAR=$rbar, ZMBAR=$zmbar)"
    if nbody6.checkpoint
        n_dumps = floor(Int, tcrit / dtadj) + 1
        dump_bytes = _checkpoint_dump_bytes(N_total, nbody6.nnbopt)
        @info "Checkpoint: a restart dump of at most $(Base.format_bytes(dump_bytes)) at every " *
              "adjustment, $(n_dumps) to TCRIT; at most $(Base.format_bytes(n_dumps * dump_bytes)) " *
              "if none is pruned"
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Internal: sample a single cluster from its spec
# ---------------------------------------------------------------------------
# `zmet` sets the ZAMS radii of the Roche guard (`binaries.roche_fill_max > 0`).
function _sample_cluster(
    spec::ClusterSpec;
    rng::AbstractRNG = Random.default_rng(),
    nmax::Int,
    zmet::Real = NaN,
)
    # Total-mass handling is delegated to the IMFSpec subtype (natural Kroupa
    # keeps the sampled sum; rescaled/equal enforce their targets).
    masses = sample_masses(spec.imf, spec.N, rng)

    # Primordial pairs: the density sampler places *systems* (pairs by their
    # centre of mass); pairs are expanded into bodies after the clusters are
    # combined and truncated. Systems are ordered pairs first, then singles.
    bins = sample_binaries(spec.binaries, masses, rng; zmet = zmet)
    n_b = length(bins.primary)
    paired = falses(spec.N)
    paired[bins.primary] .= true
    paired[bins.secondary] .= true
    singles = findall(!, paired)
    sys_mass = vcat(bins.m1 .+ bins.m2, bins.masses[singles])
    system_binary = vcat(collect(1:n_b), zeros(Int, length(singles)))
    N_sys = length(sys_mass)

    pos, vel = if spec.segregation.method == "bdk08"
        # Primordial mass segregation of the systems on the same model and scaling
        seg_pos, seg_vel, _ = segregate_systems(
            spec.profile,
            sys_mass,
            spec.rbar;
            degree = spec.segregation.degree,
            rng = rng,
            max_points = spec.segregation.max_points,
        )
        (seg_pos, seg_vel)
    elseif spec.profile isa PlummerProfile
        # Scale radius from the target half-mass radius (r_hm = 1.305 a)
        sample_plummer(N_sys, spec.rbar / _PLUMMER_RHM_OVER_A; rng = rng)
    elseif spec.profile isa KingProfile
        # Sampled with unit tidal radius; the empirical rescale below sets r_hm
        sample_king(N_sys, spec.profile.W0, 1.0; rng = rng)
    elseif spec.profile isa RotatingKingProfile
        # Model units (r_K = σ = 1), rotation about ẑ carried onto the profile axis
        model = _rotating_king_cached(spec.profile.W0, spec.profile.omega0)
        pos_m, vel_m = sample_rotating_king(N_sys, model; rng = rng)
        _rotate_to_axis!(pos_m, vel_m, spec.profile.axis)
        (pos_m, vel_m)
    else
        error("Unknown profile: $(typeof(spec.profile))")
    end

    # Scale positions so the mass-based half-mass radius equals the target
    r_hm = half_mass_radius(sys_mass, pos; centre = zeros(3))
    if r_hm > 0
        pos .*= spec.rbar / r_hm
    end

    energies = virialise!(sys_mass, pos, vel; nmax = nmax)
    binaries =
        (system_binary = system_binary, m1 = bins.m1, m2 = bins.m2, a_pc = bins.a_pc, e = bins.e)
    # Model and realised rotation of a rotating King cluster (pc, code velocity unit)
    rotation = if spec.profile isa RotatingKingProfile
        (
            moments = rotating_king_moments(
                _rotating_king_cached(spec.profile.W0, spec.profile.omega0),
            ),
            realised = rotation_moments(sys_mass, pos, vel; axis = spec.profile.axis),
        )
    else
        nothing
    end
    return (
        pos = pos,
        vel = vel,
        mass = sys_mass,
        T = energies.T,
        W = energies.W,
        binaries = binaries,
        rotation = rotation,
    )
end

"""
    _cluster_rng(seed::Int, i::Int) -> MersenneTwister

Random stream of cluster `i` under `seed_streams = "per_cluster"`: a
`MersenneTwister` seeded with the two 32-bit halves of `seed` and the cluster
index, so each cluster's stream depends on the seed and its index only.
"""
_cluster_rng(seed::Int, i::Int) = Random.MersenneTwister(
    UInt32[
        UInt32(reinterpret(UInt64, seed) & 0xffffffff),
        UInt32(reinterpret(UInt64, seed) >> 32),
        UInt32(i),
    ],
)

"""Mass-weighted rms speed [input velocity unit] of `vel` (3×N) about the mass-weighted mean velocity."""
function _rms_speed_about_mean(mass::AbstractVector{Float64}, vel::AbstractMatrix{Float64})::Float64
    M = sum(mass)
    v̄ = zeros(3)
    for i in eachindex(mass), k in 1:3
        v̄[k] += mass[i] * vel[k, i]
    end
    v̄ ./= M
    s = 0.0
    for i in eachindex(mass)
        s += mass[i] * ((vel[1, i] - v̄[1])^2 + (vel[2, i] - v̄[2])^2 + (vel[3, i] - v̄[3])^2)
    end
    return sqrt(s / M)
end

"""Largest distance of the systems `idx` of `pos` (3×N) from their mass-weighted centre; `NaN` for an empty set."""
function _outer_radius(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    idx::AbstractVector{Int},
)::Float64
    isempty(idx) && return NaN
    M = 0.0
    c = zeros(3)
    for j in idx
        M += mass[j]
        for k in 1:3
            c[k] += mass[j] * pos[k, j]
        end
    end
    c ./= M
    r2_max = 0.0
    for j in idx
        r2 = (pos[1, j] - c[1])^2 + (pos[2, j] - c[2])^2 + (pos[3, j] - c[3])^2
        r2_max = max(r2_max, r2)
    end
    return sqrt(r2_max)
end

"""
    _two_cluster_derived(cfg, cluster_data, pos, mass, system_ranges, d_apo, e_total) -> NamedTuple

[`two_cluster_orbit_quantities`](@ref) of a placed Kepler-mode pair: masses
and outer radii of the kept systems (`system_ranges` into the combined
system arrays `pos` [pc], `mass` [M☉], before pair expansion), σ₁ from the
sampled systems of cluster 1 (code-unit velocities converted to km s⁻¹) and
the configuration energy `e_total` [M☉ km² s⁻²].
"""
function _two_cluster_derived(
    cfg::MergerConfig,
    cluster_data::AbstractVector,
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    system_ranges::AbstractVector{<:AbstractVector{Int}},
    d_apo::Float64,
    e_total::Float64,
)::NamedTuple
    M1 = sum(mass[system_ranges[1]])
    M2 = sum(mass[system_ranges[2]])
    c1 = cluster_data[1]
    σ1 = _rms_speed_about_mean(c1.mass, c1.vel) * _CODE_VSTAR_KMS
    r_out = (_outer_radius(pos, mass, system_ranges[1]), _outer_radius(pos, mass, system_ranges[2]))
    return two_cluster_orbit_quantities(
        M1,
        M2,
        d_apo,
        cfg.orbit.eccentricity,
        cfg.clusters[1].rbar;
        sigma1 = σ1,
        e_total = e_total,
        r_out = r_out,
    )
end

"""Sampled minus kept system mass [M☉] of every cluster (`0` in single mode, which does not truncate)."""
function _truncated_masses(
    orbit_mode::AbstractString,
    cluster_data::AbstractVector,
    mass::AbstractVector{Float64},
    system_ranges::AbstractVector{<:AbstractVector{Int}},
)::Vector{Float64}
    orbit_mode == "single" && return zeros(length(cluster_data))
    return Float64[
        sum(cluster_data[i].mass) - sum(mass[system_ranges[i]]) for i in eachindex(cluster_data)
    ]
end

"""
    _complex_members(spec::ComplexSpec, template::ClusterSpec, rng::AbstractRNG) -> NamedTuple

Members of a cluster complex before sampling: masses `M` [M☉]
([`complex_member_masses`](@ref)), star counts `Nᵢ = round(Mᵢ / ⟨m⟩)` with the
mean stellar mass of the template IMF, half-mass radii
`r_h,i = complex_member_radius(spec, Mᵢ) 10^(size_scatter_dex gᵢ)` [pc] (`gᵢ`
standard normal, one draw per member after the masses), centres and raw
velocities ([`complex_centres`](@ref)), and the member `ClusterSpec`s (the
template with `N = Nᵢ`, `rbar = r_h,i`, `position` the centre, a zero
`velocity` placeholder). With `spec.imf_partition = "global"` one star list is
drawn and partitioned after the masses (`_complex_partition_stars`): member i
carries its share as an `ExplicitMassesIMF`, `Nᵢ` is the share's
count, and `Mᵢ` (returned as `masses` and used in the size–mass relation) its
realised sum. `target_masses` are the drawn masses in member order. Every
draw comes from `rng`. Throws an `ArgumentError` for fewer
than two members or a member below `spec.member_n_min` stars.
"""
function _complex_members(spec::ComplexSpec, template::ClusterSpec, rng::AbstractRNG)
    M_drawn = complex_member_masses(spec, rng)
    n = length(M_drawn)
    n ≥ 2 || throw(
        ArgumentError(
            "merger.complex: the mass function gave $n member(s); a complex needs at least 2",
        ),
    )
    m̄ = expected_mass(template.imf, 1)
    global_partition = spec.imf_partition == "global"
    stars = if global_partition
        _complex_partition_stars(template.imf, M_drawn, rng)
    else
        Vector{Float64}[]
    end
    N = if global_partition
        Int[length(s) for s in stars]
    else
        Int[round(Int, Mi / m̄) for Mi in M_drawn]
    end
    M::Vector{Float64} = global_partition ? Float64[sum(s) for s in stars] : M_drawn
    for i in 1:n
        N[i] ≥ spec.member_n_min || throw(
            ArgumentError(
                "merger.complex: member $i of mass $(round(M[i]; sigdigits = 4)) M☉ has " *
                "N = $(N[i]) stars (mean stellar mass $(round(m̄; sigdigits = 4)) M☉), below " *
                "member_n_min = $(spec.member_n_min)",
            ),
        )
    end
    r_h = Float64[
        complex_member_radius(spec, M[i]) * 10.0^(spec.size_scatter_dex * randn(rng)) for i in 1:n
    ]
    x, v = complex_centres(spec, n, rng)
    specs = ClusterSpec[
        ClusterSpec(;
            N = N[i],
            rbar = r_h[i],
            profile = template.profile,
            imf = global_partition ? ExplicitMassesIMF(stars[i]) : template.imf,
            position = x[:, i],
            velocity = zeros(3),
            binaries = template.binaries,
            segregation = template.segregation,
        ) for i in 1:n
    ]
    return (masses = M, target_masses = M_drawn, r_h = r_h, pos = x, vel = v, specs = specs)
end

"""
    _scale_complex_members(cs::ComplexSpec, members::NamedTuple, cluster_data) -> NamedTuple

Bulk velocities of the sampled members of a complex: the centres `members.pos`
and raw velocities `members.vel` (modified in place) are scaled by
[`complex_scale_velocities!`](@ref) with the realised member masses, softening
`cs.softening` or the mean member half-mass radius. Returns the member specs
with their centres and velocities [km s⁻¹] (`specs`), the scaling result
(`scaling`), the member masses (`masses`, M☉) and the softening (`softening`, pc).
"""
function _scale_complex_members(
    cs::ComplexSpec,
    members::NamedTuple,
    cluster_data::AbstractVector,
)::NamedTuple
    x_c, v_c = members.pos, members.vel
    m_c = Float64[sum(cd.mass) for cd in cluster_data]
    ε = cs.softening > 0 ? cs.softening : sum(members.r_h) / length(members.r_h)
    scaling = complex_scale_velocities!(
        v_c,
        x_c,
        m_c;
        virial_ratio = cs.virial_ratio,
        rotation_fraction = cs.rotation_fraction,
        softening = ε,
    )
    specs = ClusterSpec[
        ClusterSpec(;
            N = s.N,
            rbar = s.rbar,
            profile = s.profile,
            imf = s.imf,
            position = x_c[:, i],
            velocity = v_c[:, i],
            binaries = s.binaries,
            segregation = s.segregation,
        ) for (i, s) in enumerate(members.specs)
    ]
    return (specs = specs, scaling = scaling, masses = m_c, softening = ε)
end

"""
    _complex_energy(cfg::MergerConfig, members::NamedTuple, rngs) -> NamedTuple

System energies of the complex `members` as `generate_merger_ic` builds it
without a control: member i sampled from `rngs[i]`, the bulk velocities scaled
(`_scale_complex_members`) and the members combined as in explicit mode, with
the primordial pairs as point masses (before pair expansion), the footing on
which the control's own energies are taken. The internal energy of the pairs
is not matched (the control draws its own pairs). Returns `T`, `W` (code
units, G = 1) and the star count `n_stars` (systems plus pairs). Nothing is
written. Throws an `ArgumentError` when the star count exceeds
`cfg.virial_max_n`.
"""
function _complex_energy(
    cfg::MergerConfig,
    members::NamedTuple,
    rngs::AbstractVector{<:AbstractRNG},
)::NamedTuple
    cluster_data = [
        _sample_cluster(s; rng = rngs[i], nmax = cfg.virial_max_n, zmet = cfg.stellar.zmet) for
        (i, s) in enumerate(members.specs)
    ]
    specs = _scale_complex_members(cfg.complex, members, cluster_data).specs
    r_trunc_rh = cfg.output.truncate_radius_rh
    pos, vel, mass, _, kept = combine_clusters_explicit(
        cluster_data,
        specs;
        truncate_jacobi_flag = cfg.output.truncate_jacobi,
        truncation_radii = r_trunc_rh > 0 ? Float64[r_trunc_rh * c.rbar for c in specs] : Float64[],
        truncate_rule = cfg.output.truncate_rule,
    )
    # Every kept pair adds one star to its system
    n_pairs = sum(
        count(k -> cd.binaries.system_binary[k] > 0, kept[i]) for
        (i, cd) in enumerate(cluster_data)
    )
    n_stars = length(mass) + n_pairs
    n_stars ≤ cfg.virial_max_n || throw(
        ArgumentError(
            "merger.complex.control: the complex has N = $n_stars stars, above " *
            "merger.virial_max_n = $(cfg.virial_max_n); its energy is an exact O(N²) pair sum",
        ),
    )
    T, W = _kinetic_and_potential(mass, pos, vel)
    return (T = T, W = W, n_stars = n_stars)
end

"""
    _rescale_about_centre(mass, pos, vel, s, k) -> (pos, vel)

Copies of `pos` and `vel` (3 × N) scaled about the mass-weighted centre of
mass: positions by `s`, velocities by `k`.
"""
function _rescale_about_centre(
    mass::AbstractVector{Float64},
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    s::Float64,
    k::Float64,
)::Tuple{Matrix{Float64},Matrix{Float64}}
    M = sum(mass)
    c = zeros(3)
    u = zeros(3)
    for i in eachindex(mass), d in 1:3
        c[d] += mass[i] * pos[d, i]
        u[d] += mass[i] * vel[d, i]
    end
    c ./= M
    u ./= M
    pos_s = Matrix{Float64}(undef, 3, length(mass))
    vel_s = Matrix{Float64}(undef, 3, length(mass))
    for i in eachindex(mass), d in 1:3
        pos_s[d, i] = c[d] + s * (pos[d, i] - c[d])
        vel_s[d, i] = u[d] + k * (vel[d, i] - u[d])
    end
    return pos_s, vel_s
end

"""
    _control_scales(kind, T_c, W_c, E_cx, Q_cx) -> (s, k, q_target)

Length and velocity factors of an energy-matched control (code units, G = 1)
from the energies `T_c`, `W_c` of the sampled Plummer cluster and the total
energy `E_cx < 0` and virial ratio `Q_cx` of the complex. `"monolithic"`
(C-E): `s = E_c / E_cx`, `k = 1/√s`, so `E = E_cx` at the cluster's Q.
`"cold"` (C-cold): `|W| = |E_cx| / (1 − Q_cx)`, `T = Q_cx |W|`,
`s = |W_c| / |W|`, `k = √(T / T_c)`, so `E = E_cx` and `Q = Q_cx`.
"""
function _control_scales(
    kind::String,
    T_c::Float64,
    W_c::Float64,
    E_cx::Float64,
    Q_cx::Float64,
)::NTuple{3,Float64}
    if kind == "monolithic"
        s = (T_c + W_c) / E_cx
        return (s, 1 / sqrt(s), 0.5)
    end
    W_target = abs(E_cx) / (1 - Q_cx)
    T_target = Q_cx * W_target
    return (abs(W_c) / W_target, sqrt(T_target / T_c), Q_cx)
end

"""
    generate_merger_ic(cfg::MergerConfig; rng=Random.default_rng(),
                       output_dir="") -> MergerICResult

Generate complete merger initial conditions from a [`MergerConfig`](@ref).
Produces `dat.10`, `merger.inp`, and `merger_summary.txt`.

Supports four orbit modes:
- `"kepler"`: 2 clusters on a Keplerian orbit (positions/velocities auto-computed)
- `"explicit"`: N clusters with user-specified CM positions and velocities
- `"single"`: one cluster at the origin, at rest (no truncation)
- `"complex"`: the members of a cluster complex (`cfg.complex`): masses,
  radii and centres drawn from the complex's own stream `_cluster_rng(seed, 0)`,
  member `i` sampled from `_cluster_rng(seed, i)`; after sampling, the centres
  and bulk velocities are scaled to the target virial ratio of the realised
  member masses as softened point masses (`complex_scale_velocities!`,
  softening `cfg.complex.softening` or the mean member half-mass radius), and
  the members are combined as in explicit mode

With `cfg.complex.control = "monolithic"` or `"cold"` the complex is built
first as without the control (same streams, nothing written) and the energies
of its systems are evaluated (pairs as point masses; their internal energy is
not matched, the control draws its own pairs); one Plummer cluster of the union
of the members' star lists, sampled and virialised from `_cluster_rng(seed, 1)`
at r_h = 1 pc, is then rescaled about its centre of mass to the complex's total
energy (C-E, Q = 0.5) or to its total energy and virial ratio (C-cold), and is
written instead as a single cluster at the origin, at rest.

With `cfg.seed_streams = "per_cluster"` (always in complex mode) and no
external `rng`, cluster `i` is sampled and its pairs expanded from its own
stream (`_cluster_rng(seed, i)`); an external `rng` replaces every stream.

Every cluster with an [`IMBHSpec`](@ref) gets its IMBH, in cluster order,
on its systems after truncation and before pair expansion: centred on their
centre of mass (plus the offset), its velocity from the cluster's stream for
`velocity = "equipartition"`, and with `response = "jeans"` the systems'
velocities rescaled. The IMBHs are the last bodies of `dat.10`, in cluster
order (names `N_stars + 1, …`). The route follows
[`resolved_imbh_route`](@ref): under `:kz24` (one IMBH) its state is repeated
on the `IMBHINIT` line of `merger.inp` (`KZ(24) = 1`); under `:flag`
(`KZ(24) = 3`, up to 16) the lines `NIMBHS DTBH` and the IMBH names follow
`&INSCALE`. `N_total` and `M_total` include the IMBHs, the cluster ranges and
`RBAR` (the stars' half-mass radius) do not. Refused: M•/M_stars > 0.05 and
M•/⟨m⟩ < 20 (warned above 0.02, below 50 and above 10⁴).

Returns a [`MergerICResult`](@ref) for downstream plotting and inspection.
"""
function generate_merger_ic(
    cfg::MergerConfig;
    rng::Union{AbstractRNG,Nothing} = nothing,
    output_dir::AbstractString = "",
)
    out_dir = isempty(output_dir) ? cfg.output.output_dir : output_dir
    isempty(out_dir) && (out_dir = pwd())
    mkpath(out_dir)

    n_clusters = length(cfg.clusters)
    n_clusters ≥ 1 || error("Need at least 1 cluster, got $n_clusters")
    (cfg.orbit_mode == "kepler" && n_clusters != 2) &&
        error("Kepler orbit mode requires exactly 2 clusters, got $n_clusters")
    complex_mode = cfg.orbit_mode == "complex"
    (complex_mode && n_clusters != 1) &&
        error("Complex mode takes exactly one member template in cfg.clusters, got $n_clusters")
    complex_mode && _validate_complex(cfg.complex)
    for (i, spec) in enumerate(cfg.clusters)
        _validate_cluster_spec(spec, i)
    end
    _validate_nbody6(cfg.nbody6)
    _validate_stellar(cfg.stellar, cfg.output)
    _validate_binary_physics(cfg.nbody6, cfg.stellar)
    _validate_tidal(cfg.tidal)
    _validate_tidal_tolerance(cfg.tidal, cfg.nbody6)
    # IMBH route checks of the loader, for configurations built in code
    imbh_carriers = findall(c -> c.imbh !== nothing, cfg.clusters)
    imbh_violation = _imbh_route_violation(cfg.nbody6.imbh_route, imbh_carriers)
    imbh_violation === nothing ||
        error("generate_merger_ic: $(imbh_violation[1]): $(imbh_violation[2])")
    imbh_route = resolved_imbh_route(cfg)
    (!isempty(imbh_carriers) && haskey(cfg.nbody6.kz, 24)) && error(
        "generate_merger_ic: merger.nbody6.kz KZ(24) is set by the IMBH route; remove the override",
    )
    cfg.seed_streams in _SEED_STREAMS || throw(
        ArgumentError(
            "seed_streams must be \"shared\" or \"per_cluster\"; got \"$(cfg.seed_streams)\"",
        ),
    )

    # Resolve RNG and effective seed. When the caller supplies an RNG the
    # sampling is NOT reproducible from `effective_seed`; the metadata records
    # this honestly via `external_rng` (the seed still feeds Nbody6's NRAND).
    external_rng = rng !== nothing
    effective_seed::Int = cfg.seed === nothing ? Int(rand(UInt32) % typemax(Int32)) : cfg.seed
    use_rng = external_rng ? rng : Random.MersenneTwister(effective_seed)
    # Complex mode: member masses, radii and centres from the complex's own
    # stream (index 0); the members then take the place of cfg.clusters.
    members = if complex_mode
        rng0::AbstractRNG = rng === nothing ? _cluster_rng(effective_seed, 0) : rng
        _complex_members(cfg.complex, cfg.clusters[1], rng0)
    else
        nothing
    end
    # Controls (C-E, C-cold): the all-star energy of the complex as built
    # without the control, then one cluster of the union star list instead.
    control = complex_mode ? cfg.complex.control : "none"
    control_spec = nothing
    control_target = nothing
    if control != "none"
        cm = members::NamedTuple
        n_m = length(cm.specs)
        member_rngs = if external_rng
            AbstractRNG[rng::AbstractRNG for _ in 1:n_m]
        else
            AbstractRNG[_cluster_rng(effective_seed, i) for i in 1:n_m]
        end
        cx = _complex_energy(cfg, cm, member_rngs)
        E_cx = cx.T + cx.W
        (cx.W < 0 && E_cx < 0) || throw(
            ArgumentError(
                "merger.complex.control: the complex is unbound (E = $E_cx, W = $(cx.W) in " *
                "code units); an energy-matched control needs E < 0",
            ),
        )
        stars = Float64[m for s in cm.specs for m in (s.imf::ExplicitMassesIMF).masses]
        control_spec = ClusterSpec(;
            N = length(stars),
            rbar = 1.0,
            profile = PlummerProfile(),
            imf = ExplicitMassesIMF(stars),
            binaries = cfg.clusters[1].binaries,
        )
        control_target = (e = E_cx, q = cx.T / abs(cx.W))
        members = nothing
    end
    specs::Vector{ClusterSpec} = if control_spec !== nothing
        ClusterSpec[control_spec]
    elseif members === nothing
        cfg.clusters
    else
        members.specs
    end
    n_clusters = length(specs)
    # Per-cluster streams: cluster i is sampled and its pairs expanded from
    # its own stream, independent of the other clusters and of the orbit.
    per_cluster = (cfg.seed_streams == "per_cluster" || complex_mode) && !external_rng
    (external_rng && cfg.seed_streams == "per_cluster") &&
        @warn "seed_streams = \"per_cluster\" needs the configured seed; " *
              "the external RNG is used for every cluster"
    cluster_rngs = if per_cluster
        [_cluster_rng(effective_seed, i) for i in 1:n_clusters]
    else
        Random.MersenneTwister[]
    end

    seed_note = external_rng ? "external RNG (seed not reproducible)" : "seed=$effective_seed"
    @info "Generating merger ICs for $n_clusters clusters ($(cfg.orbit_mode) mode, $seed_note)..."

    # Sample each cluster independently
    cluster_data = [
        begin
            @info "  Cluster $i: $(profile_name(spec.profile)) profile, N=$(spec.N), " *
                  "imf=$(imf_name(spec.imf)), M≈$(round(expected_mass(spec.imf, spec.N); digits=1)) M☉"
            _sample_cluster(
                spec;
                rng = per_cluster ? cluster_rngs[i] : use_rng,
                nmax = cfg.virial_max_n,
                zmet = cfg.stellar.zmet,
            )
        end for (i, spec) in enumerate(specs)
    ]

    # Complex mode: centres and bulk velocities [km s⁻¹] scaled to the target
    # virial ratio of the realised member masses as softened point masses.
    complex_info = nothing
    if members !== nothing
        scaled = _scale_complex_members(cfg.complex, members, cluster_data)
        specs = scaled.specs
        scaling = scaled.scaling
        m_c = scaled.masses
        ε = scaled.softening
        M_c = sum(m_c)
        λ_c =
            sqrt(sum(abs2, scaling.l_total)) * sqrt(abs(scaling.e_total)) /
            (_G_PC_KMS2_MSUN * M_c^2.5)
        complex_info = (
            scaling = scaling,
            masses = m_c,
            target_masses = members.target_masses,
            r_h = members.r_h,
            lambda_peebles = λ_c,
        )
        @info @sprintf(
            "Complex: %d members, M = %.4g M☉, Q_cc = %.3f, λ = %.3g, ε = %.3g pc",
            length(m_c),
            M_c,
            scaling.virial_ratio,
            λ_c,
            ε
        )
    end

    # Control: the virialised Plummer cluster rescaled about its centre of
    # mass to the energy (and, for C-cold, the virial ratio) of the complex.
    control_info = nothing
    if control_target !== nothing
        sampled = cluster_data[1]
        T_c, W_c = _kinetic_and_potential(sampled.mass, sampled.pos, sampled.vel)
        s_len, k_vel, q_target =
            _control_scales(control, T_c, W_c, control_target.e, control_target.q)
        pos_s, vel_s = _rescale_about_centre(sampled.mass, sampled.pos, sampled.vel, s_len, k_vel)
        scaled_data = (
            pos = pos_s,
            vel = vel_s,
            mass = sampled.mass,
            T = k_vel^2 * T_c,
            W = W_c / s_len,
            binaries = sampled.binaries,
        )
        cluster_data = [scaled_data]
        sampled_spec = specs[1]
        control_cluster = ClusterSpec(;
            N = sampled_spec.N,
            rbar = s_len * sampled_spec.rbar,
            profile = sampled_spec.profile,
            imf = sampled_spec.imf,
            binaries = sampled_spec.binaries,
        )
        specs = ClusterSpec[control_cluster]
        control_info = (kind = control, e_target = control_target.e, q_target = q_target)
        @info @sprintf(
            "Control %s: E = %.6g, Q = %.4f matched (code units), r_h = %.3g pc",
            control,
            control_target.e,
            q_target,
            specs[1].rbar
        )
    end

    # Kepler apocentre from the configured size key; fixed truncation radii
    # (0 = Jacobi radius) in units of each cluster's input half-mass radius.
    d_apo = cfg.orbit_mode == "kepler" ? orbit_apocentre(cfg.orbit, cfg.clusters[1].rbar) : NaN
    r_trunc_rh = cfg.output.truncate_radius_rh

    # Combine clusters (as systems) according to orbit mode; both paths
    # return the combined arrays, the per-cluster index ranges of the systems
    # that survived truncation, and their indices in each sampled set.
    pos_combined, vel_combined, mass_combined, system_ranges, kept = if cfg.orbit_mode == "kepler"
        c1, c2 = cluster_data[1], cluster_data[2]
        r_trunc = (r_trunc_rh * cfg.clusters[1].rbar, r_trunc_rh * cfg.clusters[2].rbar)
        setup_two_cluster_orbit(
            c1.pos,
            c1.vel,
            c1.mass,
            c2.pos,
            c2.vel,
            c2.mass,
            d_apo,
            cfg.orbit.eccentricity;
            truncate_jacobi_flag = cfg.output.truncate_jacobi,
            truncation_radii = r_trunc,
            truncate_rule = cfg.output.truncate_rule,
            orbit_normal = cfg.orbit.orbit_normal,
            apocentre_angle = cfg.orbit.apocentre_angle,
        )
    elseif control_info !== nothing
        # The control is written as a single cluster: origin, at rest
        combine_single_cluster(cluster_data[1])
    elseif cfg.orbit_mode == "explicit" || complex_mode
        combine_clusters_explicit(
            cluster_data,
            specs;
            truncate_jacobi_flag = cfg.output.truncate_jacobi,
            truncation_radii = r_trunc_rh > 0 ? Float64[r_trunc_rh * c.rbar for c in specs] :
                               Float64[],
            truncate_rule = cfg.output.truncate_rule,
        )
    elseif cfg.orbit_mode == "single"
        combine_single_cluster(cluster_data[1])
    else
        error("Unknown orbit_mode: $(cfg.orbit_mode)")
    end

    M_total = sum(mass_combined)
    N_systems = length(mass_combined)

    # Half-mass radius of the combined system — this is the RBAR length unit
    # written to the .inp file and used for the NB-unit conversion of dat.10.
    rbar = half_mass_radius(mass_combined, pos_combined)

    # Combined-system diagnostics in code units (G = 1) on the systems (pair
    # binding energies excluded): the virial ratio of the whole
    # configuration, the resolution of the members by the length unit, and
    # the integration parameters scaled to the smallest member.
    q_virial, t_cr_config_code, e_config_code = if N_systems ≤ cfg.virial_max_n
        T, W = _kinetic_and_potential(mass_combined, pos_combined, vel_combined)
        (W < 0 ? T / abs(W) : NaN, crossing_time(M_total, T + W), T + W)
    else
        @warn "Combined virial ratio not evaluated: N = $N_systems systems exceed merger.virial_max_n = $(cfg.virial_max_n) (O(N²) pair sum)"
        (NaN, NaN, NaN)
    end

    # Mass removed by truncation and, in Kepler mode, the derived orbit
    # quantities, from the system arrays before pair expansion. The code
    # energy unit G M☉²/pc is _G_PC_KMS2_MSUN M☉ km² s⁻².
    mass_truncated = _truncated_masses(cfg.orbit_mode, cluster_data, mass_combined, system_ranges)
    mass_kept = Float64[sum(mass_combined[r]) for r in system_ranges]
    # The radius each cluster was truncated at [pc]; 0 when it was not truncated
    truncation_radii_pc = if !cfg.output.truncate_jacobi || control_info !== nothing
        zeros(n_clusters)
    elseif cfg.orbit_mode == "kepler"
        collect(
            _pair_truncation_radii(
                d_apo,
                sum(cluster_data[1].mass),
                sum(cluster_data[2].mass),
                cfg.output.truncate_rule,
                (r_trunc_rh * cfg.clusters[1].rbar, r_trunc_rh * cfg.clusters[2].rbar),
            ),
        )
    elseif cfg.orbit_mode == "explicit" || complex_mode
        _explicit_truncation_radii(
            cluster_data,
            specs;
            truncate_jacobi_flag = true,
            truncate_rule = cfg.output.truncate_rule,
            truncation_radii = r_trunc_rh > 0 ? Float64[r_trunc_rh * c.rbar for c in specs] :
                               Float64[],
        )
    else
        zeros(n_clusters)
    end
    orbit_derived = if cfg.orbit_mode == "kepler"
        _two_cluster_derived(
            cfg,
            cluster_data,
            pos_combined,
            mass_combined,
            system_ranges,
            d_apo,
            e_config_code * _G_PC_KMS2_MSUN,
        )
    else
        nothing
    end
    # A bound pair in a rotating-frame tide must stay inside its pair Jacobi radius.
    if orbit_derived !== nothing && cfg.tidal.kz14 in (1, 2)
        r_j = pair_jacobi_radius(sum(mass_combined), cfg.tidal)
        ratio = d_apo / r_j
        ratio < 1 || error(
            "merger.orbit: apocentre $(round(d_apo; sigdigits = 3)) pc reaches the pair Jacobi " *
            "radius $(round(r_j; sigdigits = 3)) pc of the tide (kz14 = $(cfg.tidal.kz14)): the pair " *
            "is not bound to itself; choose a smaller orbit",
        )
        ratio > 0.5 && @warn "Apocentre is $(round(ratio; digits = 2)) of the pair Jacobi radius " *
              "$(round(r_j; sigdigits = 3)) pc: the tide strips the pair orbit"
        orbit_derived = merge(orbit_derived, (apocentre_over_jacobi_pair = ratio,))
    end

    # Realised energy of the control on the systems, the footing of its target
    if control_info !== nothing
        T_r, W_r = _kinetic_and_potential(mass_combined, pos_combined, vel_combined)
        control_info = merge(control_info, (e_realised = T_r + W_r, q_realised = T_r / abs(W_r)))
    end

    # IMBHs on the systems of their clusters, in cluster order, before pair
    # expansion: a primordial pair counts as one system in ⟨m⟩ and in the dispersion.
    imbh_states = NamedTuple[]
    for k_imbh in findall(c -> c.imbh !== nothing, specs)
        spec_bh = specs[k_imbh].imbh::IMBHSpec
        members_bh = collect(Int, system_ranges[k_imbh])
        M_i = sum(mass_combined[members_bh])
        m_mean = M_i / length(members_bh)
        M_bh = isfinite(spec_bh.mass) ? spec_bh.mass : spec_bh.mass_fraction * M_i
        frac_bh = M_bh / M_i
        μ_bh = M_bh / m_mean
        frac_bh > 0.05 && error(
            "IMBH of cluster $k_imbh: M•/M_stars = $(round(frac_bh; sigdigits = 3)) exceeds 0.05 " *
            "(BME04a: the cusp extends to M(<r) = M• above 5 %)",
        )
        frac_bh > 0.02 &&
            @warn "IMBH of cluster $k_imbh: M•/M_stars = $(round(frac_bh; sigdigits = 3)) exceeds 0.02"
        μ_bh < 20 && error(
            "IMBH of cluster $k_imbh: not a massive particle (M•/⟨m⟩ = $(round(μ_bh; sigdigits = 3)) " *
            "< 20); raise the mass",
        )
        μ_bh < 50 &&
            @warn "IMBH of cluster $k_imbh: M•/⟨m⟩ = $(round(μ_bh; sigdigits = 3)) lies below 50, " *
                  "outside the calibrated range of the Brownian tests"
        μ_bh > 1e4 &&
            @warn "IMBH of cluster $k_imbh: M•/⟨m⟩ = $(round(μ_bh; sigdigits = 3)); the wandering " *
                  "law and the FP32 regular force are untested above μ ≈ 10⁴"
        x_c, v_c = _imbh_centre(pos_combined, vel_combined, mass_combined, members_bh)
        σ_c = _imbh_central_dispersion(
            pos_combined,
            vel_combined,
            mass_combined,
            members_bh,
            x_c,
            v_c,
        )
        x_bh, v_bh = _imbh_body(
            spec_bh,
            M_bh,
            m_mean,
            x_c,
            v_c,
            σ_c,
            per_cluster ? cluster_rngs[k_imbh] : use_rng,
        )
        jeans_factors = if spec_bh.response == "jeans"
            _imbh_jeans_response!(
                pos_combined,
                vel_combined,
                mass_combined,
                members_bh,
                v_c,
                M_bh,
                x_bh,
            )
        else
            Float64[]
        end
        push!(
            imbh_states,
            (
                cluster = k_imbh,
                spec = spec_bh,
                mass = M_bh,
                mass_fraction = frac_bh,
                m_mean = m_mean,
                σ_c = σ_c,
                pos = x_bh,
                vel = v_bh,
                jeans_factors = jeans_factors,
            ),
        )
        @info @sprintf(
            "IMBH: cluster %d, M• = %.4g M☉ (M•/M_stars = %.3g, M•/⟨m⟩ = %.3g), σ_c = %.3g km/s, velocity %s, response %s",
            k_imbh,
            M_bh,
            frac_bh,
            μ_bh,
            σ_c * _CODE_VSTAR_KMS,
            spec_bh.velocity,
            spec_bh.response
        )
    end

    # Static black hole at the origin ([merger.mbh]): the stellar velocities
    # respond to its potential as to a seeded IMBH's (isotropic Jeans rescale)
    mbh_factors = Float64[]
    if cfg.mbh !== nothing
        members_all = collect(Int, 1:size(pos_combined, 2))
        M_sys = sum(mass_combined)
        v_mean = [sum(mass_combined .* vel_combined[k, :]) / M_sys for k in 1:3]
        if cfg.mbh.response == "jeans"
            mbh_factors = _imbh_jeans_response!(
                pos_combined,
                vel_combined,
                mass_combined,
                members_all,
                v_mean,
                cfg.mbh.mass_msun,
                zeros(3),
            )
        end
        @info @sprintf(
            "Static MBH: M = %.4g M☉ (M/M_stars = %.3g), r_t mode %s, partial %s, response %s",
            cfg.mbh.mass_msun,
            cfg.mbh.mass_msun / M_sys,
            cfg.mbh.rt_mode,
            cfg.mbh.partial,
            cfg.mbh.response
        )
    end

    # Expand primordial pairs into bodies: every cluster's pairs first, then
    # the singles (the engine's convention for the first 2·NBIN0 bodies).
    expanded = expand_binaries(
        pos_combined,
        vel_combined,
        mass_combined,
        system_ranges,
        [cd.binaries for cd in cluster_data],
        kept;
        rng = per_cluster ? cluster_rngs : use_rng,
    )
    pos_combined, vel_combined, mass_combined = expanded.pos, expanded.vel, expanded.mass
    # The IMBHs are the last bodies in cluster order (particle N of IMBHINIT
    # under kz24, the flagged names under flag), in no cluster's blocks; they
    # count in N_total and M_total, while RBAR stays the half-mass radius of
    # the stars.
    for st in imbh_states
        pos_combined = hcat(pos_combined, st.pos)
        vel_combined = hcat(vel_combined, st.vel)
        mass_combined = vcat(mass_combined, st.mass)
        M_total += st.mass
    end
    # Inertial velocities into the engine's rotating frame (kz14 = 1, 2): v − Ω × r,
    # with Ω converted to the code velocity unit per pc; no rotating frame otherwise.
    if cfg.tidal.velocity_frame == "inertial" && cfg.tidal.kz14 in (1, 2)
        to_rotating_frame!(
            pos_combined,
            vel_combined,
            tidal_frame_rotation(cfg.tidal) / _CODE_VSTAR_KMS,
        )
    end
    cluster_blocks = expanded.cluster_blocks
    cluster_ranges = [_members_from_blocks(b) for b in cluster_blocks]
    n_pairs = expanded.n_pairs
    nbin0 = sum(n_pairs)
    N_total = length(mass_combined)
    zmbar = M_total / N_total
    nbin0 > 0 && @info "Primordial binaries: NBIN0 = $nbin0 pairs (" *
          join(["cluster $i: $(n_pairs[i])" for i in 1:length(n_pairs)], ", ") *
          "); hard fractions " *
          join([isnan(f) ? "-" : @sprintf("%.2f", f) for f in expanded.hard_fraction], ", ")
    # Crossing times: configuration and smallest member (pre-truncation
    # energies of the virialised members), in code units → NB and Myr.
    t_star_myr = _nbody_time_myr(M_total, rbar)
    code_to_nb = _CODE_TIME_MYR / t_star_myr
    t_cr_member_code = minimum(crossing_time(sum(c.mass), c.T + c.W) for c in cluster_data)
    r_h_min_pc = minimum(c.rbar for c in specs)
    rbar_over_rhm_min = rbar / r_h_min_pc
    nbody6 = resolve_nbody6_parameters(cfg.nbody6, specs, cluster_ranges, N_total, rbar)
    # Bound on the per-pair relative error of the FP32 regular force from the
    # coordinate offset: 4 u X / RS0, u = 2⁻²⁴, X the largest |x| [NB]
    # (the device rounds the coordinates to FP32, so the offset of a pair
    # from the origin sets the error).
    x_max_nb =
        maximum(
            sqrt(pos_combined[1, i]^2 + pos_combined[2, i]^2 + pos_combined[3, i]^2) for
            i in axes(pos_combined, 2)
        ) / rbar
    fp32_offset_bound = 4 * 2.0^-24 * x_max_nb / nbody6.rs0
    # Time parameters: physical values are converted with the realised T*.
    out = cfg.output
    tcrit_nb = out.tcrit_myr > 0 ? out.tcrit_myr / t_star_myr : out.tcrit
    # The intervals pass through the engine's decimal digit counter: write
    # dyadic values with an exact decimal expansion (engine_interval).
    dtadj_nb =
        _engine_interval_logged("DTADJ", out.dtadj_myr > 0 ? out.dtadj_myr / t_star_myr : out.dtadj)
    deltat_nb = _engine_interval_logged(
        "DELTAT",
        out.deltat_myr > 0 ? out.deltat_myr / t_star_myr : out.deltat,
    )
    dtplot_nb = _engine_interval_logged(
        "DTPLOT",
        cfg.stellar.dtplot_myr > 0 ? cfg.stellar.dtplot_myr / t_star_myr : cfg.stellar.dtplot,
    )
    dtplot_nb ≥ deltat_nb || error(
        "merger.stellar.dtplot ($(round(dtplot_nb; sigdigits = 4)) NB) must be ≥ merger.output.deltat " *
        "($(round(deltat_nb; sigdigits = 4)) NB) after conversion with T* = $(round(t_star_myr; sigdigits = 4)) Myr",
    )
    stellar = StellarSpec(;
        kz19 = cfg.stellar.kz19,
        level = cfg.stellar.level,
        zmet = cfg.stellar.zmet,
        epoch0 = cfg.stellar.epoch0,
        dtplot = dtplot_nb,
        dtplot_myr = cfg.stellar.dtplot_myr,
    )
    _check_multicluster_regime(
        q_virial,
        rbar_over_rhm_min,
        nbody6,
        r_h_min_pc / rbar;
        deltat = deltat_nb,
        t_cr_member_min_nb = t_cr_member_code * code_to_nb,
        fp32_offset_bound = fp32_offset_bound,
    )
    regime = (
        q_virial = q_virial,
        rbar_over_rhm_min = rbar_over_rhm_min,
        nbody6 = nbody6,
        cluster_blocks = cluster_blocks,
        n_pairs = n_pairs,
        nbin0 = nbin0,
        tcrit_nb = tcrit_nb,
        dtadj_nb = dtadj_nb,
        deltat_nb = deltat_nb,
        dtplot_nb = dtplot_nb,
        hard_fraction = expanded.hard_fraction,
        t_cr_config_nb = t_cr_config_code * code_to_nb,
        t_cr_config_myr = t_cr_config_code * _CODE_TIME_MYR,
        t_cr_member_min_nb = t_cr_member_code * code_to_nb,
        t_cr_member_min_myr = t_cr_member_code * _CODE_TIME_MYR,
        t_star_myr = t_star_myr,
        fp32_offset_bound = fp32_offset_bound,
        cluster_rotation = Any[get(c, :rotation, nothing) for c in cluster_data],
    )
    @info @sprintf(
        "Combined system: Q = %.3f, RBAR/r_hm,min = %.2f, t_cr = %.3g NB (%.3g Myr), smallest member t_cr = %.3g NB; RS0 = %.3g, RMIN = %.3g, DTMIN = %.3g, NNBOPT = %d",
        q_virial,
        rbar_over_rhm_min,
        regime.t_cr_config_nb,
        regime.t_cr_config_myr,
        regime.t_cr_member_min_nb,
        nbody6.rs0,
        nbody6.rmin,
        nbody6.dtmin,
        nbody6.nnbopt
    )

    # Keep physical-unit copies before conversion (sampled and Kepler
    # velocities carry the (G = 1, M☉, pc) code unit, _CODE_VSTAR_KMS km/s).
    mass_phys = copy(mass_combined)
    pos_phys = copy(pos_combined)
    vel_phys = vel_combined .* _CODE_VSTAR_KMS  # code units → km/s

    # Convert to N-body units for dat.10 (KZ(22)=2, the only supported
    # output format). `to_nbody_units!` expects velocity in
    # km/s, so convert from code units first.
    cfg.output.format == "nbody" || error(
        "Unsupported output format \"$(cfg.output.format)\"; only \"nbody\" (KZ(22)=2) is supported.",
    )
    vel_combined .*= _CODE_VSTAR_KMS  # code units → km/s (match to_nbody_units! API)
    to_nbody_units!(mass_combined, pos_combined, vel_combined, M_total, rbar)

    # IMBH input of merger.inp (DTBH with the realised T*) and the run record:
    # the scalar fields of the first IMBH, every IMBH under `bodies`
    imbh_inp = nothing
    imbh_record = nothing
    n_imbh = length(imbh_states)
    if n_imbh > 0
        (imbh_route === :flag || (imbh_route === :kz24 && n_imbh == 1)) || error(
            "generate_merger_ic: $n_imbh IMBHs do not fit the resolved route $(repr(imbh_route))",
        )
        dtbh_nb = nbody6.imbh_output_dt_myr / t_star_myr
        n_stars = N_total - n_imbh
        imbh_names = collect((n_stars + 1):N_total)
        imbh_inp = if imbh_route === :kz24
            (
                route = :kz24,
                mass_nb = mass_combined[end],
                pos_nb = pos_combined[:, end],
                vel_nb = vel_combined[:, end],
                dtbh_nb = dtbh_nb,
            )
        else
            (route = :flag, names = imbh_names, dtbh_nb = dtbh_nb)
        end
        imbh_bodies = [
            Dict{String,Any}(
                "name" => imbh_names[j],
                "cluster" => st.cluster,
                "mass_msun" => st.mass,
                "mass_fraction" => st.mass_fraction,
                "mean_stellar_mass_msun" => st.m_mean,
                "sigma_c_kms" => st.σ_c * _CODE_VSTAR_KMS,
                "placement" => st.spec.placement,
                "velocity_mode" => st.spec.velocity,
                "response" => st.spec.response,
                "jeans_factors" => collect(Float64, st.jeans_factors),
                "spin" => st.spec.spin,
                "position_pc" => pos_phys[:, imbh_names[j]],
                "velocity_kms" => vel_phys[:, imbh_names[j]],
            ) for (j, st) in enumerate(imbh_states)
        ]
        imbh_record = merge(
            copy(imbh_bodies[1]),
            Dict{String,Any}(
                "route" => String(imbh_route),
                "names" => imbh_names,
                "clusters" => Int[st.cluster for st in imbh_states],
                "bodies" => imbh_bodies,
                "dtbh_nb" => dtbh_nb,
                "rbar_note" => "length unit from the stars alone",
            ),
        )
    end

    # Write files
    write_dat10(joinpath(out_dir, "dat.10"), mass_combined, pos_combined, vel_combined)
    generate_merger_inp(
        joinpath(out_dir, "merger.inp"),
        N_total,
        rbar,
        zmbar;
        nbody6 = nbody6,
        stellar = stellar,
        tidal = cfg.tidal,
        tcrit = tcrit_nb,
        dtadj = dtadj_nb,
        deltat = deltat_nb,
        nrand = effective_seed,
        mass_bounds = (
            minimum(_imf_mass_bounds(c.imf)[1] for c in specs),
            maximum(_imf_mass_bounds(c.imf)[2] for c in specs),
        ),
        nbin0 = nbin0,
        imbh = imbh_inp,
        mbh = cfg.mbh === nothing ? nothing :
              (
            mass_msun = cfg.mbh.mass_msun,
            softening_pc = cfg.mbh.softening_pc,
            rtmode = cfg.mbh.rt_mode == "stellar" ? 1 : 2,
            r_acc_pc = cfg.mbh.r_acc_pc,
            beta_full = cfg.mbh.beta_full,
            ipart = _MBH_IPART[cfg.mbh.partial],
        ),
    )

    # Summary log
    _write_merger_summary(
        joinpath(out_dir, "merger_summary.txt"),
        cfg,
        cluster_data,
        cluster_ranges,
        N_total,
        M_total,
        rbar,
        zmbar,
        regime;
        d_apo = d_apo,
        orbit_derived = orbit_derived,
        mass_truncated = mass_truncated,
        cluster_specs = specs,
        complex_info = complex_info,
        control_info = control_info,
        imbh_record = imbh_record,
    )

    # Structured metadata for post-hoc regeneration of IC plots
    _write_merger_ic_metadata(
        joinpath(out_dir, "merger_ic.toml"),
        cfg,
        cluster_ranges,
        effective_seed,
        external_rng,
        N_total,
        M_total,
        rbar,
        zmbar,
        regime;
        d_apo = d_apo,
        orbit_derived = orbit_derived,
        mass_truncated = mass_truncated,
        mass_kept = mass_kept,
        truncation_radii_pc = truncation_radii_pc,
        cluster_specs = specs,
        complex_info = complex_info,
        control_info = control_info,
        imbh_record = imbh_record,
        mbh_record = cfg.mbh === nothing ? nothing :
                     Dict{String,Any}(
            "mass_msun" => cfg.mbh.mass_msun,
            "softening_pc" => cfg.mbh.softening_pc,
            "rt_mode" => cfg.mbh.rt_mode,
            "r_acc_pc" => cfg.mbh.r_acc_pc,
            "partial" => cfg.mbh.partial,
            "beta_full" => cfg.mbh.beta_full,
            "response" => cfg.mbh.response,
            "jeans_factors" => mbh_factors,
        ),
    )

    return MergerICResult(
        out_dir,
        N_total,
        M_total,
        rbar,
        zmbar,
        cluster_ranges,
        n_pairs,
        collect(specs),
        cfg.orbit_mode,
        cfg.orbit,
        mass_phys,
        pos_phys,
        vel_phys,
    )
end

"""
    _write_merger_ic_metadata(path, cfg, cluster_ranges, seed, external_rng,
                              N_total, M_total, rbar, zmbar, regime;
                              d_apo = NaN, orbit_derived = nothing,
                              mass_truncated = Float64[], mass_kept = Float64[],
                              truncation_radii_pc = Float64[], cluster_specs = cfg.clusters,
                              complex_info = nothing, control_info = nothing,
                              imbh_record = nothing)

Write a machine-readable TOML snapshot of the IC generation (schema version 4),
sufficient to reconstruct a [`MergerICResult`](@ref) from the `dat.10` file
later (so `plot_merger_ic` can be re-run after code fixes without
re-sampling). Cluster specs are stored in the structured form
(`profile = {type=...}`, `imf = {type=...}`) that
[`_parse_cluster_table`](@ref) also accepts, so the file round-trips.

The `[orbit]` table holds every [`OrbitSpec`](@ref) field as configured; in
Kepler mode also `apocentre_resolved` (`d_apo` [pc]) and every field of
`orbit_derived` ([`two_cluster_orbit_quantities`](@ref), pairs as 2-element
arrays). `mass_truncated[i]` [M☉], the sampled minus the kept mass,
`mass_kept[i]` (as `mass`) [M☉], the kept mass of the cluster's stars, and
`truncation_radii_pc[i]` (as `truncation_radius_pc`) [pc], the radius the
cluster was truncated at (0 when it was not), go to the table of cluster `i`
when given. The `[output]` table records `truncate_rule`. A cluster of a
[`RotatingKingProfile`](@ref) gets a `rotation` subtable with `omega0`, `axis`,
the model's `lambda_model`, `t_rot_model` and `r_h_over_r_c_model`
([`rotating_king_moments`](@ref)) and the realised `lambda_realised` and
`l_axis_realised_pc_kms` [pc km s⁻¹] of its virialised systems
([`rotation_moments`](@ref)), read from `regime.cluster_rotation`.

One `cluster<i>` table is written per entry of `cluster_specs` (the members
in complex mode). With `complex_info` (complex mode) a `[complex]` table holds
every [`ComplexSpec`](@ref) field and a `[complex.realised]` subtable the
realised `n_members`, `mass_total` [M☉], `member_target_masses` [M☉] (the
drawn member masses, in member order), `virial_ratio`, `e_total`
[M☉ km² s⁻²], `l_total` [M☉ pc km s⁻¹], `rotation_fraction`, `softening` [pc]
and `lambda_peebles` = |L| √|E| / (G M^{5/2}) of the member centres.

With `control_info` (a monolithic control) the `[complex]` table holds the
[`ComplexSpec`](@ref) fields only, and a `[control]` table the `kind`, the
targets `e_target` [M☉ km² s⁻²] and `q_target` from the systems of the complex,
the realised `e_realised` and `q_realised` of the control's systems,
`r_h_pc` [pc], `n_stars` and `mass_total` [M☉] of the control.

With `imbh_record` (seeded IMBHs) the top-level `[imbh]` table holds the
per-IMBH fields of the first IMBH: its body name `name`, `cluster`,
`mass_msun`, the realised `mass_fraction` M•/M_stars, `mean_stellar_mass_msun`
(per system), `sigma_c_kms`, `placement`, `velocity_mode`, `response`,
`jeans_factors` (empty for a rigid response), `spin` and the written state
`position_pc` and `velocity_kms` (in the engine's frame). It adds `route`
(`"kz24"` or `"flag"`), `names` and `clusters` of every IMBH in `dat.10`
order, the array of tables `bodies` with the per-IMBH fields of each,
`dtbh_nb` and `rbar_note`: `rbar` is the half-mass radius of the stars, while
`N_total`, `M_total` and `zmbar` include the IMBHs.
"""
function _write_merger_ic_metadata(
    path::AbstractString,
    cfg::MergerConfig,
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    seed::Int,
    external_rng::Bool,
    N_total::Int,
    M_total::Float64,
    rbar::Float64,
    zmbar::Float64,
    regime::NamedTuple;
    d_apo::Float64 = NaN,
    orbit_derived::Union{NamedTuple,Nothing} = nothing,
    mass_truncated::AbstractVector{Float64} = Float64[],
    mass_kept::AbstractVector{Float64} = Float64[],
    truncation_radii_pc::AbstractVector{Float64} = Float64[],
    cluster_specs::AbstractVector{ClusterSpec} = cfg.clusters,
    complex_info::Union{NamedTuple,Nothing} = nothing,
    control_info::Union{NamedTuple,Nothing} = nothing,
    imbh_record::Union{Nothing,Dict{String,Any}} = nothing,
    mbh_record::Union{Nothing,Dict{String,Any}} = nothing,
)
    orbit_table =
        Dict{String,Any}(String(k) => getfield(cfg.orbit, k) for k in fieldnames(OrbitSpec))
    if orbit_derived !== nothing
        orbit_table["apocentre_resolved"] = d_apo
        for (k, v) in pairs(orbit_derived)
            orbit_table[String(k)] = v isa Tuple ? collect(v) : v
        end
    end
    d = Dict{String,Any}(
        "meta" => Dict{String,Any}(
            "generated_at" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
            # 3: explicit-orbit velocities in km s⁻¹; 4: orbit size keys and
            # orientation, derived orbit quantities, seed streams
            "schema_version" => 4,
            "commit" => _source_stamp(_PACKAGE_ROOT),
            "seed" => seed,
            "seed_streams" => cfg.seed_streams,
            "external_rng" => external_rng,
            "N_total" => N_total,
            "M_total" => M_total,
            "rbar" => rbar,
            "zmbar" => zmbar,
            "q_virial" => regime.q_virial,
            "rbar_over_rhm_min" => regime.rbar_over_rhm_min,
            "t_cr_config_nb" => regime.t_cr_config_nb,
            "t_cr_config_myr" => regime.t_cr_config_myr,
            "t_cr_member_min_nb" => regime.t_cr_member_min_nb,
            "t_cr_member_min_myr" => regime.t_cr_member_min_myr,
            "t_star_myr" => regime.t_star_myr,
            "nbin0" => regime.nbin0,
            "fp32_offset_bound" => regime.fp32_offset_bound,
        ),
        "nbody6" => _nbody6_table(regime.nbody6),
        "stellar" => _struct_to_dict(cfg.stellar),
        "tidal" => merge(
            _struct_to_dict(cfg.tidal),
            Dict{String,Any}("omega_kms_pc" => tidal_frame_rotation(cfg.tidal)),
        ),
        "hardware" => _hardware_fingerprint(),
        "orbit_mode" => cfg.orbit_mode,
        "orbit" => orbit_table,
        "output" => Dict{String,Any}(
            "format" => cfg.output.format,
            "truncate_jacobi" => cfg.output.truncate_jacobi,
            "truncate_radius_rh" => cfg.output.truncate_radius_rh,
            "truncate_rule" => cfg.output.truncate_rule,
            "tcrit" => regime.tcrit_nb,
            "dtadj" => regime.dtadj_nb,
            "deltat" => regime.deltat_nb,
            "tcrit_myr" => regime.tcrit_nb * regime.t_star_myr,
            "dtadj_myr" => regime.dtadj_nb * regime.t_star_myr,
            "deltat_myr" => regime.deltat_nb * regime.t_star_myr,
            "dtplot" => regime.dtplot_nb,
        ),
        "cluster_blocks" =>
            [[[first(r), last(r)] for r in blocks] for blocks in regime.cluster_blocks],
    )

    if complex_info !== nothing
        sc = complex_info.scaling
        complex_table =
            Dict{String,Any}(String(k) => getfield(cfg.complex, k) for k in fieldnames(ComplexSpec))
        complex_table["realised"] = Dict{String,Any}(
            "n_members" => length(cluster_specs),
            "mass_total" => sum(complex_info.masses),
            "member_target_masses" => collect(Float64, complex_info.target_masses),
            "virial_ratio" => sc.virial_ratio,
            "e_total" => sc.e_total,
            "l_total" => collect(sc.l_total),
            "rotation_fraction" => sc.rotation_fraction,
            "softening" => sc.softening,
            "lambda_peebles" => complex_info.lambda_peebles,
        )
        d["complex"] = complex_table
    end
    if control_info !== nothing
        d["complex"] =
            Dict{String,Any}(String(k) => getfield(cfg.complex, k) for k in fieldnames(ComplexSpec))
        d["control"] = Dict{String,Any}(
            "kind" => control_info.kind,
            "e_target" => control_info.e_target * _G_PC_KMS2_MSUN,
            "q_target" => control_info.q_target,
            "e_realised" => control_info.e_realised * _G_PC_KMS2_MSUN,
            "q_realised" => control_info.q_realised,
            "r_h_pc" => cluster_specs[1].rbar,
            "n_stars" => N_total,
            "mass_total" => M_total,
        )
    end
    imbh_record === nothing || (d["imbh"] = imbh_record)
    mbh_record === nothing || (d["mbh"] = mbh_record)

    # One table per cluster so spec + post-truncation count are co-located
    for (i, spec) in enumerate(cluster_specs)
        cluster_table = Dict{String,Any}(
            "N" => spec.N,
            "rbar" => spec.rbar,
            "profile" => _profile_table(spec.profile),
            "imf" => _imf_table(spec.imf),
            "position" => spec.position,
            "velocity" => spec.velocity,
            "binaries" => _binaries_table(spec.binaries),
            "segregation" => _segregation_table(spec.segregation),
            "N_after_trunc" => length(cluster_ranges[i]),
            "N_pairs" => regime.n_pairs[i],
            "hard_fraction" => regime.hard_fraction[i],
        )
        isempty(mass_truncated) || (cluster_table["mass_truncated"] = mass_truncated[i])
        isempty(mass_kept) || (cluster_table["mass"] = mass_kept[i])
        isempty(truncation_radii_pc) ||
            (cluster_table["truncation_radius_pc"] = truncation_radii_pc[i])
        rotation = get(regime, :cluster_rotation, Any[])
        if spec.profile isa RotatingKingProfile && i ≤ length(rotation) && rotation[i] !== nothing
            mom = rotation[i].moments
            realised = rotation[i].realised
            cluster_table["rotation"] = Dict{String,Any}(
                "omega0" => spec.profile.omega0,
                "axis" => copy(spec.profile.axis),
                "lambda_model" => mom.λ,
                "t_rot_model" => mom.t_rot,
                "r_h_over_r_c_model" => mom.r_h_over_r_c,
                "lambda_realised" => realised.λ,
                # code velocity unit → km s⁻¹
                "l_axis_realised_pc_kms" => realised.l_axis * _CODE_VSTAR_KMS,
            )
        end
        d["cluster$i"] = cluster_table
    end

    _backup_existing(path)
    _atomic_write_toml(path, d)
    return nothing
end

"""
    load_merger_ic_result(dir::AbstractString) -> MergerICResult

Reconstruct a [`MergerICResult`](@ref) from the `dat.10` and `merger_ic.toml`
files in `dir`, without re-sampling. Useful when the IC plots need to be
regenerated after a code fix — e.g. bug fixes to `plot_merger_ic` or legend
rendering — without changing the underlying particle data.

The particle positions and velocities in the returned result are in
**physical units** (pc and km/s), matching what `plot_merger_ic` expects.
Conversion from the NB-unit `dat.10` uses the stored `rbar` and the mass
scaling (G = 1 in the (M☉, pc, km/s) system with the bundled scaling).
A cluster recorded with an explicit star list (`imf.type = "explicit"`) gets
an `ExplicitMassesIMF` of its masses read back from `dat.10`, with `N` their count.
"""
function load_merger_ic_result(dir::AbstractString)::MergerICResult
    meta_path = joinpath(dir, "merger_ic.toml")
    dat_path = joinpath(dir, "dat.10")
    isfile(meta_path) || error("Missing merger_ic.toml in $dir")
    isfile(dat_path) || error("Missing dat.10 in $dir")

    raw = TOML.parsefile(meta_path)
    meta = raw["meta"]::Dict
    M_total = Float64(meta["M_total"])
    rbar = Float64(meta["rbar"])
    zmbar = Float64(meta["zmbar"])
    N_total = Int(meta["N_total"])

    orbit_mode = String(raw["orbit_mode"])

    # Read cluster specs and ranges — the metadata stores the structured
    # profile/imf tables that `_parse_cluster_table` accepts directly.
    # Parse with the "metadata" mode, which applies no mode rule to
    # position/velocity: metadata always records them (empty in kepler mode).
    cluster_specs = ClusterSpec[]
    n_pairs = Int[]
    # Clusters with an explicit star list (complex members of a global
    # partition, controls): restored from the dat.10 masses below.
    explicit = Int[]
    i = 1
    while haskey(raw, "cluster$i")
        c = raw["cluster$i"]::Dict
        # The metadata table carries the specification plus the realised
        # counts; only the specification keys go to the parser.
        spec_table = Dict{String,Any}(k => v for (k, v) in c if k in _CLUSTER_TABLE_KEYS)
        imf_raw = get(spec_table, "imf", nothing)
        if imf_raw isa AbstractDict && get(imf_raw, "type", "") == "explicit"
            spec_table["imf"] = Dict{String,Any}("type" => "equal", "particle_mass" => 1.0)
            push!(explicit, i)
        end
        push!(cluster_specs, _parse_cluster_table(spec_table, i, "metadata"))
        push!(n_pairs, Int(get(c, "N_pairs", 0)))
        i += 1
    end
    cluster_ranges = [
        _members_from_blocks([Int(b[1]):Int(b[2]) for b in blocks]) for
        blocks in raw["cluster_blocks"]
    ]
    # Schema 3 files carry apocentre and eccentricity only; the keys of
    # schema 4 take their defaults when absent.
    orbit_raw = raw["orbit"]::Dict
    normal_raw = get(orbit_raw, "orbit_normal", nothing)
    orbit = OrbitSpec(;
        apocentre = Float64(orbit_raw["apocentre"]),
        eccentricity = Float64(orbit_raw["eccentricity"]),
        apocentre_rh = Float64(get(orbit_raw, "apocentre_rh", 0.0)),
        pericentre = Float64(get(orbit_raw, "pericentre", 0.0)),
        pericentre_rh = Float64(get(orbit_raw, "pericentre_rh", 0.0)),
        orbit_normal = normal_raw isa AbstractVector ? Float64[Float64(x) for x in normal_raw] :
                       [0.0, 0.0, -1.0],
        apocentre_angle = Float64(get(orbit_raw, "apocentre_angle", 0.0)),
    )

    # Read dat.10 (NB units, the only supported output format) and
    # reconstruct physical-unit copies from the stored scaling.
    format = String(raw["output"]["format"])
    format == "nbody" || error(
        "Unsupported output format \"$format\" in merger_ic.toml; " *
        "only \"nbody\" (KZ(22)=2) is supported.",
    )
    lines = readlines(dat_path)
    N = length(lines)
    N == N_total || @warn "dat.10 row count ($N) differs from metadata N_total ($N_total)"

    mass = zeros(Float64, N)
    pos = zeros(Float64, 3, N)
    vel = zeros(Float64, 3, N)
    for (k, line) in enumerate(lines)
        parts = parse.(Float64, split(line))
        mass[k] = parts[1]
        pos[:, k] = parts[2:4]
        vel[:, k] = parts[5:7]
    end

    mass_phys = mass .* M_total                              # NB mass fraction → M☉
    pos_phys = pos .* rbar                                   # NB length → pc
    vel_phys = vel .* (_CODE_VSTAR_KMS * sqrt(M_total / rbar))  # NB velocity → km/s

    for j in explicit
        s = cluster_specs[j]
        m = mass_phys[cluster_ranges[j]]
        cluster_specs[j] = ClusterSpec(
            length(m),
            s.rbar,
            s.profile,
            ExplicitMassesIMF(m),
            s.position,
            s.velocity,
            s.binaries,
            s.segregation,
        )
    end

    return MergerICResult(
        abspath(dir),
        N_total,
        M_total,
        rbar,
        zmbar,
        cluster_ranges,
        n_pairs,
        cluster_specs,
        orbit_mode,
        orbit,
        mass_phys,
        pos_phys,
        vel_phys,
    )
end

function _write_merger_summary(
    path,
    cfg,
    cluster_data,
    cluster_ranges,
    N_total,
    M_total,
    rbar,
    zmbar,
    regime;
    d_apo::Float64 = NaN,
    orbit_derived::Union{NamedTuple,Nothing} = nothing,
    mass_truncated::AbstractVector{Float64} = Float64[],
    cluster_specs::AbstractVector{ClusterSpec} = cfg.clusters,
    complex_info::Union{NamedTuple,Nothing} = nothing,
    control_info::Union{NamedTuple,Nothing} = nothing,
    imbh_record::Union{Nothing,Dict{String,Any}} = nothing,
)
    n_clusters = length(cluster_specs)
    _backup_existing(path)
    open(path, "w") do io
        println(io, "=" ^ 60)
        println(
            io,
            cfg.orbit_mode == "single" ? "  Single-cluster IC Summary" :
            "  Multi-Cluster Merger IC Summary",
        )
        println(io, "=" ^ 60)
        println(io)
        println(io, "Generated: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
        println(io, "N clusters: ", n_clusters)
        println(io, "Orbit mode: ", cfg.orbit_mode)
        @printf(
            io,
            "Time unit: T* = %.4g Myr; tcrit = %.4g NB = %.4g Myr, deltat = %.4g NB\n",
            regime.t_star_myr,
            regime.tcrit_nb,
            regime.tcrit_nb * regime.t_star_myr,
            regime.deltat_nb,
        )
        println(io)
        for (i, spec) in enumerate(cluster_specs)
            Ni = length(cluster_ranges[i])
            M_sampled = sum(cluster_data[i].mass)  # pre-truncation sampled mass
            @printf(
                io,
                "  Cluster %d: %s, imf=%s, N=%d (after trunc: %d, binaries: %d), M=%.1f M☉, r_hm=%.2f pc",
                i,
                profile_name(spec.profile),
                imf_name(spec.imf),
                spec.N,
                Ni,
                regime.n_pairs[i],
                M_sampled,
                spec.rbar
            )
            if spec.profile isa KingProfile
                @printf(io, ", W0=%.1f", spec.profile.W0)
            end
            if !isempty(spec.position)
                @printf(
                    io,
                    "\n             pos=[%.2f, %.2f, %.2f] pc",
                    spec.position[1],
                    spec.position[2],
                    spec.position[3]
                )
                @printf(
                    io,
                    "  vel=[%.3f, %.3f, %.3f] km/s",
                    spec.velocity[1],
                    spec.velocity[2],
                    spec.velocity[3]
                )
            end
            println(io)
        end
        # Not cluster lines: parse_merger_summary reads the stellar ranges only
        if imbh_record !== nothing
            for b in get(imbh_record, "bodies", [imbh_record])
                @printf(
                    io,
                    "  IMBH: name %d, M = %.4g M☉ (M•/M_stars = %.3g), cluster %d, velocity %s, response %s, route %s\n",
                    b["name"],
                    b["mass_msun"],
                    b["mass_fraction"],
                    b["cluster"],
                    b["velocity_mode"],
                    b["response"],
                    imbh_record["route"]
                )
            end
        end
        println(io)
        if complex_info !== nothing
            cs = cfg.complex
            sc = complex_info.scaling
            m_c = complex_info.masses
            r_c = complex_info.r_h
            @printf(
                io,
                "Complex: %d members, M = %.4g–%.4g M☉ (total %.4g M☉), r_h = %.3g–%.3g pc\n",
                length(m_c),
                minimum(m_c),
                maximum(m_c),
                sum(m_c),
                minimum(r_c),
                maximum(r_c)
            )
            @printf(
                io,
                "         centres %s (radius %.3g pc), mass function %s, size–mass %s\n",
                cs.centres,
                cs.radius,
                cs.mass_function,
                cs.size_mass
            )
            @printf(
                io,
                "         realised Q_cc = %.4f (target %.4f), rotation fraction %.3f, λ = %.3g, ε = %.3g pc\n",
                sc.virial_ratio,
                cs.virial_ratio,
                sc.rotation_fraction,
                complex_info.lambda_peebles,
                sc.softening
            )
            println(io)
        end
        if control_info !== nothing
            @printf(
                io,
                "Control: %s, matched E = %.6g M☉ km² s⁻² (realised %.6g), Q = %.4f (realised %.4f)\n",
                control_info.kind == "monolithic" ? "monolithic (C-E)" : "cold (C-cold)",
                control_info.e_target * _G_PC_KMS2_MSUN,
                control_info.e_realised * _G_PC_KMS2_MSUN,
                control_info.q_target,
                control_info.q_realised
            )
            println(io)
        end
        if cfg.orbit_mode == "kepler"
            rh1 = cfg.clusters[1].rbar
            ecc = cfg.orbit.eccentricity
            d_orb = isnan(d_apo) ? orbit_apocentre(cfg.orbit, rh1) : d_apo
            r_p = d_orb * (1.0 - ecc) / (1.0 + ecc)
            @printf(
                io,
                "Orbit: d_apo = %.2f pc (%.2f r_h,1), r_p = %.2f pc (%.2f r_h,1), e = %.3f\n",
                d_orb,
                d_orb / rh1,
                r_p,
                r_p / rh1,
                ecc
            )
            if orbit_derived !== nothing
                q = orbit_derived
                @printf(
                    io,
                    "       a = %.2f pc, P = %.3g Myr, j_orb = %.3f, λ = %.3g, v_p/σ₁ = %.3g\n",
                    q.a_pc,
                    q.period_myr,
                    q.j_orb,
                    q.lambda_peebles,
                    q.v_peri_over_sigma
                )
                @printf(
                    io,
                    "       Roche filling at apocentre: %.3f, %.3f\n",
                    q.roche_filling[1],
                    q.roche_filling[2]
                )
            end
            n_orb = cfg.orbit.orbit_normal
            @printf(io, "       orbit normal = [%.3g, %.3g, %.3g]\n", n_orb[1], n_orb[2], n_orb[3])
            if any(>(0), mass_truncated)
                @printf(
                    io,
                    "       truncated mass: %.1f M☉ (cluster 1), %.1f M☉ (cluster 2)\n",
                    mass_truncated[1],
                    mass_truncated[2]
                )
            end
            println(io)
        end
        if regime.nbin0 > 0
            @printf(
                io,
                "Primordial binaries: NBIN0 = %d pairs written first (KZ(8) = 2); hard fraction per cluster: %s\n",
                regime.nbin0,
                join([isnan(f) ? "-" : @sprintf("%.2f", f) for f in regime.hard_fraction], ", ")
            )
            println(io)
        end
        @printf(io, "Combined: N_total = %d, M_total = %.1f M☉\n", N_total, M_total)
        @printf(io, "          RBAR = %.4f pc, ZMBAR = %.4f M☉\n", rbar, zmbar)
        @printf(
            io,
            "          virial ratio Q = T/|W| = %.3f (0.5 = equilibrium), RBAR/r_hm,min = %.2f\n",
            regime.q_virial,
            regime.rbar_over_rhm_min
        )
        @printf(
            io,
            "          crossing time: configuration %.3g NB = %.3g Myr, smallest member %.3g NB = %.3g Myr (T* = %.3g Myr)\n",
            regime.t_cr_config_nb,
            regime.t_cr_config_myr,
            regime.t_cr_member_min_nb,
            regime.t_cr_member_min_myr,
            regime.t_star_myr
        )
        p = regime.nbody6
        @printf(
            io,
            "Integration (merger.inp, N-body units): preset=%s binary_physics=%s QE=%.1E ETAI=%.3g ETAR=%.3g NNBOPT=%d RS0=%.3g RMIN=%.3E DTMIN=%.3E KZ(16)=%d TCRTP0=%.4g Myr\n",
            p.preset,
            p.binary_physics,
            p.qe,
            p.etai,
            p.etar,
            p.nnbopt,
            p.rs0,
            p.rmin,
            p.dtmin,
            p.kz16,
            p.tcrtp0
        )
        st = cfg.stellar
        @printf(
            io,
            "Stellar evolution: KZ(19)=%d Level=%s ZMET=%.4g EPOCH0=%.3g DTPLOT=%.3g\n",
            st.kz19,
            st.level,
            st.zmet,
            st.epoch0,
            regime.dtplot_nb
        )
        td = cfg.tidal
        if td.kz14 == 0
            println(io, "External field: none (KZ(14)=0, isolated)")
        elseif td.kz14 == 1
            println(io, "External field: KZ(14)=1 solar-neighbourhood linearised tide")
        elseif td.kz14 == 2
            @printf(
                io,
                "External field: KZ(14)=2 point-mass galaxy GMG=%.4g M☉ at RG0=%.4g kpc\n",
                td.gmg,
                td.rg0
            )
        else
            @printf(
                io,
                "External field: KZ(14)=5 MWPotential2014, RG=[%.3g, %.3g, %.3g] kpc VG=[%.3g, %.3g, %.3g] km/s\n",
                td.rg[1],
                td.rg[2],
                td.rg[3],
                td.vg[1],
                td.vg[2],
                td.vg[3]
            )
        end
        println(io)
        println(io, "Output format: ", cfg.output.format, " (KZ(22)=2)")
        if cfg.orbit_mode == "single"
            println(io, "Jacobi truncation: not applicable (single cluster)")
        else
            println(io, "Jacobi truncation: ", cfg.output.truncate_jacobi)
        end
        println(io)
        println(io, "Files:")
        println(io, "  dat.10     — particle data")
        println(io, "  merger.inp — Nbody6++ input file")
        println(io, "=" ^ 60)
    end
    @info "Summary written to $path"
end
