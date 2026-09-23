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

Convert from physical units (M☉, pc, km/s) to Hénon N-body units
(G=1, M_total=1, E=-1/4) in-place.

# N-body scaling
- Mass: `m_nb = m_solar / M_total_solar`
- Length: `r_nb = r_pc / rbar_pc` (rbar = virial radius in pc)
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
solved model, with `r_h` located on the cumulative-mass profile.
"""
_central_density_contrast(::PlummerProfile) = 2.0 * _PLUMMER_RHM_OVER_A^3
function _central_density_contrast(p::KingProfile)
    rhat, ρ, m_cum = _king_cumulative(p.W0)
    i_h = findfirst(≥(0.5 * m_cum[end]), m_cum)
    r_h = rhat[i_h]
    ρ_mean = 0.5 * m_cum[end] / (4π / 3 * r_h^3)
    return ρ[1] / ρ_mean
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

Fill the zero entries of `spec` from the member clusters, in N-body units of
the combined system (length unit `rbar_pc`). With `r_h` the smallest member
half-mass radius in those units, `N_min` the smallest post-truncation
membership, and `ρ̂` the central density contrast of that cluster's profile
([`_central_density_contrast`](@ref)):

- `NNBOPT = clamp(round(√N_total), 20, 300)`
- `RS0 = 2 r_h (2 NNBOPT / N_min)^{1/3}`, twice the radius enclosing about
  `NNBOPT` stars at the mean density of the half-mass sphere, capped at
  `r_h`. The factor 2 matches the engine's own example inputs
  (`RS0 ≈ 0.5 r_h` for `NNBOPT ≈ √N`). The engine regrows an empty
  neighbour sphere itself (`fpoly0.F`), so the radius is a matter of
  start-up cost, not of correctness; the start-up hang once attributed to
  it was the interval digit counter ([`engine_interval`](@ref)).
- `RMIN = 4 r_h / (N_min ρ̂^{1/3})`, the functional form of the engine's own
  re-derivation (`adjust.F`) evaluated for the member cluster instead of the
  whole configuration
- `DTMIN = 0.04 √(ETAI/0.02) √(RMIN³ N_total)`, the engine's rule with the
  mean body mass `1/N_total`

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
    nnbopt = spec.nnbopt > 0 ? spec.nnbopt : clamp(round(Int, sqrt(N_total)), 20, 300)
    rs0 = spec.rs0 > 0 ? spec.rs0 : min(2 * r_h * cbrt(2 * nnbopt / N_min), r_h)
    ρ̂ = _central_density_contrast(clusters[i_min].profile)
    rmin = spec.rmin > 0 ? spec.rmin : 4 * r_h / (N_min * cbrt(ρ̂))
    dtmin = spec.dtmin > 0 ? spec.dtmin : 0.04 * sqrt(spec.etai / 0.02) * sqrt(rmin^3 * N_total)
    kw = Dict{Symbol,Any}(k => getfield(spec, k) for k in fieldnames(Nbody6ParameterSpec))
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
    for k in (:nnbopt, :rs0, :rmin, :dtmin)
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
diagnostics), and when the snapshot interval `deltat` exceeds the smallest
member crossing time (member dynamics undersampled).
"""
function _check_multicluster_regime(
    q_virial::Float64,
    rbar_over_rhm_min::Float64,
    nbody6::Nbody6ParameterSpec,
    r_h_min_nb::Float64;
    deltat::Float64 = 0.0,
    t_cr_member_min_nb::Float64 = NaN,
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
    if !isnan(t_cr_member_min_nb) && deltat > t_cr_member_min_nb
        @warn "Snapshot interval deltat = $deltat exceeds the smallest member crossing time " *
              "$(round(t_cr_member_min_nb; sigdigits = 3)) (N-body units): member dynamics are " *
              "undersampled by conf.3; lower merger.output.deltat"
    end
    return nothing
end

"""
    generate_merger_inp(path, N_total, rbar, zmbar; nbody6, stellar = StellarSpec(),
                        kz14 = 0, tcrit, dtadj, deltat, nrand, mass_bounds = (0.08, 100.0))

Generate an Nbody6++ `.inp` file configured for external particle
input via `dat.10` in N-body units (`KZ(22) = 2`, the only supported
input mode).

Follows the Fortran NAMELIST read order expected by `nbody6.F → start.F`:
  `&INNBODY6` → `&ININPUT` → `&INSSE` → `&INBSE` → `&INCOLL` →
  `&INDATA` → `&INSCALE` (→ `&INXTRNL0` if KZ(14)>0).

# Key settings
- `nbody6`: resolved [`Nbody6ParameterSpec`](@ref); obtain it from
  [`resolve_nbody6_parameters`](@ref) — unresolved zeros are refused. Its
  `kz` overrides are applied last.
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
)
    _assert_resolved(nbody6)
    _validate_tidal(tidal)
    nbin0 ≥ 0 || throw(ArgumentError("nbin0 must be ≥ 0; got $nbin0"))
    kz14 = tidal.kz14
    kz = zeros(Int, 50)
    kz[1] = 1
    kz[2] = -1
    kz[3] = 2
    kz[7] = 3
    kz[8] = nbin0 > 0 ? 2 : 0   # 2: primordial pairs are the first 2·NBIN0 bodies of dat.10
    kz[12] = stellar.kz19 == 0 ? 0 : 1
    kz[14] = kz14
    kz[16] = nbody6.kz16
    kz[19] = stellar.kz19
    kz[22] = 2
    kz[23] = 2
    kz[26] = 1
    kz[30] = 1
    for (i, v) in nbody6.kz
        i in (8, 14, 16, 19) &&
            @warn "merger.nbody6.kz overrides KZ($i) = $v, which also has a named setting (binaries, kz14/tidal, kz16, stellar.kz19)"
        kz[i] = v
    end

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
        println(io, "&INBSE /")
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

        # --- No binaries (NBIN0=0) or hierarchical triples ---
    end

    @info "Wrote .inp file: $path (N=$N_total, NBIN0=$nbin0, KZ(22)=2, KZ(14)=$kz14, RBAR=$rbar, ZMBAR=$zmbar)"
    return nothing
end

# ---------------------------------------------------------------------------
# Internal: sample a single cluster from its spec
# ---------------------------------------------------------------------------
function _sample_cluster(spec::ClusterSpec; rng::AbstractRNG = Random.default_rng(), nmax::Int)
    # Total-mass handling is delegated to the IMFSpec subtype (natural Kroupa
    # keeps the sampled sum; rescaled/equal enforce their targets).
    masses = sample_masses(spec.imf, spec.N, rng)

    # Primordial pairs: the density sampler places *systems* (pairs by their
    # centre of mass); pairs are expanded into bodies after the clusters are
    # combined and truncated. Systems are ordered pairs first, then singles.
    bins = sample_binaries(spec.binaries, masses, rng)
    n_b = length(bins.primary)
    paired = falses(spec.N)
    paired[bins.primary] .= true
    paired[bins.secondary] .= true
    singles = findall(!, paired)
    sys_mass = vcat(bins.m1 .+ bins.m2, bins.masses[singles])
    system_binary = vcat(collect(1:n_b), zeros(Int, length(singles)))
    N_sys = length(sys_mass)

    pos, vel = if spec.profile isa PlummerProfile
        # Scale radius from the target half-mass radius (r_hm = 1.305 a)
        sample_plummer(N_sys, spec.rbar / _PLUMMER_RHM_OVER_A; rng = rng)
    elseif spec.profile isa KingProfile
        # Sampled with unit tidal radius; the empirical rescale below sets r_hm
        sample_king(N_sys, spec.profile.W0, 1.0; rng = rng)
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
    return (
        pos = pos,
        vel = vel,
        mass = sys_mass,
        T = energies.T,
        W = energies.W,
        binaries = binaries,
    )
end

"""
    generate_merger_ic(cfg::MergerConfig; rng=Random.default_rng(),
                       output_dir="") -> MergerICResult

Generate complete merger initial conditions from a [`MergerConfig`](@ref).
Produces `dat.10`, `merger.inp`, and `merger_summary.txt`.

Supports two orbit modes:
- `"kepler"`: 2 clusters on a Keplerian orbit (positions/velocities auto-computed)
- `"explicit"`: N clusters with user-specified CM positions and velocities

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
    for (i, spec) in enumerate(cfg.clusters)
        _validate_cluster_spec(spec, i)
    end
    _validate_nbody6(cfg.nbody6)
    _validate_stellar(cfg.stellar, cfg.output)
    _validate_tidal(cfg.tidal)
    _validate_tidal_tolerance(cfg.tidal, cfg.nbody6)

    # Resolve RNG and effective seed. When the caller supplies an RNG the
    # sampling is NOT reproducible from `effective_seed`; the metadata records
    # this honestly via `external_rng` (the seed still feeds Nbody6's NRAND).
    external_rng = rng !== nothing
    effective_seed::Int = cfg.seed === nothing ? Int(rand(UInt32) % typemax(Int32)) : cfg.seed
    use_rng = external_rng ? rng : Random.MersenneTwister(effective_seed)

    seed_note = external_rng ? "external RNG (seed not reproducible)" : "seed=$effective_seed"
    @info "Generating merger ICs for $n_clusters clusters ($(cfg.orbit_mode) mode, $seed_note)..."

    # Sample each cluster independently
    cluster_data = [
        begin
            @info "  Cluster $i: $(profile_name(spec.profile)) profile, N=$(spec.N), " *
                  "imf=$(imf_name(spec.imf)), M≈$(round(expected_mass(spec.imf, spec.N); digits=1)) M☉"
            _sample_cluster(spec; rng = use_rng, nmax = cfg.virial_max_n)
        end for (i, spec) in enumerate(cfg.clusters)
    ]

    # Combine clusters (as systems) according to orbit mode; both paths
    # return the combined arrays, the per-cluster index ranges of the systems
    # that survived truncation, and their indices in each sampled set.
    pos_combined, vel_combined, mass_combined, system_ranges, kept = if cfg.orbit_mode == "kepler"
        c1, c2 = cluster_data[1], cluster_data[2]
        setup_two_cluster_orbit(
            c1.pos,
            c1.vel,
            c1.mass,
            c2.pos,
            c2.vel,
            c2.mass,
            cfg.orbit.apocentre,
            cfg.orbit.eccentricity;
            truncate_jacobi_flag = cfg.output.truncate_jacobi,
        )
    elseif cfg.orbit_mode == "explicit"
        combine_clusters_explicit(
            cluster_data,
            cfg.clusters;
            truncate_jacobi_flag = cfg.output.truncate_jacobi,
        )
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
    q_virial, t_cr_config_code = if N_systems ≤ cfg.virial_max_n
        T, W = _kinetic_and_potential(mass_combined, pos_combined, vel_combined)
        (W < 0 ? T / abs(W) : NaN, crossing_time(M_total, T + W))
    else
        @warn "Combined virial ratio not evaluated: N = $N_systems systems exceed merger.virial_max_n = $(cfg.virial_max_n) (O(N²) pair sum)"
        (NaN, NaN)
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
        rng = use_rng,
    )
    pos_combined, vel_combined, mass_combined = expanded.pos, expanded.vel, expanded.mass
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
    r_h_min_pc = minimum(c.rbar for c in cfg.clusters)
    rbar_over_rhm_min = rbar / r_h_min_pc
    nbody6 = resolve_nbody6_parameters(cfg.nbody6, cfg.clusters, cluster_ranges, N_total, rbar)
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
    # output format; decision D4). `to_nbody_units!` expects velocity in
    # km/s, so convert from code units first.
    cfg.output.format == "nbody" || error(
        "Unsupported output format \"$(cfg.output.format)\"; only \"nbody\" (KZ(22)=2) is supported.",
    )
    vel_combined .*= _CODE_VSTAR_KMS  # code units → km/s (match to_nbody_units! API)
    to_nbody_units!(mass_combined, pos_combined, vel_combined, M_total, rbar)

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
            minimum(_imf_mass_bounds(c.imf)[1] for c in cfg.clusters),
            maximum(_imf_mass_bounds(c.imf)[2] for c in cfg.clusters),
        ),
        nbin0 = nbin0,
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
        regime,
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
        regime,
    )

    return MergerICResult(
        out_dir,
        N_total,
        M_total,
        rbar,
        zmbar,
        cluster_ranges,
        n_pairs,
        collect(cfg.clusters),
        cfg.orbit_mode,
        cfg.orbit,
        mass_phys,
        pos_phys,
        vel_phys,
    )
end

"""
    _write_merger_ic_metadata(path, cfg, cluster_ranges, seed, external_rng,
                              N_total, M_total, rbar, zmbar, regime)

Write a machine-readable TOML snapshot of the IC generation (schema v2),
sufficient to reconstruct a [`MergerICResult`](@ref) from the `dat.10` file
later (so `plot_merger_ic` can be re-run after code fixes without
re-sampling). Cluster specs are stored in the structured form
(`profile = {type=...}`, `imf = {type=...}`) that
[`_parse_cluster_table`](@ref) also accepts, so the file round-trips.
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
    regime::NamedTuple,
)
    d = Dict{String,Any}(
        "meta" => Dict{String,Any}(
            "generated_at" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
            "schema_version" => 3,   # 3: explicit-orbit velocities in km s⁻¹
            "commit" => _source_stamp(_PACKAGE_ROOT),
            "seed" => seed,
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
        ),
        "nbody6" => _nbody6_table(regime.nbody6),
        "stellar" => _struct_to_dict(cfg.stellar),
        "tidal" => _struct_to_dict(cfg.tidal),
        "hardware" => _hardware_fingerprint(),
        "orbit_mode" => cfg.orbit_mode,
        "orbit" => Dict{String,Any}(
            "apocentre" => cfg.orbit.apocentre,
            "eccentricity" => cfg.orbit.eccentricity,
        ),
        "output" => Dict{String,Any}(
            "format" => cfg.output.format,
            "truncate_jacobi" => cfg.output.truncate_jacobi,
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

    # One table per cluster so spec + post-truncation count are co-located
    for (i, spec) in enumerate(cfg.clusters)
        d["cluster$i"] = Dict{String,Any}(
            "N" => spec.N,
            "rbar" => spec.rbar,
            "profile" => _profile_table(spec.profile),
            "imf" => _imf_table(spec.imf),
            "position" => spec.position,
            "velocity" => spec.velocity,
            "binaries" => _binaries_table(spec.binaries),
            "N_after_trunc" => length(cluster_ranges[i]),
            "N_pairs" => regime.n_pairs[i],
            "hard_fraction" => regime.hard_fraction[i],
        )
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
    # Parse with orbit_mode "kepler" semantics to skip position/velocity
    # validation — metadata always records them explicitly (possibly empty
    # for kepler mode).
    cluster_specs = ClusterSpec[]
    n_pairs = Int[]
    i = 1
    while haskey(raw, "cluster$i")
        c = raw["cluster$i"]::Dict
        # The metadata table carries the specification plus the realised
        # counts; only the specification keys go to the parser.
        spec_table = Dict{String,Any}(k => v for (k, v) in c if k in _CLUSTER_TABLE_KEYS)
        push!(cluster_specs, _parse_cluster_table(spec_table, i, "kepler"))
        push!(n_pairs, Int(get(c, "N_pairs", 0)))
        i += 1
    end
    cluster_ranges = [
        _members_from_blocks([Int(b[1]):Int(b[2]) for b in blocks]) for
        blocks in raw["cluster_blocks"]
    ]
    orbit_raw = raw["orbit"]::Dict
    orbit = OrbitSpec(;
        apocentre = Float64(orbit_raw["apocentre"]),
        eccentricity = Float64(orbit_raw["eccentricity"]),
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
    regime,
)
    n_clusters = length(cfg.clusters)
    _backup_existing(path)
    open(path, "w") do io
        println(io, "=" ^ 60)
        println(io, "  Multi-Cluster Merger IC Summary")
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
        for (i, spec) in enumerate(cfg.clusters)
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
        println(io)
        if cfg.orbit_mode == "kepler"
            @printf(
                io,
                "Orbit: d_apo = %.2f pc, e = %.3f\n",
                cfg.orbit.apocentre,
                cfg.orbit.eccentricity
            )
            a_orb = cfg.orbit.apocentre / (1.0 + cfg.orbit.eccentricity)
            @printf(io, "       a = %.2f pc (semi-major axis)\n", a_orb)
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
            "Integration (merger.inp, N-body units): QE=%.1E ETAI=%.3g ETAR=%.3g NNBOPT=%d RS0=%.3g RMIN=%.3E DTMIN=%.3E KZ(16)=%d TCRTP0=%.4g Myr\n",
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
            st.dtplot
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
        println(io, "Jacobi truncation: ", cfg.output.truncate_jacobi)
        println(io)
        println(io, "Files:")
        println(io, "  dat.10     — particle data")
        println(io, "  merger.inp — Nbody6++ input file")
        println(io, "=" ^ 60)
    end
    @info "Summary written to $path"
end
