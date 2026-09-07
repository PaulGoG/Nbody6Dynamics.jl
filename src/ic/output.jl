# =============================================================================
# dat.10 writer and .inp generator for Nbody6++ merger ICs
# =============================================================================

# Velocity unit of the internal (G = 1, M☉, pc) code-unit system:
# sqrt(G M☉ / pc) in km/s. Sampled and Kepler velocities carry this unit
# until the NB-unit conversion for dat.10.
const _CODE_VSTAR_KMS = 0.06557

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
    rhat, _, ρ = _solve_king(p.W0)
    n = length(rhat)
    m_cum = zeros(Float64, n)
    for i in 2:n
        f_prev = 4π * rhat[i - 1]^2 * ρ[i - 1]
        f_here = 4π * rhat[i]^2 * ρ[i]
        m_cum[i] = m_cum[i - 1] + 0.5 * (f_prev + f_here) * (rhat[i] - rhat[i - 1])
    end
    i_h = findfirst(≥(0.5 * m_cum[end]), m_cum)
    r_h = rhat[i_h]
    ρ_mean = 0.5 * m_cum[end] / (4π / 3 * r_h^3)
    return ρ[1] / ρ_mean
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
- `RS0 = r_h (2 NNBOPT / N_min)^{1/3}`, the radius enclosing about `NNBOPT`
  stars at the mean density of the half-mass sphere, capped at `r_h`
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
    cluster_ranges::Vector{UnitRange{Int}},
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
    rs0 = spec.rs0 > 0 ? spec.rs0 : min(r_h * cbrt(2 * nnbopt / N_min), r_h)
    ρ̂ = _central_density_contrast(clusters[i_min].profile)
    rmin = spec.rmin > 0 ? spec.rmin : 4 * r_h / (N_min * cbrt(ρ̂))
    dtmin = spec.dtmin > 0 ? spec.dtmin : 0.04 * sqrt(spec.etai / 0.02) * sqrt(rmin^3 * N_total)
    return Nbody6ParameterSpec(;
        qe = spec.qe,
        etai = spec.etai,
        etar = spec.etar,
        nnbopt = nnbopt,
        rs0 = rs0,
        rmin = rmin,
        dtmin = dtmin,
        kz16 = spec.kz16,
    )
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
    _check_multicluster_regime(q_virial, rbar_over_rhm_min, nbody6, r_h_min_nb)

Generation-time guards against configurations the engine mis-integrates or
mis-diagnoses: refuse a neighbour radius `rs0` wider than the smallest
member half-mass radius (every neighbour list would overflow at start-up);
warn when the combined virial ratio is below `$(_Q_COLD_COLLAPSE)`
(cold-collapse regime, global diagnostics meaningless until the remnant
forms) and when the length unit exceeds `$(_RBAR_OVER_RHM_UNRESOLVED)`
member half-mass radii (members unresolved by the single-centre
diagnostics).
"""
function _check_multicluster_regime(
    q_virial::Float64,
    rbar_over_rhm_min::Float64,
    nbody6::Nbody6ParameterSpec,
    r_h_min_nb::Float64,
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
    return nothing
end

"""
    generate_merger_inp(path, N_total, rbar, zmbar; nbody6, kz14 = 0, tcrit, dtadj, deltat, nrand)

Generate an Nbody6++ `.inp` file configured for external particle
input via `dat.10` in N-body units (`KZ(22) = 2`, the only supported
input mode).

Follows the Fortran NAMELIST read order expected by `nbody6.F → start.F`:
  `&INNBODY6` → `&ININPUT` → `&INSSE` → `&INBSE` → `&INCOLL` →
  `&INDATA` → `&INSCALE` (→ `&INXTRNL0` if KZ(14)>0).

# Key settings
- `nbody6`: resolved [`Nbody6ParameterSpec`](@ref) (`QE`, `ETAI`, `ETAR`,
  `NNBOPT`, `RS0`, `RMIN`, `DTMIN`, `KZ(16)`); obtain it from
  [`resolve_nbody6_parameters`](@ref) — unresolved zeros are refused
- `KZ(14) = kz14`: tidal field option (0 = isolated)
"""
function generate_merger_inp(
    path::AbstractString,
    N_total::Int,
    rbar::Float64,
    zmbar::Float64;
    nbody6::Nbody6ParameterSpec,
    kz14::Int = 0,
    tcrit::Float64 = 100.0,
    dtadj::Float64 = 1.0,
    deltat::Float64 = 1.0,
    nrand::Int = 10000,
)
    _assert_resolved(nbody6)
    kz = zeros(Int, 50)
    kz[1] = 1
    kz[2] = -1
    kz[3] = 2
    kz[7] = 3
    kz[12] = 1
    kz[14] = kz14
    kz[16] = nbody6.kz16
    kz[19] = 3
    kz[22] = 2
    kz[23] = 2
    kz[26] = 1
    kz[30] = 1

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
        println(io, "KSTART=1,TCOMP=1.0E8,TCRTP0=3600.0,isernb=40,iserreg=40,iserks=0 /")
        println(io)

        # --- 2. &ININPUT: main simulation parameters + KZ options ---
        println(io, "&ININPUT")
        @printf(
            io,
            "N=%d,NFIX=1,NCRIT=10,NRAND=%d,NNBOPT=%d,NRUN=1,NCOMM=10,\n",
            N_total,
            abs(nrand) % typemax(Int32),
            nbody6.nnbopt
        )
        @printf(
            io,
            "ETAI=%.4G,ETAR=%.4G,RS0=%.4G,DTADJ=%.4f,DELTAT=%.4f,TCRIT=%.2f,QE=%.3E,RBAR=%.6f,ZMBAR=%.6f,\n",
            nbody6.etai,
            nbody6.etar,
            nbody6.rs0,
            dtadj,
            deltat,
            tcrit,
            nbody6.qe,
            rbar,
            zmbar
        )
        println(io, join(kz_strs, "\n"))
        @printf(
            io,
            "DTMIN=%.3E,RMIN=%.3E,ETAU=0.1,ECLOSE=1.0,GMIN=1.0E-06,GMAX=0.01,SMAX=1.0,\n",
            nbody6.dtmin,
            nbody6.rmin
        )
        println(io, "Level='C' /")
        println(io)

        # --- 3-5. SSE/BSE/Coll: empty → use Level defaults ---
        println(io, "&INSSE /")
        println(io)
        println(io, "&INBSE /")
        println(io)
        println(io, "&INCOLL /")
        println(io)

        # --- 6. &INDATA: IMF and mass function ---
        # Inert under KZ(22)=2: masses come from dat.10, so the IMF fields
        # here are never sampled; they are namelist boilerplate only.
        println(io, "&INDATA")
        println(
            io,
            "ALPHAS=2.35,BODY1=150.0,BODYN=0.08,NBIN0=0,NHI0=0,ZMET=0.001,EPOCH0=0,DTPLOT=1.0 /",
        )
        println(io)

        # --- 7. &INSETUP: placeholder (no Fortran NAMELIST reads this) ---
        println(io, "&INSETUP SEMI=,ECC=,APO=,N2=,SCALE=,ZM1=,ZM2,ZMH,RCUT= /")
        println(io)

        # --- 8. &INSCALE: virial ratio, rotation, tidal ---
        println(io, "&INSCALE")
        println(io, "Q=0.5,VXROT=0.0,VZROT=0.0,RTIDE=0.0 /")
        println(io)

        # --- 9. External potential (only needed if KZ(14)>0) ---
        if kz14 > 0
            println(io, "&INXTRNL0")
            println(
                io,
                "GMG=1.78E11,RG0=13.3,DISK=,A=,B=,VCIRC=,RCIRC=,GMB=,AR=,GAM=,RG=,,,VG=,,,MP=,AP2=,MPDOT=,TDELAY= /",
            )
            println(io)
        end

        # --- No binaries (NBIN0=0) or hierarchical triples ---
    end

    @info "Wrote .inp file: $path (N=$N_total, KZ(22)=2, RBAR=$rbar, ZMBAR=$zmbar)"
    return nothing
end

# ---------------------------------------------------------------------------
# Internal: sample a single cluster from its spec
# ---------------------------------------------------------------------------
function _sample_cluster(spec::ClusterSpec; rng::AbstractRNG = Random.default_rng())
    # Total-mass handling is delegated to the IMFSpec subtype (natural Kroupa
    # keeps the sampled sum; rescaled/equal enforce their targets).
    masses = sample_masses(spec.imf, spec.N, rng)

    pos, vel = if spec.profile isa PlummerProfile
        # Scale radius from the target half-mass radius (r_hm = 1.305 a)
        sample_plummer(spec.N, spec.rbar / _PLUMMER_RHM_OVER_A; rng = rng)
    elseif spec.profile isa KingProfile
        # Sampled with unit tidal radius; the empirical rescale below sets r_hm
        sample_king(spec.N, spec.profile.W0, 1.0; rng = rng)
    else
        error("Unknown profile: $(typeof(spec.profile))")
    end

    # Scale positions so the mass-based half-mass radius equals the target
    r_hm = half_mass_radius(masses, pos; centre = zeros(3))
    if r_hm > 0
        pos .*= spec.rbar / r_hm
    end

    virialise!(masses, pos, vel)
    return (pos = pos, vel = vel, mass = masses)
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
    mkpath(out_dir)

    n_clusters = length(cfg.clusters)
    n_clusters ≥ 2 || error("Need at least 2 clusters, got $n_clusters")
    for (i, spec) in enumerate(cfg.clusters)
        _validate_cluster_spec(spec, i)
    end
    _validate_nbody6(cfg.nbody6)

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
            _sample_cluster(spec; rng = use_rng)
        end for (i, spec) in enumerate(cfg.clusters)
    ]

    # Combine clusters according to orbit mode; both paths return the
    # combined arrays plus per-cluster (post-truncation) index ranges.
    pos_combined, vel_combined, mass_combined, cluster_ranges = if cfg.orbit_mode == "kepler"
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

    N_total = length(mass_combined)
    M_total = sum(mass_combined)
    zmbar = M_total / N_total

    # Half-mass radius of the combined system — this is the RBAR length unit
    # written to the .inp file and used for the NB-unit conversion of dat.10.
    rbar = half_mass_radius(mass_combined, pos_combined)

    # Combined-system diagnostics in code units (G = 1): the virial ratio of
    # the whole configuration, the resolution of the members by the length
    # unit, and the integration parameters scaled to the smallest member.
    q_virial = if N_total ≤ _VIRIAL_NMAX
        T, W = _kinetic_and_potential(mass_combined, pos_combined, vel_combined)
        W < 0 ? T / abs(W) : NaN
    else
        @warn "Combined virial ratio not evaluated: N_total = $N_total exceeds $(_VIRIAL_NMAX) (O(N²) pair sum)"
        NaN
    end
    r_h_min_pc = minimum(c.rbar for c in cfg.clusters)
    rbar_over_rhm_min = rbar / r_h_min_pc
    nbody6 = resolve_nbody6_parameters(cfg.nbody6, cfg.clusters, cluster_ranges, N_total, rbar)
    _check_multicluster_regime(q_virial, rbar_over_rhm_min, nbody6, r_h_min_pc / rbar)
    regime = (q_virial = q_virial, rbar_over_rhm_min = rbar_over_rhm_min, nbody6 = nbody6)
    @info @sprintf(
        "Combined system: Q = %.3f, RBAR/r_hm,min = %.2f; RS0 = %.3g, RMIN = %.3g, DTMIN = %.3g, NNBOPT = %d (N-body units)",
        q_virial,
        rbar_over_rhm_min,
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
        tcrit = cfg.output.tcrit,
        dtadj = cfg.output.dtadj,
        deltat = cfg.output.deltat,
        nrand = effective_seed,
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
    cluster_ranges::Vector{UnitRange{Int}},
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
            "schema_version" => 2,
            "commit" => _git_commit(_PROJECT_ROOT),
            "seed" => seed,
            "external_rng" => external_rng,
            "N_total" => N_total,
            "M_total" => M_total,
            "rbar" => rbar,
            "zmbar" => zmbar,
            "q_virial" => regime.q_virial,
            "rbar_over_rhm_min" => regime.rbar_over_rhm_min,
        ),
        "nbody6" => _struct_to_dict(regime.nbody6),
        "hardware" => _hardware_fingerprint(),
        "orbit_mode" => cfg.orbit_mode,
        "orbit" => Dict{String,Any}(
            "apocentre" => cfg.orbit.apocentre,
            "eccentricity" => cfg.orbit.eccentricity,
        ),
        "output" => Dict{String,Any}(
            "format" => cfg.output.format,
            "truncate_jacobi" => cfg.output.truncate_jacobi,
            "tcrit" => cfg.output.tcrit,
            "dtadj" => cfg.output.dtadj,
            "deltat" => cfg.output.deltat,
        ),
        "cluster_ranges" => [[first(r), last(r)] for r in cluster_ranges],
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
            "N_after_trunc" => length(cluster_ranges[i]),
        )
    end

    _backup_existing(path)
    open(path, "w") do io
        TOML.print(io, d)
    end
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
    cluster_ranges = UnitRange{Int}[]
    i = 1
    while haskey(raw, "cluster$i")
        c = raw["cluster$i"]::Dict
        push!(cluster_specs, _parse_cluster_table(c, i, "kepler"))
        i += 1
    end
    for pair in raw["cluster_ranges"]
        push!(cluster_ranges, Int(pair[1]):Int(pair[2]))
    end
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
        println(io)
        for (i, spec) in enumerate(cfg.clusters)
            Ni = length(cluster_ranges[i])
            M_sampled = sum(cluster_data[i].mass)  # pre-truncation sampled mass
            @printf(
                io,
                "  Cluster %d: %s, imf=%s, N=%d (after trunc: %d), M=%.1f M☉, r_hm=%.2f pc",
                i,
                profile_name(spec.profile),
                imf_name(spec.imf),
                spec.N,
                Ni,
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
                    "  vel=[%.2f, %.2f, %.2f]",
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
        @printf(io, "Combined: N_total = %d, M_total = %.1f M☉\n", N_total, M_total)
        @printf(io, "          RBAR = %.4f pc, ZMBAR = %.4f M☉\n", rbar, zmbar)
        @printf(
            io,
            "          virial ratio Q = T/|W| = %.3f (0.5 = equilibrium), RBAR/r_hm,min = %.2f\n",
            regime.q_virial,
            regime.rbar_over_rhm_min
        )
        p = regime.nbody6
        @printf(
            io,
            "Integration (merger.inp, N-body units): QE=%.1E ETAI=%.3g ETAR=%.3g NNBOPT=%d RS0=%.3g RMIN=%.3E DTMIN=%.3E KZ(16)=%d\n",
            p.qe,
            p.etai,
            p.etar,
            p.nnbopt,
            p.rs0,
            p.rmin,
            p.dtmin,
            p.kz16
        )
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
