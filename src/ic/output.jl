# =============================================================================
# dat.10 writer and .inp generator for Nbody6++ merger ICs
# =============================================================================

"""
    write_dat10(path, mass, pos, vel)

Write a `dat.10` particle file for Nbody6++ (`KZ(22)=2`).

Each line: `MASS  X  Y  Z  VX  VY  VZ`

All values must be in N-body units (G=1, M_total=1).
"""
function write_dat10(path::AbstractString, mass::Vector{Float64},
                     pos::Matrix{Float64}, vel::Matrix{Float64})
    N = length(mass)
    size(pos) == (3, N) || error("pos must be 3×$N, got $(size(pos))")
    size(vel) == (3, N) || error("vel must be 3×$N, got $(size(vel))")

    _backup_existing(path)
    open(path, "w") do io
        for i in 1:N
            @printf(io, "%.15e  %.15e  %.15e  %.15e  %.15e  %.15e  %.15e\n",
                    mass[i], pos[1, i], pos[2, i], pos[3, i],
                    vel[1, i], vel[2, i], vel[3, i])
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
- Velocity: `v_nb = v_kms / vstar` where `vstar = 0.06557 × sqrt(M_total / rbar)` km/s

The factor 0.06557 comes from `sqrt(G M☉ / pc)` in km/s.
"""
function to_nbody_units!(mass::Vector{Float64}, pos::Matrix{Float64},
                         vel::Matrix{Float64}, M_total_solar::Float64,
                         rbar_pc::Float64)
    vstar = 0.06557 * sqrt(M_total_solar / rbar_pc)
    mass ./= M_total_solar
    pos  ./= rbar_pc
    vel  ./= vstar
    return nothing
end

"""
    generate_merger_inp(path, N_total, rbar, zmbar; kz22=2, kz14=0)

Generate an Nbody6++ `.inp` file configured for external particle
input via `dat.10`.

Follows the Fortran NAMELIST read order expected by `nbody6.F → start.F`:
  `&INNBODY6` → `&ININPUT` → `&INSSE` → `&INBSE` → `&INCOLL` →
  `&INDATA` → `&INSCALE` (→ `&INXTRNL0` if KZ(14)>0).

# Key settings
- `KZ(22) = kz22`: 2 for N-body units, 10 for astrophysical units
- `KZ(14) = kz14`: tidal field option (0 = isolated)
"""
function generate_merger_inp(path::AbstractString, N_total::Int,
                             rbar::Float64, zmbar::Float64;
                             kz22::Int = 2, kz14::Int = 0,
                             tcrit::Float64 = 100.0,
                             dtadj::Float64 = 1.0,
                             deltat::Float64 = 1.0,
                             nrand::Int = 10000)
    kz = zeros(Int, 50)
    kz[1]  = 1;  kz[2]  = -1; kz[3]  = 2;  kz[7]  = 3
    kz[12] = 1;  kz[14] = kz14; kz[19] = 3; kz[22] = kz22
    kz[23] = 2;  kz[26] = 1;  kz[30] = 1

    # Build KZ lines: KZ(1:10) = ... etc.
    kz_strs = String[]
    for row in 1:5
        s = (row - 1) * 10 + 1
        push!(kz_strs, "KZ($(s):$(s+9))=$(join(string.(kz[s:s+9]), " "))")
    end

    # Neighbour number: ~sqrt(N), clamped to [20, 300]
    nnbopt = clamp(round(Int, sqrt(N_total)), 20, 300)

    _backup_existing(path)
    open(path, "w") do io
        # --- 1. &INNBODY6: start/restart, CPU time, checkpointing ---
        println(io, "&INNBODY6")
        println(io, "KSTART=1,TCOMP=1.0E8,TCRTP0=3600.0,isernb=40,iserreg=40,iserks=0 /")
        println(io)

        # --- 2. &ININPUT: main simulation parameters + KZ options ---
        println(io, "&ININPUT")
        @printf(io, "N=%d,NFIX=1,NCRIT=10,NRAND=%d,NNBOPT=%d,NRUN=1,NCOMM=10,\n",
                N_total, abs(nrand) % typemax(Int32), nnbopt)
        @printf(io, "ETAI=0.02,ETAR=0.02,RS0=0.5,DTADJ=%.4f,DELTAT=%.4f,TCRIT=%.2f,QE=1.0,RBAR=%.6f,ZMBAR=%.6f,\n",
                dtadj, deltat, tcrit, rbar, zmbar)
        println(io, join(kz_strs, "\n"))
        println(io, "DTMIN=2.5E-6,RMIN=8.E-5,ETAU=0.1,ECLOSE=1.0,GMIN=1.0E-06,GMAX=0.01,SMAX=1.0,")
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
        println(io, "&INDATA")
        println(io, "ALPHAS=2.35,BODY1=150.0,BODYN=0.08,NBIN0=0,NHI0=0,ZMET=0.001,EPOCH0=0,DTPLOT=1.0 /")
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
            println(io, "GMG=1.78E11,RG0=13.3,DISK=,A=,B=,VCIRC=,RCIRC=,GMB=,AR=,GAM=,RG=,,,VG=,,,MP=,AP2=,MPDOT=,TDELAY= /")
            println(io)
        end

        # --- No binaries (NBIN0=0) or hierarchical triples ---
    end

    @info "Wrote .inp file: $path (N=$N_total, KZ(22)=$kz22, RBAR=$rbar, ZMBAR=$zmbar)"
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
        # r_hm = 1.305 a for a Plummer sphere → scale radius from target r_hm
        sample_plummer(spec.N, spec.rbar / 1.305; rng = rng)
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
function generate_merger_ic(cfg::MergerConfig;
                            rng::Union{AbstractRNG,Nothing} = nothing,
                            output_dir::AbstractString = "")
    out_dir = isempty(output_dir) ? cfg.output.output_dir : output_dir
    mkpath(out_dir)

    n_clusters = length(cfg.clusters)
    n_clusters ≥ 2 || error("Need at least 2 clusters, got $n_clusters")

    # Resolve RNG and effective seed. When the caller supplies an RNG the
    # sampling is NOT reproducible from `effective_seed`; the metadata records
    # this honestly via `external_rng` (the seed still feeds Nbody6's NRAND).
    external_rng = rng !== nothing
    effective_seed::Int = cfg.seed === nothing ? Int(rand(UInt32) % typemax(Int32)) :
                                                 cfg.seed
    use_rng = external_rng ? rng : Random.MersenneTwister(effective_seed)

    seed_note = external_rng ? "external RNG (seed not reproducible)" : "seed=$effective_seed"
    @info "Generating merger ICs for $n_clusters clusters ($(cfg.orbit_mode) mode, $seed_note)..."

    # Sample each cluster independently
    cluster_data = [begin
        @info "  Cluster $i: $(profile_name(spec.profile)) profile, N=$(spec.N), " *
              "imf=$(imf_name(spec.imf)), M≈$(round(expected_mass(spec.imf, spec.N); digits=1)) M☉"
        _sample_cluster(spec; rng = use_rng)
    end for (i, spec) in enumerate(cfg.clusters)]

    # Combine clusters according to orbit mode
    pos_combined, vel_combined, mass_combined, cluster_ranges = if cfg.orbit_mode == "kepler"
        c1, c2 = cluster_data[1], cluster_data[2]
        p, v, m = setup_two_cluster_orbit(
            c1.pos, c1.vel, c1.mass,
            c2.pos, c2.vel, c2.mass,
            cfg.orbit.apocentre, cfg.orbit.eccentricity;
            truncate_jacobi_flag = cfg.output.truncate_jacobi
        )
        # Reconstruct ranges from the returned combined array
        # cluster 1 occupies 1:n1, cluster 2 n1+1:end
        # But truncation may have changed counts — infer from total
        # setup_two_cluster_orbit concatenates: first N1' then N2'
        # We need to know how many survived truncation. Since cluster_data
        # is NOT modified in-place by setup_two_cluster_orbit (it copies),
        # we must infer: if truncation is on, re-truncate to count.
        if cfg.output.truncate_jacobi
            M1, M2 = sum(c1.mass), sum(c2.mass)
            rJ1 = jacobi_radius(cfg.orbit.apocentre, M1, M2)
            rJ2 = jacobi_radius(cfg.orbit.apocentre, M2, M1)
            n1 = count(j -> sqrt(c1.pos[1,j]^2+c1.pos[2,j]^2+c1.pos[3,j]^2) ≤ rJ1,
                       1:length(c1.mass))
            n2 = count(j -> sqrt(c2.pos[1,j]^2+c2.pos[2,j]^2+c2.pos[3,j]^2) ≤ rJ2,
                       1:length(c2.mass))
        else
            n1 = length(c1.mass)
            n2 = length(c2.mass)
        end
        ranges = [1:n1, (n1+1):(n1+n2)]
        (p, v, m, ranges)
    elseif cfg.orbit_mode == "explicit"
        combine_clusters_explicit(cluster_data, cfg.clusters;
            truncate_jacobi_flag = cfg.output.truncate_jacobi)
    else
        error("Unknown orbit_mode: $(cfg.orbit_mode)")
    end

    N_total = length(mass_combined)
    M_total = sum(mass_combined)
    zmbar = M_total / N_total

    # Half-mass radius of the combined system — this is the RBAR length unit
    # written to the .inp file and used for the NB-unit conversion of dat.10.
    rbar = half_mass_radius(mass_combined, pos_combined)

    # Keep physical-unit copies before conversion.
    # All sampled & Kepler velocities are in "code units" where G=1 with
    # (M_sun, pc) bases, so 1 code unit = sqrt(G M_sun / pc) = 0.06557 km/s.
    mass_phys = copy(mass_combined)
    pos_phys  = copy(pos_combined)
    vel_phys  = vel_combined .* 0.06557  # code units → km/s

    # Convert to N-body units for dat.10. `to_nbody_units!` expects velocity
    # in km/s, so we must convert from code units first.
    kz22 = if cfg.output.format == "nbody"
        vel_combined .*= 0.06557  # code units → km/s (match to_nbody_units! API)
        to_nbody_units!(mass_combined, pos_combined, vel_combined, M_total, rbar)
        2
    elseif cfg.output.format == "astro"
        # In astrophysical output we write the physical-units copies so the
        # file contains km/s, pc, M_sun — convert in-place here too.
        vel_combined .*= 0.06557
        10
    else
        error("Unknown output format: $(cfg.output.format)")
    end

    # Write files
    write_dat10(joinpath(out_dir, "dat.10"), mass_combined, pos_combined, vel_combined)
    generate_merger_inp(joinpath(out_dir, "merger.inp"), N_total, rbar, zmbar;
                        kz22 = kz22, tcrit = cfg.output.tcrit,
                        dtadj = cfg.output.dtadj, deltat = cfg.output.deltat,
                        nrand = effective_seed)

    # Summary log
    _write_merger_summary(joinpath(out_dir, "merger_summary.txt"),
                          cfg, cluster_data, cluster_ranges,
                          N_total, M_total, rbar, zmbar, kz22)

    # Structured metadata for post-hoc regeneration of IC plots
    _write_merger_ic_metadata(joinpath(out_dir, "merger_ic.toml"),
                               cfg, cluster_ranges, effective_seed, external_rng,
                               N_total, M_total, rbar, zmbar)

    return MergerICResult(
        out_dir, N_total, M_total, rbar, zmbar,
        cluster_ranges, collect(cfg.clusters), cfg.orbit_mode, cfg.orbit,
        mass_phys, pos_phys, vel_phys,
    )
end

"""
    _write_merger_ic_metadata(path, cfg, cluster_ranges, seed, external_rng,
                              N_total, M_total, rbar, zmbar)

Write a machine-readable TOML snapshot of the IC generation (schema v2),
sufficient to reconstruct a [`MergerICResult`](@ref) from the `dat.10` file
later (so `plot_merger_ic` can be re-run after code fixes without
re-sampling). Cluster specs are stored in the structured form
(`profile = {type=...}`, `imf = {type=...}`) that
[`_parse_cluster_table`](@ref) also accepts, so the file round-trips.
"""
function _write_merger_ic_metadata(path::AbstractString, cfg::MergerConfig,
                                    cluster_ranges::Vector{UnitRange{Int}},
                                    seed::Int, external_rng::Bool,
                                    N_total::Int, M_total::Float64,
                                    rbar::Float64, zmbar::Float64)
    d = Dict{String,Any}(
        "meta" => Dict{String,Any}(
            "generated_at" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
            "schema_version" => 2,
            "seed"         => seed,
            "external_rng" => external_rng,
            "N_total" => N_total,
            "M_total" => M_total,
            "rbar"    => rbar,
            "zmbar"   => zmbar,
        ),
        "orbit_mode" => cfg.orbit_mode,
        "orbit" => Dict{String,Any}(
            "apocentre"    => cfg.orbit.apocentre,
            "eccentricity" => cfg.orbit.eccentricity,
        ),
        "output" => Dict{String,Any}(
            "format"          => cfg.output.format,
            "truncate_jacobi" => cfg.output.truncate_jacobi,
            "tcrit"           => cfg.output.tcrit,
            "dtadj"           => cfg.output.dtadj,
            "deltat"          => cfg.output.deltat,
        ),
        "cluster_ranges" => [[first(r), last(r)] for r in cluster_ranges],
    )

    # One table per cluster so spec + post-truncation count are co-located
    for (i, spec) in enumerate(cfg.clusters)
        d["cluster$i"] = Dict{String,Any}(
            "N"        => spec.N,
            "rbar"     => spec.rbar,
            "profile"  => _profile_table(spec.profile),
            "imf"      => _imf_table(spec.imf),
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
    dat_path  = joinpath(dir, "dat.10")
    isfile(meta_path) || error("Missing merger_ic.toml in $dir")
    isfile(dat_path)  || error("Missing dat.10 in $dir")

    raw = TOML.parsefile(meta_path)
    meta = raw["meta"]::Dict
    M_total = Float64(meta["M_total"])
    rbar    = Float64(meta["rbar"])
    zmbar   = Float64(meta["zmbar"])
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
    orbit_raw  = raw["orbit"]::Dict
    orbit = OrbitSpec(;
        apocentre    = Float64(orbit_raw["apocentre"]),
        eccentricity = Float64(orbit_raw["eccentricity"]),
    )

    # Read dat.10 — format depends on the output format. "nbody" stores NB
    # units; "astro" stores physical. We always reconstruct physical copies.
    format = String(raw["output"]["format"])
    lines = readlines(dat_path)
    N = length(lines)
    N == N_total || @warn "dat.10 row count ($N) differs from metadata N_total ($N_total)"

    mass = zeros(Float64, N)
    pos  = zeros(Float64, 3, N)
    vel  = zeros(Float64, 3, N)
    for (k, line) in enumerate(lines)
        parts = parse.(Float64, split(line))
        mass[k]  = parts[1]
        pos[:, k] = parts[2:4]
        vel[:, k] = parts[5:7]
    end

    mass_phys = if format == "nbody"
        mass .* M_total                  # NB mass fraction → M☉
    else
        mass                             # already M☉
    end
    pos_phys = if format == "nbody"
        pos .* rbar                      # NB length → pc
    else
        pos                              # already pc
    end
    vel_phys = if format == "nbody"
        # NB velocity → km/s: v_nb × vstar_kms, vstar = 0.06557 √(M_tot/rbar)
        vstar_kms = 0.06557 * sqrt(M_total / rbar)
        vel .* vstar_kms
    else
        vel                              # already km/s
    end

    return MergerICResult(
        abspath(dir), N_total, M_total, rbar, zmbar,
        cluster_ranges, cluster_specs, orbit_mode, orbit,
        mass_phys, pos_phys, vel_phys,
    )
end

function _write_merger_summary(path, cfg, cluster_data, cluster_ranges,
                               N_total, M_total, rbar, zmbar, kz22)
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
            @printf(io, "  Cluster %d: %s, imf=%s, N=%d (after trunc: %d), M=%.1f M☉, r_hm=%.2f pc",
                    i, profile_name(spec.profile), imf_name(spec.imf),
                    spec.N, Ni, M_sampled, spec.rbar)
            if spec.profile isa KingProfile
                @printf(io, ", W0=%.1f", spec.profile.W0)
            end
            if !isempty(spec.position)
                @printf(io, "\n             pos=[%.2f, %.2f, %.2f] pc",
                        spec.position[1], spec.position[2], spec.position[3])
                @printf(io, "  vel=[%.2f, %.2f, %.2f]",
                        spec.velocity[1], spec.velocity[2], spec.velocity[3])
            end
            println(io)
        end
        println(io)
        if cfg.orbit_mode == "kepler"
            @printf(io, "Orbit: d_apo = %.2f pc, e = %.3f\n",
                    cfg.orbit.apocentre, cfg.orbit.eccentricity)
            a_orb = cfg.orbit.apocentre / (1.0 + cfg.orbit.eccentricity)
            @printf(io, "       a = %.2f pc (semi-major axis)\n", a_orb)
            println(io)
        end
        @printf(io, "Combined: N_total = %d, M_total = %.1f M☉\n", N_total, M_total)
        @printf(io, "          RBAR = %.4f pc, ZMBAR = %.4f M☉\n", rbar, zmbar)
        println(io)
        println(io, "Output format: ", cfg.output.format, " (KZ(22)=$kz22)")
        println(io, "Jacobi truncation: ", cfg.output.truncate_jacobi)
        println(io)
        println(io, "Files:")
        println(io, "  dat.10     — particle data")
        println(io, "  merger.inp — Nbody6++ input file")
        println(io, "=" ^ 60)
    end
    @info "Summary written to $path"
end
