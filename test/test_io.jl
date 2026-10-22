# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Fortran binary I/O" begin
    # Create a synthetic Fortran binary file
    fpath = joinpath(TESTDIR, "test_fortran.bin")
    open(fpath, "w") do io
        # Write a record: 3 Float32 values
        data = Float32[1.0, 2.0, 3.0]
        marker = Int32(sizeof(data))
        write(io, marker)
        write(io, data)
        write(io, marker)

        # Write a second record: 2 Int32 values
        idata = Int32[42, 99]
        marker2 = Int32(sizeof(idata))
        write(io, marker2)
        write(io, idata)
        write(io, marker2)
    end

    open(fpath, "r") do io
        rec1 = Nbody6Dynamics.read_fortran_record(io, Float32, 3)
        @test rec1 ≈ Float32[1.0, 2.0, 3.0]

        rec2 = Nbody6Dynamics.read_fortran_record(io, Int32, 2)
        @test rec2 == Int32[42, 99]
    end

    # peek_record_size
    open(fpath, "r") do io
        sz = Nbody6Dynamics.peek_record_size(io)
        @test sz == Int32(12)  # 3 × 4 bytes
    end
end

# =====================================================================

@testset "conf.3 reader" begin
    using StableRNGs
    rng = StableRNG(63)
    # Create a synthetic conf.3 file in standard format
    fpath = joinpath(TESTDIR, "conf.3_test")
    n_particles = 5

    open(fpath, "w") do io
        # Record 1: header integers
        hdr = Int32[n_particles, 1, 1, 20]
        _write_fortran_record(io, hdr)

        # Record 2: AS(1:20) parameters
        params = zeros(Float32, 20)
        params[1] = 0.5f0    # time
        params[3] = 1.0f0    # rbar
        params[4] = 0.5f0    # zmbar
        params[11] = 10.0f0  # tscale
        params[12] = 5.0f0   # vstar
        params[18] = 2.0f0   # rscale
        _write_fortran_record(io, params)

        # Particle records (standard: M, X, Y, Z, VX, VY, VZ, NAME)
        for i in 1:n_particles
            m = Float32(1.0 / n_particles)
            x, y, z = Float32.(randn(rng, 3))
            vx, vy, vz = Float32.(0.1 .* randn(rng, 3))
            name = Int32(i)
            buf = IOBuffer()
            write(buf, m, x, y, z, vx, vy, vz, name)
            _write_fortran_record(io, take!(buf))
        end
    end

    snap = read_conf3(fpath)
    @test nparticles(snap) == n_particles
    @test time_nb(snap.header) ≈ 0.5
    @test tscale(snap.header) ≈ 10.0
    @test rscale(snap.header) ≈ 2.0
    @test length(snap.mass) == n_particles
    @test size(snap.pos) == (3, n_particles)
    @test size(snap.vel) == (3, n_particles)
    @test snap.name == Int32.(1:n_particles)
    @test isempty(snap.rho)  # standard format has no density

    # Test read_all_conf3
    snaps = read_all_conf3(TESTDIR, "conf.3_")
    @test length(snaps) >= 1
end

# =====================================================================

@testset "conf.3 reader — extended format" begin
    using StableRNGs
    rng = StableRNG(108)
    fpath = joinpath(TESTDIR, "conf.3_ext")
    n_particles = 3

    open(fpath, "w") do io
        hdr = Int32[n_particles, 1, 1, 20]
        _write_fortran_record(io, hdr)

        params = zeros(Float32, 20)
        params[1] = 1.0f0
        _write_fortran_record(io, params)

        # Extended: M, RHO, XNS, X, Y, Z, VX, VY, VZ, PHI, NAME
        for i in 1:n_particles
            buf = IOBuffer()
            write(buf, Float32(0.1 * i))   # mass
            write(buf, Float32(100.0))      # rho
            write(buf, Float32(0.01))       # xns
            write(buf, Float32.(randn(rng, 3))...)  # pos
            write(buf, Float32.(0.1 .* randn(rng, 3))...)  # vel
            write(buf, Float32(-0.5))       # phi
            write(buf, Int32(i))            # name
            _write_fortran_record(io, take!(buf))
        end
    end

    snap = read_conf3(fpath)
    @test nparticles(snap) == n_particles
    @test length(snap.rho) == n_particles
    @test length(snap.phi) == n_particles
    @test all(snap.rho .≈ 100.0f0)
    @test all(snap.phi .≈ -0.5f0)
end

# =====================================================================

@testset "conf.3 reader — pair centre-of-mass entries" begin
    using StableRNGs
    # Engine layout: NTOT = N + NPAIRS, components of pair k at 2k − 1 and 2k,
    # its centre of mass at N + k (output.F).
    rng = StableRNG(11)
    n_stars, n_pairs = 7, 2
    ntot = n_stars + n_pairs
    mass = Float32.(0.05 .+ rand(rng, n_stars))
    pos = randn(rng, Float32, 3, n_stars)
    vel = randn(rng, Float32, 3, n_stars)
    cm_mass = [mass[2k - 1] + mass[2k] for k in 1:n_pairs]
    cm_pos = reduce(
        hcat,
        [
            (mass[2k - 1] .* pos[:, 2k - 1] .+ mass[2k] .* pos[:, 2k]) ./ cm_mass[k] for
            k in 1:n_pairs
        ],
    )
    cm_vel = reduce(
        hcat,
        [
            (mass[2k - 1] .* vel[:, 2k - 1] .+ mass[2k] .* vel[:, 2k]) ./ cm_mass[k] for
            k in 1:n_pairs
        ],
    )
    names = vcat(Int32.(1:n_stars), Int32.(1000 .+ (1:n_pairs)))
    function write_bulk(path, npairs_slot)
        open(path, "w") do io
            _write_fortran_record(io, Int32[ntot, 1, 1, 20])
            params = zeros(Float32, 20)
            params[1] = 2.0f0
            params[2] = npairs_slot
            body = vcat(
                params,
                mass,
                cm_mass,
                fill(1.0f0, ntot),                     # RHOS
                zeros(Float32, ntot),                  # XNS
                vec(hcat(pos, cm_pos)),
                vec(hcat(vel, cm_vel)),
                fill(-0.5f0, ntot),                    # PHI
                reinterpret(Float32, names),
            )
            _write_fortran_record(io, body)
        end
        return path
    end
    snap = read_conf3(write_bulk(joinpath(TESTDIR, "conf3_pairs"), Float32(n_pairs)))
    @test nparticles(snap) == n_stars
    @test snap.npairs == n_pairs
    @test snap.header.ntot == ntot
    @test pair_indices(snap) == [(1, 2), (3, 4)]
    @test snap.mass == mass
    @test snap.pos == pos && snap.vel == vel
    @test snap.name == Int32.(1:n_stars)
    @test snap.pair_name == Int32[1001, 1002]
    @test length(snap.rho) == n_stars && length(snap.phi) == n_stars
    @test sum(Float64, snap.mass) ≈ sum(Float64, mass)

    # An active chain's centre-of-mass body (NAME = 0) among the singles is dropped too.
    names[6] = 0
    chain = read_conf3(write_bulk(joinpath(TESTDIR, "conf3_chain"), Float32(n_pairs)))
    @test nparticles(chain) == n_stars - 1
    @test chain.name == Int32[1, 2, 3, 4, 5, 7]
    @test chain.mass == mass[[1, 2, 3, 4, 5, 7]]
    @test chain.pos == pos[:, [1, 2, 3, 4, 5, 7]]
    @test pair_indices(chain) == [(1, 2), (3, 4)]
    names[6] = 6

    # Inconsistent pair counts are refused, not truncated silently.
    @test_throws ArgumentError read_conf3(write_bulk(joinpath(TESTDIR, "conf3_pairs_bad"), 5.0f0))
    @test_throws ArgumentError read_conf3(write_bulk(joinpath(TESTDIR, "conf3_pairs_frac"), 1.5f0))
    @test_throws ArgumentError Snapshot(
        snap.header,
        snap.name,
        snap.mass,
        snap.pos,
        snap.vel,
        Float32[],
        Float32[],
        4,
    )
    @test_throws DimensionMismatch Snapshot(
        snap.header,
        snap.name,
        snap.mass[1:3],
        snap.pos,
        snap.vel,
        Float32[],
        Float32[],
    )
    @test_throws DimensionMismatch Snapshot(
        snap.header,
        snap.name,
        snap.mass,
        snap.pos,
        snap.vel,
        Float32[],
        Float32[],
        2,
        Int32[5],
    )
end

# =====================================================================

@testset "Diagnostics parser" begin
    diag_path = joinpath(TESTDIR, "out1000")
    write(
        diag_path,
        """
 Some header text...

 PHYSICAL SCALING:  R* = 1.234  M* = 5678.0  V* = 3.456  T* = 7.890
                <M> = 0.567  SU = 1.0  AU = 1.0

 ADJUST:    0.0000      0.00  1.000  0.00E+00 -0.2500  10000    500   1.234
 ADJUST:    0.5000     50.00  0.987  1.23E-06 -0.2499   9998    498   1.245
 ADJUST:    1.0000    100.00  0.995  2.45E-06 -0.2498   9990    495   1.267

 END RUN
""",
    )

    diag = read_diagnostics(diag_path)

    @test length(diag.adjust) == 3
    @test diag.adjust[1].time_nb ≈ 0.0
    @test diag.adjust[2].time_myr ≈ 50.0
    @test diag.adjust[3].n == 9990
    @test diag.adjust[2].npairs == 498

    # Physical scaling
    @test haskey(diag.physical_scaling, "R*")
    @test diag.physical_scaling["R*"] ≈ 1.234
    @test diag.physical_scaling["V*"] ≈ 3.456

    # Unit extraction
    units = extract_scaling(diag)
    @test units.rbar ≈ 1.234
    @test units.tscale ≈ 7.890
end

# =====================================================================

@testset "Lagrangian radii reader" begin
    lagr_path = joinpath(TESTDIR, "lagr.7")

    # Synthetic lagr.7 in the upstream lagr.f layout: ##/TIME header,
    # then one row per epoch (3 epochs, 5 radii columns, simplified)
    write(
        lagr_path,
        """
## Time, [Column number | Label]:    2 R_{lagr}   (synthetic)
TIME   0.01  0.05  0.20  0.50  1.00
0.0  0.01  0.05  0.20  0.50  1.00
0.5  0.012 0.055 0.22  0.52  1.05
1.0  0.015 0.060 0.25  0.55  1.10
""",
    )

    lagr = read_lagr(lagr_path)

    @test length(lagr.time) == 3
    @test lagr.time ≈ [0.0, 0.5, 1.0]
    @test size(lagr.radii) == (5, 3)
    @test lagr.radii[1, 1] ≈ 0.01
    @test lagr.radii[5, 3] ≈ 1.10
end

# =====================================================================

@testset "Escaper reader" begin
    # Real esc.11 shape (escape.F): 5 NB-unit diagnostics, then the
    # physical-unit block T[Myr] M[M*] EESC VI[km/s] K* NAME, then the
    # direction angles ANGLE PHI / ANGLE THETA (tokens 12–13).  The last
    # line is truncated to 11 tokens (no angles) and must be skipped.
    esc_path = joinpath(TESTDIR, "esc.11")
    write(
        esc_path,
        """
     TTOT         BODY         RI           VI           STEP         T[Myr]       M[M*]      EESC      VI[km/s]     K*  NAME      ANGLE PHI     ANGLE THETA
   1.00000E+00  5.00000E-04  2.00000E+01  4.00000E+01  1.95312E-03  1.23400E+00  5.00000E-01 -1.23000E-01  1.56000E+01   0       101  1.20000E+01  2.50000E+01
   2.00000E+00  3.00000E-04  2.50000E+01  4.20000E+01  1.95312E-03  2.56700E+00  3.00000E-01  4.56000E-01  2.23000E+01   1       202  1.85000E+02 -3.00000E+01
   3.00000E+00  1.20000E-03  3.00000E+01  4.40000E+01  1.95312E-03  3.89000E+00  1.20000E+00 -7.89000E-01  1.01000E+01  14       303  3.40000E+02  6.00000E+01
   4.00000E+00  1.00000E-03  3.10000E+01  4.50000E+01  1.95312E-03  4.50000E+00  1.00000E+00  1.00000E-01  1.20000E+01   1       404
""",
    )

    escs = read_escapers(esc_path)
    @test length(escs) == 3   # the 11-token line is skipped

    @test escs[1].time_myr ≈ 1.234
    @test escs[1].mass_solar ≈ 0.500
    @test escs[1].escape_energy ≈ -0.123
    @test escs[1].velocity_kms ≈ 15.6
    @test escs[1].stellar_type == 0
    @test escs[1].name == 101
    @test escs[1].phi_deg ≈ 12.0
    @test escs[1].theta_deg ≈ 25.0

    @test escs[2].phi_deg ≈ 185.0
    @test escs[2].theta_deg ≈ -30.0

    @test escs[3].stellar_type == 14   # BH (Hurley convention)
    @test escs[3].name == 303
    @test escs[3].phi_deg ≈ 340.0
    @test escs[3].theta_deg ≈ 60.0

    # Empty file
    empty_path = joinpath(TESTDIR, "esc_empty.11")
    write(empty_path, "")
    @test isempty(read_escapers(empty_path))
end

# =====================================================================

@testset "Stellar evolution reader" begin
    # Header time is TPHYS [Myr]; data-line token 1 is TTOT [NB] —
    # different clocks, both kept. 15-token v2026.07+ layout.
    sev_path = joinpath(TESTDIR, "sev.83_0")
    write(
        sev_path,
        """
  3  0.5000
   0.4315   1   101  0  1.20  0.800   0.123  -0.456  3.750  0.0  0.0  90.0  0.0  0.0  1.0e-10
   0.4315   2   202  1  0.80  1.200   1.500   0.200  4.100  0.0  0.0  5.20  0.0  0.0  1.0e-10
   0.4315   3   303  13 2.50  10.00   5.000   1.500  4.500  0.0  0.0  4.10  1.4  1.0e-5  1.0e-10
""",
    )

    sev = read_stellar_evolution(sev_path)
    @test sev.n_stars == 3
    @test sev.time_myr ≈ 0.5
    @test length(sev.records) == 3

    r1 = sev.records[1]
    @test r1.time_nb ≈ 0.4315   # per-line TTOT [NB], not the Myr header
    @test r1.index == Int32(1)
    @test r1.name == Int32(101)
    @test r1.stellar_type == Int32(0)  # MS
    @test r1.ri ≈ 1.20
    @test r1.mass_solar ≈ 0.800
    @test r1.log_luminosity ≈ 0.123
    @test r1.log_radius ≈ -0.456
    @test r1.log_teff ≈ 3.750
    @test r1.ms_lifetime_myr ≈ 90.0    # TM

    r3 = sev.records[3]
    @test r3.stellar_type == Int32(13)  # NS
    @test r3.mass_core ≈ 1.4

    # Test read_all_stellar_evolution
    sev2_path = joinpath(TESTDIR, "sev.83_1")
    write(
        sev2_path,
        """
  2  1.0000
   1.0000   1   101  2  1.50  0.750   0.500  -0.200  3.600  0.0  0.0  8.0  0.2  0.01  2.0
   1.0000   2   202  4  0.90  1.100   2.000   0.400  3.900  0.0  0.0  6.5  0.4  0.02  4.0
""",
    )

    sevs = read_all_stellar_evolution(TESTDIR, "sev.83_*")
    @test length(sevs) == 2
    @test sevs[1].time_myr < sevs[2].time_myr

    # Older 11-token layout (the one the engine manual documents): the line
    # ends at AGE EPOCH and the SSE tail is unknown, not absent data.
    # Outside TESTDIR: the sev.83_ count above must stay at two.
    alt_dir = mktempdir()
    short_path = joinpath(alt_dir, "sev.83_short")
    write(
        short_path,
        """
  2  0.5000
   0.4315   1   101  0  1.20  0.800   0.123  -0.456  3.750  0.0  0.0
   0.4315   2   202  1  0.80  1.200   1.500   0.200  4.100  0.0  0.0
""",
    )
    short = @test_logs read_stellar_evolution(short_path)   # a known layout: no warning
    @test length(short.records) == 2
    @test all(r -> isnan(r.ms_lifetime_myr), short.records)
    @test all(
        r -> isnan(r.mass_core) && isnan(r.radius_core) && isnan(r.radius_envelope),
        short.records,
    )
    @test short.records[1].log_teff ≈ 3.750

    # A line of neither layout is skipped, and the file reports its tally once.
    mixed_path = joinpath(alt_dir, "sev.83_mixed")
    write(
        mixed_path,
        """
  2  0.5000
   0.4315   1   101  0  1.20  0.800   0.123  -0.456  3.750  0.0  0.0  90.0  0.0  0.0  1.0e-10
   0.4315   2   202  1  0.80  1.200   1.500   0.200  4.100  0.0  0.0  5.20  0.0
""",
    )
    mixed = @test_logs (:warn, r"1 of 2 data lines skipped") read_stellar_evolution(mixed_path)
    @test length(mixed.records) == 1
    @test mixed.records[1].ms_lifetime_myr ≈ 90.0
end

# =====================================================================

@testset "Binary evolution reader" begin
    # hrplot.F FORMAT 5: 32 tokens per line; header NPAIRS TPHYS [Myr];
    # token 1 of every data line is TTOT [NB]. Record 1 is a verbatim
    # engine line, record 2 carries a recognisable 1..12 tail.
    bev_path = joinpath(TESTDIR, "bev.82_0")
    write(
        bev_path,
        """
   2      0.0
   0.00000E+00       1       2       1       2  0  1   0  1.34446E+00  9.21880E-01  1.20647E+00  1.30277E+00  2.21134E-01  1.98278E-01 -1.99897E+00 -2.08678E+00 -6.28737E-01 -6.64486E-01  3.57640E+00  3.57233E+00  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  4.98658E+05  6.20908E+05  0.00000E+00  0.00000E+00  0.00000E+00  0.00000E+00  1.52819E-01  1.40743E-01
   0.00000E+00       3       4       3       4  1  0   0  5.14787E-01  6.97120E-01  5.58338E+00  4.42220E+00  1.11885E+00  5.68179E-01  3.81907E-01 -9.85920E-01 -1.23430E-02 -2.97930E-01  3.86343E+00  3.70000E+00  1.0  2.0  3.0  4.0  5.0  6.0  7.0  8.0  9.0  10.0  11.0  12.0
   0.00000E+00       5       6
""",
    )
    bev = read_binary_evolution(bev_path)
    @test bev.n_pairs == 2
    @test bev.time_myr == 0.0
    @test length(bev.records) == 2        # the short line is skipped

    r1 = bev.records[1]
    @test r1.time_nb == 0.0
    @test (r1.index1, r1.index2, r1.name1, r1.name2) == (1, 2, 1, 2)
    @test (r1.stellar_type1, r1.stellar_type2, r1.stellar_type_cm) == (0, 1, 0)
    @test r1.ri ≈ 1.34446
    @test r1.eccentricity ≈ 0.92188
    @test r1.log_period_days ≈ 1.20647
    @test r1.log_semi_major_axis_rsun ≈ 1.30277
    @test r1.mass1 ≈ 0.221134
    @test r1.mass2 ≈ 0.198278
    @test r1.log_luminosity1 ≈ -1.99897
    @test r1.log_radius2 ≈ -0.664486
    @test r1.log_teff1 ≈ 3.5764
    @test r1.ms_lifetime1_myr ≈ 4.98658e5     # TM of a 0.22 M☉ dwarf
    @test r1.radius_envelope2 ≈ 0.140743

    r2 = bev.records[2]
    @test (r2.age1_myr, r2.age2_myr) == (1.0, 2.0)
    @test (r2.epoch1_myr, r2.epoch2_myr) == (3.0, 4.0)
    @test (r2.ms_lifetime1_myr, r2.ms_lifetime2_myr) == (5.0, 6.0)
    @test (r2.mass_core1, r2.mass_core2) == (7.0, 8.0)
    @test (r2.radius_core1, r2.radius_core2) == (9.0, 10.0)
    @test (r2.radius_envelope1, r2.radius_envelope2) == (11.0, 12.0)

    # The engine's period and semi-major axis columns obey Kepler III,
    # a³/P² = m₁ + m₂ in au, yr, M☉ — a check of the column map and of
    # the R☉ → pc → au conversions.
    P_yr = 10^r1.log_period_days / 365.25
    a_au = semi_major_axis_pc(r1) / Nbody6Dynamics._AU_IN_PC
    @test a_au^3 / P_yr^2 ≈ r1.mass1 + r1.mass2 rtol = 1e-2
    @test binding_energy(r1) ≈
          Nbody6Dynamics._G_PC_KMS2_MSUN * r1.mass1 * r1.mass2 / (2 * semi_major_axis_pc(r1))

    write(joinpath(TESTDIR, "bev.82_1"), "       0      1.5\n")
    bevs = read_all_binary_evolution(TESTDIR, "bev.82_*")
    @test length(bevs) == 2
    @test bevs[1].time_myr < bevs[2].time_myr
    @test bevs[2].n_pairs == 0 && isempty(bevs[2].records)
    @test isempty(read_all_binary_evolution(joinpath(TESTDIR, "absent"), "bev.82_*"))
    @test_throws ErrorException read_binary_evolution(joinpath(TESTDIR, "missing.82"))

    # Older 24-token layout (the engine manual): the line ends at
    # EPOCH(I1) EPOCH(I2), so the eight SSE fields are unknown.
    alt_dir = mktempdir()
    short_path = joinpath(alt_dir, "bev.82_short")
    write(
        short_path,
        """
   1      0.0
   0.00000E+00       3       4       3       4  1  0   0  5.14787E-01  6.97120E-01  5.58338E+00  4.42220E+00  1.11885E+00  5.68179E-01  3.81907E-01 -9.85920E-01 -1.23430E-02 -2.97930E-01  3.86343E+00  3.70000E+00  1.0  2.0  3.0  4.0
""",
    )
    short = @test_logs read_binary_evolution(short_path)   # a known layout: no warning
    @test length(short.records) == 1
    s1 = short.records[1]
    @test (s1.index1, s1.index2) == (3, 4)
    @test s1.eccentricity ≈ 0.69712
    @test (s1.age1_myr, s1.age2_myr) == (1.0, 2.0)
    @test (s1.epoch1_myr, s1.epoch2_myr) == (3.0, 4.0)
    @test all(
        isnan,
        (
            s1.ms_lifetime1_myr,
            s1.ms_lifetime2_myr,
            s1.mass_core1,
            s1.mass_core2,
            s1.radius_core1,
            s1.radius_core2,
            s1.radius_envelope1,
            s1.radius_envelope2,
        ),
    )

    # A line of neither layout is skipped, and the file reports its tally once.
    mixed_path = joinpath(alt_dir, "bev.82_mixed")
    write(
        mixed_path,
        """
   1      0.0
   0.00000E+00       3       4       3       4  1  0   0  5.14787E-01  6.97120E-01  5.58338E+00  4.42220E+00  1.11885E+00  5.68179E-01  3.81907E-01 -9.85920E-01 -1.23430E-02 -2.97930E-01  3.86343E+00  3.70000E+00  1.0  2.0  3.0  4.0
   0.00000E+00       5       6
""",
    )
    mixed = @test_logs (:warn, r"data lines skipped") read_binary_evolution(mixed_path)
    @test length(mixed.records) == 1
end

# =====================================================================

@testset "Real-output fixtures" begin
    FIXDIR = joinpath(@__DIR__, "fixtures")

    @testset "esc.11" begin
        escs = read_escapers(joinpath(FIXDIR, "esc.11"))
        @test length(escs) == 88          # every data line parses
        e1 = escs[1]
        @test e1.time_myr ≈ 8.11195 rtol = 1e-5     # T[Myr], not TTOT
        @test e1.mass_solar ≈ 18.1687 rtol = 1e-5   # M[M*], not NB mass
        @test e1.escape_energy ≈ 1017.22 rtol = 1e-5
        @test e1.velocity_kms ≈ 209.841 rtol = 1e-5
        @test e1.stellar_type == 14                 # BH escaper
        @test e1.name == 2248
        @test e1.phi_deg ≈ 22.2302 rtol = 1e-5      # ANGLE PHI [deg]
        @test e1.theta_deg ≈ 28.7010 rtol = 1e-5    # ANGLE THETA [deg]
        # escape.F angle conventions: φ ∈ [0, 360], θ ∈ [-90, 90]
        @test all(0.0 ≤ e.phi_deg ≤ 360.0 for e in escs)
        @test all(-90.0 ≤ e.theta_deg ≤ 90.0 for e in escs)
    end

    @testset "out1000 scaling + ADJUST" begin
        diag = read_diagnostics(joinpath(FIXDIR, "out1000"))
        u = extract_scaling(diag)
        @test u.rbar ≈ 5.503639 rtol = 1e-6
        @test u.zmbar ≈ 27695.934588 rtol = 1e-6   # M* (scale), not <M>
        @test u.tscale ≈ 1.15885054 rtol = 1e-6
        @test u.vstar ≈ 4.65224425 rtol = 1e-6
        @test length(diag.adjust) ≥ 6
        a0 = diag.adjust[1]
        @test a0.time_nb == 0.0
        @test a0.qvir ≈ 0.321 atol = 1e-3           # Q = T/|W|
        @test a0.e_tot ≈ -0.5408 atol = 1e-4

        # Energy bookkeeping fields of the seven fixture adjustments
        adj = diag.adjust[1:7]
        @test [a.de_rel for a in adj] == [0.0, 2.644e-6, -1.163e-6, -2.627e-6, -2.978e-5, -2.715e-5, 3.858e-4]
        @test [a.de_abs for a in adj] == [0.0, 1.424e-6, -6.220e-7, -1.388e-6, -2.251e-5, -3.203e-5, 1.821e-4]
        @test [a.detot for a in adj] == [0.0, 1.424e-6, 8.018e-7, -5.860e-7, -2.310e-5, -5.512e-5, 1.269e-4]
        cum = cumulative_energy_error(diag)
        @test cum.time_nb[1:7] == [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
        @test cum.errtot[7] ≈ 3.27724e-4 rtol = 1e-12
        @test cum.detot[1:7] == [a.detot for a in adj]
        # DETOT is the running sum of DELTA to the printed precision
        @test maximum(abs, cumsum([a.de_abs for a in adj]) .- cum.detot[1:7]) < 1e-7
        # Records without DELTA/DETOT: accumulated from de_abs, else NaN
        strip_detot(a) = AdjustRecord(
            a.time_nb,
            a.time_myr,
            a.qvir,
            a.de_rel,
            a.e_tot,
            a.n,
            a.npairs,
            a.rscale,
            a.de_abs,
            NaN,
        )
        from_abs =
            cumulative_energy_error(DiagnosticsData(strip_detot.(adj), diag.physical_scaling))
        @test from_abs.detot ≈ cumsum([a.de_abs for a in adj]) rtol = 1e-12
        bare(a) = AdjustRecord(
            a.time_nb,
            a.time_myr,
            a.qvir,
            a.de_rel,
            a.e_tot,
            a.n,
            a.npairs,
            a.rscale,
            NaN,
            NaN,
        )
        @test all(
            isnan,
            cumulative_energy_error(DiagnosticsData(bare.(adj), diag.physical_scaling)).detot,
        )
        empty_cum = cumulative_energy_error(DiagnosticsData(AdjustRecord[], diag.physical_scaling))
        @test isempty(empty_cum.time_nb) && isempty(empty_cum.errtot) && isempty(empty_cum.detot)

        # END RUN totals: absent in the fixture excerpt, parsed from the engine's line format
        @test read_energy_totals(joinpath(FIXDIR, "out1000")) === nothing
        with_end = joinpath(mktempdir(), "out1000")
        write(
            with_end,
            read(joinpath(FIXDIR, "out1000"), String) *
            "\n         END RUN    TIME[Myr] =   10.00  TOFF/TIME/TTOT=      0.00000000      8.00000000      8.00000000  CPUTOT =    0.2  ERRTOT = 1.25000D-04  DETOT = 3.00000D-05\n" *
            "\n         END RUN    TIME[Myr] =  261.31  TOFF/TIME/TTOT=      0.00000000     40.00000000     40.00000000  CPUTOT =    1.0  ERRTOT =-4.70618D-03  DETOT = 7.54813D-03\n",
        )
        tot = read_energy_totals(with_end)
        @test tot == (time_nb = 40.0, errtot = -4.70618e-3, detot = 7.54813e-3)
        @test_throws ErrorException read_energy_totals(joinpath(mktempdir(), "absent"))

        # Energy terms of the key-value line, the options echo, no tidal echo
        a1 = diag.adjust[2]
        @test all(
            isfinite,
            (a1.ekin, a1.pot, a1.etide, a1.ebin, a1.emerge, a1.esub, a1.ecoll, a1.emdot, a1.ecdot),
        )
        @test a1.ekin > 0 && a1.pot > 0 && a1.etide == 0.0
        @test a1.ekin - a1.pot ≈
              a1.e_tot - (a1.ebin + a1.emerge + a1.esub + a1.ecoll + a1.emdot + a1.ecdot) atol =
            2.0e-3
        @test length(diag.kz) == 50
        @test diag.kz[14] == 0 && diag.kz[19] == 3 && diag.kz[30] == 1 && diag.kz[2] == -1
        @test diag.tidal == (0.0, 0.0, 0.0, 0.0)
        @test tidal_mode(diag) == :none
        @test diag.nbin0 == 0 && diag.nzero > 0
        @test all(a -> a.nmerge == 0, diag.adjust)
        # A tidal echo and an options row that is not 50 integers
        tidal_log = joinpath(mktempdir(), "out1000")
        write(
            tidal_log,
            "            OPTIONS\n\n           1  2  3\n\n           1 -1  2\n\n" *
            "             TIDAL PARAMETERS:   5.34429E-03  0.00000E+00 -1.78143E-03  8.44140E-02  TSCALE =  1.56233E+00  RTIDE =  5.71866E+00\n" *
            " ADJUST:  TIME    0.00000E+00  T[Myr]   0.000E+00  Q   0.428E+00  DE   0.000E+00 DELTA   0.000E+00 DETOT   0.000E+00 E  -4.547E-01 EKIN   3.405E-01 POT   7.952E-01 ETIDE   0.000E+00 ETOT  -4.547E-01 EBIN   0.000E+00 EMERGE   0.000E+00 ESUB   0.000E+00 ECOLL   0.000E+00 EMDOT   0.000E+00 ECDOT   0.000E+00\n",
        )
        td = read_diagnostics(tidal_log)
        @test isempty(td.kz) && td.nzero == -1
        @test td.tidal == (5.34429e-3, 0.0, -1.78143e-3, 8.44140e-2)
        @test td.adjust[1].ekin == 0.3405 && td.adjust[1].pot == 0.7952
    end

    @testset "TIME[NB] NMERGE, the N and NBIN0 echoes" begin
        capture = joinpath(mktempdir(), "out1000")
        write(
            capture,
            "\n               N  NFIX  NCRIT  NRAND  NNBOPT  NRUN NCOMM\n\n" *
            "            2000     1     10   1612     100     1     1\n\n" *
            "           ALPHAS =  2.35000E+00  BODY1 =  1.00000E+00  BODYN =  8.00000E-02 ZMASS =  2.00000E+03 NBIN0 =      12 NHI0 =    0 ZMET =  2.00000E-04\n" *
            " ADJUST:  TIME    0.00000E+00  T[Myr]   0.000E+00  Q   0.362E+00  DE   0.000E+00 DELTA   0.000E+00 DETOT   0.000E+00 E  -2.934E-01\n" *
            " TIME[NB]    0.000000E+00  N     2000  <NB>   162  NPAIRS      15  NMERGE     3 MULT    0 NS     1970\n" *
            " ADJUST:  TIME    5.00000E-01  T[Myr]   0.100E+00  Q   0.500E+00  DE   1.000E-06 DELTA   1.000E-07 DETOT   1.000E-07 E  -2.934E-01\n",
        )
        d = read_diagnostics(capture)
        @test d.nbin0 == 12 && d.nzero == 2000
        @test length(d.adjust) == 2
        @test d.adjust[1].npairs == 15 && d.adjust[1].nmerge == 3
        # Forward fill of the record without a TIME[NB] line
        @test d.adjust[2].npairs == 15 && d.adjust[2].nmerge == 3
    end

    @testset "out1000 event lines" begin
        events = read_engine_events(joinpath(FIXDIR, "out1000"))
        kinds = [e.kind for e in events]
        @test count(==(:ks_rectification), kinds) == 30
        @test count(==(:kick_full), kinds) == 5
        @test count(==(:kick_zero), kinds) == 1
        @test count(==(:remnant_formation), kinds) == 1
        @test count(==(:escape), kinds) == 1
        @test count(==(:supernova_no_remnant), kinds) == 5
        @test length(events) == 43
        @test issorted([e.line for e in events]) && allunique([e.line for e in events])
        @test all(1 ≤ e.interval ≤ 8 for e in events)
        # Positions: seven ADJUST lines; the kicks of the last interval, the
        # escape line printed after the seventh ADJUST belongs to the eighth.
        esc = only(filter(e -> e.kind == :escape, events))
        @test esc.interval == 8 && esc.line == 365 && esc.time_nb == 3.0
        @test esc.count == 5 && esc.n_total == 3943
        @test esc.mass_solar == 2.48116e4 && esc.energy_nb == -0.540706
        @test esc.kstar1 == KSTAR_ABSENT && esc.name1 == 0 && isnan(esc.velocity_kms)
        form = only(filter(e -> e.kind == :remnant_formation, events))
        @test form.interval == 7 && form.line == 354
        @test form.name1 == 1598 && form.kstar1 == 5 && form.kstar2 == 14
        @test form.time_nb == 2.97986 && form.mass_solar > 0 && form.dmass_solar > 0
        zk = only(filter(e -> e.kind == :kick_zero, events))
        @test zk.line == 355 && zk.name1 == 1598 && (zk.kstar1, zk.kstar2, zk.kstar3) == (5, 14, 0)
        @test zk.velocity_kms == 0.0 && zk.time_nb == 2.979858398
        fk = filter(e -> e.kind == :kick_full, events)
        @test all(e -> e.interval == 7 && e.kstar2 == 15 && e.velocity_kms > 0, fk)
        @test [e.name1 for e in fk] == [1699, 341, 1933, 3609, 637]
        sn = filter(e -> e.kind == :supernova_no_remnant, events)
        @test [e.name1 for e in sn] == [1699, 341, 1933, 3609, 637]
        @test all(
            e -> e.kstar1 == 5 && e.kstar2 == 15 && e.mass_solar == 0.0 && e.dmass_solar > 100,
            sn,
        )
        @test sn[1].time_nb == 2.85193 && sn[1].dmass_solar == 233.833
        ks = filter(e -> e.kind == :ks_rectification, events)
        @test all(
            e ->
                e.name1 > 0 &&
                e.name2 > 0 &&
                e.name3 > e.name2 &&
                isfinite(e.energy_nb) &&
                e.mass_solar > 0,
            ks,
        )
        @test all(e -> e.kstar1 != KSTAR_ABSENT && e.kstar3 != KSTAR_ABSENT, ks)
        # The remaining writers on their real line layouts
        other = joinpath(mktempdir(), "out1000")
        write(
            other,
            " ADJUST:  TIME    0.00000E+00  T[Myr]   0.000E+00  Q   0.500E+00  DE   0.000E+00\n" *
            " NS/BH BINARY T[NB] 1.1231445312E+01 NM1,2,S=       386      3732      5306 KW1,2,S=   0  14   0 IPAIR      570 DTAU  1.148E-02 M1,2[NB]  6.718E-05  8.618E-03 R12[NB]  1.913E-04 e,a,eb[NB]=  6.7250E-01  1.1441E-04 -2.530E-03 P[d]=  1.711E+04 H -3.796E+01 GAMMA  0.000E+00 STEP(ICM)  9.766E-04 NPERT    0 NB(ICM)   82 M1,2[*]  2.017E-01  2.588E+01 RAD1,2,S[*]  2.194E-01  1.097E-04  1.385E+04 RI,VI[NB]=  1.0E+00  1.0E+00\n" *
            " MDEND KSREG TIME[NB] 1.1231445312E+01 NM1,2,S=       386      3732      5306 KW1,2,S=   0   5   0 IPAIR      546 DTAU  1.051E-02 M1,2[NB]  6.718E-05  9.580E-03 R12[NB]  1.446E-04 e,a,eb[NB]=  6.9157E-01  1.0813E-04 -2.976E-03 P[d]=  1.492E+04 H -4.461E+01 GAMMA  6.035E-07 STEP(ICM)  9.766E-04 NPERT    0 NB(ICM)   84 M1,2[*]  2.017E-01  2.877E+01 RAD1,2,S[*]  2.194E-01  1.097E-04  1.047E+04 RI,VI[NB]=  1.0E+00  1.0E+00 DM,DMR=  1.0E-03  1.0E-01\n" *
            " MDNEW KSREG TIME[NB] 1.1231445312E+01 NM1,2,S=       386      3732      3732 KW1,2,S=   0  14  14 IPAIR      546 DTAU  1.616E-01 M1,2[NB]  6.718E-05  8.618E-03 R12[NB]  3.210E-03 e,a,eb[NB]=  6.9157E-01  1.0813E-04 -2.976E-03 P[d]=  1.579E+04 H -9.366E-03 GAMMA  6.326E-07 STEP(ICM)  9.537E-07 NPERT   83 NB(ICM)   83 M1,2[*]  2.017E-01  2.588E+01 RAD1,2,S[*]  2.194E-01  1.097E-04  2.324E+05 RI,VI[NB]=  1.0E+00  1.0E+00\n" *
            " ZERO VELOCITY KICK: TIME[NB]  1.123144531E+01 I,N=      1140      3732 K*0 K* K*(ICM) KMECH ICASE    5  14   0   4   1  Time[Myr]  6.37828E+00 M0[M*]  3.50983E+01 MN[M*]  2.58799E+01 ASPN,XJSPIN[NB],SPNFAC=  0.00000E+00  0.00000E+00  2.77000E+13  VI[km/s]  1.08726-157 VK[km/s]  0.00000E+00 VF[km/s]          NaN  VK0[km/s]  0.00000E+00 FBFAC FBTOT ECS   1.00000E+00  2.71555E+01   0  X(1:3)[pc] V(1:3)[km/s]  1.0E+00  1.0E+00  1.0E+00  1.0E+00  1.0E+00  1.0E+00\n" *
            " ADJUST:  TIME    1.00000E+00  T[Myr]   0.000E+00  Q   0.500E+00  DE   0.000E+00\n" *
            " BINARY  COLLISION    NSYS =  2 TIME[NB] 3.8146972656E-06 NM1,2,S=       747       748      6557 KW1,2,S=   0   0   0 IPAIR      233 M1,2[NB]  1.8656E-04  0.0000E+00 e,a,eb[NB],P[d]=  9.9980E-01  4.6994E-07 -8.0169E-03  3.1435E+01 STEP(ICM)  3.9062E-03 NPERT,NB(ICM)=   267   243 M1,2,3[*]  4.5593E-01  1.4905E-01  5.6026E-01 RCOLL,DP,VINF=  9.4938E-11 -1.0164E-20  0.0000E+00 RI,VI[NB]=  2.6305E+00  1.0E+00\n" *
            " ENFORCED MASS LOSS    FAC RIJ VIJ EB KE   0.081    0.01 5796.37 -2.49E+47  3.77E+49\n" *
            " MIX TYPE CHANGE:  NAM K* DM    747   748   1   0  0.04\n" *
            " NEW STAR TIME[NB] 3.8146972656E-06 MIX OF NM1,2=       747       748 IQCOLL -2 K*(I1)  0 K*(I2)  0 K*(INEW)  0 M(I1)[M*]  4.5593E-01 M(I2)[M*]  1.4905E-01 M(INEW)[M*]  5.6026E-01 RS(I1)[R*]  4.0398E-01 RS(I2)[R*]  1.7041E-01 RS(INEW)[R*]  4.9591E-01 DM[M*]  4.4716E-02\n" *
            " HYPERB  COLLISION    NSYS =  2 TIME[NB] 3.7111358643E+01 NM1,2,S=     35841     28959     85841 KW1,2,S=   1   0   0 IPAIR        1 M1,2[NB]  1.8194E-03  0.0000E+00 e,a,eb[NB],P[d]=  1.0003E+00 -1.9929E-03  1.4540E-06 -4.2260E+05 STEP(ICM)  3.0518E-05 NPERT,NB(ICM)=    83    83 M1,2,3[*]  5.2579E+01  9.2376E-02  5.2672E+01 RCOLL,DP,VINF=  5.6541E-07 -2.7105E-19  1.0641E+01 RI,VI[NB]=  2.7579E-01  1.0E+00\n" *
            " BINARY COAL: IQCOLL  3 TIME[NB] 4.0000000000E+00 N1,2        11        12 KW1,2,S   1   1   0 M1,2[NB]  1.00000E-03  2.00000E-03 R12[NB]  1.00000E-04 e,a,eb,dp[NB]=  1.00000E-01  2.00000E-04 -5.00000E-03  1.00000E-06 P[d]=  3.00000E+02  M12S,DM[*]  1.00000E+01  2.00000E+01  2.95000E+01  5.00000E-01 RAD1,2[*]  1.00000E+00  2.00000E+00 RCOLL[R*]  3.00000E+00 VINF[km/s]  0.00000E+00 RI,VI[NB]=  1.00000E+00  1.00000E+00\n" *
            " NEW CHAIN   NCH NP ENERGY RSUM RGRAV TCR RMAXS NAMEC(1-NM) [ALL in NB UNIT]   3  28  -0.01234  1.0E-03  2.0E-03  3.0E-02  4.0E-02       101       102       103\n" *
            " CHTERM: TIME[NB] I3 I4 ICM1 NI3,4,CM  STEP(I3,I4,ICM1)[NB]     5.1234       201       202       203       301       302       303  1.00E-03  2.00E-03  3.00E-03\n" *
            " NEW KSREG   TIME[NB] 2.0000000000E+00 NM1,2,S=         1         2      5001 KW1,2,S=   0   0   0 IPAIR        1 DTAU  1.000E-02 M1,2[NB]  1.000E-04  2.000E-04 R12[NB]  1.000E-03 e,a,eb[NB]=  5.0000E-01  2.0000E-03 -5.000E-06 P[d]=  1.000E+03 H -1.000E-01 GAMMA  1.000E-05 STEP(ICM)  1.000E-03 NPERT    1 NB(ICM)   10 M1,2[*]  1.000E+00  2.000E+00 RAD1,2,S[*]  1.000E+00  1.000E+00  1.000E+03 RI,VI[NB]=  1.0E+00  1.0E+00\n" *
            " LOW VELOCITY KICK: TIME[NB]  4.771484375E+01 I,N=      5902     34667 K*0 K* K*(ICM) KMECH ICASE    8  14   0   4   1  Time[Myr]  4.21285E+00 M0[M*]  1.45478E+01 MN[M*]  1.27114E+01 ASPN,XJSPIN[NB],SPNFAC=  0.00000E+00  0.00000E+00  6.50342E+14  VI[km/s]  2.74219E+00 VK[km/s]  6.03652E+00 VF[km/s]  8.65831E+00  VK0[km/s]  6.03652E+00 FBFAC FBTOT ECS   9.86847E-01  1.27238E+01   0  X(1:3)[pc] V(1:3)[km/s]  1.0E+00  1.0E+00  1.0E+00  1.0E+00  1.0E+00  1.0E+00\n",
        )
        ev = read_engine_events(other)
        @test [e.kind for e in ev] == [
            :ks_new_compact,
            :ks_end_mdot,
            :ks_new_mdot,
            :kick_zero,
            :collision,
            :merged_star,
            :collision,
            :coalescence,
            :chain_start,
            :chain_end,
            :ks_new,
            :kick_low,
        ]
        @test [e.interval for e in ev] == [2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3]
        c = ev[1]
        @test c.time_nb == 11.231445312 && (c.name1, c.name2, c.name3) == (386, 3732, 5306)
        @test (c.kstar1, c.kstar2, c.kstar3) == (0, 14, 0) && c.energy_nb == -2.530e-3
        @test c.mass_solar ≈ 0.2017 + 25.88
        @test ev[2].energy_nb == -2.976e-3 && ev[2].kstar2 == 5 && ev[3].kstar2 == 14
        k = ev[4]       # glued exponent in VI and NaN in VF leave the parsed fields intact
        @test k.name1 == 3732 && k.mass_solar == 25.8799 && k.dmass_solar ≈ 35.0983 - 25.8799
        @test k.velocity_kms == 0.0
        col = ev[5]
        @test col.count == 2 && col.energy_nb == -8.0169e-3 && col.mass_solar == 0.56026
        @test (col.name1, col.name2, col.name3) == (747, 748, 6557) &&
              col.time_nb == 3.8146972656e-6
        ns = ev[6]
        @test (ns.name1, ns.name2) == (747, 748) && (ns.kstar1, ns.kstar2, ns.kstar3) == (0, 0, 0)
        @test ns.mass_solar == 0.56026 && ns.dmass_solar == 4.4716e-2
        @test ev[7].count == 2 && ev[7].kstar1 == 1 && ev[7].mass_solar == 52.672
        co = ev[8]
        @test co.count == 3 && co.time_nb == 4.0 && (co.name1, co.name2) == (11, 12)
        @test (co.kstar1, co.kstar2, co.kstar3) == (1, 1, 0) && co.energy_nb == -5.0e-3
        @test co.mass_solar == 29.5 && co.dmass_solar == 0.5
        ch = ev[9]
        @test isnan(ch.time_nb) && ch.count == 3 && ch.energy_nb == -0.01234
        @test (ch.name1, ch.name2, ch.name3) == (101, 102, 103)
        ce = ev[10]
        @test ce.time_nb == 5.1234 && (ce.name1, ce.name2, ce.name3) == (301, 302, 303)
        nk = ev[11]
        @test nk.time_nb == 2.0 &&
              nk.name3 == 5001 &&
              nk.energy_nb == -5.0e-6 &&
              nk.mass_solar == 3.0
        @test ev[12].velocity_kms == 6.03652 && ev[12].name1 == 34667 && ev[12].kstar2 == 14
        @test_throws ErrorException read_engine_events(joinpath(mktempdir(), "absent"))
        @test isempty(read_engine_events(joinpath(FIXDIR, "lagr.7")))
    end

    @testset "event.35" begin
        row =
            "   0.00000000000000000E+00" *
            repeat("           0", 13) *
            repeat("   0.00000000000000000E+00", 10) *
            "  -0.19225963191497719E+00   0.00000000000000000E+00" *
            repeat("           0", 9) *
            repeat("   0.00000000000000000E+00", 6) *
            "   0.12500000000000000E+01" *
            "       42483        7517" *
            repeat("           0", 14)
        row2 =
            "   0.44146120028899682E+00" *
            repeat("           2", 13) *
            repeat("  -0.10000000000000000E+00", 12) *
            repeat("           3", 9) *
            repeat("   0.50000000000000000E+00", 7) *
            repeat("           1", 16)
        path = joinpath(mktempdir(), "event.35")
        write(
            path,
            "TIME[Myr] NDISS  NTIDE  NSYNC NCOLL NCOAL NDD NCIRC NROCHE NRO NCE NHYP NHYPC NKICK EBIN EMERGE ECOLL EMDOT ECDOT EKICK ESESC EBESC EMESC DEGRAV EBIND MMAX NMDOT NRG NHE NRS NNH NWD NSN NBH NBS ZMRG ZMHE ZMRS ZMNH ZMWD ZMSN ZMDOT NTYPE(1:16)\n" *
            row *
            "\n\n" *
            row2 *
            "\n   1.0  2  3\n",
        )
        cs = read_event_counters(path)
        @test length(cs) == 2          # the short row is skipped
        c0 = cs[1]
        @test c0.time_myr == 0.0 && c0.ndiss == 0 && c0.nkick == 0
        @test c0.ebind == -0.19225963191497719 && c0.mmax == 0.0
        @test c0.zmdot == 1.25 && c0.zmsn == 0.0
        @test c0.ntype[1] == 42483 && c0.ntype[2] == 7517 && c0.ntype[16] == 0
        c1 = cs[2]
        @test c1.time_myr ≈ 0.44146120028899682
        @test c1.ndiss == 2 && c1.nkick == 2 && c1.ebin == -0.1 && c1.mmax == -0.1
        @test c1.nmdot == 3 && c1.nbs == 3 && c1.zmrg == 0.5 && c1.zmdot == 0.5
        @test all(==(1), c1.ntype)
        @test_throws ErrorException read_event_counters(joinpath(mktempdir(), "absent"))
    end

    @testset "lagr.7" begin
        lagr = read_lagr(joinpath(FIXDIR, "lagr.7"))
        @test length(lagr.mass_fractions) == 18
        @test lagr.time[1] == 0.0
        i50 = argmin(abs.(lagr.mass_fractions .- 0.5))
        @test lagr.radii[i50, 1] ≈ 1.2442064 rtol = 1e-6
    end

    @testset "sev.83" begin
        sev = read_stellar_evolution(joinpath(FIXDIR, "sev.83_0"))
        @test sev.n_stars == 40
        @test sev.time_myr == 0.0
        @test length(sev.records) == 40
        r1 = sev.records[1]
        @test r1.time_nb == 0.0
        @test r1.name == Int32(1)
        @test r1.stellar_type == Int32(1)            # massive MS star, Hurley K*=1
        @test r1.ri ≈ 0.694596 rtol = 1e-5           # RI [pc]
        @test r1.mass_solar ≈ 52.5373 rtol = 1e-5
        @test r1.log_teff ≈ 4.71026 rtol = 1e-5
        @test r1.ms_lifetime_myr ≈ 4.46907 rtol = 1e-5   # TM
        @test r1.mass_core == 0.0                    # MC (MS star)
        @test r1.radius_envelope ≈ 1e-10 rtol = 1e-3 # RE placeholder on MS
    end

    @testset "bev.82" begin
        # Truncated excerpt (40 of 224 pairs) of the t = 0 record of a
        # two-cluster run with 20 % primordial binaries per cluster.
        bev = read_binary_evolution(joinpath(FIXDIR, "bev.82_0"))
        @test bev.n_pairs == 40
        @test bev.time_myr == 0.0
        @test length(bev.records) == 40
        r1 = bev.records[1]
        @test (r1.name1, r1.name2) == (1, 2)          # primordial pairs lead the body list
        @test r1.eccentricity ≈ 0.92188 rtol = 1e-5
        @test r1.log_period_days ≈ 1.20647 rtol = 1e-5
        @test r1.mass1 ≈ 0.221134 rtol = 1e-5
        @test r1.ms_lifetime1_myr ≈ 4.98658e5 rtol = 1e-5
        # Every record is a bound orbit whose own columns obey Kepler III.
        for r in bev.records
            P_yr = 10^r.log_period_days / 365.25
            a_au = semi_major_axis_pc(r) / Nbody6Dynamics._AU_IN_PC
            @test a_au^3 / P_yr^2 ≈ r.mass1 + r.mass2 rtol = 2e-2
            @test 0 ≤ r.eccentricity < 1
        end
        @test all(r -> r.name1 < r.name2, bev.records)
    end
end

# =====================================================================

@testset "Types and accessors" begin
    params = zeros(Float32, 20)
    params[1] = 1.5f0   # time
    params[3] = 2.0f0   # rbar
    params[4] = 0.6f0   # zmbar
    params[11] = 8.0f0  # tscale
    params[12] = 4.0f0  # vstar

    hdr = SnapshotHeader(Int32(100), Int32(1), Int32(1), Int32(20), params)

    @test time_nb(hdr) ≈ 1.5
    @test rbar(hdr) ≈ 2.0
    @test tscale(hdr) ≈ 8.0
    @test time_myr(hdr) ≈ 12.0   # 1.5 * 8.0

    u = UnitScaling(2.0, 0.6, 8.0, 4.0)
    @test Nbody6Dynamics.to_pc(u, 1.0) ≈ 2.0
    @test Nbody6Dynamics.to_msun(u, 1.0) ≈ 0.6
    @test Nbody6Dynamics.to_myr(u, 1.0) ≈ 8.0
    @test Nbody6Dynamics.to_kms(u, 1.0) ≈ 4.0
end

# =====================================================================

@testset "Diagnostics — key-value ADJUST format" begin
    diag_path = joinpath(TESTDIR, "out1000_kv")
    write(
        diag_path,
        """
 PHYSICAL SCALING:  R* = 2.500  M* = 1000.0  V* = 4.200  T* = 12.00
                <M> = 0.500  SU = 1.0  AU = 1.0

 ADJUST:  TIME   0.0000  T[Myr]      0.00  Q  1.000  DE  0.00E+00  ETOT  -0.2500
 RMIN =    0.001 RSCALE =    1.234
 TIME[NB]    0.0000 N    10000 <NB>      0 NPAIRS    500

 ADJUST:  TIME   0.5000  T[Myr]     50.00  Q  0.987  DE  1.23E-06  ETOT  -0.2499
 RMIN =    0.002 RSCALE =    1.300
 TIME[NB]    0.5000 N     9990 <NB>      0 NPAIRS    495

 END RUN
""",
    )

    diag = read_diagnostics(diag_path)

    @test length(diag.adjust) == 2
    # First epoch
    @test diag.adjust[1].time_nb ≈ 0.0
    @test diag.adjust[1].time_myr ≈ 0.0
    @test diag.adjust[1].qvir ≈ 1.0
    @test diag.adjust[1].n == 10000
    @test diag.adjust[1].npairs == 500
    @test diag.adjust[1].rscale ≈ 1.234

    # Second epoch — TIME[NB] values override ADJUST defaults
    @test diag.adjust[2].time_nb ≈ 0.5
    @test diag.adjust[2].n == 9990
    @test diag.adjust[2].npairs == 495
    @test diag.adjust[2].rscale ≈ 1.300

    # Physical scaling
    @test diag.physical_scaling["R*"] ≈ 2.5
    @test diag.physical_scaling["T*"] ≈ 12.0
end

# =====================================================================

@testset "Diagnostics — mixed format with partial epochs" begin
    diag_path = joinpath(TESTDIR, "out1000_mixed")
    write(
        diag_path,
        """
 ADJUST:    0.0000      0.00  1.000  0.00E+00 -0.2500  10000    500   1.234
 ADJUST:  TIME   1.0000  T[Myr]    100.00  Q  0.990  DE  5.00E-06  ETOT  -0.2480
 TIME[NB]    1.0000 N     9500 <NB>      0 NPAIRS    450
""",
    )
    diag = read_diagnostics(diag_path)
    @test length(diag.adjust) == 2

    # First: positional format
    @test diag.adjust[1].time_nb ≈ 0.0
    @test diag.adjust[1].n == 10000
    @test diag.adjust[1].rscale ≈ 1.234

    # Second: key-value + TIME[NB] merge
    @test diag.adjust[2].time_nb ≈ 1.0
    @test diag.adjust[2].n == 9500
    @test diag.adjust[2].npairs == 450
end

# =====================================================================

@testset "Degenerate reader inputs" begin
    # sev.83 with zero stars: header-only file parses to an empty snapshot
    p_empty = joinpath(TESTDIR, "sev.83_empty")
    write(p_empty, "  0  0.0\n")
    sev0 = read_stellar_evolution(p_empty)
    @test sev0.n_stars == 0 && isempty(sev0.records)

    # esc.11 line with only 12 tokens (angle column truncated) is skipped
    p_trunc = joinpath(TESTDIR, "esc_trunc.11")
    write(p_trunc, "  1.0 5.0e-4 20.0 40.0 2.0e-3 1.234 0.5 -0.1 15.6 0 101 22.2\n")
    @test isempty(read_escapers(p_trunc))
end

@testset "Threaded snapshot reading" begin
    dir = mktempdir()
    write_snap(path, t, n) = open(path, "w") do io
        _write_fortran_record(io, Int32[n, 1, 1, 20])
        params = zeros(Float32, 20)
        params[1] = Float32(t)
        params[3] = 1.0f0
        params[4] = 0.5f0
        params[11] = 10.0f0
        params[12] = 5.0f0
        params[18] = 2.0f0
        _write_fortran_record(io, params)
        for i in 1:n
            buf = IOBuffer()
            write(
                buf,
                Float32(1 / n),
                Float32(i),
                Float32(-i),
                0.0f0,
                0.1f0,
                0.0f0,
                0.0f0,
                Int32(i),
            )
            _write_fortran_record(io, take!(buf))
        end
    end
    for (k, t) in enumerate((0.0, 0.5, 1.0, 1.5, 2.0))
        write_snap(joinpath(dir, "conf.3_$(t)"), t, 4 + k)
    end
    write(joinpath(dir, "conf.3_2.5"), "corrupt")
    serial =
        @test_logs (:warn, r"Skipping corrupt snapshot") (:warn, r"1 / 6 snapshot") read_all_conf3(
            dir;
            threaded = false,
        )
    threaded =
        @test_logs (:warn, r"Skipping corrupt snapshot") (:warn, r"1 / 6 snapshot") read_all_conf3(
            dir;
            threaded = true,
        )
    @test length(serial) == length(threaded) == 5
    @test [time_nb(s.header) for s in threaded] == [0.0, 0.5, 1.0, 1.5, 2.0]
    @test all(
        a.pos == b.pos && a.mass == b.mass && a.name == b.name for (a, b) in zip(serial, threaded)
    )
    @test [nparticles(s) for s in threaded] == [5, 6, 7, 8, 9]
    # The ordered reader is generic over the element type and the reader
    names = Nbody6Dynamics._read_ordered(
        String,
        p -> uppercase(basename(p)),
        dir,
        ["conf.3_0.0", "conf.3_1.0"];
        desc = "",
        what = "name",
        threaded = true,
    )
    @test names == ["CONF.3_0.0", "CONF.3_1.0"]
    @test occursin(
        "--threads=$(Threads.nthreads())",
        string(Nbody6Dynamics._sweep_worker_command(dir)),
    )
end

# =====================================================================
