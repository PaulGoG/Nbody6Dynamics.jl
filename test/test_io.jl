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
            x, y, z = Float32.(randn(3))
            vx, vy, vz = Float32.(0.1 .* randn(3))
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
            write(buf, Float32.(randn(3))...)  # pos
            write(buf, Float32.(0.1 .* randn(3))...)  # vel
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
