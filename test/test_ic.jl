# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Initial Conditions — Merger IC Generator" begin
    # StableRNGs guarantees an identical stream across Julia versions, so
    # reference values in these tests survive upgrades.
    using StableRNGs
    rng = StableRNG(12345)

    # --- Plummer sampler ---
    @testset "Plummer sampler" begin
        N = 1000
        a = 1.0
        pos, vel = sample_plummer(N, a; rng = rng)
        @test size(pos) == (3, N)
        @test size(vel) == (3, N)
        # Particles should be centred near origin
        cm = vec(sum(pos, dims = 2)) ./ N
        @test all(abs.(cm) .< 0.5)
        # Radii should be finite and positive
        radii = [sqrt(sum(pos[:, i] .^ 2)) for i in 1:N]
        @test all(radii .> 0)
        @test all(isfinite.(radii))
        # Half-mass radius should be roughly 1.305 × a for Plummer
        r_sorted = sort(radii)
        r_hm = r_sorted[N ÷ 2]
        @test 0.5 * a < r_hm < 3.0 * a
    end

    # --- King sampler ---
    @testset "King sampler" begin
        N = 1000
        W0 = 6.0
        rt = 10.0
        pos, vel = sample_king(N, W0, rt; rng = rng)
        @test size(pos) == (3, N)
        @test size(vel) == (3, N)
        # All particles should be within the tidal radius (with some tolerance
        # from the sampling/scaling)
        radii = [sqrt(sum(pos[:, i] .^ 2)) for i in 1:N]
        @test all(isfinite.(radii))
        # Invalid W0
        @test_throws ArgumentError sample_king(10, -1.0, 5.0)
        @test_throws ArgumentError sample_king(10, 6.0, -1.0)
    end

    # --- King ODE solver ---
    @testset "King ODE solver" begin
        rhat, What, rho = Nbody6Dynamics._solve_king(6.0)
        @test length(rhat) == length(What) == length(rho)
        @test What[1] ≈ 6.0
        @test What[end] ≈ 0.0 atol = 0.01
        @test rhat[1] ≈ 0.0
        @test rhat[end] > 0.0
        # Density should be monotonically decreasing
        @test all(diff(rho[1:(end - 1)]) .≤ 0.01)
    end

    # --- King concentration vs published values (physics validation) ---
    @testset "King concentration c(W0)" begin
        # c = log10(r_t/r_0) for the standard dimensionless King equation;
        # reference values from the King (1966) model tables.
        for (W0, c_ref) in
            [(3.0, 0.672), (5.0, 1.029), (6.0, 1.255), (7.0, 1.528), (9.0, 2.119), (12.0, 2.739)]
            rhat, _, _ = Nbody6Dynamics._solve_king(W0)
            c = log10(rhat[end])
            @test isapprox(c, c_ref; rtol = 0.01)
        end
    end

    # --- Kroupa IMF ---
    @testset "Kroupa IMF" begin
        masses = sample_kroupa(5000; m_low = 0.08, m_up = 100.0, rng = rng)
        @test length(masses) == 5000
        @test all(masses .≥ 0.08)
        @test all(masses .≤ 100.0)
        # Most stars should be low-mass (IMF is bottom-heavy)
        @test count(m -> m < 1.0, masses) > 3000
        # Custom mass limits
        masses_narrow = sample_kroupa(500; m_low = 1.0, m_up = 10.0, rng = rng)
        @test all(1.0 .≤ masses_narrow .≤ 10.0)
    end

    # --- Virialisation ---
    @testset "Virialise" begin
        N = 200
        mass = ones(N) ./ N
        pos = randn(rng, 3, N)
        vel = randn(rng, 3, N) .* 0.1
        virialise!(mass, pos, vel)
        # CM should be at origin
        cm_pos = vec(sum(mass' .* pos, dims = 2))
        cm_vel = vec(sum(mass' .* vel, dims = 2))
        @test all(abs.(cm_pos) .< 1e-10)
        @test all(abs.(cm_vel) .< 1e-10)
        # Virial ratio should be ~0.5
        T = 0.5 * sum(mass[i] * sum(vel[:, i] .^ 2) for i in 1:N)
        W = 0.0
        for i in 1:N, j in (i + 1):N
            dr = pos[:, i] .- pos[:, j]
            W -= mass[i] * mass[j] / sqrt(sum(dr .^ 2))
        end
        Q = T / abs(W)
        @test Q ≈ 0.5 atol = 0.01
    end

    # --- Kepler velocity ---
    @testset "Kepler velocity" begin
        # Equal mass, circular orbit
        v1, v2 = kepler_velocity(1.0, 1.0, 10.0, 0.0)
        @test v1 ≈ v2  # symmetric
        @test v1 > 0
        # Eccentric orbit has lower apocentre velocity
        v1e, v2e = kepler_velocity(1.0, 1.0, 10.0, 0.9)
        @test v1e < v1
        # Invalid inputs
        @test_throws ArgumentError kepler_velocity(1.0, 1.0, 10.0, 1.0)
        @test_throws ArgumentError kepler_velocity(1.0, 1.0, -1.0, 0.5)
    end

    # --- Jacobi radius ---
    @testset "Jacobi radius" begin
        rj = jacobi_radius(10.0, 1.0, 3.0)
        @test rj ≈ 10.0 * (1.0 / 9.0)^(1/3)
        # Equal mass: rJ = d × (1/3)^{1/3}
        rj_eq = jacobi_radius(10.0, 1.0, 1.0)
        @test rj_eq ≈ 10.0 * (1/3)^(1/3) atol = 1e-10
    end

    # --- Two-cluster orbit setup ---
    @testset "Two-cluster orbit" begin
        N1, N2 = 100, 100
        pos1 = randn(rng, 3, N1) .* 0.5
        vel1 = randn(rng, 3, N1) .* 0.01
        mass1 = ones(N1) ./ N1
        pos2 = randn(rng, 3, N2) .* 0.5
        vel2 = randn(rng, 3, N2) .* 0.01
        mass2 = ones(N2) ./ N2
        # Centre each cluster exactly (as virialise! does in the real
        # pipeline) so the COM/momentum/energy identities below are exact.
        for (p, v, m) in ((pos1, vel1, mass1), (pos2, vel2, mass2))
            p .-= sum(m' .* p, dims = 2) ./ sum(m)
            v .-= sum(m' .* v, dims = 2) ./ sum(m)
        end

        pos, vel, mass, ranges, _ = Nbody6Dynamics.setup_two_cluster_orbit(
            pos1,
            vel1,
            mass1,
            pos2,
            vel2,
            mass2,
            20.0,
            0.5;
            truncate_jacobi_flag = false,
        )
        @test length(mass) == N1 + N2
        @test size(pos, 2) == N1 + N2
        @test ranges == [1:N1, (N1 + 1):(N1 + N2)]
        # Combined CM should be near origin
        M = sum(mass)
        cm = vec(sum(mass' .* pos, dims = 2)) ./ M
        @test all(abs.(cm) .< 0.1)

        # Total momentum must vanish and the two-COM orbit must have the
        # Keplerian energy of the requested (d_apo, e) orbit.
        p = vec(sum(mass' .* vel, dims = 2))
        @test all(abs.(p) .< 1e-10)
        r1 = 1:N1
        r2 = (N1 + 1):(N1 + N2)
        M1 = sum(mass[r1])
        M2 = sum(mass[r2])
        com1 = vec(sum(mass[r1]' .* pos[:, r1], dims = 2)) ./ M1
        com2 = vec(sum(mass[r2]' .* pos[:, r2], dims = 2)) ./ M2
        vcom1 = vec(sum(mass[r1]' .* vel[:, r1], dims = 2)) ./ M1
        vcom2 = vec(sum(mass[r2]' .* vel[:, r2], dims = 2)) ./ M2
        d = sqrt(sum((com1 .- com2) .^ 2))
        @test d ≈ 20.0 rtol = 1e-10
        v_rel2 = sum((vcom1 .- vcom2) .^ 2)
        a_orb = 20.0 / (1.0 + 0.5)
        ε = 0.5 * v_rel2 - (M1 + M2) / d          # specific orbital energy, G = 1
        @test ε ≈ -(M1 + M2) / (2a_orb) rtol = 1e-10
    end

    # --- Plummer half-mass relation (physics validation) ---
    @testset "Plummer r_hm = 1.305 a" begin
        N = 20_000
        a = 1.0
        pos, _ = sample_plummer(N, a; rng = rng)
        mass = fill(1.0 / N, N)
        r_hm = Nbody6Dynamics.half_mass_radius(mass, pos; centre = zeros(3))
        @test r_hm ≈ 1.3048 rtol = 0.05    # statistical tolerance
    end

    # --- Kroupa mean mass vs analytic (physics validation) ---
    @testset "Kroupa mean mass" begin
        @test kroupa_mean_mass(0.08, 100.0) ≈ 0.58 rtol = 0.05
        N = 50_000
        m = sample_kroupa(N; m_low = 0.08, m_up = 100.0, rng = rng)
        μ = kroupa_mean_mass(0.08, 100.0)
        # Sample mean within a generous statistical band of the analytic mean
        @test abs(sum(m) / N - μ) / μ < 0.1
    end

    # --- verif_triorbit.toml Lagrange equilibrium (config regression) ---
    @testset "verif_triorbit config equilibrium" begin
        path = joinpath(@__DIR__, "..", "input_files", "verif_triorbit.toml")
        cfg = load_merger_config(path)
        @test length(cfg.clusters) == 3
        m = expected_mass(cfg.clusters[1].imf, cfg.clusters[1].N)   # 900 M☉, equal bodies
        @test m ≈ 900.0
        rc = 6.0
        ω = sqrt(Nbody6Dynamics._G_PC_KMS2_MSUN * m / (sqrt(3) * rc^3))   # km s⁻¹ pc⁻¹
        for spec in cfg.clusters
            r⃗ = spec.position
            v⃗ = spec.velocity
            @test sqrt(sum(r⃗ .^ 2)) ≈ rc rtol = 1e-3
            # Speed = ω r and velocity ⟂ radius (circular Lagrange orbit)
            @test sqrt(sum(v⃗ .^ 2)) ≈ ω * rc rtol = 1e-3
            @test abs(sum(r⃗ .* v⃗)) / (rc * ω * rc) < 1e-3
        end
        # Zero net momentum for equal-mass clusters
        vsum = sum(spec.velocity for spec in cfg.clusters)
        @test all(abs.(vsum) .< 0.05 * ω * rc)
    end

    # --- dat.10 writer ---
    @testset "dat.10 writer" begin
        N = 50
        mass = ones(N) ./ N
        pos = randn(rng, 3, N)
        vel = randn(rng, 3, N)
        dat10_path = joinpath(TESTDIR, "test_dat10.dat")
        write_dat10(dat10_path, mass, pos, vel)
        lines = readlines(dat10_path)
        @test length(lines) == N
        # Each line should have 7 columns
        cols = split(lines[1])
        @test length(cols) == 7
        # First column should be mass
        @test parse(Float64, cols[1]) ≈ 1.0 / N
    end

    # --- N-body unit conversion ---
    @testset "N-body unit conversion" begin
        mass = [0.5, 0.5]
        pos = [1.0 -1.0; 0.0 0.0; 0.0 0.0]
        vel = [0.0 0.0; 1.0 -1.0; 0.0 0.0]
        M_tot = 1e5  # M☉
        rbar = 2.0   # pc
        to_nbody_units!(mass, pos, vel, M_tot, rbar)
        @test sum(mass) ≈ 1e-5  # each mass = 0.5/1e5
        @test pos[1, 1] ≈ 0.5   # 1.0 / 2.0
    end

    # --- TOML config loading (kepler mode) ---
    @testset "Merger config TOML — kepler" begin
        cfg_path = joinpath(TESTDIR, "test_merger.toml")
        write(
            cfg_path,
            """
[merger]
n_clusters = 2
orbit_mode = "kepler"

[merger.cluster1]
model = "plummer"
N = 100
mass_total = 60.0
rbar = 1.0

[merger.cluster2]
model = "king"
N = 100
W0 = 5.0
mass_total = 60.0
rbar = 1.0

[merger.orbit]
apocentre = 10.0
eccentricity = 0.5

[merger.output]
format = "nbody"
truncate_jacobi = false
""",
        )
        cfg = load_merger_config(cfg_path)
        @test length(cfg.clusters) == 2
        @test cfg.orbit_mode == "kepler"
        @test cfg.clusters[1].profile isa PlummerProfile
        @test cfg.clusters[2].profile isa KingProfile
        @test cfg.clusters[2].profile.W0 == 5.0
        @test cfg.orbit.apocentre == 10.0
        @test cfg.orbit.eccentricity == 0.5
        @test cfg.output.format == "nbody"
        @test cfg.output.truncate_jacobi == false

        # Unknown merger keys are refused, naming the key.
        bad = joinpath(TESTDIR, "merger_unknown.toml")
        base_text = read(joinpath(@__DIR__, "..", "input_files", "merger_demo_small.toml"), String)
        write(bad, replace(base_text, "eccentricity" => "ecentricity"; count = 1))
        err = try
            load_merger_config(bad)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("ecentricity", err.msg)
        write(bad, base_text * "\n[merger.cluster3]\nN = 10\n")
        err = try
            load_merger_config(bad)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("merger.cluster3", err.msg)

        # The O(N²) virialisation limit is a config key; clusters above it are refused.
        write(bad, replace(base_text, "[merger]\n" => "[merger]\nvirial_max_n = 500\n"; count = 1))
        err = try
            load_merger_config(bad)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError &&
              occursin("merger.cluster1.N", err.msg) &&
              occursin("merger.virial_max_n", err.msg)
        write(bad, replace(base_text, "[merger]\n" => "[merger]\nvirial_max_n = 0\n"; count = 1))
        err = try
            load_merger_config(bad)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("merger.virial_max_n", err.msg)
        write(bad, replace(base_text, "[merger]\n" => "[merger]\nvirial_max_n = 5000\n"; count = 1))
        @test load_merger_config(bad).virial_max_n == 5000
        write(bad, base_text)
        @test load_merger_config(bad).virial_max_n == Nbody6Dynamics._VIRIAL_NMAX
    end

    # --- TOML config loading (explicit mode) ---
    @testset "Merger config TOML — explicit" begin
        cfg_path = joinpath(TESTDIR, "test_merger_explicit.toml")
        write(
            cfg_path,
            """
[merger]
n_clusters = 3
orbit_mode = "explicit"

[merger.cluster1]
model = "plummer"
N = 50
mass_total = 30.0
rbar = 1.0
position = [-5.0, 0.0, 0.0]
velocity = [1.0, 0.0, 0.0]

[merger.cluster2]
model = "king"
N = 50
W0 = 4.0
mass_total = 30.0
rbar = 1.0
position = [2.5, 4.33, 0.0]
velocity = [-0.5, -0.87, 0.0]

[merger.cluster3]
model = "plummer"
N = 50
mass_total = 30.0
rbar = 0.8
position = [2.5, -4.33, 0.0]
velocity = [-0.5, 0.87, 0.0]

[merger.output]
format = "nbody"
truncate_jacobi = false
""",
        )
        cfg = load_merger_config(cfg_path)
        @test length(cfg.clusters) == 3
        @test cfg.orbit_mode == "explicit"
        @test cfg.clusters[1].position == [-5.0, 0.0, 0.0]
        @test cfg.clusters[3].velocity == [-0.5, 0.87, 0.0]
    end

    # --- Fail-fast validation of merger TOML values ---
    @testset "Merger config validation" begin
        # Template with one overridable block per constrained section
        function _merger_toml(;
            cluster1::String = "model = \"king\"\nN = 100\nW0 = 5.0\nrbar = 1.0",
            orbit::String = "apocentre = 10.0\neccentricity = 0.5",
            output::String = "tcrit = 10.0\ndtadj = 1.0\ndeltat = 1.0",
        )
            """
            [merger]
            n_clusters = 2
            orbit_mode = "kepler"

            [merger.cluster1]
            $cluster1

            [merger.cluster2]
            model = "plummer"
            N = 100
            rbar = 1.0

            [merger.orbit]
            $orbit

            [merger.output]
            $output
            """
        end

        # Valid template loads cleanly — smoke test
        ok_path = joinpath(TESTDIR, "merger_val_ok.toml")
        write(ok_path, _merger_toml())
        @test load_merger_config(ok_path) isa MergerConfig

        cases = [
            (
                "N",
                _merger_toml(cluster1 = "model = \"king\"\nN = 1\nrbar = 1.0"),
                "merger.cluster1.N",
            ),
            (
                "rbar",
                _merger_toml(cluster1 = "model = \"king\"\nN = 100\nrbar = -1.0"),
                "merger.cluster1.rbar",
            ),
            (
                "W0",
                _merger_toml(cluster1 = "model = \"king\"\nN = 100\nW0 = -2.0\nrbar = 1.0"),
                "merger.cluster1.W0",
            ),
            (
                "kroupa_bounds",
                _merger_toml(
                    cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n" *
                               "imf = \"kroupa\"\nbodyn = 50.0\nbody1 = 0.1",
                ),
                "0 < bodyn < body1",
            ),
            (
                "target_mass",
                _merger_toml(
                    cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n" *
                               "imf = \"kroupa\"\nmass_total = -500.0",
                ),
                "target_mass",
            ),
            (
                "particle_mass",
                _merger_toml(
                    cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n" *
                               "imf = \"equal\"\nmass_total = 0.0",
                ),
                "particle_mass",
            ),
            (
                "apocentre",
                _merger_toml(orbit = "apocentre = -5.0\neccentricity = 0.5"),
                "merger.orbit.apocentre",
            ),
            (
                "eccentricity",
                _merger_toml(orbit = "apocentre = 10.0\neccentricity = 1.0"),
                "merger.orbit.eccentricity",
            ),
            (
                "tcrit",
                _merger_toml(output = "tcrit = 0.0\ndtadj = 1.0\ndeltat = 1.0"),
                "merger.output.tcrit",
            ),
            (
                "dtadj",
                _merger_toml(output = "tcrit = 10.0\ndtadj = -0.5\ndeltat = 1.0"),
                "merger.output.dtadj",
            ),
            (
                "deltat",
                _merger_toml(output = "tcrit = 10.0\ndtadj = 1.0\ndeltat = 0.0"),
                "merger.output.deltat",
            ),
            ("nbody6_qe", _merger_toml() * "\n[merger.nbody6]\nqe = 0.0\n", "merger.nbody6.qe"),
            ("nbody6_kz16", _merger_toml() * "\n[merger.nbody6]\nkz16 = 5\n", "merger.nbody6.kz16"),
            ("nbody6_rs0", _merger_toml() * "\n[merger.nbody6]\nrs0 = -0.1\n", "merger.nbody6.rs0"),
            (
                "nbody6_tcomp",
                _merger_toml() * "\n[merger.nbody6]\ntcomp = 0.0\n",
                "merger.nbody6.tcomp",
            ),
            (
                "nbody6_gmin_gmax",
                _merger_toml() * "\n[merger.nbody6]\ngmin = 0.1\ngmax = 0.01\n",
                "gmin < gmax",
            ),
            (
                "nbody6_ncrit",
                _merger_toml() * "\n[merger.nbody6]\nncrit = 0\n",
                "merger.nbody6.ncrit",
            ),
            (
                "nbody6_kz_index",
                _merger_toml() * "\n[merger.nbody6.kz]\n\"99\" = 1\n",
                "merger.nbody6.kz index",
            ),
            (
                "stellar_zmet",
                _merger_toml() * "\n[merger.stellar]\nzmet = 0.1\n",
                "merger.stellar.zmet",
            ),
            (
                "stellar_level",
                _merger_toml() * "\n[merger.stellar]\nlevel = \"X\"\n",
                "merger.stellar.level",
            ),
            (
                "stellar_dtplot",
                _merger_toml() * "\n[merger.stellar]\ndtplot = 0.5\n",
                "merger.stellar.dtplot",
            ),
            (
                "stellar_epoch0",
                _merger_toml() * "\n[merger.stellar]\nepoch0 = 5.0\n",
                "merger.stellar.epoch0",
            ),
            ("tidal_kz14_3", _merger_toml() * "\n[merger.tidal]\nkz14 = 3\n", "not supported"),
            (
                "tidal_kz14_bad",
                _merger_toml() * "\n[merger.tidal]\nkz14 = 7\n",
                "merger.tidal.kz14",
            ),
            (
                "tidal_qe",
                _merger_toml() * "\n[merger.tidal]\nkz14 = 2\ngmg = 1.0e11\nrg0 = 8.5\n",
                "tidal work",
            ),
            (
                "tidal_gmg",
                _merger_toml() * "\n[merger.tidal]\nkz14 = 2\nrg0 = 8.5\n",
                "merger.tidal.gmg",
            ),
            (
                "tidal_rg",
                _merger_toml() * "\n[merger.tidal]\nkz14 = 5\nvg = [0.0, 220.0, 0.0]\n",
                "merger.tidal.rg",
            ),
            (
                "binaries_fraction",
                _merger_toml(
                    cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n[merger.cluster1.binaries]\nfraction = 1.5",
                ),
                "binaries.fraction",
            ),
            (
                "imf_rescale_factor",
                _merger_toml(
                    cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n" *
                               "imf = \"kroupa\"\nmass_total = 5000.0",
                ),
                "rescale factor",
            ),
        ]
        for (label, body, expected) in cases
            path = joinpath(TESTDIR, "merger_val_bad_$label.toml")
            write(path, body)
            @test_throws expected load_merger_config(path)
        end

        # Moderate rescale factors load with a warning (×1.6 here)
        warn_path = joinpath(TESTDIR, "merger_val_warn_rescale.toml")
        write(
            warn_path,
            _merger_toml(
                cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\nimf = \"kroupa\"\n" *
                           "mass_total = $(round(1.6 * 100 * kroupa_mean_mass(0.08, 100.0); digits = 1))",
            ),
        )
        @test_logs (:warn, r"rescales the Kroupa IMF") match_mode = :any load_merger_config(
            warn_path,
        )
    end

    # --- [merger.nbody6] integration parameters and regime guards ---
    @testset "Nbody6 integration parameters" begin
        p0 = Nbody6ParameterSpec()
        @test p0.qe == 2.0e-4 && p0.nnbopt == 0 && p0.rs0 == 0.0 && p0.kz16 == 0
        ρ̂_plummer = Nbody6Dynamics._central_density_contrast(PlummerProfile())
        @test ρ̂_plummer ≈ 2 * 1.305^3 rtol = 0.01
        c6 = Nbody6Dynamics._central_density_contrast(KingProfile(W0 = 6.0))
        c3 = Nbody6Dynamics._central_density_contrast(KingProfile(W0 = 3.0))
        @test isfinite(c6) && c6 > c3 > 1        # monotonic in W0 (3.3 and 9.3)
        @test c3 < ρ̂_plummer < c6                # Plummer's extended halo sits between W0 = 3 and 6

        clusters = [
            ClusterSpec(profile = PlummerProfile(), N = 100, rbar = 1.0, imf = KroupaIMF()),
            ClusterSpec(profile = KingProfile(W0 = 6.0), N = 120, rbar = 2.0, imf = KroupaIMF()),
        ]
        ranges = [1:100, 101:220]
        r = resolve_nbody6_parameters(p0, clusters, ranges, 220, 4.0)
        @test r.nnbopt == 20                          # clamp(round(√220) = 15, 20, 300)
        @test r.rs0 ≈ min(2 * 0.25 * cbrt(2 * 20 / 100), 0.25)   # r_h = 1/4, N_min = 100, doubled rule
        @test r.rmin ≈ 4 * 0.25 / (100 * cbrt(ρ̂_plummer))
        @test r.dtmin ≈ 0.04 * sqrt(r.rmin^3 * 220)
        @test r.qe == 2.0e-4 && r.kz16 == 0
        user = Nbody6ParameterSpec(;
            qe = 1e-3,
            nnbopt = 50,
            rs0 = 0.1,
            rmin = 1e-4,
            dtmin = 1e-6,
            kz16 = 2,
        )
        ru = resolve_nbody6_parameters(user, clusters, ranges, 220, 4.0)
        @test ru.nnbopt == 50 && ru.rs0 == 0.1 && ru.rmin == 1e-4 && ru.dtmin == 1e-6
        @test ru.kz16 == 2 && ru.qe == 1e-3
        @test_throws ArgumentError Nbody6Dynamics._assert_resolved(p0)

        # Guards: oversized RS0 refused, cold collapse and unresolved members warned
        wide = Nbody6ParameterSpec(; nnbopt = 20, rs0 = 0.5, rmin = 1e-4, dtmin = 1e-6)
        @test_throws ErrorException Nbody6Dynamics._check_multicluster_regime(0.5, 2.0, wide, 0.25)
        @test_logs (:warn, r"cold-collapse") Nbody6Dynamics._check_multicluster_regime(
            0.1,
            2.0,
            r,
            0.25,
        )
        @test_logs (:warn, r"unresolved") Nbody6Dynamics._check_multicluster_regime(
            0.5,
            8.0,
            r,
            0.25,
        )
        @test_logs Nbody6Dynamics._check_multicluster_regime(0.5, 2.0, r, 0.25)
        @test_logs (:warn, r"undersampled") Nbody6Dynamics._check_multicluster_regime(
            0.5,
            2.0,
            r,
            0.25;
            deltat = 1.0,
            t_cr_member_min_nb = 0.2,
        )
        @test_logs Nbody6Dynamics._check_multicluster_regime(
            0.5,
            2.0,
            r,
            0.25;
            deltat = 0.1,
            t_cr_member_min_nb = 0.2,
        )

        # Writer: unresolved specs refused, resolved values written
        inp = joinpath(mktempdir(), "t.inp")
        @test_throws ArgumentError generate_merger_inp(inp, 220, 4.0, 0.6; nbody6 = p0)
        generate_merger_inp(inp, 220, 4.0, 0.6; nbody6 = ru, tcrit = 5.0)
        txt = read(inp, String)
        @test occursin("NNBOPT=50,", txt) && occursin("QE=1.000E-03", txt)
        @test occursin("RS0=0.1,", txt) && occursin("DTMIN=1.000E-06,RMIN=1.000E-04,", txt)
        @test occursin("KZ(11:20)=0 1 0 0 0 2 0 0 3 0", txt)

        # [merger.nbody6] round trip
        nb_path = joinpath(TESTDIR, "merger_nbody6.toml")
        write(
            nb_path,
            """
[merger]
n_clusters = 2
orbit_mode = "kepler"

[merger.cluster1]
model = "king"
N = 100
rbar = 1.0

[merger.cluster2]
model = "plummer"
N = 100
rbar = 1.0

[merger.nbody6]
qe = 1.0e-3
nnbopt = 30
kz16 = 1
tcomp = 7200.0
tcrtp0 = 500.0
isernb = 20
ncrit = 5
smax = 0.5

[merger.nbody6.kz]
"40" = 2
"19" = 4

[merger.stellar]
kz19 = 0
level = "0"
zmet = 0.02
epoch0 = -1.0
dtplot = 2.0
""",
        )
        cfg_nb = load_merger_config(nb_path)
        @test cfg_nb.nbody6.qe == 1.0e-3 && cfg_nb.nbody6.nnbopt == 30
        @test cfg_nb.nbody6.kz16 == 1 && cfg_nb.nbody6.rs0 == 0.0
        @test cfg_nb.nbody6.tcomp == 7200.0 && cfg_nb.nbody6.tcrtp0 == 500.0
        @test cfg_nb.nbody6.isernb == 20 && cfg_nb.nbody6.ncrit == 5 && cfg_nb.nbody6.smax == 0.5
        @test cfg_nb.nbody6.kz == Dict(40 => 2, 19 => 4)
        @test cfg_nb.stellar.kz19 == 0 && cfg_nb.stellar.level == "0"
        @test cfg_nb.stellar.zmet == 0.02 &&
              cfg_nb.stellar.epoch0 == -1.0 &&
              cfg_nb.stellar.dtplot == 2.0

        # Writer: run control, stellar settings, KZ overrides (an override of a named
        # index warns), and mass bounds reach the file
        rnb = resolve_nbody6_parameters(cfg_nb.nbody6, clusters, ranges, 220, 4.0)
        inp2 = joinpath(mktempdir(), "t2.inp")
        @test_logs (:warn, r"overrides KZ\(19\)") match_mode = :any generate_merger_inp(
            inp2,
            220,
            4.0,
            0.6;
            nbody6 = rnb,
            stellar = cfg_nb.stellar,
            mass_bounds = (0.5, 20.0),
        )
        txt2 = read(inp2, String)
        @test occursin("KSTART=1,TCOMP=7200,TCRTP0=500,isernb=20,iserreg=40,iserks=0 /", txt2)
        @test occursin("N=220,NFIX=1,NCRIT=5,", txt2) &&
              occursin(",NNBOPT=30,NRUN=1,NCOMM=10,", txt2)
        @test occursin("SMAX=0.5,", txt2) && occursin("Level='0' /", txt2)
        @test occursin("KZ(11:20)=0 0 0 0 0 1 0 0 4 0", txt2)   # KZ(12) off with kz19 = 0; override 19 → 4
        @test occursin("KZ(31:40)=0 0 0 0 0 0 0 0 0 2", txt2)
        @test occursin("BODY1=20,BODYN=0.5,NBIN0=0,NHI0=0,ZMET=0.02,EPOCH0=-1,DTPLOT=2 /", txt2)

        # Tidal field: point-mass and MWPotential2014 namelists, isolated by default
        @test cfg_nb.tidal.kz14 == 0
        td2 = TidalSpec(; kz14 = 2, gmg = 1.0e11, rg0 = 8.5)
        inp3 = joinpath(mktempdir(), "t3.inp")
        generate_merger_inp(inp3, 220, 4.0, 0.6; nbody6 = ru, tidal = td2)
        txt3 = read(inp3, String)
        @test occursin("KZ(11:20)=0 1 0 2 0 2 0 0 3 0", txt3)
        @test occursin("&INXTRNL0\nGMG=1E+11,RG0=8.5 /", txt3)
        td5 = TidalSpec(; kz14 = 5, rg = [8.0, 0.0, 0.0], vg = [0.0, 220.0, 0.0])
        generate_merger_inp(inp3, 220, 4.0, 0.6; nbody6 = ru, tidal = td5)
        txt5 = read(inp3, String)
        @test occursin("KZ(11:20)=0 1 0 5 0 2 0 0 3 0", txt5)
        @test occursin("&INXTRNL0\nRG=8,0,0,VG=0,220,0 /", txt5)
        generate_merger_inp(inp3, 220, 4.0, 0.6; nbody6 = ru, tidal = TidalSpec(; kz14 = 1))
        @test !occursin("&INXTRNL0", read(inp3, String))
        @test_throws ArgumentError generate_merger_inp(
            inp3,
            220,
            4.0,
            0.6;
            nbody6 = ru,
            tidal = TidalSpec(; kz14 = 4),
        )
        @test_throws ArgumentError Nbody6Dynamics._validate_tidal(TidalSpec(; kz14 = 2, gmg = 1e11))
        @test_throws ArgumentError Nbody6Dynamics._validate_tidal(
            TidalSpec(; kz14 = 5, rg = [8.0, 0, 0]),
        )
        @test_throws ArgumentError Nbody6Dynamics._validate_tidal_tolerance(
            td2,
            Nbody6ParameterSpec(),
        )
        @test Nbody6Dynamics._validate_tidal_tolerance(td2, Nbody6ParameterSpec(; qe = 0.05)) ===
              nothing
        @test Nbody6Dynamics._validate_tidal_tolerance(TidalSpec(), Nbody6ParameterSpec()) ===
              nothing

        # Crossing time: G = 1 virial system with E = -M²/(4 r_v) has t_cr = (2 r_v)^{3/2}/√M
        @test crossing_time(1.0, -0.25) ≈ 2.0^1.5
        @test isnan(crossing_time(1.0, 0.1))
        @test Nbody6Dynamics._nbody_time_myr(1.0, 1.0) ≈ 14.91 rtol = 0.01
    end

    # --- Full pipeline — kepler mode (small N) ---
    @testset "Full merger IC pipeline — kepler" begin
        out_dir = mktempdir()
        cfg = MergerConfig(
            [
                ClusterSpec(
                    profile = PlummerProfile(),
                    N = 100,
                    rbar = 1.0,
                    imf = RescaledKroupaIMF(
                        bodyn = 0.1,
                        body1 = 50.0,
                        target_mass = 100 * kroupa_mean_mass(0.1, 50.0),
                    ),
                ),
                ClusterSpec(
                    profile = KingProfile(W0 = 5.0),
                    N = 100,
                    rbar = 1.0,
                    imf = RescaledKroupaIMF(
                        bodyn = 0.1,
                        body1 = 50.0,
                        target_mass = 100 * kroupa_mean_mass(0.1, 50.0),
                    ),
                ),
            ],
            "kepler",
            OrbitSpec(apocentre = 10.0, eccentricity = 0.5),
            MergerOutputSpec(format = "nbody", truncate_jacobi = true, output_dir = out_dir),
        )
        result = generate_merger_ic(cfg; rng = rng)
        @test result isa MergerICResult
        @test isdir(result.output_dir)
        @test isfile(joinpath(result.output_dir, "dat.10"))
        @test isfile(joinpath(result.output_dir, "merger.inp"))
        @test isfile(joinpath(result.output_dir, "merger_summary.txt"))
        @test result.N_total > 0
        @test result.M_total > 0
        @test result.orbit_mode == "kepler"
        @test length(result.cluster_ranges) == 2

        # Verify dat.10 mass normalisation
        lines = readlines(joinpath(result.output_dir, "dat.10"))
        masses = [parse(Float64, split(l)[1]) for l in lines]
        @test sum(masses) ≈ 1.0 atol = 1e-10

        # Verify .inp has KZ(22)=2 and the resolved integration parameters
        inp_text = read(joinpath(result.output_dir, "merger.inp"), String)
        @test occursin("KZ(21:30)=0 2 2 0 0 1 0 0 0 1", inp_text)
        @test occursin("QE=2.000E-04", inp_text) && !occursin("RS0=0.5,", inp_text)
        ic_meta = Nbody6Dynamics.TOML.parsefile(joinpath(result.output_dir, "merger_ic.toml"))
        @test 0 < ic_meta["nbody6"]["rs0"] ≤ minimum(c.rbar for c in cfg.clusters) / result.rbar
        @test ic_meta["nbody6"]["rmin"] > 0 && ic_meta["nbody6"]["dtmin"] > 0
        @test isfinite(ic_meta["meta"]["q_virial"]) && ic_meta["meta"]["q_virial"] > 0
        @test ic_meta["meta"]["rbar_over_rhm_min"] > 1
        @test ic_meta["meta"]["t_cr_config_nb"] > ic_meta["meta"]["t_cr_member_min_nb"] > 0
        @test ic_meta["meta"]["t_star_myr"] > 0
        @test ic_meta["stellar"]["level"] == "C" && ic_meta["nbody6"]["tcrtp0"] == 3600.0
        @test ic_meta["nbody6"]["kz"] == Dict{String,Any}()
        summary_text = read(joinpath(result.output_dir, "merger_summary.txt"), String)
        @test occursin("virial ratio Q = T/|W|", summary_text)
        @test occursin("crossing time: configuration", summary_text)
        @test occursin("Integration (merger.inp", summary_text)
        @test occursin("Stellar evolution: KZ(19)=3 Level=C", summary_text)
        @test occursin("ZMET=0.001,EPOCH0=0,DTPLOT=1 /", inp_text)
        @test occursin("BODY1=50,BODYN=0.1,", inp_text)
    end

    # --- Full pipeline — explicit mode (3 clusters) ---
    @testset "Full merger IC pipeline — explicit 3-cluster" begin
        out_dir = mktempdir()
        cfg = MergerConfig(
            [
                ClusterSpec(
                    profile = PlummerProfile(),
                    N = 80,
                    rbar = 1.0,
                    imf = RescaledKroupaIMF(
                        bodyn = 0.1,
                        body1 = 50.0,
                        target_mass = 80 * kroupa_mean_mass(0.1, 50.0),
                    ),
                    position = [-5.0, 0.0, 0.0],
                    velocity = [1.0, 0.0, 0.0],
                ),
                ClusterSpec(
                    profile = KingProfile(W0 = 5.0),
                    N = 80,
                    rbar = 1.0,
                    imf = RescaledKroupaIMF(
                        bodyn = 0.1,
                        body1 = 50.0,
                        target_mass = 80 * kroupa_mean_mass(0.1, 50.0),
                    ),
                    position = [2.5, 4.33, 0.0],
                    velocity = [-0.5, -0.87, 0.0],
                ),
                ClusterSpec(
                    profile = PlummerProfile(),
                    N = 60,
                    rbar = 0.8,
                    imf = EqualMassIMF(particle_mass = 400.0 / 60),
                    position = [2.5, -4.33, 0.0],
                    velocity = [-0.5, 0.87, 0.0],
                ),
            ],
            "explicit",
            OrbitSpec(),  # ignored in explicit mode
            MergerOutputSpec(format = "nbody", truncate_jacobi = false, output_dir = out_dir),
        )
        result = generate_merger_ic(cfg; rng = rng)
        @test result isa MergerICResult
        @test result.orbit_mode == "explicit"
        @test length(result.cluster_ranges) == 3
        @test result.N_total == 80 + 80 + 60
        # Total mass should sum to 1 in N-body units
        lines = readlines(joinpath(result.output_dir, "dat.10"))
        @test sum(parse(Float64, split(l)[1]) for l in lines) ≈ 1.0 atol = 1e-10
    end

    # --- Equal-mass IMF ---
    @testset "Equal mass IMF" begin
        out_dir = mktempdir()
        cfg = MergerConfig(
            [
                ClusterSpec(
                    profile = PlummerProfile(),
                    N = 50,
                    rbar = 1.0,
                    imf = EqualMassIMF(particle_mass = 500.0 / 50),
                ),
                ClusterSpec(
                    profile = PlummerProfile(),
                    N = 50,
                    rbar = 1.0,
                    imf = EqualMassIMF(particle_mass = 500.0 / 50),
                ),
            ],
            "kepler",
            OrbitSpec(apocentre = 8.0, eccentricity = 0.3),
            MergerOutputSpec(format = "nbody", truncate_jacobi = false, output_dir = out_dir),
        )
        generate_merger_ic(cfg; rng = rng)
        lines = readlines(joinpath(out_dir, "dat.10"))
        @test length(lines) == 100
    end

    # --- Primordial binaries ---
    @testset "Primordial binaries" begin
        rng_b = StableRNG(777)
        m = sample_kroupa(1000; rng = rng_b)
        b = sample_binaries(BinarySpec(; fraction = 0.3), m, rng_b)
        n_b = round(Int, 0.3 * 1000 / 1.3)
        @test length(b.primary) == n_b == length(b.a_pc) == length(b.e)
        @test isempty(intersect(b.primary, b.secondary)) && allunique(vcat(b.primary, b.secondary))
        @test all(b.m1 .≥ b.m2)
        @test all(0 .≤ b.e .< 1) && 0.55 < sum(b.e) / n_b < 0.78            # thermal mean 2/3
        @test all(5e-8 .< b.a_pc .< 0.2)                                      # 0.01 AU … 4×10⁴ AU
        bc = sample_binaries(
            BinarySpec(;
                fraction = 0.5,
                period = "loguniform",
                a_min = 1.0,
                a_max = 10.0,
                eccentricity = "circular",
            ),
            m,
            rng_b,
        )
        au = Nbody6Dynamics._AU_PC
        @test all(bc.e .== 0) && all(au .≤ bc.a_pc .≤ 10au)
        @test length(bc.primary) == round(Int, 0.5 * 1000 / 1.5)
        bq = sample_binaries(
            BinarySpec(; fraction = 0.2, pairing = "uniform_q", q_min = 0.5),
            m,
            rng_b,
        )
        @test all(0.5 .≤ bq.m2 ./ bq.m1 .≤ 1.0)
        @test isempty(sample_binaries(BinarySpec(), m, rng_b).primary)
        @test_throws ErrorException Nbody6Dynamics._validate_binaries(
            BinarySpec(; fraction = 1.0),
            "c",
        )
        @test_throws ErrorException Nbody6Dynamics._validate_binaries(
            BinarySpec(; pairing = "x"),
            "c",
        )
        @test_throws ErrorException Nbody6Dynamics._validate_binaries(
            BinarySpec(; a_min = 5.0, a_max = 1.0),
            "c",
        )
        # Kroupa (1995) periods stay within the distribution's support
        lp = [Nbody6Dynamics._sample_log_period_kroupa1995(rng_b) for _ in 1:2000]
        @test all(1.0 .≤ lp .≤ 8.43) && 3.0 < sum(lp) / length(lp) < 6.0

        # Kepler relative orbit: energy −GM/(2a) and r within [a(1−e), a(1+e)]
        for e_test in (0.0, 0.5, 0.9)
            r, v = Nbody6Dynamics._kepler_relative_orbit(2.0, 1e-3, e_test, rng_b)
            rn = sqrt(sum(r .^ 2))
            @test isapprox(0.5 * sum(v .^ 2) - 2.0 / rn, -2.0 / (2 * 1e-3); rtol = 1e-8)
            @test 1e-3 * (1 - e_test) - 1e-12 ≤ rn ≤ 1e-3 * (1 + e_test) + 1e-12
        end
        R = Nbody6Dynamics._random_rotation(rng_b)
        @test isapprox(R * R', [1 0 0; 0 1 0; 0 0 1]; atol = 1e-12)

        # Expansion of a tiny system set: pairs first, then singles; centres of mass kept
        pos_s = [0.0 1.0 2.0 10.0 11.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]
        vel_s = zeros(3, 5)
        mass_s = [2.0, 1.0, 1.5, 1.0, 1.0]
        cb1 = (
            system_binary = [1, 0, 2],
            m1 = [1.2, 0.9],
            m2 = [0.8, 0.6],
            a_pc = [1e-3, 2e-3],
            e = [0.0, 0.5],
        )
        cb2 = (
            system_binary = [0, 0],
            m1 = Float64[],
            m2 = Float64[],
            a_pc = Float64[],
            e = Float64[],
        )
        ex = expand_binaries(
            pos_s,
            vel_s,
            mass_s,
            [1:3, 4:5],
            [cb1, cb2],
            [[1, 2, 3], [1, 2]];
            rng = rng_b,
        )
        @test ex.n_pairs == [2, 0] && length(ex.mass) == 7
        @test ex.mass[1:4] ≈ [1.2, 0.8, 0.9, 0.6] && ex.mass[5:7] ≈ [1.0, 1.0, 1.0]
        @test ex.cluster_blocks[1] == [1:4, 5:5] && ex.cluster_blocks[2] == [5:4, 6:7]
        for (k, sys) in ((1, 1), (2, 3))
            i1, i2 = 2k - 1, 2k
            M = ex.mass[i1] + ex.mass[i2]
            com = (ex.mass[i1] .* ex.pos[:, i1] .+ ex.mass[i2] .* ex.pos[:, i2]) ./ M
            @test isapprox(com, pos_s[:, sys]; atol = 1e-12)
            r_rel = ex.pos[:, i1] .- ex.pos[:, i2]
            v_rel = ex.vel[:, i1] .- ex.vel[:, i2]
            @test 0.5 * sum(v_rel .^ 2) - M / sqrt(sum(r_rel .^ 2)) < 0
        end
        @test sum(ex.mass) ≈ sum(mass_s)
        @test isnan(ex.hard_fraction[2]) && 0 ≤ ex.hard_fraction[1] ≤ 1
        @test Nbody6Dynamics._members_from_blocks([1:4, 5:5]) == [1, 2, 3, 4, 5]

        # Summary parsing with and without pair counts
        sp = joinpath(TESTDIR, "summary_bin.txt")
        write(
            sp,
            "  Cluster 1: king, imf=kroupa, N=10 (after trunc: 9, binaries: 2), M=1 M☉\n" *
            "  Cluster 2: plummer, imf=kroupa, N=6 (after trunc: 5, binaries: 1), M=1 M☉\n",
        )
        mem = parse_merger_summary(sp)
        @test mem[1] == vcat(1:4, 7:11) && mem[2] == vcat(5:6, 12:14)
        write(
            sp,
            "  Cluster 1: king, imf=kroupa, N=10 (after trunc: 9), M=1 M☉\n" *
            "  Cluster 2: plummer, imf=kroupa, N=6 (after trunc: 5), M=1 M☉\n",
        )
        mem0 = parse_merger_summary(sp)
        @test mem0[1] == collect(1:9) && mem0[2] == collect(10:14)

        # Full generation with a binary-rich cluster: ordering, input, metadata round trip
        cfg_bin = MergerConfig(
            [
                ClusterSpec(
                    profile = PlummerProfile(),
                    N = 200,
                    rbar = 1.0,
                    imf = KroupaIMF(),
                    binaries = BinarySpec(; fraction = 0.3),
                ),
                ClusterSpec(
                    profile = KingProfile(W0 = 5.0),
                    N = 200,
                    rbar = 1.0,
                    imf = KroupaIMF(),
                ),
            ],
            "kepler",
            OrbitSpec(apocentre = 10.0, eccentricity = 0.5),
            MergerOutputSpec(; output_dir = mktempdir());
            seed = 5,
        )
        res_b = generate_merger_ic(cfg_bin)
        @test res_b.n_pairs[2] == 0 && res_b.n_pairs[1] ≥ 30
        nb = sum(res_b.n_pairs)
        inp_b = read(joinpath(res_b.output_dir, "merger.inp"), String)
        @test occursin("NBIN0=$(nb),", inp_b) && occursin("KZ(1:10)=1 -1 2 0 0 0 3 2 0 0", inp_b)
        dat = [parse.(Float64, split(l)) for l in eachline(joinpath(res_b.output_dir, "dat.10"))]
        @test length(dat) == res_b.N_total
        seps = [sqrt(sum((dat[2k - 1][2:4] .- dat[2k][2:4]) .^ 2)) for k in 1:nb]
        @test maximum(seps) < 0.05                          # pairs are far smaller than a cluster
        @test length(res_b.cluster_ranges[1]) + length(res_b.cluster_ranges[2]) == res_b.N_total
        @test res_b.cluster_ranges[1][1:(2nb)] == collect(1:(2nb))
        ic_b = load_merger_ic_result(res_b.output_dir)
        @test ic_b.cluster_ranges == res_b.cluster_ranges && ic_b.n_pairs == res_b.n_pairs
        @test ic_b.cluster_specs[1].binaries.fraction == 0.3
        @test parse_merger_summary(joinpath(res_b.output_dir, "merger_summary.txt")) ==
              res_b.cluster_ranges
        @test occursin(
            "Primordial binaries: NBIN0 = $(nb)",
            read(joinpath(res_b.output_dir, "merger_summary.txt"), String),
        )
        ic_meta_b = Nbody6Dynamics.TOML.parsefile(joinpath(res_b.output_dir, "merger_ic.toml"))
        @test ic_meta_b["meta"]["nbin0"] == nb && length(ic_meta_b["cluster_blocks"][1]) == 2
        # TOML round trip of the binaries table
        bt_path = joinpath(TESTDIR, "merger_bin.toml")
        write(
            bt_path,
            """
[merger]
n_clusters = 2
orbit_mode = "kepler"

[merger.cluster1]
model = "plummer"
N = 100
rbar = 1.0

[merger.cluster1.binaries]
fraction = 0.25
pairing = "uniform_q"
period = "loguniform"
a_min = 0.5
a_max = 50.0
eccentricity = "circular"

[merger.cluster2]
model = "king"
N = 100
rbar = 1.0
""",
        )
        cfg_bt = load_merger_config(bt_path)
        @test cfg_bt.clusters[1].binaries.fraction == 0.25 &&
              cfg_bt.clusters[1].binaries.pairing == "uniform_q"
        @test cfg_bt.clusters[1].binaries.a_max == 50.0 &&
              cfg_bt.clusters[2].binaries.fraction == 0.0
    end

    # --- Summary write→parse round-trip (guards the format/regex coupling) ---
    @testset "Merger summary round-trip" begin
        out_dir = mktempdir()
        cfg = MergerConfig(
            [
                ClusterSpec(
                    profile = KingProfile(W0 = 5.0),
                    N = 150,
                    rbar = 1.0,
                    imf = RescaledKroupaIMF(
                        bodyn = 0.1,
                        body1 = 50.0,
                        target_mass = 150 * kroupa_mean_mass(0.1, 50.0),
                    ),
                ),
                ClusterSpec(
                    profile = PlummerProfile(),
                    N = 100,
                    rbar = 1.0,
                    imf = EqualMassIMF(particle_mass = 8e2 / 100),
                ),
            ],
            "kepler",
            OrbitSpec(apocentre = 10.0, eccentricity = 0.4),
            MergerOutputSpec(format = "nbody", truncate_jacobi = true, output_dir = out_dir),
        )
        result = generate_merger_ic(cfg; rng = rng)
        ranges = parse_merger_summary(joinpath(out_dir, "merger_summary.txt"))
        @test ranges == result.cluster_ranges
        @test sum(length, ranges) == result.N_total
    end

    # --- MergerPipelineConfig in Nbody6Config ---
    @testset "MergerPipelineConfig in load_config" begin
        cfg = load_config(joinpath(TESTDIR, "test_config.toml"))
        @test cfg.merger.enabled == false
        @test cfg.merger.config_file == ""
    end
end

# =====================================================================
