# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# The two-cluster orbit by pericentre or apocentre (absolute or in units of r_h
# of cluster 1), its orientation, orbit-independent truncation, the
# density-ratio key, per-cluster random streams and the derived orbit
# quantities of the IC record.

@testset "Two-cluster orbit" begin
    using StableRNGs
    nd = Nbody6Dynamics
    G = nd._G_PC_KMS2_MSUN

    @testset "Size keys" begin
        @test nd.orbit_apocentre(OrbitSpec(; apocentre = 12.0, eccentricity = 0.5), 2.0) == 12.0
        @test nd.orbit_apocentre(OrbitSpec(; apocentre_rh = 6.0, eccentricity = 0.5), 2.0) == 12.0
        @test nd.orbit_apocentre(OrbitSpec(; pericentre = 4.0, eccentricity = 0.5), 2.0) ≈ 12.0
        @test nd.orbit_apocentre(OrbitSpec(; pericentre_rh = 2.0, eccentricity = 0.5), 2.0) ≈ 12.0
        @test nd.orbit_apocentre(OrbitSpec(; eccentricity = 0.5), 2.0) == 15.0
        @test_throws ArgumentError nd.orbit_apocentre(OrbitSpec(), 0.0)
    end

    base_toml(
        orbit;
        extra_merger = "",
        extra_cluster2 = "",
        output = "truncate_jacobi = false",
    ) = """
[merger]
n_clusters = 2
orbit_mode = "kepler"
seed = 4242
$extra_merger

[merger.cluster1]
model = "plummer"
N = 300
rbar = 1.0
imf = "kroupa"
body1 = 1.0
bodyn = 0.1

[merger.cluster1.binaries]
fraction = 0.3

[merger.cluster2]
model = "plummer"
N = 300
imf = "kroupa"
body1 = 1.0
bodyn = 0.1
$extra_cluster2

[merger.cluster2.binaries]
fraction = 0.3

[merger.orbit]
$orbit

[merger.output]
format = "nbody"
$output
"""
    function load(toml)
        path = joinpath(mktempdir(), "merger.toml")
        write(path, toml)
        return load_merger_config(path)
    end

    @testset "Configuration" begin
        # One size key only
        @test_throws ArgumentError load(
            base_toml(
                "apocentre = 10.0\npericentre = 2.0\neccentricity = 0.5";
                extra_cluster2 = "rbar = 1.0",
            ),
        )
        @test_throws ArgumentError load(
            base_toml("pericentre_rh = -1.0\neccentricity = 0.5"; extra_cluster2 = "rbar = 1.0"),
        )
        @test_throws ArgumentError load(
            base_toml(
                "apocentre = 10.0\neccentricity = 0.5\norbit_normal = [0.0, 0.0, 0.0]";
                extra_cluster2 = "rbar = 1.0",
            ),
        )
        @test_throws ArgumentError load(
            base_toml(
                "apocentre = 10.0\neccentricity = 0.5";
                extra_merger = "seed_streams = \"other\"",
                extra_cluster2 = "rbar = 1.0",
            ),
        )
        @test_throws ArgumentError load(
            base_toml(
                "apocentre = 10.0\neccentricity = 0.5";
                extra_cluster2 = "rbar = 1.0",
                output = "truncate_radius_rh = -1.0",
            ),
        )
        c = load(
            base_toml("pericentre_rh = 2.0\neccentricity = 0.5"; extra_cluster2 = "rbar = 1.0"),
        )
        @test c.orbit.pericentre_rh == 2.0 && c.seed_streams == "shared"
        # Density ratio: equal masses, ρ₂/ρ₁ = 8 → r_h,2 = r_h,1 / 2
        d = load(
            base_toml(
                "apocentre = 10.0\neccentricity = 0.5";
                extra_cluster2 = "density_ratio = 8.0",
            ),
        )
        @test d.clusters[2].rbar ≈ 0.5
        @test_throws ArgumentError load(
            base_toml(
                "apocentre = 10.0\neccentricity = 0.5";
                extra_cluster2 = "rbar = 1.0\ndensity_ratio = 8.0",
            ),
        )
    end

    @testset "Placement and orientation" begin
        rng = StableRNG(613)
        p1 = randn(rng, 3, 40)
        v1 = randn(rng, 3, 40)
        m1 = rand(rng, 40) .+ 0.5
        p2 = randn(rng, 3, 30)
        v2 = randn(rng, 3, 30)
        m2 = rand(rng, 30) .+ 0.5
        # Each set about its own centre of mass, as the generator samples them
        for (p, v, m) in ((p1, v1, m1), (p2, v2, m2))
            p .-= sum(m' .* p; dims = 2) ./ sum(m)
            v .-= sum(m' .* v; dims = 2) ./ sum(m)
        end
        d_apo, e = 20.0, 0.5
        pos, vel, mass, ranges, _ = nd.setup_two_cluster_orbit(
            p1,
            v1,
            m1,
            p2,
            v2,
            m2,
            d_apo,
            e;
            truncate_jacobi_flag = false,
        )
        # The layout of earlier versions, bit for bit: cluster 1 at −d1 x̂ moving +v1 ŷ, cluster 2 at +d2 x̂ moving −v2 ŷ
        M1, M2 = sum(m1), sum(m2)
        dd1 = (M2 / (M1 + M2)) * d_apo
        dd2 = (M1 / (M1 + M2)) * d_apo
        w1, w2 = kepler_velocity(M1, M2, d_apo, e)
        # Largest deviation per component, so that a failure reports its size
        maxdev(a, b) = maximum(abs, a .- b)
        @test (
            maxdev(pos[1, 1:40], p1[1, :] .- dd1),
            maxdev(pos[2, 1:40], p1[2, :]),
            maxdev(pos[3, 1:40], p1[3, :]),
        ) == (0.0, 0.0, 0.0)
        @test (
            maxdev(vel[1, 1:40], v1[1, :]),
            maxdev(vel[2, 1:40], v1[2, :] .+ w1),
            maxdev(vel[3, 1:40], v1[3, :]),
        ) == (0.0, 0.0, 0.0)
        @test (maxdev(pos[1, 41:70], p2[1, :] .+ dd2), maxdev(vel[2, 41:70], v2[2, :] .- w2)) ==
              (0.0, 0.0)
        # Orbital angular momentum of the two centres along the orbit normal
        function centres(pos, vel, mass, ranges)
            c = [vec(sum(mass[r]' .* pos[:, r]; dims = 2)) ./ sum(mass[r]) for r in ranges]
            u = [vec(sum(mass[r]' .* vel[:, r]; dims = 2)) ./ sum(mass[r]) for r in ranges]
            return c, u
        end
        cross3(a, b) = [a[2]b[3] - a[3]b[2], a[3]b[1] - a[1]b[3], a[1]b[2] - a[2]b[1]]
        for n in ([0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 1.0])
            ps, vs, ms, rs, _ = nd.setup_two_cluster_orbit(
                p1,
                v1,
                m1,
                p2,
                v2,
                m2,
                d_apo,
                e;
                truncate_jacobi_flag = false,
                orbit_normal = n,
                apocentre_angle = 30.0,
            )
            c, u = centres(ps, vs, ms, rs)
            L = cross3(c[2] - c[1], u[2] - u[1])
            @test L ./ sqrt(sum(abs2, L)) ≈ n ./ sqrt(sum(abs2, n)) atol = 1e-10
            @test sqrt(sum(abs2, c[2] - c[1])) ≈ d_apo rtol = 1e-10
            # Apocentre: the relative velocity is perpendicular to the separation
            @test abs(sum((c[2] - c[1]) .* (u[2] - u[1]))) < 1e-9
        end
        # Apocentre angle 90° with the default normal turns the separation from +x to −y
        ps, vs, ms, rs, _ = nd.setup_two_cluster_orbit(
            p1,
            v1,
            m1,
            p2,
            v2,
            m2,
            d_apo,
            e;
            truncate_jacobi_flag = false,
            apocentre_angle = 90.0,
        )
        c, _ = centres(ps, vs, ms, rs)
        @test (c[2] - c[1]) ./ d_apo ≈ [0.0, -1.0, 0.0] atol = 1e-10
        # Fixed truncation radii replace the orbit-dependent Jacobi radius
        _, _, mt, rt, kept = nd.setup_two_cluster_orbit(
            p1,
            v1,
            m1,
            p2,
            v2,
            m2,
            d_apo,
            e;
            truncation_radii = (1.0, 1.5),
        )
        @test kept[1] == findall(i -> sqrt(sum(abs2, p1[:, i])) ≤ 1.0, 1:40)
        @test kept[2] == findall(i -> sqrt(sum(abs2, p2[:, i])) ≤ 1.5, 1:30)
    end

    @testset "Derived orbit quantities" begin
        q = nd.two_cluster_orbit_quantities(1000.0, 1000.0, 12.0, 0.5, 2.0; sigma1 = 2.0)
        @test q.a_pc ≈ 8.0 && q.pericentre_pc ≈ 4.0 && q.apocentre_pc == 12.0
        @test q.apocentre_rh ≈ 6.0 && q.pericentre_rh ≈ 2.0
        @test q.j_orb ≈ sqrt(3.0)
        @test q.period_myr ≈ 2π * sqrt(8.0^3 / (G * 2000.0)) * nd._PC_KM / nd._MYR_S
        @test q.e_orb ≈ -G * 1000.0^2 / 16.0
        @test q.l_orb ≈ 500.0 * sqrt(G * 2000.0 * 6.0)
        @test q.v_peri_kms ≈ sqrt(G * 2000.0 * 1.5 / 4.0)
        @test q.v_peri_over_sigma ≈ q.v_peri_kms / 2.0
        rl = 12.0 * 0.49 / (0.6 + log(2.0))
        @test q.roche_lobe_pc[1] ≈ rl && q.roche_lobe_pc[2] ≈ rl
        @test isnan(q.lambda_peebles) && all(isnan, q.roche_filling)
        qe = nd.two_cluster_orbit_quantities(1000.0, 1000.0, 12.0, 0.5, 2.0; e_total = -50.0)
        @test qe.lambda_peebles ≈ qe.l_orb * sqrt(50.0) / (G * 2000.0^2.5)
        @test_throws ArgumentError nd.two_cluster_orbit_quantities(1000.0, 1000.0, 12.0, 1.0, 2.0)
        @test_throws ArgumentError nd.two_cluster_orbit_quantities(0.0, 1000.0, 12.0, 0.5, 2.0)
    end

    @testset "Common random numbers across an orbit grid" begin
        function cluster1(orbit, streams, output)
            cfg = load(
                base_toml(
                    orbit;
                    extra_merger = "seed_streams = \"$streams\"",
                    extra_cluster2 = "rbar = 1.0",
                    output = output,
                ),
            )
            res = generate_merger_ic(cfg; output_dir = mktempdir())
            r = res.cluster_ranges[1]
            m = res.mass_physical[r]
            x = res.pos_physical[:, r]
            return m, x .- sum(m' .* x; dims = 2) ./ sum(m), res
        end
        trunc = "truncate_jacobi = true\ntruncate_radius_rh = 8.0"
        ma, xa, ra = cluster1("apocentre = 12.0\neccentricity = 0.5", "per_cluster", trunc)
        mb, xb, _ = cluster1("pericentre_rh = 1.0\neccentricity = 0.9", "per_cluster", trunc)
        @test ma == mb
        @test xa ≈ xb atol = 1e-12
        meta = Nbody6Dynamics.TOML.parsefile(joinpath(ra.output_dir, "merger_ic.toml"))
        @test meta["meta"]["schema_version"] == 4 && meta["meta"]["seed_streams"] == "per_cluster"
        @test meta["orbit"]["apocentre_resolved"] == 12.0
        @test meta["orbit"]["j_orb"] ≈ sqrt(12.0 / 1.5 * (1 - 0.25) / 1.0)
        @test haskey(meta["cluster1"], "mass_truncated")
        @test meta["meta"]["fp32_offset_bound"] > 0
    end

    @testset "Tidal frame and the pair Jacobi radius" begin
        t1 = TidalSpec(; kz14 = 1)
        @test nd.tidal_frame_rotation(t1) ≈ 0.0264
        @test nd.tidal_frame_rotation(TidalSpec()) == 0.0
        t2 = TidalSpec(; kz14 = 2, gmg = 1.0e11, rg0 = 8.0)
        Ω2 = sqrt(G * 1.0e11 / 8000.0^3)
        @test nd.tidal_frame_rotation(t2) ≈ Ω2
        @test nd.pair_jacobi_radius(6.0e4, t1) ≈ cbrt(G * 6.0e4 / (4 * 0.0144 * 0.0264))
        @test nd.pair_jacobi_radius(6.0e4, t1) ≈ 1.4142 * cbrt(6.0e4) rtol = 1e-4
        @test nd.pair_jacobi_radius(6.0e4, t2) ≈ cbrt(G * 6.0e4 / (3Ω2^2))
        @test nd.pair_jacobi_radius(6.0e4, TidalSpec()) == Inf
        pos = [1.0 0.0; 0.0 2.0; 0.0 0.0]
        vel = zeros(3, 2)
        nd.to_rotating_frame!(pos, vel, 0.5)
        @test vel == [0.0 1.0; -0.5 0.0; 0.0 0.0]

        tidal(frame) = "\n[merger.tidal]\nkz14 = 1\nvelocity_frame = \"$frame\"\n"
        # The frame key needs a rotating-frame tide
        @test_throws ArgumentError load(
            base_toml("apocentre = 10.0\neccentricity = 0.5"; extra_cluster2 = "rbar = 1.0") *
            "\n[merger.tidal]\nkz14 = 0\nvelocity_frame = \"inertial\"\n",
        )
        gen(frame; orbit = "apocentre = 3.0\neccentricity = 0.5") = generate_merger_ic(
            load(
                base_toml(
                    orbit;
                    extra_cluster2 = "rbar = 1.0",
                    extra_merger = "seed_streams = \"per_cluster\"",
                ) *
                tidal(frame) *
                "\n[merger.nbody6]\nqe = 0.01\n",
            );
            output_dir = mktempdir(),
        )
        rot = gen("rotating")
        ine = gen("inertial")
        Ω = nd.tidal_frame_rotation(t1)
        # Same bodies, velocities differ by −Ω × r
        @test rot.pos_physical == ine.pos_physical
        @test ine.vel_physical[1, :] ≈ rot.vel_physical[1, :] .+ Ω .* rot.pos_physical[2, :]
        @test ine.vel_physical[2, :] ≈ rot.vel_physical[2, :] .- Ω .* rot.pos_physical[1, :]
        @test ine.vel_physical[3, :] == rot.vel_physical[3, :]
        meta = Nbody6Dynamics.TOML.parsefile(joinpath(ine.output_dir, "merger_ic.toml"))
        @test meta["tidal"]["velocity_frame"] == "inertial" && meta["tidal"]["omega_kms_pc"] ≈ Ω
        @test 0 < meta["orbit"]["apocentre_over_jacobi_pair"] < 0.5
        # An apocentre beyond the pair Jacobi radius is refused
        @test_throws ErrorException gen("rotating"; orbit = "apocentre = 500.0\neccentricity = 0.5")
    end

    @testset "Control derivation with the new size keys" begin
        ca(d) = nd._control_apocentre(d, 2.0)
        @test ca(Dict("apocentre" => 12.0)) == 12.0
        @test ca(Dict("apocentre_rh" => 6.0)) == 12.0
        @test ca(Dict("pericentre" => 4.0, "eccentricity" => 0.5)) ≈ 12.0
        @test ca(Dict("pericentre_rh" => 2.0, "eccentricity" => 0.5)) ≈ 12.0
        @test ca(Dict{String,Any}()) == 15.0
    end
end
