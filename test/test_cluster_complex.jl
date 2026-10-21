# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Roadmap 6.5: multi-cluster complexes — member mass function, size–mass
# presets, centre distributions and the bulk velocity field scaled to a
# target virial ratio of the members as softened point masses.

@testset "Cluster complexes (6.5)" begin
    using StableRNGs
    nd = Nbody6Dynamics
    G = nd._G_PC_KMS2_MSUN

    @testset "Member masses" begin
        eq = ComplexSpec(; n_members = 10, mass_function = "equal", mass_total = 5000.0)
        m = nd.complex_member_masses(eq, StableRNG(1))
        @test length(m) == 10 && all(==(500.0), m)
        pl = ComplexSpec(;
            n_members = 4000,
            mass_function = "power_law",
            beta = 2.0,
            mass_lo = 100.0,
            mass_hi = 1.0e4,
        )
        m = nd.complex_member_masses(pl, StableRNG(2))
        @test length(m) == 4000 && all(x -> 100.0 ≤ x ≤ 1.0e4, m)
        # β = 2: the CDF F(M) = (1/M_lo − 1/M)/(1/M_lo − 1/M_hi); Kolmogorov distance small
        F(x) = (1 / 100 - 1 / x) / (1 / 100 - 1 / 1.0e4)
        s = sort(m)
        D = maximum(abs(F(s[k]) - k / length(s)) for k in eachindex(s))
        @test D < 1.63 / sqrt(length(s))   # 1 % level of the Kolmogorov distribution
        # β = 1 and a mass target instead of a count
        p1 = ComplexSpec(;
            n_members = 0,
            mass_total = 2.0e4,
            mass_function = "power_law",
            beta = 1.0,
            mass_lo = 100.0,
            mass_hi = 1000.0,
        )
        m1 = nd.complex_member_masses(p1, StableRNG(3))
        @test sum(m1) ≥ 2.0e4 && sum(m1) - m1[end] < 2.0e4
    end

    @testset "Size–mass relations" begin
        r(preset; kw...) = nd.complex_member_radius(ComplexSpec(; size_mass = preset, kw...), 1.0e4)
        @test r("constant"; rh = 0.7) == 0.7
        @test r("marks_kroupa") ≈ 0.10 * 1.0e4^0.13
        @test r("frost") ≈ 2.365 / 8
        @test r("dobbs") ≈ 2.26
        @test nd.complex_member_radius(ComplexSpec(; size_mass = "frost"), 1.0e5) ≈ 0.45 rtol = 0.01
    end

    @testset "Centre distributions" begin
        n = 3000
        sph = ComplexSpec(; n_members = n, centres = "uniform_sphere", radius = 4.0)
        x, v = nd.complex_centres(sph, n, StableRNG(4))
        r = vec(sqrt.(sum(abs2, x; dims = 1)))
        @test size(x) == (3, n) && size(v) == (3, n) && maximum(r) ≤ 4.0
        @test sort(r)[n ÷ 2] ≈ 4.0 * 0.5^(1 / 3) rtol = 0.05
        pl = ComplexSpec(; n_members = n, centres = "plummer", radius = 2.0, plummer_cutoff = 5.0)
        x, v = nd.complex_centres(pl, n, StableRNG(5))
        r = vec(sqrt.(sum(abs2, x; dims = 1)))
        @test maximum(r) ≤ 10.0
        # Plummer enclosed mass M(r)/M = r³/(r² + a²)^{3/2}, renormalised to the cut-off
        Mr(q) = q^3 / (q^2 + 4.0)^1.5 / (10.0^3 / (104.0)^1.5)
        rs = sort(r)
        @test maximum(abs(Mr(rs[k]) - k / n) for k in eachindex(rs)) < 1.63 / sqrt(n)
        fr = ComplexSpec(;
            n_members = 200,
            centres = "fractal",
            radius = 3.0,
            fractal_dimension = 1.6,
        )
        x, v = nd.complex_centres(fr, 200, StableRNG(6))
        @test size(x) == (3, 200) && maximum(sqrt.(sum(abs2, x; dims = 1))) ≤ 3.0
        fil = ComplexSpec(;
            n_members = 500,
            centres = "filament",
            filament_length = 20.0,
            filament_width = 1.0,
        )
        x, v = nd.complex_centres(fil, 500, StableRNG(7))
        @test all(abs.(x[1, :]) .≤ 10.0)
        y = x[2, :]
        @test sqrt(sum(abs2, y .- sum(y) / length(y)) / (length(y) - 1)) ≈ 1.0 rtol = 0.15
    end

    @testset "Velocity scaling" begin
        rng = StableRNG(8)
        n = 30
        pos = 5 .* randn(rng, 3, n)
        vel = randn(rng, 3, n)
        mass = 100 .+ 900 .* rand(rng, n)
        ε = 0.5
        out = nd.complex_scale_velocities!(
            vel,
            pos,
            mass;
            virial_ratio = 0.3,
            rotation_fraction = 0.2,
            softening = ε,
        )
        # Centred, scaled to Q, with the rotational fraction about +z
        M = sum(mass)
        @test all(abs.(vec(sum(mass' .* pos; dims = 2))) ./ M .< 1e-12)
        @test all(abs.(vec(sum(mass' .* vel; dims = 2))) ./ M .< 1e-12)
        T = 0.5 * sum(mass .* vec(sum(abs2, vel; dims = 1)))
        W =
            -sum(
                G * mass[i] * mass[j] / sqrt(sum(abs2, pos[:, i] .- pos[:, j]) + ε^2) for i in 1:n
                for j in (i + 1):n
            )
        @test T / abs(W) ≈ 0.3 rtol = 1e-10
        @test out.virial_ratio ≈ 0.3 && out.e_total ≈ T + W && out.rotation_fraction ≈ 0.2
        L = sum(mass[i] .* [pos[1, i] * vel[2, i] - pos[2, i] * vel[1, i]] for i in 1:n)[1]
        @test L > 0 && out.l_total[3] ≈ L
        # Cold start: Q = 0 gives zero velocities
        vel2 = randn(rng, 3, n)
        nd.complex_scale_velocities!(vel2, copy(pos), mass; virial_ratio = 0.0, softening = ε)
        @test all(==(0.0), vel2)
    end

    @testset "Configuration and generation" begin
        toml(complex) = """
[merger]
orbit_mode = "complex"
seed = 65

[merger.complex]
$complex

[merger.complex.member]
model = "plummer"
imf = "kroupa"
body1 = 1.0
bodyn = 0.1

[merger.output]
format = "nbody"
truncate_jacobi = false
"""
        load(t) = (p = joinpath(mktempdir(), "c.toml"); write(p, t); load_merger_config(p))
        base = "n_members = 6\nmass_function = \"equal\"\nmass_total = 600.0\nsize_mass = \"constant\"\nrh = 0.5\ncentres = \"uniform_sphere\"\nradius = 5.0\nvirial_ratio = 0.4"
        cfg = load(toml(base))
        @test cfg.orbit_mode == "complex" && cfg.complex.n_members == 6
        # A power law takes either a count or a mass target, not both
        @test_throws ArgumentError load(
            toml(replace(base, "mass_function = \"equal\"" => "mass_function = \"power_law\"")),
        )
        @test_throws ArgumentError load(toml(replace(base, "radius = 5.0" => "radius = -5.0")))
        @test_throws ArgumentError load(
            toml(replace(base, "centres = \"uniform_sphere\"" => "centres = \"torus\"")),
        )
        @test_throws ArgumentError load(toml(base * "\nrotation_fraction = 1.0"))
        res = generate_merger_ic(cfg; output_dir = mktempdir())
        @test length(res.cluster_ranges) == 6
        meta = Nbody6Dynamics.TOML.parsefile(joinpath(res.output_dir, "merger_ic.toml"))
        @test meta["orbit_mode"] == "complex"
        @test meta["complex"]["realised"]["n_members"] == 6
        @test meta["complex"]["realised"]["virial_ratio"] ≈ 0.4 rtol = 1e-8
        # The member centres of mass sit where the complex put them, the complex at the origin
        cms = [
            vec(sum(res.mass_physical[r]' .* res.pos_physical[:, r]; dims = 2)) ./
            sum(res.mass_physical[r]) for r in res.cluster_ranges
        ]
        Mtot = sum(res.mass_physical)
        @test all(abs.(vec(sum(res.mass_physical' .* res.pos_physical; dims = 2))) ./ Mtot .< 1e-8)
        @test all(sqrt(sum(abs2, c)) ≤ 10.0 for c in cms)
        # Same seed, same complex
        res2 = generate_merger_ic(cfg; output_dir = mktempdir())
        @test res2.pos_physical == res.pos_physical
    end
end
