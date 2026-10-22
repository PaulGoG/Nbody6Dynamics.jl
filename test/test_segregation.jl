# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Primordial mass segregation after Baumgardt, De Marchi & Kroupa (2008) with
# the McLuster degree S (Küpper et al. 2011, Eq. A1), acting on the systems
# before virialisation.

@testset "Primordial mass segregation" begin
    using StableRNGs
    nd = Nbody6Dynamics

    @testset "McLuster slot order" begin
        n = 20_000
        # S = 1: perfect order (the i-th most massive takes the i-th slot); S = 0: a random permutation
        @test nd.segregation_order(n, 1.0, StableRNG(1)) == collect(1:n)
        p0 = nd.segregation_order(n, 0.0, StableRNG(2))
        @test sort(p0) == collect(1:n) && p0 != collect(1:n)
        # S = 0.5: the relative slot y = j/(N − i) among the free slots has density 2(1 − y)
        y = nd._segregation_relative_slots(n, 0.5, StableRNG(3))
        F(Y) = 1 - (1 - Y)^2
        s = sort(y)
        @test maximum(abs(F(s[k]) - k / length(s)) for k in eachindex(s)) < 1.63 / sqrt(length(s))
        @test_throws ArgumentError nd.segregation_order(10, 1.5, StableRNG(1))
    end

    @testset "Energy ranking of the systems" begin
        rng = StableRNG(4)
        N = 4000
        mass = sort!(0.1 .+ 10 .* rand(rng, N) .^ 3; rev = true)   # a steep mass function
        for prof in (PlummerProfile(), KingProfile(6.0))
            pos, vel, energy =
                nd.segregate_systems(prof, mass, 1.0; degree = 1.0, rng = StableRNG(5))
            @test size(pos) == (3, N) && size(vel) == (3, N)
            # S = 1: heavier systems never sit at higher specific energy than lighter ones
            @test issorted(energy)   # energy in the order of `mass` (descending mass)
            # BDK08 keeps the mass-density profile: the mass-weighted radial distribution is the same at S = 1
            # and S = 0 (two-sample Kolmogorov distance with the effective number (Σm)²/Σm²)
            pos0, _, _ = nd.segregate_systems(prof, mass, 1.0; degree = 0.0, rng = StableRNG(6))
            r1 = vec(sqrt.(sum(abs2, pos; dims = 1)))
            r0 = vec(sqrt.(sum(abs2, pos0; dims = 1)))
            M = sum(mass)
            Mcum(r, x) = sum(mass[i] for i in eachindex(r) if r[i] ≤ x; init = 0.0) / M
            grid = sort(vcat(r0, r1))[1:10:end]
            D = maximum(abs(Mcum(r1, x) - Mcum(r0, x)) for x in grid)
            n_eff = M^2 / sum(abs2, mass)
            @test D < 1.95 * sqrt(2 / n_eff)
            # Heavy systems are concentrated at S = 1 and not at S = 0
            rh(p, sel) = sort(vec(sqrt.(sum(abs2, p[:, sel]; dims = 1))))[cld(length(sel), 2)]
            heavy = 1:200
            @test rh(pos, heavy) < 0.7 * rh(pos0, heavy)
        end
    end

    @testset "Generator integration" begin
        toml(seg) = """
[merger]
n_clusters = 1
orbit_mode = "single"
seed = 67

[merger.cluster1]
model = "plummer"
N = 3000
rbar = 1.0
imf = "kroupa"
body1 = 20.0
bodyn = 0.1
$seg

[merger.output]
format = "nbody"
"""
        load(t) = (p = joinpath(mktempdir(), "s.toml"); write(p, t); load_merger_config(p))
        base = load(toml(""))
        seg = load(toml("\n[merger.cluster1.segregation]\nmethod = \"bdk08\"\ndegree = 1.0"))
        @test seg.clusters[1].segregation.method == "bdk08" &&
              seg.clusters[1].segregation.degree == 1.0
        @test base.clusters[1].segregation.method == "none"
        @test_throws ArgumentError load(
            toml("\n[merger.cluster1.segregation]\nmethod = \"bdk08\"\ndegree = 2.0"),
        )
        @test_throws ArgumentError load(toml("\n[merger.cluster1.segregation]\nmethod = \"other\""))
        r0 = generate_merger_ic(base; output_dir = mktempdir())
        r1 = generate_merger_ic(seg; output_dir = mktempdir())
        # Same stars, heavier ones more concentrated
        @test sort(r0.mass_physical) == sort(r1.mass_physical)
        top(r) = partialsortperm(r.mass_physical, 1:30; rev = true)
        rad(r, idx) = sum(sqrt.(sum(abs2, r.pos_physical[:, idx]; dims = 1))) / length(idx)
        @test rad(r1, top(r1)) < 0.7 * rad(r0, top(r0))
    end
end
