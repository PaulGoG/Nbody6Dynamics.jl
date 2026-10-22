# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Massive binaries: the Sana et al. (2012) branch for massive primaries with
# pool or replace pairing, and the Roche-lobe guard at periastron with the
# engine's ZAMS radii (Tout et al. 1996).

@testset "Massive binaries and the Roche guard" begin
    using StableRNGs
    nd = Nbody6Dynamics

    @testset "ZAMS radii (engine coefficients)" begin
        # Reference values from the engine's own zcnsts.f / zfuncs.f (rzamsf), compiled standalone
        ref = [
            (2.0000000000000001E-04, 1.0000000000000001E-01, 1.1469037082677808E-01),
            (2.0000000000000001E-04, 5.0000000000000000E-01, 4.3581814999117341E-01),
            (2.0000000000000001E-04, 1.0000000000000000E+00, 8.3707202033584849E-01),
            (2.0000000000000001E-04, 2.0000000000000000E+00, 1.0186539110758197E+00),
            (2.0000000000000001E-04, 5.0000000000000000E+00, 1.7275602203763625E+00),
            (2.0000000000000001E-04, 1.0000000000000000E+01, 2.6504620066683291E+00),
            (2.0000000000000001E-04, 3.0000000000000000E+01, 5.1634098178386765E+00),
            (2.0000000000000001E-04, 6.0000000000000000E+01, 7.5319258661551123E+00),
            (2.0000000000000001E-04, 1.0000000000000000E+02, 9.5296096025000594E+00),
            (1.0000000000000000E-03, 1.0000000000000001E-01, 1.1138780880542040E-01),
            (1.0000000000000000E-03, 5.0000000000000000E-01, 4.4064370035247707E-01),
            (1.0000000000000000E-03, 1.0000000000000000E+00, 8.9453836505900786E-01),
            (1.0000000000000000E-03, 2.0000000000000000E+00, 1.1228456854898075E+00),
            (1.0000000000000000E-03, 5.0000000000000000E+00, 1.9095413377787960E+00),
            (1.0000000000000000E-03, 1.0000000000000000E+01, 2.9150558969563822E+00),
            (1.0000000000000000E-03, 3.0000000000000000E+01, 5.6724325652131542E+00),
            (1.0000000000000000E-03, 6.0000000000000000E+01, 8.3355154979313895E+00),
            (1.0000000000000000E-03, 1.0000000000000000E+02, 1.0652625822867837E+01),
            (2.0000000000000000E-02, 1.0000000000000001E-01, 1.3045599411646885E-01),
            (2.0000000000000000E-02, 5.0000000000000000E-01, 4.5875776207476171E-01),
            (2.0000000000000000E-02, 1.0000000000000000E+00, 8.8824945029751212E-01),
            (2.0000000000000000E-02, 2.0000000000000000E+00, 1.6102103854263528E+00),
            (2.0000000000000000E-02, 5.0000000000000000E+00, 2.6355930591542558E+00),
            (2.0000000000000000E-02, 1.0000000000000000E+01, 3.9419056202596461E+00),
            (2.0000000000000000E-02, 3.0000000000000000E+01, 7.7402760346714334E+00),
            (2.0000000000000000E-02, 6.0000000000000000E+01, 1.2171932202439724E+01),
            (2.0000000000000000E-02, 1.0000000000000000E+02, 1.7226115851833896E+01),
        ]
        for (z, m, r) in ref
            @test nd.zams_radius(m, z) ≈ r rtol = 1e-12
        end
        @test_throws ArgumentError nd.zams_radius(1.0, 0.05)
        @test_throws ArgumentError nd.zams_radius(-1.0, 0.02)
    end

    massive(; kw...) = nd.MassiveBinarySpec(; model = "sana2012", kw...)

    @testset "Sana 2012 marginals" begin
        n = 200_000
        lP, q, e = nd._sana_draw(massive(), n, StableRNG(1))
        Fl(x) = (x^0.45 - 0.15^0.45) / (3.5^0.45 - 0.15^0.45)
        Fq(x) = (x^0.9 - 0.1^0.9) / (1 - 0.1^0.9)
        Fe(x) = (x / 0.9)^0.55
        for (v, F, lo, hi) in ((lP, Fl, 0.15, 3.5), (q, Fq, 0.1, 1.0), (e, Fe, 0.0, 0.9))
            s = sort(v)
            @test all(lo .≤ s .≤ hi)
            @test maximum(abs(F(s[k]) - k / n) for k in eachindex(s)) < 1.63 / sqrt(n)
        end
    end

    @testset "Pairing and guard" begin
        rng = StableRNG(2)
        N = 6000
        m = vcat(5 .+ 45 .* rand(rng, 300), 0.1 .+ 0.9 .* rand(rng, N - 300))   # 300 massive stars
        for pairing in ("pool", "replace")
            spec = BinarySpec(;
                fraction = 0.0,
                massive = massive(; fraction = 0.69, pairing = pairing),
                roche_fill_max = 0.7,
            )
            b = nd.sample_binaries(spec, m, StableRNG(3); zmet = 0.02)
            np = length(b.primary)
            @test all(b.m1 .≥ 5.0) && all(b.m2 .≤ b.m1)
            # Binary fraction of the massive stars reached as unused primaries, binomial within 4σ
            # (a massive star taken as a partner earlier is not a trial)
            trials = 300 - count(m[b.secondary] .≥ 5.0)
            @test abs(np / trials - 0.69) < 4 * sqrt(0.69 * 0.31 / trials)
            if pairing == "pool"
                @test sort(b.masses) == sort(m)      # the IMF multiset is untouched
            else
                @test sort(b.masses) != sort(m)
            end
            # Periastron outside the Roche lobe of both components by the fill limit
            au_rsun = 215.032
            for k in 1:np
                a_rsun = b.a_pc[k] / nd._AU_PC * au_rsun
                for (mi, mo) in ((b.m1[k], b.m2[k]), (b.m2[k], b.m1[k]))
                    qi = mi / mo
                    rl = 0.49 * qi^(2 / 3) / (0.6 * qi^(2 / 3) + log(1 + cbrt(qi)))
                    @test nd.zams_radius(mi, 0.02) / (rl * a_rsun * (1 - b.e[k])) ≤ 0.7 + 1e-12
                end
            end
        end
        # Without the guard and the massive branch: the earlier sampler, draw for draw
        old = nd.sample_binaries(BinarySpec(; fraction = 0.3), m, StableRNG(4))
        new = nd.sample_binaries(BinarySpec(; fraction = 0.3), m, StableRNG(4); zmet = 0.02)
        @test old.a_pc == new.a_pc && old.primary == new.primary
    end

    @testset "Configuration" begin
        toml(b) = """
[merger]
n_clusters = 1
orbit_mode = "single"
seed = 68

[merger.cluster1]
model = "plummer"
N = 2000
rbar = 1.0
imf = "kroupa"
body1 = 50.0
bodyn = 0.1

[merger.cluster1.binaries]
fraction = 0.1
$b

[merger.stellar]
kz19 = 3
zmet = 0.001

[merger.output]
format = "nbody"
"""
        load(t) = (p = joinpath(mktempdir(), "b.toml"); write(p, t); load_merger_config(p))
        c = load(
            toml(
                "roche_fill_max = 0.7\n\n[merger.cluster1.binaries.massive]\nmodel = \"sana2012\"",
            ),
        )
        @test c.clusters[1].binaries.massive.model == "sana2012" &&
              c.clusters[1].binaries.roche_fill_max == 0.7
        @test_throws ArgumentError load(toml("roche_fill_max = 1.5"))
        @test_throws ArgumentError load(
            toml("\n[merger.cluster1.binaries.massive]\nmodel = \"other\""),
        )
        @test_throws ArgumentError load(
            toml("\n[merger.cluster1.binaries.massive]\nmodel = \"sana2012\"\necc_exponent = -1.0"),
        )
        @test_throws ArgumentError load(
            toml(
                "\n[merger.cluster1.binaries.massive]\nmodel = \"sana2012\"\nlog_period_min = 4.0",
            ),
        )
        res = generate_merger_ic(c; output_dir = mktempdir())
        @test res.n_pairs[1] > 0
    end
end
