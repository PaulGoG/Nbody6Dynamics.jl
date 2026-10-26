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

    moe(; kw...) = nd.MassiveBinarySpec(; model = "moe2017", kw...)

    @testset "Moe & Di Stefano 2017: published integrals" begin
        # Moe & Di Stefano (2017) §9.4 and Table 13; inner-binary fractions of Moe's code
        for (m₁, f) in ((1.0, 0.367), (7.0, 0.633), (28.0, 1.026))
            @test nd.moe2017_integrals(m₁).f_mult_q03 ≈ f atol = 0.01
        end
        for (m₁, f) in ((1.0, 0.500), (7.0, 1.355), (12.0, 1.604), (28.0, 2.113))
            @test nd.moe2017_integrals(m₁).f_mult_q01 ≈ f atol = 0.02
        end
        for (m₁, f) in ((1.0, 0.407), (28.0, 0.955))
            @test nd.moe2017_integrals(m₁).f_bin ≈ f atol = 0.01
        end
        for (lP, f) in ((1.0, 0.299), (3.0, 0.320), (5.0, 0.301), (7.0, 0.177))
            @test nd._moe_companion_frequency(lP, 28.0) * nd._moe_q_factor(lP, 28.0) ≈ f atol = 0.01
        end
    end

    @testset "Moe & Di Stefano 2017: conditional q and e" begin
        lP, m₁ = 3.0, 28.0
        p = nd._moe_q_parameters(lP, m₁)
        rng = StableRNG(5)
        n = 200_000
        q = [nd._moe_draw_q(rand(rng), rand(rng), p, 0.1) for _ in 1:n]
        @test all(0.1 .≤ q .≤ 1.0)
        # Expected fraction above 0.95 from the segment weights
        γl, γs, F = p.γ_large, p.γ_small, p.F_twin
        I_l = (1 - 0.3^(γl + 1)) / (γl + 1)
        W_s = (1 - F) * 0.3^γl / I_l * 0.3 * (1 - 3.0^(-(γs + 1))) / (γs + 1)
        P_hi = (F + (1 - F) * (1 - 0.95^(γl + 1)) / (γl + 1) / I_l) / (1 + W_s)
        @test abs(count(>(0.95), q) / n - P_hi) < 3 * sqrt(P_hi * (1 - P_hi) / n)
        # Slope of the density per unit q on log-uniform bins in 0.3–0.95
        edges = exp10.(range(log10(0.3), log10(0.95); length = 21))
        counts = zeros(Int, 20)
        for x in q
            k = searchsortedlast(edges, x)
            1 ≤ k ≤ 20 && (counts[k] += 1)
        end
        x = log.(edges[1:(end - 1)])
        y = log.(counts ./ diff(edges))
        x̄, ȳ = sum(x) / 20, sum(y) / 20
        slope = sum((x .- x̄) .* (y .- ȳ)) / sum((x .- x̄) .^ 2)
        @test abs(slope - γl) < 0.1

        η = nd._moe_eta(lP, m₁)
        e_max = nd._moe_e_max(lP)
        n_e = 100_000
        e = [nd._moe_draw_e(rand(rng), η, e_max, true) for _ in 1:n_e]
        @test all(0 .≤ e .≤ e_max)
        p_pow = 8 / (9 + η)
        @test abs(count(<(0.8 * e_max), e) / n_e - p_pow) < 3 * sqrt(p_pow * (1 - p_pow) / n_e)
        e = [nd._moe_draw_e(rand(rng), η, e_max, false) for _ in 1:n_e]
        @test abs(-1 - n_e / sum(log.(e ./ e_max)) - η) < 0.05
    end

    @testset "Moe & Di Stefano 2017: pool pairing, guard, determinism" begin
        m = nd.sample_kroupa(20_000; m_low = 0.08, m_up = 100.0, rng = StableRNG(6))
        spec = BinarySpec(;
            fraction = 0.0,
            massive = moe(; m_min = 5.0, pairing = "pool", log_period_max = 8.0),
        )
        b = @test_logs (:warn, r"outside 0.8–40") match_mode = :any nd.sample_binaries(
            spec,
            m,
            StableRNG(7),
        )
        @test sort(b.masses) == sort(m)
        a_au = b.a_pc ./ nd._AU_PC
        @test all(0 .≤ b.e .< 1) && all(b.m2 .≤ b.m1) && all(a_au .> 0)
        # e = 0 only where the law is circular: P ≤ 2 d (e_max = 0) or η ≤ −1
        lP = [nd._log_period_days(a_au[k], b.m1[k] + b.m2[k]) for k in eachindex(a_au)]
        @test all(
            b.e[k] > 0 for k in eachindex(a_au) if
            lP[k] > log10(2.0) + 1e-9 && nd._moe_eta(lP[k], nd._moe_clamp_m1(b.m1[k])) > -1 + 1e-9
        )
        # Binary fraction of the primaries ≥ 20 M☉ reached unpaired
        is_secondary = falses(length(m))
        is_secondary[b.secondary] .= true
        trials = [i for i in eachindex(m) if m[i] ≥ 20.0 && !is_secondary[i]]
        n_bin = count(≥(20.0), b.m1)
        m̄ = sum(m[trials]) / length(trials)
        @test abs(n_bin / length(trials) - nd.moe2017_integrals(nd._moe_clamp_m1(m̄)).f_bin) < 0.05

        # Roche guard: periastron beyond the fill-limited separation of every pair
        φ, z = 0.7, 2e-4
        spec_g = BinarySpec(;
            fraction = 0.0,
            roche_fill_max = φ,
            massive = moe(; m_min = 5.0, pairing = "pool", log_period_max = 8.0),
        )
        g = @test_logs (:warn, r"outside 0.8–40") match_mode = :any nd.sample_binaries(
            spec_g,
            m,
            StableRNG(8);
            zmet = z,
        )
        a_g = g.a_pc ./ nd._AU_PC
        @test !isempty(a_g)
        @test all(
            a_g[k] * (1 - g.e[k]) ≥ nd._roche_separation_au(g.m1[k], g.m2[k], z, φ) * (1 - 1e-9) for
            k in eachindex(a_g)
        )

        # Determinism
        b1 = @test_logs (:warn,) match_mode = :any nd.sample_binaries(spec, m, StableRNG(7))
        b2 = @test_logs (:warn,) match_mode = :any nd.sample_binaries(spec, m, StableRNG(7))
        @test b1.primary == b2.primary && b1.secondary == b2.secondary
        @test b1.a_pc == b2.a_pc && b1.e == b2.e && b1.masses == b2.masses
    end

    @testset "Moe & Di Stefano 2017: configuration" begin
        moe_toml(b) = """
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

[merger.cluster1.binaries.massive]
model = "moe2017"
$b

[merger.stellar]
kz19 = 3
zmet = 0.001

[merger.output]
format = "nbody"
"""
        moe_load(t) = (p = joinpath(mktempdir(), "b.toml"); write(p, t); load_merger_config(p))
        c = moe_load(moe_toml("ecc_turnover = false"))
        @test c.clusters[1].binaries.massive.model == "moe2017" &&
              c.clusters[1].binaries.massive.ecc_turnover == false
        refusal(t) =
            try
                moe_load(t)
                nothing
            catch err
                err
            end
        err = refusal(moe_toml("log_period_min = 0.1"))
        @test err isa ArgumentError && occursin("log_period_min", sprint(showerror, err))
        err = refusal(moe_toml("log_period_max = 8.5"))
        @test err isa ArgumentError && occursin("log_period_max", sprint(showerror, err))
        @test_logs (:warn, r"m_min") match_mode = :any moe_load(moe_toml("m_min = 0.5"))
    end
end
