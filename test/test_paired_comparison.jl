# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Paired comparison of sweep arms (roadmap 4.4): Student-t quantiles and tail
# probabilities against tabulated values, Holm's correction, the paired shift with
# its equivalence and difference tests on synthetic arms (an injected shift found,
# a null arm equivalent, the family-wise error of many null arms bounded), and the
# window mean of a difference of two series.

using StableRNGs: StableRNG

@testset "Paired comparison of arms" begin
    nd = Nbody6Dynamics

    @testset "Student t" begin
        # Tabulated two-sided 95 % and one-sided 95 % quantiles
        @test isapprox(nd._student_t_quantile(0.975, 10), 2.228138852; atol = 1e-8)
        @test isapprox(nd._student_t_quantile(0.95, 4), 2.131846786; atol = 1e-8)
        @test isapprox(nd._student_t_quantile(0.975, 1), 12.70620474; atol = 1e-6)
        @test isapprox(nd._student_t_quantile(0.025, 10), -2.228138852; atol = 1e-8)
        @test nd._student_t_quantile(0.5, 7) == 0.0
        for ν in (1, 3, 9, 30), p in (0.01, 0.2, 0.7, 0.99)
            @test isapprox(nd._student_t_cdf(nd._student_t_quantile(p, ν), ν), p; atol = 1e-10)
        end
        @test nd._student_t_cdf(0.0, 5) == 0.5
        @test isapprox(nd._student_t_cdf(-1.3, 6) + nd._student_t_cdf(1.3, 6), 1.0; atol = 1e-14)
        @test_throws ArgumentError nd._student_t_quantile(1.0, 3)
    end

    @testset "Holm's correction" begin
        @test nd.holm_adjust([0.01, 0.04, 0.03]) ≈ [0.03, 0.06, 0.06]
        @test nd.holm_adjust([0.5, NaN, 0.01]) ≈ [0.5, NaN, 0.02] nans = true
        @test nd.holm_adjust([0.9, 0.8]) == [1.0, 1.0]
    end

    @testset "paired shift on synthetic arms" begin
        rng = StableRNG(4401)
        n = 8
        σ_floor = 0.01                                   # chaos floor of the window means
        null = 0.01 .* randn(rng, n)
        shifted = 0.05 .+ 0.01 .* randn(rng, n)          # five floors
        r0 = nd.paired_shift(null; margin = 3σ_floor)
        r1 = nd.paired_shift(shifted; margin = 3σ_floor)
        @test r0.n == n && r0.p_equivalence < 0.05 && r0.p_difference > 0.05
        @test r1.p_difference < 0.05 && r1.p_equivalence > 0.05
        @test r1.ci_low > 0 && r0.ci_low < 0 < r0.ci_high
        @test isnan(nd.paired_shift([1.0]; margin = 1.0).mean)
        @test nd.paired_shift([0.0, 0.0, 0.0]; margin = 0.1).p_equivalence == 0.0
        # Family-wise error under Holm: 200 null arms with a zero margin
        p = [nd.paired_shift(0.01 .* randn(rng, n); margin = 0.0).p_difference for _ in 1:200]
        @test count(<(0.05), p) > 0                      # uncorrected tests reject some
        @test count(<(0.05), nd.holm_adjust(p)) ≤ 1      # the corrected family hardly ever
    end

    @testset "window mean of a difference" begin
        t = collect(0.0:1.0:10.0)
        a = (t, 2.0 .* t)
        b = (t .+ 0.5, 2.0 .* (t .+ 0.5) .- 1.0)       # the same line shifted down by 1
        @test isapprox(nd._window_mean_difference(a, b, 0.0), 1.0; atol = 1e-12)
        @test isapprox(nd._window_mean_difference(a, b, 4.0), 1.0; atol = 1e-12)
        @test isnan(nd._window_mean_difference(a, b, 20.0))
        @test isapprox(nd._window_mean(a, 5.0), 15.0; atol = 1e-12)
    end
end
