# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Survival statistics of binary lifetimes: Kaplan–Meier with the Greenwood
# variance and delayed entry, Nelson–Aalen, Aalen–Johansen cumulative incidence
# of competing causes, the Turnbull NPMLE, the log-rank test and Cox–Snell
# residuals, against hand values, closed forms and simulated exponential
# samples (StableRNGs, so the reference samples hold across Julia versions),
# and the event tables built from binary spans.

using StableRNGs: StableRNG

@testset "Survival statistics" begin
    nd = Nbody6Dynamics
    # Exponential deviate of the given rate
    draw(rng, rate) = -log1p(-rand(rng)) / rate
    # Index of the last jump time ≤ t
    at(times, t) = searchsortedlast(times, t)
    # An event table without delayed entry unless `entry` is given
    mktable(x, c; entry = zeros(length(x)), labels = [:event]) = nd.EventTable(;
        time = Float64.(x),
        entry = Float64.(entry),
        cause = Int.(c),
        id = collect(1:length(x)),
        cause_labels = labels,
    )
    # Exponential exit times of rate 1 with independent exponential censoring of the given rate
    function censored_sample(rng, N, rate_c)
        x = zeros(N)
        c = zeros(Int, N)
        for i in 1:N
            T = draw(rng, 1.0)
            C = draw(rng, rate_c)
            x[i] = min(T, C)
            c[i] = T <= C ? 1 : 0
        end
        return x, c
    end

    @testset "Hand example" begin
        # No censoring: times 1, 2, 2, 3
        km = nd.kaplan_meier(mktable([1, 2, 2, 3], [1, 1, 1, 1]))
        @test km.time == [1.0, 2.0, 3.0]
        @test km.survival ≈ [0.75, 0.25, 0.0]
        @test km.n_risk == [4, 3, 1] && km.n_events == [1, 2, 1]
        # Greenwood: 0.75² · 1/(4·3); 0.25² · [1/(4·3) + 2/(3·1)]; n = d at t = 3
        @test km.variance[1] ≈ 0.75^2 / 12
        @test km.variance[2] ≈ 0.25^2 * (1 / 12 + 2 / 3)
        @test km.variance[3] == Inf
        # The last subject censored at 3: two event times, Ŝ = 0.75, 0.25, 0.25 at 1, 2, 3
        kc = nd.kaplan_meier(mktable([1, 2, 2, 3], [1, 1, 1, 0]))
        @test length(kc.time) == 2 && kc.censor_times == [3.0]
        @test kc.n_censored == [0, 0]    # the censoring follows the last event time
        @test [nd.survival_at(kc, t) for t in (1.0, 2.0, 3.0)] ≈ [0.75, 0.25, 0.25]
        @test nd.survival_at(kc, 0.5) == 1.0
        @test nd.survival_at(kc, 1.5) == 0.75
        @test nd.survival_at(kc, 10.0) == 0.25
        @test kc.variance ≈ km.variance[1:2]
        # Nelson–Aalen: Ĥ = 1/4, 1/4 + 2/3; Var = 1·3/4³, + 2·1/3³
        na = nd.nelson_aalen(mktable([1, 2, 2, 3], [1, 1, 1, 0]))
        @test na.time == kc.time
        @test na.hazard ≈ [0.25, 0.25 + 2 / 3]
        @test na.variance ≈ [3 / 64, 3 / 64 + 2 / 27]
        # Intervals bracket the estimate; the normal quantile
        @test nd._survival_normal_quantile(0.975) ≈ 1.959963984540054 rtol = 1e-10
        @test nd._survival_normal_quantile(0.5) ≈ 0.0 atol = 1e-12
        lo, hi = nd.confidence_interval(kc)
        @test all(lo .< kc.survival .< hi) && all(0 .<= lo) && all(hi .<= 1)
        lo_p, hi_p = nd.confidence_interval(kc; transform = :plain)
        @test lo_p ≈ max.(kc.survival .- 1.959963984540054 .* sqrt.(kc.variance), 0.0)
        @test all(hi_p .<= 1)
        l0, h0 = nd.confidence_interval(km)
        @test isnan(l0[3]) && isnan(h0[3])    # Ŝ = 0
        @test_throws ArgumentError nd.confidence_interval(kc; level = 1.0)
        @test_throws ArgumentError nd.confidence_interval(kc; transform = :log)
        # Rows with entry ≥ time or a cause without label are refused
        @test_throws ArgumentError nd.kaplan_meier(mktable([1.0], [1]; entry = [1.0]))
        @test_throws ArgumentError nd.kaplan_meier(mktable([1.0], [2]))
        # No events: an empty curve, Ŝ = 1 everywhere
        k0 = nd.kaplan_meier(mktable([1.0, 2.0], [0, 0]))
        @test isempty(k0.time) && nd.survival_at(k0, 5.0) == 1.0 && k0.censor_times == [1.0, 2.0]
    end

    @testset "Exponential sample with censoring" begin
        x, c = censored_sample(StableRNG(1), 5000, 0.5)
        tab = mktable(x, c)
        km = nd.kaplan_meier(tab)
        lo, hi = nd.confidence_interval(km)
        na = nd.nelson_aalen(tab)
        @test sum(km.n_events) == count(==(1), c)
        @test issorted(km.time) && all(diff(km.survival) .< 0)
        for t in (0.5, 1.0, 2.0)
            j = at(km.time, t)
            S = exp(-t)
            @test abs(km.survival[j] - S) <= 3 * sqrt(km.variance[j])
            @test lo[j] <= S <= hi[j]
            @test abs(na.hazard[j] - t) <= 3 * sqrt(na.variance[j])
        end
    end

    @testset "Delayed entry" begin
        # Exponential exit times of rate 1 observed only when they exceed the entry L ~ U(0, 1)
        rng = StableRNG(2)
        x = Float64[]
        ent = Float64[]
        while length(x) < 5000
            L = rand(rng)
            T = draw(rng, 1.0)
            if T > L
                push!(x, T)
                push!(ent, L)
            end
        end
        km = nd.kaplan_meier(mktable(x, ones(Int, length(x)); entry = ent))
        # Risk sets respect the entry times
        for j in 1:20
            @test km.n_risk[j] == count(i -> ent[i] < km.time[j] <= x[i], eachindex(x))
        end
        # Ŝ(2)/Ŝ(1) estimates the survival from 1 to 2 given survival to 1, e^{−1};
        # its Greenwood variance is the sum over (1, 2]
        j1 = at(km.time, 1.0)
        j2 = at(km.time, 2.0)
        R = km.survival[j2] / km.survival[j1]
        var_R = R^2 * (km.variance[j2] / km.survival[j2]^2 - km.variance[j1] / km.survival[j1]^2)
        @test abs(R - exp(-1)) <= 3 * sqrt(var_R)
    end

    @testset "Aalen–Johansen, competing exponential causes" begin
        # Cause 1 at rate 1, cause 2 at rate 2, censoring at rate 0.5: F_k = (λ_k/3)(1 − e^{−3t})
        rng = StableRNG(3)
        N = 5000
        x = zeros(N)
        c = zeros(Int, N)
        for i in 1:N
            T1 = draw(rng, 1.0)
            T2 = draw(rng, 2.0)
            C = draw(rng, 0.5)
            x[i] = min(T1, T2, C)
            c[i] = x[i] == T1 ? 1 : x[i] == T2 ? 2 : 0
        end
        tab = mktable(x, c; labels = [:first, :second])
        ci = nd.aalen_johansen(tab)
        km = nd.kaplan_meier(tab)
        @test ci.cause_labels == [:first, :second]
        @test ci.time == km.time && ci.survival == km.survival && ci.n_risk == km.n_risk
        @test vec(sum(ci.n_events; dims = 2)) == km.n_events
        @test maximum(abs.(vec(sum(ci.incidence; dims = 2)) .- (1 .- ci.survival))) <= 1e-12
        for t in (0.3, 0.6, 1.0)
            j = at(ci.time, t)
            for (k, λ) in ((1, 1.0), (2, 2.0))
                F = λ / 3 * (1 - exp(-3t))
                @test abs(ci.incidence[j, k] - F) <= 3 * sqrt(ci.variance[j, k])
            end
        end
        # 1 − KM with cause 2 censored overestimates F_1
        naive = nd.kaplan_meier(mktable(x, ifelse.(c .== 1, 1, 0)))
        @test 1 - naive.survival[at(naive.time, 1.0)] > ci.incidence[at(ci.time, 1.0), 1]
    end

    @testset "Aalen–Johansen variance calibration" begin
        rng = StableRNG(4)
        F = 1 / 3 * (1 - exp(-1.8))
        n_rep = 300
        N = 400
        hits = 0
        x = zeros(N)
        c = zeros(Int, N)
        for _ in 1:n_rep
            for i in 1:N
                T1 = draw(rng, 1.0)
                T2 = draw(rng, 2.0)
                x[i] = min(T1, T2)
                c[i] = T1 < T2 ? 1 : 2
            end
            ci = nd.aalen_johansen(mktable(x, c; labels = [:first, :second]))
            j = at(ci.time, 0.6)
            hits += abs(ci.incidence[j, 1] - F) <= 1.96 * sqrt(ci.variance[j, 1])
        end
        @test 0.92 <= hits / n_rep <= 0.98
    end

    @testset "Turnbull NPMLE" begin
        # Exact and right-censored data: the Kaplan–Meier estimate
        x, c = censored_sample(StableRNG(5), 60, 0.5)
        upper = [c[i] == 1 ? x[i] : Inf for i in eachindex(x)]
        tb = nd.turnbull(x, upper; tol = 1e-14, max_iter = 100_000)
        km = nd.kaplan_meier(mktable(x, c))
        @test tb.time == km.time
        @test maximum(abs.(tb.survival .- km.survival)) <= 1e-10
        @test all(isnan, tb.variance) && all(iszero, tb.n_risk) && all(iszero, tb.n_events)
        @test tb.censor_times == sort(x[c .== 0])
        # Interval censoring: inspections every 0.2 from a random offset, the last at or before 2
        rng = StableRNG(6)
        N = 4000
        lower = zeros(N)
        upper = zeros(N)
        for i in 1:N
            T = draw(rng, 1.0)
            o = 0.2 * rand(rng)
            last_inspection = o + 0.2 * floor((2 - o) / 0.2)
            if T > last_inspection
                lower[i] = last_inspection
                upper[i] = Inf
            elseif T <= o
                lower[i] = 0.0
                upper[i] = o
            else
                lower[i] = o + 0.2 * floor((T - o) / 0.2)
                upper[i] = lower[i] + 0.2
            end
        end
        ti = nd.turnbull(lower, upper)
        for t in (0.5, 1.0, 1.5)
            @test abs(nd.survival_at(ti, t) - exp(-t)) <= 0.03
        end
        @test issorted(ti.time) && all(diff(ti.survival) .<= 0)
        @test_throws ArgumentError nd.turnbull([2.0], [1.0])
        @test_throws DimensionMismatch nd.turnbull([1.0, 2.0], [3.0])
    end

    @testset "Log-rank test" begin
        rng = StableRNG(7)
        x = [draw(rng, 1.0) for _ in 1:1000]
        g = [isodd(i) ? 1 : 2 for i in 1:1000]
        same = nd.logrank_test(mktable(x, ones(Int, 1000)), g)
        @test same.p_value > 0.01
        @test sum(same.observed) == 1000 && sum(same.expected) ≈ 1000
        y = [i <= 500 ? draw(rng, 1.0) : draw(rng, 2.0) for i in 1:1000]
        h = [i <= 500 ? 1 : 2 for i in 1:1000]
        differ = nd.logrank_test(mktable(y, ones(Int, 1000)), h)
        @test differ.p_value < 1e-6
        @test differ.observed[2] > differ.expected[2]    # the faster group fails more often
        @test_throws ArgumentError nd.logrank_test(mktable(x, ones(Int, 1000)), mod1.(1:1000, 3))
        @test_throws ArgumentError nd.logrank_test(mktable(x, ones(Int, 1000)), ones(Int, 1000))
    end

    @testset "Cox–Snell residuals" begin
        x, c = censored_sample(StableRNG(8), 5000, 0.5)
        cs = nd.cox_snell_residuals(mktable(x, c), t -> t)
        @test cs.residual ≈ x && cs.event == (c .== 1)
        na = nd.nelson_aalen(mktable(cs.residual, Int.(cs.event)))
        for q in (0.5, 1.0, 2.0)
            j = at(na.time, q)
            @test abs(na.hazard[j] - na.time[j]) <= 3 * sqrt(na.variance[j])
        end
        # Delayed entry: Λ(time) − Λ(entry)
        r = nd.cox_snell_residuals(mktable([1.5], [1]; entry = [0.5]), t -> t^2)
        @test r.residual ≈ [2.0] && r.event == [true]
        @test_throws ArgumentError nd.cox_snell_residuals(mktable([1.0], [1]), t -> -t)
    end

    @testset "Event tables from binary spans" begin
        # A record with the given fields and the others absent (NaN, 0, false, :none)
        function filled(::Type{T}; kw...) where {T}
            values = map(fieldnames(T), fieldtypes(T)) do f, F
                haskey(kw, f) ? kw[f] :
                F <: AbstractFloat ? F(NaN) :
                F <: Bool ? false :
                F <: Integer ? zero(F) :
                F === Symbol ? :none : error("no absent value for field $f::$F")
            end
            return T(values...)
        end
        # T* = 2 Myr
        span(t0, t1, end_class; kw...) = filled(
            nd.BinarySpan;
            name1 = 1,
            name2 = 2,
            t_form_nb = t0,
            t_end_nb = t1,
            t_form_myr = 2t0,
            t_end_myr = 2t1,
            origin = :dynamical,
            end_class = end_class,
            kw...,
        )
        spans = [
            span(0.0, 10.0, :gw_merger; left_censored = true, origin = :primordial),
            span(2.0, 5.0, :exchange),
            span(4.0, 11.5, :censored_run_end; right_censored = true),
            span(6.0, 7.0, :unresolved),
            span(1.0, 9.0, :gw_merger),
        ]
        # Age axis: the left-censored span is left out; the open span is censored at the record end
        age = @test_logs (:info, r"left-censored") nd.event_table(spans; t_end = 12.0, unit = :nb)
        @test age.id == [2, 3, 4, 5]
        @test age.time ≈ [3.0, 8.0, 1.0, 8.0]
        @test all(iszero, age.entry)
        @test age.cause_labels == [:exchange, :gw_merger]
        @test age.cause == [1, 0, 0, 2]
        age_myr = @test_logs (:info, r"left-censored") nd.event_table(spans; t_end = 24.0)
        @test age_myr.time ≈ 2 .* age.time && age_myr.cause == age.cause
        # Cluster axis: delayed entry at formation, the left-censored span enters at 0
        cl = nd.event_table(spans; time_axis = :cluster, t_end = 12.0, unit = :nb)
        @test cl.id == collect(1:5)
        @test cl.entry == [0.0, 2.0, 4.0, 6.0, 1.0]
        @test cl.time == [10.0, 5.0, 12.0, 7.0, 9.0]
        @test cl.cause_labels == [:gw_merger, :exchange]
        @test cl.cause == [1, 2, 0, 0, 1]
        # At risk at the events 5, 9, 10: spans 1, 2, 3, 5; 1, 3, 5; 1, 3
        @test nd.kaplan_meier(cl).n_risk == [4, 3, 2]
        # A user mapping: exchanges count as censorings
        m = Dict(
            :gw_merger => :gw_merger,
            :exchange => :censored,
            :unresolved => :censored,
            :censored_run_end => :censored,
        )
        cm = nd.event_table(spans; time_axis = :cluster, t_end = 12.0, unit = :nb, causes = m)
        @test cm.cause == [1, 0, 0, 0, 1] && cm.cause_labels == [:gw_merger]
        # Refused: a span ending after the record, a zero lifetime, unknown axis or class
        @test_throws ArgumentError nd.event_table(
            spans;
            time_axis = :cluster,
            t_end = 10.5,
            unit = :nb,
        )
        @test_throws ArgumentError nd.event_table(
            [span(3.0, 3.0, :exchange)];
            t_end = 12.0,
            unit = :nb,
        )
        @test_throws ArgumentError nd.event_table(spans; time_axis = :binary, t_end = 12.0)
        @test_throws ArgumentError nd.event_table(spans; t_end = 24.0, unit = :yr)
        @test_throws ArgumentError nd.event_table(
            [span(1.0, 2.0, :roche_lobe)];
            t_end = 12.0,
            unit = :nb,
        )
    end

    @testset "Writers" begin
        dir = mktempdir()
        km = nd.kaplan_meier(mktable([1, 2, 2, 3], [1, 1, 1, 0]))
        path = nd.write_survival_curve(joinpath(dir, "survival.csv"), km)
        lines = readlines(path)
        @test lines[1] == "time,survival,variance,lower,upper,n_risk,n_events,n_censored"
        @test length(lines) == 1 + length(km.time)
        @test length(split(lines[2], ",")) == 8
        ci = nd.aalen_johansen(mktable([1, 2, 2, 3], [1, 2, 1, 0]; labels = [:a, :b]))
        path2 = nd.write_cumulative_incidence(joinpath(dir, "incidence.csv"), ci)
        lines2 = readlines(path2)
        @test lines2[1] == "time,cause,incidence,variance,survival,n_risk,n_events"
        @test length(lines2) == 1 + 2 * length(ci.time)
        @test split(lines2[2], ",")[2] == "a" && split(lines2[3], ",")[2] == "b"
        # An existing file is backed up
        nd.write_survival_curve(path, km)
        @test isfile(joinpath(dir, "survival#1.csv"))
    end
end
