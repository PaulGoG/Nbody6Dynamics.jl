# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Binary census: the union per pair of the overlapping or adjacent level-E
# (engine-regularised) and level-P (bound) spans.

@testset "Binary census (levels E and P)" begin
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
    E(n1, n2, t0, t1; kw...) = filled(
        BinarySpan;
        name1 = n1,
        name2 = n2,
        t_form_nb = t0,
        t_end_nb = t1,
        origin = :dynamical,
        end_class = :unresolved,
        ks_time_fraction = 1.0,
        kw...,
    )
    P(n1, n2, lo, hi, elo, ehi; kw...) = filled(
        BoundSpan;
        name1 = n1,
        name2 = n2,
        t_form_lo = lo,
        t_form_hi = hi,
        t_end_lo = elo,
        t_end_hi = ehi,
        kw...,
    )

    # A soft pair: KS episodes 2–3 and 6–7, bound at the snapshots 2…8 (unit cadence): one census span
    c = binary_census(
        [E(5, 6, 2.2, 3.1), E(5, 6, 6.3, 7.4)],
        [P(5, 6, 1.0, 2.0, 8.0, 9.0; end_unbound = true)];
        nbin0 = 0,
    )
    s = only(c)
    @test s.levels === :EP && s.n_level_e == 2 && s.n_level_p == 1
    @test s.t_form_nb == 2.0 && s.t_end_nb == 8.0 && s.t_form_e == 2.2 && s.t_end_e == 7.4
    @test s.origin === :unresolved && s.end_class === :separated   # level P first and last
    @test s.ks_time_fraction == 1.0
    # A hard pair seen only by level E, and a wide pair only by level P
    c2 = binary_census(
        [E(1, 2, 0.0, 50.0; origin = :primordial, end_class = :gw_merger, primordial = true)],
        [P(3, 4, NaN, 0.0, 99.0, NaN; right_censored = true, left_censored = true)];
        nbin0 = 2,
    )
    @test length(c2) == 2
    a = only(filter(s -> s.name1 == 1, c2))
    b = only(filter(s -> s.name1 == 3, c2))
    @test a.levels === :E &&
          a.origin === :primordial &&
          a.end_class === :gw_merger &&
          isnan(a.t_form_p)
    @test b.levels === :P &&
          b.origin === :primordial &&
          b.end_class === :censored_run_end &&
          b.right_censored
    @test isnan(b.ks_time_fraction)
    # Separated spans of one pair stay separate; flybys are left out
    c3 = binary_census(
        [E(7, 8, 1.0, 2.0), E(7, 8, 20.0, 21.0; end_class = :flyby)],
        [P(7, 8, 9.0, 10.0, 12.0, 13.0)];
        nbin0 = 0,
    )
    @test length(c3) == 2 && c3[1].levels === :E && c3[2].levels === :P
    @test issorted(c3; by = s -> (s.t_form_nb, s.name1, s.name2))
    @test_throws ArgumentError binary_census(BinarySpan[], BoundSpan[]; nbin0 = -1)
    path = write_binary_census(joinpath(mktempdir(), "binary_census.csv"), c)
    @test countlines(path) == 2
end
