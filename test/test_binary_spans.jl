# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# WP10 10.2 increment 1: engine-regularised binary spans (level E) from the
# 10.1 event records. The implementation is to satisfy this file; the method
# is specified in the workspace notes (binary spans, increment 1). API:
#
#   binary_spans(log::BinaryEventLog;
#                listings::AbstractVector{KSPairListing} = KSPairListing[],
#                collisions::AbstractVector{CollisionRecord} = CollisionRecord[],
#                coalescences::AbstractVector{CoalescenceRecord} = CoalescenceRecord[],
#                nzero::Integer, nbin0::Integer, tstar_myr::Real,
#                t_start_nb::Real = 0.0, t_end_nb::Real = NaN, n_p::Real = 2.0)
#       -> Vector{BinarySpan}, sorted by (t_form_nb, name1, name2)
#
# BinarySpan fields used here: name1 < name2 (component names), primordial,
# origin, origin_sublabel, t_form_nb, t_end_nb, t_form_myr, t_end_myr,
# left_censored, right_censored, end_class, n_ks_episodes, n_gaps_merged,
# in_chain_count, hierarchy_flag, ecc_form, semi_form_nb, kstar1_form,
# kstar2_form, ecc_end.

@testset "Binary spans (level E)" begin
    nd = Nbody6Dynamics
    # T* = 1 Myr: a period of P N-body units is P × 3.6525e8 days
    days(p_nb) = p_nb * 3.6525e8
    line = Ref(0)
    function ks(kind, t, n1, n2; p_nb = 0.3, semi = 1.0e-3, ecc = 0.5, eb = -1.0e-4, k1 = 0, k2 = 0)
        line[] += 1
        return KSEvent(;
            line = line[],
            interval = 1,
            time_nb = t,
            kind = kind,
            name1 = n1,
            name2 = n2,
            name_cm = 10_000 + min(n1, n2),
            kstar1 = k1,
            kstar2 = k2,
            kstar_cm = 0,
            semi_nb = semi,
            ecc = ecc,
            eb_nb = eb,
            period_d = days(p_nb),
            m1_msun = 1.0,
            m2_msun = 0.5,
        )
    end
    function eventlog(;
        ks = KSEvent[],
        fewbody = FewBodyEvent[],
        exchange = ExchangeEvent[],
        gw = GWEvent[],
        escape = BinaryEscapeEvent[],
        mass_loss = MassLossEvent[],
    )
        log = BinaryEventLog()
        append!(log.ks, ks)
        append!(log.fewbody, fewbody)
        append!(log.exchange, exchange)
        append!(log.gw, gw)
        append!(log.escape, escape)
        append!(log.mass_loss, mass_loss)
        return log
    end
    spans(log; kw...) =
        binary_spans(log; nzero = 10_000, nbin0 = 2, tstar_myr = 1.0, t_end_nb = 10.0, kw...)
    # A file record with the given fields and the others absent (NaN, 0, false, :none)
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

    # 1. A soft primordial pair that leaves KS at every apocentre and re-enters within one
    #    period: one span, all episodes merged, censored at the run end
    soft = [
        ks(:ks_new_primordial, 0.0, 1, 2),
        ks(:ks_end, 0.5, 1, 2),
        ks(:ks_new, 0.6, 1, 2),
        ks(:ks_end, 1.1, 1, 2),
        ks(:ks_new, 1.2, 2, 1),   # component order swapped: the same pair
        ks(:ks_end, 1.7, 1, 2),
        ks(:ks_new, 1.8, 1, 2),
    ]
    s1 = spans(eventlog(; ks = soft))
    @test length(s1) == 1
    @test (s1[1].name1, s1[1].name2) == (1, 2)
    @test s1[1].primordial && s1[1].origin === :primordial
    @test s1[1].n_ks_episodes == 4 && s1[1].n_gaps_merged == 3
    @test s1[1].t_form_nb == 0.0 && !s1[1].left_censored
    @test s1[1].right_censored && s1[1].end_class === :censored_run_end && s1[1].t_end_nb == 10.0
    @test s1[1].t_form_myr == 0.0 && s1[1].t_end_myr == 10.0
    @test s1[1].ecc_form == 0.5 && s1[1].semi_form_nb == 1.0e-3

    # 2. A gap of three periods ends the span; the pair re-forms as primordial_reformed
    s2 = spans(
        eventlog(;
            ks = [
                ks(:ks_new_primordial, 0.0, 1, 2),
                ks(:ks_end, 0.5, 1, 2),
                ks(:ks_new, 1.4, 1, 2),
                ks(:ks_end, 1.9, 1, 2),
            ],
        ),
    )
    @test length(s2) == 2
    @test s2[1].end_class === :unresolved && s2[1].t_end_nb == 0.5 && !s2[1].right_censored
    @test s2[2].origin === :primordial_reformed && s2[2].t_form_nb == 1.4
    @test s2[2].end_class === :unresolved
    # with n_p = 4 the 0.9 gap (3 periods) is bridged
    @test length(
        spans(
            eventlog(;
                ks = [
                    ks(:ks_new_primordial, 0.0, 1, 2),
                    ks(:ks_end, 0.5, 1, 2),
                    ks(:ks_new, 1.4, 1, 2),
                ],
            );
            n_p = 4.0,
        ),
    ) == 1

    # 3. A chain holding both names bridges a gap longer than n_p periods
    chain_start = FewBodyEvent(;
        line = 50,
        interval = 1,
        time_nb = 1.0,
        kind = :chain_new,
        n_members = 3,
        names = (5, 6, 7, 0, 0, 0),
    )
    chain_end = FewBodyEvent(;
        line = 60,
        interval = 1,
        time_nb = 2.5,
        kind = :chain_end,
        n_members = 3,
        names = (5, 6, 7, 0, 0, 0),
    )
    s3 = spans(
        eventlog(;
            ks = [ks(:ks_new, 0.2, 5, 6), ks(:ks_end, 1.0, 5, 6), ks(:ks_new, 2.5, 5, 6)],
            fewbody = [chain_start, chain_end],
        ),
    )
    pair56 = filter(s -> (s.name1, s.name2) == (5, 6), s3)
    @test length(pair56) == 1 && pair56[1].n_ks_episodes == 2 && pair56[1].in_chain_count == 1
    @test pair56[1].origin === :dynamical

    # 4. An exchange: the old pair ends by exchange, the new one is of exchange origin
    xch = ExchangeEvent(;
        line = 70,
        interval = 1,
        time_nb = 3.0,
        kind = :exchange,
        name1_initial = 11,
        name2_initial = 12,
        name1_final = 11,
        name2_final = 13,
    )
    s4 = spans(
        eventlog(;
            ks = [ks(:ks_new, 1.0, 11, 12), ks(:ks_end, 3.0, 11, 12), ks(:ks_new, 3.0, 11, 13)],
            exchange = [xch],
        ),
    )
    old = only(filter(s -> (s.name1, s.name2) == (11, 12), s4))
    new = only(filter(s -> (s.name1, s.name2) == (11, 13), s4))
    @test old.end_class === :exchange && old.t_end_nb == 3.0
    @test new.origin === :exchange && new.right_censored

    # 5. A GW merger: the trigger line without names joined by time to the coll.13 record
    trigger = GWEvent(; line = 80, interval = 1, time_nb = 5.0, kind = :gw_ksint_coll)
    coll = filled(
        CollisionRecord;
        time_nb = 5.0,
        name1 = 21,
        name2 = 22,
        kstar1 = 14,
        kstar2 = 14,
        kstar_product = 14,
    )
    s5 = spans(
        eventlog(;
            ks = [
                ks(:ks_new, 2.0, 21, 22; k1 = 14, k2 = 14),
                ks(:ks_end, 5.0, 21, 22; k1 = 14, k2 = 14),
            ],
            gw = [trigger],
        );
        collisions = [coll],
    )
    @test only(s5).end_class === :gw_merger && only(s5).t_end_nb == 5.0
    @test only(s5).kstar1_form == 14 && only(s5).kstar2_form == 14
    # the same collision of two main-sequence stars without a GW trigger is a collision
    coll_ms = filled(
        CollisionRecord;
        time_nb = 5.0,
        name1 = 21,
        name2 = 22,
        kstar1 = 1,
        kstar2 = 1,
        kstar_product = 1,
    )
    s5b = spans(
        eventlog(; ks = [ks(:ks_new, 2.0, 21, 22), ks(:ks_end, 5.0, 21, 22)]);
        collisions = [coll_ms],
    )
    @test only(s5b).end_class === :collision

    # 6. A supernova disruption: the mass-loss termination is not followed by a re-regularisation
    mdend = MassLossEvent(;
        line = 90,
        interval = 1,
        time_nb = 4.0,
        kind = :ks_end_mass_loss,
        name1 = 31,
        name2 = 32,
        period_d = days(0.3),
    )
    s6 = spans(eventlog(; ks = [ks(:ks_new, 1.0, 31, 32)], mass_loss = [mdend]))
    @test only(s6).end_class === :sn_disruption && only(s6).t_end_nb == 4.0
    mdnew = MassLossEvent(;
        line = 91,
        interval = 1,
        time_nb = 4.0,
        kind = :ks_new_mass_loss,
        name1 = 31,
        name2 = 32,
    )
    s6b = spans(
        eventlog(;
            ks = [ks(:ks_new, 1.0, 31, 32), ks(:ks_new, 4.0, 31, 32)],
            mass_loss = [mdend, mdnew],
        ),
    )
    @test only(s6b).end_class === :censored_run_end && only(s6b).n_ks_episodes == 2

    # 7. An escaping bound binary
    esc = BinaryEscapeEvent(;
        line = 95,
        interval = 1,
        time_nb = 6.0,
        kind = :binary_escape,
        name1 = 41,
        name2 = 42,
    )
    s7 =
        spans(eventlog(; ks = [ks(:ks_new, 1.0, 41, 42), ks(:ks_end, 6.0, 41, 42)], escape = [esc]))
    @test only(s7).end_class === :escape_bound && !only(s7).right_censored

    # 8. The NS/BH BINARY line duplicates the NEW KSREG of the same pair and time: one episode
    s8 = spans(
        eventlog(;
            ks = [
                ks(:ks_new, 1.0, 51, 52; k1 = 14, k2 = 13),
                ks(:ks_new_compact, 1.0, 51, 52; k1 = 14, k2 = 13),
            ],
        ),
    )
    @test only(s8).n_ks_episodes == 1

    # 9. A capture that starts mid-episode: left-censored at the first listing holding the pair
    listing =
        KSPairListing(; time_nb = 2.0, records = [filled(KSPairRecord; name1 = 61, name2 = 62)])
    s9 = spans(eventlog(; ks = [ks(:ks_end, 3.0, 61, 62)]); listings = [listing])
    @test only(s9).left_censored && only(s9).t_form_nb == 2.0
    s9b = spans(eventlog(; ks = [ks(:ks_end, 3.0, 61, 62)]); t_start_nb = 1.5)
    @test only(s9b).left_censored && only(s9b).t_form_nb == 1.5

    # 10. Primordial rule: (2k-1, 2k) with 2k ≤ 2 NBIN0 only; (2, 3) is not primordial
    s10 = spans(eventlog(; ks = [ks(:ks_new, 0.5, 2, 3), ks(:ks_new_primordial, 0.0, 3, 4)]))
    @test !only(filter(s -> (s.name1, s.name2) == (2, 3), s10)).primordial
    @test only(filter(s -> (s.name1, s.name2) == (3, 4), s10)).primordial

    # 12. A hyperbolic encounter is a flyby; a bound pair that leaves KS unbound is a disruption
    hyp = (semi = -1.0e-3, ecc = 1.5, eb = 1.0e-4, p_nb = 0.3)
    s12 =
        spans(eventlog(; ks = [ks(:ks_new, 1.0, 71, 72; hyp...), ks(:ks_end, 1.2, 71, 72; hyp...)]))
    @test only(s12).end_class === :flyby && only(s12).ecc_form == 1.5
    s12b = spans(eventlog(; ks = [ks(:ks_new, 1.0, 71, 72), ks(:ks_end, 1.2, 71, 72; hyp...)]))
    @test only(s12b).end_class === :disruption

    # 13. A GW capture at the pericentre of a hyperbolic NS/BH encounter, then the GW merger
    capture = GWEvent(;
        line = 100,
        interval = 1,
        time_nb = 2.003,
        kind = :gw_capture,
        name1 = 82,
        name2 = 81,
        kstar1 = 14,
        kstar2 = 14,
    )
    trigger13 = GWEvent(; line = 101, interval = 1, time_nb = 2.004, kind = :gw_unpert_coll)
    coll13 =
        filled(CollisionRecord; time_nb = 2.004, name1 = 81, name2 = 82, kstar1 = 14, kstar2 = 14)
    s13 = spans(
        eventlog(;
            ks = [
                ks(:ks_new, 2.0, 81, 82; k1 = 14, k2 = 14, hyp...),
                ks(:ks_end, 2.004, 81, 82; k1 = 14, k2 = 14),
            ],
            gw = [capture, trigger13],
        );
        collisions = [coll13],
    )
    @test only(s13).origin === :dynamical && only(s13).origin_sublabel === :gw_capture
    @test only(s13).end_class === :gw_merger

    # 11. Bad arguments
    @test_throws ArgumentError binary_spans(BinaryEventLog(); nzero = 0, nbin0 = 0, tstar_myr = 1.0)
    @test_throws ArgumentError binary_spans(
        BinaryEventLog();
        nzero = 10,
        nbin0 = 0,
        tstar_myr = 1.0,
        n_p = 0.0,
    )
    @test isempty(binary_spans(BinaryEventLog(); nzero = 10, nbin0 = 0, tstar_myr = 1.0))
end
