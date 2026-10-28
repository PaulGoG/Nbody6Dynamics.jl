# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Gravitational-wave catalogues: pair classes, the Peters track, band
# residence against the circular closed forms, the state after GW decay, and
# the compact-binary, escaped-binary and merger rows of a synthetic event log.

using StableRNGs: StableRNG

@testset "GW catalogues" begin
    nd = Nbody6Dynamics

    @testset "Pair classes" begin
        @test nd.pair_class(14, 14) === :bhbh
        @test nd.pair_class(14, 13) === :bhns
        @test nd.pair_class(13, 14) === :bhns
        @test nd.pair_class(13, 13) === :nsns
        @test nd.pair_class(11, 14) === :bhwd
        @test nd.pair_class(13, 12) === :nswd
        @test nd.pair_class(10, 12) === :wdwd
        @test nd.pair_class(14, 1) === :bh_star
        @test nd.pair_class(4, 13) === :ns_star
        @test nd.pair_class(11, 0) === :wd_star
        @test nd.pair_class(0, 9) === :star_star
        @test nd.pair_class(KSTAR_ABSENT, 14) === :unknown
        @test nd.pair_class(14, 15) === :unknown
    end

    @testset "Peters track" begin
        g(e) = e^(12 / 19) * (1 + 121 * e^2 / 304)^(870 / 2299) / (1 - e^2)
        @test nd.peters_semi_major_axis(0.3, 0.9, 0.9) == 0.3
        @test nd.peters_semi_major_axis(0.2, 0.0, 0.0) == 0.2
        @test nd.peters_semi_major_axis(1.0, 0.9, 0.5) ≈ g(0.5) / g(0.9) rtol = 1e-12
        @test nd.peters_semi_major_axis(1.0, 0.9, 0.0) == 0.0
        @test_throws ArgumentError nd.peters_semi_major_axis(1.0, 0.5, 0.6)
        @test_throws ArgumentError nd.peters_semi_major_axis(1.0, 0.5, -0.1)
        @test_throws ArgumentError nd.peters_semi_major_axis(1.0, 1.0, 0.5)
        @test_throws ArgumentError nd.peters_semi_major_axis(0.0, 0.5, 0.1)
        @test_throws ArgumentError nd.peters_semi_major_axis(1.0, 0.0, 0.1)
    end

    @testset "Band residence, circular references" begin
        # 30 + 30 M☉, f₂₂: τ(f) = (5/256)(G𝓜/c³)^{−5/3}(πf)^{−8/3}
        b = nd.band_residence(30, 30, 0.1, 0.0, 1e-4, 1e-1)
        @test b.t_in_yr ≈ 5.905668e8 rtol = 1e-4
        @test b.t_out_yr ≈ 5.947063e8 rtol = 1e-4
        @test !b.above_lo && !b.above_hi
        # a₀ of f₂₂ = 10⁻⁴ Hz, rounded down: already in band
        @test gw_frequency(30, 30, 0.02888433, 0.0) ≈ 1e-4 rtol = 1e-6
        b = nd.band_residence(30, 30, 0.02888433, 0.0, 1e-4, 1e-1)
        @test b.t_in_yr == 0
        @test b.above_lo && !b.above_hi
        @test b.t_out_yr ≈ 4.1395e6 rtol = 1e-3
        # Coalescence times from f₂₂ = 10⁻³ and 10⁻² Hz
        @test gw_frequency(30, 30, 0.0062229403, 0.0) ≈ 1e-3 rtol = 1e-6
        @test gw_inspiral_time(30, 30, 0.0062229403, 0.0) ≈ 8918.36 rtol = 1e-4
        @test gw_frequency(30, 30, 0.0013406918, 0.0) ≈ 1e-2 rtol = 1e-6
        @test gw_inspiral_time(30, 30, 0.0013406918, 0.0) ≈ 19.214 rtol = 1e-3
        @test chirp_mass(30, 30) ≈ 26.1165 rtol = 1e-4
        @test nd._gwcat_fdot_circular(26.1165, 1e-3) ≈ 1.33242e-15 rtol = 1e-4
        # Above both edges
        b = nd.band_residence(30, 30, 1e-5, 0.0, 1e-4, 1e-1)
        @test b.above_lo && b.above_hi && b.t_in_yr == 0 && b.t_out_yr == 0
        @test_throws ArgumentError nd.band_residence(30, 30, 0.1, 0.0, 1e-1, 1e-4)
        @test_throws ArgumentError nd.band_residence(30, 30, 0.1, 0.0, 0.0, 1e-1)
        @test_throws ArgumentError nd.band_residence(
            30,
            30,
            0.1,
            0.0,
            1e-4,
            1e-1;
            definition = :none,
        )
    end

    @testset "Eccentric track and the state after GW decay" begin
        T = gw_inspiral_time(30, 30, 0.1, 0.9)
        b = nd.band_residence(30, 30, 0.1, 0.9, 1e-3, 1e-1)
        r = eccentricity_at_frequency(30, 30, 0.1, 0.9, 1e-3)
        a_f = nd.peters_semi_major_axis(0.1, 0.9, r.e)
        @test b.t_in_yr + gw_inspiral_time(30, 30, a_f, r.e) ≈ T rtol = 1e-8
        @test 0 < b.t_in_yr < b.t_out_yr < T
        @test !b.above_lo && !b.above_hi

        a, e = nd.eccentricity_after(30, 30, 0.1, 0.9, 0.3T)
        @test gw_inspiral_time(30, 30, a, e) ≈ 0.7T rtol = 1e-8
        @test a == nd.peters_semi_major_axis(0.1, 0.9, e)
        @test 0 < e < 0.9
        @test nd.eccentricity_after(30, 30, 0.1, 0.9, 0.0) == (0.1, 0.9)
        @test all(isnan, nd.eccentricity_after(30, 30, 0.1, 0.9, T))
        @test all(isnan, nd.eccentricity_after(30, 30, 0.1, 0.9, 1.01T))
        T₀ = gw_inspiral_time(30, 30, 0.1, 0.0)
        a, e = nd.eccentricity_after(30, 30, 0.1, 0.0, 0.3T₀)
        @test e == 0.0
        @test a ≈ 0.1 * 0.7^(1 / 4) rtol = 1e-12
        @test_throws ArgumentError nd.eccentricity_after(30, 30, 0.1, 0.9, -1.0)
    end

    @testset "Catalogues of a synthetic log" begin
        tstar, rbar, vstar = 0.2, 1.0, 5.0
        au_per_nb = 206264.806
        evlog = BinaryEventLog()
        # Applied Peters steps of the pair (7, 3), name 3 the heavier component
        decay(t, semi, ecc) = GWEvent(;
            line = 1,
            interval = 1,
            time_nb = t,
            kind = :gw_ksint,
            res = 1,
            name1 = 7,
            name2 = 3,
            kstar1 = 14,
            kstar2 = 14,
            m1_msun = 20.0,
            m2_msun = 30.0,
            ecc = ecc,
            semi_nb = semi,
            gamma = 0.01,
            t_orb_nb = 1e-6,
            t_gr_nb = 10.0,
        )
        push!(evlog.gw, decay(200.0, 1e-10, 0.01))
        push!(evlog.gw, decay(100.0, 1e-7, 0.5))
        # A decay record with an absent type is dropped
        push!(
            evlog.gw,
            GWEvent(;
                line = 2,
                interval = 1,
                time_nb = 150.0,
                kind = :gw_tides,
                name1 = 7,
                name2 = 3,
                kstar2 = 14,
                m1_msun = 20.0,
                m2_msun = 30.0,
                ecc = 0.4,
                semi_nb = 1e-8,
            ),
        )
        push!(
            evlog.ks,
            KSEvent(;
                line = 3,
                interval = 1,
                time_nb = 50.0,
                kind = :ks_new_compact,
                name1 = 7,
                name2 = 3,
                kstar1 = 14,
                kstar2 = 14,
                m1_msun = 20.0,
                m2_msun = 30.0,
                ecc = 0.3,
                semi_nb = 2e-7,
                gamma = 0.2,
            ),
        )
        push!(
            evlog.gw,
            GWEvent(;
                line = 4,
                interval = 1,
                time_nb = 250.0,
                kind = :gw_recoil,
                time_myr = 250.0 * tstar,
                name1 = 3,
                name2 = 7,
                kstar1 = 14,
                kstar2 = 14,
                m1_msun = 30.0,
                m2_msun = 20.0,
                spin1 = 0.7,
                spin2 = 0.7,
                v_kick_kms = 100.0,
                v_esc_kms = 50.0,
            ),
        )
        push!(
            evlog.escape,
            BinaryEscapeEvent(;
                line = 5,
                interval = 1,
                time_nb = 250.5,
                kind = :binary_escape,
                name1 = 7,
                name2 = 3,
                kstar1 = 14,
                kstar2 = 14,
                m1_msun = 20.0,
                m2_msun = 30.0,
                semi_nb = 1e-7,
                ecc = 0.5,
                eb_nb = -1.0,
                period_d = 1.0,
                v_nb = 2.0,
                r_nb = 3.0,
            ),
        )
        span = BinarySpan(;
            name1 = 3,
            name2 = 7,
            primordial = true,
            hierarchy_outer = false,
            origin = :primordial,
            origin_sublabel = :none,
            mixed = false,
            t_form_nb = 40.0,
            t_end_nb = 250.0,
            t_form_myr = 40.0 * tstar,
            t_end_myr = 250.0 * tstar,
            left_censored = false,
            right_censored = false,
            end_class = :gw_merger,
            n_ks_episodes = 1,
            n_gaps_merged = 0,
            max_gap_nb = 0.0,
            ks_time_fraction = 1.0,
            in_chain_count = 0,
            hierarchy_flag = false,
            outer_partner_name = 0,
            ecc_form = 0.3,
            semi_form_nb = 2e-7,
            period_form_d = NaN,
            eb_form_nb = NaN,
            m1_form_msun = 30.0,
            m2_form_msun = 20.0,
            kstar1_form = 14,
            kstar2_form = 14,
            kstar_cm_form = 0,
            r_cm_form_nb = NaN,
            gamma_form = 0.2,
            ecc_end = NaN,
            semi_end_nb = NaN,
            period_end_d = NaN,
            eb_end_nb = NaN,
            m1_end_msun = NaN,
            m2_end_msun = NaN,
            kstar1_end = KSTAR_ABSENT,
            kstar2_end = KSTAR_ABSENT,
            kstar_cm_end = KSTAR_ABSENT,
            r_cm_end_nb = NaN,
            gamma_end = NaN,
            line_form = 3,
            line_end = 0,
        )
        # One bev.82 record of the pair: a = 10^{log₁₀(a/R☉)}/215.032 au
        a_bev = 1e-7 * au_per_nb
        bev = BinaryEvolutionSnapshot(
            120.0 * tstar,
            1,
            [
                BinaryRecord(
                    120.0,
                    1,
                    2,
                    7,
                    3,
                    14,
                    14,
                    0,
                    0.1,
                    0.4,
                    0.0,
                    log10(a_bev * 215.032),
                    20.0,
                    30.0,
                    ntuple(_ -> 0.0, 18)...,
                ),
            ],
        )

        # Compact-binary catalogue
        rows = nd.compact_binary_catalogue(
            evlog,
            [span];
            tstar_myr = tstar,
            rbar_pc = rbar,
            binary_evo = [bev],
        )
        @test length(rows) == 4
        @test [r.time_nb for r in rows] == [50.0, 100.0, 120.0, 200.0]
        @test [r.record_source for r in rows] == [:ns_bh_binary, :gw_decay, :bev82, :gw_decay]
        for r in rows
            @test (r.name1, r.name2) == (3, 7)
            @test (r.m1_msun, r.m2_msun) == (30.0, 20.0)
            @test r.q ≤ 1
            @test r.pair_class === :bhbh
            @test r.channel === :primordial
            @test r.span_origin === :primordial
            @test r.time_myr == r.time_nb * tstar
            @test r.t_merge_abs_myr == r.time_myr + r.t_gw_peters_myr
            @test isfinite(r.t_gw_peters_myr)
        end
        @test rows[1].a_au ≈ 2e-7 * au_per_nb rtol = 1e-14
        @test rows[2].a_au ≈ 1e-7 * au_per_nb rtol = 1e-14
        @test rows[3].a_au ≈ a_bev rtol = 1e-12
        @test rows[2].mchirp_msun ≈ chirp_mass(30, 20) rtol = 1e-14
        @test rows[2].t_gw_peters_myr ≈ gw_inspiral_time(30, 20, 1e-7 * au_per_nb, 0.5) / 1e6 rtol =
            1e-14
        @test rows[2].tgr_engine_nb == 10.0
        @test isnan(rows[1].tgr_engine_nb) && isnan(rows[3].gamma_pert)
        @test rows[1].flag_perturbed_no_brake && !rows[2].flag_perturbed_no_brake
        # The last decay record lies inside 25 Schwarzschild radii: post-Newtonian regime, above the band
        @test rows[4].flag_pn_regime && rows[4].above_band_hi
        @test isnan(rows[4].e_at_fref_lo) && 0 < rows[4].e_at_fref_hi < 0.01

        # Escaped-binary catalogue
        escaped =
            nd.escaped_binary_catalogue(evlog; tstar_myr = tstar, rbar_pc = rbar, vstar_kms = vstar)
        @test length(escaped) == 1
        x = only(escaped)
        @test (x.name1, x.name2) == (3, 7)
        @test x.escape_mode === :gw_recoil
        @test x.flag_kept_in_sim
        @test x.r_esc_pc == 3.0 * rbar && x.v_esc_kms == 2.0 * vstar
        @test x.t_merge_abs_myr == x.t_esc_myr + x.t_gw_peters_myr
        @test isnan(x.e_now) && isnan(x.t_to_merge_from_now_myr)
        escapers = [EscaperRecord(250.5 * tstar, 20.0, 0.0, 10.0, 14, 7, 0.0, 0.0)]
        x = only(
            nd.escaped_binary_catalogue(
                evlog;
                tstar_myr = tstar,
                rbar_pc = rbar,
                vstar_kms = vstar,
                escapers,
            ),
        )
        @test !x.flag_kept_in_sim
        x = only(
            nd.escaped_binary_catalogue(
                evlog;
                tstar_myr = tstar,
                rbar_pc = rbar,
                vstar_kms = vstar,
                t_now_myr = 250.5 * tstar + 0.1,
            ),
        )
        @test x.e_now < x.e && x.a_now_au < x.a_au
        @test x.t_to_merge_from_now_myr ≈ x.t_gw_peters_myr - 0.1 rtol = 1e-6

        # GW-merger catalogue
        coal(t, n1, n2) = CoalescenceRecord(;
            channel = :binary,
            iqcoll = 0,
            time_nb = t,
            name1 = n1,
            name2 = n2,
            kstar1 = 14,
            kstar2 = 14,
            kstar_product = 14,
            mass1_bse_msun = 30.0,
            mass2_bse_msun = 20.0,
            separation_nb = 1e-10,
            eccentricity = 0.0,
            semi_major_axis_nb = 1e-10,
            binding_energy_nb = -1.0,
            tidal_potential_correction_nb = 0.0,
            period_d = 1e-6,
            mass1_msun = 30.0,
            mass2_msun = 20.0,
            mass_product_msun = 50.0,
            mass_lost_msun = 0.0,
            radius_survivor_rsun = 1e-5,
            collision_distance_rsun = 1e-4,
            velocity_infinity_kms = 0.0,
            distance_density_centre_nb = 0.1,
            speed_nb = 1.0,
        )
        mergers = nd.gw_merger_catalogue(
            evlog,
            [coal(400.0, 3, 11), coal(250.0, 3, 7)],
            [span];
            tstar_myr = tstar,
            rbar_pc = rbar,
            vstar_kms = vstar,
        )
        @test length(mergers) == 2
        first_merger, second_merger = mergers
        @test first_merger.t_merge_nb == 250.0 && first_merger.t_merge_myr == 250.0 * tstar
        @test first_merger.pair_class === :bhbh && first_merger.channel === :binary
        @test first_merger.a_trigger_au ≈ 1e-10 * au_per_nb rtol = 1e-14
        @test first_merger.e_trigger == 0.01
        @test first_merger.q_peri_trigger_rs < 25
        @test first_merger.trigger_kind === :q25rs
        @test first_merger.v_kick_kms == 100.0 && first_merger.spin1 == 0.7
        @test !first_merger.retained
        @test 0 < first_merger.e_at_fref_hi < 0.01
        @test first_merger.generation == 1
        @test first_merger.span_primordial
        @test second_merger.generation == 2
        # coal.f may give the remnant the second name (larger BODY0): a later merger naming it counts
        swapped = nd.gw_merger_catalogue(
            evlog,
            [coal(250.0, 3, 7), coal(600.0, 7, 12)],
            [span];
            tstar_myr = tstar,
            rbar_pc = rbar,
            vstar_kms = vstar,
        )
        @test [m.generation for m in swapped] == [1, 2]
        @test second_merger.trigger_kind === :unknown
        @test isnan(second_merger.a_trigger_au) && isnan(second_merger.v_kick_kms)
        @test !second_merger.retained && !second_merger.span_primordial

        @testset "Writers" begin
            dir = mktempdir()
            for (writer, name, data, T) in (
                (nd.write_compact_binary_catalogue, "compact.csv", rows, nd.CompactBinaryRow),
                (nd.write_escaped_binary_catalogue, "escaped.csv", escaped, nd.EscapedBinaryRow),
                (nd.write_gw_merger_catalogue, "mergers.csv", mergers, nd.GWMergerRow),
            )
                path = writer(joinpath(dir, name), data)
                @test path == joinpath(dir, name)
                lines = readlines(path)
                @test first(lines) == join(fieldnames(T), ",")
                @test length(lines) == 1 + length(data)
                @test all(l -> count(==(','), l) == fieldcount(T) - 1, lines)
            end
            lines = readlines(joinpath(dir, "compact.csv"))
            cells = split(lines[2], ",")
            header = split(lines[1], ",")
            @test cells[findfirst(==("record_source"), header)] == "ns_bh_binary"
            @test cells[findfirst(==("tgr_engine_nb"), header)] == "NaN"
            @test cells[findfirst(==("flag_perturbed_no_brake"), header)] == "true"
        end

        @testset "Argument checks" begin
            @test_throws ArgumentError nd.compact_binary_catalogue(
                evlog,
                [span];
                tstar_myr = 0.0,
                rbar_pc = rbar,
            )
            @test_throws ArgumentError nd.compact_binary_catalogue(
                evlog,
                [span];
                tstar_myr = tstar,
                rbar_pc = -1.0,
            )
            @test_throws ArgumentError nd.compact_binary_catalogue(
                evlog,
                [span];
                tstar_myr = tstar,
                rbar_pc = rbar,
                band_hz = (1e-1, 1e-4),
            )
            @test_throws ArgumentError nd.escaped_binary_catalogue(
                evlog;
                tstar_myr = tstar,
                rbar_pc = rbar,
                vstar_kms = 0.0,
            )
            @test_throws ArgumentError nd.escaped_binary_catalogue(
                evlog;
                tstar_myr = tstar,
                rbar_pc = rbar,
                vstar_kms = vstar,
                f_ref_hz = (10.0, 1e-2),
            )
            @test_throws ArgumentError nd.gw_merger_catalogue(
                evlog,
                CoalescenceRecord[],
                [span];
                tstar_myr = -1.0,
                rbar_pc = rbar,
                vstar_kms = vstar,
            )
        end
    end

    @testset "One row per application" begin
        # The decision line (GR KSINT2) and the KSTIDE line of one application share the time
        # and the state: only the KSTIDE row is kept unless decision_lines = true.
        log2 = nd.BinaryEventLog()
        common = (
            interval = 1,
            name1 = 5,
            name2 = 9,
            kstar1 = 14,
            kstar2 = 14,
            m1_msun = 30.0,
            m2_msun = 30.0,
            semi_nb = 1.0e-6,
            ecc = 0.5,
            gamma = 1.0e-9,
            t_gr_nb = 100.0,
        )
        push!(log2.gw, nd.GWEvent(; line = 1, time_nb = 2.0, kind = :gw_ksint, common...))
        push!(log2.gw, nd.GWEvent(; line = 2, time_nb = 2.0, kind = :gw_decay_step, common...))
        push!(log2.gw, nd.GWEvent(; line = 3, time_nb = 3.0, kind = :gw_ksint, common...))
        rows2 = compact_binary_catalogue(log2, BinarySpan[]; tstar_myr = 0.2, rbar_pc = 1.0)
        @test [r.record_kind for r in rows2] == [:gw_decay_step, :gw_ksint]
        rows3 = compact_binary_catalogue(
            log2,
            BinarySpan[];
            tstar_myr = 0.2,
            rbar_pc = 1.0,
            decision_lines = true,
        )
        @test length(rows3) == 3
    end

    @testset "Strain columns and LDC export" begin
        log3 = nd.BinaryEventLog()
        # Three compact pairs, two of them quasi-circular
        for (i, (n1, n2, ecc)) in enumerate(((1, 2, 0.001), (3, 4, 0.005), (5, 6, 0.5)))
            push!(
                log3.gw,
                nd.GWEvent(;
                    line = i,
                    interval = 1,
                    time_nb = 10.0 * i,
                    kind = :gw_decay_step,
                    name1 = n1,
                    name2 = n2,
                    kstar1 = 14,
                    kstar2 = 13,
                    m1_msun = 20.0,
                    m2_msun = 1.4,
                    semi_nb = 1.0e-8,
                    ecc = ecc,
                    gamma = 1.0e-9,
                    t_gr_nb = 100.0,
                ),
            )
        end
        rows = compact_binary_catalogue(log3, BinarySpan[]; tstar_myr = 0.2, rbar_pc = 1.0)
        @test length(rows) == 3
        for r in rows
            @test isfinite(r.amp_ldc22_kpc) && isfinite(r.h2_avg_kpc) && isfinite(r.hc2_kpc)
            @test r.amp_ldc22_kpc ≈ gw_amplitude_ldc(r.mchirp_msun, r.f_gw22_hz) rtol = 1e-12
            @test r.h2_avg_kpc ≈ gw_harmonic_strain(r.m1_msun, r.m2_msun, r.a_au, r.e, 2) rtol =
                1e-12
            snr_ok = r.e ≤ 0.99 && 1e-5 ≤ r.f_gw22_hz ≤ 1.0
            @test snr_ok ? r.snr_lisa_kpc ≈ gw_snr(r.m1_msun, r.m2_msun, r.a_au, r.e; n_time = 3) :
                  isnan(r.snr_lisa_kpc)
        end

        ldc_names = [
            "EclipticLatitude",
            "EclipticLongitude",
            "Amplitude",
            "Frequency",
            "FrequencyDerivative",
            "Inclination",
            "Polarization",
            "InitialPhase",
            "ObservationDuration",
            "Cadence",
        ]
        dir = mktempdir()
        path = write_ldc_galactic_binaries(
            joinpath(dir, "ldc.csv"),
            rows;
            distance_kpc = 8.0,
            rng = StableRNG(2),
        )
        @test path == joinpath(dir, "ldc.csv")
        lines = readlines(path)
        @test startswith(lines[1], "#")
        @test occursin("n_pairs=3", lines[1]) && occursin("n_exported=2", lines[1])
        @test split(lines[2], ",") == ldc_names
        body = [split(l, ",") for l in lines[3:end]]
        @test length(body) == 2
        col(name) = findfirst(==(name), ldc_names)
        circular = filter(r -> r.e ≤ 0.01, rows)
        for (cells, r) in zip(body, sort(circular; by = r -> (r.name1, r.name2)))
            value(name) = parse(Float64, cells[col(name)])
            @test value("Amplitude") ≈ r.amp_ldc22_kpc / 8 rtol = 1e-12
            @test value("Frequency") == r.f_gw22_hz
            @test value("FrequencyDerivative") ≈ nd._gwcat_fdot_circular(r.mchirp_msun, r.f_gw22_hz) rtol =
                1e-12
            @test -π / 2 ≤ value("EclipticLatitude") ≤ π / 2
            @test 0 ≤ value("EclipticLongitude") < 2π
            @test 0 ≤ value("Inclination") ≤ π
            @test 0 ≤ value("Polarization") < π
            @test 0 ≤ value("InitialPhase") < 2π
            @test value("ObservationDuration") == 4 * 365.25 * 86400
            @test value("Cadence") == 15.0
        end

        # Without a generator the angle cells stay empty
        path = write_ldc_galactic_binaries(joinpath(dir, "ldc_norng.csv"), rows; distance_kpc = 8.0)
        lines = readlines(path)
        @test length(lines) == 2 + 2
        angles =
            ("EclipticLatitude", "EclipticLongitude", "Inclination", "Polarization", "InitialPhase")
        for l in lines[3:end]
            cells = split(l, ",")
            @test all(name -> cells[col(name)] == "", angles)
            @test !isempty(cells[col("Amplitude")])
        end

        # A looser cut exports every pair
        path = write_ldc_galactic_binaries(
            joinpath(dir, "ldc_all.csv"),
            rows;
            distance_kpc = 8.0,
            e_max = 0.6,
        )
        @test length(readlines(path)) == 2 + 3

        @test_throws ArgumentError write_ldc_galactic_binaries(
            joinpath(dir, "bad.csv"),
            rows;
            distance_kpc = 0.0,
        )
        @test_throws ArgumentError write_ldc_galactic_binaries(
            joinpath(dir, "bad.csv"),
            rows;
            distance_kpc = 8.0,
            e_max = 1.0,
        )
    end
end
