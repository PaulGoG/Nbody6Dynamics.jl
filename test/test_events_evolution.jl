# Readers of the gravitational-radiation, escape, mass-loss, tidal and Roche
# event lines (src/io/events_evolution.jl). The lines are synthetic: the
# layouts of the engine samples with other numbers, every field at the width
# of its FORMAT descriptor (fields built with lpad).

@testset "Evolution event lines" begin
    E = Nbody6Dynamics
    # Fields of width w, right-justified as Fortran writes them.
    fw(w, xs...) = join(lpad(string(x), w) for x in xs)
    kind_of(s) = E._evolution_line_kind(lstrip(s))
    parse_one(s; interval = 1) = E._parse_evolution_line(kind_of(s), [s], 1, interval)

    # Driver of the stdout interface: count ADJUST lines, classify, parse, skip
    # malformed records.
    function collect_records(lines)
        records = Any[]
        n_bad = 0
        interval = 1
        k = 1
        while k <= length(lines)
            s = lstrip(lines[k])
            if startswith(s, "ADJUST:")
                interval += 1
                k += 1
                continue
            end
            kind = E._evolution_line_kind(s)
            if kind === :none
                k += 1
                continue
            end
            try
                rec, n = E._parse_evolution_line(kind, lines, k, interval)
                push!(records, rec)
                k += n
            catch err
                err isa ArgumentError || rethrow()
                n_bad += 1
                k += 1
            end
        end
        return records, n_bad
    end

    # Shared middle block of the KSTIDE / KSINT2 / TIDES lines.
    gw_core =
        " I IP NPERT" *
        fw(10, 9811) *
        fw(6, 41, 3) *
        " NM1,2,S=" *
        fw(10, 3001, 3002, 13001) *
        " KW1,2,S=" *
        fw(4, 14, 13, -25) *
        " M1,2[M*]" *
        fw(13, "3.00000E+01", "1.40000E+00") *
        " e,a,QP,RS[NB]=" *
        fw(13, "9.90000E-01", "6.00000E-04", "6.00000E-06", "5.00000E-12") *
        " H, GAMMA=" *
        fw(13, "-1.25000E+01", "5.00000E-02") *
        " A_EIN, TORB, TGR(PM)=" *
        fw(13, "8.00000E-05", "7.50000E-04", "2.00000E+05")

    kstide =
        "  GR KSTIDE T DTGW STEP" *
        fw(13, "3.12500E+02", "2.50000E+01", "4.00000E-06") *
        gw_core *
        " dsemi,decc=" *
        fw(13, "-5.00000E-08", "-4.00000E-08") *
        " de,dh,oldQP=" *
        fw(13, "-3.00000E-06", "-1.00000E-03", "6.10000E-06")

    @testset "Kinds and classifier" begin
        @test length(E.EVOLUTION_EVENT_KINDS) == 22
        @test allunique(E.EVOLUTION_EVENT_KINDS)
        @test kind_of(kstide) === :gw_decay_step
        # Lines of the other family, the escape summary and coal.f lines.
        @test E._evolution_line_kind("") === :none
        @test E._evolution_line_kind(
            "NEW KSREG   TIME[NB] 0.0000000000E+00 NM1,2,S=        13        14     20013",
        ) === :none
        @test E._evolution_line_kind("END KSREG   TIME[NB] 5.7220458984E-06") === :none
        @test E._evolution_line_kind("CHAIN: NSTEP T   1234  1.0E+00") === :none
        @test E._evolution_line_kind("CHAIN ESCAPE: TT T TC =  1.1E+00") === :none
        @test E._evolution_line_kind(
            "ESCAPE T,NESC =  1.40000E+01     1 total N,NSESC,NBESC =      9988",
        ) === :none
        @test E._evolution_line_kind("ADJUST:  TIME  1.50000E+01") === :none
        @test E._evolution_line_kind("ROCHE COAL: IQCOLL   -2") === :none
        @test E._evolution_line_kind("BINARY COAL: IQCOLL   -2") === :none
        @test E._evolution_line_kind("BRAKE3    NM T E A DTGR RCOAL") === :none
        @test E._evolution_line_kind("NEW MERGER  TIME 1.57455E+00") === :none
        @test_throws ArgumentError E._parse_evolution_line(:unknown, [kstide], 1, 1)
        @test_throws ArgumentError E._parse_evolution_line(:gw_decay_step, [kstide], 2, 1)
    end

    @testset "GR KSTIDE" begin
        r, n = parse_one(kstide; interval = 3)
        @test n == 1
        @test r isa E.GWEvent
        @test r.kind === :gw_decay_step && r.line == 1 && r.interval == 3
        @test r.time_nb == 312.5 && r.dt_gw_nb == 25.0 && r.step_nb == 4.0e-6
        @test r.index_cm == 9811 && r.ipair == 41 && r.npert == 3
        @test (r.name1, r.name2, r.name_cm) == (3001, 3002, 13001)
        # The printed -25 is the temporary Peters flag, not the c.m. type.
        @test (r.kstar1, r.kstar2, r.kstar_cm) == (14, 13, KSTAR_ABSENT)
        @test r.m1_msun == 30.0 && r.m2_msun == 1.4
        @test r.ecc == 0.99 && r.semi_nb == 6.0e-4
        @test r.qperi_nb == 6.0e-6 && r.rschw_nb == 5.0e-12
        @test r.h_nb == -12.5 && r.gamma == 0.05
        @test r.a_ein == 8.0e-5 && r.t_orb_nb == 7.5e-4 && r.t_gr_nb == 2.0e5
        @test r.da_nb == -5.0e-8 && r.de == -4.0e-8
        @test r.d_energy_nb == -3.0e-6 && r.dh_nb == -1.0e-3 && r.qperi_old_nb == 6.1e-6
        @test isnan(r.t_gr_yr) && !r.overflow

        # Overflowed I6 IPAIR and E13.5 TGR.
        ovf = replace(kstide, "9811    41" => "9811******")
        ovf = replace(ovf, "  2.00000E+05" => "*"^13)
        @test length(ovf) == length(kstide)
        r = first(parse_one(ovf))
        @test r.overflow && r.ipair == 0 && isnan(r.t_gr_nb)
        @test r.npert == 3 && r.t_orb_nb == 7.5e-4

        # Malformed: a literal and a field content.
        @test_throws ArgumentError parse_one(replace(kstide, "I IP NPERT" => "I IP NPERX"))
        @test_throws ArgumentError parse_one(replace(kstide, "3.12500E+02" => "3.12500X+02"))
    end

    @testset "GR KSINT2 and GR TIDES" begin
        ksint =
            "  GR KSINT2 RES KMIN T DTGW STEP" *
            fw(3, 1, 2) *
            fw(13, "3.12500E+02", "2.50000E+01", "3.00000E-06") *
            gw_core
        r = first(parse_one(ksint))
        @test r.kind === :gw_ksint && r.res == 1 && r.kmin == 2
        @test r.time_nb == 312.5 && r.dt_gw_nb == 25.0 && r.step_nb == 3.0e-6
        @test r.index_cm == 9811 && r.ipair == 41 && r.npert == 3
        @test r.name_cm == 13001 && r.kstar_cm == KSTAR_ABSENT && r.kstar2 == 13
        @test r.h_nb == -12.5 && r.t_gr_nb == 2.0e5 && isnan(r.da_nb)

        tides =
            "  GR TIDES RES KMIN T DTGW STEP" *
            fw(4, 0, 3) *
            fw(13, "4.00000E+02", "1.25000E+01", "2.00000E-06") *
            gw_core
        r = first(parse_one(tides))
        @test kind_of(tides) === :gw_tides
        @test r.kind === :gw_tides && r.res == 0 && r.kmin == 3
        @test r.time_nb == 400.0 && r.dt_gw_nb == 12.5 && r.semi_nb == 6.0e-4
    end

    @testset "GR KSINT COLL and COLL DUE TO GW UNPERT (list-directed)" begin
        coll =
            "  GR KSINT COLL: QP SEPMRG C  RSCH ECC SEMI    2.5000000000000000E-010" *
            "   1.2500000000000000E-010   12345.678900000000        5.0000000000000000E-012" *
            "  0.50000000000000000        6.0000000000000000E-004"
        r = first(parse_one(coll))
        @test r.kind === :gw_ksint_coll && isnan(r.time_nb) && r.name1 == 0
        @test r.qperi_nb == 2.5e-10 && r.sep_merge_nb == 1.25e-10
        @test r.c_light_nb == 12345.6789 && r.rschw_nb == 5.0e-12
        @test r.ecc == 0.5 && r.semi_nb == 6.0e-4

        unpert =
            " COLL DUE TO GW UNPERT QPERI CLIGHT RSCH   2.5000000000000000E-010" *
            "   12345.678900000000        5.0000000000000000E-012"
        r = first(parse_one(unpert))
        @test r.kind === :gw_unpert_coll && isnan(r.time_nb)
        @test r.qperi_nb == 2.5e-10 && r.c_light_nb == 12345.6789 && r.rschw_nb == 5.0e-12
        # One copy per MPI rank: each copy is its own record.
        records, n_bad = collect_records([unpert, unpert])
        @test length(records) == 2 && n_bad == 0

        @test_throws ArgumentError parse_one(
            " COLL DUE TO GW UNPERT QPERI CLIGHT RSCH   2.5000000000000000E-010",
        )
    end

    @testset "GW recoil" begin
        recoil =
            " GW recoil:    " *
            " TIME[NB] I12 N12 K12 " *
            fw(17, "3.125000000E+02") *
            fw(10, 9811, 9812, 3001, 3002) *
            fw(4, 14, 14, 2) *
            " TIME[Myr] M012[M*] M12[M*] " *
            fw(13, "4.25000E+01", "3.00000E+01", "2.50000E+01", "2.90000E+01", "2.40000E+01") *
            " ASPN12,QR,XJSPIN12,SPNFAC " *
            fw(13, "6.80000E-01", "6.80000E-01", "8.27586E-01") *
            fw(13, "1.00000E-07", "2.00000E-07", "3.00000E+05") *
            " VKICK GW[km/s] (V,vx,vy,vz)" *
            fw(13, "1.50000E+02", "1.00000E+02", "-1.00000E+02", "5.00000E+01") *
            " VESC[km/s] VDIS[km/s] R/SEMI EB R " *
            fw(13, "2.00000E+01", "1.00000E+01", "1.00000E+00", "-2.00000E-03", "1.00000E-07") *
            " X1(1:3)[pc] V1(1:3)[km/s] X2(1:3)[pc] V2(1:3)[km/s] " *
            fw(13, "1.00000E+00", "2.00000E+00", "-3.00000E+00") *
            fw(13, "4.00000E+00", "5.00000E+00", "6.00000E+00") *
            fw(13, "1.50000E+00", "2.50000E+00", "-3.50000E+00") *
            fw(13, "-4.00000E+00", "-5.00000E+00", "-6.00000E+00")
        r = first(parse_one(recoil))
        @test r.kind === :gw_recoil && r.time_nb == 312.5 && r.time_myr == 42.5
        @test (r.index1, r.index2, r.name1, r.name2) == (9811, 9812, 3001, 3002)
        @test (r.kstar1, r.kstar2, r.icase) == (14, 14, 2)
        @test r.m01_msun == 30.0 && r.m02_msun == 25.0
        @test r.m1_msun == 29.0 && r.m2_msun == 24.0
        @test r.spin1 == 0.68 && r.q_ratio == 0.827586
        @test r.j_spin2_nb == 2.0e-7 && r.spin_factor == 3.0e5
        @test r.v_kick_kms == 150.0 && r.v_kick_xyz_kms == (100.0, -100.0, 50.0)
        @test r.v_esc_kms == 20.0 && r.v_dis_kms == 10.0
        @test r.eb_nb == -2.0e-3 && r.r_nb == 1.0e-7
        @test r.x1_pc == (1.0, 2.0, -3.0) && r.v1_kms == (4.0, 5.0, 6.0)
        @test r.x2_pc == (1.5, 2.5, -3.5) && r.v2_kms == (-4.0, -5.0, -6.0)
    end

    @testset "GR BRAKE (mdot.F, roche.f) and GR BRAKE3" begin
        brake =
            "  GR BRAKE T DTGR STEP" *
            fw(13, "2.00500E+02", "5.00000E+07", "3.00000E-05") *
            " I IP NPERT" *
            fw(10, 4401) *
            fw(6, 2201, 0) *
            " NM1,2,S=" *
            fw(10, 2201, 2202, 12201) *
            " KW1,2,S=" *
            fw(4, 8, 0, 0) *
            " M1,2[M*]" *
            fw(13, "1.50000E-01", "2.50000E+00") *
            " e,a,QP,RS[NB]=" *
            fw(13, "2.00000E-03", "2.50000E-08", "2.49500E-08", "2.50000E-13") *
            " H, GAMMA=" *
            fw(13, "-9.00000E+03", "0.00000E+00") *
            " A_EIN, TORB, TGR(PM)=" *
            fw(13, "9.50000E-05", "1.25000E-09", "2.00000E+09") *
            " j,dj=" *
            fw(13, "6.00000E+03", "3.00000E-06") *
            " DTXMIN,de=" *
            fw(13, "6.50000E+03", "3.25000E-12")
        r = first(parse_one(brake))
        @test r.kind === :gr_brake && r.time_nb == 200.5
        # Names and types reordered to the masses: the type-8 star is the 2.5 M☉ one.
        @test (r.name1, r.kstar1, r.m1_msun) == (2202, 0, 0.15)
        @test (r.name2, r.kstar2, r.m2_msun) == (2201, 8, 2.5)
        @test r.name_cm == 12201 && r.kstar_cm == 0
        @test r.index_star == 4401 && r.index_cm == 0 && r.ipair == 0
        @test r.dt_gr_yr == 5.0e7 && r.step_nb == 3.0e-5 && r.h_nb == -9000.0
        @test r.t_gr_yr == 2.0e9 && isnan(r.t_gr_nb)
        @test r.j_orb == 6000.0 && r.dj_dt == 3.0e-6
        @test r.dt_bse_yr == 6500.0 && r.de_dt_per_yr == 3.25e-12

        brake_roche =
            "  GR BRAKE T DTGW STEP" *
            fw(13, "2.00500E+02", "1.25000E+00", "3.00000E-05") *
            " I IP NPERT" *
            fw(10, 12201) *
            fw(6, 2201, 0) *
            " NM1,2,S=" *
            fw(10, 2201, 2202, 12201) *
            " KW1,2,S=" *
            fw(4, 3, 1, 13) *
            " M1,2[M*]" *
            fw(13, "2.50000E+00", "1.50000E-01") *
            " e,a,QP,RSCH[NB]=" *
            fw(13, "2.00000E-03", "2.50000E-08", "2.49500E-08", "2.50000E-13") *
            " H, GAMMA=" *
            fw(13, "-9.00000E+03", "0.00000E+00") *
            " A_EIN, TORB, TGR(PM)=" *
            fw(13, "9.50000E-05", "1.25000E-09", "4.00000E+04") *
            " j,dj=" *
            fw(13, "6.00000E+03", "1.50000E-01") *
            " DTM0,de=" *
            fw(13, "5.00000E+04", "1.00000E-07")
        r = first(parse_one(brake_roche))
        @test r.kind === :gr_brake_roche && r.index_cm == 12201 && r.ipair == 2201
        @test (r.name1, r.name2, r.name_cm) == (2201, 2202, 12201)
        @test (r.kstar1, r.kstar2, r.kstar_cm) == (3, 1, 13) && r.m1_msun == 2.5
        @test r.dt_gw_nb == 1.25 && r.t_gr_steps == 4.0e4 && isnan(r.t_gr_yr)
        @test r.dj == 0.15 && r.dt_bse_yr == 5.0e4 && r.de == 1.0e-7 && isnan(r.dj_dt)

        brake3 =
            " GR BRAKE3 T M1 M2 K1 K2 SEP ECC0 JORB DJ DTGR " *
            fw(14, "2.00500E+02", "1.00000E+01", "8.00000E+00") *
            fw(4, 14, 13) *
            fw(14, "3.00000E+01", "1.00000E-01", "1.00000E+04", "2.00000E-06", "5.00000E+07")
        r = first(parse_one(brake3))
        @test r.kind === :gr_brake3 && r.time_nb == 200.5
        @test r.m1_msun == 10.0 && r.m2_msun == 8.0 && (r.kstar1, r.kstar2) == (14, 13)
        @test r.semi_rsun == 30.0 && r.ecc == 0.1 && r.j_orb == 1.0e4
        @test r.dj_dt == 2.0e-6 && r.dt_gr_yr == 5.0e7
    end

    @testset "BINARY / HIARCH / QUAD ESCAPE" begin
        escape_line(name_cm, kstar_cm) =
            " BINARY ESCAPE  IP,N,NCORR =" *
            fw(5, 312) *
            fw(10, 9975) *
            fw(5, 2) *
            "  NM1,2,S=" *
            fw(10, 2047, 2048, name_cm) *
            "  KW1,2,S,LS=" *
            fw(4, 1, 0, kstar_cm, -1) *
            "  M1,2* =" *
            fw(13, "1.20000E+00", "8.00000E-01") *
            "  a,e,eb,pb =" *
            fw(13, "1.00000E-04", "5.00000E-01", "-4.00000E-06", "9.00000E+04") *
            " R*/PM =" *
            fw(13, "8.00000E-05") *
            "  RI,VI =" *
            fw(13, "1.50000E+00", "2.50000E+01") *
            "  EI =" *
            fw(13, "7.00000E-05") *
            "  RAD1,2[*]" *
            fw(13, "1.10000E+00", "7.50000E-01")
        summary = " ESCAPE T,NESC =  1.50000E+01     2 total N,NSESC,NBESC =      9975       0     2"
        adjust = " ADJUST:  TIME  1.50000E+01"
        lines = [
            "",
            escape_line(12047, 0),
            "",
            escape_line(-12047, 33),
            summary,
            adjust,
            "",
            escape_line(12047, 0),
            adjust,
            summary,
        ]

        r, n = E._parse_evolution_line(:binary_escape, lines, 2, 1)
        @test n == 1 && r isa E.BinaryEscapeEvent
        @test r.kind === :binary_escape && r.line == 2 && r.time_nb == 15.0
        @test (r.ipair, r.n_after, r.ncorr) == (312, 9975, 2)
        @test (r.name1, r.name2, r.name_cm) == (2047, 2048, 12047)
        @test (r.kstar1, r.kstar2, r.kstar_cm, r.ls_flag) == (1, 0, 0, -1)
        @test r.m1_msun == 1.2 && r.m2_msun == 0.8
        @test r.semi_nb == 1.0e-4 && r.ecc == 0.5 && r.eb_nb == -4.0e-6
        @test r.period_d == 9.0e4 && r.rstar_over_peri == 8.0e-5
        # Speed printed under RI, distance under VI.
        @test r.v_nb == 1.5 && r.r_nb == 25.0
        @test r.e_cm_nb == 7.0e-5 && r.radius1_rsun == 1.1 && r.radius2_rsun == 0.75

        # Hierarchy: K*S carries +30.
        r = first(E._parse_evolution_line(:binary_escape, lines, 4, 1))
        @test r.name_cm == -12047 && r.kstar_cm == 3 && r.time_nb == 15.0

        # No summary before the next ADJUST line.
        r = first(E._parse_evolution_line(:binary_escape, lines, 8, 2))
        @test isnan(r.time_nb)

        # The summary is not a record of this family and is not consumed.
        @test kind_of(summary) === :none
        records, n_bad = collect_records(lines)
        @test n_bad == 0 && length(records) == 3
        @test all(x -> x isa E.BinaryEscapeEvent, records)
        @test [x.line for x in records] == [2, 4, 8]
        @test [x.interval for x in records] == [1, 1, 2]
        @test records[1].time_nb == 15.0 && records[2].time_nb == 15.0
        @test isnan(records[3].time_nb)

        hiarch =
            " HIARCH ESCAPE    NM =" *
            fw(6, 2047, 2048) *
            "  K* =" *
            fw(3, 1, 0, 0) *
            "  M =" *
            fw(5, "1.2", "0.8") *
            "  EB =" *
            fw(8, "-0.0123") *
            "  E =" *
            fw(7, "0.250") *
            "  PM/PC =" *
            fw(8, "2.5E+00") *
            "  A =" *
            fw(8, "1.0E-05") *
            "  P =" *
            fw(8, "3.0E+02")
        r = first(E._parse_evolution_line(:hierarchy_escape, [hiarch, summary], 1, 1))
        @test kind_of(hiarch) === :hierarchy_escape
        @test r.kind === :hierarchy_escape && r.time_nb == 15.0
        @test (r.name1, r.name2) == (2047, 2048) && (r.kstar1, r.kstar2, r.kstar_cm) == (1, 0, 0)
        @test r.m1_msun == 1.2 && r.m2_msun == 0.8 && r.eb_nb == -0.0123
        @test r.ecc == 0.25 && r.pmin_over_pcrit == 2.5
        @test r.semi_nb == 1.0e-5 && r.period_d == 300.0 && !r.overflow

        # F5.1 mass overflowed.
        r = first(parse_one(replace(hiarch, "  M =  1.2" => "  M =*****")))
        @test r.overflow && isnan(r.m1_msun) && r.m2_msun == 0.8 && isnan(r.time_nb)

        quad =
            " QUAD ESCAPE    KS =" *
            fw(5, 57) *
            "  NM =" *
            fw(6, 111, 222) *
            "  K* =" *
            fw(3, -1, 5) *
            "  M =" *
            fw(5, "0.9", "0.6") *
            "  EB =" *
            fw(8, "-0.0050") *
            "  E =" *
            fw(7, "0.300") *
            "  A =" *
            fw(8, "2.0E-04") *
            "  P =" *
            fw(8, "1.5E+04")
        r = first(parse_one(quad))
        @test r.kind === :quad_escape && r.ipair == 57 && isnan(r.time_nb)
        # Names and KSTAR(JM) of this branch are stale: absent.
        @test r.name1 == 0 && r.name2 == 0 && r.kstar_cm == KSTAR_ABSENT
        @test r.ls_flag == -1 && r.m1_msun == 0.9 && r.m2_msun == 0.6
        @test r.eb_nb == -0.005 && r.ecc == 0.3 && r.semi_nb == 2.0e-4 && r.period_d == 1.5e4
    end

    @testset "MDEND / MDNEW KSREG and BINARY KICK" begin
        mdot(label, name_cm, kstar_cm) =
            label *
            fw(17, "6.2500000000E+00") *
            " NM1,2,S=" *
            fw(10, 1801, 1802, name_cm) *
            " KW1,2,S=" *
            fw(4, 1, 0, kstar_cm) *
            " IPAIR" *
            fw(9, 3500) *
            " DTAU" *
            fw(11, "7.500E-02") *
            " M1,2[NB]" *
            fw(11, "2.000E-03", "4.000E-05") *
            " R12[NB]" *
            fw(11, "9.500E-04") *
            " e,a,eb[NB]=" *
            fw(12, "5.5000E-01", "1.2000E-03") *
            fw(11, "-3.500E-05") *
            " P[d]=" *
            fw(11, "2.700E+05") *
            " H" *
            fw(11, "-9.000E-01") *
            " GAMMA" *
            fw(11, "1.000E-05") *
            " STEP(ICM)" *
            fw(11, "1.953E-03") *
            " NPERT" *
            fw(5, 4) *
            " NB(ICM)" *
            fw(5, 220) *
            " M1,2[*]" *
            fw(11, "2.400E+01", "5.000E-01") *
            " RAD1,2,S[*]" *
            fw(11, "7.000E+00", "4.500E-01", "4.300E+04") *
            " RI,VI[NB]=" *
            fw(11, "4.800E-01", "4.500E-01") *
            " DM,DMR=" *
            fw(13, "2.50000E-06", "1.25000E-03")

        mdend = mdot(" MDEND KSREG TIME[NB]", 21801, 0)
        r = first(parse_one(mdend))
        @test r isa E.MassLossEvent && r.kind === :ks_end_mass_loss && r.time_nb == 6.25
        @test (r.name1, r.name2, r.name_cm) == (1801, 1802, 21801)
        @test (r.kstar1, r.kstar2, r.kstar_cm) == (1, 0, 0)
        @test r.ipair == 3500 && r.dtau_nb == 0.075
        @test r.m1_nb == 2.0e-3 && r.m2_nb == 4.0e-5 && r.r12_nb == 9.5e-4
        @test r.ecc == 0.55 && r.semi_nb == 1.2e-3 && r.eb_nb == -3.5e-5
        @test r.period_d == 2.7e5 && r.h_nb == -0.9 && r.gamma == 1.0e-5
        @test r.step_cm_nb == 1.953e-3 && r.npert == 4 && r.n_neighbours == 220
        @test r.m1_msun == 24.0 && r.m2_msun == 0.5
        @test r.radius1_rsun == 7.0 && r.radius2_rsun == 0.45
        # Third value under RAD1,2,S[*] is the separation.
        @test r.r12_rsun == 4.3e4
        @test r.r_nb == 0.48 && r.v_nb == 0.45 && r.dm_nb == 2.5e-6 && r.dmr == 1.25e-3

        mdnew = mdot(" MDNEW KSREG TIME[NB]", 1801, 1)
        r = first(parse_one(mdnew))
        @test r.kind === :ks_new_mass_loss && r.time_nb == 6.25
        @test (r.name1, r.name2, r.kstar1, r.kstar2) == (1801, 1802, 1, 0)
        @test r.m1_nb == 2.0e-3 && r.m1_msun == 24.0 && r.radius2_rsun == 0.45
        @test r.dm_nb == 2.5e-6 && r.dmr == 1.25e-3
        # Stale fields of MDNEW stored as absent.
        @test r.name_cm == 0 && r.kstar_cm == KSTAR_ABSENT && r.ipair == 0
        @test isnan(r.ecc) && isnan(r.semi_nb) && isnan(r.eb_nb) && isnan(r.period_d)
        @test isnan(r.dtau_nb) && isnan(r.r12_nb) && isnan(r.h_nb) && isnan(r.gamma)
        @test isnan(r.step_cm_nb) && r.npert == 0 && r.n_neighbours == 0
        @test isnan(r.r12_rsun) && isnan(r.r_nb) && isnan(r.v_nb)

        kick =
            " BINARY KICK:  TIME[Myr] N12 K12 " *
            fw(17, "1.234500000E+01") *
            fw(10, 501, 502) *
            fw(4, 12, 1) *
            " ASPN12,QR,XJSPN12[NB],SPNFAC " *
            fw(13, "0.00000E+00", "0.00000E+00", "2.00000E+00") *
            fw(13, "0.00000E+00", "0.00000E+00", "1.00000E+05") *
            " M12[M*] " *
            " VESC[km/s] VDIS[km/s] R12/SEMI EB R12[au] " *
            fw(13, "1.20000E+00", "6.00000E-01", "1.50000E+01", "1.40000E+01") *
            fw(13, "9.00000E-01", "-2.00000E-05", "3.00000E-04") *
            " X1(1:3)[pc] V1(1:3)[km/s] X2(1:3)[pc] V2(1:3)[km/s] " *
            fw(13, "1.00000E+00", "2.00000E+00", "-3.00000E+00") *
            fw(13, "4.00000E+00", "5.00000E+00", "6.00000E+00") *
            fw(13, "1.50000E+00", "2.50000E+00", "-3.50000E+00") *
            fw(13, "-4.00000E+00", "-5.00000E+00", "-6.00000E+00")
        r = first(parse_one(kick))
        @test r.kind === :binary_kick && isnan(r.time_nb) && r.time_myr == 12.345
        @test (r.name1, r.name2, r.kstar1, r.kstar2) == (501, 502, 12, 1)
        @test r.q_ratio == 2.0 && r.spin_factor == 1.0e5
        @test r.m1_msun == 1.2 && r.m2_msun == 0.6
        @test r.v_esc_kms == 15.0 && r.v_dis_kms == 14.0 && r.r_over_semi == 0.9
        # Labelled R12[au], printed in N-body units.
        @test r.eb_nb == -2.0e-5 && r.r12_nb == 3.0e-4
        @test r.x1_pc == (1.0, 2.0, -3.0) && r.v2_kms == (-4.0, -5.0, -6.0)
    end

    @testset "NEW / END CIRC, CIRCULARIZED, INACTIIVE PHASE, NEW CAPTURE" begin
        circ(label, kstar_cm, ecc_after, ri, vi) =
            label *
            fw(13, "7.25000E+00", "8.00000E-03") *
            " I IP NPERT" *
            fw(10, 21700) *
            fw(6, 1700, 1) *
            " NM1,2,S=" *
            fw(10, 6601, 6602, 26601) *
            " KW1,2,S=" *
            fw(4, 1, 0, kstar_cm) *
            " M1,2[M*]" *
            fw(13, "2.00000E+00", "2.50000E-01") *
            " e0,e,em,a0,a,QP=" *
            fw(13, "8.75000E-01", ecc_after, "2.00000E-03") *
            fw(13, "1.00000E-06", "2.50000E-07", "1.25000E-07") *
            " H, GAMMA=" *
            fw(13, "-9.25000E+01", "0.00000E+00") *
            " PD, R1,2,SEP,RCOLL[R*]=" *
            fw(13, "2.75000E+00", "1.50000E+00", "2.25000E-01", "4.50000E+00", "1.25000E+00") *
            " RI,VI=" *
            fw(13, ri, vi)

        new_circ = circ(" NEW CIRC T DT ", 9, "3.75000E-01", "1.12500E+00", "6.87500E-01")
        r = first(parse_one(new_circ))
        @test r isa E.TidalEvent && r.kind === :circ_new && !r.inactive_phase
        @test r.time_nb == 7.25 && r.step_nb == 8.0e-3
        @test (r.index_cm, r.ipair, r.npert) == (21700, 1700, 1)
        @test (r.name1, r.name2, r.name_cm) == (6601, 6602, 26601)
        @test (r.kstar1, r.kstar2, r.kstar_cm) == (1, 0, 9)
        @test r.m1_msun == 2.0 && r.m2_msun == 0.25
        @test r.ecc0 == 0.875 && r.ecc == 0.375 && r.ecc_max == 0.002
        @test r.semi0_nb == 1.0e-6 && r.semi_nb == 2.5e-7 && r.qperi_nb == 1.25e-7
        @test r.h0_nb == -92.5 && r.gamma == 0.0
        @test r.period_d == 2.75 && isnan(r.period0_d)
        @test r.radius1_rsun == 1.5 && r.radius2_rsun == 0.225
        @test r.r12_rsun == 4.5 && r.r_coll_rsun == 1.25
        @test r.r_nb == 1.125 && r.v_nb == 0.6875

        end_circ = circ(" END CIRC T DT ", 10, "1.99750E-03", "1.12500E+00", "6.87500E-01")
        r = first(parse_one(end_circ))
        @test r.kind === :circ_end && r.kstar_cm == 10 && r.ecc == 1.9975e-3
        @test r.period_d == 2.75

        # D exponents in RI, VI; PD is the period before the step.
        circd = circ("  CIRCULARIZED T DT ", 10, "1.99750E-03", "3.25000D-01", "8.75000D-01")
        r = first(parse_one(circd))
        @test r.kind === :circularized && !r.inactive_phase
        @test r.r_nb == 0.325 && r.v_nb == 0.875
        @test r.period0_d == 2.75 && isnan(r.period_d)

        inactive = circ("  INACTIIVE PHASE T DT ", 9, "1.50000E-03", "3.25000D-01", "8.75000D-01")
        r = first(parse_one(inactive))
        @test r.kind === :circularized && r.inactive_phase
        @test r.ecc == 1.5e-3 && r.r_nb == 0.325
        @test isnan(r.period_d) && isnan(r.period0_d)

        capture =
            " NEW CAPTURE    T NM E EF QP/R* A1 " *
            fw(9, "12.50") *
            fw(6, 701, 702) *
            fw(9, "1.0500", "0.9500") *
            fw(11, "3.500E+00", "2.000E-05")
        r = first(parse_one(capture))
        @test r.kind === :capture && r.time_nb == 12.5
        @test (r.name1, r.name2, r.name_cm) == (701, 702, 0) && r.kstar1 == KSTAR_ABSENT
        @test r.ecc0 == 1.05 && r.ecc == 0.95
        @test r.qperi_over_rstar == 3.5 && r.semi_nb == 2.0e-5
    end

    @testset "NEW ROCHE, ROCHE COAL, END ROCHE" begin
        roche_head(which, kstar2) =
            which *
            "  TIME[NB]" *
            fw(17, "1.7250000000E+01") *
            " NM1,2,S=" *
            fw(10, 1901, 3401, 11901) *
            " KW1,2,S=" *
            fw(4, 1, kstar2, 11) *
            " IPAIR" *
            fw(9, 1500) *
            " DTAU" *
            fw(13, "8.50000E-04")
        roche_body =
            " e,a,eb[NB]=" *
            fw(13, "2.00000E-03", "6.60000E-07", "-2.90000E-01") *
            " P[d]=" *
            fw(13, "2.00000E+00") *
            " H" *
            fw(13, "-1.12500E+04") *
            " GAMMA" *
            fw(13, "0.00000E+00") *
            " STEP(ICM)" *
            fw(13, "4.88281E-04") *
            " NPERT" *
            fw(5, 0) *
            " NB(ICM)" *
            fw(5, 112) *
            " M0(1:2) M(1:2) MC(1:2)[*]" *
            fw(13, "8.75000E+01", "1.50000E-01", "8.70000E+01") *
            fw(13, "1.50000E-01", "2.50000E+01", "0.00000E+00") *
            " RI,VI[NB]=" *
            fw(13, "9.50000E-02", "5.25000E-01") *
            " SEP[*]" *
            fw(13, "2.97500E+01") *
            " RAD(1:2) ROL(1:2)[*]=" *
            fw(13, "2.45000E+01", "1.75000E-01", "2.30000E+01", "1.72500E+00")

        new_roche =
            roche_head(" NEW ROCHE  ", 0) *
            " TIME,AGE[Myr] " *
            fw(13, "3.50000E+00", "4.50000E-01") *
            " M1,2[NB]" *
            fw(13, "1.50000E-02", "2.50000E-05") *
            " SPIN12,OSPIN12[NB] SPNFAC JORB,OORB[*] " *
            fw(13, "8.00000E-09", "1.00000E-16", "5.00000E-12", "2.00000E-13") *
            fw(13, "6.00000E+13", "1.50000E+05", "1.12500E+03") *
            " R12[NB]" *
            fw(13, "6.50000E-07") *
            roche_body *
            " SEMI[*],RL1[*]=" *
            fw(13, "2.98000E+01", "2.30000E+01")
        r = first(parse_one(new_roche))
        @test r isa E.RocheEvent && r.kind === :roche_new && r.time_nb == 17.25
        @test (r.name1, r.name2, r.name_cm) == (1901, 3401, 11901)
        @test (r.kstar1, r.kstar2, r.kstar_cm) == (1, 0, 11)
        @test r.ipair == 1500 && r.dtau_nb == 8.5e-4
        @test r.time_myr == 3.5 && r.age_myr == 0.45
        @test r.m1_nb == 0.015 && r.m2_nb == 2.5e-5
        @test r.spin1_nb == 8.0e-9 && r.ospin2_nb == 2.0e-13 && r.spin_factor == 6.0e13
        @test r.j_orb == 1.5e5 && r.omega_orb_per_yr == 1125.0 && r.r12_nb == 6.5e-7
        @test r.ecc == 0.002 && r.semi_nb == 6.6e-7 && r.eb_nb == -0.29
        @test r.period_d == 2.0 && r.h_nb == -11250.0 && r.step_cm_nb == 4.88281e-4
        @test r.npert == 0 && r.n_neighbours == 112
        @test r.m01_msun == 87.5 && r.m1_msun == 87.0 && r.mc1_msun == 25.0 && r.mc2_msun == 0.0
        @test r.r_nb == 0.095 && r.v_nb == 0.525
        # SEP[*] is the separation, SEMI[*] the semi-major axis.
        @test r.r12_rsun == 29.75 && r.semi_rsun == 29.8 && r.rl1_rsun == 23.0
        @test r.radius1_rsun == 24.5 && r.roche_radius1_rsun == 23.0
        @test r.roche_radius2_rsun == 1.725

        coal_body(which, kstar2) =
            roche_head(which, kstar2) *
            " T[Myr] " *
            fw(13, "3.50000E+00") *
            " M1,2[NB]" *
            fw(13, "1.50000E-02", "2.50000E-05") *
            " R12[NB]" *
            fw(13, "6.50000E-07") *
            roche_body *
            " DM1/2,DT=" *
            fw(13, "0.00000E+00", "1.00000E-03", "1.00000E-01")
        # The logical COALS glued to the preceding E13.5.
        roche_coal =
            coal_body(" ROCHE COAL ", 15) *
            " AGE COALS, IXXX=" *
            fw(13, "4.50000E-01") *
            "T" *
            fw(3, 8)
        @test endswith(roche_coal, "4.50000E-01T  8")
        r = first(parse_one(roche_coal))
        @test r.kind === :roche_coal && r.time_nb == 17.25 && r.kstar2 == 15
        @test r.time_myr == 3.5 && r.age_myr == 0.45
        @test r.coalescence && r.ixxx == 8
        @test r.dm1_msun == 0.0 && r.dm2_msun == 1.0e-3 && r.dt_myr == 0.1
        @test r.ecc == 0.002 && r.r12_rsun == 29.75 && r.m1_msun == 87.0
        # TK = 0 in ROCHE COAL; EB, RI, VI not formed on this path.
        @test isnan(r.period_d) && isnan(r.eb_nb) && isnan(r.r_nb) && isnan(r.v_nb)
        @test isnan(r.semi_rsun)

        r = first(parse_one(replace(roche_coal, "E-01T" => "E-01F")))
        @test !r.coalescence && r.ixxx == 8
        @test_throws ArgumentError parse_one(replace(roche_coal, "E-01T" => "E-01X"))

        end_roche = coal_body(" END ROCHE  ", 0)
        r = first(parse_one(end_roche))
        @test kind_of(end_roche) === :roche_end
        @test r.kind === :roche_end && r.period_d == 2.0 && isnan(r.eb_nb)
        @test isnan(r.age_myr) && !r.coalescence && r.ixxx == 0 && r.dt_myr == 0.1
    end

    @testset "Driver: malformed line skipped, empty input" begin
        bad = replace(kstide, "I IP NPERT" => "I IP NPERX")
        records, n_bad = collect_records(["", kstide, "", bad, " ADJUST: T", kstide])
        @test n_bad == 1 && length(records) == 2
        @test [r.line for r in records] == [2, 6]
        @test [r.interval for r in records] == [1, 2]

        records, n_bad = collect_records(String[])
        @test isempty(records) && n_bad == 0
        path = joinpath(mktempdir(), "out1000")
        write(path, "")
        records, n_bad = collect_records(readlines(path))
        @test isempty(records) && n_bad == 0
    end
end
