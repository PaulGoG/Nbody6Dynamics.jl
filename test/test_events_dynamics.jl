# Tests for the KS, few-body, hierarchy, exchange and termination lines of the
# engine log (src/io/events_dynamics.jl). Every line is synthetic: the layout
# of the engine FORMAT with altered values, each field at its width.

_tdyn_kind(line::AbstractString) = Nbody6Dynamics._dynamics_line_kind(lstrip(line))
_tdyn_parse(kind::Symbol, lines, i::Int) = Nbody6Dynamics._parse_dynamics_line(kind, lines, i, 3)

# kspinit.f FORMAT 60 (28 values)
_tdyn_ks_primordial() = " NEW KSREG   TIME[NB] 0.0000000000E+00 NM1,2,S=        27        28     20027 KW1,2,S=   1   0   0 IPAIR        5 DTAU  2.345E-01 M1,2[NB]  5.111E-05  3.222E-05 R12[NB]  7.654E-04 e,a,eb[NB]=  6.1234E-01  4.1111E-04 -2.222E-06 P[d]=  3.444E+05 H -9.999E-02 GAMMA  1.500E-03 STEP(ICM)  1.953E-03 NPERT    2 NB(ICM)  117 M1,2[*]  6.123E-01  3.871E-01 RAD1,2,S[*]  5.987E-01  3.996E-01  3.456E+04 RI,VI[NB]=  2.718E+00  3.141E-01"

# ksinit.F FORMAT 60 (40 values)
_tdyn_ks_new() = " NEW KSREG   TIME[NB] 1.2345678901E-01 NM1,2,S=        91      6641     20091 KW1,2,S=   0   1   0 IPAIR     4123 DTAU  1.111E-01 M1,2[NB]  2.074E-05  1.592E-04 R12[NB]  5.726E-04 e,a,eb[NB]=  1.4790E+00 -7.4950E-04  2.344E-06 P[d]= -5.111E+05 H  2.354E-01 GAMMA  6.452E-08 STEP(ICM)  9.537E-07 NPERT    3 NB(ICM)  133 M1,2[*]  2.370E-01  1.878E+00 RAD1,2,S[*]  2.741E-01  1.644E+00  3.152E+04 RCM,VCM[NB]=  7.91570E-01  5.33564E-01 I1 RXY[NB],i,alp[rad]   7.51641E-01 -4.56317E-01  1.84464E+00 I1 VXY[NB],i,alp[rad]   4.97994E-01  6.24533E-02  8.63672E-01 I2 RXY[NB],i,alp[rad]   7.51556E-01 -4.55731E-01  1.84465E+00 I2 VXY[NB],i,alp[rad]   5.25654E-01 -1.05358E+00  3.09319E+00"

# ksterm.F FORMAT 15 (40 values)
_tdyn_ks_end() = " END KSREG   TIME[NB] 6.8331569095E-06 NM1,2,S=      4212      4213     24212 KW1,2,S=   0   0   0 IPAIR      811 DTAU  1.082E-02 M1,2[NB]  3.258E-05  7.383E-06 R12[NB]  1.088E-09 e,a,eb[NB]=  9.8878E-01  5.7809E-07 -2.570E-04 P[d]=  1.582E+01 H -4.269E+01 GAMMA  0.000E+00 STEP(ICM)  9.766E-04 NPERT    0 NB(ICM)  124 M1,2[*]  3.888E-01  8.811E-02 RAD1,2,S[*]  3.652E-01  1.404E-01  5.924E-02 RCM,VCM[NB]=  5.08786E-01  1.01705E+00 I1 RXY[NB],i,alp[rad]   1.59973E-01  1.25095E+00  7.51556E-01 I1 VXY[NB],i,alp[rad]   9.28947E-01  2.90952E-01  4.92748E-01 I2 RXY[NB],i,alp[rad]   1.59972E-01  1.25095E+00  7.51559E-01 I2 VXY[NB],i,alp[rad]   1.19292E+00  2.68167E-01  5.38604E-01"

# impact.f FORMAT 20 (NEW CHAIN; the other labels are substituted)
_tdyn_impact() = " NEW CHAIN   TIME 9.40626E-02 INCM    23118 I3    12926 INPAIR     3118 N1,2,3,INCM      9184      9185     11876     29184 KW1,2,3,IN   1   0   1   0 M1,2[NB]  8.000E-05  6.000E-05 M3,TOT  1.400E-04  2.800E-04  R12  2.043E-04 H -7.146E-01 IN e,a,eb[NB]  5.6062E-01  1.2671E-04 -3.176E-05 P[d]  2.000E+04 OUT e,a,eb[NB]  1.1327E+00 -2.4202E-02  8.331E-07 P[d] -4.000E+07  G4 5.780E-03  RIN3 5.4E-04  PERIM 3.1E-04  EB1/EB0 -4.6E-02  NP   4 M1,2,3,TOT[*]  9.596E-01  7.167E-01  1.629E+00  3.294E+00 RAD1,2,3[*]  8.400E-01  6.659E-01  1.490E+00 IN,OUT Sep[*]  9.657E+03  1.921E+04  RI,VI[NB]  4.047E-01  1.426E+00"

@testset "Dynamics events — classification" begin
    kinds = Nbody6Dynamics.DYNAMICS_EVENT_KINDS
    @test length(kinds) == 36 && allunique(kinds)
    @test _tdyn_kind(_tdyn_ks_primordial()) === :ks_new_primordial
    @test _tdyn_kind(_tdyn_ks_new()) === :ks_new
    @test _tdyn_kind(_tdyn_ks_end()) === :ks_end
    @test _tdyn_kind(_tdyn_impact()) === :chain_new
    # Per-step chain diagnostics, other families and engine-event lines.
    @test _tdyn_kind(
        " CHAIN: NSTEP T DTR R NAMEC    41    0.1953   -0.0000      9184  1.148E-03      9185",
    ) === :none
    @test _tdyn_kind(
        " CHAIN T,TC[NB] K,N,M,R[*],R,V,a,e,tgr,aein[NB-CH]= 1.9531250000E-01   1.84619E-02",
    ) === :none
    @test _tdyn_kind(
        "  GR KSTIDE T DTGW STEP  1.23456E+02  2.50000E+01  4.00000E-06 I IP NPERT      8801",
    ) === :none
    @test _tdyn_kind(" ADJUST:  TIME[NB]   2.00000E+00  TIME[Myr]   2.72E-01") === :none
    @test _tdyn_kind(" MDNEW KSREG TIME[NB] 1.0000000000E+00 NM1,2,S=") === :none
    @test _tdyn_kind(" MDEND KSREG TIME[NB] 1.0000000000E+00 NM1,2,S=") === :none
    @test _tdyn_kind(" BEG KSRECT TIME[NB] 1.0000000000E+00 NM1,2,S=") === :none
    @test _tdyn_kind(" RECOIL:    NAM E E1 PM/A RD A A1 RP GB IN   9185  0.7901") === :none
    @test _tdyn_kind(" CHAIN ESCAPE: TT T TC =  2.2799820E-06  5.6496631E+00") === :none
    @test _tdyn_kind("") === :none
    # No lines at all: the parser refuses the index.
    @test_throws ArgumentError _tdyn_parse(:ks_new, String[], 1)
    @test_throws ArgumentError _tdyn_parse(:no_such_kind, [_tdyn_ks_new()], 1)
end

@testset "Dynamics events — KS regularisation" begin
    nd = Nbody6Dynamics

    ev, n = _tdyn_parse(:ks_new_primordial, [_tdyn_ks_primordial()], 1)
    @test n == 1 && ev isa nd.KSEvent
    @test ev.kind === :ks_new_primordial && ev.line == 1 && ev.interval == 3
    @test ev.time_nb == 0.0
    @test (ev.name1, ev.name2, ev.name_cm) == (27, 28, 20027)
    @test (ev.kstar1, ev.kstar2, ev.kstar_cm) == (1, 0, 0)
    @test ev.ipair == 5 && ev.dtau == 0.2345
    @test ev.m1_nb == 5.111e-5 && ev.m2_nb == 3.222e-5 && ev.r12_nb == 7.654e-4
    @test ev.ecc == 0.61234 && ev.semi_nb == 4.1111e-4 && ev.eb_nb == -2.222e-6
    @test ev.period_d == 3.444e5 && ev.h_nb == -0.09999 && ev.gamma == 1.5e-3
    @test ev.step_cm_nb == 1.953e-3 && ev.n_perturbers == 2 && ev.n_neighbours_cm == 117
    @test ev.m1_msun == 0.6123 && ev.m2_msun == 0.3871
    @test ev.radius1_rsun == 0.5987 && ev.radius2_rsun == 0.3996
    # FORMAT field 26 (label RAD1,2,S[*]) is the separation in R☉.
    @test ev.r12_rsun == 3.456e4
    @test ev.r_cm_nb == 2.718 && ev.v_cm_nb == 0.3141
    @test isnan(ev.r1_xy_nb) && isnan(ev.v2_azi) && isnan(ev.deb_nb) && !ev.overflow

    ev, _ = _tdyn_parse(:ks_new, [_tdyn_ks_new()], 1)
    @test ev.kind === :ks_new && ev.time_nb == 0.12345678901
    @test (ev.name1, ev.name2, ev.name_cm) == (91, 6641, 20091)
    @test (ev.kstar1, ev.kstar2, ev.kstar_cm) == (0, 1, 0)
    @test ev.ipair == 4123 && ev.ecc == 1.479 && ev.semi_nb == -7.495e-4
    @test ev.eb_nb == 2.344e-6 && ev.period_d == -5.111e5 && ev.gamma == 6.452e-8
    @test ev.n_perturbers == 3 && ev.n_neighbours_cm == 133
    @test ev.r12_rsun == 3.152e4 && ev.radius2_rsun == 1.644
    @test ev.r_cm_nb == 0.79157 && ev.v_cm_nb == 0.533564
    @test ev.r1_xy_nb == 0.751641 && ev.r1_lat == -0.456317 && ev.r1_azi == 1.84464
    @test ev.v1_xy_nb == 0.497994 && ev.v1_lat == 0.0624533 && ev.v1_azi == 0.863672
    @test ev.r2_xy_nb == 0.751556 && ev.r2_lat == -0.455731 && ev.r2_azi == 1.84465
    @test ev.v2_xy_nb == 0.525654 && ev.v2_lat == -1.05358 && ev.v2_azi == 3.09319
    @test isnan(ev.deb_nb) && !ev.overflow

    ev, _ = _tdyn_parse(:ks_end, [_tdyn_ks_end()], 1)
    @test ev.kind === :ks_end && ev.time_nb == 6.8331569095e-6
    @test (ev.name1, ev.name2, ev.name_cm, ev.ipair) == (4212, 4213, 24212, 811)
    @test ev.ecc == 0.98878 && ev.semi_nb == 5.7809e-7 && ev.eb_nb == -2.57e-4
    @test ev.r12_rsun == 5.924e-2 && ev.radius1_rsun == 0.3652
    @test ev.r_cm_nb == 0.508786 && ev.v2_azi == 0.538604

    # ksinit.F FORMAT 65: NS/BH label, 26 core values, then RI, VI, DEB.
    compact = replace(
        _tdyn_ks_primordial(),
        " NEW KSREG   TIME[NB]" => " NS/BH BINARY T[NB]",
        " KW1,2,S=   1   0   0" => " KW1,2,S=  14  13   0",
        " RI,VI[NB]=  2.718E+00  3.141E-01" => " RI,VI[NB],DEB=  2.718E+00  3.141E-01  1.500E-01",
    )
    @test _tdyn_kind(compact) === :ks_new_compact
    ev, _ = _tdyn_parse(:ks_new_compact, [compact], 1)
    @test ev.kind === :ks_new_compact && (ev.kstar1, ev.kstar2) == (14, 13)
    @test ev.name_cm == 20027 && ev.r12_rsun == 3.456e4
    @test ev.r_cm_nb == 2.718 && ev.v_cm_nb == 0.3141 && ev.deb_nb == 0.15
    @test isnan(ev.r1_xy_nb)

    # A ten-digit name fills its I10 field and touches the label.
    glued = replace(_tdyn_ks_new(), " NM1,2,S=        91" => " NM1,2,S=1234567890")
    ev, _ = _tdyn_parse(:ks_new, [glued], 1)
    @test ev.name1 == 1234567890 && ev.name2 == 6641 && ev.name_cm == 20091 && !ev.overflow

    # An overflowed name reads as absent and flags the record.
    over = replace(_tdyn_ks_new(), "      6641     20091" => "**********     20091")
    ev, _ = _tdyn_parse(:ks_new, [over], 1)
    @test ev.overflow && ev.name2 == 0 && ev.name1 == 91 && ev.name_cm == 20091

    # Malformed: a damaged label, and the long layout asked of a short record.
    damaged = replace(_tdyn_ks_end(), " DTAU" => " DTAX")
    @test_throws ArgumentError _tdyn_parse(:ks_end, [damaged], 1)
    @test_throws ArgumentError _tdyn_parse(:ks_new, [_tdyn_ks_primordial()], 1)
end

@testset "Dynamics events — impact decisions and new mergers" begin
    nd = Nbody6Dynamics
    line = _tdyn_impact()
    ev, n = _tdyn_parse(:chain_new, [line], 1)
    @test n == 1 && ev isa nd.FewBodyEvent && ev.kind === :chain_new
    @test ev.time_nb == 9.40626e-2 && ev.n_members == 0
    @test ev.indices == (23118, 12926, 3118, 0)
    @test ev.names == (9184, 9185, 11876, 29184, 0, 0)
    @test ev.kstar_values == (1, 0, 1, 0)
    @test ev.masses_nb == (8.0e-5, 6.0e-5, 1.4e-4, 2.8e-4)
    @test ev.r12_nb == 2.043e-4 && ev.h_nb == -0.7146
    @test ev.ecc_in == 0.56062 && ev.semi_in_nb == 1.2671e-4 && ev.eb_in_nb == -3.176e-5
    @test ev.ecc_out == 1.1327 && ev.semi_out_nb == -2.4202e-2 && ev.eb_out_nb == 8.331e-7
    # Printed periods use the total mass (inner) and the inner mass (outer).
    m_in = 8.0e-5 + 6.0e-5
    @test ev.period_in_d ≈ 2.0e4 * sqrt(2.8e-4 / m_in) rtol = 1e-12
    @test ev.period_out_d ≈ -4.0e7 * sqrt(m_in / 2.8e-4) rtol = 1e-12
    @test ev.period_in_d ≈ 2.0e4 * sqrt(2.0) rtol = 1e-9
    @test ev.g4 == 5.78e-3 && ev.rij_nb == 5.4e-4 && ev.pmin_nb == 3.1e-4
    @test ev.eb_ratio == -4.6e-2 && ev.n_perturbers == 4
    @test ev.masses_msun == (0.9596, 0.7167, 1.629, 3.294)
    @test ev.radii_rsun == (0.84, 0.6659, 1.49)
    @test ev.r12_rsun == 9.657e3 && ev.rij_rsun == 1.921e4
    @test ev.r_cm_nb == 0.4047 && ev.v_cm_nb == 1.426 && !ev.overflow

    for (label, kind, type, n_members) in (
        (" NEW TRIPLE  TIME", :triple_new, nd.FewBodyEvent, 3),
        (" NEW QUAD    TIME", :quad_new, nd.FewBodyEvent, 4),
        (" NEW MERGER  TIME", :hierarchy_new, nd.HierarchyEvent, 0),
        (" NEW MERGE2  TIME", :hierarchy_new_merge2, nd.HierarchyEvent, 0),
    )
        variant = replace(line, " NEW CHAIN   TIME" => label)
        @test _tdyn_kind(variant) === kind
        ev, _ = _tdyn_parse(kind, [variant], 1)
        @test ev isa type && ev.kind === kind && ev.names[4] == 29184
        @test ev.indices == (23118, 12926, 3118, 0) && ev.semi_out_nb == -2.4202e-2
        @test ev.period_in_d ≈ 2.0e4 * sqrt(2.8e-4 / m_in) rtol = 1e-12
        if type === nd.FewBodyEvent
            @test ev.n_members == n_members
        end
    end

    # E8.1 overflow of the pericentre.
    over = replace(line, "  PERIM 3.1E-04" => "  PERIM********")
    ev, _ = _tdyn_parse(:chain_new, [over], 1)
    @test ev.overflow && isnan(ev.pmin_nb) && ev.eb_ratio == -4.6e-2
    damaged = replace(line, " INPAIR" => " INPAIX")
    @test_throws ArgumentError _tdyn_parse(:chain_new, [damaged], 1)
end

@testset "Dynamics events — chain start, termination and exchange" begin
    nd = Nbody6Dynamics

    init = " NEW CHAIN   NCH NP ENERGY RSUM RGRAV TCR RMAXS NAMEC(1-NM) [ALL in NB UNIT]   4   1  -0.00031  1.7E-03  8.4E-04  5.3E-03  9.5E-03      9184      9185     11876     11877"
    step = " CHAIN: NSTEP T DTR R NAMEC    41    0.1953   -0.0000      9184  1.148E-03      9185  2.979E-04     11876  1.238E-04     11877"
    @test _tdyn_kind(init) === :chain_init
    ev, n = _tdyn_parse(:chain_init, [init, step], 1)
    @test n == 1 && ev isa nd.FewBodyEvent
    @test ev.n_members == 4 && ev.n_perturbers == 1
    @test ev.energy_nb == -0.00031 && ev.rsum_nb == 1.7e-3 && ev.rgrav_nb == 8.4e-4
    @test ev.tcr_nb == 5.3e-3 && ev.rmaxs_nb == 9.5e-3
    @test ev.names == (9184, 9185, 11876, 11877, 0, 0)
    # No time field: T0S + TIMEC of the following CHAIN: NSTEP line.
    @test ev.time_nb == 0.1953
    ev, _ = _tdyn_parse(:chain_init, [init], 1)
    @test isnan(ev.time_nb)

    exchange = " EXCHANGE    T NAM E0 E A0 A EB0 EB      0.33  9185 11876 11876  9184  0.920  0.790  3.8E-04  3.7E-04  -2.6E-05  -3.2E-05"
    binary = " CHAIN BINARY  NAM =     11876      9184  A =  3.6927E-04  E =  7.9015E-01  EB =  9.8288E-01  GB =  3.0391E-04  G4 =  0.0000E+00  EB1 = -0.0000E+00  E1 =  3.1799E-01  ET =  2.1756E-04  DB =  5.1741E-01"
    recoil = " RECOIL:    NAM E E1 PM/A RD A A1 RP GB IN   9185  0.7901  0.79  1.7   0.2  3.69E-04  2.10E-03  2.81E-03  2.93E-04   93.6"
    endchain = " END CHAIN NCH I1-4   3  1  2  3  1 N1-4 =      9184      9185     11876      9184 M1-4 =  8.05914E-05  6.11652E-05  1.47623E-04  8.05914E-05 RB =  2.95079E-04 R13 =  4.76219E-03 R24 =  1.00000E+10 DE =  5.37979E-13 TC =  3.47151E+00 #   712    5   3 GPERT =  4.89053E-04 DB =  3.30368E-02 EC =  2.33164E-04"
    namec = " CHTERM: NAMEC       9184      9185     11876     11877 M1-4   8.0591400E-05  6.1165200E-05  1.4762300E-04  3.3363647E-05  RIJ[NB] (1-2 1-3 2-3 2-4 3-4)   2.9507900E-04  4.7621900E-03  4.7998050E-03  4.7998050E-03  0.0000000E+00 M[*]   9.5959085E-01  7.2777547E-01  1.7401952E+00  3.8859840E-01"
    chterm = " CHTERM: TIME[NB] I3 I4 ICM1 NI3,4,CM  STEP(I3,I4,ICM1)[NB]     0.1031     12926     12926     24891     11876     11876     29184  2.64E-05  2.64E-05  3.75E-05"
    second =
        " CHTERM  SECOND BINARY(ICM2)[All in NB]: I3 I4 NI3,4 STEP STEPR R12 H GAMMA" *
        "     12001     12002      8001      8002" *
        "  1.53E-05  3.05E-05  2.10E-04 -1.25E+00  4.00E-03"
    lines = [exchange, binary, recoil, "", endchain, namec, "", _tdyn_ks_new(), chterm, second]
    @test [_tdyn_kind(s) for s in lines] == [
        :exchange,
        :chain_binary,
        :none,
        :none,
        :chain_end,
        :chterm_names,
        :none,
        :ks_new,
        :chterm,
        :chterm_second_binary,
    ]

    ev, _ = _tdyn_parse(:exchange, lines, 1)
    @test ev isa nd.ExchangeEvent && ev.kind === :exchange && ev.time_nb == 0.33
    @test (ev.name1_initial, ev.name2_initial) == (9185, 11876)
    @test (ev.name1_final, ev.name2_final) == (11876, 9184)
    @test ev.e0 == 0.92 && ev.e == 0.79 && ev.a0_nb == 3.8e-4 && ev.a_nb == 3.7e-4
    @test ev.eb0_nb == -2.6e-5 && ev.eb_nb == -3.2e-5 && !ev.overflow
    glued = replace(exchange, "      0.33  9185" => "      0.33123456")
    ev, _ = _tdyn_parse(:exchange, [glued], 1)
    @test ev.time_nb == 0.33 && ev.name1_initial == 123456 && ev.name2_initial == 11876
    over = replace(exchange, "  9185 11876 11876" => "****** 11876 11876")
    ev, _ = _tdyn_parse(:exchange, [over], 1)
    @test ev.overflow && ev.name1_initial == 0 && ev.name2_initial == 11876

    ev, _ = _tdyn_parse(:chain_binary, lines, 2)
    @test ev.kind === :chain_binary && ev.line == 2 && ev.time_nb == 0.1031
    @test ev.names == (11876, 9184, 0, 0, 0, 0)
    @test ev.semi_in_nb == 3.6927e-4 && ev.ecc_in == 0.79015
    @test ev.eb_over_energy == 0.98288 && ev.gpert == 3.0391e-4 && ev.g4 == 0.0
    @test ev.eb1_over_energy == 0.0 && ev.e1_over_eb0 == 0.31799
    @test ev.et_over_cm8 == 2.1756e-4 && ev.db == 0.51741

    ev, _ = _tdyn_parse(:chain_end, lines, 5)
    @test ev.kind === :chain_end && ev.time_nb == 0.1031 && ev.n_members == 3
    # NCH = 3: I4 = I1 repeated and R24 = 1e10 are stored absent.
    @test ev.indices == (1, 2, 3, 0)
    @test ev.names == (9184, 9185, 11876, 0, 0, 0)
    @test ev.masses_nb[1:3] == (8.05914e-5, 6.11652e-5, 1.47623e-4) && isnan(ev.masses_nb[4])
    @test ev.separations_nb[1:2] == (2.95079e-4, 4.76219e-3) && isnan(ev.separations_nb[3])
    @test ev.de_nb == 5.37979e-13 && ev.duration_crossings == 3.47151
    @test ev.n_steps == 712 && ev.n_reregularisations == 5 && ev.n_perturbers == 3
    @test ev.gpert == 4.89053e-4 && ev.db == 3.30368e-2 && ev.ec == 2.33164e-4

    ev, _ = _tdyn_parse(:chterm_names, lines, 6)
    @test ev.kind === :chterm_names && ev.time_nb == 0.1031
    # RIJ(3-4) = 0: three-body chain, the stale fourth member is dropped.
    @test ev.n_members == 3 && ev.names == (9184, 9185, 11876, 0, 0, 0)
    @test ev.masses_nb[1:3] == (8.05914e-5, 6.11652e-5, 1.47623e-4) && isnan(ev.masses_nb[4])
    @test ev.separations_nb[1:3] == (2.950790e-4, 4.762190e-3, 4.799805e-3)
    @test isnan(ev.separations_nb[4]) && isnan(ev.separations_nb[5])
    @test ev.masses_msun[1:3] == (0.95959085, 0.72777547, 1.7401952) && isnan(ev.masses_msun[4])

    ev, _ = _tdyn_parse(:chterm, lines, 9)
    @test ev.kind === :chterm && ev.time_nb == 0.1031
    @test ev.indices == (12926, 12926, 24891, 0)
    @test ev.names == (11876, 11876, 29184, 0, 0, 0)
    @test ev.steps_nb == (2.64e-5, 2.64e-5, 3.75e-5)

    ev, _ = _tdyn_parse(:chterm_second_binary, lines, 10)
    @test ev.kind === :chterm_second_binary && ev.time_nb == 0.1031
    @test ev.indices == (12001, 12002, 0, 0) && ev.names == (8001, 8002, 0, 0, 0, 0)
    @test ev.steps_nb[1:2] == (1.53e-5, 3.05e-5) && isnan(ev.steps_nb[3])
    @test ev.r12_nb == 2.1e-4 && ev.h_nb == -1.25 && ev.gamma == 4.0e-3

    # Collision termination: the time comes from the following CHEND: line.
    collision =
        " END CHAIN:   RIJ RCOLL QPERI GPERT " *
        "  1.2345678E-04  2.3456789E-03  3.4567890E-03" *
        "  1.0000000E-06  5.0000000E-07  2.5000000E-04"
    chend = " CHEND: TNB,TC[NB] K,M,R[*],R,V[NB-CH]=" * " 3.1250000000E-01" * "   1.00000E-02"
    @test _tdyn_kind(collision) === :chain_end_collision && _tdyn_kind(chend) === :none
    ev, _ = _tdyn_parse(:chain_end_collision, [collision, chend], 1)
    @test ev.time_nb == 0.3125
    @test ev.separations_nb[1:3] == (1.2345678e-4, 2.3456789e-3, 3.456789e-3)
    @test ev.rcoll_nb == 1.0e-6 && ev.qperi_nb == 5.0e-7 && ev.gpert == 2.5e-4

    # Binary–binary chain selection: time of the NEW CHAIN line IMPACT writes just before it,
    # else of the following one.
    bb = " CHAIN B-B    NAME(I1-I4)      7461      7462      9906      9907  K*(INCM)  1  K*(OCM)  0  ECC0 0.774  ECC1 2.965  SEMI0  6.517E-05  SEMI2  3.329E-04  RIN3  7.276E-04  SEMI1 -3.569E-04  PERIM  3.211E-04"
    @test _tdyn_kind(bb) === :chain_bb
    later = replace(_tdyn_impact(), "TIME 9.40626E-02" => "TIME 9.50000E-02")
    ev, _ = _tdyn_parse(:chain_bb, ["", _tdyn_impact(), bb, "", later], 3)
    @test ev.time_nb == 9.40626e-2
    ev, _ = _tdyn_parse(:chain_bb, [bb, "", _tdyn_impact()], 1)
    @test ev.time_nb == 9.40626e-2 && ev.n_members == 4
    @test ev.names == (7461, 7462, 9906, 9907, 0, 0)
    @test ev.kstar_values == (1, 0, nd.KSTAR_ABSENT, nd.KSTAR_ABSENT)
    @test ev.ecc_in == 0.774 && ev.ecc_out == 2.965
    @test ev.semi_in_nb == 6.517e-5 && ev.semi_second_nb == 3.329e-4
    @test ev.rij_nb == 7.276e-4 && ev.semi_out_nb == -3.569e-4 && ev.pmin_nb == 3.211e-4
end

@testset "Dynamics events — chain escapes" begin
    nd = Nbody6Dynamics
    single = " CHAIN ESCAPE: TT T TC =  1.0391360E-01  1.0390564E-01  9.0444133E-03 IESC NM RI VI =   3   11877  5.6049101E-03  6.7597194E-01 GPERT 2*M/R VF   4.3682444E-06  2.4141645E-01  4.2286278E+00"
    @test _tdyn_kind(single) === :chain_escape_single
    ev, n = _tdyn_parse(:chain_escape_single, [single], 1)
    @test n == 1 && ev.kind === :chain_escape_single
    @test ev.time_nb == 0.1039136 && ev.time_block_nb == 0.10390564
    @test ev.time_chain_nb == 9.0444133e-3
    @test ev.indices == (3, 0, 0, 0) && ev.names == (11877, 0, 0, 0, 0, 0)
    @test ev.r_escape_nb == 5.6049101e-3 && ev.rdot_escape_nb == 0.67597194
    @test ev.gpert == 4.3682444e-6 && ev.v2_escape_nb == 0.24141645
    @test ev.vinf_kms == 4.2286278 && isnan(ev.rb_nb)

    escape1 = " CHAIN ESCAPE: TT T TC =  2.5586819E+00  2.5586739E+00  4.0198436E-03 IESC JESC NM1/2 M1/2 RI VI =   2   3    5066    5067  2.7414361E-04  4.3689725E-05 GPERT 2*M/R RB VF =  3.3778760E-03  2.1406562E+00  4.1429941E-05  3.3924825E-01"
    continuation = " CHAIN ESCAPE: TT T TC =  2.2799820E-06  5.6496631E+00"
    @test _tdyn_kind(escape1) === :chain_escape_binary
    @test _tdyn_kind(continuation) === :none
    ev, n = _tdyn_parse(:chain_escape_binary, [escape1, continuation], 1)
    # Two physical lines, one record; values under their true meaning.
    @test n == 2 && ev.kind === :chain_escape_binary && !ev.overflow
    @test ev.time_nb == 2.5586819 && ev.time_block_nb == 2.5586739
    @test ev.time_chain_nb == 4.0198436e-3
    @test ev.indices == (2, 3, 0, 0) && ev.names == (5066, 5067, 0, 0, 0, 0)
    @test ev.masses_nb[1:2] == (2.7414361e-4, 4.3689725e-5) && isnan(ev.masses_nb[3])
    @test ev.r_escape_nb == 3.377876e-3 && ev.rdot_escape_nb == 2.1406562
    @test ev.gpert == 4.1429941e-5 && ev.v2_escape_nb == 0.33924825
    @test ev.rb_nb == 2.279982e-6 && ev.vinf_kms == 5.6496631
    # Without its continuation the record is read from the first line alone.
    ev, n = _tdyn_parse(:chain_escape_binary, [escape1], 1)
    @test n == 1 && isnan(ev.rb_nb) && isnan(ev.vinf_kms) && ev.gpert == 4.1429941e-5
end

@testset "Dynamics events — hierarchy ends and higher-order mergers" begin
    nd = Nbody6Dynamics
    endmerger = " END MERGER:  TIME 2.68621E+00 IC,JC,IP=      8640      8641      4870 N IC,JC,ICM=      5938      5939     26334 KW IC,JC,ICM=   1   0   0 M1,2[NB]  9.526E-05  5.661E-05 JX1,JX2,JP=      8640      8640      4870 N JX1,JX2,JP=      5938      5938     26334 KW JX1,JX2,JCM=   1   1   0 M3,4[NB]  9.526E-05  9.526E-05  R12  5.795E-04 H -2.724E-02 R34  5.795E-04 H2 -2.724E-02 M1,2[*]  1.115E+00  6.541E-01 RAD1,2[*]  9.037E-01  6.112E-01  M3,4[*]  1.115E+00  1.115E+00 RAD3,4[*]  9.037E-01  9.037E-01 SEPI/J  3.183E+04  3.183E+04"
    @test _tdyn_kind(endmerger) === :hierarchy_end
    ev, n = _tdyn_parse(:hierarchy_end, [endmerger], 1)
    @test n == 1 && ev isa nd.HierarchyEvent && ev.kind === :hierarchy_end
    @test ev.time_nb == 2.68621
    @test ev.indices == (8640, 8641, 0, 0)
    @test ev.names == (5938, 5939, 0, 0, 0, 0)
    @test ev.kstar_values == (1, 0, nd.KSTAR_ABSENT, nd.KSTAR_ABSENT)
    @test ev.masses_nb[1:2] == (9.526e-5, 5.661e-5) && all(isnan, ev.masses_nb[3:4])
    @test ev.masses_msun[1:2] == (1.115, 0.6541) && all(isnan, ev.masses_msun[3:4])
    @test ev.radii_rsun[1:2] == (0.9037, 0.6112) && isnan(ev.radii_rsun[3])
    # Stale slot data (IP, NAME/KSTAR of ICM, R12, H, SEPI/J) is not stored.
    @test isnan(ev.r12_nb) && isnan(ev.h_nb) && isnan(ev.r12_rsun) && !ev.overflow

    t9, t10 = lpad("12.50", 9), lpad("12.50", 10)
    zm = lpad("3.25", 6)
    names3 = lpad("10001", 6) * lpad("10002", 6) * lpad("20003", 6)
    names4 = names3 * lpad("10004", 6)
    e3 = lpad("2.50E-04", 10) * lpad("1.10E-03", 10) * lpad("3.00E-03", 10)
    e6 = lpad("4.00E-05", 10) * lpad("5.00E-05", 10) * e3 * lpad("6.00E-03", 10)
    labels_new5 = "   T MT NM1 NM2 NM3 RI E0 EB0 A0 PC G0"
    labels_end5 = "   T MT NM1 NM2 NM3 RI E1 EB1 A1 PM G1 "
    ri, ecc = lpad("0.45", 6), lpad("0.30", 6)
    new5 = t10 * zm * names3 * ri * ecc * lpad("1.20E+01", 10) * e3
    end5 = t9 * zm * names3 * ri * ecc * lpad("12.0", 6) * e3
    quartet_new = " NEW QUARTET    T MT NM1 NMG NM3 E0 A0 PC G0 " * t9 * zm * names3 * ecc * e3
    quintup2_new = " NEW QUINTUP2    T MT NM1 NMG NM3 E0 A0 PC G0" * t10 * zm * names3 * ecc * e3
    hitrip_new = " NEW HITRIP    T MT NM E0 AI AJ RI RJ PC GX " * t9 * zm * names4 * ecc * e6
    quartet_end = " END QUARTET    T MT NM1 NMG NM3 E1 A1 R1 G1 " * t9 * zm * names3 * ecc * e3
    quintup2_end = " END QUINTUP2    T MT NM1 NMG NM3 E1 A1 R1 G1" * t10 * zm * names3 * ecc * e3
    hitrip_end = " END HITRIP    T MT NM E1 AI AJ R1 A1 PM G1 " * t9 * zm * names4 * ecc * e6

    cases = (
        (" NEW QUINTUPLET" * labels_new5 * new5, :hierarchy_new_quintuplet),
        (" NEW SEXTUPLET " * labels_new5 * new5, :hierarchy_new_sextuplet),
        (quartet_new, :hierarchy_new_quartet),
        (quintup2_new, :hierarchy_new_quintup2),
        (hitrip_new, :hierarchy_new_hitrip),
        (" END QUINTUPLET" * labels_end5 * end5, :hierarchy_end_quintuplet),
        (" END SEXTUPLET " * labels_end5 * end5, :hierarchy_end_sextuplet),
        (quartet_end, :hierarchy_end_quartet),
        (quintup2_end, :hierarchy_end_quintup2),
        (hitrip_end, :hierarchy_end_hitrip),
    )
    for (line, kind) in cases
        @test _tdyn_kind(line) === kind
        ev, n = _tdyn_parse(kind, [line], 1)
        @test n == 1 && ev isa nd.HierarchyEvent && ev.kind === kind
        @test ev.time_nb == 12.5 && ev.masses_msun[4] == 3.25 && !ev.overflow
        @test ev.names[1:3] == (10001, 10002, 20003)
    end

    ev, _ = _tdyn_parse(:hierarchy_new_quintuplet, [cases[1][1]], 1)
    @test ev.r_cm_nb == 0.45 && ev.ecc_in == 0.3 && ev.eb_kinetic == 12.0
    @test ev.semi_in_nb == 2.5e-4 && ev.pcrit_nb == 1.1e-3 && ev.gamma == 3.0e-3
    ev, _ = _tdyn_parse(:hierarchy_new_quartet, [cases[3][1]], 1)
    @test ev.ecc_in == 0.3 && ev.semi_in_nb == 2.5e-4
    @test ev.pcrit_nb == 1.1e-3 && ev.gamma == 3.0e-3
    ev, _ = _tdyn_parse(:hierarchy_new_hitrip, [cases[5][1]], 1)
    @test ev.names[4] == 10004 && ev.ecc_in == 0.3
    @test ev.semi_inner1_nb == 4.0e-5 && ev.semi_inner2_nb == 5.0e-5
    @test ev.r12_nb == 2.5e-4 && ev.r34_nb == 1.1e-3
    @test ev.pcrit_nb == 3.0e-3 && ev.gamma == 6.0e-3
    ev, _ = _tdyn_parse(:hierarchy_end_sextuplet, [cases[7][1]], 1)
    @test ev.r_cm_nb == 0.45 && ev.ecc_out == 0.3 && ev.eb_kinetic == 12.0
    @test ev.semi_out_nb == 2.5e-4 && ev.pmin_nb == 1.1e-3 && ev.gamma_merged == 3.0e-3
    ev, _ = _tdyn_parse(:hierarchy_end_quintup2, [cases[9][1]], 1)
    @test ev.ecc_out == 0.3 && ev.semi_out_nb == 2.5e-4 && ev.r12_nb == 1.1e-3
    @test ev.gamma_merged == 3.0e-3
    ev, _ = _tdyn_parse(:hierarchy_end_hitrip, [cases[10][1]], 1)
    @test ev.semi_inner1_nb == 4.0e-5 && ev.r12_nb == 2.5e-4 && ev.semi_out_nb == 1.1e-3
    @test ev.pmin_nb == 3.0e-3 && ev.gamma_merged == 6.0e-3

    # I6 names overflow at 10⁶.
    over = replace(cases[8][1], lpad("10002", 6) => "******")
    ev, _ = _tdyn_parse(:hierarchy_end_quartet, [over], 1)
    @test ev.overflow && ev.names[1:3] == (10001, 0, 20003)

    # END OUTER MERGER has no time field: T of the following END MERGE2 line.
    outer =
        " END OUTER MERGER" *
        lpad("12", 5) *
        "  H =" *
        lpad("-12.5", 7) *
        "  M =" *
        lpad("0.0012", 7) *
        lpad("0.0034", 7) *
        "  E1 =" *
        lpad("0.250", 6) *
        "  EB2 =" *
        lpad("-0.125", 6) *
        "  RB2 =" *
        lpad("1.0E-04", 8) *
        "  G2 =" *
        lpad("2.0E-03", 8) *
        "  DP =" *
        lpad("-3.0E-09", 8)
    merge2 =
        " END MERGE2" *
        lpad("4", 3) *
        "  T =" *
        lpad("12.75", 8) *
        "  M =" *
        lpad("0.0056", 7) *
        lpad("0.0078", 7) *
        "  R1 =" *
        lpad("5.0E-05", 8) *
        "  A1 =" *
        lpad("6.0E-04", 8) *
        "  EB =" *
        lpad("-0.500", 6) *
        "  E1 =" *
        lpad("0.125", 6) *
        "  GB =" *
        lpad("7.0E-03", 8) *
        "  G =" *
        lpad("0.010", 6) *
        "  NB =" *
        lpad("25", 3)
    @test _tdyn_kind(outer) === :hierarchy_end_outer
    @test _tdyn_kind(merge2) === :hierarchy_end_merge2
    ev, _ = _tdyn_parse(:hierarchy_end_outer, [outer, merge2], 1)
    @test ev.time_nb == 12.75 && ev.indices == (0, 0, 12, 0) && ev.h_nb == -12.5
    @test ev.masses_nb[1:2] == (0.0012, 0.0034) && ev.energy_ratio == 0.25
    @test ev.eb_fraction == -0.125 && ev.r12_nb == 1.0e-4 && ev.gamma == 2.0e-3
    @test ev.dpot_nb == -3.0e-9
    ev, _ = _tdyn_parse(:hierarchy_end_merge2, [outer, merge2], 2)
    @test ev.time_nb == 12.75 && ev.merger_index == 4 && ev.n_neighbours == 25
    @test ev.masses_nb[1:2] == (0.0056, 0.0078) && ev.r12_nb == 5.0e-5
    @test ev.semi_out_nb == 6.0e-4 && ev.eb_fraction == -0.5 && ev.energy_ratio == 0.125
    @test ev.gamma == 7.0e-3 && ev.gamma_merged == 0.01

    # HI MERGE follows NEW MERGE2 and takes its time.
    himerge =
        " HI MERGE    NAM E E1 A A1 PM PC " *
        join(lpad(string(k), 7) for k in (101, 102, 103, 104, 105, 106)) *
        lpad("0.500", 7) *
        lpad("0.900", 7) *
        lpad("1.000E-04", 11) *
        lpad("2.000E-03", 11) *
        lpad("3.000E-04", 11) *
        lpad("4.000E-04", 11)
    newmerge2 = replace(_tdyn_impact(), " NEW CHAIN   TIME" => " NEW MERGE2  TIME")
    @test _tdyn_kind(himerge) === :hierarchy_hi_merge
    ev, _ = _tdyn_parse(:hierarchy_hi_merge, [newmerge2, himerge], 2)
    @test ev.time_nb == 9.40626e-2 && ev.line == 2
    @test ev.names == (101, 102, 103, 104, 105, 106)
    @test ev.ecc_in == 0.5 && ev.ecc_out == 0.9 && ev.semi_in_nb == 1.0e-4
    @test ev.semi_out_nb == 2.0e-3 && ev.pmin_nb == 3.0e-4 && ev.pcrit_nb == 4.0e-4

    @test_throws ArgumentError _tdyn_parse(
        :hierarchy_end,
        [replace(endmerger, " IC,JC,IP=" => " IC,JC,IQ=")],
        1,
    )
end

@testset "Dynamics events — termination diagnostics" begin
    nd = Nbody6Dynamics
    deferred =
        repeat(" ", 9) *
        "END RUN DEFERRED   TTOT =" *
        lpad("1.25000000E+02", 16) *
        "  NSUB =" *
        lpad("2", 3)
    @test _tdyn_kind(deferred) === :end_run_deferred
    ev, n = _tdyn_parse(:end_run_deferred, ["", deferred], 2)
    @test n == 1 && ev isa nd.TerminationEvent && ev.kind === :end_run_deferred
    @test ev.line == 2 && ev.time_nb == 125.0 && ev.n_subsystems == 2 && !ev.overflow

    tiny =
        " ENFORCED CHAIN TINY STEPS:  NSTEP NCH STEP R" *
        lpad("1000000", 10) *
        lpad("3", 4) *
        lpad("1.0E-15", 9) *
        lpad("2.0E-04", 9) *
        lpad("3.0E-04", 9)
    @test _tdyn_kind(tiny) === :chain_tiny_steps
    ev, _ = _tdyn_parse(:chain_tiny_steps, [tiny], 1)
    @test ev.kind === :chain_tiny_steps && isnan(ev.time_nb)
    @test ev.n_steps == 1000000 && ev.n_members == 3 && ev.step_nb == 1.0e-15
    @test ev.separations_nb[1:2] == (2.0e-4, 3.0e-4) && all(isnan, ev.separations_nb[3:8])

    @test_throws ArgumentError _tdyn_parse(:end_run_deferred, [lstrip(deferred)], 1)
end
