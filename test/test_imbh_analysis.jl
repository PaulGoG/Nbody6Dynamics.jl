# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# IMBH diagnostics: the mbh.45 reader, the trajectory about the stellar
# density centre on synthetic Gaussian clouds (StableRNGs), the wandering
# moments and the Brownian predictions, the bound-companion search on a
# circular test orbit, the event list from synthetic engine records, the IC
# record lookup and the CSV writers.

using StableRNGs: StableRNG

@testset "IMBH diagnostics" begin
    nd = Nbody6Dynamics
    G = nd._G_PC_KMS2_MSUN
    M_tot = 1.0e5                       # mass unit [M☉]
    v_unit = sqrt(G * M_tot / 1.0)      # G = 1 with RBAR = 1 pc
    vs = Float64(Float32(v_unit))       # the header holds VSTAR in single precision
    function mksnap(t, names, mass, pos, vel; npairs = 0)
        params = zeros(Float32, 20)
        params[1] = t
        params[3] = 1.0f0
        params[4] = Float32(M_tot)
        params[11] = 1.0f0
        params[12] = Float32(v_unit)
        return nd.Snapshot(
            nd.SnapshotHeader(Int32(length(mass)), Int32(1), Int32(1), Int32(20), params),
            Int32.(collect(names)),
            Float32.(mass),
            Float32.(pos),
            Float32.(vel),
            Float32[],
            Float32[],
            npairs,
        )
    end

    @testset "mbh.45 reader" begin
        mktempdir() do dir
            record(which, t, ibh, n) =
                " $which $t $ibh 0.1 -0.2 0.3 1.0 2.0 3.0 12 " * join(fill("0.5", n - 9), " ")
            path = joinpath(dir, "mbh.45")
            write(
                path,
                "     STAT         Time[Myr]                 IBH         X(1)[PC]\n" *
                record("BINARY", 1.5, 1000, 20) *
                "\n\n" *
                record("SINGLE", 2.5, 999, 20) *
                "\n",
            )
            rows = nd.read_mbh45(path)
            @test length(rows) == 2
            @test rows[1].state === :binary && rows[2].state === :single
            @test rows[1].time_myr == 1.5 && rows[2].time_myr == 2.5
            @test rows[1].index == 1000 && rows[2].index == 999
            @test rows[1].pos_pc == (0.1, -0.2, 0.3)
            @test rows[1].vel_kms == (1.0, 2.0, 3.0)
            @test rows[1].n_neighbours == 12
            @test rows[1].vsigma_km2s2 == (0.5, 0.5, 0.5)
            @test rows[1].name == 0 && rows[2].name == 0
            # fork v2026.07.11 on: the body name as a 21st number, integer or real
            named = joinpath(dir, "named.45")
            write(
                named,
                record("SINGLE", 3.5, 998, 20) *
                " 8193\n" *
                record("BINARY", 4.5, 997, 20) *
                " 8.193000000000000E+03\n",
            )
            rows_named = nd.read_mbh45(named)
            @test length(rows_named) == 2
            @test all(r -> r.name == 8193, rows_named)
            @test rows_named[1].index == 998 && rows_named[1].vsigma_km2s2 == (0.5, 0.5, 0.5)
            bad = joinpath(dir, "bad.45")
            write(bad, record("SINGLE", 1.0, 5, 20) * "\n" * record("SINGLE", 2.0, 5, 19) * "\n")
            @test_throws ArgumentError nd.read_mbh45(bad)
            @test_throws ErrorException nd.read_mbh45(joinpath(dir, "absent.45"))
        end
    end

    # Gaussian star cloud of σ = 0.5 pc with the IMBH (1 % of the mass) offset
    N = 3000
    rng = StableRNG(3)
    pos_s = 0.5 .* randn(rng, 3, N)
    vel_s = 0.1 .* randn(rng, 3, N)
    m_s = fill(0.99 / N, N)
    bh = N + 1
    pos_all = hcat(pos_s, [0.3, 0.0, 0.0])
    vel_all = hcat(vel_s, [0.0, 0.2, 0.0])
    mass_all = vcat(m_s, 0.01)
    snaps = [
        mksnap(0.0, 1:(N + 1), mass_all, pos_all, vel_all),
        mksnap(1.0, 1:(N + 1), mass_all, pos_all, vel_all),
        mksnap(2.0, 1:N, m_s, pos_s, vel_s),
    ]

    @testset "Trajectory" begin
        traj = nd.imbh_trajectory(snaps, bh)
        @test length(traj) == 2
        s = traj[1]
        @test s.index == N + 1
        @test s.time_nb == 0.0 && traj[2].time_myr == 1.0
        @test isapprox(s.mass_msun, 0.01 * M_tot; rtol = 1e-6)
        @test isapprox(s.r_pc, 0.3; atol = 0.1)   # density-centre noise of a 3000-star Gaussian cloud
        @test isapprox(s.v_kms[2], 0.2 * vs; atol = 0.02 * vs)
        @test !s.in_pair && s.partner_name == 0
        @test isapprox(s.r_half_pc, 1.538 * 0.5; rtol = 0.2)
        @test isnan(s.r_core_pc)
        @test isapprox(s.r_cm_pc, 0.3; atol = 0.05)
        # KS component: the IMBH in column 1 of a snapshot with one pair
        ks = mksnap(
            0.0,
            vcat(bh, 1:N),
            vcat(0.01, m_s),
            hcat([0.3, 0.0, 0.0], pos_s),
            hcat([0.0, 0.2, 0.0], vel_s);
            npairs = 1,
        )
        s_ks = only(nd.imbh_trajectory([ks], bh))
        @test s_ks.index == 1
        @test s_ks.in_pair
        @test s_ks.partner_name == 1
        @test_throws ArgumentError nd.imbh_trajectory(nd.Snapshot[], bh)
        # Effective mass of the stage: all stars without a core radius
        @test isapprox(nd._imbh_core_effective_mass(snaps, bh, traj), 0.99 / N * M_tot; rtol = 1e-5)
    end

    @testset "Wandering and Brownian predictions" begin
        rng5 = StableRNG(5)
        n_w = 2000
        states = map(1:n_w) do k
            x = Tuple(0.2 .* randn(rng5, 3))
            nd.IMBHState(;
                time_nb = k - 1.0,
                time_myr = k - 1.0,
                index = 1,
                mass_msun = 1.0,
                r_pc = sqrt(sum(abs2, x)),
                x_pc = x,
                v_kms = (0.0, 0.0, 0.0),
                r_cm_pc = NaN,
                r_core_pc = NaN,
                r_half_pc = NaN,
                in_pair = false,
                partner_name = 0,
            )
        end
        w = only(nd.imbh_wandering(states))
        @test w.frame === :density_centre
        @test w.n == n_w
        # centre-of-mass frame: components from x_cm_pc (half the density-centre offsets here)
        states_cm = [
            nd.IMBHState(;
                (f => getfield(st, f) for f in fieldnames(nd.IMBHState))...,
                x_cm_pc = st.x_pc ./ 2,
            ) for st in states
        ]
        w_cm = only(nd.imbh_wandering(states_cm; frame = :centre_of_mass))
        @test w_cm.frame === :centre_of_mass
        @test isapprox(w_cm.x2_mean_pc2, w.x2_mean_pc2 / 4; rtol = 1e-12)
        @test_throws ArgumentError nd.imbh_wandering(states; frame = :barycentre)
        @test w.t_start_myr == 0.0 && w.t_end_myr == n_w - 1.0
        @test isapprox(w.x2_mean_pc2, 0.04; rtol = 0.1)
        @test isapprox(w.r2_mean_pc2, 0.12; rtol = 0.1)
        @test w.r_rms_pc ≈ sqrt(w.r2_mean_pc2)
        @test isapprox(w.kurtosis, 3.0; atol = 0.3)
        @test w.v2_mean_km2s2 == 0.0
        halves = nd.imbh_wandering(states; window_myr = (n_w - 1) / 2)
        @test length(halves) == 2
        @test halves[1].n + halves[2].n ≥ n_w - 2
        @test isempty(nd.imbh_wandering(nd.IMBHState[]))
        @test isapprox(nd.brownian_rms_plummer(3π / 16, 8.0e-5, 0.01), 0.02483; rtol = 1e-3)
        @test isapprox(nd.brownian_rms_core(1.0, 1.0, 1.0, 1.0), 7.450; atol = 0.005)
        # Σm²/Σm = 6/4
        @test nd.effective_mass([1, 1, 2]) == 1.5
        @test nd.effective_mass(fill(0.3, 10)) ≈ 0.3
        @test_throws ArgumentError nd.effective_mass(Float64[])
        @test_throws ArgumentError nd.brownian_rms_core(1.0, 0.0, 1.0, 1.0)
    end

    @testset "Companions" begin
        # IMBH at rest at the origin, a test star on a circular orbit of 10⁻³ pc
        m_c = 1 / N
        v_circ = sqrt((0.01 + m_c) / 1.0e-3)
        comp_snap(v) = mksnap(
            0.0,
            1:(N + 2),
            vcat(m_s, 0.01, m_c),
            hcat(pos_s, [0.0, 0.0, 0.0], [1.0e-3, 0.0, 0.0]),
            hcat(vel_s, [0.0, 0.0, 0.0], [0.0, v, 0.0]),
        )
        c = only(nd.imbh_companions([comp_snap(v_circ)], bh))
        @test c.companion_name == N + 2
        @test isapprox(c.companion_mass_msun, m_c * M_tot; rtol = 1e-6)
        @test isapprox(c.semi_major_axis_au, 1.0e-3 * 206264.806; rtol = 1e-5)
        @test isapprox(c.separation_pc, 1.0e-3; rtol = 1e-6)
        @test c.eccentricity < 1e-5
        @test c.specific_energy_km2s2 < 0
        @test !c.in_pair
        @test isfinite(c.t_gw_myr) && c.t_gw_myr > 0
        # Twice the circular speed: unbound. The cloud lies inside the IMBH's
        # radius of influence, so its nearest stars are two-body bound; k = 1
        # restricts the candidates to the test star.
        u = only(nd.imbh_companions([comp_snap(2v_circ)], bh; k = 1))
        @test u.companion_name == 0
        @test isnan(u.semi_major_axis_au) && isnan(u.eccentricity)
        @test isnan(u.specific_energy_km2s2) && isnan(u.t_gw_myr)
        @test isempty(nd.imbh_companions([snaps[3]], bh))
        @test_throws ArgumentError nd.imbh_companions(snaps, bh; k = -1)

        bound(t, n) = nd.IMBHCompanion(t, n, 1.0, 1.0e-3, 200.0, 0.1, -1.0, false, 1.0)
        unbound(t) = nd.IMBHCompanion(t, 0, NaN, NaN, NaN, NaN, NaN, false, NaN)
        s = nd.imbh_companion_summary([bound(1.0, 7), bound(2.0, 8), unbound(3.0)])
        @test s.n_snapshots == 3
        @test s.fraction_bound ≈ 2 / 3
        @test s.n_exchanges == 1
        @test s.longest_name == 7 && s.longest_span_myr == 0.0
        s2 = nd.imbh_companion_summary([bound(1.0, 7), unbound(2.0), bound(3.0, 8), bound(5.0, 8)])
        @test s2.longest_name == 8 && s2.longest_span_myr == 2.0
        @test s2.n_exchanges == 1
        s0 = nd.imbh_companion_summary([unbound(1.0)])
        @test s0.fraction_bound == 0.0 && s0.longest_name == 0 && isnan(s0.longest_span_myr)
    end

    @testset "Events" begin
        imbh = 5001
        coll(t, n1, n2, k2) = nd.CollisionRecord(;
            time_nb = t,
            name1 = n1,
            name2 = n2,
            kstar1 = 14,
            kstar2 = k2,
            kstar_product = 14,
            mass1_msun = 1000.0,
            mass2_msun = 0.5,
            mass_product_msun = 1000.2,
            mass_lost_msun = 0.3,
            radius1_rsun = 0.0,
            radius2_rsun = 0.5,
            distance_over_core_radius = 0.1,
            separation_rsun = 1.0,
            eccentricity = 0.9,
            period_d = 1.0,
        )
        coal = nd.CoalescenceRecord(;
            channel = :binary,
            iqcoll = 0,
            time_nb = 3.0,
            name1 = 60,
            name2 = imbh,
            kstar1 = 13,
            kstar2 = 14,
            kstar_product = 14,
            mass1_bse_msun = 1.4,
            mass2_bse_msun = 1000.0,
            separation_nb = 1e-8,
            eccentricity = 0.0,
            semi_major_axis_nb = 1e-8,
            binding_energy_nb = -1.0,
            tidal_potential_correction_nb = 0.0,
            period_d = 1e-4,
            mass1_msun = 1.4,
            mass2_msun = 1000.0,
            mass_product_msun = 1001.3,
            mass_lost_msun = 0.1,
            radius_survivor_rsun = 1e-5,
            collision_distance_rsun = 1e-4,
            velocity_infinity_kms = 0.0,
            distance_density_centre_nb = 0.1,
            speed_nb = 1.0,
        )
        recoil = nd.GWEvent(;
            line = 1,
            interval = 1,
            time_nb = 4.0,
            kind = :gw_recoil,
            time_myr = 8.0,
            name1 = imbh,
            name2 = 77,
            kstar1 = 14,
            kstar2 = 14,
            m1_msun = 1000.0,
            m2_msun = 10.0,
        )
        other_gw =
            nd.GWEvent(; line = 2, interval = 1, time_nb = 4.5, kind = :gw_ksint, name1 = imbh)
        escape = nd.EscaperRecord(10.0, 1001.3, 1.0, 50.0, 14, imbh, 0.0, 0.0)
        collisions = [coll(2.0, imbh, 43, 13), coll(1.0, imbh, 42, 1), coll(0.5, 3, 4, 1)]
        kw = (
            collisions = collisions,
            coalescences = [coal],
            escapers = [escape],
            gw_events = [recoil, other_gw],
        )
        ev = nd.imbh_events(imbh; kw..., tstar_myr = 2.0)
        @test length(ev) == 5
        @test [r.kind for r in ev] == [:tde, :collision, :coalescence, :recoil, :escape]
        @test [r.time_myr for r in ev] == [2.0, 4.0, 6.0, 8.0, 10.0]
        @test ev[1].partner_name == 42 && ev[1].partner_kstar == 1 && ev[1].partner_mass_msun == 0.5
        @test ev[2].partner_name == 43 && ev[2].partner_kstar == 13
        @test ev[3].partner_name == 60 &&
              ev[3].partner_kstar == 13 &&
              ev[3].partner_mass_msun == 1.4
        @test ev[4].partner_name == 77 && ev[4].partner_mass_msun == 10.0
        @test ev[5].partner_name == 0 && ev[5].partner_kstar == nd.KSTAR_ABSENT
        @test isnan(ev[5].partner_mass_msun)
        # Without T* the N-body-timed records have no time and sort last
        ev_nan = nd.imbh_events(imbh; kw...)
        @test [r.kind for r in ev_nan][1:2] == [:recoil, :escape]
        @test count(r -> isnan(r.time_myr), ev_nan) == 3
        @test isempty(nd.imbh_events(imbh))
    end

    @testset "IC record lookup" begin
        mktempdir() do dir
            path = joinpath(dir, "merger_ic.toml")
            @test nd._imbh_name(path) == 0
            write(path, "[meta]\nnbin0 = 0\n")
            @test nd._imbh_name(path) == 0
            write(path, "[meta]\nnbin0 = 0\n\n[imbh]\nname = 3001\nmass_msun = 1000.0\n")
            @test nd._imbh_name(path) == 3001
        end
    end

    @testset "Writers" begin
        mktempdir() do dir
            traj = nd.imbh_trajectory(snaps, bh)
            wand = nd.imbh_wandering(traj)
            comp = nd.imbh_companions(snaps, bh)
            ev = nd.imbh_events(
                bh;
                escapers = [nd.EscaperRecord(3.0, 1000.0, 1.0, 50.0, 14, bh, 0.0, 0.0)],
            )
            files = [
                (
                    nd.write_imbh_trajectory(joinpath(dir, "imbh_trajectory.csv"), traj; name = bh),
                    "time_nb,",
                    2,
                ),
                (
                    nd.write_imbh_wandering(joinpath(dir, "imbh_wandering.csv"), wand; name = bh),
                    "frame,t_start_myr,",
                    1,
                ),
                (
                    nd.write_imbh_companions(joinpath(dir, "imbh_companions.csv"), comp; name = bh),
                    "time_myr,",
                    2,
                ),
                (nd.write_imbh_events(joinpath(dir, "imbh_events.csv"), ev; name = bh), "kind,", 1),
            ]
            for (path, head, n_rows) in files
                @test isfile(path)
                lines = readlines(path)
                @test startswith(lines[1], "# imbh_name=$bh")
                @test startswith(lines[2], head)
                @test length(lines) == 2 + n_rows
                @test length(split(lines[3], ",")) == length(split(lines[2], ","))
            end
            nd.write_imbh_events(joinpath(dir, "imbh_events.csv"), ev; name = bh)
            @test isfile(joinpath(dir, "imbh_events#1.csv"))
        end
    end
end
