# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# IMBH diagnostics: the mbh.45 reader, the trajectory about the stellar
# density centre on synthetic Gaussian clouds (StableRNGs), the wandering
# moments and the Brownian predictions, the bound-companion search on a
# circular test orbit and its perturbation criterion, the elements of IMBH
# pairs on a circular binary, the hardening estimator on a synthetic pair
# series, the event list from synthetic engine records, the IC record lookup
# and the CSV writers.

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
        # Core radius of the stars from their own densities: the IMBH is no neighbour
        ρ_s = nd.local_density(Float64.(Float32.(pos_s)), Float64.(Float32.(m_s)))
        @test isapprox(
            s.r_core_pc,
            nd.core_radius(Float64.(Float32.(pos_s)), ρ_s).r_core;
            rtol = 1e-10,
        )
        @test isapprox(s.r_cm_pc, 0.3; atol = 0.05)
        # A tight cusp of 30 stars about the IMBH pulls the core radius to the cusp; without the stars inside
        # the influence radius (stellar mass 2 M• about the IMBH) it is the cloud's again
        cusp = [0.3, 0.0, 0.0] .+ 0.003 .* randn(rng, 3, 30)
        pos_c = hcat(pos_s, cusp, [0.3, 0.0, 0.0])
        snap_c = mksnap(
            0.0,
            1:(N + 31),
            vcat(m_s, fill(0.99 / N, 30), 0.01),
            pos_c,
            hcat(vel_s, zeros(3, 31)),
        )
        s_c = only(nd.imbh_trajectory([snap_c], N + 31))
        @test s_c.r_core_pc < 0.5 * s.r_core_pc
        @test isapprox(s_c.r_core_out_pc, s.r_core_pc; rtol = 0.15)
        @test isfinite(s.r_core_out_pc)
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
        # Effective mass of the stage over the stars inside the core radius (equal masses)
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

        bound(t, n) = nd.IMBHCompanion(t, n, 1.0, 1.0e-3, 200.0, 0.1, -1.0, false, 1.0, 0.0)
        unbound(t) = nd.IMBHCompanion(t, 0, NaN, NaN, NaN, NaN, NaN, false, NaN, NaN)
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

    @testset "Companion perturbation criterion" begin
        # The circular test star at d = 10⁻³ pc and a perturber of mass 1/N at
        # rest at 10⁻² pc, the IMBH at rest at the origin
        m_c = 1 / N
        d_c = 1.0e-3
        v_circ = sqrt((0.01 + m_c) / d_c)
        x_p = 1.0e-2
        function pert_snap(with_perturber)
            n_p = with_perturber ? 1 : 0
            return mksnap(
                0.0,
                1:(N + 2 + n_p),
                vcat(m_s, 0.01, m_c, fill(m_c, n_p)),
                hcat(pos_s, [0.0, 0.0, 0.0], [d_c, 0.0, 0.0], repeat([x_p, 0.0, 0.0], 1, n_p)),
                hcat(vel_s, [0.0, 0.0, 0.0], [0.0, v_circ, 0.0], zeros(3, n_p)),
            )
        end
        # Expected γ from the single-precision masses and positions of the snapshot
        f32(x) = Float64(Float32(x))
        M_b, m_j, m_p, x_j = f32(0.01), f32(m_c), f32(m_c), f32(d_c)
        x_cm = m_j * x_j / (M_b + m_j)
        γ_exp = 2 * m_p * x_j^3 / ((M_b + m_j) * (f32(x_p) - x_cm)^3)
        c = only(nd.imbh_companions([pert_snap(true)], bh))
        @test c.companion_name == N + 2
        @test isapprox(c.gamma, γ_exp; rtol = 1e-6)
        r = only(nd.imbh_companions([pert_snap(true)], bh; gamma_max = 1.0e-5))
        @test r.companion_name == 0
        @test isnan(r.semi_major_axis_au)
        @test r.gamma == c.gamma
        # Without the perturber the cloud stars are far
        c0 = only(nd.imbh_companions([pert_snap(false)], bh))
        @test c0.companion_name == N + 2
        @test c0.gamma < 1.0e-5
        s = nd.imbh_companion_summary([c, r])
        @test s.fraction_bound == 0.5
    end

    # Two IMBHs (names N + 1, N + 2) on a circular orbit of separation d about
    # their common centre of mass, inside the cloud; an optional third at 1 pc
    d_pair = 0.05
    m_pair = 0.01
    v_pair = sqrt(2m_pair / d_pair)     # relative circular speed, G = 1
    function pair_snap(t, v_rel; third = false)
        names = third ? (1:(N + 3)) : (1:(N + 2))
        mass = vcat(m_s, m_pair, m_pair, third ? [m_pair] : Float64[])
        pos = hcat(
            pos_s,
            [-d_pair / 2, 0.0, 0.0],
            [d_pair / 2, 0.0, 0.0],
            third ? reshape([1.0, 0.0, 0.0], 3, 1) : zeros(3, 0),
        )
        vel = hcat(
            vel_s,
            [0.0, -v_rel / 2, 0.0],
            [0.0, v_rel / 2, 0.0],
            third ? zeros(3, 1) : zeros(3, 0),
        )
        return mksnap(t, names, mass, pos, vel)
    end

    @testset "IMBH pairs" begin
        pair_names = [N + 1, N + 2]
        rows = nd.imbh_pairs([pair_snap(t, v_pair) for t in 0.0:2.0], pair_names)
        @test length(rows) == 3
        @test all(r -> r.name1 == N + 1 && r.name2 == N + 2, rows)
        @test [r.time_myr for r in rows] == [0.0, 1.0, 2.0]
        r = rows[1]
        @test r.bound && !r.in_pair
        @test r.specific_energy_km2s2 < 0
        @test isapprox(r.separation_pc, d_pair; rtol = 1e-6)
        @test isapprox(r.semi_major_axis_au, d_pair * 1.0 * 206264.806; rtol = 1e-5)
        @test r.eccentricity < 1e-5
        M12 = 2 * Float64(Float32(m_pair)) * M_tot
        a_pc = r.semi_major_axis_au / 206264.806
        P = 2π * sqrt(a_pc^3 / (G * M12)) * (3.0856775814913673e13 / 3.15576e13)
        @test isapprox(r.period_myr, P; rtol = 1e-9)
        @test isfinite(r.t_gw_myr) && r.t_gw_myr > 0
        # Influence radius: the stars nearest the pair's centre of mass (the origin) holding 2 (M₁ + M₂)
        ps = Float64.(Float32.(pos_s))
        vs_s = Float64.(Float32.(vel_s))
        ms_s = Float64.(Float32.(m_s))
        mp = Float64(Float32(m_pair))
        dist = vec(sqrt.(sum(abs2, ps; dims = 1)))
        o = sortperm(dist)
        n_in = findfirst(≥(4mp), cumsum(ms_s[o]))
        inner = o[1:n_in]
        m_in = sum(ms_s[inner])
        v̄ = sum(ms_s[inner]' .* vs_s[:, inner]; dims = 2) ./ m_in
        σ_ref = sqrt(sum(ms_s[inner]' .* sum(abs2, vs_s[:, inner] .- v̄; dims = 1)) / (3m_in))
        @test isapprox(r.r_infl_pc, dist[o[n_in]]; rtol = 1e-6)
        @test isapprox(r.sigma_infl_kms, σ_ref * vs; rtol = 1e-5)
        @test isapprox(r.rho_infl_msun_pc3, 3m_in * M_tot / (4π * dist[o[n_in]]^3); rtol = 1e-5)
        # Twice the circular speed: unbound
        u = only(nd.imbh_pairs([pair_snap(0.0, 2v_pair)], pair_names))
        @test !u.bound && u.specific_energy_km2s2 > 0
        @test isnan(u.semi_major_axis_au) && isnan(u.eccentricity)
        @test isnan(u.period_myr) && isnan(u.t_gw_myr)
        # Three names: three pairs, sorted by the names
        rows3 = nd.imbh_pairs([pair_snap(0.0, v_pair; third = true)], [N + 3, N + 1, N + 2])
        @test length(rows3) == 3
        @test [(q.name1, q.name2) for q in rows3] == [(N + 1, N + 2), (N + 1, N + 3), (N + 2, N + 3)]
        @test rows3[1].r_infl_pc == r.r_infl_pc               # the third IMBH is no star
        @test isempty(nd.imbh_pairs([pair_snap(0.0, v_pair)], [N + 1]))
        @test isempty(nd.imbh_pairs([snaps[1]], pair_names))   # N + 2 absent
        # The two IMBHs as the components of the first KS pair
        ks = mksnap(
            0.0,
            vcat(N + 1, N + 2, 1:N),
            vcat(m_pair, m_pair, m_s),
            hcat([-d_pair / 2, 0.0, 0.0], [d_pair / 2, 0.0, 0.0], pos_s),
            hcat([0.0, -v_pair / 2, 0.0], [0.0, v_pair / 2, 0.0], vel_s);
            npairs = 1,
        )
        @test only(nd.imbh_pairs([ks], pair_names)).in_pair
        s = nd.imbh_pair_summary(rows)
        @test s.n_snapshots == 3
        @test s.fraction_bound == 1
        @test s.t_first_bound_myr == 0.0 && s.t_last_bound_myr == 2.0
        @test isapprox(s.a_min_au, d_pair * 206264.806; rtol = 1e-5)
        @test s.e_at_a_min < 1e-5
        @test_throws ArgumentError nd.imbh_pair_summary(rows3)
        s0 = nd.imbh_pair_summary(nd.IMBHPair[])
        @test s0.n_snapshots == 0 && isnan(s0.t_first_bound_myr) && isnan(s0.a_min_au)
    end

    @testset "IMBH pair hardening" begin
        # 1/a = 10⁻³ + 5×10⁻⁴ t [au⁻¹], e = 0.3 + 0.1 ln(1/a)
        function hrow(t; bound = true, name2 = N + 2)
            inv_a = 1.0e-3 + 5.0e-4 * t
            return nd.IMBHPair(;
                time_myr = t,
                name1 = N + 1,
                name2 = name2,
                separation_pc = 1.0e-3,
                specific_energy_km2s2 = bound ? -1.0 : 1.0,
                bound = bound,
                semi_major_axis_au = bound ? 1 / inv_a : NaN,
                eccentricity = bound ? 0.3 + 0.1 * log(inv_a) : NaN,
                period_myr = NaN,
                t_gw_myr = NaN,
                in_pair = false,
                sigma_infl_kms = 3.0,
                rho_infl_msun_pc3 = 2.0e3,
            )
        end
        rows = [hrow(t) for t in 0.0:1.0:20.0]
        h = only(nd.imbh_pair_hardening(rows))
        @test h.name1 == N + 1 && h.name2 == N + 2
        @test h.n == 21
        @test h.t_start_myr == 0.0 && h.t_end_myr == 20.0
        @test isapprox(h.a_mean_au, sum(r -> r.semi_major_axis_au, rows) / 21; rtol = 1e-12)
        @test isapprox(h.s_au_per_myr, 5.0e-4; rtol = 1e-9)
        @test isapprox(h.s_err_au_per_myr, 0.0; atol = 1e-12)
        @test isapprox(h.k_de_dlninva, 0.1; rtol = 1e-9)
        @test h.sigma_infl_kms == 3.0 && h.rho_infl_msun_pc3 == 2.0e3
        @test isapprox(
            h.h_dimensionless,
            nd.hardening_rate_dimensionless(h.s_au_per_myr, 3.0, 2.0e3);
            rtol = 1e-12,
        )
        hw = nd.imbh_pair_hardening(rows; window_myr = 7)
        @test length(hw) == 3
        @test [q.t_start_myr for q in hw] == [0.0, 7.0, 14.0]
        @test [q.t_end_myr for q in hw] == [7.0, 14.0, 20.0]
        @test [q.n for q in hw] == [7, 7, 7]
        @test all(q -> isapprox(q.s_au_per_myr, 5.0e-4; rtol = 1e-9), hw)
        # Unbound rows are ignored
        with_unbound =
            vcat(rows[1:5], [hrow(4.5; bound = false), hrow(10.5; bound = false)], rows[6:end])
        hu = only(nd.imbh_pair_hardening(with_unbound))
        @test hu.n == 21
        @test hu.s_au_per_myr == h.s_au_per_myr
        @test_throws ArgumentError nd.imbh_pair_hardening(vcat(rows, [hrow(0.0; name2 = N + 3)]))
        @test isempty(nd.imbh_pair_hardening(vcat(rows[1:2], [hrow(3.0; bound = false)])))
        @test isempty(nd.imbh_pair_hardening(nd.IMBHPair[]))
        @test isapprox(nd.hardening_rate_dimensionless(1.0e-3, 10.0, 1.0e4), 46.89; rtol = 1e-3)
        @test_throws ArgumentError nd.hardening_rate_dimensionless(1.0e-3, 0.0, 1.0e4)
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
            # Names of several IMBHs from the run record
            names_path = joinpath(dir, "names.toml")
            write(names_path, "[imbh]\nname = 5\nnames = [5, 9]\n")
            @test nd._imbh_names(names_path) == [5, 9]
            @test nd._imbh_name(names_path) == 5
            write(names_path, "[imbh]\nname = 5\n")
            @test nd._imbh_names(names_path) == [5]
            write(names_path, "[meta]\nnbin0 = 0\n")
            @test nd._imbh_names(names_path) == Int[]
            @test nd._imbh_names(joinpath(dir, "absent.toml")) == Int[]
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
                @test startswith(lines[2], "imbh_name," * head)
                @test length(lines) == 2 + n_rows
                @test length(split(lines[3], ",")) == length(split(lines[2], ","))
                @test first(split(lines[3], ",")) == string(bh)
            end
            nd.write_imbh_events(joinpath(dir, "imbh_events.csv"), ev; name = bh)
            @test isfile(joinpath(dir, "imbh_events#1.csv"))
            # Two IMBHs in one product: the name column separates the series
            two = nd.write_imbh_trajectory(
                joinpath(dir, "imbh_trajectory_two.csv"),
                [(bh, traj), (bh + 1, traj)],
            )
            lt = readlines(two)
            @test startswith(lt[1], "# imbh_name=$bh;$(bh + 1) n_states=4")
            @test startswith(lt[2], "imbh_name,time_nb,")
            @test length(lt) == 2 + 4
            @test [first(split(l, ",")) for l in lt[3:end]] == string.([bh, bh, bh + 1, bh + 1])
            le = readlines(
                nd.write_imbh_events(
                    joinpath(dir, "imbh_events_two.csv"),
                    [(bh, ev), (bh + 1, ev)],
                ),
            )
            @test startswith(le[1], "# imbh_name=$bh;$(bh + 1) n_events=2")
            @test length(le) == 2 + 2
            # Pair states
            prow = nd.imbh_pairs([pair_snap(t, v_pair) for t in 0.0:1.0], [N + 1, N + 2])
            lp = readlines(
                nd.write_imbh_pairs(joinpath(dir, "imbh_pairs.csv"), prow; names = [N + 1, N + 2]),
            )
            @test lp[1] == "# imbh_name=$(N + 1);$(N + 2) n_pairs=1 n_rows=2"
            @test lp[2] == join(string.(fieldnames(nd.IMBHPair)), ",")
            @test length(lp) == 2 + 2
            cells = split(lp[3], ",")
            @test length(cells) == fieldcount(nd.IMBHPair)
            @test cells[2:3] == [string(N + 1), string(N + 2)]
            @test cells[6] == "true"
            lp_default =
                readlines(nd.write_imbh_pairs(joinpath(dir, "imbh_pairs_default.csv"), prow))
            @test lp_default[1] == lp[1]
            # Hardening rows: four snapshots of the same circular pair, 1/a constant
            hrows = nd.imbh_pair_hardening(
                nd.imbh_pairs([pair_snap(t, v_pair) for t in 0.0:3.0], [N + 1, N + 2]),
            )
            lh = readlines(
                nd.write_imbh_pair_hardening(
                    joinpath(dir, "imbh_pair_hardening.csv"),
                    hrows;
                    names = [N + 1, N + 2],
                ),
            )
            @test lh[1] == "# imbh_name=$(N + 1);$(N + 2) n_pairs=1 n_rows=1"
            @test lh[2] ==
                  "name1,name2,t_start_myr,t_end_myr,n,a_mean_au,s_au_per_myr," *
                  "s_err_au_per_myr,k_de_dlninva,sigma_infl_kms,rho_infl_msun_pc3,h_dimensionless"
            @test length(lh) == 3
            hc = split(lh[3], ",")
            @test length(hc) == 12
            @test hc[1:2] == [string(N + 1), string(N + 2)]
            @test hc[5] == "4"
            @test hc[9] == ""   # K undefined for a constant 1/a
            @test parse(Float64, hc[10]) > 0 && parse(Float64, hc[11]) > 0
            @test parse(Float64, hc[12]) == 0   # no hardening
        end
    end
end

@testset "TDE events at the static black hole" begin
    nd = Nbody6Dynamics
    mktempdir() do dir
        path = joinpath(dir, "out")
        write(
            path,
            """
             ADJUST:  TIME    1.00000E-02  T[Myr]   0.100E+00
             TDE PARTIAL T[NB] NAME K* M[M*] R*[R*] Q[PC] RT[PC] BETA E[NB] DM[M*] 8.2855224609E-03         1   0  9.99975E-01  0.00000E+00  1.00007E-02  1.00000E-02  9.999E-01 -5.34077E-03  8.08630E-02
             TDE        T[NB] NAME K* M[M*] R*[R*] Q[PC] RT[PC] BETA E[NB] MBH[M*] 4.1399511229E-02       142   1  9.99975E-03  1.20000E+00  3.31501E-06  1.00000E-05  3.017E+00 -7.39236E+00  1.00000E+04
             TDE PLUNGE T[NB] NAME K* M[M*] R*[R*] Q[PC] RT[PC] BETA E[NB] MBH[M*] 6.7252578025E-02       140  13  1.40000E+00  1.00000E-05  4.36891E-10  1.00000E-11  2.289E-02 -5.57977E+00  1.00000E+04
             TDE        T[NB] NAME K* M[M*] R*[R*] Q[PC] RT[PC] BETA E[NB] MBH[M*] ****************       141   1  1.0E+00
            """,
        )
        ev = @test_logs (:warn, r"1 TDE lines") nd.read_tde_events(path; tstar_myr = 2.0)
        @test [e.kind for e in ev] == [:partial, :full, :plunge]
        p, f, g = ev
        @test p.name == 1 &&
              p.kstar == 0 &&
              p.time_nb ≈ 8.2855224609e-3 &&
              p.time_myr ≈ 2 * p.time_nb
        @test p.dm_msun ≈ 8.0863e-2 && isnan(p.mbh_msun) && p.beta ≈ 0.9999
        @test f.name == 142 && f.kstar == 1 && f.radius_rsun ≈ 1.2 && f.q_pc ≈ 3.31501e-6
        @test f.r_t_pc ≈ 1.0e-5 &&
              f.e_nb ≈ -7.39236 &&
              f.dm_msun == f.mass_msun &&
              f.mbh_msun ≈ 1.0e4
        @test g.kstar == 13 && g.dm_msun ≈ 1.4
        s = nd.tde_summary(ev)
        @test (s.n_full, s.n_partial, s.n_plunge) == (1, 1, 1)
        @test s.dm_msun ≈ 8.0863e-2 + 9.99975e-3 + 1.4
        @test s.rate_per_myr ≈ 2 / (2 * (6.7252578025e-2 - 8.2855224609e-3))
        @test isnan(nd.tde_summary(ev[1:1]).rate_per_myr)
        @test_throws ErrorException nd.read_tde_events(joinpath(dir, "absent"))
        csv = joinpath(dir, "tde_events.csv")
        nd.write_tde_events(csv, ev)
        lines = readlines(csv)
        @test lines[1] == "# n_full=1 n_partial=1 n_plunge=1 dm_msun=$(s.dm_msun)"
        @test split(lines[2], ',') == string.(collect(fieldnames(nd.TDEEvent)))
        @test length(lines) == 5
        @test split(lines[3], ',')[1] == "partial" && split(lines[3], ',')[end] == ""
        rates = nd.tde_rate_series(ev; bin_myr = 0.05)
        @test length(rates) == 3
        @test [(r.n_partial, r.n_full, r.n_plunge) for r in rates] == [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        @test all(r -> r.rate_per_myr ≈ 20.0, rates) && rates[3].dm_msun ≈ 1.4
        @test length(nd.tde_rate_series(ev; bin_myr = 0.05, t_end_myr = 1.0)) == 20
        @test_throws ArgumentError nd.tde_rate_series(ev; bin_myr = 0.0)
        @test_throws ArgumentError nd.tde_rate_series(nd.read_tde_events(path); bin_myr = 1.0)
        nd.write_tde_events(csv, ev[2:2])
        @test isfile(joinpath(dir, "tde_events#1.csv")) && length(readlines(csv)) == 3
    end
end
