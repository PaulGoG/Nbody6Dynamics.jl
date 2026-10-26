# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Collision catalogue: progenitor tags, primordial-pair and mixed-origin
# flags, product names and generations, kicks before the event, the summary,
# the CSV writer and the argument checks on synthetic coll.13 and coal.24
# records.

@testset "Collision catalogue" begin
    nd = Nbody6Dynamics
    tstar = 0.87
    rbar = 2.5
    ranges = [collect(1:100), collect(101:260)]
    n_primordial = 10

    coll(t, n1, n2; m1 = 2.0, m2 = 1.0, dm = 0.01) = nd.CollisionRecord(;
        time_nb = t,
        name1 = n1,
        name2 = n2,
        kstar1 = 1,
        kstar2 = 1,
        kstar_product = 1,
        mass1_msun = m1,
        mass2_msun = m2,
        mass_product_msun = m1 + m2 - dm,
        mass_lost_msun = dm,
        radius1_rsun = 1.5,
        radius2_rsun = 1.0,
        distance_over_core_radius = 0.4,
        separation_rsun = 2.0,
        eccentricity = 0.9,
        period_d = 0.5,
    )
    coal(t, n1, n2) = nd.CoalescenceRecord(;
        channel = :binary,
        iqcoll = 0,
        time_nb = t,
        name1 = n1,
        name2 = n2,
        kstar1 = 13,
        kstar2 = 14,
        kstar_product = 14,
        mass1_bse_msun = 1.4,
        mass2_bse_msun = 10.0,
        separation_nb = 1e-8,
        eccentricity = 0.0,
        semi_major_axis_nb = 1e-8,
        binding_energy_nb = -1.0,
        tidal_potential_correction_nb = 0.0,
        period_d = 1e-4,
        mass1_msun = 1.4,
        mass2_msun = 10.0,
        mass_product_msun = 11.2,
        mass_lost_msun = 0.2,
        radius_survivor_rsun = 1e-5,
        collision_distance_rsun = 1e-4,
        velocity_infinity_kms = 0.0,
        distance_density_centre_nb = 0.1,
        speed_nb = 1.0,
    )
    kick(t, n1, n2) = nd.MassLossEvent(;
        line = 1,
        interval = 1,
        time_nb = NaN,
        kind = :binary_kick,
        time_myr = t * tstar,
        name1 = n1,
        name2 = n2,
    )

    collisions = [
        coll(3.0, 3, 7; dm = 0.03),
        coll(1.0, 3, 4; dm = 0.01),
        coll(4.0, 4, 5; dm = 0.04),
        coll(2.0, 50, 150; dm = 0.02),
    ]
    coalescences = [coal(0.5, 200, 201)]
    build(; tstar_myr = tstar, nbin0 = n_primordial, rbar_pc = NaN, mass_loss = nothing) =
        nd.collision_catalogue(
            collisions,
            coalescences,
            ranges;
            tstar_myr = tstar_myr,
            nbin0 = nbin0,
            rbar_pc = rbar_pc,
            mass_loss = mass_loss,
        )

    @testset "Rows" begin
        rows = build(; rbar_pc = rbar)
        @test length(rows) == 5
        @test [r.time_nb for r in rows] == [0.5, 1.0, 2.0, 3.0, 4.0]
        @test all(r -> r.time_myr == r.time_nb * tstar, rows)
        c, p, m, g, q = rows
        @test c.kind === :coalescence && c.channel === :binary
        @test (c.name1, c.name2) == (200, 201)
        @test (c.progenitor1, c.progenitor2) == (2, 2)
        @test !c.mixed && !c.primordial_pair
        @test c.separation_rsun ≈ 1e-8 * rbar * 4.435e7 rtol = 1e-12
        @test c.r_centre == 0.1
        @test all(r -> r.kind === :collision && r.channel === :none, rows[2:end])
        @test (p.name1, p.name2) == (3, 4)
        @test (p.progenitor1, p.progenitor2) == (1, 1)
        @test p.primordial_pair && !p.mixed
        @test p.r_centre == 0.4 && p.separation_rsun == 2.0
        @test (m.progenitor1, m.progenitor2) == (1, 2)
        @test m.mixed && !m.primordial_pair
        @test !q.primordial_pair && !q.mixed
        @test [r.generation for r in rows] == [1, 1, 1, 2, 1]
        @test !any(r -> r.kicked_before, rows)
        @test isnan(first(build()).separation_rsun)
    end

    @testset "Product name" begin
        # coll.13: the product takes the name of the heavier body
        r = only(
            nd.collision_catalogue(
                [coll(1.0, 3, 120; m1 = 0.5, m2 = 1.5)],
                nd.CoalescenceRecord[],
                ranges;
                tstar_myr = tstar,
            ),
        )
        @test (r.name1, r.name2) == (120, 3)
        @test (r.m1_msun, r.m2_msun) == (1.5, 0.5)
        @test (r.progenitor1, r.progenitor2) == (2, 1)
        tie = only(
            nd.collision_catalogue(
                [coll(1.0, 3, 120; m1 = 1.0, m2 = 1.0)],
                nd.CoalescenceRecord[],
                ranges;
                tstar_myr = tstar,
            ),
        )
        @test (tie.name1, tie.name2) == (3, 120)
        @test nd.progenitor_of(150, ranges) == 2
        @test nd.progenitor_of(1, ranges) == 1
        @test nd.progenitor_of(999, ranges) == 0
    end

    @testset "Kicks before the event" begin
        rows = build(; mass_loss = [kick(1.5, 150, 151)])
        @test [r.kicked_before for r in rows] == [false, false, true, false, false]
        rows = build(; mass_loss = [kick(2.5, 150, 151)])
        @test !any(r -> r.kicked_before, rows)
        mdend = nd.MassLossEvent(;
            line = 2,
            interval = 1,
            time_nb = 1.0,
            kind = :ks_end_mass_loss,
            name1 = 150,
            name2 = 151,
        )
        @test !any(r -> r.kicked_before, build(; mass_loss = [mdend]))
    end

    @testset "Unknown progenitor" begin
        r = only(
            nd.collision_catalogue(
                [coll(5.0, 999, 7)],
                nd.CoalescenceRecord[],
                ranges;
                tstar_myr = tstar,
            ),
        )
        @test (r.progenitor1, r.progenitor2) == (0, 1)
        @test !r.mixed
    end

    @testset "Summary" begin
        rows = build(; mass_loss = [kick(1.5, 150, 151)])
        s = nd.collision_summary(rows)
        @test s.n_collisions == 4 && s.n_coalescences == 1
        @test s.n_mixed == 1
        @test s.n_primordial_pairs == 1
        @test s.n_kicked_before == 1
        @test s.m_lost_total_msun ≈ sum(r.m_lost_msun for r in rows) rtol = 1e-14
        @test s.m_lost_total_msun ≈ 0.3 rtol = 1e-12
        @test s.t_first_myr == 0.5 * tstar && s.t_last_myr == 4.0 * tstar
        @test s.by_progenitor_pair == Dict((2, 2) => 1, (1, 1) => 3, (1, 2) => 1)
        @test s.by_pair_class == Dict(:star_star => 4, :bhns => 1)
        @test s.max_generation == 2
        e = nd.collision_summary(nd.CollisionRow[])
        @test e.n_collisions == 0 && e.max_generation == 0
        @test isnan(e.t_first_myr) && e.m_lost_total_msun == 0
    end

    @testset "Writer" begin
        rows = build()
        path = nd.write_collision_catalogue(joinpath(mktempdir(), "collisions.csv"), rows)
        lines = readlines(path)
        @test first(lines) == join(fieldnames(nd.CollisionRow), ",")
        @test length(lines) == 1 + length(rows)
        header = split(lines[1], ",")
        cells = split(lines[2], ",")
        @test length(cells) == fieldcount(nd.CollisionRow)
        @test cells[findfirst(==("kind"), header)] == "coalescence"
        @test cells[findfirst(==("separation_rsun"), header)] == "NaN"
        @test cells[findfirst(==("kicked_before"), header)] == "false"
    end

    @testset "Argument checks" begin
        @test_throws ArgumentError build(; tstar_myr = 0.0)
        @test_throws ArgumentError build(; nbin0 = -1)
        @test_throws ArgumentError nd.collision_catalogue(
            collisions,
            coalescences,
            [collect(1:100), collect(90:200)];
            tstar_myr = tstar,
        )
    end
end
