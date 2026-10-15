# Angular-momentum budget and resolved rotation on a static rotating model: the budget closes
# when stars are removed and masses are reduced, the rotation measures sit above their floors.

using StableRNGs: StableRNG
using LinearAlgebra: ×, norm

function _am_snapshot(mass, pos, vel; t = 0.0, names = 1:length(mass))
    params = zeros(Float32, 20)
    params[1] = t
    params[3] = 1.0f0
    params[4] = 1.0f0
    params[11] = 1.0f0
    params[12] = 1.0f0
    n = length(mass)
    hdr = SnapshotHeader(Int32(n), Int32(1), Int32(1), Int32(20), params)
    return Snapshot(
        hdr,
        Int32.(collect(names)),
        Float32.(mass),
        Float32.(pos),
        Float32.(vel),
        Float32[],
        Float32[],
        0,
        Int32[],
    )
end

@testset "angular_momentum_budget closes on a static series" begin
    rng = StableRNG(2201)
    N = 6000
    pos, vel = Nbody6Dynamics.sample_plummer(N, 1.0; rng)
    lynden_bell_flip!(vel, pos, [0.0, 0.0, 1.0], 1.0; rng)
    m_solar = Nbody6Dynamics.sample_kroupa(N; m_low = 0.08, m_up = 50.0, rng)
    M_solar = sum(m_solar)
    mass = m_solar ./ M_solar
    units = UnitScaling(1.0, M_solar, 1.0, 1.0)

    r = vec(sqrt.(sum(abs2, pos; dims = 1)))
    keep = sortperm(r)[1:(N - 200)]
    gone = sortperm(r)[(N - 199):N]
    m3 = copy(mass[keep])
    shed = 1:300
    m3[shed] .*= 0.9
    snaps = [
        _am_snapshot(mass, pos, vel; t = 0.0),
        _am_snapshot(mass[keep], pos[:, keep], vel[:, keep]; t = 1.0, names = keep),
        _am_snapshot(m3, pos[:, keep], vel[:, keep]; t = 2.0, names = keep),
    ]
    bsets = [bound_set(s) for s in snaps]
    @test all(b -> length(b.stars) > 0.9 * nparticles(snaps[1]), bsets)

    settings = PostprocessConfig(; rotation_n_min = 500, rotation_n_null = 10)
    b = angular_momentum_budget(snaps, bsets, units; settings, n_min_class = 100)
    @test b isa AngularMomentumBudget
    @test b.time_nb == [0.0, 1.0, 2.0]
    L0 = sqrt(sum(abs2, b.l_total[:, 1]))
    @test L0 > 0
    @test maximum(abs, b.residual) < 1e-8 * L0

    # Removed stars: the internal angular momentum of the first snapshot with and without them;
    # shed mass: the internal angular momentum of the kept stars with the two mass vectors
    p32 = Float64.(Float32.(pos))
    v32 = Float64.(Float32.(vel))
    m32 = Float64.(Float32.(mass))
    L_int = Nbody6Dynamics._internal_angular_momentum
    L_all = L_int(p32, v32, 1:N, m32)
    L_kept = L_int(p32, v32, keep, m32[keep])
    @test b.l_total[:, 1] ≈ L_all rtol = 1e-9
    @test b.l_removed[:, 2] ≈ L_all .- L_kept rtol = 1e-9
    @test b.l_removed[:, 3] ≈ L_all .- L_kept rtol = 1e-9
    @test maximum(abs, b.l_added) < 1e-12 * L0
    @test maximum(abs, b.l_stellar[:, 1:2]) < 1e-12 * L0
    @test all(iszero, b.l_remnant)
    @test all(b.l_scale .> norm.(eachcol(b.l_total)))

    # A star turned into a kicked remnant: named by an engine event, its change is its own term
    kicked = keep[length(keep) ÷ 2]   # a star at the median radius, so the kick carries angular momentum
    m4 = copy(m3)
    v4 = copy(vel[:, keep])
    q = findfirst(==(kicked), keep)
    m4[q] *= 0.1
    v4[3, q] += 50.0
    snap4 = _am_snapshot(m4, pos[:, keep], v4; t = 3.0, names = keep)
    events =
        [EngineEvent(; time_nb = 2.5, kind = :kick_full, interval = 3, line = 1, name1 = kicked)]
    bk = angular_momentum_budget(
        [snaps; snap4],
        [bsets; bound_set(snap4)],
        units;
        events,
        settings,
        n_min_class = 100,
    )
    @test maximum(abs, bk.residual) < 1e-8 * L0
    @test norm(bk.l_remnant[:, 4]) > 1e-4 * L0 && maximum(abs, bk.l_remnant[:, 1:3]) < 1e-12 * L0
    @test bk.l_stellar[:, 4] ≈ bk.l_stellar[:, 3] rtol = 1e-12
    # Without the event the same change lands in the residual
    bn = angular_momentum_budget([snaps; snap4], [bsets; bound_set(snap4)], units; settings)
    @test norm(bn.residual[:, 4]) > 1e-4 * L0 && maximum(abs, bn.l_remnant) < 1e-12 * L0
    m3f = Float64.(Float32.(m3))
    @test b.l_stellar[:, 3] ≈ L_kept .- L_int(p32, v32, keep, m3f) rtol = 1e-9

    # Elapsed relaxation times accumulate; the spin axis is the flip axis
    @test b.n_rh[:, 1] == [0.0, 0.0]
    @test all(b.n_rh[:, 2] .> 0) && all(b.n_rh[:, 3] .> b.n_rh[:, 2])
    @test b.n_rh[1, 3] < b.n_rh[2, 3]   # γ = 0.02 gives the longer t_rh
    @test all(b.axis[3, :] .> 0.99)
    @test all(isfinite, b.r_lagr) && all(b.r_lagr[1, :] .< b.r_lagr[2, :] .< b.r_lagr[3, :])
    # Maximal rotation inside every Lagrangian radius and in every populated class
    @test all(b.lambda_r_radius .> 0.5) && all(b.lambda_r_radius .> 3 .* b.lambda_r_radius_floor)
    @test sum(b.n_class[:, 1]) == length(bsets[1].stars)
    @test sum(b.l_z_class_fraction[:, 1]) ≈ 1 atol = 1e-9
    for c in 1:size(b.n_class, 1), k in 1:3
        if b.n_class[c, k] >= 100
            @test b.lambda_r_class[c, k] > 5 * b.lambda_r_class_floor[c, k]
        else
            @test isnan(b.lambda_r_class[c, k])
        end
    end
    @test b.class_edges == [0.5, 5.0]

    # Errors
    @test_throws ArgumentError angular_momentum_budget(Snapshot[], BoundSet[], units)
    @test_throws ArgumentError angular_momentum_budget(snaps, bsets[1:2], units)
    @test_throws ArgumentError angular_momentum_budget(snaps, bsets, units; fractions = (0.5, 0.1))
    @test_throws ArgumentError angular_momentum_budget(
        snaps,
        bsets,
        units;
        class_edges_solar = (5.0, 0.5),
    )
    @test_throws ArgumentError angular_momentum_budget(snaps, bsets, units; n_min_class = 1)

    # Writer
    mktempdir() do dir
        path = joinpath(dir, "angular_momentum.csv")
        @test write_angular_momentum(path, b) == path
        lines = readlines(path)
        @test length(lines) == 5
        @test startswith(lines[1], "# fractions=0.1;0.5;0.9 class_edges_solar=0.5;5")
        header = split(lines[2], ",")
        @test length(header) == 29 + 3 * 3 + 5 * 3
        @test header[1] == "time_nb" &&
              header[30] == "r_lagr_1" &&
              header[end] == "lambda_r_class_floor_3"
        @test header[17:20] == ["l_remnant_x", "l_remnant_y", "l_remnant_z", "residual_x"]
        @test header[23] == "l_scale"
        header[26] == "r_lagr_1" && header[end] == "lambda_r_class_floor_3"
        @test all(l -> length(split(l, ",")) == length(header), lines[3:5])
    end
end
