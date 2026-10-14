# =============================================================================
# Potential provider, reduced systems, bound sets and snapshot integrals
# =============================================================================

"""Snapshot of the given stars (Float32 storage, as read from conf.3) at time `t`."""
function _integrals_snapshot(mass, pos, vel; t = 0.0, npairs = 0, pair_name = Int32[])
    params = zeros(Float32, 20)
    params[1] = t
    params[2] = npairs
    params[3] = 1.0f0
    params[4] = 1.0f0
    params[11] = 1.0f0
    params[12] = 1.0f0
    n = length(mass)
    hdr = SnapshotHeader(Int32(n + npairs), Int32(1), Int32(1), Int32(20), params)
    return Snapshot(
        hdr,
        Int32.(1:n),
        Float32.(mass),
        Float32.(pos),
        Float32.(vel),
        Float32[],
        Float32[],
        npairs,
        pair_name,
    )
end

"""Relative positions and velocities of a Kepler orbit (G = 1) at eccentric anomaly `E`."""
function _kepler_relative(M, a, e, E)
    n = sqrt(M / a^3)
    Edot = n / (1 - e * cos(E))
    x = [a * (cos(E) - e), a * sqrt(1 - e^2) * sin(E), 0.0]
    v = [-a * sin(E) * Edot, a * sqrt(1 - e^2) * cos(E) * Edot, 0.0]
    return x, v
end

@testset "Potential provider" begin
    using StableRNGs
    # Two and three bodies against closed forms.
    @test pair_potential([0.0 1.0; 0.0 0.0; 0.0 0.0], [1.0, 2.0]) ≈ [-2.0, -1.0]
    x3 = [0.0 3.0 0.0; 0.0 0.0 4.0; 0.0 0.0 0.0]
    m3 = [1.0, 2.0, 3.0]
    @test pair_potential(x3, m3) ≈ [-(2 / 3 + 3 / 4), -(1 / 3 + 3 / 5), -(1 / 4 + 2 / 5)]
    # Coincident bodies do not interact; empty inputs give empty or zero potentials.
    @test pair_potential([0.0 0.0; 0.0 0.0; 0.0 0.0], [1.0, 1.0]) == [0.0, 0.0]
    @test pair_potential(zeros(3, 0), Float64[]) == Float64[]
    @test pair_potential(x3, zeros(3, 0), Float64[]) == zeros(3)
    @test_throws DimensionMismatch pair_potential(zeros(2, 3), ones(3))
    @test_throws DimensionMismatch pair_potential(x3, x3, ones(2))

    # Agreement with the serial symmetric sum, and the subset form by subtraction.
    rng = StableRNG(21)
    pos, _ = sample_plummer(2000, 1.0; rng = rng)
    mass = 0.5 .+ rand(rng, 2000)
    Φ = pair_potential(pos, mass)
    Φref = zeros(2000)
    for i in 1:1999, j in (i + 1):2000
        r = sqrt(sum(abs2, pos[:, i] .- pos[:, j]))
        Φref[i] -= mass[j] / r
        Φref[j] -= mass[i] / r
    end
    @test maximum(abs.(Φ .- Φref) ./ abs.(Φref)) < 1e-12
    set = 1:1500
    rest = 1501:2000
    Φset = pair_potential(pos[:, set], mass[set])
    Φsub = Φ[set] .- pair_potential(pos[:, set], pos[:, rest], mass[rest])
    @test maximum(abs.(Φsub .- Φset) ./ abs.(Φset)) < 1e-12
    # Float32 input is promoted, not truncated further.
    @test pair_potential(Float32.(pos), Float32.(mass)) ≈
          pair_potential(Float64.(Float32.(pos)), Float64.(Float32.(mass)))
end

@testset "Reduced systems" begin
    m = [1.0, 3.0, 2.0, 5.0]
    x = [0.0 4.0 7.0 9.0; 0.0 0.0 1.0 2.0; 0.0 0.0 0.0 0.0]
    v = [1.0 -1.0 0.5 0.0; 0.0 0.0 2.0 1.0; 0.0 0.0 0.0 0.0]
    sys = reduce_pairs(m, x, v, [(1, 2)])
    @test sys.mass == [4.0, 2.0, 5.0]
    @test sys.pos[:, 1] ≈ [3.0, 0.0, 0.0]
    @test sys.vel[:, 1] ≈ [-0.5, 0.0, 0.0]
    @test sys.star1 == [1, 3, 4] && sys.star2 == [2, 0, 0]
    @test_throws ArgumentError reduce_pairs(m, x, v, [(1, 2), (2, 3)])
    @test_throws ArgumentError reduce_pairs(m, x, v, [(1, 5)])
    @test_throws ArgumentError reduce_pairs(m, x, v, [(2, 2)])
    @test_throws DimensionMismatch reduce_pairs(m, x[:, 1:3], v, NTuple{2,Int}[])
end

@testset "Bound set" begin
    using StableRNGs
    rng = StableRNG(5)
    N = 2000
    a = 3π / 16
    pos, vel = sample_plummer(N, a; rng = rng)
    mass = fill(1.0 / N, N)

    # Static Plummer plus 100 stars at twice the local escape speed of the model.
    n_inj = 100
    pinj, _ = sample_plummer(n_inj, a; rng = rng)
    vinj = similar(pinj)
    for i in 1:n_inj
        r = sqrt(sum(abs2, pinj[:, i]))
        u = randn(rng, 3)
        vinj[:, i] = 2 * sqrt(2 / sqrt(r^2 + a^2)) .* u ./ sqrt(sum(abs2, u))
    end
    m_all = vcat(mass, fill(1.0 / N, n_inj))
    x_all = hcat(pos, pinj)
    v_all = hcat(vel, vinj)
    b = bound_set(m_all, x_all, v_all)
    @test b.status === :converged
    @test b.passes ≥ 2
    @test isempty(intersect(b.stars, (N + 1):(N + n_inj)))   # every injected star removed
    @test length(b.stars) ≥ 0.98N                             # the model itself stays bound
    @test issorted(b.stars)
    # Incremental potential equals the exact potential of the final set.
    sb = b.stars
    Wexact = 0.5 * sum(m_all[sb] .* pair_potential(x_all[:, sb], m_all[sb]))
    @test b.potential ≈ Wexact rtol = 1e-12
    vref = [sum(m_all[sb] .* v_all[k, sb]) / sum(m_all[sb]) for k in 1:3]
    @test collect(b.v_ref) ≈ vref
    @test b.kinetic ≈ 0.5 * sum(m_all[i] * sum(abs2, v_all[:, i] .- vref) for i in sb)
    # A precomputed potential gives the same result and is not modified.
    Φ0 = pair_potential(x_all, m_all)
    Φcopy = copy(Φ0)
    @test bound_set(reduce_pairs(m_all, x_all, v_all, NTuple{2,Int}[]); potential = Φ0).stars ==
          b.stars
    @test Φ0 == Φcopy

    # An escaping hard pair: bound star by star (internal binding), removed as a centre of mass.
    mp = 5.0 / N
    r12 = 1e-4
    vorb = sqrt(2mp / r12)
    xc = [0.3, 0.0, 0.0]
    vc = [0.0, 3 * sqrt(2 / sqrt(0.09 + a^2)), 0.0]
    xp = hcat(xc .+ [0.0, 0.0, r12 / 2], xc .- [0.0, 0.0, r12 / 2])
    vp = hcat(vc .+ [vorb / 2, 0.0, 0.0], vc .- [vorb / 2, 0.0, 0.0])
    m2 = vcat(mass, [mp, mp])
    x2 = hcat(pos, xp)
    v2 = hcat(vel, vp)
    star_level = bound_set(m2, x2, v2)
    @test (N + 1) in star_level.stars && (N + 2) in star_level.stars
    pair_level = bound_set(m2, x2, v2; pairs = [(N + 1, N + 2)])
    @test !((N + 1) in pair_level.stars) && !((N + 2) in pair_level.stars)
    # A bound hard pair at rest at the centre stays, with both components.
    x3 = hcat(pos, [r12 / 2 -r12 / 2; 0.0 0.0; 0.0 0.0])
    v3 = hcat(vel, [0.0 0.0; vorb/2 -vorb/2; 0.0 0.0])
    b3 = bound_set(m2, x3, v3; pairs = [(N + 1, N + 2)])
    @test (N + 1) in b3.stars && (N + 2) in b3.stars

    # Pass limit and dissolution are reported, never silent.
    b1 = @test_logs (:warn, r"not converged") bound_set(m_all, x_all, v_all; max_passes = 1)
    @test b1.status === :pass_limit && b1.passes == 1
    hot = bound_set(mass[1:50], pos[:, 1:50], 100 .* vel[:, 1:50])
    @test hot.status === :dissolved
    @test_throws ArgumentError bound_set(m_all, x_all, v_all; max_passes = 0)
    @test_throws ArgumentError bound_set(m_all, x_all, v_all; min_members = 1)

    # Snapshot form: snapshot indices, pairs split by a subset are left out.
    snap = _integrals_snapshot(m2, x3, v3; npairs = 0)
    @test bound_set(snap).stars == bound_set(m2, x3, v3).stars
    # Pair components first, as in the engine: move the pair to stars 1 and 2.
    order = vcat([N + 1, N + 2], 1:N)
    snap_p = _integrals_snapshot(
        m2[order],
        x3[:, order],
        v3[:, order];
        npairs = 1,
        pair_name = Int32[N + 3],
    )
    bp = bound_set(snap_p)
    @test 1 in bp.stars && 2 in bp.stars && isempty(bp.excluded)
    split = bound_set(snap_p, vcat(1, 3:(N + 2)))
    @test split.excluded == [1]
    @test !(1 in split.stars)
end

@testset "Snapshot integrals" begin
    using StableRNGs
    # Kepler two-body orbit sampled at eight phases: E and L constant, closed forms.
    m1, m2 = 0.6, 0.4
    M = m1 + m2
    a, e = 1.0, 0.5
    μ = m1 * m2 / M
    E_kep = -m1 * m2 / (2a)
    L_kep = μ * sqrt(M * a * (1 - e^2))
    for E in range(0, 2π; length = 9)[1:8]
        xr, vr = _kepler_relative(M, a, e, E)
        pos = hcat(m2 / M .* xr, -m1 / M .* xr)
        vel = hcat(m2 / M .* vr, -m1 / M .* vr)
        s = snapshot_integrals(_integrals_snapshot([m1, m2], pos, vel); min_members = 2)
        @test s.kinetic + s.potential ≈ E_kep rtol = 1e-5
        @test s.angular_momentum[3] ≈ L_kep rtol = 1e-5
        @test s.pair_internal == 0.0
        # The same two stars as a regularised pair: one system, energy in the pair term.
        sp =
            snapshot_integrals(_integrals_snapshot([m1, m2], pos, vel; npairs = 1); min_members = 2)
        @test sp.potential == 0.0
        @test sp.pair_internal ≈ E_kep rtol = 1e-5
        @test sp.kinetic ≈ 0.0 atol = 1e-12
        @test sp.kinetic_stars ≈ s.kinetic_stars
    end

    # Plummer: invariance under a rotation and a Galilean boost.
    rng = StableRNG(9)
    N = 1000
    pos, vel = sample_plummer(N, 3π / 16; rng = rng)
    mass = fill(1.0 / N, N)
    s0 = snapshot_integrals(_integrals_snapshot(mass, pos, vel))
    @test s0.mass ≈ 1.0 rtol = 1e-6
    @test s0.n_stars == N && s0.n_pairs == 0 && s0.n_hierarchies == 0
    θ = 0.7
    R = [cos(θ) -sin(θ) 0.0; sin(θ) cos(θ) 0.0; 0.0 0.0 1.0]
    sr = snapshot_integrals(_integrals_snapshot(mass, R * pos, R * vel))
    @test collect(sr.angular_momentum) ≈ R * collect(s0.angular_momentum) rtol = 1e-4
    @test sr.potential ≈ s0.potential rtol = 1e-6
    @test sr.kinetic ≈ s0.kinetic rtol = 1e-6
    u = [0.3, -0.2, 0.1]
    sb = snapshot_integrals(_integrals_snapshot(mass, pos, vel .+ u))
    @test collect(sb.momentum) ≈ collect(s0.momentum) .+ s0.mass .* u rtol = 1e-5
    @test collect(sb.angular_momentum_com) ≈ collect(s0.angular_momentum_com) rtol = 1e-4 atol =
        1e-7
    @test sb.potential ≈ s0.potential rtol = 1e-6
    @test sb.bound.kinetic ≈ s0.bound.kinetic rtol = 1e-5
    @test sb.bound.stars == s0.bound.stars

    # CSV: header, one row per snapshot, L_removed relative to the first row.
    path = joinpath(TESTDIR, "snapshot_integrals.csv")
    write_snapshot_integrals(path, [s0, sr])
    lines = readlines(path)
    @test length(lines) == 3
    header = split(lines[1], ",")
    @test header[1:3] == ["time_nb", "time_myr", "n_stars"]
    @test "e_total" in header && "bound_status" in header && "lz_removed" in header
    row1 = split(lines[2], ",")
    @test parse(Float64, row1[findfirst(==("lz_removed"), header)]) == 0.0
    @test row1[findfirst(==("bound_status"), header)] == "converged"
    write_snapshot_integrals(path, [s0])
    @test isfile(joinpath(TESTDIR, "snapshot_integrals#1.csv"))
    @test length(readlines(path)) == 2
end
