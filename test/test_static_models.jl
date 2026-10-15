# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Static analytic test models with known answers (estimator calibration,
# literature note §14): every estimator of the revision is checked against a
# closed form or a quadrature on a model whose answer is known. Samples use
# StableRNGs so the reference values hold across Julia versions.

using StableRNGs: StableRNG
using Nbody6Dynamics:
    sample_plummer_om,
    lynden_bell_flip!,
    stretch_positions!,
    rank_masses_by_energy!,
    shrinking_sphere_centre,
    lagrangian_grid,
    relaxation_times,
    local_relaxation_time,
    shape_tensor,
    projected_ellipticity,
    aplus,
    two_sample_d,
    mixing,
    local_density,
    snapshot_pairs,
    reconcile_pairs,
    core_collapse_markers,
    _local_density,
    _lagrangian_radii,
    _PLUMMER_RHM_OVER_A

"""Plummer radius enclosing the fraction `f` of the mass, scale `a`."""
_plummer_rf(f, a) = a / sqrt(f^(-2 / 3) - 1)
"""Plummer enclosed mass fraction inside `r`, scale `a`."""
_plummer_mf(r, a) = r^3 / (r^2 + a^2)^1.5
"""Plummer density, total mass 1, scale `a`."""
_plummer_rho(r, a) = 3 / (4π * a^3) * (1 + (r / a)^2)^(-2.5)

"""Snapshot of stars with names `names`, masses `mass`, positions and velocities (3 × n), time `t`."""
function _sm_snapshot(names, mass, pos, vel, t)
    params = zeros(Float32, 20)
    params[1] = Float32(t)
    params[3] = 1.0f0
    params[4] = 1.0f0
    params[11] = 1.0f0
    params[12] = 1.0f0
    n = length(mass)
    hdr = SnapshotHeader(Int32(n), Int32(0), Int32(1), Int32(20), params)
    return Snapshot(
        hdr,
        Int32.(names),
        Float32.(mass),
        Float32.(pos),
        Float32.(vel),
        Float32[],
        Float32[],
    )
end

@testset "Static models — Osipkov–Merritt Plummer" begin
    rng = StableRNG(11)
    N = 60_000
    a = 1.0
    @test_throws ArgumentError sample_plummer_om(10, a, 0.5a)
    for r_a in (Inf, 1.5a, 0.75a)
        pos, vel = sample_plummer_om(N, a, r_a; rng)
        @test size(pos) == (3, N) && size(vel) == (3, N)
        m = fill(1.0 / N, N)
        c = zeros(3)
        # Density: number-weighted Lagrangian radii against the Plummer closed form
        g = lagrangian_grid(pos, m, c; weights = :number)
        for (k, f) in enumerate(g.fractions)
            g.masked[k] && continue
            @test isapprox(
                log(g.radii[k]),
                log(_plummer_rf(f, a));
                atol = 4 * g.sigma_log[k] + 0.01,
            )
        end
        # Anisotropy β(r) = r²/(r² + r_a²) on equal-membership shells
        prof = radial_profile(pos, vel, m, c; nbins = 12, n_min_beta = 300)
        for k in eachindex(prof.r)
            isnan(prof.beta[k]) && continue
            β_true = isinf(r_a) ? 0.0 : prof.r[k]^2 / (prof.r[k]^2 + r_a^2)
            @test isapprox(prof.beta[k], β_true; atol = 3 * prof.beta_err[k] + 0.02)
        end
        if isinf(r_a)
            # Isotropic radial dispersion σ_r² = GM/(6√(r² + a²))
            for k in eachindex(prof.r)
                isnan(prof.sigma_r[k]) && continue
                σ_true = sqrt(1 / (6 * sqrt(prof.r[k]^2 + a^2)))
                @test isapprox(prof.sigma_r[k], σ_true; rtol = 3 * prof.sigma_err[k] + 0.03)
            end
        end
        # Total energy of the sample: a Plummer sphere has E = −3π/64 GM²/a
        T = 0.5 * sum(m[i] * sum(abs2, @view vel[:, i]) for i in 1:N)
        W = 0.5 * sum(m .* pair_potential(pos, m))
        @test isapprox(T + W, -3π / 64; rtol = 0.03)
    end
end

@testset "Static models — Lynden-Bell rotation flip" begin
    rng = StableRNG(12)
    N = 20_000
    pos, vel = sample_plummer(N, 1.0; rng)
    m = fill(1.0 / N, N)
    speeds = [sqrt(sum(abs2, @view vel[:, i])) for i in 1:N]
    v1 = copy(vel)
    lynden_bell_flip!(v1, pos, [0.0, 0.0, 1.0], 1.0; rng)
    vφ(v, x) = (x[1] * v[2] - x[2] * v[1]) / hypot(x[1], x[2])
    @test all(vφ(view(v1, :, i), view(pos, :, i)) ≥ -1e-12 for i in 1:N)
    @test maximum(abs, [sqrt(sum(abs2, @view v1[:, i])) - speeds[i] for i in 1:N]) < 1e-10
    v0 = copy(vel)
    lynden_bell_flip!(v0, pos, [0.0, 0.0, 1.0], 0.0; rng)
    @test all(vφ(view(v0, :, i), view(pos, :, i)) ≤ 1e-12 for i in 1:N)
    @test_throws ArgumentError lynden_bell_flip!(copy(vel), pos, [0.0, 0.0, 1.0], 1.5)
    # Null (p = ½): the spin parameter sits at its noise floor; maximal rotation
    # (p = 1) stands far above it and the shell rotation equals ⟨|v_φ|⟩.
    vh = copy(vel)
    lynden_bell_flip!(vh, pos, [0.0, 0.0, 1.0], 0.5; rng = StableRNG(5))
    null = rotation_analysis(pos, vh, m; n_min = 1000, n_null = 20, seed = 3)
    @test null.lambda_peebles < 3 * null.lambda_peebles_floor_analytic
    @test isapprox(null.lambda_peebles_floor, null.lambda_peebles_floor_analytic; rtol = 0.5)
    rot = rotation_analysis(pos, v1, m; n_min = 1000, n_null = 20, seed = 3)
    @test rot.lambda_peebles > 10 * rot.lambda_peebles_floor
    @test rot.lambda_r > 5 * rot.lambda_r_floor
    @test abs(rot.axis[3]) > 0.99
    for b in eachindex(rot.profile.radius)
        @test rot.profile.v_rot[b] > 0
        @test isapprox(rot.profile.v_rot_debiased[b], rot.profile.v_rot[b]; rtol = 0.05)
        @test rot.profile.v_rot[b] > 3 * rot.profile.v_rot_floor[b]
    end
    @test rot.v_peak > 0 && 0 < rot.r_peak < 5
end

@testset "Static models — shape" begin
    rng = StableRNG(13)
    N = 40_000
    pos, _ = sample_plummer(N, 1.0; rng)
    c = zeros(3)
    sph = shape_tensor(pos, c; fractions = (0.1, 0.5, 0.9))
    @test length(sph) == 3
    @test isnan(sph[1].q)                      # 4000 points inside r_10: fewer than 7000 enclosed
    for k in 2:3
        @test sph[k].converged && sph[k].n ≥ 3000
        @test sph[k].q > 0.95 && sph[k].s > 0.93
    end
    ε0 = projected_ellipticity(pos, c, _PLUMMER_RHM_OVER_A)
    @test size(ε0) == (3, 3) && all(ε0 .< 0.05)
    # y scaled by 0.6, z by 0.8: axes a = x, b = z, c = y, so q = b/a = 0.8, s = c/a = 0.6
    stretched = stretch_positions!(copy(pos), 0.6, 0.8)
    sh = shape_tensor(stretched, c; fractions = (0.1, 0.5, 0.9))
    for k in 2:3
        @test sh[k].converged
        @test isapprox(sh[k].q, 0.8; atol = 0.05)
        @test isapprox(sh[k].s, 0.6; atol = 0.05)
        @test abs(sh[k].axes[1, 1]) > 0.98          # major axis along x
        @test abs(sh[k].axes[3, 2]) > 0.98          # intermediate axis along z
    end
    @test_throws ArgumentError stretch_positions!(copy(pos), 0.0, 1.0)
    # Projected ellipticity: a circular aperture cuts the major axis first, so the
    # value at R_h sits below the axis ratio; a wide aperture approaches it.
    ε = projected_ellipticity(stretched, c, _PLUMMER_RHM_OVER_A; apertures = (0.25, 0.5, 1.0, 4.0))
    @test size(ε) == (3, 4)
    @test 0.2 < ε[1, 3] ≤ 0.45 && ε[1, 3] > ε[2, 3]           # (x, y) at R_h: 1 − 0.6 at most
    @test isapprox(ε[1, 4], 1 - 0.6; atol = 0.07)             # (x, y) at 4 R_h
    @test isapprox(ε[2, 4], 1 - 0.8; atol = 0.07)             # (x, z)
    @test isapprox(ε[3, 4], 1 - 0.6 / 0.8; atol = 0.07)       # (y, z): axes 0.6 and 0.8
end

@testset "Static models — centre, core radius, Lagrangian grid" begin
    rng = StableRNG(14)
    N = 50_000
    a = 1.0
    pos, vel = sample_plummer(N, a; rng)
    m = fill(1.0 / N, N)
    cen = shrinking_sphere_centre(pos, m)
    @test cen.n_final == 500
    # The fixed point of the iteration wanders inside the flat Plummer core: the
    # centroid noise σ is amplified by (a/R_f)² (verification log §7.3), so the
    # test asks for the centre inside the final sphere.
    @test maximum(abs, cen.centre) < 0.5 * cen.radius
    @test cen.radius > 0 && cen.sigma ≈ cen.radius / sqrt(5 * 500)
    cm = shrinking_sphere_centre(pos, m; weights = :mass)
    @test maximum(abs, cm.centre .- cen.centre) < 0.1
    small = shrinking_sphere_centre(pos[:, 1:200], m[1:200])
    @test small.n_final == 100
    @test_throws ArgumentError shrinking_sphere_centre(pos, m; weights = :volume)
    # Two Plummer spheres: per-tag centres at ±d/2, the global finder on one peak
    for d in (0.5, 1.0, 2.0, 5.0) .* _PLUMMER_RHM_OVER_A
        p2, _ = sample_plummer(N, a; rng = StableRNG(round(Int, 100d)))
        pa = pos .+ [d / 2, 0.0, 0.0]
        pb = p2 .- [d / 2, 0.0, 0.0]
        ca = shrinking_sphere_centre(pa, m).centre
        cb = shrinking_sphere_centre(pb, m).centre
        @test isapprox(ca[1], d / 2; atol = 0.15) && isapprox(cb[1], -d / 2; atol = 0.15)
        both = shrinking_sphere_centre(hcat(pa, pb), vcat(m, m))
        d ≥ 2 * _PLUMMER_RHM_OVER_A &&
            @test min(abs(both.centre[1] - d / 2), abs(both.centre[1] + d / 2)) < 0.2
    end
    # Core radius: 0.548 a with exact densities; the sixth-neighbour estimate a few per cent below
    r = [sqrt(sum(abs2, @view pos[:, i])) for i in 1:N]
    ρ_exact = _plummer_rho.(r, a)
    @test isapprox(core_radius(pos, ρ_exact).r_core, sqrt(0.3) * a; rtol = 0.02)
    ρ_n = local_density(pos, m; weighting = :number)
    rc_n = core_radius(pos, ρ_n).r_core
    @test 0.85 * sqrt(0.3) * a < rc_n < 1.02 * sqrt(0.3) * a
    @test isapprox(
        local_density(pos[:, 1:500], m[1:500]),
        _local_density(pos[:, 1:500], m[1:500]);
        rtol = 1e-12,
    )
    @test isapprox(
        local_density(pos[:, 1:500], m[1:500]; weighting = :number),
        _local_density(pos[:, 1:500], ones(500) ./ N; weighting = :number);
        rtol = 1e-12,
    )
    # Lagrangian grid: the Plummer closed form within the noise, masks per M2
    g = lagrangian_grid(pos, m, zeros(3); weights = :number, n_min = 100)
    @test g.masked[1] && g.masked[end]                 # f = 0.001 (50 < 100) and f = 1
    @test !g.masked[findfirst(==(0.5), g.fractions)]
    for (k, f) in enumerate(g.fractions)
        g.masked[k] && continue
        @test isapprox(log(g.radii[k]), log(_plummer_rf(f, a)); atol = 4 * g.sigma_log[k] + 0.01)
    end
    gm = lagrangian_grid(pos, m, zeros(3); weights = :mass)
    @test gm.n_eff ≈ N
    @test all(isnan, g.radii[g.masked]) && all(isfinite, g.radii[.!g.masked])
    # Relaxation and crossing times from the closed forms
    rt = relaxation_times(pos, vel, m, zeros(3))
    @test isapprox(rt.r_h, _PLUMMER_RHM_OVER_A * a; rtol = 0.03)
    @test rt.t_rh[1] ≈ 0.138 * sqrt(N) * rt.r_h^1.5 / (sqrt(1 / N) * log(0.02 * N))
    @test rt.t_rh[2] ≈ 0.138 * sqrt(N) * rt.r_h^1.5 / (sqrt(1 / N) * log(0.11 * N))
    @test rt.t_cr ≈ 2 * rt.r_h / rt.sigma_3d
    @test isapprox(rt.sigma_3d, sqrt(3π / 32); rtol = 0.03)   # Plummer ⟨v²⟩ = 3π/32 GM/a
    @test local_relaxation_time(1.0, 1.0, 0.1, 1000, 0.11) ≈ 0.338 / (0.1 * log(110.0))
    @test isnan(local_relaxation_time(0.0, 1.0, 0.1, 1000, 0.11))
end

@testset "Static models — segregation and mixing" begin
    rng = StableRNG(15)
    N = 20_000
    a = 1.0
    pos, vel = sample_plummer(N, a; rng)
    masses = sample_kroupa(N; m_low = 0.08, m_up = 100.0, rng = StableRNG(2))
    m = masses ./ sum(masses)
    # Independent masses: every ratio consistent with 1
    tracer = sortperm(m; rev = true)[1:100]
    seg0 = mass_segregation(pos, m; tracer, seed = 1)
    @test seg0.n_mst_used == [10, 20, 50, 100]
    @test abs(seg0.lambda_msr - 1) < 3 * seg0.lambda_err + 0.1
    @test abs(seg0.lambda_median - 1) < 0.3
    @test 1 / seg0.gamma_err^3 < seg0.gamma_msr < seg0.gamma_err^3 || abs(seg0.gamma_msr - 1) < 0.2
    # Energy-ranked masses: Λ grows monotonically with the degree S
    e = [0.5 * sum(abs2, @view vel[:, i]) for i in 1:N] .+ pair_potential(pos, m)
    λ = Float64[]
    for S in (0.0, 0.5, 1.0)
        mS = rank_masses_by_energy!(copy(m), e, S; rng = StableRNG(3))
        @test isapprox(sum(mS), 1.0; rtol = 1e-12) && sort(mS) == sort(m)
        tr = sortperm(mS; rev = true)[1:100]
        push!(λ, mass_segregation(pos, mS; tracer = tr, seed = 1).lambda_msr)
    end
    @test λ[1] < λ[2] < λ[3] && λ[3] > 1.5
    @test_throws ArgumentError rank_masses_by_energy!(copy(m), e, 1.5)
    # A⁺ between two Plummer populations of scale a₁ = a/2 and a₂ = a, by quadrature
    pos1, _ = sample_plummer(N, 0.5a; rng = StableRNG(21))
    pos2, _ = sample_plummer(N, a; rng = StableRNG(22))
    r1 = [sqrt(sum(abs2, @view pos1[:, i])) for i in 1:N]
    r2 = [sqrt(sum(abs2, @view pos2[:, i])) for i in 1:N]
    r_h = _PLUMMER_RHM_OVER_A * a
    ap = aplus(r1, r2, r_h)
    φ(r, s) = _plummer_mf(r, s) / _plummer_mf(r_h, s)      # cumulative fraction inside r_h
    xs = range(log10(1e-3), 0.0; length = 20001)
    quad = sum((φ(10^x * r_h, 0.5a) - φ(10^x * r_h, a)) for x in xs) * step(xs)
    @test ap.aplus > 0 && isapprox(ap.aplus, quad; atol = 4 * ap.error + 0.01)
    @test aplus(r2, r2, r_h).aplus == 0.0
    # Mixing: random tags give A₂⁺ ≈ 0 and D at the null floor; concentric populations do not
    tags = [i ≤ N ÷ 2 ? 1 : 2 for i in 1:N]
    mx = mixing(pos, vel, m, tags, zeros(3), r_h)
    @test abs(mx.aplus2) < 4 * mx.aplus2_error + 0.01
    @test mx.d_radius < 2.5 * mx.floor && mx.d_angular_momentum < 2.5 * mx.floor
    @test isnan(mx.d_energy)
    two = hcat(pos1, pos2)
    mx2 = mixing(two, hcat(vel, vel), vcat(m, m), vcat(fill(1, N), fill(2, N)), zeros(3), r_h)
    @test mx2.aplus2 > 0.1 && mx2.d_radius > 5 * mx2.floor
    @test_throws ArgumentError mixing(pos, vel, m, ones(Int, N), zeros(3), r_h)
    d = two_sample_d(r1, r1)
    @test d.d == 0.0 && d.floor ≈ 0.87 * sqrt(2N / N^2)
end

@testset "Static models — snapshot binaries" begin
    rng = StableRNG(16)
    N = 4000
    pos, vel = sample_plummer(N, 1.0; rng)
    m = fill(1.0 / N, N)
    # Inject 20 pairs: replace stars 1..40 by binaries at the hosts' positions
    a_b, e_b = 2.0e-4, 0.5    # E_b = m²/(2a) = 1.6×10⁻⁴ against kT = 10⁻⁴: hard
    inj = Tuple{Int,Int}[]
    for k in 1:20
        i, j = 2k - 1, 2k
        μ = m[i] + m[j]
        # Relative orbit at apocentre in the x–y plane, c.m. at the host position with the host velocity
        r_rel = a_b * (1 + e_b)
        v_rel = sqrt(μ * (2 / r_rel - 1 / a_b))
        cx, cv = pos[:, i], vel[:, i]
        pos[:, i] = cx .+ [r_rel * m[j] / μ, 0.0, 0.0]
        pos[:, j] = cx .- [r_rel * m[i] / μ, 0.0, 0.0]
        vel[:, i] = cv .+ [0.0, v_rel * m[j] / μ, 0.0]
        vel[:, j] = cv .- [0.0, v_rel * m[i] / μ, 0.0]
        push!(inj, (i, j))
    end
    pairs = snapshot_pairs(pos, vel, m; kT = 1.0e-4)
    @test length(pairs) ≥ 20
    found = Dict((min(p.i, p.j), max(p.i, p.j)) => p for p in pairs)
    for (i, j) in inj
        @test haskey(found, (i, j))
        p = found[(i, j)]
        @test isapprox(p.a, a_b; rtol = 1e-6) && isapprox(p.e, e_b; rtol = 1e-5)
        @test p.tidal_ok && p.hardness > 1 && p.binding_energy ≈ m[i] * m[j] / (2a_b)
    end
    @test issorted([p.binding_energy for p in pairs]; rev = true)
    @test all(isnan(p.hardness) for p in snapshot_pairs(pos, vel, m))
    # Reconciliation: the injected pairs as the engine's KS list
    snap = _sm_snapshot(1:N, m, pos, vel, 0.0)
    ks_snap =
        Snapshot(snap.header, snap.name, snap.mass, snap.pos, snap.vel, snap.rho, snap.phi, 20)
    rec = reconcile_pairs(pairs, ks_snap)
    @test rec.n_ks == 20 && rec.n_both == 20 && rec.n_ks_only == 0
    @test rec.n_snapshot_only == length(pairs) - 20 && all(rec.ks_found)
end

@testset "Static models — coalescence and core-collapse markers" begin
    rng = StableRNG(17)
    N = 3000
    pos, vel = sample_plummer(N, 1.0; rng)
    p2, v2 = sample_plummer(N, 1.0; rng = StableRNG(18))
    m = fill(1.0 / N, N)
    r_h = _PLUMMER_RHM_OVER_A
    names = 1:(2N)
    ranges = [1:N, (N + 1):(2N)]
    function two(d, t)
        pa = pos .+ [d / 2, 0.0, 0.0]
        pb = p2 .- [d / 2, 0.0, 0.0]
        return _sm_snapshot(names, vcat(m, m), hcat(pa, pb), hcat(vel, v2), t)
    end
    # Separated (d = 5 r_h), overlapping (d = 0), separated again, then overlapping for good
    snaps = [two(5r_h, 0.0), two(0.0, 1.0), two(5r_h, 2.0), two(0.0, 3.0), two(0.0, 4.0)]
    co = coalescence_time(snaps, ranges)
    @test co.n_distinct == [2, 1, 2, 1, 1]
    @test co.index == 4 && co.time_nb == 3.0
    @test size(co.separation) == (1, 5) && size(co.threshold) == (1, 5)
    @test isapprox(co.separation[1, 1], 5r_h; rtol = 0.1) &&
          co.separation[1, 2] < co.threshold[1, 2]
    never = coalescence_time([two(5r_h, 0.0), two(5r_h, 1.0)], ranges)
    @test never.index === nothing && isnan(never.time_nb)
    mk = core_collapse_markers(
        [0.0, 1.0, 2.0, 3.0],
        [1.0, 0.6, 0.5, 0.9],
        [1.0, 1.5, 2.5, 1.2];
        window = 1,
    )
    @test mk.i_rc_min == 3 && mk.i_rho_max == 3 && mk.t_rc_min == 2.0
    mk3 = core_collapse_markers(
        collect(0.0:6.0),
        [1.0, 1.0, 0.2, 1.0, 0.9, 0.95, 1.0],
        ones(7);
        window = 3,
    )
    @test mk3.i_rc_min == 4                                      # the running mean moves the minimum
    @test_throws ArgumentError core_collapse_markers([0.0, 1.0], [1.0, 1.0], [1.0, 1.0]; window = 2)
end
