# Equipartition-mass fit (Bianchini et al. 2016) on synthetic dispersion–mass relations and the
# per-snapshot series on a static model (uses `_am_snapshot` of test_angular_momentum.jl).

using StableRNGs: StableRNG
using LinearAlgebra: ×

_eq_g(m, m_eq) = m <= m_eq ? exp(-m / (2 * m_eq)) : exp(-0.5) * sqrt(m_eq / m)

@testset "equipartition_fit recovers a synthetic m_eq" begin
    rng = StableRNG(1101)
    n = 40_000
    mass = 0.1 .+ 1.9 .* rand(rng, n)
    vel = 10 .* _eq_g.(mass, 0.8)' .* randn(rng, 3, n)
    fit = equipartition_fit(mass, vel)
    @test fit isa EquipartitionFit
    @test fit.n_bins == 19 && fit.n_stars == n
    @test all(fit.bin_n .>= 50) && issorted(fit.bin_mass)
    @test 0.75 < fit.m_eq < 0.85
    @test fit.m_eq_low <= fit.m_eq <= fit.m_eq_high < 1.0
    # The Δχ² = 1 interval covers the truth within one step of the 3.3 % grid
    @test fit.m_eq_low / 1.04 < 0.8 < fit.m_eq_high * 1.04
    @test abs(fit.sigma0 - 10) < 0.3
    @test fit.chi2 < fit.chi2_flat
    @test fit.chi2 < 3 * fit.n_bins
    # The bin dispersions follow the model at the bin masses within their errors
    pull = (fit.bin_sigma .- fit.sigma0 .* _eq_g.(fit.bin_mass, fit.m_eq)) ./ fit.bin_sigma_err
    @test maximum(abs, pull) < 4
end

@testset "equipartition_fit without a mass dependence" begin
    rng = StableRNG(1102)
    n = 40_000
    mass = 0.1 .+ 1.9 .* rand(rng, n)
    vel = 10 .* randn(rng, 3, n)
    flat = equipartition_fit(mass, vel)
    @test flat.m_eq > 10 * maximum(mass)
    @test isinf(flat.m_eq_high)
    @test flat.chi2_flat - flat.chi2 < 4
    @test abs(flat.sigma0 - 10) < 0.2
end

@testset "equipartition_fit edge cases and errors" begin
    rng = StableRNG(1103)
    mass = 0.1 .+ 1.9 .* rand(rng, 120)
    vel = randn(rng, 3, 120)
    few = equipartition_fit(mass, vel)
    @test few.n_bins < 3 && isnan(few.m_eq) && isnan(few.sigma0) && isnan(few.chi2)
    @test few.n_bins == 0 || !isnan(few.chi2_flat) || few.n_bins < 2
    two = equipartition_fit(mass, vel; n_min = 2, bin_width = 1.0)
    @test two.n_bins == 2 && isnan(two.m_eq) && !isnan(two.chi2_flat)
    @test_throws DimensionMismatch equipartition_fit(mass, vel[:, 1:100])
    @test_throws DimensionMismatch equipartition_fit(mass, vel[1:2, :])
    @test_throws ArgumentError equipartition_fit(mass, vel; bin_width = 0.0)
    @test_throws ArgumentError equipartition_fit(mass, vel; n_min = 1)
    @test_throws ArgumentError equipartition_fit(mass, vel; n_grid = 1)
    @test_throws ArgumentError equipartition_fit(mass, vel; m_eq_max_factor = 1.0)
    @test_throws ArgumentError equipartition_fit(Float64[], zeros(3, 0))
    bad = copy(mass)
    bad[1] = 0.0
    @test_throws ArgumentError equipartition_fit(bad, vel)
end

@testset "write_equipartition" begin
    rows = [
        EquipartitionRow(0.0, 0.0, 4000, 19, 0.8, 0.7, 0.9, 10.0, 15.0, 400.0, 0.1, 2.0),
        EquipartitionRow(1.0, 1.5, 3900, 18, Inf, 20.0, Inf, 9.9, 17.0, 17.5, 0.1, 1.9),
        EquipartitionRow(2.0, 3.0, 100, 2, NaN, NaN, NaN, NaN, NaN, 1.0, 0.1, 1.2),
    ]
    mktempdir() do dir
        path = joinpath(dir, "equipartition.csv")
        @test write_equipartition(path, rows) == path
        lines = readlines(path)
        @test length(lines) == 4
        @test lines[1] ==
              "time_nb,time_myr,n_stars,n_bins,m_eq,m_eq_low,m_eq_high,sigma0_kms,chi2,chi2_flat,mass_min,mass_max"
        c2 = split(lines[3], ",")
        @test c2[5] == "Inf" && c2[7] == "Inf" && c2[3] == "3900"
        c3 = split(lines[4], ",")
        @test c3[5] == "" && c3[6] == "" && c3[7] == "" && c3[10] == "1"
        write_equipartition(path, rows[1:1])
        @test length(readlines(path)) == 2 && any(f -> occursin("#1", f), readdir(dir))
    end
end

@testset "equipartition_series on a static model" begin
    rng = StableRNG(1104)
    N = 20_000
    pos, vel = Nbody6Dynamics.sample_plummer(N, 1.0; rng)
    m_solar = 0.1 .+ 1.9 .* rand(rng, N)
    M_solar = sum(m_solar)
    mass = m_solar ./ M_solar
    # Velocities in km/s follow σ(m) = 10 g(m; 0.8); the NB velocities are those over vstar
    vstar = 2.5
    vel_kms = 10 .* _eq_g.(m_solar, 0.8)' .* randn(rng, 3, N)
    units = UnitScaling(1.0, M_solar, 1.0, vstar)
    snaps = [
        _am_snapshot(mass, pos, vel_kms ./ vstar; t = 0.0),
        _am_snapshot(mass, pos, vel_kms ./ vstar; t = 1.0),
    ]
    bsets = [
        BoundSet(collect(1:N), Int[], 0, :converged, (0.0, 0.0, 0.0), 0.0, 0.0, Int[]) for
        _ in snaps
    ]
    rows = equipartition_series(snaps, bsets, units)
    @test length(rows) == 2 && rows[2].time_nb == 1.0
    r = rows[1]
    @test 0.6 < r.m_eq < 1.0 && r.m_eq_low <= r.m_eq <= r.m_eq_high
    @test 0.1 <= r.mass_min < 0.2 && r.mass_max > 1.9
    @test r.n_stars > 0.4 * N && r.n_bins >= 15
    # A stellar-evolution snapshot at the same time restricts the fit to K* ≤ 9
    records = [
        StellarRecord(
            0.0,
            Int32(i),
            Int32(i),
            Int32(i <= 100 ? 14 : 1),
            1.0,
            m_solar[i],
            0.0,
            0.0,
            3.7,
        ) for i in 1:N
    ]
    sev = StellarEvolutionSnapshot(0.0, N, records)
    lum = equipartition_series(snaps, bsets, units; sevs = [sev])
    @test lum[1].n_stars < r.n_stars && lum[1].n_stars >= r.n_stars - 100
    @test isnan(lum[2].m_eq) == false   # the second snapshot has no matching sev, all stars enter
    small = equipartition_series(
        snaps,
        [
            BoundSet(collect(1:50), Int[], 0, :converged, (0.0, 0.0, 0.0), 0.0, 0.0, Int[]) for
            _ in snaps
        ],
        units,
    )
    @test all(x -> isnan(x.m_eq) && x.n_stars == 0, small)
    @test_throws ArgumentError equipartition_series(snaps, bsets[1:1], units)
end
