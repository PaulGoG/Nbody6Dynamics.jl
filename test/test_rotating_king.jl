# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Rotating King model of Einsel & Spurzem (1999, DOI 10.1046/j.1365-8711.1999.02083.x):
# local velocity moments against King (1966), the self-consistent Legendre–Poisson
# solution against the published W₀ = 6 sequence (E&S99 Table 1), exact equilibrium
# properties, and the phase-space sampler. Samples use StableRNGs so the reference
# values hold across Julia versions.

using StableRNGs: StableRNG
using LinearAlgebra: Symmetric, eigvals

@testset "Rotating King model" begin
    nd = Nbody6Dynamics
    xq, wq = nd._rk_gauss_legendre(64)

    @testset "Local moments" begin
        # k = 0: the King (1966) density up to a constant factor
        ratios = [nd._rk_moments(ψ, 0.0, xq, wq).ρ / nd._king_density(ψ) for ψ in (0.5, 2.0, 6.0)]
        @test all(q -> isapprox(q, ratios[1]; rtol = 1e-10), ratios)
        # k → 0: ⟨v_φ⟩ → k⟨u²⟩/3
        k = 1e-3
        mo = nd._rk_moments(3.0, k, xq, wq)
        @test isapprox(mo.ρv_φ / mo.ρ, k * (mo.ρu² / mo.ρ) / 3; rtol = 1e-6)
        # ⟨v_φ²⟩ + 2⟨v_R²⟩ = ⟨u²⟩
        for (ψ, k) in ((0.3, 0.0), (2.0, 0.7), (6.0, 4.0))
            mo = nd._rk_moments(ψ, k, xq, wq)
            @test isapprox((mo.ρv_φ² + 2 * mo.ρv_R²) / mo.ρ, mo.ρu² / mo.ρ; rtol = 1e-12)
        end
        # the 64-node rule against a 256-node rule
        x256, w256 = nd._rk_gauss_legendre(256)
        rng = StableRNG(11)
        for _ in 1:10
            ψ = 1e-3 + 6 * rand(rng)
            k = 5 * rand(rng)
            @test isapprox(
                nd._rk_density(ψ, k, xq, wq),
                nd._rk_density(ψ, k, x256, w256);
                rtol = 1e-6,
            )
        end
    end

    # Central one-dimensional dispersion of the W₀ = 6 model (identical for every ω₀,
    # k = 0 at the centre), for the core radius r_K σ_c/σ reported alongside r_K.
    mo_c = nd._rk_moments(6.0, 0.0, xq, wq)
    σ_c = sqrt(mo_c.ρu² / mo_c.ρ / 3)

    m0 = nd.solve_rotating_king(6.0, 0.0)
    @testset "Non-rotating limit" begin
        @test maximum(maximum(m0.ψ[i, :]) - minimum(m0.ψ[i, :]) for i in axes(m0.ψ, 1)) < 1e-8
        rk, wk, _ = nd._solve_king(6.0; n_grid = 20_000)
        r_tK = rk[end]
        function king_ψ(x)
            j = clamp(searchsortedlast(rk, x), 1, length(rk) - 1)
            f = (x - rk[j]) / (rk[j + 1] - rk[j])
            return (1 - f) * wk[j] + f * wk[j + 1]
        end
        inner = findall(<(0.9 * r_tK), m0.r)
        @test maximum(abs(m0.ψ[i, 1] - king_ψ(m0.r[i])) / king_ψ(m0.r[i]) for i in inner) < 1e-5
        mom = nd.rotating_king_moments(m0)
        @test mom.λ == 0
        @info(
            "Rotating King W₀ = 6, ω₀ = 0",
            r_h_over_r_K = mom.r_h_over_r_K,
            r_t_over_r_K = mom.r_t_over_r_K,
            r_t_King = r_tK,
            σ_c = mom.σ_c,
            r_h_over_r_c = mom.r_h_over_r_c,
            r_t_over_r_c = mom.r_t_over_r_c,
            iterations = m0.iterations,
        )
        # E&S99 Table 1 with their core radius r_c = r_K σ_c/σ (King 1966, Eq. 9):
        # r_t/r_c = 18.72 is met to 0.4 %, r_h/r_c = 2.70 to 1.6 %.
        @test isapprox(mom.σ_c, σ_c; rtol = 1e-12)
        @test isapprox(mom.r_h_over_r_c, 2.70; rtol = 0.02)
        @test isapprox(mom.r_t_over_r_c, 18.72; rtol = 0.01)
    end

    # Einsel & Spurzem (1999) Table 1, W₀ = 6
    ω_table = (0.2, 0.4, 0.6, 0.8, 1.0)
    λ_table = (0.0338, 0.1123, 0.1981, 0.2718, 0.3299)
    r_h_table = (2.66, 2.55, 2.39, 2.20, 2.04)
    models = [nd.solve_rotating_king(6.0, ω) for ω in ω_table]
    moments = [nd.rotating_king_moments(m) for m in models]
    @testset "Einsel–Spurzem W₀ = 6 sequence" begin
        for q in eachindex(ω_table)
            mom = moments[q]
            @info(
                "Rotating King W₀ = 6, ω₀ = $(ω_table[q])",
                λ = mom.λ,
                t_rot = mom.t_rot,
                e_dyn = mom.e_dyn,
                r_t_equator = mom.r_t_over_r_K,
                r_t_pole = models[q].r_t_pole,
                r_h_over_r_K = mom.r_h_over_r_K,
                r_h_over_r_c = mom.r_h_over_r_c,
                iterations = models[q].iterations,
            )
            # Measured against the table: λ within 0.3 % (ω₀ = 0.2), 1.6 % (0.6), 2.7 % (0.8),
            # 3.6 % (1.0) and r_h/r_c within 1.7 %, 0.9 %, 3.1 %, 6.9 %, a drift growing with the
            # rotation that is insensitive to l_max (12 → 20), n_μ (48 → 64) and n_r (600 → 1200);
            # the 1999 (r, z)-mesh values are taken as the less converged ones at strong rotation.
            tol_λ = ω_table[q] ≤ 0.6 ? 0.03 : 0.05
            tol_r = ω_table[q] ≤ 0.6 ? 0.03 : 0.08
            @test isapprox(mom.λ, λ_table[q]; rtol = tol_λ)
            @test isapprox(mom.r_h_over_r_c, r_h_table[q]; rtol = tol_r)
            @test models[q].r_t_equator > models[q].r_t_pole
        end
        @test all(>(0), diff([mom.e_dyn for mom in moments]))
    end

    m6 = models[3]
    mom6 = moments[3]
    c6 = nd._rk_legendre_coefficients(m6)
    @testset "Equilibrium properties at ω₀ = 0.6" begin
        # rotation in the equatorial plane: positive inside, zero at the surface ψ = 0
        for R in range(0.05, 0.95; length = 19) .* m6.r_t_equator
            ψ_eq = nd._rk_psi_at(m6, c6, R, 0.0)
            mo = nd._rk_moments(ψ_eq, m6.ω0 * R, xq, wq)
            @test ψ_eq > 0 && mo.ρv_φ / mo.ρ > 0
        end
        mo_s = nd._rk_moments(1e-8, m6.ω0 * m6.r_t_equator, xq, wq)
        @test 0 < mo_s.ρv_φ / mo_s.ρ < 1e-6
        # scalar virial theorem, rotation included in T_kin
        @test abs(2 * mom6.T_kin + mom6.W) / abs(mom6.W) < 1e-3
        # Poisson consistency of the converged solution: the non-spherical part of the
        # multipole reconstruction is the fixed point of the iteration (to the tolerance);
        # its spherical part differs from the exactly integrated one by the second-order
        # error of the product-trapezoid rules (≈ 3e-3 at n_r = 600).
        _, wμ = nd._rk_gauss_legendre(length(m6.μ))
        Φ = nd._rk_multipole_potential(m6.r, m6.μ, wμ, m6.ρ, m6.l_max)
        ψ_re = m6.W0 .- (Φ .- Φ[1, 1])
        mean_μ(a) = [0.5 * sum(wμ[j] * a[i, j] for j in axes(a, 2)) for i in axes(a, 1)]
        ψh_re = ψ_re .- mean_μ(ψ_re)
        ψh = m6.ψ .- mean_μ(m6.ψ)
        @test maximum(abs.(ψh_re .- ψh)) ≤ 1e-7
        @test maximum(abs.(mean_μ(ψ_re) .- mean_μ(m6.ψ))) ≤ 5e-3
    end

    @testset "Critical rotation" begin
        @test_throws ArgumentError nd.solve_rotating_king(6.0, 3.0)
    end

    @testset "Sampling at ω₀ = 0.6" begin
        N = 20_000
        pos, vel = nd.sample_rotating_king(N, m6; rng = StableRNG(1))
        @test size(pos) == (3, N) && size(vel) == (3, N)
        r = [sqrt(pos[1, p]^2 + pos[2, p]^2 + pos[3, p]^2) for p in 1:N]
        ψ_p = [nd._rk_psi_at(m6, c6, r[p], r[p] > 0 ? pos[3, p] / r[p] : 0.0) for p in 1:N]
        u² = [vel[1, p]^2 + vel[2, p]^2 + vel[3, p]^2 for p in 1:N]
        @test all(>(0), ψ_p)
        @test all(p -> u²[p] < 2 * ψ_p[p], 1:N)

        # λ from the mean v_φ in 20 cylindrical bins of equal counts
        R = [hypot(pos[1, p], pos[2, p]) for p in 1:N]
        L_z = [pos[1, p] * vel[2, p] - pos[2, p] * vel[1, p] for p in 1:N]
        v_φ = L_z ./ R
        order = sortperm(R)
        n_bin = N ÷ 20
        T_rot = sum(
            0.5 * n_bin * (sum(v_φ[order[((b - 1) * n_bin + 1):(b * n_bin)]]) / n_bin)^2 for
            b in 1:20
        )
        T_kin = 0.5 * sum(u²)
        λ_s = T_rot / T_kin
        @test sum(L_z) / N > 0

        # axis ratio c/a of the second-moment tensor against the continuum density
        S = Symmetric(pos * transpose(pos) ./ N)
        ev = eigvals(S)
        ca_sample = sqrt(ev[1] / ev[3])
        _, wμ = nd._rk_gauss_legendre(length(m6.μ))
        I_R² = sum(
            m6.ρ[i, j] * wμ[j] * m6.r[i]^4 * (1 - m6.μ[j]^2) for
            i in axes(m6.ρ, 1), j in axes(m6.ρ, 2)
        )
        I_z² = sum(
            m6.ρ[i, j] * wμ[j] * m6.r[i]^4 * m6.μ[j]^2 for i in axes(m6.ρ, 1), j in axes(m6.ρ, 2)
        )
        ca_continuum = sqrt(I_z² / (0.5 * I_R²))
        @info("Rotating King sample, ω₀ = 0.6", λ_s, λ_continuum = mom6.λ, ca_sample, ca_continuum)
        @test abs(λ_s - mom6.λ) < 0.01
        @test abs(ca_sample - ca_continuum) < 0.02

        # velocity draws at a fixed position against the quadrature moments
        rng = StableRNG(2)
        n = 200_000
        draws = [nd._rotating_king_velocity(3.0, 0.6, rng) for _ in 1:n]
        vφ_d = [d[2] for d in draws]
        u²_d = [d[1]^2 + d[2]^2 + d[3]^2 for d in draws]
        mo = nd._rk_moments(3.0, 0.6, xq, wq)
        mean_vφ = sum(vφ_d) / n
        mean_u² = sum(u²_d) / n
        se_vφ = sqrt((sum(abs2, vφ_d) / n - mean_vφ^2) / n)
        se_u² = sqrt((sum(abs2, u²_d) / n - mean_u²^2) / n)
        @test abs(mean_vφ - mo.ρv_φ / mo.ρ) < 3 * se_vφ
        @test abs(mean_u² - mo.ρu² / mo.ρ) < 3 * se_u²
    end

    @testset "Determinism" begin
        @test nd.sample_rotating_king(2000, m6; rng = StableRNG(3)) ==
              nd.sample_rotating_king(2000, m6; rng = StableRNG(3))
    end

    # m6 was solved with the default settings of the cache; the generator tests below
    # reuse it instead of solving the model again.
    get!(nd._ROTATING_KING_CACHE, (6.0, 0.6), m6)

    @testset "Rotating King profile: parsing and validation" begin
        table = Dict{String,Any}(
            "type" => "rotating_king",
            "W0" => 6.0,
            "omega0" => 0.6,
            "axis" => [0, 0, -1],
        )
        p = nd._parse_profile_table(table, "merger.cluster1")
        @test p isa nd.RotatingKingProfile
        @test p.W0 == 6.0 && p.omega0 == 0.6 && p.axis == [0.0, 0.0, -1.0]
        @test nd.profile_name(p) == "rotating_king"
        q = nd._parse_profile_table(nd._profile_table(p), "merger.cluster1")
        @test q isa nd.RotatingKingProfile
        @test (q.W0, q.omega0, q.axis) == (p.W0, p.omega0, p.axis)
        d = nd._parse_profile_table(Dict{String,Any}("type" => "rotating_king"), "merger.cluster1")
        @test (d.W0, d.omega0, d.axis) == (6.0, 0.6, [0.0, 0.0, 1.0])
        # omega0 and axis are not keys of the King and Plummer tables
        @test_throws ArgumentError nd._parse_profile_table(
            Dict{String,Any}("type" => "king", "omega0" => 0.6),
            "merger.cluster1",
        )
        @test_throws ArgumentError nd._parse_profile_table(
            Dict{String,Any}("type" => "plummer", "axis" => [0, 0, 1]),
            "merger.cluster1",
        )
        function validation_message(profile)
            try
                nd._validate_profile(profile, "merger.cluster1")
            catch err
                err isa ArgumentError || rethrow()
                return err.msg
            end
            return ""
        end
        @test validation_message(nd.RotatingKingProfile()) == ""
        @test occursin(
            "merger.cluster1.omega0",
            validation_message(nd.RotatingKingProfile(omega0 = -1.0)),
        )
        @test occursin(
            "merger.cluster1.axis",
            validation_message(nd.RotatingKingProfile(axis = [0.0, 0.0, 0.0])),
        )
        @test occursin(
            "merger.cluster1.axis",
            validation_message(nd.RotatingKingProfile(axis = [0.0, 1.0])),
        )
    end

    @testset "Rotating King profile: central density contrast" begin
        ρ̂_rot = nd._central_density_contrast(nd.RotatingKingProfile(W0 = 6.0, omega0 = 0.0))
        ρ̂_king = nd._central_density_contrast(nd.KingProfile(W0 = 6.0))
        @test isapprox(ρ̂_rot, ρ̂_king; rtol = 0.01)
        @test nd._rotating_king_cached(6.0, 0.0) === nd._rotating_king_cached(6.0, 0.0)
    end

    @testset "Rotation onto the profile axis" begin
        cross3(a, b) =
            [a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1]]
        for axis in ([1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [1.0, 1.0, 1.0] ./ sqrt(3), [0.0, 0.0, 1.0])
            pos = [1.0 0.3; 0.0 -0.2; 0.0 0.5]
            vel = [0.0 -0.4; 1.0 0.1; 0.0 0.7]
            norms_before = [sqrt(sum(abs2, M[:, p])) for M in (pos, vel) for p in 1:2]
            nd._rotate_to_axis!(pos, vel, axis)
            @test maximum(abs.(cross3(pos[:, 1], vel[:, 1]) .- axis)) < 1e-12
            norms_after = [sqrt(sum(abs2, M[:, p])) for M in (pos, vel) for p in 1:2]
            @test maximum(abs.(norms_after .- norms_before)) < 1e-12
        end
    end

    @testset "Rotating King cluster sample" begin
        λ_model = nd.rotating_king_moments(nd._rotating_king_cached(6.0, 0.6)).λ
        sample(axis) = nd._sample_cluster(
            nd.ClusterSpec(;
                N = 20_000,
                rbar = 1.0,
                profile = nd.RotatingKingProfile(W0 = 6.0, omega0 = 0.6, axis = axis),
                imf = nd.EqualMassIMF(particle_mass = 1.0),
            );
            rng = StableRNG(21),
            nmax = nd._VIRIAL_NMAX,
            zmet = nd.StellarSpec().zmet,
        )
        c = sample([0.0, 0.0, 1.0])
        rot = nd.rotation_moments(c.mass, c.pos, c.vel)
        @info("Rotating King cluster, ω₀ = 0.6", λ_realised = rot.λ, λ_model)
        @test abs(rot.λ - λ_model) < 0.02
        @test rot.l_axis > 0
        @test c.rotation.realised == rot
        @test c.rotation.moments.λ == λ_model
        # The radius is set about the model centre, before virialise! moves the
        # centre of mass to the origin; the shift is of order r_h/√N.
        @test isapprox(nd.half_mass_radius(c.mass, c.pos), 1.0; rtol = 0.01)

        c_retro = sample([0.0, 0.0, -1.0])
        @test nd.rotation_moments(c_retro.mass, c_retro.pos, c_retro.vel).l_axis < 0
        c_x = sample([1.0, 0.0, 0.0])
        l_x = nd.rotation_moments(c_x.mass, c_x.pos, c_x.vel; axis = [1.0, 0.0, 0.0]).l_axis
        l_z = nd.rotation_moments(c_x.mass, c_x.pos, c_x.vel).l_axis
        @test l_x > 0
        @test abs(l_z) < 0.03 * l_x
    end

    @testset "Rotating King profile refuses mass segregation" begin
        @test_throws ArgumentError nd.segregate_systems(
            nd.RotatingKingProfile(),
            fill(1.0, 100),
            1.0;
            degree = 0.5,
            rng = StableRNG(5),
        )
    end

    @testset "Rotating King cluster record" begin
        out = mktempdir()
        cfg = nd.MergerConfig(
            [
                nd.ClusterSpec(;
                    N = 2000,
                    rbar = 1.0,
                    profile = nd.RotatingKingProfile(W0 = 6.0, omega0 = 0.6),
                    imf = nd.EqualMassIMF(particle_mass = 1.0),
                ),
            ],
            "single",
            nd.OrbitSpec(),
            nd.MergerOutputSpec(format = "nbody", output_dir = out);
            seed = 17,
        )
        nd.generate_merger_ic(cfg; output_dir = out)
        meta = nd.TOML.parsefile(joinpath(out, "merger_ic.toml"))
        @test meta["cluster1"]["profile"]["type"] == "rotating_king"
        rot = meta["cluster1"]["rotation"]
        @test rot["omega0"] == 0.6
        @test rot["axis"] == [0.0, 0.0, 1.0]
        @test rot["lambda_model"] > 0
        @test isfinite(rot["lambda_realised"])
        @test rot["l_axis_realised_pc_kms"] > 0
        p = nd.load_merger_ic_result(out).cluster_specs[1].profile
        @test p isa nd.RotatingKingProfile
        @test (p.W0, p.omega0, p.axis) == (6.0, 0.6, [0.0, 0.0, 1.0])
    end
end
