# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Gravitational-wave inspiral time and frequencies of a binary (Peters 1964
# orbit average), against the closed forms, the published limits, the Mandel
# (2021) fit and an independent integration of the Peters equations.

@testset "GW inspiral time and frequencies" begin
    nd = Nbody6Dynamics

    @testset "Peters–Mathews enhancement and the time factor" begin
        @test gw_enhancement(0.0) == 1.0
        @test gw_enhancement(0.5) ≈ (1 + 73 / 24 * 0.25 + 37 / 96 * 0.0625) / 0.75^3.5
        # F(e₀) = T / [T_c (1 − e₀²)^{7/2}]
        @test gw_time_factor(0.0) ≈ 1.0 rtol = 1e-12
        for (e, F) in (
            (0.5, 0.979318),
            (0.9, 1.139561),
            (0.99, 1.500922),
            (0.999, 1.696682),
            (0.9999, 1.770611),
        )
            @test gw_time_factor(e) ≈ F rtol = 1e-6
        end
        @test gw_time_factor(1 - 1e-14) ≈ 768 / 425 rtol = 1e-5
        # The two quadratures agree where both apply
        for e in (0.2, 0.4, 0.5, 0.6, 0.8)
            @test nd._gw_time_factor_x(e) ≈ nd._gw_time_factor_tau(e) rtol = 1e-8
        end
        # Low-e limit: F → [(1 + k e²)^{−870/2299}]⁴ (1 − e²)^{−7/2} (1 − e²)⁴, within 1 % at e = 0.1
        k = 121 / 304
        @test gw_time_factor(0.1) ≈ (1 - 0.01)^0.5 * (1 + k * 0.01)^(-3480 / 2299) rtol = 0.01
        @test_throws ArgumentError gw_time_factor(1.0)
        @test_throws ArgumentError gw_time_factor(-0.1)
    end

    @testset "Inspiral time" begin
        # Circular: T_c = 5 c⁵ a⁴ / (256 G³ m₁ m₂ M)
        G, c, au, yr = nd._GW_G, nd._GW_C, nd._GW_AU, nd._GW_YR
        msun = nd._GW_MSUN
        Tc = 5 * c^5 * au^4 / (256 * G^3 * (30msun)^2 * 60msun) / yr
        @test gw_inspiral_time(30.0, 30.0, 1.0, 0.0) ≈ Tc rtol = 1e-12
        # 30 + 30 M☉, 1 au, e = 0.999: ≈ 3.5×10³ yr (Mandel fit), the exact value within the fit's 3 %
        T = gw_inspiral_time(30.0, 30.0, 1.0, 0.999)
        @test 3.0e3 < T < 4.0e3
        @test gw_inspiral_time_mandel(30.0, 30.0, 1.0, 0.999) ≈ T rtol = 0.03
        for e in (0.0, 0.3, 0.6, 0.9, 0.99, 0.9999)
            @test gw_inspiral_time_mandel(10.0, 5.0, 0.1, e) ≈ gw_inspiral_time(10.0, 5.0, 0.1, e) rtol =
                0.03
        end
        # High-e limit (768/425) T_c (1 − e²)^{7/2} is 2 % high at 1 − e = 10⁻⁴
        e = 0.9999
        @test 768 / 425 * (1 - e^2)^3.5 * gw_inspiral_time(30.0, 30.0, 1.0, 0.0) /
              gw_inspiral_time(30.0, 30.0, 1.0, e) ≈ 1.02 atol = 0.01
        # Independent integration of the Peters equations to a(1 − e) ≪ a₀, plus the circular remainder
        for e0 in (0.5, 0.9)
            @test nd._gw_inspiral_time_ode(20.0, 10.0, 0.05, e0) ≈
                  gw_inspiral_time(20.0, 10.0, 0.05, e0) rtol = 1e-5
        end
        # The 1PN-recalibrated factor of Zwick et al. 2020: exp(2.5 r_S / p₀), p₀ = a(1 − e)
        rs = 2G * 60msun / c^2
        @test gw_zwick_factor(30.0, 30.0, 1.0, 0.999) ≈ exp(2.5rs / (au * 1.0e-3)) rtol = 1e-12
        @test_throws ArgumentError gw_inspiral_time(-1.0, 30.0, 1.0, 0.5)
        @test_throws ArgumentError gw_inspiral_time(30.0, 30.0, 0.0, 0.5)
    end

    @testset "Masses and frequencies" begin
        @test chirp_mass(30.0, 30.0) ≈ 30.0^1.2 / 60.0^0.2
        @test chirp_mass(30.0, 30.0) ≈ 26.12 rtol = 1e-3
        # 30 + 30 M☉, 1 au, e = 0.999: f_Wen = 12.6 mHz, f_peri = 15.5 mHz, f₂₂ = 0.49 μHz
        @test gw_frequency(30.0, 30.0, 1.0, 0.999; definition = :wen) ≈ 12.6e-3 rtol = 0.01
        @test gw_frequency(30.0, 30.0, 1.0, 0.999; definition = :peri) ≈ 15.5e-3 rtol = 0.01
        @test gw_frequency(30.0, 30.0, 1.0, 0.999; definition = :f22) ≈ 0.49e-6 rtol = 0.01
        @test gw_frequency(30.0, 30.0, 1.0, 0.0; definition = :f22) ≈
              2 * orbital_frequency(30.0, 30.0, 1.0)
        @test_throws ArgumentError gw_frequency(30.0, 30.0, 1.0, 0.5; definition = :other)
    end

    @testset "Eccentricity at a reference frequency" begin
        r = eccentricity_at_frequency(30.0, 30.0, 1.0, 0.999, 10.0; definition = :wen)
        @test r.e ≈ 4.42e-4 rtol = 0.01
        @test !r.above_at_record
        for d in (:f22, :peri)
            @test eccentricity_at_frequency(30.0, 30.0, 1.0, 0.999, 10.0; definition = d).e ≈
                  4.42e-4 rtol = 0.01
        end
        @test eccentricity_at_frequency(30.0, 30.0, 1.0, 0.999, 10.0e-3; definition = :f22).e ≈
              0.441 rtol = 0.01
        above = eccentricity_at_frequency(30.0, 30.0, 1.0, 0.999, 10.0e-3; definition = :wen)
        @test isnan(above.e) && above.above_at_record
        for (d, e) in ((:wen, 5.37e-3), (:f22, 5.34e-3), (:peri, 5.38e-3))
            @test eccentricity_at_frequency(30.0, 30.0, 0.2, 0.9, 10.0e-3; definition = d).e ≈ e rtol =
                0.01
        end
        r2 = eccentricity_at_frequency(30.0, 30.0, 0.2, 0.9, 10.0; definition = :wen)
        @test r2.e ≈ 3.6e-6 rtol = 0.05
        # The time to the reference frequency is the difference of the two inspiral times
        @test r2.time_to_ref_yr ≈
              gw_inspiral_time(30.0, 30.0, 0.2, 0.9) - gw_inspiral_time(30.0, 30.0, r2.a_au, r2.e) rtol =
            1e-8
    end

    @testset "Strain amplitudes" begin
        # Circular orbit: only n = 2 radiates
        @test peters_mathews_g(2, 0.0) ≈ 1.0 atol = 1e-12
        @test peters_mathews_g(1, 0.0) ≈ 0.0 atol = 1e-12
        @test peters_mathews_g(3, 0.0) ≈ 0.0 atol = 1e-12
        # The harmonics sum to the Peters–Mathews enhancement
        for (e, N, tol) in ((0.3, 400, 1e-8), (0.6, 400, 1e-8), (0.9, 4000, 1e-5))
            @test sum(n -> peters_mathews_g(n, e), 1:N) ≈ gw_enhancement(e) rtol = tol
        end

        # LDC amplitude of 𝓜 = 0.5 M☉ at f = 1 mHz and 1 kpc
        @test gw_amplitude_ldc(0.5, 1e-3) ≈ 1.87e-22 rtol = 1e-2

        m1 = m2 = 0.6
        a = 0.01
        f_orb = orbital_frequency(m1, m2, a)
        𝓜 = chirp_mass(m1, m2)
        @test gw_harmonic_strain(m1, m2, a, 0.0, 2) ≈ 4 / sqrt(5) * gw_amplitude_ldc(𝓜, 2f_orb) rtol =
            1e-12

        # h_{c,n}² = (fₙ²/ḟₙ) hₙ²
        e = 0.4
        for n in (1, 2, 3, 5)
            f_n = n * f_orb
            fdot_n = gw_harmonic_frequency_derivative(m1, m2, a, e, n)
            h_n = gw_harmonic_strain(m1, m2, a, e, n)
            @test gw_characteristic_strain(m1, m2, a, e, n)^2 ≈ f_n^2 / fdot_n * h_n^2 rtol = 1e-8
        end
        # The n = 2 drift at e = 0 is the circular chirp of f₂₂; ḟ_orb = ḟ₂₂ F(e)/2 in general
        @test gw_harmonic_frequency_derivative(m1, m2, a, 0.0, 2) ≈
              nd._gwcat_fdot_circular(𝓜, 2f_orb) rtol = 1e-12
        @test gw_harmonic_frequency_derivative(m1, m2, a, e, 1) ≈
              nd._gwcat_fdot_circular(𝓜, 2f_orb) * gw_enhancement(e) / 2 rtol = 1e-12

        # Every strain scales as 1/D
        @test gw_amplitude_ldc(𝓜, 2f_orb; distance_kpc = 8.0) ≈ gw_amplitude_ldc(𝓜, 2f_orb) / 8 rtol =
            1e-12
        @test gw_harmonic_strain(m1, m2, a, e, 3; distance_kpc = 2.0) ≈
              gw_harmonic_strain(m1, m2, a, e, 3; distance_kpc = 8.0) * 4 rtol = 1e-12
        @test gw_characteristic_strain(m1, m2, a, e, 3; distance_kpc = 2.0) ≈
              gw_characteristic_strain(m1, m2, a, e, 3; distance_kpc = 8.0) * 4 rtol = 1e-12

        @test_throws ArgumentError peters_mathews_g(0, 0.5)
        @test_throws ArgumentError peters_mathews_g(2, 1.0)
        @test_throws ArgumentError peters_mathews_g(2, -0.1)
        @test_throws ArgumentError gw_amplitude_ldc(0.0, 1e-3)
        @test_throws ArgumentError gw_amplitude_ldc(0.5, 0.0)
        @test_throws ArgumentError gw_amplitude_ldc(0.5, 1e-3; distance_kpc = 0.0)
        @test_throws ArgumentError gw_harmonic_strain(-0.6, m2, a, e, 2)
        @test_throws ArgumentError gw_harmonic_strain(m1, m2, 0.0, e, 2)
        @test_throws ArgumentError gw_harmonic_strain(m1, m2, a, 1.0, 2)
        @test_throws ArgumentError gw_harmonic_strain(m1, m2, a, e, 0)
        @test_throws ArgumentError gw_harmonic_strain(m1, m2, a, e, 2; distance_kpc = -1.0)
        @test_throws ArgumentError gw_characteristic_strain(m1, m2, a, e, 2; distance_kpc = 0.0)
        @test_throws ArgumentError gw_characteristic_strain(m1, 0.0, a, e, 2)
        @test_throws ArgumentError gw_harmonic_frequency_derivative(m1, m2, a, e, 0)
        @test_throws ArgumentError gw_harmonic_frequency_derivative(m1, m2, -a, e, 1)
    end
end
