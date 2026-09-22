# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Binary population diagnostics" begin
    # Synthetic snapshot: four singles and one pair with known velocities.
    # Scalings: 10 M☉ per NB mass unit, 2 km/s per NB velocity unit.
    params = zeros(Float32, 20)
    params[3] = 1.0f0
    params[4] = 10.0f0
    params[11] = 1.0f0
    params[12] = 2.0f0
    names = Int32.(1:6)
    mass = Float32[0.1, 0.1, 0.1, 0.1, 0.3, 0.3]
    pos = zeros(Float32, 3, 6)
    vel = zeros(Float32, 3, 6)
    vel[1, 1:4] .= Float32[1, -1, 1, -1]     # singles: ±1 along x
    vel[1, 5], vel[1, 6] = 3.0f0, -3.0f0     # pair (5, 6): equal masses, c.m. at rest in x
    vel[2, 5] = vel[2, 6] = 0.5f0            # c.m. moving at 0.5 along y
    hdr = SnapshotHeader(Int32(6), Int32(1), Int32(1), Int32(20), params)
    snap = Snapshot(hdr, names, mass, pos, vel, Float32[], Float32[])

    _bin_rec(n1, n2, m1, m2, loga, e) = BinaryRecord(
        0.0,
        Int32(n1),
        Int32(n2),
        Int32(n1),
        Int32(n2),
        Int32(0),
        Int32(0),
        Int32(0),
        0.0,
        e,
        0.0,
        loga,
        m1,
        m2,
        zeros(18)...,
    )
    hard_rec = _bin_rec(5, 6, 3.0, 3.0, 3.0, 0.1)   # a = 10³ R☉
    soft_rec = _bin_rec(5, 6, 3.0, 3.0, 6.0, 0.9)   # a = 10⁶ R☉
    bev0 = BinaryEvolutionSnapshot(0.0, 2, [hard_rec, soft_rec])
    bev1 = BinaryEvolutionSnapshot(1.0, 1, [hard_rec])

    # Systems: four singles (m = 0.1, v = ±x̂) and the pair's c.m.
    # (m = 0.6, v = 0.5 ŷ). M = 1, v_c = 0.3 ŷ,
    # σ² = [4 × 0.1 × (1 + 0.09) + 0.6 × 0.04] / (3 × 1) = 0.46/3.
    sc = hardness_scale(snap, bev0)
    @test sc.n_systems == 5
    @test sc.m_mean ≈ 2.0 rtol = 1e-6            # Float32 particle masses
    @test sc.sigma_kms ≈ 2 * sqrt(0.46 / 3) rtol = 1e-6
    # A pair with a component missing from the snapshot is not reduced.
    sc_missing =
        hardness_scale(snap, BinaryEvolutionSnapshot(0.0, 1, [_bin_rec(5, 99, 3.0, 3.0, 3.0, 0.0)]))
    @test sc_missing.n_systems == 6
    @test_throws ArgumentError hardness_scale(
        Snapshot(
            hdr,
            Int32[],
            Float32[],
            zeros(Float32, 3, 0),
            zeros(Float32, 3, 0),
            Float32[],
            Float32[],
        ),
        bev0,
    )

    x = binary_hardness(bev0, sc.m_mean, sc.sigma_kms)
    scale = sc.m_mean * sc.sigma_kms^2
    @test x[1] ≈ binding_energy(hard_rec) / scale
    @test x[1] > 1 && x[2] < 1
    @test x[1] / x[2] ≈ 1e3            # same masses, a ratio 10³
    @test semi_major_axis_pc(hard_rec) ≈ 1e3 * Nbody6Dynamics._RSUN_PC
    @test_throws ArgumentError binary_hardness(bev0, 0.0, 1.0)
    @test_throws ArgumentError binary_hardness(bev0, 1.0, -1.0)

    pop = binary_population([bev0, bev1]; n_stars = 6, m_mean = sc.m_mean, sigma_kms = sc.sigma_kms)
    @test pop.classified
    @test pop.time_myr == [0.0, 1.0]
    @test pop.n_pairs == [2, 1]
    @test pop.n_hard == [1, 1]
    @test pop.n_soft == [1, 0]
    @test pop.binary_fraction ≈ [2 / 4, 1 / 5]   # N_b / (N_stars − N_b)
    pop_u = binary_population([bev0, bev1])
    @test !pop_u.classified
    @test all(iszero, pop_u.n_hard) && pop_u.n_soft == pop_u.n_pairs
    @test all(isnan, pop_u.binary_fraction)
    pop_v = binary_population(
        [bev0, bev1];
        n_stars = [6, 4],
        m_mean = [2.0, 2.0],
        sigma_kms = [sc.sigma_kms, 1e-3],
    )
    @test pop_v.binary_fraction ≈ [0.5, 1 / 3]
    @test pop_v.n_hard == [1, 1]
    @test_throws DimensionMismatch binary_population([bev0, bev1]; n_stars = [6, 6, 6])
    @test_throws ArgumentError binary_population(BinaryEvolutionSnapshot[])

    scales = binary_scales([bev0, bev1], [snap])
    @test scales.n_stars == [6, 6]
    @test scales.m_mean ≈ [2.0, 2.0] rtol = 1e-6
    @test scales.sigma_kms ≈ fill(sc.sigma_kms, 2) rtol = 1e-6
    @test binary_scales([bev0], Snapshot[]) === nothing

    # Figures: time series, single epoch, classified and unclassified
    # elements, period histograms, and the no-data paths.
    vis_b = VisualizationConfig(; format = "png", column = "single", output_dir = mktempdir())
    @test isfile(plot_binary_population(pop, vis_b; filename = "pop_series"))
    @test isfile(
        plot_binary_population(
            binary_population([bev0]; n_stars = 6),
            vis_b;
            filename = "pop_single",
        ),
    )
    @test isfile(
        plot_binary_orbital_elements(
            bev0,
            vis_b;
            m_mean = sc.m_mean,
            sigma_kms = sc.sigma_kms,
            filename = "ae_classified",
        ),
    )
    @test isfile(plot_binary_orbital_elements(bev0, vis_b; filename = "ae_plain"))
    @test isfile(plot_binary_period_distribution([bev0, bev1], vis_b; filename = "period_two"))
    @test isfile(plot_binary_period_distribution([bev0], vis_b; filename = "period_one"))
    empty_bev = BinaryEvolutionSnapshot(2.0, 0, BinaryRecord[])
    @test isfile(plot_binary_orbital_elements(empty_bev, vis_b; filename = "ae_empty"))
    @test isfile(plot_binary_period_distribution([empty_bev], vis_b; filename = "period_empty"))
    @test_throws ErrorException plot_binary_period_distribution(BinaryEvolutionSnapshot[], vis_b)
    rm(vis_b.output_dir; recursive = true, force = true)
end

# =====================================================================

@testset "Remnant diagnostics" begin
    using StableRNGs
    rng_r = StableRNG(7)
    N = 2000
    pos, vel = sample_plummer(N, 1.0; rng = rng_r)
    mass = fill(1.0 / N, N)
    virialise!(mass, pos, vel)

    # Core radius: with exact Plummer densities the particle-weighted
    # Casertano–Hut estimator gives √0.3 a; the 6-neighbour estimate
    # sits a few per cent below.
    r_all = vec(sqrt.(sum(pos .^ 2; dims = 1)))
    ρ_exact = (1 .+ r_all .^ 2) .^ (-2.5)
    cr = core_radius(pos, ρ_exact)
    @test cr.r_core ≈ sqrt(0.3) rtol = 0.06
    @test all(abs.(cr.centre) .< 0.1)
    ρ6 = Nbody6Dynamics._local_density(pos, mass)
    @test length(ρ6) == N && all(≥(0), ρ6)
    @test 0.42 < core_radius(pos, ρ6).r_core < 0.62
    @test_throws ArgumentError Nbody6Dynamics._local_density(pos[:, 1:5], mass[1:5])
    @test_throws DimensionMismatch core_radius(pos, ρ6[1:10])
    @test_throws ArgumentError core_radius(pos, zeros(N))

    # Rotation: an isotropic sphere has λ_R at the noise level; adding
    # solid-body rotation about z recovers the axis and raises λ_R.
    rot0 = rotation_analysis(pos, vel, mass)
    @test rot0.lambda_r < 0.15
    @test rot0.lambda_peebles < 0.05
    @test length(rot0.profile.radius) == 8 && issorted(rot0.profile.radius)
    vrot = copy(vel)
    Ω = 0.5
    for i in 1:N
        vrot[1, i] -= Ω * pos[2, i]
        vrot[2, i] += Ω * pos[1, i]
    end
    rot1 = rotation_analysis(pos, vrot, mass)
    @test rot1.lambda_r > 0.6
    @test abs(rot1.axis[3]) > 0.99
    @test rot1.axis[3] > 0                          # L along +z for counter-clockwise motion
    @test all(>(0), rot1.profile.v_rot)              # v_φ > 0 in every shell
    @test rot1.profile.v_rot_over_sigma[end - 1] > rot1.profile.v_rot_over_sigma[1]   # Ω R grows outward
    @test rotation_analysis(pos, vrot, mass; nbins = 3).profile.radius |> length == 3
    @test_throws ArgumentError rotation_analysis(pos[:, 1:5], vrot[:, 1:5], mass[1:5])

    # Mass segregation: random masses give Λ ≈ 1; the heaviest stars
    # placed innermost give Λ ≫ 1 and a small half-mass ratio.
    m_rand = mass .* (1 .+ rand(rng_r, N))
    ms0 = mass_segregation(pos, m_rand; seed = 1)
    @test abs(ms0.lambda_msr - 1) < 3 * ms0.lambda_err + 0.3
    @test 0.7 < ms0.r_half_ratio < 1.3
    m_seg = zeros(N)
    m_seg[sortperm(r_all)] = sort(1 .+ 9 .* rand(rng_r, N); rev = true)
    ms1 = mass_segregation(pos, m_seg; seed = 1)
    @test ms1.lambda_msr > 5
    @test ms1.r_half_ratio < 0.4
    @test mass_segregation(pos, m_seg; seed = 1).lambda_msr == ms1.lambda_msr   # reproducible
    @test_throws ArgumentError mass_segregation(pos[:, 1:20], m_seg[1:20])       # ≤ n_massive
    @test_throws ArgumentError mass_segregation(pos, m_seg; n_massive = 1)
    @test Nbody6Dynamics._mst_length(pos, Int[]) == 0.0
    @test Nbody6Dynamics._mst_length([0.0 3.0 3.0; 0.0 0.0 4.0; 0.0 0.0 0.0], [1, 2, 3]) ≈ 7.0   # 3 + 4

    # Coalescence and the full time series on two synthetic snapshots:
    # separated clusters, then superposed with their orbital velocities.
    p1, v1 = sample_plummer(600, 0.3; rng = rng_r)
    m1 = fill(1.0 / 1200, 600)
    virialise!(m1, p1, v1)
    p2, v2 = sample_plummer(600, 0.3; rng = rng_r)
    m2 = fill(1.0 / 1200, 600)
    virialise!(m2, p2, v2)
    pos_o, vel_o, mass_o, ranges_o, _ = Nbody6Dynamics.setup_two_cluster_orbit(
        p1,
        v1,
        m1,
        p2,
        v2,
        m2,
        6.0,
        0.0;
        truncate_jacobi_flag = false,
    )
    function _rsnap(pos, vel, mass, t)
        params = zeros(Float32, 20)
        params[1] = t
        params[3] = 2.0f0     # rbar: 1 NB = 2 pc
        params[4] = 1.0f0
        params[11] = 3.0f0    # tscale: 1 NB = 3 Myr
        params[12] = 1.0f0
        Snapshot(
            SnapshotHeader(Int32(length(mass)), Int32(1), Int32(1), Int32(20), params),
            Int32.(1:length(mass)),
            Float32.(mass),
            Float32.(pos),
            Float32.(vel),
            Float32[],
            Float32[],
        )
    end
    s0 = _rsnap(pos_o, vel_o, mass_o, 0.0)
    # "Merged": the clusters overlap (0.5 apart, radii ≈ 0.35) and keep
    # their orbital velocities, so the pair's angular momentum survives
    # with the orbital sense; superposing them exactly would cancel it.
    pos_m = copy(pos_o)
    pos_m[:, ranges_o[1]] .= p1 .- [0.25, 0.0, 0.0]
    pos_m[:, ranges_o[2]] .= p2 .+ [0.25, 0.0, 0.0]
    s1 = _rsnap(pos_m, vel_o, mass_o, 1.0)
    coal = coalescence_time([s0, s1], ranges_o)
    @test coal.index == 2 && coal.time_nb == 1.0 && coal.time_myr == 3.0
    @test coal.n_distinct == [2, 1]
    far = coalescence_time([s0, s0], ranges_o)
    @test far.index === nothing && isnan(far.time_nb) && far.n_distinct == [2, 2]
    @test_throws ArgumentError coalescence_time([s0], ranges_o[1:1])
    L_orb = Nbody6Dynamics._orbital_angular_momentum(s0, ranges_o)
    @test abs(L_orb[3]) > 0.5 && abs(L_orb[1]) < 1e-6 && abs(L_orb[2]) < 1e-6   # orbit in the x–y plane

    diag = remnant_diagnostics([s0, s1], ranges_o; n_massive = 20, n_random = 20)
    @test diag.time == [0.0, 1.0] && diag.time_myr == [0.0, 3.0] && diag.rbar == [2.0, 2.0]
    @test diag.coalescence_time == 1.0 && diag.coalescence_time_myr == 3.0
    @test all(diag.n_bound .≥ 1100)
    @test all(0.9 .< diag.bound_mass_fraction .≤ 1.0)
    @test diag.r_core[1] > diag.r_core[2] > 0           # two separated clusters → one compact remnant
    @test diag.r_half[1] > diag.r_half[2] > 0
    @test diag.spin_alignment[1] ≈ 1.0 atol = 1e-3       # spin of the pair = orbital L at t = 0
    @test diag.spin_alignment[2] > 0.9                   # remnant keeps the orbital sense
    @test all(isfinite, diag.lambda_r) && all(isfinite, diag.lambda_peebles)
    @test all(isfinite, diag.lambda_msr) && all(isfinite, diag.segregation_ratio)
    @test isnan(diag.segregation_time)
    @test !isempty(diag.profile.radius)
    @test_throws ArgumentError remnant_diagnostics(Snapshot[], ranges_o)

    csv_path = joinpath(mktempdir(), "remnant.csv")
    @test write_remnant_diagnostics(csv_path, diag) == csv_path
    lines = readlines(csv_path)
    @test length(lines) == 4
    @test startswith(lines[1], "# coalescence_time_nb=1.0 coalescence_time_myr=3.0")
    @test lines[2] ==
          "time_nb,time_myr,rbar_pc,n_bound,bound_mass_fraction,r_core,r_half,lambda_r,lambda_peebles,spin_x,spin_y,spin_z,spin_alignment,lambda_msr,lambda_msr_err,segregation_ratio"
    @test startswith(lines[4], "1,3,2,")

    vis_r = VisualizationConfig(; format = "png", column = "single", output_dir = mktempdir())
    figs = remnant_figures(diag, vis_r)
    @test length(figs) == 4 && all(isfile, figs)
    vis_nb = VisualizationConfig(;
        format = "png",
        column = "single",
        units = "nbody",
        output_dir = vis_r.output_dir,
    )
    @test isfile(plot_remnant_structure(diag, vis_nb; filename = "structure_nb"))
    @test isfile(
        plot_mass_segregation_evolution(diag, vis_nb; lambda_threshold = 1.5, filename = "seg_nb"),
    )
    @test isfile(
        plot_rotation_profile(
            diag.profile,
            vis_r;
            rbar = 2.0,
            lambda_r = 0.3,
            filename = "profile",
        ),
    )
    @test_throws ErrorException plot_rotation_profile(
        RotationProfile(Float64[], Float64[], Float64[], Float64[]),
        vis_r,
    )
    rm(vis_r.output_dir; recursive = true, force = true)
end

# =====================================================================

@testset "Stellar classes, HR population and census" begin
    # Every K* of the SSE/BSE range belongs to exactly one class.
    @test all(k -> count(c -> k in c.kstar, STELLAR_CLASSES) == 1, -1:15)
    @test [stellar_class(k).key for k in (-1, 0, 1, 2, 3, 4, 6, 8, 11, 13, 14, 15)] == [
        :pre_main_sequence,
        :main_sequence,
        :main_sequence,
        :hertzsprung_gap,
        :red_giant,
        :core_helium_burning,
        :asymptotic_giant,
        :helium_star,
        :white_dwarf,
        :neutron_star,
        :black_hole,
        :massless_remnant,
    ]
    @test [c.key for c in STELLAR_CLASSES if !c.luminous] == [:neutron_star, :black_hole, :massless_remnant]
    @test_throws ArgumentError stellar_class_index(16)
    @test_throws ArgumentError stellar_class_index(-2)

    star(k, log_l, log_t) =
        StellarRecord(0.0, Int32(1), Int32(1), Int32(k), 1.0, 1.0, log_l, 0.0, log_t)
    pair(k1, k2, l1, l2, t1, t2) = BinaryRecord(
        0.0,
        Int32(1),
        Int32(2),
        Int32(1),
        Int32(2),
        Int32(k1),
        Int32(k2),
        Int32(0),
        0.0,
        0.1,
        0.0,
        3.0,
        1.0,
        1.0,
        l1,
        l2,
        0.0,
        0.0,
        t1,
        t2,
        zeros(12)...,
    )
    # t = 0: two MS stars, a black hole at the SSE placeholder, a neutron
    # star hot enough to pass any numeric cut; one KS pair, MS + CHeB.
    sev1 = StellarEvolutionSnapshot(
        0.0,
        4,
        [star(1, 2.0, 4.2), star(0, -1.0, 3.6), star(14, -10.0, 3.3), star(13, 0.4, 6.3)],
    )
    bev1 = BinaryEvolutionSnapshot(0.0, 1, [pair(1, 4, 1.0, 5.0, 4.0, 3.9)])
    # t = 2: an MS star, a Hertzsprung-gap star, an MS placeholder record;
    # no binary snapshot for this epoch.
    sev2 = StellarEvolutionSnapshot(
        2.0,
        3,
        [star(1, 2.0, 4.2), star(2, 4.0, 4.0), star(1, -10.0, 3.5)],
    )
    ms, hg, cheb, ns, bh = stellar_class_index.((1, 2, 4, 13, 14))

    pop = hr_population(sev1, bev1)
    @test length(pop) == 4                       # BH and NS are not on the plane
    @test pop.class == [ms, ms, ms, cheb]
    @test pop.binary_member == [false, false, true, true]
    @test pop.log_teff == [4.2, 3.6, 4.0, 3.9]
    @test length(hr_population(sev1)) == 2       # without the pair

    pops = hr_populations([sev1, sev2], [bev1])  # pairs by time, not position
    @test length(pops[1]) == 4
    @test pops[2].class == [ms, hg]              # placeholder record dropped

    census = stellar_census([sev1, sev2], [bev1])
    @test census.time_myr == [0.0, 2.0]
    @test census.single[1, [ms, ns, bh]] == [2, 1, 1]
    @test census.binary_member[1, [ms, cheb]] == [1, 1]
    @test census.single[2, [ms, hg]] == [2, 1]   # the census counts every record
    @test sum(census.binary_member[2, :]) == 0
    @test sum(class_counts(census)) == 6 + 3
    @test classes_present(census) == [ms, hg, cheb, ns, bh]

    mktempdir() do dir
        path = write_stellar_census(joinpath(dir, "stellar_census.csv"), census)
        rows = readlines(path)
        @test length(rows) == 3
        @test startswith(rows[1], "time_myr,pre_main_sequence_single,pre_main_sequence_binary,")
        @test length(split(rows[1], ',')) == 1 + 2 * length(STELLAR_CLASSES)
        @test startswith(rows[2], "0.0,0,0,2,1,")
        write_stellar_census(path, census)
        @test isfile(joinpath(dir, "stellar_census#1.csv"))
    end
end

@testset "Stellar type labels" begin
    # Standard Hurley et al. (2000) SSE/BSE table used by this fork
    @test startswith(STELLAR_TYPE_LABELS[0], "MS")
    @test startswith(STELLAR_TYPE_LABELS[1], "MS")
    @test startswith(STELLAR_TYPE_LABELS[10], "HeWD")
    @test startswith(STELLAR_TYPE_LABELS[12], "ONeWD")
    @test startswith(STELLAR_TYPE_LABELS[13], "NS")
    @test startswith(STELLAR_TYPE_LABELS[14], "BH")
    @test length(STELLAR_TYPE_LABELS) == 16   # K* = 0..15 complete
end

# =====================================================================
# Regression tests against verbatim Nbody6PPGPU-beijing output.
# The synthetic fixtures above mirror the readers' assumptions by
# construction; these fixtures were copied from a real run and would
# have caught the esc.11 column-map, M*-vs-<M>, and K*-label bugs.

@testset "Per-cluster structure" begin
    rng_s = StableRNG(2024)
    function _synthetic_cluster(N, a, centre, vcm, rng)
        pos, vel = sample_plummer(N, a; rng = rng)
        mass = fill(1.0 / N, N)
        virialise!(mass, pos, vel)
        pos .+= centre
        vel .+= vcm
        return pos, vel, mass
    end
    N1, N2 = 400, 300
    p1, v1, m1 = _synthetic_cluster(N1, 0.2, [-3.0, 0.0, 0.0], [0.0, 0.3, 0.0], rng_s)
    p2, v2, m2 = _synthetic_cluster(N2, 0.15, [3.0, 0.0, 0.0], [0.0, -0.3, 0.0], rng_s)
    # Unbind the last 5 % of cluster 1: radial speed well above escape (√(2/a) ≈ 3.2)
    n_unb = 20
    for j in (N1 - n_unb + 1):N1
        dir = p1[:, j] .- [-3.0, 0.0, 0.0]
        dir ./= max(sqrt(sum(dir .^ 2)), 1e-6)
        v1[:, j] .= [0.0, 0.3, 0.0] .+ 5.0 .* dir
    end
    pos = hcat(p1, p2)
    vel = hcat(v1, v2)
    mass = vcat(m1, m2)
    names = Int32.(1:(N1 + N2))
    function _snap(t)
        params = zeros(Float32, 20)
        params[1] = t
        params[3] = 1.0f0
        params[4] = 1.0f0
        params[11] = 1.0f0
        params[12] = 1.0f0
        Snapshot(
            SnapshotHeader(Int32(N1 + N2), Int32(1), Int32(1), Int32(20), params),
            names,
            Float32.(mass),
            Float32.(pos),
            Float32.(vel),
            Float32[],
            Float32[],
        )
    end
    snaps = [_snap(0.0), _snap(1.0)]
    ranges = [1:N1, (N1 + 1):(N1 + N2)]

    st = cluster_structure(snaps, ranges)
    @test st isa ClusterStructure
    @test st.time == [0.0, 1.0]
    @test st.n_members[:, 1] == [N1, N2]
    @test N1 - n_unb - 20 ≤ st.n_bound[1, 1] ≤ N1 - n_unb
    @test st.n_bound[2, 1] ≥ 0.9 * N2
    @test 0.88 ≤ st.bound_mass_fraction[1, 1] ≤ 0.95
    @test isapprox(st.centre[:, 1, 1], [-3.0, 0.0, 0.0]; atol = 0.05)
    @test isapprox(st.centre[:, 2, 1], [3.0, 0.0, 0.0]; atol = 0.05)
    @test isapprox(st.r_lagr[2, 1, 1], 1.305 * 0.2; rtol = 0.15)
    @test isapprox(st.r_lagr[2, 2, 1], 1.305 * 0.15; rtol = 0.15)
    @test st.r_lagr[1, 1, 1] < st.r_lagr[2, 1, 1] < st.r_lagr[3, 1, 1]
    @test 0.4 ≤ st.q_virial[2, 1] ≤ 0.6
    @test isfinite(st.sigma_1d[1, 1]) && st.sigma_1d[1, 1] > 0
    @test st.sigma_1d[1, 1] > st.sigma_1d[2, 1] * 0.5   # both of order √(M/r)

    # Unbound fast members inflate the all-member ratio; the bound selection removes them
    Q_all, _ = per_cluster_virial(snaps, ranges; bound_only = false)
    Q_bound, nm = per_cluster_virial(snaps, ranges)
    @test Q_all[1, 1] > Q_bound[1, 1]
    @test 0.4 ≤ Q_bound[1, 1] ≤ 0.65
    @test nm[1, 1] == N1
    st_all = cluster_structure(snaps, ranges; bound_only = false)
    @test st_all.n_bound[1, 1] == N1 && st_all.bound_mass_fraction[1, 1] == 1.0

    # Whole-system bound fraction: the 20 kicked stars are unbound to the pair as well
    fb = bound_fraction(snaps[1])
    @test 0.9 ≤ fb ≤ 0.98

    vis_p = VisualizationConfig(;
        format = "png",
        column = "single",
        units = "nbody",
        output_dir = mktempdir(),
    )
    # Radial profiles against the generating models
    rng_p = StableRNG(9090)
    Np = 6000
    pp, vp = sample_plummer(Np, 0.2; rng = rng_p)
    mp = fill(1.0 / Np, Np)
    virialise!(mp, pp, vp)
    prof = radial_profile(pp, vp, mp, zeros(3))
    @test prof isa RadialProfile && length(prof.r) == 12 && sum(prof.n) ≥ 0.98 * Np
    @test isapprox(prof.r_h, 1.305 * 0.2; rtol = 0.08)
    ρ_pl = model_density(PlummerProfile(), 1.0, prof.r_h)
    good = prof.n .≥ 80
    ratios = prof.rho[good] ./ [ρ_pl(x) for x in prof.r[good]]
    @test all(0.7 .≤ ratios .≤ 1.3)
    @test all(abs.(prof.beta[prof.n .≥ 150]) .< 0.25)           # isotropic Plummer
    σ_iso = prof.sigma_r[prof.n .≥ 150]
    @test all(0.3 .< σ_iso .< 3.0)
    @test_throws ArgumentError radial_profile(pp[:, 1:10], vp[:, 1:10], mp[1:10], zeros(3))
    # King W0 = 6, scaled to r_h = 0.3: the model integrates to the mass and half-mass radius
    pk, vk = sample_king(Np, 6.0, 1.0; rng = rng_p)
    mk = fill(1.0 / Np, Np)
    r_hk = Nbody6Dynamics.half_mass_radius(mk, pk; centre = zeros(3))
    pk .*= 0.3 / r_hk
    virialise!(mk, pk, vk)
    profk = radial_profile(pk, vk, mk, zeros(3))
    ρ_k = model_density(KingProfile(W0 = 6.0), 1.0, 0.3)
    goodk = profk.n .≥ 80
    ratk = profk.rho[goodk] ./ [ρ_k(x) for x in profk.r[goodk]]
    @test all(0.65 .≤ ratk .≤ 1.35)
    @test ρ_k(100.0) == 0.0                                       # beyond the tidal radius
    r_grid = exp10.(range(-3, 1; length = 400))
    m_int = sum(
        4π * r_grid[i]^2 * ρ_k(r_grid[i]) * (r_grid[i + 1] - r_grid[i]) for
        i in 1:(length(r_grid) - 1)
    )
    @test isapprox(m_int, 1.0; rtol = 0.05)
    # Per-cluster profiles on the synthetic pair and the system profile
    profs = cluster_profiles(snaps[1], ranges)
    @test length(profs) == 2 && all(!isnothing, profs)
    @test isapprox(profs[2].r_h, 1.305 * 0.15; rtol = 0.15)
    sp = system_profile(snaps[1])
    @test sp.M ≈ sum(mass) && sp.r_h > profs[1].r_h                 # the pair is wider than a member
    # Figures: with and without the generating models, physical off
    specs_s = [
        ClusterSpec(profile = PlummerProfile(), N = N1, rbar = 1.305 * 0.2, imf = KroupaIMF()),
        ClusterSpec(profile = PlummerProfile(), N = N2, rbar = 1.305 * 0.15, imf = KroupaIMF()),
    ]
    pd1 = plot_density_profiles(snaps[1], ranges, vis_p; specs = specs_s)
    @test isfile(pd1)
    pd2 = plot_density_profiles(snaps[1], ranges, vis_p; filename = "density_nomodel")
    @test isfile(pd2)
    pv = plot_velocity_dispersion(snaps[1], ranges, vis_p)
    @test isfile(pv)
    @test_throws ArgumentError plot_density_profiles(snaps[1], ranges, vis_p; specs = specs_s[1:1])

    # Figure: with and without the engine overlay (NB units, raster draft)
    vis_s = VisualizationConfig(;
        format = "png",
        column = "single",
        units = "nbody",
        output_dir = mktempdir(),
    )
    path = plot_cluster_structure(snaps, ranges, vis_s)
    @test isfile(path)
    lagr_s = LagrangianData([0.0, 1.0], [0.1, 0.5, 0.9], [0.5 0.5; 3.0 3.0; 6.0 6.0])
    path2 =
        plot_cluster_structure(snaps, ranges, vis_s; lagr = lagr_s, filename = "structure_overlay")
    @test isfile(path2)
end

# =====================================================================

@testset "half_mass_radius boundaries" begin
    # Single particle: the half-mass radius is that particle's radius
    @test Nbody6Dynamics.half_mass_radius([2.0], zeros(3, 1); centre = zeros(3)) == 0.0
    # Two equal masses: cumulative ≥ M/2 at the inner particle
    pos2 = [1.0 2.0; 0.0 0.0; 0.0 0.0]
    @test Nbody6Dynamics.half_mass_radius([0.5, 0.5], pos2; centre = zeros(3)) ≈ 1.0
end

@testset "Kepler near-parabolic limit" begin
    # Apocentre speed → 0 as e → 1 (vis-viva); must stay finite/nonneg
    v1, v2 = kepler_velocity(1.0, 1.0, 10.0, 1.0 - 1e-12)
    @test 0.0 ≤ v1 < 1e-5 && 0.0 ≤ v2 < 1e-5
end

# =====================================================================
# Static QA ships with the tests.
# =====================================================================
