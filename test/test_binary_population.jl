# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# WP10 10.3 increment 1: binary population series (design W8) — number, star
# and mass fractions per radial shell, stellar class and census level, on
# synthetic snapshots with an injected, radially dependent binary population.

@testset "Binary population series (10.3)" begin
    using StableRNGs
    using Random: randperm
    nd = Nbody6Dynamics

    # Injected population: Plummer systems, binaries inside r < 1 with probability 0.6 and outside with 0.2;
    # hard circular pairs (a = 1e-5) and bulk speeds of 50 so that no two systems are bound to each other.
    rng = StableRNG(103)
    n_sys = 2000
    xs, _ = nd.sample_plummer(n_sys, 1.0; rng = rng)
    m_sys = (0.2 .+ 0.8 .* rand(rng, n_sys)) ./ n_sys
    r_sys = vec(sqrt.(sum(abs2, xs; dims = 1)))
    is_bin = [rand(rng) < (r_sys[i] < 1.0 ? 0.6 : 0.2) for i in 1:n_sys]
    q_sys = 0.2 .+ 0.8 .* rand(rng, n_sys)
    vbulk = [50.0 .* (u -> u ./ sqrt(sum(abs2, u)))(randn(rng, 3)) for _ in 1:n_sys]
    zmbar = 1000.0
    function bodies_of(i, name0)
        x, v, m = xs[:, i], vbulk[i], m_sys[i]
        is_bin[i] || return [(name0, m, x, v)]
        m1 = m / (1 + q_sys[i])
        m2 = m - m1
        a = 1.0e-5
        vrel = sqrt(m / a)
        ehat = [1.0, 0, 0]
        uhat = [0, 1.0, 0]
        return [
            (name0, m1, x .- (m2 / m) * a .* ehat, v .- (m2 / m) * vrel .* uhat),
            (name0 + 1, m2, x .+ (m1 / m) * a .* ehat, v .+ (m1 / m) * vrel .* uhat),
        ]
    end
    bins_idx = findall(is_bin)
    sing_idx = findall(.!is_bin)
    n_b = length(bins_idx)
    # Level-E snapshot: the pairs first (KS components 2k − 1, 2k), then the singles
    bodies = vcat(
        reduce(vcat, [bodies_of(i, 2k - 1) for (k, i) in enumerate(bins_idx)]),
        [bodies_of(i, 2n_b + k)[1] for (k, i) in enumerate(sing_idx)],
    )
    function snap(t, bodies; npairs = 0)
        n = length(bodies)
        params = zeros(Float32, 20)
        params[1] = t
        params[3] = 1.0
        params[4] = zmbar
        params[11] = 1.0
        hdr = SnapshotHeader(Int32(n), Int32(0), Int32(1), Int32(20), params)
        return Snapshot(
            hdr,
            Int32[b[1] for b in bodies],
            Float32[b[2] for b in bodies],
            Float32.(reduce(hcat, [b[3] for b in bodies])),
            Float32.(reduce(hcat, [b[4] for b in bodies])),
            Float32[],
            Float32[],
            npairs,
        )
    end
    s_ks = snap(0.0, bodies; npairs = n_b)
    # The same bodies with no KS pair, in a shuffled order: level P must find the pairs itself
    s_free = snap(1.0, bodies[randperm(rng, length(bodies))])

    # Oracle: systems from the snapshot data in Float64 (centre of mass of each injected pair), r_h of the systems
    mass = Float64.(s_ks.mass)
    pos = Float64.(s_ks.pos)
    sys_m = vcat([mass[2k - 1] + mass[2k] for k in 1:n_b], mass[(2n_b + 1):end])
    sys_x = hcat(
        reduce(
            hcat,
            [(mass[2k - 1] .* pos[:, 2k - 1] .+ mass[2k] .* pos[:, 2k]) ./ sys_m[k] for k in 1:n_b],
        ),
        pos[:, (2n_b + 1):end],
    )
    sys_r = vec(sqrt.(sum(abs2, sys_x; dims = 1)))
    sys_bin = vcat(trues(n_b), falses(length(sing_idx)))
    sys_q = vcat(
        [min(mass[2k - 1], mass[2k]) / max(mass[2k - 1], mass[2k]) for k in 1:n_b],
        zeros(length(sing_idx)),
    )
    sys_m1 = vcat([max(mass[2k - 1], mass[2k]) for k in 1:n_b], mass[(2n_b + 1):end])
    order = sortperm(sys_r)
    cum = cumsum(sys_m[order]) ./ sum(sys_m)
    r_h = sys_r[order[findfirst(≥(0.5), cum)]]

    edges = [0.0, 0.5, 1.0, 2.0, Inf]
    window = (0.15, 0.35)
    rows = nd.binary_population_series(
        [s_ks, s_free];
        levels = (:ks, :bound),
        radius_frame = :rh,
        edges = edges,
        centre = :origin,
        q_min = 0.5,
        primary_window_msun = window,
    )
    @test rows isa Vector{nd.BinaryPopulationRow}
    # One row per snapshot × level × bin (0 = all radii, then 1…4) × class (only :all without stellar types)
    @test length(rows) == 2 * 2 * 5
    @test all(r -> r.class == :all && r.radius_frame == :rh, rows)
    pick(s, lev, b) = only(filter(r -> r.snapshot == s && r.level == lev && r.bin == b, rows))

    @testset "Exact recovery of the injected population" begin
        for (s, lev) in ((1, :ks), (1, :bound), (2, :bound))
            for b in 0:4
                sel =
                    b == 0 ? trues(length(sys_r)) : (edges[b] * r_h .≤ sys_r .< edges[b + 1] * r_h)
                row = pick(s, lev, b)
                @test row.n_binaries == count(sel .& sys_bin)
                @test row.n_singles == count(sel .& .!sys_bin)
                @test row.n_in_binaries == 2 * row.n_binaries
                @test row.m_binaries_msun ≈ zmbar * sum(sys_m[sel .& sys_bin]; init = 0.0) rtol =
                    1e-5
                @test row.m_singles_msun ≈ zmbar * sum(sys_m[sel .& .!sys_bin]; init = 0.0) rtol = 1e-5
                b > 0 && @test row.r_lo_nb ≈ edges[b] * r_h rtol = 1e-6
                b > 0 && @test row.edge_lo == edges[b] && row.edge_hi == edges[b + 1]
            end
        end
        # Level E of the snapshot without KS pairs is empty
        @test pick(2, :ks, 0).n_binaries == 0 && pick(2, :ks, 0).n_singles == length(bodies)
        # Shells partition the population
        for (s, lev) in ((1, :ks), (2, :bound))
            @test sum(pick(s, lev, b).n_binaries for b in 1:4) == pick(s, lev, 0).n_binaries
            @test sum(pick(s, lev, b).n_singles for b in 1:4) == pick(s, lev, 0).n_singles
        end
        # The radial dependence: more binaries per system inside r_h than outside 2 r_h
        @test pick(1, :ks, 1).f_b_number > 2 * pick(1, :ks, 4).f_b_number
    end

    @testset "Fraction definitions" begin
        row = pick(1, :ks, 0)
        f = n_b / n_sys
        @test row.f_b_number ≈ f
        @test row.f_b_stars ≈ 2f / (1 + f)
        @test row.f_b_mass ≈ sum(sys_m[sys_bin]) / sum(sys_m) rtol = 1e-5
        @test row.f_b_number_se ≈ sqrt(f * (1 - f) / n_sys)
        # q > 0.5 analogue: systems whose single or primary mass lies in the window [M☉]
        inwin = window[1] .≤ zmbar .* sys_m1 .≤ window[2]
        @test row.f_b_q ≈ count(inwin .& sys_bin .& (sys_q .> 0.5)) / count(inwin)
        # Half-mass radii of binaries and singles on the all-radii row only
        @test row.r_h_binaries_nb < row.r_h_singles_nb
        @test isnan(pick(1, :ks, 2).r_h_binaries_nb)
        @test row.t_nb == 0.0 && pick(2, :bound, 0).t_nb == 1.0
        @test row.n_hierarchies == 0
    end

    @testset "Lagrangian shells" begin
        fr = [0.0, 0.25, 0.5, 0.75, 1.0]
        lrows = nd.binary_population_series(
            [s_ks];
            levels = (:ks,),
            radius_frame = :lagrangian,
            edges = fr,
            centre = :origin,
        )
        @test length(lrows) == 5
        shell(b) = only(filter(r -> r.bin == b, lrows))
        @test sum(shell(b).n_binaries + shell(b).n_singles for b in 1:4) == n_sys
        for b in 1:4
            msh = (shell(b).m_binaries_msun + shell(b).m_singles_msun) / zmbar
            # Each shell holds a quarter of the mass to within one system
            @test abs(msh - 0.25 * sum(sys_m)) ≤ maximum(sys_m) + 1e-9
            @test shell(b).edge_lo == fr[b] && shell(b).edge_hi == fr[b + 1]
        end
        @test shell(1).r_lo_nb == 0.0 && shell(2).r_lo_nb == shell(1).r_hi_nb
    end

    @testset "Stellar classes" begin
        # Black holes: the 40 heaviest bodies of the snapshot; every other body a main-sequence star
        heavy = Set(Int.(s_ks.name[partialsortperm(s_ks.mass, 1:40; rev = true)]))
        kmap = Dict(Int(nm) => (Int(nm) in heavy ? 14 : 1) for nm in s_ks.name)
        crows = nd.binary_population_series(
            [s_ks];
            levels = (:ks,),
            radius_frame = :rh,
            edges = edges,
            centre = :origin,
            kstar = [kmap],
        )
        @test Set(r.class for r in crows) == Set([:all, :main_sequence, :black_hole])
        @test length(crows) == 5 * 3
        c(cl) = only(filter(r -> r.bin == 0 && r.class == cl, crows))
        isbh(i) = Int(s_ks.name[i]) in heavy
        # Binaries counted by the class of the primary (the heavier component)
        n_bh_primary = count(k -> isbh(mass[2k - 1] ≥ mass[2k] ? 2k - 1 : 2k), 1:n_b)
        @test c(:black_hole).n_binaries == n_bh_primary
        @test c(:main_sequence).n_binaries == n_b - n_bh_primary
        # Objects of the class in binaries (any role) over all objects of the class
        n_bh_in_bin = count(isbh, 1:(2n_b))
        @test c(:black_hole).n_in_binaries == n_bh_in_bin
        @test c(:black_hole).n_singles == count(isbh, (2n_b + 1):length(bodies))
        @test c(:black_hole).f_b_stars ≈ n_bh_in_bin / 40
        @test c(:all).n_binaries == n_b
        # A body without a stellar type is an error naming it
        bad = copy(kmap)
        delete!(bad, Int(s_ks.name[end]))
        @test_throws ArgumentError nd.binary_population_series(
            [s_ks];
            levels = (:ks,),
            centre = :origin,
            kstar = [bad],
        )
    end

    @testset "Arguments and writer" begin
        @test_throws ArgumentError nd.binary_population_series([s_ks]; levels = (:other,))
        @test_throws ArgumentError nd.binary_population_series([s_ks]; radius_frame = :other)
        @test_throws ArgumentError nd.binary_population_series(
            [s_ks];
            radius_frame = :lagrangian,
            edges = [0.0, 0.6, 0.5, 1.0],
        )
        @test_throws ArgumentError nd.binary_population_series(
            [s_ks];
            kstar = [Dict{Int,Int}(), Dict{Int,Int}()],
        )
        @test_throws ArgumentError nd.binary_population_series(Snapshot[])
        # Default frame and edges, shrinking-sphere centre: runs and partitions the systems
        drows = nd.binary_population_series([s_ks]; levels = (:ks,))
        @test drows[1].radius_frame == :lagrangian && drows[1].bin == 0
        @test sum(r.n_binaries for r in drows if r.bin > 0) == n_b
        path = joinpath(mktempdir(), "binary_population.csv")
        @test nd.write_binary_population(path, rows) == path
        lines = readlines(path)
        @test lines[1] == join(string.(fieldnames(nd.BinaryPopulationRow)), ",")
        @test length(lines) == length(rows) + 1
    end

    @testset "Stellar-type map" begin
        fx = joinpath(@__DIR__, "fixtures")
        sev = nd.read_stellar_evolution(joinpath(fx, "sev.83_0"))
        bev = nd.read_binary_evolution(joinpath(fx, "bev.82_0"))
        @test nd.stellar_type_map(sev, nothing) ==
              Dict(Int(r.name) => Int(r.stellar_type) for r in sev.records)
        m = nd.stellar_type_map(sev, bev)
        # The pair components of bev.82 take precedence over sev.83
        bev_names = Set(Iterators.flatten((Int(r.name1), Int(r.name2)) for r in bev.records))
        @test all(
            m[Int(r.name1)] == r.stellar_type1 && m[Int(r.name2)] == r.stellar_type2 for
            r in bev.records
        )
        @test all(
            m[Int(r.name)] == r.stellar_type for r in sev.records if !(Int(r.name) in bev_names)
        )
        @test keys(m) == union(bev_names, Set(Int(r.name) for r in sev.records))
    end
end
