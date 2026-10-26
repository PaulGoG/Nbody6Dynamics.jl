# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# IMBH seeding of the generator (kz24 route): the [merger.clusterK.imbh] keys
# and their guards, the IMBH as the last body of dat.10 repeated on the
# engine's IMBHINIT line (identity overwrite under the engine's own unit
# conversion), the velocity modes with the Brownian equilibrium of Merritt
# (2005, DOI 10.1086/429964), and the isotropic Jeans response against its
# virial identity. Samples use StableRNGs so the reference values hold across
# Julia versions.

using StableRNGs: StableRNG

@testset "IMBH seeding" begin
    nd = Nbody6Dynamics

    function imbh_config(;
        imbh::String = "imbh = {mass = 50.0, velocity = \"equipartition\"}",
        kz19::Int = 0,
        escapers::String = "keep",
        nbody6_lines::String = "",
    )
        text = """
        [merger]
        seed = 5
        n_clusters = 1
        orbit_mode = "single"

        [merger.cluster1]
        N = 2000
        rbar = 1.0
        model = "king"
        W0 = 7.0
        imf = "equal"
        mass_total = 1000.0
        $imbh

        [merger.stellar]
        kz19 = $kz19

        [merger.nbody6]
        escapers = "$escapers"
        $nbody6_lines
        """
        path = joinpath(mktempdir(), "imbh.toml")
        write(path, text)
        return path
    end

    # Mass-weighted mean of the columns `idx` of a 3×N matrix
    mean_of(x, m, idx) = [sum(m[i] * x[k, i] for i in idx) / sum(m[idx]) for k in 1:3]

    @testset "Parsing and guards" begin
        cfg = nd.load_merger_config(imbh_config())
        s = cfg.clusters[1].imbh
        @test s isa IMBHSpec
        @test s.mass == 50.0 && isnan(s.mass_fraction) && s.velocity == "equipartition"
        @test cfg.nbody6.imbh_route == "auto" && cfg.nbody6.imbh_output_dt_myr == 1.0
        # the eight-argument form of earlier versions carries no IMBH
        c8 = ClusterSpec(
            100,
            1.0,
            PlummerProfile(),
            KroupaIMF(),
            Float64[],
            Float64[],
            BinarySpec(),
            SegregationSpec(),
        )
        @test c8.imbh === nothing

        @test_throws ArgumentError nd.load_merger_config(
            imbh_config(; imbh = "imbh = {mass = 50.0, mass_fraction = 0.01}"),
        )
        @test_throws ArgumentError nd.load_merger_config(
            imbh_config(; imbh = "imbh = {velocity = \"rest\"}"),
        )
        @test_throws ArgumentError nd.load_merger_config(
            imbh_config(; imbh = "imbh = {mass_fraction = 0.06}"),
        )
        @test_logs (:warn, r"mass_fraction") match_mode = :any nd.load_merger_config(
            imbh_config(; imbh = "imbh = {mass_fraction = 0.03}"),
        )
        @test_throws ArgumentError nd.load_merger_config(
            imbh_config(; nbody6_lines = "imbh_route = \"flag\""),
        )
        @test_throws ArgumentError nd.load_merger_config(
            imbh_config(; imbh = "", nbody6_lines = "imbh_route = \"flag\""),
        )
        @test_throws ArgumentError nd.load_merger_config(
            imbh_config(; nbody6_lines = "kz = {24 = 1}"),
        )
        @test_throws ArgumentError nd.load_merger_config(imbh_config(; kz19 = 3))
        @test_logs (:warn, r"restart") match_mode = :any nd.load_merger_config(
            imbh_config(; escapers = "remove"),
        )
        @test_throws ArgumentError nd.load_merger_config(
            imbh_config(; nbody6_lines = "imbh_output_dt_myr = 0"),
        )
        @test_throws ArgumentError nd.load_merger_config(
            imbh_config(; imbh = "imbh = {mass = 50.0, placement = \"edge\"}"),
        )

        # two clusters with one IMBH each
        two = """
        [merger]
        seed = 5
        n_clusters = 2
        orbit_mode = "kepler"

        [merger.cluster1]
        N = 2000
        rbar = 1.0
        model = "plummer"
        imbh = {mass = 50.0}

        [merger.cluster2]
        N = 2000
        rbar = 1.0
        model = "plummer"
        imbh = {mass = 50.0}

        [merger.stellar]
        kz19 = 0

        [merger.orbit]
        apocentre = 10.0
        eccentricity = 0.5
        """
        path_two = joinpath(mktempdir(), "two.toml")
        write(path_two, two)
        @test_throws ArgumentError nd.load_merger_config(path_two)
    end

    @testset "Generation and the IMBHINIT line" begin
        cfg = nd.load_merger_config(imbh_config())
        dir = mktempdir()
        result = nd.generate_merger_ic(cfg; rng = StableRNG(7), output_dir = dir)
        @test result.N_total == 2001
        M_stars = sum(result.mass_physical[1:2000])
        @test result.M_total ≈ M_stars + 50
        @test collect(result.cluster_ranges[1]) == collect(1:2000)

        rows = [parse.(Float64, split(l)) for l in readlines(joinpath(dir, "dat.10"))]
        @test length(rows) == 2001
        last_row = rows[end]
        @test isapprox(last_row[1], 50 / result.M_total; rtol = 1e-12)

        inp = read(joinpath(dir, "merger.inp"), String)
        kz_line = match(r"KZ\(21:30\)=([0-9 -]+)", inp).captures[1]
        @test parse(Int, split(kz_line)[4]) == 1
        units = match(r"RBAR=([0-9.]+),ZMBAR=([0-9.]+)", inp)
        rbar_w = parse(Float64, units.captures[1])
        zmbar_w = parse(Float64, units.captures[2])
        zmass_w = zmbar_w * 2001
        numeric_lines = filter(readlines(joinpath(dir, "merger.inp"))) do l
            p = split(l)
            length(p) == 8 && all(x -> tryparse(Float64, x) !== nothing, p)
        end
        @test length(numeric_lines) == 1
        lines = readlines(joinpath(dir, "merger.inp"))
        i_line = findfirst(==(numeric_lines[1]), lines)
        @test i_line > findfirst(==("&INSCALE"), lines)
        v = parse.(Float64, split(numeric_lines[1]))
        vstar = nd._imbh_engine_vstar_kms(zmass_w, rbar_w)
        @test isapprox(v[1] / zmass_w, last_row[1]; rtol = 1e-12)
        @test isapprox(v[1], 50.0; rtol = 1e-5)
        for k in 1:3
            @test isapprox(v[1 + k] / rbar_w, last_row[1 + k]; rtol = 1e-12)
            @test isapprox(v[4 + k] / vstar, last_row[4 + k]; rtol = 1e-12)
        end

        meta = nd.TOML.parsefile(joinpath(dir, "merger_ic.toml"))
        t_star_myr = meta["meta"]["t_star_myr"]
        @test isapprox(v[8], 1.0 / t_star_myr; rtol = 1e-12)
        @test haskey(meta, "imbh")
        @test meta["imbh"]["name"] == 2001
        @test meta["imbh"]["route"] == "kz24"
        @test meta["imbh"]["velocity_mode"] == "equipartition"
        @test isempty(meta["imbh"]["jeans_factors"])

        summary_path = joinpath(dir, "merger_summary.txt")
        @test occursin("IMBH: name 2001", read(summary_path, String))
        ranges = nd.parse_merger_summary(summary_path)
        @test length(ranges) == 1 && collect(ranges[1]) == collect(1:2000)
    end

    @testset "Velocity modes" begin
        rest = nd.generate_merger_ic(
            nd.load_merger_config(
                imbh_config(; imbh = "imbh = {mass = 50.0, velocity = \"rest\"}"),
            );
            rng = StableRNG(7),
            output_dir = mktempdir(),
        )
        v_cm = mean_of(rest.vel_physical, rest.mass_physical, 1:2000)
        @test isapprox(rest.vel_physical[:, end], v_cm; atol = 1e-9)

        explicit = nd.generate_merger_ic(
            nd.load_merger_config(
                imbh_config(;
                    imbh = "imbh = {mass = 50.0, velocity = \"explicit\", velocity_kms = [1.0, 0.0, 0.0]}",
                ),
            );
            rng = StableRNG(7),
            output_dir = mktempdir(),
        )
        v_cm_e = mean_of(explicit.vel_physical, explicit.mass_physical, 1:2000)
        @test isapprox(explicit.vel_physical[:, end] .- v_cm_e, [1.0, 0.0, 0.0]; atol = 1e-9)

        # Brownian equilibrium: per-component variance σ_c² ⟨m⟩/M• = 1/100
        spec = IMBHSpec(; mass = 100.0, velocity = "equipartition")
        rng = StableRNG(1)
        draws = reduce(
            hcat,
            [nd._imbh_body(spec, 100.0, 1.0, zeros(3), zeros(3), 1.0, rng)[2] for _ in 1:200],
        )
        for k in 1:3
            x = draws[k, :]
            μ = sum(x) / length(x)
            s² = sum(abs2, x .- μ) / (length(x) - 1)
            @test abs(s² / 0.01 - 1) < 0.3
        end
    end

    @testset "Jeans response" begin
        N = 4000
        spec = ClusterSpec(;
            N = N,
            rbar = 1.0,
            profile = PlummerProfile(),
            imf = EqualMassIMF(; particle_mass = 1 / N),
        )
        sampled = nd._sample_cluster(spec; rng = StableRNG(3), nmax = N)
        pos, vel, mass = copy(sampled.pos), copy(sampled.vel), copy(sampled.mass)
        members = collect(1:N)
        x_c, v_c = nd._imbh_centre(pos, vel, mass, members)
        M_bh = 0.02
        kinetic(v) = 0.5 * sum(mass[i] * sum(abs2, v[:, i] .- v_c) for i in members)
        K_old = kinetic(vel)
        f = nd._imbh_jeans_response!(pos, vel, mass, members, v_c, M_bh, x_c; n_shell = 100)
        K_new = kinetic(vel)
        @test all(≥(1), f)
        @test f[1] == maximum(f)
        W_ss = 0.0
        for i in 1:(N - 1), j in (i + 1):N
            r_ij = sqrt(
                (pos[1, i] - pos[1, j])^2 + (pos[2, i] - pos[2, j])^2 + (pos[3, i] - pos[3, j])^2,
            )
            W_ss -= mass[i] * mass[j] / r_ij
        end
        W_sbh = -M_bh * sum(mass[i] / sqrt(sum(abs2, pos[:, i] .- x_c)) for i in members)
        @test abs(2K_new + W_ss + W_sbh) ≤ 0.1 * abs(W_sbh)
        @test K_new > K_old
    end

    @testset "Offset placement" begin
        res = nd.generate_merger_ic(
            nd.load_merger_config(
                imbh_config(;
                    imbh = "imbh = {mass = 50.0, placement = \"offset\", offset = [0.5, 0.0, 0.0]}",
                ),
            );
            rng = StableRNG(7),
            output_dir = mktempdir(),
        )
        x_cm = mean_of(res.pos_physical, res.mass_physical, 1:2000)
        @test isapprox(sqrt(sum(abs2, res.pos_physical[:, end] .- x_cm)), 0.5; atol = 1e-9)
    end
end
