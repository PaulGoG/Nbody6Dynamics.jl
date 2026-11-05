# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# N ladder (roadmap 4.5): rungs that differ in N alone, with equal-mass bodies
# M/N at fixed cluster mass, half-mass radius and orbit. The generated initial
# conditions carry the same mass, length and time units on every rung, up to
# the sampling scatter of the realised half-mass radius.

@testset "N ladder" begin
    nd = Nbody6Dynamics
    G = nd._G_PC_KMS2_MSUN
    dir = mktempdir()
    write(
        joinpath(dir, "pipeline.toml"),
        """
        [install]
        enabled = false
        install_dir = "engine"

        [simulation]
        input_file = "unused.inp"
        runs_dir = "runs"

        [postprocess]
        enabled = false

        [visualization]
        enabled = false

        [merger]
        enabled = false
        config_file = ""
        """,
    )
    sweep = Dict{String,Any}("pipeline_config" => "../pipeline.toml", "concurrency" => 2)
    ns = [256, 512, 1024, 2048]
    rung_seeds = [[11, 12, 13], [21, 22], [31], [41]]
    quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
    TOML = nd.TOML

    # The realised half-mass radius, the length unit RBAR, scatters about the
    # input r_h: each cluster is scaled to r_h about the model centre and then
    # recentred on its centre of mass. Over 20 seeds at N = 256 … 8192 the
    # relative scatter is ≈ 0.35 N^{-1/2} for a single King cluster and for a
    # Plummer pair; 2 N^{-1/2} bounds it by about six standard deviations, and
    # T* ∝ RBAR^{3/2} by 3 N^{-1/2}.
    tol_r(N) = 2 / sqrt(N)

    single = joinpath(dir, "king.toml")
    write(
        single,
        """
        [merger]
        seed = 5
        n_clusters = 1
        orbit_mode = "single"

        [merger.cluster1]
        model = "king"
        W0 = 6.0
        N = 1000
        rbar = 1.0
        imf = "equal"
        mass_total = 1.0e5

        [merger.output]
        tcrit = 50.0
        dtadj = 0.5
        deltat = 0.5

        [merger.stellar]
        kz19 = 0
        level = "0"
        """,
    )

    @testset "Rung files" begin
        out = joinpath(dir, "single")
        rungs = n_ladder(single, out; n_total = ns, seeds = rung_seeds, sweep = sweep)
        @test [r.n_total for r in rungs] == ns
        @test [r.n for r in rungs] == [[n] for n in ns]
        for (r, s) in zip(rungs, rung_seeds)
            @test r.seeds == s
            @test basename(r.merger_config) == "king_n$(r.n_total).toml"
            @test basename(r.sweep_config) == "king_n$(r.n_total)_sweep.toml"
            @test startswith(readline(r.merger_config), "# N ladder king, rung")
            cfg = load_merger_config(r.merger_config)
            @test cfg.clusters[1].N == r.n_total && cfg.seed == first(s)
            @test cfg.clusters[1].imf.particle_mass == 1.0e5 / r.n_total
            sw = load_sweep_config(r.sweep_config)
            @test sw.name == "king_n$(r.n_total)" && sw.seeds == s && sw.concurrency == 2
            @test sw.merger_config == r.merger_config
            @test sw.pipeline_config == joinpath(dir, "pipeline.toml")
            @test isempty(sw.grid)
        end
        # A second call keeps the first files as backups
        quiet() do
            n_ladder(single, out; n_total = ns[1:1], seeds = rung_seeds[1:1], sweep = sweep)
        end
        @test isfile(joinpath(out, "king_n256#1.toml"))
        @test isfile(joinpath(out, "king_n256_sweep#1.toml"))
        # Sweeps of the rungs prepare without a grid (seed ensembles)
        sw = load_sweep_config(rungs[1].sweep_config)
        sdir, points = quiet(() -> prepare_sweep(sw; sweep_dir = joinpath(dir, "sweep_256")))
        @test [p.seed for p in points] == rung_seeds[1]
        @test load_merger_config(joinpath(sdir, points[2].id, "merger.toml")).clusters[1].N == 256
    end

    @testset "Invariants of a single cluster" begin
        rungs = n_ladder(
            single,
            joinpath(dir, "single_ic");
            n_total = ns,
            seeds = rung_seeds,
            sweep = sweep,
        )
        t_star_model = nd._nbody_time_myr(1.0e5, 1.0)
        for r in rungs
            N = r.n_total
            res = quiet(
                () -> generate_merger_ic(
                    load_merger_config(r.merger_config);
                    output_dir = joinpath(dir, "ic_single_$N"),
                ),
            )
            meta = TOML.parsefile(joinpath(res.output_dir, "merger_ic.toml"))["meta"]
            @test res.N_total == N == meta["N_total"]
            @test all(==(1.0e5 / N), res.mass_physical)
            @test res.M_total ≈ 1.0e5 rtol = 1e-12
            @test meta["M_total"] ≈ 1.0e5 rtol = 1e-12
            @test meta["zmbar"] ≈ 1.0e5 / N rtol = 1e-12
            @test meta["q_virial"] ≈ 0.5 atol = 1e-10
            @test abs(meta["rbar"] - 1.0) ≤ tol_r(N)
            @test nd.half_mass_radius(res.mass_physical, res.pos_physical) ≈ meta["rbar"] rtol =
                1e-12
            @test meta["t_star_myr"] ≈ nd._nbody_time_myr(meta["M_total"], meta["rbar"]) rtol =
                1e-12
            @test abs(meta["t_star_myr"] / t_star_model - 1) ≤ 1.5 * tol_r(N)
            out = TOML.parsefile(joinpath(res.output_dir, "merger_ic.toml"))["output"]
            @test (out["tcrit"], out["dtadj"], out["deltat"]) == (50.0, 0.5, 0.5)
        end
    end

    @testset "Invariants of a cluster pair" begin
        pair = joinpath(dir, "pair_merger.toml")
        write(
            pair,
            """
            [merger]
            seed = 9
            n_clusters = 2
            orbit_mode = "kepler"

            [merger.cluster1]
            model = "plummer"
            N = 384
            rbar = 1.0
            imf = { type = "equal", target_mass = 7.5e4 }

            [merger.cluster2]
            model = "plummer"
            N = 128
            rbar = 0.8
            imf = "equal"
            mass_total = 2.5e4

            [merger.orbit]
            apocentre = 10.0
            eccentricity = 0.5

            [merger.output]
            truncate_jacobi = false
            tcrit = 20.0
            dtadj = 0.5
            deltat = 0.5

            [merger.stellar]
            kz19 = 0
            """,
        )
        rungs =
            n_ladder(pair, joinpath(dir, "pair"); n_total = ns, seeds = rung_seeds, sweep = sweep)
        @test [r.n for r in rungs] == [[3n ÷ 4, n ÷ 4] for n in ns]
        @test basename(rungs[1].merger_config) == "pair_n256.toml"
        masses, rbars = (7.5e4, 2.5e4), (1.0, 0.8)
        v_apo = sqrt(G * 1.0e5 * (1 - 0.5) / 10.0)   # relative speed at apocentre [km s⁻¹]
        orbit_ref = nothing
        t_star_ref = nothing
        for r in rungs
            N = r.n_total
            res = quiet(
                () -> generate_merger_ic(
                    load_merger_config(r.merger_config);
                    output_dir = joinpath(dir, "ic_pair_$N"),
                ),
            )
            ic = TOML.parsefile(joinpath(res.output_dir, "merger_ic.toml"))
            @test res.N_total == N
            @test all(==(1.0e5 / N), res.mass_physical)   # equal masses across both clusters
            centres = map(1:2) do i
                idx = res.cluster_ranges[i]
                m = res.mass_physical[idx]
                Mi = sum(m)
                @test Mi ≈ masses[i] rtol = 1e-12
                @test length(idx) == r.n[i]
                x = res.pos_physical[:, idx] * m / Mi
                v = res.vel_physical[:, idx] * m / Mi
                # Internal virial ratio and half-mass radius about the cluster's own centre
                T, W = nd._kinetic_and_potential(
                    m,
                    res.pos_physical[:, idx],
                    res.vel_physical[:, idx] .- v,
                )
                @test T / (G * abs(W)) ≈ 0.5 rtol = 1e-8
                @test abs(nd.half_mass_radius(m, res.pos_physical[:, idx]) / rbars[i] - 1) ≤
                      tol_r(r.n[i])
                (x, v)
            end
            @test sqrt(sum(abs2, centres[1][1] - centres[2][1])) ≈ 10.0 rtol = 1e-10
            @test sqrt(sum(abs2, centres[1][2] - centres[2][2])) ≈ v_apo rtol = 1e-10
            orbit = ic["orbit"]
            keys_orbit =
                ("apocentre_pc", "pericentre_pc", "period_myr", "v_peri_kms", "e_orb", "l_orb")
            if orbit_ref === nothing
                orbit_ref = orbit
                t_star_ref = ic["meta"]["t_star_myr"]
            else
                for k in keys_orbit
                    @test orbit[k] ≈ orbit_ref[k] rtol = 1e-10
                end
                # The length unit is the configuration's half-mass radius: equal within its scatter
                @test abs(ic["meta"]["t_star_myr"] / t_star_ref - 1) ≤ 3 * tol_r(256)
            end
            @test ic["meta"]["M_total"] ≈ 1.0e5 rtol = 1e-12
        end
        # The split must be exact: 3:1 takes multiples of 4
        @test_throws r"multiple of 4" n_ladder(
            pair,
            joinpath(dir, "pair_bad");
            n_total = [258],
            seeds = [[1]],
            sweep = sweep,
        )
    end

    @testset "Refusals" begin
        base_text = read(single, String)
        variant(name, from, to) =
            (p = joinpath(dir, name); write(p, replace(base_text, from => to)); p)
        ladder(path; kw...) = n_ladder(
            path,
            joinpath(dir, "refused");
            n_total = get(kw, :n_total, [256]),
            seeds = get(kw, :seeds, [[1]]),
            sweep = get(kw, :sweep, sweep),
        )
        # Only equal-mass bodies given by the cluster mass
        kroupa = variant("kroupa.toml", "imf = \"equal\"\nmass_total = 1.0e5", "imf = \"kroupa\"")
        @test_throws r"merger.cluster1.imf must be \"equal\"" ladder(kroupa)
        by_body = variant(
            "by_body.toml",
            "imf = \"equal\"\nmass_total = 1.0e5",
            "imf = { type = \"equal\", particle_mass = 1.0 }",
        )
        @test_throws r"merger.cluster1.imf.particle_mass fixes the body mass" ladder(by_body)
        # Both or neither of the body mass and the cluster mass: refused on load
        both = variant(
            "both.toml",
            "imf = \"equal\"\nmass_total = 1.0e5",
            "imf = { type = \"equal\", particle_mass = 1.0, target_mass = 1.0e5 }",
        )
        @test_throws r"particle_mass and .*target_mass are mutually exclusive" load_merger_config(
            both,
        )
        neither = variant(
            "neither.toml",
            "imf = \"equal\"\nmass_total = 1.0e5",
            "imf = { type = \"equal\" }",
        )
        @test_throws r"particle_mass or .*target_mass is required" load_merger_config(neither)
        flat_none = variant("flat_none.toml", "mass_total = 1.0e5", "")
        @test_throws r"merger.cluster1.mass_total is required for imf = \"equal\"" load_merger_config(
            flat_none,
        )
        # Stellar evolution, primordial binaries, an IMBH, the complex mode
        evolving = variant("evolving.toml", "kz19 = 0\nlevel = \"0\"", "kz19 = 3")
        @test_throws r"merger.stellar.kz19 must be 0" ladder(evolving)
        binaries = variant(
            "binaries.toml",
            "[merger.output]",
            "[merger.cluster1.binaries]\nfraction = 0.1\n\n[merger.output]",
        )
        @test_throws r"merger.cluster1.binaries must be absent" ladder(binaries)
        imbh = variant(
            "imbh.toml",
            "mass_total = 1.0e5",
            "mass_total = 1.0e5\nimbh = { mass_fraction = 0.01 }",
        )
        @test_throws r"merger.cluster1.imbh must be absent" quiet(() -> ladder(imbh))
        # A rung the loader refuses names the rung and the key
        capped = variant(
            "capped.toml",
            "orbit_mode = \"single\"",
            "orbit_mode = \"single\"\nvirial_max_n = 1000",
        )
        @test_throws r"rung n_total = 2048: .*merger.virial_max_n" ladder(
            capped;
            n_total = [512, 2048],
            seeds = [[1], [2]],
        )
        @test !isdir(joinpath(dir, "refused"))   # nothing written before every rung loads
        # Arguments
        @test_throws r"seeds holds 1 lists for 2 rungs" ladder(single; n_total = [256, 512])
        @test_throws r"n_total must be distinct" ladder(
            single;
            n_total = [256, 256],
            seeds = [[1], [2]],
        )
        @test_throws r"seeds\[1\] must be a nonempty list" ladder(single; seeds = [[1, 1]])
        @test_throws r"sweep.pipeline_config is required" ladder(single; sweep = Dict{String,Any}())
        @test_throws r"sweep.pipeline_config not found" ladder(
            single;
            sweep = Dict{String,Any}("pipeline_config" => "missing.toml"),
        )
        @test_throws r"sweep.seeds is set by the ladder" ladder(
            single;
            sweep = merge(sweep, Dict("seeds" => [1])),
        )
        @test_throws r"sweep.arms.alt.\"merger.cluster1.N\"" ladder(
            single;
            sweep = merge(sweep, Dict("arms" => Dict("alt" => Dict("merger.cluster1.N" => 100)))),
        )
    end

    @testset "Shipped V3-B set" begin
        # The rung files under input_files/campaign/validation/v3_b are this call's output
        src = joinpath(nd._PACKAGE_ROOT, "input_files", "campaign", "validation", "v3_b")
        tmp = mktempdir()
        cp(joinpath(src, "v3_b_pipeline.toml"), joinpath(tmp, "v3_b_pipeline.toml"))
        rungs = quiet(
            () -> n_ladder(
                joinpath(src, "v3_b_king.toml"),
                tmp;
                n_total = [25_000, 50_000, 100_000, 200_000, 400_000, 800_000],
                seeds = [3211:3218, 3221:3224, 3231:3234, 3241:3242, [3251], [3261]],
                sweep = Dict{String,Any}(
                    "pipeline_config" => "v3_b_pipeline.toml",
                    "runs_dir" => "../../../../runs/campaign/v3_b",
                    "concurrency" => 8,
                    "omp_threads" => 8,
                    "poll_interval" => 10.0,
                    "devices" => [0, 1],
                    "engines_per_device" => 4,
                    "arms" => Dict(
                        "base" => Dict{String,Any}(),
                        "twin" => Dict{String,Any}("merger.output.twin_amplitude" => 1.0e-10),
                    ),
                ),
            ),
        )
        for r in rungs, path in (r.merger_config, r.sweep_config)
            shipped = joinpath(src, basename(path))
            @test isfile(shipped)
            @test TOML.parsefile(shipped) == TOML.parsefile(path)
        end
        @test sum(length(r.seeds) for r in rungs) == 20
        top = load_merger_config(joinpath(src, "v3_b_king_n800000.toml"))
        @test top.clusters[1].imf.particle_mass == 0.125 && top.stellar.kz19 == 0
    end
end
