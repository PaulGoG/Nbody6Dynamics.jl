# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Block collector (roadmap 8.6): the paired comparison, the per-arm ensembles and the
# block tables over an explicit run list. A synthetic sweep with arms — a reference,
# a perturbed twin, a null arm and an arm with an injected shift, eight seeds — carries
# Lagrangian radii only (`lagr.7`) and the pairing record (`merger_ic.toml`); the run
# list of its points must reproduce the sweep's table exactly, the shift must be found
# and the null arm accepted, and the ledger front-end must resolve and refuse runs.

using StableRNGs: StableRNG

"""
    _synthetic_arm_sweep(work; seeds, shift, σ) -> (sdir, pts)

A sweep of the arms `base`, `null`, `shifted` and `twin` over `seeds` at one grid
point, every point done, with a synthetic `lagr.7` per run (t = 0…10 in steps of 0.1,
identity scaling): the reference radius of a seed, plus — for every other arm — a
per-seed offset of standard deviation `σ` and AR(1) noise (φ = 0.8, innovations
`σ/5`); the arm `shifted` adds `shift` to every radius. `merger_ic.toml` holds the
seed's `dat10_sha256`; the twin's carries `[twin]` and its own hash.
"""
function _synthetic_arm_sweep(work::AbstractString; seeds = 1:8, shift = 0.05, σ = 0.01)
    write(
        joinpath(work, "config.toml"),
        "[install]\nenabled = false\ninstall_dir = \"engine\"\n\n[simulation]\n" *
        "input_file = \"unused.inp\"\nruns_dir = \"runs\"\n\n[postprocess]\nenabled = false\n\n" *
        "[visualization]\nenabled = false\n\n[merger]\nenabled = false\nconfig_file = \"\"\n",
    )
    write(
        joinpath(work, "merger.toml"),
        "[merger]\nseed = 5\nn_clusters = 1\norbit_mode = \"single\"\n\n[merger.cluster1]\n" *
        "N = 400\nrbar = 1.0\nmodel = \"plummer\"\n\n[merger.output]\ntcrit = 10.0\n" *
        "dtadj = 1.0\ndeltat = 1.0\n",
    )
    write(
        joinpath(work, "sweep.toml"),
        "[sweep]\nname = \"blocks\"\npipeline_config = \"config.toml\"\n" *
        "merger_config = \"merger.toml\"\nseeds = $(collect(seeds))\n\n[sweep.grid]\n" *
        "\"merger.cluster1.rbar\" = [1.0]\n\n[sweep.arms.base]\n\n[sweep.arms.null]\n" *
        "\"pipeline.simulation.omp_threads\" = 4\n\n[sweep.arms.shifted]\n" *
        "\"pipeline.simulation.omp_threads\" = 2\n\n[sweep.arms.twin]\n" *
        "\"merger.output.twin_amplitude\" = 1e-10\n",
    )
    cfg = load_sweep_config(joinpath(work, "sweep.toml"))
    sdir, _ = prepare_sweep(cfg; sweep_dir = joinpath(work, "sweep"))
    pts = sweep_points(cfg)
    rng = StableRNG(860)
    t = collect(0.0:0.1:10.0)
    fractions = Nbody6Dynamics.LAGR_MASS_FRACTIONS
    reference = Dict(s => (1.0 + 0.1 * randn(rng)) .* (1.0 .+ 0.05 .* t) for s in seeds)
    status = Dict{Int,Dict{String,Any}}()
    for p in pts
        out = joinpath(sdir, p.id, "run", "output")
        mkpath(out)
        r = copy(reference[p.seed])
        if p.arm != "base"
            noise = zeros(length(t))
            for k in 2:length(t)
                noise[k] = 0.8 * noise[k - 1] + (σ / 5) * randn(rng)
            end
            r .+= σ * randn(rng) .+ noise .+ (p.arm == "shifted" ? shift : 0.0)
        end
        open(joinpath(out, "lagr.7"), "w") do io
            println(io, "TIME  ", join(string.(fractions), "  "))
            for (k, tk) in enumerate(t)
                println(io, tk, "  ", join(string.(r[k] .* (0.2 .+ 2 .* fractions)), "  "))
            end
        end
        twin = p.arm == "twin" ? "\n[twin]\nbody = 1\namplitude = 1e-10\n" : ""
        hash = p.arm == "twin" ? "t$(p.seed)" : "h$(p.seed)"
        write(joinpath(out, "merger_ic.toml"), "[meta]\ndat10_sha256 = \"$(hash)\"\n$(twin)")
        status[p.index] =
            Dict{String,Any}("status" => "done", "exit_status" => 0, "elapsed_seconds" => 1.0)
    end
    write_sweep_index(sdir, cfg, pts; status = status)
    return sdir, pts
end

@testset "Block collector" begin
    nd = Nbody6Dynamics
    work = mktempdir()
    seeds = 1:8
    sdir, pts = _synthetic_arm_sweep(work; seeds)
    obs = ["r10", "r50"]
    margins = Dict("r10" => ("chaos", 3.0), "r50" => ("chaos", 3.0))
    kw = (; observables = obs, t_settle_myr = 1.0, margins)
    grid = "merger.cluster1.rbar=1.0"
    by_hand =
        [(; arm = p.arm, seed = p.seed, run_dir = joinpath(sdir, p.id, "run"), grid) for p in pts]

    @testset "run list of a sweep: identical rows" begin
        ref_rows = paired_comparison(sdir; kw...)
        @test length(ref_rows) == 2 * 3
        @test isequal(block_comparison(by_hand; kw...), ref_rows)
        listed = sweep_runs(sdir)
        @test length(listed) == length(pts)
        @test Set((r.arm, r.seed, r.run_dir, r.grid) for r in listed) ==
              Set((r.arm, r.seed, r.run_dir, r.grid) for r in by_hand)
        @test all(r -> r.run_id in [p.id for p in pts], listed)
        @test isequal(block_comparison(listed; kw...), ref_rows)
        # The default margin of a table with a twin arm: twice the chaos floor
        @test isequal(block_comparison(by_hand), paired_comparison(sdir))
    end

    @testset "injected shift, null arm, twin" begin
        rows = block_comparison(by_hand; kw...)
        row(o, a) = only(filter(r -> r.observable == o && r.arm == a, rows))
        for o in obs
            @test row(o, "shifted").verdict == "different"
            @test row(o, "shifted").ci_low > 0
            @test row(o, "null").verdict == "equivalent"
            @test row(o, "null").ci_low < 0 < row(o, "null").ci_high
            @test row(o, "twin").n == length(seeds)
            @test isfinite(row(o, "twin").chaos_floor)
            @test row(o, "null").margin ≈ 3 * row(o, "null").chaos_floor
        end
        # r50 carries the radius times (0.2 + 2 × 0.5): the shift scales with it
        @test isapprox(row("r50", "shifted").shift, 0.05 * 1.2; rtol = 0.25)
        # Without the grid field every run belongs to the grid point ""
        no_grid = [(; arm = r.arm, seed = r.seed, run_dir = r.run_dir) for r in by_hand]
        rows0 = block_comparison(no_grid; kw...)
        @test all(r -> r.grid == "", rows0)
        @test [r.verdict for r in rows0] == [r.verdict for r in rows]
        @test [r.shift for r in rows0] == [r.shift for r in rows]
    end

    @testset "refusals" begin
        bad_dir = joinpath(work, "mismatch")
        cp(by_hand[1].run_dir, bad_dir)               # a base run of seed 1 …
        write(joinpath(bad_dir, "output", "merger_ic.toml"), "[meta]\ndat10_sha256 = \"other\"\n")
        mixed =
            [r.arm == "null" && r.seed == 1 ? merge(r, (; run_dir = bad_dir)) : r for r in by_hand]
        err = try
            block_comparison(mixed; kw...)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError &&
              occursin("dat.10 hashes differ", err.msg) &&
              occursin("null run of seed 1", err.msg)
        @test_throws ArgumentError block_comparison(NamedTuple[])
        @test_throws ArgumentError block_comparison([
            merge(by_hand[1], (; run_dir = joinpath(work, "absent"))),
        ])
        @test_throws ArgumentError block_comparison([by_hand[1], by_hand[1]])     # duplicate
        @test_throws ArgumentError block_comparison([(;
            arm = "base",
            run_dir = by_hand[1].run_dir,
        )])
        @test_throws ArgumentError block_comparison([(;
            arm = "base",
            seed = 1.5,
            run_dir = by_hand[1].run_dir,
        )])
        @test_throws ArgumentError block_comparison(by_hand; observables = ["r42"])
        @test_throws ArgumentError block_comparison(filter(r -> r.arm != "base", by_hand))
    end

    @testset "ensembles per arm" begin
        ens = block_ensembles(by_hand; observables = obs, n_grid = 50)
        @test length(ens) == 2 * 4
        @test [e.arm for e in ens] == repeat(["base", "null", "shifted", "twin"], 2)
        @test all(e -> e.stats.n == length(seeds) && length(e.stats.time) == 50, ens)
        base_r50 = only(filter(e -> e.observable == "r50" && e.arm == "base", ens)).stats
        members = [
            nd._run_series(r.run_dir, :lagrangian; fraction = 0.5) for
            r in by_hand if r.arm == "base"
        ]
        direct =
            ensemble_statistics([(Float64.(m[1]), Float64.(m[2])) for m in members]; n_grid = 50)
        @test base_r50.median == direct.median && base_r50.q975 == direct.q975
        shifted_r50 = only(filter(e -> e.observable == "r50" && e.arm == "shifted", ens)).stats
        @test isapprox(sum(shifted_r50.median .- base_r50.median) / 50, 0.06; atol = 0.03)
    end

    @testset "ledger front-end and block tables" begin
        root = joinpath(work, "campaign", "A")
        mkpath(root)
        ledger = joinpath(work, "ledger.csv")
        open(ledger, "w") do io
            println(
                io,
                "block,arm,seed,run_id,config,engine_tag,package_commit,n_bodies,span_myr,span_nb," *
                "mode,submitted,job_ids,segments,status,wall_h,device_h,disk_gb,pulled,analysed," *
                "accepted,note",
            )
            for r in by_hand
                id = "A_$(r.arm)_$(r.seed)"
                cp(r.run_dir, joinpath(root, id))
                println(
                    io,
                    "A,$(r.arm),$(r.seed),$(id),cfg.toml,v2026.07.24,v0.8.0,400,10,10,chain,,,1," *
                    "completed,,,,,,,\"note, with a comma\"",
                )
            end
            println(io, "A,base,9,A_base_9,cfg.toml,,,400,10,10,chain,,,0,running,,,,,,,")
            println(io, "B,base,1,B_base_1,cfg.toml,,,400,10,10,chain,,,1,completed,,,,,,,")
        end
        runs = block_runs(ledger, root; block = "A")
        @test length(runs) == length(by_hand)
        @test all(r -> isdir(r.run_dir) && r.run_dir == joinpath(root, r.run_id), runs)
        @test Set((r.arm, r.seed) for r in runs) == Set((r.arm, r.seed) for r in by_hand)
        rows_ledger = block_comparison(runs; kw...)
        rows_sweep = paired_comparison(sdir; kw...)
        strip_grid(rows) = [Base.structdiff(r, NamedTuple{(:grid,)}) for r in rows]
        @test isequal(strip_grid(rows_ledger), strip_grid(rows_sweep))

        # Refusals: a selected run without a directory, an empty selection, the header
        @test_throws ArgumentError block_runs(ledger, root; block = "A", statuses = ["running"])
        err = try
            block_runs(ledger, root; block = "A", statuses = ["completed", "running"])
            nothing
        catch e
            e
        end
        @test err isa ArgumentError &&
              occursin("A_base_9", err.msg) &&
              occursin("no directory", err.msg)
        @test_throws ArgumentError block_runs(ledger, root; block = "C")
        @test_throws ArgumentError block_runs(joinpath(work, "absent.csv"), root; block = "A")
        bad = joinpath(work, "bad_ledger.csv")
        write(bad, "block,arm,run_id,status\nA,base,A_base_1,completed\n")
        err = try
            block_runs(bad, root; block = "A")
            nothing
        catch e
            e
        end
        @test err isa ArgumentError && occursin("no column seed", err.msg)
        write(bad, "block,arm,seed,run_id,status\nA,base,one,A_base_1,completed\n")
        @test_throws ArgumentError block_runs(bad, root; block = "A")
        write(bad, "block,arm,seed,run_id,status\nA,base,1,A_base_1\n")
        @test_throws ArgumentError block_runs(bad, root; block = "A")

        # Block tables: comparison, ensembles and the member sidecar; #k backups
        out = joinpath(work, "tables")
        paths = write_block_tables(out, runs; kw...)
        @test all(isfile, values(paths))
        cmp_lines = readlines(paths.comparison)
        @test length(cmp_lines) == 1 + 2 * 3
        @test startswith(cmp_lines[1], "grid,observable,arm,n,shift,se")
        ens_lines = readlines(paths.ensembles)
        @test ens_lines[1] == "grid,observable,arm,n,time_myr,median,q16,q84,q025,q975"
        @test length(ens_lines) == 1 + 2 * 4 * 200
        mem_lines = readlines(paths.members)
        @test mem_lines[1] == "arm,seed,grid,run_id,run_dir"
        @test length(mem_lines) == 1 + length(runs)
        @test split(mem_lines[2], ",")[4] == runs[1].run_id
        write_block_tables(out, runs; kw...)
        @test isfile(joinpath(out, "block_comparison#1.csv")) &&
              isfile(joinpath(out, "block_ensembles#1.csv")) &&
              isfile(joinpath(out, "block_members#1.csv"))
    end
end
