# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Calibration report (roadmap 4.6): the time-resolved mean paired difference with its
# standard error, the LaTeX calibration table and the two figures per observable. The
# synthetic runs are exact: radii linear in time, so an injected constant shift must be
# recovered at every node of the common grid and the standard error must be that of
# the per-seed offsets.

"""
    _write_lagr_run(dir, t, r; hash, twin = false) -> String

A run directory with a `lagr.7` whose radius of mass fraction f is `r .* (0.2 + 2f)`
at the times `t`, and a `merger_ic.toml` with the `dat10_sha256` `hash` (and `[twin]`
when `twin`). Returns `dir`.
"""
function _write_lagr_run(dir, t, r; hash, twin = false)
    out = joinpath(dir, "output")
    mkpath(out)
    fractions = Nbody6Dynamics.LAGR_MASS_FRACTIONS
    open(joinpath(out, "lagr.7"), "w") do io
        println(io, "TIME  ", join(string.(fractions), "  "))
        for (k, tk) in enumerate(t)
            println(io, tk, "  ", join(string.(r[k] .* (0.2 .+ 2 .* fractions)), "  "))
        end
    end
    write(
        joinpath(out, "merger_ic.toml"),
        "[meta]\ndat10_sha256 = \"$(hash)\"\n" * (twin ? "\n[twin]\nbody = 1\n" : ""),
    )
    return dir
end

@testset "Calibration report" begin
    nd = Nbody6Dynamics
    work = mktempdir()
    seeds = 1:4
    Δ = 0.05                                  # injected shift of the radius
    c = [0.01, 0.03, -0.02, 0.02]             # per-seed offsets of the arm `offset`
    runs = NamedTuple[]
    for s in seeds
        t = collect(0.0:0.1:(10.0 - 0.5s))    # seeds end at 9.5, 9, 8.5, 8
        r = (1 + 0.1s) .* (1 .+ 0.05 .* t)
        for (arm, rr, twin) in (
            ("base", r, false),
            ("shift", r .+ Δ, false),
            ("offset", r .+ c[s], false),
            ("twin", r .+ 0.5 .* c[s], true),
        )
            dir = _write_lagr_run(
                joinpath(work, "$(arm)_$(s)"),
                t,
                rr;
                hash = twin ? "t$(s)" : "h$(s)",
                twin,
            )
            push!(runs, (; arm, seed = s, run_dir = dir))
        end
    end
    obs = ["r10", "r50"]

    @testset "time-resolved paired differences" begin
        series = paired_difference_series(runs; observables = obs, n_grid = 41)
        @test length(series) == 2 * 3
        @test [s.arm for s in series] == repeat(["offset", "shift", "twin"], 2)
        get_s(o, a) = only(filter(s -> s.observable == o && s.arm == a, series))
        for (o, k) in (("r10", 0.4), ("r50", 1.2))   # 0.2 + 2f at f = 0.1, 0.5
            sh = get_s(o, "shift")
            @test sh.n == length(seeds) && length(sh.time) == 41
            @test first(sh.time) == 0.0 && last(sh.time) ≈ 8.0     # the shortest pair
            @test all(isapprox.(sh.mean, k * Δ; atol = 1e-12))
            @test all(<(1e-12), sh.se)
            off = get_s(o, "offset")
            m = sum(c) / length(c)
            se = sqrt(sum(abs2, c .- m) / (length(c) - 1)) / sqrt(length(c))
            @test all(isapprox.(off.mean, k * m; atol = 1e-12))
            @test all(isapprox.(off.se, k * se; atol = 1e-12))
        end
        # A single seed: the mean without an error
        one_seed = paired_difference_series(filter(r -> r.seed == 1, runs); observables = ["r50"])
        @test all(s -> s.n == 1 && all(isnan, s.se), one_seed)
        # Arms without a series of the observable: rows with empty vectors
        @test all(
            s -> s.n == 0 && isempty(s.time),
            paired_difference_series(runs; observables = ["energy"]),
        )
        @test_throws ArgumentError paired_difference_series(runs; n_grid = 1)
        @test_throws ArgumentError paired_difference_series(runs; observables = ["r42"])
        @test_throws ArgumentError paired_difference_series(runs; reference = "absent")
    end

    @testset "LaTeX calibration table" begin
        cells = nd._tex_shift_cells
        @test cells(-3.43285e-6, 1.24515e-6, -1.12944e-5, 4.42871e-6) ==
              (raw"$(-3.4 \pm 1.2) \times 10^{-6}$", raw"$[-11.3, 4.4] \times 10^{-6}$")
        @test cells(-0.102147, 0.107495, -0.780842, 0.576548) ==
              (raw"$-0.10 \pm 0.11$", raw"$[-0.78, 0.58]$")
        @test cells(0.0, 0.0, 0.0, 0.0) == (raw"$0 \pm 0$", raw"$[0, 0]$")
        @test cells(NaN, NaN, NaN, NaN) == ("--", "--")
        # A standard error at round-off level: three significant digits of the shift
        @test cells(0.02, 2e-18, 0.02 - 1e-17, 0.02 + 1e-17) ==
              (raw"$0.0200 \pm 0.0000$", raw"$[0.0200, 0.0200]$")
        @test cells(0.0012, 0.104, -0.2, 0.2) == (raw"$0.00 \pm 0.10$", raw"$[-0.20, 0.20]$")
        @test cells(0.012, NaN, NaN, NaN) == (raw"$0.0120$", "--")
        @test nd._tex_sig(0.304041) == "0.304" && nd._tex_sig(2.12017) == "2.12"
        @test nd._tex_sig(3.52181e-6) == raw"3.52 \times 10^{-6}"
        @test nd._tex_sig(1e-5) == "10^{-5}" && nd._tex_sig(12345) == raw"1.23 \times 10^{4}"
        @test nd._latex_escape("a_b&c%~") == raw"a\_b\&c\%\textasciitilde{}"

        rows = block_comparison(runs; observables = obs)
        path = write_calibration_table(joinpath(work, "calibration.tex"), rows)
        lines = readlines(path)
        @test lines[1] == raw"\begin{tabular}{llrllll}"
        @test lines[2] == raw"\toprule" && lines[4] == raw"\midrule"
        @test occursin(raw"$90\,\%$ interval", lines[3]) &&
              occursin(raw"Margin $\varepsilon$", lines[3])
        @test lines[end - 1] == raw"\bottomrule" && lines[end] == raw"\end{tabular}"
        body = lines[5:(end - 2)]
        @test count(==(raw"\addlinespace"), body) == 1
        data = filter(l -> l != raw"\addlinespace", body)
        @test length(data) == length(rows)
        @test all(l -> endswith(l, " \\\\") && count("&", l) == 6, data)
        @test startswith(data[1], raw"$r_{10\,\%}$ [pc] & offset & 4 & ")
        @test startswith(data[2], " & shift & 4 & ")              # the observable once per group
        @test occursin("different", only(filter(l -> occursin("& shift &", l), data[1:3])))
        # Escaping, several grid points, α, backups
        odd = [merge(rows[1], (; arm = "fast_io", grid = "g=1")), merge(rows[2], (; grid = "g=2"))]
        write_calibration_table(path, odd; α = 0.025)
        lines = readlines(path)
        @test isfile(joinpath(work, "calibration#1.tex"))
        @test lines[1] == raw"\begin{tabular}{lllrllll}"
        @test occursin(raw"Grid point & Observable", lines[3]) &&
              occursin(raw"$95\,\%$ interval", lines[3])
        @test occursin(raw"fast\_io", lines[5]) && startswith(lines[5], "g=1 & ")
        @test_throws ArgumentError write_calibration_table(path, NamedTuple[])
        @test_throws ArgumentError write_calibration_table(path, rows; α = 0.5)
    end

    @testset "figures (smoke)" begin
        vis = VisualizationConfig(;
            output_dir = joinpath(TESTDIR, "calibration_plots"),
            format = "png",
            dpi = 72,
        )
        rows = block_comparison(runs; observables = obs, t_settle_myr = 1.0)
        series = paired_difference_series(runs; observables = obs)
        paths = calibration_figures(rows, series, vis; t_settle_myr = 1.0)
        @test length(paths) == 4 && all(isfile, paths)
        @test basename.(paths) == [
            "paired_difference_r10.png",
            "paired_shift_r10.png",
            "paired_difference_r50.png",
            "paired_shift_r50.png",
        ]
        @test isfile(
            plot_paired_shifts(
                rows,
                vis;
                observable = "r10",
                holm = false,
                filename = "shift_plain",
            ),
        )
        # Arms identical to the reference are named, not drawn; a table without a margin
        same =
            [(; arm = "copy", seed = r.seed, run_dir = r.run_dir) for r in runs if r.arm == "base"]
        both = vcat(filter(r -> r.arm == "base", runs), same)
        rows0 =
            block_comparison(both; observables = ["r50"], margins = Dict("r50" => ("chaos", 2.0)))
        @test only(rows0).shift == 0 && isnan(only(rows0).margin)   # no twin: no chaos floor
        @test only(rows0).verdict == "inconclusive"
        @test isnan(only(rows0).p_equivalence)                      # an undefined margin tests nothing
        # A single seed per arm: no statistics; the Holm axis of the shift figure shows no value
        rows1 = block_comparison(filter(r -> r.seed == 1, runs); observables = ["r50"])
        @test all(r -> r.n == 1 && isnan(r.p_equivalence_holm), rows1)
        @test isfile(plot_paired_shifts(rows1, vis; filename = "single_seed_shift"))
        @test isfile(
            plot_paired_differences(
                paired_difference_series(both; observables = ["r50"]),
                rows0,
                vis;
                filename = "identical",
            ),
        )
        @test isfile(plot_paired_shifts(rows0, vis; filename = "identical_shift"))
        @test_throws ArgumentError plot_paired_shifts(rows, vis; observable = "energy")
        @test_throws ArgumentError plot_paired_shifts(rows, vis; grid = "absent")
    end
end
