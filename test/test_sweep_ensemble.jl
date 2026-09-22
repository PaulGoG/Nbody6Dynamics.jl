# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Parameter sweeps" begin
    work = mktempdir()
    fixtures = joinpath(@__DIR__, "fixtures")
    base_pipeline = joinpath(work, "config.toml")
    write(
        base_pipeline,
        """
[install]
enabled = false
install_dir = "backend/Nbody6PPGPU-beijing"

[simulation]
run_test = true
input_file = "../../input_files/N1k_quick.inp"
runs_dir = "runs"
omp_threads = 8
run_id_prefix = "sw"

[postprocess]
enabled = true

[visualization]
enabled = false
format = "png"
column = "double"
output_dir = "plots"

[merger]
enabled = false
config_file = ""
""",
    )
    base_merger = joinpath(work, "merger.toml")
    write(
        base_merger,
        replace(
            read(joinpath(@__DIR__, "..", "input_files", "merger_demo_small.toml"), String),
            "[merger]\n" => "[merger]\nseed = 7\n";
            count = 1,
        ),
    )
    sweep_toml(extra = "") = """
[sweep]
name = "unit"
pipeline_config = "config.toml"
merger_config = "merger.toml"
seeds = [1, 2]
concurrency = 2
omp_threads = 3
$(extra)
[sweep.grid]
"merger.orbit.eccentricity" = [0.0, 0.5]
"merger.cluster2.N" = [300, 400]
"""
    spath = joinpath(work, "sweep.toml")
    write(spath, sweep_toml())
    scfg = load_sweep_config(spath)
    @test scfg.name == "unit"
    @test scfg.pipeline_config == base_pipeline && scfg.merger_config == base_merger
    @test first.(scfg.grid) == ["merger.cluster2.N", "merger.orbit.eccentricity"]   # key order
    @test last.(scfg.grid) == [[300, 400], [0.0, 0.5]]
    @test scfg.seeds == [1, 2] && scfg.concurrency == 2 && scfg.omp_threads == 3
    @test scfg.runs_dir == joinpath(work, "runs") && scfg.poll_interval == 2.0

    pts = sweep_points(scfg)
    @test length(pts) == 8                       # 2 × 2 grid × 2 seeds
    @test [p.index for p in pts] == 1:8
    @test pts[1].id == "001_cluster2-N=300_orbit-eccentricity=0_seed=1"
    @test pts[2].id == "002_cluster2-N=300_orbit-eccentricity=0_seed=2"
    @test pts[3].id == "003_cluster2-N=400_orbit-eccentricity=0_seed=1"   # first axis fastest
    @test pts[5].id == "005_cluster2-N=300_orbit-eccentricity=0.5_seed=1"
    @test pts[8].values == ["merger.cluster2.N" => 400, "merger.orbit.eccentricity" => 0.5]
    @test pts[8].seed == 2
    @test Nbody6Dynamics._axis_short("merger.cluster1.binaries.fraction") ==
          "cluster1-binaries-fraction"
    @test Nbody6Dynamics._format_axis_value(0.25) == "0.25"
    @test Nbody6Dynamics._format_axis_value(1000) == "1000"
    @test Nbody6Dynamics._format_axis_value(true) == "true"
    @test Nbody6Dynamics._format_axis_value("king W0=6") == "king-W0-6"

    d = Dict{String,Any}(
        "merger" => Dict{String,Any}("orbit" => Dict{String,Any}("apocentre" => 1.0)),
    )
    Nbody6Dynamics._set_nested!(d, "merger.orbit.eccentricity", 0.3)
    @test d["merger"]["orbit"]["eccentricity"] == 0.3
    @test_throws ArgumentError Nbody6Dynamics._set_nested!(d, "merger.tidal.kz14", 2)

    # Validation
    bad(extra_or_text) =
        (write(spath, extra_or_text); @test_throws ErrorException load_sweep_config(spath))
    bad(replace(sweep_toml(), "name = \"unit\"" => "name = \"unit sweep\""))
    bad(replace(sweep_toml(), "seeds = [1, 2]" => "seeds = [1, 1]"))
    bad(replace(sweep_toml(), "seeds = [1, 2]" => "seeds = []"))
    bad(replace(sweep_toml(), "concurrency = 2" => "concurrency = 0"))
    bad(replace(sweep_toml(), "omp_threads = 3" => "omp_threads = 0"))
    bad(replace(sweep_toml(), "\"merger.cluster2.N\" = [300, 400]\n" => "\"merger.seed\" = [3]\n"))
    bad(
        replace(
            sweep_toml(),
            "\"merger.cluster2.N\" = [300, 400]\n" => "\"orbit.apocentre\" = [3.0]\n",
        ),
    )
    bad(
        replace(
            sweep_toml(),
            "\"merger.cluster2.N\" = [300, 400]\n" => "\"merger.cluster2.N\" = []\n",
        ),
    )
    bad(
        replace(
            sweep_toml(),
            "merger_config = \"merger.toml\"" => "merger_config = \"absent.toml\"",
        ),
    )
    # An absent grid is a seed ensemble of the base configuration, not an error
    write(spath, replace(sweep_toml(), r"\[sweep\.grid\][\s\S]*" => ""))
    @test isempty(load_sweep_config(spath).grid)
    write(spath, sweep_toml("poll_interval = 0.5\n"))
    scfg = load_sweep_config(spath)
    @test scfg.poll_interval == 0.5

    # Preparation (dry run): derived configs, index, summary, figures
    sdir = joinpath(work, "sweep_unit")
    @test run_sweep(scfg; dry_run = true, sweep_dir = sdir) == sdir
    @test_throws ErrorException prepare_sweep(scfg; sweep_dir = sdir)   # never reuse a sweep dir
    @test isfile(joinpath(sdir, "base_config.toml")) && isfile(joinpath(sdir, "base_merger.toml"))
    m3 = load_merger_config(joinpath(sdir, pts[3].id, "merger.toml"))
    @test m3.clusters[2].N == 400 && m3.orbit.eccentricity == 0.0 && m3.seed == 1
    m8 = load_merger_config(joinpath(sdir, pts[8].id, "merger.toml"))
    @test m8.clusters[2].N == 400 && m8.orbit.eccentricity == 0.5 && m8.seed == 2
    c3 = load_config(joinpath(sdir, pts[3].id, "config.toml"))
    @test !c3.install.enabled
    @test c3.install.install_dir == joinpath(work, "backend", "Nbody6PPGPU-beijing")
    @test c3.simulation.runs_dir == joinpath(sdir, pts[3].id)
    @test c3.simulation.omp_threads == 3 && c3.simulation.run_test && !c3.simulation.monitor
    @test c3.merger.enabled && c3.merger.config_file == joinpath(sdir, pts[3].id, "merger.toml")
    @test c3.visualization.column == "double"     # untouched base settings survive
    idx = read_sweep_index(sdir)
    @test idx["sweep"]["name"] == "unit" && idx["sweep"]["n_points"] == 8
    @test idx["sweep"]["axes"] == ["merger.cluster2.N", "merger.orbit.eccentricity"]
    @test all(p -> p["status"] == "pending", idx["points"])
    @test idx["points"][5]["values"]["merger.orbit.eccentricity"] == 0.5
    @test idx["points"][5]["dir"] == joinpath(sdir, pts[5].id)

    # Mark two points done with real output excerpts and summarise
    status = Dict{Int,Dict{String,Any}}()
    for i in (1, 3)
        run_out = joinpath(sdir, pts[i].id, "run", "output")
        mkpath(run_out)
        cp(joinpath(fixtures, "out1000"), joinpath(run_out, "out1000"))
        cp(joinpath(fixtures, "lagr.7"), joinpath(run_out, "lagr.7"))
        write(
            joinpath(sdir, pts[i].id, "run", "RUN_INFO.toml"),
            "[run]\nelapsed_seconds = 12.5\n\n[[segments]]\nexit_status = 0\n",
        )
        status[i] =
            Dict{String,Any}("status" => "done", "exit_status" => 0, "elapsed_seconds" => 12.5)
    end
    status[2] = Dict{String,Any}("status" => "failed", "exit_status" => 1, "elapsed_seconds" => 3.0)
    write_sweep_index(sdir, scfg, pts; status = status)
    columns, rows = sweep_summary(sdir)
    @test columns[1:5] == ["index", "id", "kind", "control_of", "seed"]
    @test columns[6:7] == idx["sweep"]["axes"]
    @test length(rows) == 8
    @test rows[1]["status"] == "done" &&
          rows[1]["elapsed_seconds"] == 12.5 &&
          rows[1]["exit_status"] == 0
    @test rows[1]["n_final"] > 0 && isfinite(rows[1]["de_final"]) && isfinite(rows[1]["q_final"])
    @test rows[1]["t_final_myr"] > 0
    @test rows[2]["status"] == "failed" && rows[2]["exit_status"] == 1
    @test rows[4]["status"] == "pending" && isnan(rows[4]["de_final"]) && rows[4]["n_final"] == -1
    csv_path = write_sweep_summary(sdir)
    lines = readlines(csv_path)
    @test length(lines) == 9
    @test lines[1] == join(columns, ",")
    @test startswith(
        lines[2],
        "1,001_cluster2-N=300_orbit-eccentricity=0_seed=1,merger,,1,300,0,done,0,true,12.5,",
    )
    @test occursin(",pending,-1,false,,,-1,-1,,", lines[5])   # NaN → empty cells

    vis_sw = sweep_visualization(scfg, sdir)
    @test vis_sw.output_dir == joinpath(sdir, "plots") &&
          vis_sw.format == "png" &&
          vis_sw.column == "double"
    figs = sweep_figures(sdir, vis_sw)
    @test length(figs) == 4 && all(isfile, figs)   # two seeds → comparison and ensemble figures
    @test isfile(
        plot_sweep_lagrangian(
            sdir,
            vis_sw;
            axis = "merger.orbit.eccentricity",
            fraction = 0.1,
            filename = "r10",
        ),
    )
    @test_throws ArgumentError plot_sweep_energy(sdir, vis_sw; axis = "merger.orbit.apocentre")
    rm(work; recursive = true, force = true)

    # Annotation corner: the least occupied of the data's bounding box
    xs = collect(0.0:0.1:1.0)
    @test MakieExt._emptiest_corner(xs, xs) == :tl          # rising series frees the top-left
    @test MakieExt._emptiest_corner(xs, 1 .- xs) == :tr     # falling series frees the top-right
    @test MakieExt._emptiest_corner(Float64[], Float64[]) == :tl
end

# =====================================================================

@testset "Seeded ensembles" begin
    # Type-7 quantiles and linear interpolation
    q = Nbody6Dynamics._quantile_sorted
    v = [1.0, 2.0, 3.0, 4.0, 5.0]
    @test q(v, 0.5) == 3.0 && q(v, 0.0) == 1.0 && q(v, 1.0) == 5.0
    @test q(v, 0.16) ≈ 1.64 && q(v, 0.84) ≈ 4.36
    @test q(v, 0.025) ≈ 1.1 && q(v, 0.975) ≈ 4.9
    @test q([7.0], 0.3) == 7.0
    @test_throws ArgumentError q(Float64[], 0.5)
    @test_throws ArgumentError q(v, 1.5)
    interp = Nbody6Dynamics._interpolate_linear
    @test interp([0.0, 1.0, 2.0], [0.0, 10.0, 0.0], [0.0, 0.5, 1.0, 1.5, 2.0]) ==
          [0.0, 5.0, 10.0, 5.0, 0.0]
    @test interp([1.0], [3.0], [1.0]) == [3.0]
    @test_throws ArgumentError interp([0.0, 1.0], [0.0, 1.0], [1.5])
    @test_throws DimensionMismatch interp([0.0, 1.0], [0.0], [0.5])

    # Members y_i = i t on different spans: common grid = intersection
    members = [(collect(0.0:0.5:10.0), i .* collect(0.0:0.5:10.0)) for i in 1:5]
    members[2] = (collect(1.0:0.5:12.0), 2 .* collect(1.0:0.5:12.0))   # starts later, ends later
    st = ensemble_statistics(members; n_grid = 19)
    @test st.n == 5
    @test st.time[1] == 1.0 && st.time[end] == 10.0 && length(st.time) == 19
    @test st.median ≈ 3 .* st.time
    @test st.q16 ≈ 1.64 .* st.time && st.q84 ≈ 4.36 .* st.time
    @test st.q025 ≈ 1.1 .* st.time && st.q975 ≈ 4.9 .* st.time
    single = ensemble_statistics(members[1:1]; n_grid = 5)
    @test single.n == 1 &&
          single.median == single.q025 == single.q975 ≈ collect(range(0.0, 10.0; length = 5))
    @test_throws ArgumentError ensemble_statistics(typeof(members[1])[])
    @test_throws ArgumentError ensemble_statistics([
        members[1],
        (collect(20.0:1.0:30.0), zeros(11)),
    ])
    @test_throws ArgumentError ensemble_statistics([([1.0, 0.5, 2.0], [0.0, 1.0, 2.0])])
    @test_throws ArgumentError ensemble_statistics(members; n_grid = 1)

    # Series extractors on a run directory built from the fixtures
    fixtures = joinpath(@__DIR__, "fixtures")
    rdir = mktempdir()
    @test Nbody6Dynamics._run_series(rdir, :lagrangian) === nothing
    @test Nbody6Dynamics._run_series(rdir, :energy) === nothing
    mkpath(joinpath(rdir, "output"))
    cp(joinpath(fixtures, "out1000"), joinpath(rdir, "output", "out1000"))
    cp(joinpath(fixtures, "lagr.7"), joinpath(rdir, "output", "lagr.7"))
    t_l, r_l = Nbody6Dynamics._run_series(rdir, :lagrangian; fraction = 0.5)
    @test length(t_l) == length(r_l) > 1 && issorted(t_l) && all(>(0), r_l)
    t_e, de = Nbody6Dynamics._run_series(rdir, :energy)
    @test all(>(0), de) && length(t_e) == length(de)
    t_n, n = Nbody6Dynamics._run_series(rdir, :n_stars)
    @test all(>(0), n)
    t_p, np = Nbody6Dynamics._run_series(rdir, :n_pairs)
    @test length(np) == length(t_n) && all(≥(0), np)
    @test_throws ArgumentError Nbody6Dynamics._run_series(rdir, :unknown)
    @test Nbody6Dynamics._series_label(:energy, 0.5) isa AbstractString
    @test_throws ArgumentError Nbody6Dynamics._series_label(:unknown, 0.5)

    # A seeds-only sweep (no grid axes) and a one-axis sweep, marked done with fixture output
    work = mktempdir()
    write(
        joinpath(work, "config.toml"),
        "[install]\nenabled = false\ninstall_dir = \"backend\"\n\n[simulation]\nrun_test = true\nruns_dir = \"runs\"\n\n[visualization]\nformat = \"png\"\n",
    )
    write(
        joinpath(work, "merger.toml"),
        read(joinpath(@__DIR__, "..", "input_files", "merger_demo_small.toml"), String),
    )
    function _done_sweep(name, grid_text, seeds_text)
        spath = joinpath(work, "sweep_$(name).toml")
        write(
            spath,
            "[sweep]\nname = \"$(name)\"\npipeline_config = \"config.toml\"\nmerger_config = \"merger.toml\"\nseeds = $(seeds_text)\n$(grid_text)",
        )
        scfg = load_sweep_config(spath)
        sdir = joinpath(work, "sweep_$(name)")
        prepare_sweep(scfg; sweep_dir = sdir)
        pts = sweep_points(scfg)
        status = Dict{Int,Dict{String,Any}}()
        for p in pts
            out = joinpath(sdir, p.id, "run", "output")
            mkpath(out)
            cp(joinpath(fixtures, "out1000"), joinpath(out, "out1000"))
            cp(joinpath(fixtures, "lagr.7"), joinpath(out, "lagr.7"))
            status[p.index] =
                Dict{String,Any}("status" => "done", "exit_status" => 0, "elapsed_seconds" => 1.0)
        end
        write_sweep_index(sdir, scfg, pts; status = status)
        return scfg, sdir, pts
    end
    scfg0, sdir0, pts0 = _done_sweep("seeds", "", "[1, 2, 3]")
    @test isempty(scfg0.grid)
    @test [p.id for p in pts0] == ["001_seed=1", "002_seed=2", "003_seed=3"]
    ens0 = sweep_ensembles(sdir0, :lagrangian)
    @test length(ens0) == 1 && ens0[1].stats.n == 3 && isempty(ens0[1].values)
    @test ens0[1].stats.median ≈ ens0[1].stats.q975    # identical members → zero-width bands
    vis_e = VisualizationConfig(;
        format = "png",
        column = "single",
        output_dir = joinpath(work, "plots"),
    )
    @test isfile(plot_sweep_ensemble(sdir0, vis_e))
    @test isfile(plot_sweep_ensemble(sdir0, vis_e; quantity = :energy))
    @test isfile(plot_sweep_ensemble(sdir0, vis_e; quantity = :n_stars, filename = "ens_n"))
    @test isfile(plot_sweep_lagrangian(sdir0, vis_e; filename = "seeds_lagr"))   # no axes: one colour, no legend
    @test_throws ArgumentError plot_sweep_lagrangian(
        sdir0,
        vis_e;
        axis = "merger.orbit.eccentricity",
    )
    figs0 = sweep_figures(sdir0, vis_e)
    @test length(figs0) == 4 && all(isfile, figs0)

    scfg1, sdir1, pts1 = _done_sweep(
        "axes",
        "[sweep.grid]\n\"merger.orbit.eccentricity\" = [0.0, 0.5]\n\"merger.cluster2.N\" = [300, 400]\n",
        "[1, 2]",
    )
    ens1 = sweep_ensembles(sdir1, :energy)
    @test length(ens1) == 4 && all(e -> e.stats.n == 2, ens1)
    @test isfile(plot_sweep_ensemble(sdir1, vis_e; quantity = :lagrangian, filename = "ens_axes"))
    @test isfile(
        plot_sweep_ensemble(
            sdir1,
            vis_e;
            axis = "merger.orbit.eccentricity",
            fixed = Dict("merger.cluster2.N" => 400),
            filename = "ens_fixed",
        ),
    )
    @test_throws ErrorException plot_sweep_ensemble(
        sdir1,
        vis_e;
        axis = "merger.orbit.eccentricity",
        fixed = Dict("merger.cluster2.N" => 999),
    )
    @test_throws ArgumentError plot_sweep_ensemble(
        sdir1,
        vis_e;
        fixed = Dict("merger.cluster2.N" => 400),
    )   # the axis itself
    @test_throws ArgumentError plot_sweep_ensemble(
        sdir1,
        vis_e;
        fixed = Dict("merger.orbit.apocentre" => 1.0),
    )
    @test length(sweep_figures(sdir1, vis_e)) == 4
    rm(work; recursive = true, force = true)
    rm(rdir; recursive = true, force = true)
end

# =====================================================================

@testset "Control configurations" begin
    demo = joinpath(@__DIR__, "..", "input_files", "merger_demo_small.toml")
    raw = Nbody6Dynamics.TOML.parsefile(demo)
    raw["merger"]["cluster2"]["N"] = 500
    raw["merger"]["cluster2"]["rbar"] = 4.0
    raw["merger"]["cluster1"]["binaries"] = Dict{String,Any}("fraction" => 0.1)
    raw["merger"]["nbody6"] = Dict{String,Any}("qe" => 0.01)
    c = control_merger_dict(raw)
    m = c["merger"]
    @test m["n_clusters"] == 1 && m["orbit_mode"] == "explicit"
    @test !haskey(m, "orbit") && !haskey(m, "cluster2")
    @test m["cluster1"]["N"] == 1500                                  # 1000 + 500
    @test m["cluster1"]["rbar"] ≈ (1000 * 2.0 + 500 * 4.0) / 1500     # N-weighted r_h
    @test m["cluster1"]["position"] == [0.0, 0.0, 0.0] &&
          m["cluster1"]["velocity"] == [0.0, 0.0, 0.0]
    @test m["cluster1"]["model"] == "king" && m["cluster1"]["W0"] == 6.0
    @test m["cluster1"]["binaries"]["fraction"] == 0.1                # cluster 1's population kept
    @test m["nbody6"]["qe"] == 0.01 && haskey(m, "output")             # other sections verbatim
    # NB intervals scale by (RBAR_est / r_h,c)^{3/2}: apocentre 12 × 500/1500 + N-weighted r_h
    rbar_c = m["cluster1"]["rbar"]
    factor = ((12.0 * 500 / 1500 + rbar_c) / rbar_c)^1.5
    @test m["output"]["tcrit"] ≈ raw["merger"]["output"]["tcrit"] * factor
    @test m["output"]["deltat"] ≈ raw["merger"]["output"]["deltat"] * factor
    raw_myr = deepcopy(raw)
    delete!(raw_myr["merger"]["output"], "tcrit")
    raw_myr["merger"]["output"]["tcrit_myr"] = 25.0
    raw_myr["merger"]["stellar"] = Dict{String,Any}("dtplot_myr" => 5.0)
    m_myr = control_merger_dict(raw_myr)["merger"]
    @test m_myr["output"]["tcrit_myr"] == 25.0 && !haskey(m_myr["output"], "tcrit")   # physical: verbatim
    @test m_myr["output"]["deltat"] ≈ raw["merger"]["output"]["deltat"] * factor     # NB: scaled
    @test m_myr["stellar"]["dtplot_myr"] == 5.0
    raw_ex = deepcopy(raw)
    raw_ex["merger"]["orbit_mode"] = "explicit"
    delete!(raw_ex["merger"], "orbit")
    raw_ex["merger"]["cluster1"]["position"] = [-6.0, 0.0, 0.0]
    raw_ex["merger"]["cluster2"]["position"] = [12.0, 0.0, 0.0]
    f_ex = Nbody6Dynamics._control_time_factor(raw_ex["merger"], [1000, 500], [2.0, 4.0], rbar_c)
    centre = (1000 * -6.0 + 500 * 12.0) / 1500
    rms = sqrt((1000 * (-6.0 - centre)^2 + 500 * (12.0 - centre)^2) / 1500)
    @test f_ex ≈ ((rms + rbar_c) / rbar_c)^1.5
    @test raw["merger"]["n_clusters"] == 2                            # input untouched
    @test_throws ArgumentError control_merger_dict(Dict{String,Any}())
    bad = deepcopy(raw)
    delete!(bad["merger"], "cluster2")
    @test_throws ArgumentError control_merger_dict(bad)

    work = mktempdir()
    dst = joinpath(work, "control.toml")
    @test write_control_merger_config(demo, dst) == dst
    ccfg = load_merger_config(dst)
    @test length(ccfg.clusters) == 1 && ccfg.clusters[1].N == 2000 && ccfg.orbit_mode == "explicit"
    @test_throws ErrorException write_control_merger_config(joinpath(work, "absent.toml"), dst)

    # The single-cluster generator path: one member block, no truncation, at rest
    small = replace(read(dst, String), r"N = \d+" => "N = 300")
    write(joinpath(work, "control_small.toml"), small)
    res = generate_merger_ic(
        load_merger_config(joinpath(work, "control_small.toml"));
        output_dir = joinpath(work, "ic"),
    )
    @test res.N_total == 300 &&
          length(res.cluster_ranges) == 1 &&
          length(res.cluster_ranges[1]) == 300
    @test isfile(joinpath(work, "ic", "dat.10")) && isfile(joinpath(work, "ic", "merger.inp"))
    @test parse_merger_summary(joinpath(work, "ic", "merger_summary.txt")) == [collect(1:300)]
    # Physical-time intervals: converted with the realised T* at generation
    phys = replace(small, r"tcrit = [0-9.]+" => "tcrit_myr = 20.0")
    phys = replace(phys, r"deltat = [0-9.]+" => "deltat_myr = 4.0")
    write(joinpath(work, "control_phys.toml"), phys)
    cfg_p = load_merger_config(joinpath(work, "control_phys.toml"))
    @test cfg_p.output.tcrit_myr == 20.0 && cfg_p.output.deltat_myr == 4.0
    res_p = generate_merger_ic(cfg_p; output_dir = joinpath(work, "ic_phys"))
    meta_p = Nbody6Dynamics.TOML.parsefile(joinpath(work, "ic_phys", "merger_ic.toml"))
    t_star = meta_p["meta"]["t_star_myr"]
    @test t_star > 0
    @test meta_p["output"]["tcrit"] ≈ 20.0 / t_star && meta_p["output"]["tcrit_myr"] ≈ 20.0
    # The interval is the dyadic rounding of the converted value (engine digit counter)
    @test meta_p["output"]["deltat"] == engine_interval(4.0 / t_star)
    @test abs(meta_p["output"]["deltat"] - 4.0 / t_star) ≤ 4.0 / t_star / 256
    inp_p = read(joinpath(work, "ic_phys", "merger.inp"), String)
    @test occursin(Nbody6Dynamics.Printf.@sprintf("TCRIT=%.2f", 20.0 / t_star), inp_p)
    @test occursin("Time unit: T* =", read(joinpath(work, "ic_phys", "merger_summary.txt"), String))
    both = replace(small, r"tcrit = [0-9.]+" => "tcrit = 1.0\ntcrit_myr = 20.0")
    write(joinpath(work, "both.toml"), both)
    @test_throws ErrorException load_merger_config(joinpath(work, "both.toml"))
    neg = replace(small, r"tcrit = [0-9.]+" => "tcrit_myr = -1.0")
    write(joinpath(work, "neg.toml"), neg)
    @test_throws ErrorException load_merger_config(joinpath(work, "neg.toml"))
    # dtplot below deltat is caught at generation once T* is known
    # (the control file carries a [merger.stellar] table with the scaled dtplot)
    @test occursin(r"dtplot = [0-9.]+", phys)
    late = replace(phys, r"dtplot = [0-9.]+" => "dtplot_myr = 1.0")
    write(joinpath(work, "late.toml"), late)
    @test_throws ErrorException load_merger_config(joinpath(work, "late.toml"))   # both physical: config time
    late_nb = replace(phys, r"dtplot = [0-9.]+" => "dtplot = 0.01")
    write(joinpath(work, "late_nb.toml"), late_nb)
    @test_throws ErrorException generate_merger_ic(
        load_merger_config(joinpath(work, "late_nb.toml"));
        output_dir = joinpath(work, "ic_late"),
    )
    kepler_one = replace(small, "orbit_mode = \"explicit\"" => "orbit_mode = \"kepler\"")
    write(joinpath(work, "kepler_one.toml"), kepler_one)
    @test_throws ErrorException load_merger_config(joinpath(work, "kepler_one.toml"))

    # Sweep with controls: companions interleaved, derived at the point's values
    write(
        joinpath(work, "config.toml"),
        "[install]\nenabled = false\ninstall_dir = \"backend\"\n\n[simulation]\nrun_test = true\nruns_dir = \"runs\"\n\n[visualization]\nformat = \"png\"\n",
    )
    cp(demo, joinpath(work, "merger.toml"))
    write(
        joinpath(work, "sweep.toml"),
        "[sweep]\nname = \"ctrl\"\npipeline_config = \"config.toml\"\nmerger_config = \"merger.toml\"\nseeds = [1]\ncontrols = true\n\n[sweep.grid]\n\"merger.cluster2.N\" = [300, 400]\n",
    )
    scfg = load_sweep_config(joinpath(work, "sweep.toml"))
    @test scfg.controls
    pts = sweep_points(scfg)
    @test length(pts) == 4
    @test [p.kind for p in pts] == ["merger", "control", "merger", "control"]
    @test pts[2].id == "002_cluster2-N=300_seed=1_control" && pts[2].control_of == pts[1].id
    @test pts[2].values == pts[1].values && pts[2].seed == pts[1].seed
    sdir = joinpath(work, "sweep_ctrl")
    prepare_sweep(scfg; sweep_dir = sdir)
    m3 = load_merger_config(joinpath(sdir, pts[3].id, "merger.toml"))
    m4 = load_merger_config(joinpath(sdir, pts[4].id, "merger.toml"))
    @test length(m3.clusters) == 2 && m3.clusters[2].N == 400
    @test length(m4.clusters) == 1 && m4.clusters[1].N == 1400 && m4.seed == 1
    idx = read_sweep_index(sdir)
    @test idx["sweep"]["controls"] == true
    @test [p["kind"] for p in idx["points"]] == ["merger", "control", "merger", "control"]
    @test idx["points"][4]["control_of"] == pts[3].id
    columns, rows = sweep_summary(sdir)
    @test columns[1:5] == ["index", "id", "kind", "control_of", "seed"]
    @test rows[2]["kind"] == "control" && rows[2]["control_of"] == pts[1].id
    @test_throws ErrorException load_sweep_config(
        (
            write(
                joinpath(work, "bad.toml"),
                replace(
                    read(joinpath(work, "sweep.toml"), String),
                    "controls = true" => "controls = \"yes\"",
                ),
            );
            joinpath(work, "bad.toml")
        ),
    )

    # Paired figure on fixture output: mergers solid, controls dashed
    fixtures = joinpath(@__DIR__, "fixtures")
    status = Dict{Int,Dict{String,Any}}()
    for p in pts
        out = joinpath(sdir, p.id, "run", "output")
        mkpath(out)
        cp(joinpath(fixtures, "out1000"), joinpath(out, "out1000"))
        cp(joinpath(fixtures, "lagr.7"), joinpath(out, "lagr.7"))
        status[p.index] =
            Dict{String,Any}("status" => "done", "exit_status" => 0, "elapsed_seconds" => 1.0)
    end
    write_sweep_index(sdir, scfg, pts; status = status)
    _, done_m = Nbody6Dynamics._sweep_done_points(sdir)
    _, done_c = Nbody6Dynamics._sweep_done_points(sdir; kind = "control")
    _, done_all = Nbody6Dynamics._sweep_done_points(sdir; kind = "")
    @test length(done_m) == 2 && length(done_c) == 2 && length(done_all) == 4
    vis_c = VisualizationConfig(;
        format = "png",
        column = "single",
        output_dir = joinpath(work, "plots"),
    )
    @test isfile(plot_control_comparison(sdir, vis_c))
    @test isfile(plot_control_comparison(sdir, vis_c; quantity = :energy, filename = "ctrl_energy"))
    figs = sweep_figures(sdir, vis_c)
    @test length(figs) == 3 && all(isfile, figs)      # one seed: comparison + control figures
    @test length(sweep_ensembles(sdir, :lagrangian)) == 2   # mergers only
    status[2]["status"] = "failed"
    status[4]["status"] = "failed"
    write_sweep_index(sdir, scfg, pts; status = status)
    @test_throws ErrorException plot_control_comparison(sdir, vis_c; filename = "none")
    rm(work; recursive = true, force = true)
end

# =====================================================================
