# Validity record (roadmap 2.4): calibration lookup, stage detectors, sections,
# epochs and the schema of validity.toml / validity_snapshots.csv.

using StableRNGs

"""Snapshot of N-body arrays at time `t` (rbar 2 pc, tscale 3 Myr), no density or potential."""
function _val_snap(pos, vel, mass, t)
    params = zeros(Float32, 20)
    params[1] = t
    params[3] = 2.0f0
    params[4] = 1.0f0
    params[11] = 3.0f0
    params[12] = 1.0f0
    return Snapshot(
        SnapshotHeader(Int32(length(mass)), Int32(1), Int32(1), Int32(20), params),
        Int32.(1:length(mass)),
        Float32.(mass),
        Float32.(pos),
        Float32.(vel),
        Float32[],
        Float32[],
    )
end

@testset "Validity: calibration table" begin
    cal = load_calibration()
    @test cal.version == 1
    @test length(cal.rows) == 1974
    rows = sort(
        [
            r for r in cal.rows if r.model == "plummer" &&
                r.estimator == "core_radius" &&
                r.quantity == "r_core_number"
        ];
        by = r -> r.n,
    )
    r4 = rows[findfirst(r -> r.n == 10_000, rows)]
    c = calibrated(cal, "plummer", "core_radius", "r_core_number", 10_000)
    @test c.relative && c.in_range
    @test c.bias ≈ r4.bias / abs(r4.truth)
    @test c.scatter ≈ r4.scatter / abs(r4.truth)
    @test c.ok == (abs(c.bias) ≤ 0.05)
    # Linear in log10 N between the bracketing grid points
    k = findfirst(r -> r.n == 10_000, rows)
    r5 = rows[k + 1]
    nmid = sqrt(10_000 * r5.n)
    cm = calibrated(cal, "plummer", "core_radius", "r_core_number", nmid)
    @test cm.bias ≈ ((r4.bias + r5.bias) / 2) / abs((r4.truth + r5.truth) / 2)
    # Outside the grid: the end row, flagged
    c_lo = calibrated(cal, "plummer", "core_radius", "r_core_number", 100)
    @test !c_lo.in_range && c_lo.bias ≈ rows[1].bias / abs(rows[1].truth)
    # tolerance 0 masks every nonzero bias
    @test calibrated(cal, "plummer", "core_radius", "r_core_number", 10_000; tolerance = 0.0).ok ==
          (r4.bias == 0)
    # A quantity without a finite nonzero truth stays absolute and is never masked
    rl = [
        r for r in cal.rows if
        r.model == "plummer" && r.estimator == "rotation" && r.quantity == "lambda_r"
    ]
    if !isempty(rl) && !(isfinite(rl[1].truth) && rl[1].truth != 0)
        ca = calibrated(cal, "plummer", "rotation", "lambda_r", rl[1].n)
        @test !ca.relative && ca.ok
    end
    @test_throws ArgumentError calibrated(cal, "plummer", "core_radius", "no_such_quantity", 1000)
    @test_throws ArgumentError calibrated(cal, "plummer", "core_radius", "r_core_number", 0)
    @test_throws ArgumentError load_calibration(mktempdir())
end

@testset "Validity: stage detectors" begin
    holds = Nbody6Dynamics._validity_holds_from
    v = [true, false, true, true]
    @test holds(j -> v[j], 4) == 3
    @test holds(j -> true, 5) == 1
    @test holds(j -> j < 3, 4) === nothing
    @test holds(j -> true, 0) === nothing
    rate = Nbody6Dynamics._validity_rate_test
    t = collect(0.0:1.0:10.0)
    @test rate(t, fill(1.0, 11), 2.0, 0.05)
    step = [j ≤ 5 ? 0.0 : 0.1 for j in 1:11]
    @test !rate(t, step, 2.0, 0.05)
    @test rate(t, step, 2.0, 0.2)
    @test !rate(t, fill(1.0, 11), 20.0, 0.05)   # window longer than the series
end

@testset "Validity: record of a two-cluster approach" begin
    rng = StableRNG(2024)
    n_cl = 600
    p1, v1 = Nbody6Dynamics.sample_plummer(n_cl, 0.2; rng = rng)
    p2, v2 = Nbody6Dynamics.sample_plummer(n_cl, 0.2; rng = rng)
    mass = fill(1.0 / (2n_cl), 2n_cl)
    vel = hcat(v1, v2) .* sqrt(0.5)
    ranges = [1:n_cl, (n_cl + 1):(2n_cl)]
    # Separation along x: approach, a pericentre at t = 2, recession, merger from t = 4
    sep = [4.0, 2.0, 1.0, 2.0, 0.1, 0.1, 0.1]
    snaps = [
        _val_snap(hcat(p1 .- [d / 2, 0, 0], p2 .+ [d / 2, 0, 0]), vel, mass, Float64(k - 1)) for
        (k, d) in enumerate(sep)
    ]
    ints = snapshot_integrals(snaps)
    pp = PostprocessConfig(; msr_n_tracer = 20, msr_n_random = 20)
    bsets = [s.bound for s in ints]
    remn = remnant_diagnostics(snaps, ranges; settings = pp, bound_sets = bsets)
    pb = progenitor_budget(snaps, ranges, bsets; settings = pp)
    # ADJUST records: one interval (t = 2 → 3) above the 1e-4 threshold
    de = [0.0, 1e-6, 2e-6, 5e-4, 1e-6, 1e-6, 1e-6]
    detot = cumsum(de .* 0.25)
    adjust = [
        Nbody6Dynamics.AdjustRecord(
            Float64(k - 1),
            3.0 * (k - 1),
            0.5,
            de[k],
            -0.25,
            2n_cl,
            0,
            0.8,
            de[k] * 0.25,
            detot[k],
        ) for k in eachindex(de)
    ]
    diag = Nbody6Dynamics.DiagnosticsData(adjust, Dict{String,Float64}())
    results = Dict{Symbol,Any}(
        :snapshots => snaps,
        :integrals => ints,
        :diagnostics => diag,
        :remnant => remn,
        :progenitor_budget => pb,
    )
    rec = validity_record(
        results,
        pp;
        build_info = Dict{String,Any}(
            "backend_commit" => "abc1234",
            "source_amendments" => ["escape-summary-format"],
        ),
        run_info = Dict{String,Any}("run" => Dict{String,Any}("id" => "val_test")),
    )
    @test rec.schema_version == Nbody6Dynamics.VALIDITY_SCHEMA_VERSION == 1
    h = rec.sections["header"]
    @test h["run_id"] == "val_test" && h["backend_commit"] == "abc1234"
    @test h["source_amendments"] == ["escape-summary-format"]
    @test h["calibration_model"] == "plummer" && h["calibration_version"] == 1
    @test h["n_snapshots"] == 7 && h["t_start_nb"] == 0.0 && h["t_end_nb"] == 6.0
    @test h["tidal_mode"] == "none"

    e = rec.sections["energy"]
    @test e["available"] && e["engine_energy_valid"]
    @test e["n_intervals"] == 6 && e["n_above_threshold"] == 1
    @test e["max_abs_de_rel"] ≈ 5e-4
    @test e["e0"] == -0.25 && e["detot"] ≈ detot[end]
    @test e["cumulative_rel"] ≈ detot[end] / 0.25
    @test e["verdict"] == "flag" && e["reasons"] == ["1 interval above 1.0e-4"]
    @test isnan(e["budget_final_residual_rel"])

    # Crossing time of the bound set (Heggie & Mathieu 1986) and the violent-relaxation window
    s1 = ints[1]
    tcr1 = Nbody6Dynamics.crossing_time(s1.bound_mass, s1.bound.kinetic + s1.bound.potential)
    @test rec.snapshots[1].t_cr ≈ tcr1
    vr = only(filter(x -> x.kind === :violent_relaxation, rec.epochs))
    @test vr.t_start == 0.0 && vr.t_end ≈ 5 * tcr1
    th = only(filter(x -> x.kind === :energy_threshold, rec.epochs))
    @test th.t_start == 2.0 && th.t_end == 3.0
    @test issorted([x.t_start for x in rec.epochs])

    # The pericentre before coalescence and the coalescence phase
    st = rec.sections["stages"]
    @test st["coalescence_time_nb"] == remn.coalescence_time == 4.0
    @test st["pericentres"] == [2.0]
    pe = only(filter(x -> x.kind === :pericentre, rec.epochs))
    @test pe.t_start < 2.0 < pe.t_end
    cp = only(filter(x -> x.kind === :coalescence_phase, rec.epochs))
    @test cp.t_start ≈ pe.t_start && cp.t_end == 4.0
    @test st["virialised"] in ("reached", "not_reached", "undetermined")
    @test st["mixed"] in ("reached", "not_reached", "undetermined")

    # Per-snapshot rows
    @test length(rec.snapshots) == 7
    windows =
        filter(x -> x.kind in (:violent_relaxation, :pericentre, :coalescence_phase), rec.epochs)
    for (k, s) in enumerate(rec.snapshots)
        @test s.time_nb == k - 1
        @test s.n_bound == length(ints[k].bound.stars)
        @test s.n_eff ≈ s.n_bound          # equal masses
        @test s.in_window == any(w -> w.t_start ≤ s.time_nb ≤ w.t_end, windows)
        @test s.equilibrium_ok == (!s.in_window && s.bound_converged)
        @test s.lagr_f_min ≈ pp.lagr_n_min / s.n_eff && s.lagr_f_max ≈ 1 - pp.lagr_n_min / s.n_eff
        @test length(s.calibration) == length(Nbody6Dynamics.VALIDITY_TRACKED)
    end
    cs = rec.sections["calibration"]
    @test cs["version"] == 1 && cs["n_reference"] == rec.snapshots[1].n_bound
    @test all(haskey(cs, "$(a)__$(b)") for (a, b) in Nbody6Dynamics.VALIDITY_TRACKED)
    @test rec.sections["masks"]["M2_lagr_n_min"] == pp.lagr_n_min
    @test rec.sections["bound_set"]["n_snapshots"] == 7
    @test rec.sections["thresholds"]["interval_de"] == pp.validity_interval_de

    # Schema of the written files
    dir = mktempdir()
    path = write_validity(joinpath(dir, "validity.toml"), rec)
    d = Nbody6Dynamics.TOML.parsefile(path)
    @test d["schema_version"] == 1
    for sec in (
        "header",
        "thresholds",
        "energy",
        "mass",
        "angular_momentum",
        "bound_set",
        "cadence",
        "stages",
        "calibration",
        "masks",
    )
        @test d[sec] isa AbstractDict
    end
    @test d["energy"]["verdict"] == "flag"
    @test d["energy"]["budget_final_residual_rel"] isa Float64 &&
          isnan(d["energy"]["budget_final_residual_rel"])
    @test d["stages"]["pericentres"] == [2.0]
    @test length(d["epochs"]) == length(rec.epochs)
    kinds = Set([
        "violent_relaxation",
        "pericentre",
        "coalescence_phase",
        "energy_interval",
        "energy_threshold",
        "mass_budget",
        "bound_set",
    ])
    @test all(
        ep -> Set(keys(ep)) == Set(["t_start", "t_end", "kind", "detail"]) && ep["kind"] in kinds,
        d["epochs"],
    )
    @test d["calibration"]["core_radius__r_core_number"]["relative"] isa Bool
    write_validity(path, rec)
    @test isfile(joinpath(dir, "validity#1.toml"))
    csv = write_validity_snapshots(joinpath(dir, "validity_snapshots.csv"), rec)
    lines = readlines(csv)
    header = split(lines[1], ',')
    @test length(header) == 10 + 3 * length(Nbody6Dynamics.VALIDITY_TRACKED)
    @test header[1:3] == ["time_nb", "n_bound", "n_eff"] &&
          header[11] == "core_radius__r_core_number_bias"
    @test length(lines) == 8 && all(l -> length(split(l, ',')) == length(header), lines[2:end])

    # Snapshots and integrals that do not correspond
    @test_throws DimensionMismatch validity_record(
        Dict{Symbol,Any}(:snapshots => snaps, :integrals => ints[1:2]),
        pp,
    )
end

@testset "Validity: empty results and configuration bounds" begin
    rec = validity_record(Dict{Symbol,Any}(), PostprocessConfig())
    @test rec.sections["energy"]["verdict"] == "unavailable"
    @test isempty(rec.snapshots) && isempty(rec.epochs)
    @test rec.sections["stages"]["mixed"] == "not_applicable"
    dir = mktempdir()
    @test isfile(write_validity(joinpath(dir, "v.toml"), rec))
    @test length(readlines(write_validity_snapshots(joinpath(dir, "v.csv"), rec))) == 1
    rec0 = validity_record(Dict{Symbol,Any}(), PostprocessConfig(); calibration = nothing)
    @test rec0.sections["calibration"]["version"] == -1
    for (key, bad) in (
        (:validity_interval_de, 0.0),
        (:validity_cumulative_de, -1.0),
        (:validity_mass_tolerance, 0.0),
        (:validity_settle_tcr, -1.0),
        (:validity_pericentre_tcr, -0.5),
        (:validity_window_tcr, 0.0),
        (:validity_virial_tolerance, 0.0),
        (:validity_rate_tolerance, 0.0),
        (:validity_bias_tolerance, 0.0),
    )
        cfg_bad = PostprocessConfig(; (key => bad,)...)
        @test_throws ArgumentError Nbody6Dynamics._validate_estimators(cfg_bad)
        @test_throws Regex("postprocess\\.$key") Nbody6Dynamics._validate_estimators(cfg_bad)
    end
    @test Nbody6Dynamics._validate_estimators(PostprocessConfig(; validity_settle_tcr = 0.0)) ===
          nothing
end
