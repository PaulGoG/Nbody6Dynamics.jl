# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Figure extension" begin
    # This process loaded CairoMakie, so the extension is in and every
    # declared entry point has an implementation behind its fallback.
    @test plotting_available()
    @test Base.get_extension(Nbody6Dynamics, :Nbody6DynamicsMakieExt) !== nothing
    for entry_point in Nbody6Dynamics._PLOTTING_ENTRY_POINTS
        @test length(methods(getfield(Nbody6Dynamics, entry_point))) > 1
    end
    # With a backend loaded, wrong arguments stay a MethodError instead
    # of being reported as a missing backend.
    @test_throws MethodError plot_energy(:not_a_diagnostics_object)
    # The fallbacks and the methods the extension adds to the same
    # functions must not be ambiguous with one another.
    @test isempty(detect_ambiguities(Nbody6Dynamics, MakieExt))

    msg = sprint(showerror, Nbody6Dynamics.PlottingUnavailable(:plot_energy, "Remedy line."))
    @test occursin("using CairoMakie", msg)
    @test occursin("Remedy line.", msg)

    # A session without a backend: the package loads, never pulls Makie
    # in, reports the remedy from an entry point, and refuses a pipeline
    # that must produce figures before it builds or integrates anything.
    headless = """
    using Nbody6Dynamics
    println("available=", plotting_available())
    try
        plot_energy("a", "b")
    catch e
        println("entry=", nameof(typeof(e)), ":", e.entry_point)
    end
    cfg_path = tempname() * ".toml"
    write(cfg_path, join([
        "[install]", "enabled = false",
        "[simulation]", "run_test = false",
        "[postprocess]", "enabled = false",
        "[visualization]", "enabled = true",
    ], "\\n"))
    try
        run_pipeline(load_config(cfg_path))
    catch e
        println("pipeline=", nameof(typeof(e)), ":", e.entry_point)
        println("hint=", !isempty(e.hint))
    end
    rm(cfg_path; force = true)
    println("makie_loaded=", any(m -> nameof(m) === :Makie, Base.loaded_modules_array()))
    """
    out = readchomp(
        `$(Base.julia_cmd()) --project=$(Nbody6Dynamics._PROJECT_ROOT) --startup-file=no -e $headless`,
    )
    @test occursin("available=false", out)
    @test occursin("entry=PlottingUnavailable:plot_energy", out)
    @test occursin("pipeline=PlottingUnavailable:run_pipeline", out)
    @test occursin("hint=true", out)
    @test occursin("makie_loaded=false", out)
end

# =====================================================================

@testset "Plotting (smoke tests)" begin
    Nbody6Dynamics.set_publication_theme!()
    # The theme's faces must be live in this (fresh) process: a face
    # captured at precompile time has a null FreeType pointer and Makie
    # falls back to its sans default without warning.
    theme_fonts = publication_theme().fonts
    for key in (:regular, :bold, :italic)
        @test getfield(theme_fonts[key][], :ft_ptr) != C_NULL
    end
    @test theme_fonts[:regular][] === MakieExt.texfont(:text)

    vis = VisualizationConfig(;
        output_dir = joinpath(TESTDIR, "test_plots"),
        format = "png",
        dpi = 72,
    )

    # --- Snapshot plots ---
    n = 50
    params = zeros(Float32, 20)
    params[1] = 1.0f0
    hdr = SnapshotHeader(Int32(n), Int32(1), Int32(1), Int32(20), params)
    snap = Snapshot(
        hdr,
        Int32.(1:n),
        Float32.(rand(n)),
        Float32.(randn(3, n)),
        Float32.(0.1 .* randn(3, n)),
        Float32[],
        Float32[],
    )

    plot_snapshot(snap, vis; filename = "test_snap")
    @test isfile(joinpath(TESTDIR, "test_plots", "test_snap_xy.png"))
    @test isfile(joinpath(TESTDIR, "test_plots", "test_snap_xz.png"))

    # --- Energy plot ---
    adj = [
        AdjustRecord(0.0, 0.0, 1.0, 0.0, -0.25, 100, 10, 1.0),
        AdjustRecord(0.5, 50.0, 0.98, 1e-6, -0.249, 99, 9, 1.1),
        AdjustRecord(1.0, 100.0, 0.99, 2e-6, -0.248, 98, 8, 1.2),
    ]
    diag = DiagnosticsData(adj, Dict{String,Float64}())
    plot_energy(diag, vis; filename = "test_energy")
    @test isfile(joinpath(TESTDIR, "test_plots", "test_energy.png"))

    # --- Lagrangian plot ---
    lagr = LagrangianData(
        [0.0, 0.5, 1.0],
        [0.1, 0.5, 1.0],
        [0.05 0.06 0.07; 0.2 0.22 0.25; 1.0 1.05 1.1],
    )
    plot_lagrangian(lagr, vis; filename = "test_lagr", selected_fractions = [0.1, 0.5, 1.0])
    @test isfile(joinpath(TESTDIR, "test_plots", "test_lagr.png"))

    # --- HR diagram ---
    sev = StellarEvolutionSnapshot(
        0.5,
        4,
        [
            StellarRecord(0.5, Int32(1), Int32(101), Int32(0), 1.0, 0.8, 0.1, -0.5, 3.75),
            StellarRecord(0.5, Int32(2), Int32(102), Int32(1), 0.8, 1.2, 1.5, 0.2, 4.10),
            StellarRecord(0.5, Int32(3), Int32(103), Int32(2), 1.5, 0.7, 2.0, 0.8, 3.60),
            StellarRecord(0.5, Int32(4), Int32(104), Int32(13), 2.0, 10.0, 5.0, 1.5, 4.50),
        ],
    )
    plot_hr(sev, vis; filename = "test_hr")
    @test isfile(joinpath(TESTDIR, "test_plots", "test_hr.png"))

    # --- HR evolution ---
    sev2 = StellarEvolutionSnapshot(
        1.0,
        2,
        [
            StellarRecord(1.0, Int32(1), Int32(101), Int32(2), 1.5, 0.75, 0.5, -0.2, 3.60),
            StellarRecord(1.0, Int32(2), Int32(102), Int32(4), 0.9, 1.10, 2.0, 0.4, 3.90),
        ],
    )
    plot_hr_evolution([sev, sev2], vis; filename = "test_hr_evo")
    @test isfile(joinpath(TESTDIR, "test_plots", "test_hr_evo.png"))

    # --- HR figure of one epoch of a run ---
    plot_hr([sev, sev2], vis; epoch = 1, filename = "test_hr_epoch")
    @test isfile(joinpath(TESTDIR, "test_plots", "test_hr_epoch.png"))
    @test_throws ArgumentError plot_hr([sev, sev2], vis; epoch = 3)
    @test_throws ArgumentError plot_hr_evolution(StellarEvolutionSnapshot[], vis)
    @test_throws ArgumentError plot_hr_evolution([sev, sev2], vis; epochs = [0])

    # The neutron star (log Teff = 4.5) is counted, not drawn: it must
    # not set the axis shared by every HR figure of the run.
    hr_run = MakieExt._hr_run([sev, sev2], BinaryEvolutionSnapshot[])
    @test hr_run.tlims[2] < 4.5
    @test stellar_class_index(13) in hr_run.classes
    @test MakieExt._hr_has_strip(hr_run)
    @test !MakieExt._hr_has_strip(MakieExt._hr_run([sev], BinaryEvolutionSnapshot[]))
    @test MakieExt._hr_dark_note(hr_run, 1) == "Not on the plane: 1 neutron star"
    @test MakieExt._hr_dark_note(hr_run, 2) == ""
    @test MakieExt._census_axis_top.((1, 3, 7, 10, 11, 120)) == (1.0, 5.0, 10.0, 10.0, 20.0, 200.0)
    # One style per class, and no two drawn classes share colour and marker.
    @test Set(keys(MakieExt._HR_CLASS_STYLE)) == Set(c.key for c in STELLAR_CLASSES)
    styles = [(st.color, st.marker) for st in values(MakieExt._HR_CLASS_STYLE)]
    @test allunique(styles)
end

# =====================================================================

@testset "Escaper and SSE plots (smoke tests)" begin
    Nbody6Dynamics.set_publication_theme!()

    plots_dir = joinpath(TESTDIR, "test_plots_f23")
    vis = VisualizationConfig(; output_dir = plots_dir, format = "png", dpi = 72)

    # --- Escaper plots ---
    escs = [
        EscaperRecord(1.0, 0.5, -0.1, 15.0, 0, 101, 10.0, -20.0),
        EscaperRecord(2.0, 0.3, 0.2, 30.0, 1, 202, 120.0, 35.0),
        EscaperRecord(3.5, 1.4, 0.5, 250.0, 14, 303, 300.0, -60.0),
        EscaperRecord(4.0, 0.6, 0.1, 22.0, 11, 404, 200.0, 5.0),
    ]

    plot_escapers(escs, vis; filename = "test_escapers")
    @test isfile(joinpath(plots_dir, "test_escapers.png"))

    plot_escape_anisotropy(escs, vis; filename = "test_esc_aniso")
    @test isfile(joinpath(plots_dir, "test_esc_aniso.png"))

    # Empty-input guards: warn, return nothing, produce no file
    ret = @test_logs (:warn, r"No escaper records") plot_escapers(
        EscaperRecord[],
        vis;
        filename = "test_escapers_empty",
    )
    @test ret === nothing
    @test !isfile(joinpath(plots_dir, "test_escapers_empty.png"))

    ret = @test_logs (:warn, r"No escaper records") plot_escape_anisotropy(
        EscaperRecord[],
        vis;
        filename = "test_esc_aniso_empty",
    )
    @test ret === nothing
    @test !isfile(joinpath(plots_dir, "test_esc_aniso_empty.png"))

    # --- SSE plots ---
    # Mixed snapshot: two MS stars (finite TM, one past turnoff), one
    # giant with a partial core, one NS remnant (MC = M).  Records built
    # with the 9-argument convenience constructor carry NaN SSE fields
    # and must be skipped gracefully.
    sev_mix = StellarEvolutionSnapshot(
        50.0,
        5,
        [
            StellarRecord(
                1.0,
                Int32(1),
                Int32(1),
                Int32(0),
                0.5,
                0.4,
                -0.6,
                -0.3,
                3.65,
                8.0e4,
                0.0,
                0.0,
                1e-10,
            ),
            StellarRecord(
                1.0,
                Int32(2),
                Int32(2),
                Int32(1),
                1.2,
                1.0,
                0.0,
                0.0,
                3.76,
                30.0,
                0.0,
                0.0,
                1e-10,
            ),
            StellarRecord(
                1.0,
                Int32(3),
                Int32(3),
                Int32(3),
                2.0,
                1.8,
                1.8,
                1.2,
                3.68,
                40.0,
                0.25,
                0.02,
                15.0,
            ),
            StellarRecord(
                1.0,
                Int32(4),
                Int32(4),
                Int32(13),
                3.0,
                1.4,
                -5.0,
                -5.0,
                5.0,
                10.0,
                1.4,
                1e-5,
                1e-10,
            ),
            StellarRecord(1.0, Int32(5), Int32(5), Int32(1), 0.9, 0.8, -0.1, -0.1, 3.70),
        ],
    )

    plot_mass_segregation(sev_mix, vis; filename = "test_mass_seg")
    @test isfile(joinpath(plots_dir, "test_mass_seg.png"))

    plot_evolutionary_clock(sev_mix, vis; filename = "test_evo_clock")
    @test isfile(joinpath(plots_dir, "test_evo_clock.png"))

    plot_core_mass([sev_mix], vis; filename = "test_core_mass")
    @test isfile(joinpath(plots_dir, "test_core_mass.png"))

    # MS-only snapshot: no evolved stars — core-mass plot must warn,
    # return nothing, and create no file.
    sev_ms = StellarEvolutionSnapshot(
        0.5,
        2,
        [
            StellarRecord(
                0.5,
                Int32(1),
                Int32(1),
                Int32(0),
                1.0,
                0.8,
                0.1,
                -0.5,
                3.75,
                90.0,
                0.0,
                0.0,
                1e-10,
            ),
            StellarRecord(
                0.5,
                Int32(2),
                Int32(2),
                Int32(1),
                0.8,
                1.2,
                1.5,
                0.2,
                4.10,
                5.2,
                0.0,
                0.0,
                1e-10,
            ),
        ],
    )
    ret = @test_logs (:warn, r"No evolved stars") plot_core_mass(
        [sev_ms],
        vis;
        filename = "test_core_mass_ms",
    )
    @test ret === nothing
    @test !isfile(joinpath(plots_dir, "test_core_mass_ms.png"))

    # NaN-TM main-sequence snapshot (pre-v2026.07 data): the clock has
    # no valid TM and must warn + return nothing without a file.
    sev_nan = StellarEvolutionSnapshot(
        0.5,
        1,
        [StellarRecord(0.5, Int32(1), Int32(1), Int32(1), 1.0, 0.8, 0.1, -0.5, 3.75)],
    )
    ret = @test_logs (:warn, r"No main-sequence records") plot_evolutionary_clock(
        sev_nan,
        vis;
        filename = "test_evo_clock_nan",
    )
    @test ret === nothing
    @test !isfile(joinpath(plots_dir, "test_evo_clock_nan.png"))
end

# =====================================================================

@testset "Auto FPS calculation" begin
    _auto_fps = MakieExt._auto_fps

    # Few frames → clamped to min
    @test _auto_fps(5; target_duration = 12.0, min_fps = 1, max_fps = 10) == 1
    # Many frames → clamped to max
    @test _auto_fps(500; target_duration = 12.0, min_fps = 2, max_fps = 30) == 30
    # Medium range → computed value
    fps = _auto_fps(120; target_duration = 12.0, min_fps = 2, max_fps = 30)
    @test fps == 10
    # Exact target
    @test _auto_fps(24; target_duration = 12.0, min_fps = 1, max_fps = 30) == 2
end

# =====================================================================

@testset "Animations (smoke tests)" begin
    Nbody6Dynamics.set_publication_theme!()

    vis = VisualizationConfig(;
        output_dir = joinpath(TESTDIR, "test_anims"),
        format = "png",
        dpi = 72,
    )

    # Build two simple test snapshots
    n = 30
    params1 = zeros(Float32, 20)
    params1[1] = 0.0f0
    params2 = zeros(Float32, 20)
    params2[1] = 1.0f0
    hdr1 = SnapshotHeader(Int32(n), Int32(1), Int32(1), Int32(20), params1)
    hdr2 = SnapshotHeader(Int32(n), Int32(2), Int32(1), Int32(20), params2)
    snap1 = Snapshot(
        hdr1,
        Int32.(1:n),
        Float32.(rand(n)),
        Float32.(randn(3, n)),
        Float32.(0.1 .* randn(3, n)),
        Float32[],
        Float32[],
    )
    snap2 = Snapshot(
        hdr2,
        Int32.(1:n),
        Float32.(rand(n)),
        Float32.(randn(3, n) .+ 0.5),
        Float32.(0.1 .* randn(3, n)),
        Float32[],
        Float32[],
    )

    # --- Cluster animation (explicit fps) ---
    outpaths = animate_cluster([snap1, snap2], vis; filename = "test_cluster_anim", fps = 2)
    @test all(isfile, outpaths)
    @test all(p -> endswith(p, ".gif"), outpaths)

    # --- Cluster animation (auto fps) ---
    outpaths = animate_cluster([snap1, snap2], vis; filename = "test_cluster_anim_auto")
    @test all(isfile, outpaths)
    @test all(p -> endswith(p, ".gif"), outpaths)

    # --- Lagrangian animation ---
    lagr = LagrangianData(
        [0.0, 0.5, 1.0],
        [0.1, 0.5, 1.0],
        [0.05 0.06 0.07; 0.2 0.22 0.25; 1.0 1.05 1.1],
    )
    outpath = animate_lagrangian(
        lagr,
        vis;
        filename = "test_lagr_anim",
        fps = 2,
        selected_fractions = [0.1, 0.5, 1.0],
    )
    @test isfile(outpath)
    @test endswith(outpath, ".gif")

    # --- Lagrangian animation (auto fps) ---
    outpath = animate_lagrangian(
        lagr,
        vis;
        filename = "test_lagr_anim_auto",
        selected_fractions = [0.1, 0.5, 1.0],
    )
    @test isfile(outpath)
    @test endswith(outpath, ".gif")

    # --- HR animation ---
    sev1 = StellarEvolutionSnapshot(
        0.5,
        2,
        [
            StellarRecord(0.5, Int32(1), Int32(101), Int32(0), 1.0, 0.8, 0.1, -0.5, 3.75),
            StellarRecord(0.5, Int32(2), Int32(102), Int32(1), 0.8, 1.2, 1.5, 0.2, 4.10),
        ],
    )
    sev2 = StellarEvolutionSnapshot(
        1.0,
        2,
        [
            StellarRecord(1.0, Int32(1), Int32(101), Int32(2), 1.5, 0.75, 0.5, -0.2, 3.60),
            StellarRecord(1.0, Int32(2), Int32(102), Int32(4), 0.9, 1.10, 2.0, 0.4, 3.90),
        ],
    )
    outpath = animate_hr([sev1, sev2], vis; filename = "test_hr_anim", fps = 2)
    @test isfile(outpath)
    @test endswith(outpath, ".gif")

    # --- HR animation (auto fps) ---
    outpath = animate_hr([sev1, sev2], vis; filename = "test_hr_anim_auto")
    @test isfile(outpath)
    @test endswith(outpath, ".gif")
end

# =====================================================================

@testset "Merger plot suite smoke" begin
    # plot_merger_ic on a small generated IC (all five figure families)
    ic_dir = mktempdir()
    rng_plot = StableRNG(4242)
    cfg_ic = MergerConfig(
        [
            ClusterSpec(
                profile = PlummerProfile(),
                N = 60,
                rbar = 1.0,
                imf = RescaledKroupaIMF(
                    bodyn = 0.1,
                    body1 = 50.0,
                    target_mass = 60 * kroupa_mean_mass(0.1, 50.0),
                ),
            ),
            ClusterSpec(
                profile = KingProfile(W0 = 5.0),
                N = 60,
                rbar = 1.0,
                imf = EqualMassIMF(particle_mass = 5.0),
            ),
        ],
        "kepler",
        OrbitSpec(apocentre = 10.0, eccentricity = 0.5),
        MergerOutputSpec(format = "nbody", truncate_jacobi = false, output_dir = ic_dir),
    )
    ic_result = generate_merger_ic(cfg_ic; rng = rng_plot)
    plots_dir = mktempdir()
    vis_smoke = VisualizationConfig(; format = "png", output_dir = plots_dir)
    plot_merger_ic(ic_result, vis_smoke)
    for stem in (
        "merger_ic_xy",
        "merger_ic_xz",
        "merger_ic_yz",
        "merger_ic_overview",
        "merger_ic_velocity",
        "merger_ic_imf",
        "merger_ic_density",
    )
        @test isfile(joinpath(plots_dir, stem * ".png"))
    end

    # Synthetic two-cluster snapshots: converging COMs, NB units
    function _smoke_snapshot(t, d)
        n_half = 40
        params = zeros(Float32, 20)
        params[1] = Float32(t)
        pos = zeros(Float32, 3, 2 * n_half)
        vel = Float32.(0.05 .* randn(rng_plot, 3, 2 * n_half))
        pos[:, 1:n_half] .= Float32.(0.3 .* randn(rng_plot, 3, n_half))
        pos[1, 1:n_half] .-= Float32(d)
        pos[:, (n_half + 1):end] .= Float32.(0.3 .* randn(rng_plot, 3, n_half))
        pos[1, (n_half + 1):end] .+= Float32(d)
        Snapshot(
            SnapshotHeader(Int32(2 * n_half), Int32(1), Int32(1), Int32(20), params),
            Int32.(1:(2 * n_half)),
            fill(Float32(1 / (2 * n_half)), 2 * n_half),
            pos,
            vel,
            Float32[],
            Float32[],
        )
    end
    snaps = [_smoke_snapshot(t, d) for (t, d) in ((0.0, 4.0), (1.0, 2.0), (2.0, 0.5))]
    ranges = [1:40, 41:80]

    Q, n_mem = per_cluster_virial(snaps, ranges)
    @test size(Q) == (2, 3)
    @test all(n_mem .== 40)
    @test all(isfinite, Q)

    plot_cluster_separation(snaps, ranges, vis_smoke)
    @test isfile(joinpath(plots_dir, "merger_cluster_separation.png"))
    plot_cluster_virial(snaps, ranges, vis_smoke)
    @test isfile(joinpath(plots_dir, "merger_cluster_virial.png"))
    # Membership as body-index blocks (what parse_merger_summary returns)
    Q_blocks, n_blocks = per_cluster_virial(snaps, collect.(ranges))
    @test Q_blocks == Q && n_blocks == n_mem
    plot_cluster_virial(
        snaps,
        collect.(ranges),
        vis_smoke;
        filename = "merger_cluster_virial_blocks",
    )
    @test isfile(joinpath(plots_dir, "merger_cluster_virial_blocks.png"))

    # Envelope statistics helper: NaNs are skipped, all-NaN columns stay NaN
    env = [1.0 NaN 3.0; 5.0 NaN 1.0]
    lo, hi, mean_vals = MakieExt._envelope_stats(env)
    @test lo == [1.0, NaN, 1.0] || (lo[1] == 1.0 && isnan(lo[2]) && lo[3] == 1.0)
    @test hi[1] == 5.0 && isnan(hi[2]) && hi[3] == 3.0
    @test mean_vals[1] == 3.0 && isnan(mean_vals[2]) && mean_vals[3] == 2.0
end

# =====================================================================

@testset "Canvases in layout units, export at the printed width" begin
    # The canvas is fixed by the figure type; the configuration only
    # sets the width it is printed at.
    for col in ("", "single", "double")
        cfg = VisualizationConfig(; column = col)
        @test MakieExt._figsize_px(cfg) == (900, 600)
        @test MakieExt._fig_two_panel(cfg) == (900, 950)
        @test MakieExt._fig_with_colorbar(cfg) == (900 + MakieExt._COLORBAR_WIDTH, 600)
    end
    cfg = VisualizationConfig(; export_width = 5.5)
    @test MakieExt._export_width_in(cfg) == 5.5
    @test MakieExt._export_width_in(VisualizationConfig(; column = "single")) == 3.4
    @test MakieExt._export_width_in(VisualizationConfig(; column = "double")) == 7.05

    # Stacks: a single panel plus 350 per further row, plus what is reserved
    @test MakieExt._fig_multipanel(cfg, 2, 1) == (900, 950)
    @test MakieExt._fig_multipanel(cfg, 3, 1; extra_height = 180) == (900, 1480)
    # Grids: at least 1200 wide, 500 per column beyond that
    @test MakieExt._grid_canvas_width.((1, 2, 3, 4)) == (900, 1200, 1500, 2000)
    gap = MakieExt._MULTIPANEL_GAP_COMPACT
    strip = MakieExt._AXIS_PROTRUSION
    box_w = (1500 - strip - 2 * gap) / 3
    @test box_w == MakieExt._multipanel_box_width(cfg, 3; inner_ticks = false)
    w3, h3 = MakieExt._fig_multipanel(cfg, 2, 3; inner_ticks = false)
    @test w3 == 1500
    @test h3 == round(Int, 2 * box_w * 600 / 900 + gap + strip)
    @test MakieExt._fig_multipanel(cfg, 2, 3; inner_ticks = false, extra_height = 30)[2] ==
          round(Int, 2 * box_w * 600 / 900 + gap + strip + 30)
    # Square panels are taller; a reserved colorbar column narrows them
    square = MakieExt._fig_multipanel(cfg, 2, 3; inner_ticks = false, panel_aspect = 1.0)[2]
    @test square > h3
    @test MakieExt._fig_multipanel(
        cfg,
        2,
        3;
        inner_ticks = false,
        panel_aspect = 1.0,
        extra_width = MakieExt._COLORBAR_WIDTH,
    )[2] < square
    @test MakieExt._multipanel_gap(3; inner_ticks = false) == MakieExt._MULTIPANEL_GAP_COMPACT
    @test MakieExt._multipanel_gap(3; inner_ticks = true) == MakieExt._MULTIPANEL_HGAP
    @test MakieExt._multipanel_gap(1; inner_ticks = false) == MakieExt._MULTIPANEL_HGAP
    # Marker scale follows the panel width; the annotation band is a data-free strip
    @test MakieExt._multipanel_scale(cfg, 1) == 1.0
    s3 = MakieExt._multipanel_scale(cfg, 3; inner_ticks = false)
    @test s3 ≈ box_w / 900
    @test MakieExt._multipanel_scale(cfg, 3; inner_ticks = false, extra_width = 140) < s3
    @test MakieExt._multipanel_scale(cfg, 3; inner_ticks = true) < s3
    @test 0 < MakieExt._MONTAGE_BAND_FRAC < 0.5

    # The canvas width becomes the printed width: 900 units over 6.5 in
    # is 0.52 pt per unit (468 pt), over a 3.4 in column 0.272; a grid
    # 1500 wide over the same 6.5 in is 0.312. Raster output is never
    # coarser than 4 pixels per unit and follows dpi above that.
    std = VisualizationConfig(; dpi = 300)
    @test MakieExt._export_scale(std, 900).pt_per_unit ≈ 0.52
    @test MakieExt._export_scale(std, 1500).pt_per_unit ≈ 0.312
    @test MakieExt._export_scale(VisualizationConfig(; column = "single"), 900).pt_per_unit ≈ 0.272
    @test MakieExt._export_scale(std, 900).px_per_unit == 4.0
    @test MakieExt._export_scale(VisualizationConfig(; dpi = 1200), 900).px_per_unit ≈
          1200 * 6.5 / 900
    mktempdir() do dir
        fig = CairoMakie.Figure(; size = (900, 600))
        CairoMakie.Axis(fig[1, 1])
        png = MakieExt._save_fig(
            VisualizationConfig(; format = "png", dpi = 72, output_dir = dir),
            "probe",
            fig,
        )
        header = read(png)[17:24]
        @test reinterpret(UInt32, reverse(header[1:4]))[1] == 3600
        @test reinterpret(UInt32, reverse(header[5:8]))[1] == 2400
        pdf = MakieExt._save_fig(
            VisualizationConfig(; format = "pdf", output_dir = dir),
            "probe",
            fig,
        )
        @test filesize(pdf) > 0
    end
end

@testset "Figure routines scope the theme" begin
    CairoMakie.set_theme!()
    ambient = CairoMakie.Makie.to_value(CairoMakie.Makie.theme(nothing, :fontsize))
    MakieExt.@publication function _theme_probe(x)
        x < 0 && return nothing
        return CairoMakie.Makie.to_value(CairoMakie.Makie.theme(nothing, :fontsize))
    end
    @test _theme_probe(1) == MakieExt._STYLE.label == 26
    @test _theme_probe(-1) === nothing
    @test CairoMakie.Makie.to_value(CairoMakie.Makie.theme(nothing, :fontsize)) == ambient
    theme = publication_theme()
    @test theme.linewidth[] == 3
    @test theme.markersize[] == 14
    # 10 all round, 30 on the right for the overhang of the last x tick label
    @test theme.figure_padding[] == (10, 30, 10, 10)
    @test theme.Axis.xticklabelsize[] == 22
    @test theme.Legend.framevisible[] == false
end

# =====================================================================
# Edge cases for pure helpers and degenerate reader inputs
# =====================================================================

@testset "Degenerate axis ranges" begin
    # Identical stars (equal-mass, unevolved) give zero-span HR data
    @test MakieExt._padded_range(3.678, 3.678) == (3.578, 3.778)
    lo, hi = MakieExt._padded_range(3.0, 4.0)
    @test lo ≈ 2.94 && hi ≈ 4.06
    @test !isempty(MakieExt._logval_ticks(3.62, 3.74))
    @test MakieExt._logval_ticks(3.678, 3.678) == [3.678]
    @test !isempty(MakieExt._nice_ticks(1.0, 1.0001))
    @test MakieExt._nice_ticks(0.0, 10.0) == collect(0.0:1.0:10.0) ||
          !isempty(MakieExt._nice_ticks(0.0, 10.0))
end

@testset "Number formatting never uses computer notation" begin
    f = MakieExt._fmt_latex_sig
    @test f(0) == "0"
    @test f(1200.0, 3) == "1200"       # "%.3g" gives 1.2e+03
    @test f(100.0, 2) == "100"         # "%.2g" gives 1e+02
    @test f(150.0, 2) == "150"
    @test f(99.96, 3) == "100"
    @test f(12.34, 2) == "12"
    @test f(0.5, 2) == "0.5"
    @test f(0.0123, 3) == "0.0123"
    @test f(12345.0, 3) == "1.23 \\times 10^{4}"
    @test f(2.5e-4, 2) == "2.5 \\times 10^{-4}"
    @test f(1.0e5, 3) == "10^{5}"      # 1×10ⁿ collapses to 10ⁿ
    @test MakieExt._fmt_latex_sig3(47.25) == "47.2" || MakieExt._fmt_latex_sig3(47.25) == "47.3"
    samples = (3.0e-7, 0.004, 0.07, 1.0, 9.99, 100.0, 999.5, 1.0e4, 6.02e23)
    @test !any(occursin(r"[0-9]e[+-]?[0-9]", f(x, n)) for x in samples, n in (2, 3))
end

@testset "Log-tick generator edge cases" begin
    # >2 in-range decades → decades only
    vals, _ = MakieExt._log_ticks(0.05, 50.0)
    @test vals == [0.1, 1.0, 10.0]
    # ≤2 in-range decades → 2×/5× intermediates appear
    vals2, labels2 = MakieExt._log_ticks(0.5, 30.0)
    @test all(v -> v in vals2, (0.5, 1.0, 2.0, 5.0, 10.0, 20.0))
    # 10^0 renders as plain "1"
    lab1 = String(labels2[findfirst(==(1.0), vals2)])
    @test occursin("1", lab1) && !occursin("10", lab1)
    # Degenerate equal endpoints still yield ≥ 2 ticks
    vals3, _ = MakieExt._log_ticks(2.0, 2.0)
    @test length(vals3) ≥ 2
    # Plain decimals throughout on a short axis within 10⁻³–10⁴ (the
    # virial-ratio panel: 0.5 … 10, never "5 × 10⁻¹ … 10¹")
    vals4, labels4 = MakieExt._log_ticks(0.45, 13.0)
    @test vals4 == [0.5, 1.0, 2.0, 5.0, 10.0]
    @test [MakieExt._log_tick_label(v, true) for v in vals4] == ["0.5", "1", "2", "5", "10"]
    @test all(l -> !occursin("times", String(l)) && !occursin("^", String(l)), labels4)
    @test [MakieExt._log_tick_label(v, true) for v in (0.001, 0.01, 0.1, 100.0, 50000.0)] == ["0.001", "0.01", "0.1", "100", "50000"]
    # Exponent form beyond the plain window, with the mandatory collapses
    _, labels5 = MakieExt._log_ticks(1e-8, 1e-3)
    @test occursin("10^{-8}", String(labels5[1])) && occursin("10^{-3}", String(labels5[end]))
    _, labels6 = MakieExt._log_ticks(0.5, 1e6)
    @test occursin("1", String(labels6[1])) && !occursin("10", String(labels6[1]))
    @test String(labels6[2]) == "\$10\$" && occursin("10^{2}", String(labels6[3]))
    @test MakieExt._log_tick_label(2e-7, false) == "2\\times 10^{-7}"
    @test MakieExt._log_tick_label(0.2, false) == "0.2"
    @test MakieExt._log_tick_label(20.0, false) == "20"
    @test MakieExt._log_tick_label(200.0, false) == "2\\times 10^{2}"
end
