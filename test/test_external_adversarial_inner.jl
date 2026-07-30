# =============================================================================
# Adversarial tests for external post-processing (included from runtests.jl)
# =============================================================================
# Tests scan_output, postprocess_external, and run_pipeline with incomplete,
# corrupt, and edge-case Nbody6++ output data.

const TDIR = mktempdir()

# Helper: build a minimal conf.3 snapshot file (N particles, time t)
# Matches the Nbody6++ bulk-array format:
#   Record 1: NTOT, MODEL, NRUN, NK  (4 × Int32)
#   Record 2: AS(1:NK) ++ BODYS(N) ++ RHOS(N) ++ XNS(N) ++ XS(3,N) ++ VS(3,N) ++ PHI(N) ++ NAME(N)
#   All Float32 except NAME which is Int32.
function write_fake_conf3(path::String, N::Int, t::Float64)
    NK = Int32(20)
    open(path, "w") do io
        # Record 1: 4 × Int32
        _write_fortran_record(io, Int32[Int32(N), Int32(1), Int32(1), NK])

        # Record 2: bulk array
        params = zeros(Float32, NK)
        params[1] = Float32(t)
        params[2] = Float32(N)

        mass = ones(Float32, N) ./ N
        rho  = zeros(Float32, N)
        xns  = zeros(Float32, N)
        pos  = randn(Float32, 3 * N) .* 5
        vel  = randn(Float32, 3 * N) .* 0.3f0
        phi  = zeros(Float32, N)
        name = reinterpret(Float32, collect(Int32(1):Int32(N)))

        all_data = vcat(params, mass, rho, xns, pos, vel, phi, name)
        _write_fortran_record(io, all_data)
    end
end

@testset "Adversarial External Post-Processing" begin

    # =====================================================================
    @testset "scan_output edge cases" begin

        # --- Non-existent directory ---
        @test_throws ErrorException scan_output("/this/does/not/exist/at/all")

        # --- Empty directory ---
        d_empty = joinpath(TDIR, "empty")
        mkpath(d_empty)
        s = scan_output(d_empty)
        @test all(v -> v == false, values(s.available))
        @test isempty(s.conf3_files)
        @test isempty(s.hdf5_files)
        @test isempty(s.stellar_evo_files)

        # --- Only a stdout file, nothing else ---
        d_stdout_only = joinpath(TDIR, "stdout_only")
        mkpath(d_stdout_only)
        open(joinpath(d_stdout_only, "out1000"), "w") do io
            println(io, "some random text that is not ADJUST")
        end
        s = scan_output(d_stdout_only)
        @test s.available[:diagnostics] == true
        @test s.available[:snapshots_conf3] == false
        @test s.available[:lagr] == false
        @test s.available[:stellar_evo] == false

        # --- Decoy files: wrong naming convention ---
        d_decoy = joinpath(TDIR, "decoy")
        mkpath(d_decoy)
        touch(joinpath(d_decoy, "conf3_0"))       # wrong: should be conf.3_0
        touch(joinpath(d_decoy, "CONF.3_0"))       # wrong case
        touch(joinpath(d_decoy, "sev_83_0"))        # wrong: should be sev.83_0
        touch(joinpath(d_decoy, "lagr7"))           # wrong: should be lagr.7
        touch(joinpath(d_decoy, "esc11"))           # wrong: should be esc.11
        touch(joinpath(d_decoy, "output.txt"))      # not a known stdout name
        s = scan_output(d_decoy)
        @test all(v -> v == false, values(s.available))

        # --- Mixed valid + invalid ---
        d_mixed = joinpath(TDIR, "mixed")
        mkpath(d_mixed)
        write_fake_conf3(joinpath(d_mixed, "conf.3_0"), 100, 0.0)
        write_fake_conf3(joinpath(d_mixed, "conf.3_1"), 100, 1.0)
        touch(joinpath(d_mixed, "conf3_bogus"))     # decoy
        touch(joinpath(d_mixed, "lagr.7"))           # empty but exists
        s = scan_output(d_mixed)
        @test s.available[:snapshots_conf3] == true
        @test length(s.conf3_files) == 2
        @test s.available[:lagr] == true
        @test s.available[:diagnostics] == false

        # --- Bare conf.3 (no suffix) ---
        d_bare = joinpath(TDIR, "bare_conf3")
        mkpath(d_bare)
        write_fake_conf3(joinpath(d_bare, "conf.3"), 50, 0.0)
        s = scan_output(d_bare)
        @test s.available[:snapshots_conf3] == true
        @test length(s.conf3_files) == 1

        # --- HDF5 detection ---
        d_hdf5 = joinpath(TDIR, "hdf5_decoy")
        mkpath(d_hdf5)
        touch(joinpath(d_hdf5, "data.h5part"))
        touch(joinpath(d_hdf5, "backup.hdf5"))
        touch(joinpath(d_hdf5, "notahdf5.h5"))     # .h5 not matched
        s = scan_output(d_hdf5)
        @test s.available[:snapshots_hdf5] == true
        @test length(s.hdf5_files) == 2  # .h5part + .hdf5

        # --- Alternative stdout names ---
        for name in ["out1", "stdout", "output.log"]
            d_alt = joinpath(TDIR, "alt_stdout_$name")
            mkpath(d_alt)
            touch(joinpath(d_alt, name))
            s = scan_output(d_alt)
            @test s.available[:diagnostics] == true
            @test basename(s.stdout_file) == name
        end

        # --- Priority: out1000 preferred over out1 ---
        d_prio = joinpath(TDIR, "stdout_priority")
        mkpath(d_prio)
        touch(joinpath(d_prio, "out1000"))
        touch(joinpath(d_prio, "out1"))
        s = scan_output(d_prio)
        @test basename(s.stdout_file) == "out1000"
    end

    # =====================================================================
    @testset "postprocess_external with corrupt data" begin

        # --- Truncated conf.3 file ---
        d_trunc = joinpath(TDIR, "truncated_conf3")
        mkpath(d_trunc)
        write_fake_conf3(joinpath(d_trunc, "conf.3_0"), 100, 0.0)
        open(joinpath(d_trunc, "conf.3_1"), "w") do io
            write(io, rand(UInt8, 20))
        end
        results = postprocess_external(d_trunc; generate_plots=false)
        @test haskey(results, :snapshots) || !haskey(results, :snapshots)

        # --- Empty stdout file ---
        d_empty_stdout = joinpath(TDIR, "empty_stdout")
        mkpath(d_empty_stdout)
        touch(joinpath(d_empty_stdout, "out1000"))
        results = postprocess_external(d_empty_stdout; generate_plots=false)
        @test haskey(results, :diagnostics)
        diag = results[:diagnostics]
        @test isempty(diag.adjust)

        # --- Garbage in stdout ---
        d_garbage_stdout = joinpath(TDIR, "garbage_stdout")
        mkpath(d_garbage_stdout)
        open(joinpath(d_garbage_stdout, "out1000"), "w") do io
            println(io, "This is total nonsense")
            println(io, "ADJUST: not actually valid format gibberish")
            println(io, "TIME[NB] also broken N abc NPAIRS xyz")
            println(io, "RSCALE  not parseable either")
            println(io, "12345 random numbers 678.9")
        end
        results = postprocess_external(d_garbage_stdout; generate_plots=false)
        @test haskey(results, :diagnostics)

        # --- Valid ADJUST but no TIME[NB] lines (N=0 warning) ---
        d_no_time = joinpath(TDIR, "no_timenb")
        mkpath(d_no_time)
        open(joinpath(d_no_time, "out1000"), "w") do io
            for t in 0:5
                println(io, " ADJUST:    $(Float64(t))   0.0   0.50  1.0E-08  -0.250   0   0   0.0")
            end
        end
        results = postprocess_external(d_no_time; generate_plots=false)
        diag = results[:diagnostics]
        @test length(diag.adjust) == 6
        @test all(r -> r.n == 0, diag.adjust)

        # --- Empty lagr.7 ---
        d_empty_lagr = joinpath(TDIR, "empty_lagr")
        mkpath(d_empty_lagr)
        touch(joinpath(d_empty_lagr, "lagr.7"))
        results = postprocess_external(d_empty_lagr; generate_plots=false)
        if haskey(results, :lagr)
            @test isempty(results[:lagr].time)
        end

        # --- lagr.7 with only a header (no data rows) ---
        d_hdr_lagr = joinpath(TDIR, "header_only_lagr")
        mkpath(d_hdr_lagr)
        open(joinpath(d_hdr_lagr, "lagr.7"), "w") do io
            println(io, "# TIME  0.01  0.10  0.50  1.00")
        end
        results = postprocess_external(d_hdr_lagr; generate_plots=false)
        if haskey(results, :lagr)
            @test isempty(results[:lagr].time)
        end

        # --- Empty esc.11 ---
        d_empty_esc = joinpath(TDIR, "empty_esc")
        mkpath(d_empty_esc)
        touch(joinpath(d_empty_esc, "esc.11"))
        results = postprocess_external(d_empty_esc; generate_plots=false)

        # --- esc.11 with only comments ---
        d_comment_esc = joinpath(TDIR, "comment_esc")
        mkpath(d_comment_esc)
        open(joinpath(d_comment_esc, "esc.11"), "w") do io
            println(io, "# This is a comment")
            println(io, "# Another comment")
            println(io, "")
        end
        results = postprocess_external(d_comment_esc; generate_plots=false)

        # --- Stellar evolution files that are empty ---
        d_empty_sev = joinpath(TDIR, "empty_sev")
        mkpath(d_empty_sev)
        touch(joinpath(d_empty_sev, "sev.83_0"))
        touch(joinpath(d_empty_sev, "sev.83_1"))
        results = postprocess_external(d_empty_sev; generate_plots=false)

        # --- Stellar evolution with garbage content ---
        d_garbage_sev = joinpath(TDIR, "garbage_sev")
        mkpath(d_garbage_sev)
        open(joinpath(d_garbage_sev, "sev.83_0"), "w") do io
            println(io, "not a valid stellar evolution file at all!")
            println(io, "banana apple")
        end
        results = postprocess_external(d_garbage_sev; generate_plots=false)
    end

    # =====================================================================
    @testset "Sanity checks: snapshots" begin

        # --- Non-monotonic time ---
        d_nonmono = joinpath(TDIR, "nonmono_snaps")
        mkpath(d_nonmono)
        write_fake_conf3(joinpath(d_nonmono, "conf.3_0"), 100, 0.0)
        write_fake_conf3(joinpath(d_nonmono, "conf.3_1"), 100, 5.0)
        write_fake_conf3(joinpath(d_nonmono, "conf.3_2"), 100, 3.0)  # out of order!
        write_fake_conf3(joinpath(d_nonmono, "conf.3_3"), 100, 10.0)
        results = postprocess_external(d_nonmono; generate_plots=false)
        @test haskey(results, :snapshots)
        @test length(results[:snapshots]) == 4

        # --- Duplicate snapshot times ---
        d_dup = joinpath(TDIR, "dup_snaps")
        mkpath(d_dup)
        write_fake_conf3(joinpath(d_dup, "conf.3_0"), 100, 0.0)
        write_fake_conf3(joinpath(d_dup, "conf.3_1"), 100, 5.0)
        write_fake_conf3(joinpath(d_dup, "conf.3_2"), 100, 5.0)  # duplicate!
        results = postprocess_external(d_dup; generate_plots=false)
        @test haskey(results, :snapshots)
        @test length(results[:snapshots]) == 3

        # --- Massive particle loss (>50%) ---
        d_loss = joinpath(TDIR, "particle_loss")
        mkpath(d_loss)
        write_fake_conf3(joinpath(d_loss, "conf.3_0"), 1000, 0.0)
        write_fake_conf3(joinpath(d_loss, "conf.3_1"), 400, 50.0)  # 60% loss!
        results = postprocess_external(d_loss; generate_plots=false)
        @test haskey(results, :snapshots)

        # --- Single snapshot (no evolution possible) ---
        d_single = joinpath(TDIR, "single_snap")
        mkpath(d_single)
        write_fake_conf3(joinpath(d_single, "conf.3_0"), 500, 0.0)
        results = postprocess_external(d_single; generate_plots=false)
        @test haskey(results, :snapshots)
        @test length(results[:snapshots]) == 1
    end

    # =====================================================================
    @testset "Sanity checks: diagnostics" begin

        # --- Huge energy error ---
        d_huge_de = joinpath(TDIR, "huge_de")
        mkpath(d_huge_de)
        open(joinpath(d_huge_de, "out1000"), "w") do io
            println(io, " ADJUST:    1.0   0.0   0.50  0.5  -0.250   1000   0   0.0")
            println(io, " ADJUST:    2.0   0.0   0.50  0.5  -0.250   1000   0   0.0")
        end
        results = postprocess_external(d_huge_de; generate_plots=false)
        diag = results[:diagnostics]
        @test length(diag.adjust) == 2

        # --- Wild virial ratio ---
        d_wild_q = joinpath(TDIR, "wild_qvir")
        mkpath(d_wild_q)
        open(joinpath(d_wild_q, "out1000"), "w") do io
            println(io, " ADJUST:    1.0   0.0   50.0  1e-8  -0.250   1000   0   0.0")
            println(io, " ADJUST:    2.0   0.0   100.0 1e-8  -0.250   1000   0   0.0")
        end
        results = postprocess_external(d_wild_q; generate_plots=false)
    end

    # =====================================================================
    @testset "run_pipeline edge cases" begin

        # --- Non-existent data_dir ---
        pp = PostprocessConfig(; enabled=true, data_dir="/this/does/not/exist")
        sim = SimulationConfig(; run_test=false)
        cfg = Nbody6Config(InstallConfig(; enabled=false), BuildConfig(),
                           sim, pp, VisualizationConfig(; enabled=false), MergerPipelineConfig())
        @test_throws ErrorException run_pipeline(cfg)

        # --- data_dir with only diagnostics (no snapshots) ---
        d_diag_only = joinpath(TDIR, "diag_only_pipeline")
        mkpath(d_diag_only)
        open(joinpath(d_diag_only, "out1000"), "w") do io
            println(io, " ADJUST:    1.0   0.5   0.50  1e-8  -0.250   500   10   1.2")
            println(io, " ADJUST:    2.0   1.0   0.51  2e-8  -0.249   498   10   1.3")
        end
        pp2 = PostprocessConfig(; enabled=true, data_dir=d_diag_only)
        sim2 = SimulationConfig(; run_test=false)
        vis2 = VisualizationConfig(; enabled=true, output_dir=joinpath(d_diag_only, "plots"))
        cfg2 = Nbody6Config(InstallConfig(; enabled=false), BuildConfig(),
                            sim2, pp2, vis2, MergerPipelineConfig())
        results = run_pipeline(cfg2)
        @test haskey(results, :diagnostics)
        @test !haskey(results, :snapshots)
        plots_dir = joinpath(d_diag_only, "plots")
        @test isfile(joinpath(plots_dir, "energy.png"))
        @test isfile(joinpath(plots_dir, "particle_count.png"))
        @test !isfile(joinpath(plots_dir, "snapshot_final_xy.png"))
        @test !isfile(joinpath(plots_dir, "lagrangian_radii.png"))

        # --- run_test=false, no data_dir, no runs/ directory ---
        pp3 = PostprocessConfig(; enabled=true, data_dir="")
        sim3 = SimulationConfig(; run_test=false, runs_dir=joinpath(TDIR, "nonexistent_runs"))
        cfg3 = Nbody6Config(InstallConfig(; enabled=false), BuildConfig(),
                            sim3, pp3, VisualizationConfig(; enabled=false), MergerPipelineConfig())
        results = run_pipeline(cfg3)
        @test isempty(results)

        # --- Full synthetic pipeline: snapshots + diagnostics + lagr ---
        d_full = joinpath(TDIR, "full_synthetic")
        mkpath(d_full)
        write_fake_conf3(joinpath(d_full, "conf.3_0"), 200, 0.0)
        write_fake_conf3(joinpath(d_full, "conf.3_1"), 195, 10.0)
        write_fake_conf3(joinpath(d_full, "conf.3_2"), 190, 20.0)
        open(joinpath(d_full, "out1000"), "w") do io
            for t in 0:2
                println(io, " ADJUST:    $(Float64(t*10))   0.0   0.50  1e-8  -0.250   $(200-t*5)   0   1.0")
            end
        end
        open(joinpath(d_full, "lagr.7"), "w") do io
            println(io, "# TIME  0.01  0.10  0.50  1.00")
            println(io, "0.0   0.05  0.3  1.0  5.0")
            println(io, "10.0  0.05  0.3  1.1  5.5")
            println(io, "20.0  0.04  0.3  1.2  6.0")
        end
        pp4 = PostprocessConfig(; enabled=true, data_dir=d_full)
        sim4 = SimulationConfig(; run_test=false)
        vis4 = VisualizationConfig(; enabled=true,
                                     output_dir=joinpath(d_full, "plots"))
        cfg4 = Nbody6Config(InstallConfig(; enabled=false), BuildConfig(),
                            sim4, pp4, vis4, MergerPipelineConfig())
        results = run_pipeline(cfg4)
        @test haskey(results, :snapshots)
        @test haskey(results, :diagnostics)
        @test haskey(results, :lagr)
        @test length(results[:snapshots]) == 3
        pdir = joinpath(d_full, "plots")
        @test isfile(joinpath(pdir, "snapshot_final_xy.png"))
        @test isfile(joinpath(pdir, "snapshot_final_xz.png"))
        @test isfile(joinpath(pdir, "snapshot_evolution_xy.png"))
        @test isfile(joinpath(pdir, "energy.png"))
        @test isfile(joinpath(pdir, "particle_count.png"))
        @test isfile(joinpath(pdir, "lagrangian_radii.png"))
        @test isfile(joinpath(pdir, "cluster_evolution_xy.gif"))
        @test isfile(joinpath(pdir, "lagrangian_anim.gif"))
    end

    # =====================================================================
    @testset "OutputScan display" begin
        d_display = joinpath(TDIR, "display_test")
        mkpath(d_display)
        write_fake_conf3(joinpath(d_display, "conf.3_0"), 50, 0.0)
        touch(joinpath(d_display, "lagr.7"))
        open(joinpath(d_display, "out1000"), "w") do io
            println(io, " ADJUST:    1.0   0.5   0.50  1e-8  -0.25   50   0   1.0")
        end

        scan = scan_output(d_display)
        buf = IOBuffer()
        show(buf, MIME("text/plain"), scan)
        output = String(take!(buf))

        @test occursin("conf.3 snapshots", output)
        @test occursin("HDF5 snapshots", output)
        @test occursin("Diagnostics", output)
        @test occursin("Lagrangian radii", output)
        @test occursin("Escapers", output)
        @test occursin("Stellar evolution", output)
        @test occursin("Available:", output)
        @test occursin("Plots available:", output)
    end

end  # Adversarial testset
