using Test
using Nbody6Dynamics

# Temporary directory for test artifacts
const TESTDIR = mktempdir()

# ---------------------------------------------------------------------------
# Helper: write a Fortran binary record (for test data generation)
# ---------------------------------------------------------------------------
function _write_fortran_record(io::IO, data::Vector{T}) where {T}
    bytes = reinterpret(UInt8, data)
    marker = Int32(length(bytes))
    write(io, marker)
    write(io, bytes)
    write(io, marker)
end

function _write_fortran_record(io::IO, data::Vector{UInt8})
    marker = Int32(length(data))
    write(io, marker)
    write(io, data)
    write(io, marker)
end

@testset "Nbody6Dynamics.jl" begin

    # =====================================================================
    @testset "Configuration" begin
        # Write a minimal config
        cfg_path = joinpath(TESTDIR, "test_config.toml")
        write(
            cfg_path,
            """
[install]
enabled = false
install_dir = "test-nbody"

[build]
enable_mpi = false
enable_hdf5 = false
enable_gpu = false
cuda_path = "/usr/local/cuda"

[simulation]
run_test = false
runs_dir = "test-runs"
run_id_prefix = "test"

[postprocess]
enabled = false
read_escapers = true
escapers_file = "esc.11"
read_stellar_evo = true
stellar_evo_pattern = "sev.83_*"

[visualization]
enabled = false
dpi = 150
figsize = [10, 8]
""",
        )

        cfg = load_config(cfg_path)

        @test cfg.install.enabled == false
        @test cfg.install.install_dir == "test-nbody"
        @test cfg.build.enable_mpi == false
        @test cfg.build.cuda_path == "/usr/local/cuda"
        @test cfg.simulation.run_test == false
        @test cfg.simulation.run_id_prefix == "test"
        @test cfg.postprocess.read_escapers == true
        @test cfg.postprocess.escapers_file == "esc.11"
        @test cfg.postprocess.read_stellar_evo == true
        @test cfg.postprocess.stellar_evo_pattern == "sev.83_*"
        @test cfg.visualization.dpi == 150
        @test cfg.visualization.figsize == (10, 8)
        @test cfg.install.source_url == "https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing.git"
    end

    # =====================================================================
    @testset "Config round-trip (save/load)" begin
        cfg_path = joinpath(TESTDIR, "test_config_rt.toml")
        write(
            cfg_path,
            """
[install]
enabled = true

[build]
enable_gpu = true
cuda_path = "/opt/cuda"

[simulation]
run_id_prefix = "bench"

[postprocess]
enabled = true

[visualization]
dpi = 300
figsize = [12, 9]
""",
        )

        cfg = load_config(cfg_path)
        save_path = joinpath(TESTDIR, "test_config_saved.toml")
        save_config(cfg, save_path)

        cfg2 = load_config(save_path)
        @test cfg2.build.enable_gpu == cfg.build.enable_gpu
        @test cfg2.build.cuda_path == cfg.build.cuda_path
        @test cfg2.simulation.run_id_prefix == cfg.simulation.run_id_prefix
        @test cfg2.visualization.dpi == cfg.visualization.dpi
        @test cfg2.visualization.figsize == cfg.visualization.figsize
    end

    # =====================================================================
    @testset "Config validation" begin
        # Valid config (defaults only) passes validation — smoke test
        ok_path = joinpath(TESTDIR, "val_ok.toml")
        write(
            ok_path,
            """
[simulation]
mpi_ranks = 1

[visualization]
format = "pdf"
""",
        )
        @test load_config(ok_path) isa Nbody6Config

        # Each entry: (label, TOML violating exactly one rule, expected
        # substring of the error message naming the offending section.key).
        bad_cases = [
            ("install_dir", "[install]\ninstall_dir = \"\"\n", "install.install_dir"),
            ("nproc", "[build]\nnproc = -1\n", "build.nproc"),
            ("mpi_ranks", "[simulation]\nmpi_ranks = 0\n", "simulation.mpi_ranks"),
            (
                "mpi_no_mpi_build",
                "[build]\nenable_mpi = false\n\n[simulation]\nmpi_ranks = 4\n",
                "requires build.enable_mpi = true",
            ),
            ("runs_dir", "[simulation]\nruns_dir = \"\"\n", "simulation.runs_dir"),
            ("run_id_prefix", "[simulation]\nrun_id_prefix = \"\"\n", "simulation.run_id_prefix"),
            ("input_file", "[simulation]\ninput_file = \"\"\n", "simulation.input_file"),
            (
                "snapshot_format_hdf5",
                "[postprocess]\nsnapshot_format = \"hdf5\"\n",
                "postprocess.snapshot_format",
            ),
            (
                "snapshot_format_other",
                "[postprocess]\nsnapshot_format = \"parquet\"\n",
                "postprocess.snapshot_format",
            ),
            (
                "snapshot_pattern",
                "[postprocess]\nsnapshot_pattern = \"\"\n",
                "postprocess.snapshot_pattern",
            ),
            (
                "stdout_file",
                "[postprocess]\nparse_stdout = true\nstdout_file = \"\"\n",
                "postprocess.stdout_file",
            ),
            (
                "lagr_file",
                "[postprocess]\nread_lagr = true\nlagr_file = \"\"\n",
                "postprocess.lagr_file",
            ),
            (
                "escapers_file",
                "[postprocess]\nread_escapers = true\nescapers_file = \"\"\n",
                "postprocess.escapers_file",
            ),
            (
                "stellar_evo_pattern",
                "[postprocess]\nread_stellar_evo = true\nstellar_evo_pattern = \"\"\n",
                "postprocess.stellar_evo_pattern",
            ),
            ("format", "[visualization]\nformat = \"gif\"\n", "visualization.format"),
            ("column", "[visualization]\ncolumn = \"triple\"\n", "visualization.column"),
            ("units", "[visualization]\nunits = \"cgs\"\n", "visualization.units"),
            ("dpi", "[visualization]\ndpi = 50\n", "visualization.dpi"),
            ("figsize", "[visualization]\nfigsize = [0.0, 6.0]\n", "visualization.figsize"),
            (
                "marker_budget",
                "[visualization.style]\nmarker_budget = 0.0\n",
                "visualization.style.marker_budget",
            ),
            (
                "marker_bounds",
                "[visualization.style]\nmarker_min = 10.0\nmarker_max = 4.0\n",
                "0 < marker_min ≤ marker_max",
            ),
            (
                "q_log_threshold",
                "[visualization.style]\nq_log_threshold = -1.0\n",
                "visualization.style.q_log_threshold",
            ),
            ("q_floor", "[visualization.style]\nq_floor = 1.5\n", "visualization.style.q_floor"),
            (
                "zoom_frac",
                "[visualization.style]\nzoom_frac = 0.0\n",
                "visualization.style.zoom_frac",
            ),
            ("anim_fps", "[visualization.style]\nanim_fps = -5\n", "visualization.style.anim_fps"),
            (
                "anim_target_seconds",
                "[visualization.style]\nanim_target_seconds = 0.0\n",
                "visualization.style.anim_target_seconds",
            ),
            (
                "merger_config_file",
                "[merger]\nenabled = true\nconfig_file = \"\"\n",
                "merger.config_file",
            ),
        ]
        for (label, body, expected) in bad_cases
            path = joinpath(TESTDIR, "val_bad_$label.toml")
            write(path, body)
            @test_throws expected load_config(path)
        end
    end

    # =====================================================================
    @testset "Platform detection" begin
        platform = detect_platform()
        @test platform isa Symbol
        @test platform in [:fedora, :ubuntu, :debian, :unknown]

        # check_command should find basic tools
        @test Nbody6Dynamics.check_command("ls") == true
        @test Nbody6Dynamics.check_command("nonexistent_tool_xyz") == false
    end

    # =====================================================================
    @testset "CUDA detection" begin
        # detect_cuda_path returns a string (may be empty if no CUDA)
        cuda = detect_cuda_path()
        @test cuda isa String

        # cuda_env_vars returns a dict
        env = Nbody6Dynamics.cuda_env_vars("/usr/local/cuda")
        @test haskey(env, "CUDA_HOME")
        @test env["CUDA_HOME"] == "/usr/local/cuda"
        @test occursin("/usr/local/cuda/bin", env["PATH"])

        # Empty path returns empty dict
        env_empty = Nbody6Dynamics.cuda_env_vars("")
        @test isempty(env_empty)
    end

    # =====================================================================
    @testset "Run ID generation" begin
        id1 = generate_run_id()
        @test startswith(id1, "run_")
        @test length(id1) > 15

        id2 = generate_run_id("bench")
        @test startswith(id2, "bench_")

        # Two successive IDs should differ (hex suffix)
        id3 = generate_run_id()
        @test id1 != id3
    end

    # =====================================================================
    @testset "Fortran binary I/O" begin
        # Create a synthetic Fortran binary file
        fpath = joinpath(TESTDIR, "test_fortran.bin")
        open(fpath, "w") do io
            # Write a record: 3 Float32 values
            data = Float32[1.0, 2.0, 3.0]
            marker = Int32(sizeof(data))
            write(io, marker)
            write(io, data)
            write(io, marker)

            # Write a second record: 2 Int32 values
            idata = Int32[42, 99]
            marker2 = Int32(sizeof(idata))
            write(io, marker2)
            write(io, idata)
            write(io, marker2)
        end

        open(fpath, "r") do io
            rec1 = Nbody6Dynamics.read_fortran_record(io, Float32, 3)
            @test rec1 ≈ Float32[1.0, 2.0, 3.0]

            rec2 = Nbody6Dynamics.read_fortran_record(io, Int32, 2)
            @test rec2 == Int32[42, 99]
        end

        # peek_record_size
        open(fpath, "r") do io
            sz = Nbody6Dynamics.peek_record_size(io)
            @test sz == Int32(12)  # 3 × 4 bytes
        end
    end

    # =====================================================================
    @testset "conf.3 reader" begin
        # Create a synthetic conf.3 file in standard format
        fpath = joinpath(TESTDIR, "conf.3_test")
        n_particles = 5

        open(fpath, "w") do io
            # Record 1: header integers
            hdr = Int32[n_particles, 1, 1, 20]
            _write_fortran_record(io, hdr)

            # Record 2: AS(1:20) parameters
            params = zeros(Float32, 20)
            params[1] = 0.5f0    # time
            params[3] = 1.0f0    # rbar
            params[4] = 0.5f0    # zmbar
            params[11] = 10.0f0  # tscale
            params[12] = 5.0f0   # vstar
            params[18] = 2.0f0   # rscale
            _write_fortran_record(io, params)

            # Particle records (standard: M, X, Y, Z, VX, VY, VZ, NAME)
            for i in 1:n_particles
                m = Float32(1.0 / n_particles)
                x, y, z = Float32.(randn(3))
                vx, vy, vz = Float32.(0.1 .* randn(3))
                name = Int32(i)
                buf = IOBuffer()
                write(buf, m, x, y, z, vx, vy, vz, name)
                _write_fortran_record(io, take!(buf))
            end
        end

        snap = read_conf3(fpath)
        @test nparticles(snap) == n_particles
        @test time_nb(snap.header) ≈ 0.5
        @test tscale(snap.header) ≈ 10.0
        @test rscale(snap.header) ≈ 2.0
        @test length(snap.mass) == n_particles
        @test size(snap.pos) == (3, n_particles)
        @test size(snap.vel) == (3, n_particles)
        @test snap.name == Int32.(1:n_particles)
        @test isempty(snap.rho)  # standard format has no density

        # Test read_all_conf3
        snaps = read_all_conf3(TESTDIR, "conf.3_")
        @test length(snaps) >= 1
    end

    # =====================================================================
    @testset "conf.3 reader — extended format" begin
        fpath = joinpath(TESTDIR, "conf.3_ext")
        n_particles = 3

        open(fpath, "w") do io
            hdr = Int32[n_particles, 1, 1, 20]
            _write_fortran_record(io, hdr)

            params = zeros(Float32, 20)
            params[1] = 1.0f0
            _write_fortran_record(io, params)

            # Extended: M, RHO, XNS, X, Y, Z, VX, VY, VZ, PHI, NAME
            for i in 1:n_particles
                buf = IOBuffer()
                write(buf, Float32(0.1 * i))   # mass
                write(buf, Float32(100.0))      # rho
                write(buf, Float32(0.01))       # xns
                write(buf, Float32.(randn(3))...)  # pos
                write(buf, Float32.(0.1 .* randn(3))...)  # vel
                write(buf, Float32(-0.5))       # phi
                write(buf, Int32(i))            # name
                _write_fortran_record(io, take!(buf))
            end
        end

        snap = read_conf3(fpath)
        @test nparticles(snap) == n_particles
        @test length(snap.rho) == n_particles
        @test length(snap.phi) == n_particles
        @test all(snap.rho .≈ 100.0f0)
        @test all(snap.phi .≈ -0.5f0)
    end

    # =====================================================================
    @testset "Diagnostics parser" begin
        diag_path = joinpath(TESTDIR, "out1000")
        write(
            diag_path,
            """
 Some header text...

 PHYSICAL SCALING:  R* = 1.234  M* = 5678.0  V* = 3.456  T* = 7.890
                    <M> = 0.567  SU = 1.0  AU = 1.0

 ADJUST:    0.0000      0.00  1.000  0.00E+00 -0.2500  10000    500   1.234
 ADJUST:    0.5000     50.00  0.987  1.23E-06 -0.2499   9998    498   1.245
 ADJUST:    1.0000    100.00  0.995  2.45E-06 -0.2498   9990    495   1.267

 END RUN
""",
        )

        diag = read_diagnostics(diag_path)

        @test length(diag.adjust) == 3
        @test diag.adjust[1].time_nb ≈ 0.0
        @test diag.adjust[2].time_myr ≈ 50.0
        @test diag.adjust[3].n == 9990
        @test diag.adjust[2].npairs == 498

        # Physical scaling
        @test haskey(diag.physical_scaling, "R*")
        @test diag.physical_scaling["R*"] ≈ 1.234
        @test diag.physical_scaling["V*"] ≈ 3.456

        # Unit extraction
        units = extract_scaling(diag)
        @test units.rbar ≈ 1.234
        @test units.tscale ≈ 7.890
    end

    # =====================================================================
    @testset "Lagrangian radii reader" begin
        lagr_path = joinpath(TESTDIR, "lagr.7")

        # Synthetic lagr.7 in the upstream lagr.f layout: ##/TIME header,
        # then one row per epoch (3 epochs, 5 radii columns, simplified)
        write(
            lagr_path,
            """
## Time, [Column number | Label]:    2 R_{lagr}   (synthetic)
TIME   0.01  0.05  0.20  0.50  1.00
0.0  0.01  0.05  0.20  0.50  1.00
0.5  0.012 0.055 0.22  0.52  1.05
1.0  0.015 0.060 0.25  0.55  1.10
""",
        )

        lagr = read_lagr(lagr_path)

        @test length(lagr.time) == 3
        @test lagr.time ≈ [0.0, 0.5, 1.0]
        @test size(lagr.radii) == (5, 3)
        @test lagr.radii[1, 1] ≈ 0.01
        @test lagr.radii[5, 3] ≈ 1.10
    end

    # =====================================================================
    @testset "Escaper reader" begin
        # Real esc.11 shape (escape.F): 5 NB-unit diagnostics, then the
        # physical-unit block T[Myr] M[M*] EESC VI[km/s] K* NAME, then the
        # direction angles ANGLE PHI / ANGLE THETA (tokens 12–13).  The last
        # line is truncated to 11 tokens (no angles) and must be skipped.
        esc_path = joinpath(TESTDIR, "esc.11")
        write(
            esc_path,
            """
         TTOT         BODY         RI           VI           STEP         T[Myr]       M[M*]      EESC      VI[km/s]     K*  NAME      ANGLE PHI     ANGLE THETA
   1.00000E+00  5.00000E-04  2.00000E+01  4.00000E+01  1.95312E-03  1.23400E+00  5.00000E-01 -1.23000E-01  1.56000E+01   0       101  1.20000E+01  2.50000E+01
   2.00000E+00  3.00000E-04  2.50000E+01  4.20000E+01  1.95312E-03  2.56700E+00  3.00000E-01  4.56000E-01  2.23000E+01   1       202  1.85000E+02 -3.00000E+01
   3.00000E+00  1.20000E-03  3.00000E+01  4.40000E+01  1.95312E-03  3.89000E+00  1.20000E+00 -7.89000E-01  1.01000E+01  14       303  3.40000E+02  6.00000E+01
   4.00000E+00  1.00000E-03  3.10000E+01  4.50000E+01  1.95312E-03  4.50000E+00  1.00000E+00  1.00000E-01  1.20000E+01   1       404
""",
        )

        escs = read_escapers(esc_path)
        @test length(escs) == 3   # the 11-token line is skipped

        @test escs[1].time_myr ≈ 1.234
        @test escs[1].mass_solar ≈ 0.500
        @test escs[1].escape_energy ≈ -0.123
        @test escs[1].velocity_kms ≈ 15.6
        @test escs[1].stellar_type == 0
        @test escs[1].name == 101
        @test escs[1].phi_deg ≈ 12.0
        @test escs[1].theta_deg ≈ 25.0

        @test escs[2].phi_deg ≈ 185.0
        @test escs[2].theta_deg ≈ -30.0

        @test escs[3].stellar_type == 14   # BH (Hurley convention)
        @test escs[3].name == 303
        @test escs[3].phi_deg ≈ 340.0
        @test escs[3].theta_deg ≈ 60.0

        # Empty file
        empty_path = joinpath(TESTDIR, "esc_empty.11")
        write(empty_path, "")
        @test isempty(read_escapers(empty_path))
    end

    # =====================================================================
    @testset "Stellar evolution reader" begin
        # Header time is TPHYS [Myr]; data-line token 1 is TTOT [NB] —
        # different clocks, both kept. 15-token v2026.07+ layout.
        sev_path = joinpath(TESTDIR, "sev.83_0")
        write(
            sev_path,
            """
  3  0.5000
   0.4315   1   101  0  1.20  0.800   0.123  -0.456  3.750  0.0  0.0  90.0  0.0  0.0  1.0e-10
   0.4315   2   202  1  0.80  1.200   1.500   0.200  4.100  0.0  0.0  5.20  0.0  0.0  1.0e-10
   0.4315   3   303  13 2.50  10.00   5.000   1.500  4.500  0.0  0.0  4.10  1.4  1.0e-5  1.0e-10
""",
        )

        sev = read_stellar_evolution(sev_path)
        @test sev.n_stars == 3
        @test sev.time_myr ≈ 0.5
        @test length(sev.records) == 3

        r1 = sev.records[1]
        @test r1.time_nb ≈ 0.4315   # per-line TTOT [NB], not the Myr header
        @test r1.index == Int32(1)
        @test r1.name == Int32(101)
        @test r1.stellar_type == Int32(0)  # MS
        @test r1.ri ≈ 1.20
        @test r1.mass_solar ≈ 0.800
        @test r1.log_luminosity ≈ 0.123
        @test r1.log_radius ≈ -0.456
        @test r1.log_teff ≈ 3.750
        @test r1.ms_lifetime_myr ≈ 90.0    # TM

        r3 = sev.records[3]
        @test r3.stellar_type == Int32(13)  # NS
        @test r3.mass_core ≈ 1.4

        # Test read_all_stellar_evolution
        sev2_path = joinpath(TESTDIR, "sev.83_1")
        write(
            sev2_path,
            """
  2  1.0000
   1.0000   1   101  2  1.50  0.750   0.500  -0.200  3.600  0.0  0.0  8.0  0.2  0.01  2.0
   1.0000   2   202  4  0.90  1.100   2.000   0.400  3.900  0.0  0.0  6.5  0.4  0.02  4.0
""",
        )

        sevs = read_all_stellar_evolution(TESTDIR, "sev.83_*")
        @test length(sevs) == 2
        @test sevs[1].time_myr < sevs[2].time_myr
    end

    # =====================================================================
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
    @testset "Real-output fixtures" begin
        FIXDIR = joinpath(@__DIR__, "fixtures")

        @testset "esc.11" begin
            escs = read_escapers(joinpath(FIXDIR, "esc.11"))
            @test length(escs) == 88          # every data line parses
            e1 = escs[1]
            @test e1.time_myr ≈ 8.11195 rtol = 1e-5     # T[Myr], not TTOT
            @test e1.mass_solar ≈ 18.1687 rtol = 1e-5   # M[M*], not NB mass
            @test e1.escape_energy ≈ 1017.22 rtol = 1e-5
            @test e1.velocity_kms ≈ 209.841 rtol = 1e-5
            @test e1.stellar_type == 14                 # BH escaper
            @test e1.name == 2248
            @test e1.phi_deg ≈ 22.2302 rtol = 1e-5      # ANGLE PHI [deg]
            @test e1.theta_deg ≈ 28.7010 rtol = 1e-5    # ANGLE THETA [deg]
            # escape.F angle conventions: φ ∈ [0, 360], θ ∈ [-90, 90]
            @test all(0.0 ≤ e.phi_deg ≤ 360.0 for e in escs)
            @test all(-90.0 ≤ e.theta_deg ≤ 90.0 for e in escs)
        end

        @testset "out1000 scaling + ADJUST" begin
            diag = read_diagnostics(joinpath(FIXDIR, "out1000"))
            u = extract_scaling(diag)
            @test u.rbar ≈ 5.503639 rtol = 1e-6
            @test u.zmbar ≈ 27695.934588 rtol = 1e-6   # M* (scale), not <M>
            @test u.tscale ≈ 1.15885054 rtol = 1e-6
            @test u.vstar ≈ 4.65224425 rtol = 1e-6
            @test length(diag.adjust) ≥ 6
            a0 = diag.adjust[1]
            @test a0.time_nb == 0.0
            @test a0.qvir ≈ 0.321 atol = 1e-3           # Q = T/|W|
            @test a0.e_tot ≈ -0.5408 atol = 1e-4
        end

        @testset "lagr.7" begin
            lagr = read_lagr(joinpath(FIXDIR, "lagr.7"))
            @test length(lagr.mass_fractions) == 18
            @test lagr.time[1] == 0.0
            i50 = argmin(abs.(lagr.mass_fractions .- 0.5))
            @test lagr.radii[i50, 1] ≈ 1.2442064 rtol = 1e-6
        end

        @testset "sev.83" begin
            sev = read_stellar_evolution(joinpath(FIXDIR, "sev.83_0"))
            @test sev.n_stars == 40
            @test sev.time_myr == 0.0
            @test length(sev.records) == 40
            r1 = sev.records[1]
            @test r1.time_nb == 0.0
            @test r1.name == Int32(1)
            @test r1.stellar_type == Int32(1)            # massive MS star, Hurley K*=1
            @test r1.ri ≈ 0.694596 rtol = 1e-5           # RI [pc]
            @test r1.mass_solar ≈ 52.5373 rtol = 1e-5
            @test r1.log_teff ≈ 4.71026 rtol = 1e-5
            @test r1.ms_lifetime_myr ≈ 4.46907 rtol = 1e-5   # TM
            @test r1.mass_core == 0.0                    # MC (MS star)
            @test r1.radius_envelope ≈ 1e-10 rtol = 1e-3 # RE placeholder on MS
        end
    end

    # =====================================================================
    @testset "Types and accessors" begin
        params = zeros(Float32, 20)
        params[1] = 1.5f0   # time
        params[3] = 2.0f0   # rbar
        params[4] = 0.6f0   # zmbar
        params[11] = 8.0f0  # tscale
        params[12] = 4.0f0  # vstar

        hdr = SnapshotHeader(Int32(100), Int32(1), Int32(1), Int32(20), params)

        @test time_nb(hdr) ≈ 1.5
        @test rbar(hdr) ≈ 2.0
        @test tscale(hdr) ≈ 8.0
        @test time_myr(hdr) ≈ 12.0   # 1.5 * 8.0

        u = UnitScaling(2.0, 0.6, 8.0, 4.0)
        @test Nbody6Dynamics.to_pc(u, 1.0) ≈ 2.0
        @test Nbody6Dynamics.to_msun(u, 1.0) ≈ 0.6
        @test Nbody6Dynamics.to_myr(u, 1.0) ≈ 8.0
        @test Nbody6Dynamics.to_kms(u, 1.0) ≈ 4.0
    end

    # =====================================================================
    @testset "Plotting (smoke tests)" begin
        Nbody6Dynamics.set_publication_theme!()

        vis = VisualizationConfig(;
            output_dir = joinpath(TESTDIR, "test_plots"),
            format = "png",
            dpi = 72,
            figsize = (6, 4),
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
    end

    # =====================================================================
    @testset "Escaper and SSE plots (smoke tests)" begin
        Nbody6Dynamics.set_publication_theme!()

        plots_dir = joinpath(TESTDIR, "test_plots_f23")
        vis = VisualizationConfig(;
            output_dir = plots_dir,
            format = "png",
            dpi = 72,
            figsize = (6, 4),
        )

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
    @testset "Elapsed time formatting" begin
        fmt = Nbody6Dynamics._format_elapsed

        # Sub-minute
        @test fmt(0.0) == "0.0 s"
        @test fmt(1.5) == "1.5 s"
        @test fmt(59.4) == "59.4 s"

        # Minutes
        @test fmt(60.0) == "1m 00s"
        @test fmt(122.0) == "2m 02s"
        @test fmt(3599.0) == "59m 59s"

        # Hours
        @test fmt(3600.0) == "1h 00m 00s"
        @test fmt(3661.0) == "1h 01m 01s"
        @test fmt(7384.0) == "2h 03m 04s"
    end

    # =====================================================================
    @testset "Diagnostics — key-value ADJUST format" begin
        diag_path = joinpath(TESTDIR, "out1000_kv")
        write(
            diag_path,
            """
 PHYSICAL SCALING:  R* = 2.500  M* = 1000.0  V* = 4.200  T* = 12.00
                    <M> = 0.500  SU = 1.0  AU = 1.0

 ADJUST:  TIME   0.0000  T[Myr]      0.00  Q  1.000  DE  0.00E+00  ETOT  -0.2500
 RMIN =    0.001 RSCALE =    1.234
 TIME[NB]    0.0000 N    10000 <NB>      0 NPAIRS    500

 ADJUST:  TIME   0.5000  T[Myr]     50.00  Q  0.987  DE  1.23E-06  ETOT  -0.2499
 RMIN =    0.002 RSCALE =    1.300
 TIME[NB]    0.5000 N     9990 <NB>      0 NPAIRS    495

 END RUN
""",
        )

        diag = read_diagnostics(diag_path)

        @test length(diag.adjust) == 2
        # First epoch
        @test diag.adjust[1].time_nb ≈ 0.0
        @test diag.adjust[1].time_myr ≈ 0.0
        @test diag.adjust[1].qvir ≈ 1.0
        @test diag.adjust[1].n == 10000
        @test diag.adjust[1].npairs == 500
        @test diag.adjust[1].rscale ≈ 1.234

        # Second epoch — TIME[NB] values override ADJUST defaults
        @test diag.adjust[2].time_nb ≈ 0.5
        @test diag.adjust[2].n == 9990
        @test diag.adjust[2].npairs == 495
        @test diag.adjust[2].rscale ≈ 1.300

        # Physical scaling
        @test diag.physical_scaling["R*"] ≈ 2.5
        @test diag.physical_scaling["T*"] ≈ 12.0
    end

    # =====================================================================
    @testset "Diagnostics — mixed format with partial epochs" begin
        diag_path = joinpath(TESTDIR, "out1000_mixed")
        write(
            diag_path,
            """
 ADJUST:    0.0000      0.00  1.000  0.00E+00 -0.2500  10000    500   1.234
 ADJUST:  TIME   1.0000  T[Myr]    100.00  Q  0.990  DE  5.00E-06  ETOT  -0.2480
 TIME[NB]    1.0000 N     9500 <NB>      0 NPAIRS    450
""",
        )
        diag = read_diagnostics(diag_path)
        @test length(diag.adjust) == 2

        # First: positional format
        @test diag.adjust[1].time_nb ≈ 0.0
        @test diag.adjust[1].n == 10000
        @test diag.adjust[1].rscale ≈ 1.234

        # Second: key-value + TIME[NB] merge
        @test diag.adjust[2].time_nb ≈ 1.0
        @test diag.adjust[2].n == 9500
        @test diag.adjust[2].npairs == 450
    end

    # =====================================================================
    @testset "Auto FPS calculation" begin
        _auto_fps = Nbody6Dynamics._auto_fps

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
            figsize = (4, 3),
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
    @testset "External post-processing" begin
        # Create a fake output directory with some recognisable files
        ext_dir = joinpath(TESTDIR, "external_output")
        mkpath(ext_dir)

        # --- scan_output on empty directory ---
        scan_empty = scan_output(ext_dir)
        @test scan_empty.dir == abspath(ext_dir)
        @test scan_empty.available[:snapshots_conf3] == false
        @test scan_empty.available[:diagnostics] == false
        @test scan_empty.available[:lagr] == false
        @test scan_empty.available[:escapers] == false
        @test scan_empty.available[:stellar_evo] == false
        @test scan_empty.available[:snapshots_hdf5] == false

        # --- Populate with dummy files ---
        # Write a minimal valid out1000
        open(joinpath(ext_dir, "out1000"), "w") do io
            println(io, " ADJUST: TIME  1.0  T[MYR]  0.5  Q  0.5  DE  1e-8  E  -0.25")
        end
        # Write a dummy lagr.7
        open(joinpath(ext_dir, "lagr.7"), "w") do io
            println(io, "# TIME  0.01  0.10  0.50  1.00")
            println(io, "0.0  0.1  0.3  1.0  5.0")
            println(io, "1.0  0.1  0.3  1.1  5.2")
        end
        # Write a dummy esc.11
        open(joinpath(ext_dir, "esc.11"), "w") do io
            println(io, "# escaper data")
            println(io, "1  1.0  0.5  1.0 2.0 3.0  0.1 0.2 0.3")
        end

        # --- scan_output with partial data ---
        scan = scan_output(ext_dir)
        @test scan.available[:diagnostics] == true
        @test scan.available[:lagr] == true
        @test scan.available[:escapers] == true
        @test scan.available[:snapshots_conf3] == false
        @test scan.available[:stellar_evo] == false

        # --- OutputScan display ---
        buf = IOBuffer()
        show(buf, MIME("text/plain"), scan)
        output_str = String(take!(buf))
        @test occursin("Diagnostics", output_str)
        @test occursin("out1000", output_str)
        @test occursin("not found", output_str)  # for missing categories

        # --- scan_output error on non-existent directory ---
        @test_throws ErrorException scan_output("/nonexistent/path")

        # --- postprocess_external (data only, no plots) ---
        results = postprocess_external(ext_dir; make_plots = false)
        @test haskey(results, :scan)
        @test results[:scan] isa OutputScan
        @test haskey(results, :diagnostics)
        @test haskey(results, :lagr)
    end

    # =====================================================================
    # Adversarial external post-processing tests (included from separate file)
    # =====================================================================
    include("test_external_adversarial_inner.jl")

    # =====================================================================
    @testset "Initial Conditions — Merger IC Generator" begin
        # StableRNGs guarantees an identical stream across Julia versions, so
        # reference values in these tests survive upgrades (§8).
        using StableRNGs
        rng = StableRNG(12345)

        # --- Plummer sampler ---
        @testset "Plummer sampler" begin
            N = 1000
            a = 1.0
            pos, vel = sample_plummer(N, a; rng = rng)
            @test size(pos) == (3, N)
            @test size(vel) == (3, N)
            # Particles should be centred near origin
            cm = vec(sum(pos, dims = 2)) ./ N
            @test all(abs.(cm) .< 0.5)
            # Radii should be finite and positive
            radii = [sqrt(sum(pos[:, i] .^ 2)) for i in 1:N]
            @test all(radii .> 0)
            @test all(isfinite.(radii))
            # Half-mass radius should be roughly 1.305 × a for Plummer
            r_sorted = sort(radii)
            r_hm = r_sorted[N ÷ 2]
            @test 0.5 * a < r_hm < 3.0 * a
        end

        # --- King sampler ---
        @testset "King sampler" begin
            N = 1000
            W0 = 6.0
            rt = 10.0
            pos, vel = sample_king(N, W0, rt; rng = rng)
            @test size(pos) == (3, N)
            @test size(vel) == (3, N)
            # All particles should be within the tidal radius (with some tolerance
            # from the sampling/scaling)
            radii = [sqrt(sum(pos[:, i] .^ 2)) for i in 1:N]
            @test all(isfinite.(radii))
            # Invalid W0
            @test_throws ArgumentError sample_king(10, -1.0, 5.0)
            @test_throws ArgumentError sample_king(10, 6.0, -1.0)
        end

        # --- King ODE solver ---
        @testset "King ODE solver" begin
            rhat, What, rho = Nbody6Dynamics._solve_king(6.0)
            @test length(rhat) == length(What) == length(rho)
            @test What[1] ≈ 6.0
            @test What[end] ≈ 0.0 atol = 0.01
            @test rhat[1] ≈ 0.0
            @test rhat[end] > 0.0
            # Density should be monotonically decreasing
            @test all(diff(rho[1:(end - 1)]) .≤ 0.01)
        end

        # --- King concentration vs published values (physics validation) ---
        @testset "King concentration c(W0)" begin
            # c = log10(r_t/r_0) for the standard dimensionless King equation;
            # reference values from the King (1966) model tables.
            for (W0, c_ref) in [
                (3.0, 0.672),
                (5.0, 1.029),
                (6.0, 1.255),
                (7.0, 1.528),
                (9.0, 2.119),
                (12.0, 2.739),
            ]
                rhat, _, _ = Nbody6Dynamics._solve_king(W0)
                c = log10(rhat[end])
                @test isapprox(c, c_ref; rtol = 0.01)
            end
        end

        # --- Kroupa IMF ---
        @testset "Kroupa IMF" begin
            masses = sample_kroupa(5000; m_low = 0.08, m_up = 100.0, rng = rng)
            @test length(masses) == 5000
            @test all(masses .≥ 0.08)
            @test all(masses .≤ 100.0)
            # Most stars should be low-mass (IMF is bottom-heavy)
            @test count(m -> m < 1.0, masses) > 3000
            # Custom mass limits
            masses_narrow = sample_kroupa(500; m_low = 1.0, m_up = 10.0, rng = rng)
            @test all(1.0 .≤ masses_narrow .≤ 10.0)
        end

        # --- Virialisation ---
        @testset "Virialise" begin
            N = 200
            mass = ones(N) ./ N
            pos = randn(rng, 3, N)
            vel = randn(rng, 3, N) .* 0.1
            virialise!(mass, pos, vel)
            # CM should be at origin
            cm_pos = vec(sum(mass' .* pos, dims = 2))
            cm_vel = vec(sum(mass' .* vel, dims = 2))
            @test all(abs.(cm_pos) .< 1e-10)
            @test all(abs.(cm_vel) .< 1e-10)
            # Virial ratio should be ~0.5
            T = 0.5 * sum(mass[i] * sum(vel[:, i] .^ 2) for i in 1:N)
            W = 0.0
            for i in 1:N, j in (i + 1):N
                dr = pos[:, i] .- pos[:, j]
                W -= mass[i] * mass[j] / sqrt(sum(dr .^ 2))
            end
            Q = T / abs(W)
            @test Q ≈ 0.5 atol = 0.01
        end

        # --- Kepler velocity ---
        @testset "Kepler velocity" begin
            # Equal mass, circular orbit
            v1, v2 = kepler_velocity(1.0, 1.0, 10.0, 0.0)
            @test v1 ≈ v2  # symmetric
            @test v1 > 0
            # Eccentric orbit has lower apocentre velocity
            v1e, v2e = kepler_velocity(1.0, 1.0, 10.0, 0.9)
            @test v1e < v1
            # Invalid inputs
            @test_throws ArgumentError kepler_velocity(1.0, 1.0, 10.0, 1.0)
            @test_throws ArgumentError kepler_velocity(1.0, 1.0, -1.0, 0.5)
        end

        # --- Jacobi radius ---
        @testset "Jacobi radius" begin
            rj = jacobi_radius(10.0, 1.0, 3.0)
            @test rj ≈ 10.0 * (1.0 / 9.0)^(1/3)
            # Equal mass: rJ = d × (1/3)^{1/3}
            rj_eq = jacobi_radius(10.0, 1.0, 1.0)
            @test rj_eq ≈ 10.0 * (1/3)^(1/3) atol = 1e-10
        end

        # --- Two-cluster orbit setup ---
        @testset "Two-cluster orbit" begin
            N1, N2 = 100, 100
            pos1 = randn(rng, 3, N1) .* 0.5
            vel1 = randn(rng, 3, N1) .* 0.01
            mass1 = ones(N1) ./ N1
            pos2 = randn(rng, 3, N2) .* 0.5
            vel2 = randn(rng, 3, N2) .* 0.01
            mass2 = ones(N2) ./ N2
            # Centre each cluster exactly (as virialise! does in the real
            # pipeline) so the COM/momentum/energy identities below are exact.
            for (p, v, m) in ((pos1, vel1, mass1), (pos2, vel2, mass2))
                p .-= sum(m' .* p, dims = 2) ./ sum(m)
                v .-= sum(m' .* v, dims = 2) ./ sum(m)
            end

            pos, vel, mass, ranges = Nbody6Dynamics.setup_two_cluster_orbit(
                pos1,
                vel1,
                mass1,
                pos2,
                vel2,
                mass2,
                20.0,
                0.5;
                truncate_jacobi_flag = false,
            )
            @test length(mass) == N1 + N2
            @test size(pos, 2) == N1 + N2
            @test ranges == [1:N1, (N1 + 1):(N1 + N2)]
            # Combined CM should be near origin
            M = sum(mass)
            cm = vec(sum(mass' .* pos, dims = 2)) ./ M
            @test all(abs.(cm) .< 0.1)

            # Total momentum must vanish and the two-COM orbit must have the
            # Keplerian energy of the requested (d_apo, e) orbit.
            p = vec(sum(mass' .* vel, dims = 2))
            @test all(abs.(p) .< 1e-10)
            r1 = 1:N1
            r2 = (N1 + 1):(N1 + N2)
            M1 = sum(mass[r1])
            M2 = sum(mass[r2])
            com1 = vec(sum(mass[r1]' .* pos[:, r1], dims = 2)) ./ M1
            com2 = vec(sum(mass[r2]' .* pos[:, r2], dims = 2)) ./ M2
            vcom1 = vec(sum(mass[r1]' .* vel[:, r1], dims = 2)) ./ M1
            vcom2 = vec(sum(mass[r2]' .* vel[:, r2], dims = 2)) ./ M2
            d = sqrt(sum((com1 .- com2) .^ 2))
            @test d ≈ 20.0 rtol = 1e-10
            v_rel2 = sum((vcom1 .- vcom2) .^ 2)
            a_orb = 20.0 / (1.0 + 0.5)
            ε = 0.5 * v_rel2 - (M1 + M2) / d          # specific orbital energy, G = 1
            @test ε ≈ -(M1 + M2) / (2a_orb) rtol = 1e-10
        end

        # --- Plummer half-mass relation (physics validation) ---
        @testset "Plummer r_hm = 1.305 a" begin
            N = 20_000
            a = 1.0
            pos, _ = sample_plummer(N, a; rng = rng)
            mass = fill(1.0 / N, N)
            r_hm = Nbody6Dynamics.half_mass_radius(mass, pos; centre = zeros(3))
            @test r_hm ≈ 1.3048 rtol = 0.05    # statistical tolerance
        end

        # --- Kroupa mean mass vs analytic (physics validation) ---
        @testset "Kroupa mean mass" begin
            @test kroupa_mean_mass(0.08, 100.0) ≈ 0.58 rtol = 0.05
            N = 50_000
            m = sample_kroupa(N; m_low = 0.08, m_up = 100.0, rng = rng)
            μ = kroupa_mean_mass(0.08, 100.0)
            # Sample mean within a generous statistical band of the analytic mean
            @test abs(sum(m) / N - μ) / μ < 0.1
        end

        # --- verif_triorbit.toml Lagrange equilibrium (config regression) ---
        @testset "verif_triorbit config equilibrium" begin
            path = joinpath(@__DIR__, "..", "input_files", "verif_triorbit.toml")
            cfg = load_merger_config(path)
            @test length(cfg.clusters) == 3
            m = 1.5e4
            rc = 6.0
            ω = sqrt(m / (sqrt(3) * rc^3))
            for spec in cfg.clusters
                r⃗ = spec.position
                v⃗ = spec.velocity
                @test sqrt(sum(r⃗ .^ 2)) ≈ rc rtol = 1e-3
                # Speed = ω r and velocity ⟂ radius (circular Lagrange orbit)
                @test sqrt(sum(v⃗ .^ 2)) ≈ ω * rc rtol = 1e-3
                @test abs(sum(r⃗ .* v⃗)) / (rc * ω * rc) < 1e-3
            end
            # Zero net momentum for equal-mass clusters
            vsum = sum(spec.velocity for spec in cfg.clusters)
            @test all(abs.(vsum) .< 0.05 * ω * rc)
        end

        # --- dat.10 writer ---
        @testset "dat.10 writer" begin
            N = 50
            mass = ones(N) ./ N
            pos = randn(rng, 3, N)
            vel = randn(rng, 3, N)
            dat10_path = joinpath(TESTDIR, "test_dat10.dat")
            write_dat10(dat10_path, mass, pos, vel)
            lines = readlines(dat10_path)
            @test length(lines) == N
            # Each line should have 7 columns
            cols = split(lines[1])
            @test length(cols) == 7
            # First column should be mass
            @test parse(Float64, cols[1]) ≈ 1.0 / N
        end

        # --- N-body unit conversion ---
        @testset "N-body unit conversion" begin
            mass = [0.5, 0.5]
            pos = [1.0 -1.0; 0.0 0.0; 0.0 0.0]
            vel = [0.0 0.0; 1.0 -1.0; 0.0 0.0]
            M_tot = 1e5  # M☉
            rbar = 2.0   # pc
            to_nbody_units!(mass, pos, vel, M_tot, rbar)
            @test sum(mass) ≈ 1e-5  # each mass = 0.5/1e5
            @test pos[1, 1] ≈ 0.5   # 1.0 / 2.0
        end

        # --- TOML config loading (kepler mode) ---
        @testset "Merger config TOML — kepler" begin
            cfg_path = joinpath(TESTDIR, "test_merger.toml")
            write(
                cfg_path,
                """
[merger]
n_clusters = 2
orbit_mode = "kepler"

[merger.cluster1]
model = "plummer"
N = 100
mass_total = 1000.0
rbar = 1.0

[merger.cluster2]
model = "king"
N = 100
W0 = 5.0
mass_total = 1000.0
rbar = 1.0

[merger.orbit]
apocentre = 10.0
eccentricity = 0.5

[merger.output]
format = "nbody"
truncate_jacobi = false
""",
            )
            cfg = load_merger_config(cfg_path)
            @test length(cfg.clusters) == 2
            @test cfg.orbit_mode == "kepler"
            @test cfg.clusters[1].profile isa PlummerProfile
            @test cfg.clusters[2].profile isa KingProfile
            @test cfg.clusters[2].profile.W0 == 5.0
            @test cfg.orbit.apocentre == 10.0
            @test cfg.orbit.eccentricity == 0.5
            @test cfg.output.format == "nbody"
            @test cfg.output.truncate_jacobi == false
        end

        # --- TOML config loading (explicit mode) ---
        @testset "Merger config TOML — explicit" begin
            cfg_path = joinpath(TESTDIR, "test_merger_explicit.toml")
            write(
                cfg_path,
                """
[merger]
n_clusters = 3
orbit_mode = "explicit"

[merger.cluster1]
model = "plummer"
N = 50
mass_total = 500.0
rbar = 1.0
position = [-5.0, 0.0, 0.0]
velocity = [1.0, 0.0, 0.0]

[merger.cluster2]
model = "king"
N = 50
W0 = 4.0
mass_total = 500.0
rbar = 1.0
position = [2.5, 4.33, 0.0]
velocity = [-0.5, -0.87, 0.0]

[merger.cluster3]
model = "plummer"
N = 50
mass_total = 300.0
rbar = 0.8
position = [2.5, -4.33, 0.0]
velocity = [-0.5, 0.87, 0.0]

[merger.output]
format = "nbody"
truncate_jacobi = false
""",
            )
            cfg = load_merger_config(cfg_path)
            @test length(cfg.clusters) == 3
            @test cfg.orbit_mode == "explicit"
            @test cfg.clusters[1].position == [-5.0, 0.0, 0.0]
            @test cfg.clusters[3].velocity == [-0.5, 0.87, 0.0]
        end

        # --- Fail-fast validation of merger TOML values ---
        @testset "Merger config validation" begin
            # Template with one overridable block per constrained section
            function _merger_toml(;
                cluster1::String = "model = \"king\"\nN = 100\nW0 = 5.0\nrbar = 1.0",
                orbit::String = "apocentre = 10.0\neccentricity = 0.5",
                output::String = "tcrit = 10.0\ndtadj = 1.0\ndeltat = 1.0",
            )
                """
                [merger]
                n_clusters = 2
                orbit_mode = "kepler"

                [merger.cluster1]
                $cluster1

                [merger.cluster2]
                model = "plummer"
                N = 100
                rbar = 1.0

                [merger.orbit]
                $orbit

                [merger.output]
                $output
                """
            end

            # Valid template loads cleanly — smoke test
            ok_path = joinpath(TESTDIR, "merger_val_ok.toml")
            write(ok_path, _merger_toml())
            @test load_merger_config(ok_path) isa MergerConfig

            cases = [
                (
                    "N",
                    _merger_toml(cluster1 = "model = \"king\"\nN = 1\nrbar = 1.0"),
                    "merger.cluster1.N",
                ),
                (
                    "rbar",
                    _merger_toml(cluster1 = "model = \"king\"\nN = 100\nrbar = -1.0"),
                    "merger.cluster1.rbar",
                ),
                (
                    "W0",
                    _merger_toml(cluster1 = "model = \"king\"\nN = 100\nW0 = -2.0\nrbar = 1.0"),
                    "merger.cluster1.W0",
                ),
                (
                    "kroupa_bounds",
                    _merger_toml(
                        cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n" *
                                   "imf = \"kroupa\"\nbodyn = 50.0\nbody1 = 0.1",
                    ),
                    "0 < bodyn < body1",
                ),
                (
                    "target_mass",
                    _merger_toml(
                        cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n" *
                                   "imf = \"kroupa\"\nmass_total = -500.0",
                    ),
                    "target_mass",
                ),
                (
                    "particle_mass",
                    _merger_toml(
                        cluster1 = "model = \"king\"\nN = 100\nrbar = 1.0\n" *
                                   "imf = \"equal\"\nmass_total = 0.0",
                    ),
                    "particle_mass",
                ),
                (
                    "apocentre",
                    _merger_toml(orbit = "apocentre = -5.0\neccentricity = 0.5"),
                    "merger.orbit.apocentre",
                ),
                (
                    "eccentricity",
                    _merger_toml(orbit = "apocentre = 10.0\neccentricity = 1.0"),
                    "merger.orbit.eccentricity",
                ),
                (
                    "tcrit",
                    _merger_toml(output = "tcrit = 0.0\ndtadj = 1.0\ndeltat = 1.0"),
                    "merger.output.tcrit",
                ),
                (
                    "dtadj",
                    _merger_toml(output = "tcrit = 10.0\ndtadj = -0.5\ndeltat = 1.0"),
                    "merger.output.dtadj",
                ),
                (
                    "deltat",
                    _merger_toml(output = "tcrit = 10.0\ndtadj = 1.0\ndeltat = 0.0"),
                    "merger.output.deltat",
                ),
            ]
            for (label, body, expected) in cases
                path = joinpath(TESTDIR, "merger_val_bad_$label.toml")
                write(path, body)
                @test_throws expected load_merger_config(path)
            end
        end

        # --- Full pipeline — kepler mode (small N) ---
        @testset "Full merger IC pipeline — kepler" begin
            out_dir = mktempdir()
            cfg = MergerConfig(
                [
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 100,
                        rbar = 1.0,
                        imf = RescaledKroupaIMF(bodyn = 0.1, body1 = 50.0, target_mass = 1e3),
                    ),
                    ClusterSpec(
                        profile = KingProfile(W0 = 5.0),
                        N = 100,
                        rbar = 1.0,
                        imf = RescaledKroupaIMF(bodyn = 0.1, body1 = 50.0, target_mass = 1e3),
                    ),
                ],
                "kepler",
                OrbitSpec(apocentre = 10.0, eccentricity = 0.5),
                MergerOutputSpec(format = "nbody", truncate_jacobi = true, output_dir = out_dir),
            )
            result = generate_merger_ic(cfg; rng = rng)
            @test result isa MergerICResult
            @test isdir(result.output_dir)
            @test isfile(joinpath(result.output_dir, "dat.10"))
            @test isfile(joinpath(result.output_dir, "merger.inp"))
            @test isfile(joinpath(result.output_dir, "merger_summary.txt"))
            @test result.N_total > 0
            @test result.M_total > 0
            @test result.orbit_mode == "kepler"
            @test length(result.cluster_ranges) == 2

            # Verify dat.10 mass normalisation
            lines = readlines(joinpath(result.output_dir, "dat.10"))
            masses = [parse(Float64, split(l)[1]) for l in lines]
            @test sum(masses) ≈ 1.0 atol = 1e-10

            # Verify .inp has KZ(22)=2
            inp_text = read(joinpath(result.output_dir, "merger.inp"), String)
            @test occursin("KZ(21:30)", inp_text)
        end

        # --- Full pipeline — explicit mode (3 clusters) ---
        @testset "Full merger IC pipeline — explicit 3-cluster" begin
            out_dir = mktempdir()
            cfg = MergerConfig(
                [
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 80,
                        rbar = 1.0,
                        imf = RescaledKroupaIMF(bodyn = 0.1, body1 = 50.0, target_mass = 800.0),
                        position = [-5.0, 0.0, 0.0],
                        velocity = [1.0, 0.0, 0.0],
                    ),
                    ClusterSpec(
                        profile = KingProfile(W0 = 5.0),
                        N = 80,
                        rbar = 1.0,
                        imf = RescaledKroupaIMF(bodyn = 0.1, body1 = 50.0, target_mass = 800.0),
                        position = [2.5, 4.33, 0.0],
                        velocity = [-0.5, -0.87, 0.0],
                    ),
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 60,
                        rbar = 0.8,
                        imf = EqualMassIMF(particle_mass = 400.0 / 60),
                        position = [2.5, -4.33, 0.0],
                        velocity = [-0.5, 0.87, 0.0],
                    ),
                ],
                "explicit",
                OrbitSpec(),  # ignored in explicit mode
                MergerOutputSpec(format = "nbody", truncate_jacobi = false, output_dir = out_dir),
            )
            result = generate_merger_ic(cfg; rng = rng)
            @test result isa MergerICResult
            @test result.orbit_mode == "explicit"
            @test length(result.cluster_ranges) == 3
            @test result.N_total == 80 + 80 + 60
            # Total mass should sum to 1 in N-body units
            lines = readlines(joinpath(result.output_dir, "dat.10"))
            @test sum(parse(Float64, split(l)[1]) for l in lines) ≈ 1.0 atol = 1e-10
        end

        # --- Equal-mass IMF ---
        @testset "Equal mass IMF" begin
            out_dir = mktempdir()
            cfg = MergerConfig(
                [
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 50,
                        rbar = 1.0,
                        imf = EqualMassIMF(particle_mass = 500.0 / 50),
                    ),
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 50,
                        rbar = 1.0,
                        imf = EqualMassIMF(particle_mass = 500.0 / 50),
                    ),
                ],
                "kepler",
                OrbitSpec(apocentre = 8.0, eccentricity = 0.3),
                MergerOutputSpec(format = "nbody", truncate_jacobi = false, output_dir = out_dir),
            )
            generate_merger_ic(cfg; rng = rng)
            lines = readlines(joinpath(out_dir, "dat.10"))
            @test length(lines) == 100
        end

        # --- Summary write→parse round-trip (guards the format/regex coupling) ---
        @testset "Merger summary round-trip" begin
            out_dir = mktempdir()
            cfg = MergerConfig(
                [
                    ClusterSpec(
                        profile = KingProfile(W0 = 5.0),
                        N = 150,
                        rbar = 1.0,
                        imf = RescaledKroupaIMF(bodyn = 0.1, body1 = 50.0, target_mass = 1e3),
                    ),
                    ClusterSpec(
                        profile = PlummerProfile(),
                        N = 100,
                        rbar = 1.0,
                        imf = EqualMassIMF(particle_mass = 8e2 / 100),
                    ),
                ],
                "kepler",
                OrbitSpec(apocentre = 10.0, eccentricity = 0.4),
                MergerOutputSpec(format = "nbody", truncate_jacobi = true, output_dir = out_dir),
            )
            result = generate_merger_ic(cfg; rng = rng)
            ranges = parse_merger_summary(joinpath(out_dir, "merger_summary.txt"))
            @test ranges == result.cluster_ranges
            @test sum(length, ranges) == result.N_total
        end

        # --- MergerPipelineConfig in Nbody6Config ---
        @testset "MergerPipelineConfig in load_config" begin
            cfg = load_config(joinpath(TESTDIR, "test_config.toml"))
            @test cfg.merger.enabled == false
            @test cfg.merger.config_file == ""
        end
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
                    imf = RescaledKroupaIMF(bodyn = 0.1, body1 = 50.0, target_mass = 600.0),
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

        # Envelope statistics helper: NaNs are skipped, all-NaN columns stay NaN
        env = [1.0 NaN 3.0; 5.0 NaN 1.0]
        lo, hi, mean_vals = Nbody6Dynamics._envelope_stats(env)
        @test lo == [1.0, NaN, 1.0] || (lo[1] == 1.0 && isnan(lo[2]) && lo[3] == 1.0)
        @test hi[1] == 5.0 && isnan(hi[2]) && hi[3] == 3.0
        @test mean_vals[1] == 3.0 && isnan(mean_vals[2]) && mean_vals[3] == 2.0
    end

    # =====================================================================
    @testset "export_for_paper provenance" begin
        src_dir = mktempdir()
        dest = mktempdir()
        fig_path = joinpath(src_dir, "plots", "energy.pdf")
        mkpath(dirname(fig_path))
        write(fig_path, "pdfbytes")
        write(
            joinpath(src_dir, "RUN_INFO.txt"),
            "Run ID:    testrun_x\nCommit:    abc1234\nBackend:   def5678-dirty\n",
        )
        out = export_for_paper([fig_path], dest; run_dir = src_dir)
        @test length(out) == 1
        @test isfile(out[1])
        @test startswith(basename(out[1]), basename(src_dir) * "__")
        sidecar = out[1] * ".provenance.toml"
        @test isfile(sidecar)
        prov = Nbody6Dynamics.TOML.parsefile(sidecar)
        @test prov["package_commit"] == "abc1234"
        @test prov["backend_commit"] == "def5678-dirty"
        # Never-overwrite: second export backs up, both exist
        export_for_paper([fig_path], dest; run_dir = src_dir)
        @test length(filter(f -> endswith(f, ".pdf"), readdir(dest))) ≥ 2
    end

    # =====================================================================
    # Edge cases for pure helpers and degenerate reader inputs
    # =====================================================================
    @testset "Log-tick generator edge cases" begin
        # >2 in-range decades → decades only
        vals, _ = Nbody6Dynamics._log_ticks(0.05, 50.0)
        @test vals == [0.1, 1.0, 10.0]
        # ≤2 in-range decades → 2×/5× intermediates appear
        vals2, labels2 = Nbody6Dynamics._log_ticks(0.5, 30.0)
        @test all(v -> v in vals2, (0.5, 1.0, 2.0, 5.0, 10.0, 20.0))
        # 10^0 renders as plain "1"
        lab1 = String(labels2[findfirst(==(1.0), vals2)])
        @test occursin("1", lab1) && !occursin("10", lab1)
        # Degenerate equal endpoints still yield ≥ 2 ticks
        vals3, _ = Nbody6Dynamics._log_ticks(2.0, 2.0)
        @test length(vals3) ≥ 2
    end

    @testset "Degenerate reader inputs" begin
        # sev.83 with zero stars: header-only file parses to an empty snapshot
        p_empty = joinpath(TESTDIR, "sev.83_empty")
        write(p_empty, "  0  0.0\n")
        sev0 = read_stellar_evolution(p_empty)
        @test sev0.n_stars == 0 && isempty(sev0.records)

        # esc.11 line with only 12 tokens (angle column truncated) is skipped
        p_trunc = joinpath(TESTDIR, "esc_trunc.11")
        write(p_trunc, "  1.0 5.0e-4 20.0 40.0 2.0e-3 1.234 0.5 -0.1 15.6 0 101 22.2\n")
        @test isempty(read_escapers(p_trunc))
    end

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
    # Static QA (§8): ships with the tests.
    # =====================================================================
    @testset "Static QA — Aqua" begin
        using Aqua
        # Method ambiguities are checked for this package only — recursing
        # into the Makie/SciML dependency tree reports upstream noise.
        Aqua.test_all(Nbody6Dynamics; ambiguities = false, persistent_tasks = false)
        Aqua.test_ambiguities(Nbody6Dynamics)
    end

    @testset "Static QA — ExplicitImports" begin
        using ExplicitImports
        # Pragmatic subset: the package uses plain `using` for its small,
        # stable dependency surface (full explicit-import migration is
        # tracked in the roadmap). These checks catch the real hazards:
        # stale explicit imports, self-qualified names, and accesses of
        # non-owning modules.
        @test check_no_stale_explicit_imports(Nbody6Dynamics) === nothing
        @test check_no_self_qualified_accesses(Nbody6Dynamics) === nothing
    end

    @testset "Static QA — JET" begin
        using JET
        # Reports scoped to this package's own frames — Base/dependency
        # internals (e.g. @sync's sync_end, tuple broadcasting) produce
        # known false positives outside our control.
        JET.test_package(Nbody6Dynamics; target_modules = (Nbody6Dynamics,))
    end
end  # top-level testset
