# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

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
    @test scan_empty.available[:binary_evo] == false
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

    # The engine's regularised-binary record (real fixture)
    cp(joinpath(@__DIR__, "fixtures", "bev.82_0"), joinpath(ext_dir, "bev.82_0"); force = true)

    # --- scan_output with partial data ---
    scan = scan_output(ext_dir)
    @test scan.available[:diagnostics] == true
    @test scan.available[:lagr] == true
    @test scan.available[:escapers] == true
    @test scan.available[:snapshots_conf3] == false
    @test scan.available[:stellar_evo] == false
    @test scan.available[:binary_evo] == true && length(scan.binary_evo_files) == 1

    # --- OutputScan display ---
    buf = IOBuffer()
    show(buf, MIME("text/plain"), scan)
    output_str = String(take!(buf))
    @test occursin("Diagnostics", output_str)
    @test occursin("out1000", output_str)
    @test occursin("not found", output_str)  # for missing categories
    @test occursin("Regularised binaries", output_str)
    @test occursin("binary_population", output_str)

    # --- scan_output error on non-existent directory ---
    @test_throws ErrorException scan_output("/nonexistent/path")

    # --- postprocess_external (data only, no plots) ---
    results = postprocess_external(ext_dir; make_plots = false)
    @test haskey(results, :scan)
    @test results[:scan] isa OutputScan
    @test haskey(results, :diagnostics)
    @test haskey(results, :lagr)
    @test haskey(results, :binary_evo)
    bevs_ext = results[:binary_evo]
    @test bevs_ext isa Vector{BinaryEvolutionSnapshot} && length(bevs_ext) == 1
    @test bevs_ext[1].n_pairs == length(bevs_ext[1].records) > 0
    # A bev.82 record without diagnostics is still read, and an absent one is absent
    rm(joinpath(ext_dir, "bev.82_0"))
    @test !haskey(postprocess_external(ext_dir; make_plots = false), :binary_evo)
end

# =====================================================================
# Adversarial external post-processing tests (included from separate file)
# =====================================================================
include("test_external_adversarial_inner.jl")

# =====================================================================
