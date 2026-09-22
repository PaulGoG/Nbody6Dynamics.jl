using Test
using CairoMakie   # triggers the figure extension; the smoke tests need it
using Nbody6Dynamics
# Public (documented, unexported) names the tests call unqualified.
using Nbody6Dynamics:
    nparticles,
    time_nb,
    time_myr,
    rbar,
    zmbar,
    tscale,
    vstar,
    rscale,
    rc,
    to_pc,
    to_msun,
    to_myr,
    to_kms,
    detect_platform,
    check_dependencies,
    detect_cuda_path,
    detect_compute_capabilities,
    cuda_arch_from_compute_cap,
    cuda_gencode_flags,
    resolve_cuda_arch,
    nvcc_release,
    nvcc_supported_archs,
    binary_hardness,
    hardness_scale,
    binary_scales,
    semi_major_axis_pc,
    binding_energy,
    stellar_class_index,
    class_counts,
    classes_present,
    run_sweep_point,
    read_sweep_index,
    write_sweep_index,
    profile_name,
    imf_name,
    expected_mass,
    sample_masses,
    kroupa_mean_mass,
    sample_plummer,
    sample_king,
    sample_kroupa,
    sample_binaries,
    expand_binaries,
    virialise!,
    kepler_velocity,
    jacobi_radius,
    crossing_time,
    to_nbody_units!,
    write_dat10,
    generate_merger_inp,
    resolve_nbody6_parameters,
    set_publication_theme!

# The figure routines and their helpers live in the package extension; tests
# that reach for an internal of the figure layer go through this module.
const MakieExt = Base.get_extension(Nbody6Dynamics, :Nbody6DynamicsMakieExt)
MakieExt === nothing && error("the Makie extension is not loaded; the figure tests cannot run")

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
    include("test_config.jl")
    include("test_platform_build.jl")
    include("test_validation.jl")
    include("test_run.jl")
    include("test_io.jl")
    include("test_diagnostics.jl")
    include("test_sweep_ensemble.jl")
    include("test_plotting.jl")
    include("test_external.jl")
    include("test_ic.jl")
    include("test_telemetry.jl")
    include("test_engine_gated.jl")
    include("test_qa.jl")
end  # top-level testset
