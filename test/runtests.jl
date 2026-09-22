using Test
using CairoMakie   # triggers the figure extension; the smoke tests need it
using Nbody6Dynamics

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
