# =============================================================================
# Core type definitions for Nbody6Setup
# =============================================================================

# ---------------------------------------------------------------------------
# Configuration types
# ---------------------------------------------------------------------------

"""
    InstallConfig

Configuration for the Nbody6++ source installation phase.
Controls git clone, reinstall behaviour, and clean builds.
"""
Base.@kwdef struct InstallConfig
    enabled::Bool             = true
    source_url::String        = "https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing.git"
    install_dir::String       = joinpath("backend", "Nbody6PPGPU-beijing")
    reinstall::Bool           = false
    clean_build::Bool         = true
end

"""
    BuildConfig

Configuration for the Nbody6++ compilation phase.
Controls configure flags, MPI/GPU/HDF5 toggles, CUDA path, and parallel make.
"""
Base.@kwdef struct BuildConfig
    configure_flags::Vector{String} = ["--enable-mcmodel=large", "--with-par=b1m"]
    enable_mpi::Bool                = false
    enable_hdf5::Bool               = true
    enable_gpu::Bool                = false
    cuda_path::String               = ""
    nproc::Int                      = 0
end

"""
    SimulationConfig

Configuration for the simulation execution phase.
Controls input file, run directory, MPI ranks, and run ID generation.
"""
Base.@kwdef struct SimulationConfig
    run_test::Bool        = true
    input_file::String    = "examples/input_files/N10k_noDat10.inp"
    runs_dir::String      = "runs"
    binary_name::String   = "nbody6++"
    mpi_ranks::Int        = 1
    run_id_prefix::String = "run"
end

"""
    PostprocessConfig

Configuration for the post-processing phase.
Controls which output files to read (snapshots, diagnostics, Lagrangian radii,
escapers, stellar evolution) and their file patterns.
"""
Base.@kwdef struct PostprocessConfig
    enabled::Bool              = true
    data_dir::String           = ""
    snapshot_format::String    = "conf3"
    snapshot_pattern::String   = "conf.3_*"
    hdf5_file::String          = "data.40.h5part"
    parse_stdout::Bool         = true
    stdout_file::String        = "out1000"
    read_lagr::Bool            = true
    lagr_file::String          = "lagr.7"
    read_escapers::Bool        = true
    escapers_file::String      = "esc.11"
    read_stellar_evo::Bool     = true
    stellar_evo_pattern::String = "sev.83_*"
end

"""
    VisualizationConfig

Configuration for the plotting and animation phase.
Controls output format (png/pdf/svg), DPI, figure size, and output directory.
"""
Base.@kwdef struct VisualizationConfig
    enabled::Bool          = true
    format::String         = "png"
    dpi::Int               = 300
    figsize::Tuple{Int,Int} = (8, 6)
    output_dir::String     = "plots"
end

"""
    MergerPipelineConfig

Configuration for the merger IC generation phase within the main pipeline.

# Fields
- `enabled::Bool`: whether to generate merger ICs before simulation
- `config_file::String`: path to the merger cluster TOML file
"""
Base.@kwdef struct MergerPipelineConfig
    enabled::Bool      = false
    config_file::String = ""
end

"""
    Nbody6Config

Top-level configuration aggregating all subsections.
"""
struct Nbody6Config
    install::InstallConfig
    build::BuildConfig
    simulation::SimulationConfig
    postprocess::PostprocessConfig
    visualization::VisualizationConfig
    merger::MergerPipelineConfig
end

# ---------------------------------------------------------------------------
# Snapshot header — maps to conf.3 AS(1:20) array
# ---------------------------------------------------------------------------

"""
    SnapshotHeader

Header metadata from a conf.3 snapshot file.  The `params` vector stores the
full AS(1:NK) Fortran array (NK = 20 by convention).
"""
struct SnapshotHeader
    ntot::Int32
    model::Int32
    nrun::Int32
    nk::Int32
    params::Vector{Float32}   # AS(1:NK)
end

"""    time_nb(h::SnapshotHeader) -> Float64
Simulation time in N-body units."""
time_nb(h::SnapshotHeader)  = Float64(h.params[1])

"""    npairs(h::SnapshotHeader) -> Int
Number of KS regularised binary pairs."""
npairs(h::SnapshotHeader)   = round(Int, h.params[2])

"""    rbar(h::SnapshotHeader) -> Float64
Length scaling factor: 1 NB length unit = `rbar` pc."""
rbar(h::SnapshotHeader)     = Float64(h.params[3])

"""    zmbar(h::SnapshotHeader) -> Float64
Mean stellar mass in solar masses."""
zmbar(h::SnapshotHeader)    = Float64(h.params[4])

rtide(h::SnapshotHeader)    = Float64(h.params[5])
tidal4(h::SnapshotHeader)   = Float64(h.params[6])
rdens(h::SnapshotHeader)    = Float64.(h.params[7:9])
time_tcr(h::SnapshotHeader) = Float64(h.params[10])

"""    tscale(h::SnapshotHeader) -> Float64
Time scaling factor: 1 NB time unit = `tscale` Myr."""
tscale(h::SnapshotHeader)   = Float64(h.params[11])

"""    vstar(h::SnapshotHeader) -> Float64
Velocity scaling factor: 1 NB velocity unit = `vstar` km/s."""
vstar(h::SnapshotHeader)    = Float64(h.params[12])

"""    rc(h::SnapshotHeader) -> Float64
Core radius in N-body units."""
rc(h::SnapshotHeader)       = Float64(h.params[13])

nc(h::SnapshotHeader)       = round(Int, h.params[14])
vc(h::SnapshotHeader)       = Float64(h.params[15])
rhom(h::SnapshotHeader)     = Float64(h.params[16])
cmax(h::SnapshotHeader)     = Float64(h.params[17])

"""    rscale(h::SnapshotHeader) -> Float64
Half-mass radius in N-body units."""
rscale(h::SnapshotHeader)   = Float64(h.params[18])

rsmin(h::SnapshotHeader)    = Float64(h.params[19])
dmin1(h::SnapshotHeader)    = Float64(h.params[20])

"""    time_myr(h::SnapshotHeader) -> Float64
Physical time in Myr, computed as `time_nb(h) * tscale(h)`."""
time_myr(h::SnapshotHeader) = time_nb(h) * tscale(h)

# ---------------------------------------------------------------------------
# Snapshot — particle data from conf.3 or HDF5
# ---------------------------------------------------------------------------

"""
    Snapshot

Particle data from a single snapshot (conf.3 or HDF5).
Positions and velocities are stored column-major: `pos[:, i]` gives particle i.
"""
struct Snapshot
    header::SnapshotHeader
    name::Vector{Int32}
    mass::Vector{Float32}
    pos::Matrix{Float32}      # 3 × N
    vel::Matrix{Float32}      # 3 × N
    rho::Vector{Float32}      # local density  (empty if unavailable)
    phi::Vector{Float32}      # potential       (empty if unavailable)
end

"""    nparticles(s::Snapshot) -> Int
Number of particles in the snapshot."""
nparticles(s::Snapshot) = length(s.name)

# ---------------------------------------------------------------------------
# Diagnostics parsed from simulation stdout
# ---------------------------------------------------------------------------

"""
    AdjustRecord

Single ADJUST output line from the simulation log.
"""
struct AdjustRecord
    time_nb::Float64
    time_myr::Float64
    qvir::Float64            # virial ratio 2T/|W|
    de_rel::Float64          # relative energy error
    e_tot::Float64           # total energy
    n::Int
    npairs::Int
    rscale::Float64          # half-mass radius (NB)
end

"""
    DiagnosticsData

Collection of parsed diagnostics from simulation stdout.
"""
struct DiagnosticsData
    adjust::Vector{AdjustRecord}
    physical_scaling::Dict{String,Float64}
end

# ---------------------------------------------------------------------------
# Lagrangian radii
# ---------------------------------------------------------------------------

"""
    LagrangianData

Lagrangian radii evolution parsed from lagr.7.
"""
struct LagrangianData
    time::Vector{Float64}
    mass_fractions::Vector{Float64}   # e.g. [0.001, 0.003, …, 1.0]
    radii::Matrix{Float64}           # n_fractions × n_times
end

# ---------------------------------------------------------------------------
# Escaper records from esc.11
# ---------------------------------------------------------------------------

"""
    EscaperRecord

A single escaper event from esc.11.
"""
struct EscaperRecord
    time_myr::Float64         # escape time [Myr]
    mass_solar::Float64       # mass [M☉]
    escape_energy::Float64    # dimensionless escape energy
    velocity_kms::Float64     # escape velocity [km/s]
    stellar_type::Int         # K* stellar type
    name::Int                 # particle identifier
end

# ---------------------------------------------------------------------------
# Stellar evolution records from sev*.83
# ---------------------------------------------------------------------------

"""
    StellarRecord

One star's properties from a single-star evolution snapshot (sev*.83).
"""
struct StellarRecord
    time_nb::Float64          # NB time
    index::Int32              # internal index
    name::Int32               # particle identifier
    stellar_type::Int32       # K* (0=MS, 1=HG, …, 13=BH)
    ri_rc::Float64            # distance / core radius
    mass_solar::Float64       # mass [M☉]
    log_luminosity::Float64   # log10(L/L☉)
    log_radius::Float64       # log10(R/R☉)
    log_teff::Float64         # log10(Teff/K)
end

"""
    StellarEvolutionSnapshot

All single-star data from one sev*.83 file.
"""
struct StellarEvolutionSnapshot
    time_myr::Float64
    n_stars::Int
    records::Vector{StellarRecord}
end

"""
    STELLAR_TYPE_LABELS

`Dict{Int,String}` mapping BSE stellar type codes (K*) to human-readable labels
used in HR diagram legends. Covers types 0 (MS) through 15 (Unknown).
"""
const STELLAR_TYPE_LABELS = Dict{Int,String}(
    0  => "MS (Main Seq.)",
    1  => "HG (Hertzsprung Gap)",
    2  => "GB (Giant Branch)",
    3  => "CHeB (Core He Burn.)",
    4  => "AGB (Asymp. Giant)",
    5  => "EAGB (Early AGB)",
    6  => "HeStar (He Star)",
    7  => "HeHG (He Hertzsp.)",
    8  => "HeGB (He Giant)",
    9  => "HeWD (He White Dwarf)",
    10 => "COWD (CO White Dwarf)",
    11 => "ONeWD (ONe White Dwarf)",
    12 => "NS (Neutron Star)",
    13 => "BH (Black Hole)",
    14 => "MSn (Naked He MS)",
    15 => "Unknown",
)

# ---------------------------------------------------------------------------
# N-body unit conversion helpers
# ---------------------------------------------------------------------------

"""
    UnitScaling

Physical unit conversion factors extracted from simulation output.
"""
struct UnitScaling
    rbar::Float64      # NB length  → pc
    zmbar::Float64     # average particle mass → M☉
    tscale::Float64    # NB time    → Myr
    vstar::Float64     # NB velocity → km/s
end

to_pc(u::UnitScaling, r_nb)    = r_nb * u.rbar
to_msun(u::UnitScaling, m_nb)  = m_nb * u.zmbar
to_myr(u::UnitScaling, t_nb)   = t_nb * u.tscale
to_kms(u::UnitScaling, v_nb)   = v_nb * u.vstar
