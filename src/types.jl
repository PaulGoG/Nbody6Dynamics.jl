# =============================================================================
# Core type definitions for Nbody6Dynamics
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
    enabled::Bool = true
    source_url::String = "https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing.git"
    install_dir::String = joinpath("backend", "Nbody6PPGPU-beijing")
    reinstall::Bool = false
    clean_build::Bool = true
end

"""
    BuildConfig

Configuration for the Nbody6++ compilation phase.
Controls configure flags, MPI/GPU/HDF5 toggles, CUDA path, and parallel make.
"""
Base.@kwdef struct BuildConfig
    configure_flags::Vector{String} = ["--enable-mcmodel=large", "--with-par=b1m"]
    enable_mpi::Bool = false
    enable_hdf5::Bool = true
    enable_gpu::Bool = false
    cuda_path::String = ""
    nproc::Int = 0
end

"""
    SimulationConfig

Configuration for the simulation execution phase.
Controls input file, run directory, MPI ranks, backend OpenMP thread count,
run ID generation, and runtime telemetry.

# Fields
- `omp_threads`: OpenMP threads for the backend; `0` leaves the runtime
  default (an inherited `OMP_NUM_THREADS`, else every logical CPU)
- `telemetry_interval`: sampling interval of the process-tree/GPU telemetry
  [s]; `0` disables the sampler (exact CPU accounting stays on)
"""
Base.@kwdef struct SimulationConfig
    run_test::Bool = true
    input_file::String = "examples/input_files/N10k_noDat10.inp"
    runs_dir::String = "runs"
    binary_name::String = "nbody6++"
    mpi_ranks::Int = 1
    omp_threads::Int = 0
    run_id_prefix::String = "run"
    monitor::Bool = false   # opt-in live progress ticker (§9); interactive stderr only
    telemetry_interval::Float64 = 5.0
end

"""
    PostprocessConfig

Configuration for the post-processing phase.
Controls which output files to read (snapshots, diagnostics, Lagrangian radii,
escapers, stellar evolution) and their file patterns.
"""
Base.@kwdef struct PostprocessConfig
    enabled::Bool = true
    data_dir::String = ""
    snapshot_format::String = "conf3"   # "hdf5" is no longer supported
    snapshot_pattern::String = "conf.3_*"
    parse_stdout::Bool = true
    stdout_file::String = "out1000"
    read_lagr::Bool = true
    lagr_file::String = "lagr.7"
    read_escapers::Bool = true
    escapers_file::String = "esc.11"
    read_stellar_evo::Bool = true
    stellar_evo_pattern::String = "sev.83_*"
end

"""
    PlotStyle

Stylistic plotting parameters, configurable via `[visualization.style]` in
`config.toml`. These are presentation knobs only — data-validity cutoffs
(e.g. BSE placeholder filtering) remain named constants in the plot code.

# Fields
- `marker_budget`: scatter marker size is `clamp(marker_budget/N, marker_min, marker_max)`
- `marker_min`, `marker_max`: clamp bounds for the scatter marker size [pt]
- `q_log_threshold`: switch virial-ratio axes to log scale when max(Q) exceeds this
- `q_floor`: clamp floor for the virial ratio on *log-scale* axes only
- `zoom_frac`: fraction of particles defining the adaptive zoom-in radius
- `anim_fps`: animation frame rate; `0` selects automatically from frame count
- `anim_target_seconds`: target duration used by the automatic FPS selection
"""
Base.@kwdef struct PlotStyle
    marker_budget::Float64 = 18000.0
    marker_min::Float64 = 4.0
    marker_max::Float64 = 20.0
    q_log_threshold::Float64 = 10.0
    q_floor::Float64 = 1e-3
    zoom_frac::Float64 = 0.15
    anim_fps::Int = 0
    anim_target_seconds::Float64 = 12.0
end

"""
    VisualizationConfig

Configuration for the plotting and animation phase.
Controls output format (png/pdf/svg), DPI, figure size, output directory,
and the [`PlotStyle`](@ref) presentation knobs.
"""
Base.@kwdef struct VisualizationConfig
    enabled::Bool = true
    format::String = "pdf"    # vector default; "png"/"svg" available
    dpi::Int = 300      # raster resolution at FINAL print size
    column::String = "single" # "single" | "double" journal-width preset;
    # "" falls back to free-form figsize
    figsize::Tuple{Float64,Float64} = (8.0, 6.0)  # inches; used only when column = ""
    units::String = "physical"      # "physical" | "nbody" axis units
    output_dir::String = "plots"
    style::PlotStyle = PlotStyle()
end

"""
    MergerPipelineConfig

Configuration for the merger IC generation phase within the main pipeline.

# Fields
- `enabled::Bool`: whether to generate merger ICs before simulation
- `config_file::String`: path to the merger cluster TOML file
"""
Base.@kwdef struct MergerPipelineConfig
    enabled::Bool = false
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
full AS(1:NK) Fortran array (NK = 20 by convention). Slot map (accessors
exist for the starred entries):

| AS  | Quantity                              | AS  | Quantity                     |
|:----|:--------------------------------------|:----|:-----------------------------|
| 1*  | TTOT, time [NB] (`time_nb`)           | 11* | TSCALE, NB time → Myr        |
| 2   | NPAIRS, KS pair count                 | 12* | VSTAR, NB velocity → km/s    |
| 3*  | RBAR, NB length → pc                  | 13* | RC, core radius [NB]         |
| 4*  | ZMBAR, NB mass → M☉ (total-mass)      | 14  | NC, core member count        |
| 5   | RTIDE, tidal radius [NB]              | 15  | VC, core velocity dispersion |
| 6   | TIDAL(4), tidal-field coefficient     | 16  | RHOM, mean core density      |
| 7–9 | RDENS(1:3), density centre [NB]       | 17  | CMAX                         |
| 10  | TTOT/TCR, time in crossing times      | 18* | RSCALE, half-mass radius [NB]|
|     |                                       | 19  | RSMIN                        |
|     |                                       | 20  | DMIN1                        |
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
time_nb(h::SnapshotHeader) = Float64(h.params[1])

"""    rbar(h::SnapshotHeader) -> Float64
Length scaling factor: 1 NB length unit = `rbar` pc."""
rbar(h::SnapshotHeader) = Float64(h.params[3])

"""    zmbar(h::SnapshotHeader) -> Float64
Mass scaling factor: 1 NB mass unit (the total cluster mass) = `zmbar` M☉.
NOT the mean stellar mass — Nbody6++ redefines ZMBAR as the total-mass
scale factor at startup (`start.F`); the mean mass is printed separately
as `<M>` in the PHYSICAL SCALING line."""
zmbar(h::SnapshotHeader) = Float64(h.params[4])

"""    tscale(h::SnapshotHeader) -> Float64
Time scaling factor: 1 NB time unit = `tscale` Myr."""
tscale(h::SnapshotHeader) = Float64(h.params[11])

"""    vstar(h::SnapshotHeader) -> Float64
Velocity scaling factor: 1 NB velocity unit = `vstar` km/s."""
vstar(h::SnapshotHeader) = Float64(h.params[12])

"""    rc(h::SnapshotHeader) -> Float64
Core radius in N-body units."""
rc(h::SnapshotHeader) = Float64(h.params[13])

"""    rscale(h::SnapshotHeader) -> Float64
Half-mass radius in N-body units."""
rscale(h::SnapshotHeader) = Float64(h.params[18])

"""    time_myr(h::SnapshotHeader) -> Float64
Physical time in Myr, computed as `time_nb(h) * tscale(h)`."""
time_myr(h::SnapshotHeader) = time_nb(h) * tscale(h)

# ---------------------------------------------------------------------------
# Snapshot — particle data from conf.3 or HDF5
# ---------------------------------------------------------------------------

"""
    Snapshot

Particle data from a single conf.3 snapshot.
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
    qvir::Float64            # virial ratio Q = T/|W| (equilibrium at 0.5)
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

A single escaper event from esc.11.  The direction angles follow the
`escape.F` convention: `phi_deg` is the azimuth of the escape direction
measured from the x-axis in [0°, 360°]; `theta_deg` is the elevation
from the xy-plane in [-90°, 90°].
"""
struct EscaperRecord
    time_myr::Float64         # escape time [Myr]
    mass_solar::Float64       # mass [M☉]
    escape_energy::Float64    # dimensionless escape energy
    velocity_kms::Float64     # escape velocity [km/s]
    stellar_type::Int         # K* stellar type
    name::Int                 # particle identifier
    phi_deg::Float64          # escape azimuth from x-axis [deg], [0, 360]
    theta_deg::Float64        # escape elevation from xy-plane [deg], [-90, 90]
end

# ---------------------------------------------------------------------------
# Stellar evolution records from sev.83_*
# ---------------------------------------------------------------------------

"""
    StellarRecord

One star's properties from a single-star evolution snapshot (sev.83_*),
matching the upstream v2026.07+ `hrplot.F` output.
"""
struct StellarRecord
    time_nb::Float64          # NB time (TTOT, first token of each data line)
    index::Int32              # internal index
    name::Int32               # particle identifier
    stellar_type::Int32       # K* (Hurley: 0/1=MS, 2=HG, …, 13=NS, 14=BH)
    ri::Float64               # RI, distance from density centre [pc]
    mass_solar::Float64       # mass [M☉]
    log_luminosity::Float64   # log10(L/L☉)
    log_radius::Float64       # log10(R/R☉)
    log_teff::Float64         # log10(Teff/K)
    ms_lifetime_myr::Float64  # TM, main-sequence lifetime [Myr]
    mass_core::Float64        # MC, core mass [M☉]
    radius_core::Float64      # RCC, core radius [R☉]
    radius_envelope::Float64  # RE, envelope radius [R☉]
end

# Convenience constructor for synthetic records (SSE fields default to NaN).
StellarRecord(
    time_nb,
    index,
    name,
    stellar_type,
    ri,
    mass_solar,
    log_luminosity,
    log_radius,
    log_teff,
) = StellarRecord(
    time_nb,
    index,
    name,
    stellar_type,
    ri,
    mass_solar,
    log_luminosity,
    log_radius,
    log_teff,
    NaN,
    NaN,
    NaN,
    NaN,
)

"""
    StellarEvolutionSnapshot

All single-star data from one sev.83_* file.
"""
struct StellarEvolutionSnapshot
    time_myr::Float64
    n_stars::Int
    records::Vector{StellarRecord}
end

"""
    STELLAR_TYPE_LABELS

`Dict{Int,String}` mapping SSE/BSE stellar type codes (K*) to human-readable
labels used in HR diagram legends. Follows the standard Hurley et al. (2000)
convention used by this fork (`global_output.F`): 0/1 = low-/high-mass MS,
2 = HG, …, 13 = NS, 14 = BH, 15 = massless supernova remnant.
"""
const STELLAR_TYPE_LABELS = Dict{Int,String}(
    0 => "MS (low-mass, M < 0.7)",
    1 => "MS (Main Seq.)",
    2 => "HG (Hertzsprung Gap)",
    3 => "GB (Giant Branch)",
    4 => "CHeB (Core He Burn.)",
    5 => "EAGB (Early AGB)",
    6 => "TPAGB (Therm. Puls. AGB)",
    7 => "HeMS (Naked He MS)",
    8 => "HeHG (He Hertzsp. Gap)",
    9 => "HeGB (He Giant Branch)",
    10 => "HeWD (He White Dwarf)",
    11 => "COWD (CO White Dwarf)",
    12 => "ONeWD (ONe White Dwarf)",
    13 => "NS (Neutron Star)",
    14 => "BH (Black Hole)",
    15 => "SNR (Massless Remnant)",
)

# ---------------------------------------------------------------------------
# N-body unit conversion helpers
# ---------------------------------------------------------------------------

"""
    UnitScaling

Physical unit conversion factors extracted from simulation output.
"""
struct UnitScaling
    rbar::Float64      # NB length → pc (R*)
    zmbar::Float64     # NB mass → M☉ (M*, the total-mass scale factor)
    tscale::Float64    # NB time → Myr (T*)
    vstar::Float64     # NB velocity → km/s (V*)
end

to_pc(u::UnitScaling, r_nb) = r_nb * u.rbar
to_msun(u::UnitScaling, m_nb) = m_nb * u.zmbar
to_myr(u::UnitScaling, t_nb) = t_nb * u.tscale
to_kms(u::UnitScaling, v_nb) = v_nb * u.vstar
