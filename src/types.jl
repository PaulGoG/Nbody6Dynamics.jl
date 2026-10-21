# =============================================================================
# Core type definitions for Nbody6Dynamics
# =============================================================================

# ---------------------------------------------------------------------------
# Configuration types
# ---------------------------------------------------------------------------

"""
    InstallConfig

Configuration for the Nbody6++ source installation phase.
Controls git clone, the checked-out reference (`ref`: commit, tag, or
branch of the backend; the default is the validated commit), reinstall
behaviour, and clean builds.
"""
Base.@kwdef struct InstallConfig
    enabled::Bool = true
    source_url::String = "https://github.com/PaulGoG/Nbody6PPGPU-beijing.git"
    ref::String = "v2026.07.7"   # validated engine tag (the fork's main branch on upstream v2026.07); "" = clone default branch
    install_dir::String = joinpath("backend", "Nbody6PPGPU-beijing")
    reinstall::Bool = false
    clean_build::Bool = true
end

"""
    BuildConfig

Configuration for the Nbody6++ compilation phase.
Controls configure flags, MPI/GPU/HDF5 toggles, CUDA path, CUDA target
architectures, and parallel make.

# Fields
- `cuda_arch`: CUDA architectures the GPU kernels are compiled for
  (`"sm_90"` Hopper, `"sm_120"` consumer Blackwell, …): native code for
  each plus PTX for the highest. Empty = the compute capabilities of the
  devices visible to `nvidia-smi`, or the `nvcc` default target (PTX
  JIT-compiled at the first launch) when no device is visible
- `nvcc_flags`: extra `nvcc` options appended to the GPU build's `CUFLAGS`
  (e.g. `["-allow-unsupported-compiler"]` or `["-ccbin", "gcc-14"]` when
  the host compiler is newer than the toolkit supports)
"""
Base.@kwdef struct BuildConfig
    configure_flags::Vector{String} = ["--enable-mcmodel=large", "--with-par=b1m"]
    enable_mpi::Bool = false
    enable_hdf5::Bool = true
    enable_gpu::Bool = false
    cuda_path::String = ""
    cuda_arch::Vector{String} = String[]
    nvcc_flags::Vector{String} = String[]
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
- `gpu_list`: CUDA device indices the engine may use, exported as its
  `GPU_LIST` variable; empty = every visible device. At most four per
  process (the engine's `MAX_GPU`); requires `build.enable_gpu`
- `telemetry_interval`: sampling interval of the process-tree/GPU telemetry
  [s]; `0` disables the sampler (exact CPU accounting stays on)
- `startup_timeout`: wall-clock seconds within which the engine must report
  its first adjustment beyond t = 0; otherwise the run is terminated with an
  error. `0` disables
- `stall_timeout`: wall-clock seconds after the first adjustment beyond
  t = 0 during which neither the stdout capture grows nor a new restart dump
  appears; then the engine is asked to stop, terminated if it has not
  stopped the same time later, and the segment is recorded as `stalled`.
  `0` disables
- `exit_grace`: after the engine has printed `END RUN`, the wall-clock
  seconds its output directory may stay unchanged with the process still
  alive; then it is terminated and the segment is recorded as completed
  (the engine has been seen to finish its integration and never exit). A
  final COMMON dump in progress keeps the directory changing and is never
  interrupted. `0` disables
- `wall_budget`: wall-clock budget of one engine segment [s]; `0` = none.
  At `wall_budget − stop_margin` the engine is asked to stop (a `STOP` file
  in its working directory: it writes a restart dump and exits, and the
  segment is recorded as `stopped`); an engine still running at
  `wall_budget` is terminated
- `stop_margin`: wall-clock seconds the engine is given to act on a stop
  request, both before the wall budget and when the driver process itself
  is terminated (SIGTERM); must be smaller than a non-zero `wall_budget`
- `checkpoint_keep`: number of periodic restart dumps (`comm.2_<t>`) kept
  while the engine runs; older ones are deleted as new ones appear. `0`
  keeps every dump; otherwise at least 2, because the newest may be
  incomplete when the engine dies. Dumps written at the end of a run or at
  a stop request (`comm.1_<t>`) and the dump a restart started from are
  never deleted
- `live_diagnostics`: with `monitor`, print in-terminal sparklines of the
  virial ratio and energy error every `live_interval` seconds (opt-in;
  interactive terminals only, never in the log file)
"""
Base.@kwdef struct SimulationConfig
    run_test::Bool = true
    input_file::String = "input_files/engine/N1k_quick.inp"
    runs_dir::String = "runs"
    binary_name::String = "nbody6++"
    mpi_ranks::Int = 1
    omp_threads::Int = 0
    gpu_list::Vector{Int} = Int[]
    run_id_prefix::String = "run"
    monitor::Bool = false   # opt-in live progress ticker; interactive stderr only
    live_diagnostics::Bool = false
    live_interval::Float64 = 30.0
    telemetry_interval::Float64 = 5.0
    startup_timeout::Float64 = 0.0
    stall_timeout::Float64 = 0.0
    exit_grace::Float64 = 120.0
    wall_budget::Float64 = 0.0
    stop_margin::Float64 = 120.0
    checkpoint_keep::Int = 0
end

"""
    PostprocessConfig

Configuration for the post-processing phase.
Controls which output files to read (snapshots, diagnostics, Lagrangian radii,
escapers, stellar evolution, the counters of `events_file`) and their file
patterns, the snapshot integrals (`snapshot_integrals`), the limits of the
bound-set iteration (`bound_max_passes`, `bound_min_members`, see
[`bound_set`](@ref)), the largest snapshot on which the O(N²) diagnostics run
(`pair_sum_max_n`), the energy-jump localisation (`energy_events`, flagging an
adjustment interval whose |DE| reaches `max(energy_jump_threshold,
energy_jump_sigma × 1.4826 × MAD)`, see [`energy_intervals`](@ref)) and the
mass budget (`mass_budget`, see [`mass_budget`](@ref)).

The estimator settings follow. The centre keys (`centre_n_stop`,
`centre_stop_fraction`, `centre_n_floor`, `centre_shrink`) set the stop rule
and step of `shrinking_sphere_centre`. `lagr_n_min` is the effective
particle number `lagrangian_grid` requires on either side of a
Lagrangian radius. The profile keys (`profile_nbins`, `profile_n_min_sigma`,
`profile_n_min_beta`, `profile_n_min_density`) are the shell count and
membership thresholds of [`radial_profile`](@ref). `density_neighbours` is
the neighbour rank of the Casertano–Hut density estimate behind the core
radius and core density. The segregation keys (`msr_n_tracer`, `msr_scan`,
`msr_n_random`, `msr_n_random_small`, `msr_n_small`, `msr_seed`,
`msr_threshold`) fix the tracer class, the N_MST scan, the random comparison
sets and seed of [`mass_segregation`](@ref) and the Λ_MSR level of the
segregation time. The rotation keys (`rotation_n_min`, `rotation_nbins_max`,
`rotation_n_null`) are the shell membership, shell cap and null draws of
[`rotation_analysis`](@ref). The shape keys (`shape_n_min_shell`,
`shape_n_min_enclosed`) are the membership thresholds of `shape_tensor`.
The pair keys (`pair_k`, `pair_tidal_factor`) are the neighbour count and
tidal criterion of `snapshot_pairs`, which feeds the hardness series of
[`remnant_diagnostics`](@ref). `relaxation_gamma` holds the
Coulomb-logarithm arguments γ of `relaxation_times`.

The validity keys set the per-run validity record of `validity_record`:
`validity` switches the output of `validity.toml` and
`validity_snapshots.csv`; `validity_interval_de` and
`validity_cumulative_de` are the per-interval |DE| and the cumulative
|DETOT/E₀| above which the energy check flags the run;
`validity_mass_tolerance` is the relative mass-budget residual above which
an epoch is flagged; `validity_settle_tcr`, `validity_pericentre_tcr` and
`validity_window_tcr` are the violent-relaxation window after the start,
the half-width of the window about a pericentre and the window of the
stationarity tests, in crossing times; `validity_virial_tolerance` bounds
|2T/|W| − 1| of a virialised bound set, `validity_rate_tolerance` the change
per window of a stationary series, and `validity_bias_tolerance` the
relative calibrated bias above which an estimator is masked.
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
    read_binary_evo::Bool = true
    binary_evo_pattern::String = "bev.82_*"
    events_file::String = "event.35"
    pair_sum_max_n::Int = 2_000_000
    snapshot_integrals::Bool = true
    bound_max_passes::Int = 100
    bound_min_members::Int = 10
    energy_events::Bool = true
    energy_jump_threshold::Float64 = 1.0e-4
    energy_jump_sigma::Float64 = 5.0
    mass_budget::Bool = true
    centre_n_stop::Int = 1000
    centre_stop_fraction::Float64 = 0.01
    centre_n_floor::Int = 500
    centre_shrink::Float64 = 0.025
    lagr_n_min::Int = 100
    profile_nbins::Int = 12
    profile_n_min_sigma::Int = 50
    profile_n_min_beta::Int = 300
    profile_n_min_density::Int = 100
    density_neighbours::Int = 6
    msr_n_tracer::Int = 50
    msr_scan::Vector{Int} = [10, 20, 50, 100, 200]
    msr_n_random::Int = 200
    msr_n_random_small::Int = 1000
    msr_n_small::Int = 20
    msr_seed::Int = 0
    msr_threshold::Float64 = 2.0
    rotation_n_min::Int = 1000
    rotation_nbins_max::Int = 16
    rotation_n_null::Int = 20
    shape_n_min_shell::Int = 3000
    shape_n_min_enclosed::Int = 7000
    pair_k::Int = 10
    pair_tidal_factor::Float64 = 1 / 3
    relaxation_gamma::Vector{Float64} = [0.02, 0.11]
    angular_momentum::Bool = true
    mass_class_edges::Vector{Float64} = [0.5, 5.0]
    mass_class_n_min::Int = 100
    equipartition::Bool = true
    equipartition_bin_width::Float64 = 0.1
    equipartition_n_min::Int = 50
    merger_tree::Bool = true
    binary_spans::Bool = false
    span_n_periods::Float64 = 2.0
    validity::Bool = true
    validity_interval_de::Float64 = 1.0e-4
    validity_cumulative_de::Float64 = 1.0e-3
    validity_mass_tolerance::Float64 = 1.0e-6
    validity_settle_tcr::Float64 = 5.0
    validity_pericentre_tcr::Float64 = 2.0
    validity_window_tcr::Float64 = 5.0
    validity_virial_tolerance::Float64 = 0.05
    validity_rate_tolerance::Float64 = 0.05
    validity_bias_tolerance::Float64 = 0.05
end

"""
    PlotStyle

Stylistic plotting parameters, configurable via `[visualization.style]` in
`config.toml`. These are presentation knobs only — data-validity cutoffs
(e.g. BSE placeholder filtering) remain named constants in the plot code.

# Fields
- `marker_budget`: scatter marker size is `clamp(marker_budget/N, marker_min, marker_max)`
- `marker_min`, `marker_max`: clamp bounds for the scatter marker size [layout units]
- `q_log_threshold`: switch virial-ratio axes to log scale when max(Q) exceeds this
- `q_floor`: clamp floor for the virial ratio on *log-scale* axes only
- `zoom_frac`: extent-ratio threshold below which snapshot panels and animation
  frames are zoomed each on its own data
- `anim_fps`: animation frame rate; `0` selects automatically from frame count
- `anim_target_seconds`: target duration used by the automatic FPS selection
- `anim_px_per_unit`: pixels per canvas unit of animation frames (a single-panel
  canvas is 900 × 600 units)
- `snapshot_render`: particle projections as a mass-coloured scatter (`"scatter"`),
  as the projected surface mass density (`"density"`), or `"auto"`: density from
  `density_min_n` particles on, where scatter markers overlap into a flat disc
- `density_min_n`: particle count from which `"auto"` renders surface density
- `density_bins`: grid cells per axis of a surface-density panel
- `density_decades`: colour range of a surface-density map below its maximum [dex]
- `density_mass_frac`: mass fraction whose projected extent sets the limits of a
  surface-density panel, so that escapers do not set the scale
"""
Base.@kwdef struct PlotStyle
    marker_budget::Float64 = 27000.0
    marker_min::Float64 = 6.0
    marker_max::Float64 = 30.0
    q_log_threshold::Float64 = 10.0
    q_floor::Float64 = 1e-3
    zoom_frac::Float64 = 0.15
    anim_fps::Int = 0
    anim_target_seconds::Float64 = 12.0
    anim_px_per_unit::Float64 = 1.4
    snapshot_render::String = "auto"
    density_min_n::Int = 20_000
    density_bins::Int = 400
    density_decades::Float64 = 4.0
    density_mass_frac::Float64 = 0.99
end

"""
    VisualizationConfig

Configuration for the plotting and animation phase: output format
(png/pdf/svg), the printed width of an export, output directory, and the
[`PlotStyle`](@ref) presentation knobs.

Figures are composed in layout units on canvases fixed by the figure type (a
single panel is 900 × 600, stacked panels and auxiliary strips add height,
panel grids are wider) and exported so that the canvas width becomes
`export_width` inches on the page; text and line weights scale with it. A
journal target is an explicit override: `column = "single"` or `"double"`
exports at that column width instead. Either way the export enters a document
at native size.
"""
Base.@kwdef struct VisualizationConfig
    enabled::Bool = true
    format::String = "pdf"    # vector default; "png"/"svg" available
    dpi::Int = 300      # minimum raster resolution at the printed width
    export_width::Float64 = 6.5   # printed width of an export [in]
    column::String = ""           # "" | "single" | "double": journal column width instead of export_width
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
    ClusterConfig

Batch-scheduler layer of a run (`[cluster]`): one queue job per engine
segment, chained until the run reaches its end time. The job script is
generated into `<run_dir>/cluster/` at submission; nothing here is read
unless `enabled`.

# Fields
- `enabled`: submit runs as scheduler jobs (`submit_run`) instead of
  running them in the calling process
- `scheduler`: `"slurm"` (the only scheduler supported)
- `partition`, `account`, `qos`: the corresponding `sbatch` options; empty
  = not passed
- `job_name`: `--job-name`; empty = the run id
- `time_limit`: `--time`, in Slurm's time format (`minutes`,
  `minutes:seconds`, `hours:minutes:seconds`, `days-hours`,
  `days-hours:minutes`, `days-hours:minutes:seconds`)
- `partition_limit`: the partition's maximum time limit, same format; when
  nonempty, `time_limit` must not exceed it
- `nodes`, `ntasks`, `cpus_per_task`: the corresponding `sbatch` options
  (each ≥ 1)
- `mem`: `--mem`, e.g. `"64G"`; empty = the scheduler's default
- `gres`: `--gres`, e.g. `"gpu:1"`; empty = none
- `exclusive`: pass `--exclusive`
- `sbatch_options`: further `#SBATCH` options, verbatim without the
  `#SBATCH` prefix (e.g. `["--hint=nomultithread"]`)
- `prologue`: shell lines run in the job before Julia starts (module
  loads, depot paths)
- `julia`: the Julia executable the job runs; empty = the one that submits
- `julia_threads`: `JULIA_NUM_THREADS` of the job; `0` = `cpus_per_task`
- `chain`: `"resubmit"` (each job submits its successor while the run is
  short of its end time) or `"dependency"` (`max_segments` jobs submitted
  at once, each `--dependency=afterok` on the one before)
- `max_segments`: number of segments the chain may run; `0` = until the
  run completes (`resubmit` only); `dependency` requires ≥ 1
- `overhead`: wall-clock seconds of the job outside the engine's wall
  budget (start-up, initial-condition generation, post-processing);
  `wall_budget + stop_margin + overhead` must fit `time_limit`
- `sbatch`: the submission executable
"""
Base.@kwdef struct ClusterConfig
    enabled::Bool = false
    scheduler::String = "slurm"
    partition::String = ""
    account::String = ""
    qos::String = ""
    job_name::String = ""
    time_limit::String = "36:00:00"
    partition_limit::String = ""
    nodes::Int = 1
    ntasks::Int = 1
    cpus_per_task::Int = 8
    mem::String = ""
    gres::String = ""
    exclusive::Bool = false
    sbatch_options::Vector{String} = String[]
    prologue::Vector{String} = String[]
    julia::String = ""
    julia_threads::Int = 0
    chain::String = "resubmit"
    max_segments::Int = 0
    overhead::Float64 = 1800.0
    sbatch::String = "sbatch"
end

"""
    CostConfig

Cost model of a run (`[cost]`): the wall time, the number of segments and
the disk of an engine launch are estimated from a table of seconds per
N-body time unit against particle number (`estimate_cost`), logged, written
into the segment record, and compared with a budget.

# Fields
- `enabled`: estimate at every launch
- `table`: rows `[N, s_per_nb]` (≥ 2 rows, strictly increasing `N`,
  positive rates) replacing the built-in table measured on one H200 NVL
  device; empty = the built-in table
- `merger_factor`: factor applied to the rate of a two-cluster (merger)
  configuration, whose time unit is set by the whole configuration; > 0
- `budget_hours`: wall-time budget of one run [h]; ≥ 0, `0` = no check
- `on_exceed`: `"warn"` or `"refuse"` when the estimate exceeds the budget
"""
Base.@kwdef struct CostConfig
    enabled::Bool = true
    table::Vector{Vector{Float64}} = Vector{Float64}[]
    merger_factor::Float64 = 2.0
    budget_hours::Float64 = 0.0
    on_exceed::String = "warn"
end

"""
    Nbody6Config

Top-level configuration aggregating all subsections.

`config_dir` is the directory every relative path of the configuration
resolves against — `install.install_dir`, `simulation.input_file`,
`simulation.runs_dir`, `postprocess.data_dir` and `merger.config_file` — and
the default `base_dir` of the pipeline entry points. [`load_config`](@ref)
sets it to the directory of the file it read; a configuration built in memory
takes the working directory. It is not a configuration key and is not
written by [`save_config`](@ref).

The positional constructor takes the eight sections and `config_dir`; the
`cluster` and `cost` sections and `config_dir` may be omitted (defaults and
the working directory).
"""
struct Nbody6Config
    install::InstallConfig
    build::BuildConfig
    simulation::SimulationConfig
    postprocess::PostprocessConfig
    visualization::VisualizationConfig
    merger::MergerPipelineConfig
    cluster::ClusterConfig
    cost::CostConfig
    config_dir::String
end

Nbody6Config(
    install::InstallConfig,
    build::BuildConfig,
    simulation::SimulationConfig,
    postprocess::PostprocessConfig,
    visualization::VisualizationConfig,
    merger::MergerPipelineConfig,
    cluster::ClusterConfig = ClusterConfig(),
    cost::CostConfig = CostConfig(),
) = Nbody6Config(
    install,
    build,
    simulation,
    postprocess,
    visualization,
    merger,
    cluster,
    cost,
    pwd(),
)

Nbody6Config(
    install::InstallConfig,
    build::BuildConfig,
    simulation::SimulationConfig,
    postprocess::PostprocessConfig,
    visualization::VisualizationConfig,
    merger::MergerPipelineConfig,
    config_dir::AbstractString,
) = Nbody6Config(
    install,
    build,
    simulation,
    postprocess,
    visualization,
    merger,
    ClusterConfig(),
    CostConfig(),
    String(config_dir),
)

# ---------------------------------------------------------------------------
# Snapshot header — maps to conf.3 AS(1:20) array
# ---------------------------------------------------------------------------

"""
    SnapshotHeader

Header metadata from a conf.3 snapshot file. `ntot` is the engine's NTOT, the
number of entries in the file: the N stars followed by one centre-of-mass
entry per regularised pair (NTOT = N + NPAIRS). The `params` vector stores the
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
    Snapshot(header, name, mass, pos, vel, rho, phi[, npairs = 0, pair_name = Int32[]])

Stars of a single conf.3 snapshot. Positions and velocities are stored
column-major: `pos[:, i]` gives star i. The engine's centre-of-mass entries of
regularised pairs are not part of the snapshot; the components of the
`npairs` regularised pairs are the first `2npairs` stars, pair k being stars
`2k − 1` and `2k` (see [`pair_indices`](@ref)). `pair_name[k]` is the name of
pair k's centre-of-mass entry: negative for a hierarchical merger, whose
stars 2k − 1 and 2k are an inner component and the outer body rather than a
pair; empty when unknown. `rho` and `phi` are empty when the file does not
carry them.

Throws `DimensionMismatch` when the arrays disagree in length and
`ArgumentError` when `2npairs` exceeds the number of stars.
"""
struct Snapshot
    header::SnapshotHeader
    name::Vector{Int32}
    mass::Vector{Float32}
    pos::Matrix{Float32}      # 3 × N
    vel::Matrix{Float32}      # 3 × N
    rho::Vector{Float32}      # local density  (empty if unavailable)
    phi::Vector{Float32}      # potential       (empty if unavailable)
    npairs::Int               # regularised pairs; components are stars 1:2npairs
    pair_name::Vector{Int32}  # names of the pairs' centre-of-mass entries (empty if unknown)

    function Snapshot(
        header,
        name,
        mass,
        pos,
        vel,
        rho,
        phi,
        npairs::Integer = 0,
        pair_name = Int32[],
    )
        n = length(name)
        length(mass) == n ||
            throw(DimensionMismatch("mass has $(length(mass)) entries, name has $n"))
        size(pos) == (3, n) || throw(DimensionMismatch("pos is $(size(pos)), expected (3, $n)"))
        size(vel) == (3, n) || throw(DimensionMismatch("vel is $(size(vel)), expected (3, $n)"))
        for (label, a) in (("rho", rho), ("phi", phi))
            isempty(a) ||
                length(a) == n ||
                throw(DimensionMismatch("$label has $(length(a)) entries, name has $n"))
        end
        0 ≤ 2npairs ≤ n || throw(ArgumentError("npairs = $npairs is inconsistent with $n stars"))
        isempty(pair_name) ||
            length(pair_name) == npairs ||
            throw(DimensionMismatch("pair_name has $(length(pair_name)) entries for $npairs pairs"))
        return new(header, name, mass, pos, vel, rho, phi, Int(npairs), pair_name)
    end
end

"""    nparticles(s::Snapshot) -> Int
Number of stars in the snapshot (regularised-pair components counted
individually, centre-of-mass entries excluded)."""
nparticles(s::Snapshot) = length(s.name)

"""
    pair_indices(s::Snapshot) -> Vector{NTuple{2,Int}}

Star indices `(2k − 1, 2k)` of the components of each regularised (KS) pair of
the snapshot, k = 1, …, `s.npairs`.
"""
pair_indices(s::Snapshot) = [(2k - 1, 2k) for k in 1:s.npairs]

# ---------------------------------------------------------------------------
# Diagnostics parsed from simulation stdout
# ---------------------------------------------------------------------------

"""
    AdjustRecord

Single ADJUST output line from the simulation log, recording the relative energy error
`de_rel` (ΔE / max(T, |E|)), the interval energy change `de_abs` (engine DELTA, NB),
the engine running sum of ΔE `detot` (DETOT, NB), and the energy terms the line
prints after them (N-body units; NaN when a line lacks them): `ekin` (EKIN, the
kinetic energy of the centre-of-mass bodies), `pot` (POT, the positive potential
energy, W = −POT), `etide` (ETIDE), and the engine's correction channels `ebin`
(EBIN), `emerge` (EMERGE), `esub` (ESUB), `ecoll` (ECOLL), `emdot` (EMDOT),
`ecdot` (ECDOT), whose increments between consecutive lines locate the energy
changes the engine has accounted for. The ten-argument constructor leaves the
energy terms NaN. `nmerge` is the NMERGE count of the `TIME[NB]` line (stable
hierarchies held in the engine's merger table), 0 when not printed; the ten- and
nineteen-argument constructors set it to 0.
"""
struct AdjustRecord
    time_nb::Float64
    time_myr::Float64
    qvir::Float64            # virial ratio Q = T/|W| (equilibrium at 0.5)
    de_rel::Float64          # relative energy error of the adjustment interval, ΔE / max(T, |E|)
    e_tot::Float64           # total energy
    n::Int
    npairs::Int
    rscale::Float64          # half-mass radius (NB)
    de_abs::Float64          # energy change of the adjustment interval ΔE (engine DELTA, NB); NaN when not printed
    detot::Float64           # engine running sum of ΔE (DETOT, NB); NaN when not printed
    ekin::Float64            # EKIN (NB); NaN when not printed
    pot::Float64             # POT, positive (NB); NaN when not printed
    etide::Float64           # ETIDE (NB); NaN when not printed
    ebin::Float64            # EBIN (NB); NaN when not printed
    emerge::Float64          # EMERGE (NB); NaN when not printed
    esub::Float64            # ESUB (NB); NaN when not printed
    ecoll::Float64           # ECOLL (NB); NaN when not printed
    emdot::Float64           # EMDOT (NB); NaN when not printed
    ecdot::Float64           # ECDOT (NB); NaN when not printed
    nmerge::Int              # NMERGE of the TIME[NB] line (stable hierarchies in the merger table); 0 when not printed
end

function AdjustRecord(
    time_nb::Real,
    time_myr::Real,
    qvir::Real,
    de_rel::Real,
    e_tot::Real,
    n::Integer,
    npairs::Integer,
    rscale::Real,
    de_abs::Real,
    detot::Real,
)
    return AdjustRecord(
        time_nb,
        time_myr,
        qvir,
        de_rel,
        e_tot,
        n,
        npairs,
        rscale,
        de_abs,
        detot,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        NaN,
        0,
    )
end

function AdjustRecord(
    time_nb::Real,
    time_myr::Real,
    qvir::Real,
    de_rel::Real,
    e_tot::Real,
    n::Integer,
    npairs::Integer,
    rscale::Real,
    de_abs::Real,
    detot::Real,
    ekin::Real,
    pot::Real,
    etide::Real,
    ebin::Real,
    emerge::Real,
    esub::Real,
    ecoll::Real,
    emdot::Real,
    ecdot::Real,
)
    return AdjustRecord(
        time_nb,
        time_myr,
        qvir,
        de_rel,
        e_tot,
        n,
        npairs,
        rscale,
        de_abs,
        detot,
        ekin,
        pot,
        etide,
        ebin,
        emerge,
        esub,
        ecoll,
        emdot,
        ecdot,
        0,
    )
end

"""
    DiagnosticsData

Collection of parsed diagnostics from simulation stdout: the ADJUST records,
the `PHYSICAL SCALING` values, the engine options `kz` (the 50 values of the
`OPTIONS` echo; empty when the capture lacks it) and the linearised tidal
field `tidal` (`TIDAL(1:4)` of the `TIDAL PARAMETERS` echo, N-body units;
zeros when absent, i.e. no external field of `KZ(14) ≤ 2`). The two-argument
constructor leaves `kz` empty and `tidal` zero. `nbin0` is NBIN0 of the
`&INDATA` echo (`NBIN0 = …` on the line that starts with `ALPHAS`), `nzero`
N of the `&ININPUT` echo (the value row under `N  NFIX  NCRIT`), each `-1`
when absent; the two-, four- and five-argument constructors set the missing
ones to `-1`.
"""
struct DiagnosticsData
    adjust::Vector{AdjustRecord}
    physical_scaling::Dict{String,Float64}
    kz::Vector{Int}
    tidal::NTuple{4,Float64}
    nbin0::Int
    nzero::Int
end

function DiagnosticsData(adjust::Vector{AdjustRecord}, physical_scaling::Dict{String,Float64})
    return DiagnosticsData(adjust, physical_scaling, Int[], (0.0, 0.0, 0.0, 0.0), -1, -1)
end

function DiagnosticsData(
    adjust::Vector{AdjustRecord},
    physical_scaling::Dict{String,Float64},
    kz::Vector{Int},
    tidal::NTuple{4,Float64},
    nbin0::Int,
)
    return DiagnosticsData(adjust, physical_scaling, kz, tidal, nbin0, -1)
end

function DiagnosticsData(
    adjust::Vector{AdjustRecord},
    physical_scaling::Dict{String,Float64},
    kz::Vector{Int},
    tidal::NTuple{4,Float64},
)
    return DiagnosticsData(adjust, physical_scaling, kz, tidal, -1)
end

# ---------------------------------------------------------------------------
# Stamped event lines of the engine log and the counters of event.35
# ---------------------------------------------------------------------------

"""Stellar type value of an [`EngineEvent`](@ref) field the line does not carry."""
const KSTAR_ABSENT = -99

"""
    ENGINE_EVENT_KINDS

The kinds of [`EngineEvent`](@ref), in the order of the engine writers:
`:ks_rectification` (`BEG KSRECT`), `:ks_new` (`NEW KSREG`), `:ks_new_compact`
(`NS/BH BINARY`), `:ks_new_mdot` (`MDNEW KSREG`), `:ks_end_mdot` (`MDEND KSREG`),
`:remnant_formation` (`WD/NS/BH FORMATION`), `:kick_full`, `:kick_zero`,
`:kick_low` (`… VELOCITY KICK`), `:collision` (`… COLLISION`), `:merged_star`
(`NEW STAR`), `:coalescence` (`…COAL:`), `:escape` (`ESCAPE T,NESC`),
`:supernova_no_remnant` (`MDOT WARNING! SN KW=15`), `:chain_start` (`NEW CHAIN`),
`:chain_end` (`CHTERM:`).
"""
const ENGINE_EVENT_KINDS = (
    :ks_rectification,
    :ks_new,
    :ks_new_compact,
    :ks_new_mdot,
    :ks_end_mdot,
    :remnant_formation,
    :kick_full,
    :kick_zero,
    :kick_low,
    :collision,
    :merged_star,
    :coalescence,
    :escape,
    :supernova_no_remnant,
    :chain_start,
    :chain_end,
)

"""
    EngineEvent

One stamped event line of the engine log, from [`read_engine_events`](@ref).
`interval` is the index of the ADJUST record whose line follows the event
in the log — the adjustment interval the event falls in, `(t_{k−1}, t_k]`,
assigned by position because the engine prints every line of an interval
between the two ADJUST lines that bound it — and `line` the line number in
the capture. `time_nb` is the stamp of the line (NaN for the lines without
one: `NEW CHAIN`). The remaining fields hold what the kind carries, with
`0`, `KSTAR_ABSENT` and `NaN` where it does not:

| kind | `name1,2,3` | `kstar1,2,3` | `count` | `n_total` | `energy_nb` | `mass_solar` | `dmass_solar` | `velocity_kms` |
|---|---|---|---|---|---|---|---|---|
| `:ks_rectification`, `:ks_new`, `:ks_new_compact`, `:ks_new_mdot`, `:ks_end_mdot` | NM1, NM2, NMS | KW1, KW2, KWS | — | — | E_b | M1 + M2 | — | — |
| `:remnant_formation` | NAME | K* before, KW after | — | — | — | M1 after | M0 − M1 | — |
| `:kick_*` | NAME | K*0, K*, K*(ICM) | — | — | — | MN after | M0 − MN | VK |
| `:collision` | NM1, NM2, NMS | KW1, KW2, KWS | NSYS | — | E_b | M3 (product) | — | — |
| `:merged_star` | NM1, NM2 | K*(I1), K*(I2), K*(INEW) | — | — | — | M(INEW) | DM lost | — |
| `:coalescence` | N1, N2 | KW1, KW2, S | IQCOLL | — | E_b | M12S (product) | DM lost | — |
| `:escape` | — | — | NESC | N after | total E (BE(3)) | M* total after | — | — |
| `:supernova_no_remnant` | NM | K* before, 15 | — | — | — | M1 | — | — |
| `:chain_start` | NAMEC(1:3) | — | NCH | — | ENERGY | — | — | — |
| `:chain_end` | NI3, NI4, NCM | — | — | — | — | — | — | — |

Masses in M☉, energies in N-body units, velocities in km s⁻¹.
"""
Base.@kwdef struct EngineEvent
    time_nb::Float64
    kind::Symbol
    interval::Int
    line::Int
    name1::Int = 0
    name2::Int = 0
    name3::Int = 0
    kstar1::Int = KSTAR_ABSENT
    kstar2::Int = KSTAR_ABSENT
    kstar3::Int = KSTAR_ABSENT
    count::Int = 0
    n_total::Int = 0
    energy_nb::Float64 = NaN
    mass_solar::Float64 = NaN
    dmass_solar::Float64 = NaN
    velocity_kms::Float64 = NaN
end

"""
    EventCounters

One row of `event.35`, the engine's counters at an output time (`events.f`):
`time_myr`; the interaction counters `ndiss`, `ntide`, `nsync`, `ncoll`,
`ncoal`, `ndd`, `ncirc`, `nroche`, `nro`, `nce`, `nhyp`, `nhypc`, `nkick`;
the energy terms `ebin`, `emerge`, `ecoll`, `emdot`, `ecdot`, `ekick`,
`esesc`, `ebesc`, `emesc`, `degrav`, `ebind` (N-body units) and `mmax`
(the most massive body, M☉); the stellar-evolution counters `nmdot`, `nrg`,
`nhe`, `nrs`, `nnh`, `nwd`, `nsn`, `nbh`, `nbs`; the mass lost by stellar
evolution per phase `zmrg`, `zmhe`, `zmrs`, `zmnh`, `zmwd`, `zmsn` and in
total `zmdot` (M☉, cumulative); `ntype`, the number of stars of stellar
type K* = 0 … 15.
"""
struct EventCounters
    time_myr::Float64
    ndiss::Int
    ntide::Int
    nsync::Int
    ncoll::Int
    ncoal::Int
    ndd::Int
    ncirc::Int
    nroche::Int
    nro::Int
    nce::Int
    nhyp::Int
    nhypc::Int
    nkick::Int
    ebin::Float64
    emerge::Float64
    ecoll::Float64
    emdot::Float64
    ecdot::Float64
    ekick::Float64
    esesc::Float64
    ebesc::Float64
    emesc::Float64
    degrav::Float64
    ebind::Float64
    mmax::Float64
    nmdot::Int
    nrg::Int
    nhe::Int
    nrs::Int
    nnh::Int
    nwd::Int
    nsn::Int
    nbh::Int
    nbs::Int
    zmrg::Float64
    zmhe::Float64
    zmrs::Float64
    zmnh::Float64
    zmwd::Float64
    zmsn::Float64
    zmdot::Float64
    ntype::NTuple{16,Int}
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
    BinaryRecord

One KS-regularised binary from a bev.82_* file (`hrplot.F`, upstream
v2026.07: 32 tokens per line). Only regularised pairs are listed, so counts
derived from these records are lower bounds on the bound-pair population;
wide pairs beyond the regularisation distance never appear.

Fields with suffix `1`/`2` refer to the two components: masses in M☉, radii
in R☉, ages, epochs and main-sequence lifetimes in Myr;
`log_period_days` is ``\\log_{10}(P/\\mathrm{d})`` and
`log_semi_major_axis_rsun` is ``\\log_{10}(a/R_\\odot)``.
"""
struct BinaryRecord
    time_nb::Float64
    index1::Int32
    index2::Int32
    name1::Int32
    name2::Int32
    stellar_type1::Int32
    stellar_type2::Int32
    stellar_type_cm::Int32
    ri::Float64                      # c.m. distance from the density centre [pc]
    eccentricity::Float64
    log_period_days::Float64
    log_semi_major_axis_rsun::Float64
    mass1::Float64
    mass2::Float64
    log_luminosity1::Float64
    log_luminosity2::Float64
    log_radius1::Float64
    log_radius2::Float64
    log_teff1::Float64
    log_teff2::Float64
    age1_myr::Float64
    age2_myr::Float64
    epoch1_myr::Float64
    epoch2_myr::Float64
    ms_lifetime1_myr::Float64
    ms_lifetime2_myr::Float64
    mass_core1::Float64
    mass_core2::Float64
    radius_core1::Float64
    radius_core2::Float64
    radius_envelope1::Float64
    radius_envelope2::Float64
end

"""
    BinaryEvolutionSnapshot

All regularised binaries from one bev.82_* file; `time_myr` is the header
TPHYS [Myr] and `n_pairs` the header NPAIRS.
"""
struct BinaryEvolutionSnapshot
    time_myr::Float64
    n_pairs::Int
    records::Vector{BinaryRecord}
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

"""
    to_pc(u::UnitScaling, r_nb) -> Float64

Length in pc from N-body units: `r_nb * u.rbar`.
"""
to_pc(u::UnitScaling, r_nb) = r_nb * u.rbar

"""
    to_msun(u::UnitScaling, m_nb) -> Float64

Mass in M☉ from N-body units: `m_nb * u.zmbar` (the total-mass scale factor).
"""
to_msun(u::UnitScaling, m_nb) = m_nb * u.zmbar

"""
    to_myr(u::UnitScaling, t_nb) -> Float64

Time in Myr from N-body units: `t_nb * u.tscale`.
"""
to_myr(u::UnitScaling, t_nb) = t_nb * u.tscale

"""
    to_kms(u::UnitScaling, v_nb) -> Float64

Velocity in km s⁻¹ from N-body units: `v_nb * u.vstar`.
"""
to_kms(u::UnitScaling, v_nb) = v_nb * u.vstar
