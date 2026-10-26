# API Reference

```@meta
CurrentModule = Nbody6Dynamics
```

Names listed here without a `Nbody6Dynamics.` prefix are either exported or `public`; the export list in `src/Nbody6Dynamics.jl` says which. A `public` name is documented and stable but is called as `Nbody6Dynamics.name`.

## Pipeline

```@docs
example_input
run_pipeline
resume_pipeline
run_gpu_validation
postprocess
generate_plots
load_config
save_config
```

## Simulation

```@docs
setup_nbody6
run_simulation
restart_simulation
resume_run
generate_run_id
export_for_paper
```

## Cluster Jobs

```@docs
submit_run
cluster_segment
job_script
read_job_ledger
```

## Cost Model

```@docs
CostEstimate
estimate_cost
cost_rate
```

## External Post-Processing

```@docs
scan_output
postprocess_external
OutputScan
```

## Configuration Types

```@docs
Nbody6Config
InstallConfig
BuildConfig
SimulationConfig
PostprocessConfig
VisualizationConfig
PlotStyle
MergerPipelineConfig
ClusterConfig
CostConfig
```

## Data Types

```@docs
Snapshot
SnapshotHeader
DiagnosticsData
AdjustRecord
EngineEvent
ENGINE_EVENT_KINDS
KSTAR_ABSENT
EventCounters
LagrangianData
UnitScaling
EscaperRecord
StellarRecord
StellarEvolutionSnapshot
STELLAR_TYPE_LABELS
BinaryRecord
BinaryEvolutionSnapshot
BinaryPopulation
```

## Snapshot Accessors

```@docs
nparticles
pair_indices
time_nb
time_myr
rbar
zmbar
tscale
vstar
rscale
rc
```

```@docs
Nbody6Dynamics.to_pc
Nbody6Dynamics.to_msun
Nbody6Dynamics.to_myr
Nbody6Dynamics.to_kms
```

## I/O Readers

```@docs
read_conf3
read_all_conf3
read_telemetry
read_run_telemetry
read_diagnostics
extract_scaling
cumulative_energy_error
read_energy_totals
read_engine_events
read_event_counters
read_lagr
read_escapers
read_stellar_evolution
read_all_stellar_evolution
read_binary_evolution
read_all_binary_evolution
```

### Binary, few-body and gravitational-wave output

```@docs
read_binary_events
BinaryEventLog
DYNAMICS_EVENT_KINDS
EVOLUTION_EVENT_KINDS
KSEvent
FewBodyEvent
HierarchyEvent
ExchangeEvent
TerminationEvent
GWEvent
BinaryEscapeEvent
MassLossEvent
TidalEvent
RocheEvent
read_ks_pair_listing
read_all_ks_pair_listings
KSPairListing
KSPairRecord
read_wide_pair_listing
read_all_wide_pair_listings
WidePairListing
WidePairRecord
read_hierarchy_listing
HierarchyListingRecord
read_hierarchy_log
HierarchyLogRecord
read_collision_log
CollisionRecord
read_coalescence_log
CoalescenceRecord
read_cmbody_log
CMBodyRecord
read_binary_stage_log
BinaryStageRecord
read_degenerate_log
DegenerateRecord
read_escaped_binary_log
EscapedBinaryRecord
read_hardest_binary_log
HardestBinaryRecord
```

## Parameter Sweeps

```@docs
SweepConfig
SweepPoint
load_sweep_config
sweep_points
prepare_sweep
run_sweep
run_sweep_point
read_sweep_index
write_sweep_index
sweep_summary
write_sweep_summary
sweep_visualization
plot_sweep_lagrangian
plot_sweep_energy
sweep_figures
control_merger_dict
write_control_merger_config
plot_control_comparison
EnsembleStatistics
ensemble_statistics
sweep_ensembles
plot_sweep_ensemble
```

## Potential, Bound Sets and Snapshot Integrals

```@docs
pair_potential
SystemSet
reduce_pairs
BoundSet
bound_set
SnapshotIntegrals
snapshot_integrals
write_snapshot_integrals
tidal_mode
engine_energy_valid
```

## Structure Estimators

```@docs
shrinking_sphere_centre
lagrangian_grid
relaxation_times
local_relaxation_time
```

## Energy-Jump Localisation and Mass Budget

```@docs
EnergyInterval
EnergyEventAnalysis
energy_intervals
energy_events_summary
write_energy_intervals
write_engine_events
MassBudgetRow
mass_budget
mass_budget_summary
write_mass_budget
EnergyBudgetRow
energy_budget
energy_budget_summary
write_energy_budget
```

## Remnant Diagnostics

```@docs
RemnantDiagnostics
RotationProfile
remnant_diagnostics
write_remnant_diagnostics
coalescence_time
core_radius
rotation_analysis
mass_segregation
core_collapse_markers
plot_remnant_rotation
plot_rotation_profile
plot_remnant_structure
plot_mass_segregation_evolution
remnant_figures
```

## Shape, Segregation and Mixing

```@docs
ShapeShell
shape_tensor
projected_ellipticity
aplus
two_sample_d
mixing
```

## Snapshot Binaries

```@docs
SnapshotPair
local_density
snapshot_pairs
reconcile_pairs
```

## Angular-Momentum Budget

```@docs
AngularMomentumBudget
angular_momentum_budget
write_angular_momentum
```

## Equipartition

```@docs
EquipartitionFit
equipartition_fit
EquipartitionRow
equipartition_series
write_equipartition
```

## Merger Tree and Progenitor Budget

```@docs
MergerEvent
MergerTree
merger_tree
write_merger_tree
ProgenitorBudgetRow
progenitor_budget
write_progenitor_budget
```

## Binary Spans

```@docs
BinarySpan
binary_spans
write_binary_spans
```

## Binary Census

```@docs
CensusSpan
binary_census
write_binary_census
```

## Gravitational-Wave Observables

```@docs
gw_enhancement
gw_time_factor
gw_inspiral_time
gw_inspiral_time_mandel
gw_zwick_factor
chirp_mass
orbital_frequency
gw_frequency
eccentricity_at_frequency
```

## Binary Population Series

```@docs
BinaryPopulationRow
binary_population_series
write_binary_population
stellar_type_map
BinaryPairTypeRow
binary_pair_type_series
write_binary_pair_types
BinaryDistributionRow
binary_distribution_series
write_binary_distributions
```

## Survival Statistics

```@docs
EventTable
event_table
SurvivalCurve
kaplan_meier
survival_at
confidence_interval
nelson_aalen
CumulativeIncidence
aalen_johansen
turnbull
logrank_test
cox_snell_residuals
write_survival_curve
write_cumulative_incidence
```

## Collision Catalogue

```@docs
CollisionRow
collision_catalogue
collision_summary
progenitor_of
write_collision_catalogue
```

## Gravitational-Wave Catalogues

```@docs
pair_class
peters_semi_major_axis
band_residence
eccentricity_after
CompactBinaryRow
compact_binary_catalogue
write_compact_binary_catalogue
EscapedBinaryRow
escaped_binary_catalogue
write_escaped_binary_catalogue
GWMergerRow
gw_merger_catalogue
write_gw_merger_catalogue
```

## Rotating King Model

```@docs
RotatingKingModel
RotatingKingProfile
solve_rotating_king
rotating_king_moments
sample_rotating_king
rotation_moments
```

## Bound Spans

```@docs
BoundPair
BoundSpan
bound_pairs
bound_spans
write_bound_spans
```

## Validity Record

```@docs
ValidityRecord
ValidityEpoch
ValiditySnapshot
validity_record
write_validity
write_validity_snapshots
CalibrationRow
CalibrationTable
CalibratedEstimate
load_calibration
calibrated
```

## Binary Population

```@docs
binary_population
binary_hardness
hardness_scale
binary_scales
binding_energy
semi_major_axis_pc
plot_binary_population
plot_binary_orbital_elements
plot_binary_period_distribution
plot_binary_timeline
plot_binary_population_series
plot_survival_curves
plot_gw_planes
```

## Stellar Population

Stellar classes, the population of an epoch on the HR plane (single stars
and members of KS-regularised pairs together), and the class census of a
run. `run_pipeline` writes the census as `stellar_census.csv`.

```@docs
StellarClass
STELLAR_CLASSES
stellar_class
stellar_class_index
HRPopulation
hr_population
hr_populations
StellarCensus
stellar_census
class_counts
classes_present
write_stellar_census
```

## Plotting

The figure routines below are implemented by the package extension
`Nbody6DynamicsMakieExt` and require a Makie backend in the session
(`using CairoMakie`); see [Visualisation](@ref "10. Visualisation").

```@docs
plotting_available
PlottingUnavailable
set_publication_theme!
publication_theme
plot_snapshot
plot_snapshot_evolution
plot_energy
plot_particle_count
plot_lagrangian
plot_hr
plot_hr_evolution
plot_escapers
plot_escape_anisotropy
plot_mass_segregation
plot_evolutionary_clock
plot_core_mass
```

## Animations

```@docs
animate_cluster
animate_lagrangian
animate_hr
```

## Run Telemetry

```@docs
TelemetrySample
plot_telemetry
```

## Platform Utilities

```@docs
detect_platform
check_dependencies
detect_cuda_path
detect_compute_capabilities
cuda_arch_from_compute_cap
cuda_gencode_flags
resolve_cuda_arch
nvcc_release
nvcc_supported_archs
```

## Static Test Models

```@docs
sample_plummer_om
lynden_bell_flip!
stretch_positions!
rank_masses_by_energy!
```

## Merger Initial Conditions

### Density Profiles

```@docs
DensityProfile
KingProfile
PlummerProfile
profile_name
```

### Initial Mass Functions

```@docs
IMFSpec
KroupaIMF
RescaledKroupaIMF
PowerLawIMF
EqualMassIMF
imf_name
expected_mass
sample_masses
kroupa_mean_mass
```

### Configuration

```@docs
ClusterSpec
SegregationSpec
BinarySpec
MassiveBinarySpec
OrbitSpec
ComplexSpec
MergerOutputSpec
Nbody6ParameterSpec
StellarSpec
TidalSpec
MergerConfig
MergerICResult
load_merger_config
```

### IC Generation

```@docs
run_merger_pipeline
generate_merger_ic
load_merger_ic_result
plot_merger_ic
sample_plummer
sample_king
model_density
sample_kroupa
virialise!
```

### Orbital Mechanics

```@docs
kepler_velocity
jacobi_radius
Nbody6Dynamics.roche_truncation_radius
Nbody6Dynamics.orbit_apocentre
Nbody6Dynamics.two_cluster_orbit_quantities
Nbody6Dynamics.tidal_frame_rotation
Nbody6Dynamics.pair_jacobi_radius
Nbody6Dynamics.to_rotating_frame!
Nbody6Dynamics.complex_member_masses
Nbody6Dynamics.complex_member_radius
Nbody6Dynamics.complex_centres
Nbody6Dynamics.complex_scale_velocities!
```

### Output

```@docs
write_dat10
generate_merger_inp
engine_interval
sample_binaries
Nbody6Dynamics.zams_radius
expand_binaries
crossing_time
resolve_nbody6_parameters
to_nbody_units!
```

### Merger Run Diagnostics

```@docs
parse_merger_summary
per_cluster_virial
bound_fraction
plot_cluster_structure
ClusterStructure
cluster_structure
RadialProfile
radial_profile
cluster_profiles
system_profile
plot_density_profiles
plot_velocity_dispersion
plot_cluster_separation
plot_cluster_virial
```

## Internals referenced from the public docstrings

```@docs
Nbody6Dynamics._engine_digit_counter_terminates
Nbody6Dynamics._validate
Nbody6Dynamics._kinetic_and_potential
Nbody6Dynamics._central_density_contrast
Nbody6Dynamics._sweep_point_outcome
Nbody6Dynamics._run_series
Nbody6Dynamics._control_time_factor
Nbody6Dynamics._local_density
```
