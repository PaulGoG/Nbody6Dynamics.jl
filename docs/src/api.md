# API Reference

## Pipeline

```@docs
run_pipeline
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
generate_run_id
export_for_paper
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
Nbody6Dynamics.PlotStyle
MergerPipelineConfig
```

## Data Types

```@docs
Snapshot
SnapshotHeader
DiagnosticsData
AdjustRecord
LagrangianData
UnitScaling
EscaperRecord
StellarRecord
StellarEvolutionSnapshot
STELLAR_TYPE_LABELS
```

## Snapshot Accessors

```@docs
nparticles
time_nb
time_myr
rbar
zmbar
tscale
vstar
rscale
rc
```

## I/O Readers

```@docs
read_conf3
read_all_conf3
read_diagnostics
extract_scaling
read_lagr
read_escapers
read_stellar_evolution
read_all_stellar_evolution
```

## Plotting

```@docs
set_publication_theme!
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

## Platform Utilities

```@docs
detect_platform
check_dependencies
detect_cuda_path
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
EqualMassIMF
imf_name
expected_mass
sample_masses
kroupa_mean_mass
```

### Configuration

```@docs
ClusterSpec
OrbitSpec
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
```

### Output

```@docs
write_dat10
generate_merger_inp
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
Nbody6Dynamics._validate
Nbody6Dynamics._bound_members
Nbody6Dynamics._kinetic_and_potential
Nbody6Dynamics._central_density_contrast
```
