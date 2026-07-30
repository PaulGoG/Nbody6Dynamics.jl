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
generate_run_id
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
read_hdf5_snapshot
read_hdf5_snapshots
list_hdf5_steps
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

### Configuration

```@docs
MergerPipelineConfig
ClusterSpec
OrbitSpec
MergerOutputSpec
MergerConfig
MergerICResult
load_merger_config
```

### IC Generation

```@docs
run_merger_pipeline
generate_merger_ic
plot_merger_ic
sample_plummer
sample_king
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
to_nbody_units!
```
