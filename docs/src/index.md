# Nbody6Setup.jl

A Julia package for automated setup, execution, post-processing, and visualisation of the [Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing) N-body astrophysical simulation code.

## Features

- **Install & Build** -- Clone, configure, patch (HDF5), and compile Nbody6++GPU with auto-detected CUDA/MPI
- **Simulate** -- Launch simulations with unique run IDs, frozen configs, and real-time ADJUST monitoring
- **Post-process** -- Read all Nbody6++ output formats (conf.3, H5Part, out1000, lagr.7, esc.11, sev*.83)
- **Visualise** -- Generate 18 publication-quality plots and GIF animations with CairoMakie

## Quick Start

```bash
cd Nbody6Setup

# Install Julia dependencies (one-time)
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Edit configuration
$EDITOR config.toml

# Run the full pipeline
julia --project=. -e 'using Nbody6Setup; run_pipeline(load_config("config.toml"))'
```

## Usage Modes

| Mode              | Settings                                                      |
|:------------------|:--------------------------------------------------------------|
| Full pipeline     | `install.enabled=true`, `simulation.run_test=true`            |
| Simulate only     | `install.enabled=false`, `simulation.run_test=true`           |
| Postprocess only  | `simulation.run_test=false`, `postprocess.data_dir="..."`     |
| Re-plot latest    | `simulation.run_test=false`, `postprocess.data_dir=""`        |

## Programmatic API

```julia
using Nbody6Setup

# Full pipeline from config
cfg = load_config("config.toml")
results = run_pipeline(cfg)

# Post-process an external output directory directly
results = postprocess_external("/path/to/output")

# Quick scan of available output files
scan = scan_output("/path/to/output")
println(scan)
```

See the [API Reference](@ref) for the complete public interface.

## Contents

```@contents
Pages = ["manual.md", "input_files.md", "api.md"]
Depth = 2
```
