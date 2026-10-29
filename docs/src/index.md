# Nbody6Dynamics.jl

A Julia package that runs [Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing) star-cluster simulations end to end, from the engine build and the initial conditions to the diagnostics and figures.

![Projected surface density of two merging star clusters: the pair at release, the two cores at 4 Myr and the remnant at 50 Myr](assets/merger_evolution.png)

*Two King clusters of 25 000 stars each, released 10 pc apart on an orbit of eccentricity 0.5,
drawn as projected surface mass density at release, at 4 Myr and at 50 Myr; produced end to end
by this package.*

## Features

The package clones and builds the engine (by default my fork at `v2026.07.14`; the CUDA toolkit is located automatically, MPI and HDF5 are switched in `[build]`), generates initial conditions for one cluster, a Kepler pair, N clusters on explicit orbits or a cluster complex, and writes them as `dat.10` with a matching `merger.inp` (`KZ(22) = 2`). It runs the integration under wall budgets, stop requests and start-up and stall watchdogs, resumes it across queue segments or as a Slurm job chain, and records commits, hardware and cost per segment. It reads every standard output file, the binary and few-body event lines included, and derives snapshot integrals, energy and mass budgets, remnant structure, rotation and segregation, binary spans, the binary census and binary fractions, gravitational-wave observables, and a validity record per run. Figures and GIF animations come from a CairoMakie extension; sweeps over a grid with seeds and isolated controls run as concurrent workers.

## Quick Start

As a dependency of your own project (the package is not registered):

```julia
using Pkg
Pkg.add(url = "https://github.com/PaulGoG/Nbody6Dynamics.jl")
Pkg.add("CairoMakie")     # the figure backend
```

The directory of your `config.toml` is the project directory: the engine is built under its `backend/` and runs land under its `runs/`. `example_input("showcase/equal_pipeline.toml")` returns a shipped configuration to copy and edit, together with `equal_merger.toml` beside it, since its paths are relative to its own directory.

From a checkout:

```bash
git clone https://github.com/PaulGoG/Nbody6Dynamics.jl.git
cd Nbody6Dynamics.jl

# Install Julia dependencies (one-time)
julia activate.jl          # package environment
julia scripts/activate.jl  # entry scripts: the package plus its figure backend

# Edit configuration
$EDITOR config.toml

# Run the full pipeline
julia scripts/run_setup.jl config.toml
```

The figure routines are a package extension (see [Visualisation](@ref "10. Visualisation")): a
session that produces figures loads a Makie backend first, a numerics-only session needs
neither the backend nor its dependencies.

```julia
using CairoMakie
using Nbody6Dynamics
run_pipeline(load_config("config.toml"))
```

## Usage Modes

The six ways of running the pipeline (full, simulate only, post-process only, re-plot, initial conditions only, initial conditions plus simulation) are listed in the manual, [Pipeline Modes](@ref "5. Pipeline Modes").

## Programmatic API

```julia
using Nbody6Dynamics

# Full pipeline from config
cfg = load_config("config.toml")
results = run_pipeline(cfg)

# One-call merger IC generation (no config.toml needed)
result = run_merger_pipeline("input_files/mergers/merger_demo_small.toml")

# Post-process an external output directory directly
results = postprocess_external("/path/to/output")

# Quick scan of available output files
scan = scan_output("/path/to/output")
println(scan)

# Reload a previously generated merger IC (no re-sampling)
result = load_merger_ic_result("runs/merger_run_.../output")
```

See the [API Reference](@ref) for the complete public interface.

## Contents

```@contents
Pages = ["walkthrough.md", "manual.md", "input_files.md", "showcase.md", "multi_cluster_mergers.md", "api.md", "references.md"]
Depth = 2
```
