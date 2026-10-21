# Nbody6Dynamics.jl

| Documentation | Build and tests | Code quality | Release |
|:-:|:-:|:-:|:-:|
| [![Docs stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://PaulGoG.github.io/Nbody6Dynamics.jl/stable/)<br>[![Docs dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://PaulGoG.github.io/Nbody6Dynamics.jl/dev/) | [![CI](https://github.com/PaulGoG/Nbody6Dynamics.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/PaulGoG/Nbody6Dynamics.jl/actions/workflows/CI.yml)<br>[![Backend build](https://github.com/PaulGoG/Nbody6Dynamics.jl/actions/workflows/Backend.yml/badge.svg)](https://github.com/PaulGoG/Nbody6Dynamics.jl/actions/workflows/Backend.yml)<br>[![Codecov](https://codecov.io/gh/PaulGoG/Nbody6Dynamics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/PaulGoG/Nbody6Dynamics.jl) | [![Format check](https://github.com/PaulGoG/Nbody6Dynamics.jl/actions/workflows/Format.yml/badge.svg)](https://github.com/PaulGoG/Nbody6Dynamics.jl/actions/workflows/Format.yml)<br>[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)<br>[![JET](https://img.shields.io/badge/tested_with-JET.jl-233f9a.svg)](https://github.com/aviatesk/JET.jl) | [![Release](https://img.shields.io/github/v/release/PaulGoG/Nbody6Dynamics.jl?display_name=tag&sort=semver)](https://github.com/PaulGoG/Nbody6Dynamics.jl/releases)<br>[![Julia](https://img.shields.io/badge/Julia-1.13%2B-9558B2.svg?logo=julia&logoColor=white)](https://julialang.org)<br>[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)<br>[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) |

A Julia package that runs [Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing) star-cluster simulations end to end: it clones and builds the Fortran engine, generates initial conditions for single clusters and multi-cluster mergers, launches and supervises the integration (with restart dumps, wall budgets and resumption across queue segments), reads every standard output file, computes the cluster diagnostics and draws the figures. One TOML configuration drives every phase through a single entry point, `run_pipeline`.

![Projected surface density of two merging star clusters: the pair at release, the two cores at 4 Myr and the remnant at 50 Myr](docs/src/assets/merger_evolution.png)

*Two King clusters of 25 000 stars each, released 10 pc apart at the apocentre of an orbit of
eccentricity 0.5, drawn as projected surface mass density: the pair at release, the two cores
still distinct at 4 Myr, and the remnant with its halo at 50 Myr. The integration took 7.6 min
of engine time on one NVIDIA H200 NVL. Initial conditions, integration, analysis and figure were
all produced by this package; the case is `input_files/gpu/merger_50k.toml`.*

![Animation of the same merger over 50 Myr](docs/src/assets/merger_evolution.gif)

*The same run animated on fixed axes, one frame per 2 Myr. The cores fall together and merge
between 6 and 8 Myr; the remnant then relaxes while its halo of loosely bound stars spreads past
the frame, which holds 95 % of the mass at every epoch. Every run writes animations like this
one alongside its static figures.*

## Project structure

```
Nbody6Dynamics/
├── Project.toml      # Package metadata and dependencies
├── activate.jl       # Activates and instantiates the package environment
├── config.toml       # Main pipeline configuration (edit this)
├── CONTRIBUTING.md   # What a change needs
├── src/              # Library: orchestration, IC generator, readers, diagnostics
├── ext/              # Makie extension: figures and animations
├── scripts/          # CLI entry points
├── input_files/      # Shipped inputs by purpose: engine/, clusters/, mergers/, verification/, sweeps/, showcase/, gpu/
├── test/             # Unit and physics-validation suite
├── bench/            # Benchmarks and scaling scripts
├── docs/             # Documenter.jl site sources
├── assets/           # Shipped tables: the estimator calibration (assets/calibration/)
├── deps/cuda/        # NVIDIA samples headers used by the CUDA build
├── deps/engine/      # Corrections applied to the engine source before the build
├── backend/          # (gitignored) engine checkout and build
└── runs/             # (gitignored) run outputs
```

The full tree is at the end of this file.

## Installation

Requires Julia ≥ 1.13, installed through [juliaup](https://github.com/JuliaLang/juliaup) (`juliaup add release`). The package is not registered; add it to a project of your own by URL, together with a Makie backend for the figures:

```julia
using Pkg
Pkg.add(url = "https://github.com/PaulGoG/Nbody6Dynamics.jl")
Pkg.add("CairoMakie")
```

```julia
using CairoMakie, Nbody6Dynamics
cfg = load_config("config.toml")        # your project's configuration
results = run_pipeline(cfg)
```

The directory of that `config.toml` is the project directory: the engine is cloned and built under its `backend/` and the runs land under its `runs/`, never inside the package. The shipped inputs are reachable with `example_input("N1k_quick.inp")`, and `example_input("showcase/equal_pipeline.toml")` returns a complete pipeline configuration to copy and edit. Building the Fortran engine needs `git`, `gfortran` and `make`, and optionally HDF5 and CUDA (auto-detected; see `[build]` in `config.toml`).

## Environment setup of a checkout

Development and the validated hosts run Julia 1.13. Manifests are not under version control: each environment resolves from its `Project.toml` on first activation, and every run stores the manifest it resolved as `environment_manifest.toml` in its run directory.

Every environment ships an activation script that activates and instantiates it silently. Running one on a new machine performs the dependency resolution and precompilation once:

```bash
julia activate.jl          # package environment (required)
julia scripts/activate.jl  # entry scripts: the package plus its figure backend
julia docs/activate.jl     # documentation build (optional)
julia bench/activate.jl    # benchmarks (optional)
```

An interactive session starts with `julia -i activate.jl`.

The scripts under `scripts/`, `docs/` and `bench/` include their environment's activation script, so they need no project flag; those environments develop the package by a relative path and always run against the local source.

### Figures are an extension

The figure routines are a package extension triggered by a Makie backend, so a headless host installs no plotting stack. Initial conditions, engine build, integration, readers, diagnostics and benchmarks all run without them.

```julia
using CairoMakie        # brings in every plot_*, animate_* and generate_plots
using Nbody6Dynamics
```

Without a backend the figure routines report the remedy instead of failing obscurely, and a pipeline whose `[visualization] enabled = true` is refused before it builds or integrates anything rather than after. `plotting_available()` answers the question in code. The entry scripts under `scripts/` load CairoMakie for you; a numerics-only session (`julia activate.jl`, then `using Nbody6Dynamics`) does not.

## Entry points

One invocation each; details in the sections below.

| Task | Invocation |
|------|------------|
| Instantiate the package environment | `julia activate.jl` |
| Instantiate the script environment (figures) | `julia scripts/activate.jl` |
| Run the main pipeline | `julia scripts/run_setup.jl [config.toml]` |
| Continue a stopped or interrupted run | `julia scripts/run_resume.jl <run_dir>` (one segment per call; post-processing and figures once the run is complete) |
| Submit a run as a chain of Slurm jobs | `julia scripts/run_cluster.jl config.toml [--dry-run] [--run-id=<id>]` (`[cluster] enabled = true`; one job per engine segment, the job script generated into the run directory) |
| Generate a single cluster or a merger and run it | `julia scripts/run_setup.jl input_files/showcase/equal_pipeline.toml` (the showcase cases; a single-cluster case: `[merger] config_file = "input_files/clusters/imbh_runaway.toml"`) |
| Run a parameter sweep | `julia scripts/run_sweep.jl input_files/sweeps/sweep_demo.toml [--dry-run]` |
| Run the verification suite | `julia scripts/run_verif_suite.jl` |
| Regenerate the estimator calibration table | `julia -t 6 scripts/calibrate_estimators.jl` (`--quick` for a schema check; writes `assets/calibration/`) |
| Execute the test suite | `julia -e 'include("activate.jl"); Pkg.test()'` |
| Run the benchmarks | `julia bench/benchmarks.jl` |
| Validate a CUDA host | `julia scripts/run_gpu_validation.jl` (`--dry-run`, `--stages=`, `--stop-on-failure`; results under `runs/gpu_validation_<host>_<timestamp>/`) |
| Build and run on a CUDA host by hand | `julia scripts/run_setup.jl input_files/gpu/gpu_pipeline.toml` (CPU reference: `cpu_pipeline.toml`; recipe in the manual) |
| Measure the GPU speed-up | `julia bench/gpu_scaling.jl 20000,50000 4,8 "0;0,1" 0.25` (needs the CPU and the GPU binary) |
| Measure the device's FP32 peak | `julia bench/fp32_peak_probe.jl --out=runs` (`--backend=cpu` without a device) |
| Build the documentation | `julia docs/make.jl` (executes the walkthrough; the site lands in `docs/build/`) |

## Usage

`run_pipeline(cfg)` is the single entry point; `config.toml` flags select which phases run. Six modes:

| # | Mode | Config flags |
|---|------|--------------|
| 1 | Full pipeline | `install.enabled=true`, `simulation.run_test=true` |
| 2 | Simulate only | `install.enabled=false`, `simulation.run_test=true` |
| 3 | Postprocess only | `simulation.run_test=false`, `postprocess.data_dir="/path/to/output"` |
| 4 | Re-plot latest run | `simulation.run_test=false`, `postprocess.data_dir=""` |
| 5 | Merger ICs only | `merger.enabled=true`, `merger.config_file="input_files/..."` |
| 6 | Merger + simulate | `merger.enabled=true`, `simulation.run_test=true` |

```julia
using Nbody6Dynamics
cfg = load_config("config.toml")
results = run_pipeline(cfg)   # Dict with :snapshots, :diagnostics, :lagr, :escapers, :stellar_evo
```

Or via the CLI wrapper:

```bash
julia scripts/run_setup.jl [config.toml]
```

Initial conditions alone, without a pipeline configuration (a single cluster or a merger; the file decides):

```julia
result = run_merger_pipeline("input_files/clusters/gc_bh_subsystem.toml")
# writes dat.10, merger.inp, merger_summary.txt, merger_ic.toml and the IC figures
```

Post-process any directory containing Nbody6++ output, config-free:

```julia
scan = scan_output("/scratch/sim42/output")     # report available files
results = postprocess_external("/scratch/sim42/output")
```

Each run gets an isolated `runs/<run_id>/` directory (`output/`, `plots/`, frozen `config.toml`, `RUN_INFO.toml` with one entry per engine segment). Run IDs come from `generate_run_id` (prefix + timestamp + 4-hex uniqueness suffix); runs from generated initial conditions prepend `merger_` to the configured prefix.

## Status

| Component | Status | Notes |
|-----------|--------|-------|
| Install / build | Working | Clone of the engine fork at its release tag (the default `install.source_url` and `install.ref`; the upstream repository at any ref also works), source amendments for an upstream checkout (`deps/engine/amendments.toml`: the escaper summary format that aborted a run at more than 1000 escapers in one step, the tidal energy in the total), `configure`, HDF5 Makefile patch, parallel make; CUDA path auto-detection; GPU builds compiled for the visible devices' compute capabilities or an explicit `cuda_arch` list (the upstream configure emits no architecture flag), `BUILD_INFO.toml` next to the binary, CPU/GPU/MPI binary variants selected by suffix. Validated on five CUDA hosts (compute capability 7.5 to 12.0, CUDA 11.8 to 13.3); see the manual, section Validated hardware |
| Initial-condition generator | Working | Plummer and King samplers (King c(W0) validated against published concentrations); Kroupa (2001) IMF; one cluster at rest, two clusters on a Kepler orbit given by apocentre or pericentre (in pc or member half-mass radii) with its orientation and the tidal-frame convention, N clusters on explicit orbits, or a cluster complex drawn from a member mass function, a size–mass relation and a centre distribution (uniform sphere, Plummer, box fractal, filament) with the bulk velocities scaled to a target virial ratio, its energy-matched monolithic controls and one star list shared across arms (the shipped science cases and the sweep controls are single-cluster configurations); primordial binaries (Kroupa 1995 periods, thermal eccentricities, written in the engine's pair convention); Jacobi truncation; integration parameters from the engine's rules or from the baseline preset; seeded reproducibility (the TOML `seed` drives the sampler and the engine's `NRAND`); output intervals written as dyadic rationals with exact decimals, which the engine's digit counter requires. Static test models with a known answer for one estimator each: the Osipkov–Merritt Plummer sphere, the Lynden-Bell rotation flip, stretched and energy-ranked models. The engine has no multi-centre diagnostics: cluster-level results before coalescence come from the snapshot-based per-cluster tools, not from `lagr.7` or `esc.11` (merger documentation, "Feasibility and limitations") |
| Snapshot integrals and bound sets | Working | Exact Newtonian potential in Float64 (a threaded, vectorised loop over tiles that visits each pair once; a KernelAbstractions kernel as the device form); regularised pairs as centres of mass; iterative energetic unbinding with an explicit convergence status; per snapshot the mass, centre of mass, momenta, angular momentum, kinetic and potential energy of the reduced system (the engine's `EKIN` and `POT`), the pairs' internal energy and the bound set; the Jacobi integral under a linearised tide (`e_jacobi`); `snapshot_integrals.csv` per run |
| Energy-jump localisation and mass budget | Working | The stamped event lines of the engine log as records with their adjustment interval; per interval the engine's `DE` with the increments of its correction channels and the event counts, flagged by a floor and a robust scale (`energy_intervals.csv`, `engine_events.csv`); per snapshot the mass against the escapers, the stellar-evolution and the collision losses (`mass_budget.csv`) and the energy against the removed escapers and the engine's corrections (`energy_budget.csv`); the Jacobi integral of a linearised tide, under which the engine's `DE` is a drift |
| Remnant diagnostics | Working | Bound remnant of the whole system per snapshot: both Casertano–Hut core radii, from the mass density (engine densities or sixth-neighbour estimate on a k-d tree) and from the number density, core density, half-mass radius, rotation (λ_R, Peebles λ_P, spin alignment with the orbital angular momentum, v_rot/σ profile, floors from a sign-flip null, Bianchini rotation-curve fit), mass segregation of a fixed-identity tracer class with Λ_MSR (Allison et al. 2009), Λ̃_MSR and Γ_MSR over a scan of tree sizes and the segregation time, union-find coalescence time from per-tag density centres and half-mass radii with persistence, snapshot binary census (pairs, hard pairs, largest hardness), the core-collapse marker set; `remnant_diagnostics.csv` and four figures per merger run |
| Angular momentum, equipartition, merger tree | Working | Angular-momentum budget of the stars against escapers, newcomers and stellar mass loss, rotation of the bound set inside its Lagrangian radii and per initial-mass class with sign-flip floors, elapsed relaxation times (`angular_momentum.csv`); Bianchini equipartition mass of the luminous stars within r_h (`equipartition.csv`); merger tree with persistent joins and the per-progenitor mass budget with the mixing of every progenitor in the product (`merger_tree.csv`, `progenitor_budget.csv`) |
| Binary spans | Working | `binary_spans` (`[postprocess] binary_spans`, off by default): the KS episodes of every pair of component names joined across short gaps, chains and hierarchies, with the origin (primordial, reformed, exchange, dynamical with GW, tidal or three-body capture) and the end class (GW and stellar mergers, collisions, escapes, exchanges, supernova disruptions, flybys, disruptions) from the engine's event lines (`binary_spans.csv`); the physically bound pairs of every snapshot and their spans across the snapshots, KS pairs included, with interval-censored ends (`bound_spans`, `bound_spans.csv`) |
| Validity record | Working | `validity.toml` and `validity_snapshots.csv` per run: energy, mass and angular-momentum checks with a verdict, bound-set convergence, snapshot cadence against the crossing and relaxation times, the windows excluded from equilibrium estimators (violent relaxation, pericentres, coalescence), the times of virialisation, kinematic relaxation and mixing, the calibrated bias and scatter of the tracked estimators at the bound star count, and the masking rules; thresholds as `[postprocess]` keys |
| Shape, segregation and mixing | Working | Intrinsic axis ratios from the iterative unweighted shape tensor in ellipsoidal shells (Zemp et al. 2011), projected ellipticities in three projections, the A⁺ concentration index with a jackknife error, the two-sample Kolmogorov–Smirnov distance with its null floor, and the mixing of two tagged populations in radius, energy and angular momentum |
| Stellar population | Working | `STELLAR_CLASSES` (engine grouping of K* = −1…15), `hr_population` (single stars and members of KS pairs together), `stellar_census` per epoch and class, written as `stellar_census.csv`; HR figures with a run-level legend and a census strip |
| Parameter sweeps | Working | Sweep TOML → Cartesian grid over dotted merger-TOML keys × seeds; one directory per point with derived configs validated before launch; points run as concurrent worker processes (`omp_threads` per job from the cost model; on a GPU node `devices` and `engines_per_device` assign a device per point and admit k engines per device); `sweep_index.toml` kept current, `sweep_summary.csv` with final N, pairs, energy error and virial ratio; comparison figures on common axes coloured by one grid axis; seeded ensembles (a sweep without axes, or the seeds of every grid point) summarised by median and central 68/95 % bands on a common time grid; `controls = true` adds the isolated single-cluster equivalent of every point and a paired merger-versus-control figure |
| Simulation runner | Working | Launch script with `OMP_NUM_THREADS` and `GPU_LIST` control, live stdout monitoring, start-up watchdog and completion monitor (`exit_grace`: an engine that printed `END RUN` but never exits is terminated and recorded as completed), opt-in in-terminal sparklines of the diagnostics (`live_diagnostics`), run summary with exact CPU accounting, sampled CPU/memory/GPU telemetry (`telemetry.csv`), the devices the engine initialised and the build record. Long runs: a wall budget per segment (`wall_budget`, `stop_margin`) stops the engine through its `STOP` file with a restart dump, as does SIGTERM to the driver or to the launch script, which is the engine's parent and turns the signal into the stop request; `resume_run` continues an interrupted run from its last complete dump and moves aside what a killed segment wrote after it; every launch is recorded in `RUN_INFO.toml` and closed with a status (`completed`, `stopped`, `halted`, `watchdog`, `killed`, `failed`). Under Slurm, `submit_run` turns a run into a chain of jobs, one per segment, from a job script generated into the run directory, with a job ledger; a cost estimate (wall time, segments, disk from a rate table) is logged at every launch and checked against a budget. Generated runs execute inside their IC output directory so the engine finds `dat.10`; an operator interrupt asks the engine to stop, and the run record carries the commits captured at launch |
| I/O readers | Working | `conf.3` (standard + extended; the stars only, with the regularised pairs' components first and the engine's centre-of-mass entries dropped; `read_all_conf3` reads in threaded chunks when Julia has more than one thread), `out1000` diagnostics (ADJUST + physical scaling; virial ratio Q = T/\|W\|, equilibrium at 0.5) and `out1000` event lines, `event.35` event counters, `lagr.7`, and `esc.11` (incl. the ANGLE PHI / ANGLE THETA escape-direction columns) / `sev.83_*` / `bev.82_*` in the fork's real formats; `STELLAR_TYPE_LABELS` follow the Hurley convention (13 = NS, 14 = BH); `UnitScaling.zmbar` is the total-mass scale factor M*, not the mean stellar mass. `.h5part` files are detected and reported as unsupported |
| Plotting / animation | Working | Publication theme: no titles, no minor ticks, Computer Modern fonts, dashed grey low-opacity grid on line plots; presentation knobs config-driven via `[visualization.style]` (`PlotStyle`); escaper suite (cumulative mass loss, velocity classes, escape anisotropy) and SSE-quantity plots (mass segregation, t/T_MS evolutionary clock, core-mass growth); binary suite (pair counts with the Heggie hard/soft split, binary fraction, `a`–`e` diagrams, period histograms); existing figures are never overwritten (safesave-style `#1`, `#2`, … backups) |
| Load latency | Working | PrecompileTools workload over the configuration, diagnostics-reader and merger-IC paths; the first `load_config` of a session runs compiled |
| Tests | Passing | Unit and physics-validation suite (engine-dependent tests behind `NBODY6_BINARY_TESTS=1`, a GPU-gated set behind `NBODY6_GPU_TESTS=1`), incl. physics validation, adversarial external-input tests, telemetry/thread-control, per-cluster structure, restart, and tidal-field tests |

## Testing

```bash
julia -e 'include("activate.jl"); Pkg.test()'
```

The default suite needs no Fortran binary. The engine-dependent tests (build, launch, restart, tidal field, live telemetry) run when `NBODY6_BINARY_TESTS=1` is set; without `NBODY6_BACKEND_ROOT` they clone and build the backend in a temporary directory, with it they use the build under `<root>/backend/Nbody6PPGPU-beijing`:

```bash
NBODY6_BINARY_TESTS=1 NBODY6_BACKEND_ROOT=$PWD julia -e 'include("activate.jl"); Pkg.test()'
```

The suite covers config round-trips, all I/O readers (synthetic binaries plus real fork output fixtures under `test/fixtures/`), plotting/animation smoke tests, external post-processing (including adversarial malformed inputs), and the merger IC generator. Physics-validation tests check King concentration c(W0) against published values, the Plummer half-mass relation r_hm = 1.305 a, the Kroupa mean mass, virial equilibrium Q = T/|W| = 0.5 after `virialise!`, and the Keplerian orbital energy of generated two-cluster orbits.

For interactive work against the test dependencies, activate the test environment with [`TestEnv.jl`](https://github.com/JuliaTesting/TestEnv.jl): start `julia -i activate.jl` and run `using TestEnv; TestEnv.activate()` in the session.

## Verification Suite

```bash
julia scripts/run_verif_suite.jl
```

Runs three end-to-end targets through `run_pipeline` (requires a built binary in `backend/`):

1. `single` — plain N=5000 run (`N5k_medium.inp`)
2. `triorbit` — bound three-cluster Lagrange triangle (`verif_triorbit.toml`)
3. `3d5cluster` — five clusters distributed in 3D (`verif_3d5cluster.toml`)

## Documentation

The documentation site (Manual, Input Files, Showcase Cases, Cluster Mergers, API Reference, a walkthrough from a merger TOML to the generated initial conditions) is at [paulgog.github.io/Nbody6Dynamics.jl/dev](https://PaulGoG.github.io/Nbody6Dynamics.jl/dev/). It builds from `docs/`:

```bash
julia docs/make.jl   # builds to docs/build/
```

## Contributing

Issues and pull requests are welcome; `CONTRIBUTING.md` says what a change needs (tests, formatting, docstrings, a changelog entry, and the formulation behind any physics).

## Licence

The package is released under the MIT licence (`LICENSE`). Two other parties' material is involved. The CUDA helper headers under `deps/cuda/` come from NVIDIA's cuda-samples and carry their BSD-3-Clause licence (`deps/cuda/LICENSE`). The engine is not part of this package: the install phase clones it at build time, by default from my fork [PaulGoG/Nbody6PPGPU-beijing](https://github.com/PaulGoG/Nbody6PPGPU-beijing) at a tagged release of its `main` branch, and nothing of it is redistributed here. The upstream repository, [nbody6ppgpu/Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing), carries no licence file; the fork's `NOTICE.md` states the provenance, the credit due to the code's authors, the MIT terms of the earlier public repository of the code under which I treat the inherited code, and the MIT licence of my own modifications.

## How to cite

If this package contributes to published work, please cite it through the metadata in `CITATION.cff`:

```bibtex
@software{Nbody6Dynamics_jl,
  author  = {Gogîță, Paul-Adrian},
  title   = {Nbody6Dynamics.jl},
  version = {0.5.0},
  year    = {2026},
  url     = {https://github.com/PaulGoG/Nbody6Dynamics.jl},
  license = {MIT}
}
```

## Full tree

<details><summary>Full project tree</summary>

```
Nbody6Dynamics/
├── .github/                         # CI, Format, Documentation and Backend workflows; Dependabot; issue and PR templates
├── .JuliaFormatter.toml             # Formatter configuration
├── .mailmap                         # Author identities folded into one
├── README.md
├── CONTRIBUTING.md                  # What a change needs
├── LICENSE                          # MIT
├── CITATION.cff                     # Citation metadata
├── CHANGELOG.md                     # Release history (Changelog.jl conventions)
├── Project.toml                     # Package metadata & dependencies
├── activate.jl                      # Activates and instantiates the package environment
├── config.toml                      # Main pipeline configuration (edit this)
├── deps/cuda/                       # helper_cuda.h, helper_string.h from NVIDIA cuda-samples v13.0 (CUDA 13 build of the engine)
├── assets/calibration/              # estimators.csv (bias, scatter, errors of every estimator versus N on the static models) and estimators.toml (version, grid, seeds)
├── deps/engine/amendments.toml      # Source corrections applied to the engine tree before configure (escaper summary format, tidal energy in the total under a linearised tide)
├── deps/activation.jl               # Shared by the activate.jl scripts: resolves a manifest left stale by a dependency change
├── src/
│   ├── Nbody6Dynamics.jl               # Module root; run_pipeline orchestrator; exports
│   ├── types.jl                     # Config structs (incl. PlotStyle), Snapshot, records, UnitScaling
│   ├── config.jl                    # TOML loader/serialiser for Nbody6Config
│   ├── util.jl                      # safesave-style never-overwrite backups
│   ├── platform.jl                  # OS/CUDA/dependency detection, hardware fingerprint
│   ├── telemetry.jl                 # Runtime CPU/memory/GPU telemetry of the backend process tree
│   ├── install.jl                   # Clone, configure, HDF5 Makefile patch, build
│   ├── run.jl                       # Simulation launcher: run dirs, launch script, live monitoring
│   ├── resume.jl                    # Restart dump selection, the join of segments (what a killed segment wrote after its dump is moved aside)
│   ├── cluster.jl                   # Slurm layer: generated job script, submission and chaining of segments, job ledger
│   ├── cost.jl                      # Cost model: rate table, wall-time, segment and disk estimate at launch, budget check
│   ├── gpu_validation.jl            # run_gpu_validation: logged validation stages of a CUDA host under runs/
│   ├── external.jl                  # scan_output + postprocess_external for arbitrary output dirs
│   ├── potential.jl                 # Exact pair-sum potential (KernelAbstractions), pairs as centres of mass, bound sets
│   ├── integrals.jl                 # Snapshot integrals: momenta, energies, bound set; snapshot_integrals.csv
│   ├── energy_events.jl             # Energy-jump localisation: adjustment intervals, correction channels, event counts; energy_intervals.csv, engine_events.csv
│   ├── budgets.jl                   # Mass and energy budgets per snapshot: escapers, stellar-evolution and collision losses, engine corrections; mass_budget.csv, energy_budget.csv
│   ├── cluster_structure.jl         # Per-cluster structure from snapshots: bound members, centres, radii, dispersions
│   ├── binary_population.jl         # Binary diagnostics: Heggie hard/soft split, binary fraction, energy scale from snapshots
│   ├── stellar_population.jl        # Stellar classes, HR-plane population, class census
│   ├── remnant.jl                   # Remnant diagnostics: coalescence, Casertano–Hut core radius, rotation, mass segregation
│   ├── shape_mixing.jl              # Shape tensor, projected ellipticity, A⁺ index, two-sample D, mixing of two populations
│   ├── pair_finder.jl               # k-d tree local density, snapshot binary finder, reconciliation with the KS pairs
│   ├── equipartition.jl             # Bianchini equipartition-mass fit and its per-snapshot series; equipartition.csv
│   ├── angular_momentum.jl          # Angular-momentum budget, rotation by radius and mass class, elapsed relaxation times; angular_momentum.csv
│   ├── merger_tree.jl               # Merger tree from the density-centre separations, per-progenitor mass budget and mixing; merger_tree.csv, progenitor_budget.csv
│   ├── validity.jl                  # Validity record of a run and the calibration-table lookup; validity.toml, validity_snapshots.csv
│   ├── binary_spans.jl              # Engine-regularised binary spans: KS episodes joined per pair, origin, end class; binary_spans.csv
│   ├── bound_spans.jl               # Physically bound pairs at the snapshots and their spans; bound_spans.csv
│   ├── sweep.jl                     # Parameter sweeps: grid × seeds, concurrent worker processes, index and summary
│   ├── ensemble.jl                  # Seeded ensembles: series extractors, common-grid percentiles (median, 68/95 %)
│   ├── io/
│   │   ├── io.jl                    # I/O submodule includes
│   │   ├── fortran_binary.jl        # Fortran unformatted record reader
│   │   ├── conf3.jl                 # conf.3 snapshot reader (standard + extended layout; stars only)
│   │   ├── diagnostics.jl           # out1000 parser: ADJUST lines with the energy terms, PHYSICAL SCALING, options and tidal echoes
│   │   ├── lagr.jl                  # lagr.7 Lagrangian radii reader
│   │   ├── esc.jl                   # esc.11 escaper reader (fork's real column format)
│   │   ├── events.jl                # out1000 event lines, event.35 counters
│   │   ├── stellar_evolution.jl     # sev.83_* single-star evolution reader
│   │   └── binary_evolution.jl      # bev.82_* regularised-binary reader
│   ├── ic/
│   │   ├── ic.jl                    # run_merger_pipeline / generate_merger_ic entry points
│   │   ├── config.jl                # IC TOML parser: ClusterSpec, OrbitSpec, Nbody6ParameterSpec (presets), stellar and tidal specs
│   │   ├── control.jl               # Isolated single-cluster control derived from a merger TOML
│   │   ├── complex.jl               # Cluster complexes: member masses and radii, centre distributions, bulk velocity scaling
│   │   ├── models.jl                # Plummer & King (ODE-solved) density samplers
│   │   ├── static_models.jl         # Static test models: Osipkov–Merritt Plummer, Lynden-Bell flip, stretched, energy-ranked masses
│   │   ├── imf.jl                   # Kroupa (2001) IMF, rescaled & equal-mass variants
│   │   ├── orbits.jl                # Kepler two-body (size, orientation, derived quantities), explicit and single placement, tidal frame, Jacobi radii
│   │   ├── binaries.jl              # Primordial binary populations
│   │   └── output.jl                # dat.10 writer, merger.inp generator, integration-parameter resolution per preset
│   └── plotting_api.jl              # Figure interface: the public routines, implemented by the Makie extension
├── ext/
│   ├── Nbody6DynamicsMakieExt.jl    # Extension module, loaded by `using CairoMakie`
│   └── plotting/
│       ├── common.jl                # Publication theme (Computer Modern fonts), figure sizing, axis helpers
│       ├── dispatch.jl              # generate_plots: the single figure dispatcher
│       ├── snapshots.jl             # Projection figures (surface density or scatter), cluster separation & per-cluster virial
│       ├── energy.jl                # Energy error & particle count from diagnostics
│       ├── lagrangian.jl            # Lagrangian radii evolution
│       ├── hr.jl                    # HR figures: run-level legend, panels, census strip; HR animation
│       ├── escapers.jl              # Escaper analysis: cumulative mass loss, velocities, anisotropy
│       ├── sse.jl                   # SSE-quantity plots: mass segregation, t/T_MS clock, core masses
│       ├── merger.jl                # Merger-run figures: separation, per-cluster virial and structure
│       ├── merger_ic.jl             # Merger IC diagnostic figures
│       ├── binaries.jl              # Binary population, a–e diagram, period distribution
│       ├── remnant.jl               # Remnant figures: rotation parameters and profile, structure, mass segregation
│       ├── sweep.jl                 # Sweep comparison figures on common axes (Lagrangian radius, energy error)
│       ├── control.jl               # Merger against isolated control, paired series per sweep point
│       ├── ensemble.jl              # Ensemble figures: median and 68/95 % bands per grid point
│       ├── animation.jl             # GIF animations (cluster, HR, Lagrangian)
│       └── telemetry.jl             # Run telemetry figure
├── scripts/
│   ├── run_setup.jl                 # CLI wrapper: load config → run_pipeline
│   ├── run_resume.jl                # CLI wrapper: resume_pipeline on a run directory
│   ├── run_cluster.jl               # CLI wrapper: submit_run — a run as a chain of Slurm jobs
│   ├── run_cluster_segment.jl       # Job entry: cluster_segment on a run directory (called by the generated job script)
│   ├── run_sweep.jl                 # CLI wrapper: sweep TOML → run_sweep → comparison figures
│   ├── run_verif_suite.jl           # Three-target verification suite
│   ├── calibrate_estimators.jl      # Bias and variance of every estimator on the static models versus N → assets/calibration/
│   ├── run_gpu_validation.jl        # CLI wrapper: logged validation stages of a CUDA host
│   ├── activate.jl                  # Activates the script environment (package + CairoMakie)
│   └── Project.toml                 # Script environment: the package and its figure backend
├── test/
│   ├── runtests.jl                  # Suite entry: shared helpers, includes the component files below
│   ├── test_config.jl               # Configuration loading, round trip, validation
│   ├── test_platform_build.jl       # Platform and CUDA detection, engine source tree, GPU build target
│   ├── test_validation.jl           # GPU validation driver and stage verdicts
│   ├── test_run.jl                  # Run IDs, machine identity, watchdogs, restarts, provenance
│   ├── test_resume.jl               # Dump completeness, dump selection, the join on synthetic run directories
│   ├── test_cluster.jl              # Job script rendering, submission chains against a stand-in sbatch, segments and ledger
│   ├── test_cost.jl                 # Rate table, estimate arithmetic, input reader, launch-time budget check
│   ├── test_io.jl                   # Readers: conf.3, out1000, lagr.7, esc.11, sev.83, bev.82, fixtures
│   ├── test_diagnostics.jl          # Binary population, remnant, stellar classes, cluster structure
│   ├── test_energy_budget.jl        # Jacobi integral, energy-jump localisation, mass budget
│   ├── test_potential.jl            # Pair potential, reduced systems, bound sets, snapshot integrals
│   ├── test_static_models.jl        # Estimators against the static test models with known answers
│   ├── test_sweep_ensemble.jl       # Sweeps, seeded ensembles, control configurations
│   ├── test_plotting.jl             # Figure extension, smoke tests, canvases, formatting
│   ├── test_external.jl             # External post-processing (+ test_external_adversarial_inner.jl)
│   ├── test_ic.jl                   # Merger initial-condition generator and its physics checks
│   ├── test_telemetry.jl            # Hardware telemetry sampler, readers and figure
│   ├── test_engine_gated.jl         # Shipped GPU-host configurations; engine- and GPU-gated runs
│   ├── test_qa.jl                   # Static QA: Aqua, ExplicitImports, JET
│   └── fixtures/                    # Real Nbody6++ output excerpts (esc.11, lagr.7, out1000, sev.83_0, bev.82_0) with their PROVENANCE.md
├── bench/
│   ├── thread_scaling.jl            # Thread- and N-scaling of the backend from run telemetry (cost model)
│   ├── gpu_scaling.jl               # GPU-versus-CPU binary at equal N and threads, per GPU_LIST (speed-up table)
│   ├── gpu_cells.jl                 # Benchmark cells on the device: (cell, N, threads, GPU list, MPI ranks) grid, one CSV row per run
│   ├── gpu_cells_figures.jl         # Scaling (strong, weak, host threads) and cross-hardware figures from the cell CSVs; --synthetic rehearsal
│   ├── fp32_peak_probe.jl           # Measured FP32 peak of the device (GPUDiagnostics), fp32_peak_<host>.toml for the fraction-of-peak panel
│   ├── merger_case.jl               # The single- and dual-cluster cells shared by the scaling scripts
│   ├── benchmarks.jl                # BenchmarkTools suite (kept out of tests)
│   ├── activate.jl                  # Activates the bench environment (package developed by relative path)
│   └── Project.toml                 # Bench-local environment
├── docs/
│   ├── make.jl                      # Documenter.jl build script
│   ├── activate.jl                  # Activates the docs environment (package developed by relative path)
│   ├── Project.toml                 # Documentation build environment
│   ├── src/                         # index, walkthrough (Literate), manual, input_files, showcase, multi_cluster_mergers, api, references
│   ├── src/references.bib           # BibTeX of the sources cited (DocumenterCitations)
│   └── src/assets/                  # figure and animation used by the README and the docs site, with their PROVENANCE.md
├── input_files/                     # Shipped inputs, one folder per purpose (input_files/README.md)
│   ├── engine/                      # Single-cluster engine inputs with engine-generated initial conditions: N1k_quick … N100k_production
│   ├── clusters/                    # Single-cluster generator configurations: the four science cases
│   ├── mergers/                     # Merger initial-condition TOMLs: demo, equal-mass, minor Plummer, triple, 3- and 5-cluster, 27-cluster grid
│   ├── verification/                # verif_triorbit (Lagrange triangle), verif_3d5cluster (five clusters in 3D)
│   ├── sweeps/                      # sweep_demo (eccentricity × secondary size × seeds over the demo merger)
│   ├── showcase/                    # Five end-to-end cases with Myr intervals: equal, binary, tidal, flagship, sweep with controls
│   └── gpu/                         # CUDA-host recipe: GPU and CPU builds, the 2 × 25k merger, the 6 × 10⁵-body probes
├── backend/                         # (gitignored) cloned Nbody6PPGPU-beijing source + build
└── runs/                            # (gitignored) per-run output/, plots/, frozen config.toml
```

</details>
