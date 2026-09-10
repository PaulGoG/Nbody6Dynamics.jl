# Nbody6Dynamics.jl

[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Julia](https://img.shields.io/badge/Julia-1.10%2B-9558B2.svg?logo=julia&logoColor=white)](https://julialang.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Project Status: WIP](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)

A Julia package that automates the full lifecycle of [Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing) star-cluster simulations: install/build of the Fortran code, multi-cluster merger initial-condition generation, simulation execution, post-processing of all standard output files, and publication-quality visualization. Every phase is driven by a single TOML configuration and orchestrated through one entry point, `run_pipeline`.

![Two star clusters merging: the initial conditions and the remnant after 12 Myr](docs/src/assets/merger_evolution.png)

*Two King clusters of 4000 stars each, released 6 pc apart on an eccentric orbit, and the
remnant they leave 12 Myr later. Note the axes: the pair spans 12 pc, the remnant and its halo
nearly 70 pc. Initial conditions, integration, analysis and figure were all produced by this
package from one configuration file; the case is `input_files/showcase/`.*

![Animation of the same merger over 12 Myr](docs/src/assets/merger_evolution.gif)

*The same run animated on fixed axes, one frame per Myr. The two clusters fall together and
coalesce at 4 Myr, after which the remnant relaxes and sheds the halo of loosely bound stars
that fills the frame. Every run writes animations like this one alongside its static figures.*

## Project Structure

```
Nbody6Dynamics/
├── .github/workflows/               # CI matrix, static QA, Format check, Docs build, CompatHelper (parked until the repository is public)
├── .JuliaFormatter.toml             # Formatter configuration
├── .mailmap                         # Author identities folded into one
├── README.md
├── LICENSE                          # MIT
├── CITATION.cff                     # Citation metadata
├── CHANGELOG.md                     # Release history (Changelog.jl conventions)
├── Project.toml                     # Package metadata & dependencies
├── Manifest.toml                    # Version-controlled — exact dependency versions
├── activate.jl                      # Activates and instantiates the package environment
├── config.toml                      # Main pipeline configuration (edit this)
├── src/
│   ├── Nbody6Dynamics.jl               # Module root; run_pipeline orchestrator; exports
│   ├── types.jl                     # Config structs (incl. PlotStyle), Snapshot, records, UnitScaling
│   ├── config.jl                    # TOML loader/serialiser for Nbody6Config
│   ├── util.jl                      # safesave-style never-overwrite backups
│   ├── platform.jl                  # OS/CUDA/dependency detection, hardware fingerprint
│   ├── telemetry.jl                 # Runtime CPU/memory/GPU telemetry of the backend process tree
│   ├── install.jl                   # Clone, configure, HDF5 Makefile patch, build
│   ├── run.jl                       # Simulation launcher: run dirs, launch script, live monitoring
│   ├── external.jl                  # scan_output + postprocess_external for arbitrary output dirs
│   ├── cluster_structure.jl         # Per-cluster structure from snapshots: bound members, centres, radii, dispersions
│   ├── binary_population.jl         # Binary diagnostics: Heggie hard/soft split, binary fraction, energy scale from snapshots
│   ├── remnant.jl                   # Remnant diagnostics: coalescence, Casertano–Hut core radius, rotation, mass segregation
│   ├── sweep.jl                     # Parameter sweeps: grid × seeds, concurrent worker processes, index and summary
│   ├── ensemble.jl                  # Seeded ensembles: series extractors, common-grid percentiles (median, 68/95 %)
│   ├── io/
│   │   ├── io.jl                    # I/O submodule includes
│   │   ├── fortran_binary.jl        # Fortran unformatted record reader
│   │   ├── conf3.jl                 # conf.3 snapshot reader (standard + extended layout)
│   │   ├── diagnostics.jl           # out1000 parser: ADJUST lines + PHYSICAL SCALING
│   │   ├── lagr.jl                  # lagr.7 Lagrangian radii reader
│   │   ├── esc.jl                   # esc.11 escaper reader (fork's real column format)
│   │   ├── stellar_evolution.jl     # sev.83_* single-star evolution reader
│   │   └── binary_evolution.jl      # bev.82_* regularised-binary reader
│   ├── ic/
│   │   ├── ic.jl                    # run_merger_pipeline / generate_merger_ic entry points
│   │   ├── config.jl                # Merger TOML parser: ClusterSpec, OrbitSpec, seed
│   │   ├── control.jl               # Isolated single-cluster control derived from a merger TOML
│   │   ├── models.jl                # Plummer & King (ODE-solved) density samplers
│   │   ├── imf.jl                   # Kroupa (2001) IMF, rescaled & equal-mass variants
│   │   ├── orbits.jl                # Kepler two-body + explicit N-cluster orbits, Jacobi radius, virialise!
│   │   ├── output.jl                # dat.10 writer, merger.inp generator (seed → NRAND)
│   │   └── plotting.jl              # Merger IC diagnostic plots
│   └── plotting/
│       ├── plotting.jl              # Publication theme (Computer Modern fonts)
│       ├── snapshots.jl             # Projection scatter plots, cluster separation & per-cluster virial
│       ├── energy.jl                # Energy error & particle count from diagnostics
│       ├── lagrangian.jl            # Lagrangian radii evolution
│       ├── hr.jl                    # HR diagrams colour-coded by stellar type
│       ├── escapers.jl              # Escaper analysis: cumulative mass loss, velocities, anisotropy
│       ├── sse.jl                   # SSE-quantity plots: mass segregation, t/T_MS clock, core masses
│       ├── merger.jl                # Merger IC overview figures
│       ├── binaries.jl              # Binary population, a–e diagram, period distribution
│       ├── remnant.jl               # Remnant figures: rotation parameters and profile, structure, mass segregation
│       ├── sweep.jl                 # Sweep comparison figures on common axes (Lagrangian radius, energy error)
│       ├── control.jl               # Merger against isolated control, paired series per sweep point
│       ├── ensemble.jl              # Ensemble figures: median and 68/95 % bands per grid point
│       ├── animation.jl             # GIF animations (cluster, HR, Lagrangian)
│       └── telemetry.jl             # telemetry.csv readers and the run telemetry figure
├── scripts/
│   ├── run_setup.jl                 # CLI wrapper: load config → run_pipeline
│   ├── run_sweep.jl                 # CLI wrapper: sweep TOML → run_sweep → comparison figures
│   └── run_verif_suite.jl           # Three-target verification suite
├── test/
│   ├── runtests.jl                  # Full unit + physics-validation suite
│   ├── test_external_adversarial_inner.jl  # Adversarial external post-processing tests
│   └── fixtures/                    # Real Nbody6++ output excerpts (esc.11, lagr.7, out1000, sev.83_0, bev.82_0)
├── bench/
│   ├── thread_scaling.jl            # Thread- and N-scaling of the backend from run telemetry (cost model)
│   ├── gpu_scaling.jl               # GPU-versus-CPU binary at equal N and threads, per GPU_LIST (speed-up table)
│   ├── merger_case.jl               # The two-cluster case shared by the scaling scripts
│   ├── benchmarks.jl                # BenchmarkTools suite (kept out of tests)
│   ├── activate.jl                  # Activates the bench environment (package developed by relative path)
│   ├── Project.toml                 # Bench-local environment
│   └── Manifest.toml                # Version-controlled — exact benchmark dependency versions
├── docs/
│   ├── make.jl                      # Documenter.jl build script
│   ├── activate.jl                  # Activates the docs environment (package developed by relative path)
│   ├── Project.toml                 # Documentation build environment
│   ├── Manifest.toml                # Version-controlled — exact documentation dependency versions
│   ├── src/                         # index, walkthrough (Literate), manual, input_files, multi_cluster_mergers, api, references
│   ├── src/references.bib           # BibTeX of the sources cited (DocumenterCitations)
│   └── src/assets/                  # figure and animation used by the README and the docs site
├── input_files/
│   ├── N1k_quick.inp                # N=1000 smoke test (seconds)
│   ├── N5k_medium.inp               # N=5000 medium verification run
│   ├── N25k_production.inp          # N=25000 production run
│   ├── N100k_production.inp         # N=100000 production run
│   ├── gc_bh_subsystem.inp          # Globular cluster with BH subsystem
│   ├── imbh_runaway.inp             # IMBH formation via runaway collisions
│   ├── pop3_cluster.inp             # Population III near-zero-metallicity cluster
│   ├── tidal_tails.inp              # Tidal-tail formation, Galactic-centre cluster
│   ├── young_massive_binaries.inp   # Young massive cluster, high binary fraction
│   ├── merger_demo_small.toml       # Quick equal-mass King merger demo (N=1000/cluster)
│   ├── merger_equal_mass.toml       # Equal-mass King merger (q=1), eccentric orbit
│   ├── merger_minor_plummer.toml    # Minor Plummer merger (q=0.1), inspiral setup
│   ├── merger_triple_cluster.toml   # Triple cluster, explicit-position orbit mode
│   ├── merger_3cluster_small.toml   # 3-cluster equilateral triangle, small
│   ├── merger_5cluster_small.toml   # 5-cluster pentagon, small
│   ├── N10k_long.inp                # 10k single cluster, extended TCRIT (long run)
│   ├── merger_27cluster_cubic.toml  # 27 clusters on a 3×3×3 cubic grid
│   ├── verif_triorbit.toml          # Bound Lagrange-triangle verification target
│   ├── verif_3d5cluster.toml        # 5 clusters distributed in 3D (projection/COM verification)
│   ├── gpu/                         # CUDA-host recipe: GPU and CPU reference builds + a 2 × 25k merger (see the manual)
│   ├── sweep_demo.toml              # Demonstration sweep (eccentricity × secondary size × seeds)
│   └── showcase/                    # Four science cases with Myr intervals: equal-mass eccentric merger,
│                                    #   binary-rich unequal merger, tidal-field merger, sweep with controls
├── backend/                         # (gitignored) cloned Nbody6PPGPU-beijing source + build
└── runs/                            # (gitignored) per-run output/, plots/, frozen config.toml
```

## Environment Setup

Requires Julia ≥ 1.10, installed through [juliaup](https://github.com/JuliaLang/juliaup) (`juliaup add release`). The tracked Manifests were resolved with Julia 1.13, the version of the GPU hosts the package targets, and reproduce that dependency set exactly there.

Every environment ships an activation script that activates and instantiates it silently. Running one on a new machine performs the dependency resolution and precompilation once:

```bash
julia activate.jl          # package environment (required)
julia docs/activate.jl     # documentation build (optional)
julia bench/activate.jl    # benchmarks (optional)
```

The scripts under `scripts/`, `docs/` and `bench/` include their environment's activation script, so they need no `--project` flag; the docs and bench environments develop the package by a relative path and always run against the local source. Building the Fortran backend additionally needs `git`, `gfortran`/`make`, and optionally HDF5 and CUDA (auto-detected; see `[build]` in `config.toml`).

## Entry Points

One invocation each; details in the sections below.

| Task | Invocation |
|------|------------|
| Instantiate the package environment | `julia activate.jl` |
| Run the main pipeline | `julia scripts/run_setup.jl [config.toml]` |
| Run a parameter sweep | `julia scripts/run_sweep.jl input_files/sweep_demo.toml [--dry-run]` |
| Run a showcase case | `julia scripts/run_setup.jl input_files/showcase/equal_pipeline.toml` (also `binary_`, `tidal_`; the sweep via `scripts/run_sweep.jl input_files/showcase/sweep.toml`) |
| Run the verification suite | `julia scripts/run_verif_suite.jl` |
| Execute the test suite | `julia --project=. -e 'using Pkg; Pkg.test()'` |
| Run the benchmarks | `julia bench/benchmarks.jl` |
| Build and run on a CUDA host | `julia scripts/run_setup.jl input_files/gpu/gpu_pipeline.toml` (CPU reference: `cpu_pipeline.toml`; recipe in the manual) |
| Measure the GPU speed-up | `julia bench/gpu_scaling.jl 20000,50000 4,8 "0;0,1" 0.25` (needs the CPU and the GPU binary) |
| Build the documentation | `julia docs/make.jl` (also executes the walkthrough; the site lands in `docs/build/`) |

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

Standalone merger IC generation (no main config needed):

```julia
result = run_merger_pipeline("input_files/merger_equal_mass.toml")
# writes dat.10, merger.inp, merger_summary.txt + diagnostic plots
```

Post-process any directory containing Nbody6++ output, config-free:

```julia
scan = scan_output("/scratch/sim42/output")     # report available files
results = postprocess_external("/scratch/sim42/output")
```

Each run gets an isolated `runs/<run_id>/` directory (`output/`, `plots/`, frozen `config.toml`). Run IDs come from `generate_run_id` (prefix + timestamp + 4-hex uniqueness suffix); merger runs prepend `merger_` to the configured prefix.

## Component Status

| Component | Status | Notes |
|-----------|--------|-------|
| Install / build | Working | Clone, `configure`, HDF5 Makefile patch, parallel make; CUDA path auto-detection; GPU builds compiled for the visible devices' compute capabilities or an explicit `cuda_arch` list (the upstream configure emits no architecture flag), `BUILD_INFO.toml` next to the binary, CPU/GPU/MPI binary variants selected by suffix. Untested on real GPUs so far: the GPU-gated suite (`NBODY6_GPU_TESTS=1`) is the acceptance test for the first NVIDIA host |
| Merger IC generator | Working | Plummer + King samplers (King c(W0) validated against published concentrations); Kroupa (2001) IMF; Kepler two-body and explicit N-cluster orbit modes; primordial binaries (Kroupa 1995 periods, thermal eccentricities, written in the engine's pair convention); Jacobi truncation; seeded reproducibility — the TOML `seed` drives the sampler RNG and propagates to Nbody6's `NRAND`; the output intervals are written as dyadic rationals with exact decimals, because the engine's digit counter never terminates on other values. The engine itself has no multi-centre diagnostics: cluster-level results before coalescence come from the snapshot-based per-cluster tools, not from `lagr.7`/`esc.11`; see "Feasibility and limitations" in the merger documentation |
| Remnant diagnostics | Working | Bound remnant of the whole system per snapshot: Casertano–Hut core radius (engine densities or sixth-neighbour estimate), half-mass radius, rotation (λ_R, Peebles λ_P, spin alignment with the orbital angular momentum, v_rot/σ profile), Allison et al. (2009) Λ_MSR mass segregation with segregation time, union-find coalescence time; `remnant_diagnostics.csv` and four figures per merger run |
| Parameter sweeps | Working | Sweep TOML → Cartesian grid over dotted merger-TOML keys × seeds; one directory per point with derived configs validated before launch; points run as concurrent worker processes (`omp_threads` per job from the cost model); `sweep_index.toml` kept current, `sweep_summary.csv` with final N, pairs, energy error and virial ratio; comparison figures on common axes coloured by one grid axis; seeded ensembles (a sweep without axes, or the seeds of every grid point) summarised by median and central 68/95 % bands on a common time grid; `controls = true` adds the isolated single-cluster equivalent of every point and a paired merger-versus-control figure |
| Simulation runner | Working | Launch script with `OMP_NUM_THREADS` and `GPU_LIST` control, live stdout monitoring, start-up watchdog and completion monitor (`exit_grace`: an engine that printed `END RUN` but never exits is terminated and recorded as completed), opt-in in-terminal sparklines of the diagnostics (`live_diagnostics`), run summary with exact CPU accounting, sampled CPU/memory/GPU telemetry (`telemetry.csv`), the devices the engine initialised and the build record; restarts from the engine's COMMON dumps (`restart_simulation`) with per-segment bookkeeping; merger runs execute inside the IC output dir so `dat.10` is found |
| I/O readers | Working | `conf.3` (standard + extended; `read_all_conf3` reads in threaded chunks when Julia has more than one thread), `out1000` diagnostics (ADJUST + physical scaling; virial ratio Q = T/\|W\|, equilibrium at 0.5), `lagr.7`, and `esc.11` (incl. the ANGLE PHI / ANGLE THETA escape-direction columns) / `sev.83_*` / `bev.82_*` in the fork's real formats; `STELLAR_TYPE_LABELS` follow the Hurley convention (13 = NS, 14 = BH); `UnitScaling.zmbar` is the total-mass scale factor M*, not the mean stellar mass. HDF5 reader removed — the fork's KZ(46) H5Part layout was never supported; `.h5part` files are detected and warned about |
| Plotting / animation | Working | Publication theme: no titles, no minor ticks, Computer Modern fonts, dashed grey low-opacity grid on line plots; presentation knobs config-driven via `[visualization.style]` (`PlotStyle`); escaper suite (cumulative mass loss, velocity classes, escape anisotropy) and SSE-quantity plots (mass segregation, t/T_MS evolutionary clock, core-mass growth); binary suite (pair counts with the Heggie hard/soft split, binary fraction, `a`–`e` diagrams, period histograms); existing figures are never overwritten (safesave-style `#1`, `#2`, … backups) |
| Load latency | Working | PrecompileTools workload over the configuration, diagnostics-reader and merger-IC paths; the first `load_config` of a session no longer compiles |
| Tests | Passing | 1371/1371 as of this commit (plus 33 engine-dependent tests behind `NBODY6_BINARY_TESTS=1` and a GPU-gated set behind `NBODY6_GPU_TESTS=1`), incl. physics validation, adversarial external-input tests, telemetry/thread-control, per-cluster structure, restart, and tidal-field tests |

## Testing

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The default suite needs no Fortran binary. The engine-dependent tests (build, launch, restart, tidal field, live telemetry) run when `NBODY6_BINARY_TESTS=1` is set; without `NBODY6_BACKEND_ROOT` they clone and build the backend in a temporary directory (what the weekly `Backend` workflow does), with it they use the build under `<root>/backend/Nbody6PPGPU-beijing`:

```bash
NBODY6_BINARY_TESTS=1 NBODY6_BACKEND_ROOT=$PWD julia --project=. -e 'using Pkg; Pkg.test()'
```

The suite covers config round-trips, all I/O readers (synthetic binaries plus real fork output fixtures under `test/fixtures/`), plotting/animation smoke tests, external post-processing (including adversarial malformed inputs), and the merger IC generator. Physics-validation tests check King concentration c(W0) against published values, the Plummer half-mass relation r_hm = 1.305 a, the Kroupa mean mass, virial equilibrium Q = T/|W| = 0.5 after `virialise!`, and the Keplerian orbital energy of generated two-cluster orbits.

For interactive work against the test dependencies, activate the test environment with [`TestEnv.jl`](https://github.com/JuliaTesting/TestEnv.jl): `julia --project=. -e 'using TestEnv; TestEnv.activate()'` in a REPL session.

## Verification Suite

```bash
julia scripts/run_verif_suite.jl
```

Runs three end-to-end targets through `run_pipeline` (requires a built binary in `backend/`):

1. `single` — plain N=5000 run (`N5k_medium.inp`)
2. `triorbit` — bound three-cluster Lagrange triangle (`verif_triorbit.toml`)
3. `3d5cluster` — five clusters distributed in 3D (`verif_3d5cluster.toml`)

## Documentation

Documenter.jl docs live under `docs/` (Manual, Input Files, Cluster Mergers, API Reference):

```bash
julia docs/make.jl   # builds to docs/build/
```
