# Nbody6Setup.jl

A Julia package that automates the full lifecycle of [Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing) star-cluster simulations: install/build of the Fortran code, multi-cluster merger initial-condition generation, simulation execution, post-processing of all standard output files, and publication-quality visualization. Every phase is driven by a single TOML configuration and orchestrated through one entry point, `run_pipeline`.

## Project Structure

```
Nbody6Setup/
├── README.md
├── LICENSE                          # MIT
├── Project.toml                     # Package metadata & dependencies
├── Manifest.toml                    # Version-controlled — exact dependency versions
├── config.toml                      # Main pipeline configuration (edit this)
├── src/
│   ├── Nbody6Setup.jl               # Module root; run_pipeline orchestrator; exports
│   ├── types.jl                     # Config structs (incl. PlotStyle), Snapshot, records, UnitScaling
│   ├── config.jl                    # TOML loader/serialiser for Nbody6Config
│   ├── util.jl                      # safesave-style never-overwrite backups
│   ├── platform.jl                  # OS/CUDA/dependency detection
│   ├── install.jl                   # Clone, configure, HDF5 Makefile patch, build
│   ├── run.jl                       # Simulation launcher: run dirs, launch script, live monitoring
│   ├── external.jl                  # scan_output + postprocess_external for arbitrary output dirs
│   ├── io/
│   │   ├── io.jl                    # I/O submodule includes
│   │   ├── fortran_binary.jl        # Fortran unformatted record reader
│   │   ├── conf3.jl                 # conf.3 snapshot reader (standard + extended layout)
│   │   ├── diagnostics.jl           # out1000 parser: ADJUST lines + PHYSICAL SCALING
│   │   ├── lagr.jl                  # lagr.7 Lagrangian radii reader
│   │   ├── esc.jl                   # esc.11 escaper reader (fork's real column format)
│   │   └── stellar_evolution.jl     # sev.83_* single-star evolution reader
│   ├── ic/
│   │   ├── ic.jl                    # run_merger_pipeline / generate_merger_ic entry points
│   │   ├── config.jl                # Merger TOML parser: ClusterSpec, OrbitSpec, seed
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
│       ├── merger.jl                # Merger IC overview figures
│       └── animation.jl             # GIF animations (cluster, HR, Lagrangian)
├── scripts/
│   ├── run_setup.jl                 # CLI wrapper: load config → run_pipeline
│   └── run_verif_suite.jl           # Three-target verification suite
├── test/
│   ├── runtests.jl                  # Full unit + physics-validation suite
│   ├── test_external_adversarial_inner.jl  # Adversarial external post-processing tests
│   └── fixtures/                    # Real Nbody6++ output excerpts (esc.11, lagr.7, out1000, sev.83_0)
├── benchmark/
│   └── benchmarks.jl                # BenchmarkTools suite (kept out of tests)
├── docs/
│   ├── make.jl                      # Documenter.jl build script
│   └── src/                         # index, manual, input_files, multi_cluster_mergers, api
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
│   ├── merger_2cluster_medium.toml  # Two King clusters on eccentric Kepler orbit
│   ├── merger_3cluster_small.toml   # 3-cluster equilateral triangle, small
│   ├── merger_3cluster_medium.toml  # 3-cluster equilateral triangle, medium
│   ├── merger_5cluster_small.toml   # 5-cluster pentagon, small
│   ├── merger_5cluster_medium.toml  # 5-cluster pentagon, medium
│   ├── merger_27cluster_cubic.toml  # 27 clusters on a 3×3×3 cubic grid
│   ├── verif_triorbit.toml          # Bound Lagrange-triangle verification target
│   └── verif_3d5cluster.toml        # 5 clusters distributed in 3D (projection/COM verification)
├── backend/                         # (gitignored) cloned Nbody6PPGPU-beijing source + build
└── runs/                            # (gitignored) per-run output/, plots/, frozen config.toml
```

## Environment Setup

Requires Julia ≥ 1.10.

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Manifest.toml` is version-controlled, so `instantiate` reproduces the exact dependency set on any machine. Building the Fortran backend additionally needs `git`, `gfortran`/`make`, and optionally HDF5 and CUDA (auto-detected; see `[build]` in `config.toml`).

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
using Nbody6Setup
cfg = load_config("config.toml")
results = run_pipeline(cfg)   # Dict with :snapshots, :diagnostics, :lagr, :escapers, :stellar_evo
```

Or via the CLI wrapper:

```bash
julia --project=. scripts/run_setup.jl [config.toml]
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
| Install / build | Working | Clone, `configure`, HDF5 Makefile patch, parallel make; CUDA path auto-detection |
| Merger IC generator | Working | Plummer + King samplers (King c(W0) validated against published concentrations); Kroupa (2001) IMF; Kepler two-body and explicit N-cluster orbit modes; Jacobi truncation; seeded reproducibility — the TOML `seed` drives the sampler RNG and propagates to Nbody6's `NRAND` |
| Simulation runner | Working | Launch script, live stdout monitoring, run summary; merger runs execute inside the IC output dir so `dat.10` is found |
| I/O readers | Working | `conf.3` (standard + extended), `out1000` diagnostics (ADJUST + physical scaling; virial ratio Q = T/\|W\|, equilibrium at 0.5), `lagr.7`, and `esc.11` / `sev.83_*` in the fork's real formats; `STELLAR_TYPE_LABELS` follow the Hurley convention (13 = NS, 14 = BH); `UnitScaling.zmbar` is the total-mass scale factor M*, not the mean stellar mass. HDF5 reader removed — the fork's KZ(46) H5Part layout was never supported; `.h5part` files are detected and warned about |
| Plotting / animation | Working | Publication theme: no titles, no minor ticks, Computer Modern fonts, dashed grey low-opacity grid on line plots; presentation knobs config-driven via `[visualization.style]` (`PlotStyle`); existing figures are never overwritten (safesave-style `#1`, `#2`, … backups) |
| Tests | Passing | 364/364 as of this commit, incl. physics validation and adversarial external-input tests |

## Testing

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The suite covers config round-trips, all I/O readers (synthetic binaries plus real fork output fixtures under `test/fixtures/`), plotting/animation smoke tests, external post-processing (including adversarial malformed inputs), and the merger IC generator. Physics-validation tests check King concentration c(W0) against published values, the Plummer half-mass relation r_hm = 1.305 a, the Kroupa mean mass, virial equilibrium Q = T/|W| = 0.5 after `virialise!`, and the Keplerian orbital energy of generated two-cluster orbits.

## Verification Suite

```bash
julia --project=. scripts/run_verif_suite.jl
```

Runs three end-to-end targets through `run_pipeline` (requires a built binary in `backend/`):

1. `single` — plain N=5000 run (`N5k_medium.inp`)
2. `triorbit` — bound three-cluster Lagrange triangle (`verif_triorbit.toml`)
3. `3d5cluster` — five clusters distributed in 3D (`verif_3d5cluster.toml`)

## Documentation

Documenter.jl docs live under `docs/` (Manual, Input Files, Cluster Mergers, API Reference):

```bash
julia --project=docs docs/make.jl   # builds to docs/build/
```
