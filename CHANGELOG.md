# Changelog

Notable changes to Nbody6Dynamics.jl, by version. The package follows
semantic versioning; before 1.0 a minor version may change the API, and the
entries say so.

## [Unreleased]

### Added

- `[merger] virial_max_n` (default 200 000): the largest N, per cluster and
  for the combined system, for which the initial-condition generator
  evaluates the exact O(N²) potential. A cluster above it is refused when the
  configuration is loaded, naming the key, instead of failing inside the
  sampler; large initial conditions raise it deliberately.
- `[postprocess] pair_sum_max_n` (default 100 000): the largest snapshot
  particle count for which the O(N²) pair-sum diagnostics run (remnant
  diagnostics, per-cluster virial ratio and structure, bound-member
  profiles). Above it they are skipped with a warning naming the key; the
  readers and the other figures are unaffected. `generate_plots` takes the
  same limit as a keyword.
- The probes above 5 × 10⁵ bodies under `input_files/gpu/`:
  `merger_600k.toml` (two King clusters of 300 000 stars on the orbit of
  `merger_50k.toml`) and `single_600k.toml` (one King cluster of 600 000 at
  rest), each for two N-body time units, with the CUDA and AVX pipeline
  configurations that run them on the trees the recipe builds. The
  benchmark case raises its virialisation limit to the engine's `b1m`
  capacity.
- `run_gpu_validation`: the four probe pipelines as stages
  (`gpu_merger_600k`, `gpu_single_600k`, `cpu_merger_600k`,
  `cpu_single_600k`, opt-in) and `stop_on_failure` (`--stop-on-failure`),
  which skips every stage after a failed one and runs nothing needing `nvcc`
  when the host-compiler probe did not pass, so a workstation campaign is one
  gated invocation. The probe and benchmark stages inherit
  `JULIA_NUM_THREADS` or get `auto`.
- `bench/gpu_cells.jl`: benchmark cells on the device, a grid of (cell, N,
  host threads, `GPU_LIST`, MPI ranks) points run with post-processing off
  and read back into one CSV row per run (wall time, backend timing table,
  kernel throughput, GPU utilisation, power and memory, host RSS, CPU
  efficiency, machine identity, run ID). `bench/merger_case.jl` gains the
  single-cluster cell (`single_toml`) next to the two-cluster case and
  `cell_toml` to select either.

## [0.3.0] — 2026-09-22

First public release. Nbody6Dynamics.jl automates the lifecycle of
Nbody6++GPU star-cluster simulations from one TOML configuration: it clones
and builds the engine (CPU, GPU and MPI variants, CUDA architectures
resolved from the visible devices), generates multi-cluster merger initial
conditions, runs the integration with live monitoring, watchdogs, restarts
from the engine's dumps and hardware telemetry, reads every standard output
file into typed records, derives the diagnostics of the run, and draws
publication figures and animations through a Makie extension. Every run
directory carries its frozen configuration, the resolved environment, the
package and engine commits and a hardware fingerprint.

### What the release contains

- **Configuration.** One `config.toml` drives the pipeline; the merger
  generator and the parameter sweeps have their own TOML schemas. The
  parsers reject unknown keys and wrong types and raise an `ArgumentError`
  naming `section.key`; every default is stated once, in the configuration
  structs. Relative paths resolve against the directory of the file
  (`cfg.config_dir`), which is where the engine tree and the run directories
  are created, so an installation by URL behaves like a checkout;
  `example_input(name)` reaches the shipped inputs either way.
- **Engine build.** Clone at the validated commit, configure, HDF5 patch,
  parallel make; CUDA toolkit detection, `nvcc` host-compiler probe,
  `-gencode` entries per architecture, `BUILD_INFO.toml` next to the binary.
- **Merger initial conditions.** King and Plummer samplers, Kroupa IMF in
  natural, rescaled and equal-mass forms, primordial binaries, Kepler or
  explicit orbits, Jacobi truncation, external tidal fields, engine
  parameters derived from the smallest member cluster, output intervals as
  dyadic rationals the engine's digit counter can terminate on, seeds
  recorded for reproduction.
- **Runs.** Isolated `runs/<id>/` directories, a shell-safe launch script,
  OpenMP and `GPU_LIST` control, start-up and completion watchdogs
  (SIGTERM, then SIGKILL), restarts from COMMON dumps with per-segment
  bookkeeping, CPU, memory and GPU telemetry, records written atomically
  and never overwritten, a completion marker that requires the engine's
  `END RUN`.
- **Readers and diagnostics.** `conf.3`, `out1000`, `lagr.7`, `esc.11`,
  `sev.83_*` and `bev.82_*` in the engine's current and documented layouts,
  with skipped lines reported; per-cluster structure and radial profiles,
  binary population with the Heggie split, stellar classes and a class
  census, remnant diagnostics (core radius, rotation, mass segregation,
  coalescence time), written as CSV data products with or without a figure
  backend.
- **Sweeps and ensembles.** Cartesian grids over merger keys times seeds,
  run as concurrent worker processes with an index kept current, seeded
  ensembles with median and central intervals, isolated single-cluster
  controls per point.
- **Figures.** The publication layout: canvases in layout units exported at
  a printed width, Computer Modern type, one formatter for every number in
  figure text, a run-level HR layout with a census strip, montage,
  telemetry, sweep, ensemble and control figures, GIF animations. The figure
  routines are declared by the package and implemented by the
  `Nbody6DynamicsMakieExt` extension; without a backend a call is a
  `MethodError` whose message names the remedy, and the orchestration entry
  points refuse before doing any work.
- **Validation of a CUDA host.** `run_gpu_validation` runs the GPU-gated
  suite, the GPU and CPU pipelines and the scaling benchmark as logged
  stages with a retry policy, identifies the machine by hostname, hardware
  and GPU UUID, and judges each pipeline stage by the run directory it
  assigned.
- **Tooling.** A component-split test suite with static QA (Aqua, JET,
  ExplicitImports) and engine- and GPU-gated sets, a `Documenter` site with
  a Literate walkthrough and DOI-keyed references, benchmarks, and a
  `CITATION.cff`.

### Compatibility

Julia 1.13 or later; CairoMakie 0.13 or 0.15 (Makie 0.22 or 0.24) for the
figure extension, both series tested. The engine is Nbody6PPGPU-beijing at
its v2026.07 commit `618d7a4`, cloned at build time and not redistributed.
The strict parser refuses a frozen `config.toml` of an earlier run that
carries a key the schema has since dropped (`visualization.figsize`); remove
the key to load or restart such a run.

Earlier versions were private development releases; their history is in the
repository.
