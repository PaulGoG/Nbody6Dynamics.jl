# Changelog

Notable changes to Nbody6Dynamics.jl, by version. The package follows
semantic versioning; before 1.0 a minor version may change the API, and the
entries say so.

## [Unreleased]

### Changed

- The launch script owns the stop request. `_launch.sh` stays the parent of
  the engine: it starts the engine in a subshell that ignores SIGTERM, turns
  a SIGTERM to itself into the `STOP` file, waits, and exits the way the
  engine did (its exit status, or the signal re-raised). A scheduler's
  job-wide SIGTERM therefore stops the engine with a dump; a termination by
  the package is a stop request followed, 15 s later, by a SIGKILL of the
  process group. Segment records gain `engine_pid` next to `pid` (the
  script); `output/engine.pid` holds it while the engine runs. An engine
  killed at its wall budget is recorded with exit status −9 (was −15); an
  operator interrupt closes the record as `stopped` with the dump when the
  engine honoured the request (`stop_requested = "interrupt"`), else
  `killed`. A SIGTERM sent to the engine process by hand no longer ends it;
  a SIGKILL of the launch script alone orphans the engine — kill
  `engine_pid`, or the process group.

### Added

- Slurm job layer (`[cluster]`, `ClusterConfig`). `submit_run` and
  `scripts/run_cluster.jl` turn a run into a chain of queue jobs, one per
  engine segment: the job script is generated into `<run_dir>/cluster/`
  from the section's keys (partition, account, qos, counts, memory, gres,
  time limit, extra options, prologue lines) with `--signal=B:TERM` before
  the limit, and `cluster_segment` runs a segment inside the job (the
  pipeline on the first, a resume afterwards), resubmitting its successor
  (`chain = "resubmit"`) or relying on a chain submitted at once
  (`chain = "dependency"`). Job events are appended to `cluster/jobs.csv`
  and `<runs_dir>/cluster_jobs.csv` (`read_job_ledger`). A configuration
  whose segment does not fit the time limit, or whose time limit exceeds
  the partition's, is refused at load. `--dry-run` renders the script and
  the submission commands.

- Cost model (`[cost]`, `CostConfig`). `estimate_cost` derives the wall
  time, the number of wall-budget segments and the disk of a run from a
  table of seconds per N-body time unit against particle number
  (`cost_rate`; built-in rows measured on one H200 NVL device, replaceable
  through `table`), with a factor for two-cluster configurations. Every
  engine launch logs the estimate, writes it into the segment record
  (`cost_*` keys) and warns or refuses against `budget_hours`.

- Engines per device in sweeps. `[sweep] devices` and `engines_per_device`
  assign a CUDA device to each point when it starts (the least-loaded one)
  and admit at most that many engines per device, in addition to
  `concurrency`; the point's frozen `config.toml`, the index and the
  summary (`device` column) record the assignment.

- Surface-density rendering of particle projections. `plot_snapshot`,
  `plot_snapshot_evolution` and `animate_cluster` draw the projected surface
  mass density (log₁₀ Σ in M☉ pc⁻², one colour scale per figure or
  animation) instead of the mass-coloured scatter, whose markers overlap into
  a flat disc above a few thousand particles. `[visualization.style]` gains
  `snapshot_render` (`"auto"`, `"scatter"`, `"density"`), `density_min_n`,
  `density_bins`, `density_decades` and `density_mass_frac`; the limits of a
  density panel enclose that mass fraction, so escapers no longer set the
  scale.

- Energy bookkeeping of the engine in the diagnostics reader. `AdjustRecord`
  gains `de_abs` (the engine's `DELTA`, the energy change of one adjustment
  interval) and `detot` (`DETOT`, its running sum, an absolute energy);
  `cumulative_energy_error` returns the cumulative relative and absolute
  series per adjustment and `read_energy_totals` the totals of the `END RUN`
  line. Code that constructs `AdjustRecord` positionally passes two more
  arguments.

- Segment status in `RUN_INFO.toml`. Every engine launch is recorded when it
  starts (`status = "running"`, launch time, host, process id, package
  commit, start time in N-body units, stdout offset) and closed at exit with
  one of `completed`, `stopped`, `halted`, `watchdog`, `killed`, `failed`;
  `run.status` repeats the last segment's. A halt by the engine's energy
  check and a stop on request are told apart from a crash.

- Wall budget and stop requests. `simulation.wall_budget` asks the engine to
  stop `stop_margin` seconds before the budget expires, through the `STOP`
  file it polls for; it writes a restart dump and exits, and the segment is
  recorded as `stopped`. A driver that receives SIGTERM makes the same
  request for its running engines, waits for them and closes their records,
  so a job ended by a scheduler leaves a dump and no engine behind.

- Build-limit check before launch. The particle number, the number of
  primordial binaries and the neighbour number of an input are compared
  with `NMAX`, `KMAX` and `LMAX` of the configured engine tree
  (`include/params.h`), under the conditions the engine itself applies at
  start-up; an input that does not fit is refused with the limit and the
  remedy. `BUILD_INFO.toml` records the limits as `engine_limits`.

- Checkpoint mode of the generated merger input. `checkpoint = true` in
  `[merger.nbody6]` makes the engine write a restart dump at every
  adjustment (`KZ(2) = 1`); the generator logs the size of one dump and of
  the whole series.

- Retention of periodic restart dumps. `simulation.checkpoint_keep = k`
  keeps the newest `k` dumps of a checkpointed run and deletes older ones
  while the engine runs; `0`, the default, keeps all.

- `resume_run(run_dir)` continues an interrupted run to the end time of its
  input and can be called again safely. It starts from the last complete
  dump at the exact time the engine reported, and after a kill it moves
  what the dead segment wrote beyond that dump (stdout tail, later records
  of the series the package reads, later snapshots and dumps) to
  `output/discarded/`, so the series stay single-valued; nothing is deleted.

- `resume_pipeline(run_dir)` and `scripts/run_resume.jl` continue an
  interrupted run from the configuration frozen in its directory and, once
  it has reached its end time, perform post-processing and figures. The
  call can be repeated, one per queue segment.

- Single-cluster mode of the generator. `orbit_mode = "single"` with
  `n_clusters = 1` samples one cluster, places it at the origin at rest
  without truncation and writes the same files as a merger; `position`,
  `velocity` and `[merger.orbit]` are refused or ignored. The isolated
  controls of a sweep and `gpu/single_600k.toml` use it.

- Single-cluster science cases as generator configurations under
  `input_files/clusters/`: a globular cluster with a black-hole subsystem,
  runaway collisions in a dense young cluster, a young massive cluster with
  a binary fraction of 0.5, and an open cluster of Hyades-like mass with its escapers
  kept. Each has been run through its first time units on the CPU build.

- `merger.nbody6.preset`: `"engine"` (the rules of the engine's example
  inputs, the default) or `"campaign"` (the calibration baseline: η_I 0.01,
  η_R 0.02, NNBOPT 100, RMIN three times the 90° deflection scale of the
  smallest member, DTMIN 0.4 of the Aarseth close-encounter step, no
  automatic re-derivation of the integration parameters). A zero of
  `etai`, `etar`, `etau`, `nnbopt`, `rs0`, `rmin` or `dtmin` now means the
  preset's value or rule; explicit keys override it.

### Changed

- The defaults of `Nbody6ParameterSpec` for `etai`, `etar` and `etau` are
  `0.0` (resolved from the preset at generation time; the engine preset
  gives the former 0.02, 0.02, 0.1). A specification passed directly to
  `generate_merger_inp` must carry nonzero values, as it already must for
  the derived keys.

- The shipped engine inputs `N1k_quick` to `N100k_production` write a
  restart dump at every adjustment (`NCOMM = 1`, was 10: with the larger
  value the dump of an end, a stop or a halt was written one time in ten).
  The 1k and 5k inputs are isolated (`KZ(14) = 0`): they requested a
  point-mass field without its namelist, which the engine read as a galaxy
  of zero mass. The 25k and 100k inputs generate their primordial binaries
  with the engine's `BINPOP` (`KZ(8) = 1`): `KZ(8) = 2` with internally
  generated initial conditions declared the first `2·NBIN0` bodies to be
  pairs. The 10k input has no binaries and says so (`KZ(8) = 0`).

- The five science inputs `gc_bh_subsystem`, `imbh_runaway`, `pop3_cluster`,
  `tidal_tails` and `young_massive_binaries` are withdrawn: they generated
  Plummer models where their headers said King, declared binaries the
  engine never generated, quoted spans in Myr wrong by factors of three to
  fifty, and the Population III case cannot start (the engine refuses
  `ZMET` below 10⁻⁴). The configurations under `input_files/clusters/`
  replace four of them.

- With the default `snapshot_render = "auto"`, snapshots of 20 000 particles
  or more are drawn as surface density; `"scatter"` restores the former
  figures.

- The start-up watchdog and the completion monitor of a restart read only
  the part of the stdout capture written by that segment; an `END RUN` or
  an adjustment line of an earlier segment no longer counts.

- `merger.nbody6.ncomm` defaults to 1 (was 10). The engine counts every
  dump call and writes one in `ncomm`, so a larger value could skip the
  dump of a stop request or of an energy halt; inputs without periodic
  dumps are unaffected.

- A restart no longer keeps backups of `output/comm.1`: it is a working
  copy of a dump that stays under its own name, and it is removed when the
  engine exits. The restart input carries the time increment unrounded.

- `run_pipeline` ends after a simulation phase that stopped on request
  (`wall_budget`, or a signal): post-processing and figures are left to
  the completed run instead of being made from every segment.

### Fixed

- Snapshots no longer count regularised pairs twice. The engine writes
  NTOT = N + NPAIRS entries to `conf.3`, the N stars followed by one
  centre-of-mass entry per KS pair carrying the pair's mass, and
  `read_conf3` returned all of them as stars. Every sum over a snapshot
  therefore included each regularised pair's mass twice: up to 23 % excess
  mass in a 4920-star merger with 819 primordial pairs, 0.2 % in a
  48 000-body merger with one dynamically formed pair. The reader now drops the
  centre-of-mass entries, and the centre-of-mass body of an active chain
  regularisation (NAME = 0), which the engine leaves beside the chain
  members. `Snapshot` gains `npairs`, the components of pair k being stars
  2k − 1 and 2k (`pair_indices`), and `pair_name`, the names of the pairs'
  centre-of-mass entries (negative for a hierarchical merger), and validates
  its arrays on construction. Remnant diagnostics, per-cluster structure, profiles, the
  binary hardness scale and the snapshot figures of runs with regularised
  pairs change accordingly.
- `animate_cluster` fits its canvas to the square axis, as `plot_snapshot`
  does: a blank strip about a quarter of the frame wide separated the axis
  from the colourbar. Under adaptive zoom the widest tick labels of the
  animation are reserved, so the axis no longer shifts between frames.
- The binary-population figures no longer run an ungated O(N²) pair sum:
  their hard/soft scale selected the bound systems of every snapshot
  whatever `pair_sum_max_n`, forty minutes of the post-processing of the
  575 080-body probe. Above the limit the scale is now taken over every
  system, with a warning; `binary_scales` takes the limit as a keyword.
- The telemetry figure's legend is at most two columns wide. With a device
  present it carries six entries with their summary values, and three
  columns of those overran the canvas on every campaign run.
- A step that removes more than 1000 escapers no longer ends a run. The
  engine writes the escaper names of one adjustment interval with a format
  holding a thousand of them and aborted on the next (`Expected REAL … in
  formatted transfer`, exit status 2), as a 27-cluster cold collapse did at
  1442 escapers. The install phase now corrects that format in the engine
  tree before the build (`deps/engine/amendments.toml`) and records the
  amendments in effect as `source_amendments` in `BUILD_INFO.toml`; an
  engine built earlier needs a rebuild.
- GPU telemetry follows `simulation.gpu_list`. The sampler queried every
  device of the host, so a run on one of two cards recorded half its
  utilisation and the power of both; it now names the run's devices to
  `nvidia-smi`. With a device list the launch script also exports
  `CUDA_DEVICE_ORDER=PCI_BUS_ID`, so that an index means the same device
  for the engine and for `nvidia-smi`.
- The input-file reference described `TCRTP0` as a wall-clock limit in
  seconds. It is the engine's termination time in Myr; the run ends at the
  earlier of `TCRIT` and `TCRTP0`.

## [0.4.0] — 2026-09-26

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
- `bench/fp32_peak_probe.jl`: the measured FP32 peak of the device with the
  peak probe of GPUDiagnostics.jl, written as `fp32_peak_<hostname>.toml`
  for the fraction-of-peak panel of `gpu_cells_figures.jl`. The script keeps
  its own shared environment (`@nb6_fp32_peak`) and installs
  GPUDiagnostics, KernelAbstractions and CUDA there on first use.
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
- `bench/gpu_cells_figures.jl`: the scaling figure (strong, weak and
  host-thread scaling of the two-device machine) and the cross-hardware
  figure (wall time per card with the regular-force share as a tick; kernel
  rate as a fraction of the FP32 peak) from the cell CSVs, a console table,
  and a `--synthetic` rehearsal mode.

### Changed

- `input_files/` is laid out by purpose: `engine/` (the ten single-cluster
  engine inputs), `mergers/` (the seven merger configurations),
  `verification/` (the two targets of the verification suite), `sweeps/`
  (the demonstration sweep), next to the unchanged `showcase/` and `gpu/`.
  File names are unchanged; `example_input` resolves a bare name across the
  folders, so `example_input("N1k_quick.inp")` and the other names of 0.3.0
  keep working, while a name that ships in more than one place must be given
  with its folder. The shipped pipeline and sweep TOMLs point at the new
  paths; a project's own configuration that named a shipped input by its old
  flat path must add the folder. `input_files/README.md` maps the folders to
  the documentation.
- The root `config.toml` runs the two-cluster demo merger
  (`mergers/merger_demo_small.toml`, seconds) instead of the 27-cluster
  cubic grid, a cold-collapse stress test outside the engine's regime that
  is kept under `mergers/` and described as such.
- The five science engine inputs carry a short factual header; their
  science cases, expected phenomena, cost estimates and references (with
  DOIs) moved to the Input Files page, which is organised by folder. A new
  documentation page, Showcase Cases, describes the five showcase
  configurations, their reference runs of 2026-09-11 (wall times, energy
  errors, run identifiers) and what a run produces; the flagship case is
  identified with the validation merger of `gpu/merger_50k.toml`, whose
  integration is timed on every validated host. The manual's validated
  hardware table carries the H200 NVL node. The default of
  `simulation.input_file` names a file that ships
  (`input_files/engine/N1k_quick.inp`); it named a path that did not exist.
- The package's author line in `Project.toml` carries the name only, as the
  licence does.

### Fixed

- A retried validation stage starts from a clean run directory: the
  crashed attempt's `runs/<run_id>` is kept as
  `runs/<run_id>.attempt<k>.signal<N>` beside its log, where the rerun used
  to append a second engine segment to the same record and accumulate the
  elapsed time. This makes `--retry-stages` usable for the pipeline stages,
  which the host with the JIT crashes needs.
- The start-up watchdog and the completion monitor terminate the engine
  before they report it, and `_terminate` sends SIGKILL before it reports the
  escalation. The report came first, and a pipeline whose driver had died
  wrote it to a closed pipe: the watchdog task ended on `EPIPE` with the
  kill never sent, and the engine it should have ended ran for thirteen
  hours beside its rerun (2026-09-25). The console sink of a run's logger
  now absorbs a failed write as well (`_ResilientLogger`), so a pipeline
  whose launcher is gone still writes its run record and its file log.
- An operator interrupt reaches the engine. The validation driver and every
  entry script that starts a stage or an engine set
  `Base.exit_on_sigint(false)`, so SIGINT raises `InterruptException` instead
  of ending the script at once — the forwarding of `dbf7274` never ran,
  because a non-interactive Julia exits on SIGINT before any handler. Inside
  a pipeline the helper tasks (telemetry sampler, start-up watchdog,
  completion monitor) hand an interrupt that lands in them to the task
  waiting on the engine, whose handler terminates it: SIGINT is delivered to
  whichever task is current, and under a long wait that is usually a helper.
  A SIGTERM or SIGKILL of a pipeline still leaves its engine (own session)
  running; the manual names the process sweep that ends the whole tree.
- The run record's `package_commit` and `backend_commit` are captured when
  the engine is launched, not when the record is written: a `git pull` on
  the host during a long run used to stamp the record with a commit that
  had not run.
- `startup_timeout` of the four 6 × 10⁵-body probe pipelines under
  `input_files/gpu/` raised from 3600 s to 14 400 s: the AVX build needs more
  than an hour to its first adjustment on a slow host (i7-10750H, EPYC 7551P),
  and the watchdog killed the `cpu_merger_600k` probe there at 3600 s.
- `gpu_cells_figures.jl`: FP32 peak records are matched by host *and*
  device, so two cards in hosts of the same name (one workstation name, two
  machines) each take their own measured peak instead of one overwriting the
  other in the lookup; the host comparison ignores the domain part, which
  `gethostname()` includes on some machines and the cards file does not.
- Stage processes of the validation driver and the engine itself start in
  their own session (`detach`). A terminal closing above a `nohup`-ed
  driver used to kill the running stage and its engine: `nohup` shields
  the driver alone, a child spawned by Julia has libuv's default signal
  dispositions. Two campaigns of 2026-09-23/24 lost their running stage on
  every host that way. An operator interrupt (Ctrl-C) is forwarded
  explicitly: SIGINT to a stage, SIGTERM then SIGKILL to the engine.
- `bench/gpu_cells_figures.jl`: the logarithmic floor of the hardware
  figure's wall-time panel is derived from the regular-force times as well
  as the wall times, so the regular-force tick of a bar is never clipped
  below the axis (it was, for every card whose regular-force time fell under
  the first decade of the wall times).

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
