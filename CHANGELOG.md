# Changelog

All notable changes to Nbody6Dynamics.jl. Follows semantic versioning;
pre-1.0 minor versions may break APIs (private project, no-compat policy).

## [Unreleased]

### Added
- `publication_theme()`: the publication theme as a value, for `with_theme`
  scoping.
- `run_gpu_validation` and `scripts/run_gpu_validation.jl`: the acceptance
  sequence of a CUDA host (GPU-gated suite, GPU and CPU pipelines, scaling
  benchmark) as logged stages, with the host record (including the host
  compilers found for `-ccbin` and the verdict of the `nvcc` host-compiler
  probe), every log, the benchmark results and a summary under
  `runs/gpu_validation_<host>_<stamp>/`.
- GPU builds probe `nvcc` on a trivial kernel and add
  `-allow-unsupported-compiler` themselves, with a warning, when the toolkit
  rejects the host compiler; the effective options are recorded in
  `BUILD_INFO.toml`. When no host compiler works, the error quotes the
  `nvcc` output of every attempt.
- `[build] nvcc_flags`: extra `nvcc` options for the GPU build (host-compiler
  overrides such as `-allow-unsupported-compiler`), recorded in
  `BUILD_INFO.toml`.
- `deps/cuda/`: `helper_cuda.h` and `helper_string.h` from NVIDIA's
  cuda-samples (tag v13.0, BSD-3 licence alongside).
- `.mailmap` folding the earlier author identities into one.
- Activation scripts for every environment: `activate.jl` (package),
  `docs/activate.jl` and `bench/activate.jl` (auxiliary environments,
  package developed by a relative path). The scripts under `scripts/`,
  `docs/` and `bench/` include their environment's script, so none needs a
  `--project` flag; `julia activate.jl` bootstraps a new machine.
- `CITATION.cff`.
- `bench/Manifest.toml` is tracked now that the develop path is relative.
- CUDA-host recipe: `input_files/gpu/gpu_pipeline.toml` (CUDA build into
  its own tree, architectures from `nvidia-smi`, device 0, eight host
  threads), `cpu_pipeline.toml` (the AVX reference build) and
  `merger_50k.toml` (two King clusters of 25 000 stars, `qe = 0.01`), with
  the command sequence in the manual; `bench/gpu_scaling.jl` takes the CPU
  and GPU trees from `NBODY6_CPU_BACKEND`/`NBODY6_GPU_BACKEND`.
- Documentation site (P2): literature citations through DocumenterCitations
  (`docs/src/references.bib`, `@cite` markers in the pages, a References
  page in author–year style, replacing the per-page plain lists); a
  Literate walkthrough (`docs/src/walkthrough.jl`, executed at build time)
  from a merger TOML to the generated initial conditions, the engine input
  file and the IC figures; `deploydocs` for the GitHub Pages deployment,
  which publishes once the workflows run again. Docstring templates through
  DocStringExtensions were considered and not adopted: every public
  docstring already opens with its signature, and a template would repeat
  it.
- Precompile workload (P4, PrecompileTools): configuration load and
  save, the diagnostics reader, `engine_interval`, run-ID and elapsed-time
  helpers, and a small merger initial-condition generation run at package
  precompilation, so the first call of each in a session is compiled
  already; plotting is excluded (the CairoMakie precompile is its own).
- Live sparklines (F8): `[simulation] live_diagnostics` with `live_interval`
  prints, through the interactive monitor, UnicodePlots sparklines of the
  virial ratio and `log10 |ΔE/E|` against time from the ADJUST records so
  far; stderr only, never the log file. New dependency `UnicodePlots`.
- Threaded snapshot reading (F10): `read_all_conf3(dir; threaded)` reads
  the files in chunks on `Threads.@spawn` tasks and assembles them in time
  order (default on when Julia has more than one thread; identical result to
  the serial read); the ordered reader `_read_ordered` is generic. Sweep
  worker processes inherit the driver's thread count, which
  `Base.julia_cmd` does not carry.
- Run telemetry figure: `read_telemetry`/`read_run_telemetry` read the
  sampler's `telemetry*.csv` (segments of restarted runs concatenated with
  cumulative offsets) and `plot_telemetry` draws cores busy with the host
  load average on a twin axis, resident memory (RSS, high-water mark) and,
  when sampled, GPU utilisation, as stacked panels with the means and the
  peak RSS as legend entries; `generate_plots` draws it for every run
  directory holding telemetry (`telemetry` in the plot suite).
- GPU build target (F12). `[build] cuda_arch` lists the CUDA architectures
  the kernels are compiled for (`sm_90` Hopper, `sm_120` consumer
  Blackwell, …); empty means the compute capabilities `nvidia-smi` reports,
  or the `nvcc` default when no device is visible. The upstream configure
  script emits no architecture flag, so the build passes the resolved
  `-gencode` entries (native code per architecture, PTX for the highest) to
  `make` as a `CUFLAGS` override, after checking them against the
  toolkit's own `nvcc --list-gpu-arch` so a toolkit older than the device
  fails the install phase with the release and the missing architecture
  named. `[simulation] gpu_list` names the devices
  the engine may use (its `GPU_LIST` variable, at most four per process);
  the launch script exports it. The binary is chosen by its suffix tags
  (`.gpu`, `.mpi`) so CPU and GPU variants can coexist in one build tree.
  `BUILD_INFO.toml` next to the binary records date, host, backend commit,
  configure arguments, switches, CUDA path, architectures and `nvcc`
  release; each run copies it and merges it into `RUN_INFO.toml` as
  `[build]`. The run summary also records `run.gpu_list`, the devices the
  GPU library reported at initialisation (`run.gpu_devices`) and the
  kernel label of the throughput profile (`GPU Reg.F` versus `AVX Reg.F`);
  the hardware fingerprint adds the compute capability per device.
  `bench/gpu_scaling.jl` measures the GPU speed-up over the CPU binary at
  equal N and threads for a list of `GPU_LIST` values; the two-cluster
  case moved to `bench/merger_case.jl`, shared with `thread_scaling.jl`. A
  GPU-gated testset (`NBODY6_GPU_TESTS=1`) builds the engine with CUDA in a
  temporary tree and runs the 1k input on one device and, when present, on
  two. No NVIDIA device was available for this change: the GPU path is
  verified by the gated suite on the first such host.
- README and docs index carry a figure of the equal-mass merger, at the
  start and after 12 Myr, produced by the pipeline itself; the asset lives
  in `docs/src/assets/` and is shared by both. The README also carries the
  animation of the same run (281 kB), which shows the dynamics the static
  pair cannot.
- Badges for Aqua static QA, the supported Julia version, and repository
  status. A documentation badge is deliberately absent until the site is
  deployed, and a Backend-workflow badge until that workflow has run.
- The showcase sweep is rebased on a geometry the engine can start from
  (2 pc clusters on a 10 pc orbit, three seeds): with 1 pc clusters on a
  5 pc orbit the derived neighbour radius is too small and every merger
  point hangs during neighbour-list construction. (Superseded: the hang
  was the interval digit counter, see "Fixed"; the neighbour radius was
  never involved. The geometry is kept because the showcase document was
  produced with it.)
- Control runs (F14): `control_merger_dict`/`write_control_merger_config`
  derive the isolated single-cluster equivalent of a merger TOML (summed
  `N`, `N`-weighted half-mass radius, cluster 1's model, IMF and binaries,
  other sections verbatim); the generator accepts one cluster in explicit
  mode. `controls = true` in a sweep adds a control companion to every
  point (`kind`/`control_of` in the index and summary), the comparison and
  ensemble figures use merger points only, and `plot_control_comparison`
  draws each merger against its control. Merger output intervals may be
  given in Myr (`tcrit_myr`, `dtadj_myr`, `deltat_myr`, `dtplot_myr`),
  converted at generation with the realised time unit and recorded in
  `merger_ic.toml` and the summary, so merger and control span the same
  physical time.
- Remnant diagnostics (F13): `remnant_diagnostics` analyses the bound
  remnant of the whole system at every snapshot — Casertano & Hut (1985)
  core radius (`core_radius`, engine densities or a sixth-neighbour
  estimate), half-mass radius, rotation (`rotation_analysis`: angular
  momentum, spin axis, intrinsic `λ_R`, Peebles `λ_P`, `v_rot/σ` profile),
  Allison et al. (2009) `Λ_MSR` mass segregation with a segregation time
  (`mass_segregation`), and the union-find coalescence time
  (`coalescence_time`, now marked on the separation figure). Written to
  `remnant_diagnostics.csv` in the run directory and drawn by
  `remnant_figures` (rotation, rotation profile, structure, mass
  segregation).
- Seeded ensembles (F11): `sweep_ensembles` groups the completed points
  of a sweep by grid values, `ensemble_statistics` interpolates the
  members' series onto a common time grid and takes the median and the
  central 68 % and 95 % intervals (`EnsembleStatistics`), and
  `plot_sweep_ensemble` draws them as line and bands for the Lagrangian
  radius, energy error, star or pair count, coloured by one grid axis
  with the other axes held fixed. A sweep without `[sweep.grid]` is a
  seed ensemble of the base configuration; `sweep_figures` includes the
  ensemble figures whenever a sweep has more than one seed.
- Parameter sweeps (F6): `scripts/run_sweep.jl` runs a sweep TOML
  (`[sweep]` with base pipeline and merger configs, `seeds`,
  `concurrency`, `omp_threads`; `[sweep.grid]` axes as dotted merger-TOML
  keys, Cartesian product). `prepare_sweep` writes one directory per point
  with the derived `merger.toml`/`config.toml` (validated through the
  regular parsers), `run_sweep` executes the points as concurrent worker
  processes (`run_sweep_point`) while keeping `sweep_index.toml` current,
  `write_sweep_summary` collects final N, pair count, energy error and
  virial ratio into `sweep_summary.csv`, and `sweep_figures` draws the
  Lagrangian radius and |ΔE/E| of every point on common axes coloured by
  one grid axis. `run_pipeline`/`run_simulation` accept a fixed `run_id`.
- Binary diagnostics (F5b): `read_binary_evolution` parses the engine's
  `bev.82_*` records (KS-regularised pairs with orbital elements and the
  SSE state of both components; `BinaryRecord`, `BinaryEvolutionSnapshot`);
  `binary_population` reduces a run to pair counts, binary fraction and the
  Heggie hard/soft split, with the energy scale `⟨m⟩ σ²` measured on the
  systems of the nearest snapshot (`hardness_scale`, `binary_scales`,
  `binary_hardness`); `plot_binary_population`,
  `plot_binary_orbital_elements` and `plot_binary_period_distribution` join
  the pipeline, controlled by `postprocess.read_binary_evo` and
  `binary_evo_pattern`.
- Primordial binaries (F5a): `[merger.clusterN.binaries]` (`BinarySpec`)
  pairs a fraction of a cluster's stars (random or uniform-`q` pairing,
  Kroupa 1995 periods or log-uniform semi-major axes, thermal or circular
  eccentricities); systems are sampled, virialised, and truncated as
  units and expanded into Keplerian pairs written first in `dat.10` with
  `NBIN0` and `KZ(8) = 2` (`sample_binaries`, `expand_binaries`). Cluster
  membership is now a body-index vector (`cluster_blocks` in
  `merger_ic.toml`, pair counts in `merger_summary.txt`), and every
  per-cluster diagnostic accepts it. The hard fraction of each cluster's
  pairs is reported.
- Radial profiles (F4): `radial_profile`, `cluster_profiles`, and
  `system_profile` measure density, radial and tangential velocity
  dispersions, and anisotropy in log-spaced shells about a cluster's own
  centre (bound members) or the whole system; `model_density` evaluates
  the generating King or Plummer profile; `plot_density_profiles` (with a
  `ρ/ρ_model` ratio strip) and `plot_velocity_dispersion` join the merger
  plot suite, drawn for the first and last snapshots. `bench/thread_scaling.jl`
  sweeps thread count and N for the cost model of the science sweeps.
- `[simulation] startup_timeout`: a start-up watchdog that terminates a run
  which reports no adjustment beyond t = 0 within the given wall-clock time.
  Signal terminations are recorded in `RUN_INFO.toml` as negative signal
  numbers.
- `install.ref`: the backend commit, tag, or branch checked out after
  cloning (default: the validated upstream v2026.07 commit `618d7a4`),
  so a fresh install builds the version the package was verified against.
- Engine-dependent tests behind `NBODY6_BINARY_TESTS=1` (build in a
  temporary tree or reuse `NBODY6_BACKEND_ROOT`; single run, merger run
  with restart, point-mass tidal field, telemetry on the live process)
  and a weekly/manual `Backend` GitHub workflow that builds the backend
  from source and runs them.
- `[merger.tidal]` (`TidalSpec`): external galactic field for merger runs —
  `KZ(14) = 1` solar-neighbourhood tide, `2` point-mass galaxy (`gmg`,
  `rg0`), `5` MWPotential2014 (`rg`, `vg`), with the `&INXTRNL0` namelist
  written accordingly; options 3 and 4 are refused because the engine
  rescales all velocities on those paths, and a tidal configuration with
  `qe < 0.01` is refused because the engine's energy check omits the tidal
  potential energy and would halt the run. `bound_fraction(snap)`: the
  snapshot-based bound mass fraction of the whole system.
- Restarts: `restart_simulation(run_dir; tcrit_extra, dump, tcrtp0)`
  continues a run from an engine COMMON dump (`KSTART = 2`) in the same
  output directory with appended outputs; the original input is copied
  into `output/` at launch and recorded in `RUN_INFO.toml`, which now
  keeps a `segments` list (one entry per launch) with accumulated elapsed
  time and per-segment telemetry CSVs.
- `[merger.nbody6]` (`Nbody6ParameterSpec`): the integration parameters
  of `merger.inp` (`QE`, `ETAI`, `ETAR`, `NNBOPT`, `RS0`, `RMIN`, `DTMIN`,
  `KZ(16)`) are configurable; the derivable ones default to values scaled
  to the smallest member cluster (`resolve_nbody6_parameters`) instead of
  the former fixed single-cluster constants, and `QE` defaults to
  `2e-4` instead of the former `1.0` (which disabled the energy check).
- Generation-time regime diagnostics: the combined virial ratio `T/|W|`
  and `RBAR/r_hm,min` are computed, printed in `merger_summary.txt`, and
  stored in `merger_ic.toml`; cold-collapse (`Q < 0.3`) and unresolved
  members (`RBAR/r_hm,min > 5`) warn, an `rs0` wider than the smallest
  member is refused.
- Fail-fast bound on the rescaled Kroupa IMF: a rescale factor outside
  ×[0.5, 2] is refused at load time, outside ×[0.7, 1.4] warned.
- No hardcoded input: every numerical field of `merger.inp`'s `&INNBODY6`
  and `&ININPUT` namelists is a `[merger.nbody6]` key (run-time and
  termination limits, block thresholds, output multipliers, chain and KS
  parameters) with an explicit `[merger.nbody6.kz]` override table; the
  stellar-evolution settings (`KZ(19)`, `Level`, `ZMET`, `EPOCH0`,
  `DTPLOT`) form `[merger.stellar]` (`StellarSpec`), with evolution
  switchable off; `&INDATA` mass bounds follow the clusters' IMFs.
- Crossing times of the configuration and of the smallest member
  (`crossing_time`, N-body units and Myr) and the N-body time unit are
  reported in `merger_summary.txt` and `merger_ic.toml`; generation warns
  when `deltat` exceeds the smallest member crossing time.
- Per-cluster structure from snapshots (`cluster_structure`,
  `ClusterStructure`): self-consistently bound members, shrinking-sphere
  centre, 10/50/90 % radii, velocity dispersion, bound mass fraction, and
  virial ratio per initial cluster and snapshot; `plot_cluster_structure`
  (bound half-mass radius with the engine's global `r₅₀` overlaid, bound
  mass fraction) joins the merger plot suite. `per_cluster_virial` and
  `plot_cluster_virial` use the bound members by default
  (`bound_only = false` restores the all-member ratio).
- `[simulation] omp_threads`: OpenMP thread cap for the backend, exported
  as `OMP_NUM_THREADS` by the launch script (`0` = runtime default);
  oversubscription against the host's logical CPUs warns at launch. The
  effective and the backend-reported thread counts and the MPI rank count
  are recorded in `RUN_INFO.toml [run]`.
- Runtime hardware telemetry (`src/telemetry.jl`): exact child CPU
  accounting through `getrusage(RUSAGE_CHILDREN)` (user/system seconds,
  CPU efficiency against the reserved threads) for every run, plus an
  opt-out sampler (`[simulation] telemetry_interval`, default 5 s) of the
  backend process tree (resident memory, high-water mark, cores busy, load
  average) and, with `build.enable_gpu`, `nvidia-smi` utilization, memory,
  power, and temperature. Time series in `runs/<id>/telemetry.csv`,
  summary in `RUN_INFO.toml [telemetry]`, together with the backend's own
  last cumulative timing table (`[telemetry.backend_timing]`) and the
  regular-force kernel throughput from its stderr profile
  (`[telemetry.force_kernel_gflops]`).
- Documentation of what the engine can and cannot do with multi-cluster
  initial conditions ("Feasibility and limitations" in the merger docs):
  single-centre diagnostics, escape-radius behaviour and the `escape.F`
  format overflow, hardcoded integration parameters, the `QE = 1.0`
  choice, the super-particle mass caveat of the shipped small
  configurations, and the recommended 2–5 cluster regime; checked against
  the engine source and two small confirmation runs.

### Changed
- Manifests resolved with Julia 1.13, the version of the target GPU hosts;
  Julia 1.10 remains the compat floor.
- The `nvcc` and `pkg-config` probes and the live monitor's ADJUST parser
  catch only the exception types they expect and log the miss at debug
  level instead of swallowing every error.
- Comments citing section numbers of an external style guide removed.
- Multi-panel montages (`plot_snapshot_evolution`, `plot_hr_evolution`) are
  one column wide: `_fig_multipanel` keeps the preset width and divides it
  among the panels (a reserved colorbar column comes out of the panel
  area) instead of multiplying the width by the column count, so a montage
  enters a document at native size. Montage panels use three ticks per
  axis, a compact gap when the inner tick labels are hidden, a data-free
  band above the data for the time annotation, and markers scaled with the
  panel width. Two-panel stacks are unchanged. The canvas height reserves
  one axis-decoration strip in each direction (`_AXIS_PROTRUSION`), so the
  boxes keep their intended size when a single row of small panels carries
  its x labels; the initial-condition overview (`merger_ic_overview`, three
  projections in a row) follows the same layout, with one shared y label at
  the single preset, fewer colorbar ticks, and the orbit parameters in an
  annotation row above the panels.
- The derived initial neighbour radius is `2 r_h (2 NNBOPT/N_min)^{1/3}`
  (capped at the member half-mass radius): the undoubled value hung the
  engine's neighbour-list initialisation on one of three random
  realisations of the two-cluster demo (outer stars without neighbours),
  independent of the tidal field, the thread count, and the KS parameters;
  six seeds ran with the doubled value.
- Explicit-orbit `velocity` entries of the merger TOML are in km s⁻¹
  (formerly the generator's code unit, 0.0656 km s⁻¹); the conversion
  happens in `combine_clusters_explicit`, `merger_ic.toml` stores km s⁻¹
  (schema version 3), and the shipped explicit configurations were
  converted. The code velocity unit derives from an explicit
  `G = 4.30091e-3 pc (km/s)² M☉⁻¹` (0.06558 km s⁻¹ instead of the rounded
  0.06557).
- Shipped merger configurations sample the natural Kroupa IMF (no
  `mass_total`); explicit-orbit velocities were rescaled by the square
  root of the mass ratio to preserve each configuration's virial state.
  `verif_triorbit.toml` uses equal 0.6 M☉ bodies (`m = 900` M☉ per
  cluster, `|v| = 9.306` code units) so the Lagrange equilibrium is exact.
- `MergerConfig` carries an `nbody6::Nbody6ParameterSpec` field;
  `generate_merger_inp` requires a resolved `nbody6` keyword.
- `virialise!` shares its energy evaluation with the new
  `_kinetic_and_potential` helper.

Structural audit: redundancy, dead-code, and naming sweep (no physics
changes; all numerical outputs unchanged).

### Changed (audit)
- One plot dispatcher: `generate_plots(results, vis::VisualizationConfig;
  sim_dir, animations)` is the single dispatch core; the `Nbody6Config`
  method and `postprocess_external` both route through it (removes a
  ~90-line drifting copy of the plot pipeline).
- One simulation executor: `_execute_simulation` backs both
  `run_simulation` and the merger path (removes a ~35-line copy).
- `setup_two_cluster_orbit` returns per-cluster index ranges, matching
  `combine_clusters_explicit`; the duplicate truncation recount in
  `generate_merger_ic` is gone.
- `ClusterSpec` keyword constructor is structured-only (`profile::DensityProfile`,
  `imf::IMFSpec`); the flat schema remains the canonical *TOML* form,
  handled by `load_merger_config` (decision D6).
- `postprocess_external` kwargs: `make_plots`/`make_animations` (no longer
  shadow `generate_plots`), `column`/`units` presets replace the dead
  `figsize`; defaults now match the project visualization standard (PDF,
  single column). Same for `run_merger_pipeline(make_plots = ...)`.
- lagr.7 reader is modern-single-row-format only, with fail-fast layout
  validation (`rows_per_block` and the block-format path removed).
- Naming: `missing_deps` (was `missing`, shadowing `Base.missing`),
  `w_hat`/`sqrt_w`/`abs_w`/`r_jacobi*`/`integral_ab`/`n_threads`
  replace camelCase/ad-hoc locals; code-unit velocity and Plummer
  r_hm/a factors are named constants.
- Plotting layer deduplicated behind shared helpers (rendered output
  unchanged, verified by figure QA): `_annotate!`/`_no_data_note!`
  standard in-axis annotations (18 sites), `_log_color_range` mass colour
  scale (4), `_envelope_stats` min/max/mean bands (2), `_hr_limits` (3),
  `_scatter_escaper_classes!` two-class encoding (2), shared Lagrangian
  axis prologue for plot + animation, and named layout constants
  (`_COLORBAR_WIDTH`, `_COLORBAR_COLGAP`, `_TWO_PANEL_ROWGAP`) replacing
  scattered literals; IMF reference slopes read the canonical
  `_KROUPA_ALPHAS`/`_KROUPA_BREAKS`.
- Test suite: smoke coverage added for the previously untested merger
  plot functions (`plot_merger_ic`, `plot_cluster_separation`,
  `plot_cluster_virial`, `per_cluster_virial`).
- Run summary is machine-readable: `RUN_INFO.txt` replaced by
  `RUN_INFO.toml` (run identity/timing, package + backend commits, output
  inventory) extended with a hardware fingerprint — host, OS, CPU model
  and logical cores, total memory, Julia version, Julia/BLAS thread
  counts, and the GPU (nvidia-smi probe) on `enable_gpu` builds. The same
  fingerprint is stored in `merger_ic.toml`; `export_for_paper` reads the
  TOML provenance instead of scraping text. Adds the `LinearAlgebra`
  stdlib dependency (BLAS thread count).
- `docs/make.jl` self-activates its environment and develops the package
  by path (`julia docs/make.jl` with no `--project` flag); the Docs
  workflow simplified accordingly.
- README: consolidated entry-points table (instantiate, pipeline,
  verification suite, tests, benchmarks, docs build) and a `TestEnv.jl`
  note; config headers trimmed to the one-line form the config comment
  policy prescribes (run narratives and usage examples belong to docs).

### Fixed
- An interrupted `git clone` of the engine (a killed run, a host reset) left
  a directory without `configure`, and the next build took it for a finished
  checkout and died on a bare `ENOENT` spawning `./configure`. The source
  tree is now checked before use and repaired: `git checkout --force` first,
  which keeps untracked work, then a fresh clone when no usable git state
  remains.
- GPU builds with the toolkit off `PATH` failed in the engine's `configure`
  ("Cannot find CUDA compiler nvcc"): the CUDA environment reached `make`
  and the launch script but not `configure`, whose `--with-cuda` fallback
  never runs (it reuses the cached `PATH` check). `configure` now runs with
  the toolkit's `bin` on `PATH` and an explicit `--with-cuda=<path>`.
- GPU builds on hosts whose glibc (2.42 and later: Fedora 43/44, Ubuntu
  26.04, Debian 13) declares `rsqrt`/`rsqrtf` with an exception
  specification the CUDA ≤ 13.1 headers lack: `nvcc` rejected every host
  compiler with "exception specification is incompatible". The probe
  recognises the conflict and retries the host-compiler choice with
  `-U_GNU_SOURCE -D_DEFAULT_SOURCE`, keeps the override in the build flags,
  and names the header-patch alternative when that does not resolve it;
  the validation host record carries the glibc version.
- Plain-text labels (axis labels, tick labels, annotations, legends) were
  rendered in Makie's default sans font instead of Computer Modern whenever
  the package was loaded from its precompile cache: the theme held FreeType
  faces created at precompile time, serialised with null pointers, and Makie
  fell back silently. The theme is now built at call time from
  MathTeXEngine's live faces; LaTeX strings were never affected. A test
  asserts live faces in the (fresh) test process.
- Log-axis tick labels: plain decimals throughout (`0.5, 1, 2, 5, 10`
  instead of `5 × 10⁻¹ … 10¹`) on axes whose ticks lie within 10⁻³–10⁴ and
  span at most four decades; in the exponent form `10^1` and the 2×/5×
  multiples of `10^{-1}`–`10^{1}` collapse (`10`, `0.2`, `20`) as the
  typography standard prescribes.
- GPU builds on hosts whose default compiler is newer than the toolkit
  accepts (GCC 16 with CUDA 13.1): the `nvcc` probe now tries
  `-allow-unsupported-compiler`, then `-ccbin` with `CUDAHOSTCXX` and the
  versioned `g++`/`clang++` compilers on `PATH`, captures the compiler
  output through a file and reports an excerpt of every attempt; a
  configured `-ccbin` is final. The validation driver skips the benchmark stage, with the reason
  recorded, when a binary it needs is missing.
- The build dependency check accepts an `nvcc` under the configured or
  auto-detected toolkit rather than only on `PATH`; the GPU-gated tests
  query the same toolkit. On a host with the toolkit under `/usr/local/cuda`
  and nothing on `PATH`, GPU builds refused to start with "Missing
  dependencies: nvcc".
- GPU builds with a CUDA 13 toolkit: the engine's vendored `helper_cuda.h`
  reads `cudaDeviceProp.clockRate` and `.computeMode`, which CUDA 13.0
  removed; the `CUFLAGS` override now puts `deps/cuda` ahead of the engine's
  `extra_inc/cuda` and is passed for every GPU build, not only when target
  architectures were resolved.
- `scan_output`/`postprocess_external` read the engine's `bev.82_*`
  records (`:binary_evo`, category "Regularised binaries", the three
  binary figures in the plot inventory), so the external path produces the
  same binary diagnostics as the config-driven pipeline instead of
  silently dropping them.
- Start-up hang of merger runs. The engine's `string_left.f` counts the
  decimal digits of `DELTAT`, `DTADJ` and `DTPLOT` by multiplying by ten
  until the value is an integer, with a default-kind `int`; an interval
  such as `0.6302`, which the Myr conversion produces, never becomes an
  integer in binary arithmetic, overflows the conversion past 2³¹, and
  loops forever inside the first output (backtrace: `string_left_` ←
  `output_` ← `adjust_`). The generator now writes the three intervals as
  the nearest dyadic rational with an exact decimal expansion
  (`engine_interval`, change below 0.4 %, logged) instead of a
  four-decimal format, and the emulated counter is a unit test on the
  intervals of every run that hung or completed so far. The hang had been
  attributed to a small initial neighbour radius; that rule and its guard
  stay, as a matter of start-up cost only, and the documentation is
  corrected.
- A watchdog kill is now written to `RUN_INFO.toml` (segment field
  `watchdog = true`) before the error is raised; previously the run
  directory had no summary at all.
- The sweep's pre-launch binary check called the removed two-argument
  locator; it now selects the binary of the base config's build variant.
- Completion monitor (`[simulation] exit_grace`, default 120 s): an engine
  that printed `END RUN` and did not exit is terminated once its output
  directory has been idle for the grace period (a final COMMON dump in
  progress keeps it alive) and recorded as completed (`completed`,
  `terminated_after_completion` in the segment; `completed` column in the
  sweep summary), so post-processing proceeds instead of waiting on an
  external timeout.
- The main-sequence turnoff label of the evolutionary-clock figure is placed
  on whichever side of the reference line has room; it was clipped by the
  axis whenever the distribution ended just past t = T_MS.
- The rotation figure's alignment label is placed at the end of the axis
  away from the coalescence marker, and the lambda_R annotation treats the
  marker as occupied when choosing its corner.
- The IMF reference slopes are continuous at the 0.5 Msun break; the
  alpha = 2.3 segment shared the low-mass anchor and was drawn a factor of
  five below the histogram, which made a correctly sampled population look
  wrong.
- The hard/soft energy scale of the binary diagnostics is taken over the
  bound systems: escapers and kicked stellar remnants had inflated the
  mass-weighted dispersion by a factor of several at single epochs.
- The per-cluster virial figure accepted only index ranges and crashed on
  the body-index blocks written since the primordial-binary change.
- `Manifest.toml` self-entry still carried the pre-rename package name
  (`Nbody6Setup`); three shipped merger TOML headers likewise.
- README structure tree: `bench/` (was stale `benchmark/`), previously
  missing `.github/`, `.JuliaFormatter.toml`, `CHANGELOG.md`, and
  environment files added.
- `RescaledKroupaIMF` warning no longer advises a silencing mechanism
  that did not exist.

### Removed
- Dead code: unused `SnapshotHeader` accessors (slot map documented on the
  struct instead), dead `"astro"` format branches in
  `load_merger_ic_result`, the vestigial `kz22` keyword, the unreachable
  HDF5 branch in `postprocess`, legacy flat-layout run-directory
  fallbacks, `BenchmarkTools` from the package extras, stale `.gitignore`
  entries.

## [0.2.0] — 2026-08-07

First version pushed to the private remote (CI/CD shakedown).

### Added
- Fail-fast configuration validation naming the offending `section.key`.
- Static QA in the test suite (Aqua, ExplicitImports, JET package-scoped);
  StableRNGs test streams; real-output fixtures for every reader.
- Escaper analysis suite (cumulative mass loss, velocity–time, escape
  anisotropy from the new angle columns) and SSE-quantity plots
  (mass segregation via RI [pc], evolutionary clock t/T_MS, core masses).
- Journal column-width figure presets with true-print-size PDF export;
  horizontal top-outside legends; decade log ticks; Okabe–Ito semantic
  colour system; physical-units axes (pc / Myr / km s⁻¹ / M☉) by default.
- Provenance: package + backend commit stamps in run metadata;
  `export_for_paper` with provenance sidecars; structured per-run logging.
- Merger IC generator: validated King sampler (c(W0) to <1% of published
  values), Kroupa IMF modes, Kepler/explicit orbits, seeded reproducibility
  through to Nbody6's NRAND.

### Changed
- Package renamed Nbody6Setup.jl → Nbody6Dynamics.jl (2026-08-02).
- Backend refreshed to upstream v2026.07; sev.83 reader targets the
  15-token format exclusively (RI in pc; TM/MC/RCC/RE fields).
- Formatter-enforced style (`.JuliaFormatter.toml`); benchmarks moved to
  a bench-local environment.

### Removed
- Dead HDF5 reader path; "astro"/KZ(22)=10 output branch; legacy seed
  sentinel; compatibility shims (private-project no-compat policy).

## [0.1.0] — 2026-04-24

Pre-remediation baseline (Nbody6Setup.jl): install/build pipeline, conf.3
and diagnostics readers, initial plotting, first merger-IC implementation.
