# Changelog

All notable changes to Nbody6Dynamics.jl. Follows semantic versioning;
pre-1.0 minor versions may break APIs (private project, no-compat policy).

## [Unreleased]

### Added
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
