# Changelog

All notable changes to Nbody6Dynamics.jl. Follows semantic versioning;
pre-1.0 minor versions may break APIs (private project, no-compat policy).

## [Unreleased]

Structural audit: redundancy, dead-code, and naming sweep (no physics
changes; all numerical outputs unchanged).

### Changed
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
