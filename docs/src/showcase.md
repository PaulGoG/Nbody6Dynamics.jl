# Showcase Cases

`input_files/showcase/` holds five cases that take the pipeline from a merger
configuration to a finished, figured run: four small ones that complete within
minutes on a workstation, and the 2 × 25 000-star merger over 50 Myr at which
the package's CUDA path was validated. Each case is a merger TOML with a
matching pipeline TOML (paths relative to the folder), so one invocation runs
it end to end:

```bash
julia scripts/run_setup.jl input_files/showcase/equal_pipeline.toml     # also binary_, tidal_, flagship_
julia scripts/run_sweep.jl input_files/showcase/sweep.toml               # the sweep with its controls
```

Every case gives its intervals in Myr (`tcrit_myr`, `dtadj_myr`, `deltat_myr`,
`dtplot_myr`), so the merger, its isolated control and the points of the sweep
cover the same physical span whatever their N-body time unit; the pipeline
TOMLs use the AVX engine built under `backend/Nbody6PPGPU-beijing`
(`gpu/gpu_pipeline.toml` builds the CUDA variant) and write into `runs/`.

## The cases

| Case | Clusters | Orbit | Span | Tolerance `qe` | What it shows |
|---|---|---|---|---|---|
| `equal` | 2 × 4000, King `W0 = 6`, `r_h = 1` pc, Kroupa 0.1–50 M☉ | apocentre 6 pc, `e = 0.9` | 20 Myr | 0.01 | an equal pair that falls together and coalesces within the run; the remnant sheds a halo of loosely bound stars (the README figure and animation) |
| `binary` | 4000 (`W0 = 6`, 1 pc) + 1000 (`W0 = 5`, 0.7 pc), `q ≈ 0.25`, 20 % primordial binaries in both | apocentre 5 pc, `e = 0.7` | 15 Myr | 0.01 | the binary suite (pair counts with the hard/soft split, `a`–`e` diagrams, period distributions) on an unequal pair, and the relaxed tolerance a binary-rich run needs |
| `tidal` | 2 × 2500 as `equal` | apocentre 6 pc, `e = 0.9`, point-mass galaxy of 10¹¹ M☉ at 8.5 kpc (`KZ(14) = 2`) | 20 Myr | 0.05 | the same encounter inside an external field: escapers removed on the distance criterion, and an engine energy check that measures the tidal work rather than the integration error |
| `sweep` | 2 × 500, King `W0 = 6`, `r_h = 1` pc | apocentre 5 pc, `e ∈ {0.5, 0.9}`, seeds 2, 3, 4, with an isolated control per point | 20 Myr | 0.01 | a twelve-point grid: comparison figures on common axes, seed ensembles with 68/95 % bands, and each merger against its control |
| `flagship` | 2 × 25 000, King `W0 = 6`, `r_h = 2` pc, Kroupa 0.08–100 M☉ | apocentre 10 pc, `e = 0.5` | 50 Myr | 0.01 | the merger of the hardware validation: the same configuration as `gpu/merger_50k.toml`, run on every validated host on both engine builds |

The infall time from apocentre is half the Kepler period, so the eccentric
cases coalesce early in the run (`equal` at about 3 Myr) and the sweep's two
eccentricities separate in time (12 Myr at `e = 0.5` against 8 Myr at
`e = 0.9`), which is what its eccentricity axis is for.

## Reference runs

The four small cases were run on 2026-09-11 on a laptop (Intel Core Ultra 7
155H, four OpenMP threads, package commit `f89eb2e`, engine
Nbody6PPGPU-beijing at `618d7a4`, upstream v2026.07). Engine wall time is the
integration alone; the pipeline adds the initial conditions, the readers, the
diagnostics and about fifty figures and five animations per run, a few
minutes at these sizes.

| Case | Run | Bodies after truncation | Engine wall | Final `DETOT` |
|---|---|---|---|---|
| `equal` | `merger_showcase_equal_20260911_200411_f8f6` | 7832 | 72 s | +1.7×10⁻³ |
| `binary` | `merger_showcase_binary_20260911_200459_aac4` | 4895 | 124 s | +3.9×10⁻³ |
| `tidal` | `merger_showcase_tidal_20260911_200554_1584` | 4878 | 37 s | +1.5×10⁻³ (`ERRTOT` −5.3×10⁻², the tidal work) |
| `sweep` | `sweep_showcase_20260911_202212` | 941–974 per merger, 995–999 per control | 42–52 s per point, four at a time | ≤ 5×10⁻³ at every point |

`DETOT` is the engine's cumulative relative energy change at the end of the
run; the merger runs evolve stars and form binaries, so it sits well above
the 10⁻⁷ of a quiet single cluster, and inside the tolerance the
configurations set. The documentation assets (`docs/src/assets/`) are
rendered from the `equal` run; their provenance is recorded next to them.

The `flagship` configuration is the merger of the hardware validation
(`gpu/merger_50k.toml` differs from `flagship_merger.toml` in its header
comment only). Its integration over 45.8 N-body units, 50.6 Myr, at eight host
threads measured in the validation chains of 2026-09-24 to 2026-09-26:

| Host | CUDA build | AVX build |
|---|---|---|
| RTX 5090, Ryzen 9 9950X | 360 s | 636 s |
| RTX 5070 Ti, i9-13900KS | 388 s | 821 s |
| H200 NVL, EPYC 9755 | 454 s | 851 s |
| RTX 2080 Super Max-Q, i7-10750H | 954 s | 2136 s |
| Tesla T4, EPYC 7551P | 1470 s | 4192 s |

At this size the device removes about half of the work, because the
irregular force, host work in both builds, is as large as the regular force;
the gain grows with N (see the manual, section Validated hardware).

## What a run produces

`runs/<run_id>/` holds the frozen `config.toml`, `RUN_INFO.toml` (timing, the
thread layout, the hardware fingerprint, the commits, the telemetry summary),
`nbody6dynamics.log`, `telemetry.csv`, the engine output under `output/`
(`conf.3_*` snapshots, `out1000`, `lagr.7`, `esc.11`, `sev.83_*`, `bev.82_*`,
`merger.inp`, `dat.10`, `merger_summary.txt`, `merger_ic.toml`), the derived
tables (`stellar_census.csv`, `remnant_diagnostics.csv`) and `plots/`:

- initial conditions: projections, overview, velocity field, sampled IMF, density profiles (`merger_ic_*`);
- the encounter: cluster separation, per-cluster virial ratio and structure, density profiles and velocity dispersions at the first and last snapshot (`merger_*`);
- the system: energy error and particle count, Lagrangian radii, escapers and their anisotropy, snapshots (`energy`, `particle_count`, `lagrangian_radii`, `escapers`, `escape_anisotropy`, `snapshot_*`);
- the stars: HR diagrams at three epochs and their evolution, mass segregation, the evolutionary clock, core masses (`hr_*`, `mass_segregation`, `evolutionary_clock`, `core_mass_growth`);
- the binaries: population and hardness, orbital elements, periods (`binary_*`);
- the remnant: structure, rotation, rotation profile, mass segregation (`remnant_*`);
- the run itself: telemetry (`telemetry`);
- animations: the cluster in three projections, the Lagrangian radii, the HR diagram (`*_anim.gif`, `cluster_evolution_*.gif`).

A sweep adds `sweep_index.toml`, `sweep_summary.csv` and, under its `plots/`,
the comparison figures (`sweep_lagrangian`, `sweep_energy`), the ensemble
bands (`ensemble_*`) and the merger-against-control figure
(`control_lagrangian`). The figure routines are documented in the manual,
section Visualisation, and in the API reference.
