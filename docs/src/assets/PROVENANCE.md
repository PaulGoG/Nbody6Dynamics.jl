# Provenance of the documentation assets

`merger_evolution.png` and `merger_evolution.gif` show the equal-mass
showcase merger: two King clusters of 4000 stars each (`W0 = 6`, Kroupa
IMF), released 6 pc apart on an orbit of eccentricity 0.9 and followed for
20 Myr. The configuration is `input_files/showcase/equal_pipeline.toml` with
`equal_merger.toml` (seed 21).

| Asset | What | How |
|---|---|---|
| `merger_evolution.png` | the initial conditions and the remnant at 20 Myr, xy projection, one panel each | `plot_snapshot_evolution` on the first and last snapshot, downscaled to 1600 px width for the README |
| `merger_evolution.gif` | the same run animated, one frame per Myr, fixed axes | `animate_cluster`, xy projection, halved in resolution and frame-optimised for the README |

Run. `runs/merger_showcase_equal_20260911_200411_f8f6` (not shipped),
produced on 2026-09-11 by `scripts/run_setup.jl` from the configuration
above, package at commit `f89eb2e` (with local changes at the time), engine
Nbody6PPGPU-beijing at commit `618d7a4` (upstream v2026.07), 72 s of engine
time on eight threads. The figures were rendered from that run's snapshots
with the release figure code (package commit `99cbc43`); no simulation was
repeated for the release, so the run data predate the strict configuration
parser of 0.3.0 (its frozen `config.toml` carries the since-removed
`visualization.figsize` key).

Licence. Simulation output and figures produced by the package author,
released with the package under its MIT licence.
