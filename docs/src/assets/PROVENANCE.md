# Provenance of the documentation assets

`merger_evolution.png` and `merger_evolution.gif` show the merger of two King
clusters (`W0 = 6`, Kroupa IMF over 0.08–100 M☉, half-mass radius 2 pc) of
25 000 stars each, 47 999 stars and 2.75 × 10⁴ M☉ together after Jacobi
truncation, released 10 pc apart at the apocentre of an orbit of eccentricity
0.5 and followed for 50 Myr. The configuration is
`input_files/gpu/gpu_pipeline.toml` with `input_files/gpu/merger_50k.toml`
(seed 31), run as the `gpu` stage of `scripts/run_gpu_validation.jl`.

| Asset | What | How |
|---|---|---|
| `merger_evolution.png` | projected surface mass density at 0, 4.0 and 50.1 Myr, xy projection, one panel each on common axes | `plot_snapshot_evolution` on those three snapshots with `snapshot_render = "density"` and `density_mass_frac = 0.95`, other style keys at their defaults; downscaled to 1600 px width for the README |
| `merger_evolution.gif` | the same run animated, all 26 snapshots (one per 2.0 Myr), fixed axes | `animate_cluster`, xy projection, same style; halved in resolution and frame-optimised for the README |

Run. `gpu_validation_20260924_002436_gpu` (not shipped), produced on
2026-09-24 on one NVIDIA H200 NVL (AMD EPYC 9755 host, eight OpenMP threads,
CUDA build of the engine) by package commit `b376fa7`, engine
Nbody6PPGPU-beijing at commit `618d7a4` (upstream v2026.07); 454 s of engine
time for the 50 Myr. The figures were rendered from that run's snapshots with
the figure code of package commit `24740d0`; no simulation was repeated.

Licence. Simulation output and figures produced by the package author,
released with the package under its MIT licence.
