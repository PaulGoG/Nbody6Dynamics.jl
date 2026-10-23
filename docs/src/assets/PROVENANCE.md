# Provenance of the documentation assets

Every asset here comes from one of the showcase runs of 2026-10-23 listed in
the Showcase Cases page (`input_files/showcase/`), produced on a laptop
(Intel Core Ultra 7 155H, AVX build of the engine at the pin `v2026.07.9`,
engine commit `c67097a`) by package commit `063a40a` (its tree also carried the source-layout move of the following commit, which changed no code). Nothing was
edited by hand: each figure is the package's own export, reduced in size with
ImageMagick for the web (`magick in.png -resize 1200x -strip out.png`; the
README montage to 1600 px, the animation halved and frame-optimised).

| Asset | Run | What | How |
|---|---|---|---|
| `merger_evolution.png` | `merger_showcase_flagship_20261023_154446_4dd7` | projected surface mass density at 0, 4.0 and 50.2 Myr, xy projection, one panel each on common axes | `plot_snapshot_evolution` on those three snapshots with `snapshot_render = "density"`, `zoom_frac = 0.15` and `density_mass_frac = 0.99`, other style keys at their defaults; 1600 px wide |
| `merger_evolution.gif` | the same | the same run animated, all 26 snapshots (one per 2.0 Myr), fixed axes | `animate_cluster`, xy projection, same style; halved in resolution and frame-optimised |
| `showcase_equal_separation.png` | `merger_showcase_equal_20261023_154447_8ace` | separation of the two density centres against the smaller half-mass radius, coalescence epoch marked | `plot_cluster_separation` as drawn by the pipeline (`plots/merger_cluster_separation.png`) |
| `showcase_binary_orbital_elements.png` | `merger_showcase_binary_20261023_154448_eb6e` | the regularised pairs at 15 Myr in the a–e plane, hard and soft by the Heggie criterion | `plot_binary_orbital_elements` at the last snapshot (`plots/binary_orbital_elements_final.png`) |
| `showcase_tidal_escapers.png` | `merger_showcase_tidal_20261023_154446_b65d` | cumulative escaped mass and the escape velocities by stellar class | `plot_escapers` as drawn by the pipeline (`plots/escapers.png`) |
| `showcase_sweep_lagrangian.png` | `sweep_showcase_20261023_154908` | half-mass radius of the five finished mergers on common axes, colour by eccentricity | `plot_sweep_lagrangian` as drawn by the sweep driver (`plots/sweep_lagrangian.png`) |
| `showcase_sweep_ensemble.png` | the same | median half-mass radius per eccentricity with 68 and 95 % bands over the seeds | `plot_sweep_ensemble` as drawn by the sweep driver (`plots/ensemble_lagrangian.png`) |

The flagship integration (2 × 25 000 stars, 50 Myr) took 2044 s (34 min) of
engine time at eight OpenMP threads; the small cases a few minutes each (the
table of the Showcase Cases page).

Licence. Simulation output and figures produced by the package author,
released with the package under its MIT licence.
