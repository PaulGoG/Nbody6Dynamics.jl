# Provenance of the documentation assets

Every asset here comes from one of the showcase runs of 2026-11-02 listed in
the Showcase Cases page (`input_files/showcase/`), produced on a laptop
(Intel Core Ultra 7 155H, AVX build of the engine at the pin `v2026.07.23`,
engine commit `ce10277`) by the package tree of commit `916fb0a` with the
three-cluster flagship configuration of the following commit, and rendered
again on 2026-11-03 from the same engine output by the figure code of
release 0.7.0 (the merger epochs on the separation figure). Nothing was
edited by hand: each figure is the package's own export, reduced in size with
ImageMagick for the web (`magick in.png -resize 1200x -strip out.png`; the
README montage to 1600 px, the animation halved and frame-optimised), as
`utilities/render_readme_assets.jl` of the development workspace does for the
two README files.

| Asset | Run | What | How |
|---|---|---|---|
| `merger_evolution.png` | `merger_showcase_flagship_20261102_170948_a01b` | projected surface mass density at 0, 4.0 and 50.1 Myr, xy projection, one panel each on common axes | `plot_snapshot_evolution` on those three snapshots, `snapshot_render = "density"`, `zoom_frac = 0.15`, `density_mass_frac = 0.99` |
| `merger_evolution.gif` | the same | the same run animated, all 26 snapshots (one per 2.0 Myr), fixed axes | `animate_cluster`, xy projection, same style; halved in resolution and frame-optimised |
| `showcase_equal_separation.png` | `merger_showcase_equal_20261102_165008_d8ff` | separation of the two density centres against the smaller half-mass radius; contact, centre coalescence and core coalescence marked, the regimes shaded, the centre offset after the core coalescence | `plot_cluster_separation` as drawn by the pipeline (`plots/merger_cluster_separation.png`) |
| `showcase_binary_orbital_elements.png` | `merger_showcase_binary_20261102_165418_b7bf` | the regularised pairs at 15 Myr in the a–e plane, hard and soft by the Heggie criterion | `plot_binary_orbital_elements` at the last snapshot (`plots/binary_orbital_elements_final.png`) |
| `showcase_tidal_escapers.png` | `merger_showcase_tidal_20261102_165852_4a2d` | cumulative escaped mass and the escape velocities by stellar class | `plot_escapers` as drawn by the pipeline (`plots/escapers.png`) |
| `showcase_sweep_lagrangian.png` | `sweep_showcase_20261102_170236` | half-mass radius of the six mergers on common axes, colour by eccentricity | `plot_sweep_lagrangian` as drawn by the sweep driver (`plots/sweep_lagrangian.png`) |
| `showcase_sweep_ensemble.png` | the same | median half-mass radius per eccentricity with 68 and 95 % bands over the seeds | `plot_sweep_ensemble` as drawn by the sweep driver (`plots/ensemble_lagrangian.png`) |

The flagship integration (3 × 25 000 stars, 58 221 bodies after the Jacobi
truncation, 50 Myr) took 3328 s (55 min) of engine time at eight OpenMP
threads; the small cases a few minutes each (the table of the Showcase Cases
page).

Licence. Simulation output and figures produced by the package author,
released with the package under its MIT licence.
