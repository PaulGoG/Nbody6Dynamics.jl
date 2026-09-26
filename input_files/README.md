# Shipped inputs

One folder per purpose. `example_input("engine/N1k_quick.inp")` returns the
absolute path of any file here for a project outside the checkout; a bare
name (`example_input("N1k_quick.inp")`) is looked up across the folders.

```
input_files/
├── engine/          # Nbody6++ NAMELIST inputs of single clusters (simulation.input_file)
├── mergers/         # merger initial-condition TOMLs (merger.config_file, run_merger_pipeline)
├── verification/    # the two merger targets of scripts/run_verif_suite.jl
├── sweeps/          # sweep TOMLs (scripts/run_sweep.jl)
├── showcase/        # five end-to-end cases, each a merger TOML with its pipeline TOML
└── gpu/             # CUDA-host validation: GPU and CPU builds, the 2 × 25k merger, the 6 × 10⁵ probes
```

| Folder | Files | Documented in |
|---|---|---|
| `engine/` | `N1k_quick`, `N5k_medium`, `N10k_long`, `N25k_production`, `N100k_production` (single clusters of growing size); `gc_bh_subsystem`, `imbh_runaway`, `pop3_cluster`, `tidal_tails`, `young_massive_binaries` (science cases) | Input Files, "Engine inputs" and "Engine inputs: the science cases" |
| `mergers/` | `merger_demo_small`, `merger_equal_mass`, `merger_minor_plummer`, `merger_triple_cluster`, `merger_3cluster_small`, `merger_5cluster_small`, `merger_27cluster_cubic` | Input Files, "Merger configurations"; Cluster Mergers |
| `verification/` | `verif_triorbit`, `verif_3d5cluster` | Input Files, "Verification configurations" |
| `sweeps/` | `sweep_demo` (base: the root `config.toml` and `mergers/merger_demo_small.toml`) | Manual, "Parameter sweeps" |
| `showcase/` | `{equal,binary,tidal,flagship}_{merger,pipeline}.toml`, `sweep.toml` with `sweep_{merger,pipeline}.toml` | Showcase Cases |
| `gpu/` | `gpu_pipeline`, `cpu_pipeline`, `merger_50k`, `merger_600k`, `single_600k`, `{gpu,cpu}_{merger,single}_600k` | Input Files, "The GPU-host cases"; Manual, "Recipe for a CUDA host" |

Relative paths inside a TOML resolve against the file's own directory; the
pipeline TOMLs of `showcase/` and `gpu/` point at the checkout's `backend/`
and `runs/` two levels up. The root `config.toml` runs `mergers/merger_demo_small.toml`.
