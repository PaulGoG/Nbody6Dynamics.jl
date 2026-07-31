# Nbody6Setup.jl — User Manual

## Table of Contents

1. [Introduction](#1-introduction)
2. [Prerequisites](#2-prerequisites)
3. [Installation](#3-installation)
4. [Configuration Reference](#4-configuration-reference)
5. [Pipeline Modes](#5-pipeline-modes)
6. [Build Phase Details](#6-build-phase-details)
7. [Simulation Execution](#7-simulation-execution)
8. [Post-processing](#8-post-processing)
9. [Visualisation](#9-visualisation)
10. [Using the Julia API Directly](#10-using-the-julia-api-directly)
11. [Output File Reference](#11-output-file-reference)
12. [N-body Units and Conversions](#12-n-body-units-and-conversions)
13. [Troubleshooting](#13-troubleshooting)

---

## 1. Introduction

**Nbody6Setup.jl** automates the full lifecycle of N-body star cluster simulations with [Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing):

1. **Download and compile** the Fortran simulation code
2. **Generate merger initial conditions** (optional; see [Multi-Cluster Merger Simulations](@ref))
3. **Execute** the simulation with proper environment setup
4. **Parse** output files into Julia data structures
5. **Generate** publication-quality plots and GIF animations

Everything is driven by a single `config.toml`; every parsed key and its default is listed below. The single entry point is:

```julia
cfg = load_config("config.toml")
results = run_pipeline(cfg)
```

---

## 2. Prerequisites

- **Linux** (tested on Fedora and Ubuntu)
- **Julia 1.10+**
- **GCC toolchain**: `gcc`, `g++`, `gfortran`, `make`
- **Git** (for cloning the simulation code)

Optional, depending on config: MPI (`mpicc`, `mpif90`, `mpirun`), CUDA toolkit (`nvcc`), HDF5 development headers (`hdf5-devel` / `libhdf5-dev`) for the *build-side* HDF5 patch.

---

## 3. Installation

```bash
cd Nbody6Setup
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

---

## 4. Configuration Reference

Every key below is parsed by `load_config` (`src/config.jl`). Missing keys fall back to the defaults shown. `save_config` writes a frozen snapshot of the full configuration into each run directory.

### `[install]`

| Key           | Type   | Default | Description |
|---------------|--------|---------|-------------|
| `enabled`     | Bool   | `true`  | Enable the install/build phase |
| `source_url`  | String | `"https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing.git"` | Git repository to clone |
| `install_dir` | String | `"backend/Nbody6PPGPU-beijing"` | Source directory, relative to the package root |
| `reinstall`   | Bool   | `false` | Delete and re-clone if `true` |
| `clean_build` | Bool   | `true`  | Run `make clean` before building |

### `[build]`

| Key               | Type     | Default | Description |
|-------------------|----------|---------|-------------|
| `configure_flags` | [String] | `["--enable-mcmodel=large", "--with-par=b1m"]` | Flags passed to `./configure` |
| `enable_mpi`      | Bool     | `false` | Enable MPI parallelism |
| `enable_hdf5`     | Bool     | `true`  | Patch the Makefile with HDF5 build flags |
| `enable_gpu`      | Bool     | `false` | Enable GPU acceleration (requires CUDA) |
| `cuda_path`       | String   | `""`    | CUDA installation path; empty = auto-detect |
| `nproc`           | Int      | `0`     | Parallel `make` jobs; 0 = auto-detect |

### `[simulation]`

| Key             | Type   | Default | Description |
|-----------------|--------|---------|-------------|
| `run_test`      | Bool   | `true`  | Run the simulation phase |
| `input_file`    | String | `"examples/input_files/N10k_noDat10.inp"` | Path to the `.inp` input file |
| `runs_dir`      | String | `"runs"` | Base directory for run output |
| `binary_name`   | String | `"nbody6++"` | Expected binary name |
| `mpi_ranks`     | Int    | `1`     | Number of MPI ranks (when MPI enabled) |
| `run_id_prefix` | String | `"run"` | Prefix for run directory names |

### `[postprocess]`

| Key                   | Type   | Default | Description |
|-----------------------|--------|---------|-------------|
| `enabled`             | Bool   | `true`  | Enable post-processing |
| `data_dir`            | String | `""`    | External data directory. Empty → use `run_dir/output/` from the simulation, or the most recent run under `runs/`. Set → read from this directory instead (any dir with Nbody6++ output) |
| `snapshot_format`     | String | `"conf3"` | Snapshot source. **Only `"conf3"` is supported** (see note below) |
| `snapshot_pattern`    | String | `"conf.3_*"` | Glob pattern for conf.3 files |
| `hdf5_file`           | String | `""`    | **Deprecated, ignored** — kept only so old TOML files still parse |
| `parse_stdout`        | Bool   | `true`  | Parse simulation stdout for ADJUST diagnostics |
| `stdout_file`         | String | `"out1000"` | Name of the stdout capture file |
| `read_lagr`           | Bool   | `true`  | Read Lagrangian radii |
| `lagr_file`           | String | `"lagr.7"` | Lagrangian radii file name |
| `read_escapers`       | Bool   | `true`  | Read escaper data |
| `escapers_file`       | String | `"esc.11"` | Escaper file name |
| `read_stellar_evo`    | Bool   | `true`  | Read stellar evolution snapshots |
| `stellar_evo_pattern` | String | `"sev.83_*"` | Glob pattern for stellar evolution files (matches what this fork's `hrplot.F` writes: `sev.83_<time>`) |

!!! note "HDF5 snapshot support was removed"
    Setting `snapshot_format = "hdf5"` raises an error. The old HDF5 reader
    targeted dataset names this fork never writes — the fork's `KZ(46)` writer
    produces `snap.40_*.h5part` files with `Step#i` groups and numbered
    datasets, a layout the removed reader could never parse. Re-adding support
    means porting to that layout (see git history for the old reader). HDF5
    files found in an output directory are still *detected* by `scan_output`
    and reported with a warning, but `conf.3` is the only readable snapshot
    source.

### `[visualization]`

| Key          | Type    | Default | Description |
|--------------|---------|---------|-------------|
| `enabled`    | Bool    | `true`  | Enable plot generation |
| `format`     | String  | `"png"` | Static plot format: `"png"`, `"pdf"`, `"svg"` (animations are always GIF) |
| `dpi`        | Int     | `300`   | Resolution for raster formats |
| `figsize`    | [Float] | `[8.0, 6.0]` | Figure size in inches `[width, height]` |
| `output_dir` | String  | `"plots"` | Plot directory, relative to each run directory |

### `[visualization.style]`

Presentation knobs collected in the `PlotStyle` struct (`cfg.visualization.style`). These are stylistic only — data-validity cutoffs (e.g. BSE placeholder filtering in HR diagrams) remain named constants in the plot code.

| Key                   | Type  | Default   | Description |
|-----------------------|-------|-----------|-------------|
| `marker_budget`       | Float | `18000.0` | Scatter marker size is `clamp(marker_budget/N, marker_min, marker_max)` |
| `marker_min`          | Float | `4.0`     | Lower clamp bound for the scatter marker size [pt] |
| `marker_max`          | Float | `20.0`    | Upper clamp bound for the scatter marker size [pt] |
| `q_log_threshold`     | Float | `10.0`    | Switch virial-ratio axes to log scale when `max(Q)` exceeds this |
| `q_floor`             | Float | `1e-3`    | Clamp floor for the virial ratio on *log-scale* axes only |
| `zoom_frac`           | Float | `0.15`    | Extent-ratio threshold for adaptive per-panel zoom in snapshot evolution plots and cluster animations |
| `anim_fps`            | Int   | `0`       | Animation frame rate; `0` selects automatically from frame count |
| `anim_target_seconds` | Float | `12.0`    | Target GIF duration used by the automatic FPS selection |

### `[merger]`

| Key           | Type   | Default | Description |
|---------------|--------|---------|-------------|
| `enabled`     | Bool   | `false` | Generate merger ICs before the simulation phase |
| `config_file` | String | `""`    | Path to the merger cluster TOML (e.g. `"input_files/merger_demo_small.toml"`); relative paths resolve against the package root |

The merger TOML schema itself (clusters, profiles, IMFs, orbit, output, seed) is documented in [Input File Reference](@ref) and [Multi-Cluster Merger Simulations](@ref).

---

## 5. Pipeline Modes

`run_pipeline` executes up to five phases, each controlled by config flags:

| Phase              | Controlled by           | Notes |
|:-------------------|:------------------------|:------|
| 1. Install/Build   | `install.enabled`       | Clone + compile Nbody6++ |
| 1.5 Merger ICs     | `merger.enabled`        | Generate `dat.10` + `merger.inp` into a new `runs/merger_<prefix>_.../output/` directory |
| 2. Simulation      | `simulation.run_test`   | Run the N-body integration (with merger ICs when phase 1.5 ran) |
| 3. Post-process    | `postprocess.enabled`   | Read output data + sanity checks |
| 4. Plots           | `visualization.enabled` | Static plots + GIF animations |

Common configurations (also listed in the `config.toml` header):

1. **Full pipeline:** `install.enabled=true`, `simulation.run_test=true`
2. **Simulate only:** `install.enabled=false`, `simulation.run_test=true`
3. **Postprocess only:** `simulation.run_test=false`, `postprocess.data_dir="/path/to/output"`
4. **Re-plot latest run:** `simulation.run_test=false`, `postprocess.data_dir=""` — the most recent run under `runs/` (plain or merger) is post-processed
5. **Merger ICs:** `merger.enabled=true`, `merger.config_file="input_files/..."` (with `run_test=false` to only generate)
6. **Merger + simulate:** `merger.enabled=true`, `simulation.run_test=true` — the binary runs inside the merger output directory so `dat.10` is found in its working directory

For output produced entirely outside this project there is also the config-free path:

```julia
results = postprocess_external("/scratch/sim42/output")   # scan + read + plot
```

---

## 6. Build Phase Details

### HDF5 build patching

The upstream `./configure --enable-hdf5` flag is broken. Nbody6Setup patches the build with a two-file mechanism:

1. Writes `hdf5_flags.mk` in the source root with the correct `-DCONFIG_HDF5`, include paths, and library flags
2. Appends `-include ../hdf5_flags.mk` to `build/Makefile`

This survives `./configure` re-runs: the flags file persists and only the one-line include needs reapplication. (This affects only what the *binary* can write — the Julia post-processing side reads `conf.3` snapshots, not HDF5; see the note in [Configuration Reference](#4-configuration-reference).)

### CUDA auto-detection

When `enable_gpu = true` and `cuda_path` is empty, the build searches in order:

1. Environment variables: `CUDA_HOME`, `CUDA_PATH`, `CUDA_ROOT`
2. Standard paths: `/usr/local/cuda`, `/usr/local/cuda-12`, `/usr/local/cuda-11.8`, …
3. `nvcc` location on `PATH`

### Fedora h5pfc workaround

On Fedora the HDF5 parallel Fortran wrapper `h5pfc` may carry an erroneous `/openmpi-x86_64` suffix in its `includedir`, preventing `hdf5.mod` from being found. Nbody6Setup detects this and prints the fix:

```bash
sudo sed -i 's|/openmpi-x86_64||g' /usr/lib64/openmpi/bin/h5pfc
```

---

## 7. Simulation Execution

### Run ID system

Each run gets a unique ID `{prefix}_YYYYMMDD_HHMMSS_{4hex}` (e.g. `run_20260325_143022_a1f3`); merger runs are prefixed `merger_` (e.g. `merger_run_20260325_143022_a1f3`). Run directory layout:

```
runs/run_20260325_143022_a1f3/
├── config.toml          # frozen snapshot of the configuration used
├── RUN_INFO.txt         # summary with timestamps, hostname, file listing
├── output/              # all simulation artefacts
│   ├── nbody6++         # binary copy (reproducibility)
│   ├── _launch.sh       # generated bash launch script
│   ├── out1000          # captured stdout
│   ├── err1000          # captured stderr
│   ├── conf.3_*         # particle snapshots
│   ├── lagr.7           # Lagrangian radii
│   ├── esc.11           # escaper events
│   └── sev.83_*         # stellar evolution snapshots
└── plots/               # post-processing plots & GIF animations
```

Merger runs additionally contain `dat.10`, `merger.inp`, `merger_summary.txt`, and `merger_ic.toml` in `output/`.

### Launch script

The generated `_launch.sh` sets `ulimit -s unlimited` (Fortran stack), `OMP_STACKSIZE=4096M`, CUDA environment variables (if GPU enabled), and `stdbuf -oL` for line-buffered output where available.

### Real-time monitoring

During execution, ADJUST summaries are echoed live:

```
[ Info:   t_NB=0.0500  t_Myr=0.4  N=9998  |ΔE/E|=1.23e-06  Q_vir=0.987
```

---

## 8. Post-processing

### Snapshots (conf.3, Fortran binary)

```julia
snap = read_conf3("path/to/conf.3_0")
snap.header        # SnapshotHeader with AS(1:20) parameters
snap.pos           # 3 × N matrix (Float32)
snap.vel           # 3 × N matrix (Float32)
snap.mass          # N-element vector (Float32)
snap.name          # particle identifiers (Int32)
snap.rho, snap.phi # local density / potential (empty if unavailable)
nparticles(snap)

snaps = read_all_conf3(dir, "conf.3_*")   # numerically sorted by time suffix
```

The reader auto-detects the bulk-array format (one record with AS params + all particle arrays) and the legacy per-particle record format (32- or 44-byte records; the extended form carries `rho`/`phi`). Corrupt files are skipped with a warning.

### Diagnostics (stdout)

```julia
diag = read_diagnostics("out1000")
diag.adjust              # Vector{AdjustRecord}: time, Q, ΔE/E, E_tot, N, N_pairs, R_scale
diag.physical_scaling    # Dict with R*, M*, V*, T*, <M>, …

units = extract_scaling(diag)   # → UnitScaling
```

Each ADJUST epoch is merged from up to three stdout lines (`ADJUST:`, `RMIN/RSCALE`, `TIME[NB]`); both positional and key-value ADJUST formats are handled, and particle counts are forward-filled across epochs where `TIME[NB]` lines are sparse.

### Lagrangian radii

```julia
lagr = read_lagr("lagr.7")
lagr.time                # time values [NB]
lagr.mass_fractions      # [0.001, 0.003, …, 1.0] (18 standard fractions)
lagr.radii               # n_fractions × n_times matrix
```

Auto-detects the modern single-row format (with `##`/`TIME` headers) and the legacy block format.

### Escapers

```julia
escs = read_escapers("esc.11")
escs[1].time_myr         # escape time [Myr]
escs[1].mass_solar       # mass [M☉]
escs[1].velocity_kms     # escape velocity [km/s]
escs[1].stellar_type     # K* type (Hurley codes, see below)
```

The physical-unit quantities are taken from tokens 6–11 of each line (the first five columns are NB-unit diagnostics).

### Stellar evolution

```julia
sev = read_stellar_evolution("sev.83_0")
sev.time_myr             # header TPHYS [Myr]
sev.records              # Vector{StellarRecord}
sev.records[1].log_teff
sev.records[1].log_luminosity

sevs = read_all_stellar_evolution(run_out, "sev.83_*")
```

Note the header time is TPHYS in **Myr** while each data line's first token is TTOT in **NB units**; both clocks are kept (`StellarEvolutionSnapshot.time_myr`, `StellarRecord.time_nb`).

---

## 9. Visualisation

All plots use the built-in publication theme, activated globally with `set_publication_theme!()` (called automatically by the pipeline).

### Theme rules

- **No plot titles** — contextual information (time, projection label, orbit parameters) is placed as in-axis annotations at the top-left corner
- **No minor ticks**; major ticks face inward; fully boxed axes
- **Grid**: dashed grey at very low opacity — on *line plots only*. Dense scatter plots (cluster projections, HR diagrams) disable the grid locally
- **Fonts**: Computer Modern (NewComputerModern, bundled with MathTeXEngine) for text and math
- **Never-overwrite policy**: before saving, any existing file at the target path is moved to the first free `name#k.ext` sibling (DrWatson-`safesave` style); this applies to plots, animations, `dat.10`, `merger.inp`, and metadata files alike
- Legends appear only when ≥ 2 items are plotted; axis limits snap to nice round values with a visual buffer

### Plot inventory

Static plots take the extension from `visualization.format`; animations are always `.gif`.

| Function | Description | Output files |
|----------|-------------|--------------|
| `plot_snapshot(snap, vis)` | Scatter projections, viridis mass colouring `log₁₀(m/M_tot)` | `snapshot_final_{xy,xz,yz}` |
| `plot_snapshot_evolution(snaps, vis)` | Up to 6 panels, shared colorbar; adaptive per-panel zoom when extents vary by more than `1/zoom_frac` | `snapshot_evolution_{xy,xz,yz}` |
| `plot_energy(diag, vis)` | Two panels: log `\|ΔE/E\|` and virial ratio with `Q = 0.5` reference | `energy` |
| `plot_particle_count(diag, vis)` | Two panels: bound N and KS pairs, integer ticks | `particle_count` |
| `plot_lagrangian(lagr, vis)` | Selected mass-fraction radii vs time, log y | `lagrangian_radii` |
| `plot_hr(sev, vis)` | HR diagram coloured by stellar type K*, reversed Teff axis | `hr_diagram_{early,mid,final}` |
| `plot_hr_evolution(sevs, vis)` | HR panel grid across up to 6 epochs | `hr_evolution` |
| `plot_cluster_separation(snaps, ranges, vis)` | Pairwise COM separations of the initial clusters; per-pair lines for ≤ 5 clusters, min/max/mean envelope + surviving-cluster staircase and coalescence marker otherwise | `merger_cluster_separation` |
| `plot_cluster_virial(snaps, ranges, vis)` | Internal virial ratio `Q_i(t)` per initial cluster (COM-subtracted, self-gravity only); log axis when `max(Q) > q_log_threshold` | `merger_cluster_virial` |
| `plot_merger_ic(result, vis)` | IC diagnostics: projections, 3-panel overview, velocity quiver by cluster, IMF histogram with Kroupa reference slopes, per-cluster radial density | `merger_ic_{xy,xz,yz}`, `merger_ic_overview`, `merger_ic_velocity`, `merger_ic_imf`, `merger_ic_density` |
| `animate_cluster(snaps, vis)` | Animated scatter per projection, global or adaptive limits | `cluster_evolution_{xy,xz,yz}.gif` |
| `animate_lagrangian(lagr, vis)` | Progressive line draw with ghost background and time cursor | `lagrangian_anim.gif` |
| `animate_hr(sevs, vis)` | Animated HR diagram | `hr_evolution_anim.gif` |

The merger plots (`plot_cluster_separation`, `plot_cluster_virial`) are generated automatically by `generate_plots`/`postprocess_external` when a `merger_summary.txt` is found next to the snapshots; the cluster index ranges come from `parse_merger_summary`. `plot_merger_ic` runs whenever a merger IC was generated in the same pipeline invocation (or can be re-run later via `load_merger_ic_result`).

### Customising plots

All plot functions accept a `VisualizationConfig` and a `filename` keyword:

```julia
vis = VisualizationConfig(format = "pdf", dpi = 600, figsize = (10.0, 8.0),
                          output_dir = "my_plots", style = Nbody6Setup.PlotStyle())
plot_snapshot(snap, vis; filename = "cluster_final")
```

---

## 10. Using the Julia API Directly

```julia
using Nbody6Setup

cfg = load_config("config.toml")

run_out = "runs/run_20260325_143022_a1f3/output"
snap = read_conf3(joinpath(run_out, "conf.3_0"))
diag = read_diagnostics(joinpath(run_out, "out1000"))
lagr = read_lagr(joinpath(run_out, "lagr.7"))
escs = read_escapers(joinpath(run_out, "esc.11"))
sevs = read_all_stellar_evolution(run_out)

units = extract_scaling(diag)
pos_pc = snap.pos .* units.rbar          # positions in parsecs

set_publication_theme!()
plot_hr(sevs[end], cfg.visualization; filename = "my_hr")

# Full post-processing for a specific run
run_dir = "runs/run_20260325_143022_a1f3"
results = postprocess(cfg; run_dir)
generate_plots(results, cfg; run_dir)
```

### Key exported types

| Type | Description |
|------|-------------|
| `Nbody6Config` | Top-level config (install, build, simulation, postprocess, visualization, merger) |
| `Snapshot` / `SnapshotHeader` | Particle data from conf.3 + AS(1:20) header |
| `DiagnosticsData` / `AdjustRecord` | Parsed ADJUST records + physical scaling |
| `LagrangianData` | Lagrangian radii evolution |
| `EscaperRecord` | Single escaper event |
| `StellarRecord` / `StellarEvolutionSnapshot` | Stellar evolution data from one sev.83 file |
| `UnitScaling` | Physical unit conversion factors |
| `MergerConfig` / `MergerICResult` | Merger IC specification and structured result |
| `OutputScan` | Result of `scan_output` on an external directory |

---

## 11. Output File Reference

Files produced by the Nbody6++ simulation and read by this package:

| File | Description | Reader |
|------|-------------|--------|
| `conf.3_*` | Particle snapshots (Fortran binary) | `read_conf3` / `read_all_conf3` |
| `out1000` | Simulation stdout (ADJUST lines, scaling) | `read_diagnostics` |
| `err1000` | Simulation stderr | — |
| `lagr.7` | Lagrangian radii (ASCII) | `read_lagr` |
| `esc.11` | Escaper events (ASCII) | `read_escapers` |
| `sev.83_*` | Single-star evolution snapshots (ASCII) | `read_stellar_evolution` |
| `snap.40_*.h5part` | HDF5 snapshots (KZ(46)) | **not readable** — detected and warned about only |

### Stellar type codes (K*)

Hurley et al. (2000) convention as used by this fork (`STELLAR_TYPE_LABELS`):

| K* | Type | K* | Type |
|----|------|----|------|
| 0 | MS (low-mass, M < 0.7) | 8 | HeHG (He Hertzsprung Gap) |
| 1 | MS (Main Sequence) | 9 | HeGB (He Giant Branch) |
| 2 | HG (Hertzsprung Gap) | 10 | HeWD (He White Dwarf) |
| 3 | GB (Giant Branch) | 11 | COWD (CO White Dwarf) |
| 4 | CHeB (Core He Burning) | 12 | ONeWD (ONe White Dwarf) |
| 5 | EAGB (Early AGB) | 13 | NS (Neutron Star) |
| 6 | TPAGB (Thermally Pulsing AGB) | 14 | BH (Black Hole) |
| 7 | HeMS (Naked He MS) | 15 | SNR (Massless Remnant) |

---

## 12. N-body Units and Conversions

Nbody6++ uses Hénon N-body units internally: `G = 1`, `M_total = 1`, `E_total = -1/4`. The virial ratio is `Q = T/|W|` with equilibrium at `Q = 0.5`.

Physical conversions come from the scaling factors printed in the `PHYSICAL SCALING` stdout block and collected into a `UnitScaling`:

| Field    | Stdout symbol | Converts | Physical unit |
|----------|---------------|----------|---------------|
| `rbar`   | `R*`          | NB length → parsecs | pc |
| `zmbar`  | `M*`          | NB mass → solar masses | M☉ |
| `tscale` | `T*`          | NB time → megayears | Myr |
| `vstar`  | `V*`          | NB velocity → km/s | km/s |

!!! warning "ZMBAR is the total-mass scale, not the mean mass"
    Nbody6++ redefines `ZMBAR` at startup (`start.F`) as the **total-mass**
    scale factor — 1 NB mass unit (the whole cluster) = `zmbar` M☉ — and
    `units.f` converts masses as `ZMBAR*M`. The mean stellar mass is printed
    separately as `<M>` in the `PHYSICAL SCALING` line. `extract_scaling`
    therefore reads `M*`, not `<M>`.

```julia
units = extract_scaling(diag)
r_pc  = r_nb * units.rbar
t_myr = t_nb * units.tscale
v_kms = v_nb * units.vstar
m_sun = m_nb * units.zmbar
```

The same factors are available per snapshot from the conf.3 header via the accessors `rbar(h)`, `zmbar(h)`, `tscale(h)`, `vstar(h)`, plus `time_nb(h)`, `time_myr(h)`, `rscale(h)` (half-mass radius, NB), and `rc(h)` (core radius, NB).

---

## 13. Troubleshooting

### Build fails with "hdf5.mod not found"

**Fedora**: the `h5pfc` wrapper likely has a broken includedir; see [Build Phase Details](#6-build-phase-details).
**Ubuntu**: ensure `libhdf5-dev` (or `libhdf5-openmpi-dev` with MPI) is installed.

### Simulation segfaults immediately

Fortran code requires a large stack. The launch script sets `ulimit -s unlimited`; if running manually, set this in your shell first.

### "CUDA not found" during build

Set `cuda_path` explicitly in `[build]`, or export one of `CUDA_HOME`, `CUDA_PATH`, `CUDA_ROOT`.

### Post-processing finds no files

Check `simulation.runs_dir` and, for external data, `postprocess.data_dir`. For a specific run:

```julia
results = postprocess(cfg; run_dir = "runs/run_20260325_143022_a1f3")
generate_plots(results, cfg; run_dir = "runs/run_20260325_143022_a1f3")
```

### HDF5 snapshots present but nothing is read

Expected: the `KZ(46)` H5Part layout is not supported (see the note in [Configuration Reference](#4-configuration-reference)). Enable `conf.3` output in the `.inp` file (`KZ(3)` ≥ 1) and use `snapshot_format = "conf3"`.
