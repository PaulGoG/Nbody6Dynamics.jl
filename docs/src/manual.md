# Nbody6Setup.jl — User Manual

## Table of Contents

1. [Introduction](#1-introduction)
2. [Prerequisites](#2-prerequisites)
3. [Installation](#3-installation)
4. [Configuration Reference](#4-configuration-reference)
5. [Running the Pipeline](#5-running-the-pipeline)
6. [Build Phase Details](#6-build-phase-details)
7. [Simulation Execution](#7-simulation-execution)
8. [Post-processing](#8-post-processing)
9. [Visualisation](#9-visualisation)
10. [Using the Julia API Directly](#10-using-the-julia-api-directly)
11. [Output File Reference](#11-output-file-reference)
12. [N-body Units and Conversions](#12-n-body-units-and-conversions)
13. [Troubleshooting](#13-troubleshooting)
14. [FAQ](#14-faq)

---

## 1. Introduction

**Nbody6Setup.jl** automates the full lifecycle of running an N-body star cluster simulation with [Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing):

1. **Download and compile** the Fortran/C++ simulation code
2. **Execute** the simulation with proper environment setup
3. **Parse** all output files into Julia data structures
4. **Generate** publication-quality plots

Everything is driven by a single `config.toml` file, making it easy to reproduce runs and share configurations.

---

## 2. Prerequisites

### System requirements

- **Linux** (tested on Fedora 42+ and Ubuntu 20.04+)
- **Julia 1.10+**
- **GCC toolchain**: `gcc`, `g++`, `gfortran`, `make`
- **Git** (for cloning the simulation code)

### Optional (depending on config)

| Feature | Requires |
|---------|----------|
| MPI support | `mpicc`, `mpif90`, `mpirun` (OpenMPI or MPICH) |
| GPU support | CUDA toolkit with `nvcc` |
| HDF5 output | `libhdf5-dev` (Ubuntu) or `hdf5-devel` (Fedora) |

### Installing system dependencies

**Fedora:**
```bash
sudo dnf install gcc gcc-c++ gcc-gfortran make git
sudo dnf install hdf5-devel                        # HDF5 support
sudo dnf install openmpi-devel                      # MPI support (optional)
```

**Ubuntu/Debian:**
```bash
sudo apt install gcc g++ gfortran make git
sudo apt install libhdf5-dev                        # HDF5 support
sudo apt install libopenmpi-dev openmpi-bin          # MPI support (optional)
```

---

## 3. Installation

```bash
# Clone or navigate to the project
cd Nbody6Setup

# Install Julia package dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Verify with tests
julia --project=. -e 'using Pkg; Pkg.test()'
```

---

## 4. Configuration Reference

The `config.toml` file controls all phases. Here is a complete reference:

### `[install]`

| Key           | Type   | Default | Description |
|---------------|--------|---------|-------------|
| `enabled`     | Bool   | `true`  | Enable the install/build phase |
| `source_url`  | String | GitHub URL | Git repository to clone |
| `install_dir` | String | `"backend/Nbody6PPGPU-beijing"` | Local directory for the source (inside `backend/`) |
| `reinstall`   | Bool   | `false` | Delete and re-clone if `true` |
| `clean_build` | Bool   | `true`  | Run `make clean` before building |

### `[build]`

| Key               | Type     | Default | Description |
|-------------------|----------|---------|-------------|
| `configure_flags` | [String] | `["--enable-mcmodel=large", "--with-par=b1m"]` | Flags passed to `./configure` |
| `enable_mpi`      | Bool     | `false` | Enable MPI parallelism |
| `enable_hdf5`     | Bool     | `true`  | Patch Makefile for HDF5 output |
| `enable_gpu`      | Bool     | `false` | Enable GPU acceleration (requires CUDA) |
| `cuda_path`       | String   | `""`    | CUDA installation path; empty = auto-detect |
| `nproc`           | Int      | `0`     | Parallel `make` jobs; 0 = auto-detect |

### `[simulation]`

| Key             | Type   | Default | Description |
|-----------------|--------|---------|-------------|
| `run_test`      | Bool   | `true`  | Run a test simulation |
| `input_file`    | String | `"examples/input_files/N10k_noDat10.inp"` | Path to input file (relative to source dir) |
| `runs_dir`      | String | `"runs"` | Base directory for run output |
| `binary_name`   | String | `"nbody6++"` | Expected binary name |
| `mpi_ranks`     | Int    | `1`     | Number of MPI ranks (when MPI enabled) |
| `run_id_prefix` | String | `"run"` | Prefix for run directory names |

### `[postprocess]`

| Key                   | Type   | Default | Description |
|-----------------------|--------|---------|-------------|
| `enabled`             | Bool   | `true`  | Enable post-processing |
| `snapshot_format`     | String | `"conf3"` | `"conf3"` or `"hdf5"` |
| `snapshot_pattern`    | String | `"conf.3_*"` | Glob pattern for conf.3 files |
| `hdf5_file`           | String | `"data.40.h5part"` | HDF5 snapshot file name |
| `parse_stdout`        | Bool   | `true`  | Parse simulation stdout |
| `stdout_file`         | String | `"out1000"` | Name of the stdout capture file |
| `read_lagr`           | Bool   | `true`  | Read Lagrangian radii |
| `lagr_file`           | String | `"lagr.7"` | Lagrangian radii file name |
| `read_escapers`       | Bool   | `true`  | Read escaper data |
| `escapers_file`       | String | `"esc.11"` | Escaper file name |
| `read_stellar_evo`    | Bool   | `true`  | Read stellar evolution snapshots |
| `stellar_evo_pattern` | String | `"sev*.83"` | Glob pattern for sev files |

### `[visualization]`

| Key          | Type   | Default | Description |
|--------------|--------|---------|-------------|
| `enabled`    | Bool   | `true`  | Enable plot generation |
| `format`     | String | `"png"` | Output format: `"png"`, `"pdf"`, `"svg"` |
| `dpi`        | Int    | `300`   | Resolution for raster formats |
| `figsize`    | [Int]  | `[10, 8]` | Figure size in inches `[width, height]` |
| `output_dir` | String | `"plots"` | Directory for saved plots (relative to each run) |

---

## 5. Running the Pipeline

### Full pipeline via script

```bash
julia --project=. scripts/run_setup.jl [config.toml]
```

The script executes four phases in sequence:
1. **Install & Build** (if `install.enabled = true`)
2. **Simulation** (if `simulation.run_test = true`)
3. **Post-processing** (if `postprocess.enabled = true`)
4. **Visualisation** (if `visualization.enabled = true` and data exists)

### Skipping phases

Set `enabled = false` or `run_test = false` in the relevant section. Common patterns:

**Build only (no simulation):**
```toml
[simulation]
run_test = false

[postprocess]
enabled = false
```

**Post-process existing output (no build, no simulation):**
```toml
[install]
enabled = false

[simulation]
run_test = false
```

---

## 6. Build Phase Details

### How HDF5 patching works

The upstream `./configure --enable-hdf5` flag is broken. Nbody6Setup patches the build using a robust two-file mechanism:

1. Writes `hdf5_flags.mk` in the source root with the correct `-DCONFIG_HDF5`, include paths, and library flags
2. Appends `-include ../hdf5_flags.mk` to `build/Makefile`

This design survives `./configure` re-runs: the flags file persists, and only the one-line include injection needs reapplication.

### CUDA auto-detection

When `enable_gpu = true` and `cuda_path` is empty, the build phase searches for CUDA in order:
1. Environment variables: `CUDA_HOME`, `CUDA_PATH`, `CUDA_ROOT`
2. Standard paths: `/usr/local/cuda`, `/usr/local/cuda-12`, `/usr/local/cuda-11.8`, etc.
3. `nvcc` location on PATH

### Fedora h5pfc workaround

On Fedora, the HDF5 parallel Fortran wrapper `h5pfc` has a known bug where the `includedir` contains an erroneous `/openmpi-x86_64` suffix, preventing `hdf5.mod` from being found. Nbody6Setup detects this and prints a warning with the fix:

```bash
sudo sed -i 's|/openmpi-x86_64||g' /usr/lib64/openmpi/bin/h5pfc
```

---

## 7. Simulation Execution

### Run ID system

Each simulation run gets a unique ID: `{prefix}_YYYYMMDD_HHMMSS_{4hex}`, e.g. `run_20260325_143022_a1f3`.

The run directory (`runs/run_20260325_143022_a1f3/`) has the following structure:

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
│   └── sev*.83          # stellar evolution snapshots
└── plots/               # post-processing plots & GIF animations
```

### Launch script

The generated `_launch.sh` script sets:
- `ulimit -s unlimited` (required for Fortran stack-heavy code)
- `OMP_STACKSIZE=4096M`
- CUDA environment variables (if GPU enabled)
- `stdbuf -oL` for line-buffered output (if available)

### Real-time monitoring

During execution, the terminal displays ADJUST summaries in real time:

```
[ Info:   t_NB=0.0500  t_Myr=0.4  N=9998  |ΔE/E|=1.23e-06  Q_vir=0.987
[ Info:   t_NB=0.1000  t_Myr=0.8  N=9995  |ΔE/E|=2.45e-06  Q_vir=0.991
```

---

## 8. Post-processing

### Snapshot readers

**conf.3 (Fortran binary):**
```julia
snap = read_conf3("path/to/conf.3_0")
snap.header        # SnapshotHeader with 20 parameters
snap.pos           # 3 × N matrix (Float32)
snap.vel           # 3 × N matrix (Float32)
snap.mass          # N-element vector (Float32)
snap.name          # particle identifiers (Int32)
nparticles(snap)   # number of particles
```

Auto-detects standard (32-byte) vs extended (44-byte) per-particle records. Extended format includes density (`snap.rho`) and potential (`snap.phi`).

**HDF5/H5Part:**
```julia
snaps = read_hdf5_snapshots("data.40.h5part")
```

Reads `/Step#N/` groups with flexible field name detection (handles case variations).

### Diagnostics

```julia
diag = read_diagnostics("out1000")
diag.adjust              # Vector{AdjustRecord}
diag.physical_scaling    # Dict with R*, T*, V*, <M>, etc.

units = extract_scaling(diag)
units.rbar               # pc per NB length
units.tscale             # Myr per NB time
```

### Lagrangian radii

```julia
lagr = read_lagr("lagr.7")
lagr.time                # time values
lagr.mass_fractions      # [0.001, 0.003, ..., 1.0]
lagr.radii               # n_fractions × n_times matrix
```

### Escapers

```julia
escs = read_escapers("esc.11")
escs[1].time_myr         # escape time [Myr]
escs[1].mass_solar       # mass [M☉]
escs[1].velocity_kms     # escape velocity [km/s]
escs[1].stellar_type     # K* type (0=MS, ..., 13=BH)
```

### Stellar evolution

```julia
sev = read_stellar_evolution("sev001.83")
sev.time_myr             # NB time of snapshot
sev.records              # Vector{StellarRecord}
sev.records[1].log_teff  # log10(Teff/K)
sev.records[1].log_luminosity  # log10(L/L☉)

# Batch read all files
sevs = read_all_stellar_evolution("run_dir/", "sev*.83")
```

---

## 9. Visualisation

All plots use the built-in publication theme (fully-boxed axes, appropriate font sizes, 300 DPI default). Activate it with:

```julia
set_publication_theme!()
```

### Available plots

| Function | Description | Output |
|----------|-------------|--------|
| `plot_snapshot(snap, vis)` | XY/XZ/YZ scatter plots with mass colorbar | `snapshot_final_{xy,xz,yz}.png` |
| `plot_snapshot_evolution(snaps, vis)` | Multi-panel time sequence with mass colorbar | `snapshot_evolution_{xy,xz,yz}.png` |
| `plot_energy(diag, vis)` | Energy error + virial ratio | `energy.png` |
| `plot_particle_count(diag, vis)` | N(t) evolution with integer ticks | `particle_count.png` |
| `plot_lagrangian(lagr, vis)` | Lagrangian radii vs time | `lagrangian_radii.png` |
| `plot_hr(sev, vis)` | HR diagram (early/mid/final epochs) | `hr_diagram_{early,mid,final}.png` |
| `plot_hr_evolution(sevs, vis)` | HR diagram panels across 6 epochs | `hr_evolution.png` |
| `animate_cluster(snaps, vis)` | Animated cluster scatter in XY/XZ/YZ | `cluster_evolution_{xy,xz,yz}.gif` |
| `animate_lagrangian(lagr, vis)` | Lagrangian radii progressive draw | `lagrangian_anim.gif` |
| `animate_hr(sevs, vis)` | Animated HR diagram evolution | `hr_evolution_anim.gif` |

### Plot legends and visual encoding

Each plot uses specific legends, colorbars, and visual encodings. This section documents what every symbol, colour, and label means across all generated plots.

#### Snapshot projections (`plot_snapshot`)

- **Scatter markers**: each point is one particle (star or compact object)
- **Colorbar**: `log₁₀(m / M_tot)` — logarithmic particle mass normalised to total cluster mass. Colour map: *viridis* (dark purple = low mass, yellow = high mass)
- **Axes**: spatial coordinates in N-body units `[NB]`
- **Title**: particle count `N` and simulation time `t [NB]`
- **Marker size**: auto-scaled as `clamp(18000/N, 4, 20)` — visible but non-overlapping
- **Projections**: XY, XZ, and YZ (one file per projection)

#### Snapshot evolution (`plot_snapshot_evolution`)

- Up to 6 panels showing the cluster at evenly-spaced epochs (one file per projection: XY, XZ, YZ)
- **Scatter markers**: coloured by `log₁₀(m / M_tot)` using global viridis colormap with shared colorbar
- **Consistent axis limits** across all panels for direct visual comparison (nice round values ± visual buffer)
- **Panel titles**: simulation time `t [NB]`

#### Energy diagnostics (`plot_energy`)

Two vertically-stacked panels with linked time axes:

**Top panel — Energy conservation:**
- **Blue line** (`Relative energy error`): absolute relative energy error `|ΔE/E|` on a log₁₀ y-axis. Values below ~10⁻⁶ indicate excellent energy conservation; values approaching 10⁻² signal numerical problems.

**Bottom panel — Virial ratio:**
- **Red line** (`Q = T/|W| (virial ratio)`): ratio of kinetic to absolute potential energy as reported by Nbody6++. For a virialised system, `Q ≈ 0.5`.
- **Grey dashed line** (`Q = 0.5 (virial equilibrium)`): reference at the virial equilibrium value. Persistent deviations indicate the cluster is not in virial equilibrium (e.g., during core collapse, tidal stripping, or energetic binary ejections).

#### Particle count (`plot_particle_count`)

- **Blue line** (`N (bound particles)`): total number of gravitationally bound particles vs time. Decreases as stars escape the cluster.
- **Orange line** (`N_pairs (KS regularised binaries)`): number of active Kustaanheimo–Stiefel regularised binary pairs. Reflects binary formation, hardening, and disruption activity.

#### Lagrangian radii (`plot_lagrangian`)

- **Coloured lines**: each line tracks the radius enclosing a specific mass fraction of the cluster. Legend labels are of the form `M(r)/M_tot = X%`, e.g.:
  - `M(r)/M_tot = 1%` — core radius (innermost 1% of mass)
  - `M(r)/M_tot = 10%` — inner halo
  - `M(r)/M_tot = 50%` — half-mass radius
  - `M(r)/M_tot = 90%` — outer envelope
  - `M(r)/M_tot = 100%` — total cluster extent (tidal radius)
- **y-axis**: log₁₀ scale in N-body length units `[NB]`
- **Colours**: Wong colour palette (colourblind-safe)
- **Physical interpretation**: contraction of inner radii signals core collapse; expansion of outer radii indicates tidal mass loss.

#### HR diagram (`plot_hr`)

- **Scatter markers**: each point is one star; x-axis is reversed (hot stars on the left)
- **Colour by stellar type (K\*)**: each BSE stellar type gets a distinct colour:

| Colour | Stellar type |
|--------|-------------|
| Royal blue | MS (Main Sequence, K\*=0) |
| Dodger blue | HG (Hertzsprung Gap, K\*=1) |
| Orange | GB (Giant Branch, K\*=2) |
| Gold | CHeB (Core He Burning, K\*=3) |
| Orange-red | AGB (Asymptotic Giant Branch, K\*=4) |
| Red | EAGB (Early AGB, K\*=5) |
| Cyan | HeStar (Helium Star, K\*=6) |
| Teal | HeHG (He Hertzsprung Gap, K\*=7) |
| Olive | HeGB (He Giant Branch, K\*=8) |
| Light grey | HeWD (Helium White Dwarf, K\*=9) |
| Silver | COWD (CO White Dwarf, K\*=10) |
| Grey | ONeWD (ONe White Dwarf, K\*=11) |
| Purple | NS (Neutron Star, K\*=12) |
| Black | BH (Black Hole, K\*=13) |

- **Legend**: shown when ≤12 stellar types are present (arranged in 2 columns)
- **Marker size**: 10 px for all types

#### HR evolution (`plot_hr_evolution`)

- Same colour encoding as `plot_hr`, applied across up to 6 time-spaced panels
- **Consistent axis limits** across all panels
- **No per-panel legend** (colour mapping is identical to single HR diagram above)

#### Lagrangian animation (`animate_lagrangian`)

- **Ghost lines** (light grey): full time extent of each Lagrangian radius, drawn as background
- **Coloured lines**: progressively revealed from left to right, matching the static Lagrangian plot colours and labels
- **Dashed vertical cursor**: marks the current animation time
- **Frame rate**: auto-calculated for ~12 s total duration (2–30 fps), or manually overridden

#### Cluster animation (`animate_cluster`)

- Same encoding as `plot_snapshot`: viridis mass colouring with global colour range
- **Global axis limits**: computed across all frames for stable viewport
- **Title**: updates per frame with `N` and `t [NB]`
- **Frame rate**: auto-calculated for ~12 s total duration (1–10 fps)

#### HR animation (`animate_hr`)

- Same stellar-type colour encoding as `plot_hr`
- **Global axis limits**: computed across all epochs
- **Title**: updates per frame with `t_NB` and `N_star`
- **Frame rate**: auto-calculated for ~15 s total duration (1–8 fps); slower pace because HR frames are information-dense

### Customising plots

All plot functions accept a `VisualizationConfig` and a `filename` keyword:

```julia
vis = VisualizationConfig(format="pdf", dpi=600, figsize=(10, 8), output_dir="my_plots")
plot_snapshot(snap, vis; filename="cluster_final")
```

---

## 10. Using the Julia API Directly

For interactive analysis or custom workflows, use the package directly from the Julia REPL:

```julia
using Nbody6Setup

# Load config
cfg = load_config("config.toml")

# Read specific files from a run's output directory
run_out = "runs/run_20260325_143022_a1f3/output"
snap = read_conf3(joinpath(run_out, "conf.3_0"))
diag = read_diagnostics(joinpath(run_out, "out1000"))
lagr = read_lagr(joinpath(run_out, "lagr.7"))
escs = read_escapers(joinpath(run_out, "esc.11"))
sevs = read_all_stellar_evolution(run_out)

# Extract unit conversions
units = extract_scaling(diag)
pos_pc = snap.pos .* units.rbar  # positions in parsecs

# Generate specific plots
set_publication_theme!()
vis = cfg.visualization
plot_hr(sevs[end], vis; filename="my_hr")
plot_energy(diag, vis; filename="my_energy")

# Or run full post-processing for a specific run
run_dir = "runs/run_20260325_143022_a1f3"
results = postprocess(cfg; run_dir, base_dir=pwd())
generate_plots(results, cfg; run_dir)
```

### Key exported types

| Type | Description |
|------|-------------|
| `Nbody6Config` | Top-level config (install, build, simulation, postprocess, visualization) |
| `Snapshot` | Particle data from conf.3 or HDF5 |
| `SnapshotHeader` | Header metadata with 20 AS parameters |
| `DiagnosticsData` | Parsed ADJUST records + physical scaling |
| `LagrangianData` | Lagrangian radii evolution |
| `EscaperRecord` | Single escaper event |
| `StellarRecord` | Single star properties at one epoch |
| `StellarEvolutionSnapshot` | All stars from one sev*.83 file |
| `UnitScaling` | Physical unit conversion factors |

---

## 11. Output File Reference

Files produced by the Nbody6++ simulation:

| File | Description | Reader |
|------|-------------|--------|
| `conf.3_*` | Particle snapshots (Fortran binary) | `read_conf3` |
| `data.40.h5part` | Particle snapshots (HDF5) | `read_hdf5_snapshots` |
| `out1000` | Simulation stdout (ADJUST lines, scaling) | `read_diagnostics` |
| `err1000` | Simulation stderr | — |
| `lagr.7` | Lagrangian radii (ASCII blocks) | `read_lagr` |
| `esc.11` | Escaper events (ASCII) | `read_escapers` |
| `sev*.83` | Single-star evolution snapshots (ASCII) | `read_stellar_evolution` |

### Stellar type codes (K*)

| K* | Type | K* | Type |
|----|------|----|------|
| 0 | Main Sequence | 7 | He Hertzsprung Gap |
| 1 | Hertzsprung Gap | 8 | He Giant Branch |
| 2 | Giant Branch | 9 | He White Dwarf |
| 3 | Core He Burning | 10 | CO White Dwarf |
| 4 | AGB | 11 | ONe White Dwarf |
| 5 | Early AGB | 12 | Neutron Star |
| 6 | He Star | 13 | Black Hole |

---

## 12. N-body Units and Conversions

Nbody6++ uses Heggie & Mathieu (1986) N-body units where G = 1, M_total = 1, and E_total = -1/4. Physical conversions are given by scaling factors printed during initialisation:

| Symbol | Converts | Physical unit |
|--------|----------|---------------|
| `RBAR` (R*) | NB length → parsecs | pc |
| `ZMBAR` (<M>) | NB mass → solar masses | M☉ |
| `TSCALE` (T*) | NB time → megayears | Myr |
| `VSTAR` (V*) | NB velocity → km/s | km/s |

These are automatically extracted by `extract_scaling(diag)` from the `PHYSICAL SCALING` block in stdout.

```julia
units = extract_scaling(diag)
r_physical = r_nb * units.rbar    # NB → pc
t_physical = t_nb * units.tscale  # NB → Myr
v_physical = v_nb * units.vstar   # NB → km/s
m_physical = m_nb * units.zmbar   # NB → M☉
```

Convenience functions:
```julia
to_pc(units, r_nb)
to_myr(units, t_nb)
to_kms(units, v_nb)
to_msun(units, m_nb)
```

---

## 13. Troubleshooting

### Build fails with "hdf5.mod not found"

**Fedora**: The `h5pfc` wrapper likely has a broken includedir. Run:
```bash
sudo sed -i 's|/openmpi-x86_64||g' /usr/lib64/openmpi/bin/h5pfc
```

**Ubuntu**: Ensure `libhdf5-dev` is installed. If using MPI, install `libhdf5-openmpi-dev`.

### Simulation segfaults immediately

Fortran code requires a large stack. The launch script sets `ulimit -s unlimited`, but if running manually ensure this is set in your shell before launching.

### "CUDA not found" warning during build

Set `cuda_path` explicitly in `config.toml`:
```toml
[build]
cuda_path = "/usr/local/cuda-12"
```

Or ensure one of `CUDA_HOME`, `CUDA_PATH`, or `CUDA_ROOT` is set in your environment.

### Post-processing finds no files

Ensure `runs_dir` in `[simulation]` points to the correct location. When using the script pipeline, all paths are relative to the directory containing `config.toml`.

For post-processing a specific run, pass the `run_dir` keyword:
```julia
results = postprocess(cfg; run_dir="runs/run_20260325_143022_a1f3")
generate_plots(results, cfg; run_dir="runs/run_20260325_143022_a1f3")
```

### MPI: "all ports busy" or binding errors

The launch script uses `mpirun --bind-to none`. If you see port conflicts, reduce `mpi_ranks` or ensure OpenMPI is properly configured:
```bash
# Fedora: load the MPI module
module load mpi/openmpi-x86_64
```

---

## 14. FAQ

**Q: Can I use this without GPU support?**
A: Yes. Set `enable_gpu = false` in `[build]`. This is the default.

**Q: Can I post-process output from a simulation run elsewhere?**
A: Yes. Set `install.enabled = false` and `simulation.run_test = false`, then use `postprocess(cfg; run_dir="runs/<run_id>")` to point at the directory containing the output files.

**Q: How do I change the number of particles?**
A: Edit the input file specified by `simulation.input_file`. The default `N10k_noDat10.inp` runs 10,000 particles. Smaller test files (e.g., 1000 particles) are available in `examples/input_files/`.

**Q: Can I produce PDF plots for a paper?**
A: Set `format = "pdf"` in `[visualization]`. SVG is also available.

**Q: How do I read just one output file interactively?**
A: Use the reader functions directly:
```julia
using Nbody6Setup
snap = read_conf3("conf.3_0")
```
