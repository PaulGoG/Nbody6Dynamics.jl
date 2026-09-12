# Nbody6Dynamics.jl — User Manual

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

**Nbody6Dynamics.jl** automates the full lifecycle of N-body star cluster simulations with [Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing):

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
- **Julia 1.10+** through [juliaup](https://github.com/JuliaLang/juliaup); the tracked Manifests were resolved with Julia 1.13
- **GCC toolchain**: `gcc`, `g++`, `gfortran`, `make`
- **Git** (for cloning the simulation code)

Optional, depending on config: MPI (`mpicc`, `mpif90`, `mpirun`), CUDA toolkit (`nvcc`), HDF5 development headers (`hdf5-devel` / `libhdf5-dev`) for the *build-side* HDF5 patch.

---

## 3. Installation

```bash
cd Nbody6Dynamics
julia activate.jl                                   # resolves, instantiates and precompiles the package environment
julia --project=. -e 'using Pkg; Pkg.test()'
```

`activate.jl` activates and instantiates the package environment silently; `docs/activate.jl` and `bench/activate.jl` do the same for the documentation and benchmark environments, developing the package by a relative path. The scripts under `scripts/`, `docs/` and `bench/` include their environment's activation script, so they run without a `--project` flag. The first `using Nbody6Dynamics` after an install or a source change precompiles the package together with a small workload (configuration parsing, the diagnostics reader, a merger initial-condition generation), so those paths run compiled in every later session.

---

## 4. Configuration Reference

Every key below is parsed by `load_config` (`src/config.jl`). Missing keys fall back to the defaults shown. All constraints listed below are enforced fail-fast at load time: `load_config` raises an error naming the offending `section.key` and the actual value, so a pipeline cannot start from a configuration it cannot honor. `save_config` writes a frozen snapshot of the full configuration into each run directory.

### `[install]`

| Key           | Type   | Default | Description |
|---------------|--------|---------|-------------|
| `enabled`     | Bool   | `true`  | Enable the install/build phase |
| `source_url`  | String | `"https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing.git"` | Git repository to clone |
| `ref`         | String | `"618d7a4"` | Commit, tag, or branch checked out after cloning (the default is the validated upstream v2026.07 commit); `""` keeps the default branch. An existing source directory is left as is |
| `install_dir` | String | `"backend/Nbody6PPGPU-beijing"` | Source directory, relative to the package root; must be nonempty |
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
| `cuda_arch`       | [String] | `[]`    | CUDA architectures compiled into the GPU kernels, `sm_<major><minor>` (`"sm_90"` H100/H200, `"sm_120"` RTX 50 series); native code for each plus PTX for the highest. Empty = the compute capabilities `nvidia-smi` reports, or the `nvcc` default target when no device is visible (see [GPU builds](#gpu-builds-and-target-architectures)) |
| `nvcc_flags`      | [String] | `[]`    | Extra `nvcc` options appended to the GPU build flags, e.g. `["-allow-unsupported-compiler"]` or `["-ccbin", "gcc-14"]` when the host compiler is newer than the toolkit supports; entries must be nonempty. The build probes `nvcc` on a trivial kernel first: when the toolkit rejects the host compiler it tries `-allow-unsupported-compiler`, then `-ccbin` with `CUDAHOSTCXX` and the versioned `g++-15` … `g++-12`, `clang++` found on `PATH`, and records what it settled on; a configured `-ccbin` is used as is. When the host's glibc (≥ 2.42) declares `rsqrt`/`rsqrtf` with an exception specification the CUDA headers lack, the probe retries the host-compiler choice with `-U_GNU_SOURCE -D_DEFAULT_SOURCE` and keeps it. When no combination works the error quotes the `nvcc` output of every attempt |
| `nproc`           | Int      | `0`     | Parallel `make` jobs; 0 = auto-detect; must be ≥ 0 |

### `[simulation]`

| Key             | Type   | Default | Description |
|-----------------|--------|---------|-------------|
| `run_test`      | Bool   | `true`  | Run the simulation phase |
| `input_file`    | String | `"examples/input_files/N10k_noDat10.inp"` | Path to the `.inp` input file; must be nonempty |
| `runs_dir`      | String | `"runs"` | Base directory for run output; must be nonempty |
| `binary_name`   | String | `"nbody6++"` | Expected binary name |
| `mpi_ranks`     | Int    | `1`     | Number of MPI ranks; must be ≥ 1, and > 1 requires `build.enable_mpi = true` |
| `omp_threads`   | Int    | `0`     | OpenMP threads for the backend, exported as `OMP_NUM_THREADS`; must be ≥ 0. `0` leaves the OpenMP runtime default: an inherited `OMP_NUM_THREADS`, else every logical CPU. Oversubscription (`omp_threads × mpi_ranks` above the host's logical CPUs) warns at launch |
| `gpu_list`      | [Int]  | `[]`    | CUDA device indices the engine may use, exported as its `GPU_LIST` variable; empty = every visible device. Entries ≥ 0 and distinct, at most 4 per process (the engine's `MAX_GPU`); nonempty requires `build.enable_gpu = true` |
| `run_id_prefix` | String | `"run"` | Prefix for run directory names; must be nonempty |
| `monitor`       | Bool   | `false` | Live ADJUST ticker on stderr; interactive terminals only |
| `live_diagnostics` | Bool | `false` | With `monitor`: print in-terminal sparklines of the virial ratio and `log10 |ΔE/E|` against time every `live_interval` seconds (UnicodePlots; interactive terminals only, never in the log file) |
| `live_interval` | Float | `30.0` | Period of the sparkline panel [s]; must be ≥ 1 |
| `telemetry_interval` | Float | `5.0` | Sampling interval of the process-tree/GPU telemetry [s]; must be ≥ 0. `0` disables the sampler; the exact CPU accounting stays on |
| `startup_timeout` | Float | `0.0` | Start-up watchdog [s]; must be ≥ 0. When set, a run that reports no adjustment beyond t = 0 within this wall-clock time is terminated with an error, and the kill is recorded in `RUN_INFO.toml` (`segments[].watchdog`). `0` disables; large-N runs need a generous value. The known cause of such a hang is fixed in the generator; see [Troubleshooting](#run-never-advances-past-t-0) for hand-written inputs |
| `exit_grace` | Float | `120.0` | Completion monitor [s]; must be ≥ 0. Once the stdout shows the engine's `END RUN` line, a process whose output directory has not changed for this long is terminated and the segment is recorded as `completed` with `terminated_after_completion = true`; a final COMMON dump still being written keeps it alive. `0` disables |

### `[postprocess]`

| Key                   | Type   | Default | Description |
|-----------------------|--------|---------|-------------|
| `enabled`             | Bool   | `true`  | Enable post-processing |
| `data_dir`            | String | `""`    | External data directory. Empty → use `run_dir/output/` from the simulation, or the most recent run under `runs/`. Set → read from this directory instead (any dir with Nbody6++ output) |
| `snapshot_format`     | String | `"conf3"` | Snapshot source. **Only `"conf3"` is supported** (see note below) |
| `snapshot_pattern`    | String | `"conf.3_*"` | Glob pattern for conf.3 files; must be nonempty |
| `parse_stdout`        | Bool   | `true`  | Parse simulation stdout for ADJUST diagnostics |
| `stdout_file`         | String | `"out1000"` | Name of the stdout capture file; must be nonempty when `parse_stdout = true` |
| `read_lagr`           | Bool   | `true`  | Read Lagrangian radii |
| `lagr_file`           | String | `"lagr.7"` | Lagrangian radii file name; must be nonempty when `read_lagr = true` |
| `read_escapers`       | Bool   | `true`  | Read escaper data |
| `escapers_file`       | String | `"esc.11"` | Escaper file name; must be nonempty when `read_escapers = true` |
| `read_stellar_evo`    | Bool   | `true`  | Read stellar evolution snapshots |
| `stellar_evo_pattern` | String | `"sev.83_*"` | Glob pattern for stellar evolution files (matches what this fork's `hrplot.F` writes: `sev.83_<time>`); must be nonempty when `read_stellar_evo = true` |
| `read_binary_evo`     | Bool   | `true`  | Read the regularised-binary snapshots |
| `binary_evo_pattern`  | String | `"bev.82_*"` | Glob pattern for the binary files (`bev.82_<time>`, written by `hrplot.F` alongside `sev.83`); must be nonempty when `read_binary_evo = true` |

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
| `format`     | String  | `"pdf"` | Static plot format; one of `"pdf"`, `"svg"`, `"png"` (animations are always GIF) |
| `dpi`        | Int     | `300`   | Resolution for raster formats; must be ≥ 72 |
| `column`     | String  | `"single"` | Journal-width preset; one of `"single"`, `"double"`, `""` (empty = free-form `figsize`) |
| `figsize`    | [Float] | `[8.0, 6.0]` | Figure size in inches `[width, height]`; both entries must be > 0 |
| `units`      | String  | `"physical"` | Axis units; one of `"physical"`, `"nbody"` |
| `output_dir` | String  | `"plots"` | Plot directory, relative to each run directory |

### `[visualization.style]`

Presentation knobs collected in the `PlotStyle` struct (`cfg.visualization.style`). These are stylistic only — data-validity cutoffs (e.g. BSE placeholder filtering in HR diagrams) remain named constants in the plot code.

| Key                   | Type  | Default   | Description |
|-----------------------|-------|-----------|-------------|
| `marker_budget`       | Float | `18000.0` | Scatter marker size is `clamp(marker_budget/N, marker_min, marker_max)`; must be > 0 |
| `marker_min`          | Float | `4.0`     | Lower clamp bound for the scatter marker size [pt]; must satisfy `0 < marker_min ≤ marker_max` |
| `marker_max`          | Float | `20.0`    | Upper clamp bound for the scatter marker size [pt]; must be ≥ `marker_min` |
| `q_log_threshold`     | Float | `10.0`    | Switch virial-ratio axes to log scale when `max(Q)` exceeds this; must be > 0 |
| `q_floor`             | Float | `1e-3`    | Clamp floor for the virial ratio on *log-scale* axes only; must satisfy `0 < q_floor < 1` |
| `zoom_frac`           | Float | `0.15`    | Extent-ratio threshold for adaptive per-panel zoom in snapshot evolution plots and cluster animations; must satisfy `0 < zoom_frac ≤ 1` |
| `anim_fps`            | Int   | `0`       | Animation frame rate; `0` selects automatically from frame count; must be ≥ 0 |
| `anim_target_seconds` | Float | `12.0`    | Target GIF duration used by the automatic FPS selection; must be > 0 |

### `[merger]`

| Key           | Type   | Default | Description |
|---------------|--------|---------|-------------|
| `enabled`     | Bool   | `false` | Generate merger ICs before the simulation phase |
| `config_file` | String | `""`    | Path to the merger cluster TOML (e.g. `"input_files/merger_demo_small.toml"`); relative paths resolve against the package root; must be nonempty when `enabled = true` |

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

The upstream `./configure --enable-hdf5` flag is broken. Nbody6Dynamics patches the build with a two-file mechanism:

1. Writes `hdf5_flags.mk` in the source root with the correct `-DCONFIG_HDF5`, include paths, and library flags
2. Appends `-include ../hdf5_flags.mk` to `build/Makefile`

This survives `./configure` re-runs: the flags file persists and only the one-line include needs reapplication. (This affects only what the *binary* can write — the Julia post-processing side reads `conf.3` snapshots, not HDF5; see the note in [Configuration Reference](#4-configuration-reference).)

### CUDA auto-detection

When `enable_gpu = true` and `cuda_path` is empty, the build searches in order:

1. Environment variables: `CUDA_HOME`, `CUDA_PATH`, `CUDA_ROOT`
2. Standard paths: `/usr/local/cuda`, `/usr/local/cuda-12`, `/usr/local/cuda-11.8`, …
3. `nvcc` location on `PATH`

### GPU builds and target architectures

The engine's `./configure` finds `nvcc` but emits no architecture flag, so on its own `nvcc` compiles the kernels for its default target and the driver JIT-compiles the embedded PTX at the first launch on any newer device. The build phase therefore resolves the targets itself and passes them to `make` as a `CUFLAGS` override (the configure-generated value plus one `-gencode` entry per architecture and PTX for the highest):

| Device | Compute capability | `cuda_arch` entry | Toolkit |
|--------|--------------------|-------------------|---------|
| RTX 5070 Ti, RTX 5090 (consumer Blackwell) | 12.0 | `"sm_120"` | CUDA ≥ 12.8 |
| H100, H200 (Hopper) | 9.0 | `"sm_90"` | CUDA ≥ 11.8 |
| A100 (Ampere) | 8.0 | `"sm_80"` | CUDA ≥ 11.0 |

With `cuda_arch` empty the build queries `nvidia-smi --query-gpu=compute_cap` and compiles for every distinct capability it reports; a build host without a visible device keeps the `nvcc` default (a warning says so), which still runs by JIT. A binary meant for several machines lists all of their architectures. The device capability can only come from the driver (`nvidia-smi`); `nvcc` knows the toolkit, not the hardware, so the install phase also reads `nvcc --list-gpu-arch` and stops before `make`, naming the release and the missing architecture, when the toolkit cannot compile for a target.

Other facts of the GPU build:

- The binary is named after the configure options: `nbody6++.avx.gpu` (and `nbody6++.avx.mpi.gpu` with MPI). CPU and GPU binaries may coexist in `build/`; the launcher picks the variant whose suffix tags match `enable_gpu`/`enable_mpi`. Switching one source tree between CPU and GPU builds requires `clean_build = true`, because the Fortran objects are compiled with `-D GPU` in one case only and `make` does not track the flag change.
- `BUILD_INFO.toml` is written next to the binary (date, host, backend commit, configure arguments, switches, CUDA path, compiled architectures, `nvcc` release). Every run copies it into its output directory and merges it into `RUN_INFO.toml` as the `[build]` table.
- One engine process drives at most four devices (`MAX_GPU` in the kernel source): with `gpu_list` empty it takes every visible device, otherwise those listed, and splits the regular-force j-range evenly across them by one OpenMP thread per device. The kernels compute in single precision on every architecture.
- The GPU library reports the devices it initialised on stderr; the run summary records them as `run.gpu_devices`, and the kernel throughput profile carries the label `GPU Reg.F` in `telemetry.force_kernel_gflops.kernel`, so a run that silently fell back to the CPU path is visible.
- The upstream authors advise the GPU build only above roughly 5×10⁴ bodies: below that the regular force is a minor share of the work and host–device transfers can make the run slower. `bench/gpu_scaling.jl` measures the crossover on the machine at hand.
- HDF5 output is unnecessary for the Julia side (it reads `conf.3`), so a GPU build can use `enable_hdf5 = false`.

### Validated hardware

The GPU path has been validated end to end on four hosts spanning four NVIDIA generations and two toolchain eras:

| GPU | Compute capability | CUDA | Host CPU | OS, glibc, GCC |
|---|---|---|---|---|
| RTX 5090, 32 GiB | 12.0 | 13.1 | Ryzen 9 9950X | Fedora 44, 2.43, 16.2.1 |
| RTX 5070 Ti, 16 GiB | 12.0 | 13.1 | i9-13900KS | Fedora 44, 2.43, 16.2.1 |
| RTX 2080 Super Max-Q, 8 GiB | 7.5 | 13.1 | i7-10750H | Fedora 44, 2.43, 16.2.1 |
| Tesla T4, 15 GiB | 7.5 | 11.8 | EPYC 7551P | Ubuntu 20.04, 2.31, 9.4 |

No host needed a package change. The `nvcc` host-compiler probe selected the plain toolchain on the Ubuntu 20.04 host and `-ccbin g++-15 -U_GNU_SOURCE -D_DEFAULT_SOURCE` on the Fedora ones on its own, so the same configuration file builds on a 2020 userspace with CUDA 11.8 and on a 2026 one with CUDA 13.1.

Measured against the AVX build of the same engine on the same host, for two King clusters merging on a Kepler orbit, the CUDA build is ahead at every size measured:

| N | Host threads | RTX 5090 | RTX 5070 Ti | RTX 2080 Super Max-Q |
|---|---|---|---|---|
| 19 638 | 4 | 2.00 | 1.95 | 1.84 |
| 19 638 | 8 | 1.46 | 1.39 | 1.69 |
| 49 226 | 4 | 3.38 | 3.38 | 2.93 |
| 49 226 | 8 | 2.62 | 2.54 | 2.75 |
| 98 426 | 4 | 4.55 | 4.61 | 4.02 |
| 98 426 | 8 | 3.26 | 3.52 | 3.77 |

There is no size below which the CPU build is the better choice; any crossover lies under N ≈ 2×10⁴, where a run costs seconds either way. Beyond that, the fleet shows what actually governs the gain. Force-kernel throughput spans a factor 5.7 across these three cards at N ≈ 10⁵ — 4.2, 15.2 and 23.9 TFLOP s⁻¹ — while the end-to-end speed-up spans only 4.02 to 4.61, and the 2019 mobile card reaches 87 % of what the RTX 5090 delivers. The pipeline is bound by the irregular force, which stays on the host: it takes 47 % to 89 % of the engine's accounted time in the 50 Myr merger runs, where mean device utilisation is 6 % to 16 %; nowhere in the campaign did it exceed 22 %. Peak device memory is under 0.6 GiB at N = 10⁵, so an 8 GiB card is not the constraint either. Below N ≈ 10⁵, a faster card is not how this pipeline gets faster; more concurrent jobs per card and faster host cores are.

Halving the host threads costs almost nothing on the GPU path. A four-thread CUDA job reaches 79 %, 81 % and 98 % of the corresponding eight-thread wall time on the three hosts above, and the speed-up over the AVX build is always larger at four threads. Four threads per job is the better unit for parameter sweeps.

The device computes the regular force in single precision and returns it as `double`. Over a 50 Myr merger of 2 × 25 000 stars, the fleet's cumulative energy errors run from −3.4×10⁻³ to +2.5×10⁻³ and the largest per-step excursion from 4.7×10⁻⁴ to 1.6×10⁻³, all well within what `qe = 0.01` tolerates, with no systematic separation between the CUDA and AVX builds. The same case ends with 47 023 to 47 119 of 47 999 bodies bound across all hosts, a 0.20 % spread that reflects chaotic divergence under a different summation order rather than a numerical defect. Results are reproducible to the last digit for a fixed host, binary and thread count, and only statistically so across hosts, so an energy error should always be quoted with the machine that produced it. A study needing tighter energy conservation should use the CPU build or a shorter regular-force timestep.

### Recipe for a CUDA host

Shipped under `input_files/gpu/`: `gpu_pipeline.toml` clones and builds the engine with CUDA into `backend/Nbody6PPGPU-beijing-gpu` (architectures from `nvidia-smi`, HDF5 off) and runs `merger_50k.toml`, two King clusters of 25 000 stars, on device 0 with eight host threads; `cpu_pipeline.toml` builds the AVX engine into the default tree and runs the same merger, so the two are comparable binary against binary. The order on a fresh host:

```bash
git clone git@github.com:PaulGoG/Nbody6Dynamics.jl.git Nbody6Dynamics && cd Nbody6Dynamics
julia activate.jl
julia scripts/run_gpu_validation.jl        # under tmux or nohup: the CPU reference of the 2 × 25k merger is the long stage
```

`run_gpu_validation` runs four logged stages as separate Julia processes and collects everything to pull back under `runs/gpu_validation_<machine>_<timestamp>/`: `HOST_INFO.toml` (GPU, driver, compute capabilities, CUDA path, `nvcc`, `gcc`, `gfortran`, glibc, the host compilers found for `-ccbin` and the verdict of the host-compiler probe, Julia, package commit), `suite.log` (the test suite with `NBODY6_GPU_TESTS=1`: CUDA build in a temporary tree, N = 1000 on one device, then on two when present), `gpu.log` (`gpu_pipeline.toml`), `cpu.log` (`cpu_pipeline.toml`), `bench.log` with the benchmark CSV and the `RUN_INFO.toml` and `telemetry.csv` of every benchmark run, and `VALIDATION.toml` (status, exit code and duration per stage, the run directories created). A failed stage does not stop the later ones; the benchmark stage is skipped, with the reason recorded, when a build tree it needs is missing.

A stage is judged by what it produced, not only by how it exited. A pipeline stage that exits zero is recorded `incomplete`, with the reason, unless it also left a run directory carrying the `[pipeline] completed` marker: the run summary is written the moment the engine exits, so a process killed during post-processing or plotting leaves a directory that would otherwise pass inspection. A stage killed by a signal is recorded with that signal, its output kept as `<stage>.signal<N>.log`, and retried once.

`<machine>` is the identity of the machine rather than its hostname: the hostname followed by compact CPU and GPU tags, as `[hardware] machine` in every run summary and host record. Hostnames are frequently not unique across a cloned workstation deployment, which would otherwise make the returned datasets indistinguishable. `--dry-run` writes the host record and the planned commands only; `--stages=suite,gpu` selects stages; `--n=`, `--threads=`, `--gpus="0;0,1"` and `--tcrit=` set the benchmark grid (the GPU lists default to device 0, plus devices 0 and 1 when two are visible). The same stages by hand:

```bash
NBODY6_GPU_TESTS=1 julia --project=. -e 'using Pkg; Pkg.test()'
julia scripts/run_setup.jl input_files/gpu/gpu_pipeline.toml
julia scripts/run_setup.jl input_files/gpu/cpu_pipeline.toml
NBODY6_GPU_BACKEND=backend/Nbody6PPGPU-beijing-gpu julia bench/gpu_scaling.jl 20000,50000,100000 4,8 "0" 0.25   # "0;0,1" with two devices
```

The toolkit need not be on `PATH`: the build looks at `CUDA_HOME`, `/usr/local/cuda` and the other standard locations, or at `cuda_path`, and exports the toolkit's `bin` and `lib64` to `configure`, `make` and the launch script; `configure` also receives `--with-cuda=<path>`. The engine's `configure` finds `nvcc` only on `PATH` (its `--with-cuda` fallback reuses the cached result of the `PATH` check), which is why the exported `PATH` matters. CUDA 13 toolkits work: the engine's own copy of `helper_cuda.h` reads `cudaDeviceProp` fields that CUDA 13.0 removed, so the build places the package's `deps/cuda` (NVIDIA's current samples headers) ahead of it in the `nvcc` include path. A host compiler newer than the toolkit supports (CUDA 13.x accepts GCC up to 15; Fedora 44 ships GCC 16) needs a compatible one for the CUDA host code only: install the distribution's side-by-side compiler (`sudo dnf install gcc15 gcc15-c++` on Fedora, which provides `g++-15`; `gcc14-c++` likewise) and the build picks it up through `-ccbin` on its own, or point `nvcc_flags = ["-ccbin", "<path>"]` at any supported `g++`. Fortran and the rest of the engine still compile with the system toolchain. A glibc newer than the toolkit supports is the other obstacle: glibc 2.42 and later (Fedora 43/44, Ubuntu 26.04, Debian 13) declare `rsqrt` and `rsqrtf` with an exception specification that CUDA 13.1's `crt/math_functions.h` lacks, and `nvcc` then rejects every host compiler with "exception specification is incompatible". The probe recognises the message and adds `-U_GNU_SOURCE -D_DEFAULT_SOURCE`, which keeps glibc's GNU-extension declarations out of the CUDA sources; the engine's GPU sources use none of them. The alternative is to patch `<toolkit>/targets/x86_64-linux/include/crt/math_functions.h`, adding `noexcept(true)` to the two declarations, which needs root and is undone by a toolkit update.

Each run's `RUN_INFO.toml` carries the build record (`[build]`: architectures, `nvcc` release), the devices the engine initialised (`run.gpu_devices`), the GPU telemetry summary and the kernel label of the throughput profile; the benchmark prints the GPU speed-up over the CPU binary at equal N and threads and writes `bench/results/gpu_scaling_<timestamp>.csv`. To use one binary on several machines set `cuda_arch = ["sm_90", "sm_120"]` explicitly; `gpu_list = [0, 1]` puts two devices in one process.

### Fedora h5pfc workaround

On Fedora the HDF5 parallel Fortran wrapper `h5pfc` may carry an erroneous `/openmpi-x86_64` suffix in its `includedir`, preventing `hdf5.mod` from being found. Nbody6Dynamics detects this and prints the fix:

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
├── RUN_INFO.toml        # run summary: identity/timing/thread layout, commits, hardware fingerprint, telemetry summary, file inventory, pipeline completion
├── telemetry.csv        # hardware telemetry time series (when telemetry_interval > 0)
├── nbody6dynamics.log   # teed pipeline log
├── output/              # all simulation artefacts
│   ├── nbody6++         # binary copy (reproducibility)
│   ├── _launch.sh       # generated bash launch script
│   ├── out1000          # captured stdout
│   ├── err1000          # captured stderr
│   ├── conf.3_*         # particle snapshots
│   ├── lagr.7           # Lagrangian radii
│   ├── esc.11           # escaper events
│   ├── sev.83_*         # stellar evolution snapshots
│   └── bev.82_*         # regularised-binary snapshots
└── plots/               # post-processing plots & GIF animations
```

Merger runs additionally contain `dat.10`, `merger.inp`, `merger_summary.txt`, and `merger_ic.toml` in `output/`.

### Launch script

The generated `_launch.sh` sets `ulimit -s unlimited` (Fortran stack), `OMP_STACKSIZE=4096M`, `OMP_NUM_THREADS` when `simulation.omp_threads > 0`, `GPU_LIST` when `simulation.gpu_list` is nonempty, CUDA environment variables (if GPU enabled), and `stdbuf -oL` for line-buffered output where available. The backend takes its thread count from the OpenMP runtime alone (there is no input parameter for it) and echoes it at start-up; that echoed value is recorded as `run.omp_threads_reported` in `RUN_INFO.toml` next to the configured `run.omp_threads`, `run.mpi_ranks` and `run.gpu_list`; the devices the GPU library initialised appear as `run.gpu_devices`, and the binary's build record as the `[build]` table.

### Choosing the thread count

The backend's OpenMP parallelism saturates early for the particle numbers a workstation handles: on a 22-thread machine the two-cluster benchmark (`bench/thread_scaling.jl`) reaches its shortest wall time at four threads for N ≤ 2×10⁴, while more threads only add CPU time (efficiency 0.38 at 22 threads). Parameter sweeps and ensembles are therefore best run as several four-thread jobs in parallel; the benchmark script sweeps thread count and N on your hardware and reports the fitted cost, using the telemetry of each run.

### Choosing the GPU devices

`bench/gpu_scaling.jl` runs the same two-cluster case with the CPU and the GPU binary over a grid of N, thread counts and `GPU_LIST` values, and prints the GPU speed-up at equal N and threads together with the backend's regular-force share and the kernel throughput. The regular force grows as N² against N⟨N_nb⟩ for the neighbour force, so its share, and with it the gain from the GPU, rises with N; at N = 2×10⁴ it is two thirds of the backend CPU time on the reference workstation, which bounds the gain there at about 3×. Two devices in one process (`gpu_list = [0, 1]`) halve only the regular-force term, so they pay off later in N than the first device. Under MPI every rank on a host applies the same `gpu_list`; per-rank device slicing is not provided.

### Real-time monitoring

During execution, ADJUST summaries are echoed live:

```
[ Info:   t_NB=0.0500  t_Myr=0.4  N=9998  |ΔE/E|=1.23e-06  Q_vir=0.987
```

### Completion detection

The engine prints `END RUN` (with the final timing tables) when its termination criterion is met, and normally exits right after. In a tidal-field run it has been seen to print everything and then never exit, which blocked the pipeline until an external timeout. With `simulation.exit_grace > 0` a monitor watches the last 64 KiB of the stdout capture for that line and, once no file of the output directory has changed for the grace period with the process still alive, terminates it. The engine writes its final COMMON dump after `END RUN` when `KZ(1) > 0`; that write keeps the directory changing, so a dump of any size is never interrupted. The segment then carries `completed = true`, `terminated_after_completion = true` and the signal exit status, post-processing proceeds, and the sweep summary counts the point as completed (`completed` column) even though its exit status is nonzero. A run that never printed `END RUN` is never touched by this monitor.

### Live sparklines

With `simulation.live_diagnostics = true` (and `monitor = true` on an interactive terminal) the monitor prints, every `live_interval` seconds, two in-terminal sparklines built with UnicodePlots from the ADJUST records so far: the virial ratio `Q = T/|W|` and `log10 |ΔE/E|` against time (Myr when the scaling is known). The panel is written to stderr below the log lines, the spinner resumes underneath, and nothing of it reaches `nbody6dynamics.log`. Runs that print fewer than two adjustments show no panel.

### Restarting a run

The engine writes a COMMON dump every `ncomm × deltat` N-body time units (`output/comm.1_<t>`, `output/comm.2_<t>`, alternating). `restart_simulation(run_dir; tcrit_extra = 5.0)` continues the run from the latest dump (or a chosen one, `dump = "comm.2_20.0"`) for `tcrit_extra` more N-body time units: the dump is copied to `output/comm.1`, which `KSTART = 2` reads; the original input recorded in `RUN_INFO.toml` (`run.input_file`, copied into `output/` at launch) supplies the `&ININPUT` block of the restart input, with `TCRIT` set to the increment the engine adds to the saved time; the termination time in Myr can be raised with `tcrtp0`. The engine runs in the same output directory with stdout and stderr appended, so the diagnostics, Lagrangian-radii, and escaper files continue and the time-stamped snapshot and stellar-evolution files carry on; post-processing reads the concatenated run as one. `RUN_INFO.toml` keeps one entry per launch in `segments` (kind, dump, extra time, elapsed, exit status), `run.elapsed_seconds` accumulates, and each segment's telemetry goes to `telemetry_<k>.csv`.

```julia
restart_simulation("runs/merger_demo_20260907_181355_2b7f"; tcrit_extra = 5.0)
```

### Hardware telemetry

Every run records the exact CPU consumption of the backend process tree from `getrusage(RUSAGE_CHILDREN)` deltas: `cpu_user_s`, `cpu_system_s`, `threads_total` (effective OpenMP threads × MPI ranks), and `cpu_efficiency = (user + system) / (elapsed × threads_total)`, the fraction of the reserved CPU capacity the integration actually used. These land in the `[telemetry]` table of `RUN_INFO.toml`.

With `simulation.telemetry_interval > 0` (default 5 s) an asynchronous sampler additionally writes `runs/<run_id>/telemetry.csv`, one row per interval, with the columns of `TelemetrySample`: elapsed time, process count, resident memory and high-water mark summed over the tree (from `/proc/<pid>/status`), cumulative CPU time and the derived cores-busy rate (from `/proc/<pid>/stat`), the 1-minute load average, and — when `build.enable_gpu` — one `nvidia-smi --query-gpu` sample per interval (utilization, memory utilization, memory used, power, temperature; mean, sum, or maximum over devices). The summary adds `peak_rss_mib`, `mean_cores_busy`, `peak_cores_busy`, `peak_load_1min`, and the GPU means and peaks when GPU samples exist. Sampling is Linux-only; a failing `nvidia-smi` disables GPU sampling for the rest of the run after one warning, and a sampler failure never aborts the run.

The backend's own performance report is captured as well. The last timing table it prints to stdout (one per `DTADJ`, cumulative CPU seconds per code section: regular and irregular force, prediction, KS, adjust, output, communication, …) becomes the `[telemetry.backend_timing]` sub-table with lower-case keys (`total`, `reg`, `irr`, `ks`, `reg_gpu_s`, …), and the `Perf.(Gflops)` lines of the AVX/SSE or GPU regular-force profiles on stderr are summarised as `[telemetry.force_kernel_gflops]` (`samples`, `mean`, `peak`).

---

### Parameter sweeps

A sweep runs the merger pipeline over a Cartesian grid of merger-TOML values times a set of seeds, one directory per point:

```toml
[sweep]
name = "demo"                              # letters, digits, "_", "-"
pipeline_config = "../config.toml"         # base pipeline config, relative to this file
merger_config = "merger_demo_small.toml"   # base merger config, relative to this file
seeds = [11, 12]                           # merger.seed at every grid point; distinct integers
concurrency = 5                            # simultaneous jobs; ≥ 1
omp_threads = 4                            # OpenMP threads per job; ≥ 1
runs_dir = "../runs"                       # sweep root: runs_dir/sweep_<name>_<timestamp>/
poll_interval = 2.0                        # seconds between job checks; > 0
controls = false                           # isolated single-cluster control companion per point

[sweep.grid]                               # dotted merger-TOML keys → arrays; Cartesian product
"merger.orbit.eccentricity" = [0.0, 0.6]
"merger.cluster2.N" = [500, 1000]
```

```bash
julia scripts/run_sweep.jl input_files/sweep_demo.toml            # run
julia scripts/run_sweep.jl input_files/sweep_demo.toml --dry-run  # prepare only
```

Axes address the merger TOML from its root; the parent table must exist in the base file (add `[merger.tidal]` with `kz14 = 0` to sweep the tidal field), and `merger.seed` is reserved for `seeds`. Axes are processed in key order with the first varying fastest; point directories read `<index>_<axis>=<value>_…_seed=<seed>`.

`prepare_sweep` writes each point's `merger.toml` and `config.toml` (the base pipeline config with the install phase disabled, an absolute backend path, the point directory as run root, the derived merger file and the sweep's `omp_threads`) and loads both back through the regular parsers, so an invalid point fails before anything runs. `run_sweep` then executes the points as separate worker processes (`run_sweep_point`, at most `concurrency` at a time, each logging to `<point>/sweep_point.log`, the run itself in `<point>/run/`), rewrites `sweep_index.toml` on every state change (`pending`, `running`, `done`, `failed` with exit status and elapsed time), and writes `sweep_summary.csv`: index, id, seed, the axis values, status, elapsed time, exit status, and the final time, star count, pair count, energy error and virial ratio from the last ADJUST record. `sweep_figures` draws the half-mass Lagrangian radius and |ΔE/E| of every completed run on common physical axes, coloured by the value of one grid axis (the first by default; seeds share the colour), into `<sweep>/plots` with the base config's `[visualization]` settings. Four threads per job saturate the backend for N ≲ 2×10⁴ on the reference workstation, so five concurrent four-thread jobs use it fully.

A sweep without `[sweep.grid]` is a seed ensemble of the base merger configuration (points `001_seed=11`, …). Whenever a sweep has more than one seed, `sweep_ensembles` groups the completed points by grid values and `ensemble_statistics` interpolates every member's series onto a common time grid (the interval covered by all members) and takes the median and the central 68 % and 95 % intervals per node (type-7 quantiles). `plot_sweep_ensemble` draws the median line with the two bands for `:lagrangian`, `:energy`, `:n_stars` or `:n_pairs`, one colour per value of the chosen axis with the other axes held at `fixed` values (their first value by default, annotated); `sweep_figures` adds the Lagrangian-radius and energy-error ensembles automatically.

With `controls = true` every point gets a companion `<index>_…_control` whose merger TOML is the isolated equivalent of the point's configuration (`control_merger_dict`): one cluster at rest with the summed `N`, the `N`-weighted mean half-mass radius, and cluster 1's model, IMF and binary population, with the `output`, `nbody6`, `stellar` and `tidal` sections kept verbatim. The index and summary carry `kind` and `control_of`; the comparison and ensemble figures use the merger points only, and `plot_control_comparison` draws every merger (solid) against its control (dashed). `write_control_merger_config` derives a control file from any merger TOML outside a sweep; the generator accepts a single cluster in explicit mode for this purpose.

## 8. Post-processing

### Snapshots (conf.3, Fortran binary)

```julia
snap = read_conf3("path/to/conf.3_0")
snap.header        # SnapshotHeader with AS(1:20) parameters
snap.pos           # 3 × N matrix (Float32)
snap.vel           # 3 × N matrix (Float32)
snap.mass          # N-element vector (Float32)
snap.name          # particle identifiers (Int32)

snaps = read_all_conf3(run_out, "conf.3_*")              # ordered by the time suffix
snaps = read_all_conf3(run_out; threaded = true)         # chunked over Threads.@spawn
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

### Regularised binaries

```julia
bev = read_binary_evolution("bev.82_0")
bev.n_pairs                       # header NPAIRS
bev.records[1].eccentricity
bev.records[1].log_period_days    # log10(P / d)
semi_major_axis_pc(bev.records[1])

bevs = read_all_binary_evolution(run_out, "bev.82_*")
snaps = read_all_conf3(run_out)
scales = binary_scales(bevs, snaps)          # ⟨m⟩, σ and N per epoch from the nearest snapshot
pop = binary_population(bevs; n_stars = scales.n_stars,
                        m_mean = scales.m_mean, sigma_kms = scales.sigma_kms)
pop.n_hard, pop.n_soft, pop.binary_fraction
```

`bev.82` lists the KS-regularised pairs only (32 columns per line, the same header convention as `sev.83`): component indices, names and stellar types, the distance of the centre of mass from the density centre, eccentricity, `log10(P/d)`, `log10(a/R☉)`, and the SSE quantities of both components. Pairs wider than the regularisation distance are absent, so every count derived from it is a lower bound on the bound-pair population. Hardness follows [Heggie1975](@cite): a pair is hard when its binding energy `G m₁ m₂ / (2a)` exceeds `⟨m⟩ σ²`, with the mean system mass and the one-dimensional, mass-weighted dispersion of the systems (singles plus pair centres of mass) taken from the snapshot nearest in time (`hardness_scale`). The pipeline draws `binary_population` (pair counts and hard/soft split above the binary and hard fractions), `binary_period_distribution` (first against last epoch), and `binary_orbital_elements_initial`/`_final` (semi-major axis against eccentricity by hardness class, with the boundary of a pair of mean component-mass product). `scan_output` detects the files and `postprocess_external` reads and plots them like the config-driven pipeline, so run directories produced elsewhere get the same binary diagnostics.

---

## 9. Visualisation

### The figure backend is an extension

Every figure routine — `plot_*`, `animate_*`, `generate_plots`, `sweep_figures`, `remnant_figures`, `publication_theme` — lives in the package extension `Nbody6DynamicsMakieExt`, which Julia loads as soon as a Makie backend is in the session:

```julia
using CairoMakie
using Nbody6Dynamics
```

The package itself declares no plotting dependency. A headless host therefore installs 109 packages instead of 278 and never fetches or precompiles Cairo, Pango, GLib, HarfBuzz or the font artifacts: initial conditions, engine build, integration, readers, diagnostics, sweeps and benchmarks are all backend-free. The distinction is not cosmetic — on a compute server whose system GLib is older than the artifact one, the plotting stack is exactly what fails to load, and before the split it took the numerical work down with it.

Consequences to know:

- Without a backend, a figure routine throws `PlottingUnavailable` naming the remedy; it never fails with an `UndefVarError`.
- `run_pipeline` with `[visualization] enabled = true`, `run_merger_pipeline(...; make_plots = true)` and `postprocess_external(...; make_plots = true)` check for the backend **before** doing any work, so a missing backend costs nothing instead of costing an integration. Set the flag to `false` for a numerics-only run.
- `plotting_available()` reports whether the extension is loaded.
- The entry scripts under `scripts/` run in their own environment (`scripts/Project.toml`) which carries the backend, so `julia scripts/run_setup.jl config.toml` plots as before. `bench/` has no backend, by design.
- A sweep worker inherits the driver's environment and loads the backend only if the driver has one, so sweeps behave like the session that launched them.

All plots use the built-in publication theme, activated globally with `set_publication_theme!()` (called automatically by the pipeline); `publication_theme()` returns it as a value for `with_theme` scoping. The theme is built at call time so that its Computer Modern faces come from MathTeXEngine's live registry: a face captured while the package precompiles carries a null FreeType pointer, and Makie would then render every plain-text label in its default sans font without warning.

### Theme rules

- **No plot titles** — contextual information (time, projection label, orbit parameters) is placed as in-axis annotations at the top-left corner
- **No minor ticks**; major ticks face inward; fully boxed axes
- **Grid**: dashed grey at very low opacity — on *line plots only*. Dense scatter plots (cluster projections, HR diagrams) disable the grid locally
- **Fonts**: Computer Modern (NewComputerModern, bundled with MathTeXEngine) for text and math
- **Never-overwrite policy**: before saving, any existing file at the target path is moved to the first free `name#k.ext` sibling (DrWatson-`safesave` style); this applies to plots, animations, `dat.10`, `merger.inp`, and metadata files alike
- Legends appear only when ≥ 2 items are plotted; axis limits snap to nice round values with a visual buffer

Multi-panel montages (`snapshot_evolution_*`, `hr_evolution`) are one column wide like every other figure: the panels share the preset width (a shared colorbar column is taken from the panel area), use three ticks per axis, keep a data-free band above the data for the in-axis time annotation, and scale their markers with the panel width; inner tick labels are hidden and the panel gap is compact unless the adaptive zoom shows every panel's ticks. A three-column grid at the `single` preset has 25 mm panels, so `double` is the preset for montages placed in a manuscript.

### Plot inventory

Static plots take the extension from `visualization.format`; animations are always `.gif`.

| Function | Description | Output files |
|----------|-------------|--------------|
| `plot_snapshot(snap, vis)` | Scatter projections, viridis mass colouring `log₁₀(m/M_tot)` | `snapshot_final_{xy,xz,yz}` |
| `plot_snapshot_evolution(snaps, vis)` | Up to 6 panels, shared colorbar; adaptive per-panel zoom when extents vary by more than `1/zoom_frac` | `snapshot_evolution_{xy,xz,yz}` |
| `plot_energy(diag, vis)` | Two panels: log `\|ΔE/E\|` and virial ratio with `Q = 0.5` reference | `energy` |
| `plot_particle_count(diag, vis)` | Two panels: bound N and KS pairs, integer ticks | `particle_count` |
| `plot_lagrangian(lagr, vis)` | Selected mass-fraction radii vs time, log y | `lagrangian_radii` |
| `plot_cluster_separation(snaps, ranges, vis)` | Pairwise centre-of-mass separations of the initial clusters (merger runs) | `merger_cluster_separation` |
| `plot_cluster_virial(snaps, ranges, vis)` | Virial ratio of each initial cluster from its bound members, `Q = 0.5` reference | `merger_cluster_virial` |
| `plot_cluster_structure(snaps, ranges, vis)` | Two panels: bound half-mass radius per cluster with the engine's global `r₅₀` overlaid, and bound mass fraction | `merger_cluster_structure` |
| `plot_density_profiles(snap, ranges, vis; specs)` | Density profiles of each cluster about its own centre (log–log) with the generating King/Plummer model dashed and a `ρ/ρ_model` ratio strip with unity guide | `merger_density_profiles_{initial,final}` |
| `plot_velocity_dispersion(snap, ranges, vis)` | Radial (solid) and tangential (dashed) velocity dispersion profiles per cluster, and the anisotropy `β(r)` with the isotropic guide | `merger_velocity_dispersion` |
| `plot_hr(sev, vis)` | HR diagram coloured by stellar type K*, reversed Teff axis | `hr_diagram_{early,mid,final}` |
| `plot_hr_evolution(sevs, vis)` | HR panel grid across up to 6 epochs | `hr_evolution` |
| `plot_escapers(escs, vis)` | Two panels: cumulative escaped mass step curve with totals annotation, and escape velocity vs time (log y) split into luminous / compact-remnant classes | `escapers` |
| `plot_escape_anisotropy(escs, vis)` | Sky projection of escape directions φ ∈ [0°, 360°], θ ∈ [-90°, 90°] by stellar class | `escape_anisotropy` |
| `plot_mass_segregation(sev, vis)` | Distance from density centre RI vs stellar mass, log–log, epoch annotated | `mass_segregation` |
| `plot_evolutionary_clock(sev, vis)` | Histogram of the MS age fraction t/T_MS (K* ≤ 1), turnoff boundary and past-turnoff fraction annotated | `evolutionary_clock` |
| `plot_core_mass(sevs, vis)` | Core mass MC vs total mass for evolved stars (K* ≥ 2) at the last epoch, MC = M identity guide | `core_mass_growth` |
| `plot_cluster_separation(snaps, ranges, vis)` | Pairwise COM separations of the initial clusters; per-pair lines for ≤ 5 clusters, min/max/mean envelope + surviving-cluster staircase and coalescence marker otherwise | `merger_cluster_separation` |
| `plot_cluster_virial(snaps, ranges, vis)` | Internal virial ratio `Q_i(t)` per initial cluster (COM-subtracted, self-gravity only); log axis when `max(Q) > q_log_threshold` | `merger_cluster_virial` |
| `plot_merger_ic(result, vis)` | IC diagnostics: projections, 3-panel overview, velocity quiver by cluster, IMF histogram with Kroupa reference slopes, per-cluster radial density | `merger_ic_{xy,xz,yz}`, `merger_ic_overview`, `merger_ic_velocity`, `merger_ic_imf`, `merger_ic_density` |
| `animate_cluster(snaps, vis)` | Animated scatter per projection, global or adaptive limits | `cluster_evolution_{xy,xz,yz}.gif` |
| `animate_lagrangian(lagr, vis)` | Progressive line draw with ghost background and time cursor | `lagrangian_anim.gif` |
| `animate_hr(sevs, vis)` | Animated HR diagram | `hr_evolution_anim.gif` |
| `plot_telemetry(samples, vis)` | Run telemetry against wall-clock time: cores busy with the host load average on a twin axis, resident memory (RSS, high-water mark), GPU utilisation when sampled; the means and the peak RSS are legend entries. Drawn for every run directory that holds the sampler's `telemetry*.csv` (`read_run_telemetry` concatenates the segments of restarted runs) | `telemetry` |

The merger plots (`plot_cluster_separation`, `plot_cluster_virial`) are generated automatically by `generate_plots`/`postprocess_external` when a `merger_summary.txt` is found next to the snapshots; the cluster index ranges come from `parse_merger_summary`. `plot_merger_ic` runs whenever a merger IC was generated in the same pipeline invocation (or can be re-run later via `load_merger_ic_result`).

### Customising plots

All plot functions accept a `VisualizationConfig` and a `filename` keyword:

```julia
vis = VisualizationConfig(format = "pdf", dpi = 600, column = "double",
                          output_dir = "my_plots", style = Nbody6Dynamics.PlotStyle())
plot_snapshot(snap, vis; filename = "cluster_final")
```

---

## 10. Using the Julia API Directly

```julia
using Nbody6Dynamics

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
| `BinaryRecord` / `BinaryEvolutionSnapshot` | Regularised binaries from one bev.82 file |
| `BinaryPopulation` | Pair counts, hard/soft split and binary fraction against time |
| `RemnantDiagnostics` / `RotationProfile` | Bound-remnant core radius, rotation, mass segregation and coalescence time against snapshots |
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
| `bev.82_*` | Regularised-binary snapshots (ASCII) | `read_binary_evolution` |
| `snap.40_*.h5part` | HDF5 snapshots (KZ(46)) | **not readable** — detected and warned about only |

### Stellar type codes (K*)

[Hurley2000](@cite) convention as used by this fork (`STELLAR_TYPE_LABELS`):

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

Nbody6++ uses Hénon N-body units [Henon1971](@cite) internally: `G = 1`, `M_total = 1`, `E_total = -1/4`. The virial ratio is `Q = T/|W|` with equilibrium at `Q = 0.5`.

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

### Run never advances past t = 0

The engine prints the t = 0 adjustment, then sits at 100 % of one core forever; the start-up watchdog ends it. The cause is the engine's `string_left.f`, which counts the decimal digits of `DELTAT`, `DTADJ` and `DTPLOT` by multiplying by ten until the value is an integer, with a default-kind `int`. A value such as `0.6302` never becomes an integer in binary arithmetic; past 2³¹ the conversion overflows and the loop never ends, inside the very first output. The generator therefore writes these three intervals as dyadic rationals with an exact decimal expansion (`engine_interval`, change below 0.4 %) and logs the rounding. For a hand-written input, choose intervals such as `0.5`, `0.25`, `0.125`, `0.0625` or any other `m/2ᵏ` with `k ≤ 9`, written out in full (exactly ten decimals trip a second defect of the same routine, an invalid `I1` format, which ends the run with a runtime error at the next output); `engine_interval(x)` gives the nearest such value. This hang was for a time attributed to a small initial neighbour radius; the neighbour radius has no part in it.

### Simulation segfaults immediately

Fortran code requires a large stack. The launch script sets `ulimit -s unlimited`; if running manually, set this in your shell first.

### Simulation exits with status 2 after a burst of escapers

`err1000` ends with `Fortran runtime error: Expected REAL for item ... in formatted transfer` at `escape.F`. The engine's escaper summary `WRITE` has a fixed item list, and a single adjustment interval that flags on the order of a thousand escapers overflows it. This happens in isolated multi-cluster runs, where the escape radius is `2 × 10 RSCALE` about the global density centre and ejecta from an early collapse cross it together many crossing times later. Shorten `DTADJ`, run in an external tidal field, or treat escapers in post-processing from the snapshots; the output written before the failure is intact.

### "CUDA not found" during build

Set `cuda_path` explicitly in `[build]`, or export one of `CUDA_HOME`, `CUDA_PATH`, `CUDA_ROOT`.

### "nvcc fatal: Unsupported gpu architecture 'compute_120'"

The toolkit predates the device: consumer Blackwell (`sm_120`) needs CUDA 12.8 or later, Hopper (`sm_90`) CUDA 11.8 or later. Install a current toolkit and point `cuda_path` at it, or drop the architecture from `cuda_arch`.

### "no kernel image is available for execution on the device"

The binary carries native code for other architectures and no PTX the driver can compile for this one. Rebuild with `cuda_arch` empty (the visible device's capability is detected) or including the device's `sm_<major><minor>`; a binary shared between machines lists all of their architectures.

### GPU build runs, but `RUN_INFO.toml` shows no `run.gpu_devices`

The engine used the CPU path: the launcher picked a binary without the `.gpu` suffix (check `build.binary` in the run summary and `enable_gpu` in the config), or the GPU library found no device (`GPU_LIST` names an index that does not exist, or the driver is not loaded — compare `hardware.gpu`).

### Post-processing finds no files

Check `simulation.runs_dir` and, for external data, `postprocess.data_dir`. For a specific run:

```julia
results = postprocess(cfg; run_dir = "runs/run_20260325_143022_a1f3")
generate_plots(results, cfg; run_dir = "runs/run_20260325_143022_a1f3")
```

### HDF5 snapshots present but nothing is read

Expected: the `KZ(46)` H5Part layout is not supported (see the note in [Configuration Reference](#4-configuration-reference)). Enable `conf.3` output in the `.inp` file (`KZ(3)` ≥ 1) and use `snapshot_format = "conf3"`.
