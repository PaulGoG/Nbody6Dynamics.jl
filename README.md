# Nbody6Setup.jl

A Julia project for automated setup, execution, post-processing, and visualisation of the [Nbody6PPGPU-beijing](https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing) N-body astrophysical simulation code.

## Project Structure

```
Nbody6Setup/
├── Project.toml                     # Julia package metadata & dependencies
├── LICENSE                          # MIT license
├── config.toml                      # User configuration (edit this)
├── .gitignore
├── README.md
├── src/
│   ├── Nbody6Setup.jl               # Main module & orchestration
│   ├── types.jl                     # Core type definitions
│   ├── config.jl                    # TOML configuration loader & serialiser
│   ├── platform.jl                  # Platform detection, CUDA, HDF5, dependencies
│   ├── install.jl                   # Git clone, configure, HDF5 patch, build
│   ├── run.jl                       # Simulation execution with run ID & monitoring
│   ├── external.jl                  # scan_output, postprocess_external entry points
│   ├── io/
│   │   ├── io.jl                    # I/O submodule includes
│   │   ├── fortran_binary.jl        # Fortran unformatted binary record reader
│   │   ├── conf3.jl                 # conf.3 snapshot reader (standard + extended)
│   │   ├── hdf5_reader.jl           # HDF5/H5Part snapshot reader
│   │   ├── diagnostics.jl           # stdout (out1000) parser — ADJUST lines
│   │   ├── lagr.jl                  # lagr.7 Lagrangian radii reader
│   │   ├── esc.jl                   # esc.11 escaper reader
│   │   └── stellar_evolution.jl     # sev*.83 stellar evolution reader
│   ├── ic/
│   │   ├── ic.jl                    # IC submodule includes
│   │   ├── config.jl                # MergerConfig types & TOML parser
│   │   ├── models.jl                # Plummer & King model samplers
│   │   ├── imf.jl                   # Kroupa (2001) broken power-law IMF
│   │   ├── orbits.jl                # Kepler orbits, Jacobi truncation
│   │   └── output.jl               # dat.10 writer & .inp generator
│   └── plotting/
│       ├── plotting.jl              # Publication theme (CM fonts) & helpers
│       ├── snapshots.jl             # 2D projection scatter plots
│       ├── lagrangian.jl            # Lagrangian radii evolution
│       ├── energy.jl                # Energy error & virial ratio plots
│       ├── hr.jl                    # Hertzsprung-Russell diagram
│       └── animation.jl             # GIF animations (cluster, HR, Lagrangian)
├── scripts/
│   └── run_setup.jl                 # Main entry-point script
├── test/
│   ├── runtests.jl                  # Full test suite (311 tests)
│   └── test_external_adversarial_inner.jl  # Adversarial edge-case tests
├── benchmark/
│   └── benchmarks.jl                # I/O performance benchmarks
├── docs/
│   ├── make.jl                      # Documenter.jl build script
│   ├── Project.toml                 # Docs environment dependencies
│   └── src/                         # Documenter source pages
│       ├── index.md                 # Landing page
│       ├── api.md                   # API reference (autodocs)
│       ├── manual.md                # User manual
│       └── input_files.md           # .inp file format reference & design guide
├── input_files/                     # Custom simulation input files
│   ├── N25k_production.inp          #   25k particles, standard King model
│   ├── N100k_production.inp         #   100k particles, production run
│   ├── imbh_runaway.inp             #   IMBH formation via runaway collisions
│   ├── gc_bh_subsystem.inp          #   Globular cluster with BH subsystem
│   ├── tidal_tails.inp              #   Tidal stripping in Milky Way potential
│   ├── young_massive_binaries.inp   #   Binary-rich young massive cluster
│   ├── pop3_cluster.inp             #   Population III near-zero metallicity
│   ├── merger_equal_mass.toml       #   Equal-mass King merger, kepler mode (q=1)
│   ├── merger_minor_plummer.toml    #   Minor Plummer merger, kepler mode (q=0.1)
│   └── merger_triple_cluster.toml   #   3-cluster triangular merger, explicit mode
├── backend/                         # ← cloned N-body source code (black-box)
│   └── Nbody6PPGPU-beijing/         #   upstream Fortran/C++ simulation code
└── runs/                            # ← simulation outputs, grouped by run ID
    └── run_YYYYMMDD_HHMMSS_XXXX/    #   one directory per run
        ├── config.toml              #   frozen config snapshot
        ├── RUN_INFO.txt             #   run metadata summary
        ├── output/                  #   simulation output files
        │   ├── nbody6++             #   binary copy (reproducibility)
        │   ├── _launch.sh           #   generated launch script
        │   ├── out1000              #   captured stdout
        │   ├── err1000              #   captured stderr
        │   ├── conf.3_*             #   particle snapshots
        │   ├── lagr.7               #   Lagrangian radii
        │   ├── esc.11               #   escaper events
        │   └── sev*.83              #   stellar evolution snapshots
        └── plots/                   #   post-processing plots & GIFs (18 files)
```

## Quick Start

```bash
cd Nbody6Setup

# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Edit configuration
$EDITOR config.toml

# Run the full pipeline (install → build → simulate → postprocess → plot)
julia --project=. -e 'using Nbody6Setup; run_pipeline(load_config("config.toml"))'

# Or run with a custom config
julia --project=. -e 'using Nbody6Setup; run_pipeline(load_config("/path/to/my_config.toml"))'
```

### Usage Modes

The pipeline supports four modes, controlled by `config.toml`:

| Mode              | Settings                                                      |
|-------------------|---------------------------------------------------------------|
| Full pipeline     | `install.enabled=true`, `simulation.run_test=true`            |
| Simulate only     | `install.enabled=false`, `simulation.run_test=true`           |
| Postprocess only  | `simulation.run_test=false`, `postprocess.data_dir="..."`     |
| Re-plot latest    | `simulation.run_test=false`, `postprocess.data_dir=""`        |

### Programmatic API

```julia
using Nbody6Setup

# Full pipeline from config
cfg = load_config("config.toml")
results = run_pipeline(cfg)

# Post-process an external output directory directly
results = postprocess_external("/path/to/output")

# Quick scan of available output files
scan = scan_output("/path/to/output")
println(scan)

# Generate multi-cluster merger initial conditions
cfg = load_merger_config("input_files/merger_equal_mass.toml")
generate_merger_ic(cfg)  # → dat.10 + merger.inp
```

## Merger Initial Conditions

Generate initial conditions for multi-cluster merger simulations.

### Two orbit modes

| Mode | Clusters | Positions | Use case |
|------|----------|-----------|----------|
| `"kepler"` | Exactly 2 | Auto-computed from apocentre + eccentricity | Binary cluster mergers |
| `"explicit"` | Any N ≥ 2 | User-specified per-cluster `position`/`velocity` | Hierarchical or multi-body encounters |

### Standalone usage

```julia
using Nbody6Setup

# One-liner (loads TOML, generates ICs + plots)
result = run_merger_pipeline("input_files/merger_equal_mass.toml")

# 3-cluster explicit mode
result = run_merger_pipeline("input_files/merger_triple_cluster.toml")
```

### Integrated via config.toml

```toml
# In config.toml:
[merger]
enabled     = true
config_file = "input_files/merger_equal_mass.toml"
```
```julia
cfg = load_config("config.toml")
results = run_pipeline(cfg)  # Phase 1.5 generates ICs, Phase 4 plots them
```

### Programmatic API

```julia
cfg = MergerConfig(
    [
        ClusterSpec(model="king", N=50000, W0=6.0, mass_total=1e5, rbar=2.0),
        ClusterSpec(model="king", N=50000, W0=4.0, mass_total=5e4, rbar=3.0),
    ],
    "kepler",
    OrbitSpec(apocentre=15.0, eccentricity=0.7),
    MergerOutputSpec(format="nbody", truncate_jacobi=true)
)
result = generate_merger_ic(cfg)
```

Supported density profiles: **King** (recommended, finite tidal radius) and **Plummer**. Masses drawn from a Kroupa (2001) IMF. Optional nearest-neighbour Jacobi truncation. Output: `dat.10` + `.inp` for Nbody6++ with `KZ(22)=2`.

See [docs/src/multi_cluster_mergers.md](docs/src/multi_cluster_mergers.md) for the full physics background and parameter guide.

## Configuration

All behaviour is controlled by `config.toml`.  Key sections:

| Section           | Purpose                                                    |
|-------------------|------------------------------------------------------------|
| `[install]`       | Source URL, backend install directory, reinstall toggle     |
| `[build]`         | Configure flags, MPI/GPU/HDF5 toggles, CUDA path, nproc   |
| `[simulation]`    | Input file, runs directory, MPI ranks, run ID prefix       |
| `[postprocess]`   | Output readers: conf.3, HDF5, stdout, lagr, escapers, sev  |
| `[visualization]` | Plot format (png/pdf/svg), DPI, figure size                |

Set `enabled = false` or `run_test = false` in any section to skip that phase.

## Phases

### 1. Install & Build

- Clones the repository into `backend/` (or skips if already present)
- Checks all build dependencies (compilers, MPI, CUDA)
- Runs `./configure` with user-specified flags
- Patches `build/Makefile` for HDF5 via a robust `-include` mechanism that survives `./configure` reruns
- Auto-detects CUDA installation path from environment variables and standard locations
- Detects and warns about the Fedora `h5pfc` broken includedir bug
- Compiles with `make -j$(nproc)`

### 2. Simulation

- Creates an isolated run directory: `runs/<run_id>/output/`
- Saves a frozen copy of `config.toml` and the binary for reproducibility
- Generates a bash launch script with `ulimit -s unlimited`, `OMP_STACKSIZE=4096M`, and CUDA environment
- Uses `stdbuf -oL` for line-buffered output when available
- Monitors stdout in real-time, printing ADJUST summaries to the terminal
- Writes `RUN_INFO.txt` with run metadata after completion

### 3. Post-processing

Reads simulation output files:

| Reader                      | File             | Format                                        |
|-----------------------------|------------------|-----------------------------------------------|
| `read_conf3`                | `conf.3_*`       | Fortran unformatted binary (standard/extended) |
| `read_hdf5_snapshot`        | `data.40.h5part`  | HDF5/H5Part                                  |
| `read_diagnostics`          | `out1000`        | ASCII (ADJUST lines, PHYSICAL SCALING)         |
| `read_lagr`                 | `lagr.7`         | ASCII (block-structured Lagrangian radii)      |
| `read_escapers`             | `esc.11`         | ASCII (escaper events)                         |
| `read_stellar_evolution`    | `sev*.83`        | ASCII (stellar properties per epoch)           |

### 4. Visualisation

Publication-quality plots generated with CairoMakie, using **Computer Modern (LaTeX) fonts** via `LaTeXStrings.jl` and the `NewComputerModern` OTF family bundled with MathTeXEngine:

**Static plots (PNG/PDF/SVG):**

- **Snapshot projections** (×3) — XY/XZ/YZ scatter plots coloured by mass with colorbar
- **Snapshot evolution** (×3) — Multi-panel time sequence with mass colormap and shared colorbar (XY/XZ/YZ)
- **Energy diagnostics** — Relative energy error |ΔE/E| and virial ratio vs time
- **Particle count** — N(t) and N_pairs(t) evolution with integer ticks
- **Lagrangian radii** — Selected mass-fraction radii vs time (log scale)
- **HR diagrams** (×3) — Hertzsprung-Russell diagram at early, mid, and final epochs, coloured by stellar type
- **HR evolution** — Multi-panel HR diagram across 6 epochs

**Animations (GIF):**

- **Cluster evolution** (×3) — Animated scatter plot in XY/XZ/YZ with consistent axis limits and colour scale
- **Lagrangian radii** — Progressive line drawing with time cursor over ghost background
- **HR evolution** — Animated HR diagram with per-frame stellar type colouring

Total: **18 output files** per run (13 static + 5 GIF animations).

All plots use a fully-boxed publication theme with Computer Modern fonts, minor ticks, and LaTeX-rendered axis labels.

**Legend and colour encoding summary:**
- **Snapshot scatter**: viridis colormap → `log₁₀(m/M_tot)` (particle mass)
- **Energy plot**: blue line = `|ΔE/E|`, red line = `Q = T/|W|`, dashed grey = virial equilibrium (Q=0.5)
- **Particle count**: blue = `N` (bound particles), orange = `N_pairs` (KS binaries)
- **Lagrangian radii**: coloured lines labelled `M(r)/M_tot = X%` (1%, 10%, 50%, 90%, 100%)
- **HR diagrams**: colour-coded by BSE stellar type K* (MS=royal blue, GB=orange, BH=black, etc.)

For detailed legend documentation and colour tables, see [docs/src/manual.md § 9](docs/src/manual.md).

## N-body Units

The simulation uses Heggie & Mathieu (1986) N-body units:

| Quantity  | N-body unit      | Physical conversion         |
|-----------|------------------|-----------------------------|
| G         | 1                | —                           |
| M_total   | 1                | `ZMBAR × N` M☉             |
| E_total   | −1/4             | —                           |
| Length    | 1 (virial radius) | `RBAR` pc                  |
| Time      | 1                | `TSCALE` Myr               |
| Velocity  | 1                | `VSTAR` km/s               |

Use `extract_scaling(diagnostics)` to get a `UnitScaling` struct for conversions.

## Testing

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Benchmarks

```bash
julia --project=. benchmark/benchmarks.jl
```

## Dependencies

| Package        | Purpose                              |
|----------------|--------------------------------------|
| CairoMakie     | Publication-quality plotting & GIFs  |
| LaTeXStrings   | LaTeX-rendered axis labels & titles  |
| HDF5           | H5Part snapshot reading              |
| ProgressMeter  | Progress bars for batch reads        |
| SpecialFunctions | Error function for King model       |
| TOML           | Configuration parsing (stdlib)       |

## Platform Support

Tested on:
- **Fedora 42+** (DNF 5, environment modules for MPI)
- **Ubuntu 20.04+** (apt, standard MPI/HDF5 packages)

The install phase auto-detects the platform and adjusts HDF5 linking flags accordingly.

## Input Files

Custom simulation input files live in `input_files/` at the project root (separate from the backend's `examples/` directory). These use the Nbody6++GPU Fortran NAMELIST format with blocks `&INNBODY6`, `&ININPUT`, `&INDATA`, `&INSTAR`.

See [docs/src/input_files.md](docs/src/input_files.md) for a comprehensive reference on the `.inp` file format, KZ option flags, physical scaling, and a guide for designing new simulations.

## Documentation

Built with [Documenter.jl](https://documenter.juliadocs.org/). To build locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate(io=devnull)'
julia --project=docs docs/make.jl
```

Output: `docs/build/index.html`

Source documents:
- [docs/src/manual.md](docs/src/manual.md) — Full user manual
- [docs/src/input_files.md](docs/src/input_files.md) — Input file format reference & simulation design guide

## License

[MIT](LICENSE)
