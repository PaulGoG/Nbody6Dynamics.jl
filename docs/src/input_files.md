# Nbody6++GPU Input File Reference

This document describes the Fortran NAMELIST input format used by
Nbody6PPGPU-beijing and how to design new `.inp` files for custom
simulations.

## File structure

An `.inp` file is a sequence of Fortran NAMELIST blocks.  Each block
starts with `&BLOCKNAME` and ends with `/`.  Parameters are
comma-separated `KEY=VALUE` pairs.  Omitted parameters keep their
compiled defaults.  Fortran comments (`!`) are allowed outside blocks.

The blocks, **in order**, are:

| Block | Purpose |
|---|---|
| `&INNBODY6` | Run control: start mode, CPU limits, output intervals |
| `&ININPUT` | Core physics: N, timesteps, KZ options, tolerances |
| `&INSSE` | Single-star evolution overrides (usually empty) |
| `&INBSE` | Binary-star evolution overrides (usually empty) |
| `&INCOLL` | Collision parameters (usually empty) |
| `&INDATA` | IMF, metallicity, binary count, snapshot interval |
| `&INSETUP` | External initial conditions (for KZ(22)>=2) |
| `&INSCALE` | Virial ratio, rotation, tidal radius |
| `&INXTRNL0` | External tidal field parameters |
| `&INBINPOP` | Primordial hard binary population |
| `&INHIPOP` | Primordial hierarchical (wide) binary population |

Empty blocks are written as `&BLOCKNAME /` and are required for
format compliance even when no overrides are needed.


## Block-by-block reference

### `&INNBODY6` — Run control

| Parameter | Type | Description |
|---|---|---|
| `KSTART` | Int | Start mode: 1 = new run, 2 = restart from COMMON, 3/4/5 = restart variants |
| `TCOMP` | Real | Maximum CPU time (seconds) before graceful termination |
| `TCRTP0` | Real | Wall-clock time limit (seconds); the code checks this periodically |
| `isernb` | Int | Block-step parameter for irregular force (typically 40) |
| `iserreg` | Int | Block-step parameter for regular force (typically 40) |
| `iserks` | Int | KS block-step parameter (0 = use default) |

**Practical notes:**
- Set `TCOMP=1.0E8` (effectively infinite) and control runtime via `TCRTP0`.
- `TCRTP0` is wall-clock seconds: 3600 = 1 hour, 86400 = 1 day.
- For test runs use `TCRTP0=3600`; for production, 86400-259200 (1-3 days).


### `&ININPUT` — Core physics

#### Particle and run parameters

| Parameter | Type | Description |
|---|---|---|
| `N` | Int | Total number of particles (single stars + binary centers of mass) |
| `NFIX` | Int | Output frequency control (1 = standard) |
| `NCRIT` | Int | Neighbour list membership criterion (typically 10) |
| `NRAND` | Int | Random seed (change between runs for different realisations) |
| `NNBOPT` | Int | Optimal neighbour number (80-100 for N>10k) |
| `NRUN` | Int | Run number (always 1 for first run) |
| `NCOMM` | Int | COMMON save frequency (10 = every 10 DTADJ) |

#### Integration parameters

| Parameter | Type | Description |
|---|---|---|
| `ETAI` | Real | Irregular timestep accuracy (0.01-0.02; smaller = more accurate but slower) |
| `ETAR` | Real | Regular timestep accuracy (0.01-0.02) |
| `RS0` | Real | Initial radius of neighbour sphere in NB-units (0.08-0.15) |
| `DTADJ` | Real | ADJUST output interval in NB-units (energy, particle counts, etc.) |
| `DELTAT` | Real | Snapshot output interval in NB-units (conf.3 files) |
| `TCRIT` | Real | Final integration time in NB-units |
| `QE` | Real | Energy tolerance for restart (1.0 = no restart on energy error) |
| `RBAR` | Real | Half-mass radius in parsecs (sets physical length scale) |
| `ZMBAR` | Real | Mean mass in solar masses (sets physical mass scale) |

**Units:** Nbody6 uses Henon N-body units internally.  The physical
scaling is set by `RBAR` (pc) and `ZMBAR` (Msun).  The time scale
`TSCALE` (in Myr per NB-unit) is derived from these.  A cluster with
`RBAR=1.0, ZMBAR=0.6` has `TSCALE ~ 1 Myr/NB-unit`.

**Output intervals:**
- `DTADJ`: Controls how often diagnostics (energy, virial ratio, N,
  N_pairs) are written to stdout.  Use 0.5-2.0.
- `DELTAT`: Controls conf.3 snapshot output.  Use 1.0-5.0 for dense
  output, 10-50 for long integrations.  `DELTAT=DTPLOT` by convention.

#### Tolerances

| Parameter | Type | Description |
|---|---|---|
| `DTMIN` | Real | Minimum timestep (2.5E-6 typical) |
| `RMIN` | Real | Minimum two-body distance for KS regularisation (8E-5 typical) |
| `ETAU` | Real | KS regularisation accuracy (0.1) |
| `ECLOSE` | Real | Binding energy criterion for KS termination (1.0) |
| `GMIN` | Real | Minimum relative perturbation for KS (1E-6) |
| `GMAX` | Real | Maximum relative perturbation for unperturbed KS (0.01) |
| `SMAX` | Real | Maximum step for search of next apocentre (0.5-1.0) |

#### Stellar evolution level

| Parameter | Type | Description |
|---|---|---|
| `Level` | Char | Stellar evolution prescription level: `'A'`, `'B'`, or `'C'` |

- **Level A:** Standard Hurley SSE/BSE (2000/2002).
- **Level B:** Updated prescriptions (Banerjee et al. 2020).
- **Level C:** Kamlah et al. (2022) — most up-to-date.  Includes
  updated remnant mass prescriptions, pair-instability supernovae,
  pulsational pair-instability, and fallback-modulated natal kicks.
  **Recommended for all new simulations.**


### KZ option flags

The `KZ(1:50)` array controls physics modules and output options.  Set
them as `KZ(1:10)= 1 1 1 0 1 ...` (space-separated within each group
of 10).

**Most important flags:**

| KZ | Values | Description |
|---|---|---|
| KZ(1) | 0/1/2 | COMMON save: 0=off, 1=end only, 2=also at TCOMP |
| KZ(2) | 0/1/2 | COMMON save + restart on energy error |
| KZ(3) | 0/1/2 | Basic data output to file: 1=conf.3, 2=conf.3+extended |
| KZ(5) | 0/1 | Initial density profile: 0=uniform sphere, 1=Plummer |
| KZ(7) | 0-4 | Lagrangian radii output: 3=output, 4=+center averaging |
| KZ(8) | 0/1/2/3 | Primordial binaries: 0=none, 1/2=eigenevolution, 3=Kroupa period |
| KZ(9) | 0/1/2/3 | Stellar density centre: 3=density-weighted |
| KZ(12) | 0/1 | HR diagnostics output (sev.83 files): 1=enabled |
| KZ(14) | 0/1/2/3 | External tidal field: 0=none, 2=point-mass galaxy, 3=disk+halo |
| KZ(15) | 0/1/2 | Triple, quad, chain regularisation: 2=full |
| KZ(19) | 0/1/2/3 | Stellar evolution: 0=off, 3=Hurley SSE/BSE |
| KZ(20) | 0-7 | IMF: 6=Kroupa (2001) random pairing, 7=Kroupa corrected |
| KZ(22) | 0/1/2 | Initial conditions: 0=King model internal, 2=read dat.10 |
| KZ(23) | 0/1/2 | Escaper removal: 1=remove at 2*r_t, 2=with tidal tail tracking |
| KZ(25) | 0/1/2 | WD kicks: 2=enabled |
| KZ(26) | 0/1/2 | Slow-down regularisation for hard binaries: 2=enabled |
| KZ(27) | 0/1/2/3 | Tidal circularisation: 2=enabled, 3=Hurley |
| KZ(30) | 0/1/2 | Chain regularisation: 2=full |
| KZ(46) | 0/1/2/3/4 | Output format: 1=HDF5, 2=CSV, 3=bdat.9, 4=conf.3+bdat.9 |

**Standard KZ template for production runs:**
```
KZ(1:10)= 1 1 1 0 1 0 3 2 3 2
KZ(11:20)=0 1 0 2 2 0 0 0 3 6
KZ(21:30)=1 1 2 0 2 2 3 2 0 2
KZ(31:40)=1 0 2 2 1 0 1 1 2 1
KZ(41:50)=0 0 0 0 0 4 3 0 3 0
```

Key differences from the template for specific science cases:
- **No tidal field (isolated):** Set `KZ(14)=0`.
- **No primordial binaries:** Set `KZ(8)=0` and `NBIN0=0`.
- **HDF5 output:** Set `KZ(46)=1` (requires HDF5-enabled build).


### `&INDATA` — IMF and stellar population

| Parameter | Type | Description |
|---|---|---|
| `ALPHAS` | Real | IMF power-law slope: 2.35 = Salpeter/Kroupa upper end, 1.0 = flat (top-heavy) |
| `BODY1` | Real | Maximum stellar mass in Msun (100-300) |
| `BODYN` | Real | Minimum stellar mass in Msun (0.08 = hydrogen-burning limit) |
| `NBIN0` | Int | Number of primordial binaries (binary fraction ~ 2*NBIN0/N) |
| `NHI0` | Int | Number of primordial hierarchical triples (usually 0) |
| `ZMET` | Real | Metallicity Z (solar = 0.02, metal-poor GC = 0.0002, Pop III ~ 1E-8) |
| `EPOCH0` | Real | Initial age of stars in Myr (0 = zero-age main sequence) |
| `DTPLOT` | Real | Stellar evolution snapshot interval in NB-units (sev.83 files) |

**IMF notes:**
- `KZ(20)=6` with `ALPHAS=2.35` gives the standard Kroupa (2001) IMF
  with slopes 0.3 (0.08-0.5 Msun) and 2.35 (0.5-BODY1 Msun).
- For a top-heavy Pop III IMF, use `ALPHAS=1.0` with `BODYN=8.0`.
- `BODY1=150.0` is standard; use `BODY1=300.0` for very massive stars.

**Binary fraction:**
The binary fraction is approximately `f_b = 2*NBIN0/(N + NBIN0)`.
For example: N=50000, NBIN0=12500 gives f_b ~ 40%.  N=100000,
NBIN0=2500 gives f_b ~ 5%.


### `&INSCALE` — Virial ratio and rotation

| Parameter | Type | Description |
|---|---|---|
| `Q` | Real | Initial virial ratio: 0.5 = virial equilibrium, <0.5 = cold (collapsing), >0.5 = supervirial (expanding) |
| `VXROT` | Real | Rotation angular velocity about x-axis (0 = none) |
| `VZROT` | Real | Rotation angular velocity about z-axis (0 = none) |
| `RTIDE` | Real | Tidal radius override (0 = computed from tidal field) |


### `&INXTRNL0` — Tidal field

Active only when `KZ(14)>=2`.

**Point-mass galaxy (`KZ(14)=2`):**

| Parameter | Type | Description |
|---|---|---|
| `GMG` | Real | Galaxy mass in Msun (e.g., 1.78E11 for Milky Way-like) |
| `RG0` | Real | Galactocentric distance in kpc (circular orbit radius) |

**Disk + halo model (`KZ(14)=3`):**

| Parameter | Type | Description |
|---|---|---|
| `GMG` | Real | Bulge mass (Msun) |
| `DISK` | Real | Disk mass (Msun) |
| `A`, `B` | Real | Miyamoto-Nagai disk scale lengths (kpc) |
| `VCIRC` | Real | Circular velocity (km/s) |
| `RCIRC` | Real | Reference radius for VCIRC (kpc) |

For an isolated cluster (no tidal field), set `KZ(14)=0` and
`GMG=0.0, RG0=0.0`.

**Typical values for Milky Way-like galaxy:**
- Solar neighbourhood: `GMG=1.78E11, RG0=8.5`
- Inner disk (strong tide): `GMG=1.78E11, RG0=2.0`
- Outer halo (weak tide): `GMG=1.0E11, RG0=50.0`


### `&INBINPOP` — Primordial binary parameters

Active only when `NBIN0 > 0`.

| Parameter | Type | Description |
|---|---|---|
| `SEMI0` | Real | Minimum semi-major axis in AU (0.0005-0.001 typical) |
| `ECC0` | Real | Eccentricity distribution: -1.0 = thermal f(e)=2e |
| `RATIO` | Real | Mass ratio range parameter (1.0 = uniform q distribution) |
| `RANGE` | Real | Semi-major axis range: log10(a_max/a_min) = RANGE (5.0 typical) |
| `NSKIP` | Int | Frequency of binary assignment: 1 = every star can be binary, 3 = every 3rd, etc. |
| `IDORM` | Int | Dormant binaries: 0 = all active |

**Interpretation:** With `SEMI0=0.001, RANGE=5.0`, the semi-major axis
distribution is uniform in log(a) from 0.001 AU to 100 AU.


## Physical scaling

Nbody6 uses Henon N-body units internally where G=1, M_total=1,
E=-1/4.  The physical scaling is:

| Quantity | Formula | Example (RBAR=1.0, ZMBAR=0.6, N=100k) |
|---|---|---|
| Length | 1 NB = RBAR pc | 1.0 pc |
| Mass | 1 NB = N * ZMBAR Msun | 60,000 Msun |
| Time | TSCALE = sqrt(RBAR^3 / (G * N * ZMBAR)) | ~0.8 Myr |
| Velocity | VSTAR = RBAR / TSCALE | ~1.2 km/s |

The time scale determines the physical duration of your simulation.
With `TCRIT=100` and `TSCALE~0.8 Myr`, the physical time is ~80 Myr.


## Designing a new simulation

### Checklist

1. **Choose N** based on available compute:
   - N=1k-10k: minutes to hours (testing)
   - N=25k-50k: hours to 1 day (moderate production)
   - N=100k: 6-24 hours (our standard production)
   - N=500k+: days to weeks (heavy production)

2. **Set physical scales** (`RBAR`, `ZMBAR`):
   - Open cluster: RBAR=1-3 pc, ZMBAR=0.5-0.7
   - Globular cluster: RBAR=3-10 pc, ZMBAR=0.6
   - Ultra-dense (IMBH runs): RBAR=0.1-0.5 pc

3. **Choose IMF** (`ALPHAS`, `BODY1`, `BODYN`):
   - Standard Kroupa: ALPHAS=2.35, BODY1=100-150, BODYN=0.08
   - Top-heavy Pop III: ALPHAS=1.0, BODY1=300, BODYN=8.0

4. **Set metallicity** (`ZMET`):
   - Solar: 0.02
   - LMC/SMC: 0.004-0.008
   - Metal-poor GC: 0.0002
   - Pop III: 1E-8

5. **Configure binaries** (`NBIN0`, `SEMI0`, `ECC0`, `RANGE`):
   - No binaries: NBIN0=0
   - 5% hard: NBIN0=N/40
   - 50% total: NBIN0=N/4

6. **Set tidal field** (`KZ(14)`, `GMG`, `RG0`):
   - Isolated: KZ(14)=0
   - Point-mass: KZ(14)=2 with GMG, RG0

7. **Choose output intervals:**
   - `DTADJ`: 0.5-2.0 (diagnostics frequency)
   - `DELTAT` and `DTPLOT`: 1.0-20.0 (snapshot frequency)
   - Smaller = more output files; bigger = fewer snapshots
   - For animations, aim for 20-50 snapshots total

8. **Set time limit:**
   - `TCRIT`: integration time in NB-units
   - `TCRTP0`: wall-clock limit in seconds (safety cutoff)


## Custom input files

The custom simulation files for this project live in `input_files/`
at the project root (not in the backend's `examples/` directory).

| File | Science case | N | Z | Binaries | Tidal | Time |
|---|---|---|---|---|---|---|
| `N25k_production.inp` | Standard open cluster | 25k | 0.001 | 200 (1.6%) | Yes, 13.3 kpc | 200 NB |
| `N100k_production.inp` | Standard open cluster | 100k | 0.001 | 500 (1%) | Yes, 13.3 kpc | 100 NB |
| `imbh_runaway.inp` | IMBH via runaway collisions | 100k | 0.001 | None | None | 20 NB |
| `gc_bh_subsystem.inp` | GC black hole subsystem | 100k | 0.0002 | 2500 (5%) | Yes, 8 kpc | 2000 NB |
| `tidal_tails.inp` | Tidal stripping near GC | 50k | 0.02 | 2500 (10%) | Yes, 2 kpc | 500 NB |
| `young_massive_binaries.inp` | Binary-rich young cluster | 50k | 0.02 | 12500 (50%) | None | 100 NB |
| `pop3_cluster.inp` | Population III (primordial) | 50k | 1E-8 | 250 (1%) | Yes, 13.3 kpc | 200 NB |

To use any of these, set in `config.toml`:
```toml
[simulation]
run_test   = true
input_file = "../../input_files/imbh_runaway.inp"
```

Then run:
```
julia --project=. -e 'using Nbody6Setup; run_pipeline(load_config("config.toml"))'
```


## References

- Aarseth, S.J. (2003). *Gravitational N-Body Simulations*. Cambridge University Press.
- Wang, L. et al. (2015). MNRAS 450, 4070. — Nbody6++GPU code paper.
- Kamlah, A.W.H. et al. (2022). MNRAS 511, 4060. — Level C stellar evolution, DRAGON-II.
- Vergara, M.C. et al. (2025). A&A. — IMBH via runaway collisions.
- Arca Sedda, M. et al. (2026). arXiv:2603.29657. — Pop III cluster simulations.
- Banerjee, S. (2021). MNRAS 503, 3371. — Binary-rich young massive clusters.
- Khalisi, E. & Spurzem, R. — Nbody6++ manual (Heidelberg).
