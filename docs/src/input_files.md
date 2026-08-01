# Input File Reference

The project's input files live in `input_files/` at the package root. Two kinds exist:

- **`.inp`** — Fortran NAMELIST input files for single-cluster Nbody6++ runs (used via `simulation.input_file`)
- **`.toml`** — merger IC configurations for the Julia IC generator (used via `merger.config_file` or `run_merger_pipeline`)

## Directory contents

| File | Purpose |
|---|---|
| `N1k_quick.inp` | 1k-particle smoke test; isolated, no binaries, TCRIT=5 NB, Level C |
| `N5k_medium.inp` | 5k-particle medium test; isolated, no binaries, TCRIT=10 NB |
| `N25k_production.inp` | 25k production open cluster; Z=0.001, 200 binaries, point-mass tide at 13.3 kpc, TCRIT=200 NB |
| `N100k_production.inp` | 100k production open cluster; Z=0.001, 500 binaries, tidal field, TCRIT=100 NB |
| `imbh_runaway.inp` | IMBH formation via runaway collisions; 100k, ultra-dense (RBAR=0.5 pc), no binaries, isolated, TCRIT=20 NB |
| `gc_bh_subsystem.inp` | Globular-cluster BH subsystem; 100k, Z=0.0002, 2500 binaries, tide at 8 kpc, TCRIT=2000 NB |
| `tidal_tails.inp` | Tidal stripping near the Galactic centre; 50k, Z=0.02, tide at 2 kpc, TCRIT=500 NB |
| `young_massive_binaries.inp` | Binary-rich young massive cluster; 50k, 50% binaries (NBIN0=12500), isolated, TCRIT=100 NB |
| `pop3_cluster.inp` | Population III cluster; 50k, Z=1e-8, top-heavy IMF (ALPHAS=1.0, 8–300 M☉), TCRIT=200 NB |
| `merger_demo_small.toml` | 2×1000 King clusters, Kepler orbit; runs in seconds — full-pipeline demo, TCRIT=5 NB |
| `merger_equal_mass.toml` | 2×50000 King clusters, q=1 production merger (core-merger / IMBH science case) |
| `merger_minor_plummer.toml` | Plummer minor merger, q=0.1 (100k primary + 10k satellite); dynamical-friction inspiral |
| `merger_triple_cluster.toml` | 3×30000 King clusters, explicit triangular infall configuration |
| `merger_3cluster_small.toml` | 3×800 clusters (2 King + 1 Plummer), explicit triangle — small demo |
| `merger_5cluster_small.toml` | 5×500 clusters, pentagon layout with inward velocities — small demo |
| `N10k_long.inp` | 10k-body single cluster, extended TCRIT (long-run variant of the upstream example) |
| `merger_27cluster_cubic.toml` | 27×1000 clusters on a 3×3×3 cubic grid, inward velocities; stress test for the many-cluster plot paths |
| `verif_triorbit.toml` | Verification: 3 equal clusters on a rotating Lagrange-equilibrium triangle (seed 7; see below) |
| `verif_3d5cluster.toml` | Verification: 5 clusters distributed out of the z=0 plane; exercises xz/yz projections and 3D COM tracking (seed 13) |

To use an `.inp` file, set in `config.toml`:

```toml
[simulation]
run_test   = true
input_file = "../../input_files/N25k_production.inp"
```

To use a merger TOML:

```toml
[merger]
enabled     = true
config_file = "input_files/merger_demo_small.toml"
```

or directly `run_merger_pipeline("input_files/merger_demo_small.toml")`.

---

## Single-cluster `.inp` format (Fortran NAMELIST)

An `.inp` file is a sequence of NAMELIST blocks, each starting with `&BLOCKNAME` and ending with `/`, read in the order expected by `nbody6.F → start.F`:

| Block | Purpose |
|---|---|
| `&INNBODY6` | Run control: `KSTART` (1 = new run), `TCOMP` (CPU limit), `TCRTP0` (wall-clock limit, s) |
| `&ININPUT` | Core physics: `N`, `NRAND` (seed), `NNBOPT`, timestep accuracies `ETAI`/`ETAR`, output intervals `DTADJ`/`DELTAT`, end time `TCRIT`, physical scales `RBAR`/`ZMBAR`, the `KZ(1:50)` option array, tolerances, and the stellar-evolution `Level` (`'C'` = Kamlah et al. 2022, recommended) |
| `&INSSE` / `&INBSE` / `&INCOLL` | SSE/BSE/collision overrides (usually empty — Level defaults apply) |
| `&INDATA` | IMF (`ALPHAS`, `BODY1`, `BODYN`), binaries (`NBIN0`), metallicity `ZMET`, `DTPLOT` (sev.83 interval) |
| `&INSETUP` | External-IC placeholder block |
| `&INSCALE` | Initial virial ratio `Q`, rotation, tidal radius override |
| `&INXTRNL0` | External tidal field (only read when `KZ(14) > 0`) |
| `&INBINPOP` / `&INHIPOP` | Primordial binary / hierarchy populations (when `NBIN0`/`NHI0` > 0) |

KZ flags most relevant to this project: `KZ(3)` conf.3 snapshot output, `KZ(7)=3` Lagrangian radii (`lagr.7`), `KZ(12)=1` HR diagnostics (`sev.83_*`), `KZ(14)` tidal field (0 = isolated, 2 = point-mass galaxy), `KZ(19)=3` stellar evolution, `KZ(22)` initial conditions (0 = internal model, 2 = read `dat.10` in NB units, 10 = `dat.10` in astrophysical units), `KZ(23)` escaper removal (`esc.11`), `KZ(46)` HDF5 output (produces `snap.40_*.h5part` — **not readable** by this package; keep conf.3 output enabled).

For the full option catalogue see the Nbody6++ manual (Khalisi & Spurzem) and Wang et al. (2015).

---

## Merger TOML schema

`load_merger_config` accepts **two interchangeable schemas** for each `[merger.clusterN]` table: the flat legacy form and the structured form. Both may appear in the same file (but not mixed within one cluster table).

### Top-level `[merger]` keys

| Key | Type | Default | Description |
|---|---|---|---|
| `n_clusters` | Int | `2` | Number of `[merger.clusterN]` tables to read (N = 1…n_clusters, all required) |
| `orbit_mode` | String | `"kepler"` | `"kepler"` (exactly 2 clusters, placement auto-computed) or `"explicit"` (≥ 2 clusters, per-cluster `position`/`velocity` required) |
| `seed` | Int | *absent* | RNG seed; see [Seed semantics](#seed-semantics) |

### `[merger.orbit]` (Kepler mode only)

| Key | Type | Default | Description |
|---|---|---|---|
| `apocentre` | Float | `15.0` | Apocentre separation [pc] |
| `eccentricity` | Float | `0.7` | Orbital eccentricity, `0 ≤ e < 1` |

In `"explicit"` mode this section is ignored (a warning is emitted if present).

### `[merger.output]`

| Key | Type | Default | Description |
|---|---|---|---|
| `format` | String | `"nbody"` | `dat.10` in N-body units, `KZ(22)=2` (only supported value) |
| `truncate_jacobi` | Bool | `true` | Truncate each cluster at its nearest-neighbour Jacobi radius before combining |
| `output_dir` | String | `"."` | Output directory (overridden by the pipeline, which writes into `runs/merger_.../output/`) |
| `tcrit` | Float | `100.0` | Simulation end time (NB units for `"nbody"`) |
| `dtadj` | Float | `1.0` | ADJUST diagnostic interval |
| `deltat` | Float | `1.0` | Snapshot (conf.3) interval |

### Flat legacy cluster form

Fully annotated example (Kepler mode):

```toml
[merger]
n_clusters = 2
orbit_mode = "kepler"       # 2 clusters placed at apocentre automatically
seed       = 42             # optional; omit for a random seed

[merger.cluster1]
model      = "king"         # density profile: "king" or "plummer"
N          = 5000           # bodies to sample
W0         = 6.0            # King central potential (ignored for "plummer", with a warning)
rbar       = 2.0            # target half-mass radius [pc]
imf        = "kroupa"       # "kroupa", "kroupa_rescaled", or "equal"
bodyn      = 0.08           # lower IMF bound [M☉]
body1      = 100.0          # upper IMF bound [M☉]
mass_total = 5e4            # optional — see IMF resolution below

[merger.cluster2]
model      = "king"
N          = 5000
W0         = 6.0
rbar       = 2.0
imf        = "kroupa"
bodyn      = 0.08
body1      = 100.0
mass_total = 5e4

[merger.orbit]
apocentre    = 12.0         # [pc]
eccentricity = 0.6

[merger.output]
format          = "nbody"
truncate_jacobi = true
tcrit           = 10.0
dtadj           = 0.5
deltat          = 0.5
```

**IMF resolution in the flat form:**

- `imf = "kroupa"` **without** `mass_total` → natural Kroupa sampling in `[bodyn, body1]`; the total mass is an *output* (the sum of the samples)
- `imf = "kroupa"` **with** `mass_total` → rescaled Kroupa ("super-particle" mode): samples are uniformly rescaled so the sum equals `mass_total`. A warning is emitted at sample time when the rescale factor falls outside ×[0.7, 1.4], since individual body masses then no longer correspond to real stars and stellar-evolution output is non-physical
- `imf = "kroupa_rescaled"` → same as above but explicit (silences nothing by itself, but documents intent); requires `mass_total`
- `imf = "equal"` → all bodies get mass `mass_total / N`; requires `mass_total`

In `"explicit"` orbit mode each cluster table additionally requires:

```toml
position = [x, y, z]        # centre-of-mass position [pc]
velocity = [vx, vy, vz]     # centre-of-mass velocity in code units (G = 1 with
                            # M☉ and pc bases; 1 code unit ≈ 0.0656 km/s)
```

### Structured cluster form (preferred for new configs)

Profile and IMF become inline TOML tables with a `type` discriminator:

```toml
[merger]
n_clusters = 2
orbit_mode = "kepler"
seed       = 42

[merger.cluster1]
N       = 1500
rbar    = 1.0                                        # half-mass radius [pc]
profile = { type = "king", W0 = 6.0 }                # or { type = "plummer" }
imf     = { type = "kroupa", bounds = [0.08, 100.0] }

[merger.cluster2]
N       = 1500
rbar    = 1.0
profile = { type = "plummer" }
# rescaled Kroupa: explicit target mass (super-particle mode)
imf     = { type = "kroupa_rescaled", bounds = [0.08, 100.0], target_mass = 1.5e4 }
# equal-mass alternative:
# imf   = { type = "equal", particle_mass = 10.0 }   # or target_mass (→ /N)

[merger.orbit]
apocentre    = 12.0
eccentricity = 0.6

[merger.output]
format          = "nbody"
truncate_jacobi = true
tcrit           = 10.0
```

Accepted table keys:

- `profile`: `type = "king"` (+ `W0`, default 6.0) or `type = "plummer"`
- `imf`: `type = "kroupa" | "kroupa_rescaled" | "equal"`; mass bounds either as `bounds = [lo, hi]` or as separate `bodyn`/`body1` keys (defaults 0.08 / 100.0); `kroupa_rescaled` requires `target_mass`; `equal` requires `particle_mass` **or** `target_mass` (divided by `N`)

The metadata file `merger_ic.toml` written next to `dat.10` stores cluster specs in this structured form, and the same parser reads it back — see `load_merger_ic_result`.

### Seed semantics

- `seed` **omitted** (or the legacy sentinel `0`) → a random seed is drawn at generation time. The IC is still reproducible after the fact: the effective seed is recorded in `merger_ic.toml` (`meta.seed`)
- `seed = n` (nonzero) → deterministic sampling with `MersenneTwister(n)`
- The effective seed also feeds the `NRAND` parameter of the generated `merger.inp`, so the Nbody6++ run inherits it
- If a caller passes an explicit `rng` to `generate_merger_ic`, sampling is **not** reproducible from the recorded seed; the metadata records this honestly via `meta.external_rng = true` (the seed still feeds `NRAND`)

---

## `verif_triorbit.toml` — Lagrange-equilibrium triangle

Three equal-mass clusters (`m = 1.5×10⁴` M☉ each) on an equilateral triangle at radius `r = 6` pc from the centre of mass, given tangential (counter-clockwise) velocities so the configuration rotates rigidly as a bound Lagrange solution instead of falling inward.

For three equal masses `m` on an equilateral triangle of circumradius `r`, the net force on each body points at the centre with magnitude `G m² / (√3 r²)`; setting this equal to `m ω² r` gives the equilibrium angular velocity

```math
\omega^2 = \frac{G m}{\sqrt{3}\, r^3},
\qquad
\mathbf{v}_i = \omega\, (-y_i,\ x_i,\ 0).
```

With `G = 1` (code units), `m = 1.5e4`, `r = 6`: `ω ≈ 6.332` and `|v| = ω r ≈ 37.99` code units (≈ 2.49 km/s) per cluster — exactly the `velocity` entries in the file.

Expected behaviour: the three cluster COMs trace a slowly rotating triangle while each cluster relaxes internally. This verifies the COM-trajectory tracking and the per-cluster virial diagnostic (`per_cluster_virial`) for *separated, stable* clusters; the config's equilibrium condition is also asserted in the test suite.

---

## References

- Aarseth, S.J. (2003). *Gravitational N-Body Simulations*. Cambridge University Press.
- Wang, L. et al. (2015). MNRAS 450, 4070. — Nbody6++GPU code paper.
- Kamlah, A.W.H. et al. (2022). MNRAS 511, 4060. — Level C stellar evolution.
- Kroupa, P. (2001). MNRAS 322, 231. — IMF.
- King, I.R. (1966). AJ 71, 64. — King models.
- Khalisi, E. & Spurzem, R. — Nbody6++ manual (Heidelberg).
