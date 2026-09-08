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

### `[merger.nbody6]`

Integration parameters written to `merger.inp`, in N-body units of the combined system (length unit `RBAR`, mass unit `M_total`). A zero for a derivable key means "derive from the member clusters at generation time": with `r_h` the smallest member half-mass radius in those units, `N_min` the smallest post-truncation membership, and `ρ̂` the central density contrast of that member's profile, the rules are those of the engine's own `adjust.F` evaluated for the member cluster instead of the whole configuration.

| Key | Type | Default | Description |
|---|---|---|---|
| `qe` | Float | `2.0e-4` | Energy-error tolerance per adjustment interval (`QE`); must be > 0 |
| `etai` | Float | `0.02` | Irregular time-step factor; must be > 0 |
| `etar` | Float | `0.02` | Regular time-step factor; must be > 0 |
| `nnbopt` | Int | `0` | Target neighbour number; `0` = `clamp(round(√N_total), 20, 300)`; must be ≥ 0 |
| `rs0` | Float | `0.0` | Initial neighbour-sphere radius; `0` = `2 r_h (2 NNBOPT / N_min)^{1/3}` capped at `r_h` (the factor 2 follows the engine's example inputs; a smaller radius left outer stars without neighbours and hung the engine's start-up on one of three realisations); must be ≥ 0 and, when set, no larger than the smallest member half-mass radius |
| `rmin` | Float | `0.0` | KS regularisation distance; `0` = `4 r_h / (N_min ρ̂^{1/3})`; must be ≥ 0 |
| `dtmin` | Float | `0.0` | KS time-step threshold; `0` = `0.04 √(ETAI/0.02) √(RMIN³ N_total)`; must be ≥ 0 |
| `kz16` | Int | `0` | `KZ(16)`: the engine's re-derivation of `RMIN`, `DTMIN`, and `ECLOSE` from its global scale radius and core density every `DTADJ`; `0` keeps the written values (recommended for multi-cluster systems); one of 0, 1, 2, 3 |
| `etau` | Float | `0.1` | Regularised time-step factor; must be > 0 |
| `eclose` | Float | `1.0` | Binding energy per unit mass of a hard binary; must be > 0 |
| `gmin` | Float | `1.0e-6` | Relative perturbation for unperturbed KS motion; `0 < gmin < gmax` |
| `gmax` | Float | `0.01` | Termination parameter for soft KS binaries; `> gmin` |
| `smax` | Float | `1.0` | Maximum time step, a power of two commensurate with 1; must be > 0 |
| `tcomp` | Float | `1.0e8` | Run-time limit [s]; must be > 0 |
| `tcrtp0` | Float | `3600.0` | Termination time [Myr]; must be > 0 |
| `isernb`, `iserreg`, `iserks` | Int | `40`, `40`, `0` | MPI block-size thresholds below which irregular, regular, and KS blocks run serially; must be ≥ 0 |
| `nfix` | Int | `1` | Multiplier of `deltat` for `conf.3` and binary output; must be ≥ 1 |
| `ncrit` | Int | `10` | Minimum particle number, alternative termination criterion; must be ≥ 1 |
| `nrun` | Int | `1` | Run identification index; must be ≥ 1 |
| `ncomm` | Int | `10` | Multiplier of `deltat` for the restart (`COMMON`) dump interval; must be ≥ 1 |

`[merger.nbody6.kz]` holds explicit `KZ(i) = v` overrides with string keys `"1"`–`"50"` and integer values, applied after every named option (an override of an index that also has a named key, 14, 16, or 19, warns):

```toml
[merger.nbody6.kz]
"8" = 2      # primordial binaries from dat.10
"47" = 1
```

The resolved values appear in `merger_summary.txt` and in `merger_ic.toml` (`[nbody6]`, `[stellar]`), next to the combined virial ratio `Q = T/|W|`, the ratio `RBAR/r_hm,min`, the crossing times of the configuration and of the smallest member in N-body units and Myr, and the N-body time unit `T*` (`[meta]`). Generation warns when `Q < 0.3` (cold-collapse regime: global infall dominates and the engine's global diagnostics are meaningless until the remnant forms), when `RBAR/r_hm,min > 5` (members unresolved by the single-centre diagnostics), and when `deltat` exceeds the smallest member crossing time (member dynamics undersampled), and refuses an `rs0` wider than the smallest member.

### `[merger.stellar]`

| Key | Type | Default | Description |
|---|---|---|---|
| `kz19` | Int | `3` | `KZ(19)`, stellar evolution and mass-loss scheme: `0` off (HR diagnostics `KZ(12)` are then switched off too), `1`–`2` supernova schemes, `≥ 3` Eggleton–Tout–Hurley; must be ≥ 0 |
| `level` | String | `"C"` | SSE/BSE parameter level (Kamlah et al. 2022); one of `"A"`, `"B"`, `"C"`, `"0"` (no level: the engine's independent defaults) |
| `zmet` | Float | `0.001` | Metal abundance; `0.0001 ≤ zmet ≤ 0.03` (the engine's own bounds) |
| `epoch0` | Float | `0.0` | Formation time of the population [Myr]; must be ≤ 0 (the age at start is `−epoch0`) |
| `dtplot` | Float | `1.0` | Interval of the stellar-evolution diagnostics (`sev.83_*`) [NB]; must be > 0 and ≥ `deltat` |

### `[merger.tidal]`

| Key | Type | Default | Description |
|---|---|---|---|
| `kz14` | Int | `0` | `KZ(14)`: `0` isolated; `1` standard solar-neighbourhood linearised tide (no parameters); `2` point-mass galaxy on a circular orbit; `5` MWPotential2014 with the configuration on a galactocentric orbit. `3` and `4` are refused (see below); one of 0, 1, 2, 5 |
| `gmg` | Float | `0.0` | Galaxy mass [M☉] for `kz14 = 2`; must be > 0 |
| `rg0` | Float | `0.0` | Galactocentric distance of the circular orbit [kpc] for `kz14 = 2`; must be > 0 |
| `rg` | [Float] | `[0, 0, 0]` | Galactocentric position of the configuration's centre of mass [kpc] for `kz14 = 5`; non-zero |
| `vg` | [Float] | `[0, 0, 0]` | Galactocentric velocity [km/s] for `kz14 = 5`; non-zero |

Options `3` (point mass + Miyamoto–Nagai disk + logarithmic halo + bulge) and `4` (Plummer potential) are refused: on those paths the engine rescales every velocity to the `&INSCALE` virial ratio including the external potential (`xtrnl0.F`), which destroys the prescribed orbital kinematics of a multi-cluster configuration. The `&INSCALE` tidal radius stays `0` so the engine derives it from the field with the generator's `RBAR`; a non-zero value would override `RBAR`. In a tidal field the engine removes escapers on the distance criterion alone, and because it does not evaluate the tidal potential energy its energy check measures the tidal work: a tidal configuration requires `merger.nbody6.qe ≥ 0.01` (see the merger documentation).

The `&INDATA` mass bounds `BODY1`/`BODYN` are written from the clusters' IMF specifications (inert under `KZ(22) = 2`, where masses come from `dat.10`); `ALPHAS` is inert for the same reason. `NBIN0` is the number of primordial pairs (see `[merger.clusterN.binaries]`); `NHI0 = 0`, no hierarchies are generated.

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
- `imf = "kroupa"` **with** `mass_total` → rescaled Kroupa: samples are uniformly rescaled so the sum equals `mass_total`. The rescale factor `mass_total / (N ⟨m⟩)` is checked when the file is loaded: outside ×[0.7, 1.4] the loader warns, outside ×[0.5, 2] it refuses the configuration, because the bodies would no longer be stars while stellar evolution is always active in merger runs. Choose `N` and `mass_total` consistent with the IMF mean (about 0.58 M☉ over 0.08–100 M☉), drop `mass_total`, or use `imf = "equal"` for a collisionless super-particle model
- `imf = "kroupa_rescaled"` → same as above but explicit (documents intent; the same bounds apply); requires `mass_total`
- `imf = "equal"` → all bodies get mass `mass_total / N`; requires `mass_total`

In `"explicit"` orbit mode each cluster table additionally requires:

```toml
position = [x, y, z]        # centre-of-mass position [pc]
velocity = [vx, vy, vz]     # centre-of-mass velocity [km/s]
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

### `[merger.clusterN.binaries]` — primordial binaries

| Key | Type | Default | Description |
|---|---|---|---|
| `fraction` | Float | `0.0` | Binary fraction by systems, `N_b / (N_s + N_b)`; `0 ≤ fraction < 1` |
| `pairing` | String | `"random"` | `"random"`: both components drawn from the IMF, the more massive is the primary; `"uniform_q"`: the secondary's mass is `q m₁` with `q` uniform in `[q_min, 1]` (it replaces the drawn mass of the partner star) |
| `period` | String | `"kroupa1995"` | `"kroupa1995"`: the Kroupa (1995) birth period distribution `f(log P) ∝ (log P − 1)/(45 + (log P − 1)²)`, `1 ≤ log(P/d) ≤ 8.43`, converted to a semi-major axis with Kepler's third law; `"loguniform"`: semi-major axis log-uniform in `[a_min, a_max]` |
| `a_min`, `a_max` | Float | `0.01`, `100.0` | Semi-major axis bounds [AU] for `"loguniform"`; `0 < a_min < a_max` |
| `q_min` | Float | `0.1` | Lower mass-ratio bound for `"uniform_q"`; `0 < q_min < 1` |
| `eccentricity` | String | `"thermal"` | `"thermal"` (`f(e) = 2e`) or `"circular"` |

The density sampler places *systems* (a pair by its centre of mass), which are virialised and truncated as such; each pair is then expanded into two bodies on a Keplerian orbit with a random phase and orientation. The engine reads primordial pairs as bodies `2i − 1, 2i` for `i ≤ NBIN0`, so the pairs of all clusters are written first, cluster by cluster, followed by every cluster's singles; `NBIN0` and `KZ(8) = 2` are set in `merger.inp`, and a cluster's members are then two contiguous blocks recorded in `merger_ic.toml` (`cluster_blocks`) and recoverable from `merger_summary.txt` (`binaries: N` per cluster). The summary also reports the hard fraction of each cluster's pairs (`G m₁ m₂ / 2a > ⟨m⟩ σ²` of the cluster's systems). Binary-rich runs carry a larger energy error per adjustment interval (regularised pairs, chains): the two-cluster demo with 20 % binaries reached 4×10⁻³ over its first time unit, so set `merger.nbody6.qe` accordingly (10⁻³ to 10⁻²) or the engine halts at the first adjustment that exceeds the default 2×10⁻⁴.

### Seed semantics

- `seed` **omitted** (or the legacy sentinel `0`) → a random seed is drawn at generation time. The IC is still reproducible after the fact: the effective seed is recorded in `merger_ic.toml` (`meta.seed`)
- `seed = n` (nonzero) → deterministic sampling with `MersenneTwister(n)`
- The effective seed also feeds the `NRAND` parameter of the generated `merger.inp`, so the Nbody6++ run inherits it
- If a caller passes an explicit `rng` to `generate_merger_ic`, sampling is **not** reproducible from the recorded seed; the metadata records this honestly via `meta.external_rng = true` (the seed still feeds `NRAND`)

---

## `verif_triorbit.toml` — Lagrange-equilibrium triangle

Three equal-mass clusters (`m = 900` M☉ each: 1500 equal bodies of 0.6 M☉, so the equilibrium is exact and no massive star evolves within the run) on an equilateral triangle at radius `r = 6` pc from the centre of mass, given tangential (counter-clockwise) velocities so the configuration rotates rigidly as a bound Lagrange solution instead of falling inward.

For three equal masses `m` on an equilateral triangle of circumradius `r`, the net force on each body points at the centre with magnitude `G m² / (√3 r²)`; setting this equal to `m ω² r` gives the equilibrium angular velocity

```math
\omega^2 = \frac{G m}{\sqrt{3}\, r^3},
\qquad
\mathbf{v}_i = \omega\, (-y_i,\ x_i,\ 0).
```

With `G = 4.30091×10⁻³ pc (km s⁻¹)² M☉⁻¹`, `m = 900 M☉`, `r = 6 pc`: `ω ≈ 0.1017 km s⁻¹ pc⁻¹` and `|v| = ω r ≈ 0.6103 km s⁻¹` per cluster — exactly the `velocity` entries in the file.

Expected behaviour: the three cluster COMs trace a slowly rotating triangle while each cluster relaxes internally. This verifies the COM-trajectory tracking and the per-cluster virial diagnostic (`per_cluster_virial`) for *separated, stable* clusters; the config's equilibrium condition is also asserted in the test suite.

---

## References

- Aarseth, S.J. (2003). *Gravitational N-Body Simulations*. Cambridge University Press.
- Wang, L. et al. (2015). MNRAS 450, 4070. — Nbody6++GPU code paper.
- Kamlah, A.W.H. et al. (2022). MNRAS 511, 4060. — Level C stellar evolution.
- Kroupa, P. (2001). MNRAS 322, 231. — IMF.
- King, I.R. (1966). AJ 71, 64. — King models.
- Khalisi, E. & Spurzem, R. — Nbody6++ manual (Heidelberg).
