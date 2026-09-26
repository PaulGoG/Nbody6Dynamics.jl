# Input File Reference

The inputs ship under `input_files/` in the package tree, one folder per purpose. `example_input(name)` returns the absolute path of one for a project outside the checkout (`example_input("engine/N1k_quick.inp")`; a bare name such as `example_input("N1k_quick.inp")` is looked up across the folders). Four kinds of file exist:

- **engine inputs** (`.inp`, under `engine/`) — Fortran NAMELIST inputs of single-cluster Nbody6++ runs, used through `simulation.input_file`;
- **merger configurations** (`.toml`, under `mergers/` and `verification/`, and the merger TOMLs of `showcase/` and `gpu/`) — initial conditions for the Julia generator, used through `merger.config_file` or `run_merger_pipeline`;
- **pipeline configurations** (`*_pipeline.toml` under `showcase/` and `gpu/`) — complete configurations for `scripts/run_setup.jl`, with every path relative to their own folder;
- **sweep configurations** (`sweeps/`, `showcase/sweep.toml`) — a grid over a merger TOML for `scripts/run_sweep.jl`.

```
input_files/
├── engine/          # single-cluster engine inputs
├── mergers/         # merger initial-condition configurations
├── verification/    # the two targets of scripts/run_verif_suite.jl
├── sweeps/          # sweep configurations
├── showcase/        # five end-to-end cases (Showcase Cases)
└── gpu/             # CUDA-host validation and the probes above 5 × 10⁵ bodies
```

## Engine inputs

| File | Purpose |
|---|---|
| `engine/N1k_quick.inp` | 1k-particle smoke test; isolated, no binaries, TCRIT=5 NB, Level C |
| `engine/N5k_medium.inp` | 5k-particle medium test; isolated, no binaries, TCRIT=10 NB; the `single` target of the verification suite |
| `engine/N10k_long.inp` | 10k-body single cluster, extended TCRIT=50 NB (long-run variant of the upstream example) |
| `engine/N25k_production.inp` | 25k production open cluster; Z=0.001, 200 binaries, point-mass tide at 13.3 kpc, TCRIT=200 NB |
| `engine/N100k_production.inp` | 100k production open cluster; Z=0.001, 500 binaries, tidal field, TCRIT=100 NB |
| `engine/gc_bh_subsystem.inp` | Globular-cluster black-hole subsystem; 100k, Z=0.0002, 2500 binaries, tide at 8 kpc, TCRIT=2000 NB — [science case](@ref "Engine inputs: the science cases") |
| `engine/imbh_runaway.inp` | IMBH formation via runaway collisions; 100k, ultra-dense (RBAR=0.5 pc), no binaries, isolated, TCRIT=20 NB — science case |
| `engine/pop3_cluster.inp` | Population III cluster; 50k, Z=1e-8, top-heavy IMF (ALPHAS=1.0, 8–300 M☉), TCRIT=200 NB — science case |
| `engine/tidal_tails.inp` | Tidal stripping near the Galactic centre; 50k, Z=0.02, tide at 2 kpc, TCRIT=500 NB — science case |
| `engine/young_massive_binaries.inp` | Binary-rich young massive cluster; 50k, 50 % binaries (NBIN0=12500), isolated, TCRIT=100 NB — science case |

Every shipped `.inp` sets `QE = 1.0E-02`: with `KZ(2) = 1` the engine halts when the relative energy change over one adjustment interval exceeds `5 QE`, and these runs evolve stars (`KZ(19) = 3`), carry binaries, and mostly sit in a point-mass tidal field whose work the engine's energy check charges to the error, so a tolerance of the order of `1e-4` stops them on ordinary events (a KS termination halted the 1k smoke run at t = 3 of 5) while `1.0` would disable the check. At `1e-2` the check still catches an integration that goes wrong by 5 % per interval, which is what it is for.

To use an `.inp` file, point `simulation.input_file` at it, relative to the directory of the configuration (from the root `config.toml` of a checkout):

```toml
[simulation]
run_test   = true
input_file = "input_files/engine/N25k_production.inp"
```

## Merger configurations

| File | Purpose |
|---|---|
| `mergers/merger_demo_small.toml` | 2×1000 King clusters, Kepler orbit; runs in seconds — the full-pipeline demo of the root `config.toml`, TCRIT=5 NB |
| `mergers/merger_equal_mass.toml` | 2×50000 King clusters, q=1 production merger (core-merger / IMBH science case) |
| `mergers/merger_minor_plummer.toml` | Plummer minor merger, q=0.1 (100k primary + 10k satellite); dynamical-friction inspiral |
| `mergers/merger_triple_cluster.toml` | 3×30000 King clusters, explicit triangular infall configuration |
| `mergers/merger_3cluster_small.toml` | 3×800 clusters (2 King + 1 Plummer), explicit triangle — small demo |
| `mergers/merger_5cluster_small.toml` | 5×500 clusters, pentagon layout with inward velocities — small demo |
| `mergers/merger_27cluster_cubic.toml` | 27×1000 clusters on a 3×3×3 cubic grid, inward velocities; a cold-collapse stress test of the many-cluster figure paths, not a merger in the engine's regime (see Cluster Mergers, "Feasibility and limitations") |

To use a merger TOML:

```toml
[merger]
enabled     = true
config_file = "input_files/mergers/merger_demo_small.toml"
```

or directly `run_merger_pipeline("input_files/mergers/merger_demo_small.toml")`. The schema is documented below ([Merger TOML schema](@ref)).

## Verification configurations

| File | Purpose |
|---|---|
| `verification/verif_triorbit.toml` | Three equal clusters on a rotating Lagrange-equilibrium triangle (seed 7); the derivation is [below](@ref "`verif_triorbit.toml` — Lagrange-equilibrium triangle") |
| `verification/verif_3d5cluster.toml` | Five clusters distributed out of the z=0 plane; exercises the xz/yz projections and the three-dimensional COM tracking (seed 13) |

Both run through `scripts/run_verif_suite.jl` together with `engine/N5k_medium.inp`.

## Sweep configurations

`sweeps/sweep_demo.toml` is a small demonstration grid (orbital eccentricity × secondary size × two seeds) over `mergers/merger_demo_small.toml` with the root `config.toml` as its pipeline base; `showcase/sweep.toml` is the showcase sweep with isolated controls. The file format is in the manual, section Parameter sweeps.

## Showcase and GPU-host cases

The five cases of `showcase/`, their reference runs and the figures they produce are on the page [Showcase Cases](@ref). The files of `gpu/` (the CUDA and AVX builds, the 2 × 25 000-star validation merger and the probes above 5 × 10⁵ bodies) are described [below](@ref "The GPU-host and verification cases"); the recipe for a CUDA host is in the manual.

## Engine inputs: the science cases

Five of the engine inputs are science cases rather than tests: each sets up a cluster in which a specific dynamical process dominates, with the parameters that process needs. They are shipped as starting points, not as validated results: none of them has been run by the package's author, and the cost figures are estimates for the AVX build from the N-scaling of the shipped runs, with the wall-clock limit `TCRTP0` of each file set accordingly. The stellar-evolution level of every file is `C` [Kamlah2022](@cite).

### `gc_bh_subsystem.inp` — globular cluster with a black-hole subsystem

Metal-poor globular clusters retain a population of stellar-mass black holes that segregates to the core and forms a dynamically decoupled subsystem: it delays core collapse by acting as a central energy source, hardens black-hole binaries in three-body encounters and produces mergers, and is depleted by gravitational-wave recoil and dynamical ejection; once it is exhausted, the cluster undergoes a late core collapse [ArcaSedda2024](@cite), [Banerjee2022](@cite). Setup: N = 100 000, King `W0 = 6`, `Z = 0.0002`, 2500 primordial binaries, Kroupa IMF 0.08–100 M☉, point-mass galaxy at 8 kpc on a circular orbit, `TCRIT = 2000` NB (about 2 Gyr): the formation of the subsystem and its early depletion, not the 12 Gyr life of the cluster. Expected: mass segregation of the black holes within about 100 Myr, binary hardening and ejected mergers, gradual depletion, tidal stripping. Cost: days on 16–24 cores.

### `imbh_runaway.inp` — IMBH formation through runaway collisions

In ultra-dense young clusters (`r_h` of 0.1–0.5 pc) the most massive stars sink to the core within a Myr and collide repeatedly, building a very massive star that can collapse into an intermediate-mass black hole [Vergara2025](@cite), [ArcaSedda2023](@cite). Setup: N = 100 000, King `W0 = 9`, `r_h = 0.5` pc, `Z = 0.001`, no primordial binaries, Kroupa IMF 0.08–150 M☉, isolated, `TCRIT = 20` NB (about 20 Myr, the runaway phase). Expected: mass segregation of the O and B stars within 1–2 Myr, runaway mergers in the core, a very massive star and its collapse at a few Myr, exchange encounters that form an IMBH binary, hypervelocity ejections. Cost: hours on 16–24 cores.

### `pop3_cluster.inp` — Population III cluster

The first clusters formed from primordial gas with a top-heavy mass function, because the lack of metals suppresses cooling and fragmentation; with negligible winds, massive stars keep their mass to core collapse, which produces very massive black holes by direct collapse, pair-instability supernovae that leave no remnant, and pulsational pair instability [Wu2026](@cite). Setup: N = 50 000, King `W0 = 6`, `Z = 10⁻⁸`, flat IMF (`ALPHAS = 1.0`, equal mass per logarithmic bin) over 8–300 M☉, 250 primordial binaries, point-mass galaxy at 13.3 kpc, `TCRIT = 200` NB (about 500 Myr). Expected: rapid segregation, pair-instability explosions removing the most massive stars, a black-hole subsystem within a few Myr, black-hole mergers, dissolution as the supernova mass loss unbinds the system. Caveat: `Z = 10⁻⁸` lies below the metallicity range of the engine's stellar-evolution fits (10⁻⁴ to 0.03, [Hurley2000](@cite)), so the stars evolve at the edge of the tables. Cost: one to two days on 16–24 cores.

### `tidal_tails.inp` — tidal tails near the Galactic centre

Clusters in strong tidal fields lose stars through the Lagrange points and develop tails whose morphology, with its epicyclic overdensities, records the orbit and the mass-loss history [Kupper2010](@cite); clusters near the Galactic centre dissolve fastest [Park2018](@cite), and Palomar 5 is the archetype of a tail-dominated halo cluster [Gieles2021](@cite). Setup: N = 50 000, King `W0 = 5`, `Z = 0.02`, 2500 primordial binaries, Kroupa IMF 0.08–100 M☉, point-mass galaxy at 2 kpc on a circular orbit, `TCRIT = 500` NB (about 1 Gyr). Expected: stripping from the start, tails within about 50 Myr, preferential loss of low-mass stars, more than half the mass lost, core contraction as the tidal boundary shrinks. The engine's point-mass field (`KZ(14) = 2`) charges the tidal work to its energy check (see Cluster Mergers, "Escapers"); `QE = 0.01` accommodates it. Cost: one to two days on 16–24 cores.

### `young_massive_binaries.inp` — young massive cluster with a high binary fraction

More than seventy per cent of massive O stars are born in binaries [Sana2012](@cite); in a dense young cluster these drive mass transfer, stripped-envelope supernovae, X-ray binaries, double compact objects and runaway stars [Banerjee2021](@cite). Setup: N = 50 000 as 25 000 singles and 12 500 pairs (`NBIN0 = 12500`), King `W0 = 6`, `Z = 0.02`, Kroupa IMF 0.08–150 M☉, isolated, `TCRIT = 100` NB (about 50 Myr). Expected: early supernovae ejecting neutron stars, binary disruptions and runaway OB stars, blue stragglers from mass transfer, X-ray binaries, double compact objects within about 30 Myr, a diverse escaper population. Cost: half a day to a day on 16–24 cores.

---

## Single-cluster `.inp` format (Fortran NAMELIST)

An `.inp` file is a sequence of NAMELIST blocks, each starting with `&BLOCKNAME` and ending with `/`, read in the order expected by `nbody6.F → start.F`:

| Block | Purpose |
|---|---|
| `&INNBODY6` | Run control: `KSTART` (1 = new run), `TCOMP` (CPU limit), `TCRTP0` (wall-clock limit, s) |
| `&ININPUT` | Core physics: `N`, `NRAND` (seed), `NNBOPT`, timestep accuracies `ETAI`/`ETAR`, output intervals `DTADJ`/`DELTAT`, end time `TCRIT`, physical scales `RBAR`/`ZMBAR`, the `KZ(1:50)` option array, tolerances, and the stellar-evolution `Level` (`'C'` = [Kamlah2022](@cite), recommended) |
| `&INSSE` / `&INBSE` / `&INCOLL` | SSE/BSE/collision overrides (usually empty — Level defaults apply) |
| `&INDATA` | IMF (`ALPHAS`, `BODY1`, `BODYN`), binaries (`NBIN0`), metallicity `ZMET`, `DTPLOT` (sev.83 interval) |
| `&INSETUP` | External-IC placeholder block |
| `&INSCALE` | Initial virial ratio `Q`, rotation, tidal radius override |
| `&INXTRNL0` | External tidal field (only read when `KZ(14) > 0`) |
| `&INBINPOP` / `&INHIPOP` | Primordial binary / hierarchy populations (when `NBIN0`/`NHI0` > 0) |

KZ flags most relevant to this project: `KZ(3)` conf.3 snapshot output, `KZ(7)=3` Lagrangian radii (`lagr.7`), `KZ(12)=1` HR diagnostics (`sev.83_*`), `KZ(14)` tidal field (0 = isolated, 2 = point-mass galaxy), `KZ(19)=3` stellar evolution, `KZ(22)` initial conditions (0 = internal model, 2 = read `dat.10` in NB units, 10 = `dat.10` in astrophysical units), `KZ(23)` escaper removal (`esc.11`), `KZ(46)` HDF5 output (produces `snap.40_*.h5part` — **not readable** by this package; keep conf.3 output enabled).

For the full option catalogue see the Nbody6++ manual (Khalisi & Spurzem, Heidelberg, unpublished), the code papers [Spurzem1999](@cite), [NitadoriAarseth2012](@cite) and [Wang2015](@cite), the stellar-evolution levels of [Kamlah2022](@cite), and the method review of [SpurzemKamlah2023](@cite).

---

## Merger TOML schema

`load_merger_config` accepts **two interchangeable schemas** for each `[merger.clusterN]` table: the flat form and the structured form. Both may appear in the same file (but not mixed within one cluster table).

### Top-level `[merger]` keys

| Key | Type | Default | Description |
|---|---|---|---|
| `n_clusters` | Int | `2` | Number of `[merger.clusterN]` tables to read (N = 1…n_clusters, all required) |
| `orbit_mode` | String | `"kepler"` | `"kepler"` (exactly 2 clusters, placement auto-computed) or `"explicit"` (≥ 2 clusters, per-cluster `position`/`velocity` required) |
| `seed` | Int | *absent* | RNG seed; see [Seed semantics](#seed-semantics) |
| `virial_max_n` | Int | `200000` | Largest N, per cluster and for the combined system, for which the exact O(N²) potential is evaluated (virialisation of each cluster, combined virial ratio). A cluster above it is refused at load, naming the key; raise the value deliberately for large initial conditions. The pair sum is threaded and costs ≈ N²/2 evaluations: 11 s at N = 2×10⁵ on twelve threads, quadratic from there |

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
| `output_dir` | String | `""` | Output directory; a relative path resolves against the TOML's own directory. Empty leaves the choice to the caller: `run_merger_pipeline` creates `runs/merger_<timestamp>/` under the working directory, and the main pipeline always writes into `runs/merger_.../output/` |
| `tcrit` | Float | `100.0` | Simulation end time (NB units for `"nbody"`) |
| `dtadj` | Float | `1.0` | ADJUST diagnostic interval |
| `deltat` | Float | `1.0` | Snapshot (conf.3) interval |
| `tcrit_myr`, `dtadj_myr`, `deltat_myr` | Float | `0.0` | The same three in Myr; a positive value replaces the NB one and is converted at generation with the realised time unit `T*` (both forms of one key: error). `dtadj`, `deltat` and the stellar `dtplot` are written as the nearest dyadic rational with an exact decimal expansion (`engine_interval`, change below 0.4 %), because the engine counts their decimal digits with a loop that never terminates on other values |

### `[merger.nbody6]`

Integration parameters written to `merger.inp`, in N-body units of the combined system (length unit `RBAR`, mass unit `M_total`). A zero for a derivable key means "derive from the member clusters at generation time": with `r_h` the smallest member half-mass radius in those units, `N_min` the smallest post-truncation membership, and `ρ̂` the central density contrast of that member's profile, the rules are those of the engine's own `adjust.F` evaluated for the member cluster instead of the whole configuration.

| Key | Type | Default | Description |
|---|---|---|---|
| `qe` | Float | `2.0e-4` | Energy-error tolerance per adjustment interval (`QE`); must be > 0 |
| `etai` | Float | `0.02` | Irregular time-step factor; must be > 0 |
| `etar` | Float | `0.02` | Regular time-step factor; must be > 0 |
| `nnbopt` | Int | `0` | Target neighbour number; `0` = `clamp(round(√N_total), 20, 300)`; must be ≥ 0 |
| `rs0` | Float | `0.0` | Initial neighbour-sphere radius; `0` = `2 r_h (2 NNBOPT / N_min)^{1/3}` capped at `r_h` (the factor 2 follows the engine's example inputs; the engine regrows an empty sphere itself, so the value sets start-up cost, not correctness); must be ≥ 0 and, when set, no larger than the smallest member half-mass radius |
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

Runs with stellar evolution need the relaxed tolerance `qe ≈ 1e-2`: supernova kicks change the energy budget between adjustments, and the default tolerance halts the run.

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
| `kz19` | Int | `3` | `KZ(19)`, stellar evolution and mass-loss scheme: `0` off (HR diagnostics `KZ(12)` are then switched off too), `1`–`2` supernova schemes, `≥ 3` Eggleton–Tout–Hurley [Hurley2000](@cite), with binaries after [Hurley2002](@cite); must be ≥ 0 |
| `level` | String | `"C"` | SSE/BSE parameter level [Kamlah2022](@cite); one of `"A"`, `"B"`, `"C"`, `"0"` (no level: the engine's independent defaults) |
| `zmet` | Float | `0.001` | Metal abundance; `0.0001 ≤ zmet ≤ 0.03` (the engine's own bounds) |
| `epoch0` | Float | `0.0` | Formation time of the population [Myr]; must be ≤ 0 (the age at start is `−epoch0`) |
| `dtplot` | Float | `1.0` | Interval of the stellar-evolution diagnostics (`sev.83_*`) [NB]; must be > 0 and ≥ `deltat` |
| `dtplot_myr` | Float | `0.0` | The same in Myr (positive replaces `dtplot`); the ordering against `deltat` is checked at generation when either side is physical |

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

### Flat cluster form

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
| `period` | String | `"kroupa1995"` | `"kroupa1995"`: the [Kroupa1995](@cite) birth period distribution `f(log P) ∝ (log P − 1)/(45 + (log P − 1)²)`, `1 ≤ log(P/d) ≤ 8.43`, converted to a semi-major axis with Kepler's third law; `"loguniform"`: semi-major axis log-uniform in `[a_min, a_max]` |
| `a_min`, `a_max` | Float | `0.01`, `100.0` | Semi-major axis bounds [AU] for `"loguniform"`; `0 < a_min < a_max` |
| `q_min` | Float | `0.1` | Lower mass-ratio bound for `"uniform_q"`; `0 < q_min < 1` |
| `eccentricity` | String | `"thermal"` | `"thermal"` (`f(e) = 2e`) or `"circular"` |

The density sampler places *systems* (a pair by its centre of mass), which are virialised and truncated as such; each pair is then expanded into two bodies on a Keplerian orbit with a random phase and orientation. The engine reads primordial pairs as bodies `2i − 1, 2i` for `i ≤ NBIN0`, so the pairs of all clusters are written first, cluster by cluster, followed by every cluster's singles; `NBIN0` and `KZ(8) = 2` are set in `merger.inp`, and a cluster's members are then two contiguous blocks recorded in `merger_ic.toml` (`cluster_blocks`) and recoverable from `merger_summary.txt` (`binaries: N` per cluster). The summary also reports the hard fraction of each cluster's pairs (`G m₁ m₂ / 2a > ⟨m⟩ σ²` of the cluster's systems). Binary-rich runs carry a larger energy error per adjustment interval (regularised pairs, chains): the two-cluster demo with 20 % binaries reached 4×10⁻³ over its first time unit, so set `merger.nbody6.qe` accordingly (10⁻³ to 10⁻²) or the engine halts at the first adjustment that exceeds the default 2×10⁻⁴.

### Seed semantics

- `seed` **omitted** → a random seed is drawn at generation time (`0` is an ordinary seed). The IC is still reproducible after the fact: the effective seed is recorded in `merger_ic.toml` (`meta.seed`)
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

With stellar evolution active (the Level C defaults), BSE mass loss drives a slow symmetric adiabatic expansion of the triangle: all three separations grow together, by about a factor two over 20 NB. The pairwise distances must stay equal — asymmetric growth, or a shrinking triangle, indicates broken equilibrium velocities.

---

## The GPU-host and verification cases

### `verif_3d5cluster.toml`

Five clusters of 800 bodies placed well out of the z = 0 plane, which exercises the xz- and yz-projection panels and the three-dimensional COM tracking. The inward radial velocities are mild, so the clusters interact over `tcrit` without falling together completely.

### `gpu/merger_50k.toml`

Two equal King clusters of 25 000 stars on an eccentric orbit, the case of the hardware validation runs. The CUDA build is faster than the AVX build at every N measured, from N ≈ 2×10⁴ upwards, and the gain grows with N as the regular force scales as N²; the figures are in the manual, section Validated hardware.

### `gpu/gpu_pipeline.toml`

Clones and builds the engine with CUDA into its own tree — CPU and GPU objects differ — and runs `merger_50k.toml` on the first device. The host needs `nvcc` (CUDA ≥ 12.8 for the RTX 50 series, ≥ 11.8 for H100/H200), `nvidia-smi`, `gfortran`, `gcc`/`g++`, `make` and `git`.

### `gpu/cpu_pipeline.toml`

Clones and builds the AVX engine without CUDA into the default tree on the same host and runs the same merger as `gpu_pipeline.toml`, so the two runs and `bench/gpu_scaling.jl` are comparable binary against binary.

### `gpu/merger_600k.toml` and `gpu/single_600k.toml`

The probes above 5 × 10⁵ bodies. The merger is the case of `merger_50k.toml` at twelve times the membership: two King clusters of 300 000 stars, `W0 = 6`, half-mass radius 2 pc, Kroupa masses between 0.08 and 100 M☉, released at a 10 pc apocentre with `e = 0.5`; Jacobi truncation leaves 575 080 bodies, the combined half-mass radius is 5.3 pc and the N-body time unit 0.32 Myr. The single cluster has the same structure with 600 000 stars at rest (`n_clusters = 1`, explicit mode, no truncation; unit 0.07 Myr). Both are integrated for two time units with adjustments every 0.25 and outputs every 0.5, enough for the engine's cost per time unit to be measured well past start-up while the CPU reference stays within a few hours on a workstation. Both raise `virial_max_n` to 10⁶, since the generator otherwise refuses to virialise a cluster above 2 × 10⁵ stars; the pair sum takes about three minutes per case on twelve threads.

### `gpu/gpu_merger_600k.toml`, `gpu/cpu_merger_600k.toml`, `gpu/gpu_single_600k.toml`, `gpu/cpu_single_600k.toml`

The pipeline configurations of the two probes on the CUDA and the AVX binary: device 0, eight host threads, a start-up watchdog of one hour (the initial force and neighbour lists of 6 × 10⁵ bodies take minutes on the host), post-processing and figures on. Their install phase is off; they use the trees that `gpu_pipeline.toml` and `cpu_pipeline.toml` build. `postprocess.pair_sum_max_n` stays at its default, so the remnant and per-cluster diagnostics are skipped with a warning at this size and the binaries are compared on the engine time recorded in `RUN_INFO.toml`.

## Sources

The sources cited on this page are listed under [References](@ref); the King models are those of [King1966](@cite).
