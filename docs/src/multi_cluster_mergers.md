# Multi-Cluster Merger Simulations

Nbody6Dynamics ships a merger initial-condition generator (`src/ic/`) that produces `dat.10` particle files and matching `.inp` files for Nbody6++ external-IC runs (`KZ(22)=2`). It supports any number of clusters ≥ 2, King or Plummer density profiles per cluster, three IMF modes, automatic Kepler placement for cluster pairs, and Jacobi truncation.

## Entry points

```julia
# One-call: TOML → dat.10 + merger.inp + summary + metadata + diagnostic plots
result = run_merger_pipeline("input_files/merger_demo_small.toml")

# Config-driven: generate ICs, run the simulation, post-process, plot
# (config.toml: merger.enabled = true, merger.config_file = "input_files/...")
results = run_pipeline(load_config("config.toml"))

# Programmatic
cfg    = load_merger_config("input_files/merger_demo_small.toml")
result = generate_merger_ic(cfg; output_dir = "my_ics")

# Reload a previously generated IC from disk (no re-sampling)
result = load_merger_ic_result("runs/merger_run_.../output")
```

The TOML schema (flat and structured forms) is documented in [Input File Reference](@ref).

## Scientific context

Cluster merger simulations model the gravitational encounter and coalescence of star clusters, relevant to nuclear star cluster assembly, hierarchical formation of young massive clusters (R136, Westerlund 1), suspected merger remnants (NGC 1851, Terzan 5), and IMBH formation via runaway collisions during core mergers.

## Density profile samplers

Each cluster's spatial sampler is selected by its `DensityProfile` tag: `KingProfile(W0)` or `PlummerProfile()`.

### King (1966) [King1966](@cite) — `sample_king`

The King model is solved from its Poisson equation in standard dimensionless form, `r̂ = r/r₀` with King radius `r₀² = 9σ²/(4πGρ₀)`:

```math
\frac{d^2\hat{W}}{d\hat{r}^2} + \frac{2}{\hat{r}}\frac{d\hat{W}}{d\hat{r}}
    = -9\,\frac{\hat{\rho}(\hat{W})}{\hat{\rho}(W_0)},
\qquad
\hat{\rho}(\hat{W}) = e^{\hat{W}} \operatorname{erf}\!\big(\sqrt{\hat{W}}\big)
    - \sqrt{\tfrac{4\hat{W}}{\pi}} \left(1 + \tfrac{2\hat{W}}{3}\right)
```

The density on the right-hand side is **normalised to the central value** — with the unnormalised density the radial unit is compressed by `√ρ̂₀` and the concentration comes out wrong (a bug fixed in the current rewrite). The ODE is integrated with `Tsit5` from OrdinaryDiffEqTsit5 at `abstol = reltol = 1e-12`, starting slightly off-centre with the Taylor expansion `Ŵ ≈ W0 − (3/2) r̂²`, and terminated exactly at the tidal radius by a continuous callback on the `Ŵ = 0` crossing. The solution is evaluated on a log-spaced grid (dense in the core, resolved out to the edge for any concentration).

**Validation:** the concentration `c = log₁₀(r̂_t)` of the solved profile matches published King-model values (e.g. `c ≈ 1.25` for `W0 = 6`) to better than 1%; this is asserted in the test suite across several `W0`.

Particle positions are drawn by inverse-CDF sampling of `ρ(r̂) r̂²`; at each radius, speeds are rejection-sampled from the lowered Maxwellian `f(v) ∝ v² [e^{Ŵ − v²/2} − 1]` for `v < v_esc = √(2Ŵ)`.

### Plummer (1911) [Plummer1911](@cite) — `sample_plummer`

Inverse-CDF sampling for radius (`r = a/√(X^{-2/3} − 1)`) and von Neumann rejection for the velocity distribution `f(q) ∝ q²(1−q²)^{7/2}` [AarsethHenonWielen1974](@cite). Useful relations: half-mass radius `r_hm ≈ 1.305 a`, virial radius `r_v ≈ 1.70 a`. The generator converts the requested half-mass radius to the scale radius via `a = rbar/1.305`.

## IMF sampling

Masses are drawn from the [Kroupa2001](@cite) continuous broken power law (`ξ(m) ∝ m^{-α}` with `α = 0.3, 1.3, 2.3` across breaks at 0.01, 0.08, 0.5 M☉) by exact inverse-CDF sampling over the segments overlapping the requested `[bodyn, body1]` window — O(N), no rejection. `kroupa_mean_mass(m_low, m_up)` gives the analytic mean (`≈ 0.58` M☉ for [0.08, 100]).

Three modes, selected by the `IMFSpec` tag on each cluster:

- **`KroupaIMF`** — natural sampling in `[bodyn, body1]`. The total cluster mass is an *output* (the sum of the N samples), not a parameter. Use this when the scientific intent is a realistic stellar population: Nbody6++'s stellar-evolution prescriptions are well-defined here.
- **`RescaledKroupaIMF`** — "super-particle" mode: sample from Kroupa, then apply a uniform linear rescale so the sum equals `target_mass`. The IMF *shape* is preserved but the effective mass range shifts by the rescale factor. A warning is emitted at sample time when the factor falls outside ×[0.7, 1.4]: in that regime individual body masses no longer correspond to real stars and downstream stellar-evolution output (SEV/BEV files, HR diagrams) is non-physical. The flat form `imf = "kroupa"` with `mass_total` selects it as well; `imf = "kroupa_rescaled"` states the intent.
- **`EqualMassIMF`** — every body gets the same `particle_mass`.

`expected_mass(imf, N)` returns the expected total per cluster — exact for the rescaled/equal modes, the analytic expectation `N⟨m⟩` for natural Kroupa.

## Per-cluster scaling and virialisation

For each cluster the generator:

1. Samples `N` masses from the IMF and positions/velocities from the profile
2. Rescales positions so the **mass-weighted half-mass radius** equals the target `rbar` (not the count-median radius — with an IMF the two differ by sampling noise, and only the mass-based definition matches RBAR semantics)
3. Calls `virialise!`: shifts to the centre-of-mass frame and scales velocities so `Q = T/|W| = 0.5` exactly, using the exact O(N²) pairwise potential (threaded over strided rows). A guard refuses above `nmax = 200_000` particles per cluster unless raised explicitly.

Two independently virialised clusters are *not* in virial equilibrium as a combined system (the mutual potential is unaccounted), which is precisely why the output targets `KZ(22)=2`: that path sets `LSCALE=.FALSE.` in Nbody6++, so no centre-of-mass correction or velocity rescaling is applied to the combined ICs.

## Primordial binaries

`[merger.clusterN.binaries]` gives a cluster a primordial binary population: a fraction of its stars are paired (random pairing from the IMF or a uniform mass-ratio distribution), periods follow the [Kroupa1995](@cite) birth distribution or a log-uniform semi-major axis, eccentricities are thermal or circular. Pairs travel through sampling, virialisation, orbit placement, and Jacobi truncation as single systems and are expanded into two bodies on Keplerian orbits at the end, so the cluster's structure and virial state are those of the systems. The engine's convention puts every primordial pair at the head of `dat.10`; membership is therefore recorded as index blocks (`cluster_blocks`, `parse_merger_summary`) and every per-cluster diagnostic accepts body-index vectors. The hard fraction of each cluster's pairs is reported at generation.

The engine reports its binaries in `bev.82_<t>` (one file per `dtplot` interval, KS-regularised pairs only). `read_binary_evolution` parses them, `binary_population` reduces a run to pair counts, the binary fraction and the hard/soft split, and the pipeline draws the population against time, the initial and final `a`–`e` diagrams by hardness class, and the log-period histograms. The hard/soft scale `⟨m⟩ σ²` is measured on the systems of the nearest snapshot exactly as the generator defines it, so the initial hard fraction of the *regularised* pairs can be compared with the generated one: in the two-cluster demonstration with 20 % binaries per cluster, 58–60 % of the 324 generated pairs are hard, and 224 of them are regularised at `t = 0`, of which 77 % are hard, because the engine regularises the tight pairs preferentially while the wide, soft ones remain unresolved by these diagnostics. Counts from `bev.82` are therefore lower bounds on the bound-pair population.

## Orbit placement

### Kepler mode (2 clusters)

`setup_two_cluster_orbit` places the pair at apocentre on the x-axis. With semi-major axis `a = d_apo/(1+e)`, the vis-viva relation at `r = d_apo` gives the purely tangential apocentre speed

```math
v_{\rm apo} = \sqrt{\frac{G\,(M_1+M_2)\,(1-e)}{a\,(1+e)}},
```

decomposed into centre-of-mass frame speeds `v₁ = (M₂/M) v_apo`, `v₂ = (M₁/M) v_apo` (`kepler_velocity`). Cluster 1 sits at `(−d₁, 0, 0)` with velocity `(0, +v₁, 0)` and cluster 2 at `(+d₂, 0, 0)` with `(0, −v₂, 0)`, where `d₁ = (M₂/M) d_apo`, `d₂ = (M₁/M) d_apo`.

### Explicit mode (N ≥ 2 clusters)

Each cluster specifies its own COM `position` [pc] and `velocity` [km s⁻¹]. `combine_clusters_explicit` converts the velocities to the generator's code units, applies the offsets, and shifts the combined system to its centre-of-mass frame.

## Jacobi truncation

With `output.truncate_jacobi = true` (default), each cluster is truncated **before combining** at its instantaneous Jacobi radius with respect to its nearest neighbour at separation `d`:

```math
r_J = d \left(\frac{M_{\rm self}}{3\,M_{\rm other}}\right)^{1/3}
```

(`jacobi_radius`). Stars beyond `r_J` would be immediately unbound at the initial separation; removing them avoids transient mass loss and energy errors at simulation start. The number removed per cluster is logged, and the post-truncation counts are recorded in `merger_summary.txt` and `merger_ic.toml`.

## Units and output files

Sampling and orbit placement happen in internal code units: `G = 1` with masses in M☉ and lengths in pc, so one code velocity unit is `√(G M☉/pc) = 0.06558` km s⁻¹ (with `G = 4.30091×10⁻³ pc (km s⁻¹)² M☉⁻¹`) and the code time unit is 14.91 Myr. User-facing quantities (positions, velocities, masses, times of the summary) are physical.

`generate_merger_ic` writes four files (all protected by the never-overwrite `name#k.ext` backup policy):

- **`dat.10`** — one line per particle, `MASS X Y Z VX VY VZ` at full precision. The arrays are converted with `to_nbody_units!` to Hénon units (`G = 1`, `M_total = 1`).
- **`merger.inp`** — matching NAMELIST input file: `N` = post-truncation total, `KZ(22) = 2` (N-body-unit `dat.10` input, the only supported mode), `KZ(14) = 0` (isolated), `NRAND` = effective seed, `TCRIT`/`DTADJ`/`DELTAT` from `[merger.output]`, every `&INNBODY6`/`&ININPUT` numerical field from `[merger.nbody6]` (integration parameters derived from the smallest member cluster unless set, run-time limits, block thresholds, output multipliers, explicit `KZ` overrides) and the stellar-evolution settings from `[merger.stellar]` (see [Input File Reference](@ref)), and:
  - **`RBAR` = the combined system's mass-weighted half-mass radius [pc]** — this is the NB length unit used for the `dat.10` conversion, so the physical scaling in Nbody6++ is self-consistent
  - **`ZMBAR` = mean particle mass `M_total/N_total` [M☉]**
- **`merger_summary.txt`** — human-readable summary: per-cluster profile/IMF/N (before and after truncation)/mass/`r_hm`/`W0`/COM state, orbit parameters, combined totals, output format. Parsed later by `parse_merger_summary` to recover the per-cluster particle index ranges.
- **`merger_ic.toml`** — machine-readable metadata (schema v2): generation timestamp, effective seed, `external_rng` flag, combined totals (`N_total`, `M_total`, `rbar`, `zmbar`), orbit mode and parameters, output spec, cluster index ranges, and one structured table per cluster (`profile = {type = ...}`, `imf = {type = ...}`, position/velocity, `N_after_trunc`). The cluster tables use the same structured schema the config parser accepts, so the file round-trips.

The returned `MergerICResult` carries everything downstream plotting needs without re-reading files: output dir, totals, cluster ranges and specs, orbit info, and **physical-unit** copies of the particle arrays (M☉, pc, km/s).

## Reloading ICs — `load_merger_ic_result`

`load_merger_ic_result(dir)` reconstructs a `MergerICResult` from `dat.10` + `merger_ic.toml` without re-sampling — e.g. to regenerate `plot_merger_ic` output after a plotting fix without touching the particle data. NB-unit files are converted back to physical units using the stored `rbar` and mass scale (`vstar = 0.06558 √(M_total/rbar)` km s⁻¹).

## Running and diagnosing merger simulations

When `run_pipeline` runs with `merger.enabled = true` and `simulation.run_test = true`, the binary executes *inside* the merger output directory so `dat.10` is found in its working directory, using `merger.inp` as input.

Post-processing then adds two merger-specific plots whenever `merger_summary.txt` sits next to the snapshots (this also works in `postprocess_external`):

- **`plot_cluster_separation`** — pairwise COM separations of the initial clusters over time (tracked by particle NAME ranges). For > 5 clusters: min/max envelope + mean, a union-find staircase counting spatially distinct surviving clusters, and a heuristic coalescence-time marker.
- **`plot_cluster_virial`** — internal virial ratio `Q_i(t)` of each initial cluster from its bound members (COM-velocity subtracted, self-gravity only; `per_cluster_virial` for the raw matrix, `bound_only = false` for all present members). Meaningful *before* coalescence; after merging, the bound subset of an ID-group is whatever remains self-bound in the remnant.
- **`plot_density_profiles`** — density profile of each initial cluster about its own centre from its bound members (`cluster_profiles`, log-spaced shells), with the generating King or Plummer model (`model_density`, evaluated at the cluster's initial half-mass radius and its bound mass) dashed and a ratio strip `ρ/ρ_model` beneath. Drawn for the first snapshot (a check of the sampling and virialisation) and the last; `system_profile` gives the remnant's profile about the global shrinking-sphere centre.
- **`plot_velocity_dispersion`** — radial and one-dimensional tangential velocity dispersions per cluster and the anisotropy `β(r) = 1 − σ_t²/σ_r²` (0 isotropic), from the same shells.
- **`plot_binary_population`**, **`plot_binary_orbital_elements`**, **`plot_binary_period_distribution`** — drawn whenever `bev.82_*` files are present: regularised pair counts with the hard/soft split above the binary and hard fractions; `a`–`e` diagrams at the first and last epoch, classed by the Heggie criterion [Heggie1975](@cite) with the boundary of a pair of mean component-mass product marked; and the `log10(P/d)` histograms of the first against the last epoch (`binary_population`, `binary_hardness`, `hardness_scale`).
- **Remnant diagnostics** (`remnant_diagnostics`, written to `remnant_diagnostics.csv` in the run directory) — at every snapshot the self-consistent bound set of the whole system is analysed: the Casertano–Hut core radius [CasertanoHut1985](@cite) (from the engine's local densities, or a sixth-neighbour estimate) and half-mass radius (`plot_remnant_structure`, with `r_h/r_c`); the rotation about the spin axis — the ordered-motion parameter `λ_R` [Emsellem2007](@cite), the Peebles spin `λ_P` [Peebles1969](@cite), and the alignment of the spin axis with the orbital angular momentum of the initial clusters (`plot_remnant_rotation`), plus the `v_rot/σ` profile in cylindrical shells at the last snapshot (`plot_rotation_profile`); the [Allison2009](@cite) mass-segregation ratio `Λ_MSR` with its error and the half-mass ratio of the most massive stars (`plot_mass_segregation_evolution`, with a configurable threshold defining the segregation time); and the coalescence time, the first snapshot at which every initial cluster overlaps its neighbours in the union-find sense (`coalescence_time`, also marked on the separation figure). Before coalescence `λ_R` measures the orbital motion of the pair, so the rotation parameters are meaningful for the remnant only after `t_coalesce`. Cost O(N²) per snapshot.
- **`plot_cluster_structure`** — two stacked panels: the half-mass radius of each initial cluster about its own centre from its bound members, with the engine's global 50 % Lagrangian radius overlaid for contrast, and the bound mass fraction. `cluster_structure` returns the underlying `ClusterStructure` (members, bound members, bound mass fraction, shrinking-sphere centre, 10/50/90 % radii, one-dimensional velocity dispersion, virial ratio; all per cluster and snapshot).

`plot_merger_ic` documents the ICs themselves: spatial projections and a 3-panel overview (viridis mass colouring, orbit annotation), a velocity quiver coloured by cluster membership, the sampled IMF histogram against the Kroupa reference slopes (`α = 1.3`, `2.3`), and per-cluster radial density profiles.

## Parameter sweeps

A sweep TOML (`input_files/sweep_demo.toml`) names a base pipeline config and a base merger config, the grid axes as dotted keys into the merger TOML (`"merger.orbit.eccentricity" = [0.0, 0.6]`, `"merger.cluster2.N" = [500, 1000]`, Cartesian product; the parent table must exist in the base file) and the seeds replicated at every grid point. `scripts/run_sweep.jl` (or `run_sweep`) writes one directory per point with the derived `merger.toml` and `config.toml`, runs the points as concurrent worker processes with the configured `omp_threads` (four threads per job and five jobs at a time follow the cost model of the reference workstation), keeps `sweep_index.toml` current, writes `sweep_summary.csv` (final time, star and pair counts, energy error, virial ratio per point) and draws the comparison figures: the half-mass Lagrangian radius and the energy error of every run on common axes, coloured by the value of one grid axis with seeds sharing the colour. `--dry-run` prepares everything without launching. With several seeds the completed points of each grid point form an ensemble: `plot_sweep_ensemble` draws the median and the central 68 % and 95 % bands of the Lagrangian radius, the energy error, or the star and pair counts on a common time grid, and a sweep without grid axes is the plain seed ensemble of one configuration. The mass ratio, the orbital eccentricity, the apocentre, the tidal field (`merger.tidal.kz14`, with a `[merger.tidal]` table in the base file) and the binary fraction are all reachable this way; see the manual for the file format.

## Control runs

Every merger case has an isolated equivalent: the same total number of stars, structural model, IMF, binary population and integration, stellar-evolution and tidal settings in one cluster at rest, with the `N`-weighted mean of the progenitors' half-mass radii as its half-mass radius. `write_control_merger_config` derives it from any merger TOML at the file level (no sampling involved). A control's N-body time unit is shorter than the merger's (one compact cluster against a wide configuration), so intervals given in N-body units are scaled by the ratio of the estimated time units, `(RBAR_est / r_h,control)^{3/2}`; intervals given in Myr (`tcrit_myr`, `dtadj_myr`, `deltat_myr`, `dtplot_myr`, converted by the generator with each run's realised `T*`) match the merger's physical span exactly and are the recommended form for controlled comparisons. `controls = true` in a sweep adds the control as a companion point of every grid point and seed, drawn against its merger by `plot_control_comparison`. The control isolates the effect of the merger geometry on the quantities the pipeline tracks (energy error, Lagrangian radii, counts, binary population, remnant structure); since the engine evaluates no tidal energy, a tidal merger's accuracy is judged against its control rather than against the energy tolerance.

## Feasibility and limitations

Nbody6++GPU integrates a multi-cluster system correctly as a set of point masses, but every diagnostic and adaptive control inside the engine assumes one centrally concentrated cluster. Knowing which numbers to trust before the clusters coalesce decides whether a merger run produces science or artefacts. The statements below were checked against the engine source (upstream v2026.07) and confirmed with two small runs of the shipped configurations (`merger_demo_small.toml`, two clusters of N = 1000 on a Kepler orbit; `merger_3cluster_small.toml`, three clusters of N = 600–800 on explicit orbits).

**What the engine does with the initial conditions.** The `KZ(22) = 2` route disables the engine's own centring, mass normalisation, and virial rescaling (`scale.F`, `LSCALE = .FALSE.`), so the particles are integrated exactly as generated; `virialise!` and the orbit setup fix the dynamical state, and the `Q = 0.5` line in `&INSCALE` is inert. The length unit is the half-mass radius of the *combined* configuration, so a member cluster of `r_hm = 2 pc` in a system whose `RBAR` is 6–12 pc occupies only 0.1–0.3 N-body length units.

**One centre, one scale radius.** The engine locates a single density centre `RDENS` (`core.f`, a density-weighted mean inside one sphere that is inflated until it holds N/5 particles) and one scale radius `RSCALE = 0.5 M²/|W|` from the global potential. Before coalescence both describe the configuration, not its members. In the two-cluster run at t = 0 the engine's core radius was 0.80 length units (5 pc) against member half-mass radii of 0.31 (2 pc); its density centre sat 0.50 and 1.34 units from the two cluster centres, which were 1.84 units apart, and drifted between them as the orbit evolved. Consequently:

- `lagr.7` measures the mass distribution of the configuration about `RDENS` (its 50 % radius was 1.05 units, 6.7 pc, at t = 0), then that of the merged remnant. It never describes a member cluster.
- The ADJUST line's `RSCALE`, `RC`, and `RTIDE`, and the escape angles and radii in `esc.11`, are all taken about the same centre.
- The regular time-step control, the crossing-time counter, and the neighbour-number stabilisation (`regcor_gpu.F` reduces the target membership outside the global half-mass radius) treat every member cluster as the outskirts of one big system.
- The Ahmad–Cohen neighbour scheme [AhmadCohen1973](@cite): the initial neighbour radii are set from the coordinate origin as `RS0 √(1 + r²)` (`fpoly0.F`). A fixed `RS0 = 0.5`, several member half-mass radii wide, makes every neighbour list overflow at start-up, and the engine's repair loop then shrinks them: in such runs `RSMIN` falls from 0.50 to 0.03 within the first adjustment interval. `RS0`, `RMIN`, and `DTMIN` are therefore derived from the smallest member cluster by default (`[merger.nbody6]`), the generator refuses an `RS0` wider than a member, and `KZ(16) = 0` keeps the engine from re-deriving them from its global quantities.

**Escapers.** For an isolated system (`KZ(14) = 0`) the engine sets `RTIDE = 10 RSCALE` at every adjustment and removes particles beyond `2 RTIDE` from `RDENS` with positive energy (`adjust.F`, `escape.F`); the `RTIDE` given in the input is overwritten, so the escape radius cannot be chosen by the user. For the two-cluster run this radius was 13 length units (84 pc). Stars ejected early therefore stay in the system for many crossing times and distort the global virial ratio and energy bookkeeping until they cross that sphere. In a 27-cluster grid run this delay produced 1442 escapers in a single step, which overflowed the fixed item list of the summary `WRITE` in `escape.F` and terminated the run with a Fortran runtime error (exit status 2); any large ejection event can trigger the same failure. In an external field (`[merger.tidal]`: the solar-neighbourhood tide, a point-mass galaxy, or MWPotential2014) the engine removes escapers on the distance criterion alone, and its energy bookkeeping omits the tidal potential energy altogether (the evaluation is commented out in `energy.F`), so the relative energy change it checks against `QE` is the work done by the field: on the two-cluster demo it scales with the tidal strength as Ω², from 6×10⁻⁵ per adjustment interval at 30 kpc to 2×10⁻³ at 8.5 kpc around a 10¹¹ M☉ point mass, against 2×10⁻⁶ isolated. The generator therefore refuses a tidal configuration with `qe < 0.01`, and the integration accuracy of a tidal run has to be judged against its isolated control. `bound_fraction(snap)` gives a snapshot-based escaper measure independent of the engine's escape sphere: the mass fraction bound to the whole system in its own frame.

**What to trust instead.** Cluster-level results before coalescence must come from the `conf.3` snapshots by original membership. `merger_ic.toml` records the particle index range of every cluster (`cluster_ranges`); `cluster_structure` measures every cluster about its own shrinking-sphere centre from the members it still binds (centre, 10/50/90 % radii, velocity dispersion, bound mass fraction, virial ratio), `plot_cluster_structure` and `plot_cluster_virial` draw these, and `plot_cluster_separation` tracks the centres of mass. Per-cluster centres and half-mass radii computed this way agree with the generated clusters at t = 0 and follow their tidal growth through pericentre, where the engine's global numbers do not. The bound-member selection (a self-consistent negative-energy criterion in each cluster's own frame) removes stripped stars and kicked stellar remnants, which otherwise dominate the velocity dispersion and the virial ratio after the first supernovae; once a cluster has dissolved into the remnant, its bound subset is whatever remains self-bound and the diagnostics stop describing a cluster.

**Regime.** The engine has no native multi-cluster mode: `KZ(5) = 2` (two orbiting Plummer models) is dead code in this version, since it reads list-directed input that the NAMELIST format no longer supplies. Two clusters on a Kepler orbit and three to five clusters on explicit orbits are the configurations the package handles well. A large grid such as the 27-cluster example is a cold-collapse experiment of a cluster of clusters: its combined virial ratio at t = 0 is 0.13, every member is disrupted within about two crossing times of the configuration (about ten internal crossing times), and none of the engine's diagnostics is meaningful before the remnant forms. The generator computes the combined virial ratio `Q = T/|W|` and the ratio of the length unit to the smallest member half-mass radius, prints both in `merger_summary.txt`, stores them in `merger_ic.toml`, and warns below `Q = 0.3` and above `RBAR/r_hm,min = 5`.

**Integration accuracy.** The generated input uses `QE = 2×10⁻⁴` by default (`[merger.nbody6]`); the shipped single-cluster inputs, which evolve stars in a tidal field, use `1×10⁻²` (see the input-file reference). The small runs without binaries conserved energy to better than 10⁻⁴ per adjustment interval; with 20 % primordial binaries the same demo reached 4×10⁻³ over its first time unit, so binary-rich configurations need `qe` of order 10⁻³ to 10⁻².

**Stellar masses.** `mass_total` together with `N` fixes the mean stellar mass. Earlier versions of the shipped configurations requested 10⁴ M☉ for 1000 stars, a mean of 10 M☉ that the rescaled Kroupa sampler reached by scaling every mass up by a factor of about 17, so they contained stars of several hundred solar masses: supernova kicks began after 3 Myr, mass loss took a quarter of the CPU time, and the membership-based virial ratios blew up when the kicked remnants left. The shipped files now sample the natural Kroupa IMF (explicit-orbit velocities were rescaled by the square root of the mass ratio, which preserves the virial state), and the loader refuses a rescale factor outside ×[0.5, 2] and warns outside ×[0.7, 1.4]. Science configurations need `N` and `mass_total` consistent with the IMF mean (about 0.58 M☉ for Kroupa over 0.08–100 M☉), or no `mass_total` at all.

**Execution.** The engine takes its thread count from the OpenMP runtime alone and reports it at start-up; `simulation.omp_threads` controls it, and `RUN_INFO.toml` records the configured and reported values with the CPU accounting of the run. GPU devices are selected through the engine's `GPU_LIST` environment variable, not `CUDA_VISIBLE_DEVICES`.

## Showcase configurations

`input_files/showcase/` holds four cases that exercise the pipeline end to end within minutes on a workstation, each as a merger TOML with a matching pipeline TOML (paths relative to the folder): an equal-mass King pair on an eccentric orbit that coalesces within the run (`equal_*`), an unequal pair with 20 % primordial binaries in both clusters and the relaxed energy tolerance such runs need (`binary_*`), the equal pair inside a point-mass galactic potential (`tidal_*`), and an eccentricity × seed sweep with isolated controls (`sweep.toml`). All intervals are given in Myr, so merger, control and sweep points cover the same physical span.

## Verification configs

Two configs in `input_files/` exercise the machinery end-to-end: `verif_triorbit.toml` (three clusters on a rotating Lagrange-equilibrium triangle, `ω² = Gm/(√3 r³)` — see [Input File Reference](@ref) for the derivation) and `verif_3d5cluster.toml` (five clusters distributed in 3D). The test suite additionally validates the King concentration `c(W0)`, the Plummer `r_hm = 1.305 a` relation, the Kroupa mean mass, `virialise!` reaching `Q = 0.5`, and the Kepler/Jacobi relations.

## Related work

Cluster assembly through mergers of sub-clusters is the subject of [Fujii2012](@cite). The initial-condition conventions follow the McLuster generator [Kupper2011](@cite) and the engine's own documentation [Aarseth2003](@cite); the direct N-body method and its regularisations are set out in [HeggieHut2003](@cite). Every source cited in these pages is listed under [References](@ref).
- Wang, L. et al. (2015). MNRAS 450, 4070. — Nbody6++GPU code paper.
