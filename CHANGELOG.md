# Changelog

Notable changes to Nbody6Dynamics.jl, by version. The package follows
semantic versioning; before 1.0 a minor version may change the API, and the
entries say so.

## [Unreleased]

### Changed

- Engine pin `v2026.07.10`: the Peters decay of compact KS pairs is accumulated
  at every check from the elements of that interval and applied as sums
  (`GWACCUM`; an encounter that changes e no longer charges the interval before
  it at the new rate: on a planted fly-by raising e from 0.50 to 0.66 the first
  application removed 3.85 times the correct Δa with the previous tag and 1.00
  with this one), and the block-step list stays consistent at an escaper's
  removal (a duplicated index after the removal of a freshly kicked remnant
  corrected one particle twice with a zero step and spun the run in `CORE`
  without output, the stall of the showcase sweep point of 2026-10-23).
  Restart dumps written by earlier tags cannot be restarted with it.
- Source layout: the orchestration modules (platform, install, telemetry,
  run, resume, cost, cluster, GPU validation, sweeps, external
  post-processing) live under `src/pipeline/` and the diagnostics derived
  from the output (potential and integrals, budgets, structure, binaries,
  remnant, validity, gravitational-wave observables, ensembles) under
  `src/analysis/`; the core types, configuration, utilities and the figure
  interface stay at `src/`. No name or interface changed.

### Added

- Rotating King profile of the generator: `profile = { type = "rotating_king",
  W0, omega0, axis }` samples the Einsel & Spurzem (1999) model
  (`sample_rotating_king`, solved once per `(W0, omega0)`), rotates its
  angular momentum onto `axis` (`[0, 0, 1]` prograde under `KZ(14) = 1, 2`),
  scales and virialises it like the King profile, writes the model's and the
  realised λ and specific angular momentum to `merger_ic.toml`
  (`clusterN.rotation`, `rotation_moments`) and refuses primordial mass
  segregation; `model_density` gives its angle-averaged profile.
- `[merger.nbody6] gw_recoil` (`"off"`, `"fit_q"`, `"nr_fits"`,
  `"test_1e4"` → `&INBSE kicktype`), `bh_spin` (`"zero"` … `"mesa"` →
  `bhspinfl`, the unseeded draws warned) and `escapers` (`"remove"` or
  `"keep"` → `KZ(23)` 2 or 0); any recoil needs `stellar.kz19 ≥ 3`.
- Gravitational-wave catalogues (`[postprocess] gw_catalogues`, default `false`):
  `compact_binary_catalogue` (one row per Peters-application line, NS/BH BINARY
  line or `bev.82` record of a pair of white dwarfs, neutron stars or black
  holes: masses, chirp mass, elements, GW frequencies, Peters, Mandel and
  circular coalescence times, eccentricities and times at 10 mHz and 10 Hz,
  band residence times, channel from the covering span),
  `escaped_binary_catalogue` (one row per `BINARY ESCAPE` of a compact pair,
  with the escape mode, the Peters evolution and the state at an observation
  epoch) and `gw_merger_catalogue` (one row per compact coalescence of
  `coal.24` joined with the last decay line and the recoil line: trigger
  elements, kick, retention, generation). `band_residence` and
  `eccentricity_after` evolve a state along the Peters track
  (`peters_semi_major_axis`); keys `gw_band_hz`, `gw_t_now_myr`.
- Binary pair types, hardness and element distributions:
  `binary_pair_type_series` (`binary_pair_types.csv`: MS–MS, MS–evolved,
  MS–WD, WD–WD, NS–any, BH–BH, BH–star, other and untyped pairs, the hard/soft
  split against the ⟨m⟩σ² of the systems inside the half-mass radius, median
  a, e and q per shell and census level) and `binary_distribution_series`
  (`binary_distributions.csv`: histograms of log a, e, q and log |E_b|/⟨m⟩σ²
  per snapshot and level); stages `binary_pair_types`, `binary_distributions`.
  The population frames gain `"projected_rh"` and `"milone"` (the three
  projected regions of Milone et al. 2012 with the core radius from the
  projected number-density profile), a `population_projection_axis` and a
  tidal filter `population_r_tide_nb`.
- `[merger.clusterN.binaries.massive] model = "moe2017"`: the joint
  period–mass-ratio–eccentricity distribution of Moe & Di Stefano (2017) with
  the inner-binary sampler of Moe's code (period table per primary mass, twin
  excess, small- and large-q power laws, eccentricity law with the turnover
  of `ecc_turnover`), the primary mass clamped to 0.8–40 M☉ with a count of the
  clamped primaries; `moe2017_integrals` reproduces the published companion
  frequencies and the 0.41/0.96 inner-binary fractions.
- `PowerLawIMF` (`imf = { type = "powerlaw", alpha = 2.35, bounds = [lo, hi] }`):
  a single power-law mass function sampled by the inverse cumulative
  distribution, `powerlaw_mean_mass`; the Salpeter populations of the
  published two-cluster replicas.
- `solve_rotating_king`, `sample_rotating_king`, `rotating_king_moments`,
  `RotatingKingModel`: the rotating King model of Einsel & Spurzem (1999),
  `f ∝ (e^{−βE} − 1) e^{βΩ₀J_z}`, solved self-consistently by a Legendre
  expansion and under-relaxed iteration, with exact velocity sampling at every
  position (speed by rejection against the King law, azimuthal angle by the
  closed-form inverse); model units `r_K = σ = 1`.
- Replica configurations of published two-cluster mergers under
  `input_files/verification/replicas/`: Priyatikanto et al. (2016) models
  notide, draco1, draco2 and draco5 (2 × 2048, Plummer, Salpeter 0.15–50 M☉,
  10 pc circular orbit in the solar-neighbourhood tide) and de la Fuente
  Marcos & de la Fuente Marcos (2010) q = 1 pairs (2 × 4096, S₀ = 10 pc,
  e₀ = 0 and 0.6) with their single-cluster control of 8192 stars.
- `simulation.stall_timeout` (default `0`, off): a stall watchdog armed after
  the first adjustment beyond t = 0. When neither the stdout capture grows nor
  a new restart dump appears for that wall-clock time, the engine is asked to
  stop, terminated if it has not stopped the same time later, and the segment
  is recorded as `stalled`; `resume_run` refuses such a run and a cluster job
  chain ends on it. An engine whose block time stalls on a collapsed step
  otherwise runs without output until its wall budget is spent.
- Readers of the engine's binary, few-body and gravitational-wave output.
  `read_binary_events(path)` reads `out1000` in one pass into a
  `BinaryEventLog`: KS episodes (`KSEvent`: NEW/END KSREG, NS/BH BINARY),
  chain, triple and quad regularisations (`FewBodyEvent`), stable hierarchies
  (`HierarchyEvent`), exchanges, the termination diagnostics of the fork
  (`TerminationEvent`), Peters decay steps, GW mergers and recoil, BSE
  braking (`GWEvent`), escapes of binaries and hierarchies
  (`BinaryEscapeEvent`), mass-loss terminations and kicks (`MassLossEvent`),
  tidal circularisation and capture (`TidalEvent`) and Roche-lobe overflow
  (`RocheEvent`), each record with its line, adjustment interval and time.
  File readers: `read_ks_pair_listing` (`bdat.9`), `read_wide_pair_listing`
  (`bwdat.19`), `read_hierarchy_listing` (`merger.84`), `read_hierarchy_log`
  (`hiarch.12`), `read_collision_log` (`coll.13`), `read_coalescence_log`
  (`coal.24`), `read_cmbody_log` (`cmbody.26`), `read_binary_stage_log`
  (`binev.17`), `read_degenerate_log` (`degen.4`), `read_escaped_binary_log`
  (`escbin.31`), `read_hardest_binary_log` (`hbin.39`). Formatted records are
  read by their FORMAT's field widths (`read_fixed`), since the engine glues
  fields together; values the engine prints under a wrong label or unit are
  stored under their true meaning, as each record type documents. The
  `GR CAPTURE` and `GR FLYBY` lines of hyperbolic NS/BH passages (engine fork
  v2026.07.8) are `GWEvent` kinds `:gw_capture` and `:gw_flyby`.
- Physically bound pairs and their spans from snapshots (`bound_pairs`, `bound_spans`, `BoundPair`,
  `BoundSpan`, `write_bound_spans`): mutual most-bound pairs whose apocentre perturbation by the nearest bodies is
  below 0.1, or whose dominant perturber is the bound outer body of a Mardling & Aarseth (2001) stable triple,
  plus every KS pair of the snapshot; spans over the snapshot sequence bridge single misses of a pair that stays
  bound, with interval-censored formation and end times. The stage `[postprocess] bound_spans` (default `false`; `bound_gamma_max`, `bound_n_miss`) writes
  `bound_spans.csv`.
- Gravitational-wave observables of a binary: `gw_inspiral_time` (Peters 1964 orbit average by singularity-free
  quadratures of the time factor, `gw_time_factor`, exact to 10⁻⁶ against the reference table and to 768/425 as
  e → 1, checked against an integration of the Peters equations), `gw_inspiral_time_mandel` (Mandel 2021 fit),
  `gw_zwick_factor` (the 1PN-recalibrated factor of Zwick et al. 2020, reported separately), `gw_enhancement`,
  `chirp_mass`, `orbital_frequency`, `gw_frequency` (f₂₂, Wen 2003 peak, pericentre) and
  `eccentricity_at_frequency` (the eccentricity and time at a reference frequency along the Peters track).
  QuadGK is a new direct dependency (it was already in the dependency tree).
- Binary census (`binary_census`, `CensusSpan`, `write_binary_census`): the union per pair of the overlapping or
  adjacent engine-regularised and bound spans, with both levels' limits and counts, the origin of the earlier and
  the end class of the later level; written as `binary_census.csv` when both span stages run.
- Binary population series (`binary_population_series`, `BinaryPopulationRow`, `write_binary_population`,
  `stellar_type_map`; `[postprocess] binary_population`, off by default, writes `binary_population.csv`): at every
  snapshot and for both census levels (KS pairs; bound pairs), the binary fractions by systems, stars and mass, their
  binomial error and a q-restricted fraction in a primary-mass window, over all radii and per Lagrangian or r/r_h
  shell about the shrinking-sphere centre, for all stars and per stellar class (binaries by the class of the
  primary, star and mass fractions per object) with the types of `sev.83`/`bev.82`.
- Two-cluster orbits (`orbit_mode = "kepler"`) by pericentre or apocentre, absolute or in units of the
  half-mass radius of cluster 1 (`[merger.orbit]` `pericentre`, `pericentre_rh`, `apocentre_rh`; one size key
  per config), oriented by `orbit_normal` (the direction of the orbital angular momentum; default [0, 0, −1],
  the earlier layout) and `apocentre_angle`. `[merger.output] truncate_radius_rh` truncates every cluster at a
  fixed multiple of its own half-mass radius instead of the orbit-dependent Jacobi radius; the cluster key
  `density_ratio` derives the half-mass radius of cluster N ≥ 2 from its mean density relative to cluster 1;
  `[merger] seed_streams = "per_cluster"` gives every cluster its own random stream, so progenitors are
  identical across an orbit grid. `merger_ic.toml` (schema 4) records the resolved apocentre, a, r_p, the
  period, orbital energy and angular momentum, j_orb, the Peebles spin parameter, v_p/σ₁, the Roche-lobe
  filling at apocentre and the truncated mass per cluster (`two_cluster_orbit_quantities`). Existing
  configurations generate the same `dat.10` and `.inp` as before.
- Cluster complexes (`orbit_mode = "complex"`, `[merger.complex]` with a member template
  `[merger.complex.member]`, `ComplexSpec`): member masses equal or from dN/dM ∝ M^−β, half-mass radii from a
  constant or the Marks & Kroupa (2012), FROST-CLUSTERS (Rantala et al. 2024) or Dobbs et al. (2022) size–mass
  relation with optional log-normal scatter, centres in a uniform sphere, a truncated Plummer distribution with
  Aarseth–Hénon–Wielen velocities, a Goodwin & Whitworth (2004) box fractal with inherited velocities or a
  filament, and the bulk velocities scaled to a target virial ratio of the members as softened point masses, with
  an optional solid-body rotation about +z. Every member has its own random stream; the record holds the realised
  Q, E, L and Peebles λ of the centres. The derived control of a merger configuration (`control_merger_dict`)
  resolves the new orbit size keys and `density_ratio`.
- Primordial mass segregation (`[merger.clusterN.segregation] method = "bdk08"`, `degree` S ∈ [0, 1];
  `SegregationSpec`): the scheme of Baumgardt, De Marchi & Kroupa (2008) on the systems before virialisation —
  energy-ranked phase-space slices of the unsegregated Plummer or King model handed out by cumulative mass — with
  the partial ordering of McLuster (Küpper et al. 2011, Eq. A1); the mass-density profile is preserved. Configurations
  without the table generate the same files as before.
- Massive primordial binaries (`[merger.clusterN.binaries.massive] model = "sana2012"`, `MassiveBinarySpec`): the
  stars above `m_min` are paired first with the period, mass-ratio and eccentricity power laws of Sana et al. (2012),
  the partner taken from the star list (`pairing = "pool"`, the IMF multiset kept) or given the mass q m₁
  (`"replace"`); the low-mass keys then apply to the remaining stars below `m_min`. `roche_fill_max` truncates the
  period and eccentricity laws of every pair so that neither component fills more than that fraction of its Roche
  lobe (Eggleton 1983) at periastron, with the engine's ZAMS radii at `merger.stellar.zmet` (`zams_radius`, the
  Tout et al. 1996 fit of SSE). Configurations without the new keys generate the same files as before.
- Complexes from one star list and their controls: `[merger.complex] imf_partition = "global"` draws one star
  list from the member IMF and deals it out to the members (`ExplicitMassesIMF` carries a member's stars), so the
  IMF realisation is the same across arms; `control = "monolithic"` (C-E: one virial Plummer cluster of the same
  stars with the complex's total energy) and `control = "cold"` (C-cold: the same energy and global virial ratio)
  replace the complex by its control for the same seed; the record holds the matched E and Q. With primordial
  binaries the control's energy differs from the target by the internal energy of its re-drawn pairs.
- Tidal frame of generated configurations: under `kz14 = 1, 2` the engine integrates in the frame rotating
  with the circular Galactic orbit (Ω along +z) and takes the `dat.10` velocities as velocities in that frame.
  `[merger.tidal] velocity_frame = "inertial"` (the default) writes v − Ω × r, `"rotating"` keeps the sampled
  velocities as they are. A Kepler pair whose apocentre reaches its pair Jacobi radius is refused (warned
  above half of it); `tidal_frame_rotation`, `pair_jacobi_radius`, `to_rotating_frame!`. The IC record gains Ω,
  the apocentre over the pair Jacobi radius and the bound 4u·max|x|/RS0 on the FP32 regular-force error from the
  coordinate offset, warned above 10⁻⁴.
- Binary spans from the event records (`binary_spans`, `BinarySpan`,
  `write_binary_spans`): the KS episodes of each pair of component names,
  joined when the gap is at most `n_p` periods or a chain or hierarchy holding
  both components covers it, with the origin (primordial by the name-parity
  rule, reformed, exchange, dynamical with GW or tidal capture and three-body
  sublabels), the end class (GW merger, stellar merger, collision, escape,
  exchange, supernova disruption, flyby, disruption, unresolved, censored)
  and the elements of the opening and closing records.
- `binary_physics = "inherited" | "full"` in `[merger.nbody6]`: the
  few-body and binary options of the generated input. `"inherited"` (the
  default) keeps the block written so far, without the chain, stable
  hierarchies or the gravitational-wave decay of compact pairs.
  `"full"` writes the set adopted for the binary programme: chain and
  hierarchies (`KZ(15) = 2`, `KZ(30) = 2`), slow-down `KZ(26) = 2`, Peters
  decay with sequential tides `KZ(27) = 3`, `KZ(28) = 1`, Roche transfer
  `KZ(34) = 2`, the hierarchy census `KZ(18) = 1`, binary records also
  without primordial pairs `KZ(8) = 2`, the binary event files `KZ(9) = 3`
  and KS event lines `KZ(10) = 2`; it requires stellar evolution. Explicit
  `kz` keys override either set; `KZ(18) ≥ 2` (the engine's primordial-triple
  generator) is refused and `KZ(46) > 0` (per-binary custom output) warns.
- `capacity` section of the validity record (schema version 2): the
  headroom of the engine's SIMD irregular-force library (KS pairs beyond
  `NBIN0 + 100` overrun it) and of the merger table of stable hierarchies
  (`MMAX − 1` entries), with a verdict. `AdjustRecord` gains `nmerge`,
  `DiagnosticsData` gains `nbin0` (from the `&INDATA` echo), and the build
  record's `engine_limits` gains `MMAX`.
- Validity record of a run (`validity_record`, `validity.toml`,
  `validity_snapshots.csv`): energy, mass and angular-momentum checks with a
  verdict, bound-set convergence, the snapshot cadence against the crossing
  and relaxation times, the windows that exclude snapshots from equilibrium
  estimators (violent relaxation, pericentres, coalescence phase), the times
  of virialisation, kinematic relaxation and mixing, the calibrated bias and
  scatter of the tracked estimators at the bound star count, and the masking
  rules; ten `[postprocess]` keys (`validity`, `validity_*`) hold the
  thresholds. `load_calibration` and `calibrated` read the shipped
  calibration table.
- Wall time of every post-processing stage in `results[:timings]`, the log
  and `RUN_INFO.toml` (`[pipeline.postprocess_seconds]`).
- `RemnantDiagnostics` keeps the per-pair density-centre separations and
  thresholds of the coalescence test (`coalescence_pairs`,
  `coalescence_separation`, `coalescence_threshold`).
- `rotation_analysis(...; peebles = false)` skips λ_P and its pair sum.
- Angular-momentum budget (`angular_momentum_budget`, `angular_momentum.csv`):
  the internal angular momentum of the stars (about their centre of mass,
  invariant under the engine's centre-of-mass corrections) closed against
  escapers, newcomers, stellar mass loss and the kicked remnants named by
  the engine's events, by exact set differences, with the angular-momentum
  content `l_scale` as the reference of the residual, the rotation of the bound set inside its Lagrangian
  radii and per initial-mass class about one spin axis (`spin_axis` of
  `rotation_analysis`), and the elapsed relaxation times `∫dt/t_rh`;
  `[postprocess]` keys `angular_momentum`, `mass_class_edges`,
  `mass_class_n_min`.
- Equipartition mass of Bianchini et al. (2016) (`equipartition_fit`,
  `equipartition_series`, `equipartition.csv`); keys `equipartition`,
  `equipartition_bin_width`, `equipartition_n_min`.
- Merger tree and per-progenitor budget for multi-cluster runs
  (`merger_tree`, `progenitor_budget`, `merger_tree.csv`,
  `progenitor_budget.csv`): persistent pairwise joins with union–find, the
  exact mass split of every progenitor into present, removed and shed mass,
  and the mixing of its bound members in the product; `coalescence_time`
  returns the pair list and the presence matrix; key `merger_tree`.
- Engine source amendments `tidal-energy-total` and `tidal-energy-total-mpi`:
  the whole-system call of `xtrnlv.F` is restored in the energy evaluation
  for a linearised tide (`KZ(14) = 1, 2`), so that the potential energy of
  the field enters `ETIDE`, the engine's total and its energy check, and the
  external virial term enters `Q`, as in the original NBODY6; the
  post-processing then treats the engine's `DE` of such a run as an error
  measure (`engine_energy_valid(mode; tidal_energy_amended)`, from the
  build record copied into the run).
- `scripts/calibrate_estimators.jl`: bias and variance of every estimator
  versus N on the static models, seeded with `StableRNGs`, written as the
  versioned table `assets/calibration/estimators.csv` with its sidecar
  `estimators.toml`.
- Static test models, each with a known answer for one estimator:
  `sample_plummer_om` (the Osipkov–Merritt Plummer sphere,
  `β(r) = r²/(r² + r_a²)`, `r_a ≥ 0.75 a`), `lynden_bell_flip!` (rotation
  by sign flips, `⟨v_φ⟩ = (2p − 1)⟨|v_φ|⟩`), `stretch_positions!` (axis
  ratios `q`, `s`) and `rank_masses_by_energy!` (energy-ranked mass
  segregation of degree `S`).
- Centre: `shrinking_sphere_centre` returns the centre with the radius and
  membership of the final set and the centroid noise
  `radius/√(5 n_eff)`.
- Core radius from the number density next to the mass density
  (`r_core_number`), and the core density (`rho_core`).
- `lagrangian_grid`: the 18 engine fractions by mass or number, masked
  where fewer than `n_min` effective particles lie on either side, with the
  binomial noise of every radius.
- `radial_profile` errors on the dispersions and the anisotropy, and
  logarithmic density shells merged to a least membership (`r_log`,
  `rho_log`, `n_log`).
- `relaxation_times` (half-mass and core relaxation times for two
  Coulomb-logarithm factors, crossing time) and `local_relaxation_time`.
- Mass segregation: the median-edge ratio `Λ̃_MSR` and the geometric-mean
  ratio `Γ_MSR` next to `Λ_MSR`, over a scan of tree sizes.
- Rotation: sign-flip null floors of `v_rot`, `λ_R` and `λ_P`, the
  debiased shell amplitude, the analytic floor `0.295/√N_eff` and the
  Bianchini rotation-curve fit (`v_peak`, `r_peak`).
- Shape, segregation and mixing (`ShapeShell`, `shape_tensor`,
  `projected_ellipticity`, `aplus`, `two_sample_d`, `mixing`): axis ratios
  from the iterative unweighted shape tensor, projected ellipticities, the
  A⁺ index with a jackknife error, the two-sample D with its null floor,
  and the mixing of two tagged populations.
- `core_collapse_markers` (smoothed minimum of the core radius, maximum of
  the core density) and the pair-based markers of the remnant diagnostics.
- `local_density`, the Casertano–Hut density on a k-d tree, and the
  snapshot binary finder (`snapshot_pairs`, `SnapshotPair`,
  `reconcile_pairs`): mutually bound neighbours joined most bound first,
  with a tidal filter, hardness and the comparison with the KS pairs.
- Twenty-five `[postprocess]` keys for the thresholds of these estimators:
  `centre_n_stop`, `centre_stop_fraction`, `centre_n_floor`,
  `centre_shrink`, `lagr_n_min`, `profile_nbins`, `profile_n_min_sigma`,
  `profile_n_min_beta`, `profile_n_min_density`, `density_neighbours`,
  `msr_n_tracer`, `msr_scan`, `msr_n_random`, `msr_n_random_small`,
  `msr_n_small`, `msr_seed`, `msr_threshold`, `rotation_n_min`,
  `rotation_nbins_max`, `rotation_n_null`, `shape_n_min_shell`,
  `shape_n_min_enclosed`, `pair_k`, `pair_tidal_factor`,
  `relaxation_gamma`.

### Changed

- The restart input carries the `&INSSE`, `&INBSE` and `&INCOLL` blocks of
  the original input after `&ININPUT`: the engine reads them again at a
  restart, and without them non-default stellar, binary-evolution and
  collision parameters reverted to the Level defaults.
- `resume_run` and `restart_simulation` read the subsystem counter `NSUB` from
  the header of a dump and do not restart from one written while a triple,
  quad or chain was active: such a dump does not hold the subsystem, and the
  engine crashes on it. `resume_run` takes the previous restartable dump with
  a warning; `restart_simulation` does the same by default and refuses such a
  `dump` named explicitly. Builds before `v2026.07.6` wrote one at a
  run's end and at an energy-check halt.
- The exact potential runs on the CPU as a threaded, vectorised loop over
  tiles; the potential of a set on itself visits every pair once. At 10⁵
  bodies on 8 threads it takes 1.9 s instead of 9.2 s, and it agrees with the
  previous kernel to 10⁻¹⁵ (summation order). The KernelAbstractions kernel
  remains as the device form.
- The angular-momentum budget no longer computes λ_P for every Lagrangian
  fraction and mass class, a pair sum per call that it did not use.
- The engine is cloned from the package's fork,
  [PaulGoG/Nbody6PPGPU-beijing](https://github.com/PaulGoG/Nbody6PPGPU-beijing),
  at the tag `v2026.07.9` by default (`install.source_url`,
  `install.ref`; the shipped pipeline configurations pin the same tag). The
  fork's branch `main` (named `campaign` until its tags `v2026.07-campaign.1`
  to `.6`; later tags are `v2026.07.N`) carries the three source amendments and the CUDA 13
  helper headers as commits on upstream v2026.07, and no longer tracks the
  outputs of `configure`, so a build from it records the amendments as found
  in place and its `backend_commit` is clean (an upstream checkout reports
  `-dirty` after every `configure`). `deps/engine/amendments.toml` and
  `deps/cuda/` remain for a build from an upstream checkout. The fork also
  corrects engine defects the amendments did not cover: an output interval
  without a short exact binary expansion no longer hangs the run at its first
  output; the end-of-run, STOP-file and halt dumps are written whatever
  `NCOMM`, on a fixed unit; the `COMMON SAVED` line prints its values under
  their labels; `configure` finds `nvcc` under `--with-cuda`; more than 32
  OpenMP threads no longer abort the AVX build; and the supernova-ghost and
  `lagr.7` fixes of upstream `dev` are included. From `campaign.4` the
  Milky-Way potential mode (`KZ(14) = 5`) runs: escapers and mass changes no
  longer carry the absolute Galactic potential into the energy bookkeeping,
  and the NFW and bulge jerks are corrected. From `campaign.5` the
  gravitational-wave decay of compact KS pairs under `KZ(27) > 2` follows the
  Peters equations: decay below a relative change of 10⁻⁴ per step is
  accumulated instead of discarded, the KS slow-down is included, the decay
  is integrated over the accumulated interval and applied at any orbital
  phase (an unperturbed eccentric pair gained eccentricity at every
  application, merged early or reached e > 1), the stellar-evolution look-up
  no longer applies it a second time, and a decaying pair is no longer
  treated as a chaotic binary; the speed of light is 299 792.458 km/s instead
  of 3×10⁵ km/s. Runs with `KZ(27) ≤ 2` are unchanged, and restart dumps of
  earlier builds cannot be read by `campaign.5`. From `campaign.6` chain
  regularization (`KZ(15) > 0`) runs: five FORMAT statements that gfortran
  rejects at run time are repaired (the first binary–binary chain ended every
  run); a run reaching its end while a triple, quad or chain is active
  terminates the subsystem and ends at the next adjustment (`END RUN
  DEFERRED`) instead of writing a dump that crashes the restart; a dump
  skipped during a chain is written at its end also under `KZ(2) = 1`; the
  limit on tiny chain steps acts per chain and terminates the chain, not the
  run; `KZ(26)` survives the absorption of an inert binary; and `MIX` no
  longer overwrites the `&INSSE` `nsflag`. From `v2026.07.7` a body inside
  a perturbed KS pair whose step has collapsed no longer stalls the block
  time (the pair is terminated), and every mass change of a KS component
  corrects the KS orbit and books its energy: below 10⁻⁸ N-body mass units the
  binding energy of a hard pair changed without a counterpart, which in a
  metal-poor cluster of 10⁴ stars cost up to 7×10⁻³ of the cluster energy per
  adjustment. From `v2026.07.8` a hyperbolic passage of two neutron stars or
  black holes radiates the quadrupole energy of the passage at its pericentre
  and is captured when that exceeds the orbital energy (`GR CAPTURE`,
  `GR FLYBY` lines), and small mass changes in KS pairs of merged hierarchies
  are booked as well. From `v2026.07.9` a merged hierarchy terminated by `SEARCH` is restored by `RESET`; the
  time offset of `KZ(35) > 0` shifts the reference epochs of the gravitational-wave decay of KS pairs; the KS
  slow-down factor applies once after a re-initialisation by that decay; `ROCHE` no longer applies a second
  decay to compact pairs under `KZ(27) > 2`; a FORMAT statement of a chain diagnostic is repaired; and the second
  component of a pair in `bev.82` evolves with its own main-sequence time. The fork's `NOTICE.md` lists
  the corrections by tag and states the provenance, the credit and the
  licensing. An existing
  engine checkout under `install.install_dir` is reused as it is; set
  `install.reinstall = true` once to replace an upstream checkout with the fork.
- `remnant_diagnostics` takes `settings::PostprocessConfig` instead of
  `n_massive`, `n_random`, `seed`, `lambda_threshold`, `nbins` and
  `overlap_factor`.
- `RadialProfile` has new trailing fields, and its kinematic shells are of
  equal membership; the figures are unchanged.
- `RotationProfile` and `RemnantDiagnostics` have new trailing fields, and
  `remnant_diagnostics.csv` has new columns and header scalars (a
  minor-version API change).
- `coalescence_time` uses per-tag bound sets with their density centres and
  half-mass radii and requires the single group to persist to the last
  snapshot; it returns the `separation` and `threshold` matrices and no
  longer depends on `overlap_factor`.
- The shrinking-sphere centre is number-weighted with the stop rule of
  Power et al. (2003) (was mass-weighted, stopping at 20 particles or 1 %).
- `mass_segregation` takes a tracer class and a scan of tree sizes (was the
  50 most massive stars by present mass against 50 random sets).
- The sixth-neighbour density runs on a k-d tree (was O(N²)).
- New dependency: NearestNeighbors 0.4.
- Generated initial conditions. Every cluster is truncated at its Roche lobe with respect to its nearest
  neighbour (Eggleton 1983, `roche_truncation_radius`) by default, since the Hill radius exceeds the inner Lagrange point at
  comparable masses; `[merger.output] truncate_rule = "hill"` restores the earlier radius. Under `kz14 = 1, 2` the
  sampled inertial velocities are converted to the engine's rotating frame by default (`velocity_frame =
  "inertial"`); `"rotating"` restores the earlier output. A truncated Kepler pair is shifted to its centre-of-mass
  frame. The global IMF partition of a complex hits the target mass of every member. Fractal centre velocities
  are scaled by 1/2 per generation (Goodwin & Whitworth 2004). Cluster tables are read strictly: a value of the
  wrong kind and a mix of flat and structured keys are refused, and errors in the member template name
  `[merger.complex.member]`. An invalid binary specification raises an `ArgumentError`. The new key
  `[merger.clusterN.segregation] max_points` bounds the number of model points of the segregation scheme, whose
  memory grows with the mass range of the IMF.
- Binary analysis. The q-restricted fraction of the population series takes the primary window 0.47–0.76 M☉ of
  Milone et al. (2012) by default and, when stellar types are available, main-sequence pairs only. A merged
  hierarchy (a negative centre-of-mass name in `conf.3`) enters both census levels, the KS pairs and the bound
  pairs, as its inner binary with the outer body as a single, and the energy budget no longer adds the `bev.82` energy of its inner
  binary a second time. The end class `:exchange` needs an exchange within the join tolerance of the span's close;
  an open chain episode closes at the next opening; the census keeps apart bound spans that a miss separates.
  Records that a restart writes twice are dropped from `coll.13` and `coal.24`. The stellar types of a snapshot
  fall back to the latest earlier `sev.83`/`bev.82` epoch when no epoch matches its time. The binary-span stage
  needs `KZ(10) ≥ 2` and `parse_stdout`; NBIN0 falls back to `merger_ic.toml` when the stdout capture lacks it.
  Unparseable ADJUST records are skipped with a warning, and a missing T* or M* gives `NaN`.
- Run infrastructure. A fresh run into a directory that holds a `RUN_INFO.toml` is refused; such a run continues
  with `resume_run` or `resume_pipeline`. A resubmit chain stops after a failed or killed segment without progress,
  or after two failed segments in a row, and records `chain_stopped` in the job ledger. Any exception while the
  engine runs terminates it and closes the segment. A sweep point counts as done only when its run completed, and
  ensembles use completed runs only. The launch takes the CUDA toolkit recorded at build time. An existing
  checkout at a commit other than `install.ref` is reported. Every segment record carries its own backend commit
  and hardware. Configurations that cannot be honoured are refused at load: `simulation.run_test` together with
  `postprocess.data_dir`, `postprocess.binary_spans` without `parse_stdout`, and an empty `binary_name`.

### Deprecated

- `preset = "campaign"` in `[merger.nbody6]` is now `preset = "baseline"`;
  the former name is read as `"baseline"` with a deprecation warning. The
  same holds for `binary_physics = "campaign"`, read as `"full"`.

## [0.5.0] — 2026-10-15

### Added

- Exact potential and bound sets of snapshots. `pair_potential` is the
  Newtonian pair sum in Float64, written as one KernelAbstractions kernel
  and run on the CPU threads; a subset form gives the potential of given
  sources at given points. `reduce_pairs` turns the regularised pairs into
  their centres of mass (`SystemSet`). `bound_set` replaces the bound-member
  selection of the remnant, per-cluster and binary diagnostics: iterative
  energetic unbinding with the survivors' potential updated by subtraction,
  pairs as centres of mass (with the companion's potential included, the
  internal binding of a hard pair kept an escaping pair in the bound set),
  and an explicit status (`converged`, `pass_limit` with a warning,
  `dissolved`) instead of a silent stop after eight passes. Its `T` and `W`
  feed the per-cluster virial ratio, the dispersions and the Peebles spin
  parameter, which no longer recompute the pair sum. Configuration keys
  `[postprocess] bound_max_passes` and `bound_min_members`.
- Snapshot integrals (`snapshot_integrals`, `SnapshotIntegrals`,
  `write_snapshot_integrals`): mass, centre of mass, linear and angular
  momentum, the stars' kinetic energy, the kinetic and potential energy of
  the reduced system (the engine's `EKIN` and `POT`), the regularised pairs'
  internal energy (`EBIN`) and the bound set, per snapshot; written by the
  pipeline as `snapshot_integrals.csv` (`[postprocess] snapshot_integrals`).
- Readers of the engine's events. `read_engine_events` returns the stamped
  event lines of `out1000` as `EngineEvent` records, sixteen kinds
  (`ENGINE_EVENT_KINDS`: KS regularisations, remnant formation, kicks,
  collisions, mergers, coalescences, escapes, supernovae without a remnant,
  chains), and
  `read_event_counters` the rows of `event.35` as `EventCounters`. The
  ADJUST routine evaluates the energy, removes the escapers and then writes
  its output, so every event line between two ADJUST lines is assigned to the
  interval that ends at the second.
- Energy-jump localisation and mass and energy budgets. `energy_intervals`
  returns an `EnergyEventAnalysis` with one `EnergyInterval` per adjustment
  interval: the engine's `DE`, the increments of its correction channels and
  the counts of the events in it, flagged when |DE| reaches the larger of a
  floor and a multiple of the robust scale of the run's DE;
  `energy_events_summary` counts the flagged intervals the log does not
  explain. The pipeline writes `energy_intervals.csv` and
  `engine_events.csv` (`[postprocess] energy_events`,
  `energy_jump_threshold`, `energy_jump_sigma`). `mass_budget`
  (`MassBudgetRow`, `mass_budget_summary`) sets the mass of every snapshot
  against the escapers, the stellar-evolution loss of `event.35` and the
  collision and coalescence losses, `energy_budget` (`EnergyBudgetRow`,
  `energy_budget_summary`) its Jacobi integral against the energy of the
  removed escapers and the engine's corrections, both from `event.35` at
  full precision when present, with the pairs' internal energy from
  `bev.82` where it covers the snapshot time and the engine's `DETOT` for
  comparison; written as `mass_budget.csv` and `energy_budget.csv`
  (`[postprocess] mass_budget`, `events_file`).
- Tidal handling. `AdjustRecord` keeps the energy terms of the ADJUST line
  (`EKIN`, `POT`, `ETIDE` and the correction channels), `DiagnosticsData`
  the options echo `kz` and the linearised tidal field `tidal`.
  `snapshot_integrals` takes `tidal` and `snapshot_integrals.csv` gains
  `w_tidal` and the Jacobi integral `e_jacobi`; `tidal_mode` and
  `engine_energy_valid` classify `KZ(14)`. The engine never adds the
  potential energy of a linearised field to its total, so its `DE` is a
  drift there and the pipeline flags no interval.
- The activation scripts resolve a manifest that a change of the project, or
  of the package an auxiliary environment takes by path, has left stale
  (`deps/activation.jl`); the manifests are not tracked, so a pulled
  dependency no longer fails at precompilation.
- Slurm job layer (`[cluster]`, `ClusterConfig`). `submit_run` and
  `scripts/run_cluster.jl` turn a run into a chain of queue jobs, one per
  engine segment: the job script is generated into `<run_dir>/cluster/`
  from the section's keys (partition, account, qos, counts, memory, gres,
  time limit, extra options, prologue lines) with `--signal=B:TERM` before
  the limit, and `cluster_segment` runs a segment inside the job (the
  pipeline on the first, a resume afterwards), resubmitting its successor
  (`chain = "resubmit"`) or relying on a chain submitted at once
  (`chain = "dependency"`). Job events are appended to `cluster/jobs.csv`
  and `<runs_dir>/cluster_jobs.csv` (`read_job_ledger`). A configuration
  whose segment does not fit the time limit, or whose time limit exceeds
  the partition's, is refused at load. `--dry-run` renders the script and
  the submission commands.

- Cost model (`[cost]`, `CostConfig`). `estimate_cost` derives the wall
  time, the number of wall-budget segments and the disk of a run from a
  table of seconds per N-body time unit against particle number
  (`cost_rate`; built-in rows measured on one H200 NVL device, replaceable
  through `table`), with a factor for two-cluster configurations. Every
  engine launch logs the estimate, writes it into the segment record
  (`cost_*` keys) and warns or refuses against `budget_hours`.

- Engines per device in sweeps. `[sweep] devices` and `engines_per_device`
  assign a CUDA device to each point when it starts (the least-loaded one)
  and admit at most that many engines per device, in addition to
  `concurrency`; the point's frozen `config.toml`, the index and the
  summary (`device` column) record the assignment.

- Surface-density rendering of particle projections. `plot_snapshot`,
  `plot_snapshot_evolution` and `animate_cluster` draw the projected surface
  mass density (log₁₀ Σ in M☉ pc⁻², one colour scale per figure or
  animation) instead of the mass-coloured scatter, whose markers overlap into
  a flat disc above a few thousand particles. `[visualization.style]` gains
  `snapshot_render` (`"auto"`, `"scatter"`, `"density"`), `density_min_n`,
  `density_bins`, `density_decades` and `density_mass_frac`; the limits of a
  density panel enclose that mass fraction, so escapers no longer set the
  scale.

- Energy bookkeeping of the engine in the diagnostics reader. `AdjustRecord`
  gains `de_abs` (the engine's `DELTA`, the energy change of one adjustment
  interval) and `detot` (`DETOT`, its running sum, an absolute energy);
  `cumulative_energy_error` returns the cumulative relative and absolute
  series per adjustment and `read_energy_totals` the totals of the `END RUN`
  line. Code that constructs `AdjustRecord` positionally passes two more
  arguments.

- Segment status in `RUN_INFO.toml`. Every engine launch is recorded when it
  starts (`status = "running"`, launch time, host, process id, package
  commit, start time in N-body units, stdout offset) and closed at exit with
  one of `completed`, `stopped`, `halted`, `watchdog`, `killed`, `failed`;
  `run.status` repeats the last segment's. A halt by the engine's energy
  check and a stop on request are told apart from a crash.

- Wall budget and stop requests. `simulation.wall_budget` asks the engine to
  stop `stop_margin` seconds before the budget expires, through the `STOP`
  file it polls for; it writes a restart dump and exits, and the segment is
  recorded as `stopped`. A driver that receives SIGTERM makes the same
  request for its running engines, waits for them and closes their records,
  so a job ended by a scheduler leaves a dump and no engine behind.

- Build-limit check before launch. The particle number, the number of
  primordial binaries and the neighbour number of an input are compared
  with `NMAX`, `KMAX` and `LMAX` of the configured engine tree
  (`include/params.h`), under the conditions the engine itself applies at
  start-up; an input that does not fit is refused with the limit and the
  remedy. `BUILD_INFO.toml` records the limits as `engine_limits`.

- Checkpoint mode of the generated merger input. `checkpoint = true` in
  `[merger.nbody6]` makes the engine write a restart dump at every
  adjustment (`KZ(2) = 1`); the generator logs the size of one dump and of
  the whole series.

- Retention of periodic restart dumps. `simulation.checkpoint_keep = k`
  keeps the newest `k` dumps of a checkpointed run and deletes older ones
  while the engine runs; `0`, the default, keeps all.

- `resume_run(run_dir)` continues an interrupted run to the end time of its
  input and can be called again safely. It starts from the last complete
  dump at the exact time the engine reported, and after a kill it moves
  what the dead segment wrote beyond that dump (stdout tail, later records
  of the series the package reads, later snapshots and dumps) to
  `output/discarded/`, so the series stay single-valued; nothing is deleted.

- `resume_pipeline(run_dir)` and `scripts/run_resume.jl` continue an
  interrupted run from the configuration frozen in its directory and, once
  it has reached its end time, perform post-processing and figures. The
  call can be repeated, one per queue segment.

- Single-cluster mode of the generator. `orbit_mode = "single"` with
  `n_clusters = 1` samples one cluster, places it at the origin at rest
  without truncation and writes the same files as a merger; `position`,
  `velocity` and `[merger.orbit]` are refused or ignored. The isolated
  controls of a sweep and `gpu/single_600k.toml` use it.

- Single-cluster science cases as generator configurations under
  `input_files/clusters/`: a globular cluster with a black-hole subsystem,
  runaway collisions in a dense young cluster, a young massive cluster with
  a binary fraction of 0.5, and an open cluster of Hyades-like mass with its escapers
  kept. Each has been run through its first time units on the CPU build.

- `merger.nbody6.preset`: `"engine"` (the rules of the engine's example
  inputs, the default) or `"campaign"` (the calibration baseline: η_I 0.01,
  η_R 0.02, NNBOPT 100, RMIN three times the 90° deflection scale of the
  smallest member, DTMIN 0.4 of the Aarseth close-encounter step, no
  automatic re-derivation of the integration parameters). A zero of
  `etai`, `etar`, `etau`, `nnbopt`, `rs0`, `rmin` or `dtmin` now means the
  preset's value or rule; explicit keys override it.

### Changed

- The static-QA dependency JET is 0.12 (was 0.9–0.10): the General registry
  no longer offers 0.10 on Julia 1.13, which left the test environment
  unresolvable. JET 0.12 reports a call on a regex capture group as a
  possible error on its `Nothing` branch, so the readers and the launcher
  narrow their captures with `something` where the group cannot be absent;
  `pair_potential` types its `backend` keyword as a KernelAbstractions
  backend. No behaviour changes.
- `[postprocess] pair_sum_max_n` defaults to 2 000 000 (was 100 000): the
  pair sum is threaded and computed once per snapshot. The four 6 × 10⁵-body
  probe pipelines set 100 000 and keep their diagnostics off.
- New dependency: KernelAbstractions 0.9.
- `AdjustRecord` and `DiagnosticsData` gained fields; the previous
  positional constructors remain. `SnapshotIntegrals` gained
  `tidal_potential` and `snapshot_integrals.csv` the columns `w_tidal` and
  `e_jacobi` (a minor-version API change).
- The launch script owns the stop request. `_launch.sh` stays the parent of
  the engine: it starts the engine in a subshell that ignores SIGTERM, turns
  a SIGTERM to itself into the `STOP` file, waits, and exits the way the
  engine did (its exit status, or the signal re-raised). A scheduler's
  job-wide SIGTERM therefore stops the engine with a dump; a termination by
  the package is a stop request followed, 15 s later, by a SIGKILL of the
  process group. Segment records gain `engine_pid` next to `pid` (the
  script); `output/engine.pid` holds it while the engine runs. An engine
  killed at its wall budget is recorded with exit status −9 (was −15); an
  operator interrupt closes the record as `stopped` with the dump when the
  engine honoured the request (`stop_requested = "interrupt"`), else
  `killed`. A SIGTERM sent to the engine process by hand no longer ends it;
  a SIGKILL of the launch script alone orphans the engine — kill
  `engine_pid`, or the process group.

- The defaults of `Nbody6ParameterSpec` for `etai`, `etar` and `etau` are
  `0.0` (resolved from the preset at generation time; the engine preset
  gives the former 0.02, 0.02, 0.1). A specification passed directly to
  `generate_merger_inp` must carry nonzero values, as it already must for
  the derived keys.

- The shipped engine inputs `N1k_quick` to `N100k_production` write a
  restart dump at every adjustment (`NCOMM = 1`, was 10: with the larger
  value the dump of an end, a stop or a halt was written one time in ten).
  The 1k and 5k inputs are isolated (`KZ(14) = 0`): they requested a
  point-mass field without its namelist, which the engine read as a galaxy
  of zero mass. The 25k and 100k inputs generate their primordial binaries
  with the engine's `BINPOP` (`KZ(8) = 1`): `KZ(8) = 2` with internally
  generated initial conditions declared the first `2·NBIN0` bodies to be
  pairs. The 10k input has no binaries and says so (`KZ(8) = 0`).

- The five science inputs `gc_bh_subsystem`, `imbh_runaway`, `pop3_cluster`,
  `tidal_tails` and `young_massive_binaries` are withdrawn: they generated
  Plummer models where their headers said King, declared binaries the
  engine never generated, quoted spans in Myr wrong by factors of three to
  fifty, and the Population III case cannot start (the engine refuses
  `ZMET` below 10⁻⁴). The configurations under `input_files/clusters/`
  replace four of them.

- With the default `snapshot_render = "auto"`, snapshots of 20 000 particles
  or more are drawn as surface density; `"scatter"` restores the former
  figures.

- The start-up watchdog and the completion monitor of a restart read only
  the part of the stdout capture written by that segment; an `END RUN` or
  an adjustment line of an earlier segment no longer counts.

- `merger.nbody6.ncomm` defaults to 1 (was 10). The engine counts every
  dump call and writes one in `ncomm`, so a larger value could skip the
  dump of a stop request or of an energy halt; inputs without periodic
  dumps are unaffected.

- A restart no longer keeps backups of `output/comm.1`: it is a working
  copy of a dump that stays under its own name, and it is removed when the
  engine exits. The restart input carries the time increment unrounded.

- `run_pipeline` ends after a simulation phase that stopped on request
  (`wall_budget`, or a signal): post-processing and figures are left to
  the completed run instead of being made from every segment.

### Fixed

- Snapshots no longer count regularised pairs twice. The engine writes
  NTOT = N + NPAIRS entries to `conf.3`, the N stars followed by one
  centre-of-mass entry per KS pair carrying the pair's mass, and
  `read_conf3` returned all of them as stars. Every sum over a snapshot
  therefore included each regularised pair's mass twice: up to 23 % excess
  mass in a 4920-star merger with 819 primordial pairs, 0.2 % in a
  48 000-body merger with one dynamically formed pair. The reader now drops the
  centre-of-mass entries, and the centre-of-mass body of an active chain
  regularisation (NAME = 0), which the engine leaves beside the chain
  members. `Snapshot` gains `npairs`, the components of pair k being stars
  2k − 1 and 2k (`pair_indices`), and `pair_name`, the names of the pairs'
  centre-of-mass entries (negative for a hierarchical merger), and validates
  its arrays on construction. Remnant diagnostics, per-cluster structure, profiles, the
  binary hardness scale and the snapshot figures of runs with regularised
  pairs change accordingly.
- `animate_cluster` fits its canvas to the square axis, as `plot_snapshot`
  does: a blank strip about a quarter of the frame wide separated the axis
  from the colourbar. Under adaptive zoom the widest tick labels of the
  animation are reserved, so the axis no longer shifts between frames.
- The binary-population figures no longer run an ungated O(N²) pair sum:
  their hard/soft scale selected the bound systems of every snapshot
  whatever `pair_sum_max_n`, forty minutes of the post-processing of the
  575 080-body probe. Above the limit the scale is now taken over every
  system, with a warning; `binary_scales` takes the limit as a keyword.
- The telemetry figure's legend is at most two columns wide. With a device
  present it carries six entries with their summary values, and three
  columns of those overran the canvas on every campaign run.
- A step that removes more than 1000 escapers no longer ends a run. The
  engine writes the escaper names of one adjustment interval with a format
  holding a thousand of them and aborted on the next (`Expected REAL … in
  formatted transfer`, exit status 2), as a 27-cluster cold collapse did at
  1442 escapers. The install phase now corrects that format in the engine
  tree before the build (`deps/engine/amendments.toml`) and records the
  amendments in effect as `source_amendments` in `BUILD_INFO.toml`; an
  engine built earlier needs a rebuild.
- GPU telemetry follows `simulation.gpu_list`. The sampler queried every
  device of the host, so a run on one of two cards recorded half its
  utilisation and the power of both; it now names the run's devices to
  `nvidia-smi`. With a device list the launch script also exports
  `CUDA_DEVICE_ORDER=PCI_BUS_ID`, so that an index means the same device
  for the engine and for `nvidia-smi`.
- The input-file reference described `TCRTP0` as a wall-clock limit in
  seconds. It is the engine's termination time in Myr; the run ends at the
  earlier of `TCRIT` and `TCRTP0`.

## [0.4.0] — 2026-09-26

### Added

- `[merger] virial_max_n` (default 200 000): the largest N, per cluster and
  for the combined system, for which the initial-condition generator
  evaluates the exact O(N²) potential. A cluster above it is refused when the
  configuration is loaded, naming the key, instead of failing inside the
  sampler; large initial conditions raise it deliberately.
- `[postprocess] pair_sum_max_n` (default 100 000): the largest snapshot
  particle count for which the O(N²) pair-sum diagnostics run (remnant
  diagnostics, per-cluster virial ratio and structure, bound-member
  profiles). Above it they are skipped with a warning naming the key; the
  readers and the other figures are unaffected. `generate_plots` takes the
  same limit as a keyword.
- `bench/fp32_peak_probe.jl`: the measured FP32 peak of the device with the
  peak probe of GPUDiagnostics.jl, written as `fp32_peak_<hostname>.toml`
  for the fraction-of-peak panel of `gpu_cells_figures.jl`. The script keeps
  its own shared environment (`@nb6_fp32_peak`) and installs
  GPUDiagnostics, KernelAbstractions and CUDA there on first use.
- The probes above 5 × 10⁵ bodies under `input_files/gpu/`:
  `merger_600k.toml` (two King clusters of 300 000 stars on the orbit of
  `merger_50k.toml`) and `single_600k.toml` (one King cluster of 600 000 at
  rest), each for two N-body time units, with the CUDA and AVX pipeline
  configurations that run them on the trees the recipe builds. The
  benchmark case raises its virialisation limit to the engine's `b1m`
  capacity.
- `run_gpu_validation`: the four probe pipelines as stages
  (`gpu_merger_600k`, `gpu_single_600k`, `cpu_merger_600k`,
  `cpu_single_600k`, opt-in) and `stop_on_failure` (`--stop-on-failure`),
  which skips every stage after a failed one and runs nothing needing `nvcc`
  when the host-compiler probe did not pass, so a workstation campaign is one
  gated invocation. The probe and benchmark stages inherit
  `JULIA_NUM_THREADS` or get `auto`.
- `bench/gpu_cells.jl`: benchmark cells on the device, a grid of (cell, N,
  host threads, `GPU_LIST`, MPI ranks) points run with post-processing off
  and read back into one CSV row per run (wall time, backend timing table,
  kernel throughput, GPU utilisation, power and memory, host RSS, CPU
  efficiency, machine identity, run ID). `bench/merger_case.jl` gains the
  single-cluster cell (`single_toml`) next to the two-cluster case and
  `cell_toml` to select either.
- `bench/gpu_cells_figures.jl`: the scaling figure (strong, weak and
  host-thread scaling of the two-device machine) and the cross-hardware
  figure (wall time per card with the regular-force share as a tick; kernel
  rate as a fraction of the FP32 peak) from the cell CSVs, a console table,
  and a `--synthetic` rehearsal mode.

### Changed

- `input_files/` is laid out by purpose: `engine/` (the ten single-cluster
  engine inputs), `mergers/` (the seven merger configurations),
  `verification/` (the two targets of the verification suite), `sweeps/`
  (the demonstration sweep), next to the unchanged `showcase/` and `gpu/`.
  File names are unchanged; `example_input` resolves a bare name across the
  folders, so `example_input("N1k_quick.inp")` and the other names of 0.3.0
  keep working, while a name that ships in more than one place must be given
  with its folder. The shipped pipeline and sweep TOMLs point at the new
  paths; a project's own configuration that named a shipped input by its old
  flat path must add the folder. `input_files/README.md` maps the folders to
  the documentation.
- The root `config.toml` runs the two-cluster demo merger
  (`mergers/merger_demo_small.toml`, seconds) instead of the 27-cluster
  cubic grid, a cold-collapse stress test outside the engine's regime that
  is kept under `mergers/` and described as such.
- The five science engine inputs carry a short factual header; their
  science cases, expected phenomena, cost estimates and references (with
  DOIs) moved to the Input Files page, which is organised by folder. A new
  documentation page, Showcase Cases, describes the five showcase
  configurations, their reference runs of 2026-09-11 (wall times, energy
  errors, run identifiers) and what a run produces; the flagship case is
  identified with the validation merger of `gpu/merger_50k.toml`, whose
  integration is timed on every validated host. The manual's validated
  hardware table carries the H200 NVL node. The default of
  `simulation.input_file` names a file that ships
  (`input_files/engine/N1k_quick.inp`); it named a path that did not exist.
- The package's author line in `Project.toml` carries the name only, as the
  licence does.

### Fixed

- A retried validation stage starts from a clean run directory: the
  crashed attempt's `runs/<run_id>` is kept as
  `runs/<run_id>.attempt<k>.signal<N>` beside its log, where the rerun used
  to append a second engine segment to the same record and accumulate the
  elapsed time. This makes `--retry-stages` usable for the pipeline stages,
  which the host with the JIT crashes needs.
- The start-up watchdog and the completion monitor terminate the engine
  before they report it, and `_terminate` sends SIGKILL before it reports the
  escalation. The report came first, and a pipeline whose driver had died
  wrote it to a closed pipe: the watchdog task ended on `EPIPE` with the
  kill never sent, and the engine it should have ended ran for thirteen
  hours beside its rerun (2026-09-25). The console sink of a run's logger
  now absorbs a failed write as well (`_ResilientLogger`), so a pipeline
  whose launcher is gone still writes its run record and its file log.
- An operator interrupt reaches the engine. The validation driver and every
  entry script that starts a stage or an engine set
  `Base.exit_on_sigint(false)`, so SIGINT raises `InterruptException` instead
  of ending the script at once — the forwarding of `dbf7274` never ran,
  because a non-interactive Julia exits on SIGINT before any handler. Inside
  a pipeline the helper tasks (telemetry sampler, start-up watchdog,
  completion monitor) hand an interrupt that lands in them to the task
  waiting on the engine, whose handler terminates it: SIGINT is delivered to
  whichever task is current, and under a long wait that is usually a helper.
  A SIGTERM or SIGKILL of a pipeline still leaves its engine (own session)
  running; the manual names the process sweep that ends the whole tree.
- The run record's `package_commit` and `backend_commit` are captured when
  the engine is launched, not when the record is written: a `git pull` on
  the host during a long run used to stamp the record with a commit that
  had not run.
- `startup_timeout` of the four 6 × 10⁵-body probe pipelines under
  `input_files/gpu/` raised from 3600 s to 14 400 s: the AVX build needs more
  than an hour to its first adjustment on a slow host (i7-10750H, EPYC 7551P),
  and the watchdog killed the `cpu_merger_600k` probe there at 3600 s.
- `gpu_cells_figures.jl`: FP32 peak records are matched by host *and*
  device, so two cards in hosts of the same name (one workstation name, two
  machines) each take their own measured peak instead of one overwriting the
  other in the lookup; the host comparison ignores the domain part, which
  `gethostname()` includes on some machines and the cards file does not.
- Stage processes of the validation driver and the engine itself start in
  their own session (`detach`). A terminal closing above a `nohup`-ed
  driver used to kill the running stage and its engine: `nohup` shields
  the driver alone, a child spawned by Julia has libuv's default signal
  dispositions. Two campaigns of 2026-09-23/24 lost their running stage on
  every host that way. An operator interrupt (Ctrl-C) is forwarded
  explicitly: SIGINT to a stage, SIGTERM then SIGKILL to the engine.
- `bench/gpu_cells_figures.jl`: the logarithmic floor of the hardware
  figure's wall-time panel is derived from the regular-force times as well
  as the wall times, so the regular-force tick of a bar is never clipped
  below the axis (it was, for every card whose regular-force time fell under
  the first decade of the wall times).

## [0.3.0] — 2026-09-22

First public release. Nbody6Dynamics.jl automates the lifecycle of
Nbody6++GPU star-cluster simulations from one TOML configuration: it clones
and builds the engine (CPU, GPU and MPI variants, CUDA architectures
resolved from the visible devices), generates multi-cluster merger initial
conditions, runs the integration with live monitoring, watchdogs, restarts
from the engine's dumps and hardware telemetry, reads every standard output
file into typed records, derives the diagnostics of the run, and draws
publication figures and animations through a Makie extension. Every run
directory carries its frozen configuration, the resolved environment, the
package and engine commits and a hardware fingerprint.

### What the release contains

- **Configuration.** One `config.toml` drives the pipeline; the merger
  generator and the parameter sweeps have their own TOML schemas. The
  parsers reject unknown keys and wrong types and raise an `ArgumentError`
  naming `section.key`; every default is stated once, in the configuration
  structs. Relative paths resolve against the directory of the file
  (`cfg.config_dir`), which is where the engine tree and the run directories
  are created, so an installation by URL behaves like a checkout;
  `example_input(name)` reaches the shipped inputs either way.
- **Engine build.** Clone at the validated commit, configure, HDF5 patch,
  parallel make; CUDA toolkit detection, `nvcc` host-compiler probe,
  `-gencode` entries per architecture, `BUILD_INFO.toml` next to the binary.
- **Merger initial conditions.** King and Plummer samplers, Kroupa IMF in
  natural, rescaled and equal-mass forms, primordial binaries, Kepler or
  explicit orbits, Jacobi truncation, external tidal fields, engine
  parameters derived from the smallest member cluster, output intervals as
  dyadic rationals the engine's digit counter can terminate on, seeds
  recorded for reproduction.
- **Runs.** Isolated `runs/<id>/` directories, a shell-safe launch script,
  OpenMP and `GPU_LIST` control, start-up and completion watchdogs
  (SIGTERM, then SIGKILL), restarts from COMMON dumps with per-segment
  bookkeeping, CPU, memory and GPU telemetry, records written atomically
  and never overwritten, a completion marker that requires the engine's
  `END RUN`.
- **Readers and diagnostics.** `conf.3`, `out1000`, `lagr.7`, `esc.11`,
  `sev.83_*` and `bev.82_*` in the engine's current and documented layouts,
  with skipped lines reported; per-cluster structure and radial profiles,
  binary population with the Heggie split, stellar classes and a class
  census, remnant diagnostics (core radius, rotation, mass segregation,
  coalescence time), written as CSV data products with or without a figure
  backend.
- **Sweeps and ensembles.** Cartesian grids over merger keys times seeds,
  run as concurrent worker processes with an index kept current, seeded
  ensembles with median and central intervals, isolated single-cluster
  controls per point.
- **Figures.** The publication layout: canvases in layout units exported at
  a printed width, Computer Modern type, one formatter for every number in
  figure text, a run-level HR layout with a census strip, montage,
  telemetry, sweep, ensemble and control figures, GIF animations. The figure
  routines are declared by the package and implemented by the
  `Nbody6DynamicsMakieExt` extension; without a backend a call is a
  `MethodError` whose message names the remedy, and the orchestration entry
  points refuse before doing any work.
- **Validation of a CUDA host.** `run_gpu_validation` runs the GPU-gated
  suite, the GPU and CPU pipelines and the scaling benchmark as logged
  stages with a retry policy, identifies the machine by hostname, hardware
  and GPU UUID, and judges each pipeline stage by the run directory it
  assigned.
- **Tooling.** A component-split test suite with static QA (Aqua, JET,
  ExplicitImports) and engine- and GPU-gated sets, a `Documenter` site with
  a Literate walkthrough and DOI-keyed references, benchmarks, and a
  `CITATION.cff`.

### Compatibility

Julia 1.13 or later; CairoMakie 0.13 or 0.15 (Makie 0.22 or 0.24) for the
figure extension, both series tested. The engine is Nbody6PPGPU-beijing at
its v2026.07 commit `618d7a4`, cloned at build time and not redistributed.
The strict parser refuses a frozen `config.toml` of an earlier run that
carries a key the schema has since dropped (`visualization.figsize`); remove
the key to load or restart such a run.

Earlier versions were private development releases; their history is in the
repository.
