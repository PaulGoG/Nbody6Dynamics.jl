# =============================================================================
# Gravitational-radiation, escape, mass-loss, tidal and Roche event lines of
# the engine log (out1000)
# =============================================================================
#
#   kind                  line label (after the leading blanks)    writer
#   :gw_decay_step        GR KSTIDE T DTGW STEP                    kstide.f
#   :gw_ksint             GR KSINT2 RES KMIN T DTGW STEP           ksint.f
#   :gw_ksint_coll        GR KSINT COLL: (list-directed)           ksint.f
#   :gw_tides             GR TIDES RES KMIN T DTGW STEP            unpert.f
#   :gw_unpert_coll       COLL DUE TO GW UNPERT (list-directed)    unpert.f
#   :gw_recoil            GW recoil:                               kickGW.f
#   :gr_brake             GR BRAKE T DTGR STEP                     mdot.F
#   :gr_brake_roche       GR BRAKE T DTGW STEP                     roche.f
#   :gr_brake3            GR BRAKE3 T M1 M2 K1 K2 ...              brake3.f
#   :gw_capture           GR CAPTURE T[NB]                         ksint.f (GWCAPT)
#   :gw_flyby             GR FLYBY   T[NB]                         ksint.f (GWCAPT)
#   :binary_escape        BINARY ESCAPE  IP,N,NCORR =              escape.F
#   :hierarchy_escape     HIARCH ESCAPE    NM =                    escape.F
#   :quad_escape          QUAD ESCAPE    KS =                      escape.F
#   :ks_end_mass_loss     MDEND KSREG TIME[NB]                     mdot.F
#   :ks_new_mass_loss     MDNEW KSREG TIME[NB]                     mdot.F
#   :binary_kick          BINARY KICK:  TIME[Myr]                  kick.F
#   :circ_new, :circ_end  NEW CIRC T DT, END CIRC T DT             kstide.f
#   :circularized         CIRCULARIZED T DT, INACTIIVE PHASE T DT  kstide.f
#   :capture              NEW CAPTURE    T NM                      kstide.f
#   :roche_new            NEW ROCHE    TIME[NB]                    roche.f
#   :roche_coal           ROCHE COAL   TIME[NB]                    roche.f
#   :roche_end            END ROCHE    TIME[NB]                    roche.f
#
# Every line is one physical record; the empty line that a leading `/` of a
# FORMAT writes before it is not part of the record. Formatted lines are read
# descriptor by descriptor with read_fixed, the two list-directed lines are
# split on blanks. coal.f's `  ROCHE COAL: IQCOLL`, the ` ESCAPE T,NESC =`
# summary and the second BRAKE3 line are not lines of this family.

"""
    EVOLUTION_EVENT_KINDS

The kinds returned by `_evolution_line_kind`: the gravitational-radiation
kinds of [`GWEvent`](@ref) (`:gw_decay_step`, `:gw_ksint`, `:gw_ksint_coll`,
`:gw_tides`, `:gw_unpert_coll`, `:gw_recoil`, `:gr_brake`, `:gr_brake_roche`,
`:gr_brake3`, `:gw_capture`, `:gw_flyby`), the escapes of [`BinaryEscapeEvent`](@ref) (`:binary_escape`,
`:hierarchy_escape`, `:quad_escape`), the mass-loss lines of
[`MassLossEvent`](@ref) (`:ks_end_mass_loss`, `:ks_new_mass_loss`,
`:binary_kick`), the tidal lines of [`TidalEvent`](@ref) (`:circ_new`,
`:circ_end`, `:circularized`, `:capture`) and the Roche lines of
[`RocheEvent`](@ref) (`:roche_new`, `:roche_coal`, `:roche_end`).
"""
const EVOLUTION_EVENT_KINDS = (
    :gw_decay_step,
    :gw_ksint,
    :gw_ksint_coll,
    :gw_tides,
    :gw_unpert_coll,
    :gw_recoil,
    :gr_brake,
    :gr_brake_roche,
    :gr_brake3,
    :gw_capture,
    :gw_flyby,
    :binary_escape,
    :hierarchy_escape,
    :quad_escape,
    :ks_end_mass_loss,
    :ks_new_mass_loss,
    :binary_kick,
    :circ_new,
    :circ_end,
    :circularized,
    :capture,
    :roche_new,
    :roche_coal,
    :roche_end,
)

const _EVO_GW_KINDS = EVOLUTION_EVENT_KINDS[1:11]
const _EVO_ESCAPE_KINDS = EVOLUTION_EVENT_KINDS[12:14]
const _EVO_MASS_LOSS_KINDS = EVOLUTION_EVENT_KINDS[15:17]
const _EVO_TIDAL_KINDS = EVOLUTION_EVENT_KINDS[18:21]
const _EVO_ROCHE_KINDS = EVOLUTION_EVENT_KINDS[22:24]

"""Absent Cartesian triple."""
const _EVO_NAN3 = (NaN, NaN, NaN)

# ---------------------------------------------------------------------------
# Record types
# ---------------------------------------------------------------------------

"""
    GWEvent

One gravitational-radiation or GR-braking line of the engine log, built by
`_parse_evolution_line`. `line` is the line number in the capture,
`interval` the adjustment interval (1 + the ADJUST lines above), `kind` one of:

| kind | writer | line |
|---|---|---|
| `:gw_decay_step` | `kstide.f` | `GR KSTIDE T DTGW STEP`, one per applied Peters step |
| `:gw_ksint` | `ksint.f` | `GR KSINT2 RES KMIN …`, after the GR KSTIDE line of the step |
| `:gw_ksint_coll` | `ksint.f` | `GR KSINT COLL:`, list-directed merger trigger |
| `:gw_tides` | `unpert.f` | `GR TIDES RES KMIN …`, unperturbed twin of `:gw_ksint` |
| `:gw_unpert_coll` | `unpert.f` | `COLL DUE TO GW UNPERT`, list-directed merger trigger |
| `:gw_recoil` | `kickGW.f` | `GW recoil:`, from `mix.f`/`coal.f` for two NS/BH |
| `:gw_capture` | `ksint.f` (GWCAPT) | `GR CAPTURE T[NB]`, a hyperbolic NS/BH passage that ends bound |
| `:gw_flyby` | `ksint.f` (GWCAPT) | `GR FLYBY   T[NB]`, the same passage left unbound |
| `:gr_brake` | `mdot.F` | `GR BRAKE T DTGR STEP`, BSE look-up, DJGR > 10⁻⁶ |
| `:gr_brake_roche` | `roche.f` | `GR BRAKE T DTGW STEP`, inside a Roche stage |
| `:gr_brake3` | `brake3.f` | `GR BRAKE3 T M1 M2 …`, hierarchy with an NS/BH inner pair |

Fields, with the kinds that carry them (absent: `NaN`, `0`, `KSTAR_ABSENT`):

| field | unit | kinds and meaning |
|---|---|---|
| `time_nb` | NB | TTOT; `NaN` for the two COLL lines (no time printed) |
| `time_myr` | Myr | `:gw_recoil` TTOT·T* |
| `res`, `kmin` | — | `:gw_ksint`, `:gw_tides`: RES (1 = decay applied), KMIN |
| `index_cm` | — | global c.m. index I: Peters kinds, `:gr_brake_roche` |
| `index_star` | — | `:gr_brake`: index I of the star being evolved |
| `index1`, `index2` | — | `:gw_recoil`: global indices I1, I2 |
| `ipair` | — | Peters kinds, `:gr_brake_roche` |
| `npert` | — | perturber count LIST(1,I1) or LIST(1,J1) |
| `name1`, `name2`, `name_cm` | — | component and c.m. names |
| `kstar1`, `kstar2`, `kstar_cm` | — | stellar types; `:gr_brake3` KW1, KW2 of the inner pair |
| `icase` | — | `:gw_recoil` ICASE (2) |
| `m1_msun`, `m2_msun` | M☉ | current masses |
| `m01_msun`, `m02_msun` | M☉ | `:gw_recoil` initial masses BODY0 |
| `ecc`, `semi_nb` | —, NB | e and a (Peters kinds: before the step) |
| `semi_rsun` | R☉ | `:gr_brake3` SEP, the BSE semi-major axis |
| `qperi_nb` | NB | pericentre; Peters kinds: after the step |
| `qperi_old_nb` | NB | `:gw_decay_step` a(1 − e) before the step (label oldQP) |
| `rschw_nb` | NB | Schwarzschild radius 2M/c² of the pair |
| `sep_merge_nb` | NB | `:gw_ksint_coll` SEPMRG = 25 RSCHW |
| `c_light_nb` | NB | the two COLL lines: CLIGHT |
| `h_nb` | NB | H: `:gw_decay_step` before the update, `:gw_ksint`/`:gw_tides` after |
| `gamma` | — | relative perturbation |
| `a_ein` | rad | 1PN periastron advance per orbit 6πM/(a c²(1 − e²)) |
| `t_orb_nb` | NB | orbital period 2π(a³/M)^½ |
| `t_gr_nb` | NB | Peters kinds TGR = abs(Δt·a/Δa), i.e. a/abs(ȧ), not the merger time |
| `t_gr_yr` | yr | `:gr_brake` TGR = JORB/DJGR |
| `t_gr_steps` | — | `:gr_brake_roche` TGR = JORB/(DJGR·DTM0), in DTM0 steps |
| `dt_gw_nb` | NB | DTGW(IPAIR): Peters kinds, `:gr_brake_roche` |
| `dt_gr_yr` | yr | DTGR: `:gr_brake` 0.02·JORB/abs(J̇), `:gr_brake3` |
| `dt_bse_yr` | yr | `:gr_brake` DTXMIN, `:gr_brake_roche` DTM0 |
| `step_nb` | NB | STEP(I1) or STEP(J1) |
| `da_nb`, `de` | NB, — | `:gw_decay_step` Δa, Δe; `de` also `:gr_brake_roche` DELET1 |
| `d_energy_nb`, `dh_nb` | NB | `:gw_decay_step` orbital-energy change μΔH and ΔH |
| `j_orb` | M☉ R☉² yr⁻¹ | orbital angular momentum: `:gr_brake*` |
| `dj_dt` | M☉ R☉² yr⁻² | GR loss rate DJGR: `:gr_brake`, `:gr_brake3` |
| `dj` | M☉ R☉² yr⁻¹ | `:gr_brake_roche` loss over the step DJGR·DTM0 |
| `de_dt_per_yr` | yr⁻¹ | `:gr_brake` GR de/dt DELET |
| `spin1`, `spin2` | — | `:gw_recoil` ASPN after the kick model |
| `q_ratio` | — | `:gw_recoil` mass ratio ≤ 1 |
| `j_spin1_nb`, `j_spin2_nb`, `spin_factor` | NB, — | `:gw_recoil` XJSPIN, SPNFAC |
| `v_kick_kms`, `v_kick_xyz_kms` | km s⁻¹ | `:gw_recoil` kick speed and components |
| `v_esc_kms`, `v_dis_kms`, `r_over_semi`, `eb_nb`, `r_nb` | | `:gw_recoil`, see below |
| `x1_pc`, `v1_kms`, `x2_pc`, `v2_kms` | pc, km s⁻¹ | `:gw_recoil`, absolute frame, pre-kick |
| `v_inf_kms` | km s⁻¹ | capture kinds: v∞ = √(2H₀)·V* before the passage |
| `ecc_minus_1` | — | capture kinds: e − 1 of the incoming hyperbola |
| `r_capture_nb` | NB | capture kinds: Quinlan–Shapiro capture pericentre (RPMAX) |
| `h_after_nb` | NB | capture kinds: H after the passage (`h_nb` before it) |
| `one_minus_e_after` | — | capture kinds: 1 − e = r_p/a after the passage (< 0: still unbound) |
| `overflow` | — | a field was printed as asterisks (read as absent) |

Corrections and hazards. The K*S field of `:gw_decay_step`, `:gw_ksint` and
`:gw_tides` holds the temporary KSTAR(ICM) = −25 that selects the Peters
branch; `kstar_cm` is absent for these kinds. `mdot.F` prints the names and
types of `:gr_brake` in the reverse order of the masses; they are stored in
the order of the masses (`name1`, `kstar1` belong to `m1_msun`). The IPAIR of
`:gr_brake` is a local variable set only on the Roche branches of `mdot.F`
and is stored as absent. The DTGW label of `:gr_brake_roche` holds DTGW(IPAIR)
in NB, its TGR is dimensionless; its DELET1 is read as the change of e over
DTM0 (not settled from the code). `:gw_unpert_coll` has no rank test in
`unpert.f`: an MPI run writes one copy per rank, each read as its own record.
The COLL lines carry neither names nor time; the merger they trigger is the
following `CMBODY` collision (IQCOLL = −2). The capture kinds are written at
the pericentre (`time_nb` is the pericentre time, `qperi_nb` the pericentre
distance, `d_energy_nb` the radiated energy of the passage, `npert` the
perturbers after it); a passage that changes the orbital energy by less than
10⁻⁴ is not written. The `:gw_recoil` VESC, VDIS,
R/SEMI, EB and R use IPAIR = KVEC(I1) after `KSTERM` and probably refer to a
terminated pair: they are stored as printed and are not reliable.
"""
Base.@kwdef struct GWEvent
    line::Int
    interval::Int
    time_nb::Float64
    kind::Symbol
    time_myr::Float64 = NaN
    res::Int = 0
    kmin::Int = 0
    index_cm::Int = 0
    index_star::Int = 0
    index1::Int = 0
    index2::Int = 0
    ipair::Int = 0
    npert::Int = 0
    name1::Int = 0
    name2::Int = 0
    name_cm::Int = 0
    kstar1::Int = KSTAR_ABSENT
    kstar2::Int = KSTAR_ABSENT
    kstar_cm::Int = KSTAR_ABSENT
    icase::Int = 0
    m1_msun::Float64 = NaN
    m2_msun::Float64 = NaN
    m01_msun::Float64 = NaN
    m02_msun::Float64 = NaN
    ecc::Float64 = NaN
    semi_nb::Float64 = NaN
    semi_rsun::Float64 = NaN
    qperi_nb::Float64 = NaN
    qperi_old_nb::Float64 = NaN
    rschw_nb::Float64 = NaN
    sep_merge_nb::Float64 = NaN
    c_light_nb::Float64 = NaN
    h_nb::Float64 = NaN
    gamma::Float64 = NaN
    a_ein::Float64 = NaN
    t_orb_nb::Float64 = NaN
    t_gr_nb::Float64 = NaN
    t_gr_yr::Float64 = NaN
    t_gr_steps::Float64 = NaN
    dt_gw_nb::Float64 = NaN
    dt_gr_yr::Float64 = NaN
    dt_bse_yr::Float64 = NaN
    step_nb::Float64 = NaN
    da_nb::Float64 = NaN
    de::Float64 = NaN
    d_energy_nb::Float64 = NaN
    dh_nb::Float64 = NaN
    j_orb::Float64 = NaN
    dj_dt::Float64 = NaN
    dj::Float64 = NaN
    de_dt_per_yr::Float64 = NaN
    spin1::Float64 = NaN
    spin2::Float64 = NaN
    q_ratio::Float64 = NaN
    j_spin1_nb::Float64 = NaN
    j_spin2_nb::Float64 = NaN
    spin_factor::Float64 = NaN
    v_kick_kms::Float64 = NaN
    v_kick_xyz_kms::NTuple{3,Float64} = _EVO_NAN3
    v_esc_kms::Float64 = NaN
    v_dis_kms::Float64 = NaN
    r_over_semi::Float64 = NaN
    eb_nb::Float64 = NaN
    r_nb::Float64 = NaN
    v_inf_kms::Float64 = NaN
    ecc_minus_1::Float64 = NaN
    r_capture_nb::Float64 = NaN
    h_after_nb::Float64 = NaN
    one_minus_e_after::Float64 = NaN
    x1_pc::NTuple{3,Float64} = _EVO_NAN3
    v1_kms::NTuple{3,Float64} = _EVO_NAN3
    x2_pc::NTuple{3,Float64} = _EVO_NAN3
    v2_kms::NTuple{3,Float64} = _EVO_NAN3
    overflow::Bool = false
end

"""
    BinaryEscapeEvent

One escape line of a KS pair, a hierarchy or a quadruple ghost binary, written
by `ESCAPE` (`escape.F`, called from `ADJUST`) and built by
`_parse_evolution_line`. Kinds: `:binary_escape` (`BINARY ESCAPE`, a real KS
pair or the outer orbit of a merged hierarchy), `:hierarchy_escape`
(`HIARCH ESCAPE`, the inner binary of that hierarchy, after its BINARY
ESCAPE line), `:quad_escape` (`QUAD ESCAPE`, the ghost binary of a
quadruple). No escape variant beyond these three has a FORMAT in the
layout specification.

| field | unit | kinds and meaning |
|---|---|---|
| `time_nb` | NB | TTOT of the following ` ESCAPE T,NESC =` summary (see below) |
| `ipair` | — | KS pair index: `:binary_escape`, `:quad_escape` |
| `n_after` | — | `:binary_escape` particle number N after the removal |
| `ncorr` | — | `:binary_escape` escaper counter of this call |
| `name1`, `name2` | — | components; `:hierarchy_escape` inner primary, ghost |
| `name_cm` | — | `:binary_escape` c.m. name JLIST(NCORR), < 0 for a hierarchy |
| `kstar1`, `kstar2` | — | component types (`:hierarchy_escape`: inner c.m. body, ghost) |
| `kstar_cm` | — | c.m. type; `:hierarchy_escape` KSTARM, inner c.m. type at merger |
| `ls_flag` | — | LIST(2,·), −1 for a primordial-flagged pair |
| `m1_msun`, `m2_msun` | M☉ | masses; `:quad_escape` CM(3:4) of the ghost binary |
| `semi_nb`, `ecc` | NB, — | a, e (hierarchy line: inner; binary line of a hierarchy: outer) |
| `eb_nb` | NB | binding energy (`:quad_escape` EB1, `:hierarchy_escape` μ·HM) |
| `period_d` | d | period |
| `rstar_over_peri` | — | `:binary_escape` max(R*)/pericentre (label R*/PM) |
| `v_nb`, `r_nb` | NB | `:binary_escape` c.m. speed, distance from the density centre |
| `e_cm_nb` | NB | `:binary_escape` c.m. energy in the cluster field EI |
| `radius1_rsun`, `radius2_rsun` | R☉ | `:binary_escape` stellar radii |
| `pmin_over_pcrit` | — | `:hierarchy_escape` outer pericentre over its stability limit |
| `overflow` | — | a field was printed as asterisks (read as absent) |

Corrections. The lines print no time: `time_nb` is the TTOT of the next
` ESCAPE T,NESC =` line, which closes the same `ESCAPE` call; the summary is
looked up without being consumed, and `time_nb` is `NaN` when an ADJUST line
or the end of the capture comes first. `BINARY ESCAPE` prints the c.m. speed
under the label RI and the distance under VI; they are stored as `v_nb` and
`r_nb`. Its K*S field carries +30 when the c.m. is a merged hierarchy
(NAME < 0); `kstar_cm` holds KSTAR(c.m.) with the offset removed. On the
`QUAD ESCAPE` branch NAME1, NAME2 are not assigned (left from a previous
escaper) and KSTAR(JM) indexes the particle arrays with a merger-table index,
so `name1`, `name2` and `kstar_cm` are absent for `:quad_escape` (reading of
the specification, not settled from the code). The F5.1 masses of
`:hierarchy_escape`/`:quad_escape` have 0.1 M☉ resolution and F8.4 energies
keep four decimals.
"""
Base.@kwdef struct BinaryEscapeEvent
    line::Int
    interval::Int
    time_nb::Float64
    kind::Symbol
    ipair::Int = 0
    n_after::Int = 0
    ncorr::Int = 0
    name1::Int = 0
    name2::Int = 0
    name_cm::Int = 0
    kstar1::Int = KSTAR_ABSENT
    kstar2::Int = KSTAR_ABSENT
    kstar_cm::Int = KSTAR_ABSENT
    ls_flag::Int = 0
    m1_msun::Float64 = NaN
    m2_msun::Float64 = NaN
    semi_nb::Float64 = NaN
    ecc::Float64 = NaN
    eb_nb::Float64 = NaN
    period_d::Float64 = NaN
    rstar_over_peri::Float64 = NaN
    v_nb::Float64 = NaN
    r_nb::Float64 = NaN
    e_cm_nb::Float64 = NaN
    radius1_rsun::Float64 = NaN
    radius2_rsun::Float64 = NaN
    pmin_over_pcrit::Float64 = NaN
    overflow::Bool = false
end

"""
    MassLossEvent

One KS termination or re-regularisation by stellar mass loss (`mdot.F`) or a
binary kick line (`kick.F`), built by `_parse_evolution_line`. Kinds:
`:ks_end_mass_loss` (`MDEND KSREG`, written before `KSTERM` for a large
mass loss, a pair made soft by it, or a type change to a compact remnant),
`:ks_new_mass_loss` (`MDNEW KSREG`, the former components re-regularised,
written before `KSREG`), `:binary_kick` (`BINARY KICK:`, a KS pair with a
member receiving a natal kick).

| field | unit | kinds and meaning |
|---|---|---|
| `time_nb` | NB | TIME + TOFF; `NaN` for `:binary_kick`, which prints TTOT·T* only |
| `time_myr` | Myr | `:binary_kick` TTOT·T* |
| `name1`, `name2`, `name_cm` | — | components, c.m. name |
| `kstar1`, `kstar2`, `kstar_cm` | — | stellar types |
| `ipair`, `dtau_nb` | —, NB | pair index, regularised KS step |
| `m1_nb`, `m2_nb` | NB | masses |
| `r12_nb` | NB | separation; `:binary_kick` R(IPAIR) |
| `ecc`, `semi_nb`, `eb_nb` | —, NB, NB | e, a, binding energy (`:binary_kick` EB = μH) |
| `period_d` | d | period |
| `h_nb`, `gamma` | NB, — | specific binding energy, relative perturbation |
| `step_cm_nb` | NB | c.m. step |
| `npert`, `n_neighbours` | — | KS perturbers, c.m. neighbours |
| `m1_msun`, `m2_msun` | M☉ | masses (MDEND/MDNEW: before this step's loss) |
| `radius1_rsun`, `radius2_rsun` | R☉ | stellar radii |
| `r12_rsun` | R☉ | separation (label RAD1,2,S[*]) |
| `r_nb`, `v_nb` | NB | c.m. distance from the density centre, c.m. speed |
| `dm_nb`, `dmr` | NB, — | mass lost in this look-up step, relative to the mass |
| `spin1`, `spin2`, `q_ratio` | — | `:binary_kick` ASPN, QR = m₁/m₂ (not ≤ 1) |
| `j_spin1_nb`, `j_spin2_nb`, `spin_factor` | NB, — | `:binary_kick` XJSPIN, SPNFAC |
| `v_esc_kms`, `v_dis_kms` | km s⁻¹ | `:binary_kick` (2M/R + 4)^½·V*, (2M/R)^½·V* |
| `r_over_semi` | — | `:binary_kick` R(IPAIR)/SEMI |
| `x1_pc`, `v1_kms`, `x2_pc`, `v2_kms` | pc, km s⁻¹ | `:binary_kick`, absolute frame, pre-kick |
| `overflow` | — | a field was printed as asterisks (read as absent) |

Corrections. Field 26 of the KSREG layout, labelled RAD1,2,S[*], is the
separation in R☉ (`r12_rsun`). `BINARY KICK` labels the separation R12[au]
but prints R(IPAIR) in N-body units (`r12_nb`). In `MDNEW KSREG` only the
time, the component names and types, the masses [NB, M☉], the radii and DM,
DMR are current: its e, a, EB repeat the preceding MDEND line, IPAIR, DTAU,
R12, H, GAMMA and the derived period and R12[*] belong to whatever pair now
occupies the old index, the S name and type are those of the evolving star,
and STEP, NPERT, NB(ICM), RI, VI refer to that star; all of these are stored
as absent. The orbit of the re-regularised pair is in the `NEW KSREG` line
that follows. In `MDEND KSREG` the type of the evolving component may
already be negated (random-phase flag of a new NS/BH, `mdot.F`); it is
stored as printed. An MDEND line not followed by MDNEW/NEW KSREG of the same
names is a termination, not necessarily a disruption.
"""
Base.@kwdef struct MassLossEvent
    line::Int
    interval::Int
    time_nb::Float64
    kind::Symbol
    time_myr::Float64 = NaN
    name1::Int = 0
    name2::Int = 0
    name_cm::Int = 0
    kstar1::Int = KSTAR_ABSENT
    kstar2::Int = KSTAR_ABSENT
    kstar_cm::Int = KSTAR_ABSENT
    ipair::Int = 0
    dtau_nb::Float64 = NaN
    m1_nb::Float64 = NaN
    m2_nb::Float64 = NaN
    r12_nb::Float64 = NaN
    ecc::Float64 = NaN
    semi_nb::Float64 = NaN
    eb_nb::Float64 = NaN
    period_d::Float64 = NaN
    h_nb::Float64 = NaN
    gamma::Float64 = NaN
    step_cm_nb::Float64 = NaN
    npert::Int = 0
    n_neighbours::Int = 0
    m1_msun::Float64 = NaN
    m2_msun::Float64 = NaN
    radius1_rsun::Float64 = NaN
    radius2_rsun::Float64 = NaN
    r12_rsun::Float64 = NaN
    r_nb::Float64 = NaN
    v_nb::Float64 = NaN
    dm_nb::Float64 = NaN
    dmr::Float64 = NaN
    spin1::Float64 = NaN
    spin2::Float64 = NaN
    q_ratio::Float64 = NaN
    j_spin1_nb::Float64 = NaN
    j_spin2_nb::Float64 = NaN
    spin_factor::Float64 = NaN
    v_esc_kms::Float64 = NaN
    v_dis_kms::Float64 = NaN
    r_over_semi::Float64 = NaN
    x1_pc::NTuple{3,Float64} = _EVO_NAN3
    v1_kms::NTuple{3,Float64} = _EVO_NAN3
    x2_pc::NTuple{3,Float64} = _EVO_NAN3
    v2_kms::NTuple{3,Float64} = _EVO_NAN3
    overflow::Bool = false
end

"""
    TidalEvent

One tidal-circularisation or capture line of `kstide.f` (KZ(27) > 0, called
from `KSINT`/`UNPERT` at pericentre), built by `_parse_evolution_line`.
Kinds: `:circ_new` (`NEW CIRC`, start of sequential circularisation, KSTAR(ICM)
set to 9), `:circ_end` (`END CIRC`, e below ECCM, KSTAR(ICM) set to 10),
`:circularized` (`CIRCULARIZED`, and its twin `INACTIIVE PHASE`, label
misspelt in the source, flagged by `inactive_phase`), `:capture`
(`NEW CAPTURE`, a hyperbolic pair made bound).

| field | unit | kinds and meaning |
|---|---|---|
| `time_nb` | NB | TTOT (`:capture` TIME + TOFF in F9.2, 0.01 NB resolution) |
| `inactive_phase` | — | `true` for the `INACTIIVE PHASE` line |
| `step_nb` | NB | KS step STEP(I1) (label DT, not an interval between events) |
| `index_cm`, `ipair`, `npert` | — | c.m. global index, pair index, perturbers |
| `name1`, `name2`, `name_cm` | — | names (`:capture`: components only) |
| `kstar1`, `kstar2`, `kstar_cm` | — | types; `kstar_cm` already updated (9 or 10) |
| `m1_msun`, `m2_msun` | M☉ | masses |
| `ecc0`, `ecc` | — | e before and after the tidal step |
| `ecc_max` | — | ECCM = 0.002 |
| `semi0_nb`, `semi_nb` | NB | a before and after (angular momentum conserved) |
| `qperi_nb` | NB | pericentre passed to `KSTIDE` |
| `qperi_over_rstar` | — | `:capture` QPERI/max(R*) |
| `h0_nb`, `gamma` | NB, — | H before the step, relative perturbation |
| `period_d` | d | `:circ_new`, `:circ_end`: period after the step (from a) |
| `period0_d` | d | `CIRCULARIZED`: period before the step (from a0) |
| `radius1_rsun`, `radius2_rsun` | R☉ | stellar radii |
| `r12_rsun`, `r_coll_rsun` | R☉ | separation, collision distance 0.75(R₁ + R₂) |
| `r_nb`, `v_nb` | NB | abs(X), abs(Ẋ) of the c.m., not relative to the density centre |
| `overflow` | — | a field was printed as asterisks (read as absent) |

Corrections. The PD field of NEW/END CIRC is the period after the step, that
of CIRCULARIZED the period before it; each is stored under its meaning. The
`INACTIIVE PHASE` path does not assign PD (value of an earlier call, reading
of the specification), so neither period is stored for it; its leading
blanks are not fixed by the specification and the line is read from its
first non-blank character. RI, VI of CIRCULARIZED and INACTIIVE PHASE are
written with D exponents and read as such.
"""
Base.@kwdef struct TidalEvent
    line::Int
    interval::Int
    time_nb::Float64
    kind::Symbol
    inactive_phase::Bool = false
    step_nb::Float64 = NaN
    index_cm::Int = 0
    ipair::Int = 0
    npert::Int = 0
    name1::Int = 0
    name2::Int = 0
    name_cm::Int = 0
    kstar1::Int = KSTAR_ABSENT
    kstar2::Int = KSTAR_ABSENT
    kstar_cm::Int = KSTAR_ABSENT
    m1_msun::Float64 = NaN
    m2_msun::Float64 = NaN
    ecc0::Float64 = NaN
    ecc::Float64 = NaN
    ecc_max::Float64 = NaN
    semi0_nb::Float64 = NaN
    semi_nb::Float64 = NaN
    qperi_nb::Float64 = NaN
    qperi_over_rstar::Float64 = NaN
    h0_nb::Float64 = NaN
    gamma::Float64 = NaN
    period_d::Float64 = NaN
    period0_d::Float64 = NaN
    radius1_rsun::Float64 = NaN
    radius2_rsun::Float64 = NaN
    r12_rsun::Float64 = NaN
    r_coll_rsun::Float64 = NaN
    r_nb::Float64 = NaN
    v_nb::Float64 = NaN
    overflow::Bool = false
end

"""
    RocheEvent

One Roche-lobe overflow line of `roche.f` (KZ(34) > 0), built by
`_parse_evolution_line`. Kinds: `:roche_new` (`NEW ROCHE`, start of an
overflow stage, KSTAR(ICM) incremented to an odd value), `:roche_coal`
(`ROCHE COAL   TIME[NB]`, coalescence decided in the Roche stage; not coal.f's
`ROCHE COAL: IQCOLL` line), `:roche_end` (`END ROCHE`, the donor back inside
its lobe).

| field | unit | kinds and meaning |
|---|---|---|
| `time_nb` | NB | TTOT |
| `name1`, `name2`, `name_cm` | — | NAME(J1), NAME(J2) in donor order, c.m. name |
| `kstar1`, `kstar2`, `kstar_cm` | — | KW1, KW2, KSTAR(ICM) (odd: active stage) |
| `ipair`, `dtau_nb` | —, NB | pair index, regularised KS step |
| `time_myr` | Myr | TPHYS = TTOT·T* |
| `age_myr` | Myr | TPHYS − EPOCH: `:roche_new`, `:roche_coal` |
| `m1_nb`, `m2_nb` | NB | BODY(J1), BODY(J2) |
| `spin1_nb`, `spin2_nb` | NB | `:roche_new` spin angular momenta |
| `ospin1_nb`, `ospin2_nb`, `spin_factor` | NB, — | `:roche_new` OSPIN/SPNFAC, SPNFAC |
| `j_orb` | M☉ R☉² yr⁻¹ | `:roche_new` orbital angular momentum |
| `omega_orb_per_yr` | yr⁻¹ | `:roche_new` orbital frequency 2π/TB |
| `r12_nb` | NB | separation R(IPAIR) |
| `ecc`, `semi_nb`, `eb_nb` | —, NB, NB | e, a = SEP/SU, binding energy (`:roche_new` only) |
| `period_d` | d | TK (absent for `:roche_coal`) |
| `h_nb`, `gamma`, `step_cm_nb` | NB, —, NB | H, relative perturbation, c.m. step |
| `npert`, `n_neighbours` | — | perturbers LIST(1,J1), c.m. neighbours |
| `m01_msun`, `m02_msun` | M☉ | initial masses MASS0 |
| `m1_msun`, `m2_msun` | M☉ | current masses MASS |
| `mc1_msun`, `mc2_msun` | M☉ | core masses MASSC |
| `r_nb`, `v_nb` | NB | `:roche_new` c.m. distance from the density centre, c.m. speed |
| `r12_rsun` | R☉ | instantaneous separation (label SEP[*]) |
| `radius1_rsun`, `radius2_rsun` | R☉ | stellar radii |
| `roche_radius1_rsun`, `roche_radius2_rsun` | R☉ | Roche-lobe radii |
| `semi_rsun`, `rl1_rsun` | R☉ | `:roche_new` semi-major axis SEP (label SEMI[*]), RL1 |
| `dm1_msun`, `dm2_msun`, `dt_myr` | M☉, Myr | `:roche_coal`, `:roche_end` DM1, DM2, DTM |
| `coalescence` | — | `:roche_coal` logical COALS |
| `ixxx` | — | `:roche_coal` IXXX |
| `overflow` | — | a field was printed as asterisks (read as absent) |

Corrections. The label SEP[*] holds the instantaneous separation R·SU and
SEMI[*] the semi-major axis SEP; they are stored as `r12_rsun` and
`semi_rsun`. ROCHE COAL prints TK = 0 (`roche.f`), stored as an absent
period. RI, VI and EB are formed only on the NEW ROCHE path; in ROCHE COAL and
END ROCHE they repeat an earlier value or are undefined (reading of the
specification) and are stored as absent. The logical COALS is glued to the
preceding E13.5 field and is read from its single column.
"""
Base.@kwdef struct RocheEvent
    line::Int
    interval::Int
    time_nb::Float64
    kind::Symbol
    name1::Int = 0
    name2::Int = 0
    name_cm::Int = 0
    kstar1::Int = KSTAR_ABSENT
    kstar2::Int = KSTAR_ABSENT
    kstar_cm::Int = KSTAR_ABSENT
    ipair::Int = 0
    dtau_nb::Float64 = NaN
    time_myr::Float64 = NaN
    age_myr::Float64 = NaN
    m1_nb::Float64 = NaN
    m2_nb::Float64 = NaN
    spin1_nb::Float64 = NaN
    spin2_nb::Float64 = NaN
    ospin1_nb::Float64 = NaN
    ospin2_nb::Float64 = NaN
    spin_factor::Float64 = NaN
    j_orb::Float64 = NaN
    omega_orb_per_yr::Float64 = NaN
    r12_nb::Float64 = NaN
    ecc::Float64 = NaN
    semi_nb::Float64 = NaN
    eb_nb::Float64 = NaN
    period_d::Float64 = NaN
    h_nb::Float64 = NaN
    gamma::Float64 = NaN
    step_cm_nb::Float64 = NaN
    npert::Int = 0
    n_neighbours::Int = 0
    m01_msun::Float64 = NaN
    m02_msun::Float64 = NaN
    m1_msun::Float64 = NaN
    m2_msun::Float64 = NaN
    mc1_msun::Float64 = NaN
    mc2_msun::Float64 = NaN
    r_nb::Float64 = NaN
    v_nb::Float64 = NaN
    r12_rsun::Float64 = NaN
    radius1_rsun::Float64 = NaN
    radius2_rsun::Float64 = NaN
    roche_radius1_rsun::Float64 = NaN
    roche_radius2_rsun::Float64 = NaN
    semi_rsun::Float64 = NaN
    rl1_rsun::Float64 = NaN
    dm1_msun::Float64 = NaN
    dm2_msun::Float64 = NaN
    dt_myr::Float64 = NaN
    coalescence::Bool = false
    ixxx::Int = 0
    overflow::Bool = false
end

"""Record types returned by `_parse_evolution_line`."""
const _EvolutionEventRecord = Union{GWEvent,BinaryEscapeEvent,MassLossEvent,TidalEvent,RocheEvent}

# ---------------------------------------------------------------------------
# FORMATs, transcribed descriptor by descriptor
# ---------------------------------------------------------------------------

"""
Descriptors shared by the GR lines of `kstide.f`, `ksint.f`, `unpert.f`, `mdot.F`
and `roche.f`, from ` I IP NPERT` through `TGR(PM)=`: `I10,2I6`, `3I10`, `3I4`,
`2E13.5`, `4E13.5`, `2E13.5`, `3E13.5`; `qp_label` is the writer's e,a,QP,RS label.
"""
function _evo_gw_core_format(qp_label::String)::Vector{FortranItem}
    return vcat(
        fmt_lit(" I IP NPERT"),
        fmt_int(1, 10),
        fmt_int(2, 6),
        fmt_lit(" NM1,2,S="),
        fmt_int(3, 10),
        fmt_lit(" KW1,2,S="),
        fmt_int(3, 4),
        fmt_lit(" M1,2[M*]"),
        fmt_real(2, 13),
        fmt_lit(qp_label),
        fmt_real(4, 13),
        fmt_lit(" H, GAMMA="),
        fmt_real(2, 13),
        fmt_lit(" A_EIN, TORB, TGR(PM)="),
        fmt_real(3, 13),
    )
end

# kstide.f 665: 1X,' GR KSTIDE T DTGW STEP',1P,3E13.5,<core>,' dsemi,decc=',2E13.5,
# ' de,dh,oldQP=',3E13.5
const _EVO_GR_KSTIDE_FORMAT = vcat(
    fmt_lit("  GR KSTIDE T DTGW STEP"),
    fmt_real(3, 13),
    _evo_gw_core_format(" e,a,QP,RS[NB]="),
    fmt_lit(" dsemi,decc="),
    fmt_real(2, 13),
    fmt_lit(" de,dh,oldQP="),
    fmt_real(3, 13),
)

# ksint.f 665: 1X,' GR KSINT2 RES KMIN T DTGW STEP',2I3,1P,3E13.5,<core>
const _EVO_GR_KSINT2_FORMAT = vcat(
    fmt_lit("  GR KSINT2 RES KMIN T DTGW STEP"),
    fmt_int(2, 3),
    fmt_real(3, 13),
    _evo_gw_core_format(" e,a,QP,RS[NB]="),
)

# unpert.f 665: 1X,' GR TIDES RES KMIN T DTGW STEP',2I4,1P,3E13.5,<core>
const _EVO_GR_TIDES_FORMAT = vcat(
    fmt_lit("  GR TIDES RES KMIN T DTGW STEP"),
    fmt_int(2, 4),
    fmt_real(3, 13),
    _evo_gw_core_format(" e,a,QP,RS[NB]="),
)

# mdot.F 665: 1X,' GR BRAKE T DTGR STEP',1P,3E13.5,<core>,' j,dj=',2E13.5,
# ' DTXMIN,de=',2E13.5
const _EVO_GR_BRAKE_FORMAT = vcat(
    fmt_lit("  GR BRAKE T DTGR STEP"),
    fmt_real(3, 13),
    _evo_gw_core_format(" e,a,QP,RS[NB]="),
    fmt_lit(" j,dj="),
    fmt_real(2, 13),
    fmt_lit(" DTXMIN,de="),
    fmt_real(2, 13),
)

# roche.f 665: 1X,' GR BRAKE T DTGW STEP',1P,3E13.5,<core with RSCH>,' j,dj=',2E13.5,
# ' DTM0,de=',2E13.5
const _EVO_GR_BRAKE_ROCHE_FORMAT = vcat(
    fmt_lit("  GR BRAKE T DTGW STEP"),
    fmt_real(3, 13),
    _evo_gw_core_format(" e,a,QP,RSCH[NB]="),
    fmt_lit(" j,dj="),
    fmt_real(2, 13),
    fmt_lit(" DTM0,de="),
    fmt_real(2, 13),
)

# brake3.f 45: ' GR BRAKE3 T M1 M2 K1 K2 SEP ECC0 JORB DJ DTGR ',1P,3E14.5,2I4,5E14.5
const _EVO_GR_BRAKE3_FORMAT = vcat(
    fmt_lit(" GR BRAKE3 T M1 M2 K1 K2 SEP ECC0 JORB DJ DTGR "),
    fmt_real(3, 14),
    fmt_int(2, 4),
    fmt_real(5, 14),
)

# ksint.f GWCAPT 20: 1X,A10 (WHICH = 'GR CAPTURE' or 'GR FLYBY'),' T[NB]',1P,E17.10,
# ' NM1,2=',2I10,' KW1,2=',2I4,' M1,2[M*]',2E12.4,' VINF[km/s]',E12.4,' E-1',E15.7,
# ' QP,RCAP[NB]',E15.7,E12.4,' DE[NB]',E15.7,' H0,H1[NB]',2E15.7,' 1-E1',E15.7,' NP',I6
function _evo_gr_capture_format(which::String)::Vector{FortranItem}
    return vcat(
        fmt_lit(" " * rpad(which, 10)),
        fmt_lit(" T[NB]"),
        fmt_real(1, 17),
        fmt_lit(" NM1,2="),
        fmt_int(2, 10),
        fmt_lit(" KW1,2="),
        fmt_int(2, 4),
        fmt_lit(" M1,2[M*]"),
        fmt_real(2, 12),
        fmt_lit(" VINF[km/s]"),
        fmt_real(1, 12),
        fmt_lit(" E-1"),
        fmt_real(1, 15),
        fmt_lit(" QP,RCAP[NB]"),
        fmt_real(1, 15),
        fmt_real(1, 12),
        fmt_lit(" DE[NB]"),
        fmt_real(1, 15),
        fmt_lit(" H0,H1[NB]"),
        fmt_real(2, 15),
        fmt_lit(" 1-E1"),
        fmt_real(1, 15),
        fmt_lit(" NP"),
        fmt_int(1, 6),
    )
end
const _EVO_GR_CAPTURE_FORMAT = _evo_gr_capture_format("GR CAPTURE")
const _EVO_GR_FLYBY_FORMAT = _evo_gr_capture_format("GR FLYBY")

# kickGW.f 310: A15 (WHICH = ' GW recoil:    '),' TIME[NB] I12 N12 K12 ',1P,E17.9,4I10,3I4,
# ' TIME[Myr] M012[M*] M12[M*] ',5E13.5,' ASPN12,QR,XJSPIN12,SPNFAC ',6E13.5,
# ' VKICK GW[km/s] (V,vx,vy,vz)',4E13.5,' VESC[km/s] VDIS[km/s] R/SEMI EB R ',5E13.5,
# ' X1(1:3)[pc] V1(1:3)[km/s] X2(1:3)[pc] V2(1:3)[km/s] ',12E13.5
const _EVO_GW_RECOIL_FORMAT = vcat(
    fmt_lit(" GW recoil:    "),
    fmt_lit(" TIME[NB] I12 N12 K12 "),
    fmt_real(1, 17),
    fmt_int(4, 10),
    fmt_int(3, 4),
    fmt_lit(" TIME[Myr] M012[M*] M12[M*] "),
    fmt_real(5, 13),
    fmt_lit(" ASPN12,QR,XJSPIN12,SPNFAC "),
    fmt_real(6, 13),
    fmt_lit(" VKICK GW[km/s] (V,vx,vy,vz)"),
    fmt_real(4, 13),
    fmt_lit(" VESC[km/s] VDIS[km/s] R/SEMI EB R "),
    fmt_real(5, 13),
    fmt_lit(" X1(1:3)[pc] V1(1:3)[km/s] X2(1:3)[pc] V2(1:3)[km/s] "),
    fmt_real(12, 13),
)

# escape.F 200: /,' BINARY ESCAPE  IP,N,NCORR =',I5,I10,I5,'  NM1,2,S=',3I10,
# '  KW1,2,S,LS=',4I4,'  M1,2* =',1P,2E13.5,'  a,e,eb,pb =',4E13.5,' R*/PM =',E13.5,
# '  RI,VI =',2E13.5,'  EI =',E13.5,'  RAD1,2[*]',2E13.5
const _EVO_BINARY_ESCAPE_FORMAT = vcat(
    fmt_lit(" BINARY ESCAPE  IP,N,NCORR ="),
    fmt_int(1, 5),
    fmt_int(1, 10),
    fmt_int(1, 5),
    fmt_lit("  NM1,2,S="),
    fmt_int(3, 10),
    fmt_lit("  KW1,2,S,LS="),
    fmt_int(4, 4),
    fmt_lit("  M1,2* ="),
    fmt_real(2, 13),
    fmt_lit("  a,e,eb,pb ="),
    fmt_real(4, 13),
    fmt_lit(" R*/PM ="),
    fmt_real(1, 13),
    fmt_lit("  RI,VI ="),
    fmt_real(2, 13),
    fmt_lit("  EI ="),
    fmt_real(1, 13),
    fmt_lit("  RAD1,2[*]"),
    fmt_real(2, 13),
)

# escape.F 240: /,' HIARCH ESCAPE    NM =',2I6,'  K* =',3I3,'  M =',2F5.1,'  EB =',F8.4,
# '  E =',F7.3,'  PM/PC =',1P,E8.1,'  A =',E8.1,'  P =',E8.1
const _EVO_HIARCH_ESCAPE_FORMAT = vcat(
    fmt_lit(" HIARCH ESCAPE    NM ="),
    fmt_int(2, 6),
    fmt_lit("  K* ="),
    fmt_int(3, 3),
    fmt_lit("  M ="),
    fmt_real(2, 5),
    fmt_lit("  EB ="),
    fmt_real(1, 8),
    fmt_lit("  E ="),
    fmt_real(1, 7),
    fmt_lit("  PM/PC ="),
    fmt_real(1, 8),
    fmt_lit("  A ="),
    fmt_real(1, 8),
    fmt_lit("  P ="),
    fmt_real(1, 8),
)

# escape.F 202: /,' QUAD ESCAPE    KS =',I5,'  NM =',2I6,'  K* =',2I3,'  M =',2F5.1,
# '  EB =',F8.4,'  E =',F7.3,'  A =',1P,E8.1,'  P =',E8.1
const _EVO_QUAD_ESCAPE_FORMAT = vcat(
    fmt_lit(" QUAD ESCAPE    KS ="),
    fmt_int(1, 5),
    fmt_lit("  NM ="),
    fmt_int(2, 6),
    fmt_lit("  K* ="),
    fmt_int(2, 3),
    fmt_lit("  M ="),
    fmt_real(2, 5),
    fmt_lit("  EB ="),
    fmt_real(1, 8),
    fmt_lit("  E ="),
    fmt_real(1, 7),
    fmt_lit("  A ="),
    fmt_real(1, 8),
    fmt_lit("  P ="),
    fmt_real(1, 8),
)

"""
`mdot.F` FORMATs 16 (MDEND) and 60 (MDNEW): `label`, `1P,E17.10`, `3I10`, `3I4`,
`I9`, `E11.3`, `2E11.3`, `E11.3`, `2E12.4,E11.3`, five `E11.3`, `I5`, `I5`,
`2E11.3`, `3E11.3`, `2E11.3`, `2E13.5`.
"""
function _evo_mdot_format(label::String)::Vector{FortranItem}
    return vcat(
        fmt_lit(label),
        fmt_real(1, 17),
        fmt_lit(" NM1,2,S="),
        fmt_int(3, 10),
        fmt_lit(" KW1,2,S="),
        fmt_int(3, 4),
        fmt_lit(" IPAIR"),
        fmt_int(1, 9),
        fmt_lit(" DTAU"),
        fmt_real(1, 11),
        fmt_lit(" M1,2[NB]"),
        fmt_real(2, 11),
        fmt_lit(" R12[NB]"),
        fmt_real(1, 11),
        fmt_lit(" e,a,eb[NB]="),
        fmt_real(2, 12),
        fmt_real(1, 11),
        fmt_lit(" P[d]="),
        fmt_real(1, 11),
        fmt_lit(" H"),
        fmt_real(1, 11),
        fmt_lit(" GAMMA"),
        fmt_real(1, 11),
        fmt_lit(" STEP(ICM)"),
        fmt_real(1, 11),
        fmt_lit(" NPERT"),
        fmt_int(1, 5),
        fmt_lit(" NB(ICM)"),
        fmt_int(1, 5),
        fmt_lit(" M1,2[*]"),
        fmt_real(2, 11),
        fmt_lit(" RAD1,2,S[*]"),
        fmt_real(3, 11),
        fmt_lit(" RI,VI[NB]="),
        fmt_real(2, 11),
        fmt_lit(" DM,DMR="),
        fmt_real(2, 13),
    )
end

const _EVO_MDEND_FORMAT = _evo_mdot_format(" MDEND KSREG TIME[NB]")
const _EVO_MDNEW_FORMAT = _evo_mdot_format(" MDNEW KSREG TIME[NB]")

# kick.F 1: ' BINARY KICK:  TIME[Myr] N12 K12 ',1P,E17.9,2I10,2I4,
# ' ASPN12,QR,XJSPN12[NB],SPNFAC ',6E13.5,' M12[M*] ',
# ' VESC[km/s] VDIS[km/s] R12/SEMI EB R12[au] ',7E13.5,
# ' X1(1:3)[pc] V1(1:3)[km/s] X2(1:3)[pc] V2(1:3)[km/s] ',12E13.5
const _EVO_BINARY_KICK_FORMAT = vcat(
    fmt_lit(" BINARY KICK:  TIME[Myr] N12 K12 "),
    fmt_real(1, 17),
    fmt_int(2, 10),
    fmt_int(2, 4),
    fmt_lit(" ASPN12,QR,XJSPN12[NB],SPNFAC "),
    fmt_real(6, 13),
    fmt_lit(" M12[M*] "),
    fmt_lit(" VESC[km/s] VDIS[km/s] R12/SEMI EB R12[au] "),
    fmt_real(7, 13),
    fmt_lit(" X1(1:3)[pc] V1(1:3)[km/s] X2(1:3)[pc] V2(1:3)[km/s] "),
    fmt_real(12, 13),
)

"""
`kstide.f` FORMATs 608 (NEW/END CIRC, `1X,A8,' T DT '`) and 638 (CIRCULARIZED):
`label`, then `1P,2E13.5`, `I10,2I6`, `3I10`, `3I4`, `2E13.5`, `6E13.5`, `2E13.5`,
`5E13.5`, `2E13.5` (`2D13.5` in 638, same width).
"""
function _evo_circ_format(label::String)::Vector{FortranItem}
    return vcat(
        fmt_lit(label),
        fmt_real(2, 13),
        fmt_lit(" I IP NPERT"),
        fmt_int(1, 10),
        fmt_int(2, 6),
        fmt_lit(" NM1,2,S="),
        fmt_int(3, 10),
        fmt_lit(" KW1,2,S="),
        fmt_int(3, 4),
        fmt_lit(" M1,2[M*]"),
        fmt_real(2, 13),
        fmt_lit(" e0,e,em,a0,a,QP="),
        fmt_real(6, 13),
        fmt_lit(" H, GAMMA="),
        fmt_real(2, 13),
        fmt_lit(" PD, R1,2,SEP,RCOLL[R*]="),
        fmt_real(5, 13),
        fmt_lit(" RI,VI="),
        fmt_real(2, 13),
    )
end

const _EVO_NEW_CIRC_FORMAT = _evo_circ_format(" NEW CIRC T DT ")
const _EVO_END_CIRC_FORMAT = _evo_circ_format(" END CIRC T DT ")
const _EVO_CIRCULARIZED_FORMAT = _evo_circ_format("  CIRCULARIZED T DT ")
# The leading blanks of the INACTIIVE PHASE line are not given verbatim
# [UNSURE]: the label is matched from the first non-blank character.
const _EVO_INACTIVE_PHASE_FORMAT = _evo_circ_format("INACTIIVE PHASE T DT ")

# kstide.f 35: ' NEW CAPTURE    T NM E EF QP/R* A1 ',F9.2,2I6,2F9.4,1P,2E11.3
const _EVO_CAPTURE_FORMAT = vcat(
    fmt_lit(" NEW CAPTURE    T NM E EF QP/R* A1 "),
    fmt_real(1, 9),
    fmt_int(2, 6),
    fmt_real(2, 9),
    fmt_real(2, 11),
)

"""
Opening of `roche.f` FORMATs 8, 76, 77: `A12` (`which`), `'  TIME[NB]'`,
`1P,E17.10`, `3I10`, `3I4`, `I9`, `E13.5`.
"""
function _evo_roche_head_format(which::String)::Vector{FortranItem}
    return vcat(
        fmt_lit(which),
        fmt_lit("  TIME[NB]"),
        fmt_real(1, 17),
        fmt_lit(" NM1,2,S="),
        fmt_int(3, 10),
        fmt_lit(" KW1,2,S="),
        fmt_int(3, 4),
        fmt_lit(" IPAIR"),
        fmt_int(1, 9),
        fmt_lit(" DTAU"),
        fmt_real(1, 13),
    )
end

"""
Middle of `roche.f` FORMATs 8, 76, 77, from `e,a,eb[NB]=` through
`RAD(1:2) ROL(1:2)[*]=`: `3E13.5`, four `E13.5`, `I5`, `I5`, `6E13.5`, `2E13.5`,
`E13.5`, `4E13.5`.
"""
function _evo_roche_body_format()::Vector{FortranItem}
    return vcat(
        fmt_lit(" e,a,eb[NB]="),
        fmt_real(3, 13),
        fmt_lit(" P[d]="),
        fmt_real(1, 13),
        fmt_lit(" H"),
        fmt_real(1, 13),
        fmt_lit(" GAMMA"),
        fmt_real(1, 13),
        fmt_lit(" STEP(ICM)"),
        fmt_real(1, 13),
        fmt_lit(" NPERT"),
        fmt_int(1, 5),
        fmt_lit(" NB(ICM)"),
        fmt_int(1, 5),
        fmt_lit(" M0(1:2) M(1:2) MC(1:2)[*]"),
        fmt_real(6, 13),
        fmt_lit(" RI,VI[NB]="),
        fmt_real(2, 13),
        fmt_lit(" SEP[*]"),
        fmt_real(1, 13),
        fmt_lit(" RAD(1:2) ROL(1:2)[*]="),
        fmt_real(4, 13),
    )
end

# roche.f 8: A12,'  TIME[NB]',…,' TIME,AGE[Myr] ',2E13.5,' M1,2[NB]',2E13.5,
# ' SPIN12,OSPIN12[NB] SPNFAC JORB,OORB[*] ',7E13.5,' R12[NB]',E13.5,<body>,
# ' SEMI[*],RL1[*]=',2E13.5
const _EVO_NEW_ROCHE_FORMAT = vcat(
    _evo_roche_head_format(" NEW ROCHE  "),
    fmt_lit(" TIME,AGE[Myr] "),
    fmt_real(2, 13),
    fmt_lit(" M1,2[NB]"),
    fmt_real(2, 13),
    fmt_lit(" SPIN12,OSPIN12[NB] SPNFAC JORB,OORB[*] "),
    fmt_real(7, 13),
    fmt_lit(" R12[NB]"),
    fmt_real(1, 13),
    _evo_roche_body_format(),
    fmt_lit(" SEMI[*],RL1[*]="),
    fmt_real(2, 13),
)

# roche.f 76: A12,'  TIME[NB]',…,' T[Myr] ',E13.5,' M1,2[NB]',2E13.5,' R12[NB]',E13.5,
# <body>,' DM1/2,DT=',3E13.5
const _EVO_END_ROCHE_FORMAT = vcat(
    _evo_roche_head_format(" END ROCHE  "),
    fmt_lit(" T[Myr] "),
    fmt_real(1, 13),
    fmt_lit(" M1,2[NB]"),
    fmt_real(2, 13),
    fmt_lit(" R12[NB]"),
    fmt_real(1, 13),
    _evo_roche_body_format(),
    fmt_lit(" DM1/2,DT="),
    fmt_real(3, 13),
)

# roche.f 77: as 76, then ' AGE COALS, IXXX=',E13.5,L1,I3; the L1 is skipped here
# and read from its column by _evo_logical.
const _EVO_ROCHE_COAL_FORMAT = vcat(
    _evo_roche_head_format(" ROCHE COAL "),
    fmt_lit(" T[Myr] "),
    fmt_real(1, 13),
    fmt_lit(" M1,2[NB]"),
    fmt_real(2, 13),
    fmt_lit(" R12[NB]"),
    fmt_real(1, 13),
    _evo_roche_body_format(),
    fmt_lit(" DM1/2,DT="),
    fmt_real(3, 13),
    fmt_lit(" AGE COALS, IXXX="),
    fmt_real(1, 13),
    fmt_skip(1),
    fmt_int(1, 3),
)

"""Column of the `L1` field COALS in a ROCHE COAL line."""
const _EVO_ROCHE_COAL_L1_COLUMN =
    sum(item -> item.width, @view(_EVO_ROCHE_COAL_FORMAT[1:(end - 2)])) + 1

# ---------------------------------------------------------------------------
# Classifier and parser (driver interface)
# ---------------------------------------------------------------------------

"""
    _evolution_line_kind(s::AbstractString) -> Symbol

Kind of the log line `s` (already `lstrip`ped) among
[`EVOLUTION_EVENT_KINDS`](@ref), identified by its full label; `:none` for
any other line, in particular coal.f's `ROCHE COAL: IQCOLL`, the
`ESCAPE T,NESC =` summary, the second `BRAKE3` line, `brake.f`'s `BRAKE T STEP`
and every line of the regularisation, chain and merger family.
"""
function _evolution_line_kind(s::AbstractString)::Symbol
    startswith(s, "GR KSTIDE T DTGW STEP") && return :gw_decay_step
    startswith(s, "GR KSINT2 RES KMIN T DTGW STEP") && return :gw_ksint
    startswith(s, "GR KSINT COLL:") && return :gw_ksint_coll
    startswith(s, "GR TIDES RES KMIN T DTGW STEP") && return :gw_tides
    startswith(s, "COLL DUE TO GW UNPERT") && return :gw_unpert_coll
    startswith(s, "GW recoil:") && return :gw_recoil
    startswith(s, "GR BRAKE T DTGR STEP") && return :gr_brake
    startswith(s, "GR BRAKE T DTGW STEP") && return :gr_brake_roche
    startswith(s, "GR BRAKE3 T M1 M2 K1 K2") && return :gr_brake3
    startswith(s, "GR CAPTURE T[NB]") && return :gw_capture
    startswith(s, "GR FLYBY   T[NB]") && return :gw_flyby
    startswith(s, "BINARY ESCAPE  IP,N,NCORR =") && return :binary_escape
    startswith(s, "HIARCH ESCAPE    NM =") && return :hierarchy_escape
    startswith(s, "QUAD ESCAPE    KS =") && return :quad_escape
    startswith(s, "MDEND KSREG TIME[NB]") && return :ks_end_mass_loss
    startswith(s, "MDNEW KSREG TIME[NB]") && return :ks_new_mass_loss
    startswith(s, "BINARY KICK:  TIME[Myr]") && return :binary_kick
    startswith(s, "NEW CIRC T DT") && return :circ_new
    startswith(s, "END CIRC T DT") && return :circ_end
    startswith(s, "CIRCULARIZED T DT") && return :circularized
    startswith(s, "INACTIIVE PHASE T DT") && return :circularized
    startswith(s, "NEW CAPTURE    T NM") && return :capture
    startswith(s, "NEW ROCHE    TIME[NB]") && return :roche_new
    startswith(s, "ROCHE COAL   TIME[NB]") && return :roche_coal
    startswith(s, "END ROCHE    TIME[NB]") && return :roche_end
    return :none
end

"""
    _parse_evolution_line(kind, lines, i, interval) -> (record, n_consumed)

Parse `lines[i]`, a line of `kind` (from `_evolution_line_kind`), into a
[`GWEvent`](@ref), [`BinaryEscapeEvent`](@ref), [`MassLossEvent`](@ref),
[`TidalEvent`](@ref) or [`RocheEvent`](@ref) with `line = i` and the given
adjustment `interval`. Every line of the family is one physical record, so
`n_consumed == 1`; the escape kinds look ahead for the time of the following
` ESCAPE T,NESC =` summary without consuming it. A line that does not match
its FORMAT (literal, field content, or a list-directed line with too few
values) throws `ArgumentError`, as does an unknown `kind`.
"""
function _parse_evolution_line(
    kind::Symbol,
    lines::AbstractVector{<:AbstractString},
    i::Int,
    interval::Int,
)::Tuple{_EvolutionEventRecord,Int}
    checkbounds(Bool, lines, i) ||
        throw(ArgumentError("Line index $i outside the $(length(lines)) lines given"))
    line = lines[i]
    record = if kind in _EVO_GW_KINDS
        _evo_parse_gw(kind, line, i, interval)
    elseif kind in _EVO_ESCAPE_KINDS
        _evo_parse_escape(kind, lines, i, interval)
    elseif kind in _EVO_MASS_LOSS_KINDS
        _evo_parse_mass_loss(kind, line, i, interval)
    elseif kind in _EVO_TIDAL_KINDS
        _evo_parse_tidal(kind, line, i, interval)
    elseif kind in _EVO_ROCHE_KINDS
        _evo_parse_roche(kind, line, i, interval)
    else
        throw(ArgumentError("Unknown evolution event kind $kind"))
    end
    return (record, 1)
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""Name or counter of value `k`; 0 when empty or overflowed."""
_evo_int(v::Vector{Float64}, k::Int)::Int = _fixed_int(v[k], 0)

"""Stellar type of value `k`; `KSTAR_ABSENT` when empty or overflowed."""
_evo_kstar(v::Vector{Float64}, k::Int)::Int = _fixed_int(v[k], KSTAR_ABSENT)

"""Values `k`, `k + 1`, `k + 2` as a Cartesian triple."""
_evo_xyz(v::Vector{Float64}, k::Int)::NTuple{3,Float64} = (v[k], v[k + 1], v[k + 2])

"""
The `n` reals that follow the text `label` in a list-directed line `s`;
`ArgumentError` when the label is absent or fewer than `n` values follow.
"""
function _evo_list_directed(s::AbstractString, label::AbstractString, n::Int)::Vector{Float64}
    r = findfirst(label, s)
    isnothing(r) && throw(ArgumentError("Label $(repr(label)) absent from list-directed line"))
    tokens = split(SubString(s, nextind(s, last(r))))
    length(tokens) >= n ||
        throw(ArgumentError("$(length(tokens)) of $n values follow $(repr(label))"))
    return [_parse_fortran_real(tokens[k]) for k in 1:n]
end

"""
Logical `L1` field at character `column` of `line`: `true` for `T`, `false`
for `F`, `ArgumentError` otherwise.
"""
function _evo_logical(line::AbstractString, column::Int)::Bool
    chars = collect(line)
    column <= length(chars) || throw(ArgumentError("Logical field at column $column absent"))
    c = chars[column]
    c == 'T' && return true
    c == 'F' && return false
    throw(ArgumentError("Logical field at column $column is $(repr(c)), not T or F"))
end

"""Label of the escape summary line written by `escape.F` after the escapers of a call."""
const _EVO_ESCAPE_SUMMARY_LABEL = "ESCAPE T,NESC ="

"""
TTOT of the first ` ESCAPE T,NESC =` line after line `i`, read as the first
token after the label; `NaN` when an `ADJUST:` line or the end of `lines` comes
first, or when the token is not a real.
"""
function _evo_escape_summary_time(lines::AbstractVector{<:AbstractString}, i::Int)::Float64
    for j in (i + 1):lastindex(lines)
        s = lstrip(lines[j])
        startswith(s, "ADJUST:") && return NaN
        startswith(s, _EVO_ESCAPE_SUMMARY_LABEL) || continue
        tokens = split(SubString(s, ncodeunits(_EVO_ESCAPE_SUMMARY_LABEL) + 1))
        isempty(tokens) && return NaN
        return something(tryparse(Float64, replace(tokens[1], r"[dD]" => "E")), NaN)
    end
    return NaN
end

"""
Fields read the same way from the GR lines that share `_evo_gw_core_format`:
`o` is the number of values before TTOT (2 for RES, KMIN, else 0). Names and
types are left to the caller (their order differs in `mdot.F`).
"""
function _evo_gw_core_fields(v::Vector{Float64}, o::Int)
    return (;
        time_nb = v[o + 1],
        step_nb = v[o + 3],
        npert = _evo_int(v, o + 6),
        m1_msun = v[o + 13],
        m2_msun = v[o + 14],
        ecc = v[o + 15],
        semi_nb = v[o + 16],
        qperi_nb = v[o + 17],
        rschw_nb = v[o + 18],
        h_nb = v[o + 19],
        gamma = v[o + 20],
        a_ein = v[o + 21],
        t_orb_nb = v[o + 22],
    )
end

# ---------------------------------------------------------------------------
# Parsers per record type
# ---------------------------------------------------------------------------

"""Build the [`GWEvent`](@ref) of a line of a gravitational-radiation `kind`."""
function _evo_parse_gw(kind::Symbol, line::AbstractString, i::Int, interval::Int)::GWEvent
    base = (; line = i, interval = interval, kind = kind)
    if kind === :gw_decay_step
        r = read_fixed(line, _EVO_GR_KSTIDE_FORMAT)
        v = r.values
        # v[12] is the temporary KSTAR(ICM) = -25, not the c.m. type.
        return GWEvent(;
            base...,
            _evo_gw_core_fields(v, 0)...,
            dt_gw_nb = v[2],
            index_cm = _evo_int(v, 4),
            ipair = _evo_int(v, 5),
            name1 = _evo_int(v, 7),
            name2 = _evo_int(v, 8),
            name_cm = _evo_int(v, 9),
            kstar1 = _evo_kstar(v, 10),
            kstar2 = _evo_kstar(v, 11),
            t_gr_nb = v[23],
            da_nb = v[24],
            de = v[25],
            d_energy_nb = v[26],
            dh_nb = v[27],
            qperi_old_nb = v[28],
            overflow = r.overflow,
        )
    elseif kind === :gw_ksint || kind === :gw_tides
        fmt = kind === :gw_ksint ? _EVO_GR_KSINT2_FORMAT : _EVO_GR_TIDES_FORMAT
        r = read_fixed(line, fmt)
        v = r.values
        # v[14] is the temporary KSTAR(ICM) = -25, not the c.m. type.
        return GWEvent(;
            base...,
            _evo_gw_core_fields(v, 2)...,
            res = _evo_int(v, 1),
            kmin = _evo_int(v, 2),
            dt_gw_nb = v[4],
            index_cm = _evo_int(v, 6),
            ipair = _evo_int(v, 7),
            name1 = _evo_int(v, 9),
            name2 = _evo_int(v, 10),
            name_cm = _evo_int(v, 11),
            kstar1 = _evo_kstar(v, 12),
            kstar2 = _evo_kstar(v, 13),
            t_gr_nb = v[25],
            overflow = r.overflow,
        )
    elseif kind === :gw_ksint_coll
        # PRINT *: QPERI, SEPMRG, CLIGHT, RSCHW, ECC, SEMI; no time, no names.
        x = _evo_list_directed(line, "RSCH ECC SEMI", 6)
        return GWEvent(;
            base...,
            time_nb = NaN,
            qperi_nb = x[1],
            sep_merge_nb = x[2],
            c_light_nb = x[3],
            rschw_nb = x[4],
            ecc = x[5],
            semi_nb = x[6],
        )
    elseif kind === :gw_unpert_coll
        # PRINT * on every MPI rank: QPERI, CLIGHT, RSCHW; no time, no names.
        x = _evo_list_directed(line, "QPERI CLIGHT RSCH", 3)
        return GWEvent(;
            base...,
            time_nb = NaN,
            qperi_nb = x[1],
            c_light_nb = x[2],
            rschw_nb = x[3],
        )
    elseif kind === :gw_recoil
        r = read_fixed(line, _EVO_GW_RECOIL_FORMAT)
        v = r.values
        return GWEvent(;
            base...,
            time_nb = v[1],
            index1 = _evo_int(v, 2),
            index2 = _evo_int(v, 3),
            name1 = _evo_int(v, 4),
            name2 = _evo_int(v, 5),
            kstar1 = _evo_kstar(v, 6),
            kstar2 = _evo_kstar(v, 7),
            icase = _evo_int(v, 8),
            time_myr = v[9],
            m01_msun = v[10],
            m02_msun = v[11],
            m1_msun = v[12],
            m2_msun = v[13],
            spin1 = v[14],
            spin2 = v[15],
            q_ratio = v[16],
            j_spin1_nb = v[17],
            j_spin2_nb = v[18],
            spin_factor = v[19],
            v_kick_kms = v[20],
            v_kick_xyz_kms = _evo_xyz(v, 21),
            v_esc_kms = v[24],
            v_dis_kms = v[25],
            r_over_semi = v[26],
            eb_nb = v[27],
            r_nb = v[28],
            x1_pc = _evo_xyz(v, 29),
            v1_kms = _evo_xyz(v, 32),
            x2_pc = _evo_xyz(v, 35),
            v2_kms = _evo_xyz(v, 38),
            overflow = r.overflow,
        )
    elseif kind === :gr_brake
        r = read_fixed(line, _EVO_GR_BRAKE_FORMAT)
        v = r.values
        # Names and types are printed as JX(2), JX(1), the masses as JX(1), JX(2):
        # reorder to the masses. IPAIR (v[5]) is a stale local variable.
        return GWEvent(;
            base...,
            _evo_gw_core_fields(v, 0)...,
            dt_gr_yr = v[2],
            index_star = _evo_int(v, 4),
            name1 = _evo_int(v, 8),
            name2 = _evo_int(v, 7),
            name_cm = _evo_int(v, 9),
            kstar1 = _evo_kstar(v, 11),
            kstar2 = _evo_kstar(v, 10),
            kstar_cm = _evo_kstar(v, 12),
            t_gr_yr = v[23],
            j_orb = v[24],
            dj_dt = v[25],
            dt_bse_yr = v[26],
            de_dt_per_yr = v[27],
            overflow = r.overflow,
        )
    elseif kind === :gr_brake_roche
        r = read_fixed(line, _EVO_GR_BRAKE_ROCHE_FORMAT)
        v = r.values
        # DELET1 (v[27]) read as the change of e over DTM0 [UNSURE].
        return GWEvent(;
            base...,
            _evo_gw_core_fields(v, 0)...,
            dt_gw_nb = v[2],
            index_cm = _evo_int(v, 4),
            ipair = _evo_int(v, 5),
            name1 = _evo_int(v, 7),
            name2 = _evo_int(v, 8),
            name_cm = _evo_int(v, 9),
            kstar1 = _evo_kstar(v, 10),
            kstar2 = _evo_kstar(v, 11),
            kstar_cm = _evo_kstar(v, 12),
            t_gr_steps = v[23],
            j_orb = v[24],
            dj = v[25],
            dt_bse_yr = v[26],
            de = v[27],
            overflow = r.overflow,
        )
    elseif kind === :gr_brake3
        r = read_fixed(line, _EVO_GR_BRAKE3_FORMAT)
        v = r.values
        return GWEvent(;
            base...,
            time_nb = v[1],
            m1_msun = v[2],
            m2_msun = v[3],
            kstar1 = _evo_kstar(v, 4),
            kstar2 = _evo_kstar(v, 5),
            semi_rsun = v[6],
            ecc = v[7],
            j_orb = v[8],
            dj_dt = v[9],
            dt_gr_yr = v[10],
            overflow = r.overflow,
        )
    elseif kind === :gw_capture || kind === :gw_flyby
        r = read_fixed(line, kind === :gw_capture ? _EVO_GR_CAPTURE_FORMAT : _EVO_GR_FLYBY_FORMAT)
        v = r.values
        return GWEvent(;
            base...,
            time_nb = v[1],
            name1 = _evo_int(v, 2),
            name2 = _evo_int(v, 3),
            kstar1 = _evo_kstar(v, 4),
            kstar2 = _evo_kstar(v, 5),
            m1_msun = v[6],
            m2_msun = v[7],
            v_inf_kms = v[8],
            ecc_minus_1 = v[9],
            qperi_nb = v[10],
            r_capture_nb = v[11],
            d_energy_nb = v[12],
            h_nb = v[13],
            h_after_nb = v[14],
            one_minus_e_after = v[15],
            npert = _evo_int(v, 16),
            overflow = r.overflow,
        )
    end
    throw(ArgumentError("Unknown gravitational-radiation event kind $kind"))
end

"""Build the [`BinaryEscapeEvent`](@ref) of `lines[i]`, timed by the following summary."""
function _evo_parse_escape(
    kind::Symbol,
    lines::AbstractVector{<:AbstractString},
    i::Int,
    interval::Int,
)::BinaryEscapeEvent
    line = lines[i]
    if kind === :binary_escape
        r = read_fixed(line, _EVO_BINARY_ESCAPE_FORMAT)
        v = r.values
        name_cm = _evo_int(v, 6)
        kstar_cm = _evo_kstar(v, 9)
        # KSTARI carries +30 for a merged hierarchy (NAME(c.m.) < 0).
        name_cm < 0 && kstar_cm != KSTAR_ABSENT && (kstar_cm -= 30)
        return BinaryEscapeEvent(;
            line = i,
            interval = interval,
            time_nb = _evo_escape_summary_time(lines, i),
            kind = kind,
            ipair = _evo_int(v, 1),
            n_after = _evo_int(v, 2),
            ncorr = _evo_int(v, 3),
            name1 = _evo_int(v, 4),
            name2 = _evo_int(v, 5),
            name_cm = name_cm,
            kstar1 = _evo_kstar(v, 7),
            kstar2 = _evo_kstar(v, 8),
            kstar_cm = kstar_cm,
            ls_flag = _evo_int(v, 10),
            m1_msun = v[11],
            m2_msun = v[12],
            semi_nb = v[13],
            ecc = v[14],
            eb_nb = v[15],
            period_d = v[16],
            rstar_over_peri = v[17],
            # Printed under "RI,VI" as VI, RI.
            v_nb = v[18],
            r_nb = v[19],
            e_cm_nb = v[20],
            radius1_rsun = v[21],
            radius2_rsun = v[22],
            overflow = r.overflow,
        )
    elseif kind === :hierarchy_escape
        r = read_fixed(line, _EVO_HIARCH_ESCAPE_FORMAT)
        v = r.values
        return BinaryEscapeEvent(;
            line = i,
            interval = interval,
            time_nb = _evo_escape_summary_time(lines, i),
            kind = kind,
            name1 = _evo_int(v, 1),
            name2 = _evo_int(v, 2),
            kstar1 = _evo_kstar(v, 3),
            kstar2 = _evo_kstar(v, 4),
            kstar_cm = _evo_kstar(v, 5),
            m1_msun = v[6],
            m2_msun = v[7],
            eb_nb = v[8],
            ecc = v[9],
            pmin_over_pcrit = v[10],
            semi_nb = v[11],
            period_d = v[12],
            overflow = r.overflow,
        )
    elseif kind === :quad_escape
        r = read_fixed(line, _EVO_QUAD_ESCAPE_FORMAT)
        v = r.values
        # NAME1, NAME2 (v[2:3]) are not assigned on this branch and KSTAR(JM)
        # (v[5]) uses a merger-table index as particle index [UNSURE]: absent.
        return BinaryEscapeEvent(;
            line = i,
            interval = interval,
            time_nb = _evo_escape_summary_time(lines, i),
            kind = kind,
            ipair = _evo_int(v, 1),
            ls_flag = _evo_int(v, 4),
            m1_msun = v[6],
            m2_msun = v[7],
            eb_nb = v[8],
            ecc = v[9],
            semi_nb = v[10],
            period_d = v[11],
            overflow = r.overflow,
        )
    end
    throw(ArgumentError("Unknown escape event kind $kind"))
end

"""Build the [`MassLossEvent`](@ref) of a line of a mass-loss `kind`."""
function _evo_parse_mass_loss(
    kind::Symbol,
    line::AbstractString,
    i::Int,
    interval::Int,
)::MassLossEvent
    if kind === :ks_end_mass_loss
        r = read_fixed(line, _EVO_MDEND_FORMAT)
        v = r.values
        return MassLossEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            name1 = _evo_int(v, 2),
            name2 = _evo_int(v, 3),
            name_cm = _evo_int(v, 4),
            kstar1 = _evo_kstar(v, 5),
            kstar2 = _evo_kstar(v, 6),
            kstar_cm = _evo_kstar(v, 7),
            ipair = _evo_int(v, 8),
            dtau_nb = v[9],
            m1_nb = v[10],
            m2_nb = v[11],
            r12_nb = v[12],
            ecc = v[13],
            semi_nb = v[14],
            eb_nb = v[15],
            period_d = v[16],
            h_nb = v[17],
            gamma = v[18],
            step_cm_nb = v[19],
            npert = _evo_int(v, 20),
            n_neighbours = _evo_int(v, 21),
            m1_msun = v[22],
            m2_msun = v[23],
            radius1_rsun = v[24],
            radius2_rsun = v[25],
            # Labelled RAD1,2,S[*]: the third value is the separation.
            r12_rsun = v[26],
            r_nb = v[27],
            v_nb = v[28],
            dm_nb = v[29],
            dmr = v[30],
            overflow = r.overflow,
        )
    elseif kind === :ks_new_mass_loss
        r = read_fixed(line, _EVO_MDNEW_FORMAT)
        v = r.values
        # Only time, names, types, masses, radii and DM, DMR are current; the
        # orbit, IPAIR, DTAU, R12, H, GAMMA, P, STEP, NPERT, NB, RI, VI and the
        # S name/type are stale or refer to the evolving star.
        return MassLossEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            name1 = _evo_int(v, 2),
            name2 = _evo_int(v, 3),
            kstar1 = _evo_kstar(v, 5),
            kstar2 = _evo_kstar(v, 6),
            m1_nb = v[10],
            m2_nb = v[11],
            m1_msun = v[22],
            m2_msun = v[23],
            radius1_rsun = v[24],
            radius2_rsun = v[25],
            dm_nb = v[29],
            dmr = v[30],
            overflow = r.overflow,
        )
    elseif kind === :binary_kick
        r = read_fixed(line, _EVO_BINARY_KICK_FORMAT)
        v = r.values
        return MassLossEvent(;
            line = i,
            interval = interval,
            time_nb = NaN,
            kind = kind,
            time_myr = v[1],
            name1 = _evo_int(v, 2),
            name2 = _evo_int(v, 3),
            kstar1 = _evo_kstar(v, 4),
            kstar2 = _evo_kstar(v, 5),
            spin1 = v[6],
            spin2 = v[7],
            q_ratio = v[8],
            j_spin1_nb = v[9],
            j_spin2_nb = v[10],
            spin_factor = v[11],
            m1_msun = v[12],
            m2_msun = v[13],
            v_esc_kms = v[14],
            v_dis_kms = v[15],
            r_over_semi = v[16],
            eb_nb = v[17],
            # Labelled R12[au]: R(IPAIR) in N-body units.
            r12_nb = v[18],
            x1_pc = _evo_xyz(v, 19),
            v1_kms = _evo_xyz(v, 22),
            x2_pc = _evo_xyz(v, 25),
            v2_kms = _evo_xyz(v, 28),
            overflow = r.overflow,
        )
    end
    throw(ArgumentError("Unknown mass-loss event kind $kind"))
end

"""Build the [`TidalEvent`](@ref) of a line of a tidal `kind`."""
function _evo_parse_tidal(kind::Symbol, line::AbstractString, i::Int, interval::Int)::TidalEvent
    if kind === :capture
        r = read_fixed(line, _EVO_CAPTURE_FORMAT)
        v = r.values
        return TidalEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            name1 = _evo_int(v, 2),
            name2 = _evo_int(v, 3),
            ecc0 = v[4],
            ecc = v[5],
            qperi_over_rstar = v[6],
            semi_nb = v[7],
            overflow = r.overflow,
        )
    end
    s = lstrip(line)
    inactive = kind === :circularized && startswith(s, "INACTIIVE PHASE T DT")
    r = if kind === :circ_new
        read_fixed(line, _EVO_NEW_CIRC_FORMAT)
    elseif kind === :circ_end
        read_fixed(line, _EVO_END_CIRC_FORMAT)
    elseif kind === :circularized && inactive
        start = something(findfirst(!isspace, collect(line)), 1)
        read_fixed(line, _EVO_INACTIVE_PHASE_FORMAT; start = start)
    elseif kind === :circularized
        read_fixed(line, _EVO_CIRCULARIZED_FORMAT)
    else
        throw(ArgumentError("Unknown tidal event kind $kind"))
    end
    v = r.values
    # PD: after the step (from a) in NEW/END CIRC, before it (from a0) in
    # CIRCULARIZED, unassigned on the INACTIIVE PHASE path [UNSURE].
    period_after = kind === :circularized ? NaN : v[22]
    period_before = kind === :circularized && !inactive ? v[22] : NaN
    return TidalEvent(;
        line = i,
        interval = interval,
        time_nb = v[1],
        kind = kind,
        inactive_phase = inactive,
        step_nb = v[2],
        index_cm = _evo_int(v, 3),
        ipair = _evo_int(v, 4),
        npert = _evo_int(v, 5),
        name1 = _evo_int(v, 6),
        name2 = _evo_int(v, 7),
        name_cm = _evo_int(v, 8),
        kstar1 = _evo_kstar(v, 9),
        kstar2 = _evo_kstar(v, 10),
        kstar_cm = _evo_kstar(v, 11),
        m1_msun = v[12],
        m2_msun = v[13],
        ecc0 = v[14],
        ecc = v[15],
        ecc_max = v[16],
        semi0_nb = v[17],
        semi_nb = v[18],
        qperi_nb = v[19],
        h0_nb = v[20],
        gamma = v[21],
        period_d = period_after,
        period0_d = period_before,
        radius1_rsun = v[23],
        radius2_rsun = v[24],
        r12_rsun = v[25],
        r_coll_rsun = v[26],
        r_nb = v[27],
        v_nb = v[28],
        overflow = r.overflow,
    )
end

"""
Fields of the `roche.f` body block (`_evo_roche_body_format`) whose first value,
e, is value `m + 1`.
"""
function _evo_roche_body_fields(v::Vector{Float64}, m::Int)
    return (;
        ecc = v[m + 1],
        semi_nb = v[m + 2],
        eb_nb = v[m + 3],
        period_d = v[m + 4],
        h_nb = v[m + 5],
        gamma = v[m + 6],
        step_cm_nb = v[m + 7],
        npert = _evo_int(v, m + 8),
        n_neighbours = _evo_int(v, m + 9),
        m01_msun = v[m + 10],
        m02_msun = v[m + 11],
        m1_msun = v[m + 12],
        m2_msun = v[m + 13],
        mc1_msun = v[m + 14],
        mc2_msun = v[m + 15],
        r_nb = v[m + 16],
        v_nb = v[m + 17],
        # Labelled SEP[*]: R(IPAIR)*SU, the instantaneous separation.
        r12_rsun = v[m + 18],
        radius1_rsun = v[m + 19],
        radius2_rsun = v[m + 20],
        roche_radius1_rsun = v[m + 21],
        roche_radius2_rsun = v[m + 22],
    )
end

"""Build the [`RocheEvent`](@ref) of a line of a Roche `kind`."""
function _evo_parse_roche(kind::Symbol, line::AbstractString, i::Int, interval::Int)::RocheEvent
    if kind === :roche_new
        r = read_fixed(line, _EVO_NEW_ROCHE_FORMAT)
        v = r.values
        return RocheEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            name1 = _evo_int(v, 2),
            name2 = _evo_int(v, 3),
            name_cm = _evo_int(v, 4),
            kstar1 = _evo_kstar(v, 5),
            kstar2 = _evo_kstar(v, 6),
            kstar_cm = _evo_kstar(v, 7),
            ipair = _evo_int(v, 8),
            dtau_nb = v[9],
            time_myr = v[10],
            age_myr = v[11],
            m1_nb = v[12],
            m2_nb = v[13],
            spin1_nb = v[14],
            spin2_nb = v[15],
            ospin1_nb = v[16],
            ospin2_nb = v[17],
            spin_factor = v[18],
            j_orb = v[19],
            omega_orb_per_yr = v[20],
            r12_nb = v[21],
            _evo_roche_body_fields(v, 21)...,
            # Labelled SEMI[*]: SEP, the semi-major axis.
            semi_rsun = v[44],
            rl1_rsun = v[45],
            overflow = r.overflow,
        )
    end
    coal = kind === :roche_coal
    coal || kind === :roche_end || throw(ArgumentError("Unknown Roche event kind $kind"))
    r = read_fixed(line, coal ? _EVO_ROCHE_COAL_FORMAT : _EVO_END_ROCHE_FORMAT)
    v = r.values
    body = _evo_roche_body_fields(v, 13)
    # EB, RI, VI are formed on the NEW ROCHE path only [UNSURE]; ROCHE COAL
    # prints TK = 0.
    body = merge(body, (; eb_nb = NaN, r_nb = NaN, v_nb = NaN))
    coal && (body = merge(body, (; period_d = NaN)))
    return RocheEvent(;
        line = i,
        interval = interval,
        time_nb = v[1],
        kind = kind,
        name1 = _evo_int(v, 2),
        name2 = _evo_int(v, 3),
        name_cm = _evo_int(v, 4),
        kstar1 = _evo_kstar(v, 5),
        kstar2 = _evo_kstar(v, 6),
        kstar_cm = _evo_kstar(v, 7),
        ipair = _evo_int(v, 8),
        dtau_nb = v[9],
        time_myr = v[10],
        m1_nb = v[11],
        m2_nb = v[12],
        r12_nb = v[13],
        body...,
        dm1_msun = v[36],
        dm2_msun = v[37],
        dt_myr = v[38],
        age_myr = coal ? v[39] : NaN,
        coalescence = coal ? _evo_logical(line, _EVO_ROCHE_COAL_L1_COLUMN) : false,
        ixxx = coal ? _evo_int(v, 40) : 0,
        overflow = r.overflow,
    )
end
