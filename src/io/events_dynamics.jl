# =============================================================================
# KS, few-body, hierarchy and exchange event lines of the engine log (out1000)
# =============================================================================
#
# Each kind, the label of its line (leading blanks removed) and the writer:
#
#   kind                       label                                     writer
#   :ks_new                    NEW KSREG   TIME[NB] … RCM,VCM[NB]=       ksinit.F
#   :ks_new_primordial         NEW KSREG   TIME[NB] … RI,VI[NB]=         kspinit.f
#   :ks_end                    END KSREG   TIME[NB]                      ksterm.F
#   :ks_new_compact            NS/BH BINARY T[NB]                        ksinit.F
#   :triple_new                NEW TRIPLE  TIME                          impact.f
#   :quad_new                  NEW QUAD    TIME                          impact.f
#   :chain_new                 NEW CHAIN   TIME                          impact.f
#   :chain_init                NEW CHAIN   NCH NP ENERGY …               chain.f
#   :chain_end                 END CHAIN NCH I1-4                        recoil.f
#   :chain_end_collision       END CHAIN:   RIJ RCOLL QPERI GPERT        chain.f
#   :chterm_names              CHTERM: NAMEC                             chterm.f
#   :chterm                    CHTERM: TIME[NB]                          chterm.f
#   :chterm_second_binary      CHTERM  SECOND BINARY(ICM2)               chterm.f
#   :chain_bb                  CHAIN B-B    NAME(I1-I4)                  impact.f
#   :chain_escape_single       CHAIN ESCAPE: TT T TC = … IESC NM         chmod_chain.f
#   :chain_escape_binary       CHAIN ESCAPE: TT T TC = … IESC JESC       chmod_chain.f
#   :chain_binary              CHAIN BINARY  NAM =                       recoil.f
#   :hierarchy_new             NEW MERGER  TIME                          impact.f
#   :hierarchy_new_merge2      NEW MERGE2  TIME                          impact.f
#   :hierarchy_hi_merge        HI MERGE    NAM E E1 A A1 PM PC           impact.f
#   :hierarchy_end             END MERGER:  TIME                         reset.f
#   :hierarchy_new_quintuplet  NEW QUINTUPLET   T MT …                   merge2.f
#   :hierarchy_new_sextuplet   NEW SEXTUPLET    T MT …                   merge2.f
#   :hierarchy_new_quartet     NEW QUARTET    T MT …                     merge2.f
#   :hierarchy_new_quintup2    NEW QUINTUP2    T MT …                    merge2.f
#   :hierarchy_new_hitrip      NEW HITRIP    T MT …                      merge2.f
#   :hierarchy_end_quintuplet  END QUINTUPLET   T MT …                   reset2.f
#   :hierarchy_end_sextuplet   END SEXTUPLET    T MT …                   reset2.f
#   :hierarchy_end_quartet     END QUARTET    T MT …                     reset2.f
#   :hierarchy_end_quintup2    END QUINTUP2    T MT …                    reset2.f
#   :hierarchy_end_hitrip      END HITRIP    T MT …                      reset2.f
#   :hierarchy_end_outer       END OUTER MERGER                          reset2.f
#   :hierarchy_end_merge2      END MERGE2                                reset2.f
#   :exchange                  EXCHANGE    T NAM E0 E A0 A EB0 EB        recoil.f
#   :end_run_deferred          END RUN DEFERRED   TTOT =                 adjust.F
#   :chain_tiny_steps          ENFORCED CHAIN TINY STEPS:                chain.f
#
# Every record is read through the FORMAT of its writer (fixed widths): names
# fill their I10/I6 fields and glue to the labels, negative 1PEw.d values with
# w = d + 7 glue to their neighbours. The per-step chain diagnostics
# (`CHAIN: NSTEP …`, `CHAIN T,TC[NB] …`) are not events; the continuation
# record of the binary CHAIN ESCAPE is consumed with its first line.

"""
All event kinds returned by `_dynamics_line_kind`: KS regularisation
([`KSEvent`](@ref)), chain and few-body regularisation ([`FewBodyEvent`](@ref)),
hierarchical mergers ([`HierarchyEvent`](@ref)), chain exchanges
([`ExchangeEvent`](@ref)) and termination diagnostics ([`TerminationEvent`](@ref)).
"""
const DYNAMICS_EVENT_KINDS = (
    :ks_new,
    :ks_new_primordial,
    :ks_end,
    :ks_new_compact,
    :triple_new,
    :quad_new,
    :chain_new,
    :chain_init,
    :chain_end,
    :chain_end_collision,
    :chterm_names,
    :chterm,
    :chterm_second_binary,
    :chain_bb,
    :chain_escape_single,
    :chain_escape_binary,
    :chain_binary,
    :hierarchy_new,
    :hierarchy_new_merge2,
    :hierarchy_hi_merge,
    :hierarchy_end,
    :hierarchy_new_quintuplet,
    :hierarchy_new_sextuplet,
    :hierarchy_new_quartet,
    :hierarchy_new_quintup2,
    :hierarchy_new_hitrip,
    :hierarchy_end_quintuplet,
    :hierarchy_end_sextuplet,
    :hierarchy_end_quartet,
    :hierarchy_end_quintup2,
    :hierarchy_end_hitrip,
    :hierarchy_end_outer,
    :hierarchy_end_merge2,
    :exchange,
    :end_run_deferred,
    :chain_tiny_steps,
)

const _DYN_NO_NAMES = (0, 0, 0, 0, 0, 0)
const _DYN_NO_INDICES = (0, 0, 0, 0)
const _DYN_NO_KSTARS = (KSTAR_ABSENT, KSTAR_ABSENT, KSTAR_ABSENT, KSTAR_ABSENT)
const _DYN_NAN3 = (NaN, NaN, NaN)
const _DYN_NAN4 = (NaN, NaN, NaN, NaN)
const _DYN_NAN5 = (NaN, NaN, NaN, NaN, NaN)
const _DYN_NAN8 = (NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)

"""Lines searched for the time stamp of a chain record that carries none."""
const _DYN_TIME_WINDOW = 12

# ---------------------------------------------------------------------------
# Record types
# ---------------------------------------------------------------------------

"""
    KSEvent

One KS-regularisation record of the engine log, one physical line (the empty
line the leading `/` of the FORMAT writes before it is not part of it).

| kind | writer, gate | label | values |
|---|---|---|---|
| `:ks_new` | `ksinit.F` FORMAT 60, `KZ(10) > 0` | `NEW KSREG   TIME[NB] … RCM,VCM[NB]=` | 40 |
| `:ks_new_primordial` | `kspinit.f` FORMAT 60, start with `KZ(8) > 0` | `NEW KSREG   TIME[NB] … RI,VI[NB]=` | 28 |
| `:ks_end` | `ksterm.F` FORMAT 15, `KZ(10) > 1` | `END KSREG   TIME[NB]` | 40 |
| `:ks_new_compact` | `ksinit.F` FORMAT 65, ungated | `NS/BH BINARY T[NB]` | 29 |

| field | meaning | unit |
|---|---|---|
| `line`, `interval`, `time_nb`, `kind` | line number, adjustment interval, TIME + TOFF, kind | —, —, NB, — |
| `name1`, `name2`, `name_cm` | component names, c.m. name (< 0 for a merged hierarchy) | — |
| `kstar1`, `kstar2`, `kstar_cm` | stellar types; c.m. type (0, −1 chaos, 9/10 circularisation, ≥ 11 Roche) | — |
| `ipair` | KS pair index | — |
| `dtau` | regularised KS step | KS time |
| `m1_nb`, `m2_nb` | component masses | NB |
| `r12_nb` | separation | NB |
| `ecc`, `semi_nb` | eccentricity (> 1 hyperbolic), semi-major axis (< 0 hyperbolic) | —, NB |
| `eb_nb` | binding energy: μH (`ksinit.F`, `kspinit.f`), −½m₁m₂/a (`ksterm.F`); > 0 hyperbolic | NB |
| `period_d` | period 2πa√(abs(a)/M)·T*; < 0 hyperbolic | d |
| `h_nb`, `gamma` | specific binding energy, relative perturbation | NB, — |
| `step_cm_nb` | c.m. time step | NB |
| `n_perturbers`, `n_neighbours_cm` | KS perturbers (label NPERT), c.m. neighbours (label NB(ICM)) | — |
| `m1_msun`, `m2_msun` | component masses | M☉ |
| `radius1_rsun`, `radius2_rsun` | stellar radii | R☉ |
| `r12_rsun` | separation (FORMAT field 26, label `RAD1,2,S[*]`) | R☉ |
| `r_cm_nb`, `v_cm_nb` | c.m. distance from the density centre, c.m. speed | NB |
| `deb_nb` | `:ks_new_compact` only: abs(ΔEB/EB) against the previous NS/BH record | — |
| `r1_xy_nb`, `r1_lat`, `r1_azi` | component 1 position w.r.t. the density centre: cylindrical radius, latitude atan(z/R_xy), azimuth acos(x/R_xy) | NB, rad, rad |
| `v1_xy_nb`, `v1_lat`, `v1_azi` | component 1 velocity, same decomposition, absolute frame | NB, rad, rad |
| `r2_…`, `v2_…` | the same for component 2 | NB, rad, rad |
| `overflow` | a field was printed as asterisks (read as absent) | — |

Corrections and hazards: FORMAT field 26 is the separation in R☉, not a radius,
and is stored as `r12_rsun`. The 28- and 29-value variants carry no orientation
blocks (NaN); in `:ks_new_compact` `r_cm_nb` and `v_cm_nb` are ‖X(cm)‖ and
‖Ẋ(cm)‖, not relative to the density centre, and DEB is the dimensionless
relative change of EB (the `_nb` suffix of `deb_nb` notwithstanding). The
azimuths come from ACOS and lose the sign of y. END KSREG is written for every
termination (chain exchange, merger, collision, mass loss, escape), not only
disruptions; for a hierarchical outer pair NM1 is the inner c.m. and M1 its total
mass.
"""
Base.@kwdef struct KSEvent
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
    dtau::Float64 = NaN
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
    n_perturbers::Int = 0
    n_neighbours_cm::Int = 0
    m1_msun::Float64 = NaN
    m2_msun::Float64 = NaN
    radius1_rsun::Float64 = NaN
    radius2_rsun::Float64 = NaN
    r12_rsun::Float64 = NaN
    r_cm_nb::Float64 = NaN
    v_cm_nb::Float64 = NaN
    deb_nb::Float64 = NaN
    r1_xy_nb::Float64 = NaN
    r1_lat::Float64 = NaN
    r1_azi::Float64 = NaN
    v1_xy_nb::Float64 = NaN
    v1_lat::Float64 = NaN
    v1_azi::Float64 = NaN
    r2_xy_nb::Float64 = NaN
    r2_lat::Float64 = NaN
    r2_azi::Float64 = NaN
    v2_xy_nb::Float64 = NaN
    v2_lat::Float64 = NaN
    v2_azi::Float64 = NaN
    overflow::Bool = false
end

"""
    FewBodyEvent

One chain or few-body regularisation record of the engine log. The fields a
kind does not carry hold `NaN`, `0` or [`KSTAR_ABSENT`](@ref).

| kind | writer, gate | time source |
|---|---|---|
| `:triple_new`, `:quad_new`, `:chain_new` | `impact.f` FORMAT 20, `KZ(15) > 1` or `KZ(30) > 1` | own TTOT |
| `:chain_init` | `chain.f` FORMAT 15, `KZ(30) > 1` | following `CHAIN: NSTEP` line (F10.4) |
| `:chain_end` | `recoil.f` FORMAT 30, `KZ(30) > 1` | following `CHTERM: TIME[NB]` line (F10.4) |
| `:chain_end_collision` | `chain.f` FORMAT 80, `KZ(30) > 1`, QPERI < 1 | following `CHEND:` line (E17.10) |
| `:chterm_names` | `chterm.f` FORMAT 4, `KZ(30) > 1` | following `CHTERM: TIME[NB]` line |
| `:chterm` | `chterm.f` FORMAT 75, `KZ(30) > 1` | own TIME + TOFF (F10.4) |
| `:chterm_second_binary` | `chterm.f` FORMAT 90, `KZ(30) > 1` | preceding `CHTERM: TIME[NB]` line |
| `:chain_bb` | `impact.f` FORMAT 22, binary–binary chain start | preceding `NEW CHAIN   TIME` line (written by the same IMPACT call), else the following one |
| `:chain_escape_single` | `chmod_chain.f` FORMAT 35 | own TTOT |
| `:chain_escape_binary` | `chmod_chain.f` FORMAT 28 (two physical lines) | own TTOT |
| `:chain_binary` | `recoil.f` FORMAT 20, DB > 0.1, a > 0 | following `CHTERM: TIME[NB]` line |

The neighbouring stamp is searched within $(_DYN_TIME_WINDOW) lines; `time_nb`
is `NaN` when none is found.

Fields per kind (`impact` = the three FORMAT 20 kinds):

| field | meaning, unit | kinds |
|---|---|---|
| `n_members` | NCH (`chain_init`, `chain_end`); 3 for `triple_new` and for `chterm_names` of a three-body chain; 4 for `quad_new`, `chain_bb`; 0 otherwise | |
| `n_perturbers` | perturbers (NP, NPERT) | impact, `chain_init`, `chain_end` |
| `names` | impact: inner components, intruder, inner c.m.; `chain_init`: NAMEC(1…NCH); `chain_end`: closest pair I1, I2, second pair I3, I4; `chterm_names`: NAMEC(1…4); `chterm`: NAME(I3), NAME(I4), NAME(NTOT); `chterm_second_binary`: NAME(I3), NAME(I4); `chain_bb`: pair IPAIR, pair JPAIR; escapes: escaper(s); `chain_binary`: closest pair | all |
| `kstar_values` | impact: inner components, intruder, inner c.m.; `chain_bb`: c.m. of IPAIR, c.m. of JPAIR | impact, `chain_bb` |
| `indices` | impact: global I (inner c.m.), JCOMP (intruder), KS pair IPAIR; `chain_end`: chain-local I1…I4; `chterm`: global I3, I4, NTOT; `chterm_second_binary`: I3, I4; escapes: chain-local IESC, JESC | |
| `masses_nb` | impact: m₁, m₂, m₃, total; `chain_end`: M(I1…I4); `chterm_names`: M(1…4); `chain_escape_binary`: escaping components [NB] | |
| `masses_msun` | impact: m₁, m₂, m₃, total; `chterm_names`: M(1…4) [M☉] | |
| `radii_rsun` | impact: radii of the components and the intruder (0 for a c.m. intruder) [R☉] | impact |
| `separations_nb` | `chain_end`: RB, R13, R24; `chain_end_collision`: RIJ(1,2), RIJ(1,3), RIJ(2,3); `chterm_names`: RIJ of the sorted pairs 1-2, 1-3, 2-3, 2-4, 3-4 [NB] | |
| `steps_nb` | `chterm`: STEP(I3), STEP(I4), STEP(NTOT); `chterm_second_binary`: STEP(NTOT), STEPR(NTOT) [NB] | |
| `r12_nb`, `h_nb` | inner separation and specific energy (impact); R, H of the second binary (`chterm_second_binary`) [NB] | |
| `gamma` | perturbation of the second binary | `chterm_second_binary` |
| `ecc_in`, `semi_in_nb`, `eb_in_nb`, `period_in_d` | inner orbit (impact); IPAIR orbit (`chain_bb`: ECC0, SEMI0); final binary (`chain_binary`: E, A) | |
| `ecc_out`, `semi_out_nb`, `eb_out_nb`, `period_out_d` | relative orbit of inner c.m. and intruder (impact); relative c.m. orbit (`chain_bb`: ECC1, SEMI1) | |
| `semi_second_nb` | semi-major axis of JPAIR (SEMI2) [NB] | `chain_bb` |
| `g4` | tidal perturbation PERT4 (impact); perturbation by the second binary (`chain_binary`) | |
| `rij_nb`, `pmin_nb` | intruder distance and outer pericentre (impact, E8.1); c.m. separation and pericentre (`chain_bb`) [NB] | |
| `eb_ratio` | EB1/EB (impact, E8.1) | impact |
| `r12_rsun`, `rij_rsun` | inner and outer separation [R☉] | impact |
| `r_cm_nb`, `v_cm_nb` | inner c.m. distance from the density centre, c.m. speed [NB] | impact |
| `energy_nb`, `rsum_nb`, `rgrav_nb`, `tcr_nb`, `rmaxs_nb` | chain energy (F10.5), sum of separations, gravitational radius, crossing time, maximum unperturbed size [NB] | `chain_init` |
| `de_nb` | accumulated abs(ΔE) of the chain [NB] | `chain_end` |
| `duration_crossings` | chain duration in crossing times (TC) | `chain_end` |
| `n_steps`, `n_reregularisations` | chain steps, re-regularisations | `chain_end` |
| `gpert` | perturbation: GPERT (`chain_end`, `chain_end_collision`, escapes), GB (`chain_binary`) | |
| `db` | relative binary-energy change | `chain_end`, `chain_binary` |
| `ec` | ENERGY/CM(8) | `chain_end` |
| `rcoll_nb`, `qperi_nb` | collision distance, pericentre [NB] | `chain_end_collision` |
| `time_block_nb`, `time_chain_nb` | block time TIME, chain time TIMEC [NB] | escapes |
| `r_escape_nb`, `rdot_escape_nb` | escaper distance from the inner part (FAC-scaled), radial velocity [NB] | escapes |
| `v2_escape_nb` | 2·M(ICH)/RI, escape-speed² proxy [NB] | escapes |
| `rb_nb` | separation of the escaping binary [NB] | `chain_escape_binary` |
| `vinf_kms` | asymptotic speed (0 if bound) [km/s] | escapes |
| `eb_over_energy`, `eb1_over_energy` | EB/ENERGY, EB1/ENERGY of the final and the second binary | `chain_binary` |
| `e1_over_eb0`, `et_over_cm8` | (ENERGY − EB − EB1)/EB0, −m₁m₂/(2a·CM(8)) | `chain_binary` |
| `overflow` | a field was printed as asterisks (read as absent) | all |

Corrections and hazards. Impact periods: the engine computes the inner "period"
with the total mass and the outer one with the inner c.m. mass; they are stored
rescaled to the masses of their own orbits, P_in·√(M_tot/M_in) and
P_out·√(M_in/M_tot) with M_in = m₁ + m₂ (NaN when a mass is absent). For
`:quad_new` written by the regularisation branch, `semi_in_nb` is a(IPAIR) +
a(JPAIR) as printed and the corrected inner period refers to that sum. The label
`NEW QUAD    TIME` is also written by the quadruple-merger branch of `impact.f`
with the same layout; the two are told apart only by the record that follows
(`END QUAD` of the regularisation, `END MERGER` of the merger). Impact fields 3–5
are indices, not names. `chain_end` with NCH ≤ 3 repeats member I1 as I4 and
prints R24 = 10¹⁰: the fourth index, name and mass and R24 are stored absent.
`chterm_names` with RIJ(3-4) = 0 belongs to a three-body chain (I4 = I3): the
stale fourth name and masses and the 2-4, 3-4 distances are stored absent. The
binary CHAIN ESCAPE prints 15 items through 13 descriptors: the labels after
`NM1/2` are shifted by two values (masses under `RI VI`, RI, RDOT, GPERT, 2M/R
under `GPERT 2*M/R RB VF`) and RB, VINF form a second record under `TT T`; the
values are stored under their true meaning and both lines are consumed. The
F10.4 times have a resolution of 10⁻⁴ NB; E8.1 fields carry two digits.
"""
Base.@kwdef struct FewBodyEvent
    line::Int
    interval::Int
    time_nb::Float64
    kind::Symbol
    n_members::Int = 0
    n_perturbers::Int = 0
    names::NTuple{6,Int} = _DYN_NO_NAMES
    kstar_values::NTuple{4,Int} = _DYN_NO_KSTARS
    indices::NTuple{4,Int} = _DYN_NO_INDICES
    masses_nb::NTuple{4,Float64} = _DYN_NAN4
    masses_msun::NTuple{4,Float64} = _DYN_NAN4
    radii_rsun::NTuple{3,Float64} = _DYN_NAN3
    separations_nb::NTuple{5,Float64} = _DYN_NAN5
    steps_nb::NTuple{3,Float64} = _DYN_NAN3
    r12_nb::Float64 = NaN
    h_nb::Float64 = NaN
    gamma::Float64 = NaN
    ecc_in::Float64 = NaN
    semi_in_nb::Float64 = NaN
    eb_in_nb::Float64 = NaN
    period_in_d::Float64 = NaN
    ecc_out::Float64 = NaN
    semi_out_nb::Float64 = NaN
    eb_out_nb::Float64 = NaN
    period_out_d::Float64 = NaN
    semi_second_nb::Float64 = NaN
    g4::Float64 = NaN
    rij_nb::Float64 = NaN
    pmin_nb::Float64 = NaN
    eb_ratio::Float64 = NaN
    r12_rsun::Float64 = NaN
    rij_rsun::Float64 = NaN
    r_cm_nb::Float64 = NaN
    v_cm_nb::Float64 = NaN
    energy_nb::Float64 = NaN
    rsum_nb::Float64 = NaN
    rgrav_nb::Float64 = NaN
    tcr_nb::Float64 = NaN
    rmaxs_nb::Float64 = NaN
    de_nb::Float64 = NaN
    duration_crossings::Float64 = NaN
    n_steps::Int = 0
    n_reregularisations::Int = 0
    gpert::Float64 = NaN
    db::Float64 = NaN
    ec::Float64 = NaN
    rcoll_nb::Float64 = NaN
    qperi_nb::Float64 = NaN
    time_block_nb::Float64 = NaN
    time_chain_nb::Float64 = NaN
    r_escape_nb::Float64 = NaN
    rdot_escape_nb::Float64 = NaN
    v2_escape_nb::Float64 = NaN
    rb_nb::Float64 = NaN
    vinf_kms::Float64 = NaN
    eb_over_energy::Float64 = NaN
    eb1_over_energy::Float64 = NaN
    e1_over_eb0::Float64 = NaN
    et_over_cm8::Float64 = NaN
    overflow::Bool = false
end

"""
    HierarchyEvent

One hierarchical-merger record of the engine log.

| kind | writer, gate | time source |
|---|---|---|
| `:hierarchy_new`, `:hierarchy_new_merge2` | `impact.f` FORMAT 20 (NEW MERGER, NEW MERGE2), `KZ(15) > 1` | own TTOT |
| `:hierarchy_hi_merge` | `impact.f` FORMAT 60, two hierarchies merging | preceding `NEW MERGE2  TIME` line |
| `:hierarchy_end` | `reset.f` FORMAT 21 (END MERGER), `KZ(15) > 1` | own TTOT |
| `:hierarchy_new_quintuplet`, `…_sextuplet` | `merge2.f` FORMAT 5 | own (F10.2) |
| `:hierarchy_new_quartet`, `…_quintup2` | `merge2.f` FORMAT 7, 8 | own (F9.2, F10.2) |
| `:hierarchy_new_hitrip` | `merge2.f` FORMAT 10 | own (F9.2) |
| `:hierarchy_end_quintuplet`, `…_sextuplet` | `reset2.f` FORMAT 5 | own (F9.2) |
| `:hierarchy_end_quartet`, `…_quintup2` | `reset2.f` FORMAT 9, 11 | own (F9.2, F10.2) |
| `:hierarchy_end_hitrip` | `reset2.f` FORMAT 18 | own (F9.2) |
| `:hierarchy_end_outer` | `reset2.f` FORMAT 48, `KZ(15) > 1` | following `END MERGE2` line (F8.2) |
| `:hierarchy_end_merge2` | `reset2.f` FORMAT 65, `KZ(15) > 1` | own (F8.2) |

The neighbouring stamp is searched within $(_DYN_TIME_WINDOW) lines (`NaN` when
none is found). The quadruple merger of `impact.f` prints `NEW QUAD    TIME` and
is returned as `:quad_new` ([`FewBodyEvent`](@ref)).

Fields (`impact` = the two FORMAT 20 kinds, `merge2`/`reset2` = the FORMAT 5–18
kinds of those routines):

| field | meaning, unit | kinds |
|---|---|---|
| `merger_index` | merger table index IMERGE | `end_merge2` |
| `n_perturbers` | inner perturbers (NP) | impact |
| `n_neighbours` | neighbours NNB − 1 | `end_merge2` |
| `names` | impact: inner components, outer body, inner c.m.; `hi_merge`: NAME(I1), ghost, NAME(I1+1), NAME(J1), ghost, NAME(J1+1); `end`: restored inner components; merge2/reset2: the names as printed (NM1 NM2/NMG NM3, HITRIP four names) | |
| `kstar_values` | impact: inner components, outer body, inner c.m.; `end`: inner components | |
| `indices` | impact: global I, JCOMP, KS pair IPAIR; `end`: global ICOMP, JCOMP of the restored pair; `end_outer`: pair JPAIR in slot 3 | |
| `masses_nb` | impact: m₁, m₂, m₃, total; `end`: inner components; `end_outer`: outer-binary components; `end_merge2`: restored inner components [NB] | |
| `masses_msun` | impact: m₁, m₂, m₃, total; `end`: inner components; merge2/reset2: total ZM in slot 4 [M☉] | |
| `radii_rsun` | impact: components and outer body; `end`: inner components [R☉] | |
| `r12_nb`, `h_nb` | impact: inner separation and specific energy; HITRIP: R(IPAIR); END QUARTET/QUINTUP2: R1; `end_outer`: R, H of JPAIR; `end_merge2`: R1 [NB] | |
| `r34_nb` | R(JPAIR) [NB] | `new_hitrip` |
| `ecc_in`, `semi_in_nb`, `eb_in_nb`, `period_in_d` | impact: inner orbit; merge2 and `hi_merge`: e, a of KS pair IPAIR (E0, A0) | |
| `ecc_out`, `semi_out_nb`, `eb_out_nb`, `period_out_d` | impact: outer orbit; reset2: e, a of the dissolved merged pair (E1, A1); `hi_merge`: ECC1, SEMI1; `end_merge2`: A1 | |
| `semi_inner1_nb`, `semi_inner2_nb` | inner semi-major axes AI, AJ [NB] | HITRIP |
| `g4`, `rij_nb`, `pmin_nb`, `eb_ratio` | impact: PERT4, outer distance, outer pericentre, EB1/EB; reset2, `hi_merge`: PMIN | |
| `r12_rsun`, `rij_rsun`, `r_cm_nb`, `v_cm_nb` | impact: inner, outer separation [R☉]; c.m. distance from the density centre, speed [NB]; quintuplet/sextuplet: RI | |
| `pcrit_nb` | stability pericentre PCRIT [NB] | merge2, `hi_merge` |
| `gamma` | GAMMA(IPAIR) (G0), max γ (GX), GB, G2 | merge2, `end_outer`, `end_merge2` |
| `gamma_merged` | perturbation G1 of the merged pair before restoration | reset2, `end_merge2` |
| `eb_kinetic` | ½m₁m₂/a in units of the mean kinetic energy (positive) | quintuplet/sextuplet |
| `eb_fraction` | binding energy over BE(3) | `end_outer`, `end_merge2` |
| `energy_ratio` | merged-pair over outer-binary energy (E1/EB2), merged-pair over inner energy (E1/EB) | `end_outer`, `end_merge2` |
| `dpot_nb` | tidal potential change POT4 − POT3 [NB] | `end_outer` |
| `overflow` | a field was printed as asterisks (read as absent) | all |

Corrections and hazards. Impact periods are rescaled from the swapped masses
as in [`FewBodyEvent`](@ref). END MERGER is written after the restoration: IP,
NAME(ICM), KSTAR(ICM), R12, H, R34, H2 and SEPI/J belong to whatever pair now
occupies the old slot, and the JX block repeats the inner component 1; only the
time, ICOMP, JCOMP and the inner components' names, types, masses and radii are
stored, everything else stays absent (the orbit of the restored pair is in the
adjacent NEW KSREG line, the outer body in the END KSREG line of the merged
pair). merge2/reset2 names are I6 and overflow at 10⁶; ZM is F6.2.
"""
Base.@kwdef struct HierarchyEvent
    line::Int
    interval::Int
    time_nb::Float64
    kind::Symbol
    merger_index::Int = 0
    n_perturbers::Int = 0
    n_neighbours::Int = 0
    names::NTuple{6,Int} = _DYN_NO_NAMES
    kstar_values::NTuple{4,Int} = _DYN_NO_KSTARS
    indices::NTuple{4,Int} = _DYN_NO_INDICES
    masses_nb::NTuple{4,Float64} = _DYN_NAN4
    masses_msun::NTuple{4,Float64} = _DYN_NAN4
    radii_rsun::NTuple{3,Float64} = _DYN_NAN3
    r12_nb::Float64 = NaN
    h_nb::Float64 = NaN
    r34_nb::Float64 = NaN
    ecc_in::Float64 = NaN
    semi_in_nb::Float64 = NaN
    eb_in_nb::Float64 = NaN
    period_in_d::Float64 = NaN
    ecc_out::Float64 = NaN
    semi_out_nb::Float64 = NaN
    eb_out_nb::Float64 = NaN
    period_out_d::Float64 = NaN
    semi_inner1_nb::Float64 = NaN
    semi_inner2_nb::Float64 = NaN
    g4::Float64 = NaN
    rij_nb::Float64 = NaN
    pmin_nb::Float64 = NaN
    eb_ratio::Float64 = NaN
    r12_rsun::Float64 = NaN
    rij_rsun::Float64 = NaN
    r_cm_nb::Float64 = NaN
    v_cm_nb::Float64 = NaN
    pcrit_nb::Float64 = NaN
    gamma::Float64 = NaN
    gamma_merged::Float64 = NaN
    eb_kinetic::Float64 = NaN
    eb_fraction::Float64 = NaN
    energy_ratio::Float64 = NaN
    dpot_nb::Float64 = NaN
    overflow::Bool = false
end

"""
    ExchangeEvent

An exchange in a chain (`recoil.f` FORMAT 15, kind `:exchange`), written at the
chain termination when the closest pair differs from both initial pairs and is
bound.

| field | meaning | unit |
|---|---|---|
| `line`, `interval`, `kind` | line number, adjustment interval, `:exchange` | — |
| `time_nb` | T0S(ISUB) + TIMEC (F9.2) | NB |
| `name1_initial`, `name2_initial` | closest pair at chain start | — |
| `name1_final`, `name2_final` | closest pair at termination | — |
| `e0`, `e` | initial, final eccentricity (F7.3) | — |
| `a0_nb`, `a_nb` | initial (more bound initial pair), final semi-major axis (E9.1) | NB |
| `eb0_nb`, `eb_nb` | initial, final binding energy (E10.1) | NB |
| `overflow` | a field was printed as asterisks (read as absent) | — |

Names are I6: they glue at 10⁵ and overflow at 10⁶. Semi-major axes and
energies carry two significant digits; the time has a resolution of 0.01 NB.
"""
Base.@kwdef struct ExchangeEvent
    line::Int
    interval::Int
    time_nb::Float64
    kind::Symbol = :exchange
    name1_initial::Int = 0
    name2_initial::Int = 0
    name1_final::Int = 0
    name2_final::Int = 0
    e0::Float64 = NaN
    e::Float64 = NaN
    a0_nb::Float64 = NaN
    a_nb::Float64 = NaN
    eb0_nb::Float64 = NaN
    eb_nb::Float64 = NaN
    overflow::Bool = false
end

"""
    TerminationEvent

A termination diagnostic of the engine log.

| kind | writer | fields |
|---|---|---|
| `:end_run_deferred` | `adjust.F` FORMAT 64: TCRIT or NCRIT reached while subsystems are active, repeated at every ADJUST until they end | `time_nb` = TTOT [NB], `n_subsystems` = NSUB |
| `:chain_tiny_steps` | `chain.f` FORMAT 403: chain terminated after 10⁶ steps with abs(STEP) < 10⁻¹⁴ | `n_steps` = NSTEP1, `n_members` = NCH, `step_nb` = STEP [NB], `separations_nb` = 1/RINV(1…NCH−1) [NB] |

`:chain_tiny_steps` has no time field (`time_nb = NaN`). The line holds STEP and
at most eight separations (1P 9E9.1, two digits); a ten-member chain sends its
ninth separation to a continuation record, which is not read.
"""
Base.@kwdef struct TerminationEvent
    line::Int
    interval::Int
    time_nb::Float64
    kind::Symbol
    n_subsystems::Int = 0
    n_steps::Int = 0
    n_members::Int = 0
    step_nb::Float64 = NaN
    separations_nb::NTuple{8,Float64} = _DYN_NAN8
    overflow::Bool = false
end

# ---------------------------------------------------------------------------
# FORMATs, transcribed descriptor by descriptor (`1P`/`0P` and the leading `/`
# do not change widths and are omitted)
# ---------------------------------------------------------------------------

# ksinit.F FORMAT 60 / ksterm.F FORMAT 15 / kspinit.f FORMAT 60 / ksinit.F
# FORMAT 65, everything after the leading label up to RAD1,2,S[*] 3E11.3.
const _KS_CORE_FORMAT = vcat(
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
)

const _KS_ORIENTATION_FORMAT = vcat(
    fmt_lit(" RCM,VCM[NB]="),
    fmt_real(2, 13),
    fmt_lit(" I1 RXY[NB],i,alp[rad] "),
    fmt_real(3, 13),
    fmt_lit(" I1 VXY[NB],i,alp[rad] "),
    fmt_real(3, 13),
    fmt_lit(" I2 RXY[NB],i,alp[rad] "),
    fmt_real(3, 13),
    fmt_lit(" I2 VXY[NB],i,alp[rad] "),
    fmt_real(3, 13),
)

const _KS_NEW_FORMAT =
    vcat(fmt_lit(" NEW KSREG   TIME[NB]"), _KS_CORE_FORMAT, _KS_ORIENTATION_FORMAT)
const _KS_END_FORMAT =
    vcat(fmt_lit(" END KSREG   TIME[NB]"), _KS_CORE_FORMAT, _KS_ORIENTATION_FORMAT)
const _KS_PRIMORDIAL_FORMAT =
    vcat(fmt_lit(" NEW KSREG   TIME[NB]"), _KS_CORE_FORMAT, fmt_lit(" RI,VI[NB]="), fmt_real(2, 11))
# ksinit.F:275-281 (transcribed from the engine source; the dossier gives the
# layout in prose only).
const _KS_COMPACT_FORMAT = vcat(
    fmt_lit(" NS/BH BINARY T[NB]"),
    _KS_CORE_FORMAT,
    fmt_lit(" RI,VI[NB],DEB="),
    fmt_real(3, 11),
)

# impact.f FORMAT 20 (NEW TRIPLE/QUAD/CHAIN/MERGER/MERGE2); A8 is skipped. In
# the quadruple-merger branch EB may refer to the binary before an IPAIR swap
# (impact.f:764-779) [UNSURE in the dossier]; it is stored as printed.
const _IMPACT_FORMAT = vcat(
    fmt_lit(" NEW"),
    fmt_skip(8),
    fmt_lit(" TIME"),
    fmt_real(1, 12),
    fmt_lit(" INCM"),
    fmt_int(1, 9),
    fmt_lit(" I3"),
    fmt_int(1, 9),
    fmt_lit(" INPAIR"),
    fmt_int(1, 9),
    fmt_lit(" N1,2,3,INCM"),
    fmt_int(4, 10),
    fmt_lit(" KW1,2,3,IN"),
    fmt_int(4, 4),
    fmt_lit(" M1,2[NB]"),
    fmt_real(2, 11),
    fmt_lit(" M3,TOT"),
    fmt_real(2, 11),
    fmt_lit("  R12"),
    fmt_real(1, 11),
    fmt_lit(" H"),
    fmt_real(1, 11),
    fmt_lit(" IN e,a,eb[NB]"),
    fmt_real(2, 12),
    fmt_real(1, 11),
    fmt_lit(" P[d]"),
    fmt_real(1, 11),
    fmt_lit(" OUT e,a,eb[NB]"),
    fmt_real(2, 12),
    fmt_real(1, 11),
    fmt_lit(" P[d]"),
    fmt_real(1, 11),
    fmt_lit("  G4"),
    fmt_real(1, 10),
    fmt_lit("  RIN3"),
    fmt_real(1, 8),
    fmt_lit("  PERIM"),
    fmt_real(1, 8),
    fmt_lit("  EB1/EB0 "),
    fmt_real(1, 8),
    fmt_lit("  NP"),
    fmt_int(1, 4),
    fmt_lit(" M1,2,3,TOT[*]"),
    fmt_real(4, 11),
    fmt_lit(" RAD1,2,3[*]"),
    fmt_real(3, 11),
    fmt_lit(" IN,OUT Sep[*]"),
    fmt_real(2, 11),
    fmt_lit("  RI,VI[NB]"),
    fmt_real(2, 11),
)

# The time of the impact line, for records that take their stamp from it.
const _CHAIN_NEW_TIME_FORMAT = vcat(fmt_lit(" NEW CHAIN   TIME"), fmt_real(1, 12))
const _MERGE2_NEW_TIME_FORMAT = vcat(fmt_lit(" NEW MERGE2  TIME"), fmt_real(1, 12))

# chain.f FORMAT 15. NAMEC is dimensioned 6, so the 6I10 block holds all names;
# a reverted second record for NCH > 6 [UNSURE in the dossier] is not read.
const _CHAIN_INIT_FORMAT = vcat(
    fmt_lit(" NEW CHAIN   NCH NP ENERGY RSUM RGRAV TCR RMAXS NAMEC(1-NM) [ALL in NB UNIT]"),
    fmt_int(2, 4),
    fmt_real(1, 10),
    fmt_real(4, 9),
    fmt_int(6, 10),
)

# chain.f FORMAT 55, NSTEP1 and T0S + TIMEC only (time stamp of `:chain_init`).
const _CHAIN_STEP_FORMAT =
    vcat(fmt_lit(" CHAIN: NSTEP T DTR R NAMEC "), fmt_int(1, 5), fmt_real(1, 10))

# recoil.f FORMAT 30. EC = ENERGY/CM(8) with CM(8) = E(3) at CHINIT; the
# physical meaning of E(3) there is [UNSURE in the dossier], EC is stored as
# printed.
const _CHAIN_END_FORMAT = vcat(
    fmt_lit(" END CHAIN NCH I1-4 "),
    fmt_int(5, 3),
    fmt_lit(" N1-4 ="),
    fmt_int(4, 10),
    fmt_lit(" M1-4 ="),
    fmt_real(4, 13),
    fmt_lit(" RB ="),
    fmt_real(1, 13),
    fmt_lit(" R13 ="),
    fmt_real(1, 13),
    fmt_lit(" R24 ="),
    fmt_real(1, 13),
    fmt_lit(" DE ="),
    fmt_real(1, 13),
    fmt_lit(" TC ="),
    fmt_real(1, 13),
    fmt_lit(" #"),
    fmt_int(1, 6),
    fmt_int(1, 5),
    fmt_int(1, 4),
    fmt_lit(" GPERT ="),
    fmt_real(1, 13),
    fmt_lit(" DB ="),
    fmt_real(1, 13),
    fmt_lit(" EC ="),
    fmt_real(1, 13),
)

# chain.f FORMAT 80 and the time of the following FORMAT 303 (CHEND) line;
# transcribed from the engine source (chain.f:847-857).
const _CHAIN_COLLISION_FORMAT =
    vcat(fmt_lit(" END CHAIN:   RIJ RCOLL QPERI GPERT "), fmt_real(6, 15))
const _CHEND_TIME_FORMAT = vcat(fmt_lit(" CHEND: TNB,TC[NB] K,M,R[*],R,V[NB-CH]="), fmt_real(1, 17))

# chterm.f FORMAT 4
const _CHTERM_NAMES_FORMAT = vcat(
    fmt_lit(" CHTERM: NAMEC "),
    fmt_int(4, 10),
    fmt_lit(" M1-4 "),
    fmt_real(4, 15),
    fmt_lit("  RIJ[NB] (1-2 1-3 2-3 2-4 3-4) "),
    fmt_real(5, 15),
    fmt_lit(" M[*] "),
    fmt_real(4, 15),
)

# chterm.f FORMAT 75
const _CHTERM_FORMAT = vcat(
    fmt_lit(" CHTERM: TIME[NB] I3 I4 ICM1 NI3,4,CM  STEP(I3,I4,ICM1)[NB] "),
    fmt_real(1, 10),
    fmt_int(6, 10),
    fmt_real(3, 10),
)

# chterm.f FORMAT 90 (transcribed from the engine source, chterm.f:519-520)
const _CHTERM_SECOND_FORMAT = vcat(
    fmt_lit(" CHTERM  SECOND BINARY(ICM2)[All in NB]: I3 I4 NI3,4 STEP STEPR R12 H GAMMA"),
    fmt_int(4, 10),
    fmt_real(5, 10),
)

# impact.f FORMAT 22
const _CHAIN_BB_FORMAT = vcat(
    fmt_lit(" CHAIN B-B    NAME(I1-I4)"),
    fmt_int(4, 10),
    fmt_lit("  K*(INCM)"),
    fmt_int(1, 3),
    fmt_lit("  K*(OCM)"),
    fmt_int(1, 3),
    fmt_lit("  ECC0"),
    fmt_real(1, 6),
    fmt_lit("  ECC1"),
    fmt_real(1, 6),
    fmt_lit("  SEMI0"),
    fmt_real(1, 11),
    fmt_lit("  SEMI2"),
    fmt_real(1, 11),
    fmt_lit("  RIN3"),
    fmt_real(1, 11),
    fmt_lit("  SEMI1"),
    fmt_real(1, 11),
    fmt_lit("  PERIM"),
    fmt_real(1, 11),
)

# chmod_chain.f FORMAT 28 (binary escaper; the second record is the reverted
# format holding RB and VINF) and FORMAT 35 (single escaper). At the first
# FORMAT 28 call site (chmod_chain.f:419-425) VINF may be stale [UNSURE in the
# dossier]; it is stored as printed.
const _ESCAPE_BINARY_FORMAT = vcat(
    fmt_lit(" CHAIN ESCAPE: TT T TC ="),
    fmt_real(3, 15),
    fmt_lit(" IESC JESC NM1/2 M1/2 RI VI ="),
    fmt_int(2, 4),
    fmt_int(2, 8),
    fmt_real(2, 15),
    fmt_lit(" GPERT 2*M/R RB VF ="),
    fmt_real(4, 15),
)
const _ESCAPE_CONTINUATION_FORMAT = vcat(fmt_lit(" CHAIN ESCAPE: TT T TC ="), fmt_real(2, 15))
const _ESCAPE_SINGLE_FORMAT = vcat(
    fmt_lit(" CHAIN ESCAPE: TT T TC ="),
    fmt_real(3, 15),
    fmt_lit(" IESC NM RI VI ="),
    fmt_int(1, 4),
    fmt_int(1, 8),
    fmt_real(2, 15),
    fmt_lit(" GPERT 2*M/R VF "),
    fmt_real(3, 15),
)

# recoil.f FORMAT 20
const _CHAIN_BINARY_FORMAT = vcat(
    fmt_lit(" CHAIN BINARY  NAM ="),
    fmt_int(2, 10),
    fmt_lit("  A ="),
    fmt_real(1, 12),
    fmt_lit("  E ="),
    fmt_real(1, 12),
    fmt_lit("  EB ="),
    fmt_real(1, 12),
    fmt_lit("  GB ="),
    fmt_real(1, 12),
    fmt_lit("  G4 ="),
    fmt_real(1, 12),
    fmt_lit("  EB1 ="),
    fmt_real(1, 12),
    fmt_lit("  E1 ="),
    fmt_real(1, 12),
    fmt_lit("  ET ="),
    fmt_real(1, 12),
    fmt_lit("  DB ="),
    fmt_real(1, 12),
)

# reset.f FORMAT 21
const _END_MERGER_FORMAT = vcat(
    fmt_lit(" END MERGER:  TIME"),
    fmt_real(1, 12),
    fmt_lit(" IC,JC,IP="),
    fmt_int(3, 10),
    fmt_lit(" N IC,JC,ICM="),
    fmt_int(3, 10),
    fmt_lit(" KW IC,JC,ICM="),
    fmt_int(3, 4),
    fmt_lit(" M1,2[NB]"),
    fmt_real(2, 11),
    fmt_lit(" JX1,JX2,JP="),
    fmt_int(3, 10),
    fmt_lit(" N JX1,JX2,JP="),
    fmt_int(3, 10),
    fmt_lit(" KW JX1,JX2,JCM="),
    fmt_int(3, 4),
    fmt_lit(" M3,4[NB]"),
    fmt_real(2, 11),
    fmt_lit("  R12"),
    fmt_real(1, 11),
    fmt_lit(" H"),
    fmt_real(1, 11),
    fmt_lit(" R34"),
    fmt_real(1, 11),
    fmt_lit(" H2"),
    fmt_real(1, 11),
    fmt_lit(" M1,2[*]"),
    fmt_real(2, 11),
    fmt_lit(" RAD1,2[*]"),
    fmt_real(2, 11),
    fmt_lit("  M3,4[*]"),
    fmt_real(2, 11),
    fmt_lit(" RAD3,4[*]"),
    fmt_real(2, 11),
    fmt_lit(" SEPI/J"),
    fmt_real(2, 11),
)

# impact.f FORMAT 60 (transcribed from the engine source, impact.f:1038-1039)
const _HI_MERGE_FORMAT = vcat(
    fmt_lit(" HI MERGE    NAM E E1 A A1 PM PC "),
    fmt_int(6, 7),
    fmt_real(2, 7),
    fmt_real(4, 11),
)

# merge2.f FORMAT 5 (A11 skipped), 7, 8, 10
const _MERGE2_MULTIPLE_FORMAT = vcat(
    fmt_lit(" NEW"),
    fmt_skip(11),
    fmt_lit("   T MT NM1 NM2 NM3 RI E0 EB0 A0 PC G0"),
    fmt_real(1, 10),
    fmt_real(1, 6),
    fmt_int(3, 6),
    fmt_real(2, 6),
    fmt_real(4, 10),
)
const _MERGE2_QUARTET_FORMAT = vcat(
    fmt_lit(" NEW QUARTET    T MT NM1 NMG NM3 E0 A0 PC G0 "),
    fmt_real(1, 9),
    fmt_real(1, 6),
    fmt_int(3, 6),
    fmt_real(1, 6),
    fmt_real(3, 10),
)
const _MERGE2_QUINTUP2_FORMAT = vcat(
    fmt_lit(" NEW QUINTUP2    T MT NM1 NMG NM3 E0 A0 PC G0"),
    fmt_real(1, 10),
    fmt_real(1, 6),
    fmt_int(3, 6),
    fmt_real(1, 6),
    fmt_real(3, 10),
)
const _MERGE2_HITRIP_FORMAT = vcat(
    fmt_lit(" NEW HITRIP    T MT NM E0 AI AJ RI RJ PC GX "),
    fmt_real(1, 9),
    fmt_real(1, 6),
    fmt_int(4, 6),
    fmt_real(1, 6),
    fmt_real(6, 10),
)

# reset2.f FORMAT 5 (A11 skipped), 9, 11, 18, 48, 65
const _RESET2_MULTIPLE_FORMAT = vcat(
    fmt_lit(" END"),
    fmt_skip(11),
    fmt_lit("   T MT NM1 NM2 NM3 RI E1 EB1 A1 PM G1 "),
    fmt_real(1, 9),
    fmt_real(1, 6),
    fmt_int(3, 6),
    fmt_real(2, 6),
    fmt_real(1, 6),
    fmt_real(3, 10),
)
const _RESET2_QUARTET_FORMAT = vcat(
    fmt_lit(" END QUARTET    T MT NM1 NMG NM3 E1 A1 R1 G1 "),
    fmt_real(1, 9),
    fmt_real(1, 6),
    fmt_int(3, 6),
    fmt_real(1, 6),
    fmt_real(3, 10),
)
const _RESET2_QUINTUP2_FORMAT = vcat(
    fmt_lit(" END QUINTUP2    T MT NM1 NMG NM3 E1 A1 R1 G1"),
    fmt_real(1, 10),
    fmt_real(1, 6),
    fmt_int(3, 6),
    fmt_real(1, 6),
    fmt_real(3, 10),
)
const _RESET2_HITRIP_FORMAT = vcat(
    fmt_lit(" END HITRIP    T MT NM E1 AI AJ R1 A1 PM G1 "),
    fmt_real(1, 9),
    fmt_real(1, 6),
    fmt_int(4, 6),
    fmt_real(1, 6),
    fmt_real(6, 10),
)
const _RESET2_OUTER_FORMAT = vcat(
    fmt_lit(" END OUTER MERGER"),
    fmt_int(1, 5),
    fmt_lit("  H ="),
    fmt_real(1, 7),
    fmt_lit("  M ="),
    fmt_real(2, 7),
    fmt_lit("  E1 ="),
    fmt_real(1, 6),
    fmt_lit("  EB2 ="),
    fmt_real(1, 6),
    fmt_lit("  RB2 ="),
    fmt_real(1, 8),
    fmt_lit("  G2 ="),
    fmt_real(1, 8),
    fmt_lit("  DP ="),
    fmt_real(1, 8),
)
const _RESET2_MERGE2_FORMAT = vcat(
    fmt_lit(" END MERGE2"),
    fmt_int(1, 3),
    fmt_lit("  T ="),
    fmt_real(1, 8),
    fmt_lit("  M ="),
    fmt_real(2, 7),
    fmt_lit("  R1 ="),
    fmt_real(1, 8),
    fmt_lit("  A1 ="),
    fmt_real(1, 8),
    fmt_lit("  EB ="),
    fmt_real(1, 6),
    fmt_lit("  E1 ="),
    fmt_real(1, 6),
    fmt_lit("  GB ="),
    fmt_real(1, 8),
    fmt_lit("  G ="),
    fmt_real(1, 6),
    fmt_lit("  NB ="),
    fmt_int(1, 3),
)

# recoil.f FORMAT 15
const _EXCHANGE_FORMAT = vcat(
    fmt_lit(" EXCHANGE    T NAM E0 E A0 A EB0 EB "),
    fmt_real(1, 9),
    fmt_int(4, 6),
    fmt_real(2, 7),
    fmt_real(2, 9),
    fmt_real(2, 10),
)

# adjust.F FORMAT 64 (9X written as nine blanks) and chain.f FORMAT 403
const _END_RUN_DEFERRED_FORMAT = vcat(
    fmt_lit(repeat(" ", 9) * "END RUN DEFERRED   TTOT ="),
    fmt_real(1, 16),
    fmt_lit("  NSUB ="),
    fmt_int(1, 3),
)
const _TINY_STEPS_FORMAT = vcat(
    fmt_lit(" ENFORCED CHAIN TINY STEPS:  NSTEP NCH STEP R"),
    fmt_int(1, 10),
    fmt_int(1, 4),
    fmt_real(9, 9),
)

# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

# Full labels (leading blanks removed) and their kinds; NEW KSREG and CHAIN
# ESCAPE are told apart by later labels in `_dynamics_line_kind`. [UNSURE in
# the format dossier] `NEW QUAD    TIME` is written by both the regularisation
# and the merger branch of impact.f with one layout; it is classified as
# `:quad_new`, and only the following END QUAD / END MERGER tells them apart.
const _DYNAMICS_LABELS = (
    "END KSREG   TIME[NB]" => :ks_end,
    "NS/BH BINARY T[NB]" => :ks_new_compact,
    "NEW TRIPLE  TIME" => :triple_new,
    "NEW QUAD    TIME" => :quad_new,
    "NEW CHAIN   TIME" => :chain_new,
    "NEW CHAIN   NCH NP ENERGY RSUM RGRAV TCR RMAXS" => :chain_init,
    "END CHAIN NCH I1-4" => :chain_end,
    "END CHAIN:   RIJ RCOLL QPERI GPERT" => :chain_end_collision,
    "CHTERM: NAMEC" => :chterm_names,
    "CHTERM: TIME[NB]" => :chterm,
    "CHTERM  SECOND BINARY(ICM2)" => :chterm_second_binary,
    "CHAIN B-B    NAME(I1-I4)" => :chain_bb,
    "CHAIN BINARY  NAM =" => :chain_binary,
    "NEW MERGER  TIME" => :hierarchy_new,
    "NEW MERGE2  TIME" => :hierarchy_new_merge2,
    "HI MERGE    NAM E E1 A A1 PM PC" => :hierarchy_hi_merge,
    "END MERGER:  TIME" => :hierarchy_end,
    "NEW QUINTUPLET   T MT NM1 NM2 NM3" => :hierarchy_new_quintuplet,
    "NEW SEXTUPLET    T MT NM1 NM2 NM3" => :hierarchy_new_sextuplet,
    "NEW QUARTET    T MT NM1 NMG NM3" => :hierarchy_new_quartet,
    "NEW QUINTUP2    T MT NM1 NMG NM3" => :hierarchy_new_quintup2,
    "NEW HITRIP    T MT NM" => :hierarchy_new_hitrip,
    "END QUINTUPLET   T MT NM1 NM2 NM3" => :hierarchy_end_quintuplet,
    "END SEXTUPLET    T MT NM1 NM2 NM3" => :hierarchy_end_sextuplet,
    "END QUARTET    T MT NM1 NMG NM3" => :hierarchy_end_quartet,
    "END QUINTUP2    T MT NM1 NMG NM3" => :hierarchy_end_quintup2,
    "END HITRIP    T MT NM" => :hierarchy_end_hitrip,
    "END OUTER MERGER" => :hierarchy_end_outer,
    "END MERGE2" => :hierarchy_end_merge2,
    "EXCHANGE    T NAM E0 E A0 A EB0 EB" => :exchange,
    "END RUN DEFERRED   TTOT =" => :end_run_deferred,
    "ENFORCED CHAIN TINY STEPS:" => :chain_tiny_steps,
)

"""
    _dynamics_line_kind(s::AbstractString) -> Symbol

Kind of a log line `s` (leading blanks removed) of the KS, few-body, hierarchy,
exchange and termination family, or `:none`. Lines are identified by their full
label: `NEW KSREG` by `RCM,VCM[NB]=` (`:ks_new`, `ksinit.F`) or `RI,VI[NB]=`
(`:ks_new_primordial`, `kspinit.f`), `CHAIN ESCAPE: TT T TC =` by `IESC JESC`
(`:chain_escape_binary`) or `IESC NM` (`:chain_escape_single`). The continuation
record of a binary CHAIN ESCAPE, the per-step chain lines (`CHAIN: NSTEP`,
`CHAIN T,TC[NB]`), `BEG KSRECT`, `MDNEW`/`MDEND KSREG` and the GR, escape, kick,
circularisation and Roche lines return `:none`. All kinds are listed in
[`DYNAMICS_EVENT_KINDS`](@ref).
"""
function _dynamics_line_kind(s::AbstractString)::Symbol
    if startswith(s, "NEW KSREG   TIME[NB]")
        return occursin(" RI,VI[NB]=", s) ? :ks_new_primordial : :ks_new
    end
    if startswith(s, "CHAIN ESCAPE: TT T TC =")
        occursin(" IESC JESC NM1/2", s) && return :chain_escape_binary
        occursin(" IESC NM RI VI =", s) && return :chain_escape_single
        return :none
    end
    for (label, kind) in _DYNAMICS_LABELS
        startswith(s, label) && return kind
    end
    return :none
end

# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

_dyn_name(v::Float64)::Int = _fixed_int(v, 0)
_dyn_kstar(v::Float64)::Int = _fixed_int(v, KSTAR_ABSENT)

"""Names `v[r]`, padded with 0 to six entries."""
_dyn_names6(v::Vector{Float64}, r::AbstractUnitRange{Int})::NTuple{6,Int} =
    ntuple(k -> k <= length(r) ? _dyn_name(v[r[k]]) : 0, Val(6))

"""Stellar types `v[r]`, padded with `KSTAR_ABSENT` to four entries."""
_dyn_kstars4(v::Vector{Float64}, r::AbstractUnitRange{Int})::NTuple{4,Int} =
    ntuple(k -> k <= length(r) ? _dyn_kstar(v[r[k]]) : KSTAR_ABSENT, Val(4))

"""Indices `v[r]`, padded with 0 to four entries."""
_dyn_indices4(v::Vector{Float64}, r::AbstractUnitRange{Int})::NTuple{4,Int} =
    ntuple(k -> k <= length(r) ? _dyn_name(v[r[k]]) : 0, Val(4))

"""Reals `v[r]`, padded with NaN to four entries."""
_dyn_reals4(v::Vector{Float64}, r::AbstractUnitRange{Int})::NTuple{4,Float64} =
    ntuple(k -> k <= length(r) ? v[r[k]] : NaN, Val(4))

"""
Value `k` of the first line within `window` lines of `i` (forward for `step =
1`, backward for `step = -1`) whose label is `prefix`, read by `items`; NaN when
there is none or it does not parse.
"""
function _dyn_nearby_value(
    lines::AbstractVector{<:AbstractString},
    i::Int,
    step::Int,
    window::Int,
    prefix::AbstractString,
    items::AbstractVector{FortranItem},
    k::Int,
)::Float64
    j = i + step
    for _ in 1:window
        checkbounds(Bool, lines, j) || break
        if startswith(lstrip(lines[j]), prefix)
            rec = try
                read_fixed(lines[j], items)
            catch e
                e isa ArgumentError || rethrow()
                return NaN
            end
            return length(rec.values) >= k ? rec.values[k] : NaN
        end
        j += step
    end
    return NaN
end

"""Stamp of the `CHTERM: TIME[NB]` line following (`step = 1`) or preceding (`-1`) line `i`."""
_dyn_chterm_time(lines, i::Int, step::Int)::Float64 =
    _dyn_nearby_value(lines, i, step, _DYN_TIME_WINDOW, "CHTERM: TIME[NB]", _CHTERM_FORMAT, 1)

"""Read a KS record (`ksinit.F`, `kspinit.f`, `ksterm.F`)."""
function _dyn_parse_ks(kind::Symbol, line::AbstractString, i::Int, interval::Int)::KSEvent
    items = if kind === :ks_new
        _KS_NEW_FORMAT
    elseif kind === :ks_end
        _KS_END_FORMAT
    elseif kind === :ks_new_primordial
        _KS_PRIMORDIAL_FORMAT
    else
        _KS_COMPACT_FORMAT
    end
    rec = read_fixed(line, items)
    v = rec.values
    long = kind === :ks_new || kind === :ks_end
    o(k) = long ? v[28 + k] : NaN
    return KSEvent(;
        line = i,
        interval = interval,
        time_nb = v[1],
        kind = kind,
        name1 = _dyn_name(v[2]),
        name2 = _dyn_name(v[3]),
        name_cm = _dyn_name(v[4]),
        kstar1 = _dyn_kstar(v[5]),
        kstar2 = _dyn_kstar(v[6]),
        kstar_cm = _dyn_kstar(v[7]),
        ipair = _dyn_name(v[8]),
        dtau = v[9],
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
        n_perturbers = _dyn_name(v[20]),
        n_neighbours_cm = _dyn_name(v[21]),
        m1_msun = v[22],
        m2_msun = v[23],
        radius1_rsun = v[24],
        radius2_rsun = v[25],
        # FORMAT field 26 (label RAD1,2,S[*]) is R(IPAIR)*SU, the separation.
        r12_rsun = v[26],
        r_cm_nb = v[27],
        v_cm_nb = v[28],
        deb_nb = kind === :ks_new_compact ? v[29] : NaN,
        r1_xy_nb = o(1),
        r1_lat = o(2),
        r1_azi = o(3),
        v1_xy_nb = o(4),
        v1_lat = o(5),
        v1_azi = o(6),
        r2_xy_nb = o(7),
        r2_lat = o(8),
        r2_azi = o(9),
        v2_xy_nb = o(10),
        v2_lat = o(11),
        v2_azi = o(12),
        overflow = rec.overflow,
    )
end

"""
Inner and outer impact periods rescaled to the masses of their own orbits: the
engine uses M_tot for the inner and M_in for the outer orbit.
"""
function _dyn_impact_periods(pd::Float64, pd1::Float64, m_inner::Float64, m_total::Float64)
    ok = isfinite(m_inner) && isfinite(m_total) && m_inner > 0 && m_total > 0
    ok || return (NaN, NaN)
    s = sqrt(m_total / m_inner)
    return (pd * s, pd1 / s)
end

"""Fields of an `impact.f` FORMAT 20 record, shared by `FewBodyEvent` and `HierarchyEvent`."""
function _dyn_impact_fields(line::AbstractString)
    rec = read_fixed(line, _IMPACT_FORMAT)
    v = rec.values
    period_in, period_out = _dyn_impact_periods(v[22], v[26], v[13] + v[14], v[16])
    return (;
        time_nb = v[1],
        indices = _dyn_indices4(v, 2:4),
        names = _dyn_names6(v, 5:8),
        kstar_values = _dyn_kstars4(v, 9:12),
        masses_nb = _dyn_reals4(v, 13:16),
        r12_nb = v[17],
        h_nb = v[18],
        ecc_in = v[19],
        semi_in_nb = v[20],
        eb_in_nb = v[21],
        period_in_d = period_in,
        ecc_out = v[23],
        semi_out_nb = v[24],
        eb_out_nb = v[25],
        period_out_d = period_out,
        g4 = v[27],
        rij_nb = v[28],
        pmin_nb = v[29],
        eb_ratio = v[30],
        n_perturbers = _dyn_name(v[31]),
        masses_msun = _dyn_reals4(v, 32:35),
        radii_rsun = (v[36], v[37], v[38]),
        r12_rsun = v[39],
        rij_rsun = v[40],
        r_cm_nb = v[41],
        v_cm_nb = v[42],
        overflow = rec.overflow,
    )
end

"""Read a chain or few-body record other than the impact lines and the escapes."""
function _dyn_parse_chain(
    kind::Symbol,
    lines::AbstractVector{<:AbstractString},
    i::Int,
    interval::Int,
)::FewBodyEvent
    line = lines[i]
    if kind === :chain_init
        rec = read_fixed(line, _CHAIN_INIT_FORMAT)
        v = rec.values
        t = _dyn_nearby_value(lines, i, 1, 3, "CHAIN: NSTEP", _CHAIN_STEP_FORMAT, 2)
        return FewBodyEvent(;
            line = i,
            interval = interval,
            time_nb = t,
            kind = kind,
            n_members = _dyn_name(v[1]),
            n_perturbers = _dyn_name(v[2]),
            energy_nb = v[3],
            rsum_nb = v[4],
            rgrav_nb = v[5],
            tcr_nb = v[6],
            rmaxs_nb = v[7],
            names = _dyn_names6(v, 8:13),
            overflow = rec.overflow,
        )
    elseif kind === :chain_end
        rec = read_fixed(line, _CHAIN_END_FORMAT)
        v = rec.values
        n = _dyn_name(v[1])
        # NCH ≤ 3: I4 = I1 repeats the closest member, R24 = 1e10 (recoil.f).
        # The CHTERM: TIME[NB] line of the same termination follows END CHAIN
        # (order observed in the logs: END CHAIN, CHTERM: NAMEC, NEW KSREG,
        # CHTERM: TIME[NB]).
        short = 0 < n <= 3
        return FewBodyEvent(;
            line = i,
            interval = interval,
            time_nb = _dyn_chterm_time(lines, i, 1),
            kind = kind,
            n_members = n,
            indices = (
                _dyn_name(v[2]),
                _dyn_name(v[3]),
                _dyn_name(v[4]),
                short ? 0 : _dyn_name(v[5]),
            ),
            names = (
                _dyn_name(v[6]),
                _dyn_name(v[7]),
                _dyn_name(v[8]),
                short ? 0 : _dyn_name(v[9]),
                0,
                0,
            ),
            masses_nb = (v[10], v[11], v[12], short ? NaN : v[13]),
            separations_nb = (v[14], v[15], short ? NaN : v[16], NaN, NaN),
            de_nb = v[17],
            duration_crossings = v[18],
            n_steps = _dyn_name(v[19]),
            n_reregularisations = _dyn_name(v[20]),
            n_perturbers = _dyn_name(v[21]),
            gpert = v[22],
            db = v[23],
            ec = v[24],
            overflow = rec.overflow,
        )
    elseif kind === :chain_end_collision
        rec = read_fixed(line, _CHAIN_COLLISION_FORMAT)
        v = rec.values
        t = _dyn_nearby_value(lines, i, 1, 3, "CHEND:", _CHEND_TIME_FORMAT, 1)
        return FewBodyEvent(;
            line = i,
            interval = interval,
            time_nb = t,
            kind = kind,
            separations_nb = (v[1], v[2], v[3], NaN, NaN),
            rcoll_nb = v[4],
            qperi_nb = v[5],
            gpert = v[6],
            overflow = rec.overflow,
        )
    elseif kind === :chterm_names
        rec = read_fixed(line, _CHTERM_NAMES_FORMAT)
        v = rec.values
        # RIJ(3-4) = 0 means I4 = I3: three-body chain, NAMEC(4) and M(4) stale.
        three = v[13] == 0.0
        return FewBodyEvent(;
            line = i,
            interval = interval,
            time_nb = _dyn_chterm_time(lines, i, 1),
            kind = kind,
            n_members = three ? 3 : 0,
            names = (
                _dyn_name(v[1]),
                _dyn_name(v[2]),
                _dyn_name(v[3]),
                three ? 0 : _dyn_name(v[4]),
                0,
                0,
            ),
            masses_nb = (v[5], v[6], v[7], three ? NaN : v[8]),
            separations_nb = (v[9], v[10], v[11], three ? NaN : v[12], three ? NaN : v[13]),
            masses_msun = (v[14], v[15], v[16], three ? NaN : v[17]),
            overflow = rec.overflow,
        )
    elseif kind === :chterm
        rec = read_fixed(line, _CHTERM_FORMAT)
        v = rec.values
        return FewBodyEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            indices = _dyn_indices4(v, 2:4),
            names = _dyn_names6(v, 5:7),
            steps_nb = (v[8], v[9], v[10]),
            overflow = rec.overflow,
        )
    elseif kind === :chterm_second_binary
        rec = read_fixed(line, _CHTERM_SECOND_FORMAT)
        v = rec.values
        return FewBodyEvent(;
            line = i,
            interval = interval,
            time_nb = _dyn_chterm_time(lines, i, -1),
            kind = kind,
            indices = _dyn_indices4(v, 1:2),
            names = _dyn_names6(v, 3:4),
            steps_nb = (v[5], v[6], NaN),
            r12_nb = v[7],
            h_nb = v[8],
            gamma = v[9],
            overflow = rec.overflow,
        )
    elseif kind === :chain_bb
        rec = read_fixed(line, _CHAIN_BB_FORMAT)
        v = rec.values
        # IMPACT writes its NEW CHAIN line just before this one; a later one
        # belongs to the next chain and is the fallback only.
        t = _dyn_nearby_value(lines, i, -1, 3, "NEW CHAIN   TIME", _CHAIN_NEW_TIME_FORMAT, 1)
        isnan(t) &&
            (t = _dyn_nearby_value(lines, i, 1, 3, "NEW CHAIN   TIME", _CHAIN_NEW_TIME_FORMAT, 1))
        return FewBodyEvent(;
            line = i,
            interval = interval,
            time_nb = t,
            kind = kind,
            n_members = 4,
            names = _dyn_names6(v, 1:4),
            kstar_values = _dyn_kstars4(v, 5:6),
            ecc_in = v[7],
            ecc_out = v[8],
            semi_in_nb = v[9],
            semi_second_nb = v[10],
            rij_nb = v[11],
            semi_out_nb = v[12],
            pmin_nb = v[13],
            overflow = rec.overflow,
        )
    elseif kind === :chain_escape_single
        rec = read_fixed(line, _ESCAPE_SINGLE_FORMAT)
        v = rec.values
        return FewBodyEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            time_block_nb = v[2],
            time_chain_nb = v[3],
            indices = _dyn_indices4(v, 4:4),
            names = _dyn_names6(v, 5:5),
            r_escape_nb = v[6],
            rdot_escape_nb = v[7],
            gpert = v[8],
            v2_escape_nb = v[9],
            vinf_kms = v[10],
            overflow = rec.overflow,
        )
    elseif kind === :chain_binary
        rec = read_fixed(line, _CHAIN_BINARY_FORMAT)
        v = rec.values
        return FewBodyEvent(;
            line = i,
            interval = interval,
            time_nb = _dyn_chterm_time(lines, i, 1),
            kind = kind,
            names = _dyn_names6(v, 1:2),
            semi_in_nb = v[3],
            ecc_in = v[4],
            eb_over_energy = v[5],
            gpert = v[6],
            g4 = v[7],
            eb1_over_energy = v[8],
            e1_over_eb0 = v[9],
            et_over_cm8 = v[10],
            db = v[11],
            overflow = rec.overflow,
        )
    end
    throw(ArgumentError("Unknown chain event kind $kind"))
end

"""
Read the binary CHAIN ESCAPE record: line `i` and, when line `i + 1` is its
reverted continuation (label `CHAIN ESCAPE: TT T TC =` without `IESC`), the
latter; returns the event and the number of lines consumed.
"""
function _dyn_parse_escape_binary(
    lines::AbstractVector{<:AbstractString},
    i::Int,
    interval::Int,
)::Tuple{FewBodyEvent,Int}
    rec = read_fixed(lines[i], _ESCAPE_BINARY_FORMAT)
    v = rec.values
    rb, vinf, overflow, n = NaN, NaN, rec.overflow, 1
    if checkbounds(Bool, lines, i + 1)
        s = lstrip(lines[i + 1])
        if startswith(s, "CHAIN ESCAPE: TT T TC =") && !occursin("IESC", s)
            cont = read_fixed(lines[i + 1], _ESCAPE_CONTINUATION_FORMAT)
            rb, vinf = cont.values[1], cont.values[2]
            overflow |= cont.overflow
            n = 2
        end
    end
    # Line 1 labels after NM1/2 are shifted by two values: masses sit under
    # "RI VI", RI, RDOT, GPERT, 2M/R under "GPERT 2*M/R RB VF".
    event = FewBodyEvent(;
        line = i,
        interval = interval,
        time_nb = v[1],
        kind = :chain_escape_binary,
        time_block_nb = v[2],
        time_chain_nb = v[3],
        indices = _dyn_indices4(v, 4:5),
        names = _dyn_names6(v, 6:7),
        masses_nb = _dyn_reals4(v, 8:9),
        r_escape_nb = v[10],
        rdot_escape_nb = v[11],
        gpert = v[12],
        v2_escape_nb = v[13],
        rb_nb = rb,
        vinf_kms = vinf,
        overflow = overflow,
    )
    return (event, n)
end

"""Read a hierarchical-merger record."""
function _dyn_parse_hierarchy(
    kind::Symbol,
    lines::AbstractVector{<:AbstractString},
    i::Int,
    interval::Int,
)::HierarchyEvent
    line = lines[i]
    if kind === :hierarchy_new || kind === :hierarchy_new_merge2
        f = _dyn_impact_fields(line)
        return HierarchyEvent(; line = i, interval = interval, kind = kind, f...)
    elseif kind === :hierarchy_end
        rec = read_fixed(line, _END_MERGER_FORMAT)
        v = rec.values
        # Fields 4, 7, 10 and 13–27, 32–37 describe the pair now in the old
        # slot or repeat component 1 (reset.f:279-306); not stored. The
        # JCOMP > N branch of the JX block is held unreachable [UNSURE in the
        # dossier], which does not matter here since the block is dropped.
        return HierarchyEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            indices = _dyn_indices4(v, 2:3),
            names = _dyn_names6(v, 5:6),
            kstar_values = _dyn_kstars4(v, 8:9),
            masses_nb = _dyn_reals4(v, 11:12),
            masses_msun = _dyn_reals4(v, 28:29),
            radii_rsun = (v[30], v[31], NaN),
            overflow = rec.overflow,
        )
    elseif kind === :hierarchy_hi_merge
        rec = read_fixed(line, _HI_MERGE_FORMAT)
        v = rec.values
        t = _dyn_nearby_value(lines, i, -1, 3, "NEW MERGE2  TIME", _MERGE2_NEW_TIME_FORMAT, 1)
        return HierarchyEvent(;
            line = i,
            interval = interval,
            time_nb = t,
            kind = kind,
            names = _dyn_names6(v, 1:6),
            ecc_in = v[7],
            ecc_out = v[8],
            semi_in_nb = v[9],
            semi_out_nb = v[10],
            pmin_nb = v[11],
            pcrit_nb = v[12],
            overflow = rec.overflow,
        )
    elseif kind === :hierarchy_new_quintuplet || kind === :hierarchy_new_sextuplet
        rec = read_fixed(line, _MERGE2_MULTIPLE_FORMAT)
        v = rec.values
        return HierarchyEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            masses_msun = (NaN, NaN, NaN, v[2]),
            names = _dyn_names6(v, 3:5),
            r_cm_nb = v[6],
            ecc_in = v[7],
            eb_kinetic = v[8],
            semi_in_nb = v[9],
            pcrit_nb = v[10],
            gamma = v[11],
            overflow = rec.overflow,
        )
    elseif kind === :hierarchy_new_quartet || kind === :hierarchy_new_quintup2
        items = kind === :hierarchy_new_quartet ? _MERGE2_QUARTET_FORMAT : _MERGE2_QUINTUP2_FORMAT
        rec = read_fixed(line, items)
        v = rec.values
        return HierarchyEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            masses_msun = (NaN, NaN, NaN, v[2]),
            names = _dyn_names6(v, 3:5),
            ecc_in = v[6],
            semi_in_nb = v[7],
            pcrit_nb = v[8],
            gamma = v[9],
            overflow = rec.overflow,
        )
    elseif kind === :hierarchy_new_hitrip
        rec = read_fixed(line, _MERGE2_HITRIP_FORMAT)
        v = rec.values
        return HierarchyEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            masses_msun = (NaN, NaN, NaN, v[2]),
            names = _dyn_names6(v, 3:6),
            ecc_in = v[7],
            semi_inner1_nb = v[8],
            semi_inner2_nb = v[9],
            r12_nb = v[10],
            r34_nb = v[11],
            pcrit_nb = v[12],
            gamma = v[13],
            overflow = rec.overflow,
        )
    elseif kind === :hierarchy_end_quintuplet || kind === :hierarchy_end_sextuplet
        rec = read_fixed(line, _RESET2_MULTIPLE_FORMAT)
        v = rec.values
        return HierarchyEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            masses_msun = (NaN, NaN, NaN, v[2]),
            names = _dyn_names6(v, 3:5),
            r_cm_nb = v[6],
            ecc_out = v[7],
            eb_kinetic = v[8],
            semi_out_nb = v[9],
            pmin_nb = v[10],
            gamma_merged = v[11],
            overflow = rec.overflow,
        )
    elseif kind === :hierarchy_end_quartet || kind === :hierarchy_end_quintup2
        items = kind === :hierarchy_end_quartet ? _RESET2_QUARTET_FORMAT : _RESET2_QUINTUP2_FORMAT
        rec = read_fixed(line, items)
        v = rec.values
        return HierarchyEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            masses_msun = (NaN, NaN, NaN, v[2]),
            names = _dyn_names6(v, 3:5),
            ecc_out = v[6],
            semi_out_nb = v[7],
            r12_nb = v[8],
            gamma_merged = v[9],
            overflow = rec.overflow,
        )
    elseif kind === :hierarchy_end_hitrip
        rec = read_fixed(line, _RESET2_HITRIP_FORMAT)
        v = rec.values
        return HierarchyEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            masses_msun = (NaN, NaN, NaN, v[2]),
            names = _dyn_names6(v, 3:6),
            ecc_out = v[7],
            semi_inner1_nb = v[8],
            semi_inner2_nb = v[9],
            r12_nb = v[10],
            semi_out_nb = v[11],
            pmin_nb = v[12],
            gamma_merged = v[13],
            overflow = rec.overflow,
        )
    elseif kind === :hierarchy_end_outer
        rec = read_fixed(line, _RESET2_OUTER_FORMAT)
        v = rec.values
        t = _dyn_nearby_value(lines, i, 1, 3, "END MERGE2", _RESET2_MERGE2_FORMAT, 2)
        return HierarchyEvent(;
            line = i,
            interval = interval,
            time_nb = t,
            kind = kind,
            indices = (0, 0, _dyn_name(v[1]), 0),
            h_nb = v[2],
            masses_nb = _dyn_reals4(v, 3:4),
            energy_ratio = v[5],
            eb_fraction = v[6],
            r12_nb = v[7],
            gamma = v[8],
            dpot_nb = v[9],
            overflow = rec.overflow,
        )
    elseif kind === :hierarchy_end_merge2
        rec = read_fixed(line, _RESET2_MERGE2_FORMAT)
        v = rec.values
        return HierarchyEvent(;
            line = i,
            interval = interval,
            time_nb = v[2],
            kind = kind,
            merger_index = _dyn_name(v[1]),
            masses_nb = _dyn_reals4(v, 3:4),
            r12_nb = v[5],
            semi_out_nb = v[6],
            eb_fraction = v[7],
            energy_ratio = v[8],
            gamma = v[9],
            gamma_merged = v[10],
            n_neighbours = _dyn_name(v[11]),
            overflow = rec.overflow,
        )
    end
    throw(ArgumentError("Unknown hierarchy event kind $kind"))
end

"""Read an EXCHANGE record (`recoil.f` FORMAT 15)."""
function _dyn_parse_exchange(line::AbstractString, i::Int, interval::Int)::ExchangeEvent
    rec = read_fixed(line, _EXCHANGE_FORMAT)
    v = rec.values
    return ExchangeEvent(;
        line = i,
        interval = interval,
        time_nb = v[1],
        name1_initial = _dyn_name(v[2]),
        name2_initial = _dyn_name(v[3]),
        name1_final = _dyn_name(v[4]),
        name2_final = _dyn_name(v[5]),
        e0 = v[6],
        e = v[7],
        a0_nb = v[8],
        a_nb = v[9],
        eb0_nb = v[10],
        eb_nb = v[11],
        overflow = rec.overflow,
    )
end

"""Read an END RUN DEFERRED or ENFORCED CHAIN TINY STEPS record."""
function _dyn_parse_termination(
    kind::Symbol,
    line::AbstractString,
    i::Int,
    interval::Int,
)::TerminationEvent
    if kind === :end_run_deferred
        rec = read_fixed(line, _END_RUN_DEFERRED_FORMAT)
        v = rec.values
        return TerminationEvent(;
            line = i,
            interval = interval,
            time_nb = v[1],
            kind = kind,
            n_subsystems = _dyn_name(v[2]),
            overflow = rec.overflow,
        )
    elseif kind === :chain_tiny_steps
        rec = read_fixed(line, _TINY_STEPS_FORMAT)
        v = rec.values
        n = _dyn_name(v[2])
        # NCH − 1 separations follow STEP; the line holds at most eight.
        seps = ntuple(k -> k <= n - 1 ? v[3 + k] : NaN, Val(8))
        return TerminationEvent(;
            line = i,
            interval = interval,
            time_nb = NaN,
            kind = kind,
            n_steps = _dyn_name(v[1]),
            n_members = n,
            step_nb = v[3],
            separations_nb = seps,
            overflow = rec.overflow,
        )
    end
    throw(ArgumentError("Unknown termination event kind $kind"))
end

"""
    _parse_dynamics_line(kind, lines, i, interval) -> (record, n_consumed)

Read the record of `kind` (from [`_dynamics_line_kind`](@ref)) that starts at
`lines[i]`, a line of the engine log with its leading blanks intact, and return
it with the number of physical lines it occupies (2 for a binary CHAIN ESCAPE
followed by its continuation, 1 otherwise). `record` is a [`KSEvent`](@ref),
[`FewBodyEvent`](@ref), [`HierarchyEvent`](@ref), [`ExchangeEvent`](@ref) or
[`TerminationEvent`](@ref) with `line = i` and the given adjustment `interval`.
Records without a time field take it from a neighbouring line (see the type
docstrings), which may be read but is not consumed. A line that does not match
the FORMAT of its kind throws `ArgumentError`.
"""
function _parse_dynamics_line(
    kind::Symbol,
    lines::AbstractVector{<:AbstractString},
    i::Int,
    interval::Int,
)::Tuple{Union{KSEvent,FewBodyEvent,HierarchyEvent,ExchangeEvent,TerminationEvent},Int}
    checkbounds(Bool, lines, i) || throw(ArgumentError("Line index $i outside 1:$(length(lines))"))
    line = lines[i]
    if kind === :ks_new ||
       kind === :ks_new_primordial ||
       kind === :ks_end ||
       kind === :ks_new_compact
        return (_dyn_parse_ks(kind, line, i, interval), 1)
    elseif kind === :triple_new || kind === :quad_new || kind === :chain_new
        n = kind === :triple_new ? 3 : kind === :quad_new ? 4 : 0
        f = _dyn_impact_fields(line)
        event = FewBodyEvent(; line = i, interval = interval, kind = kind, n_members = n, f...)
        return (event, 1)
    elseif kind === :chain_escape_binary
        return _dyn_parse_escape_binary(lines, i, interval)
    elseif kind === :chain_init ||
           kind === :chain_end ||
           kind === :chain_end_collision ||
           kind === :chterm_names ||
           kind === :chterm ||
           kind === :chterm_second_binary ||
           kind === :chain_bb ||
           kind === :chain_escape_single ||
           kind === :chain_binary
        return (_dyn_parse_chain(kind, lines, i, interval), 1)
    elseif kind === :exchange
        return (_dyn_parse_exchange(line, i, interval), 1)
    elseif kind === :end_run_deferred || kind === :chain_tiny_steps
        return (_dyn_parse_termination(kind, line, i, interval), 1)
    elseif kind in DYNAMICS_EVENT_KINDS
        return (_dyn_parse_hierarchy(kind, lines, i, interval), 1)
    end
    throw(ArgumentError("Unknown dynamics event kind $kind"))
end
