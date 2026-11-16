# =============================================================================
# Physical constants of the engine (units.f)
# =============================================================================
#
# The engine converts its N-body units with these values, and everything the
# package writes for it or reads from it in physical units uses them too: a
# time given in Myr reaches the engine as the number of N-body time units the
# engine itself converts back to that Myr, and a length in solar radii is
# decoded with the solar radius it was encoded with.

"""Gravitational parameter G M☉ [cm³ s⁻²] of `units.f` (`GM = 6.6743D-08*1.9884D+33`)."""
const _ENGINE_GM_CGS = 6.6743e-8 * 1.9884e33

"""Parsec [cm] of `units.f` (`PC = 3.0856776D+18`)."""
const _ENGINE_PC_CM = 3.0856776e18

"""Solar radius [cm] of `units.f` (`SU = PC/(AU*6.955D+10)*RBAR*AU`)."""
const _ENGINE_RSUN_CM = 6.955e10

"""Year [s] of `units.f`, which converts the time unit to Myr with 3.15×10¹³ s per Myr."""
const _ENGINE_YR_S = 3.15e7

"""
Time unit [Myr] of the engine for a length unit of 1 pc and a mass unit of
1 M☉, `T* = √(PC³/GM)/(3.15×10¹³ s)` = 14.937 Myr; the time unit of a system
of mass `M` [M☉] and length unit `R` [pc] is this value times `√(R³/M)`.
"""
const _ENGINE_TIME_MYR = sqrt(_ENGINE_PC_CM^3 / _ENGINE_GM_CGS) / (_ENGINE_YR_S * 1.0e6)
