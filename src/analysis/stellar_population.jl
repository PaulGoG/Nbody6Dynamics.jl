# =============================================================================
# Stellar classes, the HR-plane population of an epoch, and the class census
# of a run
# =============================================================================
#
# The engine writes the stars of one epoch to two files (`hrplot.F`): single
# stars to sev.83 and the members of KS-regularised pairs to bev.82. A star
# moves between the two whenever its pair enters or leaves regularisation, so
# either file alone gives a population whose membership changes for reasons
# that have nothing to do with stellar evolution. Everything here works on the
# union.

"""
    StellarClass

A group of SSE/BSE stellar types K* (Hurley, Pols & Tout 2000,
doi:10.1046/j.1365-8711.2000.03426.x) treated as one population: `key`
identifies it, `label` is its legend entry, `kstar` the range of types it
holds. `luminous` states whether its members have a
photosphere that places them on the optical HR plane; neutron stars, black
holes and massless remnants do not (the engine assigns them cooling-curve or
placeholder luminosities), so they are counted and never drawn.
"""
struct StellarClass
    key::Symbol
    label::String
    kstar::UnitRange{Int}
    luminous::Bool
end

"""
    STELLAR_CLASSES

The stellar classes in evolutionary order. The grouping is the engine's own
("main stellar types" of `global_output.F`, used for its per-type Lagrangian
radii): AGB = K* 5–6, helium stars = 7–9, white dwarfs = 10–12, with the two
main-sequence types (K* 0 and 1, split at 0.7 M☉) taken together.
"""
const STELLAR_CLASSES = (
    StellarClass(:pre_main_sequence, "Pre-main sequence", -1:-1, true),
    StellarClass(:main_sequence, "Main sequence", 0:1, true),
    StellarClass(:hertzsprung_gap, "Hertzsprung gap", 2:2, true),
    StellarClass(:red_giant, "Red giant", 3:3, true),
    StellarClass(:core_helium_burning, "Core He burning", 4:4, true),
    StellarClass(:asymptotic_giant, "AGB", 5:6, true),
    StellarClass(:helium_star, "He star", 7:9, true),
    StellarClass(:white_dwarf, "White dwarf", 10:12, true),
    StellarClass(:neutron_star, "Neutron star", 13:13, false),
    StellarClass(:black_hole, "Black hole", 14:14, false),
    StellarClass(:massless_remnant, "Massless remnant", 15:15, false),
)

const _KSTAR_MIN = -1
const _KSTAR_MAX = 15

# Class index of every K*, addressed as k - _KSTAR_MIN + 1.
const _KSTAR_CLASS_INDEX = ntuple(
    i -> findfirst(c -> (i + _KSTAR_MIN - 1) in c.kstar, STELLAR_CLASSES)::Int,
    _KSTAR_MAX - _KSTAR_MIN + 1,
)

"""
    stellar_class_index(kstar) -> Int

Index into [`STELLAR_CLASSES`](@ref) of the class holding stellar type
`kstar`. Throws an `ArgumentError` outside the SSE/BSE range −1…15.
"""
function stellar_class_index(kstar::Integer)::Int
    _KSTAR_MIN ≤ kstar ≤ _KSTAR_MAX ||
        throw(ArgumentError("stellar type K* = $kstar is outside the SSE/BSE range -1…15"))
    return _KSTAR_CLASS_INDEX[kstar - _KSTAR_MIN + 1]
end

"""
    stellar_class(kstar) -> StellarClass

The [`StellarClass`](@ref) of stellar type `kstar`.

```julia
stellar_class(6).label      # "AGB"
stellar_class(14).luminous  # false
```
"""
stellar_class(kstar::Integer)::StellarClass = STELLAR_CLASSES[stellar_class_index(kstar)]

# Values below these are SSE/BSE placeholders for an undefined luminosity or
# temperature, not photospheres.
const _HR_MIN_LOG_L = -5.0
const _HR_MIN_LOG_TEFF = 3.0

"""
    HRPopulation

The stars of one epoch that lie on the HR plane, as columns: `log_teff`
(log₁₀ T_eff/K), `log_luminosity` (log₁₀ L/L☉), `class` (index into
[`STELLAR_CLASSES`](@ref)) and `binary_member` (the star is a member of a
KS-regularised pair). `time_myr` is the epoch.
"""
struct HRPopulation
    time_myr::Float64
    log_teff::Vector{Float64}
    log_luminosity::Vector{Float64}
    class::Vector{Int}
    binary_member::Vector{Bool}
end

Base.length(p::HRPopulation) = length(p.class)

"""
    hr_population(sev, bev = nothing) -> HRPopulation

The HR-plane population of one epoch: the single stars of `sev` together
with both members of every pair of `bev`, the latter flagged as binary
members. Stars of a non-luminous class ([`StellarClass`](@ref)) and records
carrying placeholder values (log L ≤ −5 or log T_eff ≤ 3) are left out.
"""
function hr_population(
    sev::StellarEvolutionSnapshot,
    bev::Union{Nothing,BinaryEvolutionSnapshot} = nothing,
)::HRPopulation
    n = length(sev.records) + (bev === nothing ? 0 : 2 * length(bev.records))
    pop = HRPopulation(sev.time_myr, Float64[], Float64[], Int[], Bool[])
    foreach(v -> sizehint!(v, n), (pop.log_teff, pop.log_luminosity, pop.class, pop.binary_member))
    for r in sev.records
        _push_star!(pop, r.stellar_type, r.log_teff, r.log_luminosity, false)
    end
    if bev !== nothing
        for r in bev.records
            _push_star!(pop, r.stellar_type1, r.log_teff1, r.log_luminosity1, true)
            _push_star!(pop, r.stellar_type2, r.log_teff2, r.log_luminosity2, true)
        end
    end
    return pop
end

function _push_star!(pop::HRPopulation, kstar::Integer, log_teff, log_l, member::Bool)
    k = stellar_class_index(kstar)
    STELLAR_CLASSES[k].luminous || return nothing
    (log_l > _HR_MIN_LOG_L && log_teff > _HR_MIN_LOG_TEFF) || return nothing
    push!(pop.log_teff, log_teff)
    push!(pop.log_luminosity, log_l)
    push!(pop.class, k)
    push!(pop.binary_member, member)
    return nothing
end

"""
    hr_populations(sevs, bevs = BinaryEvolutionSnapshot[]) -> Vector{HRPopulation}

[`hr_population`](@ref) of every epoch of `sevs`, each paired with the
binary snapshot of the same epoch when `bevs` holds one.
"""
function hr_populations(
    sevs::AbstractVector{StellarEvolutionSnapshot},
    bevs::AbstractVector{BinaryEvolutionSnapshot} = BinaryEvolutionSnapshot[],
)::Vector{HRPopulation}
    return [hr_population(s, _matching_bev(bevs, sevs, i)) for (i, s) in enumerate(sevs)]
end

# The engine writes both files in the same call, so their header times agree
# exactly. Equal-length lists pair by position first: the header carries one
# decimal, and epochs closer than 0.1 Myr share a time.
function _matching_bev(bevs, sevs, i::Int)
    isempty(bevs) && return nothing
    t = sevs[i].time_myr
    length(bevs) == length(sevs) && bevs[i].time_myr == t && return bevs[i]
    j = findfirst(b -> b.time_myr == t, bevs)
    return j === nothing ? nothing : bevs[j]
end

"""
    StellarCensus

Number of stars of every class at every stellar-evolution epoch of a run.
`single[i, k]` counts the single stars of class `k`
([`STELLAR_CLASSES`](@ref)) at `time_myr[i]`, `binary_member[i, k]` the
members of KS-regularised pairs. Non-luminous classes are included: this is
where neutron stars and black holes are accounted for.
"""
struct StellarCensus
    time_myr::Vector{Float64}
    single::Matrix{Int}
    binary_member::Matrix{Int}
end

"""
    stellar_census(sevs, bevs = BinaryEvolutionSnapshot[]) -> StellarCensus

Class census of a run from its stellar-evolution snapshots and, when given,
the binary snapshots of the same epochs.

```julia
census = stellar_census(results[:stellar_evo], results[:binary_evo])
class_counts(census)[end, stellar_class_index(14)]   # black holes at the last epoch
```
"""
function stellar_census(
    sevs::AbstractVector{StellarEvolutionSnapshot},
    bevs::AbstractVector{BinaryEvolutionSnapshot} = BinaryEvolutionSnapshot[],
)::StellarCensus
    nk = length(STELLAR_CLASSES)
    single = zeros(Int, length(sevs), nk)
    member = zeros(Int, length(sevs), nk)
    for (i, s) in enumerate(sevs)
        for r in s.records
            single[i, stellar_class_index(r.stellar_type)] += 1
        end
        b = _matching_bev(bevs, sevs, i)
        b === nothing && continue
        for r in b.records
            member[i, stellar_class_index(r.stellar_type1)] += 1
            member[i, stellar_class_index(r.stellar_type2)] += 1
        end
    end
    return StellarCensus([s.time_myr for s in sevs], single, member)
end

"""
    class_counts(census) -> Matrix{Int}

Stars per epoch (rows) and class (columns), single stars and binary members
together.
"""
class_counts(census::StellarCensus)::Matrix{Int} = census.single .+ census.binary_member

"""
    classes_present(census) -> Vector{Int}

Indices of the classes with at least one star at any epoch, in the order of
[`STELLAR_CLASSES`](@ref). Figures of a run take their legend from this set,
so every figure of the run lists the same classes whatever its epoch shows.
"""
function classes_present(census::StellarCensus)::Vector{Int}
    counts = class_counts(census)
    return [k for k in axes(counts, 2) if any(>(0), view(counts, :, k))]
end

"""
    write_stellar_census(path, census) -> String

Write `census` as CSV: `time_myr`, then `<class>_single` and
`<class>_binary` for every class. An existing file is backed up, never
overwritten. Returns `path`.
"""
function write_stellar_census(path::AbstractString, census::StellarCensus)::String
    _backup_existing(path)
    open(path, "w") do io
        header = ["time_myr"]
        for c in STELLAR_CLASSES
            push!(header, "$(c.key)_single", "$(c.key)_binary")
        end
        println(io, join(header, ','))
        for i in eachindex(census.time_myr)
            row = [string(census.time_myr[i])]
            for k in eachindex(STELLAR_CLASSES)
                push!(row, string(census.single[i, k]), string(census.binary_member[i, k]))
            end
            println(io, join(row, ','))
        end
    end
    return String(path)
end
