# =============================================================================
# IMBH diagnostics: trajectory about the stellar density centre, Brownian
# wandering, bound companions and the event list of a seeded IMBH, and the
# two-body elements of every pair of seeded IMBHs
# =============================================================================
#
# The IMBH is found by its name at every snapshot; the stars are every other
# body. Physical units throughout: lengths × RBAR [pc], velocities × VSTAR
# [km s⁻¹], masses × ZMBAR [M☉] (the total-mass scale of the header), times
# × TSCALE [Myr]; G = _G_PC_KMS2_MSUN [pc km² s⁻² M☉⁻¹].

"""Astronomical units per parsec."""
const _IMBH_AU_PER_PC = 206264.806

"""
Kinds of [`GWEvent`](@ref) that record a gravitational-wave merger: the two
merger triggers `GR KSINT COLL:` and `COLL DUE TO GW UNPERT`.
"""
const _IMBH_GW_MERGER_KINDS = (:gw_ksint_coll, :gw_unpert_coll)

"""Row type of [`imbh_wandering`](@ref)."""
const _IMBHWanderingRow = NamedTuple{
    (
        :frame,
        :t_start_myr,
        :t_end_myr,
        :n,
        :x2_mean_pc2,
        :r2_mean_pc2,
        :v2_mean_km2s2,
        :r_rms_pc,
        :kurtosis,
    ),
    Tuple{Symbol,Float64,Float64,Int,Float64,Float64,Float64,Float64,Float64},
}

"""Row type of [`imbh_events`](@ref)."""
const _IMBHEventRow = NamedTuple{
    (:kind, :time_myr, :partner_name, :partner_kstar, :partner_mass_msun),
    Tuple{Symbol,Float64,Int,Int,Float64},
}

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

"""
    _imbh_names(path) -> Vector{Int}

Body names of the seeded IMBHs from the IC record `merger_ic.toml` at
`path`: the vector `imbh.names` when present (its positive integer entries,
in record order), else `[imbh.name]` when that scalar is a positive integer,
else empty. Empty as well when the file or the `[imbh]` table is absent or
the file is not valid TOML.
"""
function _imbh_names(path::AbstractString)::Vector{Int}
    isfile(path) || return Int[]
    meta = try
        TOML.parsefile(path)
    catch e
        e isa Base.TOML.ParserError || rethrow()
        return Int[]
    end
    table = get(meta, "imbh", nothing)
    table isa AbstractDict || return Int[]
    vs = get(table, "names", nothing)
    vs isa AbstractVector && return Int[Int(v) for v in vs if v isa Integer && v > 0]
    v = get(table, "name", nothing)
    return v isa Integer && v > 0 ? [Int(v)] : Int[]
end

"""
    _imbh_name(path) -> Int

Body name of the first seeded IMBH of the IC record `merger_ic.toml` at
`path`, the first entry of `_imbh_names`; 0 when there is none.
"""
function _imbh_name(path::AbstractString)::Int
    names = _imbh_names(path)
    return isempty(names) ? 0 : names[1]
end

"""Column of the body `name` in `snap`, 0 when absent."""
_imbh_index(snap::Snapshot, name::Integer)::Int = something(findfirst(==(name), snap.name), 0)

"""Column of the KS partner of star `i` (`i ≤ 2 npairs`), 0 when `i` is no KS component."""
function _imbh_partner_index(snap::Snapshot, i::Int)::Int
    i ≤ 2 * snap.npairs || return 0
    return isodd(i) ? i + 1 : i - 1
end

# ---------------------------------------------------------------------------
# Trajectory
# ---------------------------------------------------------------------------

"""
    IMBHState

State of the IMBH at one snapshot, from [`imbh_trajectory`](@ref). The
stars are every body of the snapshot but the IMBH.

| field | meaning | unit |
|---|---|---|
| `time_nb`, `time_myr` | snapshot time | NB, Myr |
| `index` | column of the IMBH in the snapshot | — |
| `mass_msun` | IMBH mass | M☉ |
| `r_pc` | distance to the stellar density centre (shrinking sphere) | pc |
| `x_pc` | offset from the stellar density centre | pc |
| `v_kms` | velocity relative to the stellar centre-of-mass velocity | km s⁻¹ |
| `r_cm_pc` | distance to the stellar centre of mass | pc |
| `x_cm_pc` | offset from the stellar centre of mass | pc |
| `r_core_pc` | Casertano–Hut core radius of the stars ([`core_radius`](@ref)) from their sixth-neighbour mass densities among themselves, the IMBH excluded; `NaN` with ≤ 6 stars | pc |
| `r_half_pc` | half-mass radius of the stars about the density centre | pc |
| `in_pair` | the IMBH is a KS component (`index ≤ 2 npairs`) | — |
| `partner_name` | name of its KS partner, 0 otherwise | — |
"""
Base.@kwdef struct IMBHState
    time_nb::Float64
    time_myr::Float64
    index::Int
    mass_msun::Float64
    r_pc::Float64
    x_pc::NTuple{3,Float64}
    v_kms::NTuple{3,Float64}
    r_cm_pc::Float64
    x_cm_pc::NTuple{3,Float64} = (NaN, NaN, NaN)
    r_core_pc::Float64
    r_half_pc::Float64
    in_pair::Bool
    partner_name::Int
end

"""[`IMBHState`](@ref) of the body in column `i` of `snap`."""
function _imbh_snapshot_state(snap::Snapshot, i::Int; centre_kwargs...)::IMBHState
    h = snap.header
    L = rbar(h)
    V = vstar(h)
    Mu = zmbar(h)
    N = nparticles(snap)
    stars = [j for j in 1:N if j != i]
    isempty(stars) && throw(
        ArgumentError("imbh_trajectory: the snapshot at t = $(time_nb(h)) holds the IMBH alone"),
    )
    pos = Float64.(snap.pos[:, stars])
    vel = Float64.(snap.vel[:, stars])
    mass = Float64.(snap.mass[stars])
    M = sum(mass)
    M > 0 || throw(
        ArgumentError(
            "imbh_trajectory: the stars of the snapshot at t = $(time_nb(h)) are massless",
        ),
    )
    c = shrinking_sphere_centre(pos, mass; centre_kwargs...).centre
    x_cm = zeros(3)
    v_cm = zeros(3)
    for j in eachindex(mass), k in 1:3
        x_cm[k] += mass[j] * pos[k, j]
        v_cm[k] += mass[j] * vel[k, j]
    end
    x_cm ./= M
    v_cm ./= M
    xb = Float64.(snap.pos[:, i])
    vb = Float64.(snap.vel[:, i])
    x = ntuple(k -> (xb[k] - c[k]) * L, 3)
    v = ntuple(k -> (vb[k] - v_cm[k]) * V, 3)
    x_about_cm = ntuple(k -> (xb[k] - x_cm[k]) * L, 3)
    r_cm = sqrt(x_about_cm[1]^2 + x_about_cm[2]^2 + x_about_cm[3]^2)
    r_half = _lagrangian_radii(pos, mass, c, (0.5,))[1] * L
    # Densities of the stars among themselves: the engine's column counts the IMBH
    # among the neighbours and lets the cusp around it dominate the weights
    r_core = NaN
    if length(stars) > _DENSITY_NEIGHBOURS
        ρ = local_density(pos, mass; k = _DENSITY_NEIGHBOURS)
        sum(ρ) > 0 && (r_core = core_radius(pos, ρ).r_core * L)
    end
    j = _imbh_partner_index(snap, i)
    return IMBHState(
        time_nb(h),
        time_myr(h),
        i,
        Float64(snap.mass[i]) * Mu,
        sqrt(x[1]^2 + x[2]^2 + x[3]^2),
        x,
        v,
        r_cm,
        x_about_cm,
        r_core,
        r_half,
        j > 0,
        j > 0 ? Int(snap.name[j]) : 0,
    )
end

"""
    imbh_trajectory(snaps::AbstractVector{Snapshot}, name::Integer; centre_kwargs...)
        -> Vector{IMBHState}

State of the IMBH of body name `name` at every snapshot that holds it, one
[`IMBHState`](@ref) each, about the stars (every other body). The density
centre is [`shrinking_sphere_centre`](@ref) of the stars with
`centre_kwargs` (the centre keys of [`PostprocessConfig`](@ref)), the
reference velocity the mass-weighted mean velocity of the stars; the
half-mass radius of the stars is measured about the density centre and
their core radius by [`core_radius`](@ref) from the snapshot's density
column (`NaN` without one). A KS component (`index ≤ 2 npairs`) has its
partner at `index + 1` for odd and `index − 1` for even `index`.

Snapshots without the IMBH (escaped, merged, or inside a chain or merged
hierarchy, whose components `conf.3` does not list) are skipped and counted
in one `@info`. `ArgumentError` for an empty `snaps`.

# Example
```julia
traj = imbh_trajectory(read_all_conf3(dir), imbh_name)
r_rms = sqrt(sum(s -> s.r_pc^2, traj) / length(traj))
```
"""
function imbh_trajectory(
    snaps::AbstractVector{Snapshot},
    name::Integer;
    centre_kwargs...,
)::Vector{IMBHState}
    isempty(snaps) && throw(ArgumentError("imbh_trajectory: no snapshots"))
    out = IMBHState[]
    n_skipped = 0
    for snap in snaps
        i = _imbh_index(snap, name)
        if i == 0
            n_skipped += 1
            continue
        end
        push!(out, _imbh_snapshot_state(snap, i; centre_kwargs...))
    end
    n_skipped > 0 &&
        @info "IMBH trajectory: $n_skipped of $(length(snaps)) snapshots do not hold body $name"
    return out
end

# ---------------------------------------------------------------------------
# Wandering and the Brownian predictions
# ---------------------------------------------------------------------------

"""Wandering statistics of the states `states` over the window `[t_start, t_end]`."""
function _imbh_window_row(
    states::AbstractVector{IMBHState},
    t_start::Float64,
    t_end::Float64,
    frame::Symbol,
)::_IMBHWanderingRow
    n = length(states)
    s2 = 0.0
    s4 = 0.0
    r2 = 0.0
    v2 = 0.0
    for s in states
        x_s = frame === :centre_of_mass ? s.x_cm_pc : s.x_pc
        for k in 1:3
            x2 = x_s[k]^2
            s2 += x2
            s4 += x2^2
        end
        r2 += x_s[1]^2 + x_s[2]^2 + x_s[3]^2
        v2 += s.v_kms[1]^2 + s.v_kms[2]^2 + s.v_kms[3]^2
    end
    x2_mean = s2 / (3n)
    x4_mean = s4 / (3n)
    r2_mean = r2 / n
    return (
        frame = frame,
        t_start_myr = t_start,
        t_end_myr = t_end,
        n = n,
        x2_mean_pc2 = x2_mean,
        r2_mean_pc2 = r2_mean,
        v2_mean_km2s2 = v2 / n,
        r_rms_pc = sqrt(r2_mean),
        kurtosis = x4_mean / x2_mean^2,
    )
end

"""
    imbh_wandering(traj::AbstractVector{IMBHState}; window_myr::Real = 0.0,
                   frame::Symbol = :density_centre) -> Vector{NamedTuple}

Wandering statistics of the IMBH offsets from the stellar density centre
(`frame = :density_centre`, the shrinking-sphere centre, `x_pc`) or from the
stellar centre of mass (`frame = :centre_of_mass`, `x_cm_pc`, the frame of
Chatterjee, Hernquist & Loeb 2002), one row per time window: the whole span for `window_myr ≤ 0`, else
consecutive windows `[t₀ + k w, t₀ + (k + 1) w)` of length `w = window_myr`
from the first time `t₀`, empty windows left out; the last window takes the
remainder up to the last time and is kept when it holds at least three
states. States with a non-finite time are ignored. Rows
`(frame, t_start_myr, t_end_myr, n, x2_mean_pc2, r2_mean_pc2, v2_mean_km2s2,
r_rms_pc, kurtosis)`:

- `x2_mean_pc2`: ⟨x²⟩ per coordinate, the squared components pooled over the
  three axes [pc²], the quantity of [`brownian_rms_plummer`](@ref) and
  [`brownian_rms_core`](@ref) squared;
- `r2_mean_pc2 = ⟨r²⟩` [pc²] (≈ 3⟨x²⟩ for isotropic wandering) and
  `r_rms_pc = √⟨r²⟩`;
- `v2_mean_km2s2 = ⟨|v|²⟩` relative to the stellar centre-of-mass velocity
  [km² s⁻²];
- `kurtosis = ⟨x⁴⟩/⟨x²⟩²` pooled over the components, 3 for a Gaussian.

The moments are taken about zero (the chosen centre), not about the sample
mean. Returns an empty vector for an empty trajectory. When the wandering
amplitude is below the radius of the final shrinking sphere (≈ the radius
enclosing `n_floor` stars), the density centre follows the cusp bound to the
IMBH and the density-centre frame underestimates the wandering; the
centre-of-mass frame is the one to compare with the Brownian predictions.

# Example
```julia
w = imbh_wandering(traj; window_myr = 50.0)
ratio = [sqrt(r.x2_mean_pc2) for r in w] ./ brownian_rms_core(σ_c, ρ_c, m_eff, M_bh)
```
"""
function imbh_wandering(
    traj::AbstractVector{IMBHState};
    window_myr::Real = 0.0,
    frame::Symbol = :density_centre,
)
    frame in (:density_centre, :centre_of_mass) || throw(
        ArgumentError(
            "imbh_wandering: frame must be :density_centre or :centre_of_mass; got :$frame",
        ),
    )
    rows = _IMBHWanderingRow[]
    states = sort!([s for s in traj if isfinite(s.time_myr)]; by = s -> s.time_myr)
    isempty(states) && return rows
    t0 = states[1].time_myr
    t1 = states[end].time_myr
    if window_myr ≤ 0
        push!(rows, _imbh_window_row(states, t0, t1, frame))
        return rows
    end
    w = Float64(window_myr)
    n_full = floor(Int, (t1 - t0) / w)
    j = 1
    for k in 0:n_full
        lo = t0 + k * w
        hi = t0 + (k + 1) * w
        first_j = j
        while j ≤ length(states) && (k == n_full || states[j].time_myr < hi)
            j += 1
        end
        sel = view(states, first_j:(j - 1))
        if k < n_full
            isempty(sel) || push!(rows, _imbh_window_row(sel, lo, hi, frame))
        else
            length(sel) ≥ 3 && push!(rows, _imbh_window_row(sel, lo, t1, frame))
        end
    end
    return rows
end

"""
    effective_mass(mass::AbstractVector{<:Real}) -> Float64

Effective perturber mass `m̃ = ⟨m²⟩/⟨m⟩ = Σm²/Σm` of a stellar population
[unit of `mass`], the mass that drives the Brownian motion of a massive body
in a field of stars with a mass spectrum (the diffusion coefficients scale
with `Σ n m²`, the dynamical friction with `Σ n m`). For equal masses it is
the stellar mass. `ArgumentError` for an empty vector or a non-positive total.

# Example
```julia
effective_mass([1.0, 1.0, 2.0])   # 1.5
```
"""
function effective_mass(mass::AbstractVector{<:Real})::Float64
    isempty(mass) && throw(ArgumentError("effective_mass: no masses"))
    M = sum(Float64, mass)
    M > 0 || throw(ArgumentError("effective_mass: the total mass must be > 0; got $M"))
    return sum(m -> Float64(m)^2, mass) / M
end

"""
    brownian_rms_plummer(a_pc, m_eff_msun, m_imbh_msun) -> Float64

Equilibrium rms offset per coordinate [pc] of a massive body `M•` wandering
at the centre of a Plummer sphere of scale length `a` [pc] made of stars of
effective mass `m̃` ([`effective_mass`](@ref)),

    ⟨x²⟩ = (2/9) a² m̃ / M•

(Chatterjee, Hernquist & Loeb 2002, DOI 10.1086/339174). Valid for
`m̃ ≪ M• ≪ M_cluster`; the offsets are those about the stellar density
centre. `ArgumentError` unless `a_pc`, `m_eff_msun` and `m_imbh_msun` are
positive.

# Example
```julia
brownian_rms_plummer(3π / 16, 8.0e-5, 0.01)   # ≈ 0.0248
```
"""
function brownian_rms_plummer(a_pc::Real, m_eff_msun::Real, m_imbh_msun::Real)::Float64
    (a_pc > 0 && m_eff_msun > 0 && m_imbh_msun > 0) || throw(
        ArgumentError(
            "brownian_rms_plummer: a, m_eff and M• must be > 0; got $a_pc, $m_eff_msun, $m_imbh_msun",
        ),
    )
    return sqrt(2 / 9 * a_pc^2 * m_eff_msun / m_imbh_msun)
end

"""
    brownian_rms_core(sigma_c_kms, rho_c_msun_pc3, m_eff_msun, m_imbh_msun) -> Float64

Equilibrium rms offset per coordinate [pc] of a massive body `M•` in a core
of constant density `ρ_c` [M☉ pc⁻³] and one-dimensional velocity dispersion
`σ_c` [km s⁻¹], where the body moves in the harmonic potential of the core,

    ⟨x²⟩ = σ_c² (m̃/M•) / ω₀²,   ω₀² = (4π/3) G ρ_c,

with `m̃` the effective stellar mass ([`effective_mass`](@ref)) and
`G = 4.30091×10⁻³ pc km² s⁻² M☉⁻¹` (Lin & Tremaine 1980,
DOI 10.1086/157803; Merritt 2005, DOI 10.1086/429964). The equipartition
`⟨v²⟩ = σ_c² m̃/M•` per coordinate and the harmonic restoring force give the
offset; it holds while the offset stays inside the core. `ArgumentError`
unless all arguments are positive.

# Example
```julia
brownian_rms_core(10.0, 1.0e4, 0.5, 1000.0)
```
"""
function brownian_rms_core(
    sigma_c_kms::Real,
    rho_c_msun_pc3::Real,
    m_eff_msun::Real,
    m_imbh_msun::Real,
)::Float64
    (sigma_c_kms > 0 && rho_c_msun_pc3 > 0 && m_eff_msun > 0 && m_imbh_msun > 0) || throw(
        ArgumentError(
            "brownian_rms_core: σ_c, ρ_c, m_eff and M• must be > 0; got " *
            "$sigma_c_kms, $rho_c_msun_pc3, $m_eff_msun, $m_imbh_msun",
        ),
    )
    ω₀² = (4π / 3) * _G_PC_KMS2_MSUN * rho_c_msun_pc3
    return sqrt(sigma_c_kms^2 * (m_eff_msun / m_imbh_msun) / ω₀²)
end

"""
    _imbh_core_effective_mass(snaps, name, traj) -> Float64

[`effective_mass`](@ref) [M☉] of the stars inside `r_core_pc` of the
density centre at the first state of `traj`, in the snapshot of that state;
all stars when `r_core_pc` is `NaN` or fewer than 20 stars lie inside. `NaN`
for an empty trajectory.
"""
function _imbh_core_effective_mass(
    snaps::AbstractVector{Snapshot},
    name::Integer,
    traj::AbstractVector{IMBHState},
)::Float64
    isempty(traj) && return NaN
    s1 = traj[1]
    k = findfirst(s -> _imbh_index(s, name) > 0, snaps)
    k === nothing && return NaN
    snap = snaps[k]
    i = _imbh_index(snap, name)
    h = snap.header
    L = rbar(h)
    Mu = zmbar(h)
    centre = ntuple(c -> Float64(snap.pos[c, i]) * L - s1.x_pc[c], 3)
    all_m = Float64[]
    core_m = Float64[]
    for j in 1:nparticles(snap)
        j == i && continue
        m = Float64(snap.mass[j]) * Mu
        push!(all_m, m)
        if isfinite(s1.r_core_pc)
            d2 = sum(c -> (Float64(snap.pos[c, j]) * L - centre[c])^2, 1:3)
            d2 ≤ s1.r_core_pc^2 && push!(core_m, m)
        end
    end
    isempty(all_m) && return NaN
    return effective_mass(length(core_m) ≥ 20 ? core_m : all_m)
end

# ---------------------------------------------------------------------------
# Bound companions
# ---------------------------------------------------------------------------

"""
    IMBHCompanion

Most bound companion of the IMBH at one snapshot, from
[`imbh_companions`](@ref). Without an accepted companion (no bound
candidate, or the most bound one perturbed beyond `gamma_max`)
`companion_name` is 0 and the mass, separation, elements, energy and
`t_gw_myr` are `NaN`; `gamma` keeps the perturbation of a rejected
candidate and is `NaN` without a bound candidate.

| field | meaning | unit |
|---|---|---|
| `time_myr` | snapshot time | Myr |
| `companion_name` | name of the companion, 0 when none is accepted | — |
| `companion_mass_msun` | companion mass `m` | M☉ |
| `separation_pc` | separation `r` | pc |
| `semi_major_axis_au` | `a = −G (M• + m)/(2ε)` | au |
| `eccentricity` | `e = √(1 − |r × v|²/(G (M• + m) a))` | — |
| `specific_energy_km2s2` | `ε = ½ v_rel² − G (M• + m)/r` | km² s⁻² |
| `in_pair` | the companion is the IMBH's KS partner | — |
| `t_gw_myr` | Peters coalescence time [`gw_inspiral_time`](@ref) | Myr |
| `gamma` | perturbation `γ` of the most bound candidate by the other candidates | — |
"""
Base.@kwdef struct IMBHCompanion
    time_myr::Float64
    companion_name::Int
    companion_mass_msun::Float64
    separation_pc::Float64
    semi_major_axis_au::Float64
    eccentricity::Float64
    specific_energy_km2s2::Float64
    in_pair::Bool
    t_gw_myr::Float64
    gamma::Float64 = NaN
end

"""
    imbh_companions(snaps::AbstractVector{Snapshot}, name::Integer; k::Integer = 10,
                    gamma_max::Real = Inf) -> Vector{IMBHCompanion}

Most bound companion of the IMBH of body name `name` at every snapshot that
holds it, one [`IMBHCompanion`](@ref) each. The candidates are the `k`
nearest stars to the IMBH and, when the IMBH is a KS component, its KS
partner; for each the two-body specific energy
`ε = ½ |v_j − v•|² − G (M• + m_j)/|x_j − x•|` is formed in physical units,
and the most bound candidate `j` is the one of most negative `ε < 0`, with
the elements `a = −G (M• + m)/(2ε)` and
`e = √max(0, 1 − |r × v|²/(G (M• + m) a))` and the Peters coalescence time
of the pair.

The most bound candidate, at separation `d` from the IMBH, is then tested
against the other candidates `p` (the KS partner included when it is one)
by the perturbation parameter

    γ = max_p 2 m_p d³ / ((M• + m_j) r_p³),

`r_p` the distance of `p` from the centre of mass of the pair (IMBH, `j`),
and `γ = 0` without another candidate. `γ` is the KS perturbation measure of
the engine: the tidal acceleration of the perturber over the pair's
internal acceleration, up to the factor 2 of a perturber at the pair's
distance. The candidate is accepted as companion when `γ < gamma_max` (in
the pipeline `PostprocessConfig.bound_gamma_max`); otherwise the row has
`companion_name = 0` and `NaN` elements as for no bound candidate, with the
rejected candidate's `γ` in `gamma`. A snapshot without a bound candidate
gives a row with `companion_name = 0` and `gamma = NaN`. Without the
perturbation limit (`gamma_max = Inf`) the two-body test over the ten
nearest stars acts as a cusp detector around a massive body: the stars
inside the IMBH's radius of influence are two-body bound by construction.

`conf.3` positions are single precision, so elements of pairs whose
separation lies below ≈ 10⁻⁵ of the cluster's distance from the coordinate
origin are noise (IMBH study §4.6). The companion's stellar type is not in
`conf.3`; join by name to `sev.83`/`bev.82` for it. `ArgumentError` for
`k < 0`.

# Example
```julia
comp = imbh_companions(read_all_conf3(dir), imbh_name; gamma_max = 0.1)
imbh_companion_summary(comp).fraction_bound
```
"""
function imbh_companions(
    snaps::AbstractVector{Snapshot},
    name::Integer;
    k::Integer = 10,
    gamma_max::Real = Inf,
)::Vector{IMBHCompanion}
    k ≥ 0 || throw(ArgumentError("imbh_companions: k must be ≥ 0; got $k"))
    γ_max = Float64(gamma_max)
    G = _G_PC_KMS2_MSUN
    rows = IMBHCompanion[]
    for snap in snaps
        i = _imbh_index(snap, name)
        i == 0 && continue
        h = snap.header
        L = rbar(h)
        V = vstar(h)
        Mu = zmbar(h)
        t = time_myr(h)
        M_bh = Float64(snap.mass[i]) * Mu
        xb = ntuple(c -> Float64(snap.pos[c, i]), 3)
        vb = ntuple(c -> Float64(snap.vel[c, i]), 3)
        others = [j for j in 1:nparticles(snap) if j != i]
        d2 = [sum(c -> (Float64(snap.pos[c, j]) - xb[c])^2, 1:3) for j in others]
        n_k = min(Int(k), length(others))
        candidates = n_k > 0 ? others[partialsortperm(d2, 1:n_k)] : Int[]
        partner = _imbh_partner_index(snap, i)
        partner > 0 && !(partner in candidates) && push!(candidates, partner)
        best = 0
        ε_best = 0.0
        r_best = (NaN, NaN, NaN)
        v_best = (NaN, NaN, NaN)
        m_best = NaN
        for j in candidates
            r = ntuple(c -> (Float64(snap.pos[c, j]) - xb[c]) * L, 3)
            v = ntuple(c -> (Float64(snap.vel[c, j]) - vb[c]) * V, 3)
            m = Float64(snap.mass[j]) * Mu
            d = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
            d > 0 || continue
            ε = 0.5 * (v[1]^2 + v[2]^2 + v[3]^2) - G * (M_bh + m) / d
            if ε < ε_best
                best = j
                ε_best = ε
                r_best = r
                v_best = v
                m_best = m
            end
        end
        if best == 0
            push!(rows, IMBHCompanion(t, 0, NaN, NaN, NaN, NaN, NaN, false, NaN, NaN))
            continue
        end
        # Perturbation by the other candidates, in N-body units (γ is dimensionless)
        jb = best   # a fixed copy: `best` is reassigned in the loop above and must not be captured
        mb = Float64(snap.mass[i])
        mj = Float64(snap.mass[jb])
        xj = ntuple(c -> Float64(snap.pos[c, jb]), 3)
        x_cm = ntuple(c -> (mb * xb[c] + mj * xj[c]) / (mb + mj), 3)
        d³ = sqrt(sum(c -> (xj[c] - xb[c])^2, 1:3))^3
        γ = 0.0
        for p in candidates
            p == jb && continue
            r_p = sqrt(sum(c -> (Float64(snap.pos[c, p]) - x_cm[c])^2, 1:3))
            γ = max(γ, 2 * Float64(snap.mass[p]) * d³ / ((mb + mj) * r_p^3))
        end
        if !(γ < γ_max)
            push!(rows, IMBHCompanion(t, 0, NaN, NaN, NaN, NaN, NaN, false, NaN, γ))
            continue
        end
        μ = G * (M_bh + m_best)
        a_pc = -μ / (2ε_best)
        r, v = r_best, v_best
        hx = r[2] * v[3] - r[3] * v[2]
        hy = r[3] * v[1] - r[1] * v[3]
        hz = r[1] * v[2] - r[2] * v[1]
        e = sqrt(max(0.0, 1 - (hx^2 + hy^2 + hz^2) / (μ * a_pc)))
        a_au = a_pc * _IMBH_AU_PER_PC
        t_gw =
            (M_bh > 0 && m_best > 0 && e < 1) ? gw_inspiral_time(M_bh, m_best, a_au, e) / 1.0e6 :
            NaN
        push!(
            rows,
            IMBHCompanion(
                t,
                Int(snap.name[best]),
                m_best,
                sqrt(r[1]^2 + r[2]^2 + r[3]^2),
                a_au,
                e,
                ε_best,
                best == partner,
                t_gw,
                γ,
            ),
        )
    end
    return rows
end

"""
    imbh_companion_summary(rows::AbstractVector{IMBHCompanion})
        -> NamedTuple{(:n_snapshots, :fraction_bound, :n_exchanges, :longest_name, :longest_span_myr)}

Summary of the companion series of [`imbh_companions`](@ref):
`n_snapshots` rows, `fraction_bound` the fraction of rows with a bound
companion (`NaN` without rows), `n_exchanges` the changes of
`companion_name` between consecutive bound rows (unbound rows in between
skipped), and the longest uninterrupted run of one companion — consecutive
rows with the same name, broken by an unbound row or another name —
`longest_name` and `longest_span_myr` (last minus first time of the run;
the earlier run on ties). Without a bound row `longest_name = 0` and
`longest_span_myr = NaN`.
"""
function imbh_companion_summary(rows::AbstractVector{IMBHCompanion})
    n = length(rows)
    bound = [r for r in rows if r.companion_name != 0]
    n_exchanges =
        count(j -> bound[j].companion_name != bound[j - 1].companion_name, 2:length(bound))
    longest_name = 0
    longest_span = NaN
    j = 1
    while j ≤ n
        nm = rows[j].companion_name
        last_j = j
        while last_j < n && nm != 0 && rows[last_j + 1].companion_name == nm
            last_j += 1
        end
        if nm != 0
            span = rows[last_j].time_myr - rows[j].time_myr
            if isnan(longest_span) || span > longest_span
                longest_name = nm
                longest_span = span
            end
        end
        j = last_j + 1
    end
    return (
        n_snapshots = n,
        fraction_bound = n > 0 ? length(bound) / n : NaN,
        n_exchanges = n_exchanges,
        longest_name = longest_name,
        longest_span_myr = longest_span,
    )
end

# ---------------------------------------------------------------------------
# IMBH pairs
# ---------------------------------------------------------------------------

"""Myr per pc/(km s⁻¹): kilometres per parsec over the seconds of a megayear."""
const _IMBH_MYR_PER_PC_KMS = 3.0856775814913673e13 / 3.15576e13

"""
    IMBHPair

Two-body state of a pair of seeded IMBHs at one snapshot, from
[`imbh_pairs`](@ref); `name1 < name2`. The elements, the period and
`t_gw_myr` are `NaN` for an unbound pair (`ε ≥ 0`).

| field | meaning | unit |
|---|---|---|
| `time_myr` | snapshot time | Myr |
| `name1`, `name2` | body names of the two IMBHs | — |
| `separation_pc` | separation `r` | pc |
| `specific_energy_km2s2` | `ε = ½ v_rel² − G (M₁ + M₂)/r` | km² s⁻² |
| `bound` | `ε < 0` | — |
| `semi_major_axis_au` | `a = −G (M₁ + M₂)/(2ε)` | au |
| `eccentricity` | `e = √(1 − |r × v|²/(G (M₁ + M₂) a))` | — |
| `period_myr` | Kepler period `2π √(a³/(G (M₁ + M₂)))` | Myr |
| `t_gw_myr` | Peters coalescence time [`gw_inspiral_time`](@ref) | Myr |
| `in_pair` | the two are the components of one KS pair | — |
| `r_infl_pc` | influence radius: the radius about the pair's centre of mass enclosing a stellar mass `2 (M₁ + M₂)` | pc |
| `sigma_infl_kms` | one-dimensional velocity dispersion of the stars inside `r_infl_pc` about their mean velocity | km s⁻¹ |
| `rho_infl_msun_pc3` | mean stellar density inside `r_infl_pc` | M☉ pc⁻³ |

The stars are all bodies but the seeded IMBHs; the three influence
quantities are `NaN` when the stars hold less than `2 (M₁ + M₂)`.
"""
Base.@kwdef struct IMBHPair
    time_myr::Float64
    name1::Int
    name2::Int
    separation_pc::Float64
    specific_energy_km2s2::Float64
    bound::Bool
    semi_major_axis_au::Float64
    eccentricity::Float64
    period_myr::Float64
    t_gw_myr::Float64
    in_pair::Bool
    r_infl_pc::Float64 = NaN
    sigma_infl_kms::Float64 = NaN
    rho_infl_msun_pc3::Float64 = NaN
end

"""
    _imbh_influence(snap, i, j, skip) -> (r_pc, sigma_kms, rho_msun_pc3)

Influence radius of the bodies in columns `i` and `j` of `snap`, the radius
about their centre of mass enclosing a stellar mass of twice theirs, with
the one-dimensional velocity dispersion and the mean density of the stars
inside it (Sesana & Khan 2015, DOI 10.1093/mnrasl/slv131, measure `σ` and
`ρ` of the hardening rate there). The stars are the columns not in `skip`;
`NaN` for all three when their mass falls short.
"""
function _imbh_influence(snap::Snapshot, i::Int, j::Int, skip::Set{Int})
    h = snap.header
    L = rbar(h)
    V = vstar(h)
    Mu = zmbar(h)
    mi = Float64(snap.mass[i])
    mj = Float64(snap.mass[j])
    xc = ntuple(c -> (mi * snap.pos[c, i] + mj * snap.pos[c, j]) / (mi + mj), 3)
    stars = [k for k in 1:nparticles(snap) if k ∉ skip]
    target = 2 * (mi + mj)
    sum(k -> Float64(snap.mass[k]), stars; init = 0.0) ≥ target || return (NaN, NaN, NaN)
    d2 = [sum(c -> (snap.pos[c, k] - xc[c])^2, 1:3) for k in stars]
    order = sortperm(d2)
    m_in = 0.0
    n_in = 0
    while m_in < target
        n_in += 1
        m_in += Float64(snap.mass[stars[order[n_in]]])
    end
    inner = stars[order[1:n_in]]
    v̄ = ntuple(c -> sum(k -> Float64(snap.mass[k]) * snap.vel[c, k], inner) / m_in, 3)
    s2 =
        sum(k -> Float64(snap.mass[k]) * sum(c -> (snap.vel[c, k] - v̄[c])^2, 1:3), inner) /
        (3 * m_in)
    r = sqrt(d2[order[n_in]]) * L
    return (r, sqrt(s2) * V, 3 * m_in * Mu / (4π * r^3))
end

"""
[`IMBHPair`](@ref) of the bodies in columns `i` (name `n1`) and `j` (name
`n2`) of `snap`; the columns `skip` (the seeded IMBHs) are no stars.
"""
function _imbh_pair_row(
    snap::Snapshot,
    i::Int,
    j::Int,
    n1::Int,
    n2::Int,
    skip::Set{Int} = Set((i, j)),
)::IMBHPair
    G = _G_PC_KMS2_MSUN
    h = snap.header
    L = rbar(h)
    V = vstar(h)
    Mu = zmbar(h)
    t = time_myr(h)
    M1 = Float64(snap.mass[i]) * Mu
    M2 = Float64(snap.mass[j]) * Mu
    r = ntuple(c -> (Float64(snap.pos[c, j]) - Float64(snap.pos[c, i])) * L, 3)
    v = ntuple(c -> (Float64(snap.vel[c, j]) - Float64(snap.vel[c, i])) * V, 3)
    d = sqrt(r[1]^2 + r[2]^2 + r[3]^2)
    μ = G * (M1 + M2)
    ε = d > 0 ? 0.5 * (v[1]^2 + v[2]^2 + v[3]^2) - μ / d : NaN
    in_pair = _imbh_partner_index(snap, i) == j
    r_inf, σ_inf, ρ_inf = _imbh_influence(snap, i, j, skip)
    ε < 0 ||
        return IMBHPair(t, n1, n2, d, ε, false, NaN, NaN, NaN, NaN, in_pair, r_inf, σ_inf, ρ_inf)
    a_pc = -μ / (2ε)
    hx = r[2] * v[3] - r[3] * v[2]
    hy = r[3] * v[1] - r[1] * v[3]
    hz = r[1] * v[2] - r[2] * v[1]
    e = sqrt(max(0.0, 1 - (hx^2 + hy^2 + hz^2) / (μ * a_pc)))
    a_au = a_pc * _IMBH_AU_PER_PC
    period = 2π * sqrt(a_pc^3 / μ) * _IMBH_MYR_PER_PC_KMS
    t_gw = (M1 > 0 && M2 > 0 && e < 1) ? gw_inspiral_time(M1, M2, a_au, e) / 1.0e6 : NaN
    return IMBHPair(t, n1, n2, d, ε, true, a_au, e, period, t_gw, in_pair, r_inf, σ_inf, ρ_inf)
end

"""
    imbh_pairs(snaps::AbstractVector{Snapshot}, names::AbstractVector{<:Integer})
        -> Vector{IMBHPair}

Two-body state of every unordered pair of distinct IMBH names (`name1 <
name2`) at every snapshot that holds both, one [`IMBHPair`](@ref) each, in
the physical units of [`imbh_companions`](@ref): the specific energy
`ε = ½ |v₂ − v₁|² − G (M₁ + M₂)/|x₂ − x₁|`, and for `ε < 0` the elements
`a = −G (M₁ + M₂)/(2ε)`, `e = √max(0, 1 − |r × v|²/(G (M₁ + M₂) a))`, the
Kepler period `2π √(a³/(G (M₁ + M₂)))` converted from pc/(km s⁻¹) to Myr,
and the Peters coalescence time. `in_pair` marks two IMBHs that are the
components `2k − 1` and `2k` (`k ≤ npairs`) of one KS pair.

The criterion is two-body: the stars about the pair are ignored, so a
bound pair may be strongly perturbed. Repeated names count once. Rows are
sorted by `(name1, name2, time_myr)`; fewer than two distinct names give an
empty vector.

# Example
```julia
rows = imbh_pairs(read_all_conf3(dir), [5001, 5002])
imbh_pair_summary(rows).t_first_bound_myr
```
"""
function imbh_pairs(
    snaps::AbstractVector{Snapshot},
    names::AbstractVector{<:Integer},
)::Vector{IMBHPair}
    ns = sort!(unique(Int[Int(n) for n in names]))
    rows = IMBHPair[]
    n_bh = length(ns)
    n_bh < 2 && return rows
    for snap in snaps
        idx = [_imbh_index(snap, n) for n in ns]
        skip = Set(k for k in idx if k > 0)
        for a in 1:(n_bh - 1), b in (a + 1):n_bh
            (idx[a] > 0 && idx[b] > 0) || continue
            push!(rows, _imbh_pair_row(snap, idx[a], idx[b], ns[a], ns[b], skip))
        end
    end
    sort!(rows; by = r -> (r.name1, r.name2, r.time_myr))
    return rows
end

"""
    imbh_pair_summary(rows::AbstractVector{IMBHPair})
        -> NamedTuple{(:n_snapshots, :t_first_bound_myr, :t_last_bound_myr, :fraction_bound,
                       :a_min_au, :e_at_a_min)}

Summary of the rows of one IMBH pair from [`imbh_pairs`](@ref):
`n_snapshots` rows, the first and last time of a bound row
(`t_first_bound_myr`, `t_last_bound_myr`), the fraction of bound rows
`fraction_bound`, the smallest semi-major axis of a bound row `a_min_au`
[au] and the eccentricity of that row `e_at_a_min`. Without a bound row the
times, `a_min_au` and `e_at_a_min` are `NaN`; empty rows give
`n_snapshots = 0` and `fraction_bound = 0`. `ArgumentError` when the rows
mix several pairs.

# Example
```julia
rows = imbh_pairs(snaps, [5001, 5002])
s = imbh_pair_summary(rows)
s.fraction_bound, s.a_min_au
```
"""
function imbh_pair_summary(rows::AbstractVector{IMBHPair})
    n = length(rows)
    if n > 0
        p = (rows[1].name1, rows[1].name2)
        all(r -> (r.name1, r.name2) == p, rows) || throw(
            ArgumentError(
                "imbh_pair_summary: the rows mix several pairs; select one (name1, name2) first",
            ),
        )
    end
    bound = [r for r in rows if r.bound]
    t_first = NaN
    t_last = NaN
    a_min = NaN
    e_min = NaN
    if !isempty(bound)
        t_first = minimum(r -> r.time_myr, bound)
        t_last = maximum(r -> r.time_myr, bound)
        k = argmin([r.semi_major_axis_au for r in bound])
        a_min = bound[k].semi_major_axis_au
        e_min = bound[k].eccentricity
    end
    return (
        n_snapshots = n,
        t_first_bound_myr = t_first,
        t_last_bound_myr = t_last,
        fraction_bound = n > 0 ? length(bound) / n : 0.0,
        a_min_au = a_min,
        e_at_a_min = e_min,
    )
end

# ---------------------------------------------------------------------------
# Hardening of IMBH pairs
# ---------------------------------------------------------------------------

"""Parsec per megayear in one km s⁻¹."""
const _IMBH_PC_MYR_PER_KMS = 1.0227121650537077

"""Row type of [`imbh_pair_hardening`](@ref)."""
const _IMBHHardeningRow = NamedTuple{
    (
        :name1,
        :name2,
        :t_start_myr,
        :t_end_myr,
        :n,
        :a_mean_au,
        :s_au_per_myr,
        :s_err_au_per_myr,
        :k_de_dlninva,
        :sigma_infl_kms,
        :rho_infl_msun_pc3,
        :h_dimensionless,
    ),
    Tuple{Int,Int,Float64,Float64,Int,Float64,Float64,Float64,Float64,Float64,Float64,Float64},
}

"""Least-squares slope of `y` against `x`; `NaN` when `x` does not vary."""
function _imbh_ls_slope(x::AbstractVector{Float64}, y::AbstractVector{Float64})::Float64
    x̄ = sum(x) / length(x)
    ȳ = sum(y) / length(y)
    Sxx = sum(v -> (v - x̄)^2, x)
    Sxx > 0 || return NaN
    return sum((x[k] - x̄) * (y[k] - ȳ) for k in eachindex(x, y)) / Sxx
end

"""Hardening row of the bound rows `sel` of one pair over `[t_start, t_end]`."""
function _imbh_hardening_row(
    sel::AbstractVector{IMBHPair},
    t_start::Float64,
    t_end::Float64,
)::_IMBHHardeningRow
    n = length(sel)
    t = Float64[r.time_myr for r in sel]
    a = Float64[r.semi_major_axis_au for r in sel]
    inv_a = 1 ./ a
    s = _imbh_ls_slope(t, inv_a)
    s_err = NaN
    if n ≥ 4 && isfinite(s)
        t̄ = sum(t) / n
        Stt = sum(v -> (v - t̄)^2, t)
        b = sum(inv_a) / n - s * t̄
        res_var = sum((inv_a[k] - b - s * t[k])^2 for k in 1:n) / (n - 2)
        s_err = sqrt(res_var / Stt)
    end
    k_e = _imbh_ls_slope(log.(inv_a), Float64[r.eccentricity for r in sel])
    σs = Float64[r.sigma_infl_kms for r in sel if isfinite(r.sigma_infl_kms)]
    ρs = Float64[r.rho_infl_msun_pc3 for r in sel if isfinite(r.rho_infl_msun_pc3)]
    σ̄ = isempty(σs) ? NaN : sum(σs) / length(σs)
    ρ̄ = isempty(ρs) ? NaN : sum(ρs) / length(ρs)
    H = (isfinite(s) && σ̄ > 0 && ρ̄ > 0) ? hardening_rate_dimensionless(s, σ̄, ρ̄) : NaN
    return (
        name1 = sel[1].name1,
        name2 = sel[1].name2,
        t_start_myr = t_start,
        t_end_myr = t_end,
        n = n,
        a_mean_au = sum(a) / n,
        s_au_per_myr = s,
        s_err_au_per_myr = s_err,
        k_de_dlninva = k_e,
        sigma_infl_kms = σ̄,
        rho_infl_msun_pc3 = ρ̄,
        h_dimensionless = H,
    )
end

"""
    imbh_pair_hardening(rows::AbstractVector{IMBHPair}; window_myr::Real = 0.0)
        -> Vector{NamedTuple}

Hardening rate of one IMBH pair from its rows of [`imbh_pairs`](@ref), one
row per time window over the bound rows (unbound rows and rows with a
non-finite time are ignored): the whole span for `window_myr ≤ 0`, else
consecutive windows `[t₀ + k w, t₀ + (k + 1) w)` of length `w = window_myr`
from the first bound time `t₀`, the last window taking the remainder up to
the last bound time, as in [`imbh_wandering`](@ref). A window with fewer
than three bound rows is left out. Rows `(name1, name2, t_start_myr,
t_end_myr, n, a_mean_au, s_au_per_myr, s_err_au_per_myr, k_de_dlninva,
sigma_infl_kms, rho_infl_msun_pc3, h_dimensionless)`:

- `n` the bound rows of the window and `a_mean_au` their mean `a` [au];
- `s_au_per_myr`: the least-squares slope `s = d(1/a)/dt` of `1/a` [au⁻¹]
  against `t` [Myr], the hardening rate of Quinlan (1996,
  DOI 10.1016/S1384-1076(96)00003-6) and Sesana & Khan (2015,
  DOI 10.1093/mnrasl/slv131), whose dimensionless form `H = s σ/(G ρ)`
  ([`hardening_rate_dimensionless`](@ref)) is ≈ 15–20 for a hard binary;
- `s_err_au_per_myr`: its standard error
  `√(residual variance / Σ(t − t̄)²)`, the residual variance with `n − 2`
  degrees of freedom; `NaN` for fewer than four rows;
- `k_de_dlninva`: the least-squares slope `K = de/d ln(1/a)` of `e` against
  `ln(1/a)`, the eccentricity growth rate of Quinlan (1996); `NaN` when
  `1/a` does not vary;
- `sigma_infl_kms`, `rho_infl_msun_pc3`: the means over the window's rows of
  the dispersion and density inside the influence radius
  ([`IMBHPair`](@ref)), and `h_dimensionless` the rate `H = s σ/(G ρ)` at
  those means; `NaN` when the rows carry no influence quantities.

`ArgumentError` when the rows mix several pairs; an empty vector without
three bound rows.

# Example
```julia
rows = imbh_pairs(read_all_conf3(dir), [5001, 5002])
h = only(imbh_pair_hardening(rows))
h.h_dimensionless   # H at the influence radius
```
"""
function imbh_pair_hardening(rows::AbstractVector{IMBHPair}; window_myr::Real = 0.0)
    out = _IMBHHardeningRow[]
    isempty(rows) && return out
    p = (rows[1].name1, rows[1].name2)
    all(r -> (r.name1, r.name2) == p, rows) || throw(
        ArgumentError(
            "imbh_pair_hardening: the rows mix several pairs; select one (name1, name2) first",
        ),
    )
    bound = sort!([r for r in rows if r.bound && isfinite(r.time_myr)]; by = r -> r.time_myr)
    length(bound) ≥ 3 || return out
    t0 = bound[1].time_myr
    t1 = bound[end].time_myr
    if window_myr ≤ 0
        push!(out, _imbh_hardening_row(bound, t0, t1))
        return out
    end
    w = Float64(window_myr)
    n_full = floor(Int, (t1 - t0) / w)
    j = 1
    for k in 0:n_full
        lo = t0 + k * w
        hi = t0 + (k + 1) * w
        first_j = j
        while j ≤ length(bound) && (k == n_full || bound[j].time_myr < hi)
            j += 1
        end
        sel = view(bound, first_j:(j - 1))
        length(sel) ≥ 3 && push!(out, _imbh_hardening_row(sel, lo, k < n_full ? hi : t1))
    end
    return out
end

"""
    hardening_rate_dimensionless(s_au_per_myr, sigma_kms, rho_msun_pc3) -> Float64

Dimensionless hardening rate `H = s σ/(G ρ)` of a massive binary (Quinlan
1996, DOI 10.1016/S1384-1076(96)00003-6; Sesana & Khan 2015,
DOI 10.1093/mnrasl/slv131) from the rate `s = d(1/a)/dt` [au⁻¹ Myr⁻¹] of
[`imbh_pair_hardening`](@ref), the one-dimensional velocity dispersion `σ`
[km s⁻¹] and the stellar density `ρ` [M☉ pc⁻³] about the binary. `s` is
converted to pc⁻¹ Myr⁻¹ (× 206264.806 au pc⁻¹), `σ` to pc Myr⁻¹
(× 1.0227121650537077, one km s⁻¹ in pc Myr⁻¹), and
`G = 4.30091×10⁻³ × 1.0227121650537077²` pc³ M☉⁻¹ Myr⁻². A hard binary in
a fixed stellar background has `H ≈ 15–20`. `ArgumentError` unless `σ` and
`ρ` are positive.

# Example
```julia
hardening_rate_dimensionless(1.0e-3, 10.0, 1.0e4)   # ≈ 46.9
```
"""
function hardening_rate_dimensionless(
    s_au_per_myr::Real,
    sigma_kms::Real,
    rho_msun_pc3::Real,
)::Float64
    (sigma_kms > 0 && rho_msun_pc3 > 0) || throw(
        ArgumentError(
            "hardening_rate_dimensionless: σ and ρ must be > 0; got $sigma_kms, $rho_msun_pc3",
        ),
    )
    s = Float64(s_au_per_myr) * _IMBH_AU_PER_PC
    σ = Float64(sigma_kms) * _IMBH_PC_MYR_PER_KMS
    G = _G_PC_KMS2_MSUN * _IMBH_PC_MYR_PER_KMS^2
    return s * σ / (G * Float64(rho_msun_pc3))
end

# ---------------------------------------------------------------------------
# Events
# ---------------------------------------------------------------------------

"""
    imbh_events(name::Integer; collisions = CollisionRecord[], coalescences = CoalescenceRecord[],
                escapers = EscaperRecord[], gw_events = GWEvent[], tstar_myr::Real = NaN)
        -> Vector{NamedTuple{(:kind, :time_myr, :partner_name, :partner_kstar, :partner_mass_msun)}}

Events of the IMBH of body name `name`, from the engine's records:

| kind | source |
|---|---|
| `:tde` | a `coll.13` record ([`CollisionRecord`](@ref)) naming the IMBH whose partner has K* ≤ 12 (a star or white dwarf) |
| `:collision` | other `coll.13` records naming the IMBH |
| `:coalescence` | a `coal.24` record ([`CoalescenceRecord`](@ref)) naming the IMBH |
| `:gw_merger` | a [`GWEvent`](@ref) of the merger kinds (`:gw_ksint_coll`, `:gw_unpert_coll`) with `name1` or `name2` the IMBH |
| `:recoil` | a `:gw_recoil` [`GWEvent`](@ref) naming the IMBH |
| `:escape` | an `esc.11` record ([`EscaperRecord`](@ref)) of the IMBH |

Times [Myr]: `time_nb × tstar_myr` for the records timed in N-body units
(`NaN` when `tstar_myr` is `NaN`), the record's own `time_myr` where it has
one (`esc.11`, `:gw_recoil`). The partner fields hold the other body's name,
stellar type and mass before the event; they are 0, `KSTAR_ABSENT` and `NaN`
where they do not apply (`:escape`). The two GW merger lines of the engine
carry neither names nor a time, so `:gw_merger` rows appear only for records
that do. Sorted by time, `NaN` times last.

# Example
```julia
ev = imbh_events(imbh_name; collisions = read_collision_log(joinpath(dir, "coll.13")),
                 tstar_myr = extract_scaling(read_diagnostics(joinpath(dir, "out1000"))).tscale)
n_tde = count(r -> r.kind === :tde, ev)
```
"""
function imbh_events(
    name::Integer;
    collisions::AbstractVector{CollisionRecord} = CollisionRecord[],
    coalescences::AbstractVector{CoalescenceRecord} = CoalescenceRecord[],
    escapers::AbstractVector{EscaperRecord} = EscaperRecord[],
    gw_events::AbstractVector{GWEvent} = GWEvent[],
    tstar_myr::Real = NaN,
)::Vector{_IMBHEventRow}
    nm = Int(name)
    ts = Float64(tstar_myr)
    nb_myr(t) = isnan(ts) ? NaN : t * ts
    rows = _IMBHEventRow[]
    for c in collisions
        (c.name1 == nm || c.name2 == nm) || continue
        pname, pk, pm =
            c.name1 == nm ? (c.name2, c.kstar2, c.mass2_msun) : (c.name1, c.kstar1, c.mass1_msun)
        push!(
            rows,
            (
                kind = pk ≤ 12 ? :tde : :collision,
                time_myr = nb_myr(c.time_nb),
                partner_name = pname,
                partner_kstar = pk,
                partner_mass_msun = pm,
            ),
        )
    end
    for c in coalescences
        (c.name1 == nm || c.name2 == nm) || continue
        pname, pk, pm =
            c.name1 == nm ? (c.name2, c.kstar2, c.mass2_msun) : (c.name1, c.kstar1, c.mass1_msun)
        push!(
            rows,
            (
                kind = :coalescence,
                time_myr = nb_myr(c.time_nb),
                partner_name = pname,
                partner_kstar = pk,
                partner_mass_msun = pm,
            ),
        )
    end
    for g in gw_events
        kind =
            g.kind in _IMBH_GW_MERGER_KINDS ? :gw_merger : g.kind === :gw_recoil ? :recoil : :none
        kind === :none && continue
        (g.name1 == nm || g.name2 == nm) || continue
        pname, pk, pm =
            g.name1 == nm ? (g.name2, g.kstar2, g.m2_msun) : (g.name1, g.kstar1, g.m1_msun)
        t = isfinite(g.time_myr) ? g.time_myr : nb_myr(g.time_nb)
        push!(
            rows,
            (
                kind = kind,
                time_myr = t,
                partner_name = pname,
                partner_kstar = pk,
                partner_mass_msun = pm,
            ),
        )
    end
    for e in escapers
        e.name == nm || continue
        push!(
            rows,
            (
                kind = :escape,
                time_myr = e.time_myr,
                partner_name = 0,
                partner_kstar = KSTAR_ABSENT,
                partner_mass_msun = NaN,
            ),
        )
    end
    sort!(rows; by = r -> isnan(r.time_myr) ? (1, 0.0) : (0, r.time_myr))
    return rows
end

# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

"""Write a CSV file: one `#` metadata line, the header, one line of cells per row."""
function _imbh_write_csv(
    path::AbstractString,
    meta::AbstractString,
    header::AbstractVector{<:AbstractString},
    rows,
)::String
    _backup_existing(path)
    open(path, "w") do io
        println(io, "# ", meta)
        println(io, join(header, ","))
        for cells in rows
            println(io, join(map(_csv_precise, cells), ","))
        end
    end
    return String(path)
end

"""
Write the series of several IMBHs as one CSV: the `#` line lists the names
(`imbh_name=n₁;n₂;…`) before `meta_tail`, the header is `imbh_name`
followed by `header`, and every element `r` of every series contributes one
row `(name, cells(r)...)`, in the order of `series`.
"""
function _imbh_write_series(
    path::AbstractString,
    meta_tail::AbstractString,
    header::AbstractVector{<:AbstractString},
    series::AbstractVector{<:Tuple{Integer,Any}},
    cells,
)::String
    names = Int[Int(p[1]) for p in series]
    rows = ((Int(n), cells(r)...) for (n, s) in series for r in s)
    return _imbh_write_csv(
        path,
        "imbh_name=$(join(names, ";")) $meta_tail",
        vcat("imbh_name", collect(String, header)),
        rows,
    )
end

"""Columns of [`write_imbh_trajectory`](@ref) after `imbh_name`."""
const _IMBH_TRAJECTORY_HEADER = [
    "time_nb",
    "time_myr",
    "index",
    "mass_msun",
    "r_pc",
    "x_pc",
    "y_pc",
    "z_pc",
    "vx_kms",
    "vy_kms",
    "vz_kms",
    "r_cm_pc",
    "x_cm_pc",
    "y_cm_pc",
    "z_cm_pc",
    "r_core_pc",
    "r_half_pc",
    "in_pair",
    "partner_name",
]

"""Cells of one [`IMBHState`](@ref), the vectors split into components."""
_imbh_state_cells(s::IMBHState) = (
    s.time_nb,
    s.time_myr,
    s.index,
    s.mass_msun,
    s.r_pc,
    s.x_pc...,
    s.v_kms...,
    s.r_cm_pc,
    s.x_cm_pc...,
    s.r_core_pc,
    s.r_half_pc,
    s.in_pair,
    s.partner_name,
)

"""
    write_imbh_trajectory(path, series::AbstractVector{<:Tuple{Integer,Any}}) -> String
    write_imbh_trajectory(path, traj::AbstractVector{IMBHState}; name::Integer = 0) -> String

Write the trajectories of one or several IMBHs as CSV. `series` holds
`(name, traj)` pairs, `traj` a vector of [`IMBHState`](@ref) from
[`imbh_trajectory`](@ref); the single-series method writes `[(name, traj)]`.
A `#` line lists the IMBH names (`imbh_name=n₁;n₂;…`) and the total number
of states, then one row per state with the IMBH name first (`imbh_name`)
and the fields in declaration order, the vectors split into components
(`x_pc`, `y_pc`, `z_pc`, `vx_kms`, `vy_kms`, `vz_kms`); non-finite values
as empty cells. An existing file is backed up, never overwritten. Returns
the path.

# Example
```julia
write_imbh_trajectory("imbh_trajectory.csv", [(n, imbh_trajectory(snaps, n)) for n in names])
```
"""
function write_imbh_trajectory(
    path::AbstractString,
    series::AbstractVector{<:Tuple{Integer,Any}},
)::String
    n_states = sum(p -> length(p[2]), series; init = 0)
    return _imbh_write_series(
        path,
        "n_states=$n_states",
        _IMBH_TRAJECTORY_HEADER,
        series,
        _imbh_state_cells,
    )
end

function write_imbh_trajectory(
    path::AbstractString,
    traj::AbstractVector{IMBHState};
    name::Integer = 0,
)::String
    return write_imbh_trajectory(path, [(Int(name), traj)])
end

"""
    write_imbh_wandering(path, series::AbstractVector{<:Tuple{Integer,Any}}) -> String
    write_imbh_wandering(path, rows; name::Integer = 0) -> String

Write the wandering statistics of [`imbh_wandering`](@ref) of one or several
IMBHs as CSV. `series` holds `(name, rows)` pairs; the single-series method
writes `[(name, rows)]`. A `#` line lists the IMBH names, the total number
of windows and of states, then one row per window with the IMBH name first
(`imbh_name`) and the columns of the row tuple. An existing file is backed
up, never overwritten. Returns the path.
"""
function write_imbh_wandering(
    path::AbstractString,
    series::AbstractVector{<:Tuple{Integer,Any}},
)::String
    n_windows = sum(p -> length(p[2]), series; init = 0)
    n_states = sum(p -> sum(r -> r.n, p[2]; init = 0), series; init = 0)
    return _imbh_write_series(
        path,
        "n_windows=$n_windows n_states=$n_states",
        String.(collect(fieldnames(_IMBHWanderingRow))),
        series,
        Tuple,
    )
end

function write_imbh_wandering(path::AbstractString, rows::AbstractVector; name::Integer = 0)::String
    return write_imbh_wandering(path, [(Int(name), rows)])
end

"""
    write_imbh_companions(path, series::AbstractVector{<:Tuple{Integer,Any}}) -> String
    write_imbh_companions(path, rows::AbstractVector{IMBHCompanion}; name::Integer = 0) -> String

Write the companion series of [`imbh_companions`](@ref) of one or several
IMBHs as CSV. `series` holds `(name, rows)` pairs; the single-series method
writes `[(name, rows)]`. A `#` line lists the IMBH names, the total number
of snapshots and of bound rows, then one row per [`IMBHCompanion`](@ref)
with the IMBH name first (`imbh_name`) and the fields in declaration order.
An existing file is backed up, never overwritten. Returns the path.
"""
function write_imbh_companions(
    path::AbstractString,
    series::AbstractVector{<:Tuple{Integer,Any}},
)::String
    n_rows = sum(p -> length(p[2]), series; init = 0)
    n_bound = sum(p -> count(r -> r.companion_name != 0, p[2]), series; init = 0)
    n_fields = fieldcount(IMBHCompanion)
    return _imbh_write_series(
        path,
        "n_snapshots=$n_rows n_bound=$n_bound",
        String.(collect(fieldnames(IMBHCompanion))),
        series,
        r -> ntuple(f -> getfield(r, f), n_fields),
    )
end

function write_imbh_companions(
    path::AbstractString,
    rows::AbstractVector{IMBHCompanion};
    name::Integer = 0,
)::String
    return write_imbh_companions(path, [(Int(name), rows)])
end

"""
    write_imbh_events(path, series::AbstractVector{<:Tuple{Integer,Any}}) -> String
    write_imbh_events(path, rows; name::Integer = 0) -> String

Write the event lists of [`imbh_events`](@ref) of one or several IMBHs as
CSV. `series` holds `(name, rows)` pairs; the single-series method writes
`[(name, rows)]`. A `#` line lists the IMBH names and the total number of
events, then one row per event with the IMBH name first (`imbh_name`,
`kind`, `time_myr`, `partner_name`, `partner_kstar`, `partner_mass_msun`).
An existing file is backed up, never overwritten. Returns the path.
"""
function write_imbh_events(
    path::AbstractString,
    series::AbstractVector{<:Tuple{Integer,Any}},
)::String
    n_events = sum(p -> length(p[2]), series; init = 0)
    return _imbh_write_series(
        path,
        "n_events=$n_events",
        String.(collect(fieldnames(_IMBHEventRow))),
        series,
        Tuple,
    )
end

function write_imbh_events(path::AbstractString, rows::AbstractVector; name::Integer = 0)::String
    return write_imbh_events(path, [(Int(name), rows)])
end

"""
    write_imbh_pairs(path, rows::AbstractVector{IMBHPair}; names = Int[]) -> String

Write the pair states of [`imbh_pairs`](@ref) as CSV: a `#` line with the
IMBH names (`imbh_name=n₁;n₂;…`, the names met in `rows` when `names` is
empty), the number of distinct pairs and of rows, then one row per
[`IMBHPair`](@ref) with the fields in declaration order; non-finite values
as empty cells. An existing file is backed up, never overwritten. Returns
the path.

# Example
```julia
write_imbh_pairs("imbh_pairs.csv", imbh_pairs(snaps, names); names = names)
```
"""
function write_imbh_pairs(
    path::AbstractString,
    rows::AbstractVector{IMBHPair};
    names::AbstractVector{<:Integer} = Int[],
)::String
    listed = if isempty(names)
        sort!(unique(vcat(Int[r.name1 for r in rows], Int[r.name2 for r in rows])))
    else
        Int[Int(n) for n in names]
    end
    n_pairs = length(unique((r.name1, r.name2) for r in rows))
    n_fields = fieldcount(IMBHPair)
    return _imbh_write_csv(
        path,
        "imbh_name=$(join(listed, ";")) n_pairs=$n_pairs n_rows=$(length(rows))",
        String.(collect(fieldnames(IMBHPair))),
        (ntuple(f -> getfield(r, f), n_fields) for r in rows),
    )
end

"""
    write_imbh_pair_hardening(path, rows; names = Int[]) -> String

Write the hardening rows of [`imbh_pair_hardening`](@ref) of one or several
IMBH pairs as CSV: a `#` line with the IMBH names (`imbh_name=n₁;n₂;…`, the
names met in `rows` when `names` is empty), the number of distinct pairs
and of rows, then one row per window with the columns of the row tuple;
non-finite values as empty cells. An existing file is backed up, never
overwritten. Returns the path.

# Example
```julia
write_imbh_pair_hardening("imbh_pair_hardening.csv", imbh_pair_hardening(rows); names = names)
```
"""
function write_imbh_pair_hardening(
    path::AbstractString,
    rows::AbstractVector;
    names::AbstractVector{<:Integer} = Int[],
)::String
    listed = if isempty(names)
        sort!(unique(vcat(Int[r.name1 for r in rows], Int[r.name2 for r in rows])))
    else
        Int[Int(n) for n in names]
    end
    n_pairs = length(unique((r.name1, r.name2) for r in rows))
    return _imbh_write_csv(
        path,
        "imbh_name=$(join(listed, ";")) n_pairs=$n_pairs n_rows=$(length(rows))",
        String.(collect(fieldnames(_IMBHHardeningRow))),
        (Tuple(r) for r in rows),
    )
end
