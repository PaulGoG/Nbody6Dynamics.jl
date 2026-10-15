# =============================================================================
# Per-cluster structure from snapshots by original membership
# =============================================================================
# The engine keeps one density centre and one scale radius, so before
# coalescence its lagr.7, RC, and RTIDE describe the configuration, not the
# member clusters. These routines measure every initial cluster about its own
# centre, from the members it still binds.

"""
    ClusterStructure

Per-cluster structural quantities at every snapshot, all by original
membership (`cluster_ranges` of `merger_ic.toml`) and in N-body units.
Arrays are `n_clusters × n_times`; `NaN` marks clusters with fewer than
`_MIN_MEMBERS` bound members.

# Fields
- `time`: snapshot times [NB]
- `n_members`: members present in the snapshot
- `n_bound`: members bound to the cluster (negative energy in its own frame)
- `bound_mass_fraction`: bound mass over present member mass
- `centre`: `3 × n_clusters × n_times` shrinking-sphere centre of the bound members
- `r_lagr`: `3 × n_clusters × n_times` radii enclosing 10 %, 50 %, 90 % of the bound mass
- `sigma_1d`: one-dimensional velocity dispersion of the bound members [NB]
- `q_virial`: `T/|W|` of the bound members (COM velocity subtracted)
"""
struct ClusterStructure
    time::Vector{Float64}
    n_members::Matrix{Int}
    n_bound::Matrix{Int}
    bound_mass_fraction::Matrix{Float64}
    centre::Array{Float64,3}
    r_lagr::Array{Float64,3}
    sigma_1d::Matrix{Float64}
    q_virial::Matrix{Float64}
end

const _LAGR_FRACTIONS_MEMBER = (0.1, 0.5, 0.9)
const _MIN_MEMBERS = 10

"""
    _member_indices(snap, members) -> Vector{Int}

Indices in `snap` of the particles whose original name (body index in
`dat.10`) belongs to `members`, a contiguous range or, with primordial
binaries, a vector of body indices.
"""
function _member_indices(snap::Snapshot, members::UnitRange{Int})::Vector{Int}
    return [i for (i, n) in enumerate(snap.name) if first(members) ≤ n ≤ last(members)]
end
function _member_indices(snap::Snapshot, members::AbstractVector{Int})::Vector{Int}
    isempty(members) && return Int[]
    flag = falses(maximum(members))
    flag[members] .= true
    return [i for (i, n) in enumerate(snap.name) if 1 ≤ n ≤ length(flag) && flag[n]]
end

"""
    _set_centroid!(c, pos, mass, idx, mass_weighted) -> c

Centroid of the particles `idx`, plain (`mass_weighted = false`) or
mass-weighted, written into `c`.
"""
function _set_centroid!(
    c::Vector{Float64},
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    idx::AbstractVector{Int},
    mass_weighted::Bool,
)::Vector{Float64}
    fill!(c, 0.0)
    w_sum = 0.0
    @inbounds for i in idx
        w = mass_weighted ? mass[i] : 1.0
        c[1] += w * pos[1, i]
        c[2] += w * pos[2, i]
        c[3] += w * pos[3, i]
        w_sum += w
    end
    c ./= w_sum
    return c
end

"""
    shrinking_sphere_centre(pos, mass; weights = :number, shrink = 0.025, n_stop = 1000,
                            stop_fraction = 0.01, n_floor = 500)
        -> (; centre::Vector{Float64}, radius::Float64, n_final::Int, sigma::Float64)

Centre of a particle set by the shrinking-sphere iteration of Power et al.
(2003, DOI 10.1046/j.1365-8711.2003.05925.x). Starting from all particles,
the centre is the centroid of the current set; the set is then reduced to
the `max(n_f, round(Int, (1 − shrink) n))` particles closest to that centre
(at least one particle is dropped per step) and the centroid is recomputed,
until at most `n_f` particles remain.

The stop count follows Power et al.: 1000 particles or 1 % of the set,
whichever is smaller, with a floor of 500 particles or half the set,
whichever is smaller,

    n_f = max(min(n_stop, ⌈stop_fraction N⌉), min(n_floor, ⌈N/2⌉)).

Sets of at most two particles return their centroid directly.

# Arguments
- `pos`: `3 × N` positions
- `mass`: `N` masses
- `weights`: `:number` (plain centroid) or `:mass` (mass-weighted centroid)
- `shrink`: fractional reduction of the set per step, in `(0, 1)`
- `n_stop`, `stop_fraction`, `n_floor`: stop rule above

# Returns
- `centre`: final centre
- `radius`: distance from the centre to the farthest particle of the final set
- `n_final`: particles in the final set
- `sigma`: per-axis centre noise `radius / √(5 n_eff)`, the standard error
  of the centroid of a uniform sphere, with `n_eff = n_final` for `:number`
  and `n_eff = (Σm)² / Σm²` over the final set for `:mass`

Number weighting is the default: with a broad mass spectrum a mass-weighted
centre locks onto the few massive stellar remnants instead of the stellar
density maximum.

# Example
```julia
c = shrinking_sphere_centre(pos, mass).centre
```
"""
function shrinking_sphere_centre(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    weights::Symbol = :number,
    shrink::Real = 0.025,
    n_stop::Integer = 1000,
    stop_fraction::Real = 0.01,
    n_floor::Integer = 500,
)
    weights in (:number, :mass) || throw(
        ArgumentError("shrinking_sphere_centre: weights must be :number or :mass; got :$weights"),
    )
    0 < shrink < 1 ||
        throw(ArgumentError("shrinking_sphere_centre: shrink must lie in (0, 1); got $shrink"))
    n_stop ≥ 1 || throw(ArgumentError("shrinking_sphere_centre: n_stop must be ≥ 1; got $n_stop"))
    0 < stop_fraction ≤ 1 || throw(
        ArgumentError(
            "shrinking_sphere_centre: stop_fraction must lie in (0, 1]; got $stop_fraction",
        ),
    )
    n_floor ≥ 1 ||
        throw(ArgumentError("shrinking_sphere_centre: n_floor must be ≥ 1; got $n_floor"))
    N = length(mass)
    size(pos, 2) == N ||
        throw(DimensionMismatch("shrinking_sphere_centre: $(size(pos, 2)) positions for $N masses"))
    N ≥ 1 || throw(ArgumentError("shrinking_sphere_centre: empty particle set"))
    mass_weighted = weights === :mass
    n_f = max(min(Int(n_stop), ceil(Int, stop_fraction * N)), min(Int(n_floor), cld(N, 2)))
    idx = collect(1:N)
    c = zeros(3)
    while true
        _set_centroid!(c, pos, mass, idx, mass_weighted)
        n = length(idx)
        (N ≤ 2 || n ≤ n_f) && break
        d2 = [(pos[1, i] - c[1])^2 + (pos[2, i] - c[2])^2 + (pos[3, i] - c[3])^2 for i in idx]
        # Guarantee progress where round((1 − shrink) n) == n for small sets
        n_keep = max(n_f, min(n - 1, round(Int, (1 - shrink) * n)))
        idx = idx[partialsortperm(d2, 1:n_keep)]
    end
    d2_max = 0.0
    m_sum = 0.0
    m2_sum = 0.0
    @inbounds for i in idx
        d2_max = max(d2_max, (pos[1, i] - c[1])^2 + (pos[2, i] - c[2])^2 + (pos[3, i] - c[3])^2)
        m_sum += mass[i]
        m2_sum += mass[i]^2
    end
    radius = sqrt(d2_max)
    n_final = length(idx)
    n_eff = mass_weighted ? m_sum^2 / m2_sum : Float64(n_final)
    return (; centre = c, radius = radius, n_final = n_final, sigma = radius / sqrt(5 * n_eff))
end

"""
    _shrinking_sphere_centre(pos, mass; kwargs...) -> Vector{Float64}

Centre of [`shrinking_sphere_centre`](@ref) alone, with the same keywords;
kept for the internal callers that need only the position.
"""
function _shrinking_sphere_centre(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    kwargs...,
)::Vector{Float64}
    return shrinking_sphere_centre(pos, mass; kwargs...).centre
end

"""
    _centre_keywords(settings::PostprocessConfig) -> NamedTuple

Stop rule and step of `shrinking_sphere_centre` from the centre keys of
`settings` (`centre_shrink`, `centre_n_stop`, `centre_stop_fraction`,
`centre_n_floor`), as keywords.
"""
_centre_keywords(settings::PostprocessConfig) = (;
    shrink = settings.centre_shrink,
    n_stop = settings.centre_n_stop,
    stop_fraction = settings.centre_stop_fraction,
    n_floor = settings.centre_n_floor,
)

"""
    _virial_energies(snap, idx) -> (T, W, M)

Kinetic energy about the mass-weighted mean velocity, potential energy and
mass of the systems formed by the snapshot stars `idx` (regularised pairs at
their centres of mass, components of pairs split by `idx` left out).
"""
function _virial_energies(snap::Snapshot, idx::AbstractVector{Int})
    sys, _ = _snapshot_systems(snap, collect(idx))
    isempty(sys.mass) && return (NaN, NaN, 0.0)
    vref = _mean_velocity(sys, eachindex(sys.mass))
    Φ = pair_potential(sys.pos, sys.mass)
    T = 0.0
    W = 0.0
    for s in eachindex(sys.mass)
        T +=
            0.5 *
            sys.mass[s] *
            (
                (sys.vel[1, s] - vref[1])^2 +
                (sys.vel[2, s] - vref[2])^2 +
                (sys.vel[3, s] - vref[3])^2
            )
        W += 0.5 * sys.mass[s] * Φ[s]
    end
    return (T, W, sum(sys.mass))
end

"""
    _lagrangian_radii(pos, mass, centre, fractions) -> Vector{Float64}

Radii about `centre` enclosing the given fractions of the total mass.
"""
function _lagrangian_radii(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    centre::AbstractVector{Float64},
    fractions,
)::Vector{Float64}
    r = [
        sqrt((pos[1, i] - centre[1])^2 + (pos[2, i] - centre[2])^2 + (pos[3, i] - centre[3])^2)
        for i in eachindex(mass)
    ]
    order = sortperm(r)
    cum = cumsum(mass[order])
    M = cum[end]
    return [r[order[something(findfirst(≥(f * M), cum), length(cum))]] for f in fractions]
end

"""
    lagrangian_grid(pos, mass, centre; fractions = LAGR_MASS_FRACTIONS, weights = :mass, n_min = 100)
        -> (; fractions::Vector{Float64}, radii::Vector{Float64}, masked::BitVector,
              sigma_log::Vector{Float64}, n_eff::Float64)

Lagrangian radii about `centre` on a grid of enclosed fractions, with the
entries that the particle number cannot resolve masked.

With `weights = :mass` the radius `r_f` encloses the fraction `f` of the
total mass and the effective particle number is `N_eff = (Σm)² / Σm²`; with
`weights = :number` it encloses the fraction `f` of the particle count and
`N_eff = N`.

# Masking
An entry is masked (`masked[k] = true`, radius `NaN`) when

    f N_eff < n_min  or  (1 − f) N_eff < n_min  or  f ≥ 1,

i.e. when fewer than `n_min` effective particles lie inside or outside the
radius. The 100 % radius is the maximum order statistic of the radial
distribution and is never usable.

# Noise
The enclosed fraction is binomial, `δ ln M = √((1 − f) / (f N_eff))`
(Giersz & Heggie 1994, DOI 10.1093/mnras/268.1.257), and propagates to the
radius through the local logarithmic slope of the mass profile,

    σ(ln r_f) = √((1 − f) / (f N_eff)) / (d ln M / d ln r).

The slope is the central difference `ln(f₊/f₋) / ln(r₊/r₋)` between the
nearest unmasked neighbouring fractions, one-sided at the ends of the
unmasked run; `sigma_log` is `NaN` where an entry is masked or the slope
cannot be formed or vanishes.

# Arguments
- `pos`: `3 × N` positions; `mass`: `N` masses; `centre`: reference point
- `fractions`: enclosed fractions in `(0, 1]`, ascending (default: the 18
  fractions of the engine's `lagr.7`)
- `weights`: `:mass` or `:number`
- `n_min`: minimum effective particle number on either side of a radius

# Example
```julia
g = lagrangian_grid(pos, mass, shrinking_sphere_centre(pos, mass).centre)
r_h = g.radii[findfirst(==(0.5), g.fractions)]
```
"""
function lagrangian_grid(
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    centre::AbstractVector{Float64};
    fractions = LAGR_MASS_FRACTIONS,
    weights::Symbol = :mass,
    n_min::Integer = 100,
)
    weights in (:mass, :number) ||
        throw(ArgumentError("lagrangian_grid: weights must be :mass or :number; got :$weights"))
    n_min ≥ 1 || throw(ArgumentError("lagrangian_grid: n_min must be ≥ 1; got $n_min"))
    f = collect(Float64, fractions)
    all(x -> 0 < x ≤ 1, f) ||
        throw(ArgumentError("lagrangian_grid: fractions must lie in (0, 1]; got $f"))
    issorted(f) || throw(ArgumentError("lagrangian_grid: fractions must be sorted ascending"))
    N = length(mass)
    size(pos, 2) == N ||
        throw(DimensionMismatch("lagrangian_grid: $(size(pos, 2)) positions for $N masses"))
    N ≥ 1 || throw(ArgumentError("lagrangian_grid: empty particle set"))
    mass_weighted = weights === :mass

    d = [
        sqrt((pos[1, i] - centre[1])^2 + (pos[2, i] - centre[2])^2 + (pos[3, i] - centre[3])^2)
        for i in 1:N
    ]
    order = sortperm(d)
    cum = mass_weighted ? cumsum(mass[order]) : collect(1.0:N)
    W = cum[end]
    n_eff = mass_weighted ? sum(mass)^2 / sum(abs2, mass) : Float64(N)

    nf = length(f)
    masked = BitVector(undef, nf)
    radii = fill(NaN, nf)
    for k in 1:nf
        masked[k] = f[k] * n_eff < n_min || (1 - f[k]) * n_eff < n_min || f[k] ≥ 1
        masked[k] && continue
        radii[k] = d[order[something(findfirst(≥(f[k] * W), cum), N)]]
    end

    sigma_log = fill(NaN, nf)
    for k in 1:nf
        masked[k] && continue
        lo = k - 1
        while lo ≥ 1 && masked[lo]
            lo -= 1
        end
        hi = k + 1
        while hi ≤ nf && masked[hi]
            hi += 1
        end
        a, b = if lo ≥ 1 && hi ≤ nf
            (hi, lo)
        elseif hi ≤ nf
            (hi, k)
        elseif lo ≥ 1
            (k, lo)
        else
            continue
        end
        slope = log(f[a] / f[b]) / log(radii[a] / radii[b])
        (isfinite(slope) && slope != 0) || continue
        sigma_log[k] = sqrt((1 - f[k]) / (f[k] * n_eff)) / slope
    end
    return (; fractions = f, radii = radii, masked = masked, sigma_log = sigma_log, n_eff = n_eff)
end

"""
    RadialProfile

Radial structure of one particle set about a centre on two shell sets.

Kinematic shells, of equal membership (counts differing by at most one):
- `r`: member-median radius of each shell
- `r_edges`: shell boundaries (`length(r) + 1`), midpoints between
  consecutive shell members, from 0 to the outer radius
- `rho`: shell density `ΔM / (4π/3 (r_{k+1}³ − r_k³))`
- `sigma_r`, `sigma_t`: radial and one-dimensional tangential velocity
  dispersions of the shell about its mean flow, with velocities taken
  relative to the set's centre-of-mass velocity (`σ_t² = ½ ⟨|v_t − ⟨v_t⟩|²⟩`);
  `NaN` in shells with fewer than `n_min_sigma` members
- `beta`: anisotropy `1 − σ_t²/σ_r²` (0 isotropic, 1 radial); `NaN` in
  shells with fewer than `n_min_beta` members
- `n`: members per shell
- `M`, `r_h`: total mass and half-mass radius of the set
- `sigma_err`: relative error of `σ_r`, `(2(n − 1))^{-1/2}`; `NaN` where `σ_r` is
- `beta_err`: standard error of `β`, `(1 − β)(3/n)^{1/2}`; `NaN` where `β` is

Logarithmic density shells:
- `r_log`: mass-weighted mean radius of each shell
- `r_log_edges`: shell boundaries (`length(r_log) + 1`)
- `rho_log`: shell density, shells merged outward until each holds at least
  `n_min_density` members
- `n_log`: members per shell
"""
struct RadialProfile
    r::Vector{Float64}
    r_edges::Vector{Float64}
    rho::Vector{Float64}
    sigma_r::Vector{Float64}
    sigma_t::Vector{Float64}
    beta::Vector{Float64}
    n::Vector{Int}
    M::Float64
    r_h::Float64
    sigma_err::Vector{Float64}
    beta_err::Vector{Float64}
    r_log::Vector{Float64}
    r_log_edges::Vector{Float64}
    rho_log::Vector{Float64}
    n_log::Vector{Int}
end

"""
    _sorted_median(x, lo, hi) -> Float64

Median of the ascending-sorted entries `x[lo:hi]`.
"""
function _sorted_median(x::AbstractVector{Float64}, lo::Int, hi::Int)::Float64
    n = hi - lo + 1
    m = lo + (n - 1) ÷ 2
    return isodd(n) ? x[m] : 0.5 * (x[m] + x[m + 1])
end

"""
    radial_profile(pos, vel, mass, centre; nbins = 12, n_min_inner = 20,
                   mass_fraction_max = 0.99, n_min_sigma = 50, n_min_beta = 300,
                   n_min_density = 100) -> RadialProfile

Shell profile of the particles about `centre` (N-body or any consistent
units), out to the radius `r_max` enclosing `mass_fraction_max` of the mass.

Kinematic shells: the particles with `d ≤ r_max` are split into
`max(1, min(nbins, N_used ÷ n_min_sigma))` shells of equal membership, so
that the statistical error is the same in every shell. Velocity dispersions
carry the relative error `δσ/σ = (2(n − 1))^{−1/2}` and the anisotropy the
standard error `σ_β = (1 − β)(3/n)^{1/2}`; about 300 members give `β` to
±0.1 near isotropy, hence the default `n_min_beta`. `sigma_r` and `sigma_t`
are `NaN` below `n_min_sigma` members, `beta` below `n_min_beta`.

Logarithmic density shells: `nbins` shells whose outer edges are log-spaced
between the radius enclosing `n_min_inner` particles and `r_max` (the
innermost shell is the sphere inside the first edge). From the inside out,
each shell is merged with its outer neighbour until it holds at least
`n_min_density` members; a short outermost remainder joins the shell inside
it, and a set with fewer than `n_min_density` particles inside `r_max` has a
single shell.

# Arguments
- `pos`, `vel`: `3 × N` positions and velocities; `mass`: `N` masses
- `centre`: reference point, e.g. [`shrinking_sphere_centre`](@ref)
- `nbins ≥ 2`: maximum number of kinematic shells, number of logarithmic shells
  before merging
- `n_min_inner`: particles inside the first logarithmic edge (`N ≥ 2 n_min_inner`)
- `mass_fraction_max`: enclosed mass fraction defining `r_max`
- `n_min_sigma ≥ 2`, `n_min_beta ≥ n_min_sigma`, `n_min_density ≥ 1`:
  membership thresholds above

# Example
```julia
c = shrinking_sphere_centre(pos, mass).centre
p = radial_profile(pos, vel, mass, c)
lines(p.r, p.beta)
```
"""
function radial_profile(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    centre::AbstractVector{Float64};
    nbins::Int = 12,
    n_min_inner::Int = 20,
    mass_fraction_max::Float64 = 0.99,
    n_min_sigma::Int = 50,
    n_min_beta::Int = 300,
    n_min_density::Int = 100,
)::RadialProfile
    N = length(mass)
    N ≥ 2 * n_min_inner ||
        throw(ArgumentError("radial_profile: need at least $(2 * n_min_inner) particles; got $N"))
    nbins ≥ 2 || throw(ArgumentError("radial_profile: nbins must be ≥ 2"))
    n_min_sigma ≥ 2 ||
        throw(ArgumentError("radial_profile: n_min_sigma must be ≥ 2; got $n_min_sigma"))
    n_min_beta ≥ n_min_sigma || throw(
        ArgumentError(
            "radial_profile: n_min_beta must be ≥ n_min_sigma = $n_min_sigma; got $n_min_beta",
        ),
    )
    n_min_density ≥ 1 ||
        throw(ArgumentError("radial_profile: n_min_density must be ≥ 1; got $n_min_density"))
    d = [
        sqrt((pos[1, i] - centre[1])^2 + (pos[2, i] - centre[2])^2 + (pos[3, i] - centre[3])^2)
        for i in 1:N
    ]
    order = sortperm(d)
    ds = d[order]
    cum = cumsum(mass[order])
    M = cum[end]
    r_h = ds[findfirst(≥(0.5 * M), cum)]
    r_min = max(ds[n_min_inner], eps())
    r_max = ds[findfirst(≥(mass_fraction_max * M), cum)]
    r_max > r_min || (r_max = r_min * 10)
    N_used = searchsortedlast(ds, r_max)::Int
    vc = vec(sum(vel .* mass'; dims = 2)) ./ M

    # Kinematic shells: sorted positions bounds[s]+1 : bounds[s+1]
    n_k = max(1, min(nbins, N_used ÷ n_min_sigma))
    base, extra = divrem(N_used, n_k)
    bounds = zeros(Int, n_k + 1)
    for s in 1:n_k
        bounds[s + 1] = bounds[s] + base + (s ≤ extra ? 1 : 0)
    end
    edges = Vector{Float64}(undef, n_k + 1)
    edges[1] = 0.0
    edges[end] = r_max
    for s in 1:(n_k - 1)
        j = bounds[s + 1]
        edges[s + 1] = 0.5 * (ds[j] + ds[j + 1])
    end
    r_med = fill(NaN, n_k)
    rho = fill(NaN, n_k)
    σr = fill(NaN, n_k)
    σt = fill(NaN, n_k)
    β = fill(NaN, n_k)
    σ_err = fill(NaN, n_k)
    β_err = fill(NaN, n_k)
    n = zeros(Int, n_k)
    n_buf = base + (extra > 0 ? 1 : 0)
    vr = Vector{Float64}(undef, n_buf)
    vt = Matrix{Float64}(undef, 3, n_buf)
    for s in 1:n_k
        lo = bounds[s] + 1
        hi = bounds[s + 1]
        ns = hi - lo + 1
        n[s] = ns
        ΔM = 0.0
        for j in lo:hi
            ΔM += mass[order[j]]
        end
        rho[s] = ΔM / (4π / 3 * (edges[s + 1]^3 - edges[s]^3))
        r_med[s] = _sorted_median(ds, lo, hi)
        ns < n_min_sigma && continue
        mean_vr = 0.0
        mvx, mvy, mvz = 0.0, 0.0, 0.0
        @inbounds for jj in 1:ns
            i = order[lo + jj - 1]
            rn = max(d[i], eps())
            ux = (pos[1, i] - centre[1]) / rn
            uy = (pos[2, i] - centre[2]) / rn
            uz = (pos[3, i] - centre[3]) / rn
            dvx, dvy, dvz = vel[1, i] - vc[1], vel[2, i] - vc[2], vel[3, i] - vc[3]
            v = dvx * ux + dvy * uy + dvz * uz
            vr[jj] = v
            vt[1, jj] = dvx - v * ux
            vt[2, jj] = dvy - v * uy
            vt[3, jj] = dvz - v * uz
            w = mass[i] / ΔM
            mean_vr += w * v
            mvx += w * vt[1, jj]
            mvy += w * vt[2, jj]
            mvz += w * vt[3, jj]
        end
        σr2 = 0.0
        σt2 = 0.0
        @inbounds for jj in 1:ns
            w = mass[order[lo + jj - 1]] / ΔM
            σr2 += w * (vr[jj] - mean_vr)^2
            σt2 += w * ((vt[1, jj] - mvx)^2 + (vt[2, jj] - mvy)^2 + (vt[3, jj] - mvz)^2)
        end
        σt2 *= 0.5
        σr[s] = sqrt(σr2)
        σt[s] = sqrt(σt2)
        σ_err[s] = 1 / sqrt(2 * (ns - 1))
        (ns ≥ n_min_beta && σr2 > 0) || continue
        β[s] = 1 - σt2 / σr2
        β_err[s] = (1 - β[s]) * sqrt(3 / ns)
    end

    # Logarithmic density shells: raw shell k holds sorted positions
    # raw[k] : raw[k+1]-1, the outermost one closed at r_max
    log_edges = vcat(0.0, exp10.(range(log10(r_min), log10(r_max); length = nbins)))
    log_edges[end] = r_max
    raw = Vector{Int}(undef, nbins + 1)
    raw[1] = 1
    for k in 2:nbins
        raw[k] = searchsortedfirst(ds, log_edges[k])
    end
    raw[nbins + 1] = N_used + 1
    first_pos = Int[]
    last_pos = Int[]
    r_log_edges = [0.0]
    start = 1
    for k in 1:nbins
        if raw[k + 1] - start ≥ n_min_density || k == nbins
            push!(first_pos, start)
            push!(last_pos, raw[k + 1] - 1)
            push!(r_log_edges, log_edges[k + 1])
            start = raw[k + 1]
        end
    end
    # A short outermost remainder joins the shell inside it (the pop first:
    # `end` is evaluated after the array has shrunk).
    if length(first_pos) ≥ 2 && last_pos[end] - first_pos[end] + 1 < n_min_density
        outermost = pop!(last_pos)
        last_pos[end] = outermost
        pop!(first_pos)
        deleteat!(r_log_edges, length(r_log_edges) - 1)
    end
    n_l = length(first_pos)
    r_log = fill(NaN, n_l)
    rho_log = zeros(n_l)
    n_log = zeros(Int, n_l)
    for s in 1:n_l
        n_log[s] = last_pos[s] - first_pos[s] + 1
        ΔM = 0.0
        mr = 0.0
        for j in first_pos[s]:last_pos[s]
            m = mass[order[j]]
            ΔM += m
            mr += m * ds[j]
        end
        rho_log[s] = ΔM / (4π / 3 * (r_log_edges[s + 1]^3 - r_log_edges[s]^3))
        n_log[s] > 0 && (r_log[s] = mr / ΔM)
    end
    return RadialProfile(
        r_med,
        edges,
        rho,
        σr,
        σt,
        β,
        n,
        M,
        r_h,
        σ_err,
        β_err,
        r_log,
        r_log_edges,
        rho_log,
        n_log,
    )
end

"""
    cluster_profiles(snap, cluster_ranges; settings = PostprocessConfig(), bound_only = true,
                     nbins = settings.profile_nbins, n_min_sigma = settings.profile_n_min_sigma,
                     n_min_beta = settings.profile_n_min_beta,
                     n_min_density = settings.profile_n_min_density,
                     max_passes = 100, min_members = 10)
        -> Vector{Union{Nothing,RadialProfile}}

Radial profile of every initial cluster about its own shrinking-sphere
centre from its bound members ([`cluster_structure`](@ref) conventions;
`max_passes` and `min_members` go to [`bound_set`](@ref), the shell
thresholds to [`radial_profile`](@ref)); `nothing` for clusters with too
few members. The centre takes its stop rule and step from the centre keys of
`settings`, a [`PostprocessConfig`](@ref), whose profile keys are also the
defaults of the shell thresholds; an explicit keyword overrides them.
"""
function cluster_profiles(
    snap::Snapshot,
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    settings::PostprocessConfig = PostprocessConfig(),
    bound_only::Bool = true,
    nbins::Int = settings.profile_nbins,
    n_min_sigma::Int = settings.profile_n_min_sigma,
    n_min_beta::Int = settings.profile_n_min_beta,
    n_min_density::Int = settings.profile_n_min_density,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)
    pos_all = Float64.(snap.pos)
    vel_all = Float64.(snap.vel)
    mass_all = Float64.(snap.mass)
    out = Vector{Union{Nothing,RadialProfile}}(nothing, length(cluster_ranges))
    for (i, rng) in enumerate(cluster_ranges)
        idx = _member_indices(snap, rng)
        length(idx) < 2 * _MIN_MEMBERS && continue
        sel = bound_only ? bound_set(snap, idx; max_passes, min_members).stars : idx
        length(sel) < 2 * _MIN_MEMBERS && continue
        p = pos_all[:, sel]
        c = _shrinking_sphere_centre(p, mass_all[sel]; _centre_keywords(settings)...)
        out[i] = radial_profile(
            p,
            vel_all[:, sel],
            mass_all[sel],
            c;
            nbins,
            n_min_sigma,
            n_min_beta,
            n_min_density,
        )
    end
    return out
end

"""
    system_profile(snap; settings = PostprocessConfig(), nbins = settings.profile_nbins,
                   n_min_sigma = settings.profile_n_min_sigma,
                   n_min_beta = settings.profile_n_min_beta,
                   n_min_density = settings.profile_n_min_density) -> RadialProfile

Radial profile of the whole snapshot about its shrinking-sphere centre
(the remnant after coalescence); the shell keywords go to
[`radial_profile`](@ref). The centre takes its stop rule and step from the
centre keys of `settings`, a [`PostprocessConfig`](@ref), whose profile keys
are also the defaults of the shell keywords; an explicit keyword overrides
them.
"""
function system_profile(
    snap::Snapshot;
    settings::PostprocessConfig = PostprocessConfig(),
    nbins::Int = settings.profile_nbins,
    n_min_sigma::Int = settings.profile_n_min_sigma,
    n_min_beta::Int = settings.profile_n_min_beta,
    n_min_density::Int = settings.profile_n_min_density,
)::RadialProfile
    pos = Float64.(snap.pos)
    vel = Float64.(snap.vel)
    mass = Float64.(snap.mass)
    c = _shrinking_sphere_centre(pos, mass; _centre_keywords(settings)...)
    return radial_profile(pos, vel, mass, c; nbins, n_min_sigma, n_min_beta, n_min_density)
end

"""
    local_relaxation_time(rho, sigma_1d, m_mean, N, gamma) -> Float64

Local two-body relaxation time of a shell (`G = 1`),

    t_r = 0.338 σ_1D³ / (m̄ ρ ln(γ N)),

from the mass density `rho`, the one-dimensional velocity dispersion
`sigma_1d` and the mean stellar mass `m_mean` (Giersz & Heggie 1994,
DOI 10.1093/mnras/268.1.257). `NaN` for non-positive or non-finite inputs
and for `γ N ≤ 1`, where the Coulomb logarithm is not positive.

# Example
```julia
p = radial_profile(pos, vel, mass, c)
m̄ = sum(mass) / length(mass)
t_r = [local_relaxation_time(p.rho[k], p.sigma_r[k], m̄, length(mass), 0.02)
       for k in eachindex(p.rho)]
```
"""
function local_relaxation_time(
    rho::Real,
    sigma_1d::Real,
    m_mean::Real,
    N::Integer,
    gamma::Real,
)::Float64
    (rho > 0 && sigma_1d > 0 && m_mean > 0 && N > 0 && gamma > 0) || return NaN
    (isfinite(rho) && isfinite(sigma_1d) && isfinite(m_mean) && isfinite(gamma)) || return NaN
    γN = gamma * N
    γN > 1 || return NaN
    return 0.338 * sigma_1d^3 / (m_mean * rho * log(γN))
end

"""
    relaxation_times(pos, vel, mass, centre; gamma = (0.02, 0.11), rho_core = NaN, sigma_core = NaN)
        -> (; t_rh::NTuple{2,Float64}, t_rc::NTuple{2,Float64}, t_cr::Float64, r_h::Float64,
             sigma_3d::Float64, m_mean::Float64, n::Int)

Half-mass and core relaxation times and the crossing time of a particle set
in N-body units (`G = 1`).

- Half-mass relaxation time (Spitzer 1987, *Dynamical Evolution of Globular
  Clusters*, Princeton University Press; the form quoted by Baumgardt 2001):

      t_rh = 0.138 √N r_h^{3/2} / (√m̄ ln(γ N)),

  with `r_h` the half-mass radius about `centre` and `m̄ = M/N`.
- Core relaxation time (Giersz & Heggie 1994, DOI 10.1093/mnras/268.1.257),
  from the core density `rho_core` and the one-dimensional core dispersion
  `sigma_core`:

      t_rc = 0.338 σ_c³ / (m̄ ρ_c ln(γ N));

  `NaN` when either is `NaN` or non-positive.
- Crossing time `t_cr = 2 r_h / σ_3D`, with `σ_3D² = Σ m |v − v̄|² / M` about
  the mass-weighted mean velocity; the energy-based definition
  `G M^{5/2} / (2|E|)^{3/2}` differs from it by a factor of order unity.

Both times are reported for the two Coulomb-logarithm factors in `gamma`:
`γ = 0.11` is calibrated for equal masses (Giersz & Heggie 1994), `γ ≈ 0.02`
is adopted for a mass spectrum (e.g. Vesperini et al. 2013,
DOI 10.1093/mnras/sts434). Neither is calibrated for a Kroupa spectrum with
stellar remnants, so the pair brackets the uncertainty.

# Arguments
- `pos`, `vel`: `3 × N` positions and velocities; `mass`: `N` masses
- `centre`: reference point for `r_h`
- `gamma`: two Coulomb-logarithm factors, each `> 0` with `γ N > 1`
- `rho_core`, `sigma_core`: core density and one-dimensional core dispersion

# Example
```julia
c = shrinking_sphere_centre(pos, mass).centre
t = relaxation_times(pos, vel, mass, c)
t.t_rh[1] / t.t_cr
```
"""
function relaxation_times(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    centre::AbstractVector{Float64};
    gamma::NTuple{2,Real} = (0.02, 0.11),
    rho_core::Real = NaN,
    sigma_core::Real = NaN,
)
    N = length(mass)
    N ≥ 1 || throw(ArgumentError("relaxation_times: empty particle set"))
    (size(pos, 2) == N && size(vel, 2) == N) || throw(
        DimensionMismatch(
            "relaxation_times: $(size(pos, 2)) positions and $(size(vel, 2)) velocities for $N masses",
        ),
    )
    for γ in gamma
        (γ > 0 && γ * N > 1) || throw(
            ArgumentError(
                "relaxation_times: gamma entries must be > 0 with gamma * N > 1 (N = $N); got $γ",
            ),
        )
    end
    M = sum(mass)
    m_mean = M / N
    r_h = _lagrangian_radii(pos, mass, centre, (0.5,))[1]
    vx, vy, vz = 0.0, 0.0, 0.0
    @inbounds for i in 1:N
        vx += mass[i] * vel[1, i]
        vy += mass[i] * vel[2, i]
        vz += mass[i] * vel[3, i]
    end
    vx /= M
    vy /= M
    vz /= M
    s2 = 0.0
    @inbounds for i in 1:N
        s2 += mass[i] * ((vel[1, i] - vx)^2 + (vel[2, i] - vy)^2 + (vel[3, i] - vz)^2)
    end
    sigma_3d = sqrt(s2 / M)
    t_rh = ntuple(k -> 0.138 * sqrt(N) * r_h^1.5 / (sqrt(m_mean) * log(gamma[k] * N)), 2)
    t_rc = ntuple(k -> local_relaxation_time(rho_core, sigma_core, m_mean, N, gamma[k]), 2)
    t_cr = 2 * r_h / sigma_3d
    return (;
        t_rh = t_rh,
        t_rc = t_rc,
        t_cr = t_cr,
        r_h = r_h,
        sigma_3d = sigma_3d,
        m_mean = m_mean,
        n = N,
    )
end

"""
    bound_fraction(snap::Snapshot; max_passes = 100, min_members = 10) -> Float64

Mass fraction of the stars bound to the whole system in its own frame
(iterative energetic unbinding with regularised pairs at their centres of
mass, [`bound_set`](@ref)): a snapshot-based escaper measure independent of
the engine's escape sphere. One O(N²) pair sum plus O(N·N_removed) per pass.
"""
function bound_fraction(
    snap::Snapshot;
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)::Float64
    mass = Float64.(snap.mass)
    isempty(mass) && return NaN
    sel = bound_set(snap; max_passes, min_members).stars
    return sum(mass[sel]) / sum(mass)
end

"""
    cluster_structure(snaps, cluster_ranges; settings = PostprocessConfig(), bound_only = true,
                      max_passes = 100, min_members = 10) -> ClusterStructure

Structure of every initial cluster at every snapshot, measured about the
cluster's own shrinking-sphere centre from the members it binds
(`bound_only = false` uses all present members and reproduces the plain
membership statistics). The bound selection removes tidally stripped
stars and kicked stellar remnants, which otherwise dominate the velocity
dispersion and the virial ratio after the first supernovae.

Regularised pairs enter the bound selection, `sigma_1d` and `q_virial` as
their centres of mass; `max_passes` and `min_members` go to
[`bound_set`](@ref). The centre takes its stop rule and step from the centre
keys of `settings`, a [`PostprocessConfig`](@ref). One O(N_i²) pair sum per
cluster and snapshot.
"""
function cluster_structure(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    settings::PostprocessConfig = PostprocessConfig(),
    bound_only::Bool = true,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)::ClusterStructure
    n_cl = length(cluster_ranges)
    n_t = length(snaps)
    time = [time_nb(s.header) for s in snaps]
    n_members = zeros(Int, n_cl, n_t)
    n_bound = zeros(Int, n_cl, n_t)
    fbound = fill(NaN, n_cl, n_t)
    centre = fill(NaN, 3, n_cl, n_t)
    r_lagr = fill(NaN, 3, n_cl, n_t)
    σ = fill(NaN, n_cl, n_t)
    q = fill(NaN, n_cl, n_t)
    for (k, snap) in enumerate(snaps)
        pos_all = Float64.(snap.pos)
        mass_all = Float64.(snap.mass)
        for (i, rng) in enumerate(cluster_ranges)
            idx = _member_indices(snap, rng)
            n_members[i, k] = length(idx)
            length(idx) < _MIN_MEMBERS && continue
            if bound_only
                b = bound_set(snap, idx; max_passes, min_members)
                sel = b.stars
                T, W = b.kinetic, b.potential
            else
                sel = idx
                T, W, _ = _virial_energies(snap, idx)
            end
            n_bound[i, k] = length(sel)
            fbound[i, k] = sum(mass_all[sel]) / sum(mass_all[idx])
            length(sel) < _MIN_MEMBERS && continue
            p = pos_all[:, sel]
            m = mass_all[sel]
            M = sum(m)
            c = _shrinking_sphere_centre(p, m; _centre_keywords(settings)...)
            centre[:, i, k] = c
            r_lagr[:, i, k] = _lagrangian_radii(p, m, c, _LAGR_FRACTIONS_MEMBER)
            σ[i, k] = sqrt(2 * T / (3 * M))
            q[i, k] = W < 0 ? T / abs(W) : NaN
        end
    end
    return ClusterStructure(time, n_members, n_bound, fbound, centre, r_lagr, σ, q)
end

# ---------------------------------------------------------------------------
# Original-membership diagnostics of a merger run
# ---------------------------------------------------------------------------

"""
    parse_merger_summary(path::AbstractString) -> Vector{Vector{Int}}

Read `merger_summary.txt` and return the body indices (`dat.10` order,
1-based) of every initial cluster from the post-truncation body counts
and, when present, the pair counts: pairs of all clusters come first,
cluster by cluster, then the singles, so a cluster with binaries owns two
contiguous blocks. Returns an empty vector if the file cannot be parsed.
"""
function parse_merger_summary(path::AbstractString)::Vector{Vector{Int}}
    isfile(path) || return Vector{Int}[]

    counts = Tuple{Int,Int}[]   # (bodies after truncation, pairs)
    # Match e.g. "Cluster 5: plummer, imf=kroupa, N=1000 (after trunc: 966, binaries: 12),"
    # Tolerant of extra comma-separated fields between the profile name and
    # the N=… count (the summary format has grown fields before; the
    # write→parse round-trip test in runtests.jl guards this coupling).
    pattern = r"Cluster\s+\d+:\s+.*?\bN=\d+\s+\(after trunc:\s+(\d+)(?:,\s*binaries:\s+(\d+))?\)"
    for line in eachline(path)
        m = match(pattern, line)
        m === nothing && continue
        n_b_capture = m.captures[2]
        n_b = n_b_capture === nothing ? 0 : parse(Int, n_b_capture)
        push!(counts, (parse(Int, something(m.captures[1])), n_b))
    end
    isempty(counts) && return Vector{Int}[]

    nbin0 = sum(last, counts)
    members = Vector{Vector{Int}}(undef, length(counts))
    pair_offset = 0
    single_offset = 2nbin0
    for (i, (n, n_b)) in enumerate(counts)
        n_s = n - 2n_b
        blocks = [(pair_offset + 1):(pair_offset + 2n_b), (single_offset + 1):(single_offset + n_s)]
        members[i] = _members_from_blocks(blocks)
        pair_offset += 2n_b
        single_offset += n_s
    end
    return members
end

"""
    _cluster_virial_snapshot(snap, rng; bound_only = true, max_passes, min_members) -> (Q, n_mem)

Virial ratio `Q = T/|W|` of the members of an initial cluster in one
snapshot, about the members' mass-weighted mean velocity and with only the
self-gravity among the selected members in `W`; regularised pairs count as
their centres of mass. With `bound_only` the selection is the bound subset
([`bound_set`](@ref)), which removes tidally stripped stars and kicked stellar
remnants; otherwise every present member counts. `n_mem` is the number of
present members. Returns `(NaN, n_mem)` with fewer than 3 selected members or
a non-negative `W`. N-body units (`G = 1`).
"""
function _cluster_virial_snapshot(
    snap::Snapshot,
    rng::AbstractVector{Int};
    bound_only::Bool = true,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)
    idx = _member_indices(snap, rng)
    n_mem = length(idx)
    n_mem < 3 && return (NaN, n_mem)
    if bound_only
        b = bound_set(snap, idx; max_passes, min_members)
        length(b.stars) < 3 && return (NaN, n_mem)
        T, W = b.kinetic, b.potential
    else
        T, W, _ = _virial_energies(snap, idx)
    end
    W < 0 || return (NaN, n_mem)
    return (T / abs(W), n_mem)
end

"""
    per_cluster_virial(snaps, cluster_ranges; bound_only = true)
        -> (Q::Matrix{Float64}, n_mem::Matrix{Int})

Internal virial ratio of every initial cluster at every snapshot, `Q[i, k]`
for cluster `i` at snapshot `k` (`NaN` with fewer than 3 selected members).
By default only the members bound to the cluster enter `T` and `W`
([`bound_set`](@ref), which takes `max_passes` and `min_members`);
`bound_only = false` uses every present member and is dominated by kicked
remnants after the first supernovae. One O(N_i²) pair sum per cluster and
snapshot.
"""
function per_cluster_virial(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    bound_only::Bool = true,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)
    n_cl = length(cluster_ranges)
    n_t = length(snaps)
    Q = fill(NaN, n_cl, n_t)
    n_mem = zeros(Int, n_cl, n_t)

    for (k, snap) in enumerate(snaps)
        for (i, rng) in enumerate(cluster_ranges)
            q, nm = _cluster_virial_snapshot(snap, rng; bound_only, max_passes, min_members)
            Q[i, k] = q
            n_mem[i, k] = nm
        end
    end
    return Q, n_mem
end
