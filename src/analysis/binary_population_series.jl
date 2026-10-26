# =============================================================================
# Binary population series: number, star and mass binary fractions
# per radial shell, stellar class and census level of a snapshot sequence
# =============================================================================

"""
    BinaryPopulationRow

Binary fractions of one snapshot, census level, radial bin and stellar class
from [`binary_population_series`](@ref). Radii are in the N-body units of the
snapshot, masses in M☉.

| field | meaning | unit |
|---|---|---|
| `snapshot` | index of the snapshot in the input vector | — |
| `t_nb`, `t_myr` | snapshot time | NB, Myr |
| `level` | census level: `:ks` (E, KS pairs) or `:bound` (P, bound pairs) | — |
| `radius_frame` | `:lagrangian` (enclosed-mass bins), `:rh` (bins in r/r_h), `:projected_rh` (bins in R/R_h) or `:milone` (core, half-mass and outer regions in projection) | — |
| `bin` | 0 for all systems, else the shell index 1, 2, … | — |
| `edge_lo`, `edge_hi` | shell edges: mass fractions (`:lagrangian`), multiples of r_h (`:rh`) or of R_h (`:projected_rh`, `:milone`) | — |
| `r_lo_nb`, `r_hi_nb` | shell radii; `r_hi_nb` is `NaN` for an empty Lagrangian shell, whose `r_lo_nb` is the outer radius of the last non-empty shell inside it | NB |
| `class` | `:all` or the `key` of a [`STELLAR_CLASSES`](@ref) entry | — |
| `n_singles` | single stars (of the class) | — |
| `n_binaries` | binaries (with a primary of the class) | — |
| `n_in_binaries` | binary components (of the class) | — |
| `n_hierarchies` | binaries of `n_binaries` in a hierarchy: at `:bound` the inner pair of a stable triple, at `:ks` the inner binary of a merged hierarchy (KS slot with a negative `pair_name`) | — |
| `m_singles_msun` | mass of the singles of `n_singles` | M☉ |
| `m_binaries_msun` | mass of the binaries (`:all`) or of the components of the class in binaries | M☉ |
| `f_b_number` | ``f_b = N_b / (N_s + N_b)`` | — |
| `f_b_stars` | ``f_⋆ = N_\\mathrm{in} / (N_s + N_\\mathrm{in})`` | — |
| `f_b_mass` | ``f_M = M_b / (M_s + M_b)`` | — |
| `f_b_number_se` | binomial standard error ``\\sqrt{f_b (1 - f_b) / (N_s + N_b)}`` | — |
| `f_b_q` | fraction of binaries with ``q > q_\\mathrm{min}`` (and, with stellar types, two main-sequence components) among the systems with single or primary mass in the window | — |
| `r_h_binaries_nb`, `r_h_singles_nb` | half-mass radii of the binaries and singles of the row; `NaN` on shells | NB |

Fractions with a vanishing denominator are `NaN`.
"""
struct BinaryPopulationRow
    snapshot::Int
    t_nb::Float64
    t_myr::Float64
    level::Symbol
    radius_frame::Symbol
    bin::Int
    edge_lo::Float64
    edge_hi::Float64
    r_lo_nb::Float64
    r_hi_nb::Float64
    class::Symbol
    n_singles::Int
    n_binaries::Int
    n_in_binaries::Int
    n_hierarchies::Int
    m_singles_msun::Float64
    m_binaries_msun::Float64
    f_b_number::Float64
    f_b_stars::Float64
    f_b_mass::Float64
    f_b_number_se::Float64
    f_b_q::Float64
    r_h_binaries_nb::Float64
    r_h_singles_nb::Float64
end

const _BPOP_LEVELS = (:ks, :bound)
const _BPOP_FRAMES = (:lagrangian, :rh, :projected_rh, :milone)
const _BPOP_EDGES_LAGRANGIAN = (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)
const _BPOP_EDGES_RH = (0.0, 0.5, 1.0, 2.0, 4.0, Inf)
# Milone regions in units of R_h; the inner edge R_c/R_h is set per snapshot
const _BPOP_EDGES_MILONE = (0.0, NaN, 1.0, 2.0)
# Systems per annulus of the projected number-density profile: 1 % of the
# systems, at least this many
const _BPOP_CORE_MIN_SYSTEMS = 20
# Relative tolerance of an epoch match in N-body time: the sev.83/bev.82 data
# lines carry six significant digits, conf.3 a Float32.
const _BPOP_EPOCH_RTOL = 1.0e-5

# Systems of a snapshot at one census level: the pairs first, then the singles.
# `secondary` is 0 for a single, whose `primary` is the star itself. `radius` is
# the radius of the frame (projected for :projected_rh and :milone), `r_3d` the
# three-dimensional distance from the centre; `a_nb`, `e`, `eb_nb` are the
# relative-orbit elements and binding energy m₁m₂/(2a) (NB, G = 1), NaN for
# singles and unbound pairs; `velocity` is 3 × n, centre-of-mass velocities.
struct _BpopSystems
    mass::Vector{Float64}
    radius::Vector{Float64}
    primary::Vector{Int}
    secondary::Vector{Int}
    hierarchy::Vector{Bool}
    a_nb::Vector{Float64}
    e::Vector{Float64}
    eb_nb::Vector{Float64}
    r_3d::Vector{Float64}
    velocity::Matrix{Float64}
end

# Counts and masses per (bin + 1, class column); column 1 is `:all`.
struct _BpopTally
    n_singles::Matrix{Int}
    n_binaries::Matrix{Int}
    n_in_binaries::Matrix{Int}
    n_hierarchies::Matrix{Int}
    m_singles::Matrix{Float64}
    m_binaries::Matrix{Float64}
    n_window::Matrix{Int}
    n_window_q::Matrix{Int}
end

_BpopTally(nr::Int, nc::Int) = _BpopTally(
    zeros(Int, nr, nc),
    zeros(Int, nr, nc),
    zeros(Int, nr, nc),
    zeros(Int, nr, nc),
    zeros(nr, nc),
    zeros(nr, nc),
    zeros(Int, nr, nc),
    zeros(Int, nr, nc),
)

"""
    binary_population_series(snapshots::AbstractVector{Snapshot}; levels = (:ks, :bound),
                             radius_frame::Symbol = :lagrangian, edges = nothing,
                             centre = :shrinking_sphere, centre_kwargs::NamedTuple = (;),
                             kstar = nothing, q_min::Real = 0.5,
                             primary_window_msun::Tuple{<:Real,<:Real} = (0.47, 0.76),
                             k::Integer = 10, n_perturbers::Integer = 4,
                             gamma_max::Real = 0.1, max_untyped_fraction::Real = 0.01,
                             projection_axis::Integer = 3, r_tide_nb::Real = NaN)
        -> Vector{BinaryPopulationRow}

Binary fractions of a snapshot sequence per census level, radial bin and
stellar class, one [`BinaryPopulationRow`](@ref) per snapshot × level × bin ×
class, in that order.

# Arguments
- `snapshots`: the snapshots
- `levels`: census levels, any of `:ks` (level E) and `:bound` (level P)
- `radius_frame`: `:lagrangian` (shells in enclosed mass fraction), `:rh`
  (shells in units of the half-mass radius), `:projected_rh` (shells in units
  of the projected half-mass radius) or `:milone` (the three projected regions
  of Milone et al. 2012)
- `edges`: shell edges, strictly increasing from 0; mass fractions ending at 1
  for `:lagrangian` (default `[0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]`), multiples
  of r_h for `:rh` and of R_h for `:projected_rh` (default
  `[0, 0.5, 1, 2, 4, Inf]`); `:milone` takes no edges
- `centre`: `:shrinking_sphere` ([`shrinking_sphere_centre`](@ref) of the
  stars of each snapshot, with `centre_kwargs`), `:origin`, or a fixed point
  (length-3 vector, NB)
- `kstar`: `nothing` (class `:all` only) or a vector aligned with `snapshots`
  of name → K* maps, e.g. from [`stellar_type_map`](@ref). A star without a
  type counts under `:all` only: the engine lists neither the outer body of a
  merged hierarchy (a KS pair with a negative `pair_name`) nor the body that
  represents a chain subsystem in sev.83 or bev.82
- `q_min`, `primary_window_msun`: mass-ratio threshold and inclusive primary
  mass window [M☉] of `f_b_q`; the default window is that of Milone et al.
  (2012, DOI 10.1051/0004-6361/201016384, §5), 0.47 < M < 0.76 M☉
- `k`, `n_perturbers`, `gamma_max`: pair membership at level P, see
  [`bound_pairs`](@ref)
- `max_untyped_fraction`: largest fraction of the massive stars of a snapshot
  without a stellar type; the untyped stars up to it count under `:all` only,
  more indicates mismatched epochs and raises an `ArgumentError`
- `projection_axis`: line of sight of the projected frames, 1, 2 or 3 (x, y,
  z); the projected radius is the distance from the centre in the plane normal
  to it
- `r_tide_nb`: tidal radius [NB]; when finite, the systems farther than it
  from the centre (three-dimensional distance) are removed before any count or
  radius, and the total number removed is logged once at `@info` level; `NaN`
  (default) or `Inf` keeps every system

Massless bodies are no systems: an unpaired massless star (the ghost of a
merged or chain subsystem) is left out of every count, and a pair with a
massless component (in particular one of summed mass 0) is not formed, its
massive component, if any, counting as a single.

# Method
At a level the stars of a snapshot form systems. Level `:ks` takes the KS
pairs (stars 2k − 1 and 2k) whose centre-of-mass name is positive. A KS slot
with a negative name is a merged hierarchy: star 2k − 1 is the first component
of the inner binary, star 2k the outer body, and the second inner component is
written among the single stars; the binary is then (2k − 1, j), j the star
nearest to 2k − 1 (Euclidean distance, ties to the lower index) among the
massive stars that are not KS components, counted in `n_hierarchies`, and the
outer body 2k is a single. Level `:bound` takes the pairs of
[`bound_pairs`](@ref); the inner pair of a stable triple is a binary and its
outer body a single. Every other star is a single. A binary has the mass m₁ + m₂ and the position
of its centre of mass; its primary is the heavier component and
q = m₂/m₁ ≤ 1. Masses are converted to M☉ with [`zmbar`](@ref).

Radii are measured from the centre of each snapshot, computed once for all
levels. The half-mass radius r_h of a set of systems is the radius of the
first system, in radius order, at which the enclosed mass reaches half the
total. With `:rh` shell b holds the systems with
``e_b r_h ≤ r < e_{b+1} r_h``; with `:lagrangian` those whose enclosed mass
fraction lies in ``(e_b, e_{b+1}]``. Bin 0 holds every system.

The projected frames take the projected radius R, the distance from the
centre in the plane normal to `projection_axis`, and the projected half-mass
radius R_h, defined from R as r_h is from r. With `:projected_rh` shell b holds
the systems with ``e_b R_h ≤ R < e_{b+1} R_h``. `:milone` has the three regions
of Milone et al. (2012, DOI 10.1051/0004-6361/201016384, §3.2): bin 1
``R < R_c``, bin 2 ``R_c ≤ R < R_h`` and bin 3 ``R_h ≤ R < 2R_h``, with
`edge_lo`, `edge_hi` in units of R_h. R_c is the projected radius at which the
projected number density of the systems falls to half its central value. With
``n_a = \\max(20, ⌈0.01 N⌉)`` systems, the central density is that inside the
projected radius of the ``n_a``-th system in R order, and the profile is
measured in consecutive annuli of ``n_a`` systems each, placed at their
area-weighted mean radius ``\\sqrt{(R_\\mathrm{in}^2 + R_\\mathrm{out}^2)/2}``
(the central disc at ``R_{n_a}/\\sqrt{2}``); R_c is interpolated linearly in R
between the last point above half the central density and the first at or
below it. It is `NaN`, and bins 1 and 2 are empty, when the profile does not
fall to half, or with fewer than ``2n_a`` systems.

With a finite `r_tide_nb` the systems beyond the tidal radius are removed
first; every count, r_h, R_h and R_c refer to the systems inside it.

With ``N_s`` singles and ``N_b`` binaries of total masses ``M_s`` and ``M_b``,
the number fraction is ``f_b = N_b/(N_s + N_b)``, with binomial standard error
``\\sqrt{f_b(1 - f_b)/(N_s + N_b)}``, the fraction of stars in binaries
``f_⋆ = 2N_b/(N_s + 2N_b) = 2f_b/(1 + f_b)`` and the mass fraction
``f_M = M_b/(M_s + M_b)``. A class row counts the singles of the class and the
binaries whose primary belongs to it; its `n_in_binaries` and
`m_binaries_msun` count the components of the class in any binary, so that its
``f_⋆`` is the fraction of the members of the class bound in pairs, whatever
their partner. `f_b_q` is the fraction of binaries with q > `q_min` among the
systems whose single or primary mass lies in `primary_window_msun`; with
stellar types a binary enters its numerator only when both components are of
class `:main_sequence`, since a photometric binary needs a luminous
main-sequence companion (without types every binary with q > `q_min` does).
On `:main_sequence` rows it is the analogue of the q-restricted fractions
measured photometrically in globular clusters (Milone et al. 2012,
DOI 10.1051/0004-6361/201016384). Their mass window is set relative to the
turn-off of each cluster, so the fixed default suits old globular clusters
only; younger populations need a window of their own.

A fraction is meaningful only together with the period range of its level:
level E is restricted to the pairs the integrator regularised, level P extends
to every pair bound against its nearest perturbers, and observed fractions
refer to their own period windows (Hurley et al. 2005,
DOI 10.1111/j.1365-2966.2005.09448.x). Photometric fractions of unresolved
systems, which include chance superpositions, are not computed.

Throws `ArgumentError` for an empty `snapshots`, empty or unknown `levels`, an
unknown `radius_frame`, invalid `edges` (any `edges` with `:milone`), an
unknown `centre`, a `kstar` not aligned with `snapshots` or leaving more than
`max_untyped_fraction` of the massive stars of a snapshot untyped, `q_min`
outside (0, 1), a window with lo ≥ hi or lo < 0, a `projection_axis` other
than 1, 2 or 3, `r_tide_nb ≤ 0`, and for a star in two pairs.

# Example
```julia
snaps = read_all_conf3(dir)
rows = binary_population_series(snaps; levels = (:ks,), radius_frame = :rh)
write_binary_population(joinpath(dir, "binary_population.csv"), rows)
```
"""
function binary_population_series(
    snapshots::AbstractVector{Snapshot};
    levels = (:ks, :bound),
    radius_frame::Symbol = :lagrangian,
    edges = nothing,
    centre = :shrinking_sphere,
    centre_kwargs::NamedTuple = (;),
    kstar = nothing,
    q_min::Real = 0.5,
    primary_window_msun::Tuple{<:Real,<:Real} = (0.47, 0.76),
    k::Integer = 10,
    n_perturbers::Integer = 4,
    gamma_max::Real = 0.1,
    max_untyped_fraction::Real = 0.01,
    projection_axis::Integer = 3,
    r_tide_nb::Real = NaN,
)::Vector{BinaryPopulationRow}
    args = _bpop_arguments(
        snapshots;
        levels,
        radius_frame,
        edges,
        centre,
        kstar,
        k,
        n_perturbers,
        gamma_max,
        max_untyped_fraction,
        projection_axis,
        r_tide_nb,
    )
    kmaps = args.kmaps
    0 < q_min < 1 ||
        throw(ArgumentError("binary_population_series: q_min must lie in (0, 1); got $q_min"))
    window = (Float64(primary_window_msun[1]), Float64(primary_window_msun[2]))
    (window[1] ≥ 0 && window[1] < window[2]) || throw(
        ArgumentError(
            "binary_population_series: primary_window_msun must satisfy 0 ≤ lo < hi; got $primary_window_msun",
        ),
    )
    q_thr = Float64(q_min)
    untyped = args.untyped
    # The classes present anywhere in the sequence, so that every snapshot has
    # the same rows (a class that dies out keeps its zero rows)
    classes = Int[]
    if kmaps !== nothing
        present = falses(length(STELLAR_CLASSES))
        for (j, snap) in enumerate(snapshots), c in _bpop_classes(snap, kmaps[j], j, untyped)
            c > 0 && (present[c] = true)
        end
        classes = findall(present)
    end

    rows = BinaryPopulationRow[]
    n_removed = 0
    for (j, snap) in enumerate(snapshots)
        mass = Float64.(snap.mass)
        pos = Float64.(snap.pos)
        vel = Float64.(snap.vel)
        c = _bpop_centre(pos, mass, args.c_fixed, centre_kwargs)
        kcls = kmaps === nothing ? Int[] : _bpop_classes(snap, kmaps[j], j, untyped)
        t_nb = time_nb(snap.header)
        t_myr = time_myr(snap.header)
        m_unit = zmbar(snap.header)
        for level in args.levels
            sys, n_out = _bpop_level(snap, mass, pos, vel, c, level, j, args)
            n_removed += n_out
            shells = _bpop_shells(sys.radius, sys.mass, args.edges, args.frame)
            _bpop_rows!(
                rows,
                sys,
                shells,
                mass,
                kcls,
                classes,
                j,
                t_nb,
                t_myr,
                level,
                args.frame,
                m_unit,
                window,
                q_thr,
            )
        end
    end
    _bpop_log_tidal("Binary population", args.r_tide, n_removed)
    return rows
end

# Arguments shared by the binary population, pair-type and distribution series,
# validated; `c_fixed` is `nothing` for a shrinking-sphere centre per snapshot.
function _bpop_arguments(
    snapshots::AbstractVector{Snapshot};
    levels,
    radius_frame::Symbol = :lagrangian,
    edges = nothing,
    centre,
    kstar = nothing,
    k::Integer,
    n_perturbers::Integer,
    gamma_max::Real,
    max_untyped_fraction::Real = 0.01,
    projection_axis::Integer = 3,
    r_tide_nb::Real = NaN,
)
    isempty(snapshots) && throw(ArgumentError("binary_population_series: no snapshots"))
    levs = _bpop_levels(levels)
    radius_frame in _BPOP_FRAMES || throw(
        ArgumentError(
            "binary_population_series: radius_frame must be :lagrangian, :rh, :projected_rh or :milone; got :$radius_frame",
        ),
    )
    shell_edges = _bpop_edges(edges, radius_frame)
    c_fixed = _bpop_fixed_centre(centre)
    kmaps = _bpop_kstar_maps(kstar, length(snapshots))
    0 ≤ max_untyped_fraction ≤ 1 || throw(
        ArgumentError(
            "binary_population_series: max_untyped_fraction must lie in [0, 1]; got $max_untyped_fraction",
        ),
    )
    projection_axis in (1, 2, 3) || throw(
        ArgumentError(
            "binary_population_series: projection_axis must be 1, 2 or 3; got $projection_axis",
        ),
    )
    (isnan(r_tide_nb) || r_tide_nb > 0) || throw(
        ArgumentError("binary_population_series: r_tide_nb must be > 0 or NaN; got $r_tide_nb"),
    )
    return (;
        levels = levs,
        frame = radius_frame,
        edges = shell_edges,
        c_fixed,
        kmaps,
        untyped = Float64(max_untyped_fraction),
        axis = Int(projection_axis),
        r_tide = Float64(r_tide_nb),
        k,
        n_perturbers,
        gamma_max,
    )
end

function _bpop_centre(
    pos::Matrix{Float64},
    mass::Vector{Float64},
    c_fixed::Union{Nothing,Vector{Float64}},
    centre_kwargs::NamedTuple,
)::Vector{Float64}
    return c_fixed === nothing ? shrinking_sphere_centre(pos, mass; centre_kwargs...).centre :
           c_fixed
end

# Systems of one level of snapshot j inside the tidal radius, and the number of
# systems removed beyond it
function _bpop_level(
    snap::Snapshot,
    mass::Vector{Float64},
    pos::Matrix{Float64},
    vel::Matrix{Float64},
    c::Vector{Float64},
    level::Symbol,
    j::Int,
    args::NamedTuple,
)::Tuple{_BpopSystems,Int}
    comps, hierarchy = _bpop_components(snap, level, args.k, args.n_perturbers, args.gamma_max)
    sys = _bpop_systems(snap, mass, pos, vel, c, comps, hierarchy, j, args.frame, args.axis)
    return _bpop_tidal_filter(sys, args.r_tide)
end

function _bpop_tidal_filter(sys::_BpopSystems, r_tide::Float64)::Tuple{_BpopSystems,Int}
    isfinite(r_tide) || return (sys, 0)
    keep = sys.r_3d .≤ r_tide
    n_out = count(!, keep)
    n_out == 0 && return (sys, 0)
    return (_bpop_select(sys, keep), n_out)
end

_bpop_select(sys::_BpopSystems, keep::BitVector)::_BpopSystems = _BpopSystems(
    sys.mass[keep],
    sys.radius[keep],
    sys.primary[keep],
    sys.secondary[keep],
    sys.hierarchy[keep],
    sys.a_nb[keep],
    sys.e[keep],
    sys.eb_nb[keep],
    sys.r_3d[keep],
    sys.velocity[:, keep],
)

function _bpop_log_tidal(label::String, r_tide::Float64, n_removed::Int)
    isfinite(r_tide) || return nothing
    @info "$label: tidal filter r_tide = $r_tide NB removed $n_removed systems in total " *
          "over all snapshots and levels"
    return nothing
end

# Hard/soft boundary |E_b| = ⟨m⟩σ² (NB) of a set of systems: ⟨m⟩ the mean mass
# and σ the one-dimensional velocity dispersion of the systems inside their
# half-mass radius (three-dimensional); NaN without systems.
function _bpop_hard_soft_energy(sys::_BpopSystems)::Float64
    r_h = _bpop_half_mass_radius(sys.r_3d, sys.mass)
    isnan(r_h) && return NaN
    inside = findall(≤(r_h), sys.r_3d)
    n_in = length(inside)
    m_sum = 0.0
    v̄ = zeros(3)
    @inbounds for s in inside
        m_sum += sys.mass[s]
        for d in 1:3
            v̄[d] += sys.velocity[d, s]
        end
    end
    v̄ ./= n_in
    ss = 0.0
    @inbounds for s in inside, d in 1:3
        ss += (sys.velocity[d, s] - v̄[d])^2
    end
    σ² = ss / (3n_in)
    return m_sum / n_in * σ²
end

function _bpop_levels(levels)::Vector{Symbol}
    levs = Symbol[]
    for l in (levels isa Symbol ? (levels,) : levels)
        (l isa Symbol && l in _BPOP_LEVELS && !(l in levs)) || throw(
            ArgumentError(
                "binary_population_series: levels must be distinct entries of (:ks, :bound); got $(repr(levels))",
            ),
        )
        push!(levs, l)
    end
    isempty(levs) && throw(ArgumentError("binary_population_series: levels is empty"))
    return levs
end

function _bpop_edges(edges, frame::Symbol)::Vector{Float64}
    if frame === :milone
        edges === nothing || throw(
            ArgumentError(
                "binary_population_series: the :milone frame has fixed regions and takes no edges; got $(repr(edges))",
            ),
        )
        return collect(Float64, _BPOP_EDGES_MILONE)
    end
    e = if edges === nothing
        collect(Float64, frame === :lagrangian ? _BPOP_EDGES_LAGRANGIAN : _BPOP_EDGES_RH)
    else
        collect(Float64, edges)
    end
    length(e) ≥ 2 ||
        throw(ArgumentError("binary_population_series: edges needs at least two entries; got $e"))
    e[1] == 0 || throw(ArgumentError("binary_population_series: edges must start at 0; got $e"))
    all(i -> e[i + 1] > e[i], 1:(length(e) - 1)) ||
        throw(ArgumentError("binary_population_series: edges must increase strictly; got $e"))
    frame === :lagrangian &&
        e[end] != 1 &&
        throw(
            ArgumentError(
                "binary_population_series: Lagrangian edges are mass fractions ending at 1; got $e",
            ),
        )
    return e
end

# The centre shared by all snapshots, or `nothing` for a shrinking-sphere centre
# computed per snapshot.
function _bpop_fixed_centre(centre)::Union{Nothing,Vector{Float64}}
    centre === :origin && return zeros(3)
    centre === :shrinking_sphere && return nothing
    (centre isa AbstractVector{<:Real} && length(centre) == 3) && return Vector{Float64}(centre)
    throw(
        ArgumentError(
            "binary_population_series: centre must be :origin, :shrinking_sphere or a length-3 vector; got $(repr(centre))",
        ),
    )
end

function _bpop_kstar_maps(kstar, n::Int)
    kstar === nothing && return nothing
    (kstar isa AbstractVector && length(kstar) == n) || throw(
        ArgumentError(
            "binary_population_series: kstar must be a vector of $n name → K* maps, one per snapshot",
        ),
    )
    return collect(kstar)
end

# Class index of every star, 0 for a star without a stellar type. The engine
# types neither the outer body of a merged hierarchy nor the representative of
# a chain subsystem in sev.83/bev.82; such stars count under `:all` only. More
# untyped stars than `max_untyped` (a fraction of the stars) means mismatched
# epochs and is an error.
function _bpop_classes(snap::Snapshot, kmap, j::Int, max_untyped::Float64)::Vector{Int}
    kmap isa AbstractDict ||
        throw(ArgumentError("binary_population_series: kstar[$j] must be a name → K* dictionary"))
    cls = Vector{Int}(undef, length(snap.name))
    n_untyped = 0
    n_massive = 0
    first_untyped = 0
    for i in eachindex(snap.name)
        nm = Int(snap.name[i])
        kt = get(kmap, nm, nothing)
        # Massless ghosts are no systems: neither typed nor counted
        snap.mass[i] == 0 || (n_massive += 1)
        if kt === nothing
            cls[i] = 0
            snap.mass[i] == 0 && continue
            n_untyped += 1
            first_untyped == 0 && (first_untyped = nm)
            continue
        end
        cls[i] = stellar_class_index(kt)
    end
    n_untyped > max_untyped * n_massive && throw(
        ArgumentError(
            "binary_population_series: $n_untyped of the $n_massive stars of snapshot $j " *
            "(e.g. $first_untyped) have no stellar type in kstar, more than the fraction " *
            "max_untyped_fraction = $max_untyped",
        ),
    )
    return cls
end

function _bpop_components(
    snap::Snapshot,
    level::Symbol,
    k::Integer,
    n_perturbers::Integer,
    gamma_max::Real,
)::Tuple{Vector{NTuple{2,Int}},Vector{Bool}}
    comps = NTuple{2,Int}[]
    hierarchy = Bool[]
    # A pair with a massless component is no binary
    massive(a, b) = snap.mass[a] != 0 && snap.mass[b] != 0
    if level === :ks
        n_ks = 2 * snap.npairs
        claimed = falses(length(snap.mass))
        for p in 1:(snap.npairs)
            a, b = 2p - 1, 2p
            if isempty(snap.pair_name) || snap.pair_name[p] > 0
                massive(a, b) || continue
                push!(comps, (a, b))
                push!(hierarchy, false)
                continue
            end
            # Merged hierarchy: 2k − 1 is the first inner component, 2k the
            # outer body (a single); the second inner component is the nearest
            # massive star outside the KS slots
            j = _bpop_nearest_single(snap, a, n_ks, claimed)
            (j == 0 || !massive(a, j)) && continue
            claimed[j] = true
            push!(comps, (a, j))
            push!(hierarchy, true)
        end
        return comps, hierarchy
    end
    pairs = bound_pairs(snap; k, n_perturbers, gamma_max)
    index = Dict{Int,Int}(Int(snap.name[i]) => i for i in eachindex(snap.name))
    for p in pairs
        a, b = index[p.name1], index[p.name2]
        massive(a, b) || continue
        push!(comps, (a, b))
        push!(hierarchy, p.hierarchy)
    end
    return comps, hierarchy
end

# Index of the massive, unclaimed star beyond the KS components (index > n_ks)
# nearest to star `a` in the Float64 positions, ties to the lower index; 0 when none.
function _bpop_nearest_single(snap::Snapshot, a::Int, n_ks::Int, claimed::BitVector)::Int
    best = 0
    d_best = Inf
    @inbounds for j in (n_ks + 1):length(snap.mass)
        (snap.mass[j] == 0 || claimed[j]) && continue
        d2 = 0.0
        for d in 1:3
            d2 += (Float64(snap.pos[d, j]) - Float64(snap.pos[d, a]))^2
        end
        if d2 < d_best
            best = j
            d_best = d2
        end
    end
    return best
end

function _bpop_systems(
    snap::Snapshot,
    mass::Vector{Float64},
    pos::Matrix{Float64},
    vel::Matrix{Float64},
    c::Vector{Float64},
    comps::Vector{NTuple{2,Int}},
    hierarchy::Vector{Bool},
    j::Int,
    frame::Symbol,
    axis::Int,
)::_BpopSystems
    n = length(mass)
    paired = falses(n)
    for (a, b) in comps
        for i in (a, b)
            paired[i] && throw(
                ArgumentError(
                    "binary_population_series: star $(snap.name[i]) of snapshot $j is a component of two pairs",
                ),
            )
            paired[i] = true
        end
    end
    n_pairs = length(comps)
    # Massless singles are ghosts of merged or chain subsystems, not systems
    n_ghosts = count(i -> !paired[i] && mass[i] == 0, 1:n)
    n_sys = n - n_pairs - n_ghosts
    m = Vector{Float64}(undef, n_sys)
    r = Vector{Float64}(undef, n_sys)
    primary = Vector{Int}(undef, n_sys)
    secondary = zeros(Int, n_sys)
    hier = fill(false, n_sys)
    projected = frame === :projected_rh || frame === :milone
    # Without projection the frame radius is the three-dimensional distance
    r3 = projected ? Vector{Float64}(undef, n_sys) : r
    v = Matrix{Float64}(undef, 3, n_sys)
    a_nb = fill(NaN, n_sys)
    ecc = fill(NaN, n_sys)
    eb = fill(NaN, n_sys)
    @inbounds for (s, (a, b)) in enumerate(comps)
        m_s = mass[a] + mass[b]
        d2 = 0.0
        d2p = 0.0
        for d in 1:3
            x = (mass[a] * pos[d, a] + mass[b] * pos[d, b]) / m_s
            δ² = (x - c[d])^2
            d2 += δ²
            d == axis || (d2p += δ²)
            v[d, s] = (mass[a] * vel[d, a] + mass[b] * vel[d, b]) / m_s
        end
        m[s] = m_s
        r3[s] = sqrt(d2)
        projected && (r[s] = sqrt(d2p))
        primary[s], secondary[s] = mass[b] > mass[a] ? (b, a) : (a, b)
        hier[s] = hierarchy[s]
        _, _, a_rel, e_rel = _bound_orbit(pos, vel, mass, a, b)
        if a_rel > 0
            a_nb[s] = a_rel
            ecc[s] = e_rel
            eb[s] = mass[a] * mass[b] / (2a_rel)
        end
    end
    s = n_pairs
    @inbounds for i in 1:n
        (paired[i] || mass[i] == 0) && continue
        s += 1
        d2 = 0.0
        d2p = 0.0
        for d in 1:3
            δ² = (pos[d, i] - c[d])^2
            d2 += δ²
            d == axis || (d2p += δ²)
            v[d, s] = vel[d, i]
        end
        m[s] = mass[i]
        r3[s] = sqrt(d2)
        projected && (r[s] = sqrt(d2p))
        primary[s] = i
    end
    return _BpopSystems(m, r, primary, secondary, hier, a_nb, ecc, eb, r3, v)
end

# Projected radius at which the projected number density of the systems falls
# to half its central value, from annuli of max(20, ⌈N/100⌉) systems; NaN when
# it does not fall to half or with too few systems.
function _bpop_core_radius(radius::Vector{Float64})::Float64
    n = length(radius)
    n_a = max(_BPOP_CORE_MIN_SYSTEMS, ceil(Int, 0.01n))
    n < 2n_a && return NaN
    R = sort(radius)
    R[n_a] > 0 || return NaN
    Σ_c = n_a / (π * R[n_a]^2)
    half = Σ_c / 2
    R_prev = R[n_a] / sqrt(2)
    Σ_prev = Σ_c
    lo = n_a
    while lo + n_a ≤ n
        hi = lo + n_a
        area = π * (R[hi]^2 - R[lo]^2)
        Σ = area > 0 ? n_a / area : Inf
        R_mid = sqrt((R[lo]^2 + R[hi]^2) / 2)
        if Σ ≤ half
            isfinite(Σ_prev) || return R_mid
            return R_prev + (Σ_prev - half) / (Σ_prev - Σ) * (R_mid - R_prev)
        end
        R_prev = R_mid
        Σ_prev = Σ
        lo = hi
    end
    return NaN
end

function _bpop_half_mass_radius(radius::Vector{Float64}, mass::Vector{Float64})::Float64
    isempty(radius) && return NaN
    order = sortperm(radius)
    cum = cumsum(mass[order]) ./ sum(mass)
    i = findfirst(≥(0.5), cum)
    return i === nothing ? NaN : radius[order[i]]
end

_bpop_scaled(e::Float64, r_h::Float64)::Float64 = isinf(e) ? e : e * r_h

# Shell of every system (0 when outside all shells) and the edges and radii of
# bin 0 and the shells, indexed bin + 1. The projected frames take projected
# radii in `radius`; `:milone` replaces its NaN inner edge by R_c/R_h.
function _bpop_shells(
    radius::Vector{Float64},
    mass::Vector{Float64},
    edges::Vector{Float64},
    frame::Symbol,
)
    n = length(radius)
    nb = length(edges) - 1
    shell = zeros(Int, n)
    milone = frame === :milone
    r_h = frame === :lagrangian ? NaN : _bpop_half_mass_radius(radius, mass)
    r_c = milone ? _bpop_core_radius(radius) : NaN
    e = milone ? [0.0; r_c / r_h; edges[3:end]] : edges
    edge_lo = [0.0; e[1:nb]]
    edge_hi = [frame === :lagrangian ? 1.0 : Inf; e[2:end]]
    r_lo = zeros(nb + 1)
    r_hi = fill(NaN, nb + 1)
    r_hi[1] = n == 0 ? NaN : maximum(radius)
    if frame !== :lagrangian
        for b in 1:nb
            r_lo[b + 1] = _bpop_scaled(e[b], r_h)
            r_hi[b + 1] = _bpop_scaled(e[b + 1], r_h)
        end
        if milone
            r_hi[2] = r_c
            r_lo[3] = r_c
        end
        for s in 1:n
            r = radius[s]
            for b in 1:nb
                if r_lo[b + 1] ≤ r < r_hi[b + 1]
                    shell[s] = b
                    break
                end
            end
        end
    else
        order = sortperm(radius)
        cum = cumsum(mass[order]) ./ sum(mass)
        # e_b < cum ≤ e_{b+1}; round-off beyond the last edge joins the last shell
        for (p, s) in enumerate(order)
            b = clamp(searchsortedfirst(edges, cum[p]) - 1, 1, nb)
            shell[s] = b
            r_hi[b + 1] = radius[s]
        end
        prev = 0.0
        for b in 1:nb
            r_lo[b + 1] = prev
            isnan(r_hi[b + 1]) || (prev = r_hi[b + 1])
        end
    end
    return (; shell, edge_lo, edge_hi, r_lo, r_hi)
end

function _bpop_add_single!(t::_BpopTally, row::Int, col::Int, m::Float64, in_window::Bool)
    t.n_singles[row, col] += 1
    t.m_singles[row, col] += m
    in_window && (t.n_window[row, col] += 1)
    return nothing
end

function _bpop_add_binary!(
    t::_BpopTally,
    row::Int,
    col::Int,
    hierarchy::Bool,
    in_window::Bool,
    q_above::Bool,
)
    t.n_binaries[row, col] += 1
    hierarchy && (t.n_hierarchies[row, col] += 1)
    if in_window
        t.n_window[row, col] += 1
        q_above && (t.n_window_q[row, col] += 1)
    end
    return nothing
end

function _bpop_tally!(
    t::_BpopTally,
    row::Int,
    sys::_BpopSystems,
    s::Int,
    mass::Vector{Float64},
    kcls::Vector{Int},
    colof::Vector{Int},
    m_unit::Float64,
    window::NTuple{2,Float64},
    q_min::Float64,
    ms::Int,
)
    p = sys.primary[s]
    sec = sys.secondary[s]
    in_window = window[1] ≤ m_unit * mass[p] ≤ window[2]
    if sec == 0
        _bpop_add_single!(t, row, 1, mass[p], in_window)
        col = _bpop_column(colof, kcls, p)
        col > 0 && _bpop_add_single!(t, row, col, mass[p], in_window)
        return nothing
    end
    # With stellar types a photometric binary needs two main-sequence components
    q_above = mass[sec] / mass[p] > q_min && (isempty(kcls) || (kcls[p] == ms && kcls[sec] == ms))
    _bpop_add_binary!(t, row, 1, sys.hierarchy[s], in_window, q_above)
    t.n_in_binaries[row, 1] += 2
    t.m_binaries[row, 1] += sys.mass[s]
    col = _bpop_column(colof, kcls, p)
    col > 0 && _bpop_add_binary!(t, row, col, sys.hierarchy[s], in_window, q_above)
    for i in (p, sec)
        col = _bpop_column(colof, kcls, i)
        col > 0 || continue
        t.n_in_binaries[row, col] += 1
        t.m_binaries[row, col] += mass[i]
    end
    return nothing
end

# Class column of star `i`; 0 without classes or for an unclassified star
_bpop_column(colof::Vector{Int}, kcls::Vector{Int}, i::Int)::Int =
    (isempty(kcls) || kcls[i] == 0) ? 0 : colof[kcls[i]]

_bpop_ratio(a::Real, d::Real)::Float64 = d == 0 ? NaN : Float64(a) / Float64(d)

function _bpop_rows!(
    rows::Vector{BinaryPopulationRow},
    sys::_BpopSystems,
    shells,
    mass::Vector{Float64},
    kcls::Vector{Int},
    classes::Vector{Int},
    j::Int,
    t_nb::Float64,
    t_myr::Float64,
    level::Symbol,
    frame::Symbol,
    m_unit::Float64,
    window::NTuple{2,Float64},
    q_min::Float64,
)
    n_cls = length(STELLAR_CLASSES)
    cols = isempty(kcls) ? [0] : [0; classes]
    colof = zeros(Int, n_cls)
    for (q, c) in enumerate(cols)
        c > 0 && (colof[c] = q)
    end
    n_rows = length(shells.edge_lo)
    n_cols = length(cols)
    t = _BpopTally(n_rows, n_cols)
    singles = [Int[] for _ in 1:n_cols]
    binaries = [Int[] for _ in 1:n_cols]
    ms = something(findfirst(c -> c.key === :main_sequence, STELLAR_CLASSES), 0)
    for s in eachindex(sys.mass)
        b = shells.shell[s]
        _bpop_tally!(t, 1, sys, s, mass, kcls, colof, m_unit, window, q_min, ms)
        b > 0 && _bpop_tally!(t, b + 1, sys, s, mass, kcls, colof, m_unit, window, q_min, ms)
        members = sys.secondary[s] == 0 ? singles : binaries
        push!(members[1], s)
        col = _bpop_column(colof, kcls, sys.primary[s])
        col > 0 && push!(members[col], s)
    end
    r_h_bin = [_bpop_half_mass_radius(sys.radius[l], sys.mass[l]) for l in binaries]
    r_h_sin = [_bpop_half_mass_radius(sys.radius[l], sys.mass[l]) for l in singles]
    for row in 1:n_rows, q in 1:n_cols
        key = cols[q] == 0 ? :all : STELLAR_CLASSES[cols[q]].key
        n_s = t.n_singles[row, q]
        n_b = t.n_binaries[row, q]
        n_in = t.n_in_binaries[row, q]
        m_s = t.m_singles[row, q] * m_unit
        m_b = t.m_binaries[row, q] * m_unit
        f = _bpop_ratio(n_b, n_s + n_b)
        push!(
            rows,
            BinaryPopulationRow(
                j,
                t_nb,
                t_myr,
                level,
                frame,
                row - 1,
                shells.edge_lo[row],
                shells.edge_hi[row],
                shells.r_lo[row],
                shells.r_hi[row],
                key,
                n_s,
                n_b,
                n_in,
                t.n_hierarchies[row, q],
                m_s,
                m_b,
                f,
                _bpop_ratio(n_in, n_s + n_in),
                _bpop_ratio(m_b, m_s + m_b),
                n_s + n_b == 0 ? NaN : sqrt(f * (1 - f) / (n_s + n_b)),
                _bpop_ratio(t.n_window_q[row, q], t.n_window[row, q]),
                row == 1 ? r_h_bin[q] : NaN,
                row == 1 ? r_h_sin[q] : NaN,
            ),
        )
    end
    return rows
end

"""
    write_binary_population(path, rows::AbstractVector{BinaryPopulationRow}) -> String

Write the rows of [`binary_population_series`](@ref) as CSV, one column per
field of [`BinaryPopulationRow`](@ref) in declaration order (symbols as their
names, non-finite values as empty cells). An existing file is backed up, never
overwritten. Returns the path.
"""
function write_binary_population(
    path::AbstractString,
    rows::AbstractVector{BinaryPopulationRow},
)::String
    _backup_existing(path)
    n_fields = fieldcount(BinaryPopulationRow)
    open(path, "w") do io
        println(io, join(fieldnames(BinaryPopulationRow), ","))
        for r in rows
            println(io, join(ntuple(f -> _bpop_cell(getfield(r, f)), n_fields), ","))
        end
    end
    return String(path)
end

_bpop_cell(v::Symbol) = String(v)
_bpop_cell(v) = _csv_precise(v)

function _bpop_write_rows(path::AbstractString, rows::AbstractVector{T})::String where {T}
    _backup_existing(path)
    n_fields = fieldcount(T)
    open(path, "w") do io
        println(io, join(fieldnames(T), ","))
        for r in rows
            println(io, join(ntuple(f -> _bpop_cell(getfield(r, f)), n_fields), ","))
        end
    end
    return String(path)
end

_bpop_median(v::Vector{Float64})::Float64 = isempty(v) ? NaN : _median(v)

"""
    BinaryPairTypeRow

Pair-type composition and hardness of the binaries of one snapshot, census
level and radial bin from [`binary_pair_type_series`](@ref). Lengths, masses
and energies are in the N-body units of the snapshot (G = 1).

| field | meaning | unit |
|---|---|---|
| `snapshot` | index of the snapshot in the input vector | — |
| `t_nb`, `t_myr` | snapshot time | NB, Myr |
| `level` | census level: `:ks` (E, KS pairs) or `:bound` (P, bound pairs) | — |
| `radius_frame` | `:lagrangian`, `:rh`, `:projected_rh` or `:milone`, see [`binary_population_series`](@ref) | — |
| `bin` | 0 for all systems, else the shell index 1, 2, … | — |
| `edge_lo`, `edge_hi` | shell edges, as in [`BinaryPopulationRow`](@ref) | — |
| `r_lo_nb`, `r_hi_nb` | shell radii, as in [`BinaryPopulationRow`](@ref) | NB |
| `n_binaries` | binaries of the row | — |
| `n_ms_ms`, `n_ms_evolved`, `n_ms_wd`, `n_wd_wd` | main sequence–main sequence, main sequence–evolved star, main sequence–white dwarf, white dwarf–white dwarf | — |
| `n_ns_any` | pairs with at least one neutron star | — |
| `n_bh_bh`, `n_bh_star` | black hole–black hole; black hole with a main-sequence or evolved star | — |
| `n_other` | every other typed pair (evolved–evolved, evolved–white dwarf, black hole–white dwarf) | — |
| `n_untyped` | pairs with a component without a stellar type | — |
| `n_hard`, `n_soft` | binaries with ``\\lvert E_b \\rvert`` above and at or below `e_hs_nb`; pairs without orbital elements count in neither | — |
| `e_hs_nb` | hard/soft boundary ``⟨m⟩σ²`` of the snapshot and level | NB |
| `median_a_nb`, `median_e`, `median_q` | medians of the semi-major axis, eccentricity and mass ratio over the binaries of the row with a finite value; `NaN` without any | NB, —, — |
| `f_hard` | `n_hard / n_binaries` | — |

The nine type columns partition `n_binaries`. Medians and fractions of an
empty row are `NaN`.
"""
struct BinaryPairTypeRow
    snapshot::Int
    t_nb::Float64
    t_myr::Float64
    level::Symbol
    radius_frame::Symbol
    bin::Int
    edge_lo::Float64
    edge_hi::Float64
    r_lo_nb::Float64
    r_hi_nb::Float64
    n_binaries::Int
    n_ms_ms::Int
    n_ms_evolved::Int
    n_ms_wd::Int
    n_wd_wd::Int
    n_ns_any::Int
    n_bh_bh::Int
    n_bh_star::Int
    n_other::Int
    n_untyped::Int
    n_hard::Int
    n_soft::Int
    e_hs_nb::Float64
    median_a_nb::Float64
    median_e::Float64
    median_q::Float64
    f_hard::Float64
end

# Component kinds by K*: untyped, MS (−1…1), evolved (2…9), WD (10…12), NS (13),
# BH (14), massless remnant (15)
const _BPOP_KIND_UNTYPED = 0
const _BPOP_KIND_MS = 1
const _BPOP_KIND_EVOLVED = 2
const _BPOP_KIND_WD = 3
const _BPOP_KIND_NS = 4
const _BPOP_KIND_BH = 5
const _BPOP_KIND_REMNANT = 6

# Pair-type columns of BinaryPairTypeRow, in field order
const _BPOP_PAIR_MS_MS = 1
const _BPOP_PAIR_MS_EVOLVED = 2
const _BPOP_PAIR_MS_WD = 3
const _BPOP_PAIR_WD_WD = 4
const _BPOP_PAIR_NS_ANY = 5
const _BPOP_PAIR_BH_BH = 6
const _BPOP_PAIR_BH_STAR = 7
const _BPOP_PAIR_OTHER = 8
const _BPOP_PAIR_UNTYPED = 9
const _BPOP_N_PAIR_TYPES = 9

# Kind of a star of class index `c` (0: untyped) of STELLAR_CLASSES
function _bpop_kind(c::Int)::Int
    c == 0 && return _BPOP_KIND_UNTYPED
    kt = first(STELLAR_CLASSES[c].kstar)
    kt ≤ 1 && return _BPOP_KIND_MS
    kt ≤ 9 && return _BPOP_KIND_EVOLVED
    kt ≤ 12 && return _BPOP_KIND_WD
    kt == 13 && return _BPOP_KIND_NS
    kt == 14 && return _BPOP_KIND_BH
    return _BPOP_KIND_REMNANT
end

# Pair-type column of a pair of component kinds; the rules apply in this order
function _bpop_pair_type(k1::Int, k2::Int)::Int
    (k1 == _BPOP_KIND_UNTYPED || k2 == _BPOP_KIND_UNTYPED) && return _BPOP_PAIR_UNTYPED
    (k1 == _BPOP_KIND_NS || k2 == _BPOP_KIND_NS) && return _BPOP_PAIR_NS_ANY
    (k1 == _BPOP_KIND_BH && k2 == _BPOP_KIND_BH) && return _BPOP_PAIR_BH_BH
    if k1 == _BPOP_KIND_BH || k2 == _BPOP_KIND_BH
        other = k1 == _BPOP_KIND_BH ? k2 : k1
        return other in (_BPOP_KIND_MS, _BPOP_KIND_EVOLVED) ? _BPOP_PAIR_BH_STAR : _BPOP_PAIR_OTHER
    end
    (k1 == _BPOP_KIND_WD && k2 == _BPOP_KIND_WD) && return _BPOP_PAIR_WD_WD
    lo, hi = minmax(k1, k2)
    lo == _BPOP_KIND_MS || return _BPOP_PAIR_OTHER
    hi == _BPOP_KIND_WD && return _BPOP_PAIR_MS_WD
    hi == _BPOP_KIND_MS && return _BPOP_PAIR_MS_MS
    hi == _BPOP_KIND_EVOLVED && return _BPOP_PAIR_MS_EVOLVED
    return _BPOP_PAIR_OTHER
end

"""
    binary_pair_type_series(snapshots::AbstractVector{Snapshot}; levels = (:ks, :bound),
                            radius_frame::Symbol = :lagrangian, edges = nothing,
                            centre = :shrinking_sphere, centre_kwargs::NamedTuple = (;),
                            kstar = nothing, k::Integer = 10, n_perturbers::Integer = 4,
                            gamma_max::Real = 0.1, max_untyped_fraction::Real = 0.01,
                            r_tide_nb::Real = NaN, projection_axis::Integer = 3)
        -> Vector{BinaryPairTypeRow}

Pair-type composition, hardness and median orbital elements of the binaries
of a snapshot sequence, one [`BinaryPairTypeRow`](@ref) per snapshot × level ×
bin, in that order, bin 0 holding every binary.

# Arguments
The keywords, the census levels, the systems, the radial frames, the centre
and the tidal filter are those of [`binary_population_series`](@ref).
Without `kstar` every binary counts in `n_untyped`.

# Method
Each component is classified by its stellar type K* with the ranges of
[`STELLAR_CLASSES`](@ref): main sequence (MS) K* ∈ {−1, 0, 1}, evolved 2–9,
white dwarf (WD) 10–12, neutron star (NS) 13, black hole (BH) 14. A pair
(k₁, k₂) takes the first matching rule: a component without a type →
`n_untyped`; any NS → `n_ns_any`; two BHs → `n_bh_bh`; one BH with an MS or
evolved star → `n_bh_star`; BH with WD → `n_other`; two WDs → `n_wd_wd`; MS
with WD → `n_ms_wd`; two MS stars → `n_ms_ms`; MS with an evolved star →
`n_ms_evolved`; any other pair → `n_other`.

The semi-major axis a, eccentricity e and binding energy
``E_b = m_1 m_2/(2a)`` follow from the relative orbit of the two components
in the snapshot (G = 1); pairs that are not bound have no elements. A binary
is hard when ``\\lvert E_b \\rvert > ⟨m⟩σ²`` (Heggie 1975,
DOI 10.1093/mnras/173.3.729), with ⟨m⟩ the mean mass and σ the
one-dimensional velocity dispersion,
``σ² = \\sum_i \\lvert \\mathbf{v}_i - \\bar{\\mathbf{v}} \\rvert^2 / (3n)``,
of the n systems of the level within their three-dimensional half-mass
radius, per snapshot and level and after the tidal filter; `e_hs_nb` is this
boundary. The mass ratio is q = m₂/m₁ ≤ 1.

Throws `ArgumentError` as [`binary_population_series`](@ref).

# Example
```julia
snaps = read_all_conf3(dir)
rows = binary_pair_type_series(snaps; levels = (:ks,), radius_frame = :milone,
                               kstar = kmaps)
write_binary_pair_types(joinpath(dir, "binary_pair_types.csv"), rows)
```
"""
function binary_pair_type_series(
    snapshots::AbstractVector{Snapshot};
    levels = (:ks, :bound),
    radius_frame::Symbol = :lagrangian,
    edges = nothing,
    centre = :shrinking_sphere,
    centre_kwargs::NamedTuple = (;),
    kstar = nothing,
    k::Integer = 10,
    n_perturbers::Integer = 4,
    gamma_max::Real = 0.1,
    max_untyped_fraction::Real = 0.01,
    r_tide_nb::Real = NaN,
    projection_axis::Integer = 3,
)::Vector{BinaryPairTypeRow}
    args = _bpop_arguments(
        snapshots;
        levels,
        radius_frame,
        edges,
        centre,
        kstar,
        k,
        n_perturbers,
        gamma_max,
        max_untyped_fraction,
        projection_axis,
        r_tide_nb,
    )
    rows = BinaryPairTypeRow[]
    n_removed = 0
    for (j, snap) in enumerate(snapshots)
        mass = Float64.(snap.mass)
        pos = Float64.(snap.pos)
        vel = Float64.(snap.vel)
        c = _bpop_centre(pos, mass, args.c_fixed, centre_kwargs)
        kind =
            args.kmaps === nothing ? Int[] :
            map(_bpop_kind, _bpop_classes(snap, args.kmaps[j], j, args.untyped))
        t_nb = time_nb(snap.header)
        t_myr = time_myr(snap.header)
        for level in args.levels
            sys, n_out = _bpop_level(snap, mass, pos, vel, c, level, j, args)
            n_removed += n_out
            shells = _bpop_shells(sys.radius, sys.mass, args.edges, args.frame)
            e_hs = _bpop_hard_soft_energy(sys)
            _bpop_pair_type_rows!(
                rows,
                sys,
                shells,
                mass,
                kind,
                e_hs,
                j,
                t_nb,
                t_myr,
                level,
                args.frame,
            )
        end
    end
    _bpop_log_tidal("Binary pair types", args.r_tide, n_removed)
    return rows
end

function _bpop_pair_type_rows!(
    rows::Vector{BinaryPairTypeRow},
    sys::_BpopSystems,
    shells,
    mass::Vector{Float64},
    kind::Vector{Int},
    e_hs::Float64,
    j::Int,
    t_nb::Float64,
    t_myr::Float64,
    level::Symbol,
    frame::Symbol,
)
    n_rows = length(shells.edge_lo)
    n_type = zeros(Int, n_rows, _BPOP_N_PAIR_TYPES)
    n_hard = zeros(Int, n_rows)
    n_soft = zeros(Int, n_rows)
    a_row = [Float64[] for _ in 1:n_rows]
    e_row = [Float64[] for _ in 1:n_rows]
    q_row = [Float64[] for _ in 1:n_rows]
    for s in eachindex(sys.mass)
        sec = sys.secondary[s]
        sec == 0 && continue
        p = sys.primary[s]
        col = isempty(kind) ? _BPOP_PAIR_UNTYPED : _bpop_pair_type(kind[p], kind[sec])
        eb = sys.eb_nb[s]
        b = shells.shell[s]
        # Row 1 is bin 0; a binary outside every shell enters it only
        for row in (1, b + 1)
            n_type[row, col] += 1
            if isfinite(eb) && isfinite(e_hs)
                eb > e_hs ? (n_hard[row] += 1) : (n_soft[row] += 1)
            end
            isfinite(sys.a_nb[s]) && push!(a_row[row], sys.a_nb[s])
            isfinite(sys.e[s]) && push!(e_row[row], sys.e[s])
            push!(q_row[row], mass[sec] / mass[p])
            b == 0 && break
        end
    end
    for row in 1:n_rows
        n_bin = sum(@view n_type[row, :])
        push!(
            rows,
            BinaryPairTypeRow(
                j,
                t_nb,
                t_myr,
                level,
                frame,
                row - 1,
                shells.edge_lo[row],
                shells.edge_hi[row],
                shells.r_lo[row],
                shells.r_hi[row],
                n_bin,
                n_type[row, _BPOP_PAIR_MS_MS],
                n_type[row, _BPOP_PAIR_MS_EVOLVED],
                n_type[row, _BPOP_PAIR_MS_WD],
                n_type[row, _BPOP_PAIR_WD_WD],
                n_type[row, _BPOP_PAIR_NS_ANY],
                n_type[row, _BPOP_PAIR_BH_BH],
                n_type[row, _BPOP_PAIR_BH_STAR],
                n_type[row, _BPOP_PAIR_OTHER],
                n_type[row, _BPOP_PAIR_UNTYPED],
                n_hard[row],
                n_soft[row],
                e_hs,
                _bpop_median(a_row[row]),
                _bpop_median(e_row[row]),
                _bpop_median(q_row[row]),
                _bpop_ratio(n_hard[row], n_bin),
            ),
        )
    end
    return rows
end

"""
    write_binary_pair_types(path, rows::AbstractVector{BinaryPairTypeRow}) -> String

Write the rows of [`binary_pair_type_series`](@ref) as CSV, one column per
field of [`BinaryPairTypeRow`](@ref) in declaration order (symbols as their
names, non-finite values as empty cells). An existing file is backed up, never
overwritten. Returns the path.
"""
write_binary_pair_types(path::AbstractString, rows::AbstractVector{BinaryPairTypeRow})::String =
    _bpop_write_rows(path, rows)

"""
    BinaryDistributionRow

One histogram bin of an orbital-element distribution of the binaries of one
snapshot and census level from [`binary_distribution_series`](@ref).

| field | meaning | unit |
|---|---|---|
| `snapshot` | index of the snapshot in the input vector | — |
| `t_nb`, `t_myr` | snapshot time | NB, Myr |
| `level` | census level: `:ks` (E, KS pairs) or `:bound` (P, bound pairs) | — |
| `quantity` | `:log_a_au` (log₁₀ a/au), `:e` (eccentricity), `:q` (mass ratio m₂/m₁) or `:log_eb` (log₁₀ of ``\\lvert E_b \\rvert`` over the hard/soft boundary ``⟨m⟩σ²``) | — |
| `bin` | histogram bin 1, 2, …; 0 for the values outside the edges and the binaries without the quantity | — |
| `edge_lo`, `edge_hi` | bin edges in the unit of `quantity`; `NaN` for bin 0 | — |
| `count` | binaries in the bin | — |
"""
struct BinaryDistributionRow
    snapshot::Int
    t_nb::Float64
    t_myr::Float64
    level::Symbol
    quantity::Symbol
    bin::Int
    edge_lo::Float64
    edge_hi::Float64
    count::Int
end

"""
    binary_distribution_series(snapshots::AbstractVector{Snapshot}; levels = (:ks, :bound),
                               centre = :shrinking_sphere, centre_kwargs::NamedTuple = (;),
                               k::Integer = 10, n_perturbers::Integer = 4,
                               gamma_max::Real = 0.1, r_tide_nb::Real = NaN,
                               edges_log_a_au = -3.0:0.25:3.0, edges_e = 0.0:0.1:1.0,
                               edges_q = 0.0:0.1:1.0, edges_log_eb = -2.0:0.25:4.0)
        -> Vector{BinaryDistributionRow}

Histograms of the semi-major axis, eccentricity, mass ratio and binding
energy of the binaries of a snapshot sequence, one
[`BinaryDistributionRow`](@ref) per snapshot × level × quantity × bin, in that
order, with every bin present (zeros included) and bin 0 first.

# Arguments
- `snapshots`, `levels`, `centre`, `centre_kwargs`, `k`, `n_perturbers`,
  `gamma_max`, `r_tide_nb`: as in [`binary_population_series`](@ref)
- `edges_log_a_au`, `edges_e`, `edges_q`, `edges_log_eb`: histogram edges,
  finite and strictly increasing, of log₁₀(a/au), e, q and
  log₁₀(``\\lvert E_b \\rvert``/⟨m⟩σ²)

# Method
The binaries and their elements a, e and ``E_b = m_1 m_2/(2a)`` are those of
[`binary_pair_type_series`](@ref); a is converted to au with
[`rbar`](@ref) and 1 pc = 206264.806 au, q = m₂/m₁ ≤ 1, and ``\\lvert E_b \\rvert`` is
measured in units of the hard/soft boundary ⟨m⟩σ² of the snapshot and level
(Heggie 1975, DOI 10.1093/mnras/173.3.729), so that `:log_eb` = 0 is the
boundary. Bin b holds the values ``e_b ≤ x < e_{b+1}``, the last bin also
``x = e_\\mathrm{last}``; values outside the edges and binaries without
elements (pairs not bound in the snapshot) are counted in bin 0. The counts of
a snapshot, level and quantity therefore sum to its number of binaries.

Throws `ArgumentError` for the arguments rejected by
[`binary_population_series`](@ref) and for invalid edges.

# Example
```julia
rows = binary_distribution_series(read_all_conf3(dir); levels = (:ks,))
write_binary_distributions(joinpath(dir, "binary_distributions.csv"), rows)
```
"""
function binary_distribution_series(
    snapshots::AbstractVector{Snapshot};
    levels = (:ks, :bound),
    centre = :shrinking_sphere,
    centre_kwargs::NamedTuple = (;),
    k::Integer = 10,
    n_perturbers::Integer = 4,
    gamma_max::Real = 0.1,
    r_tide_nb::Real = NaN,
    edges_log_a_au = -3.0:0.25:3.0,
    edges_e = 0.0:0.1:1.0,
    edges_q = 0.0:0.1:1.0,
    edges_log_eb = -2.0:0.25:4.0,
)::Vector{BinaryDistributionRow}
    args = _bpop_arguments(snapshots; levels, centre, k, n_perturbers, gamma_max, r_tide_nb)
    histograms = (
        (:log_a_au, _bpop_histogram_edges(edges_log_a_au, :edges_log_a_au)),
        (:e, _bpop_histogram_edges(edges_e, :edges_e)),
        (:q, _bpop_histogram_edges(edges_q, :edges_q)),
        (:log_eb, _bpop_histogram_edges(edges_log_eb, :edges_log_eb)),
    )
    rows = BinaryDistributionRow[]
    n_removed = 0
    for (j, snap) in enumerate(snapshots)
        mass = Float64.(snap.mass)
        pos = Float64.(snap.pos)
        vel = Float64.(snap.vel)
        c = _bpop_centre(pos, mass, args.c_fixed, centre_kwargs)
        t_nb = time_nb(snap.header)
        t_myr = time_myr(snap.header)
        au_per_nb = rbar(snap.header) / _AU_PC
        for level in args.levels
            sys, n_out = _bpop_level(snap, mass, pos, vel, c, level, j, args)
            n_removed += n_out
            e_hs = _bpop_hard_soft_energy(sys)
            binaries = findall(>(0), sys.secondary)
            for (quantity, edges) in histograms
                values = Float64[
                    _bpop_quantity(quantity, sys, s, mass, au_per_nb, e_hs) for s in binaries
                ]
                _bpop_histogram_rows!(rows, values, edges, j, t_nb, t_myr, level, quantity)
            end
        end
    end
    _bpop_log_tidal("Binary distributions", args.r_tide, n_removed)
    return rows
end

function _bpop_histogram_edges(edges, name::Symbol)::Vector{Float64}
    e = collect(Float64, edges)
    (length(e) ≥ 2 && all(isfinite, e) && all(i -> e[i + 1] > e[i], 1:(length(e) - 1))) || throw(
        ArgumentError(
            "binary_distribution_series: $name must hold at least two finite, strictly increasing edges; got $e",
        ),
    )
    return e
end

function _bpop_quantity(
    quantity::Symbol,
    sys::_BpopSystems,
    s::Int,
    mass::Vector{Float64},
    au_per_nb::Float64,
    e_hs::Float64,
)::Float64
    quantity === :log_a_au && return log10(sys.a_nb[s] * au_per_nb)
    quantity === :e && return sys.e[s]
    quantity === :q && return mass[sys.secondary[s]] / mass[sys.primary[s]]
    return log10(sys.eb_nb[s] / e_hs)
end

# Bin of `x` for the edges `e`: e_b ≤ x < e_{b+1}, the last bin closed; 0
# outside the edges or for NaN
function _bpop_histogram_bin(x::Float64, e::Vector{Float64})::Int
    (isnan(x) || x < e[1] || x > e[end]) && return 0
    x == e[end] && return length(e) - 1
    return searchsortedlast(e, x)
end

function _bpop_histogram_rows!(
    rows::Vector{BinaryDistributionRow},
    values::Vector{Float64},
    edges::Vector{Float64},
    j::Int,
    t_nb::Float64,
    t_myr::Float64,
    level::Symbol,
    quantity::Symbol,
)
    nb = length(edges) - 1
    counts = zeros(Int, nb + 1)
    for x in values
        counts[_bpop_histogram_bin(x, edges) + 1] += 1
    end
    push!(rows, BinaryDistributionRow(j, t_nb, t_myr, level, quantity, 0, NaN, NaN, counts[1]))
    for b in 1:nb
        push!(
            rows,
            BinaryDistributionRow(
                j,
                t_nb,
                t_myr,
                level,
                quantity,
                b,
                edges[b],
                edges[b + 1],
                counts[b + 1],
            ),
        )
    end
    return rows
end

"""
    write_binary_distributions(path, rows::AbstractVector{BinaryDistributionRow}) -> String

Write the rows of [`binary_distribution_series`](@ref) as CSV, one column per
field of [`BinaryDistributionRow`](@ref) in declaration order (symbols as
their names, non-finite values as empty cells). An existing file is backed up,
never overwritten. Returns the path.
"""
write_binary_distributions(
    path::AbstractString,
    rows::AbstractVector{BinaryDistributionRow},
)::String = _bpop_write_rows(path, rows)

"""
    stellar_type_map(sev::StellarEvolutionSnapshot,
                     bev::Union{Nothing,BinaryEvolutionSnapshot}) -> Dict{Int,Int}

Stellar type K* of every star of one epoch, keyed by name: the single stars of
the sev.83 records, then both components of every bev.82 record, which take
precedence (sev.83 omits the members of KS pairs). The `kstar` input of
[`binary_population_series`](@ref).

# Example
```julia
kmap = stellar_type_map(read_stellar_evolution(sev_path), read_binary_evolution(bev_path))
```
"""
function stellar_type_map(
    sev::StellarEvolutionSnapshot,
    bev::Union{Nothing,BinaryEvolutionSnapshot},
)::Dict{Int,Int}
    kmap = Dict{Int,Int}()
    for r in sev.records
        kmap[Int(r.name)] = Int(r.stellar_type)
    end
    if bev !== nothing
        for r in bev.records
            kmap[Int(r.name1)] = Int(r.stellar_type1)
            kmap[Int(r.name2)] = Int(r.stellar_type2)
        end
    end
    return kmap
end

# Epoch for the snapshot time `t` [NB] and its lag t − t_epoch [NB]: the epoch
# nearest to `t` within the matching tolerance (lag 0), else the latest epoch at
# or before `t`, since stellar types change on evolutionary time scales;
# `(nothing, NaN)` when there is none. Epochs are timed by their records' N-body
# time: the file headers carry the time in Myr with one decimal only. Epochs
# without records are skipped.
function _bpop_nearest_epoch(
    epochs::AbstractVector{<:Union{StellarEvolutionSnapshot,BinaryEvolutionSnapshot}},
    t::Float64,
)::Tuple{Union{Nothing,Int},Float64}
    best = 0
    d_best = Inf
    before = 0
    t_before = -Inf
    for (i, e) in enumerate(epochs)
        isempty(e.records) && continue
        te = e.records[1].time_nb
        d = abs(te - t)
        if d < d_best
            best = i
            d_best = d
        end
        if te ≤ t && te > t_before
            before = i
            t_before = te
        end
    end
    best > 0 && d_best ≤ _BPOP_EPOCH_RTOL * max(1.0, abs(t)) && return (best, 0.0)
    before > 0 && return (before, t - t_before)
    return (nothing, NaN)
end

# Name → K* maps aligned with `snapshots` from the sev.83 and bev.82 epochs of
# the same time or, failing that, the latest epochs before it (one @info with
# the largest lag); `nothing` (warned) when a snapshot has no sev.83 epoch at
# or before its time or, with KS pairs, no such bev.82 epoch (sev.83 omits the
# KS components).
function _bpop_pipeline_kstar(
    snapshots::AbstractVector{Snapshot},
    sevs::AbstractVector{StellarEvolutionSnapshot},
    bevs::AbstractVector{BinaryEvolutionSnapshot},
)::Union{Nothing,Vector{Dict{Int,Int}}}
    maps = Dict{Int,Int}[]
    n_lagged = 0
    lag_max = 0.0
    for (j, snap) in enumerate(snapshots)
        t = time_nb(snap.header)
        i, lag = _bpop_nearest_epoch(sevs, t)
        if i === nothing
            @warn "Binary population: snapshot $j (t = $t NB) has no sev.83 epoch at or before " *
                  "its time; stellar classes omitted"
            return nothing
        end
        ib, lag_b = _bpop_nearest_epoch(bevs, t)
        if ib === nothing && snap.npairs > 0
            @warn "Binary population: snapshot $j (t = $t NB) has KS pairs but no bev.82 epoch at " *
                  "or before its time; stellar classes omitted"
            return nothing
        end
        lagged = lag > 0 || (ib !== nothing && lag_b > 0)
        n_lagged += lagged
        lagged && (lag_max = max(lag_max, lag, ib === nothing ? 0.0 : lag_b))
        push!(maps, stellar_type_map(sevs[i], ib === nothing ? nothing : bevs[ib]))
    end
    n_lagged > 0 &&
        @info "Binary population: $n_lagged of $(length(snapshots)) snapshots take their stellar " *
              "types from the latest earlier sev.83/bev.82 epoch; largest lag $(lag_max) NB"
    return maps
end
