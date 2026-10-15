# =============================================================================
# Angular-momentum budget and resolved rotation of a snapshot series
# =============================================================================
# Closes the angular momentum of the stars against escapers, newcomers and stellar mass loss,
# and resolves the rotation of the bound set by radius and by initial-mass class.

using LinearAlgebra: dot, norm

"""
    AngularMomentumBudget

Angular-momentum bookkeeping of a snapshot series (N-body units unless noted), one column
per snapshot: `time_nb`, `time_myr`; the elapsed relaxation times `n_rh` (2 × n_t,
∫dt/t_rh of the bound set for the two Coulomb factors); the internal angular momentum of
all stars, about their centre of mass and mean velocity, `l_total` (3 × n_t) and the
cumulative `l_removed`, `l_added`, `l_stellar`, `l_remnant` (3 × n_t: the internal angular
momentum carried off by the stars that left the snapshots, brought in by those that
appeared, shed with the mass lost by evolving stars, and changed by the stars that became
remnants in the interval, kick and mass loss together); the budget `residual` (3 × n_t)
and the scale `l_scale` (n_t, `Σ m |r − c| |v − v̄|` of all stars about their centre of mass
and mean velocity, the angular-momentum content a residual is compared with);
the bound set's angular momentum about its centre of mass `l_bound` (3 × n_t), its
`axis` (3 × n_t, unit vector), the Lagrangian radii `r_lagr` (n_f × n_t) of the bound set
with `fractions`, `lambda_r_radius` and `lambda_r_radius_floor` (n_f × n_t) inside them;
the initial-mass classes `class_edges` [M☉] with, per class (n_c × n_t), the members in
the bound set `n_class`, the angular momentum along the axis `l_z_class`, its fraction
`l_z_class_fraction` of the bound set's, and `lambda_r_class`, `lambda_r_class_floor`.
"""
struct AngularMomentumBudget
    time_nb::Vector{Float64}
    time_myr::Vector{Float64}
    n_rh::Matrix{Float64}
    l_total::Matrix{Float64}
    l_removed::Matrix{Float64}
    l_added::Matrix{Float64}
    l_stellar::Matrix{Float64}
    l_remnant::Matrix{Float64}
    residual::Matrix{Float64}
    l_scale::Vector{Float64}
    l_bound::Matrix{Float64}
    axis::Matrix{Float64}
    fractions::Vector{Float64}
    r_lagr::Matrix{Float64}
    lambda_r_radius::Matrix{Float64}
    lambda_r_radius_floor::Matrix{Float64}
    class_edges::Vector{Float64}
    n_class::Matrix{Int}
    l_z_class::Matrix{Float64}
    l_z_class_fraction::Matrix{Float64}
    lambda_r_class::Matrix{Float64}
    lambda_r_class_floor::Matrix{Float64}
end

"""
    _internal_angular_momentum(pos, vel, idx, m) -> Vector{Float64}

Angular momentum of the columns `idx` of `pos`, `vel` with the masses `m` (aligned with
`idx`) about their own centre of mass and mass-weighted mean velocity: the internal
angular momentum, invariant under translations and boosts of the frame. Zero for an
empty set or a vanishing total mass.
"""
function _internal_angular_momentum(
    pos::Matrix{Float64},
    vel::Matrix{Float64},
    idx::AbstractVector{Int},
    m::AbstractVector{Float64},
)::Vector{Float64}
    length(m) == length(idx) || throw(DimensionMismatch("masses and columns disagree in length"))
    M = 0.0
    cx = cy = cz = 0.0
    vx = vy = vz = 0.0
    @inbounds for (s, i) in enumerate(idx)
        M += m[s]
        cx += m[s] * pos[1, i]
        cy += m[s] * pos[2, i]
        cz += m[s] * pos[3, i]
        vx += m[s] * vel[1, i]
        vy += m[s] * vel[2, i]
        vz += m[s] * vel[3, i]
    end
    M > 0 || return zeros(3)
    cx /= M
    cy /= M
    cz /= M
    vx /= M
    vy /= M
    vz /= M
    lx = ly = lz = 0.0
    @inbounds for (s, i) in enumerate(idx)
        rx, ry, rz = pos[1, i] - cx, pos[2, i] - cy, pos[3, i] - cz
        ux, uy, uz = vel[1, i] - vx, vel[2, i] - vy, vel[3, i] - vz
        lx += m[s] * (ry * uz - rz * uy)
        ly += m[s] * (rz * ux - rx * uz)
        lz += m[s] * (rx * uy - ry * ux)
    end
    return [lx, ly, lz]
end

"""
    _angular_momentum_scale(pos, vel, mass) -> Float64

`Σ mᵢ |rᵢ − c| |vᵢ − v̄|` about the centre of mass `c` and the mass-weighted mean velocity
`v̄`: the angular-momentum content of the set, an upper bound of `|L|` that is reached only
when every star moves tangentially in one sense.
"""
function _angular_momentum_scale(
    pos::Matrix{Float64},
    vel::Matrix{Float64},
    mass::Vector{Float64},
)::Float64
    M = sum(mass)
    M > 0 || return 0.0
    c = vec(sum(pos .* mass'; dims = 2)) ./ M
    v̄ = vec(sum(vel .* mass'; dims = 2)) ./ M
    s = 0.0
    @inbounds for i in eachindex(mass)
        r = sqrt((pos[1, i] - c[1])^2 + (pos[2, i] - c[2])^2 + (pos[3, i] - c[3])^2)
        u = sqrt((vel[1, i] - v̄[1])^2 + (vel[2, i] - v̄[2])^2 + (vel[3, i] - v̄[3])^2)
        s += mass[i] * r * u
    end
    return s
end

"""
    angular_momentum_budget(snaps, bound_sets, units; events = EngineEvent[],
                            fractions = (0.1, 0.5, 0.9), class_edges_solar = (0.5, 5.0),
                            settings = PostprocessConfig(), n_min_class = 100)
        -> AngularMomentumBudget

Angular-momentum budget of a snapshot series and the rotation of its bound sets, in
N-body units (see [`AngularMomentumBudget`](@ref)). Stars are identified across
snapshots by `name`.

The budget follows the internal angular momentum of all stars, about their centre of mass
and mass-weighted mean velocity (`l_total`), between consecutive snapshots `k − 1` and `k`.
The internal angular momentum is invariant under translations and boosts of the frame, so
the engine's centre-of-mass corrections after escaper removals do not enter it. The stars
absent from `k` carry off the difference between the internal angular momentum of
snapshot `k − 1` with and without them (`l_removed`); the stars new in `k` bring the same
difference at `k` (`l_added`); the stars common to both snapshots shed, with their mass
loss, the difference between their internal angular momentum with the earlier and with the
later masses, evaluated at either snapshot's phase space and averaged (`l_stellar`). The
stars that became remnants in the interval, named by the engine's `remnant_formation` and
`kick_*` events of `events` (`name1`, event time in `(t_{k−1}, t_k]`), are left out of the
shed term; the change of the internal angular momentum they cause, the difference between
the common set with and without them at `k` minus the same at `k − 1`, is their own term
`l_remnant` (the natal kick and the mass loss of the remnant together). With these
cumulative terms

    residual = l_total + l_removed − l_added + l_stellar − l_remnant − l_total(t₁),

which vanishes exactly for stars removed, masses changed or remnants kicked at fixed
phase space of the others, and otherwise measures the torque of the departing stars on
the rest between the last snapshot that holds them and their removal, together with the
integration error (isolated run) or the tidal torque (external field). A configuration
with little net rotation has `|l_total|` far below its angular-momentum content, so the
residual is judged against `l_scale`, not against `|l_total|`.

The bound set of snapshot `k` (`bound_sets[k].stars`) defines the spin `axis`, the unit
vector of its angular momentum `l_bound` about its centre of mass and mean velocity
(`(0, 0, 1)` when it vanishes). Its shrinking-sphere density centre (centre keys of
`settings`) gives the Lagrangian radii `r_lagr` for `fractions`, and its half-mass
relaxation times for the two `settings.relaxation_gamma` factors give the elapsed
relaxation times `n_rh = ∫dt/t_rh` (trapezoid, zero at the first snapshot; `NaN` where
`γ N ≤ 1`). The ordered-to-total motion parameter λ_R and its sign-flip floor
([`rotation_analysis`](@ref), about `axis`) are measured inside every Lagrangian radius
and for every initial-mass class, the class of a star being fixed by its mass in M☉ in
the first snapshot where it appears against `class_edges_solar`; a subset with fewer
than `max(10, n_min_class)` stars gets `NaN`. The axial angular momentum of each class
about the bound set's centre of mass, `l_z_class`, and its fraction of the bound set's
complete the class decomposition. An empty bound set leaves all bound-set quantities of
its snapshot `NaN`.

# Arguments
- `snaps`: snapshot series in time order
- `bound_sets`: one [`BoundSet`](@ref) per snapshot
- `units`: unit scaling; `zmbar` converts masses to M☉ for the classes
- `events`: the engine's stamped events ([`read_engine_events`](@ref)); the remnant
  formations and kicks name the stars of the `l_remnant` term
- `fractions`: enclosed mass fractions, strictly increasing in `(0, 1]`
- `class_edges_solar`: class edges [M☉], positive and strictly increasing
- `settings`: [`PostprocessConfig`](@ref) (centre, rotation and relaxation keys)
- `n_min_class`: shell membership of the class profiles and minimum subset size, `≥ 2`

# Returns
An [`AngularMomentumBudget`](@ref).

Throws `ArgumentError` without snapshots, for a bound-set count other than the snapshot
count, or for invalid `fractions`, `class_edges_solar` or `n_min_class`.

# Example
```julia
bsets = [bound_set(s) for s in snaps]
b = angular_momentum_budget(snaps, bsets, units)
maximum(abs, b.residual) / maximum(abs, b.l_total[:, 1])
```
"""
function angular_momentum_budget(
    snaps::Vector{Snapshot},
    bound_sets::AbstractVector{BoundSet},
    units::UnitScaling;
    events::AbstractVector{EngineEvent} = EngineEvent[],
    fractions = (0.1, 0.5, 0.9),
    class_edges_solar = (0.5, 5.0),
    settings::PostprocessConfig = PostprocessConfig(),
    n_min_class::Integer = 100,
)::AngularMomentumBudget
    n_t = length(snaps)
    n_t ≥ 1 || throw(ArgumentError("angular_momentum_budget needs at least one snapshot"))
    length(bound_sets) == n_t || throw(
        ArgumentError(
            "angular_momentum_budget: $(length(bound_sets)) bound sets for $n_t snapshots",
        ),
    )
    f = collect(Float64, fractions)
    (all(x -> 0 < x ≤ 1, f) && all(>(0), diff(f))) || throw(
        ArgumentError(
            "angular_momentum_budget: fractions must be strictly increasing in (0, 1]; got $f",
        ),
    )
    edges = collect(Float64, class_edges_solar)
    (all(>(0), edges) && all(>(0), diff(edges))) || throw(
        ArgumentError(
            "angular_momentum_budget: class_edges_solar must be positive and strictly increasing; got $edges",
        ),
    )
    n_min_class ≥ 2 ||
        throw(ArgumentError("angular_momentum_budget: n_min_class must be ≥ 2; got $n_min_class"))

    n_f = length(f)
    n_c = length(edges) + 1
    n_req = max(_MIN_MEMBERS, Int(n_min_class))
    γ = (settings.relaxation_gamma[1], settings.relaxation_gamma[2])
    j_half = findfirst(==(0.5), f)
    centre_kw = _centre_keywords(settings)

    t_nb = zeros(n_t)
    t_myr = zeros(n_t)
    n_rh = fill(NaN, 2, n_t)
    l_total = zeros(3, n_t)
    l_removed = zeros(3, n_t)
    l_added = zeros(3, n_t)
    l_stellar = zeros(3, n_t)
    l_remnant = zeros(3, n_t)
    residual = zeros(3, n_t)
    l_scale = zeros(n_t)
    l_bound = fill(NaN, 3, n_t)
    axis = fill(NaN, 3, n_t)
    r_lagr = fill(NaN, n_f, n_t)
    λr_radius = fill(NaN, n_f, n_t)
    λr_radius_floor = fill(NaN, n_f, n_t)
    n_class = zeros(Int, n_c, n_t)
    l_z_class = fill(NaN, n_c, n_t)
    l_z_class_fraction = fill(NaN, n_c, n_t)
    λr_class = fill(NaN, n_c, n_t)
    λr_class_floor = fill(NaN, n_c, n_t)

    # Initial-mass class of every name, from its first appearance
    class_of = Dict{Int32,Int}()
    for snap in snaps
        for i in eachindex(snap.name)
            nm = snap.name[i]
            haskey(class_of, nm) && continue
            m₀ = Float64(snap.mass[i]) * units.zmbar
            class_of[nm] = 1 + count(e -> m₀ ≥ e, edges)
        end
    end

    remnant_kinds = (:remnant_formation, :kick_full, :kick_zero, :kick_low)
    remnant_events = [e for e in events if e.kind in remnant_kinds && e.name1 != 0]
    pos_prev = zeros(3, 0)
    vel_prev = zeros(3, 0)
    mass_prev = Float64[]
    t_rh_prev = (NaN, NaN)
    Lc = zeros(3, n_c)
    for (k, snap) in enumerate(snaps)
        pos = Matrix{Float64}(snap.pos)
        vel = Matrix{Float64}(snap.vel)
        mass = Vector{Float64}(snap.mass)
        N = length(mass)
        slot = Dict{Int32,Int}()
        sizehint!(slot, N)
        @inbounds for i in 1:N
            slot[snap.name[i]] = i
        end
        t_nb[k] = time_nb(snap.header)
        t_myr[k] = time_myr(snap.header)

        l_total[:, k] .= _internal_angular_momentum(pos, vel, 1:N, mass)
        l_scale[k] = _angular_momentum_scale(pos, vel, mass)

        if k ≥ 2
            # The stars common to k − 1 and k: columns `kept` of k − 1 and `old` of k, in one order
            name_prev = snaps[k - 1].name
            n_prev = length(mass_prev)
            kept = Int[]
            old = Int[]
            sizehint!(kept, n_prev)
            sizehint!(old, n_prev)
            @inbounds for i in 1:n_prev
                q = get(slot, name_prev[i], 0)
                q == 0 && continue
                push!(kept, i)
                push!(old, q)
            end
            m_prev_c = mass_prev[kept]
            m_now_c = mass[old]
            # Departed stars: the internal angular momentum of k − 1 with and without them
            L_prev_kept = _internal_angular_momentum(pos_prev, vel_prev, kept, m_prev_c)
            removed = l_total[:, k - 1] .- L_prev_kept
            # Stars new at k: the same at k
            L_now_old_now = _internal_angular_momentum(pos, vel, old, m_now_c)
            added = l_total[:, k] .- L_now_old_now
            # Stars that became remnants in the interval: left out of the shed term, their own
            # change of the common set's internal angular momentum between the snapshots
            t_lo, t_hi = t_nb[k - 1], t_nb[k]
            transformed =
                Set{Int32}(Int32(e.name1) for e in remnant_events if t_lo < e.time_nb ≤ t_hi)
            quiet = [q for q in eachindex(kept) if !(name_prev[kept[q]] in transformed)]
            kept_q = kept[quiet]
            old_q = old[quiet]
            m_prev_q = m_prev_c[quiet]
            m_now_q = m_now_c[quiet]
            L_prev_q = _internal_angular_momentum(pos_prev, vel_prev, kept_q, m_prev_q)
            L_now_q = _internal_angular_momentum(pos, vel, old_q, m_now_q)
            remnant = (L_now_old_now .- L_now_q) .- (L_prev_kept .- L_prev_q)
            # Mass shed by the quiet common stars: their internal angular momentum with the
            # earlier masses minus that with the later ones, at either phase space, averaged
            shed_prev = L_prev_q .- _internal_angular_momentum(pos_prev, vel_prev, kept_q, m_now_q)
            shed_now = _internal_angular_momentum(pos, vel, old_q, m_prev_q) .- L_now_q
            l_removed[:, k] .= l_removed[:, k - 1] .+ removed
            l_added[:, k] .= l_added[:, k - 1] .+ added
            l_stellar[:, k] .= l_stellar[:, k - 1] .+ 0.5 .* (shed_prev .+ shed_now)
            l_remnant[:, k] .= l_remnant[:, k - 1] .+ remnant
        end

        b = bound_sets[k].stars
        t_rh = (NaN, NaN)
        if !isempty(b)
            n_b = length(b)
            pos_b = pos[:, b]
            vel_b = vel[:, b]
            mass_b = mass[b]

            # Centre of mass and mean velocity of the bound set
            M_b = 0.0
            cx = cy = cz = 0.0
            vx = vy = vz = 0.0
            @inbounds for i in 1:n_b
                m = mass_b[i]
                M_b += m
                cx += m * pos_b[1, i]
                cy += m * pos_b[2, i]
                cz += m * pos_b[3, i]
                vx += m * vel_b[1, i]
                vy += m * vel_b[2, i]
                vz += m * vel_b[3, i]
            end
            cx /= M_b
            cy /= M_b
            cz /= M_b
            vx /= M_b
            vy /= M_b
            vz /= M_b

            # Angular momentum about (c_b, v_b), in total and per class
            fill!(Lc, 0.0)
            cls = Vector{Int}(undef, n_b)
            hx_b = hy_b = hz_b = 0.0
            @inbounds for i in 1:n_b
                c = class_of[snap.name[b[i]]]
                cls[i] = c
                m = mass_b[i]
                px, py, pz = pos_b[1, i] - cx, pos_b[2, i] - cy, pos_b[3, i] - cz
                ux, uy, uz = vel_b[1, i] - vx, vel_b[2, i] - vy, vel_b[3, i] - vz
                hx = m * (py * uz - pz * uy)
                hy = m * (pz * ux - px * uz)
                hz = m * (px * uy - py * ux)
                hx_b += hx
                hy_b += hy
                hz_b += hz
                Lc[1, c] += hx
                Lc[2, c] += hy
                Lc[3, c] += hz
            end
            lb = [hx_b, hy_b, hz_b]
            l_bound[:, k] .= lb
            lb_norm = norm(lb)
            ax = lb_norm == 0 ? [0.0, 0.0, 1.0] : lb ./ lb_norm
            axis[:, k] .= ax

            centre = shrinking_sphere_centre(pos_b, mass_b; centre_kw...).centre
            r_lagr[:, k] .= _lagrangian_radii(pos_b, mass_b, centre, Tuple(f))

            # relaxation_times requires γ N > 1; smaller sets leave the clock NaN
            if !any(g -> g > 0 && g * n_b ≤ 1, γ)
                rt = relaxation_times(pos_b, vel_b, mass_b, centre; gamma = γ)
                t_rh = (Float64(rt.t_rh[1]), Float64(rt.t_rh[2]))
            end

            r_h = j_half === nothing ? NaN : r_lagr[j_half, k]

            # Rotation inside the Lagrangian radii
            d = Vector{Float64}(undef, n_b)
            @inbounds for i in 1:n_b
                d[i] = sqrt(
                    (pos_b[1, i] - centre[1])^2 +
                    (pos_b[2, i] - centre[2])^2 +
                    (pos_b[3, i] - centre[3])^2,
                )
            end
            for j in 1:n_f
                sub = findall(≤(r_lagr[j, k]), d)
                length(sub) ≥ n_req || continue
                rot = rotation_analysis(
                    pos_b[:, sub],
                    vel_b[:, sub],
                    mass_b[sub];
                    spin_axis = ax,
                    n_min = settings.rotation_n_min,
                    nbins_max = settings.rotation_nbins_max,
                    n_null = settings.rotation_n_null,
                    seed = settings.msr_seed,
                    r_h = r_h,
                )
                λr_radius[j, k] = rot.lambda_r
                λr_radius_floor[j, k] = rot.lambda_r_floor
            end

            # Rotation per initial-mass class
            lz_b = dot(ax, lb)
            for c in 1:n_c
                members = findall(==(c), cls)
                n_class[c, k] = length(members)
                lz_c = ax[1] * Lc[1, c] + ax[2] * Lc[2, c] + ax[3] * Lc[3, c]
                l_z_class[c, k] = lz_c
                l_z_class_fraction[c, k] = lz_b == 0 ? NaN : lz_c / lz_b
                length(members) ≥ n_req || continue
                rot = rotation_analysis(
                    pos_b[:, members],
                    vel_b[:, members],
                    mass_b[members];
                    spin_axis = ax,
                    n_min = Int(n_min_class),
                    nbins_max = settings.rotation_nbins_max,
                    n_null = settings.rotation_n_null,
                    seed = settings.msr_seed,
                    r_h = r_h,
                )
                λr_class[c, k] = rot.lambda_r
                λr_class_floor[c, k] = rot.lambda_r_floor
            end
        end

        # Elapsed relaxation times, trapezoid on 1/t_rh
        for g in 1:2
            n_rh[g, k] =
                k == 1 ? 0.0 :
                n_rh[g, k - 1] + (t_nb[k] - t_nb[k - 1]) * 0.5 * (1 / t_rh_prev[g] + 1 / t_rh[g])
        end

        pos_prev = pos
        vel_prev = vel
        mass_prev = mass
        t_rh_prev = t_rh
    end

    @inbounds for k in 1:n_t, d in 1:3
        residual[d, k] =
            l_total[d, k] + l_removed[d, k] - l_added[d, k] + l_stellar[d, k] - l_remnant[d, k] -
            l_total[d, 1]
    end

    return AngularMomentumBudget(
        t_nb,
        t_myr,
        n_rh,
        l_total,
        l_removed,
        l_added,
        l_stellar,
        l_remnant,
        residual,
        l_scale,
        l_bound,
        axis,
        f,
        r_lagr,
        λr_radius,
        λr_radius_floor,
        edges,
        n_class,
        l_z_class,
        l_z_class_fraction,
        λr_class,
        λr_class_floor,
    )
end

"""
    write_angular_momentum(path, b::AngularMomentumBudget) -> String

Write `b` as CSV, one line per snapshot (N-body units plus `time_myr`). A leading
comment line records `fractions` and `class_edges_solar`; the columns are the times,
`n_rh_1`, `n_rh_2`, the components of `l_total` (`lx`, `ly`, `lz`), `l_removed`,
`l_added`, `l_stellar`, `l_remnant`, `residual`, `l_scale`, `l_bound` and `axis`, then `r_lagr_j`, `lambda_r_j`,
`lambda_r_floor_j` per fraction `j` and `n_class_c`, `l_z_class_c`, `l_z_fraction_c`,
`lambda_r_class_c`, `lambda_r_class_floor_c` per class `c`. Floats carry 12 significant
digits, non-finite values are empty cells. An existing file is backed up, never
overwritten. Returns `path`.

# Example
```julia
write_angular_momentum(joinpath(out_dir, "angular_momentum.csv"), b)
```
"""
function write_angular_momentum(path::AbstractString, b::AngularMomentumBudget)::String
    _backup_existing(path)
    n_f = length(b.fractions)
    n_c = size(b.n_class, 1)
    header = String[
        "time_nb",
        "time_myr",
        "n_rh_1",
        "n_rh_2",
        "lx",
        "ly",
        "lz",
        "l_removed_x",
        "l_removed_y",
        "l_removed_z",
        "l_added_x",
        "l_added_y",
        "l_added_z",
        "l_stellar_x",
        "l_stellar_y",
        "l_stellar_z",
        "l_remnant_x",
        "l_remnant_y",
        "l_remnant_z",
        "residual_x",
        "residual_y",
        "residual_z",
        "l_scale",
        "l_bound_x",
        "l_bound_y",
        "l_bound_z",
        "axis_x",
        "axis_y",
        "axis_z",
    ]
    for j in 1:n_f
        push!(header, "r_lagr_$j", "lambda_r_$j", "lambda_r_floor_$j")
    end
    for c in 1:n_c
        push!(
            header,
            "n_class_$c",
            "l_z_class_$c",
            "l_z_fraction_$c",
            "lambda_r_class_$c",
            "lambda_r_class_floor_$c",
        )
    end
    open(path, "w") do io
        println(
            io,
            "# fractions=",
            join(_csv_precise.(b.fractions), ";"),
            " class_edges_solar=",
            join(_csv_precise.(b.class_edges), ";"),
        )
        println(io, join(header, ","))
        row = Vector{String}()
        sizehint!(row, length(header))
        for k in eachindex(b.time_nb)
            empty!(row)
            push!(row, _csv_precise(b.time_nb[k]), _csv_precise(b.time_myr[k]))
            push!(row, _csv_precise(b.n_rh[1, k]), _csv_precise(b.n_rh[2, k]))
            for m in (b.l_total, b.l_removed, b.l_added, b.l_stellar, b.l_remnant, b.residual)
                for d in 1:3
                    push!(row, _csv_precise(m[d, k]))
                end
            end
            push!(row, _csv_precise(b.l_scale[k]))
            for m in (b.l_bound, b.axis)
                for d in 1:3
                    push!(row, _csv_precise(m[d, k]))
                end
            end
            for j in 1:n_f
                push!(
                    row,
                    _csv_precise(b.r_lagr[j, k]),
                    _csv_precise(b.lambda_r_radius[j, k]),
                    _csv_precise(b.lambda_r_radius_floor[j, k]),
                )
            end
            for c in 1:n_c
                push!(
                    row,
                    string(b.n_class[c, k]),
                    _csv_precise(b.l_z_class[c, k]),
                    _csv_precise(b.l_z_class_fraction[c, k]),
                    _csv_precise(b.lambda_r_class[c, k]),
                    _csv_precise(b.lambda_r_class_floor[c, k]),
                )
            end
            println(io, join(row, ","))
        end
    end
    return String(path)
end
