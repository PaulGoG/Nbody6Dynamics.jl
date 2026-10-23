# =============================================================================
# Gravitational potential, reduced systems and bound sets of snapshots
# =============================================================================
# The potential is the exact Newtonian pair sum (G = 1, no softening) in
# Float64. On the CPU it is a threaded loop over tiles of contiguous coordinate
# vectors, vectorised by @simd; the potential of a set on itself visits every
# pair once. The KernelAbstractions kernel below evaluates the same sum with
# one work-item per target, the form for a device backend; the tests check it
# against the CPU loop. The engine's own PHI is not used: it excludes the
# companion of a regularised component, is zero for ghosts and single precision
# (notes of 2026-10-14). Regularised pairs enter the energetics as their centre
# of mass: with the companion's potential included, a component's internal
# binding outweighs the energy of the pair's centre of mass, and an escaping
# hard pair would stay in the bound set.

"""Work-group size of the potential kernel (CPU backend: 256 matches a plain threaded loop)."""
const _POTENTIAL_WORKGROUP = 256
"""Targets per task of the one-sided CPU sum."""
const _POTENTIAL_TILE = 64
"""Sources per block of the one-sided CPU sum: four Float64 vectors of this length stay in the L2 cache."""
const _POTENTIAL_BLOCK = 2048
"""Bodies per tile of the symmetric CPU sum."""
const _POTENTIAL_SYM_TILE = 512

@kernel function _potential_kernel!(Φ, @Const(xt), @Const(xs), @Const(ms), ns)
    i = @index(Global)
    @inbounds begin
        xi = xt[1, i]
        yi = xt[2, i]
        zi = xt[3, i]
        acc = zero(eltype(Φ))
        for j in 1:ns
            dx = xs[1, j] - xi
            dy = xs[2, j] - yi
            dz = xs[3, j] - zi
            r2 = dx * dx + dy * dy + dz * dz
            acc += ifelse(r2 > 0, ms[j] / sqrt(r2), zero(acc))
        end
        Φ[i] = -acc
    end
end

"""
    pair_potential(pos, mass; backend = KernelAbstractions.CPU()) -> Vector{Float64}
    pair_potential(targets, sources, source_mass; backend = KernelAbstractions.CPU()) -> Vector{Float64}

Newtonian potential ``Φ_i = −Σ_j m_j / r_{ij}`` (`G = 1`, no softening). The
first form gives the potential of every body of `pos` (3 × n) from all the
others; the second the potential at the points `targets` (3 × n_t) of the
bodies at `sources` (3 × n_s) with masses `source_mass`, the form that updates
a potential by subtraction when bodies leave a set. Bodies at identical
positions do not interact, which also excludes the self-term. The sum is exact,
O(n_t n_s), in Float64. On the CPU it is a threaded loop over tiles,
vectorised by `@simd`; the first form visits each pair once and holds one
accumulator of n values per task (two tasks per thread), the second streams
blocks of sources through the cache for tiles of targets. The summation order
differs from a sequential sum, so results agree with it to round-off
(≈ 10⁻¹⁵ relative), not bit for bit.

# Example
```julia
Φ = pair_potential([0.0 1.0; 0.0 0.0; 0.0 0.0], [1.0, 2.0])   # [-2.0, -1.0]
```
"""
function pair_potential(
    pos::AbstractMatrix{<:Real},
    mass::AbstractVector{<:Real};
    backend::CPU = CPU(),
)::Vector{Float64}
    size(pos, 1) == 3 || throw(DimensionMismatch("targets must be 3 × n, got $(size(pos))"))
    size(pos, 2) == length(mass) ||
        throw(DimensionMismatch("$(size(pos, 2)) sources but $(length(mass)) source masses"))
    isempty(mass) && return Float64[]
    x, y, z = _coordinate_vectors(pos)
    return _pair_potential_symmetric(x, y, z, Vector{Float64}(mass))
end

function pair_potential(
    targets::AbstractMatrix{<:Real},
    sources::AbstractMatrix{<:Real},
    source_mass::AbstractVector{<:Real};
    backend::CPU = CPU(),
)::Vector{Float64}
    _check_potential_arguments(targets, sources, source_mass)
    nt = size(targets, 2)
    (nt == 0 || isempty(source_mass)) && return zeros(Float64, nt)
    xt, yt, zt = _coordinate_vectors(targets)
    xs, ys, zs = _coordinate_vectors(sources)
    return _pair_potential_tiled(xt, yt, zt, xs, ys, zs, Vector{Float64}(source_mass))
end

"""Throw `DimensionMismatch` unless `targets`, `sources` are 3 × n and the masses match the sources."""
function _check_potential_arguments(targets, sources, source_mass)::Nothing
    size(targets, 1) == 3 || throw(DimensionMismatch("targets must be 3 × n, got $(size(targets))"))
    size(sources, 1) == 3 || throw(DimensionMismatch("sources must be 3 × n, got $(size(sources))"))
    size(sources, 2) == length(source_mass) || throw(
        DimensionMismatch("$(size(sources, 2)) sources but $(length(source_mass)) source masses"),
    )
    return nothing
end

"""Rows of a 3 × n coordinate matrix as three contiguous Float64 vectors."""
function _coordinate_vectors(a::AbstractMatrix{<:Real})::NTuple{3,Vector{Float64}}
    return (Vector{Float64}(a[1, :]), Vector{Float64}(a[2, :]), Vector{Float64}(a[3, :]))
end

"""
    _pair_potential_tiled(xt, yt, zt, xs, ys, zs, ms) -> Vector{Float64}

One-sided CPU sum of [`pair_potential`](@ref): tasks of `_POTENTIAL_TILE`
targets, each sweeping the sources in blocks of `_POTENTIAL_BLOCK`, the inner
loop over a block vectorised.
"""
function _pair_potential_tiled(
    xt::Vector{Float64},
    yt::Vector{Float64},
    zt::Vector{Float64},
    xs::Vector{Float64},
    ys::Vector{Float64},
    zs::Vector{Float64},
    ms::Vector{Float64},
)::Vector{Float64}
    nt = length(xt)
    ns = length(ms)
    Φ = Vector{Float64}(undef, nt)
    Threads.@threads :dynamic for k in 1:cld(nt, _POTENTIAL_TILE)
        i0 = (k - 1) * _POTENTIAL_TILE + 1
        i1 = min(k * _POTENTIAL_TILE, nt)
        acc = zeros(Float64, i1 - i0 + 1)
        for j0 in 1:_POTENTIAL_BLOCK:ns
            j1 = min(j0 + _POTENTIAL_BLOCK - 1, ns)
            @inbounds for i in i0:i1
                xi = xt[i]
                yi = yt[i]
                zi = zt[i]
                a = 0.0
                @simd for j in j0:j1
                    dx = xs[j] - xi
                    dy = ys[j] - yi
                    dz = zs[j] - zi
                    r2 = muladd(dx, dx, muladd(dy, dy, dz * dz))
                    a += ifelse(r2 > 0.0, ms[j] / sqrt(r2), 0.0)
                end
                acc[i - i0 + 1] += a
            end
        end
        @inbounds for i in i0:i1
            Φ[i] = -acc[i - i0 + 1]
        end
    end
    return Φ
end

"""
    _pair_potential_symmetric(x, y, z, m) -> Vector{Float64}

Potential of a set on itself, visiting every pair once: the tile pairs
`I ≤ J` of `_POTENTIAL_SYM_TILE` bodies are dealt round-robin to two tasks per
thread, each accumulating both members' terms in its own vector of n values;
the vectors are summed at the end.
"""
function _pair_potential_symmetric(
    x::Vector{Float64},
    y::Vector{Float64},
    z::Vector{Float64},
    m::Vector{Float64},
)::Vector{Float64}
    n = length(m)
    nb = cld(n, _POTENTIAL_SYM_TILE)
    ntasks = min(2 * Threads.nthreads(), nb * (nb + 1) ÷ 2)
    parts = Vector{Vector{Float64}}(undef, ntasks)
    Threads.@threads :dynamic for c in 1:ntasks
        acc = zeros(Float64, n)
        q = 0
        for I in 1:nb, J in I:nb
            q += 1
            mod1(q, ntasks) == c && _symmetric_tile!(acc, x, y, z, m, I, J, n)
        end
        parts[c] = acc
    end
    Φ = zeros(Float64, n)
    for p in parts
        Φ .-= p
    end
    return Φ
end

"""Add the pair terms of tiles `I` and `J` (pairs `i < j` within a diagonal tile) to `acc`."""
function _symmetric_tile!(
    acc::Vector{Float64},
    x::Vector{Float64},
    y::Vector{Float64},
    z::Vector{Float64},
    m::Vector{Float64},
    I::Int,
    J::Int,
    n::Int,
)::Nothing
    i0 = (I - 1) * _POTENTIAL_SYM_TILE + 1
    i1 = min(I * _POTENTIAL_SYM_TILE, n)
    j0 = (J - 1) * _POTENTIAL_SYM_TILE + 1
    j1 = min(J * _POTENTIAL_SYM_TILE, n)
    @inbounds for i in i0:i1
        xi = x[i]
        yi = y[i]
        zi = z[i]
        mi = m[i]
        a = 0.0
        @simd for j in (I == J ? i + 1 : j0):j1
            dx = x[j] - xi
            dy = y[j] - yi
            dz = z[j] - zi
            r2 = muladd(dx, dx, muladd(dy, dy, dz * dz))
            w = ifelse(r2 > 0.0, 1.0 / sqrt(r2), 0.0)
            a += m[j] * w
            acc[j] += mi * w
        end
        acc[i] += a
    end
    return nothing
end

"""
    _pair_potential_kernel(backend, targets, sources, source_mass) -> Vector{Float64}

The one-sided sum of [`pair_potential`](@ref) evaluated by the
KernelAbstractions kernel `_potential_kernel!` on `backend` (one work-item per
target); the form a device backend runs.
"""
function _pair_potential_kernel(
    backend::CPU,
    targets::AbstractMatrix{<:Real},
    sources::AbstractMatrix{<:Real},
    source_mass::AbstractVector{<:Real},
)::Vector{Float64}
    _check_potential_arguments(targets, sources, source_mass)
    nt = size(targets, 2)
    ns = length(source_mass)
    (nt == 0 || ns == 0) && return zeros(Float64, nt)
    xt = _device_copy(backend, targets)
    xs = _device_copy(backend, sources)
    ms = _device_copy(backend, source_mass)
    Φ = KernelAbstractions.zeros(backend, Float64, nt)
    _potential_kernel!(backend, _POTENTIAL_WORKGROUP)(Φ, xt, xs, ms, ns; ndrange = nt)
    KernelAbstractions.synchronize(backend)
    return Array(Φ)
end

"""Float64 copy of `a` in the memory of `backend`."""
function _device_copy(backend, a::AbstractArray{<:Real})
    d = KernelAbstractions.allocate(backend, Float64, size(a))
    copyto!(d, Float64.(a))
    return d
end

# ---------------------------------------------------------------------------
# Reduced systems
# ---------------------------------------------------------------------------

"""
    SystemSet

Stars reduced to dynamical systems: single stars, and regularised pairs at
their centre of mass (mass-weighted position and velocity). All arrays are
Float64; `pos` and `vel` are 3 × n. `star1[s]` and `star2[s]` are the star
indices (in the arrays the set was built from) of system `s`, `star2[s] = 0`
for a single star. Pairs come first, in the order given, then the singles in
index order.
"""
struct SystemSet
    mass::Vector{Float64}
    pos::Matrix{Float64}
    vel::Matrix{Float64}
    star1::Vector{Int}
    star2::Vector{Int}
end

"""
    reduce_pairs(mass, pos, vel, pairs) -> SystemSet

Reduce the stars `mass`, `pos`, `vel` (3 × n) to a [`SystemSet`](@ref) in which
each `(i, j)` of `pairs` is one body at the pair's centre of mass. Throws
`ArgumentError` for a star index out of range or a star in two pairs.
"""
function reduce_pairs(
    mass::AbstractVector{<:Real},
    pos::AbstractMatrix{<:Real},
    vel::AbstractMatrix{<:Real},
    pairs::AbstractVector{<:NTuple{2,Integer}},
)::SystemSet
    n = length(mass)
    size(pos) == (3, n) || throw(DimensionMismatch("pos is $(size(pos)), expected (3, $n)"))
    size(vel) == (3, n) || throw(DimensionMismatch("vel is $(size(vel)), expected (3, $n)"))
    in_pair = falses(n)
    for (i, j) in pairs
        (1 ≤ i ≤ n && 1 ≤ j ≤ n && i != j) ||
            throw(ArgumentError("pair ($i, $j) out of range 1:$n"))
        (in_pair[i] || in_pair[j]) &&
            throw(ArgumentError("star $(in_pair[i] ? i : j) belongs to two pairs"))
        in_pair[i] = in_pair[j] = true
    end
    np = length(pairs)
    ns = np + (n - 2np)
    m = Vector{Float64}(undef, ns)
    x = Matrix{Float64}(undef, 3, ns)
    v = Matrix{Float64}(undef, 3, ns)
    s1 = Vector{Int}(undef, ns)
    s2 = zeros(Int, ns)
    for (s, (i, j)) in enumerate(pairs)
        mi, mj = Float64(mass[i]), Float64(mass[j])
        M = mi + mj
        m[s] = M
        for k in 1:3
            x[k, s] = (mi * pos[k, i] + mj * pos[k, j]) / M
            v[k, s] = (mi * vel[k, i] + mj * vel[k, j]) / M
        end
        s1[s] = i
        s2[s] = j
    end
    s = np
    for i in 1:n
        in_pair[i] && continue
        s += 1
        m[s] = mass[i]
        for k in 1:3
            x[k, s] = pos[k, i]
            v[k, s] = vel[k, i]
        end
        s1[s] = i
    end
    return SystemSet(m, x, v, s1, s2)
end

"""Sorted star indices of the systems `systems` of `sys`."""
function _system_stars(sys::SystemSet, systems::AbstractVector{Int})::Vector{Int}
    out = Int[]
    sizehint!(out, 2length(systems))
    for s in systems
        push!(out, sys.star1[s])
        sys.star2[s] > 0 && push!(out, sys.star2[s])
    end
    return sort!(out)
end

# ---------------------------------------------------------------------------
# Bound set
# ---------------------------------------------------------------------------

"""
    BoundSet

Result of the iterative energetic unbinding of [`bound_set`](@ref):

- `stars`: sorted star indices of the bound systems
- `systems`: indices of the bound systems in the [`SystemSet`](@ref)
- `passes`: energy evaluations performed
- `status`: `:converged` (a pass removed nothing), `:pass_limit` (the last of
  `max_passes` passes still removed systems; the set is not verified) or
  `:dissolved` (fewer than `min_members` systems remained)
- `v_ref`: mass-weighted mean velocity of the bound systems
- `kinetic`: ``T = ½ Σ M |V − v_ref|²`` of the bound systems
- `potential`: ``W = ½ Σ M Φ`` with Φ from the bound systems only
- `excluded`: stars left out before the iteration (components of pairs split
  by a subset selection)
"""
struct BoundSet
    stars::Vector{Int}
    systems::Vector{Int}
    passes::Int
    status::Symbol
    v_ref::NTuple{3,Float64}
    kinetic::Float64
    potential::Float64
    excluded::Vector{Int}
end

"""Default limit on unbinding passes (a safety bound; convergence takes a few passes)."""
const _BOUND_MAX_PASSES = 100
"""Default number of systems below which a set counts as dissolved."""
const _BOUND_MIN_MEMBERS = 10

"""
    bound_set(sys::SystemSet; potential = nothing, max_passes = 100, min_members = 10,
              backend = KernelAbstractions.CPU()) -> BoundSet
    bound_set(snap::Snapshot, idx = 1:nparticles(snap); kwargs...) -> BoundSet

Self-bound subset of the systems by iterative energetic unbinding: at every
pass the reference velocity is the mass-weighted mean velocity of the current
set and systems with ``½|V − v_\\mathrm{ref}|² + Φ ≥ 0`` are removed, Φ being the
potential of the current set; the potential of the survivors is then updated
by subtracting the removed systems' contribution, O(n·n_removed). The
iteration stops when a pass removes nothing (Springel et al. 2001,
doi:10.1046/j.1365-8711.2001.04912.x). `potential` may pass the potential of
every system from all the others when it is already known (it is not
modified).

The snapshot form reduces the regularised pairs of the snapshot among the
stars `idx` to their centres of mass first; a pair with only one component in
`idx` is left out and listed in `excluded`. `stars` are then snapshot indices.

Throws `ArgumentError` for `max_passes < 1` or `min_members < 2`. A pass limit
reached without convergence is reported by `status` and a warning.
"""
function bound_set(
    sys::SystemSet;
    potential::Union{Nothing,AbstractVector{<:Real}} = nothing,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
    backend = CPU(),
)::BoundSet
    max_passes ≥ 1 || throw(ArgumentError("max_passes must be ≥ 1, got $max_passes"))
    min_members ≥ 2 || throw(ArgumentError("min_members must be ≥ 2, got $min_members"))
    n = length(sys.mass)
    Φ =
        potential === nothing ? pair_potential(sys.pos, sys.mass; backend = backend) :
        Vector{Float64}(potential)
    length(Φ) == n || throw(DimensionMismatch("potential has $(length(Φ)) entries for $n systems"))
    set = collect(1:n)
    passes = 0
    status = :pass_limit
    while passes < max_passes
        if length(set) < min_members
            status = :dissolved
            break
        end
        passes += 1
        vref = _mean_velocity(sys, set)
        keep = BitVector(undef, length(set))
        @inbounds for (q, s) in enumerate(set)
            v2 =
                (sys.vel[1, s] - vref[1])^2 +
                (sys.vel[2, s] - vref[2])^2 +
                (sys.vel[3, s] - vref[3])^2
            keep[q] = 0.5 * v2 + Φ[s] < 0
        end
        if all(keep)
            status = :converged
            break
        end
        removed = set[.!keep]
        set = set[keep]
        isempty(set) && continue
        ΔΦ = pair_potential(
            sys.pos[:, set],
            sys.pos[:, removed],
            sys.mass[removed];
            backend = backend,
        )
        @inbounds for (q, s) in enumerate(set)
            Φ[s] -= ΔΦ[q]
        end
    end
    status === :pass_limit &&
        @warn "Bound set not converged after $passes passes ($(length(set)) systems left); raise [postprocess] bound_max_passes"
    vref = isempty(set) ? (NaN, NaN, NaN) : _mean_velocity(sys, set)
    T = 0.0
    W = 0.0
    @inbounds for s in set
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
    return BoundSet(_system_stars(sys, set), set, passes, status, vref, T, W, Int[])
end

function bound_set(
    snap::Snapshot,
    idx::AbstractVector{<:Integer} = 1:nparticles(snap);
    kwargs...,
)::BoundSet
    sub = collect(Int, idx)
    sys, excluded = _snapshot_systems(snap, sub)
    b = bound_set(sys; kwargs...)
    return BoundSet(
        sort!([sub[i] for i in b.stars]),
        b.systems,
        b.passes,
        b.status,
        b.v_ref,
        b.kinetic,
        b.potential,
        sort!([sub[i] for i in excluded]),
    )
end

"""
    bound_set(mass, pos, vel; pairs = NTuple{2,Int}[], kwargs...) -> BoundSet

Array form: the stars `mass`, `pos`, `vel` (3 × n) with the regularised `pairs`
reduced to their centres of mass ([`reduce_pairs`](@ref)); `stars` index the
arrays.
"""
function bound_set(
    mass::AbstractVector{<:Real},
    pos::AbstractMatrix{<:Real},
    vel::AbstractMatrix{<:Real};
    pairs::AbstractVector{<:NTuple{2,Integer}} = NTuple{2,Int}[],
    kwargs...,
)::BoundSet
    return bound_set(reduce_pairs(mass, pos, vel, pairs); kwargs...)
end

"""Mass-weighted mean velocity of the systems `set`."""
function _mean_velocity(sys::SystemSet, set::AbstractVector{Int})::NTuple{3,Float64}
    M = 0.0
    px = py = pz = 0.0
    @inbounds for s in set
        m = sys.mass[s]
        M += m
        px += m * sys.vel[1, s]
        py += m * sys.vel[2, s]
        pz += m * sys.vel[3, s]
    end
    return (px / M, py / M, pz / M)
end

"""
    _snapshot_systems(snap, sub) -> (SystemSet, excluded)

Systems of the snapshot stars `sub`: the regularised pairs with both
components in `sub` reduced to their centres of mass, the components of pairs
split by `sub` left out (`excluded`, positions in `sub`). The stars of a
hierarchical merger (negative `pair_name`) stay single: its stars 2k − 1 and 2k
are an inner component and the outer body, not a pair. Star indices of the
`SystemSet` are positions in `sub` without the excluded ones re-numbered, i.e.
`sub[star1[s]]` is a snapshot index.
"""
function _snapshot_systems(snap::Snapshot, sub::Vector{Int})
    n = nparticles(snap)
    local_of = zeros(Int, n)
    for (q, i) in enumerate(sub)
        (1 ≤ i ≤ n) || throw(ArgumentError("star index $i out of range 1:$n"))
        local_of[i] = q
    end
    pairs = NTuple{2,Int}[]
    drop = falses(length(sub))
    for (k, (a, b)) in enumerate(pair_indices(snap))
        !isempty(snap.pair_name) && snap.pair_name[k] < 0 && continue
        qa, qb = local_of[a], local_of[b]
        if qa > 0 && qb > 0
            push!(pairs, (qa, qb))
        elseif qa > 0 || qb > 0
            drop[max(qa, qb)] = true
        end
    end
    excluded = findall(drop)
    mass = Float64.(@view snap.mass[sub])
    pos = Float64.(@view snap.pos[:, sub])
    vel = Float64.(@view snap.vel[:, sub])
    if isempty(excluded)
        return reduce_pairs(mass, pos, vel, pairs), excluded
    end
    keep = findall(!, drop)
    renumber = zeros(Int, length(sub))
    renumber[keep] = 1:length(keep)
    sys = reduce_pairs(
        mass[keep],
        pos[:, keep],
        vel[:, keep],
        [(renumber[a], renumber[b]) for (a, b) in pairs],
    )
    # Map the star indices back to positions in `sub`.
    return SystemSet(
        sys.mass,
        sys.pos,
        sys.vel,
        keep[sys.star1],
        [s == 0 ? 0 : keep[s] for s in sys.star2],
    ),
    excluded
end
