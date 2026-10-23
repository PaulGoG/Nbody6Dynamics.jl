# =============================================================================
# Merger tree and per-progenitor mass budget of multi-cluster runs
# =============================================================================

# ---------------------------------------------------------------------------
# Merger tree
# ---------------------------------------------------------------------------

"""
    MergerEvent

One join of the merger tree: the snapshot `index` from which two groups of progenitors stay
within the coalescence threshold, its `time_nb` and `time_myr`, the tags of the two groups
`group_a` and `group_b` (sorted, `group_a` holds the smallest tag of the two), the
density-centre `separation` and `threshold` of the pair whose join connects the two groups at
that snapshot.
"""
struct MergerEvent
    index::Int
    time_nb::Float64
    time_myr::Float64
    group_a::Vector{Int}
    group_b::Vector{Int}
    separation::Float64
    threshold::Float64
end

"""
    MergerTree

Merger tree of the initial clusters: the `events` in time order, `n_distinct` groups per
snapshot, the `pairs`, their `separation` and `threshold` matrices (`n_pairs × n_t`), the
`present` matrix (`n_cl × n_t`), the snapshot index `joined` of every pair (`n_pairs`; 0 for a
pair that never joins) and the dissolution index `dissolved` of every cluster (`n_cl`; the
first snapshot after which it is never present again, 0 when present at the last snapshot or
never present).
"""
struct MergerTree
    events::Vector{MergerEvent}
    n_distinct::Vector{Int}
    pairs::Vector{Tuple{Int,Int}}
    separation::Matrix{Float64}
    threshold::Matrix{Float64}
    present::BitMatrix
    joined::Vector{Int}
    dissolved::Vector{Int}
end

"""
    merger_tree(snaps, cluster_ranges; bound_sets = nothing, max_passes = 100, min_members = 10)
        -> MergerTree

Merger tree of the initial clusters, derived from one call of [`coalescence_time`](@ref)
with the same arguments.

The pair `q = (i, j)` is joined from snapshot `k_q`, the smallest `k` at which both clusters
are present, followed by at least one later snapshot with both present (or `k` is the last
snapshot), such that `separation[q, k′] < threshold[q, k′]` at every `k′ ≥ k` where both are
present. A first pericentre passage that separates again therefore does not count. Pairs
that never satisfy this get `joined[q] = 0`.

The joined pairs are processed in increasing `k_q`, ties in pair order, with union–find over
the tags; every pair whose two tags still belong to different groups records a
[`MergerEvent`](@ref) at `k_q` with the separation and threshold of that pair, and the groups
are then united. `dissolved[i]` is the snapshot after the last one at which cluster `i` is
present, or 0 when it is present at the last snapshot or never present.

# Arguments
- `snaps`: snapshots in time order
- `cluster_ranges`: names of every initial cluster (ranges or vectors of names)
- `bound_sets`, `max_passes`, `min_members`: passed to `coalescence_time`

Throws `ArgumentError` when `snaps` is empty or fewer than two clusters are given.

# Example
```julia
tree = merger_tree(snaps, cluster_ranges; bound_sets)
[(e.index, e.group_a, e.group_b) for e in tree.events]
```
"""
function merger_tree(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}};
    bound_sets::Union{Nothing,AbstractVector{BoundSet}} = nothing,
    max_passes::Integer = _BOUND_MAX_PASSES,
    min_members::Integer = _BOUND_MIN_MEMBERS,
)::MergerTree
    co = coalescence_time(snaps, cluster_ranges; bound_sets, max_passes, min_members)
    separation = co.separation
    threshold = co.threshold
    pairs = co.pairs
    present = co.present
    n_cl, n_t = size(present)
    n_pairs = length(pairs)

    # Persistence per pair: scan back from the last snapshot until the threshold is violated
    joined = zeros(Int, n_pairs)
    for (q, (i, j)) in enumerate(pairs)
        later = false
        @inbounds for k in n_t:-1:1
            (present[i, k] && present[j, k]) || continue
            separation[q, k] < threshold[q, k] || break
            (later || k == n_t) && (joined[q] = k)
            later = true
        end
    end

    order = sort!(findall(>(0), joined); by = q -> (joined[q], q))
    parent = collect(1:n_cl)
    events = MergerEvent[]
    for q in order
        i, j = pairs[q]
        ri = _union_root!(parent, i)
        rj = _union_root!(parent, j)
        ri == rj && continue
        group_a = [t for t in 1:n_cl if _union_root!(parent, t) == ri]
        group_b = [t for t in 1:n_cl if _union_root!(parent, t) == rj]
        first(group_a) < first(group_b) || ((group_a, group_b) = (group_b, group_a))
        k = joined[q]
        h = snaps[k].header
        push!(
            events,
            MergerEvent(
                k,
                time_nb(h),
                time_myr(h),
                group_a,
                group_b,
                separation[q, k],
                threshold[q, k],
            ),
        )
        parent[rj] = ri
    end

    dissolved = zeros(Int, n_cl)
    for i in 1:n_cl
        k_last = findlast(view(present, i, :))
        (k_last === nothing || k_last == n_t) || (dissolved[i] = k_last + 1)
    end

    return MergerTree(
        events,
        co.n_distinct,
        pairs,
        separation,
        threshold,
        present,
        joined,
        dissolved,
    )
end

# ---------------------------------------------------------------------------
# Per-progenitor mass budget and mixing
# ---------------------------------------------------------------------------

"""
    ProgenitorBudgetRow

One progenitor at one snapshot, N-body mass units: `time_nb`, `time_myr`, `tag`; `n_initial`
and `m_initial` of the names of the tag in the first snapshot; `n_present`, `m_present` of
those names in this snapshot; `n_bound`, `m_bound` of them inside the global bound set;
`m_removed` of the names absent from this snapshot, at their last mass before leaving;
`m_stellar_loss` = Σ (m_initial − m_now) over the present names plus Σ (m_initial − m_last)
over the absent names, their loss before departure, so that
`m_initial = m_present + m_removed + m_stellar_loss` holds exactly for the initial names;
`r_h` of the tag's bound members about their own density centre; `centre_offset`, the
distance of that centre from the density centre of the whole bound set; and the mixing of the
tag's bound members against the other bound members about the global centre within 2 r_h of
the whole bound set: `aplus2`, `aplus2_err`, `d_radius`, `d_energy`, `d_angular_momentum`,
`d_floor`.
"""
struct ProgenitorBudgetRow
    time_nb::Float64
    time_myr::Float64
    tag::Int
    n_initial::Int
    m_initial::Float64
    n_present::Int
    m_present::Float64
    n_bound::Int
    m_bound::Float64
    m_removed::Float64
    m_stellar_loss::Float64
    r_h::Float64
    centre_offset::Float64
    aplus2::Float64
    aplus2_err::Float64
    d_radius::Float64
    d_energy::Float64
    d_angular_momentum::Float64
    d_floor::Float64
end

"""
    progenitor_budget(snaps, cluster_ranges, bound_sets; settings = PostprocessConfig(),
                      min_members = 20) -> Vector{ProgenitorBudgetRow}

Mass budget of every initial cluster (tag of `cluster_ranges`) at every snapshot, with the
structure and radial mixing of its bound members in the product. One row per
(snapshot, tag), snapshots outer, tags inner.

The names of a tag are those of the first snapshot lying in its range; mass columns follow
these names only. A name absent from a snapshot counts in `m_removed` at its last known mass,
and its loss before departure in `m_stellar_loss`.

The global bound set is `bound_sets[k].stars`. When it holds at least `min_members` stars,
its shrinking-sphere density centre `c_g` (keywords from `settings`) and half-mass radius
`r_h,g` about it are computed once per snapshot. A tag with at least `min_members` bound
members then gets its own density centre, its half-mass radius `r_h` about that centre, the
offset `|c_i − c_g|`, and [`mixing`](@ref) of its bound members (population 1) against the
other bound stars (population 2) about `c_g` within `2 r_h,g`; the specific energies
`½|v − v_ref|² + Φ` enter when the snapshot carries per-star potentials. `r_h`,
`centre_offset` and the mixing fields are `NaN` otherwise, and the mixing fields are also
`NaN` when the tag holds every bound star.

# Arguments
- `snaps`: snapshots in time order
- `cluster_ranges`: names of every initial cluster
- `bound_sets`: global bound set of every snapshot
- `settings`: postprocessing settings (shrinking-sphere keys)
- `min_members`: least bound members for centres and mixing, `≥ 2`

Throws `ArgumentError` when `snaps` is empty, the bound sets do not match the snapshots one
to one, or `min_members < 2`.

# Example
```julia
rows = progenitor_budget(snaps, cluster_ranges, bound_sets; settings)
write_progenitor_budget("progenitor_budget.csv", rows)
```
"""
function progenitor_budget(
    snaps::Vector{Snapshot},
    cluster_ranges::AbstractVector{<:AbstractVector{Int}},
    bound_sets::AbstractVector{BoundSet};
    settings::PostprocessConfig = PostprocessConfig(),
    min_members::Integer = 20,
)::Vector{ProgenitorBudgetRow}
    isempty(snaps) && throw(ArgumentError("progenitor_budget needs at least one snapshot"))
    length(bound_sets) == length(snaps) || throw(
        ArgumentError(
            "progenitor_budget: $(length(bound_sets)) bound sets for $(length(snaps)) snapshots",
        ),
    )
    min_members ≥ 2 ||
        throw(ArgumentError("progenitor_budget: min_members must be ≥ 2; got $min_members"))
    n_cl = length(cluster_ranges)
    n_t = length(snaps)
    keywords = _centre_keywords(settings)

    # Initial names: slot = column of the first snapshot
    first_snap = snaps[1]
    n_1 = nparticles(first_snap)
    slot_of_name = Dict{Int32,Int}()
    sizehint!(slot_of_name, n_1)
    for (c, name) in enumerate(first_snap.name)
        slot_of_name[name] = c
    end
    m_1 = Float64.(first_snap.mass)
    m_last = copy(m_1)
    initial_slots = [_member_indices(first_snap, rng) for rng in cluster_ranges]
    n_initial = length.(initial_slots)
    m_initial = [sum(view(m_1, s); init = 0.0) for s in initial_slots]

    rows = Vector{ProgenitorBudgetRow}(undef, n_t * n_cl)
    present_slot = falses(n_1)
    slot_of_col = Int[]
    in_bound = falses(0)
    col_in_bound = Int[]
    tags01 = Int[]
    bound_cols = Int[]

    for (k, snap) in enumerate(snaps)
        n = nparticles(snap)
        h = snap.header
        t_nb = time_nb(h)
        t_myr = time_myr(h)

        resize!(slot_of_col, n)
        fill!(present_slot, false)
        @inbounds for c in 1:n
            s = get(slot_of_name, snap.name[c], 0)
            slot_of_col[c] = s
            if s != 0
                present_slot[s] = true
                m_last[s] = Float64(snap.mass[c])
            end
        end

        b = bound_sets[k].stars
        n_b = length(b)
        resize!(in_bound, n)
        fill!(in_bound, false)
        resize!(col_in_bound, n)
        fill!(col_in_bound, 0)
        @inbounds for (pos_in_b, c) in enumerate(b)
            in_bound[c] = true
            col_in_bound[c] = pos_in_b
        end
        resize!(tags01, n_b)

        global_ok = n_b ≥ min_members
        pos_b = Float64.(snap.pos[:, b])
        vel_b = Float64.(snap.vel[:, b])
        mass_b = Float64.(snap.mass[b])
        c_g = fill(NaN, 3)
        r_hg = NaN
        energy = nothing
        if global_ok
            c_g = shrinking_sphere_centre(pos_b, mass_b; keywords...).centre
            r_hg = _lagrangian_radii(pos_b, mass_b, c_g, (0.5,))[1]
            if length(snap.phi) == n
                v_ref = bound_sets[k].v_ref
                energy = Vector{Float64}(undef, n_b)
                @inbounds for (s, c) in enumerate(b)
                    energy[s] =
                        0.5 * (
                            (vel_b[1, s] - v_ref[1])^2 +
                            (vel_b[2, s] - v_ref[2])^2 +
                            (vel_b[3, s] - v_ref[3])^2
                        ) + Float64(snap.phi[c])
                end
            end
        end

        for i in 1:n_cl
            n_present = 0
            m_present = 0.0
            n_bound = 0
            m_bound = 0.0
            m_stellar_loss = 0.0
            empty!(bound_cols)
            @inbounds for c in _member_indices(snap, cluster_ranges[i])
                s = slot_of_col[c]
                s == 0 && continue
                m = Float64(snap.mass[c])
                n_present += 1
                m_present += m
                m_stellar_loss += m_1[s] - m
                if in_bound[c]
                    n_bound += 1
                    m_bound += m
                    push!(bound_cols, c)
                end
            end
            m_removed = 0.0
            @inbounds for s in initial_slots[i]
                present_slot[s] && continue
                m_removed += m_last[s]
                m_stellar_loss += m_1[s] - m_last[s]
            end

            r_h = NaN
            centre_offset = NaN
            aplus2 = NaN
            aplus2_err = NaN
            d_radius = NaN
            d_energy = NaN
            d_angular_momentum = NaN
            d_floor = NaN
            if global_ok && n_bound ≥ min_members
                pos_i = Float64.(snap.pos[:, bound_cols])
                mass_i = Float64.(snap.mass[bound_cols])
                c_i = shrinking_sphere_centre(pos_i, mass_i; keywords...).centre
                r_h = _lagrangian_radii(pos_i, mass_i, c_i, (0.5,))[1]
                centre_offset =
                    sqrt((c_i[1] - c_g[1])^2 + (c_i[2] - c_g[2])^2 + (c_i[3] - c_g[3])^2)
                if n_bound < n_b
                    fill!(tags01, 2)
                    @inbounds for c in bound_cols
                        tags01[col_in_bound[c]] = 1
                    end
                    mx = mixing(pos_b, vel_b, mass_b, tags01, c_g, r_hg; energy)
                    aplus2 = mx.aplus2
                    aplus2_err = mx.aplus2_error
                    d_radius = mx.d_radius
                    d_energy = mx.d_energy
                    d_angular_momentum = mx.d_angular_momentum
                    d_floor = mx.floor
                end
            end

            rows[(k - 1) * n_cl + i] = ProgenitorBudgetRow(
                t_nb,
                t_myr,
                i,
                n_initial[i],
                m_initial[i],
                n_present,
                m_present,
                n_bound,
                m_bound,
                m_removed,
                m_stellar_loss,
                r_h,
                centre_offset,
                aplus2,
                aplus2_err,
                d_radius,
                d_energy,
                d_angular_momentum,
                d_floor,
            )
        end
    end
    return rows
end

# ---------------------------------------------------------------------------
# CSV writers
# ---------------------------------------------------------------------------

"""
    write_merger_tree(path, tree::MergerTree) -> String

Write the events of [`merger_tree`](@ref) as CSV. A leading comment line records the number
of clusters, the number of snapshots and the dissolution index of every cluster
(`# n_clusters=3 n_snapshots=40 dissolved=0;0;27`); the columns are
`index,time_nb,time_myr,group_a,group_b,separation,threshold`, the groups as `;`-joined tags.
Non-finite values become empty cells. An existing file is backed up, never overwritten.
Returns the path.
"""
function write_merger_tree(path::AbstractString, tree::MergerTree)::String
    _backup_existing(path)
    n_cl, n_t = size(tree.present)
    open(path, "w") do io
        println(io, "# n_clusters=$(n_cl) n_snapshots=$(n_t) dissolved=", join(tree.dissolved, ";"))
        println(io, "index,time_nb,time_myr,group_a,group_b,separation,threshold")
        for e in tree.events
            cells = (
                _csv_precise(e.index),
                _csv_precise(e.time_nb),
                _csv_precise(e.time_myr),
                join(e.group_a, ";"),
                join(e.group_b, ";"),
                _csv_precise(e.separation),
                _csv_precise(e.threshold),
            )
            println(io, join(cells, ","))
        end
    end
    return String(path)
end

"""
    write_progenitor_budget(path, rows::AbstractVector{ProgenitorBudgetRow}) -> String

Write the rows of [`progenitor_budget`](@ref) as CSV, one column per field of
[`ProgenitorBudgetRow`](@ref) in declaration order (N-body units plus `time_myr`; non-finite
values as empty cells). An existing file is backed up, never overwritten. Returns the path.
"""
function write_progenitor_budget(
    path::AbstractString,
    rows::AbstractVector{ProgenitorBudgetRow},
)::String
    _backup_existing(path)
    n_fields = fieldcount(ProgenitorBudgetRow)
    open(path, "w") do io
        println(io, join(fieldnames(ProgenitorBudgetRow), ","))
        for r in rows
            println(io, join(ntuple(f -> _csv_precise(getfield(r, f)), n_fields), ","))
        end
    end
    return String(path)
end
