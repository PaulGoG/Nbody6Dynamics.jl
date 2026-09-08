# =============================================================================
# Seeded ensembles: percentile bands over the seeds of a sweep point
# =============================================================================
#
# A sweep replicates every grid point over its seeds. Grouping the completed
# points by grid values gives one ensemble per point; every member's time
# series is interpolated onto a common grid (the interval covered by all
# members) and summarised by its median and the central 68 % and 95 %
# intervals. Time series come from the same extractors as the sweep
# comparison figures (`_run_series`).

"""
    EnsembleStatistics

Percentile summary of a seeded ensemble on a common time grid.

# Fields
- `time`: common grid [Myr]
- `median`, `q16`, `q84`, `q025`, `q975`: median and the 16th, 84th, 2.5th
  and 97.5th percentiles per node (central 68 % and 95 % intervals)
- `n`: number of members
"""
struct EnsembleStatistics
    time::Vector{Float64}
    median::Vector{Float64}
    q16::Vector{Float64}
    q84::Vector{Float64}
    q025::Vector{Float64}
    q975::Vector{Float64}
    n::Int
end

"""Nodes of the common time grid of an ensemble."""
const _ENSEMBLE_GRID_NODES = 200

"""
    _quantile_sorted(v, p) -> Float64

Quantile `p ∈ [0, 1]` of the sorted vector `v` with linear interpolation
between order statistics (Hyndman & Fan type 7, the default of R and
NumPy).
"""
function _quantile_sorted(v::AbstractVector{<:Real}, p::Real)
    n = length(v)
    n > 0 || throw(ArgumentError("quantile of an empty vector"))
    0 ≤ p ≤ 1 || throw(ArgumentError("quantile level must lie in [0, 1], got $p"))
    n == 1 && return Float64(v[1])
    h = (n - 1) * p + 1
    lo = floor(Int, h)
    hi = min(lo + 1, n)
    return Float64(v[lo]) + (h - lo) * (Float64(v[hi]) - Float64(v[lo]))
end

"""
    _interpolate_linear(x, y, xq) -> Vector{Float64}

Piecewise-linear interpolation of the samples `(x, y)` (`x` strictly
increasing) at the query points `xq`, which must lie within `[x[1], x[end]]`.
"""
function _interpolate_linear(
    x::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    xq::AbstractVector{<:Real},
)
    n = length(x)
    n == length(y) || throw(DimensionMismatch("x and y differ in length"))
    n ≥ 1 || throw(ArgumentError("interpolation needs at least one sample"))
    out = Vector{Float64}(undef, length(xq))
    for (k, q) in enumerate(xq)
        (x[1] - 1e-12 * abs(x[1]) ≤ q ≤ x[end] + 1e-12 * abs(x[end])) ||
            throw(ArgumentError("query $q outside the sampled interval [$(x[1]), $(x[end])]"))
        if n == 1 || q ≤ x[1]
            out[k] = y[1]
        elseif q ≥ x[end]
            out[k] = y[end]
        else
            j = searchsortedlast(x, q)
            w = (q - x[j]) / (x[j + 1] - x[j])
            out[k] = (1 - w) * y[j] + w * y[j + 1]
        end
    end
    return out
end

"""
    ensemble_statistics(members; n_grid = 200) -> EnsembleStatistics

Percentile summary of `members`, a vector of `(t, y)` series (each with
strictly increasing `t`): the common grid spans the interval covered by
every member with `n_grid` nodes, each member is interpolated linearly onto
it, and the median and the 2.5, 16, 84 and 97.5 % quantiles are taken per
node. With a single member every quantile equals the member.
"""
function ensemble_statistics(
    members::AbstractVector{<:Tuple{<:AbstractVector{<:Real},<:AbstractVector{<:Real}}};
    n_grid::Int = _ENSEMBLE_GRID_NODES,
)
    isempty(members) && throw(ArgumentError("ensemble_statistics: no members"))
    n_grid ≥ 2 || throw(ArgumentError("ensemble_statistics: n_grid must be ≥ 2"))
    for (t, _) in members
        length(t) ≥ 1 || throw(ArgumentError("ensemble_statistics: empty member series"))
        issorted(t; lt = <) && allunique(t) ||
            throw(ArgumentError("ensemble_statistics: member times must be strictly increasing"))
    end
    t_lo = maximum(first(t) for (t, _) in members)
    t_hi = minimum(last(t) for (t, _) in members)
    t_hi ≥ t_lo || throw(ArgumentError("ensemble_statistics: the members share no time interval"))
    grid = t_hi > t_lo ? collect(range(t_lo, t_hi; length = n_grid)) : [t_lo]
    n = length(members)
    values = Matrix{Float64}(undef, n, length(grid))
    for (i, (t, y)) in enumerate(members)
        values[i, :] = _interpolate_linear(t, y, grid)
    end
    med = similar(grid)
    q16 = similar(grid)
    q84 = similar(grid)
    q025 = similar(grid)
    q975 = similar(grid)
    col = Vector{Float64}(undef, n)
    for k in eachindex(grid)
        col .= @view values[:, k]
        sort!(col)
        med[k] = _quantile_sorted(col, 0.5)
        q16[k] = _quantile_sorted(col, 0.16)
        q84[k] = _quantile_sorted(col, 0.84)
        q025[k] = _quantile_sorted(col, 0.025)
        q975[k] = _quantile_sorted(col, 0.975)
    end
    return EnsembleStatistics(grid, med, q16, q84, q025, q975, n)
end

"""Physical unit scaling of a run from its stdout; identity when absent."""
function _run_scaling(out_dir::AbstractString)
    path = joinpath(out_dir, "out1000")
    isfile(path) || return UnitScaling(1.0, 1.0, 1.0, 1.0)
    return extract_scaling(read_diagnostics(path))
end

"""
    _run_series(run_dir, quantity; fraction = 0.5) -> Union{Nothing, Tuple}

Time series `(t_myr, y)` of one run for `quantity`:
`:lagrangian` (radius of mass fraction `fraction` [pc]), `:energy`
(|ΔE/E| per adjustment, zeros dropped), `:n_stars` and `:n_pairs` (from
the ADJUST records). `nothing` when the source file is absent or empty.
"""
function _run_series(run_dir::AbstractString, quantity::Symbol; fraction::Real = 0.5)
    out = joinpath(run_dir, "output")
    if quantity === :lagrangian
        path = joinpath(out, "lagr.7")
        isfile(path) || return nothing
        lagr = read_lagr(path)
        isempty(lagr.time) && return nothing
        k = argmin(abs.(lagr.mass_fractions .- fraction))
        u = _run_scaling(out)
        return (to_myr(u, lagr.time), to_pc(u, lagr.radii[k, :]))
    end
    path = joinpath(out, "out1000")
    isfile(path) || return nothing
    adj = read_diagnostics(path).adjust
    isempty(adj) && return nothing
    t = [a.time_myr for a in adj]
    if quantity === :energy
        keep = [isfinite(a.de_rel) && a.de_rel != 0 for a in adj]
        any(keep) || return nothing
        return (t[keep], abs.([a.de_rel for a in adj[keep]]))
    elseif quantity === :n_stars
        return (t, Float64[a.n for a in adj])
    elseif quantity === :n_pairs
        return (t, Float64[a.npairs for a in adj])
    end
    throw(
        ArgumentError(
            "unknown series quantity $quantity; use :lagrangian, :energy, :n_stars or :n_pairs",
        ),
    )
end

"""Axis label of a series quantity."""
function _series_label(quantity::Symbol, fraction::Real)
    quantity === :lagrangian &&
        return latexstring("r_{$(round(Int, 100 * fraction))\\,\\%} \\; [\\mathrm{pc}]")
    quantity === :energy && return L"|\Delta E / E|"
    quantity === :n_stars && return L"N"
    quantity === :n_pairs && return L"N_\mathrm{pairs}"
    throw(ArgumentError("unknown series quantity $quantity"))
end

"""
    sweep_ensembles(sweep_dir, quantity; fraction = 0.5)
        -> Vector{NamedTuple{(:values, :stats)}}

Group the completed points of a sweep by their grid values (one ensemble
per grid point, its members the seeds) and summarise the series of
`quantity` ([`_run_series`](@ref)) with [`ensemble_statistics`](@ref).
Groups whose members carry no series are skipped.
"""
function sweep_ensembles(sweep_dir::AbstractString, quantity::Symbol; fraction::Real = 0.5)
    _, pts = _sweep_done_points(sweep_dir)
    groups = Dict{Vector{Pair{String,Any}},Vector{Tuple{Vector{Float64},Vector{Float64}}}}()
    order = Vector{Vector{Pair{String,Any}}}()
    for (p, run_dir) in pts
        key = sort!([String(k) => v for (k, v) in p["values"]]; by = first)
        s = _run_series(run_dir, quantity; fraction = fraction)
        s === nothing && continue
        haskey(groups, key) ||
            (groups[key] = Tuple{Vector{Float64},Vector{Float64}}[]; push!(order, key))
        push!(groups[key], (Float64.(s[1]), Float64.(s[2])))
    end
    return [(values = key, stats = ensemble_statistics(groups[key])) for key in order]
end
