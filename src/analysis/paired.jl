# =============================================================================
# Paired comparison of sweep arms: chaos floor, paired shifts, equivalence verdicts
# =============================================================================
#
# The arms of a sweep (4.3) integrate the same initial conditions with different numerical settings. Every arm is
# compared with a reference arm seed by seed: per seed the window mean of the difference of an observable is one
# sample, the seeds are the independent realisations, and the arm's shift is tested for equivalence to zero within
# a declared margin (two one-sided tests) and for a difference (paired t test), with Holm's correction over the
# table. The perturbed-twin arm measures the chaos floor that sets the default margin. Design:
# notes of 2026-11-04 (paired comparison); principle: Goodman, Heggie & Hut (1993), DOI 10.1086/173196;
# Portegies Zwart & Boekholt (2014), DOI 10.1088/2041-8205/785/1/L3.

"""Student-t distribution function at `x` with `ν > 0` degrees of freedom (regularised incomplete beta)."""
function _student_t_cdf(x::Real, ν::Real)::Float64
    ν > 0 || throw(ArgumentError("degrees of freedom must be > 0; got $ν"))
    isinf(x) && return x > 0 ? 1.0 : 0.0
    tail = 0.5 * first(SpecialFunctions.beta_inc(ν / 2, 0.5, ν / (ν + Float64(x)^2)))
    return x ≥ 0 ? 1 - tail : tail
end

"""Quantile of the Student-t distribution with `ν > 0` degrees of freedom at probability `0 < p < 1`."""
function _student_t_quantile(p::Real, ν::Real)::Float64
    0 < p < 1 || throw(ArgumentError("probability must lie in (0, 1); got $p"))
    ν > 0 || throw(ArgumentError("degrees of freedom must be > 0; got $ν"))
    p == 0.5 && return 0.0
    q = min(p, 1 - p)
    z = first(SpecialFunctions.beta_inc_inv(ν / 2, 0.5, 2q))
    t = sqrt(ν * (1 - z) / z)
    return p > 0.5 ? t : -t
end

"""
    holm_adjust(p) -> Vector{Float64}

Holm–Bonferroni adjusted p-values of the family `p` (non-finite entries are left out of
the family and returned as `NaN`): the k-th smallest of m finite values is multiplied by
m − k + 1, made monotone and capped at 1. A test of the family is rejected at family
level α when its adjusted value is below α.
"""
function holm_adjust(p::AbstractVector{<:Real})::Vector{Float64}
    out = fill(NaN, length(p))
    idx = findall(isfinite, p)
    m = length(idx)
    running = 0.0
    for (k, i) in enumerate(idx[sortperm(p[idx])])
        running = max(running, min(1.0, (m - k + 1) * p[i]))
        out[i] = running
    end
    return out
end

"""
    paired_shift(d; margin, α = 0.05) -> NamedTuple

Paired analysis of the per-seed differences `d` (one value per seed, arm minus reference):
`(; n, mean, se, ci_low, ci_high, p_equivalence, p_difference)`. `ci_low, ci_high` is the
(1 − 2α) interval of the mean shift with n − 1 degrees of freedom; `p_equivalence` the
larger p-value of the two one-sided tests of |shift| ≥ `margin` (equivalence when it lies
below α); `p_difference` the two-sided p-value of a zero shift. Fewer than two finite
samples give `NaN` statistics.

# Example
```julia
paired_shift([0.01, -0.02, 0.005, 0.0]; margin = 0.1).p_equivalence < 0.05
```
"""
function paired_shift(d::AbstractVector{<:Real}; margin::Real, α::Real = 0.05)
    margin ≥ 0 || throw(ArgumentError("margin must be ≥ 0; got $margin"))
    0 < α < 0.5 || throw(ArgumentError("α must lie in (0, 0.5); got $α"))
    x = Float64[v for v in d if isfinite(v)]
    n = length(x)
    nan = (;
        n,
        mean = NaN,
        se = NaN,
        ci_low = NaN,
        ci_high = NaN,
        p_equivalence = NaN,
        p_difference = NaN,
    )
    n ≥ 2 || return nan
    m = sum(x) / n
    se = sqrt(sum(abs2, x .- m) / (n - 1)) / sqrt(n)
    ν = n - 1
    tq = _student_t_quantile(1 - α, ν)
    if se == 0
        return (;
            n,
            mean = m,
            se,
            ci_low = m,
            ci_high = m,
            p_equivalence = abs(m) < margin ? 0.0 : 1.0,
            p_difference = m == 0 ? 1.0 : 0.0,
        )
    end
    p_lower = 1 - _student_t_cdf((m + margin) / se, ν)   # H₀: shift ≤ −margin
    p_upper = _student_t_cdf((m - margin) / se, ν)       # H₀: shift ≥ +margin
    p_difference = 2 * (1 - _student_t_cdf(abs(m) / se, ν))
    return (;
        n,
        mean = m,
        se,
        ci_low = m - tq * se,
        ci_high = m + tq * se,
        p_equivalence = max(p_lower, p_upper),
        p_difference,
    )
end

"""Observables of [`paired_comparison`](@ref) by name: the series quantity of `_run_series` and its mass fraction."""
const _PAIRED_OBSERVABLES = Dict(
    "r10" => (:lagrangian, 0.1),
    "r50" => (:lagrangian, 0.5),
    "r90" => (:lagrangian, 0.9),
    "energy" => (:energy, 0.5),
    "n_stars" => (:n_stars, 0.5),
    "n_pairs" => (:n_pairs, 0.5),
)

"""
    _window_mean_difference(a, ref, t_settle; nodes = 200) -> Float64

Mean over the window of the difference of two `(t, y)` series, `a − ref`: the window is
the interval both cover, from `t_settle` on; both are interpolated linearly on `nodes`
equally spaced points and the difference is averaged by the trapezoidal rule. `NaN` when
the window is empty.
"""
function _window_mean_difference(a, ref, t_settle::Real; nodes::Integer = 200)::Float64
    t0 = max(first(a[1]), first(ref[1]), t_settle)
    t1 = min(last(a[1]), last(ref[1]))
    t1 > t0 || return NaN
    t = collect(range(t0, t1; length = nodes))
    d = _interpolate_linear(a[1], a[2], t) .- _interpolate_linear(ref[1], ref[2], t)
    return (sum(d) - (first(d) + last(d)) / 2) / (nodes - 1)
end

"""Mean of a `(t, y)` series over the window from `t_settle` to its end (trapezoidal rule on `nodes` points); `NaN` when the window is empty."""
function _window_mean(s, t_settle::Real; nodes::Integer = 200)::Float64
    t0 = max(first(s[1]), t_settle)
    t1 = last(s[1])
    t1 > t0 || return NaN
    y = _interpolate_linear(s[1], s[2], collect(range(t0, t1; length = nodes)))
    return (sum(y) - (first(y) + last(y)) / 2) / (nodes - 1)
end

"""
    paired_comparison(sweep_dir; reference = "base", observables = ["r10", "r50", "r90",
                      "energy", "n_stars", "n_pairs"], t_settle_myr = 0.0,
                      margins = Dict(), α = 0.05) -> Vector{NamedTuple}

Calibration table of the arms of a sweep (`[sweep.arms]`) against the arm `reference`:
one row per grid point, observable and arm, from the completed points only. Every point is
paired with the reference point of the same grid point and seed; the pairing is verified
by equal `meta.dat10_sha256` of the two initial conditions, except for points carrying a
perturbed twin (`[twin]` in `merger_ic.toml`), and refused otherwise. Per pair the
difference of the observable's series is averaged over the window from `t_settle_myr` to
the end of the common interval; the seeds are the samples of [`paired_shift`](@ref).

Columns: `grid` (the grid values as text), `observable`, `arm`, `n`, `shift` (mean
difference), `se`, `ci_low`, `ci_high` (the (1 − 2α) interval), `margin`, `margin_kind`,
`chaos_floor` (standard deviation over seeds of the window-mean difference of the twin
arm, NaN without one), `seed_scatter` (standard deviation over seeds of the reference
window mean), `p_equivalence`, `p_equivalence_holm`, `p_difference`,
`p_difference_holm`, `verdict` (`"equivalent"`, `"different"`, `"inconclusive"`).

`margins` maps an observable to `("chaos", κ)` (κ times the chaos floor; the default
with κ = 2 when the sweep has a twin arm), `("relative", f)` (f times the absolute
reference window mean; the default with f = 0.01 otherwise) or `("absolute", ε)`. Holm's
correction runs over all rows of the table at family level `α`.
"""
function paired_comparison(
    sweep_dir::AbstractString;
    reference::AbstractString = "base",
    observables::AbstractVector{<:AbstractString} = [
        "r10",
        "r50",
        "r90",
        "energy",
        "n_stars",
        "n_pairs",
    ],
    t_settle_myr::Real = 0.0,
    margins::AbstractDict = Dict{String,Tuple{String,Float64}}(),
    α::Real = 0.05,
)::Vector{NamedTuple}
    for o in observables
        haskey(_PAIRED_OBSERVABLES, o) || throw(
            ArgumentError(
                "paired_comparison: unknown observable \"$o\"; one of " *
                join(sort!(collect(keys(_PAIRED_OBSERVABLES))), ", "),
            ),
        )
    end
    _, pts = _sweep_done_points(sweep_dir)
    groups = Dict{Tuple{String,Int},Dict{String,String}}()   # (grid, seed) => arm => run dir
    grid_of = Dict{String,String}()
    for (p, run_dir) in pts
        grid = join(sort!(["$(k)=$(v)" for (k, v) in p["values"]]), ", ")
        get!(groups, (grid, Int(p["seed"])), Dict{String,String}())[String(get(p, "arm", ""))] =
            run_dir
    end
    arms = sort!(unique(a for g in values(groups) for a in keys(g)))
    reference in arms ||
        throw(ArgumentError("paired_comparison: no completed point of the arm \"$reference\""))
    meta(run_dir) = (
        path = joinpath(run_dir, "output", "merger_ic.toml");
        isfile(path) ? TOML.parsefile(path) : Dict{String,Any}()
    )
    # Pairing: identical initial conditions unless the point is a twin.
    for ((grid, seed), byarm) in groups
        haskey(byarm, reference) || continue
        h_ref = get(get(meta(byarm[reference]), "meta", Dict{String,Any}()), "dat10_sha256", "")
        for (arm, run_dir) in byarm
            arm == reference && continue
            m = meta(run_dir)
            haskey(m, "twin") && continue
            h = get(get(m, "meta", Dict{String,Any}()), "dat10_sha256", "")
            (isempty(h) || isempty(h_ref) || h == h_ref) || throw(
                ArgumentError(
                    "paired_comparison: the $(arm) point of seed $(seed) at $(grid) does not share " *
                    "the initial conditions of the $(reference) point (dat.10 hashes differ)",
                ),
            )
        end
    end
    twin_arms = Set(a for g in values(groups) for (a, d) in g if haskey(meta(d), "twin"))
    rows = NamedTuple[]
    grids = sort!(unique(first.(collect(keys(groups)))))
    for grid in grids, o in observables
        quantity, fraction = _PAIRED_OBSERVABLES[o]
        series = Dict{Tuple{String,Int},Any}()
        seeds = sort!([s for (g, s) in keys(groups) if g == grid])
        for s in seeds, (arm, run_dir) in groups[(grid, s)]
            series[(arm, s)] = _run_series(run_dir, quantity; fraction = fraction)
        end
        diffs(arm) = Float64[
            (
                haskey(series, (arm, s)) &&
                haskey(series, (reference, s)) &&
                series[(arm, s)] !== nothing &&
                series[(reference, s)] !== nothing
            ) ?
            _window_mean_difference(series[(arm, s)], series[(reference, s)], t_settle_myr) :
            NaN for s in seeds
        ]
        sd(x) = (
            y = filter(isfinite, x);
            length(y) < 2 ? NaN : sqrt(sum(abs2, y .- sum(y) / length(y)) / (length(y) - 1))
        )
        floors = [sd(diffs(a)) for a in twin_arms]
        chaos_floor = isempty(floors) ? NaN : maximum(floors)
        ref_means = Float64[
            series[(reference, s)] === nothing ? NaN :
            _window_mean(series[(reference, s)], t_settle_myr) for
            s in seeds if haskey(series, (reference, s))
        ]
        seed_scatter = sd(ref_means)
        kind, value = get(margins, o, isfinite(chaos_floor) ? ("chaos", 2.0) : ("relative", 0.01))
        margin =
            kind == "chaos" ? value * chaos_floor :
            kind == "relative" ?
            value * abs(sum(filter(isfinite, ref_means)) / max(1, count(isfinite, ref_means))) :
            kind == "absolute" ? Float64(value) :
            throw(
                ArgumentError(
                    "paired_comparison: margin kind of $o must be \"chaos\", \"relative\" or \"absolute\"; got \"$kind\"",
                ),
            )
        for arm in arms
            arm == reference && continue
            r = paired_shift(diffs(arm); margin = isfinite(margin) ? margin : 0.0, α = α)
            push!(
                rows,
                (;
                    grid,
                    observable = o,
                    arm,
                    n = r.n,
                    shift = r.mean,
                    se = r.se,
                    ci_low = r.ci_low,
                    ci_high = r.ci_high,
                    margin,
                    margin_kind = String(kind),
                    chaos_floor,
                    seed_scatter,
                    p_equivalence = r.p_equivalence,
                    p_equivalence_holm = NaN,
                    p_difference = r.p_difference,
                    p_difference_holm = NaN,
                    verdict = "",
                ),
            )
        end
    end
    pe = holm_adjust([r.p_equivalence for r in rows])
    pd = holm_adjust([r.p_difference for r in rows])
    return [
        merge(
            r,
            (;
                p_equivalence_holm = pe[k],
                p_difference_holm = pd[k],
                verdict = pe[k] < α ? "equivalent" : pd[k] < α ? "different" : "inconclusive",
            ),
        ) for (k, r) in enumerate(rows)
    ]
end

"""
    write_paired_comparison(path, rows) -> String

Write the rows of [`paired_comparison`](@ref) as CSV (non-finite numbers as empty cells);
an existing file is kept as a `#k` backup. Returns `path`.
"""
function write_paired_comparison(path::AbstractString, rows::AbstractVector)::String
    isempty(rows) && throw(ArgumentError("write_paired_comparison: no rows"))
    cols = collect(keys(first(rows)))
    _backup_existing(path)
    open(path, "w") do io
        println(io, join(String.(cols), ","))
        for r in rows
            println(io, join([_csv_cell(r[c]) for c in cols], ","))
        end
    end
    return path
end
