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

"""Default observables of the paired analysis, in table order."""
const _PAIRED_DEFAULT_OBSERVABLES = ["r10", "r50", "r90", "energy", "n_stars", "n_pairs"]

"""Refuse an observable that [`paired_comparison`](@ref) does not know, naming `caller`."""
function _check_paired_observables(observables::AbstractVector{<:AbstractString}, caller::String)
    for o in observables
        haskey(_PAIRED_OBSERVABLES, o) || throw(
            ArgumentError(
                "$(caller): unknown observable \"$o\"; one of " *
                join(sort!(collect(keys(_PAIRED_OBSERVABLES))), ", "),
            ),
        )
    end
    return nothing
end

"""Grid values of a sweep point as text: `"key=value"` pairs sorted and joined by commas; `""` without a grid."""
_grid_text(values::AbstractDict)::String =
    join(sort!(String["$(k)=$(v)" for (k, v) in values]), ", ")::String

"""
    _PairingGroups

The runs of a paired analysis grouped by grid point and seed:
`(grid, seed) => (arm => run_dir)`, the grid point as text ([`_grid_text`](@ref)).
"""
const _PairingGroups = Dict{Tuple{String,Int},Dict{String,String}}

"""Pairing groups of the completed points of a sweep (the arm `""` for a sweep without arms)."""
function _sweep_pairing_groups(sweep_dir::AbstractString)::_PairingGroups
    _, pts = _sweep_done_points(sweep_dir)
    groups = _PairingGroups()
    for (p, run_dir) in pts
        get!(groups, (_grid_text(p["values"]), Int(p["seed"])), Dict{String,String}())[string(
            get(p, "arm", ""),
        )::String] = run_dir
    end
    return groups
end

"""`merger_ic.toml` of a run's output as a dictionary; empty when absent."""
function _paired_run_meta(run_dir::AbstractString)::Dict{String,Any}
    path = joinpath(run_dir, "output", "merger_ic.toml")
    return isfile(path) ? TOML.parsefile(path) : Dict{String,Any}()
end

"""
    _verify_pairing(groups, reference, caller, unit)

Refuse a run of a non-reference arm whose `meta.dat10_sha256` differs from that of
the reference run of the same grid point and seed, unless it carries a perturbed
twin (`[twin]`); runs without a recorded hash pass. `unit` names a member in the
message (`"point"` of a sweep, `"run"` of a run list).
"""
function _verify_pairing(groups::_PairingGroups, reference::String, caller::String, unit::String)
    hash_of(m) = string(get(get(m, "meta", Dict{String,Any}()), "dat10_sha256", ""))::String
    for ((grid, seed), byarm) in groups
        haskey(byarm, reference) || continue
        h_ref = hash_of(_paired_run_meta(byarm[reference]))
        for (arm, run_dir) in byarm
            arm == reference && continue
            m = _paired_run_meta(run_dir)
            haskey(m, "twin") && continue
            h = hash_of(m)
            (isempty(h) || isempty(h_ref) || h == h_ref) || throw(
                ArgumentError(
                    "$(caller): the $(arm) $(unit) of seed $(seed)" *
                    (isempty(grid) ? "" : " at $(grid)") *
                    " does not share the initial conditions of the $(reference) $(unit) " *
                    "(dat.10 hashes differ)",
                ),
            )
        end
    end
    return nothing
end

"""Arms of the groups in name order; refuses a table without a run of `reference`."""
function _paired_arms(groups::_PairingGroups, reference::String, caller::String, unit::String)
    arms = sort!(unique(a for g in values(groups) for a in keys(g)))
    reference in arms ||
        throw(ArgumentError("$(caller): no completed $(unit) of the arm \"$reference\""))
    return arms
end

"""
    _arm_series(groups, grid, quantity, fraction) -> (seeds, series)

Seeds of the grid point `grid` in ascending order and the series of `quantity`
([`_run_series`](@ref)) of every run there, keyed by `(arm, seed)` (`nothing` for a
run without the source file).
"""
function _arm_series(groups::_PairingGroups, grid::String, quantity::Symbol, fraction::Real)
    seeds = sort!([s for (g, s) in keys(groups) if g == grid])
    series = Dict{Tuple{String,Int},Any}()
    for s in seeds, (arm, run_dir) in groups[(grid, s)]
        series[(arm, s)] = _run_series(run_dir, quantity; fraction = fraction)
    end
    return seeds, series
end

"""Both series of the pair `(arm, s)`–`(reference, s)` exist."""
_has_pair(series, arm, reference, s) =
    haskey(series, (arm, s)) &&
    haskey(series, (reference, s)) &&
    series[(arm, s)] !== nothing &&
    series[(reference, s)] !== nothing

"""
    _paired_rows(groups; reference, observables, t_settle_myr, margins, α, caller, unit)

The statistical core of [`paired_comparison`](@ref) and [`block_comparison`](@ref) on
explicit pairing groups: pairing check, per-seed window means of the differences,
chaos floor, seed scatter, margins, the paired tests and Holm's correction over all
rows. `caller` and `unit` name the entry point and a member in the messages.
"""
function _paired_rows(
    groups::_PairingGroups;
    reference::String,
    observables::AbstractVector{<:AbstractString},
    t_settle_myr::Real,
    margins::AbstractDict,
    α::Real,
    caller::String,
    unit::String,
)::Vector{NamedTuple}
    arms = _paired_arms(groups, reference, caller, unit)
    # Pairing: identical initial conditions unless the member is a twin.
    _verify_pairing(groups, reference, caller, unit)
    twin_arms =
        Set(a for g in values(groups) for (a, d) in g if haskey(_paired_run_meta(d), "twin"))
    rows = NamedTuple[]
    grids = sort!(unique(first.(collect(keys(groups)))))
    for grid in grids, o in observables
        quantity, fraction = _PAIRED_OBSERVABLES[o]
        seeds, series = _arm_series(groups, grid, quantity, fraction)
        diffs(arm) = Float64[
            _has_pair(series, arm, reference, s) ?
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
                    "$(caller): margin kind of $o must be \"chaos\", \"relative\" or \"absolute\"; got \"$kind\"",
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
correction runs over all rows of the table at family level `α`. The same table over an
explicit list of run directories is [`block_comparison`](@ref).
"""
function paired_comparison(
    sweep_dir::AbstractString;
    reference::AbstractString = "base",
    observables::AbstractVector{<:AbstractString} = _PAIRED_DEFAULT_OBSERVABLES,
    t_settle_myr::Real = 0.0,
    margins::AbstractDict = Dict{String,Tuple{String,Float64}}(),
    α::Real = 0.05,
)::Vector{NamedTuple}
    _check_paired_observables(observables, "paired_comparison")
    return _paired_rows(
        _sweep_pairing_groups(sweep_dir);
        reference = String(reference),
        observables,
        t_settle_myr,
        margins,
        α,
        caller = "paired_comparison",
        unit = "point",
    )
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
