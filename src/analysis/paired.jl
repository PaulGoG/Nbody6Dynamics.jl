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

"""Refuse a series whose times are not strictly increasing (the linear interpolation assumes it)."""
function _check_series_times(t::AbstractVector{<:Real})
    (issorted(t; lt = <) && allunique(t)) || throw(
        ArgumentError(
            "paired comparison: series times must be strictly increasing (a restart joined " *
            "outside the package, or an unsorted series)",
        ),
    )
    return nothing
end

"""
    _window_mean_difference(a, ref, t_settle; nodes = 200) -> Float64

Mean over the window of the difference of two `(t, y)` series, `a − ref`: the window is
the interval both cover, from `t_settle` on; both are interpolated linearly on `nodes`
equally spaced points and the difference is averaged by the trapezoidal rule. `NaN` when
the window is empty.
"""
function _window_mean_difference(a, ref, t_settle::Real; nodes::Integer = 200)::Float64
    _check_series_times(a[1])
    _check_series_times(ref[1])
    t0 = max(first(a[1]), first(ref[1]), t_settle)
    t1 = min(last(a[1]), last(ref[1]))
    t1 > t0 || return NaN
    t = collect(range(t0, t1; length = nodes))
    d = _interpolate_linear(a[1], a[2], t) .- _interpolate_linear(ref[1], ref[2], t)
    return (sum(d) - (first(d) + last(d)) / 2) / (nodes - 1)
end

"""Mean of a `(t, y)` series over the window from `t_settle` to its end (trapezoidal rule on `nodes` points); `NaN` when the window is empty."""
function _window_mean(s, t_settle::Real; nodes::Integer = 200)::Float64
    _check_series_times(s[1])
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
        if !haskey(margins, o) && !isempty(twin_arms) && !isfinite(chaos_floor)
            @warn "$(caller): chaos floor of $o undefined (fewer than two finite twin " *
                  "differences); the margin falls back to 1 % of the reference window mean"
        end
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
        if !isfinite(margin)
            @warn "$(caller): margin of $o is undefined ($(kind) margin from a non-finite " *
                  "value); p_equivalence is NaN for its rows"
        end
        for arm in arms
            arm == reference && continue
            r = paired_shift(diffs(arm); margin = isfinite(margin) ? margin : 0.0, α = α)
            # An undefined margin tests nothing: no equivalence p-value.
            isfinite(margin) || (r = merge(r, (; p_equivalence = NaN)))
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

# -----------------------------------------------------------------------------
# Time-resolved paired differences and the LaTeX calibration table (4.6)
# -----------------------------------------------------------------------------

"""
    _difference_series(groups; reference, observables, n_grid, caller, unit) -> Vector{NamedTuple}

The core of [`paired_difference_series`](@ref) on explicit pairing groups: per grid
point, observable and arm the differences arm − reference of every seed with both
series, interpolated onto `n_grid` nodes of the interval all those pairs cover, and
their mean and standard error per node.
"""
function _difference_series(
    groups::_PairingGroups;
    reference::String,
    observables::AbstractVector{<:AbstractString},
    n_grid::Int,
    caller::String,
    unit::String,
)::Vector{NamedTuple}
    n_grid ≥ 2 || throw(ArgumentError("$(caller): n_grid must be ≥ 2; got $(n_grid)"))
    arms = _paired_arms(groups, reference, caller, unit)
    _verify_pairing(groups, reference, caller, unit)
    out = NamedTuple[]
    for grid in sort!(unique(first.(collect(keys(groups))))), o in observables
        quantity, fraction = _PAIRED_OBSERVABLES[o]
        seeds, series = _arm_series(groups, grid, quantity, fraction)
        for arm in arms
            arm == reference && continue
            pairs = [
                (series[(arm, s)], series[(reference, s)]) for
                s in seeds if _has_pair(series, arm, reference, s)
            ]
            n = length(pairs)
            t0 = n == 0 ? NaN : maximum(max(first(a[1]), first(r[1])) for (a, r) in pairs)
            t1 = n == 0 ? NaN : minimum(min(last(a[1]), last(r[1])) for (a, r) in pairs)
            if !(t1 > t0)
                push!(
                    out,
                    (;
                        grid,
                        observable = o,
                        arm,
                        n,
                        time = Float64[],
                        mean = Float64[],
                        se = Float64[],
                    ),
                )
                continue
            end
            t = collect(range(t0, t1; length = n_grid))
            d = Matrix{Float64}(undef, n, n_grid)
            for (i, (a, r)) in enumerate(pairs)
                d[i, :] = _interpolate_linear(a[1], a[2], t) .- _interpolate_linear(r[1], r[2], t)
            end
            m = vec(sum(d; dims = 1)) ./ n
            se =
                n ≥ 2 ? sqrt.(vec(sum(abs2, d .- m'; dims = 1)) ./ (n - 1)) ./ sqrt(n) :
                fill(NaN, n_grid)
            push!(out, (; grid, observable = o, arm, n, time = t, mean = m, se))
        end
    end
    return out
end

"""
    paired_difference_series(sweep_dir; reference = "base", observables = [...], n_grid = 200)
    paired_difference_series(runs; reference = "base", observables = [...], n_grid = 200)
        -> Vector{NamedTuple}

Time-resolved mean paired difference of every arm against the arm `reference`, the
companion of [`paired_comparison`](@ref) (a sweep directory) and
[`block_comparison`](@ref) (a run list `(; arm, seed, run_dir[, grid])`), with the same
pairing rules. Per grid point, observable and arm, the series of every seed with both
runs are differenced, arm minus reference, on `n_grid` equally spaced nodes of the
interval that all those pairs cover (from its start, the settling interval included),
and averaged over the seeds. Returns one
`(; grid, observable, arm, n, time, mean, se)` per grid point, observable and arm:
`time` [Myr], `mean` the mean difference and `se` its standard error (standard
deviation over the `n` seeds over √n; `NaN` for n < 2) per node. A row without
pairs or without a common interval has empty vectors.
"""
function paired_difference_series(
    sweep_dir::AbstractString;
    reference::AbstractString = "base",
    observables::AbstractVector{<:AbstractString} = _PAIRED_DEFAULT_OBSERVABLES,
    n_grid::Int = _ENSEMBLE_GRID_NODES,
)::Vector{NamedTuple}
    _check_paired_observables(observables, "paired_difference_series")
    return _difference_series(
        _sweep_pairing_groups(sweep_dir);
        reference = String(reference),
        observables,
        n_grid,
        caller = "paired_difference_series",
        unit = "point",
    )
end

function paired_difference_series(
    runs::AbstractVector{<:NamedTuple};
    reference::AbstractString = "base",
    observables::AbstractVector{<:AbstractString} = _PAIRED_DEFAULT_OBSERVABLES,
    n_grid::Int = _ENSEMBLE_GRID_NODES,
)::Vector{NamedTuple}
    _check_paired_observables(observables, "paired_difference_series")
    return _difference_series(
        _run_list_groups(runs, "paired_difference_series");
        reference = String(reference),
        observables,
        n_grid,
        caller = "paired_difference_series",
        unit = "run",
    )
end

"""Observables of the paired analysis as table entries: the symbol in math mode and the unit."""
const _OBSERVABLE_TEX = Dict(
    "r10" => raw"$r_{10\,\%}$ [pc]",
    "r50" => raw"$r_{50\,\%}$ [pc]",
    "r90" => raw"$r_{90\,\%}$ [pc]",
    "energy" => raw"$|\Delta E / E|$",
    "n_stars" => raw"$N$",
    "n_pairs" => raw"$N_\mathrm{pairs}$",
)

"""`s` with the LaTeX special characters escaped for text mode."""
function _latex_escape(s::AbstractString)::String
    io = IOBuffer()
    for c in s
        if c in ('&', '%', '$', '#', '_', '{', '}')
            print(io, '\\', c)
        elseif c == '\\'
            print(io, "\\textbackslash{}")
        elseif c == '~'
            print(io, "\\textasciitilde{}")
        elseif c == '^'
            print(io, "\\textasciicircum{}")
        else
            print(io, c)
        end
    end
    return String(take!(io))
end

"""
    _tex_fixed(x, digits) -> String

`x` with `digits` decimals as math-mode text (`"-0.10"`), without a negative zero.
"""
function _tex_fixed(x::Real, digits::Int)::String
    s = @sprintf("%.*f", max(digits, 0), Float64(x))
    return all(c -> c in ('-', '0', '.'), s) ? lstrip(s, '-') : s
end

"""
    _tex_sig(x, n = 3) -> String

`x` to `n` significant digits as math-mode text: a plain decimal for
10⁻² ≤ |x| < 10⁴ (and 0), `m \\times 10^{e}` otherwise (`10^{e}` for m = 1).
"""
function _tex_sig(x::Real, n::Int = 3)::String
    isfinite(x) || return "--"
    x == 0 && return "0"
    v = round(Float64(x); sigdigits = n)
    e = floor(Int, log10(abs(v)))
    if -2 ≤ e ≤ 3
        return _tex_fixed(v, max(0, n - 1 - e))
    end
    m = round(v / 10.0^e; sigdigits = n)
    abs(m) ≥ 10 && ((m, e) = (m / 10, e + 1))         # log10 rounding at a decade
    abs(m) == 1 && return (m < 0 ? "-" : "") * "10^{$(e)}"
    return _tex_fixed(m, n - 1) * " \\times 10^{$(e)}"
end

"""
    _tex_shift_cells(shift, se, lo, hi) -> (shift_cell, interval_cell)

The shift with its standard error and the interval as math-mode table cells: the
standard error to two significant digits and the shift and the interval bounds to the same
decimal place, but the shift to three significant digits at most (a standard error below
that resolution shows as zero), and a common power of ten factored out when the larger
of |shift| and the standard error lies outside 10⁻²–10⁴. `--` for a non-finite shift.
"""
function _tex_shift_cells(shift::Real, se::Real, lo::Real, hi::Real)
    isfinite(shift) || return ("--", "--")
    scale_of = max(abs(shift), isfinite(se) ? se : 0.0)
    e = scale_of == 0 ? 0 : floor(Int, log10(scale_of))
    e = -2 ≤ e ≤ 3 ? 0 : e
    f = 10.0^(-e)
    s, u, a, b = shift * f, se * f, lo * f, hi * f
    # Two significant digits of the standard error, at most three of the shift: a
    # standard error far below the shift (round-off of a constant difference) shows as 0.
    shift_digits = s == 0 ? 0 : max(0, 2 - floor(Int, log10(abs(s))))
    digits = if isfinite(u) && u > 0
        se_digits = max(0, 1 - floor(Int, log10(round(u; sigdigits = 2))))
        s == 0 ? se_digits : min(se_digits, shift_digits)
    else
        shift_digits
    end
    pm = isfinite(u) ? " \\pm " * _tex_fixed(u, digits) : ""
    interval =
        (isfinite(a) && isfinite(b)) ?
        "[" * _tex_fixed(a, digits) * ", " * _tex_fixed(b, digits) * "]" : "--"
    factor = e == 0 ? "" : " \\times 10^{$(e)}"
    shift_cell =
        e == 0 ? "\$" * _tex_fixed(s, digits) * pm * "\$" :
        "\$(" * _tex_fixed(s, digits) * pm * ")" * factor * "\$"
    interval_cell = interval == "--" ? "--" : "\$" * interval * factor * "\$"
    return (shift_cell, interval_cell)
end

"""
    write_calibration_table(path, rows; α = 0.05) -> String

Write the rows of [`paired_comparison`](@ref) or [`block_comparison`](@ref) as a LaTeX
`tabular` in booktabs style (`\\toprule`, `\\midrule`, `\\bottomrule`,
`\\addlinespace` between observables; the document needs `\\usepackage{booktabs}`):
observable (symbol and unit, named once per group), arm, number of seeds n, shift ±
standard error, the (1 − 2α) interval, the margin ε and the verdict; a grid-point column
is added when the rows span several. Numbers carry two (standard error, and the shift
and interval to the same decimal place) or three (margin) significant digits, powers of
ten as `\\times 10^{e}`; arm names and grid values are escaped. `α` must be the level
the rows were computed at. An existing file is kept as a `#k` backup. Returns `path`.
"""
function write_calibration_table(path::AbstractString, rows::AbstractVector; α::Real = 0.05)::String
    isempty(rows) && throw(ArgumentError("write_calibration_table: no rows"))
    0 < α < 0.5 || throw(ArgumentError("write_calibration_table: α must lie in (0, 0.5); got $α"))
    for r in rows
        haskey(_OBSERVABLE_TEX, r.observable) ||
            throw(ArgumentError("write_calibration_table: unknown observable \"$(r.observable)\""))
    end
    grids = unique(String(r.grid) for r in rows)
    with_grid = length(grids) > 1
    level = 100 * (1 - 2α)
    level_text = isinteger(round(level; digits = 6)) ? string(round(Int, level)) : string(level)
    _backup_existing(path)
    open(path, "w") do io
        println(io, "\\begin{tabular}{", with_grid ? "l" : "", "llrllll}")
        println(io, "\\toprule")
        println(
            io,
            with_grid ? "Grid point & " : "",
            "Observable & Arm & \$n\$ & Shift \$\\pm\$ SE & \$$(level_text)\\,\\%\$ interval & ",
            "Margin \$\\varepsilon\$ & Verdict \\\\",
        )
        println(io, "\\midrule")
        previous = nothing
        for r in rows
            group = (String(r.grid), String(r.observable))
            previous !== nothing && group != previous && println(io, "\\addlinespace")
            first_of_group = group != previous
            previous = group
            shift_cell, interval_cell = _tex_shift_cells(r.shift, r.se, r.ci_low, r.ci_high)
            cells = String[]
            with_grid && push!(cells, first_of_group ? _latex_escape(r.grid) : "")
            push!(cells, first_of_group ? _OBSERVABLE_TEX[r.observable] : "")
            push!(cells, _latex_escape(r.arm), string(r.n), shift_cell, interval_cell)
            push!(cells, isfinite(r.margin) ? "\$" * _tex_sig(r.margin) * "\$" : "--")
            push!(cells, _latex_escape(r.verdict))
            println(io, join(cells, " & "), " \\\\")
        end
        println(io, "\\bottomrule")
        println(io, "\\end{tabular}")
    end
    return path
end
