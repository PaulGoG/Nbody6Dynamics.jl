# =============================================================================
# Per-run validity record
# =============================================================================
# One record per post-processed run: the conservation checks (energy, mass,
# angular momentum), the state of the bound-set iteration, the snapshot cadence
# against the dynamical time scales, the stages of the evolution (violent
# relaxation, pericentres, virialisation, kinematic relaxation, mixing), the
# calibrated bias and scatter of the tracked estimators at the run's particle
# number, and the masking rules the estimators apply. Times in N-body units.

"""Schema version of `validity.toml` and `validity_snapshots.csv`."""
const VALIDITY_SCHEMA_VERSION = 2

"""Directory of the shipped estimator calibration (`estimators.csv`, `estimators.toml`)."""
const _CALIBRATION_DIR = joinpath(_PACKAGE_ROOT, "assets", "calibration")

"""Estimator quantities followed snapshot by snapshot in validity_snapshots.csv."""
const VALIDITY_TRACKED = [
    ("core_radius", "r_core_number"),
    ("lagrangian_grid", "ln_r_0.1_estimated_centre"),
    ("lagrangian_grid", "ln_r_0.5_estimated_centre"),
    ("rotation", "lambda_r"),
    ("shape_tensor", "q_shell2"),
    ("mixing", "aplus2"),
    ("centre", "offset"),
]

"""Epoch kinds that exclude a snapshot from the equilibrium estimators."""
const _VALIDITY_WINDOW_KINDS = (:violent_relaxation, :pericentre, :coalescence_phase)

"""Sections of a `ValidityRecord`, in the order `write_validity` writes them."""
const _VALIDITY_SECTIONS = (
    "header",
    "thresholds",
    "energy",
    "mass",
    "angular_momentum",
    "bound_set",
    "capacity",
    "cadence",
    "stages",
    "calibration",
    "masks",
)

# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

"""
    CalibrationRow

One row of the estimator calibration table `estimators.csv`: the static
`model` (`"plummer"`, `"plummer_kroupa"`, …), the `estimator` and the
`quantity` it reports, the particle number `n`, the known answer `truth`
(NaN without one), the replica `mean` and `scatter` (standard deviation over
the replicas), `bias = mean − truth` with its error `bias_err`, and the mean
`reported_error` of the estimator's own error estimate (NaN when it gives
none).
"""
struct CalibrationRow
    model::String
    estimator::String
    quantity::String
    n::Int
    truth::Float64
    mean::Float64
    scatter::Float64
    bias::Float64
    bias_err::Float64
    reported_error::Float64
end

"""
    CalibrationTable

The estimator calibration read by `load_calibration`: its `version`
(the `version` key of `estimators.toml`), the directory `dir` it was read
from and its `rows`, one `CalibrationRow` per line of
`estimators.csv`.
"""
struct CalibrationTable
    version::Int
    dir::String
    rows::Vector{CalibrationRow}
end

"""
    CalibratedEstimate

Calibrated bias and scatter of one estimator at a particle number, from
`calibrated`:

- `bias`, `scatter`: relative to |truth| when `relative`, otherwise absolute
  (a noise floor of a quantity whose true value is zero or unknown)
- `error_ratio`: `reported_error / scatter`; NaN when either is NaN or the
  scatter is zero
- `relative`: the interpolated truth is finite and nonzero
- `in_range`: the particle number lies inside the calibrated grid
- `ok`: `isfinite(bias) && abs(bias) ≤ tolerance` when `relative`, `true`
  otherwise
"""
struct CalibratedEstimate
    bias::Float64          # relative to |truth| when `relative`, else absolute (a noise floor)
    scatter::Float64       # same convention
    error_ratio::Float64   # reported_error / scatter; NaN if either is NaN or scatter == 0
    relative::Bool         # truth finite and nonzero
    in_range::Bool         # n inside [smallest N, largest N] of the rows
    ok::Bool               # relative ? (isfinite(bias) && abs(bias) ≤ tolerance) : true
end

"""Estimate of a snapshot without calibration or without bound stars."""
const _VALIDITY_NO_ESTIMATE = CalibratedEstimate(NaN, NaN, NaN, false, false, false)

"""
    ValidityEpoch

A time interval `[t_start, t_end]` (N-body time) of a run in which a check
of `validity_record` applies: `kind` is one of
`:violent_relaxation`, `:pericentre`, `:coalescence_phase` (the windows that
exclude snapshots from the equilibrium estimators), `:energy_interval` (an
adjustment interval flagged by `energy_intervals`), `:energy_threshold` (an
adjustment interval whose |DE| exceeds `validity_interval_de`),
`:mass_budget` (a mass-budget residual above `validity_mass_tolerance`) or
`:bound_set` (a bound set that did not converge); `detail` describes it.
"""
struct ValidityEpoch
    t_start::Float64
    t_end::Float64   # N-body time
    kind::Symbol     # one of :violent_relaxation, :pericentre, :coalescence_phase, :energy_interval, :energy_threshold, :mass_budget, :bound_set
    detail::String
end

"""
    ValiditySnapshot

The validity state of one snapshot (one entry of the snapshot integrals),
from `validity_record`:

- `time_nb`; `n_bound` bound stars and their effective number
  `n_eff = (Σm)² / Σm²`
- `t_cr`: crossing time `M^{5/2} / (−2E)^{3/2}` of the bound set (`G = 1`;
  Heggie & Mathieu 1986), `M` its mass and `E = T + W` its energy
- `q_virial = 2T / |W|` of the bound set
- `lagr_f_min = lagr_n_min / n_eff` and `lagr_f_max = 1 − lagr_n_min / n_eff`,
  the mass-fraction range a Lagrangian radius can resolve; NaN when
  `n_eff < 2 lagr_n_min`
- `bound_converged`: the bound-set iteration converged
- `in_window`: the snapshot lies inside a `:violent_relaxation`,
  `:pericentre` or `:coalescence_phase` epoch
- `equilibrium_ok = !in_window && bound_converged`
- `calibration`: one `CalibratedEstimate` per entry of
  `VALIDITY_TRACKED`, in that order, at `n_bound`
"""
struct ValiditySnapshot
    time_nb::Float64
    n_bound::Int
    n_eff::Float64            # (Σm)² / Σm² over the bound stars
    t_cr::Float64
    q_virial::Float64         # 2 T_b / |W_b|
    lagr_f_min::Float64       # lagr_n_min / n_eff; NaN when n_eff < 2 lagr_n_min
    lagr_f_max::Float64       # 1 − lagr_n_min / n_eff; NaN likewise
    bound_converged::Bool
    in_window::Bool           # inside a :violent_relaxation, :pericentre or :coalescence_phase epoch
    equilibrium_ok::Bool      # !in_window && bound_converged
    calibration::Vector{CalibratedEstimate}   # one per VALIDITY_TRACKED entry, in that order
end

"""
    ValidityRecord

The validity record of a run, from `validity_record`:
`schema_version` (`VALIDITY_SCHEMA_VERSION`); `sections`, the tables
`header`, `thresholds`, `energy`, `mass`, `angular_momentum`, `bound_set`,
`capacity`, `cadence`, `stages`, `calibration` and `masks`; the `epochs`,
sorted by `t_start`; and one `ValiditySnapshot` per snapshot integral.
"""
struct ValidityRecord
    schema_version::Int
    sections::Dict{String,Dict{String,Any}}   # header, thresholds, energy, mass, angular_momentum, bound_set, capacity, cadence, stages, calibration, masks
    epochs::Vector{ValidityEpoch}
    snapshots::Vector{ValiditySnapshot}
end

# ---------------------------------------------------------------------------
# Calibration table
# ---------------------------------------------------------------------------

"""
    load_calibration(dir = _CALIBRATION_DIR) -> CalibrationTable

Read the estimator calibration of `dir`: the integer `version` of
`estimators.toml` and every data line of `estimators.csv` (columns
`model,estimator,quantity,N,n_rep,truth,mean,scatter,bias,bias_err,reported_error,predicted_scatter`;
`NaN` literals mark unavailable values; `n_rep` and `predicted_scatter` are
not kept). The default is the calibration shipped with the package.

Throws an `ArgumentError` when either file is missing, when `version` is not
an integer, or for a line whose field count is not 12 (naming the line).

# Example
```julia
cal = load_calibration()
calibrated(cal, "plummer", "rotation", "lambda_r", 20_000)
```
"""
function load_calibration(dir::AbstractString = _CALIBRATION_DIR)::CalibrationTable
    meta_path = joinpath(dir, "estimators.toml")
    table_path = joinpath(dir, "estimators.csv")
    isfile(meta_path) || throw(ArgumentError("calibration metadata not found: $meta_path"))
    isfile(table_path) || throw(ArgumentError("calibration table not found: $table_path"))
    version = get(TOML.parsefile(meta_path), "version", nothing)
    version isa Int ||
        throw(ArgumentError("$meta_path: version must be an integer; got $(repr(version))"))
    rows = CalibrationRow[]
    for (lineno, line) in enumerate(eachline(table_path))
        lineno == 1 && continue
        text = strip(line)
        isempty(text) && continue
        f = split(text, ',')
        length(f) == 12 ||
            throw(ArgumentError("$table_path, line $lineno: expected 12 fields, got $(length(f))"))
        push!(
            rows,
            CalibrationRow(
                String(f[1]),
                String(f[2]),
                String(f[3]),
                parse(Int, f[4]),
                parse(Float64, f[6]),
                parse(Float64, f[7]),
                parse(Float64, f[8]),
                parse(Float64, f[9]),
                parse(Float64, f[10]),
                parse(Float64, f[11]),
            ),
        )
    end
    return CalibrationTable(version, String(abspath(dir)), rows)
end

"""
    _validity_calibration_rows(cal, model, estimator, quantity) -> Vector{CalibrationRow}

Rows of `cal` for the triple, sorted by particle number; throws an
`ArgumentError` naming the triple when there is none.
"""
function _validity_calibration_rows(
    cal::CalibrationTable,
    model::AbstractString,
    estimator::AbstractString,
    quantity::AbstractString,
)::Vector{CalibrationRow}
    rows = CalibrationRow[
        r for
        r in cal.rows if r.model == model && r.estimator == estimator && r.quantity == quantity
    ]
    isempty(rows) && throw(
        ArgumentError(
            "calibration version $(cal.version) has no rows for model \"$model\", " *
            "estimator \"$estimator\", quantity \"$quantity\"",
        ),
    )
    sort!(rows; by = r -> r.n)
    return rows
end

"""
    _validity_interpolate(rows, n, tolerance) -> CalibratedEstimate

The estimate of `calibrated` from the rows of one triple, sorted by
particle number.
"""
function _validity_interpolate(
    rows::AbstractVector{CalibrationRow},
    n::Real,
    tolerance::Real,
)::CalibratedEstimate
    (isfinite(n) && n > 0) || throw(ArgumentError("calibrated: n must be finite and > 0; got $n"))
    (isfinite(tolerance) && tolerance ≥ 0) ||
        throw(ArgumentError("calibrated: tolerance must be finite and ≥ 0; got $tolerance"))
    x = log10(Float64(n))
    xs = Float64[log10(Float64(r.n)) for r in rows]
    in_range = xs[1] ≤ x ≤ xs[end]
    i = searchsortedlast(xs, x)::Int
    lo, hi, w = if i == 0
        1, 1, 0.0
    elseif i == length(xs) || xs[i] == x
        i, i, 0.0
    else
        i, i + 1, (x - xs[i]) / (xs[i + 1] - xs[i])
    end
    interp(f) = lo == hi ? f(rows[lo]) : (1 - w) * f(rows[lo]) + w * f(rows[hi])
    truth = interp(r -> r.truth)
    bias = interp(r -> r.bias)
    scatter_abs = interp(r -> r.scatter)
    reported = interp(r -> r.reported_error)
    relative = isfinite(truth) && truth != 0
    scale = relative ? abs(truth) : 1.0
    error_ratio =
        (isnan(reported) || isnan(scatter_abs) || scatter_abs == 0) ? NaN : reported / scatter_abs
    rel_bias = bias / scale
    ok = relative ? (isfinite(rel_bias) && abs(rel_bias) ≤ tolerance) : true
    return CalibratedEstimate(rel_bias, scatter_abs / scale, error_ratio, relative, in_range, ok)
end

"""
    calibrated(cal::CalibrationTable, model, estimator, quantity, n; tolerance = 0.05)
        -> CalibratedEstimate

Bias and scatter of an estimator at particle number `n`, from the
calibration rows of `model`, `estimator` and `quantity`. With `x = log10(n)`
and the grid `xᵢ = log10(Nᵢ)` of the rows, `bias`, `scatter`,
`reported_error` and `truth` are interpolated linearly between the
bracketing grid points; outside the grid the end row is taken and
`in_range` is false. When the interpolated truth is finite and nonzero the
bias and the scatter are divided by |truth| (`relative`), otherwise they stay
absolute. A NaN entry of a row used by the interpolation propagates. `ok`
holds when a relative bias is finite and at most `tolerance` in modulus, and
always for an absolute one.

# Arguments
- `cal`: table of `load_calibration`
- `model`, `estimator`, `quantity`: the triple of the rows
- `n`: particle number, finite and > 0
- `tolerance`: largest admissible |relative bias|, ≥ 0

# Returns
A `CalibratedEstimate`. Throws an `ArgumentError` when the table has
no row for the triple, or for an invalid `n` or `tolerance`.

# Example
```julia
c = calibrated(load_calibration(), "plummer", "core_radius", "r_core_number", 50_000)
c.ok
```
"""
function calibrated(
    cal::CalibrationTable,
    model::AbstractString,
    estimator::AbstractString,
    quantity::AbstractString,
    n::Real;
    tolerance::Real = 0.05,
)::CalibratedEstimate
    rows = _validity_calibration_rows(cal, model, estimator, quantity)
    return _validity_interpolate(rows, n, tolerance)
end

"""
    _validity_estimates(tracked, n, tolerance) -> Vector{CalibratedEstimate}

One estimate per entry of `VALIDITY_TRACKED` at `n` bound stars from the
pre-selected rows `tracked`; the empty estimate for every entry when
`tracked` is `nothing` or `n == 0`.
"""
function _validity_estimates(
    tracked::Union{Nothing,Vector{Vector{CalibrationRow}}},
    n::Integer,
    tolerance::Real,
)::Vector{CalibratedEstimate}
    (tracked === nothing || n ≤ 0) && return fill(_VALIDITY_NO_ESTIMATE, length(VALIDITY_TRACKED))
    return CalibratedEstimate[_validity_interpolate(rows, n, tolerance) for rows in tracked]
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

"""Largest finite element of `v`; NaN when there is none."""
function _validity_finite_max(v)::Float64
    m = -Inf
    found = false
    for x in v
        isfinite(x) || continue
        found = true
        m = max(m, Float64(x))
    end
    return found ? m : NaN
end

"""Median of the finite elements of `v`; NaN when there is none."""
function _validity_finite_median(v)::Float64
    f = Float64[x for x in v if isfinite(x)]
    return isempty(f) ? NaN : _median(f)
end

"""Number with three significant digits, for the `detail` strings."""
_validity_g(x::Real) = @sprintf("%.3g", x)

"""Number in the form `2.1e-3` (one decimal, exponent without padding), for the `reasons`."""
_validity_sci(x::Real) = replace(@sprintf("%.1e", x), r"e\+?(-?)0*(\d+)" => s"e\1\2")

"""
    _validity_holds_from(pred, n) -> Union{Nothing,Int}

First index `k` of `1:n` from which `pred(j)` holds for every `j ≥ k`;
`nothing` when `pred(n)` fails or `n == 0`.
"""
function _validity_holds_from(pred, n::Integer)::Union{Nothing,Int}
    k = n + 1
    while k > 1 && pred(k - 1)
        k -= 1
    end
    return k ≤ n ? k : nothing
end

"""
    _validity_rate_test(t, y, window, tolerance) -> Bool

Stationarity of the series `y` at times `t` over the time `window`: with
`m = max(1, round(Int, window / median(diff(t))))`, true when
`|y[j + m] − y[j]| ≤ tolerance` for every `j` with `j + m ≤ length(y)`;
false when `length(y) ≤ m` or `m` is undefined (fewer than two times, a
non-finite window or a non-positive median step).
"""
function _validity_rate_test(
    t::AbstractVector{Float64},
    y::AbstractVector{Float64},
    window::Float64,
    tolerance::Float64,
)::Bool
    length(t) ≥ 2 || return false
    dt = _validity_finite_median(diff(t))
    (isfinite(window) && isfinite(dt) && dt > 0) || return false
    ratio = window / dt
    ratio ≥ length(y) && return false
    m = max(1, round(Int, ratio))
    length(y) ≤ m && return false
    for j in 1:(length(y) - m)
        abs(y[j + m] - y[j]) ≤ tolerance || return false
    end
    return true
end

"""
    _validity_lagr_stationary(lagr, t0, window, tolerance) -> Bool

Rate test of `ln r_10`, the Lagrangian radius of mass fraction 0.1 in
`lagr.7`, on the times `t ≥ t0`; true when `lagr` is `nothing` or has no 0.1
fraction (the test does not apply).
"""
function _validity_lagr_stationary(
    lagr::Union{Nothing,LagrangianData},
    t0::Float64,
    window::Float64,
    tolerance::Float64,
)::Bool
    lagr === nothing && return true
    f = findfirst(x -> isapprox(x, 0.1; atol = 1e-12), lagr.mass_fractions)
    (f === nothing || f > size(lagr.radii, 1)) && return true
    sel = findall(≥(t0), lagr.time)
    y = Float64[r > 0 ? log(r) : NaN for r in lagr.radii[f, sel]]
    return _validity_rate_test(lagr.time[sel], y, window, tolerance)
end

"""
    _validity_lagr_drift(lagr, t0, t_cr) -> Float64

Median of `|Δ ln r_10| / (Δt / t_cr)` over consecutive `lagr.7` records at
`t ≥ t0`: the relative change of the 10 % Lagrangian radius per crossing time.
NaN without `lagr`, without a 0.1 fraction, with fewer than two records or a
non-finite `t_cr`.
"""
function _validity_lagr_drift(
    lagr::Union{Nothing,LagrangianData},
    t0::Float64,
    t_cr::Float64,
)::Float64
    (lagr === nothing || !(isfinite(t_cr) && t_cr > 0)) && return NaN
    f = findfirst(x -> isapprox(x, 0.1; atol = 1e-12), lagr.mass_fractions)
    (f === nothing || f > size(lagr.radii, 1)) && return NaN
    sel = findall(≥(t0), lagr.time)
    rates = Float64[]
    for q in 2:length(sel)
        i, j = sel[q - 1], sel[q]
        dt = lagr.time[j] - lagr.time[i]
        (dt > 0 && lagr.radii[f, i] > 0 && lagr.radii[f, j] > 0) || continue
        push!(rates, abs(log(lagr.radii[f, j] / lagr.radii[f, i])) / (dt / t_cr))
    end
    return _validity_finite_median(rates)
end

"""Calibration model of the run: `"plummer"` for equal masses (to 1e-12 relative), `"plummer_kroupa"` otherwise, `""` without stars."""
function _validity_model(snaps::AbstractVector{Snapshot})::String
    isempty(snaps) && return ""
    m = snaps[1].mass
    isempty(m) && return ""
    m0 = Float64(m[1])
    return all(x -> abs(Float64(x) - m0) ≤ 1e-12 * abs(m0), m) ? "plummer" : "plummer_kroupa"
end

"""
    _validity_initial_relaxation_time(snaps, ints, settings) -> Float64

Half-mass relaxation time of the bound stars of the first snapshot for the
first Coulomb factor of `settings.relaxation_gamma` (`relaxation_times`
about their shrinking-sphere centre); NaN when it cannot be computed.
"""
function _validity_initial_relaxation_time(
    snaps::AbstractVector{Snapshot},
    ints::AbstractVector{SnapshotIntegrals},
    settings::PostprocessConfig,
)::Float64
    (isempty(snaps) || isempty(ints) || isempty(settings.relaxation_gamma)) && return NaN
    stars = ints[1].bound.stars
    n = length(stars)
    γ = (settings.relaxation_gamma[1], settings.relaxation_gamma[end])
    (n ≥ 1 && all(g -> g > 0 && g * n > 1, γ)) || return NaN
    s = snaps[1]
    pos = Float64.(s.pos[:, stars])
    vel = Float64.(s.vel[:, stars])
    mass = Float64.(s.mass[stars])
    sum(mass) > 0 || return NaN
    centre = shrinking_sphere_centre(pos, mass).centre
    return relaxation_times(pos, vel, mass, centre; gamma = γ).t_rh[1]
end

"""
    _validity_dense_crossing_time(pb, tags, t, fallback) -> Float64

Crossing time `2 r_h / σ` of the denser progenitor of the pair `tags` at
the latest progenitor-budget time ≤ `t`, with `σ² = 0.4 m_bound / r_h`
(Spitzer 1987, `r_h ≈ 0.4 G M / ⟨v²⟩`) and the density measured by
`m_bound / r_h³`; `fallback` when the budget lacks either tag there or its
`m_bound`, `r_h` are not positive and finite.
"""
function _validity_dense_crossing_time(
    pb::AbstractVector{ProgenitorBudgetRow},
    tags::Tuple{Int,Int},
    t::Float64,
    fallback::Float64,
)::Float64
    best_density = -Inf
    best_tcr = fallback
    for tag in tags
        idx = 0
        for (k, row) in enumerate(pb)
            (row.tag == tag && row.time_nb ≤ t) || continue
            (idx == 0 || row.time_nb ≥ pb[idx].time_nb) && (idx = k)
        end
        idx == 0 && return fallback
        r = pb[idx]
        (isfinite(r.r_h) && r.r_h > 0 && isfinite(r.m_bound) && r.m_bound > 0) || return fallback
        density = r.m_bound / r.r_h^3
        if density > best_density
            best_density = density
            σ = sqrt(0.4 * r.m_bound / r.r_h)
            best_tcr = 2 * r.r_h / σ
        end
    end
    return best_tcr
end

# ---------------------------------------------------------------------------
# Sections
# ---------------------------------------------------------------------------

"""The `energy` section of `validity_record`."""
function _validity_energy(
    diag::Union{Nothing,DiagnosticsData},
    ev::Union{Nothing,EnergyEventAnalysis},
    eb::AbstractVector{EnergyBudgetRow},
    amendments::Vector{String},
    settings::PostprocessConfig,
)::Dict{String,Any}
    adj = diag === nothing ? AdjustRecord[] : diag.adjust
    mode = diag === nothing ? :none : tidal_mode(diag)
    valid = engine_energy_valid(mode; tidal_energy_amended = "tidal-energy-total" in amendments)
    available = !isempty(adj)
    de = Float64[abs(a.de_rel) for a in adj[2:end] if isfinite(a.de_rel)]
    n_above = count(>(settings.validity_interval_de), de)
    n_flagged = ev === nothing ? 0 : count(iv -> iv.flagged, ev.intervals)
    # Relative error with the change of the pairs' binding energy removed
    excl = Float64[]
    if ev !== nothing
        for iv in ev.intervals
            usable =
                isfinite(iv.de_rel) &&
                iv.de_rel != 0 &&
                isfinite(iv.de_abs) &&
                iv.de_abs != 0 &&
                isfinite(iv.d_ebin)
            usable && push!(excl, abs((iv.de_abs - iv.d_ebin) * iv.de_rel / iv.de_abs))
        end
    end
    excl_pairs = _validity_finite_max(excl)
    e0 = available ? adj[1].e_tot : NaN
    detot = available ? adj[end].detot : NaN
    cumulative = abs(detot) / abs(e0)
    verdict, reasons = if !valid
        "not_applicable", String[]
    elseif !available
        "unavailable", String[]
    else
        msgs = String[]
        n_above > 0 && push!(
            msgs,
            "$n_above interval$(n_above == 1 ? "" : "s") above " *
            _validity_sci(settings.validity_interval_de),
        )
        cumulative > settings.validity_cumulative_de && push!(
            msgs,
            "cumulative |DETOT/E0| = $(_validity_sci(cumulative)) > " *
            _validity_sci(settings.validity_cumulative_de),
        )
        (isempty(msgs) ? "ok" : "flag"), msgs
    end
    return Dict{String,Any}(
        "available" => available,
        "engine_energy_valid" => valid,
        "n_intervals" => max(length(adj) - 1, 0),
        "max_abs_de_rel" => _validity_finite_max(de),
        "median_abs_de_rel" => _validity_finite_median(de),
        "n_above_threshold" => n_above,
        "n_flagged" => n_flagged,
        "max_abs_de_rel_excl_pairs" => excl_pairs,
        "e0" => e0,
        "detot" => detot,
        "cumulative_rel" => cumulative,
        "budget_max_abs_residual_rel" => _validity_finite_max(abs(r.residual_rel) for r in eb),
        "budget_final_residual_rel" => isempty(eb) ? NaN : eb[end].residual_rel,
        "verdict" => verdict,
        "reasons" => reasons,
    )
end

"""The `mass` section of `validity_record`."""
function _validity_mass(
    mb::AbstractVector{MassBudgetRow},
    settings::PostprocessConfig,
)::Dict{String,Any}
    return Dict{String,Any}(
        "available" => !isempty(mb),
        "n_rows" => length(mb),
        "max_abs_residual_rel" => _validity_finite_max(abs(r.residual_rel) for r in mb),
        "final_residual_rel" => isempty(mb) ? NaN : mb[end].residual_rel,
        "n_flagged" => count(r -> abs(r.residual_rel) > settings.validity_mass_tolerance, mb),
    )
end

"""The `angular_momentum` section of `validity_record`: `|residual| / l_scale` per snapshot."""
function _validity_angular_momentum(am::Union{Nothing,AngularMomentumBudget})::Dict{String,Any}
    if am === nothing || isempty(am.l_scale)
        return Dict{String,Any}(
            "available" => false,
            "max_residual_over_scale" => NaN,
            "final_residual_over_scale" => NaN,
        )
    end
    ratio = Float64[
        (isfinite(am.l_scale[j]) && am.l_scale[j] != 0) ?
        sqrt(am.residual[1, j]^2 + am.residual[2, j]^2 + am.residual[3, j]^2) / am.l_scale[j] :
        NaN for j in eachindex(am.l_scale)
    ]
    return Dict{String,Any}(
        "available" => true,
        "max_residual_over_scale" => _validity_finite_max(ratio),
        "final_residual_over_scale" => ratio[end],
    )
end

"""The `bound_set` section of `validity_record`."""
function _validity_bound_section(ints::AbstractVector{SnapshotIntegrals})::Dict{String,Any}
    return Dict{String,Any}(
        "n_snapshots" => length(ints),
        "n_converged" => count(s -> s.bound.status === :converged, ints),
        "n_pass_limit" => count(s -> s.bound.status === :pass_limit, ints),
        "n_dissolved" => count(s -> s.bound.status === :dissolved, ints),
        "max_passes" => isempty(ints) ? 0 : maximum(s -> s.bound.passes, ints),
    )
end

"""
    _validity_capacity(diag, build_info) -> Dict{String,Any}

Fixed-size engine arrays the run came close to. The SIMD irregular-force library is opened for
N + NBIN0 − NPAIRS + 100 bodies at each start or restart (`intgrt.F`) and indexed without a bounds check,
so KS pairs beyond NBIN0 + 100 overrun it (500 more slots at a cold start; the test ignores them); the
merger table of stable hierarchies holds MMAX − 1 entries (`impact.f`). From the adjustment records:
`npairs_max`, `nmerge_max`; `nbin0` from the log (`-1` unknown); `simd_margin` = 100,
`simd_headroom` = 100 − (npairs_max − nbin0) (NaN when nbin0 is unknown); `mmax` from
`build_info["engine_limits"]["MMAX"]` (`-1` unknown), `merge_table_headroom` = mmax − 1 − nmerge_max (`-1`
when mmax is unknown); `verdict` `"ok"`, `"flag"` (a headroom below zero for the SIMD library, or at or
below zero for the merger table) or `"unavailable"` (no records); `reasons` names each flag.
"""
function _validity_capacity(
    diag::Union{Nothing,DiagnosticsData},
    build_info::AbstractDict,
)::Dict{String,Any}
    adj = diag === nothing ? AdjustRecord[] : diag.adjust
    nbin0 = diag === nothing ? -1 : diag.nbin0
    npairs_max = isempty(adj) ? 0 : maximum(a.npairs for a in adj)
    nmerge_max = isempty(adj) ? 0 : maximum(a.nmerge for a in adj)
    limits = get(build_info, "engine_limits", nothing)
    mmax = limits isa AbstractDict ? Int(get(limits, "MMAX", -1)) : -1
    mmax ≤ 0 && (mmax = -1)
    simd_headroom = nbin0 ≥ 0 ? Float64(100 - (npairs_max - nbin0)) : NaN
    merge_table_headroom = mmax > 0 ? mmax - 1 - nmerge_max : -1
    reasons = String[]
    if nbin0 ≥ 0 && simd_headroom < 0
        excess = npairs_max - nbin0
        push!(reasons, "KS pairs reached NBIN0 + $(excess) > NBIN0 + 100 (SIMD library overrun)")
    end
    if mmax > 0 && merge_table_headroom ≤ 0
        push!(reasons, "merger table full: NMERGE = $(nmerge_max) of MMAX − 1 = $(mmax - 1)")
    end
    verdict = isempty(adj) ? "unavailable" : (isempty(reasons) ? "ok" : "flag")
    return Dict{String,Any}(
        "available" => !isempty(adj),
        "nbin0" => nbin0,
        "npairs_max" => npairs_max,
        "simd_margin" => 100,
        "simd_headroom" => simd_headroom,
        "nmerge_max" => nmerge_max,
        "mmax" => mmax,
        "merge_table_headroom" => merge_table_headroom,
        "verdict" => verdict,
        "reasons" => reasons,
    )
end

"""
The `masks` section of `validity_record`: the masking rules M1–M7
with the thresholds of `settings`. M6 restricts the equilibrium estimators
of a tidal run to 0.5 r_J (Küpper et al. 2010, DOI
10.1111/j.1365-2966.2009.16173.x).
"""
function _validity_masks(settings::PostprocessConfig)::Dict{String,Any}
    return Dict{String,Any}(
        "M1" => "profile radii r ≥ max(3σ_x, r enclosing 100 bodies)",
        "M2_lagr_n_min" => settings.lagr_n_min,
        "M3_n_min_density" => settings.profile_n_min_density,
        "M3_n_min_sigma" => settings.profile_n_min_sigma,
        "M3_n_min_beta" => settings.profile_n_min_beta,
        "M3_rotation_n_min" => settings.rotation_n_min,
        "M4_shape_n_min_shell" => settings.shape_n_min_shell,
        "M4_shape_n_min_enclosed" => settings.shape_n_min_enclosed,
        "M5_msr_scan" => Float64.(settings.msr_scan),
        "M5_msr_n_random" => settings.msr_n_random,
        "M6" => "tidal runs: equilibrium estimators inside 0.5 r_J; isolated runs: bound set inside r_99%",
        "M7_settle_tcr" => settings.validity_settle_tcr,
        "M7_pericentre_tcr" => settings.validity_pericentre_tcr,
    )
end

# ---------------------------------------------------------------------------
# Record
# ---------------------------------------------------------------------------

"""
    validity_record(results::AbstractDict{Symbol,Any}, settings::PostprocessConfig;
                    calibration = load_calibration(), build_info = Dict{String,Any}(),
                    run_info = Dict{String,Any}()) -> ValidityRecord

Validity record of a post-processed run from the dictionary of `postprocess`
(any subset of `:snapshots`, `:diagnostics`, `:integrals`, `:energy_events`,
`:mass_budget`, `:energy_budget`, `:angular_momentum`, `:remnant`,
`:progenitor_budget`, `:lagr`; a missing entry leaves its checks unavailable
or NaN) and the thresholds `validity_*` of `settings`.

The crossing time of a snapshot is the standard crossing time of Heggie &
Mathieu (1986), `t_cr = M^{5/2} / (−2E)^{3/2}` in `G = 1` units, of its bound
set (`crossing_time`; 2√2 for `M = 1`, `E = −1/4`). Around a pericentre the
time scale is the crossing time `2 r_h / σ` of the denser progenitor, with
`σ² = 0.4 m / r_h` from the relation `r_h ≈ 0.4 G M / ⟨v²⟩` (Spitzer 1987).

The sections:
- `header`: run id (`run_info["run"]["id"]`), creation time, package
  version, `backend_commit` and `source_amendments` of `build_info`, the
  calibration version and model (`"plummer"` for equal masses,
  `"plummer_kroupa"` otherwise), tidal mode, `KZ(14)`, initial star count,
  snapshot count and time span
- `thresholds`: the validity settings, keys without the `validity_` prefix
- `energy`: per-interval and cumulative energy errors of the ADJUST records,
  the flagged intervals, the error without the pair-energy channel, the
  energy-budget residuals and a `verdict` (`"ok"`, `"flag"`, `"unavailable"`
  or `"not_applicable"` when the engine's DE is not an error measure) with
  its `reasons`
- `mass`, `angular_momentum`: budget residuals
- `bound_set`: convergence of the bound-set iteration
- `capacity`: the SIMD-library and merger-table headroom (`_validity_capacity`)
- `cadence`: snapshot spacing against the crossing and the initial
  half-mass relaxation times
- `stages`: end of violent relaxation, pericentres, coalescence, and the
  times of virialisation (2T/|W| within `virial_tolerance` from then on; the
  drift of `ln r_10` per crossing time and its stationarity are reported
  beside it, since stellar mass loss and relaxation move the inner radii of a
  virialised cluster), kinematic relaxation (λ_R stationary) and mixing
  (|a₊₂| ≤ 2σ for every progenitor from then on), each `"reached"`,
  `"not_reached"`, `"undetermined"` (or `"not_applicable"` for mixing with
  fewer than two progenitors), also in initial half-mass relaxation times
- `calibration`: bias and scatter of the `VALIDITY_TRACKED` estimators at
  the bound star count of the first snapshot
- `masks`: the masking rules as applied (M6 after Küpper et al. 2010, DOI
  10.1111/j.1365-2966.2009.16173.x)

The window epochs (`:violent_relaxation`, `:pericentre`,
`:coalescence_phase`) exclude snapshots from the equilibrium estimators;
the other epochs mark energy, mass and bound-set problems.

# Arguments
- `results`: the dictionary returned by `postprocess`
- `settings`: the `[postprocess]` configuration
- `calibration`: estimator calibration; `nothing` leaves the estimates NaN
- `build_info`: the engine build record (`BUILD_INFO.toml`)
- `run_info`: the run record (`RUN_INFO.toml`)

# Returns
A `ValidityRecord`. Throws a `DimensionMismatch` when snapshots and
integrals are both given but differ in number.

# Example
```julia
results = postprocess(cfg; run_dir)
rec = validity_record(results, cfg.postprocess)
write_validity(joinpath(run_dir, "validity.toml"), rec)
```
"""
function validity_record(
    results::AbstractDict{Symbol,Any},
    settings::PostprocessConfig;
    calibration::Union{Nothing,CalibrationTable} = load_calibration(),
    build_info::AbstractDict = Dict{String,Any}(),
    run_info::AbstractDict = Dict{String,Any}(),
)::ValidityRecord
    snaps = get(results, :snapshots, Snapshot[])::AbstractVector{Snapshot}
    diag = get(results, :diagnostics, nothing)::Union{Nothing,DiagnosticsData}
    ints = get(results, :integrals, SnapshotIntegrals[])::AbstractVector{SnapshotIntegrals}
    ev = get(results, :energy_events, nothing)::Union{Nothing,EnergyEventAnalysis}
    mb = get(results, :mass_budget, MassBudgetRow[])::AbstractVector{MassBudgetRow}
    eb = get(results, :energy_budget, EnergyBudgetRow[])::AbstractVector{EnergyBudgetRow}
    am = get(results, :angular_momentum, nothing)::Union{Nothing,AngularMomentumBudget}
    remn = get(results, :remnant, nothing)::Union{Nothing,RemnantDiagnostics}
    pb =
        get(results, :progenitor_budget, ProgenitorBudgetRow[])::AbstractVector{ProgenitorBudgetRow}
    lagr = get(results, :lagr, nothing)::Union{Nothing,LagrangianData}
    (isempty(snaps) || isempty(ints) || length(snaps) == length(ints)) || throw(
        DimensionMismatch(
            "validity_record: $(length(snaps)) snapshots but $(length(ints)) snapshot integrals",
        ),
    )

    # Header
    amendments = String[string(a) for a in get(build_info, "source_amendments", String[])]
    run_table = get(run_info, "run", nothing)
    run_id = run_table isa AbstractDict ? string(get(run_table, "id", "")) : ""
    model = _validity_model(snaps)
    mode = diag === nothing ? :none : tidal_mode(diag)
    header = Dict{String,Any}(
        "run_id" => run_id,
        "created" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        "package_version" => string(pkgversion(@__MODULE__)),
        "backend_commit" => string(get(build_info, "backend_commit", "")),
        "source_amendments" => amendments,
        "calibration_version" => calibration === nothing ? -1 : calibration.version,
        "calibration_model" => model,
        "tidal_mode" => String(mode),
        "kz14" => (diag !== nothing && length(diag.kz) ≥ 14) ? diag.kz[14] : -1,
        "n_initial" => isempty(snaps) ? 0 : nparticles(snaps[1]),
        "n_snapshots" => length(snaps),
        "t_start_nb" => isempty(snaps) ? NaN : time_nb(snaps[1].header),
        "t_end_nb" => isempty(snaps) ? NaN : time_nb(snaps[end].header),
    )
    thresholds = Dict{String,Any}(
        "validity" => settings.validity,
        "interval_de" => settings.validity_interval_de,
        "cumulative_de" => settings.validity_cumulative_de,
        "mass_tolerance" => settings.validity_mass_tolerance,
        "settle_tcr" => settings.validity_settle_tcr,
        "pericentre_tcr" => settings.validity_pericentre_tcr,
        "window_tcr" => settings.validity_window_tcr,
        "virial_tolerance" => settings.validity_virial_tolerance,
        "rate_tolerance" => settings.validity_rate_tolerance,
        "bias_tolerance" => settings.validity_bias_tolerance,
    )

    # Snapshot series: the bound set of every snapshot integral
    times =
        isempty(snaps) ? Float64[s.time_nb for s in ints] :
        Float64[time_nb(s.header) for s in snaps]
    t_start = isempty(times) ? NaN : times[1]
    t_end = isempty(times) ? NaN : times[end]
    nrow = length(ints)
    t_row = Float64[s.time_nb for s in ints]
    n_bound = zeros(Int, nrow)
    n_eff = fill(NaN, nrow)
    t_cr = fill(NaN, nrow)
    q_virial = fill(NaN, nrow)
    converged = falses(nrow)
    for (k, s) in enumerate(ints)
        b = s.bound
        n_bound[k] = length(b.stars)
        if !isempty(snaps)
            m = snaps[k].mass
            Σm = 0.0
            Σm² = 0.0
            for i in b.stars
                mi = Float64(m[i])
                Σm += mi
                Σm² += mi^2
            end
            n_eff[k] = Σm² > 0 ? Σm^2 / Σm² : NaN
        end
        t_cr[k] = s.bound_mass > 0 ? crossing_time(s.bound_mass, b.kinetic + b.potential) : NaN
        q_virial[k] = b.potential != 0 ? 2 * b.kinetic / abs(b.potential) : NaN
        converged[k] = b.status === :converged
    end
    t_cr_initial = nrow ≥ 1 ? t_cr[1] : NaN
    t_cr_median = _validity_finite_median(t_cr)

    # Cadence
    dts = length(times) ≥ 2 ? diff(times) : Float64[]
    dt_max = _validity_finite_max(dts)
    t_rh_initial = _validity_initial_relaxation_time(snaps, ints, settings)
    cadence = Dict{String,Any}(
        "dt_snapshot_median" => _validity_finite_median(dts),
        "dt_snapshot_max" => dt_max,
        "t_cr_initial" => t_cr_initial,
        "t_cr_median" => t_cr_median,
        "dt_over_tcr_max" =>
            _validity_finite_max((t_row[k] - t_row[k - 1]) / t_cr[k - 1] for k in 2:nrow),
        "t_rh_initial" => t_rh_initial,
        "core_collapse_resolved" => dt_max ≤ 0.1 * t_rh_initial,
        "coalescence_time_resolution" =>
            (remn !== nothing && !isempty(remn.coalescence_pairs)) ? dt_max : NaN,
    )

    # Epochs
    epochs = ValidityEpoch[]
    t_violent_end = t_start + settings.validity_settle_tcr * t_cr_initial
    isfinite(t_violent_end) && push!(
        epochs,
        ValidityEpoch(
            t_start,
            t_violent_end,
            :violent_relaxation,
            "first $(_validity_g(settings.validity_settle_tcr)) t_cr",
        ),
    )
    t_coal = remn === nothing ? NaN : remn.coalescence_time
    pericentres = Tuple{Float64,Float64}[]   # (time, half-width)
    if remn !== nothing
        S = remn.coalescence_separation
        for p in axes(S, 1)
            i, j = remn.coalescence_pairs[p]
            for k in 2:(size(S, 2) - 1)
                a, s, c = S[p, k - 1], S[p, k], S[p, k + 1]
                (isfinite(a) && isfinite(s) && isfinite(c)) || continue
                (s < a && s ≤ c) || continue
                tk = remn.time[k]
                (isnan(t_coal) || tk < t_coal) || continue
                w =
                    settings.validity_pericentre_tcr *
                    _validity_dense_crossing_time(pb, (i, j), tk, t_cr_initial)
                push!(pericentres, (tk, w))
                push!(
                    epochs,
                    ValidityEpoch(
                        tk - w,
                        tk + w,
                        :pericentre,
                        "clusters $(i)–$(j), separation $(_validity_g(s))",
                    ),
                )
            end
        end
    end
    if !isempty(pericentres)
        tp, wp = pericentres[argmin(first.(pericentres))]
        if isfinite(t_coal)
            push!(epochs, ValidityEpoch(tp - wp, t_coal, :coalescence_phase, "until coalescence"))
        else
            push!(
                epochs,
                ValidityEpoch(tp - wp, t_end, :coalescence_phase, "not coalesced by the end"),
            )
        end
    elseif isfinite(t_coal)
        # Coalescence without a resolved pericentre: the whole approach is out of equilibrium
        push!(
            epochs,
            ValidityEpoch(
                t_start,
                t_coal,
                :coalescence_phase,
                "approach without a pericentre resolved at this cadence",
            ),
        )
    end
    cadence["pericentre_resolved"] =
        (remn === nothing || isempty(remn.coalescence_pairs)) ? "not_applicable" :
        (isempty(pericentres) ? "no" : "yes")
    if ev !== nothing
        for iv in ev.intervals
            iv.flagged && push!(
                epochs,
                ValidityEpoch(
                    iv.t_start,
                    iv.t_end,
                    :energy_interval,
                    "score $(_validity_g(iv.score)), $(iv.n_events) events",
                ),
            )
        end
    end
    # The engine's DE is no error measure where the energy section's
    # `engine_energy_valid` is false (verdict `not_applicable`): no threshold epochs.
    de_valid = engine_energy_valid(mode; tidal_energy_amended = "tidal-energy-total" in amendments)
    adj = diag === nothing || !de_valid ? AdjustRecord[] : diag.adjust
    for k in 2:length(adj)
        de_k = adj[k].de_rel
        (isfinite(de_k) && abs(de_k) > settings.validity_interval_de) || continue
        push!(
            epochs,
            ValidityEpoch(
                adj[k - 1].time_nb,
                adj[k].time_nb,
                :energy_threshold,
                "DE = $(_validity_g(de_k))",
            ),
        )
    end
    for j in 2:length(mb)
        r = mb[j].residual_rel
        abs(r) > settings.validity_mass_tolerance || continue
        push!(
            epochs,
            ValidityEpoch(
                mb[j - 1].time_nb,
                mb[j].time_nb,
                :mass_budget,
                "relative residual $(_validity_g(r))",
            ),
        )
    end
    for s in ints
        b = s.bound
        b.status === :converged && continue
        push!(
            epochs,
            ValidityEpoch(
                s.time_nb,
                s.time_nb,
                :bound_set,
                "$(String(b.status)), $(b.passes) passes",
            ),
        )
    end
    sort!(epochs; by = e -> e.t_start)
    windows = filter(e -> e.kind in _VALIDITY_WINDOW_KINDS, epochs)
    t_w_end = _validity_finite_max(e.t_end for e in windows)
    T_w = isnan(t_w_end) ? t_start : t_w_end

    # Stages
    window = settings.validity_window_tcr * t_cr_median
    rate_tol = settings.validity_rate_tolerance
    after = findall(≥(T_w), t_row)
    t_virialised = NaN
    virialised = "undetermined"
    if length(after) ≥ 3
        virialised = "not_reached"
        k_q = _validity_holds_from(
            j -> abs(q_virial[j] - 1) ≤ settings.validity_virial_tolerance,
            nrow,
        )
        if k_q !== nothing
            virialised = "reached"
            t_virialised = t_row[max(k_q, first(after))]
        end
    end
    # Secular evolution of the inner Lagrangian radius from virialisation on
    # (stellar mass loss, relaxation): reported beside the virial state, not gating it
    t_drift = isnan(t_virialised) ? T_w : t_virialised
    r10_drift = _validity_lagr_drift(lagr, t_drift, t_cr_median)
    r10_stationary =
        isnan(r10_drift) ? "not_applicable" :
        (_validity_lagr_stationary(lagr, t_drift, window, rate_tol) ? "yes" : "no")
    t_kinematic = NaN
    kinematic = "undetermined"
    if remn !== nothing
        sel = findall(≥(T_w), remn.time)
        if length(sel) ≥ 3
            kinematic = "not_reached"
            floor_median = _validity_finite_median(remn.lambda_r_floor)
            tol = isfinite(floor_median) ? max(rate_tol, 2 * floor_median) : rate_tol
            t = remn.time[sel]
            y = remn.lambda_r[sel]
            for k in eachindex(t)
                if _validity_rate_test(view(t, k:length(t)), view(y, k:length(y)), window, tol)
                    kinematic = "reached"
                    t_kinematic = t[k]
                    break
                end
            end
        end
    end
    t_mixed = NaN
    mixed = "not_applicable"
    if length(unique(r.tag for r in pb)) ≥ 2
        if any(r -> isfinite(r.aplus2), pb)
            ts = sort!(unique(Float64[r.time_nb for r in pb]))
            mixed_at = trues(length(ts))
            for r in pb
                isfinite(r.aplus2) || continue
                abs(r.aplus2) ≤ 2 * r.aplus2_err && continue
                mixed_at[searchsortedfirst(ts, r.time_nb)] = false
            end
            k_m = _validity_holds_from(i -> mixed_at[i], length(ts))
            if k_m === nothing
                mixed = "not_reached"
            else
                mixed = "reached"
                t_mixed = ts[k_m]
            end
        else
            mixed = "undetermined"
        end
    end
    peri_times = sort!(Float64[first(p) for p in pericentres])
    peri_unique = Float64[]
    for x in peri_times
        (isempty(peri_unique) || abs(x - peri_unique[end]) > 1e-12) && push!(peri_unique, x)
    end
    stages = Dict{String,Any}(
        "t_violent_end" => t_violent_end,
        "pericentres" => peri_unique,
        "coalescence_time_nb" => t_coal,
        "t_virialised" => t_virialised,
        "virialised" => virialised,
        "r10_drift_per_tcr" => r10_drift,
        "r10_stationary" => r10_stationary,
        "t_kinematic" => t_kinematic,
        "kinematic" => kinematic,
        "t_mixed" => t_mixed,
        "mixed" => mixed,
        "nh_virialised" => t_virialised / t_rh_initial,
        "nh_kinematic" => t_kinematic / t_rh_initial,
        "nh_mixed" => t_mixed / t_rh_initial,
    )

    # Calibration
    tracked =
        (calibration === nothing || isempty(model)) ? nothing :
        Vector{CalibrationRow}[
            _validity_calibration_rows(calibration, model, e, q) for (e, q) in VALIDITY_TRACKED
        ]
    btol = settings.validity_bias_tolerance
    calibration_section = if calibration === nothing
        Dict{String,Any}("version" => -1)
    else
        n_ref = nrow ≥ 1 ? n_bound[1] : 0
        cal_table = Dict{String,Any}(
            "version" => calibration.version,
            "model" => model,
            "n_reference" => n_ref,
        )
        for ((e, q), c) in zip(VALIDITY_TRACKED, _validity_estimates(tracked, n_ref, btol))
            cal_table["$(e)__$(q)"] = Dict{String,Any}(
                "bias" => c.bias,
                "scatter" => c.scatter,
                "error_ratio" => c.error_ratio,
                "relative" => c.relative,
                "in_range" => c.in_range,
                "ok" => c.ok,
            )
        end
        cal_table
    end

    # Per-snapshot rows
    lagr_n_min = settings.lagr_n_min
    snapshots = Vector{ValiditySnapshot}(undef, nrow)
    for k in 1:nrow
        resolvable = n_eff[k] ≥ 2 * lagr_n_min
        in_window = any(e -> e.t_start ≤ t_row[k] ≤ e.t_end, windows)
        snapshots[k] = ValiditySnapshot(
            t_row[k],
            n_bound[k],
            n_eff[k],
            t_cr[k],
            q_virial[k],
            resolvable ? lagr_n_min / n_eff[k] : NaN,
            resolvable ? 1 - lagr_n_min / n_eff[k] : NaN,
            converged[k],
            in_window,
            !in_window && converged[k],
            _validity_estimates(tracked, n_bound[k], btol),
        )
    end

    sections = Dict{String,Dict{String,Any}}(
        "header" => header,
        "thresholds" => thresholds,
        "energy" => _validity_energy(diag, ev, eb, amendments, settings),
        "mass" => _validity_mass(mb, settings),
        "angular_momentum" => _validity_angular_momentum(am),
        "bound_set" => _validity_bound_section(ints),
        "capacity" => _validity_capacity(diag, build_info),
        "cadence" => cadence,
        "stages" => stages,
        "calibration" => calibration_section,
        "masks" => _validity_masks(settings),
    )
    return ValidityRecord(VALIDITY_SCHEMA_VERSION, sections, epochs, snapshots)
end

# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

"""
    write_validity(path, rec::ValidityRecord) -> String

Write the record as TOML: the top-level `schema_version`, one table per
section in the order `header`, `thresholds`, `energy`, `mass`,
`angular_momentum`, `bound_set`, `capacity`, `cadence`, `stages`, `calibration` (with a
nested `[calibration.<estimator>__<quantity>]` table per tracked
estimator) and `masks`, keys sorted within a table, then the epochs as
`[[epochs]]` with `t_start`, `t_end`, `kind` and `detail`. NaN is written as
the TOML float `nan`. An existing file is backed up, never overwritten.
Returns the path.
"""
function write_validity(path::AbstractString, rec::ValidityRecord)::String
    _backup_existing(path)
    open(path, "w") do io
        TOML.print(io, Dict{String,Any}("schema_version" => rec.schema_version))
        for name in _VALIDITY_SECTIONS
            haskey(rec.sections, name) || continue
            println(io)
            TOML.print(io, Dict{String,Any}(name => rec.sections[name]); sorted = true)
        end
        println(io)
        epochs = Dict{String,Any}[
            Dict{String,Any}(
                "t_start" => e.t_start,
                "t_end" => e.t_end,
                "kind" => String(e.kind),
                "detail" => e.detail,
            ) for e in rec.epochs
        ]
        TOML.print(io, Dict{String,Any}("epochs" => epochs); sorted = true)
    end
    return String(path)
end

"""
    write_validity_snapshots(path, rec::ValidityRecord) -> String

Write the per-snapshot rows of the record as CSV: `time_nb, n_bound, n_eff,
t_cr, q_virial, lagr_f_min, lagr_f_max, bound_converged, in_window,
equilibrium_ok`, then for every entry `<estimator>__<quantity>` of
`VALIDITY_TRACKED` the columns `_bias`, `_scatter` and `_ok`. Booleans as
`true`/`false`, numbers with `string` (NaN as `NaN`). An existing file is
backed up, never overwritten. Returns the path.
"""
function write_validity_snapshots(path::AbstractString, rec::ValidityRecord)::String
    _backup_existing(path)
    header = [
        "time_nb",
        "n_bound",
        "n_eff",
        "t_cr",
        "q_virial",
        "lagr_f_min",
        "lagr_f_max",
        "bound_converged",
        "in_window",
        "equilibrium_ok",
    ]
    for (e, q) in VALIDITY_TRACKED
        push!(header, "$(e)__$(q)_bias", "$(e)__$(q)_scatter", "$(e)__$(q)_ok")
    end
    open(path, "w") do io
        println(io, join(header, ","))
        for s in rec.snapshots
            cells = String[
                string(s.time_nb),
                string(s.n_bound),
                string(s.n_eff),
                string(s.t_cr),
                string(s.q_virial),
                string(s.lagr_f_min),
                string(s.lagr_f_max),
                string(s.bound_converged),
                string(s.in_window),
                string(s.equilibrium_ok),
            ]
            for c in s.calibration
                push!(cells, string(c.bias), string(c.scatter), string(c.ok))
            end
            println(io, join(cells, ","))
        end
    end
    return String(path)
end
