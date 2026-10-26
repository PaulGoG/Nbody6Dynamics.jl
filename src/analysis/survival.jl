# =============================================================================
# Survival statistics of binary lifetimes
# =============================================================================
#
# Event tables from the binary spans on two time axes (binary age, cluster time
# with delayed entry at formation), and the nonparametric estimators on them:
# Kaplan–Meier with the Greenwood variance, Nelson–Aalen, the Aalen–Johansen
# cumulative incidence of competing end causes, the Turnbull NPMLE for
# interval-censored ends, the log-rank two-sample test and Cox–Snell residuals
# for the check of a predicted hazard.

# ---------------------------------------------------------------------------
# Event tables
# ---------------------------------------------------------------------------

"""
    EventTable

Exit times, entry times and causes of a sample of subjects, the input of the
estimators of this file ([`kaplan_meier`](@ref), [`nelson_aalen`](@ref),
[`aalen_johansen`](@ref), [`logrank_test`](@ref), [`cox_snell_residuals`](@ref)).

| field | meaning |
|---|---|
| `time` | exit time: event or censoring |
| `entry` | entry time (delayed entry, left truncation); 0 when none |
| `cause` | 0 for a censored exit, k ≥ 1 for an event of cause k |
| `id` | label of the subject (the row index of its span in [`event_table`](@ref)) |
| `cause_labels` | `cause_labels[k]` names cause k |

A subject is at risk at t when `entry < t ≤ time`. The estimators require
finite times with `entry < time` and causes in `0:length(cause_labels)`, and
throw `ArgumentError` naming the first row that violates this.
"""
Base.@kwdef struct EventTable
    time::Vector{Float64}
    entry::Vector{Float64}
    cause::Vector{Int}
    id::Vector{Int}
    cause_labels::Vector{Symbol}
end

# End class of a binary span → cause label; `:censored` gives cause 0
const _SURVIVAL_DEFAULT_CAUSES = Dict{Symbol,Symbol}(
    :gw_merger => :gw_merger,
    :stellar_merger => :stellar_merger,
    :collision => :collision,
    :escape_bound => :escape_bound,
    :exchange => :exchange,
    :sn_disruption => :sn_disruption,
    :flyby => :flyby,
    :disruption => :disruption,
    :unresolved => :censored,
    :censored_run_end => :censored,
)

"""
    event_table(spans::AbstractVector{BinarySpan}; time_axis::Symbol = :age, t_end::Real,
                unit::Symbol = :myr, causes = _SURVIVAL_DEFAULT_CAUSES) -> EventTable

The lifetimes of the binary spans of [`binary_spans`](@ref) as an
[`EventTable`](@ref), on one of two time axes.

# Arguments
- `spans`: the spans; row i of the table carries `id = i` (its index in `spans`)
- `time_axis`: `:age` for the binary age τ = t − t_form (entry 0), or `:cluster`
  for the cluster time t with delayed entry at t_form
- `t_end`: end of the record, the right-censoring time; it must be ≥ every span end
- `unit`: `:myr` (fields `t_form_myr`, `t_end_myr`) or `:nb` (`t_form_nb`, `t_end_nb`)
- `causes`: `Dict{Symbol,Symbol}` from the `end_class` of a span to a cause label;
  a label `:censored` gives cause 0. The default maps every end class to itself
  and `:unresolved` and `:censored_run_end` to `:censored`

# Method
A right-censored span (open at the end of the record) exits at `t_end` with
cause 0, whatever its end class; every other span exits at its end with the
cause its end class maps to. On the age axis a left-censored span (formed before
the first record, `left_censored = true`) has an unknown age and is left out;
the count of spans dropped is logged at the `Info` level. On the cluster axis a
left-censored span enters at 0, every other span at its formation time, so the
risk set at t holds the binaries that exist at t (delayed entry; the
product-limit estimate is then conditional on survival to entry). The
`cause_labels` are the distinct labels other than `:censored` in order of first
appearance among the rows kept.

The time axes and the treatment of escape (a censoring for the binary lifetime,
a competing cause for the in-cluster lifetime: map `:escape_bound` to
`:censored` in `causes` for the former) follow the recommended method for
binary lifetimes in N-body clusters; the astronomical use of these estimators
goes back to Feigelson & Nelson (1985, DOI 10.1086/163225).

Throws `ArgumentError` for an unknown `time_axis` or `unit`, a non-finite
`t_end`, a span with a non-finite time, a span that ends after `t_end`, a span
with `t_end ≤ t_form` (the span builders never produce one) and an end class
missing from `causes`.

# Example
```julia
spans = binary_spans(log; nzero, nbin0, tstar_myr)
age = event_table(spans; t_end = t_last_myr)
km = kaplan_meier(age)
in_cluster = aalen_johansen(event_table(spans; time_axis = :cluster, t_end = t_last_myr))
```
"""
function event_table(
    spans::AbstractVector{BinarySpan};
    time_axis::Symbol = :age,
    t_end::Real,
    unit::Symbol = :myr,
    causes::AbstractDict{Symbol,Symbol} = _SURVIVAL_DEFAULT_CAUSES,
)::EventTable
    time_axis === :age ||
        time_axis === :cluster ||
        throw(ArgumentError("event_table: time_axis must be :age or :cluster; got :$time_axis"))
    unit === :myr ||
        unit === :nb ||
        throw(ArgumentError("event_table: unit must be :myr or :nb; got :$unit"))
    isfinite(t_end) || throw(ArgumentError("event_table: t_end must be finite; got $t_end"))
    t_record = Float64(t_end)
    time = Float64[]
    entry = Float64[]
    cause = Int[]
    id = Int[]
    labels = Symbol[]
    n_left = 0
    for (i, s) in enumerate(spans)
        t0 = unit === :myr ? s.t_form_myr : s.t_form_nb
        t1 = unit === :myr ? s.t_end_myr : s.t_end_nb
        (isfinite(t0) && isfinite(t1)) || throw(
            ArgumentError("event_table: span $i has a non-finite time (t_form = $t0, t_end = $t1)"),
        )
        t1 > t_record && throw(
            ArgumentError(
                "event_table: span $i ends at $t1, after the end of the record t_end = $t_record",
            ),
        )
        t1 > t0 || throw(
            ArgumentError("event_table: span $i has a lifetime ≤ 0 (t_form = $t0, t_end = $t1)"),
        )
        if time_axis === :age && s.left_censored
            n_left += 1
            continue
        end
        label = if s.right_censored
            :censored
        else
            haskey(causes, s.end_class) || throw(
                ArgumentError(
                    "event_table: end class :$(s.end_class) of span $i is not in `causes`",
                ),
            )
            causes[s.end_class]
        end
        k = 0
        if label !== :censored
            k = something(findfirst(==(label), labels), 0)
            if k == 0
                push!(labels, label)
                k = length(labels)
            end
        end
        t_exit = s.right_censored ? t_record : t1
        if time_axis === :age
            push!(entry, 0.0)
            push!(time, t_exit - t0)
        else
            push!(entry, s.left_censored ? 0.0 : t0)
            push!(time, t_exit)
        end
        push!(cause, k)
        push!(id, i)
    end
    n_left > 0 &&
        @info "event_table: $n_left left-censored span(s) left out of the age axis (age unknown)"
    return EventTable(; time = time, entry = entry, cause = cause, id = id, cause_labels = labels)
end

# ---------------------------------------------------------------------------
# Counting
# ---------------------------------------------------------------------------

function _survival_check(table::EventTable, caller::AbstractString)::Nothing
    n = length(table.time)
    (length(table.entry) == n && length(table.cause) == n && length(table.id) == n) ||
        throw(DimensionMismatch("$caller: the columns of the event table differ in length"))
    n_causes = length(table.cause_labels)
    for i in 1:n
        x = table.time[i]
        e = table.entry[i]
        (isfinite(x) && isfinite(e) && e < x) || throw(
            ArgumentError(
                "$caller: row $i needs finite times with entry < time; got entry = $e, time = $x",
            ),
        )
        0 <= table.cause[i] <= n_causes ||
            throw(ArgumentError("$caller: row $i has cause $(table.cause[i]) outside 0:$n_causes"))
    end
    return nothing
end

# Number at risk at t, #{entry < t ≤ time}, from the sorted exit and entry times:
# #{time ≥ t} − #{entry ≥ t}, since entry ≥ t implies time > t
_survival_at_risk(sorted_time::Vector{Float64}, sorted_entry::Vector{Float64}, t::Float64)::Int =
    searchsortedfirst(sorted_entry, t) - searchsortedfirst(sorted_time, t)

# Distinct event times with the numbers at risk, of events, of events per cause
# and of censorings in (t_(j−1), t_(j)], and the sorted censoring times
struct _SurvivalCounts
    time::Vector{Float64}
    n_risk::Vector{Int}
    n_events::Vector{Int}
    n_cause::Matrix{Int}
    n_censored::Vector{Int}
    censor_times::Vector{Float64}
end

function _survival_counts(table::EventTable)::_SurvivalCounts
    times = Float64[table.time[i] for i in eachindex(table.time) if table.cause[i] >= 1]
    sort!(unique!(times))
    m = length(times)
    n_causes = length(table.cause_labels)
    d = zeros(Int, m)
    dk = zeros(Int, m, n_causes)
    c = zeros(Int, m)
    censor_times = Float64[]
    for i in eachindex(table.time)
        x = table.time[i]
        k = table.cause[i]
        j = searchsortedfirst(times, x)
        if k >= 1
            d[j] += 1
            dk[j, k] += 1
        else
            push!(censor_times, x)
            j <= m && (c[j] += 1)
        end
    end
    sort!(censor_times)
    sorted_time = sort(table.time)
    sorted_entry = sort(table.entry)
    n_risk = Int[_survival_at_risk(sorted_time, sorted_entry, t) for t in times]
    return _SurvivalCounts(times, n_risk, d, dk, c, censor_times)
end

# ---------------------------------------------------------------------------
# Kaplan–Meier
# ---------------------------------------------------------------------------

"""
    SurvivalCurve

A product-limit survival curve from [`kaplan_meier`](@ref) or
[`turnbull`](@ref), evaluated at its jump times.

| field | meaning |
|---|---|
| `time` | distinct event times t_(j) (Turnbull: right ends of the innermost intervals) |
| `survival` | Ŝ(t_(j)) |
| `variance` | Greenwood variance of Ŝ(t_(j)); `NaN` when not estimated (Turnbull) |
| `n_risk` | number at risk n_j at t_(j) |
| `n_events` | number of events d_j at t_(j) |
| `n_censored` | censorings in (t_(j−1), t_(j)]; those after the last event time are in `censor_times` only |
| `censor_times` | every censoring time, sorted (tick marks of a plot) |

Evaluate with [`survival_at`](@ref); pointwise intervals from
[`confidence_interval`](@ref).
"""
Base.@kwdef struct SurvivalCurve
    time::Vector{Float64}
    survival::Vector{Float64}
    variance::Vector{Float64}
    n_risk::Vector{Int}
    n_events::Vector{Int}
    n_censored::Vector{Int}
    censor_times::Vector{Float64}
end

"""
    kaplan_meier(table::EventTable) -> SurvivalCurve

Kaplan–Meier product-limit estimate of the survival function of the exit time
for any cause, with the Greenwood variance, under right censoring and delayed
entry.

# Method
At the distinct event times t_(1) < t_(2) < … with n_j subjects at risk
(`entry < t_(j) ≤ time`) and d_j events (cause ≥ 1),

    Ŝ(t) = ∏_{t_(j) ≤ t} (1 − d_j/n_j),
    Var Ŝ(t) = Ŝ(t)² Σ_{t_(j) ≤ t} d_j / [n_j (n_j − d_j)].

Tied events enter through d_j; a censoring tied with an event time is at risk
at that time. A term with n_j = d_j (the last subject at risk fails) makes the
sum infinite, and the variance is stored as `Inf` from there on. With delayed
entry the estimate is conditional on survival to the earliest entry.

References: Kaplan & Meier (1958, DOI 10.1080/01621459.1958.10501452);
Greenwood (1926, Rep. Public Health Med. Subj. 33); Feigelson & Nelson (1985,
DOI 10.1086/163225) eqs. 8–9 for the astronomical use.

# Example
```julia
km = kaplan_meier(event_table(spans; t_end = t_last_myr))
survival_at(km, 10.0)               # Ŝ(10 Myr)
lower, upper = confidence_interval(km)
```
"""
function kaplan_meier(table::EventTable)::SurvivalCurve
    _survival_check(table, "kaplan_meier")
    c = _survival_counts(table)
    m = length(c.time)
    survival = Vector{Float64}(undef, m)
    variance = Vector{Float64}(undef, m)
    s = 1.0
    g = 0.0
    for j in 1:m
        n = c.n_risk[j]
        d = c.n_events[j]
        s *= 1 - d / n
        g += n > d ? d / (Float64(n) * (n - d)) : Inf
        survival[j] = s
        variance[j] = isinf(g) ? Inf : s^2 * g
    end
    return SurvivalCurve(;
        time = c.time,
        survival = survival,
        variance = variance,
        n_risk = c.n_risk,
        n_events = c.n_events,
        n_censored = c.n_censored,
        censor_times = c.censor_times,
    )
end

"""
    survival_at(curve::SurvivalCurve, t::Real) -> Float64

The step function Ŝ of `curve` at `t`: right-continuous, 1 before the first
jump time and the last value after the last.

# Example
```julia
km = kaplan_meier(table)
survival_at(km, 0.0)    # 1.0
```
"""
function survival_at(curve::SurvivalCurve, t::Real)::Float64
    j = searchsortedlast(curve.time, Float64(t))
    return j == 0 ? 1.0 : curve.survival[j]
end

"""
    _survival_normal_quantile(p) -> Float64

Quantile x of the standard normal distribution, Φ(x) = p with
Φ(x) = ½ erfc(−x/√2), by bisection on [−40, 40] to the floating-point
resolution. `ArgumentError` for p outside (0, 1).
"""
function _survival_normal_quantile(p::Real)::Float64
    0 < p < 1 || throw(ArgumentError("_survival_normal_quantile: p must lie in (0, 1); got $p"))
    target = Float64(p)
    lo = -40.0
    hi = 40.0
    for _ in 1:200
        mid = (lo + hi) / 2
        (mid == lo || mid == hi) && break
        if 0.5 * erfc(-mid / sqrt(2.0)) < target
            lo = mid
        else
            hi = mid
        end
    end
    return (lo + hi) / 2
end

"""
    confidence_interval(curve::SurvivalCurve; level::Real = 0.95, transform::Symbol = :loglog)
        -> (lower::Vector{Float64}, upper::Vector{Float64})

Pointwise confidence limits of Ŝ(t_(j)) from the variance stored in `curve`,
at confidence `level`.

# Method
With z the (1 + level)/2 quantile of the standard normal distribution and
σ² = Var Ŝ:

- `:loglog` (default): the interval of θ = log(−log Ŝ), whose delta-method
  standard error is σ/(Ŝ |log Ŝ|), transformed back,
  `exp(−exp(log(−log Ŝ) ± z σ/(Ŝ log Ŝ)))`; the limits lie in [0, 1] by
  construction. At Ŝ = 1 both limits are 1; at Ŝ = 0 they are `NaN` (undefined).
- `:plain`: Ŝ ∓ zσ clamped to [0, 1].

A `NaN` variance (not estimated, [`turnbull`](@ref)) gives `NaN` limits.
Reference: Feigelson & Nelson (1985, DOI 10.1086/163225) eqs. 11–12.
`ArgumentError` for a level outside (0, 1) or an unknown transform.

# Example
```julia
lower, upper = confidence_interval(kaplan_meier(table); level = 0.9)
```
"""
function confidence_interval(
    curve::SurvivalCurve;
    level::Real = 0.95,
    transform::Symbol = :loglog,
)::Tuple{Vector{Float64},Vector{Float64}}
    0 < level < 1 ||
        throw(ArgumentError("confidence_interval: level must lie in (0, 1); got $level"))
    transform === :loglog ||
        transform === :plain ||
        throw(
            ArgumentError(
                "confidence_interval: transform must be :loglog or :plain; got :$transform",
            ),
        )
    z = _survival_normal_quantile((1 + Float64(level)) / 2)
    m = length(curve.survival)
    lower = Vector{Float64}(undef, m)
    upper = Vector{Float64}(undef, m)
    for j in 1:m
        s = curve.survival[j]
        σ = sqrt(curve.variance[j])
        if transform === :plain
            lower[j] = clamp(s - z * σ, 0.0, 1.0)
            upper[j] = clamp(s + z * σ, 0.0, 1.0)
        elseif isnan(σ) || s <= 0
            lower[j] = NaN
            upper[j] = NaN
        elseif s >= 1
            lower[j] = 1.0
            upper[j] = 1.0
        else
            θ = log(-log(s))
            se = σ / (s * abs(log(s)))
            lower[j] = exp(-exp(θ + z * se))
            upper[j] = exp(-exp(θ - z * se))
        end
    end
    return lower, upper
end

# ---------------------------------------------------------------------------
# Nelson–Aalen
# ---------------------------------------------------------------------------

"""
    nelson_aalen(table::EventTable)
        -> NamedTuple{(:time, :hazard, :variance)}

Nelson–Aalen estimate of the cumulative hazard of the exit time for any cause,
under right censoring and delayed entry.

# Method
At the distinct event times t_(j) with n_j at risk and d_j events,

    Ĥ(t) = Σ_{t_(j) ≤ t} d_j/n_j,
    Var Ĥ(t) = Σ_{t_(j) ≤ t} d_j (n_j − d_j)/n_j³,

the variance being the one that accounts for tied events. A cause-specific
hazard follows from a table in which the other causes are recoded as
censorings.

References: Nelson (1972, Technometrics 14, 945); Aalen (1978,
DOI 10.1214/aos/1176344247).

# Example
```julia
na = nelson_aalen(table)
na.hazard[end]      # Ĥ at the last event time
```
"""
function nelson_aalen(
    table::EventTable,
)::NamedTuple{(:time, :hazard, :variance),Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}}
    _survival_check(table, "nelson_aalen")
    c = _survival_counts(table)
    m = length(c.time)
    hazard = Vector{Float64}(undef, m)
    variance = Vector{Float64}(undef, m)
    h = 0.0
    v = 0.0
    for j in 1:m
        n = Float64(c.n_risk[j])
        d = Float64(c.n_events[j])
        h += d / n
        v += d * (n - d) / n^3
        hazard[j] = h
        variance[j] = v
    end
    return (time = c.time, hazard = hazard, variance = variance)
end

# ---------------------------------------------------------------------------
# Aalen–Johansen
# ---------------------------------------------------------------------------

"""
    CumulativeIncidence

Cumulative incidence functions of competing causes from
[`aalen_johansen`](@ref), at the distinct event times of any cause.

| field | meaning |
|---|---|
| `time` | distinct event times t_(j) |
| `incidence` | n_times × n_causes, F̂_k(t_(j)) |
| `variance` | n_times × n_causes, delta-method variance of F̂_k(t_(j)) |
| `survival` | the all-cause Kaplan–Meier Ŝ(t_(j)) |
| `n_risk` | number at risk n_j |
| `n_events` | n_times × n_causes, events d_kj of cause k at t_(j) |
| `cause_labels` | the labels of the causes (columns) |
"""
Base.@kwdef struct CumulativeIncidence
    time::Vector{Float64}
    incidence::Matrix{Float64}
    variance::Matrix{Float64}
    survival::Vector{Float64}
    n_risk::Vector{Int}
    n_events::Matrix{Int}
    cause_labels::Vector{Symbol}
end

"""
    aalen_johansen(table::EventTable) -> CumulativeIncidence

Aalen–Johansen estimate of the cumulative incidence F_k(t) = P(T ≤ t, cause k)
of every cause of the table, under right censoring and delayed entry.

# Method
With n_j at risk, d_j events of any cause and d_kj of cause k at the distinct
event times t_(j), and Ŝ(t_(j)−) the all-cause Kaplan–Meier estimate just
before t_(j),

    F̂_k(t) = Σ_{t_(j) ≤ t} Ŝ(t_(j)−) d_kj/n_j,

so that Σ_k F̂_k = 1 − Ŝ. The variance is the delta-method estimate
(Marubini & Valsecchi 1995, eq. 10.11; Aalen 1978)

    Var F̂_k(t) = Σ_{t_(j) ≤ t} [F̂_k(t) − F̂_k(t_(j))]² d_j/[n_j (n_j − d_j)]
               + Σ_{t_(j) ≤ t} Ŝ(t_(j)−)² (n_j − d_kj) d_kj/n_j³
               − 2 Σ_{t_(j) ≤ t} [F̂_k(t) − F̂_k(t_(j))] Ŝ(t_(j)−) d_kj/n_j²,

accumulated in one pass over the event times through the recursions of the
three sums in the increments F̂_k(t_(i)) − F̂_k(t_(i−1)) ≥ 0, which avoids the
cancellation of the expanded squares. A term with n_j = d_j has a vanishing
prefactor (F̂_k is constant after it) and contributes 0. The complement
1 − Ŝ_k of the Kaplan–Meier estimate with the other causes censored
overestimates F_k; it is not a cumulative incidence under competing causes.

References: Aalen & Johansen (1978, Scand. J. Stat. 5, 141); Aalen (1978,
DOI 10.1214/aos/1176344247); Marubini & Valsecchi (1995, Analysing Survival
Data from Clinical Trials and Observational Studies, Wiley).

# Example
```julia
ci = aalen_johansen(event_table(spans; time_axis = :cluster, t_end = t_last_myr))
k = findfirst(==(:gw_merger), ci.cause_labels)
ci.incidence[end, k]        # fraction of binaries ended by a GW merger
```
"""
function aalen_johansen(table::EventTable)::CumulativeIncidence
    _survival_check(table, "aalen_johansen")
    c = _survival_counts(table)
    m = length(c.time)
    n_causes = length(table.cause_labels)
    incidence = zeros(m, n_causes)
    variance = zeros(m, n_causes)
    survival = Vector{Float64}(undef, m)
    s_before = Vector{Float64}(undef, m)
    w = Vector{Float64}(undef, m)      # d_j/[n_j (n_j − d_j)]
    s = 1.0
    for j in 1:m
        n = Float64(c.n_risk[j])
        d = c.n_events[j]
        s_before[j] = s
        for k in 1:n_causes
            prev = j > 1 ? incidence[j - 1, k] : 0.0
            incidence[j, k] = prev + s * c.n_cause[j, k] / n
        end
        s *= 1 - d / n
        survival[j] = s
        w[j] = c.n_risk[j] > d ? d / (n * (n - d)) : Inf
    end
    for k in 1:n_causes
        # T1 = Σ (F_i − F_j)² w_j, G = Σ (F_i − F_j) w_j, A = Σ w_j,
        # T3 = Σ (F_i − F_j) h_j, H = Σ h_j with h_j = Ŝ(t_j−) d_kj/n_j², T2 = Σ second terms
        t1 = 0.0
        g = 0.0
        a = 0.0
        t2 = 0.0
        t3 = 0.0
        h = 0.0
        f_prev = 0.0
        for i in 1:m
            f = incidence[i, k]
            δ = f - f_prev
            if δ != 0
                t1 += 2δ * g + δ^2 * a
                g += δ * a
                t3 += δ * h
            end
            n = Float64(c.n_risk[i])
            dk = Float64(c.n_cause[i, k])
            a += w[i]
            h += s_before[i] * dk / n^2
            t2 += s_before[i]^2 * (n - dk) * dk / n^3
            variance[i, k] = t1 + t2 - 2t3
            f_prev = f
        end
    end
    return CumulativeIncidence(;
        time = c.time,
        incidence = incidence,
        variance = variance,
        survival = survival,
        n_risk = c.n_risk,
        n_events = c.n_cause,
        cause_labels = copy(table.cause_labels),
    )
end

# ---------------------------------------------------------------------------
# Turnbull
# ---------------------------------------------------------------------------

"""
    turnbull(lower::AbstractVector{<:Real}, upper::AbstractVector{<:Real};
             tol::Real = 1e-8, max_iter::Integer = 10_000) -> SurvivalCurve

Turnbull nonparametric maximum-likelihood estimate of the survival function
from interval-censored exit times T_i ∈ (L_i, U_i].

# Arguments
- `lower`, `upper`: L_i and U_i; an exact time has L_i = U_i (the interval
  (U_i − 0, U_i]), a right-censored one U_i = `Inf`
- `tol`: the iteration stops when max_m |Δp_m| < `tol`
- `max_iter`: largest number of iterations; reaching it is logged as a warning

# Method
The left ends L_i and the right ends U_i are sorted together, a right end
before a left end of the same value ((a, x] and (x, b] are disjoint) and the
left end of an exact time before both. The innermost intervals (l, u] are the
pairs of a left end immediately followed by a right end in that order; the
NPMLE places its mass on them only. Every observation covers a contiguous block
of innermost intervals. The masses p_m, uniform at the start, follow the
self-consistency (EM) iteration

    p_m ← (1/N) Σ_i α_im p_m / Σ_m' α_im' p_m',

α_im = 1 when innermost interval m lies in (L_i, U_i]. The curve is
Ŝ(u_m) = 1 − Σ_{m' ≤ m} p_m' at the finite right ends u_m; the mass of an
innermost interval with u = `Inf` (beyond the last finite end) lowers no point
of it. For exact and right-censored data the estimate equals the Kaplan–Meier
estimate. The variance is not estimated and is `NaN`; `n_risk`, `n_events` and
`n_censored` are 0, and `censor_times` holds the L_i of the right-censored
observations.

Reference: Turnbull (1976, DOI 10.1111/j.2517-6161.1976.tb01597.x).

Throws `DimensionMismatch` for vectors of different length and `ArgumentError`
for a non-finite or `NaN` bound other than U_i = `Inf`, L_i > U_i, `tol ≤ 0` or
`max_iter < 1`.

# Example
```julia
# Ends known to the snapshot cadence: formation at the age 0, end in (t_j, t_{j+1}]
curve = turnbull(age_last_seen, age_first_missing)
survival_at(curve, 5.0)
```
"""
function turnbull(
    lower::AbstractVector{<:Real},
    upper::AbstractVector{<:Real};
    tol::Real = 1e-8,
    max_iter::Integer = 10_000,
)::SurvivalCurve
    n = length(lower)
    length(upper) == n || throw(
        DimensionMismatch("turnbull: lower and upper differ in length ($n vs $(length(upper)))"),
    )
    (isfinite(tol) && tol > 0) ||
        throw(ArgumentError("turnbull: tol must be finite and > 0; got $tol"))
    max_iter >= 1 || throw(ArgumentError("turnbull: max_iter must be ≥ 1; got $max_iter"))
    L = Float64[x for x in lower]
    U = Float64[x for x in upper]
    for i in 1:n
        (isfinite(L[i]) && !isnan(U[i]) && U[i] > -Inf) ||
            throw(ArgumentError("turnbull: observation $i has the bounds ($(L[i]), $(U[i])]"))
        L[i] <= U[i] ||
            throw(ArgumentError("turnbull: observation $i has lower > upper ($(L[i]) > $(U[i]))"))
    end
    empty_curve = SurvivalCurve(;
        time = Float64[],
        survival = Float64[],
        variance = Float64[],
        n_risk = Int[],
        n_events = Int[],
        n_censored = Int[],
        censor_times = Float64[],
    )
    n == 0 && return empty_curve

    # End points as (value, rank): 0 left end of an exact time, 1 right end, 2 left end
    left_key = Tuple{Float64,Int}[(L[i], L[i] == U[i] ? 0 : 2) for i in 1:n]
    ends = vcat(left_key, Tuple{Float64,Int}[(U[i], 1) for i in 1:n])
    sort!(ends)
    inner_left = Tuple{Float64,Int}[]
    inner_right = Float64[]
    for q in 1:(length(ends) - 1)
        if ends[q][2] != 1 && ends[q + 1][2] == 1
            push!(inner_left, ends[q])
            push!(inner_right, ends[q + 1][1])
        end
    end
    M = length(inner_right)
    lo = Vector{Int}(undef, n)
    hi = Vector{Int}(undef, n)
    for i in 1:n
        lo[i] = searchsortedfirst(inner_left, left_key[i])
        hi[i] = searchsortedlast(inner_right, U[i])
        lo[i] <= hi[i] || error("turnbull: observation $i covers no innermost interval")
    end

    # Self-consistency iteration
    p = fill(1.0 / M, M)
    cum = zeros(M + 1)
    acc = zeros(M + 1)
    converged = false
    for _ in 1:max_iter
        for q in 1:M
            cum[q + 1] = cum[q] + p[q]
        end
        fill!(acc, 0.0)
        for i in 1:n
            den = cum[hi[i] + 1] - cum[lo[i]]
            wt = den > 0 ? 1 / den : 0.0
            acc[lo[i]] += wt
            acc[hi[i] + 1] -= wt
        end
        running = 0.0
        δ = 0.0
        for q in 1:M
            running += acc[q]
            p_new = p[q] * running / n
            δ = max(δ, abs(p_new - p[q]))
            p[q] = p_new
        end
        if δ < tol
            converged = true
            break
        end
    end
    converged ||
        @warn "turnbull: the self-consistency iteration did not reach tol = $tol within $max_iter iterations"
    p ./= sum(p)

    time = Float64[]
    survival = Float64[]
    total = 0.0
    for q in 1:M
        total += p[q]
        isfinite(inner_right[q]) || continue
        push!(time, inner_right[q])
        push!(survival, max(0.0, 1 - total))
    end
    k = length(time)
    return SurvivalCurve(;
        time = time,
        survival = survival,
        variance = fill(NaN, k),
        n_risk = zeros(Int, k),
        n_events = zeros(Int, k),
        n_censored = zeros(Int, k),
        censor_times = sort!(Float64[L[i] for i in 1:n if U[i] == Inf]),
    )
end

# ---------------------------------------------------------------------------
# Log-rank test
# ---------------------------------------------------------------------------

"""
    logrank_test(table::EventTable, group::AbstractVector{<:Integer})
        -> NamedTuple{(:statistic, :p_value, :observed, :expected)}

Log-rank (Mantel–Haenszel) test of equal hazards of two samples, events of any
cause, under right censoring and delayed entry.

# Arguments
- `table`: the pooled sample
- `group`: group label of each row, 1 or 2; both must occur

# Method
At every distinct event time t_(j) of the pooled sample with n_j at risk, n_1j
of them in group 1, and d_j events, the hypergeometric expectation and variance
of the events of group 1 are

    E_1j = d_j n_1j/n_j,
    V_j = d_j (n_1j/n_j)(1 − n_1j/n_j)(n_j − d_j)/(n_j − 1)   (0 for n_j = 1).

The statistic (O_1 − E_1)²/V with O_1 = Σ d_1j, E_1 = Σ E_1j and V = Σ V_j is
χ² with one degree of freedom under the null hypothesis, so
p = erfc(√(statistic/2)). `observed` and `expected` are the vectors
[O_1, O_2] and [E_1, E_2]; E_1 + E_2 = O_1 + O_2. Without variance (no events)
the statistic and p are `NaN`.

References: Mantel (1966, Cancer Chemother. Rep. 50, 163); Peto & Peto (1972,
J. R. Stat. Soc. A 135, 185); Feigelson & Nelson (1985, DOI 10.1086/163225).

Throws `DimensionMismatch` when `group` and the table differ in length and
`ArgumentError` unless the labels are exactly {1, 2}.

# Example
```julia
age = event_table(spans; t_end = t_last_myr)
group = [spans[i].origin === :primordial ? 1 : 2 for i in age.id]
logrank_test(age, group).p_value
```
"""
function logrank_test(
    table::EventTable,
    group::AbstractVector{<:Integer},
)::NamedTuple{
    (:statistic, :p_value, :observed, :expected),
    Tuple{Float64,Float64,Vector{Int},Vector{Float64}},
}
    _survival_check(table, "logrank_test")
    length(group) == length(table.time) || throw(
        DimensionMismatch(
            "logrank_test: group has $(length(group)) labels for $(length(table.time)) rows",
        ),
    )
    labels = sort!(unique(Int[g for g in group]))
    labels == [1, 2] || throw(
        ArgumentError(
            "logrank_test: needs exactly two groups labelled 1 and 2; got labels $labels",
        ),
    )
    c = _survival_counts(table)
    m = length(c.time)
    rows1 = findall(==(1), group)
    sorted_time1 = sort(table.time[rows1])
    sorted_entry1 = sort(table.entry[rows1])
    d1 = zeros(Int, m)
    for i in rows1
        table.cause[i] >= 1 && (d1[searchsortedfirst(c.time, table.time[i])] += 1)
    end
    o1 = 0
    e1 = 0.0
    v = 0.0
    for j in 1:m
        n = Float64(c.n_risk[j])
        d = Float64(c.n_events[j])
        n1 = Float64(_survival_at_risk(sorted_time1, sorted_entry1, c.time[j]))
        e1 += d * n1 / n
        o1 += d1[j]
        n > 1 && (v += d * (n1 / n) * (1 - n1 / n) * (n - d) / (n - 1))
    end
    total = sum(c.n_events)
    statistic = v > 0 ? (o1 - e1)^2 / v : NaN
    p_value = isnan(statistic) ? NaN : erfc(sqrt(statistic / 2))
    return (
        statistic = statistic,
        p_value = p_value,
        observed = [o1, total - o1],
        expected = [e1, total - e1],
    )
end

# ---------------------------------------------------------------------------
# Cox–Snell residuals
# ---------------------------------------------------------------------------

"""
    cox_snell_residuals(table::EventTable, cumulative_hazard)
        -> NamedTuple{(:residual, :event)}

Cox–Snell residuals of a predicted cumulative hazard Λ.

# Arguments
- `table`: the observed sample
- `cumulative_hazard`: a function t ↦ Λ(t), non-decreasing, returning a real

# Method
r_i = Λ(time_i) − Λ(entry_i), the integrated hazard over the time at risk of
subject i, with the event indicator `cause_i ≥ 1`. If Λ is the true cumulative
hazard, the r_i are a right-censored sample of the unit exponential
distribution, so the Nelson–Aalen estimate of an event table built from them
(time r_i, entry 0, the same indicator) lies on the identity line; its slope
measures the ratio of the observed to the predicted hazard.

Reference: Cox & Snell (1968, J. R. Stat. Soc. B 30, 248).

Throws `ArgumentError` for a residual that is not finite or is negative (Λ
decreasing).

# Example
```julia
cs = cox_snell_residuals(age, τ -> τ / t_evap)
check = nelson_aalen(EventTable(; time = cs.residual, entry = zeros(length(cs.residual)),
    cause = Int.(cs.event), id = age.id, cause_labels = [:end]))
```
"""
function cox_snell_residuals(
    table::EventTable,
    cumulative_hazard,
)::NamedTuple{(:residual, :event),Tuple{Vector{Float64},Vector{Bool}}}
    _survival_check(table, "cox_snell_residuals")
    n = length(table.time)
    residual = Vector{Float64}(undef, n)
    event = Vector{Bool}(undef, n)
    for i in 1:n
        λ_exit = Float64(cumulative_hazard(table.time[i])::Real)
        λ_entry = Float64(cumulative_hazard(table.entry[i])::Real)
        r = λ_exit - λ_entry
        (isfinite(r) && r >= 0) || throw(
            ArgumentError(
                "cox_snell_residuals: row $i has the residual $r; Λ must be finite and non-decreasing",
            ),
        )
        residual[i] = r
        event[i] = table.cause[i] >= 1
    end
    return (residual = residual, event = event)
end

# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

"""
    write_survival_curve(path, curve::SurvivalCurve) -> String

Write a survival curve as CSV with the columns `time`, `survival`, `variance`,
`lower`, `upper` (95 % log–log limits of [`confidence_interval`](@ref)),
`n_risk`, `n_events`, `n_censored`, one row per jump time; non-finite values as
empty cells. An existing file is backed up, never overwritten. Returns the path.

# Example
```julia
write_survival_curve(joinpath(dir, "survival_age.csv"), kaplan_meier(age))
```
"""
function write_survival_curve(path::AbstractString, curve::SurvivalCurve)::String
    lower, upper = confidence_interval(curve)
    _backup_existing(path)
    open(path, "w") do io
        println(io, "time,survival,variance,lower,upper,n_risk,n_events,n_censored")
        for j in eachindex(curve.time)
            cells = (
                _csv_precise(curve.time[j]),
                _csv_precise(curve.survival[j]),
                _csv_precise(curve.variance[j]),
                _csv_precise(lower[j]),
                _csv_precise(upper[j]),
                _csv_precise(curve.n_risk[j]),
                _csv_precise(curve.n_events[j]),
                _csv_precise(curve.n_censored[j]),
            )
            println(io, join(cells, ","))
        end
    end
    return String(path)
end

"""
    write_cumulative_incidence(path, ci::CumulativeIncidence) -> String

Write cumulative incidence functions as CSV in long format, one row per event
time and cause, with the columns `time`, `cause` (label), `incidence`,
`variance`, `survival` (all-cause Kaplan–Meier), `n_risk`, `n_events` (events
of that cause); non-finite values as empty cells. An existing file is backed
up, never overwritten. Returns the path.

# Example
```julia
write_cumulative_incidence(joinpath(dir, "incidence_cluster.csv"), aalen_johansen(table))
```
"""
function write_cumulative_incidence(path::AbstractString, ci::CumulativeIncidence)::String
    _backup_existing(path)
    open(path, "w") do io
        println(io, "time,cause,incidence,variance,survival,n_risk,n_events")
        for j in eachindex(ci.time), k in eachindex(ci.cause_labels)
            cells = (
                _csv_precise(ci.time[j]),
                _csv_precise(String(ci.cause_labels[k])),
                _csv_precise(ci.incidence[j, k]),
                _csv_precise(ci.variance[j, k]),
                _csv_precise(ci.survival[j]),
                _csv_precise(ci.n_risk[j]),
                _csv_precise(ci.n_events[j, k]),
            )
            println(io, join(cells, ","))
        end
    end
    return String(path)
end
