# =============================================================================
# Cost model: wall time, segments and disk of an engine launch
# =============================================================================

"""
    _COST_TABLE_H200

Wall seconds per N-body time unit of the GPU build on one NVIDIA H200 NVL
with 8 host threads, for one single cluster at rest (King model), measured
with the `bench/gpu_cells.jl` cells and `bench/gpu_scaling.jl` in September
2026 on host `issaf-0-6`:

| N | s per N-body unit | source |
|:--|:--|:--|
| 1.0×10⁵ | 64 | `gpu_scaling_20260924_014138.csv`, 0.25 time units (start-up included; an upper bound) |
| 5.0×10⁵ | 281 | `gpu_cells_h200_weak_20260924_184842.csv`, single cell, 2.0 time units |
| 1.0×10⁶ | 1087 | `gpu_cells_h200_strong_20260925_040734.csv`, single cell, 2.0 time units |

The two-cluster (dual) cells cost 1.85–1.99 times the single cell at equal
N (559 against 281 s at 5×10⁵, 2016 against 1087 s at 10⁶), because the time
unit of a merger is set by the whole configuration; hence the default
`cost.merger_factor = 2.0`. The device pilot (roadmap P3) replaces these
rows through `[cost] table`.
"""
const _COST_TABLE_H200 = [(1.0e5, 64.0), (5.0e5, 281.0), (1.0e6, 1087.0)]

"""
    _CONF3_BYTES_PER_BODY

Bytes per body of one `conf.3` snapshot: name, mass, three positions, three
velocities, density and potential, each a 4-byte `Int32` or `Float32`.
"""
const _CONF3_BYTES_PER_BODY = 40

"""
    _TEXT_BYTES_PER_BODY

Bytes per star of the stellar- and binary-evolution text output written at
every snapshot (`sev.83` and `bev.82` lines together, approximately).
"""
const _TEXT_BYTES_PER_BODY = 130

"""
    _cost_row(row) -> Tuple{Float64,Float64}

One `(N, s_per_nb)` row of a cost table, from a tuple or from a two-element
vector (the form of `cost.table` in the configuration). Raises an
`ArgumentError` for a vector of another length.
"""
_cost_row(row::Tuple{Real,Real})::Tuple{Float64,Float64} = (Float64(row[1]), Float64(row[2]))

function _cost_row(row::AbstractVector{<:Real})::Tuple{Float64,Float64}
    length(row) == 2 ||
        throw(ArgumentError("cost table rows must be [N, s_per_nb] pairs; got $row"))
    return (Float64(row[1]), Float64(row[2]))
end

"""
    _cost_rows(table) -> Vector{Tuple{Float64,Float64}}

The rows of a cost table as `(N, s_per_nb)` tuples, checked: at least two
rows, positive entries, strictly increasing `N`. Raises an `ArgumentError`
otherwise.
"""
function _cost_rows(table::AbstractVector)::Vector{Tuple{Float64,Float64}}
    rows = Tuple{Float64,Float64}[_cost_row(r) for r in table]
    length(rows) ≥ 2 || throw(
        ArgumentError(
            "a cost table needs at least two rows (N, s_per_nb) to interpolate; " *
            "got $(length(rows))",
        ),
    )
    all(r -> r[1] > 0 && r[2] > 0, rows) ||
        throw(ArgumentError("cost table entries must be > 0; got $rows"))
    all(k -> rows[k + 1][1] > rows[k][1], 1:(length(rows) - 1)) || throw(
        ArgumentError("cost table rows must be in strictly increasing N; got N = $(first.(rows))"),
    )
    return rows
end

"""
    _cost_n_label(n) -> String

Particle number of a table row as an integer string, for the basis of an
estimate.
"""
_cost_n_label(n::Real)::String = string(round(Int, n))

"""
    _cost_rate_basis(n, table) -> Tuple{Float64,String}

Seconds per N-body unit at `n` bodies from `table` (see [`cost_rate`](@ref))
together with the basis of the value: `"interpolated between N = … and N =
…"`, `"extrapolated below the table (N < …)"` or `"extrapolated above the
table (N > …)"`. A table row is reproduced exactly.
"""
function _cost_rate_basis(n::Real, table::AbstractVector)::Tuple{Float64,String}
    rows = _cost_rows(table)
    x = Float64(n)
    (isfinite(x) && x > 0) || throw(ArgumentError("cost: N must be > 0; got $n"))
    m = length(rows)
    k, basis = if x < rows[1][1]
        (1, "extrapolated below the table (N < $(_cost_n_label(rows[1][1])))")
    elseif x > rows[m][1]
        (m - 1, "extrapolated above the table (N > $(_cost_n_label(rows[m][1])))")
    else
        j = min(searchsortedlast(first.(rows), x), m - 1)
        (
            j,
            "interpolated between N = $(_cost_n_label(rows[j][1])) and " *
            "N = $(_cost_n_label(rows[j+1][1]))",
        )
    end
    i = findfirst(r -> r[1] == x, rows)
    i === nothing || return (rows[i][2], basis)
    (n₁, r₁), (n₂, r₂) = rows[k], rows[k + 1]
    α = log(r₂ / r₁) / log(n₂ / n₁)
    return (r₁ * (x / n₁)^α, basis)
end

"""
    cost_rate(n; table = _COST_TABLE_H200) -> Float64

Wall seconds per N-body time unit of an engine run with `n` bodies:
log–log linear interpolation between the two table rows that bracket `n`;
outside the table, extrapolation with the exponent of the nearest segment.

# Arguments
- `n`: particle number, > 0
- `table`: rows `(N, s_per_nb)` with strictly increasing `N` and positive
  rates, at least two; tuples, or the two-element vectors of `cost.table`.
  The default is the built-in H200 table (`_COST_TABLE_H200`)

Raises an `ArgumentError` for a table with fewer than two rows, unordered
`N` or non-positive entries.

# Example
```julia
cost_rate(5e5)                                     # 281.0
cost_rate(2e5; table = [[1e5, 60.0], [1e6, 1000.0]])
```
"""
cost_rate(n::Real; table::AbstractVector = _COST_TABLE_H200)::Float64 =
    first(_cost_rate_basis(n, table))

"""
    CostEstimate

Estimated cost of one engine launch, from [`estimate_cost`](@ref).

# Fields
- `n`: particle number
- `span_nb`: integrated span [N-body time units]
- `rate_s_per_nb`: wall seconds per N-body time unit, merger factor included
- `wall_s`: estimated wall time [s]
- `segments`: engine segments the wall time needs under
  `simulation.wall_budget` (1 without a budget)
- `dump_bytes`: restart dumps [bytes]
- `snapshot_bytes`: `conf.3` snapshots [bytes]
- `text_bytes`: `sev.83` and `bev.82` output [bytes]
- `disk_bytes`: the sum of the three [bytes]
- `budget_hours`: the configured budget (`cost.budget_hours`) [h]; `0` = none
- `exceeds`: whether `wall_s` exceeds a non-zero budget
- `basis`: how the rate was obtained from the table (interpolated or
  extrapolated, with the merger factor when applied)
"""
struct CostEstimate
    n::Int
    span_nb::Float64
    rate_s_per_nb::Float64
    wall_s::Float64
    segments::Int
    dump_bytes::Int
    snapshot_bytes::Int
    text_bytes::Int
    disk_bytes::Int
    budget_hours::Float64
    exceeds::Bool
    basis::String
end

"""
    show(io, ::MIME"text/plain", e::CostEstimate)

Four-line summary of an estimate: wall time (with the budget when one is
set), segments, disk with its three parts, and the basis of the rate.
"""
function Base.show(io::IO, ::MIME"text/plain", e::CostEstimate)
    budget = if e.budget_hours == 0
        ""
    elseif e.exceeds
        "; exceeds the budget of $(e.budget_hours) h"
    else
        "; within the budget of $(e.budget_hours) h"
    end
    println(
        io,
        "Wall time  ",
        _format_elapsed(e.wall_s),
        " for $(e.span_nb) N-body units at N = $(e.n) ",
        "($(round(e.rate_s_per_nb; sigdigits = 3)) s per unit)",
        budget,
    )
    println(io, "Segments   ", e.segments)
    println(
        io,
        "Disk       ",
        Base.format_bytes(e.disk_bytes),
        " (dumps ",
        Base.format_bytes(e.dump_bytes),
        ", snapshots ",
        Base.format_bytes(e.snapshot_bytes),
        ", text ",
        Base.format_bytes(e.text_bytes),
        ")",
    )
    print(io, "Basis      ", e.basis)
    return nothing
end

"""
    estimate_cost(cfg; n, span_nb, nnbopt = 0, dtadj = 0.0, deltat = 0.0,
                  merger = cfg.merger.enabled) -> CostEstimate

Wall time, segments and disk of an engine run of `n` bodies over `span_nb`
N-body time units.

The rate is [`cost_rate`](@ref) on `cfg.cost.table` (the built-in H200
table when empty), multiplied by `cfg.cost.merger_factor` for a two-cluster
configuration (`merger`); the wall time is the rate times `span_nb`. With a
non-zero `simulation.wall_budget` the number of segments is the wall time
over `wall_budget − stop_margin`, rounded up; otherwise 1.

The disk is the sum of
- restart dumps: `checkpoint_keep` of them when it is non-zero, otherwise
  one per adjustment interval `dtadj` (`⌊span_nb/dtadj⌋ + 1`), at least one
  terminal dump; each of the size of `_checkpoint_dump_bytes(n, nnbopt)`;
- `conf.3` snapshots, 40 bytes per body at each of `⌊span_nb/deltat⌋ + 1`
  outputs (none when `deltat = 0`);
- `sev.83` and `bev.82` text, about 130 bytes per star at each output.

`exceeds` is set when `cfg.cost.budget_hours > 0` and the wall time is
longer than the budget.

# Arguments
- `cfg`: the configuration (`[cost]`, `simulation.wall_budget`,
  `simulation.stop_margin`, `simulation.checkpoint_keep`)
- `n`: particle number (an integral value)
- `span_nb`: integrated span [N-body time units], > 0
- `nnbopt`: target neighbour number of the dump estimate, ≥ 0
- `dtadj`, `deltat`: adjustment and snapshot intervals [N-body time units],
  ≥ 0; `0` = not known
- `merger`: apply `cost.merger_factor`

# Example
```julia
e = estimate_cost(cfg; n = 500_000, span_nb = 100, nnbopt = 300, dtadj = 1.0, deltat = 10.0)
e.wall_s        # 28100.0 with the built-in table
```
"""
function estimate_cost(
    cfg::Nbody6Config;
    n::Real,
    span_nb::Real,
    nnbopt::Integer = 0,
    dtadj::Real = 0.0,
    deltat::Real = 0.0,
    merger::Bool = cfg.merger.enabled,
)::CostEstimate
    (n ≥ 1 && isinteger(n)) || throw(ArgumentError("cost: n must be an integer ≥ 1; got $n"))
    (isfinite(span_nb) && span_nb > 0) ||
        throw(ArgumentError("cost: span_nb must be > 0; got $span_nb"))
    nnbopt ≥ 0 || throw(ArgumentError("cost: nnbopt must be ≥ 0; got $nnbopt"))
    dtadj ≥ 0 || throw(ArgumentError("cost: dtadj must be ≥ 0; got $dtadj"))
    deltat ≥ 0 || throw(ArgumentError("cost: deltat must be ≥ 0; got $deltat"))
    c = cfg.cost
    sim = cfg.simulation
    nbodies = Int(n)
    span = Float64(span_nb)

    table = isempty(c.table) ? _COST_TABLE_H200 : c.table
    rate, basis = _cost_rate_basis(nbodies, table)
    if merger
        rate *= c.merger_factor
        basis *= " × merger factor $(c.merger_factor)"
    end
    wall_s = rate * span

    segments = if sim.wall_budget > 0
        usable = sim.wall_budget - sim.stop_margin
        usable > 0 || throw(
            ArgumentError(
                "cost: simulation.wall_budget $(sim.wall_budget) s leaves no time after " *
                "simulation.stop_margin $(sim.stop_margin) s",
            ),
        )
        max(1, ceil(Int, wall_s / usable))
    else
        1
    end

    n_dumps = if sim.checkpoint_keep > 0
        sim.checkpoint_keep
    elseif dtadj > 0
        floor(Int, span / dtadj) + 1
    else
        1
    end
    n_outputs = deltat > 0 ? floor(Int, span / deltat) + 1 : 0
    dump_bytes = _checkpoint_dump_bytes(nbodies, nnbopt) * n_dumps
    snapshot_bytes = _CONF3_BYTES_PER_BODY * nbodies * n_outputs
    text_bytes = _TEXT_BYTES_PER_BODY * nbodies * n_outputs

    exceeds = c.budget_hours > 0 && wall_s > 3600 * c.budget_hours
    return CostEstimate(
        nbodies,
        span,
        rate,
        wall_s,
        segments,
        dump_bytes,
        snapshot_bytes,
        text_bytes,
        dump_bytes + snapshot_bytes + text_bytes,
        c.budget_hours,
        exceeds,
        basis,
    )
end

"""
    estimate_cost(cfg, inp_path; span_nb = nothing, merger = cfg.merger.enabled)
        -> CostEstimate

[`estimate_cost`](@ref) for the engine input file `inp_path`: `N` and
`NNBOPT` from the file, the span from its `TCRIT` unless `span_nb` is given
(a restart passes the remaining span), and the intervals `DTADJ` and
`DELTAT` from the file when present. Raises an `ArgumentError` when `N` or
`NNBOPT` cannot be read, an error when the span is not given and the file
has no `TCRIT`.

# Example
```julia
estimate_cost(cfg, "input_files/engine/N1k_quick.inp")
```
"""
function estimate_cost(
    cfg::Nbody6Config,
    inp_path::AbstractString;
    span_nb::Union{Nothing,Real} = nothing,
    merger::Bool = cfg.merger.enabled,
)::CostEstimate
    sizes = _input_sizes(inp_path)
    sizes === nothing && throw(ArgumentError("cost: N and NNBOPT could not be read from $inp_path"))
    span = span_nb === nothing ? _tcrit_of(inp_path) : Float64(span_nb)
    return estimate_cost(
        cfg;
        n = sizes.n,
        span_nb = span,
        nnbopt = sizes.nnbopt,
        dtadj = something(_input_value(inp_path, "DTADJ"), 0.0),
        deltat = something(_input_value(inp_path, "DELTAT"), 0.0),
        merger = merger,
    )
end

"""
    _input_value(inp_path, name) -> Union{Nothing,Float64}

Value of the first `name=` entry of an engine input file outside comment
lines (first non-blank character `!`), Fortran `D` exponents accepted, as
[`_tcrit_of`](@ref) reads `TCRIT`; `nothing` when the file has no such
entry. `name` is a Fortran identifier (`"DTADJ"`, `"DELTAT"`).
"""
function _input_value(inp_path::AbstractString, name::AbstractString)::Union{Nothing,Float64}
    occursin(r"^[A-Za-z][A-Za-z0-9_]*$", name) ||
        throw(ArgumentError("not an input-file parameter name: \"$name\""))
    pattern = Regex("\\b" * name * "\\s*=\\s*([0-9.eEdD+-]+)")
    return open(inp_path) do io
        for line in eachline(io)
            startswith(lstrip(line), '!') && continue
            m = match(pattern, line)
            m === nothing || return _parse_fortran_float(something(m.captures[1]))
        end
        return nothing
    end
end

"""
    _launch_cost_estimate(cfg, input_path, span) -> Union{Nothing,CostEstimate}

The cost estimate of an engine launch from `input_path` over `span` N-body
time units (`nothing`: the input's `TCRIT`), logged in one line. `nothing`
when `cfg.cost.enabled` is false, and when the estimate cannot be made from
the input (logged as a warning; a launch is never aborted for it). An
estimate beyond `cost.budget_hours` raises an `ArgumentError` with
`cost.on_exceed = "refuse"` and is logged as a warning otherwise.
"""
function _launch_cost_estimate(
    cfg::Nbody6Config,
    input_path::AbstractString,
    span::Union{Nothing,Real},
)::Union{Nothing,CostEstimate}
    cfg.cost.enabled || return nothing
    estimate = try
        estimate_cost(cfg, input_path; span_nb = span)
    catch e
        e isa Union{ArgumentError,ErrorException,SystemError,Base.IOError} || rethrow()
        @warn "cost: no estimate for $(basename(input_path)): $(sprint(showerror, e))"
        return nothing
    end
    @info "Cost estimate: wall time $(_format_elapsed(estimate.wall_s)), " *
          "$(estimate.segments) segment(s), disk $(Base.format_bytes(estimate.disk_bytes)) " *
          "($(estimate.basis))"
    if estimate.exceeds
        sentence =
            "cost: the estimated wall time $(_format_elapsed(estimate.wall_s)) exceeds the " *
            "budget of $(cfg.cost.budget_hours) h (cost.budget_hours)"
        cfg.cost.on_exceed == "refuse" &&
            throw(ArgumentError(sentence * "; raise the budget or set cost.on_exceed = \"warn\""))
        @warn sentence
    end
    return estimate
end

"""
    _cost_fields(e::CostEstimate) -> Dict{String,Any}

The entries of an estimate in a segment record: `cost_wall_s` [s, to 0.1],
`cost_segments`, `cost_disk_bytes`, `cost_rate_s_per_nb` [s, to 0.01] and
`cost_basis`.
"""
function _cost_fields(e::CostEstimate)::Dict{String,Any}
    return Dict{String,Any}(
        "cost_wall_s" => round(e.wall_s; digits = 1),
        "cost_segments" => e.segments,
        "cost_disk_bytes" => e.disk_bytes,
        "cost_rate_s_per_nb" => round(e.rate_s_per_nb; digits = 2),
        "cost_basis" => e.basis,
    )
end
