# =============================================================================
# Block collector: ensembles, paired comparisons and block tables over a run list
# =============================================================================
#
# A block of the campaign is a set of runs submitted as independent chains (arms ×
# seeds), not a sweep directory. The collector takes the block as an explicit list of
# run directories, from the campaign ledger or from a sweep, and computes the products
# of the ensemble layer on it: percentile ensembles per arm, the paired comparison of
# every arm with the reference arm and the equivalence verdicts. It reads only the
# tier-1 files of a run (`output/lagr.7`, the engine's stdout capture,
# `output/merger_ic.toml`, the frozen `config.toml`), never dumps or snapshots.

"""
    sweep_runs(sweep_dir) -> Vector{NamedTuple}

The completed points of a sweep as a run list for [`block_comparison`](@ref) and
[`block_ensembles`](@ref): one `(; arm, seed, grid, run_id, run_dir)` per point, `arm`
`""` for a sweep without arms, `grid` the point's grid values as text (`"key=value"`
pairs sorted and joined by `", "`, `""` without a grid), `run_id` the point id.
"""
function sweep_runs(sweep_dir::AbstractString)::Vector{NamedTuple}
    _, pts = _sweep_done_points(sweep_dir)
    return NamedTuple[
        (;
            arm = string(get(p, "arm", ""))::String,
            seed = Int(p["seed"]),
            grid = _grid_text(p["values"]),
            run_id = string(p["id"])::String,
            run_dir,
        ) for (p, run_dir) in pts
    ]
end

"""
    _run_list_groups(runs, caller) -> _PairingGroups

Pairing groups of a run list: every run a named tuple with the fields `arm` (a
non-empty name), `seed` (an integer), `run_dir` (an existing directory) and optionally
`grid` (text, `""` when absent). Refuses an empty list, a missing field or directory
and two runs of one arm, seed and grid point, naming `caller`.
"""
function _run_list_groups(runs::AbstractVector{<:NamedTuple}, caller::String)::_PairingGroups
    isempty(runs) && throw(ArgumentError("$(caller): the run list is empty"))
    groups = _PairingGroups()
    for (k, r) in enumerate(runs)
        for f in (:arm, :seed, :run_dir)
            hasproperty(r, f) || throw(
                ArgumentError(
                    "$(caller): run $k has no field `$f`; a run is a named tuple " *
                    "(arm, seed, run_dir[, grid])",
                ),
            )
        end
        arm = string(r.arm)::String
        isempty(arm) && throw(ArgumentError("$(caller): run $k has an empty arm name"))
        r.seed isa Integer || throw(
            ArgumentError("$(caller): run $k has the seed $(repr(r.seed)); an integer is needed"),
        )
        seed = Int(r.seed)
        run_dir = string(r.run_dir)::String
        grid = hasproperty(r, :grid) ? string(r.grid)::String : ""
        where = isempty(grid) ? "" : " at $(grid)"
        isdir(run_dir) || throw(
            ArgumentError(
                "$(caller): the run directory $(run_dir) of arm $(arm), seed $(seed)$(where) does not exist",
            ),
        )
        byarm = get!(groups, (grid, seed), Dict{String,String}())
        haskey(byarm, arm) && throw(
            ArgumentError(
                "$(caller): two runs of arm $(arm), seed $(seed)$(where): $(byarm[arm]) and $(run_dir)",
            ),
        )
        byarm[arm] = run_dir
    end
    return groups
end

"""
    block_comparison(runs; reference = "base", observables = ["r10", "r50", "r90",
                     "energy", "n_stars", "n_pairs"], t_settle_myr = 0.0,
                     margins = Dict(), α = 0.05) -> Vector{NamedTuple}

The calibration table of [`paired_comparison`](@ref) over an explicit run list instead
of a sweep directory: `runs` is a vector of named tuples `(; arm, seed, run_dir)`, with
an optional `grid` (text naming the grid point; runs of different grid points are
never paired) — the members of a campaign block from [`block_runs`](@ref), or the
points of a sweep from [`sweep_runs`](@ref). Every run is paired with the run of the
arm `reference` of the same grid point and seed; the pairing rules (equal
`meta.dat10_sha256`, perturbed twins excepted), the statistics, the columns and Holm's
correction over the table are those of `paired_comparison`, which gives identical rows
for the run list of a sweep. Only the tier-1 files of a run are read. Refuses an empty
list, a missing field or directory and two runs of one arm, seed and grid point.

# Example
```julia
runs = block_runs("ledger.csv", "runs/campaign/A"; block = "A")
rows = block_comparison(runs; t_settle_myr = 1.0)
filter(r -> r.verdict != "equivalent", rows)
```
"""
function block_comparison(
    runs::AbstractVector{<:NamedTuple};
    reference::AbstractString = "base",
    observables::AbstractVector{<:AbstractString} = _PAIRED_DEFAULT_OBSERVABLES,
    t_settle_myr::Real = 0.0,
    margins::AbstractDict = Dict{String,Tuple{String,Float64}}(),
    α::Real = 0.05,
)::Vector{NamedTuple}
    _check_paired_observables(observables, "block_comparison")
    return _paired_rows(
        _run_list_groups(runs, "block_comparison");
        reference = String(reference),
        observables,
        t_settle_myr,
        margins,
        α,
        caller = "block_comparison",
        unit = "run",
    )
end

"""Ledger columns [`block_runs`](@ref) needs."""
const _BLOCK_LEDGER_COLUMNS = ("block", "arm", "seed", "run_id", "status")

"""
    block_runs(ledger_csv, runs_root; block, statuses = ["completed"]) -> Vector{NamedTuple}

The members of the campaign block `block` as a run list for
[`block_comparison`](@ref) and [`block_ensembles`](@ref): the rows of the ledger CSV
`ledger_csv` (columns `block`, `arm`, `seed`, `run_id`, `status` at least; quoted
fields may hold commas) whose `block` is `block` and whose `status` is one of
`statuses`, each resolved to the run directory `joinpath(runs_root, run_id)`. Returns
one `(; arm, seed, run_id, run_dir)` per row, in ledger order. Refuses an absent
ledger, a header without the needed columns, a row with another number of cells or a
seed that is not an integer, a selection without rows, two rows of one arm and seed
and a run directory that does not exist, naming the row.

# Example
```julia
runs = block_runs("artefacts/runs/h200-campaign/ledger.csv", "runs/campaign/A"; block = "A")
```
"""
function block_runs(
    ledger_csv::AbstractString,
    runs_root::AbstractString;
    block::AbstractString,
    statuses::AbstractVector{<:AbstractString} = ["completed"],
)::Vector{NamedTuple}
    isfile(ledger_csv) || throw(ArgumentError("block_runs: no ledger at $(ledger_csv)"))
    isempty(statuses) && throw(ArgumentError("block_runs: no status selected"))
    records = _csv_records(read(ledger_csv, String), ledger_csv)
    isempty(records) && throw(ArgumentError("block_runs: the ledger $(ledger_csv) is empty"))
    header = String.(strip.(records[1]))
    col = Dict(c => findfirst(==(c), header) for c in _BLOCK_LEDGER_COLUMNS)
    missing_cols = [c for c in _BLOCK_LEDGER_COLUMNS if col[c] === nothing]
    isempty(missing_cols) || throw(
        ArgumentError(
            "block_runs: the ledger $(ledger_csv) has no column " * join(missing_cols, ", "),
        ),
    )
    runs = NamedTuple[]
    seen = Set{Tuple{String,Int}}()
    for (k, cells) in enumerate(records[2:end])
        line = k + 1
        length(cells) == length(header) || throw(
            ArgumentError(
                "block_runs: record $(line) of $(ledger_csv) has $(length(cells)) cells for " *
                "$(length(header)) columns",
            ),
        )
        cell(c) = String(strip(cells[something(col[c])]))
        (cell("block") == block && cell("status") in statuses) || continue
        arm = cell("arm")
        run_id = cell("run_id")
        seed = tryparse(Int, cell("seed"))
        seed === nothing && throw(
            ArgumentError(
                "block_runs: record $(line) of $(ledger_csv) (run $(run_id)) has the seed " *
                "\"$(cell("seed"))\"; an integer is needed",
            ),
        )
        (isempty(arm) || isempty(run_id)) && throw(
            ArgumentError("block_runs: record $(line) of $(ledger_csv) has an empty arm or run_id"),
        )
        (arm, seed) in seen && throw(
            ArgumentError(
                "block_runs: block $(block) lists arm $(arm), seed $(seed) twice (record $(line))",
            ),
        )
        push!(seen, (arm, seed))
        run_dir = joinpath(runs_root, run_id)
        isdir(run_dir) || throw(
            ArgumentError(
                "block_runs: run $(run_id) of block $(block) (arm $(arm), seed $(seed); " *
                "record $(line) of $(ledger_csv)) has no directory $(run_dir)",
            ),
        )
        push!(runs, (; arm, seed, run_id, run_dir))
    end
    isempty(runs) && throw(
        ArgumentError(
            "block_runs: no run of block \"$(block)\" with status " *
            join(statuses, " or ") *
            " in $(ledger_csv)",
        ),
    )
    return runs
end

"""
    block_ensembles(runs; observables = ["r10", "r50", "r90", "energy", "n_stars",
                    "n_pairs"], n_grid = 200) -> Vector{NamedTuple}

Seeded ensembles per arm over a run list (`(; arm, seed, run_dir[, grid])`, as for
[`block_comparison`](@ref)): the runs of one arm and grid point are the members, and
the series of each observable ([`_run_series`](@ref), the observables of
[`paired_comparison`](@ref)) is summarised by [`ensemble_statistics`](@ref) — the
product of [`sweep_ensembles`](@ref) for a sweep, here kept apart per arm. Returns one
`(; grid, observable, arm, stats)` per grid point, observable and arm with at least one
series, in that order (grid points and arms by name).
"""
function block_ensembles(
    runs::AbstractVector{<:NamedTuple};
    observables::AbstractVector{<:AbstractString} = _PAIRED_DEFAULT_OBSERVABLES,
    n_grid::Int = _ENSEMBLE_GRID_NODES,
)::Vector{NamedTuple}
    _check_paired_observables(observables, "block_ensembles")
    groups = _run_list_groups(runs, "block_ensembles")
    arms = sort!(unique(a for g in values(groups) for a in keys(g)))
    out = NamedTuple[]
    for grid in sort!(unique(first.(collect(keys(groups))))), o in observables
        quantity, fraction = _PAIRED_OBSERVABLES[o]
        seeds, series = _arm_series(groups, grid, quantity, fraction)
        for arm in arms
            members = Tuple{Vector{Float64},Vector{Float64}}[
                (Float64.(series[(arm, s)][1]), Float64.(series[(arm, s)][2])) for
                s in seeds if haskey(series, (arm, s)) && series[(arm, s)] !== nothing
            ]
            isempty(members) && continue
            push!(out, (; grid, observable = o, arm, stats = ensemble_statistics(members; n_grid)))
        end
    end
    return out
end

"""
    write_block_ensembles(path, ensembles) -> String

Write the result of [`block_ensembles`](@ref) as CSV in long form — one line per grid
point, observable, arm and time node, columns `grid`, `observable`, `arm`, `n`,
`time_myr`, `median`, `q16`, `q84`, `q025`, `q975` — with the cell formatting of
[`write_paired_comparison`](@ref); an existing file is kept as a `#k` backup. Returns
`path`.
"""
function write_block_ensembles(path::AbstractString, ensembles::AbstractVector)::String
    isempty(ensembles) && throw(ArgumentError("write_block_ensembles: no ensembles"))
    _backup_existing(path)
    open(path, "w") do io
        println(io, "grid,observable,arm,n,time_myr,median,q16,q84,q025,q975")
        for e in ensembles
            s = e.stats::EnsembleStatistics
            for k in eachindex(s.time)
                cells = (
                    _csv_cell(e.grid),
                    _csv_cell(e.observable),
                    _csv_cell(e.arm),
                    _csv_cell(s.n),
                    _csv_cell(s.time[k]),
                    _csv_cell(s.median[k]),
                    _csv_cell(s.q16[k]),
                    _csv_cell(s.q84[k]),
                    _csv_cell(s.q025[k]),
                    _csv_cell(s.q975[k]),
                )
                println(io, join(cells, ","))
            end
        end
    end
    return path
end

"""
    write_block_tables(out_dir, runs; reference = "base", observables = [...],
                       t_settle_myr = 0.0, margins = Dict(), α = 0.05) -> NamedTuple

The block tables of a run list in `out_dir` (created when absent; existing files kept
as `#k` backups): `block_comparison.csv` ([`block_comparison`](@ref) through
[`write_paired_comparison`](@ref)), `block_ensembles.csv`
([`write_block_ensembles`](@ref)) and the member sidecar `block_members.csv` (`arm`,
`seed`, `grid`, `run_id`, `run_dir` of every run, in list order; `run_id` defaults to
the directory name). Returns the paths as `(; comparison, ensembles, members)`.
"""
function write_block_tables(
    out_dir::AbstractString,
    runs::AbstractVector{<:NamedTuple};
    reference::AbstractString = "base",
    observables::AbstractVector{<:AbstractString} = _PAIRED_DEFAULT_OBSERVABLES,
    t_settle_myr::Real = 0.0,
    margins::AbstractDict = Dict{String,Tuple{String,Float64}}(),
    α::Real = 0.05,
)
    rows = block_comparison(runs; reference, observables, t_settle_myr, margins, α)
    ensembles = block_ensembles(runs; observables)
    mkpath(out_dir)
    comparison = write_paired_comparison(joinpath(out_dir, "block_comparison.csv"), rows)
    ens_path = write_block_ensembles(joinpath(out_dir, "block_ensembles.csv"), ensembles)
    members = joinpath(out_dir, "block_members.csv")
    _backup_existing(members)
    open(members, "w") do io
        println(io, "arm,seed,grid,run_id,run_dir")
        for r in runs
            run_dir = string(r.run_dir)::String
            cells = (
                _csv_cell(string(r.arm)::String),
                _csv_cell(Int(r.seed)),
                _csv_cell(hasproperty(r, :grid) ? string(r.grid)::String : ""),
                _csv_cell(
                    hasproperty(r, :run_id) ? string(r.run_id)::String :
                    basename(normpath(run_dir)),
                ),
                _csv_cell(run_dir),
            )
            println(io, join(cells, ","))
        end
    end
    return (; comparison, ensembles = ens_path, members)
end
