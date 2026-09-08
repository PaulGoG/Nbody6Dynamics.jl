# =============================================================================
# Parameter sweeps: Cartesian grid × seeds over a merger configuration,
# one pipeline run per point, executed as concurrent worker processes
# =============================================================================
#
# A sweep TOML names a base pipeline config and a base merger config, a set
# of grid axes (dotted keys into the merger TOML, Cartesian product) and the
# seeds replicated at every grid point. Every point gets its own directory
# with the derived merger and pipeline configs, a worker process runs the
# full pipeline there (`run_sweep_point`), and the driver keeps
# `sweep_index.toml` current while the jobs run. Jobs are separate Julia
# processes so a crash in one point cannot take the sweep down, and the
# OpenMP thread count per job follows the cost model (4 threads saturate
# the backend for N ≲ 2×10⁴; five concurrent jobs on the reference
# workstation).

"""
    SweepConfig

Specification of a parameter sweep, parsed by [`load_sweep_config`](@ref).

# Fields
- `name`: sweep name (letters, digits, `_`, `-`); the sweep directory is
  `<runs_dir>/sweep_<name>_<timestamp>`
- `pipeline_config`, `merger_config`: absolute paths of the base pipeline and
  merger TOML files
- `grid`: axes `dotted.key => values` in key order (Cartesian product, first
  axis varying fastest); keys address the merger TOML from its root
  (`"merger.orbit.eccentricity"`) and their parent table must exist there
- `seeds`: `merger.seed` values replicated at every grid point
- `concurrency`: number of simultaneous jobs
- `omp_threads`: OpenMP threads per job
- `runs_dir`: absolute sweep root
- `poll_interval`: seconds between checks of the running jobs
"""
struct SweepConfig
    name::String
    pipeline_config::String
    merger_config::String
    grid::Vector{Pair{String,Vector{Any}}}
    seeds::Vector{Int}
    concurrency::Int
    omp_threads::Int
    runs_dir::String
    poll_interval::Float64
end

"""
    SweepPoint

One point of a sweep: its 1-based `index`, directory `id`, the axis values
in grid order, and the seed.
"""
struct SweepPoint
    index::Int
    id::String
    values::Vector{Pair{String,Any}}
    seed::Int
end

const _SWEEP_NAME_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9_-]*$"
const _SWEEP_INDEX_FILE = "sweep_index.toml"
const _SWEEP_SUMMARY_FILE = "sweep_summary.csv"
const _SWEEP_RUN_ID = "run"

"""
    load_sweep_config(path::AbstractString) -> SweepConfig

Parse and validate a sweep TOML:

```toml
[sweep]
name = "demo"
pipeline_config = "../config.toml"           # relative to this file
merger_config = "merger_demo_small.toml"
seeds = [11, 12]
concurrency = 5       # ≥ 1
omp_threads = 4       # ≥ 1
runs_dir = "../runs"  # sweep root
[sweep.grid]
"merger.orbit.eccentricity" = [0.0, 0.6]
"merger.cluster2.N" = [500, 1000]
```

Axes must address the merger TOML from its root and may not be
`merger.seed`; seeds must be distinct integers.
"""
function load_sweep_config(path::AbstractString)::SweepConfig
    isfile(path) || error("Sweep configuration not found: $path")
    raw = TOML.parsefile(path)
    haskey(raw, "sweep") || error("sweep config: missing [sweep] table in $path")
    s = raw["sweep"]
    base = dirname(abspath(path))
    resolve(p) = isabspath(p) ? String(p) : normpath(joinpath(base, p))

    name = String(get(s, "name", ""))
    occursin(_SWEEP_NAME_PATTERN, name) ||
        error("sweep config: name must match $(_SWEEP_NAME_PATTERN.pattern), got \"$name\"")
    pipeline_config = resolve(String(get(s, "pipeline_config", "")))
    isfile(pipeline_config) || error("sweep config: pipeline_config not found: $pipeline_config")
    merger_config = resolve(String(get(s, "merger_config", "")))
    isfile(merger_config) || error("sweep config: merger_config not found: $merger_config")

    grid_raw = get(s, "grid", Dict{String,Any}())
    grid_raw isa AbstractDict ||
        error("sweep config: [sweep.grid] must be a table of key = [values]")
    isempty(grid_raw) && error("sweep config: [sweep.grid] must define at least one axis")
    grid = Pair{String,Vector{Any}}[]
    for key in sort!(collect(String, keys(grid_raw)))
        vals = grid_raw[key]
        (vals isa AbstractVector && !isempty(vals)) ||
            error("sweep config: grid axis \"$key\" must be a nonempty array")
        all(v -> v isa Real || v isa AbstractString, vals) ||
            error("sweep config: grid axis \"$key\" must hold numbers, booleans or strings")
        startswith(key, "merger.") || error(
            "sweep config: grid axis \"$key\" must address the merger TOML from its root (\"merger.…\")",
        )
        key == "merger.seed" &&
            error("sweep config: seeds are given by [sweep] seeds, not as a grid axis")
        push!(grid, key => Any[vals...])
    end

    seeds_raw = get(s, "seeds", Any[])
    (seeds_raw isa AbstractVector && !isempty(seeds_raw) && all(v -> v isa Integer, seeds_raw)) ||
        error("sweep config: seeds must be a nonempty array of integers")
    seeds = Int[seeds_raw...]
    allunique(seeds) || error("sweep config: seeds must be distinct, got $seeds")

    concurrency = Int(get(s, "concurrency", 5))
    concurrency ≥ 1 || error("sweep config: concurrency must be ≥ 1, got $concurrency")
    omp_threads = Int(get(s, "omp_threads", 4))
    omp_threads ≥ 1 || error("sweep config: omp_threads must be ≥ 1, got $omp_threads")
    poll_interval = Float64(get(s, "poll_interval", 2.0))
    poll_interval > 0 || error("sweep config: poll_interval must be positive, got $poll_interval")
    runs_dir = resolve(String(get(s, "runs_dir", "runs")))

    return SweepConfig(
        name,
        pipeline_config,
        merger_config,
        grid,
        seeds,
        concurrency,
        omp_threads,
        runs_dir,
        poll_interval,
    )
end

"""Axis label used in point ids and legends: the dotted key without its
`merger.` root, components joined by `-` (`cluster2-N`)."""
function _axis_short(key::AbstractString)
    parts = split(key, '.')
    parts[1] == "merger" && length(parts) > 1 && (parts = parts[2:end])
    return join(parts, "-")
end

"""Value rendering for point ids: integers verbatim, reals with `%g`,
strings with every character outside `[A-Za-z0-9.+-]` replaced by `-`."""
function _format_axis_value(v)
    v isa Bool && return string(v)
    v isa Integer && return string(v)
    v isa Real && return @sprintf("%g", v)
    return replace(String(v), r"[^A-Za-z0-9.+-]" => "-")
end

"""
    sweep_points(cfg::SweepConfig) -> Vector{SweepPoint}

Enumerate the Cartesian product of the grid axes (first axis varying
fastest) times the seeds. Ids read
`<index>_<axis>=<value>_…_seed=<seed>`.
"""
function sweep_points(cfg::SweepConfig)::Vector{SweepPoint}
    axis_keys = first.(cfg.grid)
    value_lists = last.(cfg.grid)
    points = SweepPoint[]
    index = 0
    for combo in Iterators.product(value_lists...)
        values = Pair{String,Any}[k => v for (k, v) in zip(axis_keys, combo)]
        label = join(["$(_axis_short(k))=$(_format_axis_value(v))" for (k, v) in values], "_")
        for seed in cfg.seeds
            index += 1
            push!(
                points,
                SweepPoint(index, @sprintf("%03d_%s_seed=%d", index, label, seed), values, seed),
            )
        end
    end
    return points
end

"""
    _set_nested!(d, path, value)

Assign `value` at the dotted `path` of the nested dictionary `d`. Every
table on the path except the leaf must already exist, so a mistyped axis
fails before anything runs.
"""
function _set_nested!(d::AbstractDict, path::AbstractString, value)
    parts = split(path, '.')
    node = d
    for i in 1:(length(parts) - 1)
        p = parts[i]
        (haskey(node, p) && node[p] isa AbstractDict) || throw(
            ArgumentError(
                "sweep axis \"$path\": table \"$(join(parts[1:i], '.'))\" is not present in the base configuration",
            ),
        )
        node = node[p]
    end
    node[parts[end]] = value
    return d
end

"""Pipeline-config overrides of a sweep point: no install phase, absolute
backend path, the point directory as run root, the derived merger TOML,
and the sweep's thread count."""
function _sweep_pipeline_overrides!(
    c::AbstractDict,
    cfg::SweepConfig,
    pipeline_dir::AbstractString,
    point_dir::AbstractString,
    merger_path::AbstractString,
)
    inst = get!(c, "install", Dict{String,Any}())
    inst["enabled"] = false
    install_dir = String(get(inst, "install_dir", "backend/Nbody6PPGPU-beijing"))
    inst["install_dir"] =
        isabspath(install_dir) ? install_dir : abspath(joinpath(pipeline_dir, install_dir))
    sim = get!(c, "simulation", Dict{String,Any}())
    sim["run_test"] = true
    sim["runs_dir"] = abspath(point_dir)
    sim["omp_threads"] = cfg.omp_threads
    sim["monitor"] = false
    pp = get!(c, "postprocess", Dict{String,Any}())
    pp["data_dir"] = ""
    mg = get!(c, "merger", Dict{String,Any}())
    mg["enabled"] = true
    mg["config_file"] = abspath(merger_path)
    return c
end

"""
    prepare_sweep(cfg::SweepConfig; sweep_dir = "") -> (sweep_dir, points)

Create the sweep directory (`<runs_dir>/sweep_<name>_<timestamp>` unless
given; an existing directory is never reused) and, per point, a
subdirectory holding the derived `merger.toml` and `config.toml`. Both
derived files are loaded back through the regular parsers, so an invalid
point fails here rather than in a worker. Writes the initial
`sweep_index.toml` with every point `pending`.
"""
function prepare_sweep(cfg::SweepConfig; sweep_dir::AbstractString = "")
    if isempty(sweep_dir)
        stamp = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
        sweep_dir = joinpath(cfg.runs_dir, "sweep_$(cfg.name)_$(stamp)")
    end
    sweep_dir = abspath(sweep_dir)
    isdir(sweep_dir) && error("sweep directory already exists: $sweep_dir")
    mkpath(sweep_dir)

    merger_base = TOML.parsefile(cfg.merger_config)
    pipeline_base = TOML.parsefile(cfg.pipeline_config)
    pipeline_dir = dirname(cfg.pipeline_config)
    points = sweep_points(cfg)

    for p in points
        pdir = joinpath(sweep_dir, p.id)
        mkpath(pdir)
        m = deepcopy(merger_base)
        for (k, v) in p.values
            _set_nested!(m, k, v)
        end
        _set_nested!(m, "merger.seed", p.seed)
        merger_path = joinpath(pdir, "merger.toml")
        open(io -> TOML.print(io, m), merger_path, "w")
        load_merger_config(merger_path)

        c = deepcopy(pipeline_base)
        _sweep_pipeline_overrides!(c, cfg, pipeline_dir, pdir, merger_path)
        config_path = joinpath(pdir, "config.toml")
        open(io -> TOML.print(io, c), config_path, "w")
        load_config(config_path)
    end

    cp(cfg.pipeline_config, joinpath(sweep_dir, "base_config.toml"))
    cp(cfg.merger_config, joinpath(sweep_dir, "base_merger.toml"))
    write_sweep_index(sweep_dir, cfg, points)
    return sweep_dir, points
end

"""
    write_sweep_index(sweep_dir, cfg, points; status = Dict())

Write `sweep_index.toml`: the sweep header (name, base configs, axes,
seeds, job settings) and one `[[points]]` entry per point with its id,
index, seed, axis values, directory and status (`pending`, `running`,
`done`, `failed`; finished points carry `exit_status` and
`elapsed_seconds`). `status` maps point indices to the fields that
override the pending defaults.
"""
function write_sweep_index(
    sweep_dir::AbstractString,
    cfg::SweepConfig,
    points::Vector{SweepPoint};
    status::Dict{Int,Dict{String,Any}} = Dict{Int,Dict{String,Any}}(),
)
    header = Dict{String,Any}(
        "name" => cfg.name,
        "created" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        "pipeline_config" => cfg.pipeline_config,
        "merger_config" => cfg.merger_config,
        "axes" => first.(cfg.grid),
        "seeds" => cfg.seeds,
        "concurrency" => cfg.concurrency,
        "omp_threads" => cfg.omp_threads,
        "n_points" => length(points),
    )
    entries = map(points) do p
        e = Dict{String,Any}(
            "index" => p.index,
            "id" => p.id,
            "seed" => p.seed,
            "dir" => joinpath(sweep_dir, p.id),
            "values" => Dict{String,Any}(k => v for (k, v) in p.values),
            "status" => "pending",
        )
        merge!(e, get(status, p.index, Dict{String,Any}()))
        e
    end
    d = Dict{String,Any}("sweep" => header, "points" => entries)
    open(io -> TOML.print(io, d), joinpath(sweep_dir, _SWEEP_INDEX_FILE), "w")
    return nothing
end

"""
    read_sweep_index(sweep_dir) -> Dict{String,Any}

Parse `sweep_index.toml` (`"sweep"` header and `"points"` list, sorted by
index).
"""
function read_sweep_index(sweep_dir::AbstractString)
    path = joinpath(sweep_dir, _SWEEP_INDEX_FILE)
    isfile(path) || error("sweep index not found: $path")
    idx = TOML.parsefile(path)
    haskey(idx, "sweep") && haskey(idx, "points") || error("malformed sweep index: $path")
    sort!(idx["points"]; by = p -> p["index"])
    return idx
end

"""
    run_sweep_point(point_dir::AbstractString)

Worker entry: run the full pipeline of one sweep point from its
`config.toml`, with the run directory fixed to `<point_dir>/run`.
"""
function run_sweep_point(point_dir::AbstractString)
    cfg = load_config(joinpath(point_dir, "config.toml"))
    run_pipeline(cfg; base_dir = point_dir, run_id = _SWEEP_RUN_ID)
    return nothing
end

"""Command of a worker process for `point_dir`: the current Julia with the
package project, calling [`run_sweep_point`](@ref)."""
function _sweep_worker_command(point_dir::AbstractString)
    julia = Base.julia_cmd()
    return `$julia --project=$(_PROJECT_ROOT) --startup-file=no -e "using Nbody6Dynamics; run_sweep_point(ARGS[1])" $point_dir`
end

"""Fail before launching anything when the backend binary of the base
pipeline config cannot be found."""
function _check_sweep_binary(cfg::SweepConfig)
    base = load_config(cfg.pipeline_config)
    install_dir = base.install.install_dir
    src_dir =
        isabspath(install_dir) ? install_dir : joinpath(dirname(cfg.pipeline_config), install_dir)
    return _find_binary(src_dir, base.simulation.binary_name)
end

"""
    _launch_sweep(sweep_dir, cfg, points) -> Dict{Int,Dict{String,Any}}

Run the points as worker processes, at most `cfg.concurrency` at a time,
logging each to `<point>/sweep_point.log` and rewriting the index on every
state change. Returns the per-point status table.
"""
function _launch_sweep(sweep_dir::AbstractString, cfg::SweepConfig, points::Vector{SweepPoint})
    status = Dict{Int,Dict{String,Any}}()
    queue = copy(points)
    running = Dict{Int,Tuple{Base.Process,IOStream,Float64,SweepPoint}}()
    n_done = 0
    n_failed = 0
    while !isempty(queue) || !isempty(running)
        while length(running) < cfg.concurrency && !isempty(queue)
            p = popfirst!(queue)
            pdir = joinpath(sweep_dir, p.id)
            log = open(joinpath(pdir, "sweep_point.log"), "w")
            proc =
                run(pipeline(_sweep_worker_command(pdir); stdout = log, stderr = log); wait = false)
            running[p.index] = (proc, log, time(), p)
            status[p.index] = Dict{String,Any}(
                "status" => "running",
                "started" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
            )
            write_sweep_index(sweep_dir, cfg, points; status = status)
            @info "sweep: started $(p.id) ($(length(running)) running, $(length(queue)) queued)"
        end
        sleep(cfg.poll_interval)
        for (i, (proc, log, t0, p)) in collect(running)
            process_running(proc) && continue
            close(log)
            code = _exit_status(proc)
            ok = code == 0
            ok ? (n_done += 1) : (n_failed += 1)
            status[i] = Dict{String,Any}(
                "status" => ok ? "done" : "failed",
                "started" => status[i]["started"],
                "exit_status" => code,
                "elapsed_seconds" => round(time() - t0; digits = 1),
            )
            delete!(running, i)
            write_sweep_index(sweep_dir, cfg, points; status = status)
            if ok
                @info "sweep: finished $(p.id) in $(_format_elapsed(time() - t0))"
            else
                @warn "sweep: $(p.id) failed (exit status $code); see $(joinpath(sweep_dir, p.id, "sweep_point.log"))"
            end
        end
    end
    @info "sweep: $(n_done) done, $(n_failed) failed, index at $(joinpath(sweep_dir, _SWEEP_INDEX_FILE))"
    return status
end

"""
    run_sweep(cfg::SweepConfig; dry_run = false, sweep_dir = "") -> String

Prepare the sweep ([`prepare_sweep`](@ref)), run every point as a worker
process with at most `cfg.concurrency` jobs at a time, then write
`sweep_summary.csv` ([`write_sweep_summary`](@ref)). With `dry_run = true`
the directories, derived configs and index are written and nothing is
launched. Returns the sweep directory.
"""
function run_sweep(cfg::SweepConfig; dry_run::Bool = false, sweep_dir::AbstractString = "")
    dry_run || _check_sweep_binary(cfg)
    sweep_dir, points = prepare_sweep(cfg; sweep_dir = sweep_dir)
    @info "sweep \"$(cfg.name)\": $(length(points)) points in $sweep_dir" *
          (dry_run ? " (dry run, nothing launched)" : "")
    dry_run && return sweep_dir
    _launch_sweep(sweep_dir, cfg, points)
    write_sweep_summary(sweep_dir)
    return sweep_dir
end

"""Outcome fields of one finished point read from its run directory:
elapsed time and exit status from `RUN_INFO.toml`, final time, star and
pair counts, energy error and virial ratio from the last ADJUST record.
Missing files leave `NaN`/`-1` entries."""
function _sweep_point_outcome(run_dir::AbstractString)
    out = Dict{String,Any}(
        "elapsed_seconds" => NaN,
        "exit_status" => -1,
        "t_final_myr" => NaN,
        "n_final" => -1,
        "npairs_final" => -1,
        "de_final" => NaN,
        "q_final" => NaN,
    )
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    if isfile(info_path)
        info = TOML.parsefile(info_path)
        out["elapsed_seconds"] = Float64(get(get(info, "run", Dict()), "elapsed_seconds", NaN))
        segs = get(info, "segments", Any[])
        isempty(segs) || (out["exit_status"] = Int(get(segs[end], "exit_status", -1)))
    end
    stdout_path = joinpath(run_dir, "output", "out1000")
    if isfile(stdout_path)
        diag = read_diagnostics(stdout_path)
        if !isempty(diag.adjust)
            a = diag.adjust[end]
            out["t_final_myr"] = a.time_myr
            out["n_final"] = a.n
            out["npairs_final"] = a.npairs
            out["de_final"] = a.de_rel
            out["q_final"] = a.qvir
        end
    end
    return out
end

"""
    sweep_summary(sweep_dir) -> (columns, rows)

One row per point: index, id, seed, the axis values, status, and the
outcome fields of [`_sweep_point_outcome`](@ref) for finished points.
`columns` gives the column order used by [`write_sweep_summary`](@ref).
"""
function sweep_summary(sweep_dir::AbstractString)
    idx = read_sweep_index(sweep_dir)
    axes = String[idx["sweep"]["axes"]...]
    columns = vcat(
        ["index", "id", "seed"],
        axes,
        [
            "status",
            "exit_status",
            "elapsed_seconds",
            "t_final_myr",
            "n_final",
            "npairs_final",
            "de_final",
            "q_final",
        ],
    )
    rows = Vector{Dict{String,Any}}()
    for p in idx["points"]
        row = Dict{String,Any}(
            "index" => p["index"],
            "id" => p["id"],
            "seed" => p["seed"],
            "status" => p["status"],
        )
        for a in axes
            row[a] = p["values"][a]
        end
        outcome = _sweep_point_outcome(joinpath(p["dir"], _SWEEP_RUN_ID))
        merge!(row, outcome)
        haskey(p, "elapsed_seconds") && (row["elapsed_seconds"] = p["elapsed_seconds"])
        haskey(p, "exit_status") && (row["exit_status"] = p["exit_status"])
        push!(rows, row)
    end
    return columns, rows
end

_csv_cell(v::AbstractString) =
    occursin(r"[,\"\n]", v) ? "\"" * replace(v, "\"" => "\"\"") * "\"" : v
_csv_cell(v::Bool) = string(v)
_csv_cell(v::Integer) = string(v)
_csv_cell(v::Real) = isfinite(v) ? @sprintf("%.6g", v) : ""
_csv_cell(v) = _csv_cell(string(v))

"""
    write_sweep_summary(sweep_dir) -> String

Write [`sweep_summary`](@ref) to `sweep_summary.csv` in `sweep_dir`
(non-finite numbers as empty cells) and return the path.
"""
function write_sweep_summary(sweep_dir::AbstractString)
    columns, rows = sweep_summary(sweep_dir)
    path = joinpath(sweep_dir, _SWEEP_SUMMARY_FILE)
    open(path, "w") do io
        println(io, join(_csv_cell.(columns), ","))
        for r in rows
            println(io, join([_csv_cell(get(r, c, "")) for c in columns], ","))
        end
    end
    @info "Sweep summary written: $path"
    return path
end

"""
    sweep_visualization(cfg::SweepConfig, sweep_dir) -> VisualizationConfig

The base pipeline config's `[visualization]` settings with `output_dir`
pointing at `<sweep_dir>/plots` (every other field forwarded unchanged).
"""
function sweep_visualization(cfg::SweepConfig, sweep_dir::AbstractString)
    base = load_config(cfg.pipeline_config).visualization
    fields = Dict{Symbol,Any}(k => getfield(base, k) for k in fieldnames(VisualizationConfig))
    fields[:output_dir] = joinpath(sweep_dir, "plots")
    return VisualizationConfig(; fields...)
end
