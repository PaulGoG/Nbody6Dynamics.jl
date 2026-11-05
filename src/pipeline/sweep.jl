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
# workstation). With `devices` set, every point is pinned to one CUDA device
# (`simulation.gpu_list` of its config) and at most `engines_per_device`
# engines share a device, in addition to the `concurrency` cap: the engine
# is host-bound at large N, so several engines per device are the norm.

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
- `controls`: whether every point gets an isolated single-cluster control
  companion ([`control_merger_dict`](@ref))
- `devices`: CUDA device indices (as `nvidia-smi` numbers them) the points
  are distributed over, one device per point; empty for no device
  assignment (the base config's `simulation.gpu_list` is then untouched)
- `engines_per_device`: engines admitted on one device at a time (≥ 1)
- `arms`: named bundles of overrides, `name => [dotted key => value, …]`, in
  name order; keys address the merger TOML (`merger.…`) or the point's
  pipeline TOML (`pipeline.…`). Every arm runs every grid point and seed with
  the same initial conditions (the same seed), so the arms are paired: of the
  merger keys only `merger.nbody6.*`, `merger.stellar.*`, `merger.mbh.*`, the
  output intervals and the twin keys may be set, and no structural pipeline
  key (`simulation.runs_dir`, `merger.config_file`, `merger.enabled`). Empty
  for a sweep without arms
- `max_segments`: engine segments a point may hold on record before it is
  left to the next entry: a point stopped by its wall budget, with progress
  and a worker that exited cleanly, is launched again within the entry
  ([`continue_sweep_point`](@ref)) until it completes or its record holds this
  many segments (≥ 1; 1 leaves a stopped point to the next entry, as without
  the key)
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
    controls::Bool
    devices::Vector{Int}
    engines_per_device::Int
    arms::Vector{Pair{String,Vector{Pair{String,Any}}}}
    max_segments::Int
end

# The twelve-field form of earlier versions (before arms and max_segments).
SweepConfig(
    name,
    pipeline_config,
    merger_config,
    grid,
    seeds,
    concurrency,
    omp_threads,
    runs_dir,
    poll_interval,
    controls,
    devices,
    engines_per_device,
) = SweepConfig(
    name,
    pipeline_config,
    merger_config,
    grid,
    seeds,
    concurrency,
    omp_threads,
    runs_dir,
    poll_interval,
    controls,
    devices,
    engines_per_device,
    Pair{String,Vector{Pair{String,Any}}}[],
    1,
)

"""
    SweepPoint

One point of a sweep: its 1-based `index`, directory `id`, the axis values
in grid order, the seed, its `kind` (`"merger"` or `"control"`), for a
control the id of the merger point it belongs to (`control_of`), and the
name of its arm (`""` without arms).
"""
struct SweepPoint
    index::Int
    id::String
    values::Vector{Pair{String,Any}}
    seed::Int
    kind::String
    control_of::String
    arm::String
end
SweepPoint(index, id, values, seed) = SweepPoint(index, id, values, seed, "merger", "", "")
SweepPoint(index, id, values, seed, kind, control_of) =
    SweepPoint(index, id, values, seed, kind, control_of, "")

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
controls = false      # isolated single-cluster control per point
devices = [0, 1]      # CUDA device indices, distinct, ≥ 0; [] for none
engines_per_device = 4  # ≥ 1
[sweep.grid]
"merger.orbit.eccentricity" = [0.0, 0.6]
"merger.cluster2.N" = [500, 1000]
```

Axes must address the merger TOML from its root and may not be
`merger.seed`; seeds must be distinct integers. A sweep without grid axes
is a seed ensemble of the base merger configuration. A nonempty `devices`
requires `build.enable_gpu = true` in the base pipeline config.
"""
function load_sweep_config(path::AbstractString)::SweepConfig
    isfile(path) || error("Sweep configuration not found: $path")
    raw = TOML.parsefile(path)
    _reject_unknown(raw, ("sweep",), basename(path))
    haskey(raw, "sweep") || _config_error("sweep", "table is missing in $path")
    s = raw["sweep"]
    _reject_unknown(
        s,
        (
            "name",
            "pipeline_config",
            "merger_config",
            "grid",
            "seeds",
            "concurrency",
            "omp_threads",
            "poll_interval",
            "runs_dir",
            "controls",
            "devices",
            "engines_per_device",
            "arms",
            "max_segments",
        ),
        "sweep",
    )
    base = dirname(abspath(path))
    resolve(p) = isabspath(p) ? String(p) : normpath(joinpath(base, p))

    name = _typed(get(s, "name", ""), String, "sweep.name")
    occursin(_SWEEP_NAME_PATTERN, name) ||
        _config_error("sweep.name", "must match $(_SWEEP_NAME_PATTERN.pattern), got \"$name\"")
    pipeline_config =
        resolve(_typed(get(s, "pipeline_config", ""), String, "sweep.pipeline_config"))
    isfile(pipeline_config) || _config_error("sweep.pipeline_config", "not found: $pipeline_config")
    merger_config = resolve(_typed(get(s, "merger_config", ""), String, "sweep.merger_config"))
    isfile(merger_config) || _config_error("sweep.merger_config", "not found: $merger_config")

    grid_raw = get(s, "grid", Dict{String,Any}())
    grid_raw isa AbstractDict || _config_error("sweep.grid", "must be a table of key = [values]")
    grid = Pair{String,Vector{Any}}[]
    for key in sort!(collect(String, keys(grid_raw)))
        vals = grid_raw[key]
        (vals isa AbstractVector && !isempty(vals)) ||
            _config_error("sweep.grid.$key", "must be a nonempty array")
        all(v -> v isa Real || v isa AbstractString, vals) ||
            _config_error("sweep.grid.$key", "must hold numbers, booleans or strings")
        startswith(key, "merger.") || _config_error(
            "sweep.grid.$key",
            "must address the merger TOML from its root (\"merger.…\")",
        )
        key == "merger.seed" &&
            _config_error("sweep.grid.$key", "is not an axis: seeds are given by [sweep] seeds")
        push!(grid, key => Any[vals...])
    end

    seeds_raw = get(s, "seeds", Any[])
    (seeds_raw isa AbstractVector && !isempty(seeds_raw) && all(v -> v isa Integer, seeds_raw)) ||
        _config_error("sweep.seeds", "must be a nonempty array of integers")
    seeds = Int[seeds_raw...]
    allunique(seeds) || _config_error("sweep.seeds", "must be distinct, got $seeds")

    concurrency = _typed(get(s, "concurrency", 5), Int, "sweep.concurrency")
    concurrency ≥ 1 || _config_error("sweep.concurrency", "must be ≥ 1, got $concurrency")
    omp_threads = _typed(get(s, "omp_threads", 4), Int, "sweep.omp_threads")
    omp_threads ≥ 1 || _config_error("sweep.omp_threads", "must be ≥ 1, got $omp_threads")
    poll_interval = _typed(get(s, "poll_interval", 2.0), Float64, "sweep.poll_interval")
    poll_interval > 0 ||
        _config_error("sweep.poll_interval", "must be positive, got $poll_interval")
    runs_dir = resolve(_typed(get(s, "runs_dir", "runs"), String, "sweep.runs_dir"))
    controls = _typed(get(s, "controls", false), Bool, "sweep.controls")

    devices = _typed(get(s, "devices", Any[]), Vector{Int}, "sweep.devices")
    all(≥(0), devices) ||
        _config_error("sweep.devices", "entries must be ≥ 0 (CUDA device indices), got $devices")
    allunique(devices) || _config_error("sweep.devices", "must be distinct, got $devices")
    engines_per_device = _typed(get(s, "engines_per_device", 1), Int, "sweep.engines_per_device")
    engines_per_device ≥ 1 ||
        _config_error("sweep.engines_per_device", "must be ≥ 1, got $engines_per_device")
    if !isempty(devices)
        pipeline_raw = TOML.parsefile(pipeline_config)
        get(get(pipeline_raw, "build", Dict{String,Any}()), "enable_gpu", false) === true ||
            _config_error(
                "sweep.devices",
                "requires build.enable_gpu = true in the pipeline config $pipeline_config: " *
                "the CPU binary ignores GPU_LIST",
            )
    end

    arms = _parse_sweep_arms(get(s, "arms", Dict{String,Any}()), Set(first.(grid)))
    max_segments = _typed(get(s, "max_segments", 1), Int, "sweep.max_segments")
    max_segments ≥ 1 || _config_error("sweep.max_segments", "must be ≥ 1, got $max_segments")

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
        controls,
        devices,
        engines_per_device,
        arms,
        max_segments,
    )
end

"""
    _parse_sweep_arms(raw, grid_keys) -> Vector{Pair{String,Vector{Pair{String,Any}}}}

The `[sweep.arms.<name>]` tables in name order. Every key is a quoted dotted key that
addresses the merger TOML (`"merger.…"`) or the pipeline TOML (`"pipeline.…"`). Refused: the seed
and the grid axes (an arm changes the configuration at every grid point and seed alike),
every merger key the generator reads (`_ARM_MERGER_PREFIXES`, `_ARM_MERGER_KEYS`: the arms
share their initial conditions), the structural pipeline keys (`_ARM_PIPELINE_REFUSED`) and
an unquoted dotted key, which TOML turns into nested tables.
"""
# An arm shares the initial conditions with the other arms, so it may set only merger keys
# the generator does not read: the engine parameters, the stellar-evolution options, the
# static massive black hole, the output intervals and the twin keys.
const _ARM_MERGER_PREFIXES = ("merger.nbody6.", "merger.stellar.", "merger.mbh.")
const _ARM_MERGER_KEYS = (
    "merger.output.tcrit",
    "merger.output.dtadj",
    "merger.output.deltat",
    "merger.output.tcrit_myr",
    "merger.output.dtadj_myr",
    "merger.output.deltat_myr",
    "merger.output.twin_amplitude",
    "merger.output.twin_body",
)
# Structural keys of the sweep: the runs directory, the merger file and the merger switch.
const _ARM_PIPELINE_REFUSED =
    ("pipeline.simulation.runs_dir", "pipeline.merger.config_file", "pipeline.merger.enabled")

function _parse_sweep_arms(raw, grid_keys::AbstractSet{<:AbstractString})
    raw isa AbstractDict ||
        _config_error("sweep.arms", "must be a table of [sweep.arms.<name>] tables")
    arms = Pair{String,Vector{Pair{String,Any}}}[]
    for name in sort!(collect(String, keys(raw)))
        occursin(_SWEEP_NAME_PATTERN, name) || _config_error(
            "sweep.arms.$name",
            "the arm name must match $(_SWEEP_NAME_PATTERN.pattern)",
        )
        table = raw[name]
        table isa AbstractDict ||
            _config_error("sweep.arms.$name", "must be a table of \"<dotted key>\" = value")
        items = Pair{String,Any}[]
        for key in sort!(collect(String, keys(table)))
            value = table[key]
            value isa AbstractDict && _config_error(
                "sweep.arms.$name.$key",
                "is a table: quote the dotted key (\"merger.….…\" = value)",
            )
            (startswith(key, "merger.") || startswith(key, "pipeline.")) || _config_error(
                "sweep.arms.$name.$key",
                "must address the merger TOML (\"merger.…\") or the pipeline TOML (\"pipeline.…\")",
            )
            key == "merger.seed" && _config_error(
                "sweep.arms.$name.$key",
                "is not an arm key: seeds are given by [sweep] seeds",
            )
            key in grid_keys &&
                _config_error("sweep.arms.$name.$key", "is also a grid axis; an arm may not set it")
            if startswith(key, "merger.")
                any(p -> startswith(key, p), _ARM_MERGER_PREFIXES) ||
                    key in _ARM_MERGER_KEYS ||
                    _config_error(
                        "sweep.arms.$name.$key",
                        "changes the initial conditions, which every arm shares; an arm may set " *
                        "merger.nbody6.*, merger.stellar.*, merger.mbh.*, the output intervals " *
                        "(merger.output.tcrit, dtadj, deltat and their _myr forms) and the twin keys",
                    )
            elseif key in _ARM_PIPELINE_REFUSED
                _config_error(
                    "sweep.arms.$name.$key",
                    "is a structural key of the sweep (the runs directory, the merger file, the " *
                    "merger switch); an arm may not set it",
                )
            end
            push!(items, key => value)
        end
        push!(arms, name => items)
    end
    return arms
end

"""Overrides of the arm `name` of `cfg` (empty for `""`, a sweep without arms)."""
function _arm_items(cfg::SweepConfig, name::AbstractString)::Vector{Pair{String,Any}}
    isempty(name) && return Pair{String,Any}[]
    k = findfirst(a -> first(a) == name, cfg.arms)
    k === nothing && throw(ArgumentError("sweep: no arm named \"$name\""))
    return last(cfg.arms[k])
end

"""
    _set_nested_create!(d, path, value)

Assign `value` at the dotted `path` of `d`, creating the tables on the way (an arm may
introduce a table the base file lacks, such as `[merger.output]` keys; the regular
parsers that load the derived file refuse a mistyped key).
"""
function _set_nested_create!(d::AbstractDict, path::AbstractString, value)
    parts = split(path, '.')
    node = d
    for p in parts[1:(end - 1)]
        child = get!(node, String(p), Dict{String,Any}())
        child isa AbstractDict ||
            throw(ArgumentError("sweep arm key \"$path\": \"$p\" is a value, not a table"))
        node = child
    end
    node[String(parts[end])] = value
    return d
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
fastest) times the seeds. Ids read `<index>_<axis>=<value>_…_seed=<seed>`.
With `cfg.controls` every merger point is followed by its control companion
(`<index>_<…>_control`, `kind = "control"`, `control_of` the merger id),
sharing the axis values and seed.
"""
function sweep_points(cfg::SweepConfig)::Vector{SweepPoint}
    axis_keys = first.(cfg.grid)
    value_lists = last.(cfg.grid)
    points = SweepPoint[]
    index = 0
    arm_names = isempty(cfg.arms) ? [""] : first.(cfg.arms)
    point_id(index, arm, label, seed) =
        @sprintf("%03d_", index) *
        join(filter(!isempty, ["$(isempty(arm) ? "" : "arm=" * arm)", label, "seed=$seed"]), "_")
    for arm in arm_names, combo in Iterators.product(value_lists...)
        values = Pair{String,Any}[k => v for (k, v) in zip(axis_keys, combo)]
        label = join(["$(_axis_short(k))=$(_format_axis_value(v))" for (k, v) in values], "_")
        for seed in cfg.seeds
            index += 1
            id = point_id(index, arm, label, seed)
            push!(points, SweepPoint(index, id, values, seed, "merger", "", arm))
            if cfg.controls
                index += 1
                cid = point_id(index, arm, label, seed) * "_control"
                push!(points, SweepPoint(index, cid, values, seed, "control", id, arm))
            end
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
backend and input paths, the point directory as run root, the derived
merger TOML, and the sweep's thread count."""
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
    # The point's config lives in its own directory: paths the base config
    # states relative to its location are made absolute here.
    haskey(sim, "input_file") &&
        (sim["input_file"] = _resolve_path(pipeline_dir, sim["input_file"]))
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
        arm_items = _arm_items(cfg, p.arm)
        for (k, v) in arm_items
            startswith(k, "merger.") && _set_nested_create!(m, k, v)
        end
        p.kind == "control" && (m = control_merger_dict(m))
        merger_path = joinpath(pdir, "merger.toml")
        _atomic_write_toml(merger_path, m)
        load_merger_config(merger_path)

        c = deepcopy(pipeline_base)
        _sweep_pipeline_overrides!(c, cfg, pipeline_dir, pdir, merger_path)
        for (k, v) in arm_items
            startswith(k, "pipeline.") || continue
            key = k[(length("pipeline.") + 1):end]
            # Paths of an arm resolve like those of the base pipeline config.
            key in ("install.install_dir", "simulation.input_file") &&
                v isa AbstractString &&
                (v = _resolve_path(pipeline_dir, v))
            _set_nested_create!(c, key, v)
        end
        config_path = joinpath(pdir, "config.toml")
        _atomic_write_toml(config_path, c)
        load_config(config_path)
    end

    cp(cfg.pipeline_config, joinpath(sweep_dir, "base_config.toml"))
    cp(cfg.merger_config, joinpath(sweep_dir, "base_merger.toml"))
    write_sweep_index(sweep_dir, cfg, points)
    return sweep_dir, points
end

"""
    write_sweep_index(sweep_dir, cfg, points; status = Dict(), entry = 1, created = "")

Write `sweep_index.toml`: the sweep header (name, base configs, axes,
seeds, job settings, the number of driver entries `entries` and the time of
the first, `created`) and one `[[points]]` entry per point with its id,
index, seed, axis values, directory and status (`pending`, `running`,
`done`, `failed`, or the run's own status such as `stopped` or `halted`;
finished points carry `exit_status` and `elapsed_seconds`, points pinned to
a CUDA device carry `device`, points launched by a driver entry carry
`entry`). `status` maps point indices to the fields that override the
pending defaults; `created` empty means now.
"""
function write_sweep_index(
    sweep_dir::AbstractString,
    cfg::SweepConfig,
    points::Vector{SweepPoint};
    status::Dict{Int,Dict{String,Any}} = Dict{Int,Dict{String,Any}}(),
    entry::Integer = 1,
    created::AbstractString = "",
)
    header = Dict{String,Any}(
        "name" => cfg.name,
        "created" =>
            isempty(created) ? Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS") : String(created),
        "entries" => Int(entry),
        "pipeline_config" => cfg.pipeline_config,
        "merger_config" => cfg.merger_config,
        "axes" => first.(cfg.grid),
        "seeds" => cfg.seeds,
        "concurrency" => cfg.concurrency,
        "omp_threads" => cfg.omp_threads,
        "controls" => cfg.controls,
        "devices" => cfg.devices,
        "engines_per_device" => cfg.engines_per_device,
        "arms" => first.(cfg.arms),
        "max_segments" => cfg.max_segments,
        "n_points" => length(points),
    )
    entries = map(points) do p
        e = Dict{String,Any}(
            "index" => p.index,
            "id" => p.id,
            "seed" => p.seed,
            "dir" => joinpath(sweep_dir, p.id),
            "values" => Dict{String,Any}(k => v for (k, v) in p.values),
            "kind" => p.kind,
            "control_of" => p.control_of,
            "arm" => p.arm,
            "status" => "pending",
        )
        merge!(e, get(status, p.index, Dict{String,Any}()))
        e
    end
    d = Dict{String,Any}("sweep" => header, "points" => entries)
    _atomic_write_toml(joinpath(sweep_dir, _SWEEP_INDEX_FILE), d)
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

"""
    continue_sweep_point(point_dir::AbstractString)

Worker entry of every sweep point: continue the point's run with
[`resume_pipeline`](@ref) when its run directory holds a segment record,
otherwise start it with [`run_sweep_point`](@ref). A complete,
post-processed point is left as it is; a stopped, killed or failed run gets
one more segment from its last complete dump; a completed run without the
pipeline's completion marker is post-processed. A record without any segment
(a pipeline interrupted before the engine started) is kept as a `#k` backup
and the point is started afresh.
"""
function continue_sweep_point(point_dir::AbstractString)
    run_dir = joinpath(point_dir, _SWEEP_RUN_ID)
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    if isfile(info_path) && _segment_count(info_path) > 0
        resume_pipeline(run_dir; base_dir = point_dir)
    else
        isfile(info_path) && _backup_existing(info_path)
        run_sweep_point(point_dir)
    end
    return nothing
end

"""
    _sweep_resume_state(run_dir) -> Tuple{String,String}

State of a sweep point's run at a re-entry, read from the run's own record
(the index of a driver that was cut is not trusted): `("done", "completed")`
for a completed run whose pipeline finished; `("final", status)` for a
`halted`, `watchdog` or `stalled` run, which a resume refuses (an explicit
restart continues it); `("continue", status)` otherwise: no record,
`stopped`, `killed`, `failed`, a record left `running`, or a completed run
that was not post-processed.
"""
function _sweep_resume_state(run_dir::AbstractString)::Tuple{String,String}
    status = _run_status(run_dir)
    status == "completed" && _pipeline_completed(run_dir) && return ("done", status)
    status in ("halted", "watchdog", "stalled") && return ("final", status)
    return ("continue", status)
end

"""Command of a worker process for `point_dir`: the current Julia in the
driver's own environment, with the driver's thread count (`Base.julia_cmd`
does not carry `--threads`), calling [`continue_sweep_point`](@ref), which
starts a new point and continues an interrupted one.

The worker reproduces the driver's session rather than the package
environment, because each point runs a full pipeline and therefore needs the
same figure backend the driver has: without it `run_pipeline` refuses a point
whose config asks for figures. The backend is loaded only when the driver
itself has one, so a sweep of a headless numerics-only study stays free of
the plotting stack."""
function _sweep_worker_command(point_dir::AbstractString)
    julia = Base.julia_cmd()
    project = dirname(something(Base.active_project()))
    preamble = plotting_available() ? "using CairoMakie; " : ""
    return `$julia --threads=$(Threads.nthreads()) --project=$project --startup-file=no -e "$(preamble)using Nbody6Dynamics; Nbody6Dynamics.continue_sweep_point(ARGS[1])" $point_dir`
end

"""
    _arm_pipeline_sections(cfg, name) -> (; install, build, simulation)

The install, build and simulation sections of the base pipeline configuration of `cfg`
with the `pipeline.…` overrides of the arm `name` applied (`""`: the base alone);
`install.install_dir` and `simulation.input_file` of the arm resolve against the base
file's directory, as in [`prepare_sweep`](@ref). The sections the driver's own checks
read, without writing a point configuration.
"""
function _arm_pipeline_sections(cfg::SweepConfig, name::AbstractString)
    c = TOML.parsefile(cfg.pipeline_config)
    pipeline_dir = dirname(cfg.pipeline_config)
    for (k, v) in _arm_items(cfg, name)
        startswith(k, "pipeline.") || continue
        key = k[(length("pipeline.") + 1):end]
        key in ("install.install_dir", "simulation.input_file") &&
            v isa AbstractString &&
            (v = _resolve_path(pipeline_dir, v))
        _set_nested_create!(c, key, v)
    end
    return (
        install = _parse_install(get(c, "install", Dict{String,Any}())),
        build = _parse_build(get(c, "build", Dict{String,Any}())),
        simulation = _parse_simulation(get(c, "simulation", Dict{String,Any}())),
    )
end

"""The configurations the driver checks: the base (`""`) of a sweep without arms, else every arm."""
_sweep_variants(cfg::SweepConfig) = isempty(cfg.arms) ? [""] : [first(a) for a in cfg.arms]

"""Fail before launching anything when the backend binary of the base
pipeline config, or of any arm (an arm may change the build or the binary
name), cannot be found. Returns the binary of the first variant."""
function _check_sweep_binary(cfg::SweepConfig)
    found = ""
    for name in _sweep_variants(cfg)
        s = _arm_pipeline_sections(cfg, name)
        install_dir = s.install.install_dir
        src_dir =
            isabspath(install_dir) ? install_dir :
            joinpath(dirname(cfg.pipeline_config), install_dir)
        path = _find_binary(src_dir, s.simulation.binary_name, s.build)
        isempty(found) && (found = path)
    end
    return found
end

"""
    _assign_device(load, devices, k) -> Union{Nothing,Int}

Device for the next engine: the entry of `devices` with the fewest engines
in `load` (device => running engines, a missing key counting as 0) among
those with fewer than `k`; ties go to the lowest device index. Returns
`nothing` when every device is full or `devices` is empty.

```julia
_assign_device(Dict(0 => 1), [0, 1], 1)   # 1
```
"""
function _assign_device(
    load::AbstractDict{Int,Int},
    devices::AbstractVector{Int},
    k::Integer,
)::Union{Nothing,Int}
    best = -1
    best_load = typemax(Int)
    for d in devices
        n = get(load, d, 0)
        n < k || continue
        if n < best_load || (n == best_load && d < best)
            best, best_load = d, n
        end
    end
    return best < 0 ? nothing : best
end

"""
A sweep worker process of this driver and the time its engine may need to
stop with a dump (the base configuration's `simulation.stop_margin`).
"""
struct _ActiveWorker
    pid::Int32
    stop_margin::Float64
end

"""
Sweep workers started by this process and not yet reaped; replaced, never
mutated, as [`_ACTIVE_ENGINES`](@ref).
"""
const _ACTIVE_WORKERS = Ref{Vector{_ActiveWorker}}(_ActiveWorker[])
"""Set when the driver loop dies from an exception, so that the exit hook ends the process with status 1 rather than 143 (a termination by SIGTERM)."""
const _DRIVER_FAILED = Ref(false)

"""Status the exit hook ends the driver with: 1 after an exception of the driver loop, else 143 (SIGTERM)."""
_driver_exit_status()::Cint = _DRIVER_FAILED[] ? Cint(1) : Cint(143)

_register_worker(w::_ActiveWorker) = (_ACTIVE_WORKERS[] = vcat(_ACTIVE_WORKERS[], [w]); nothing)
_deregister_worker(pid::Integer) =
    (_ACTIVE_WORKERS[] = filter(w -> w.pid != pid, _ACTIVE_WORKERS[]); nothing)

"""
    _terminate_workers() -> Bool

Send SIGTERM to every worker of [`_ACTIVE_WORKERS`](@ref) still running and
wait for them up to the largest `stop_margin` plus twice
`_KILL_GRACE_SECONDS`. A worker is a pipeline process, whose own exit hook
([`_stop_engines_at_exit`](@ref)) asks its engine to stop with a dump and
closes the segment record. Workers still running at the deadline are left to
the scheduler, since killing one would orphan its engine. Synchronous
operations only, as in [`_stop_engines_at_exit`](@ref). Returns whether any
worker was registered.
"""
function _terminate_workers()::Bool
    workers = _ACTIVE_WORKERS[]
    isempty(workers) && return false
    alive = Dict{Int32,Bool}()
    for w in workers
        running, _ = _child_state(w.pid)
        alive[w.pid] = running
        running && ccall(:kill, Cint, (Cint, Cint), w.pid, 15)
    end
    deadline = time() + maximum(w.stop_margin for w in workers) + 2 * _KILL_GRACE_SECONDS
    while any(values(alive)) && time() ≤ deadline
        Libc.systemsleep(0.2)
        for w in workers
            alive[w.pid] || continue
            running, _ = _child_state(w.pid)
            running || (alive[w.pid] = false)
        end
    end
    _ACTIVE_WORKERS[] = _ActiveWorker[]
    return true
end

"""
    _stop_workers_at_exit()

Process-exit hook of a sweep driver: hand a termination on to the workers
([`_terminate_workers`](@ref)). A batch scheduler sends its early warning
signal (`--signal=B:TERM@…`) to the batch shell only, which in a sweep job is
the driver; without the hand-over the workers would see a signal only at the
time limit itself, too late for a dump.

When it has handed a termination on, the hook ends the process itself with
`_exit(143)`, the status of a termination by SIGTERM (1 when the driver loop
died from an exception, [`_driver_exit_status`](@ref)), after stopping any
engine of the driver's own ([`_stop_engines_at_exit`](@ref)). Julia runs the
exit of a SIGTERM on the main thread at its next safepoint, which in the
driver's poll loop can lie inside a libuv timer callback that holds the I/O
lock; the handle finalizers that follow the hooks then wait on that lock
forever (observed 2026-11-04: workers stopped cleanly, the driver spinning at
50 % CPU until killed).
"""
function _stop_workers_at_exit()
    _terminate_workers() || return nothing
    _stop_engines_at_exit()
    ccall(:_exit, Cvoid, (Cint,), _driver_exit_status())
    return nothing
end

"""Stop margin of the sweep, the time a worker's engine may take to stop with a dump: the largest over the base pipeline configuration of `cfg` and its arms (an arm may raise it); the default of [`SimulationConfig`](@ref) when the configuration cannot be read."""
function _sweep_stop_margin(cfg::SweepConfig)::Float64
    return try
        maximum(
            _arm_pipeline_sections(cfg, name).simulation.stop_margin for
            name in _sweep_variants(cfg)
        )
    catch e
        e isa Union{ArgumentError,SystemError,TOML.ParserError} || rethrow()
        SimulationConfig().stop_margin
    end
end

"""Whether the point in `point_dir` runs the CUDA build (`build.enable_gpu` of its config), and so takes a device slot."""
function _point_uses_gpu(point_dir::AbstractString)::Bool
    c = TOML.parsefile(joinpath(point_dir, "config.toml"))
    return get(get(c, "build", Dict{String,Any}()), "enable_gpu", false) === true
end

"""
    _continue_in_entry(run_dir, max_segments) -> Bool

Whether a point stopped by its wall budget runs its next segment within the same
driver entry: fewer than `max_segments` segments on record and progress in the last
one (`t_end` beyond `t_start`).
"""
function _continue_in_entry(run_dir::AbstractString, max_segments::Integer)::Bool
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    isfile(info_path) || return false
    segments = get(TOML.parsefile(info_path), "segments", Any[])
    (isempty(segments) || length(segments) ≥ max_segments) && return false
    last_segment = segments[end]
    t0 = get(last_segment, "t_start", 0.0)
    t1 = get(last_segment, "t_end", nothing)
    return t1 isa Real && t0 isa Real && t1 > t0
end

"""
    _launch_sweep(sweep_dir, cfg, points; queued = points, status = Dict(), entry = 1,
                  created = "", worker_command = _sweep_worker_command)
        -> Dict{Int,Dict{String,Any}}

Run the points of `queued` (all of `points` by default) as worker processes,
at most `cfg.concurrency` at a time, logging each to
`<point>/sweep_point.log` (appended) and rewriting the index of all `points`
on every state change, starting from the per-point fields of `status`.
Returns the per-point status table. `entry` and `created` go to the index
header ([`write_sweep_index`](@ref)); every launched point records `entry`.
Each worker is known to the driver's exit hook until it is reaped
([`_stop_workers_at_exit`](@ref)).

With `cfg.devices` nonempty a queued point is admitted only when
[`_assign_device`](@ref) finds a device holding fewer than
`cfg.engines_per_device` engines; its `config.toml` then receives
`simulation.gpu_list = [device]` before the launch, so the frozen config
records the device that ran it, and its status carries `device`.
`worker_command` maps a point directory to the worker `Cmd`.
"""
function _launch_sweep(args...; kwargs...)
    try
        return _launch_sweep_loop(args...; kwargs...)
    catch
        _DRIVER_FAILED[] = true
        rethrow()
    end
end

function _launch_sweep_loop(
    sweep_dir::AbstractString,
    cfg::SweepConfig,
    points::Vector{SweepPoint};
    queued::Vector{SweepPoint} = points,
    status::Dict{Int,Dict{String,Any}} = Dict{Int,Dict{String,Any}}(),
    entry::Integer = 1,
    created::AbstractString = "",
    worker_command::Function = _sweep_worker_command,
)
    status = deepcopy(status)
    created = isempty(created) ? Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS") : String(created)
    write_index() =
        write_sweep_index(sweep_dir, cfg, points; status = status, entry = entry, created = created)
    stop_margin = _sweep_stop_margin(cfg)
    # The entry is on record even when it has nothing to run.
    write_index()
    queue = copy(queued)
    running = Dict{Int,Tuple{Base.Process,IOStream,Float64,SweepPoint,Int,Int32}}()
    load = Dict{Int,Int}()
    n_done = 0
    n_failed = 0
    while !isempty(queue) || !isempty(running)
        while length(running) < cfg.concurrency && !isempty(queue)
            device = -1
            p = first(queue)
            pdir = joinpath(sweep_dir, p.id)
            if !isempty(cfg.devices) && _point_uses_gpu(pdir)
                d = _assign_device(load, cfg.devices, cfg.engines_per_device)
                d === nothing && break
                device = d
            end
            popfirst!(queue)
            if device ≥ 0
                config_path = joinpath(pdir, "config.toml")
                c = TOML.parsefile(config_path)
                _set_nested!(c, "simulation.gpu_list", [device])
                _atomic_write_toml(config_path, c)
            end
            # Appended: a re-entry keeps the log of the earlier entries.
            log = open(joinpath(pdir, "sweep_point.log"), "a")
            proc = run(pipeline(worker_command(pdir); stdout = log, stderr = log); wait = false)
            # Taken once: `getpid` throws for a process that has already exited.
            wpid = try
                Int32(getpid(proc))
            catch e
                e isa Base.IOError || rethrow()
                Int32(0)
            end
            wpid > 0 && _register_worker(_ActiveWorker(wpid, stop_margin))
            running[p.index] = (proc, log, time(), p, device, wpid)
            status[p.index] = Dict{String,Any}(
                "status" => "running",
                "started" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
                "entry" => Int(entry),
            )
            if device ≥ 0
                status[p.index]["device"] = device
                load[device] = get(load, device, 0) + 1
            end
            write_index()
            on_device = device ≥ 0 ? " on device $device" : ""
            @info "sweep: started $(p.id)$(on_device) " *
                  "($(length(running)) running, $(length(queue)) queued)"
        end
        sleep(cfg.poll_interval)
        for (i, (proc, log, t0, p, device, wpid)) in collect(running)
            process_running(proc) && continue
            close(log)
            _deregister_worker(wpid)
            code = _exit_status(proc)
            point_status = _sweep_point_status(code, joinpath(sweep_dir, p.id, _SWEEP_RUN_ID))
            ok = point_status == "done"
            # A clean exit only: a worker that failed before it opened a segment leaves the
            # record as it was, and the point would be launched again without end.
            again =
                point_status == "stopped" &&
                code == 0 &&
                _continue_in_entry(joinpath(sweep_dir, p.id, _SWEEP_RUN_ID), cfg.max_segments)
            again && push!(queue, p)
            ok ? (n_done += 1) : again || (n_failed += 1)
            status[i] = Dict{String,Any}(
                "status" => point_status,
                "started" => status[i]["started"],
                "entry" => Int(entry),
                "exit_status" => code,
                "elapsed_seconds" => round(time() - t0; digits = 1),
            )
            if device ≥ 0
                status[i]["device"] = device
                load[device] -= 1
            end
            delete!(running, i)
            write_index()
            if ok
                @info "sweep: finished $(p.id) in $(_format_elapsed(time() - t0))"
            elseif again
                @info "sweep: $(p.id) stopped at its wall budget; queued for its next segment"
            else
                @warn "sweep: $(p.id) ended as $(point_status) (exit status $code); see $(joinpath(sweep_dir, p.id, "sweep_point.log"))"
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

"""
    resume_sweep(cfg::SweepConfig, sweep_dir; worker_command = _sweep_worker_command) -> String

Re-enter the sweep in `sweep_dir`, prepared earlier from the same sweep
configuration `cfg`, after its driver was ended (a job cut at its time
limit, a cancelled job, a crash): every point is classified from its own run
record, not from the index. Complete,
post-processed points are skipped; `halted`, `watchdog` and `stalled` points
keep their status (a restart continues them explicitly); all others are
launched again ([`continue_sweep_point`](@ref): stopped and killed runs
resume from their last dump, unstarted points start). The launch follows the
`concurrency`, `devices` and `engines_per_device` of `cfg`, so a re-entry may
run with other resources than the first entry; the index records the entry
count. `sweep_summary.csv` is rewritten at the end.

The configuration must describe the indexed points: same number, ids, seeds
and kinds, else an `ArgumentError` names the first difference (a sweep is
not extended or changed by a re-entry). Returns `sweep_dir`.

# Example
```julia
cfg = load_sweep_config("sweep.toml")
resume_sweep(cfg, "runs/sweep_pilot_20261104_120000")
```
"""
function resume_sweep(
    cfg::SweepConfig,
    sweep_dir::AbstractString;
    worker_command::Function = _sweep_worker_command,
)::String
    sweep_dir = abspath(sweep_dir)
    idx = read_sweep_index(sweep_dir)
    points = sweep_points(cfg)
    indexed = idx["points"]
    length(indexed) == length(points) || throw(
        ArgumentError(
            "resume_sweep: $(sweep_dir) indexes $(length(indexed)) points, the configuration " *
            "describes $(length(points)); a re-entry runs the sweep that was prepared",
        ),
    )
    for (e, p) in zip(indexed, points)
        (
            e["id"] == p.id &&
            Int(e["seed"]) == p.seed &&
            get(e, "kind", "merger") == p.kind &&
            get(e, "arm", "") == p.arm
        ) || throw(
            ArgumentError(
                "resume_sweep: point $(p.index) of the configuration ($(p.id), seed " *
                "$(p.seed), $(p.kind)) differs from the index ($(e["id"]), seed " *
                "$(e["seed"]), $(get(e, "kind", "merger")))",
            ),
        )
    end
    _check_sweep_binary(cfg)
    previous = Dict{Int,Dict{String,Any}}(Int(e["index"]) => e for e in indexed)
    status = Dict{Int,Dict{String,Any}}()
    queued = SweepPoint[]
    counts = Dict("done" => 0, "final" => 0, "continue" => 0)
    for p in points
        state, run_status = _sweep_resume_state(joinpath(sweep_dir, p.id, _SWEEP_RUN_ID))
        counts[state] += 1
        if state == "continue"
            push!(queued, p)
            continue
        end
        e = Dict{String,Any}("status" => state == "done" ? "done" : run_status)
        for k in ("started", "entry", "exit_status", "elapsed_seconds", "device")
            haskey(previous[p.index], k) && (e[k] = previous[p.index][k])
        end
        status[p.index] = e
    end
    entry = Int(get(idx["sweep"], "entries", 1)) + 1
    created = String(get(idx["sweep"], "created", ""))
    @info "sweep \"$(cfg.name)\": entry $(entry) of $(sweep_dir): $(counts["done"]) done, " *
          "$(counts["final"]) ended for good (halted, watchdog, stalled), " *
          "$(counts["continue"]) to run"
    _launch_sweep(
        sweep_dir,
        cfg,
        points;
        queued = queued,
        status = status,
        entry = entry,
        created = created,
        worker_command = worker_command,
    )
    write_sweep_summary(sweep_dir)
    return sweep_dir
end

"""
    _sweep_point_status(exit_status, run_dir) -> String

Index status of a finished sweep point from the worker's `exit_status` and
the run status of its run directory ([`_run_status`](@ref)): `"done"` when the
worker exited with 0 and the run is `completed`; otherwise the run status
(`"halted"`, `"stopped"`, `"failed"`, …), and `"failed"` for a nonzero exit of
a run recorded as completed (the post-processing failed). A run without a
status (no record, or a record that predates segment statuses) is judged by
the exit status alone.
"""
function _sweep_point_status(exit_status::Integer, run_dir::AbstractString)::String
    run_status = _run_status(run_dir)
    if isempty(run_status)
        return exit_status == 0 ? "done" : "failed"
    elseif run_status == "completed"
        return exit_status == 0 ? "done" : "failed"
    end
    return run_status
end

"""Outcome fields of one finished point read from its run directory:
elapsed time, exit status and completion (`END RUN` seen, whether or not
the engine exited by itself) from `RUN_INFO.toml`, final time, star and
pair counts, energy error and virial ratio from the last ADJUST record.
Missing files leave `NaN`/`-1`/`false` entries."""
function _sweep_point_outcome(run_dir::AbstractString)
    out = Dict{String,Any}(
        "elapsed_seconds" => NaN,
        "exit_status" => -1,
        "completed" => false,
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
        if !isempty(segs)
            out["exit_status"] = Int(get(segs[end], "exit_status", -1))
            out["completed"] = Bool(get(segs[end], "completed", out["exit_status"] == 0))
        end
    end
    stdout_path = _run_stdout_path(run_dir)
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

One row per point: index, id, kind, control id, arm (only when the sweep has arms), seed,
the axis values, status, the
assigned CUDA device (`-1` when none), and the outcome fields of
[`_sweep_point_outcome`](@ref) for finished points.
`columns` gives the column order used by [`write_sweep_summary`](@ref).
"""
function sweep_summary(sweep_dir::AbstractString)
    idx = read_sweep_index(sweep_dir)
    axes = String[idx["sweep"]["axes"]...]
    columns = vcat(
        ["index", "id", "kind", "control_of"],
        isempty(get(idx["sweep"], "arms", String[])) ? String[] : ["arm"],
        ["seed"],
        axes,
        [
            "status",
            "device",
            "exit_status",
            "completed",
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
            "kind" => get(p, "kind", "merger"),
            "control_of" => get(p, "control_of", ""),
            "arm" => get(p, "arm", ""),
            "seed" => p["seed"],
            "status" => p["status"],
            "device" => get(p, "device", -1),
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
    _backup_existing(path)
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

"""Completed points of a sweep of the given `kind` (`"merger"` by default,
`"control"`, or `""` for all) as `(entry, run_dir)` pairs: index status
`done` and a run that is `completed` ([`_run_status`](@ref); a run without a
recorded status is kept). Points marked done whose run ended otherwise (an
index written before the run status was consulted) are left out with a
warning that counts them."""
function _sweep_done_points(sweep_dir::AbstractString; kind::AbstractString = "merger")
    idx = read_sweep_index(sweep_dir)
    pts = Tuple{Dict{String,Any},String}[]
    n_excluded = 0
    for p in idx["points"]
        (p["status"] == "done" && (isempty(kind) || get(p, "kind", "merger") == kind)) || continue
        run_dir = joinpath(p["dir"], _SWEEP_RUN_ID)
        if _run_status(run_dir) in ("", "completed")
            push!(pts, (p, run_dir))
        else
            n_excluded += 1
        end
    end
    n_excluded > 0 &&
        @warn "sweep: $(n_excluded) point(s) marked done in $(sweep_dir) did not complete their " *
              "run and are excluded from the ensembles"
    return idx, pts
end
