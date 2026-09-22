# =============================================================================
# Configuration loading and serialization
# =============================================================================

# ---------------------------------------------------------------------------
# Shared parser helpers — one error form, one unknown-key check, one typed
# getter, used by every configuration parser of the package.
# ---------------------------------------------------------------------------

"""
    _config_error(key, msg)

Throw the `ArgumentError` every configuration parser raises: `"config: <key> <msg>"`,
where `key` is the dotted `section.key` of the offending entry.
"""
_config_error(key::AbstractString, msg::AbstractString) = throw(ArgumentError("config: $key $msg"))

"""
    _reject_unknown(table, allowed, section)

Raise `_config_error` for the first key of `table` (a parsed TOML table) that is
not in `allowed`, naming it as `<section>.<key>` and listing the allowed keys
(sorted, comma-separated). `section` is the dotted table name (`"simulation"`,
`"merger.cluster1.imf"`). Returns `nothing`.
"""
function _reject_unknown(table::AbstractDict, allowed, section::AbstractString)
    for k in sort!(collect(String, keys(table)))
        k in allowed || _config_error(
            "$section.$k",
            "is not a recognised key; allowed keys of [$section]: " *
            join(sort!(collect(String, allowed)), ", "),
        )
    end
    return nothing
end

_typed_error(key::AbstractString, what::AbstractString, value) =
    _config_error(key, "must be a $what; got $value ($(typeof(value)))")

"""
    _typed(value, T, key) -> T

Coerce a TOML `value` to the field type `T` of a configuration struct or raise
`_config_error(key, ...)`. Accepted: a `Bool` for `Bool`; a non-`Bool` `Integer`
for an `Integer` type; a non-`Bool` `Real` for an `AbstractFloat` type
(converted with `Float64`); an `AbstractString` for `String`; an
`AbstractVector` whose every element passes the element rule for
`Vector{String}`, `Vector{Int}` and `Vector{Float64}`. Anything else raises
`"must be a <description>; got <value> (<typeof(value)>)"` with the descriptions
`"boolean"`, `"integer"`, `"number"`, `"string"`, `"array of strings"`,
`"array of integers"`, `"array of numbers"`.
"""
function _typed end

_typed(value, ::Type{Bool}, key::AbstractString) =
    value isa Bool ? value : _typed_error(key, "boolean", value)

_typed(value, ::Type{T}, key::AbstractString) where {T<:Integer} =
    (value isa Integer && !(value isa Bool)) ? T(value) : _typed_error(key, "integer", value)

_typed(value, ::Type{T}, key::AbstractString) where {T<:AbstractFloat} =
    (value isa Real && !(value isa Bool)) ? Float64(value) : _typed_error(key, "number", value)

_typed(value, ::Type{String}, key::AbstractString) =
    value isa AbstractString ? String(value) : _typed_error(key, "string", value)

function _typed(value, ::Type{Vector{String}}, key::AbstractString)
    (value isa AbstractVector && all(v -> v isa AbstractString, value)) ||
        _typed_error(key, "array of strings", value)
    return String[String(v) for v in value]
end

function _typed(value, ::Type{Vector{Int}}, key::AbstractString)
    (value isa AbstractVector && all(v -> v isa Integer && !(v isa Bool), value)) ||
        _typed_error(key, "array of integers", value)
    return Int[Int(v) for v in value]
end

function _typed(value, ::Type{Vector{Float64}}, key::AbstractString)
    (value isa AbstractVector && all(v -> v isa Real && !(v isa Bool), value)) ||
        _typed_error(key, "array of numbers", value)
    return Float64[Float64(v) for v in value]
end

"""
    _parse_section(T, table, section; skip = ()) -> T

Build the `@kwdef` struct `T` from the TOML `table` of `[section]`: unknown keys
are rejected against `fieldnames(T)` minus `skip`, every present key is coerced
with `_typed` to its field type, and absent keys take the struct's own default —
the parser states no default of its own. Keys in `skip` are neither rejected
nor read (the caller handles them).
"""
function _parse_section(
    ::Type{T},
    table::AbstractDict,
    section::AbstractString;
    skip = (),
) where {T}
    fields = setdiff(String.(fieldnames(T)), String.(collect(skip)))
    _reject_unknown(table, vcat(fields, String.(collect(skip))), section)
    kwargs = Dict{Symbol,Any}()
    for k in fields
        haskey(table, k) || continue
        kwargs[Symbol(k)] = _typed(table[k], fieldtype(T, Symbol(k)), "$section.$k")
    end
    return T(; kwargs...)
end

const _CONFIG_SECTIONS =
    ("install", "build", "simulation", "postprocess", "visualization", "merger")

"""
    load_config(path::AbstractString) -> Nbody6Config

Parse a TOML configuration file and return a fully-typed `Nbody6Config`.
Missing sections or keys fall back to defaults defined in the `@kwdef` structs;
keys the configuration schema does not define are rejected, so a typo cannot
pass silently. The parsed configuration is validated fail-fast by
[`_validate`](@ref): enumerated choices, numerical bounds, and required
non-empty fields raise an `ArgumentError` naming the offending `section.key`
before any pipeline phase starts.
The directory of `path` is recorded as `config_dir`: relative paths in the
configuration resolve against it, and it is the default `base_dir` of the
entry points, so a project directory holding a `config.toml` receives its
`backend/` and `runs/` there, wherever the package is installed.
"""
function load_config(path::AbstractString)::Nbody6Config
    raw = TOML.parsefile(path)
    _reject_unknown(raw, _CONFIG_SECTIONS, basename(path))
    for s in _CONFIG_SECTIONS
        haskey(raw, s) && (raw[s] isa AbstractDict || _config_error(s, "must be a table"))
    end
    cfg = Nbody6Config(
        _parse_install(get(raw, "install", Dict{String,Any}())),
        _parse_build(get(raw, "build", Dict{String,Any}())),
        _parse_simulation(get(raw, "simulation", Dict{String,Any}())),
        _parse_postprocess(get(raw, "postprocess", Dict{String,Any}())),
        _parse_visualization(get(raw, "visualization", Dict{String,Any}())),
        _parse_merger_pipeline(get(raw, "merger", Dict{String,Any}())),
        dirname(abspath(path)),
    )
    _validate(cfg)
    return cfg
end

"""
    save_config(cfg::Nbody6Config, path::AbstractString)

Serialize a config struct back to a TOML file (for frozen run snapshots).
`config_dir` is not written: it is a property of where the file lives. An
existing file is backed up, never overwritten, and the write is atomic.
"""
function save_config(cfg::Nbody6Config, path::AbstractString)
    d = Dict(
        "install" => _struct_to_dict(cfg.install),
        "build" => _struct_to_dict(cfg.build),
        "simulation" => _struct_to_dict(cfg.simulation),
        "postprocess" => _struct_to_dict(cfg.postprocess),
        "visualization" => _struct_to_dict(cfg.visualization),
        "merger" => _struct_to_dict(cfg.merger),
    )
    # TOML cannot serialize nested structs — convert explicitly
    d["visualization"]["style"] = _struct_to_dict(cfg.visualization.style)
    _backup_existing(path)
    _atomic_write_toml(path, d)
    return nothing
end

# ---------------------------------------------------------------------------
# Section parsers — convert Dict{String,Any} → typed struct
# ---------------------------------------------------------------------------

_parse_install(d::Dict) = _parse_section(InstallConfig, d, "install")

_parse_build(d::Dict) = _parse_section(BuildConfig, d, "build")

_parse_simulation(d::Dict) = _parse_section(SimulationConfig, d, "simulation")

_parse_postprocess(d::Dict) = _parse_section(PostprocessConfig, d, "postprocess")

_parse_merger_pipeline(d::Dict) = _parse_section(MergerPipelineConfig, d, "merger")

function _parse_visualization(d::Dict)
    st = get(d, "style", Dict{String,Any}())
    st isa AbstractDict || _config_error("visualization.style", "must be a table")
    style = _parse_section(PlotStyle, st, "visualization.style")
    vis = _parse_section(VisualizationConfig, d, "visualization"; skip = ("style",))
    return VisualizationConfig(;
        (k => getfield(vis, k) for k in fieldnames(VisualizationConfig) if k != :style)...,
        style = style,
    )
end

# ---------------------------------------------------------------------------
# Fail-fast validation — the parser enforces exactly the constraints the
# config comments document (types, enumerated choices, numerical bounds).
# ---------------------------------------------------------------------------

"""
    _validate(cfg::Nbody6Config)

Validate a parsed configuration against the constraints documented in
`config.toml` and the manual: enumerated choices, numerical bounds, and
required non-empty fields. Raises an `ArgumentError` naming the offending
`section.key` and the actual value; returns `nothing` on success. Called
unconditionally at the end of [`load_config`](@ref) so a pipeline cannot
start from a configuration it cannot honor.

File-existence checks (input files, merger config paths) remain at runtime:
those paths may be created by earlier pipeline phases.
"""
function _validate(cfg::Nbody6Config)
    inst = cfg.install
    bld = cfg.build
    sim = cfg.simulation
    pp = cfg.postprocess
    vis = cfg.visualization
    st = vis.style

    # [install] / [build]
    isempty(inst.install_dir) && _config_error("install.install_dir", "must be nonempty")
    bld.nproc ≥ 0 || _config_error("build.nproc", "must be ≥ 0; got $(bld.nproc)")
    for arch in bld.cuda_arch
        occursin(_CUDA_ARCH_PATTERN, arch) || _config_error(
            "build.cuda_arch",
            "entries must be CUDA architecture names of the form " *
            "\"sm_<major><minor>\" (e.g. \"sm_90\", \"sm_120\"); got \"$arch\"",
        )
    end
    for flag in bld.nvcc_flags
        isempty(strip(flag)) && _config_error("build.nvcc_flags", "entries must be nonempty")
    end

    # [simulation]
    sim.mpi_ranks ≥ 1 || _config_error("simulation.mpi_ranks", "must be ≥ 1; got $(sim.mpi_ranks)")
    if sim.mpi_ranks > 1 && !bld.enable_mpi
        _config_error(
            "simulation.mpi_ranks",
            "= $(sim.mpi_ranks) requires build.enable_mpi = true — " *
            "a multi-rank launch needs an MPI-enabled binary",
        )
    end
    sim.omp_threads ≥ 0 ||
        _config_error("simulation.omp_threads", "must be ≥ 0; got $(sim.omp_threads)")
    if !isempty(sim.gpu_list)
        bld.enable_gpu || _config_error(
            "simulation.gpu_list",
            "= $(sim.gpu_list) requires build.enable_gpu = true — " *
            "the CPU binary ignores GPU_LIST",
        )
        all(≥(0), sim.gpu_list) || _config_error(
            "simulation.gpu_list",
            "entries must be ≥ 0 (CUDA device indices); got $(sim.gpu_list)",
        )
        allunique(sim.gpu_list) ||
            _config_error("simulation.gpu_list", "entries must be distinct; got $(sim.gpu_list)")
        length(sim.gpu_list) ≤ _MAX_GPU_PER_PROCESS || _config_error(
            "simulation.gpu_list",
            "may name at most $(_MAX_GPU_PER_PROCESS) devices per " *
            "process (the engine's MAX_GPU); got $(length(sim.gpu_list))",
        )
    end
    sim.telemetry_interval ≥ 0 || _config_error(
        "simulation.telemetry_interval",
        "must be ≥ 0 [s]; got $(sim.telemetry_interval)",
    )
    sim.startup_timeout ≥ 0 ||
        _config_error("simulation.startup_timeout", "must be ≥ 0 [s]; got $(sim.startup_timeout)")
    sim.exit_grace ≥ 0 ||
        _config_error("simulation.exit_grace", "must be ≥ 0 [s]; got $(sim.exit_grace)")
    sim.live_interval ≥ 1 ||
        _config_error("simulation.live_interval", "must be ≥ 1 [s]; got $(sim.live_interval)")
    isempty(sim.runs_dir) && _config_error("simulation.runs_dir", "must be nonempty")
    isempty(sim.run_id_prefix) && _config_error("simulation.run_id_prefix", "must be nonempty")
    isempty(sim.input_file) && _config_error("simulation.input_file", "must be nonempty")

    # [postprocess]
    if pp.snapshot_format != "conf3"
        if pp.snapshot_format == "hdf5"
            _config_error(
                "postprocess.snapshot_format",
                "= \"hdf5\" is no longer supported; " *
                "use \"conf3\". The fork's KZ(46) H5Part layout was never readable by the old code.",
            )
        end
        _config_error(
            "postprocess.snapshot_format",
            "must be \"conf3\" (the only supported value); " * "got \"$(pp.snapshot_format)\"",
        )
    end
    isempty(pp.snapshot_pattern) &&
        _config_error("postprocess.snapshot_pattern", "must be nonempty")
    pp.parse_stdout &&
        isempty(pp.stdout_file) &&
        _config_error("postprocess.stdout_file", "must be nonempty when parse_stdout = true")
    pp.read_lagr &&
        isempty(pp.lagr_file) &&
        _config_error("postprocess.lagr_file", "must be nonempty when read_lagr = true")
    pp.read_escapers &&
        isempty(pp.escapers_file) &&
        _config_error("postprocess.escapers_file", "must be nonempty when read_escapers = true")
    pp.read_stellar_evo &&
        isempty(pp.stellar_evo_pattern) &&
        _config_error(
            "postprocess.stellar_evo_pattern",
            "must be nonempty when read_stellar_evo = true",
        )
    pp.read_binary_evo &&
        isempty(pp.binary_evo_pattern) &&
        _config_error(
            "postprocess.binary_evo_pattern",
            "must be nonempty when read_binary_evo = true",
        )

    # [visualization]
    vis.format in ("pdf", "svg", "png") || _config_error(
        "visualization.format",
        "must be one of \"pdf\", \"svg\", \"png\"; got \"$(vis.format)\"",
    )
    vis.column in ("single", "double", "") || _config_error(
        "visualization.column",
        "must be one of \"single\", \"double\", \"\"; got \"$(vis.column)\"",
    )
    vis.units in ("physical", "nbody") || _config_error(
        "visualization.units",
        "must be \"physical\" or \"nbody\"; got \"$(vis.units)\"",
    )
    vis.dpi ≥ 72 || _config_error("visualization.dpi", "must be ≥ 72; got $(vis.dpi)")
    vis.export_width > 0 ||
        _config_error("visualization.export_width", "must be > 0; got $(vis.export_width)")

    # [visualization.style]
    st.marker_budget > 0 ||
        _config_error("visualization.style.marker_budget", "must be > 0; got $(st.marker_budget)")
    (0 < st.marker_min ≤ st.marker_max) || _config_error(
        "visualization.style",
        "must satisfy 0 < marker_min ≤ marker_max; " *
        "got marker_min = $(st.marker_min), marker_max = $(st.marker_max)",
    )
    st.q_log_threshold > 0 || _config_error(
        "visualization.style.q_log_threshold",
        "must be > 0; got $(st.q_log_threshold)",
    )
    (0 < st.q_floor < 1) || _config_error(
        "visualization.style.q_floor",
        "must satisfy 0 < q_floor < 1; got $(st.q_floor)",
    )
    (0 < st.zoom_frac ≤ 1) || _config_error(
        "visualization.style.zoom_frac",
        "must satisfy 0 < zoom_frac ≤ 1; got $(st.zoom_frac)",
    )
    st.anim_fps ≥ 0 ||
        _config_error("visualization.style.anim_fps", "must be ≥ 0; got $(st.anim_fps)")
    st.anim_px_per_unit > 0 || _config_error(
        "visualization.style.anim_px_per_unit",
        "must be > 0; got $(st.anim_px_per_unit)",
    )
    st.anim_target_seconds > 0 || _config_error(
        "visualization.style.anim_target_seconds",
        "must be > 0; got $(st.anim_target_seconds)",
    )

    # [merger]
    cfg.merger.enabled &&
        isempty(cfg.merger.config_file) &&
        _config_error("merger.config_file", "must be nonempty when merger.enabled = true")

    return nothing
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _struct_to_dict(s)
    Dict(String(k) => getfield(s, k) for k in fieldnames(typeof(s)))
end

"""
    _resolve_path(base_dir, path) -> String

`path` made absolute: unchanged when already absolute, otherwise joined onto
`base_dir` and normalised. An empty `path` stays empty.
"""
function _resolve_path(base_dir::AbstractString, path::AbstractString)::String
    isempty(path) && return ""
    return isabspath(path) ? String(path) : normpath(joinpath(abspath(base_dir), path))
end

"""
    _with_absolute_paths(cfg, base_dir) -> Nbody6Config

A copy of `cfg` whose path-valued keys (`install.install_dir`,
`simulation.input_file`, `simulation.runs_dir`, `postprocess.data_dir`,
`merger.config_file`) are absolute, resolved against `base_dir`, and whose
`config_dir` is `base_dir`. The frozen `config.toml` of a run is written from
this form so that it records the paths the run used and can be loaded from
anywhere.
"""
function _with_absolute_paths(cfg::Nbody6Config, base_dir::AbstractString)::Nbody6Config
    rebuild(x; kw...) = typeof(x)(; (k => getfield(x, k) for k in fieldnames(typeof(x)))..., kw...)
    return Nbody6Config(
        rebuild(cfg.install; install_dir = _resolve_path(base_dir, cfg.install.install_dir)),
        cfg.build,
        rebuild(
            cfg.simulation;
            input_file = _resolve_path(base_dir, cfg.simulation.input_file),
            runs_dir = _resolve_path(base_dir, cfg.simulation.runs_dir),
        ),
        rebuild(cfg.postprocess; data_dir = _resolve_path(base_dir, cfg.postprocess.data_dir)),
        cfg.visualization,
        rebuild(cfg.merger; config_file = _resolve_path(base_dir, cfg.merger.config_file)),
        String(abspath(base_dir)),
    )
end
