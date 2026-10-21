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
    sec = String(section)::String
    names = sort!(collect(String, allowed))::Vector{String}
    for k in sort!(collect(String, keys(table)))
        k in names || _config_error(
            "$sec.$k",
            "is not a recognised key; allowed keys of [$sec]: " * join(names, ", "),
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

function _typed(value, ::Type{Vector{Vector{Float64}}}, key::AbstractString)
    (
        value isa AbstractVector &&
        all(row -> row isa AbstractVector && all(v -> v isa Real && !(v isa Bool), row), value)
    ) || _typed_error(key, "array of arrays of numbers", value)
    return Vector{Float64}[Float64[Float64(v) for v in row] for row in value]
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
    ("install", "build", "simulation", "postprocess", "visualization", "merger", "cluster", "cost")

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
        _parse_cluster(get(raw, "cluster", Dict{String,Any}())),
        _parse_cost(get(raw, "cost", Dict{String,Any}())),
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
        "cluster" => _struct_to_dict(cfg.cluster),
        "cost" => _struct_to_dict(cfg.cost),
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

_parse_cluster(d::Dict) = _parse_section(ClusterConfig, d, "cluster")

_parse_cost(d::Dict) = _parse_section(CostConfig, d, "cost")

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
    sim.stall_timeout ≥ 0 ||
        _config_error("simulation.stall_timeout", "must be ≥ 0 [s]; got $(sim.stall_timeout)")
    sim.exit_grace ≥ 0 ||
        _config_error("simulation.exit_grace", "must be ≥ 0 [s]; got $(sim.exit_grace)")
    sim.wall_budget ≥ 0 ||
        _config_error("simulation.wall_budget", "must be ≥ 0 [s]; got $(sim.wall_budget)")
    sim.stop_margin ≥ 0 ||
        _config_error("simulation.stop_margin", "must be ≥ 0 [s]; got $(sim.stop_margin)")
    if sim.wall_budget > 0
        sim.stop_margin < sim.wall_budget || _config_error(
            "simulation.stop_margin",
            "must be smaller than simulation.wall_budget = $(sim.wall_budget) s; " *
            "got $(sim.stop_margin)",
        )
    end
    (sim.checkpoint_keep == 0 || sim.checkpoint_keep ≥ 2) || _config_error(
        "simulation.checkpoint_keep",
        "must be 0 (keep every dump) or ≥ 2 (the newest dump may be incomplete); " *
        "got $(sim.checkpoint_keep)",
    )
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
    pp.pair_sum_max_n ≥ 0 || _config_error(
        "postprocess.pair_sum_max_n",
        "must be ≥ 0 (0 disables the pair-sum diagnostics); got $(pp.pair_sum_max_n)",
    )
    pp.bound_max_passes ≥ 1 ||
        _config_error("postprocess.bound_max_passes", "must be ≥ 1; got $(pp.bound_max_passes)")
    pp.bound_min_members ≥ 2 ||
        _config_error("postprocess.bound_min_members", "must be ≥ 2; got $(pp.bound_min_members)")
    pp.energy_jump_threshold ≥ 0 || _config_error(
        "postprocess.energy_jump_threshold",
        "must be ≥ 0; got $(pp.energy_jump_threshold)",
    )
    pp.energy_jump_sigma ≥ 0 ||
        _config_error("postprocess.energy_jump_sigma", "must be ≥ 0; got $(pp.energy_jump_sigma)")
    pp.mass_budget &&
        isempty(pp.events_file) &&
        _config_error("postprocess.events_file", "must be nonempty when mass_budget = true")
    _validate_estimators(pp)

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
    st.snapshot_render in ("auto", "scatter", "density") || _config_error(
        "visualization.style.snapshot_render",
        "must be one of \"auto\", \"scatter\", \"density\"; got \"$(st.snapshot_render)\"",
    )
    st.density_min_n ≥ 1 ||
        _config_error("visualization.style.density_min_n", "must be ≥ 1; got $(st.density_min_n)")
    (16 ≤ st.density_bins ≤ 4096) || _config_error(
        "visualization.style.density_bins",
        "must satisfy 16 ≤ density_bins ≤ 4096; got $(st.density_bins)",
    )
    st.density_decades > 0 || _config_error(
        "visualization.style.density_decades",
        "must be > 0; got $(st.density_decades)",
    )
    (0 < st.density_mass_frac ≤ 1) || _config_error(
        "visualization.style.density_mass_frac",
        "must satisfy 0 < density_mass_frac ≤ 1; got $(st.density_mass_frac)",
    )

    # [merger]
    cfg.merger.enabled &&
        isempty(cfg.merger.config_file) &&
        _config_error("merger.config_file", "must be nonempty when merger.enabled = true")

    # [cluster] / [cost]
    _validate_cluster(cfg)
    _validate_cost(cfg)

    return nothing
end

"""
    _validate_estimators(pp::PostprocessConfig)

The estimator checks of [`_validate`](@ref) for the `[postprocess]` keys
that set the shrinking-sphere centre, the Lagrangian grid, the radial
profiles, the local density, the mass segregation, the rotation analysis,
the shape tensor, the snapshot pair finder, the relaxation times and the
thresholds of the validity record. Raises
an `ArgumentError` naming the offending `postprocess.<key>`.
"""
function _validate_estimators(pp::PostprocessConfig)
    pp.centre_n_stop ≥ 1 ||
        _config_error("postprocess.centre_n_stop", "must be ≥ 1; got $(pp.centre_n_stop)")
    (0 < pp.centre_stop_fraction ≤ 1) || _config_error(
        "postprocess.centre_stop_fraction",
        "must satisfy 0 < centre_stop_fraction ≤ 1; got $(pp.centre_stop_fraction)",
    )
    pp.centre_n_floor ≥ 1 ||
        _config_error("postprocess.centre_n_floor", "must be ≥ 1; got $(pp.centre_n_floor)")
    (0 < pp.centre_shrink < 1) || _config_error(
        "postprocess.centre_shrink",
        "must satisfy 0 < centre_shrink < 1; got $(pp.centre_shrink)",
    )
    pp.lagr_n_min ≥ 1 ||
        _config_error("postprocess.lagr_n_min", "must be ≥ 1; got $(pp.lagr_n_min)")
    pp.profile_nbins ≥ 2 ||
        _config_error("postprocess.profile_nbins", "must be ≥ 2; got $(pp.profile_nbins)")
    pp.profile_n_min_sigma ≥ 2 || _config_error(
        "postprocess.profile_n_min_sigma",
        "must be ≥ 2; got $(pp.profile_n_min_sigma)",
    )
    pp.profile_n_min_beta ≥ pp.profile_n_min_sigma || _config_error(
        "postprocess.profile_n_min_beta",
        "must be ≥ profile_n_min_sigma = $(pp.profile_n_min_sigma); got $(pp.profile_n_min_beta)",
    )
    pp.profile_n_min_density ≥ 1 || _config_error(
        "postprocess.profile_n_min_density",
        "must be ≥ 1; got $(pp.profile_n_min_density)",
    )
    pp.density_neighbours ≥ 2 ||
        _config_error("postprocess.density_neighbours", "must be ≥ 2; got $(pp.density_neighbours)")
    pp.msr_n_tracer ≥ 2 ||
        _config_error("postprocess.msr_n_tracer", "must be ≥ 2; got $(pp.msr_n_tracer)")
    isempty(pp.msr_scan) && _config_error("postprocess.msr_scan", "must be nonempty")
    all(≥(2), pp.msr_scan) ||
        _config_error("postprocess.msr_scan", "entries must be ≥ 2; got $(pp.msr_scan)")
    all(diff(pp.msr_scan) .> 0) ||
        _config_error("postprocess.msr_scan", "must be strictly increasing; got $(pp.msr_scan)")
    pp.msr_n_random ≥ 2 ||
        _config_error("postprocess.msr_n_random", "must be ≥ 2; got $(pp.msr_n_random)")
    pp.msr_n_random_small ≥ 2 ||
        _config_error("postprocess.msr_n_random_small", "must be ≥ 2; got $(pp.msr_n_random_small)")
    pp.msr_n_small ≥ 0 ||
        _config_error("postprocess.msr_n_small", "must be ≥ 0; got $(pp.msr_n_small)")
    pp.msr_seed ≥ 0 || _config_error("postprocess.msr_seed", "must be ≥ 0; got $(pp.msr_seed)")
    pp.msr_threshold > 0 ||
        _config_error("postprocess.msr_threshold", "must be > 0; got $(pp.msr_threshold)")
    pp.rotation_n_min ≥ 2 ||
        _config_error("postprocess.rotation_n_min", "must be ≥ 2; got $(pp.rotation_n_min)")
    pp.rotation_nbins_max ≥ 1 ||
        _config_error("postprocess.rotation_nbins_max", "must be ≥ 1; got $(pp.rotation_nbins_max)")
    pp.rotation_n_null ≥ 0 ||
        _config_error("postprocess.rotation_n_null", "must be ≥ 0; got $(pp.rotation_n_null)")
    pp.shape_n_min_shell ≥ 1 ||
        _config_error("postprocess.shape_n_min_shell", "must be ≥ 1; got $(pp.shape_n_min_shell)")
    pp.shape_n_min_enclosed ≥ 1 || _config_error(
        "postprocess.shape_n_min_enclosed",
        "must be ≥ 1; got $(pp.shape_n_min_enclosed)",
    )
    pp.pair_k ≥ 1 || _config_error("postprocess.pair_k", "must be ≥ 1; got $(pp.pair_k)")
    pp.pair_tidal_factor > 0 ||
        _config_error("postprocess.pair_tidal_factor", "must be > 0; got $(pp.pair_tidal_factor)")
    isempty(pp.relaxation_gamma) &&
        _config_error("postprocess.relaxation_gamma", "must be nonempty")
    all(>(0), pp.relaxation_gamma) || _config_error(
        "postprocess.relaxation_gamma",
        "entries must be > 0; got $(pp.relaxation_gamma)",
    )
    isempty(pp.mass_class_edges) &&
        _config_error("postprocess.mass_class_edges", "must be nonempty")
    (all(>(0), pp.mass_class_edges) && all(diff(pp.mass_class_edges) .> 0)) || _config_error(
        "postprocess.mass_class_edges",
        "entries must be > 0 and strictly increasing; got $(pp.mass_class_edges)",
    )
    pp.mass_class_n_min ≥ 2 ||
        _config_error("postprocess.mass_class_n_min", "must be ≥ 2; got $(pp.mass_class_n_min)")
    pp.equipartition_bin_width > 0 || _config_error(
        "postprocess.equipartition_bin_width",
        "must be > 0; got $(pp.equipartition_bin_width)",
    )
    pp.equipartition_n_min ≥ 2 || _config_error(
        "postprocess.equipartition_n_min",
        "must be ≥ 2; got $(pp.equipartition_n_min)",
    )
    (isfinite(pp.span_n_periods) && pp.span_n_periods > 0) || _config_error(
        "postprocess.span_n_periods",
        "must be finite and > 0; got $(pp.span_n_periods)",
    )
    (isfinite(pp.bound_gamma_max) && pp.bound_gamma_max > 0) || _config_error(
        "postprocess.bound_gamma_max",
        "must be finite and > 0; got $(pp.bound_gamma_max)",
    )
    pp.bound_n_miss ≥ 0 ||
        _config_error("postprocess.bound_n_miss", "must be ≥ 0; got $(pp.bound_n_miss)")
    pp.population_frame in ("lagrangian", "rh") || _config_error(
        "postprocess.population_frame",
        "must be one of \"lagrangian\", \"rh\"; got \"$(pp.population_frame)\"",
    )
    0 < pp.population_q_min < 1 || _config_error(
        "postprocess.population_q_min",
        "must lie in (0, 1); got $(pp.population_q_min)",
    )
    0 ≤ pp.population_m1_min ||
        _config_error("postprocess.population_m1_min", "must be ≥ 0; got $(pp.population_m1_min)")
    pp.population_m1_min < pp.population_m1_max || _config_error(
        "postprocess.population_m1_max",
        "must exceed population_m1_min = $(pp.population_m1_min); got $(pp.population_m1_max)",
    )
    pp.validity_interval_de > 0 || _config_error(
        "postprocess.validity_interval_de",
        "must be > 0; got $(pp.validity_interval_de)",
    )
    pp.validity_cumulative_de > 0 || _config_error(
        "postprocess.validity_cumulative_de",
        "must be > 0; got $(pp.validity_cumulative_de)",
    )
    pp.validity_mass_tolerance > 0 || _config_error(
        "postprocess.validity_mass_tolerance",
        "must be > 0; got $(pp.validity_mass_tolerance)",
    )
    pp.validity_settle_tcr ≥ 0 || _config_error(
        "postprocess.validity_settle_tcr",
        "must be ≥ 0; got $(pp.validity_settle_tcr)",
    )
    pp.validity_pericentre_tcr ≥ 0 || _config_error(
        "postprocess.validity_pericentre_tcr",
        "must be ≥ 0; got $(pp.validity_pericentre_tcr)",
    )
    pp.validity_window_tcr > 0 || _config_error(
        "postprocess.validity_window_tcr",
        "must be > 0; got $(pp.validity_window_tcr)",
    )
    pp.validity_virial_tolerance > 0 || _config_error(
        "postprocess.validity_virial_tolerance",
        "must be > 0; got $(pp.validity_virial_tolerance)",
    )
    pp.validity_rate_tolerance > 0 || _config_error(
        "postprocess.validity_rate_tolerance",
        "must be > 0; got $(pp.validity_rate_tolerance)",
    )
    pp.validity_bias_tolerance > 0 || _config_error(
        "postprocess.validity_bias_tolerance",
        "must be > 0; got $(pp.validity_bias_tolerance)",
    )
    return nothing
end

"""
    _slurm_seconds(s) -> Int

Seconds of a Slurm time specification: `minutes`, `minutes:seconds`,
`hours:minutes:seconds`, `days-hours`, `days-hours:minutes` or
`days-hours:minutes:seconds` (the forms `sbatch --time` accepts). Raises an
`ArgumentError` naming the string for anything else.
"""
function _slurm_seconds(s::AbstractString)::Int
    text = strip(s)
    bad() = throw(
        ArgumentError(
            "\"$s\" is not a Slurm time (minutes, mm:ss, hh:mm:ss, d-hh, d-hh:mm, d-hh:mm:ss)",
        ),
    )
    field(x) = occursin(r"^\d+$", x) ? parse(Int, x) : bad()
    if occursin('-', text)
        days_s, rest = split(text, '-'; limit = 2)
        days = field(days_s)
        parts = split(rest, ':')
        1 ≤ length(parts) ≤ 3 || bad()
        hours = field(parts[1])
        minutes = length(parts) ≥ 2 ? field(parts[2]) : 0
        seconds = length(parts) == 3 ? field(parts[3]) : 0
        return ((days * 24 + hours) * 60 + minutes) * 60 + seconds
    end
    parts = split(text, ':')
    length(parts) == 1 && return field(parts[1]) * 60
    length(parts) == 2 && return field(parts[1]) * 60 + field(parts[2])
    length(parts) == 3 && return (field(parts[1]) * 60 + field(parts[2])) * 60 + field(parts[3])
    return bad()
end

"""
    _slurm_time(seconds) -> String

A duration in Slurm's `days-hh:mm:ss` form (`hh:mm:ss` below one day), the
inverse of [`_slurm_seconds`](@ref) for messages. `seconds` must be ≥ 0.
"""
function _slurm_time(seconds::Integer)::String
    seconds ≥ 0 || throw(ArgumentError("a duration cannot be negative; got $seconds s"))
    days, rem_d = divrem(seconds, 86_400)
    h, rem_h = divrem(rem_d, 3600)
    m, s = divrem(rem_h, 60)
    hms = @sprintf("%02d:%02d:%02d", h, m, s)
    return days == 0 ? hms : "$(days)-$(hms)"
end

const _SLURM_MEM_PATTERN = r"^\d+[KMGT]?$"

"""
    _validate_cluster(cfg::Nbody6Config)

The `[cluster]` checks of [`_validate`](@ref), applied only when
`cluster.enabled`: the scheduler, the time formats, the counts, the chain
mode, and the fit of one segment — `simulation.wall_budget` must be set and
`wall_budget + stop_margin + overhead` must not exceed `time_limit`, which in
turn must not exceed `partition_limit` when that is given. Raises an
`ArgumentError` naming the offending `cluster.<key>`.
"""
function _validate_cluster(cfg::Nbody6Config)
    cl = cfg.cluster
    cl.enabled || return nothing
    sim = cfg.simulation
    cl.scheduler == "slurm" ||
        _config_error("cluster.scheduler", "must be \"slurm\"; got \"$(cl.scheduler)\"")
    limit = try
        _slurm_seconds(cl.time_limit)
    catch e
        e isa ArgumentError || rethrow()
        _config_error("cluster.time_limit", e.msg)
    end
    if !isempty(cl.partition_limit)
        plimit = try
            _slurm_seconds(cl.partition_limit)
        catch e
            e isa ArgumentError || rethrow()
            _config_error("cluster.partition_limit", e.msg)
        end
        limit ≤ plimit || _config_error(
            "cluster.time_limit",
            "= $(cl.time_limit) exceeds the partition limit $(cl.partition_limit)",
        )
    end
    cl.nodes ≥ 1 || _config_error("cluster.nodes", "must be ≥ 1; got $(cl.nodes)")
    cl.ntasks ≥ 1 || _config_error("cluster.ntasks", "must be ≥ 1; got $(cl.ntasks)")
    cl.cpus_per_task ≥ 1 ||
        _config_error("cluster.cpus_per_task", "must be ≥ 1; got $(cl.cpus_per_task)")
    cl.julia_threads ≥ 0 ||
        _config_error("cluster.julia_threads", "must be ≥ 0; got $(cl.julia_threads)")
    cl.overhead ≥ 0 || _config_error("cluster.overhead", "must be ≥ 0; got $(cl.overhead)")
    cl.max_segments ≥ 0 ||
        _config_error("cluster.max_segments", "must be ≥ 0; got $(cl.max_segments)")
    cl.chain in ("resubmit", "dependency") || _config_error(
        "cluster.chain",
        "must be one of \"resubmit\", \"dependency\"; got \"$(cl.chain)\"",
    )
    cl.chain == "dependency" &&
        cl.max_segments < 1 &&
        _config_error(
            "cluster.max_segments",
            "must be ≥ 1 with chain = \"dependency\" (the chain length is fixed at submission)",
        )
    isempty(cl.mem) ||
        occursin(_SLURM_MEM_PATTERN, cl.mem) ||
        _config_error(
            "cluster.mem",
            "must be a size with an optional K, M, G or T suffix (e.g. \"64G\"); got \"$(cl.mem)\"",
        )
    for opt in cl.sbatch_options
        isempty(strip(opt)) && _config_error("cluster.sbatch_options", "entries must be nonempty")
        startswith(strip(opt), '#') && _config_error(
            "cluster.sbatch_options",
            "entries are sbatch options without the \"#SBATCH\" prefix; got \"$opt\"",
        )
    end
    isempty(cl.sbatch) && _config_error("cluster.sbatch", "must be nonempty")
    sim.wall_budget > 0 || _config_error(
        "cluster.enabled",
        "requires simulation.wall_budget > 0: a cluster job stops its engine segment through the " *
        "wall budget, below the job's time limit",
    )
    need = sim.wall_budget + sim.stop_margin + cl.overhead
    need ≤ limit || _config_error(
        "cluster.time_limit",
        "= $(cl.time_limit) ($(limit) s) is too short for one segment: wall_budget " *
        "$(sim.wall_budget) s + stop_margin $(sim.stop_margin) s + overhead $(cl.overhead) s = " *
        "$(need) s; lower simulation.wall_budget or raise the limit",
    )
    return nothing
end

"""
    _validate_cost(cfg::Nbody6Config)

The `[cost]` checks of [`_validate`](@ref): a replacement table of ≥ 2 rows
`[N, s_per_nb]` with strictly increasing `N` and positive values, a positive
`merger_factor`, a non-negative `budget_hours` and an `on_exceed` of
`"warn"` or `"refuse"`. Raises an `ArgumentError` naming the offending
`cost.<key>`.
"""
function _validate_cost(cfg::Nbody6Config)
    c = cfg.cost
    if !isempty(c.table)
        length(c.table) ≥ 2 || _config_error(
            "cost.table",
            "needs at least two rows [N, s_per_nb] to interpolate; got $(length(c.table))",
        )
        for row in c.table
            length(row) == 2 ||
                _config_error("cost.table", "rows must be [N, s_per_nb] pairs; got $row")
            all(>(0), row) || _config_error("cost.table", "entries must be > 0; got $row")
        end
        ns = first.(c.table)
        all(diff(ns) .> 0) ||
            _config_error("cost.table", "rows must be in strictly increasing N; got N = $ns")
    end
    c.merger_factor > 0 ||
        _config_error("cost.merger_factor", "must be > 0; got $(c.merger_factor)")
    c.budget_hours ≥ 0 || _config_error("cost.budget_hours", "must be ≥ 0; got $(c.budget_hours)")
    c.on_exceed in ("warn", "refuse") || _config_error(
        "cost.on_exceed",
        "must be one of \"warn\", \"refuse\"; got \"$(c.on_exceed)\"",
    )
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
        cfg.cluster,
        cfg.cost,
        String(abspath(base_dir)),
    )
end
