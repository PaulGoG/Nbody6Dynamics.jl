# =============================================================================
# Configuration loading and serialization
# =============================================================================

"""
    load_config(path::AbstractString) -> Nbody6Config

Parse a TOML configuration file and return a fully-typed `Nbody6Config`.
Missing sections or keys fall back to defaults defined in the `@kwdef` structs.
The parsed configuration is validated fail-fast by [`_validate`](@ref):
enumerated choices, numerical bounds, and required non-empty fields raise an
error naming the offending `section.key` before any pipeline phase starts.
The directory of `path` is recorded as `config_dir`: relative paths in the
configuration resolve against it, and it is the default `base_dir` of the
entry points, so a project directory holding a `config.toml` receives its
`backend/` and `runs/` there, wherever the package is installed.
"""
function load_config(path::AbstractString)::Nbody6Config
    raw = TOML.parsefile(path)
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
    open(path, "w") do io
        TOML.print(io, d)
    end
end

# ---------------------------------------------------------------------------
# Section parsers — convert Dict{String,Any} → typed struct
# ---------------------------------------------------------------------------

function _parse_install(d::Dict)
    InstallConfig(;
        enabled = get(d, "enabled", true),
        source_url = get(d, "source_url", "https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing.git"),
        ref = get(d, "ref", "618d7a4"),
        install_dir = get(d, "install_dir", joinpath("backend", "Nbody6PPGPU-beijing")),
        reinstall = get(d, "reinstall", false),
        clean_build = get(d, "clean_build", true),
    )
end

function _parse_build(d::Dict)
    flags = get(d, "configure_flags", ["--enable-mcmodel=large", "--with-par=b1m"])
    BuildConfig(;
        configure_flags = String.(flags),
        enable_mpi = get(d, "enable_mpi", false),
        enable_hdf5 = get(d, "enable_hdf5", true),
        enable_gpu = get(d, "enable_gpu", false),
        cuda_path = get(d, "cuda_path", ""),
        cuda_arch = Vector{String}(get(d, "cuda_arch", String[])),
        nvcc_flags = Vector{String}(get(d, "nvcc_flags", String[])),
        nproc = get(d, "nproc", 0),
    )
end

function _parse_simulation(d::Dict)
    SimulationConfig(;
        run_test = get(d, "run_test", true),
        input_file = get(d, "input_file", "examples/input_files/N10k_noDat10.inp"),
        runs_dir = get(d, "runs_dir", "runs"),
        binary_name = get(d, "binary_name", "nbody6++"),
        mpi_ranks = get(d, "mpi_ranks", 1),
        omp_threads = get(d, "omp_threads", 0),
        gpu_list = Vector{Int}(get(d, "gpu_list", Int[])),
        run_id_prefix = get(d, "run_id_prefix", "run"),
        monitor = get(d, "monitor", false),
        live_diagnostics = get(d, "live_diagnostics", false),
        live_interval = Float64(get(d, "live_interval", 30.0)),
        telemetry_interval = Float64(get(d, "telemetry_interval", 5.0)),
        startup_timeout = Float64(get(d, "startup_timeout", 0.0)),
        exit_grace = Float64(get(d, "exit_grace", 120.0)),
    )
end

function _parse_postprocess(d::Dict)
    PostprocessConfig(;
        enabled = get(d, "enabled", true),
        data_dir = get(d, "data_dir", ""),
        snapshot_format = get(d, "snapshot_format", "conf3"),
        snapshot_pattern = get(d, "snapshot_pattern", "conf.3_*"),
        parse_stdout = get(d, "parse_stdout", true),
        stdout_file = get(d, "stdout_file", "out1000"),
        read_lagr = get(d, "read_lagr", true),
        lagr_file = get(d, "lagr_file", "lagr.7"),
        read_escapers = get(d, "read_escapers", true),
        escapers_file = get(d, "escapers_file", "esc.11"),
        read_stellar_evo = get(d, "read_stellar_evo", true),
        stellar_evo_pattern = get(d, "stellar_evo_pattern", "sev.83_*"),
        read_binary_evo = get(d, "read_binary_evo", true),
        binary_evo_pattern = get(d, "binary_evo_pattern", "bev.82_*"),
    )
end

function _parse_merger_pipeline(d::Dict)
    MergerPipelineConfig(;
        enabled = get(d, "enabled", false),
        config_file = get(d, "config_file", ""),
    )
end

function _parse_visualization(d::Dict)
    st = get(d, "style", Dict{String,Any}())
    style = PlotStyle(;
        marker_budget = Float64(get(st, "marker_budget", 27000.0)),
        marker_min = Float64(get(st, "marker_min", 6.0)),
        marker_max = Float64(get(st, "marker_max", 30.0)),
        q_log_threshold = Float64(get(st, "q_log_threshold", 10.0)),
        q_floor = Float64(get(st, "q_floor", 1e-3)),
        zoom_frac = Float64(get(st, "zoom_frac", 0.15)),
        anim_fps = Int(get(st, "anim_fps", 0)),
        anim_target_seconds = Float64(get(st, "anim_target_seconds", 12.0)),
        anim_px_per_unit = Float64(get(st, "anim_px_per_unit", 1.4)),
    )
    VisualizationConfig(;
        enabled = get(d, "enabled", true),
        format = get(d, "format", "pdf"),
        dpi = get(d, "dpi", 300),
        export_width = Float64(get(d, "export_width", 6.5)),
        column = get(d, "column", ""),
        units = get(d, "units", "physical"),
        output_dir = get(d, "output_dir", "plots"),
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
required non-empty fields. Raises an `ErrorException` naming the offending
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
    isempty(inst.install_dir) && error("config: install.install_dir must be nonempty")
    bld.nproc ≥ 0 || error("config: build.nproc must be ≥ 0; got $(bld.nproc)")
    for arch in bld.cuda_arch
        occursin(_CUDA_ARCH_PATTERN, arch) || error(
            "config: build.cuda_arch entries must be CUDA architecture names of the form " *
            "\"sm_<major><minor>\" (e.g. \"sm_90\", \"sm_120\"); got \"$arch\"",
        )
    end
    for flag in bld.nvcc_flags
        isempty(strip(flag)) && error("config: build.nvcc_flags entries must be nonempty")
    end

    # [simulation]
    sim.mpi_ranks ≥ 1 || error("config: simulation.mpi_ranks must be ≥ 1; got $(sim.mpi_ranks)")
    if sim.mpi_ranks > 1 && !bld.enable_mpi
        error(
            "config: simulation.mpi_ranks = $(sim.mpi_ranks) requires build.enable_mpi = true — " *
            "a multi-rank launch needs an MPI-enabled binary",
        )
    end
    sim.omp_threads ≥ 0 ||
        error("config: simulation.omp_threads must be ≥ 0; got $(sim.omp_threads)")
    if !isempty(sim.gpu_list)
        bld.enable_gpu || error(
            "config: simulation.gpu_list = $(sim.gpu_list) requires build.enable_gpu = true — " *
            "the CPU binary ignores GPU_LIST",
        )
        all(≥(0), sim.gpu_list) || error(
            "config: simulation.gpu_list entries must be ≥ 0 (CUDA device indices); got $(sim.gpu_list)",
        )
        allunique(sim.gpu_list) ||
            error("config: simulation.gpu_list entries must be distinct; got $(sim.gpu_list)")
        length(sim.gpu_list) ≤ _MAX_GPU_PER_PROCESS || error(
            "config: simulation.gpu_list may name at most $(_MAX_GPU_PER_PROCESS) devices per " *
            "process (the engine's MAX_GPU); got $(length(sim.gpu_list))",
        )
    end
    sim.telemetry_interval ≥ 0 || error(
        "config: simulation.telemetry_interval must be ≥ 0 [s]; got $(sim.telemetry_interval)",
    )
    sim.startup_timeout ≥ 0 ||
        error("config: simulation.startup_timeout must be ≥ 0 [s]; got $(sim.startup_timeout)")
    sim.exit_grace ≥ 0 ||
        error("config: simulation.exit_grace must be ≥ 0 [s]; got $(sim.exit_grace)")
    sim.live_interval ≥ 1 ||
        error("config: simulation.live_interval must be ≥ 1 [s]; got $(sim.live_interval)")
    isempty(sim.runs_dir) && error("config: simulation.runs_dir must be nonempty")
    isempty(sim.run_id_prefix) && error("config: simulation.run_id_prefix must be nonempty")
    isempty(sim.input_file) && error("config: simulation.input_file must be nonempty")

    # [postprocess]
    if pp.snapshot_format != "conf3"
        if pp.snapshot_format == "hdf5"
            error(
                "config: postprocess.snapshot_format = \"hdf5\" is no longer supported; " *
                "use \"conf3\". The fork's KZ(46) H5Part layout was never readable by the old code.",
            )
        end
        error(
            "config: postprocess.snapshot_format must be \"conf3\" (the only supported value); " *
            "got \"$(pp.snapshot_format)\"",
        )
    end
    isempty(pp.snapshot_pattern) && error("config: postprocess.snapshot_pattern must be nonempty")
    pp.parse_stdout &&
        isempty(pp.stdout_file) &&
        error("config: postprocess.stdout_file must be nonempty when parse_stdout = true")
    pp.read_lagr &&
        isempty(pp.lagr_file) &&
        error("config: postprocess.lagr_file must be nonempty when read_lagr = true")
    pp.read_escapers &&
        isempty(pp.escapers_file) &&
        error("config: postprocess.escapers_file must be nonempty when read_escapers = true")
    pp.read_stellar_evo &&
        isempty(pp.stellar_evo_pattern) &&
        error(
            "config: postprocess.stellar_evo_pattern must be nonempty when read_stellar_evo = true",
        )
    pp.read_binary_evo &&
        isempty(pp.binary_evo_pattern) &&
        error("config: postprocess.binary_evo_pattern must be nonempty when read_binary_evo = true")

    # [visualization]
    vis.format in ("pdf", "svg", "png") || error(
        "config: visualization.format must be one of \"pdf\", \"svg\", \"png\"; got \"$(vis.format)\"",
    )
    vis.column in ("single", "double", "") || error(
        "config: visualization.column must be one of \"single\", \"double\", \"\"; got \"$(vis.column)\"",
    )
    vis.units in ("physical", "nbody") ||
        error("config: visualization.units must be \"physical\" or \"nbody\"; got \"$(vis.units)\"")
    vis.dpi ≥ 72 || error("config: visualization.dpi must be ≥ 72; got $(vis.dpi)")
    vis.export_width > 0 ||
        error("config: visualization.export_width must be > 0; got $(vis.export_width)")

    # [visualization.style]
    st.marker_budget > 0 ||
        error("config: visualization.style.marker_budget must be > 0; got $(st.marker_budget)")
    (0 < st.marker_min ≤ st.marker_max) || error(
        "config: visualization.style must satisfy 0 < marker_min ≤ marker_max; " *
        "got marker_min = $(st.marker_min), marker_max = $(st.marker_max)",
    )
    st.q_log_threshold > 0 ||
        error("config: visualization.style.q_log_threshold must be > 0; got $(st.q_log_threshold)")
    (0 < st.q_floor < 1) ||
        error("config: visualization.style.q_floor must satisfy 0 < q_floor < 1; got $(st.q_floor)")
    (0 < st.zoom_frac ≤ 1) || error(
        "config: visualization.style.zoom_frac must satisfy 0 < zoom_frac ≤ 1; got $(st.zoom_frac)",
    )
    st.anim_fps ≥ 0 || error("config: visualization.style.anim_fps must be ≥ 0; got $(st.anim_fps)")
    st.anim_px_per_unit > 0 || error(
        "config: visualization.style.anim_px_per_unit must be > 0; got $(st.anim_px_per_unit)",
    )
    st.anim_target_seconds > 0 || error(
        "config: visualization.style.anim_target_seconds must be > 0; got $(st.anim_target_seconds)",
    )

    # [merger]
    cfg.merger.enabled &&
        isempty(cfg.merger.config_file) &&
        error("config: merger.config_file must be nonempty when merger.enabled = true")

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
