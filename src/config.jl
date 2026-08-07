# =============================================================================
# Configuration loading and serialization
# =============================================================================

"""
    load_config(path::AbstractString) -> Nbody6Config

Parse a TOML configuration file and return a fully-typed `Nbody6Config`.
Missing sections or keys fall back to defaults defined in the `@kwdef` structs.
"""
function load_config(path::AbstractString)::Nbody6Config
    raw = TOML.parsefile(path)
    return Nbody6Config(
        _parse_install(get(raw, "install", Dict{String,Any}())),
        _parse_build(get(raw, "build", Dict{String,Any}())),
        _parse_simulation(get(raw, "simulation", Dict{String,Any}())),
        _parse_postprocess(get(raw, "postprocess", Dict{String,Any}())),
        _parse_visualization(get(raw, "visualization", Dict{String,Any}())),
        _parse_merger_pipeline(get(raw, "merger", Dict{String,Any}())),
    )
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
    # TOML cannot serialize Tuples or nested structs — convert explicitly
    d["visualization"]["figsize"] = collect(d["visualization"]["figsize"])
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
        run_id_prefix = get(d, "run_id_prefix", "run"),
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
    )
end

function _parse_merger_pipeline(d::Dict)
    MergerPipelineConfig(;
        enabled = get(d, "enabled", false),
        config_file = get(d, "config_file", ""),
    )
end

function _parse_visualization(d::Dict)
    fs = get(d, "figsize", [8.0, 6.0])
    st = get(d, "style", Dict{String,Any}())
    style = PlotStyle(;
        marker_budget = Float64(get(st, "marker_budget", 18000.0)),
        marker_min = Float64(get(st, "marker_min", 4.0)),
        marker_max = Float64(get(st, "marker_max", 20.0)),
        q_log_threshold = Float64(get(st, "q_log_threshold", 10.0)),
        q_floor = Float64(get(st, "q_floor", 1e-3)),
        zoom_frac = Float64(get(st, "zoom_frac", 0.15)),
        anim_fps = Int(get(st, "anim_fps", 0)),
        anim_target_seconds = Float64(get(st, "anim_target_seconds", 12.0)),
    )
    VisualizationConfig(;
        enabled = get(d, "enabled", true),
        format = get(d, "format", "pdf"),
        dpi = get(d, "dpi", 300),
        column = get(d, "column", "single"),
        figsize = (Float64(fs[1]), Float64(fs[2])),
        units = get(d, "units", "physical"),
        output_dir = get(d, "output_dir", "plots"),
        style = style,
    )
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _struct_to_dict(s)
    Dict(String(k) => getfield(s, k) for k in fieldnames(typeof(s)))
end
