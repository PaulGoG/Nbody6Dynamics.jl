#!/usr/bin/env julia
# =============================================================================
# Verification suite runner — launches three targets in-process via
# run_pipeline(cfg), the only entry point that handles merger mode correctly.
#
#   julia --project=. scripts/run_verif_suite.jl
# =============================================================================

const PROJ = normpath(joinpath(@__DIR__, ".."))

using Pkg
Pkg.activate(PROJ; io = devnull)

using Nbody6Dynamics, TOML, Dates

cd(PROJ)

function make_base_cfg()
    Dict(
        "install" => Dict(
            "enabled" => false,
            "reinstall" => false,
            "clean_build" => true,
            "source_url" => "https://github.com/nbody6ppgpu/Nbody6PPGPU-beijing.git",
            "install_dir" => "backend/Nbody6PPGPU-beijing",
        ),
        "build" => Dict(
            "enable_gpu" => false,
            "cuda_path" => "",
            "configure_flags" => ["--enable-mcmodel=large", "--with-par=b1m"],
            "nproc" => 0,
            "enable_hdf5" => true,
            "enable_mpi" => false,
        ),
        "simulation" => Dict(
            "run_test" => true,
            "runs_dir" => "runs",
            "binary_name" => "nbody6++",
            "mpi_ranks" => 1,
            "input_file" => "../../input_files/N100k_production.inp",
            "run_id_prefix" => "verif",
        ),
        "postprocess" => Dict(
            "enabled" => true,
            "data_dir" => "",
            "parse_stdout" => true,
            "stdout_file" => "out1000",
            "read_lagr" => true,
            "lagr_file" => "lagr.7",
            "read_escapers" => true,
            "escapers_file" => "esc.11",
            "read_stellar_evo" => true,
            "stellar_evo_pattern" => "sev.83_*",
            "snapshot_format" => "conf3",
            "snapshot_pattern" => "conf.3_*",
        ),
        "visualization" => Dict(
            "enabled" => true,
            "format" => "png",
            "dpi" => 300,
            "output_dir" => "plots",
            "figsize" => [8, 6],
        ),
    )
end

function cfg_for_single(inp_path::String, prefix::String)
    cfg = make_base_cfg()
    cfg["merger"] = Dict("enabled" => false, "config_file" => "")
    cfg["simulation"]["input_file"] = inp_path
    cfg["simulation"]["run_id_prefix"] = prefix
    cfg
end

function cfg_for_merger(merger_toml::String, prefix::String)
    cfg = make_base_cfg()
    cfg["merger"] = Dict("enabled" => true, "config_file" => merger_toml)
    cfg["simulation"]["run_id_prefix"] = prefix
    cfg
end

const SUITE = [
    (name = "single", raw = cfg_for_single("../../input_files/N5k_medium.inp", "verif_single")),
    (name = "triorbit", raw = cfg_for_merger("input_files/verif_triorbit.toml", "verif_triorbit")),
    (name = "3d5cluster", raw = cfg_for_merger("input_files/verif_3d5cluster.toml", "verif_3d5")),
]

# Intermediate configs go to a scratch dir (the TOML round-trip through
# load_config is deliberate — it exercises the parser). The authoritative
# frozen config for each run is written into its run directory by the
# pipeline itself; nothing is written to the project root.
mktempdir() do scratch
    for entry in SUITE
        @info "─────────────────────────────────────────────────────"
        @info "Running verification target: $(entry.name)"
        cfg_path = joinpath(scratch, "config.verif_$(entry.name).toml")
        open(cfg_path, "w") do io
            TOML.print(io, entry.raw)
        end
        t0 = time()
        try
            cfg = Nbody6Dynamics.load_config(cfg_path)
            Nbody6Dynamics.run_pipeline(cfg; base_dir = PROJ)
            @info "✓ $(entry.name) completed in $(round(time() - t0, digits=1))s"
        catch e
            @error "✗ $(entry.name) failed" exception = (e, catch_backtrace())
        end
    end
end

@info "Verification suite done."
