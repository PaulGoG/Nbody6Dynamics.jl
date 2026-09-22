#!/usr/bin/env julia
# =============================================================================
# Nbody6Dynamics — Main entry point
# =============================================================================
#
# Usage:
#   julia scripts/run_setup.jl [config.toml]
#
# If no config path is given, it defaults to config.toml in the project root.
#
# This script is a thin wrapper around `run_pipeline`, which is the single
# orchestrator handling every phase (install → merger ICs → simulation →
# post-processing → plots) from the config flags. Phase logic stays there, so
# every entry point selects the same phases from the same flags.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie   # loads the figure routines (package extension)
using Nbody6Dynamics

function main()
    config_path = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "config.toml")

    if !isfile(config_path)
        @error "Configuration file not found" path = config_path
        exit(1)
    end

    @info "Loading configuration from: $config_path"
    cfg = load_config(config_path)

    # Relative paths of the configuration resolve against its own directory
    # (cfg.config_dir), which is where backend/ and runs/ are created.
    t0 = time()
    results = run_pipeline(cfg)

    for (k, v) in results
        n = v isa AbstractVector ? length(v) : v isa DiagnosticsData ? length(v.adjust) : 1
        @info "  Loaded: $k ($n records)"
    end
    @info "Done. Total pipeline time: $(Nbody6Dynamics._format_elapsed(time() - t0))"
end

main()
