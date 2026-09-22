#!/usr/bin/env julia
# =============================================================================
# Nbody6Dynamics — Main entry point
# =============================================================================
#
# Usage:
#   julia scripts/run_setup.jl [config.toml] [--run-id=<id>]
#
# If no config path is given, it defaults to config.toml in the project root.
# --run-id pins the name of the run directory instead of generating one, so a
# caller (the GPU validation driver) can find the artefacts it asked for.
#
# This script is a thin wrapper around `run_pipeline`, which is the single
# orchestrator handling every phase (install → merger ICs → simulation →
# post-processing → plots) from the config flags. Phase logic stays there, so
# every entry point selects the same phases from the same flags.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie   # loads the figure routines (package extension)
using Nbody6Dynamics

function main()
    options = filter(startswith("--"), ARGS)
    positional = filter(!startswith("--"), ARGS)
    run_id = ""
    for opt in options
        if startswith(opt, "--run-id=")
            run_id = opt[(length("--run-id=") + 1):end]
        else
            println(stderr, "usage: run_setup.jl [config.toml] [--run-id=<id>]")
            exit(2)
        end
    end
    config_path = isempty(positional) ? joinpath(@__DIR__, "..", "config.toml") : first(positional)

    if !isfile(config_path)
        @error "Configuration file not found" path = config_path
        exit(1)
    end

    @info "Loading configuration from: $config_path"
    cfg = load_config(config_path)

    # Relative paths of the configuration resolve against its own directory
    # (cfg.config_dir), which is where backend/ and runs/ are created.
    t0 = time()
    results = run_pipeline(cfg; run_id = run_id)

    for (k, v) in results
        n = v isa AbstractVector ? length(v) : v isa DiagnosticsData ? length(v.adjust) : 1
        @info "  Loaded: $k ($n records)"
    end
    @info "Done. Total pipeline time: $(Nbody6Dynamics._format_elapsed(time() - t0))"
end

main()
