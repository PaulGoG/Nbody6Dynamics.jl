#!/usr/bin/env julia
# =============================================================================
# Nbody6Dynamics — Main entry point
# =============================================================================
#
# Usage:
#   julia --project=path/to/Nbody6Dynamics scripts/run_setup.jl [config.toml]
#
# If no config path is given, it defaults to config.toml in the project root.
#
# This script is a thin wrapper around `run_pipeline`, which is the single
# orchestrator handling every phase (install → merger ICs → simulation →
# post-processing → plots) from the config flags. Do not re-implement phase
# logic here — a previous copy of this script did, silently skipping the
# merger phase.

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

    t0 = time()
    results = run_pipeline(cfg; base_dir = dirname(abspath(config_path)))

    for (k, v) in results
        n = v isa AbstractVector ? length(v) : v isa DiagnosticsData ? length(v.adjust) : 1
        @info "  Loaded: $k ($n records)"
    end
    @info "Done. Total pipeline time: $(Nbody6Dynamics._format_elapsed(time() - t0))"
end

main()
