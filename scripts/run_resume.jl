#!/usr/bin/env julia
# =============================================================================
# Nbody6Dynamics — Continuation of an interrupted run
# =============================================================================
#
# Usage:
#   julia scripts/run_resume.jl <run_dir>
#
# Continues an interrupted run from the configuration frozen in its directory
# and performs post-processing and figures once the run is complete. The
# script can be submitted once per queue segment.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie   # loads the figure routines (package extension)
using Nbody6Dynamics

# SIGINT raises `InterruptException`, which the pipeline turns into the
# termination of the engine (a plain script would exit at once).
Base.exit_on_sigint(false)

function main()
    if length(ARGS) != 1
        println(stderr, "usage: julia scripts/run_resume.jl <run_dir>")
        exit(2)
    end
    resume_pipeline(ARGS[1])
    return nothing
end

main()
