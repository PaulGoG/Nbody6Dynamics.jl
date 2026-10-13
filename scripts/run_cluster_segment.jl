#!/usr/bin/env julia
# =============================================================================
# Nbody6Dynamics — One engine segment of a run, within a Slurm job
# =============================================================================
#
# Usage:
#   julia scripts/run_cluster_segment.jl <run_dir>
#
# Run by the job script that scripts/run_cluster.jl (submit_run) generates:
# starts or resumes the run in <run_dir>, records the job in the ledgers, and
# with chain = "resubmit" submits the job of the next segment.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie   # loads the figure routines (package extension)
using Nbody6Dynamics

# SIGINT raises `InterruptException`, which the pipeline turns into the
# termination of the engine (a plain script would exit at once).
Base.exit_on_sigint(false)

function main()
    if length(ARGS) != 1
        println(stderr, "usage: julia scripts/run_cluster_segment.jl <run_dir>")
        exit(2)
    end
    cluster_segment(ARGS[1])
    return nothing
end

main()
