#!/usr/bin/env julia
# =============================================================================
# Nbody6Dynamics — Submission of a run as a chain of Slurm jobs
# =============================================================================
#
# Usage:
#   julia scripts/run_cluster.jl config.toml [--dry-run] [--run-id=<id>]
#
# Creates the run directory, freezes the configuration and the job script
# into <run_dir>/cluster/, and submits the job chain the [cluster] section
# describes (cluster.enabled must be true). --dry-run writes the run's files
# and prints the submission commands without submitting; --run-id pins the
# name of the run directory.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie   # loads the figure routines (package extension)
using Nbody6Dynamics

function main()
    usage = "usage: julia scripts/run_cluster.jl config.toml [--dry-run] [--run-id=<id>]"
    options = filter(startswith("--"), ARGS)
    positional = filter(!startswith("--"), ARGS)
    dry_run = false
    run_id = ""
    for opt in options
        if opt == "--dry-run"
            dry_run = true
        elseif startswith(opt, "--run-id=")
            run_id = opt[(length("--run-id=") + 1):end]
        else
            println(stderr, usage)
            exit(2)
        end
    end
    if length(positional) != 1
        println(stderr, usage)
        exit(2)
    end
    config_path = only(positional)
    if !isfile(config_path)
        @error "Configuration file not found" path = config_path
        exit(1)
    end

    @info "Loading configuration from: $config_path"
    cfg = load_config(config_path)
    sub = submit_run(cfg; run_id = run_id, dry_run = dry_run)

    println("Run directory: ", sub.run_dir)
    println("Job script:    ", sub.script)
    if dry_run
        for (k, command) in enumerate(sub.commands)
            println("Job $k (not submitted): ", command)
        end
    else
        for (id, command) in zip(sub.job_ids, sub.commands)
            println("Job $id: ", command)
        end
    end
    return nothing
end

main()
