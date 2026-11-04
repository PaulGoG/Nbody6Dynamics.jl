#!/usr/bin/env julia
# =============================================================================
# Nbody6Dynamics — Submission of a run as a chain of Slurm jobs
# =============================================================================
#
# Usage:
#   julia scripts/run_cluster.jl config.toml [--dry-run] [--run-id=<id>]
#   julia scripts/run_cluster.jl <run_dir> --restart [--dump=<file>]
#         [--t-end=<NB> | --t-end-myr=<Myr> | --tcrit-extra=<NB>] [--kz=<i>:<v>,...]
#         [--engine=backend] [--dry-run]
#
# First form: creates the run directory, freezes the configuration and the job
# script into <run_dir>/cluster/, and submits the job chain the [cluster]
# section describes (cluster.enabled must be true). --dry-run writes the run's
# files and prints the submission commands without submitting; --run-id pins
# the name of the run directory.
#
# Second form (restart job): submits a chain whose first segment restarts the
# existing run from a dump (`submit_restart`): an extension of a completed run
# to a later end time, or the continuation of a halted or stalled run from an
# earlier dump, with options changed by --kz. Without --dump the dump a resume
# would take; without an end time the run's own end time. --dry-run writes the
# request and the restart input as cluster/restart_preview.{toml,inp}.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie   # loads the figure routines (package extension)
using Nbody6Dynamics

const USAGE = """
    usage: julia scripts/run_cluster.jl config.toml [--dry-run] [--run-id=<id>]
           julia scripts/run_cluster.jl <run_dir> --restart [--dump=<file>]
                 [--t-end=<NB> | --t-end-myr=<Myr> | --tcrit-extra=<NB>] [--kz=<i>:<v>,...]
                 [--engine=backend] [--dry-run]"""

function _option_value(opt::AbstractString, name::AbstractString)
    return opt[(length(name) + 2):end]   # "--name=" prefix
end

function _parse_kz(text::AbstractString)::Dict{Int,Int}
    kz = Dict{Int,Int}()
    for item in split(text, ',')
        parts = split(item, ':')
        length(parts) == 2 || error("--kz: expected <index>:<value> pairs, got \"$item\"")
        kz[parse(Int, parts[1])] = parse(Int, parts[2])
    end
    return kz
end

function main()
    options = filter(startswith("--"), ARGS)
    positional = filter(!startswith("--"), ARGS)
    dry_run = false
    restart = false
    run_id = ""
    dump = nothing
    t_end = nothing
    t_end_myr = nothing
    tcrit_extra = nothing
    kz = nothing
    engine = :run
    for opt in options
        if opt == "--dry-run"
            dry_run = true
        elseif opt == "--restart"
            restart = true
        elseif startswith(opt, "--run-id=")
            run_id = _option_value(opt, "run-id")
        elseif startswith(opt, "--dump=")
            dump = _option_value(opt, "dump")
        elseif startswith(opt, "--t-end=")
            t_end = parse(Float64, _option_value(opt, "t-end"))
        elseif startswith(opt, "--t-end-myr=")
            t_end_myr = parse(Float64, _option_value(opt, "t-end-myr"))
        elseif startswith(opt, "--tcrit-extra=")
            tcrit_extra = parse(Float64, _option_value(opt, "tcrit-extra"))
        elseif startswith(opt, "--kz=")
            kz = _parse_kz(_option_value(opt, "kz"))
        elseif opt == "--engine=backend"
            engine = :backend
        else
            println(stderr, USAGE)
            exit(2)
        end
    end
    restart_only = (dump, t_end, t_end_myr, tcrit_extra, kz)
    if length(positional) != 1 ||
       (restart && !isempty(run_id)) ||
       (!restart && (any(!isnothing, restart_only) || engine != :run))
        println(stderr, USAGE)
        exit(2)
    end
    target = only(positional)

    if restart
        isdir(target) || (@error "Run directory not found" path = target; exit(1))
        sub = submit_restart(
            target;
            t_end = t_end,
            t_end_myr = t_end_myr,
            tcrit_extra = tcrit_extra,
            dump = dump,
            kz = kz,
            engine = engine,
            dry_run = dry_run,
        )
        r = sub.request
        println("Run directory: ", sub.run_dir)
        println("Restart:       from $(r["dump"]) (t = $(r["t_start"]) NB) to t = $(r["t_end"]) NB")
        println("Job script:    ", sub.script)
        println(dry_run ? "Job (not submitted): " : "Job $(sub.job_id): ", sub.command)
        return nothing
    end

    if !isfile(target)
        @error "Configuration file not found" path = target
        exit(1)
    end
    @info "Loading configuration from: $target"
    cfg = load_config(target)
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
