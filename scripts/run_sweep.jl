#!/usr/bin/env julia
# =============================================================================
# Nbody6Dynamics — parameter sweep entry point
# =============================================================================
#
# Usage:
#   julia scripts/run_sweep.jl sweep.toml [--dry-run] [--sweep-dir=<dir>] [--resume]
#
# Prepares one directory per grid point × seed, runs the points as
# concurrent worker processes (`run_sweep`), writes sweep_index.toml and
# sweep_summary.csv, and draws the comparison figures into <sweep>/plots.
# With --dry-run the directories and derived configs are written only.
# --sweep-dir names the sweep directory instead of a time-stamped one. With
# --resume an existing sweep directory is re-entered (`resume_sweep`: done
# points skipped, stopped and unstarted ones run), and a missing one is
# prepared, so one job script serves the first job of a sweep and every job
# that continues it.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie   # loads the figure routines (package extension)
using Nbody6Dynamics

# SIGINT raises `InterruptException` (a plain script would exit at once): the
# pipeline handler then terminates the engine, which runs in its own session.
Base.exit_on_sigint(false)

function main()
    usage = "usage: run_sweep.jl <sweep.toml> [--dry-run] [--sweep-dir=<dir>] [--resume]"
    args = filter(a -> !startswith(a, "--"), ARGS)
    dry_run = false
    resume = false
    sweep_dir = ""
    for opt in filter(startswith("--"), ARGS)
        if opt == "--dry-run"
            dry_run = true
        elseif opt == "--resume"
            resume = true
        elseif startswith(opt, "--sweep-dir=")
            sweep_dir = opt[(length("--sweep-dir=") + 1):end]
        else
            println(stderr, usage)
            exit(2)
        end
    end
    if length(args) != 1 || (resume && (isempty(sweep_dir) || dry_run))
        println(stderr, usage)
        resume && println(stderr, "--resume needs --sweep-dir and excludes --dry-run")
        exit(2)
    end
    path = args[1]
    isfile(path) || (println(stderr, "sweep configuration not found: $path"); exit(1))

    cfg = load_sweep_config(path)
    t0 = time()
    sweep_dir =
        resume && isfile(joinpath(sweep_dir, "sweep_index.toml")) ? resume_sweep(cfg, sweep_dir) :
        run_sweep(cfg; dry_run = dry_run, sweep_dir = sweep_dir)
    if dry_run
        @info "Sweep prepared (dry run): $sweep_dir"
        return
    end
    figures = sweep_figures(sweep_dir, sweep_visualization(cfg, sweep_dir))
    points = Nbody6Dynamics.read_sweep_index(sweep_dir)["points"]
    n_done = count(p -> p["status"] == "done", points)
    @info "Sweep entry finished in $(Nbody6Dynamics._format_elapsed(time() - t0)): " *
          "$(n_done) of $(length(points)) points done; $sweep_dir" figures
end

main()
