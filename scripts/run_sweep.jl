#!/usr/bin/env julia
# =============================================================================
# Nbody6Dynamics — parameter sweep entry point
# =============================================================================
#
# Usage:
#   julia scripts/run_sweep.jl sweep.toml [--dry-run]
#
# Prepares one directory per grid point × seed, runs the points as
# concurrent worker processes (`run_sweep`), writes sweep_index.toml and
# sweep_summary.csv, and draws the comparison figures into <sweep>/plots.
# With --dry-run the directories and derived configs are written only.

include(joinpath(@__DIR__, "activate.jl"))

using CairoMakie   # loads the figure routines (package extension)
using Nbody6Dynamics

function main()
    args = filter(a -> !startswith(a, "--"), ARGS)
    dry_run = "--dry-run" in ARGS
    if length(args) != 1
        println(stderr, "usage: run_sweep.jl <sweep.toml> [--dry-run]")
        exit(2)
    end
    path = args[1]
    isfile(path) || (println(stderr, "sweep configuration not found: $path"); exit(1))

    cfg = load_sweep_config(path)
    t0 = time()
    sweep_dir = run_sweep(cfg; dry_run = dry_run)
    if dry_run
        @info "Sweep prepared (dry run): $sweep_dir"
        return
    end
    figures = sweep_figures(sweep_dir, sweep_visualization(cfg, sweep_dir))
    @info "Sweep complete in $(Nbody6Dynamics._format_elapsed(time() - t0)): $sweep_dir" figures
end

main()
