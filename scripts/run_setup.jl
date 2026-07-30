#!/usr/bin/env julia
# =============================================================================
# Nbody6Setup — Main entry point
# =============================================================================
#
# Usage:
#   julia --project=path/to/Nbody6Setup scripts/run_setup.jl [config.toml]
#
# If no config path is given, it defaults to config.toml in the project root.

using Nbody6Setup

function main()
    # Resolve config path
    config_path = if length(ARGS) >= 1
        ARGS[1]
    else
        joinpath(@__DIR__, "..", "config.toml")
    end

    if !isfile(config_path)
        @error "Configuration file not found" path = config_path
        exit(1)
    end

    @info "Loading configuration from: $config_path"
    cfg = load_config(config_path)

    base_dir = dirname(abspath(config_path))
    t_pipeline = time()

    # ---- Phase 1: Install / Build ----
    if cfg.install.enabled
        t0 = time()
        @info "═══ Phase 1: Install & Build ═══"
        setup_nbody6(cfg; base_dir)
        @info "Phase 1 complete ($(Nbody6Setup._format_elapsed(time() - t0)))"
    else
        @info "Install phase disabled — skipping."
    end

    # ---- Phase 2: Test Run ----
    run_dir = ""
    if cfg.simulation.run_test
        @info "═══ Phase 2: Simulation ═══"
        run_dir = run_simulation(cfg; base_dir)
    else
        @info "Simulation disabled — skipping."
    end

    # ---- Phase 3: Post-processing ----
    results = Dict{Symbol,Any}()
    if cfg.postprocess.enabled
        t0 = time()
        @info "═══ Phase 3: Post-processing ═══"
        results = postprocess(cfg; run_dir, base_dir)
        for (k, v) in results
            n = if v isa AbstractVector
                length(v)
            elseif v isa DiagnosticsData
                length(v.adjust)
            else
                1
            end
            @info "  Loaded: $k ($n records)"
        end
        @info "Phase 3 complete ($(Nbody6Setup._format_elapsed(time() - t0)))"
    else
        @info "Post-processing disabled — skipping."
    end

    # ---- Phase 4: Visualisation ----
    if cfg.visualization.enabled && !isempty(results)
        t0 = time()
        @info "═══ Phase 4: Visualisation ═══"
        generate_plots(results, cfg; run_dir)
        @info "Phase 4 complete ($(Nbody6Setup._format_elapsed(time() - t0)))"
    elseif cfg.visualization.enabled
        @info "Visualisation enabled but no data to plot."
    end

    total = Nbody6Setup._format_elapsed(time() - t_pipeline)
    @info "Done. Total pipeline time: $total"
end

main()
