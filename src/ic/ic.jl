# =============================================================================
# Initial Conditions submodule
# =============================================================================
# Multi-cluster merger IC generation for Nbody6++GPU.
# Produces dat.10 (external particle input) + matching .inp files.

using Random
using SpecialFunctions: erf

include("binaries.jl")
include("config.jl")
include("models.jl")
include("imf.jl")
include("orbits.jl")
include("output.jl")
include("plotting.jl")

# ---------------------------------------------------------------------------
# Top-level merger pipeline
# ---------------------------------------------------------------------------

"""
    run_merger_pipeline(config_path::AbstractString; rng=nothing,
                        output_dir="", make_plots=true) -> MergerICResult

One-call entry point for multi-cluster merger IC generation.

1. Loads the merger TOML configuration from `config_path`
2. Generates particle ICs (Plummer/King + Kroupa IMF + Kepler orbit)
3. Writes `dat.10` + `merger.inp` + `merger_summary.txt`
4. Generates diagnostic plots (projections, velocity field, IMF, density)

Returns a [`MergerICResult`](@ref) for programmatic inspection.

# Example
```julia
using Nbody6Dynamics
result = run_merger_pipeline("input_files/merger_equal_mass.toml")
```
"""
function run_merger_pipeline(
    config_path::AbstractString;
    rng::Union{AbstractRNG,Nothing} = nothing,
    output_dir::AbstractString = "",
    make_plots::Bool = true,
)
    @info "═══ Merger IC Pipeline ═══"
    @info "Loading config: $config_path"
    cfg = load_merger_config(config_path)

    # Resolve output directory
    out = if !isempty(output_dir)
        output_dir
    elseif cfg.output.output_dir != "."
        cfg.output.output_dir
    else
        # Auto-generate a unique run directory under the project's runs/.
        # The suffix RNG is cosmetic (directory uniqueness), so the global
        # RNG is fine here; the physics RNG is resolved in generate_merger_ic.
        joinpath(_PROJECT_ROOT, "runs", generate_run_id("merger"))
    end

    @info "Output directory: $out"
    result = generate_merger_ic(cfg; rng = rng, output_dir = out)

    if make_plots
        @info "Generating diagnostic plots..."
        # Project visualization defaults (PDF vector, single-column preset)
        vis = VisualizationConfig(; output_dir = joinpath(out, "plots"))
        plot_merger_ic(result, vis)
    end

    @info "═══ Merger IC Pipeline Complete ═══"
    @info "  dat.10:     $(joinpath(out, "dat.10"))"
    @info "  merger.inp: $(joinpath(out, "merger.inp"))"
    make_plots && @info "  plots:      $(joinpath(out, "plots/"))"

    return result
end
