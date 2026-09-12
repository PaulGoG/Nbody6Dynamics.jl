# =============================================================================
# Figure interface — implemented by the Makie extension
# =============================================================================
# The figure routines live in ext/Nbody6DynamicsMakieExt.jl and reach the
# session only when a Makie backend is loaded. Here the core package declares
# the public functions, so that the names exist, dispatch works from
# orchestration code compiled without a backend, and a call made without one
# reports the remedy instead of an `UndefVarError`.

"""
    PlottingUnavailable(entry_point, hint = "")

No Makie backend is loaded, so the figure routine `entry_point` has no
methods. Thrown by every public plotting entry point of the package and by
the pipeline phases that would produce figures; `hint` carries the
context-specific remedy.
"""
struct PlottingUnavailable <: Exception
    entry_point::Symbol
    hint::String
end

PlottingUnavailable(entry_point::Symbol) = PlottingUnavailable(entry_point, "")

function Base.showerror(io::IO, e::PlottingUnavailable)
    print(
        io,
        "PlottingUnavailable: `",
        e.entry_point,
        "` needs a Makie backend. Load one first:\n",
        "    using CairoMakie\n",
        "which triggers the package extension Nbody6DynamicsMakieExt carrying ",
        "every figure routine.",
    )
    isempty(e.hint) || print(io, "\n", e.hint)
    return nothing
end

"""
    plotting_available() -> Bool

Whether the Makie extension is loaded, i.e. whether the figure routines have
methods. `false` on a headless installation that never loaded a backend.

# Example

```julia
plotting_available() || @info "Running without figures"
```
"""
plotting_available() = Base.get_extension(@__MODULE__, :Nbody6DynamicsMakieExt) !== nothing

"""
    _require_plotting(entry_point, hint = "")

Fail before doing any work when a phase that must produce figures has no
backend. Used by the orchestration entry points, where the alternative is
discovering the missing backend after an integration has run for hours.
"""
function _require_plotting(entry_point::Symbol, hint::AbstractString = "")
    plotting_available() && return nothing
    throw(PlottingUnavailable(entry_point, String(hint)))
end

"""
The public figure routines. Each is declared here without methods and
implemented by the Makie extension; the fallback below answers a call made
without a backend.
"""
const _PLOTTING_ENTRY_POINTS = (
    # theme
    :publication_theme,
    :set_publication_theme!,
    # dispatcher
    :generate_plots,
    # snapshots, structure, initial conditions
    :plot_snapshot,
    :plot_snapshot_evolution,
    :plot_cluster_structure,
    :plot_cluster_separation,
    :plot_cluster_virial,
    :plot_density_profiles,
    :plot_velocity_dispersion,
    :plot_merger_ic,
    # integration diagnostics
    :plot_lagrangian,
    :plot_energy,
    :plot_particle_count,
    :plot_escapers,
    :plot_escape_anisotropy,
    # stellar populations
    :plot_hr,
    :plot_hr_evolution,
    :plot_mass_segregation,
    :plot_evolutionary_clock,
    :plot_core_mass,
    # binaries
    :plot_binary_population,
    :plot_binary_orbital_elements,
    :plot_binary_period_distribution,
    # remnant of a merger
    :plot_remnant_rotation,
    :plot_rotation_profile,
    :plot_remnant_structure,
    :plot_mass_segregation_evolution,
    :remnant_figures,
    # sweeps, ensembles, controls
    :plot_sweep_lagrangian,
    :plot_sweep_energy,
    :plot_sweep_ensemble,
    :plot_control_comparison,
    :sweep_figures,
    # telemetry
    :plot_telemetry,
    # animations
    :animate_cluster,
    :animate_hr,
    :animate_lagrangian,
)

"""
    _no_plotting_backend(entry_point, args)

Answer a call to a figure routine that found no method: a `MethodError` when
the extension is loaded (the arguments are simply wrong) and a
[`PlottingUnavailable`](@ref) when it is not.
"""
function _no_plotting_backend(entry_point::Symbol, args::Tuple)
    plotting_available() && throw(MethodError(getfield(@__MODULE__, entry_point), args))
    throw(PlottingUnavailable(entry_point))
end

for entry_point in _PLOTTING_ENTRY_POINTS
    @eval begin
        function $entry_point end
        $entry_point(args...; kwargs...) = _no_plotting_backend($(QuoteNode(entry_point)), args)
    end
end
