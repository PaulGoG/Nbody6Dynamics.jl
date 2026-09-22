# =============================================================================
# Publication figures and animations of Nbody6Dynamics.jl
# =============================================================================
# Loaded automatically once a Makie backend is in the session
# (`using CairoMakie`). The core package declares no plotting dependency, so a
# headless host installs and runs the pipeline — initial conditions, engine
# build, integration, readers, diagnostics, benchmarks — without the
# Cairo/Pango/GLib/HarfBuzz stack and without precompiling it.
#
# MathTeXEngine is the second trigger because the theme resolves the Computer
# Modern faces through it; Makie depends on it, so the pair is satisfied
# together by `using CairoMakie` alone.
#
# Every public figure routine is defined as a method of the core package's
# function (`function Nbody6Dynamics.plot_energy(...)`), whose fallback in
# src/plotting_api.jl reports the missing backend when this extension is not
# loaded. Helpers stay private to this module.

module Nbody6DynamicsMakieExt

using CairoMakie
using LaTeXStrings
using MathTeXEngine: texfont
using Printf: @sprintf

using Nbody6Dynamics
using Nbody6Dynamics:
    _AU_IN_PC,
    _G_PC_KMS2_MSUN,
    _KROUPA_ALPHAS,
    _KROUPA_BREAKS,
    _MSR_THRESHOLD,
    _axis_short,
    _backup_existing,
    _cluster_com_trajectories,
    _count_spatial_clusters,
    _mean,
    _run_series,
    _series_label,
    _sweep_done_points
# Public (unexported) names of the package the figure layer calls unqualified.
using Nbody6Dynamics:
    binary_hardness,
    binary_scales,
    class_counts,
    classes_present,
    nparticles,
    profile_name,
    rbar,
    read_sweep_index,
    semi_major_axis_pc,
    stellar_class_index,
    time_myr,
    time_nb,
    tscale,
    vstar,
    zmbar

# Shared theme, figure sizing, axis and annotation helpers
include("plotting/common.jl")

# Figure families
include("plotting/snapshots.jl")
include("plotting/lagrangian.jl")
include("plotting/energy.jl")
include("plotting/hr.jl")
include("plotting/escapers.jl")
include("plotting/sse.jl")
include("plotting/animation.jl")
include("plotting/merger.jl")
include("plotting/binaries.jl")
include("plotting/sweep.jl")
include("plotting/ensemble.jl")
include("plotting/remnant.jl")
include("plotting/control.jl")
include("plotting/telemetry.jl")
include("plotting/merger_ic.jl")

# The dispatcher every caller reaches
include("plotting/dispatch.jl")

end # module
