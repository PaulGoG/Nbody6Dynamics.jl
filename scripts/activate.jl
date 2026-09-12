# Activates and instantiates the entry-script environment without console
# output. It is the package plus a Makie backend: the figure routines are a
# package extension, so the pipeline scripts — which all end in figures — name
# CairoMakie here rather than in the package itself. The package is developed
# by a path relative to this directory so that the tracked Manifest stays
# portable across machines. The scripts under scripts/ include this file; on a
# new machine `julia scripts/activate.jl` bootstraps the environment on its
# own (expect a few minutes the first time: the plotting stack precompiles).
using Pkg
Pkg.activate(@__DIR__; io = devnull)
cd(@__DIR__) do
    Pkg.develop(; path = "..", io = devnull)
end
Pkg.instantiate(; io = devnull)
if abspath(PROGRAM_FILE) == @__FILE__
    println("Environment ready: ", Base.active_project(), " (Julia ", VERSION, ")")
end
