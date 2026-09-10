# Activates and instantiates the package environment without console output.
# The entry scripts under scripts/ include this file; on a new machine
# `julia activate.jl` performs the one-time dependency resolution and
# precompilation of the root environment.
using Pkg
Pkg.activate(@__DIR__; io = devnull)
Pkg.instantiate(; io = devnull)
if abspath(PROGRAM_FILE) == @__FILE__
    println("Environment ready: ", Base.active_project(), " (Julia ", VERSION, ")")
end
