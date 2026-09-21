# Activates and instantiates the documentation environment without console
# output. The package is developed by a path relative to this directory, so
# the environment resolves against the local source on any machine. make.jl includes
# this file; on a new machine `julia docs/activate.jl` bootstraps the
# environment on its own.
using Pkg
Pkg.activate(@__DIR__; io = devnull)
cd(@__DIR__) do
    Pkg.develop(; path = "..", io = devnull)
end
Pkg.instantiate(; io = devnull)
if abspath(PROGRAM_FILE) == @__FILE__
    println("Environment ready: ", Base.active_project(), " (Julia ", VERSION, ")")
end
