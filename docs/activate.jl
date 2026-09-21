# Activates and instantiates the documentation environment without console
# output. The environment takes the package by relative path (`[sources]`),
# so it always runs against the local source. make.jl includes
# this file; on a new machine `julia docs/activate.jl` bootstraps the
# environment on its own.
using Pkg
Pkg.activate(@__DIR__; io = devnull)
Pkg.instantiate(; io = devnull)
if abspath(PROGRAM_FILE) == @__FILE__
    println("Environment ready: ", Base.active_project(), " (Julia ", VERSION, ")")
end
