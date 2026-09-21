# Activates and instantiates the benchmark environment without console
# output. The environment takes the package by relative path (`[sources]`),
# so it always runs against the local source. The benchmark
# scripts include this file; on a new machine `julia bench/activate.jl`
# bootstraps the environment on its own.
using Pkg
Pkg.activate(@__DIR__; io = devnull)
Pkg.instantiate(; io = devnull)
if abspath(PROGRAM_FILE) == @__FILE__
    println("Environment ready: ", Base.active_project(), " (Julia ", VERSION, ")")
end
