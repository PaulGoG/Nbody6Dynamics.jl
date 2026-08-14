# Self-activating (§3): the docs environment consumes the package by path,
# so the build always runs against the local source.
using Pkg
Pkg.activate(@__DIR__; io = devnull)
Pkg.develop(; path = joinpath(@__DIR__, ".."), io = devnull)
Pkg.instantiate(; io = devnull)

using Documenter
using Nbody6Dynamics

makedocs(;
    modules = [Nbody6Dynamics],
    sitename = "Nbody6Dynamics.jl",
    authors = "Paul-Adrian Gogîță",
    repo = Documenter.Remotes.GitHub("PaulGoG", "Nbody6Dynamics.jl"),
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true", assets = String[]),
    pages = [
        "Home" => "index.md",
        "Manual" => "manual.md",
        "Input Files" => "input_files.md",
        "Cluster Mergers" => "multi_cluster_mergers.md",
        "API Reference" => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references],
)
