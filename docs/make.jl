# Self-activating (§3): the docs environment consumes the package by path,
# so the build always runs against the local source.
using Pkg
Pkg.activate(@__DIR__; io = devnull)
# Develop with a path relative to this directory so the tracked Manifest
# records `path = ".."` and stays portable across machines.
cd(@__DIR__) do
    Pkg.develop(; path = "..", io = devnull)
end
Pkg.instantiate(; io = devnull)

using Documenter
using DocumenterCitations
using Literate
using Nbody6Dynamics

# The walkthrough is a Literate script; its markdown is generated at build
# time (and not tracked) with Documenter @example blocks, which execute.
Literate.markdown(
    joinpath(@__DIR__, "src", "walkthrough.jl"),
    joinpath(@__DIR__, "src");
    documenter = true,
    execute = false,
    credit = false,
)

bib = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"); style = :authoryear)

makedocs(;
    modules = [Nbody6Dynamics],
    sitename = "Nbody6Dynamics.jl",
    authors = "Paul-Adrian Gogîță",
    repo = Documenter.Remotes.GitHub("PaulGoG", "Nbody6Dynamics.jl"),
    plugins = [bib],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = String[],
        # The API reference lists every public docstring on one page; exempt
        # it from Documenter's per-page size threshold.
        size_threshold_ignore = ["api.md"],
    ),
    pages = [
        "Home" => "index.md",
        "Walkthrough" => "walkthrough.md",
        "Manual" => "manual.md",
        "Input Files" => "input_files.md",
        "Cluster Mergers" => "multi_cluster_mergers.md",
        "API Reference" => "api.md",
        "References" => "references.md",
    ],
    warnonly = [:missing_docs],
)

# Publishes from CI only (a GITHUB_* environment); a local build stops here.
deploydocs(; repo = "github.com/PaulGoG/Nbody6Dynamics.jl", devbranch = "main", push_preview = false)
