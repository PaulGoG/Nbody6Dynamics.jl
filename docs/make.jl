using Documenter
using Nbody6Dynamics

makedocs(;
    modules  = [Nbody6Dynamics],
    sitename = "Nbody6Dynamics.jl",
    authors  = "paulgog",
    remotes  = nothing,
    format   = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets     = String[],
    ),
    pages = [
        "Home"              => "index.md",
        "Manual"            => "manual.md",
        "Input Files"       => "input_files.md",
        "Cluster Mergers"   => "multi_cluster_mergers.md",
        "API Reference"     => "api.md",
    ],
    warnonly = [:missing_docs],
)
