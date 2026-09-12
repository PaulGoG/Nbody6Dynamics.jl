# # Walkthrough: a two-cluster merger from configuration to initial conditions
#
# This page is generated from `docs/src/walkthrough.jl` with Literate.jl and
# executed when the documentation is built. It covers what a merger study
# needs before the engine runs: the configuration, the generated initial
# conditions, the derived integration parameters, and the initial-condition
# figures. The engine itself is not executed here; the
# [user manual](@ref "Nbody6Dynamics.jl — User Manual") covers `run_pipeline`
# and the post-processing of a finished run.

using CairoMakie #  the figure routines are a package extension; loading a
#                   Makie backend brings them into the session
using Nbody6Dynamics

# ## The merger configuration
#
# Two Plummer clusters of 300 stars with a Kroupa mass function, released on
# an eccentric Kepler orbit. Output intervals given in Myr are converted with
# the realised N-body time unit at generation.

dir = mktempdir()
toml = joinpath(dir, "merger.toml")
write(
    toml,
    """
    [merger]
    seed = 7
    n_clusters = 2
    orbit_mode = "kepler"

    [merger.cluster1]
    model = "plummer"
    N = 300
    rbar = 1.0
    imf = "kroupa"

    [merger.cluster2]
    model = "plummer"
    N = 300
    rbar = 1.0
    imf = "kroupa"

    [merger.orbit]
    apocentre = 6.0
    eccentricity = 0.5

    [merger.output]
    format = "nbody"
    truncate_jacobi = true
    output_dir = "$(dir)"
    tcrit_myr = 10.0
    dtadj_myr = 0.5
    deltat_myr = 1.0
    """,
)
cfg = load_merger_config(toml)
cfg.clusters[1]

# ## Generating the initial conditions
#
# `generate_merger_ic` samples each cluster, virialises it, places the pair on
# the orbit, truncates the members at the Jacobi radius, converts to N-body
# units and writes `dat.10`, `merger.inp`, `merger_ic.toml` and
# `merger_summary.txt` into the output directory.

result = generate_merger_ic(cfg)
(N = result.N_total, M_total = result.M_total, rbar = result.rbar, zmbar = result.zmbar)

# Cluster membership is recorded as blocks of body indices, which every
# per-cluster diagnostic accepts:

length.(result.cluster_ranges)

# ## The engine's input file
#
# The integration parameters are derived from the member clusters (see the
# [Input File Reference](@ref)). The output intervals are dyadic rationals with
# an exact decimal expansion, because the engine counts their decimal digits
# with a loop that never terminates on other values.

print(read(joinpath(dir, "merger.inp"), String))

# ## Initial-condition figures
#
# The IC plot suite draws the three projections, an overview, the velocity
# field by cluster, the mass function against the Kroupa reference slopes and
# the per-cluster density profiles. Here the figures are written next to this
# page.

vis = VisualizationConfig(;
    enabled = true,
    format = "png",
    dpi = 150,
    column = "single",
    output_dir = pwd(),
)
plot_merger_ic(result, vis)
filter(f -> startswith(f, "merger_ic"), readdir(pwd()))

# ![Initial conditions: overview of the two clusters and their orbit](merger_ic_overview.png)
#
# ![Initial conditions: sampled mass function against the Kroupa slopes](merger_ic_imf.png)

# ## From here
#
# `run_pipeline(load_config("config.toml"))` with `[merger] enabled = true` and
# `config_file` pointing at this TOML runs the engine on the generated `dat.10`,
# post-processes every output file and draws the figure suite;
# `postprocess_external(dir)` does the same for output produced elsewhere.
