# =============================================================================
# The benchmark cells shared by the scaling scripts. The dual cell is two
# King clusters of N_total/2 stars each on an eccentric Kepler orbit, natural
# Kroupa IMF, Jacobi truncation; the single cell is one King cluster of
# N_total stars at rest with the same structural parameters. Both integrate
# for `tcrit` N-body time units.
# =============================================================================

"""
    merger_toml(N_total, tcrit) -> String

Merger TOML of the benchmark case for `N_total` stars in two equal King
clusters, integrated for `tcrit` N-body time units with adjustments every
0.25 and outputs every 0.5. The virialisation limit is raised to 2×10⁶ so
that the benchmark reaches the engine's `b1m` capacity (NMAX = 1 572 864):
a benchmark is a deliberate cost.
"""
function merger_toml(N_total::Integer, tcrit::Real)
    n = N_total ÷ 2
    """
    [merger]
    n_clusters = 2
    orbit_mode = "kepler"
    seed = 11
    virial_max_n = 2000000

    [merger.cluster1]
    model = "king"
    N = $n
    W0 = 6.0
    rbar = 2.0
    imf = "kroupa"

    [merger.cluster2]
    model = "king"
    N = $n
    W0 = 6.0
    rbar = 2.0
    imf = "kroupa"

    [merger.orbit]
    apocentre = 12.0
    eccentricity = 0.6

    [merger.output]
    truncate_jacobi = true
    tcrit = $tcrit
    dtadj = 0.25
    deltat = 0.5
    """
end

"""
    single_toml(N_total, tcrit) -> String

Merger TOML of the single-cluster cell: one King cluster of `N_total` stars
at rest (`n_clusters = 1`, explicit mode, no truncation) with the structural
parameters of [`merger_toml`](@ref), integrated for `tcrit` N-body time
units with adjustments every 0.25 and outputs every 0.5.
"""
function single_toml(N_total::Integer, tcrit::Real)
    """
    [merger]
    n_clusters = 1
    orbit_mode = "explicit"
    seed = 11
    virial_max_n = 2000000

    [merger.cluster1]
    model = "king"
    N = $N_total
    W0 = 6.0
    rbar = 2.0
    imf = "kroupa"
    position = [0.0, 0.0, 0.0]
    velocity = [0.0, 0.0, 0.0]

    [merger.output]
    truncate_jacobi = false
    tcrit = $tcrit
    dtadj = 0.25
    deltat = 0.5
    """
end

"""
    cell_toml(cell, N_total, tcrit) -> String

The merger TOML of benchmark cell `cell`: `"dual"` ([`merger_toml`](@ref))
or `"single"` ([`single_toml`](@ref)).
"""
function cell_toml(cell::AbstractString, N_total::Integer, tcrit::Real)
    cell == "dual" && return merger_toml(N_total, tcrit)
    cell == "single" && return single_toml(N_total, tcrit)
    throw(ArgumentError("unknown benchmark cell \"$cell\"; choose \"single\" or \"dual\""))
end
