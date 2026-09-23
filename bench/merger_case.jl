# =============================================================================
# The benchmark case shared by the scaling scripts: two King clusters of
# N_total/2 stars each on an eccentric Kepler orbit, natural Kroupa IMF,
# Jacobi truncation, integrated for `tcrit` N-body time units.
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
