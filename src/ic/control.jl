# =============================================================================
# Control configurations: the isolated single-cluster equivalent of a merger
# =============================================================================
#
# A control run holds everything but the merger geometry fixed: the same
# total number of stars, the structural model, IMF, binary population,
# stellar-evolution, integration and tidal settings of the merger
# configuration, in one cluster at rest. Its half-mass radius is the
# N-weighted mean of the progenitors' half-mass radii (equal to the
# mass-weighted mean in expectation when the clusters share an IMF), so the
# control has the progenitors' internal scale rather than the remnant's.
# The derivation works on the TOML tables, so any merger configuration file
# yields its control without sampling. Time intervals given in Myr
# (`tcrit_myr`, `dtadj_myr`, `deltat_myr`, `dtplot_myr`) carry over
# verbatim and are converted by the generator with each run's own time
# unit, so merger and control then span the same physical time exactly.
# Intervals given in N-body units are scaled by the ratio of the estimated
# time units, (RBAR_est / r_h,control)^{3/2}, with RBAR_est the half-mass
# radius of the initial configuration estimated from the separations and
# the members' radii — an approximation, stated in the docs.

"""
    _control_time_factor(m::AbstractDict, Ns, rbars, rbar_control) -> Float64

Ratio of the estimated N-body time units of the merger configuration and
of its control, `(RBAR_est / rbar_control)^{3/2}`: the total masses are
equal, so only the length units differ. `RBAR_est` adds the `N`-weighted
mean member radius to the half-mass radius of the centres of mass — for a
Kepler orbit `apocentre × N_min / N_total`, for explicit positions their
`N`-weighted RMS distance from the `N`-weighted centre.
"""
function _control_time_factor(
    m::AbstractDict,
    Ns::AbstractVector{Int},
    rbars::AbstractVector{Float64},
    rbar_control::Real,
)
    N_total = sum(Ns)
    rbar_mean = sum(Ns .* rbars) / N_total
    mode = String(get(m, "orbit_mode", "kepler"))
    r_centres = if mode == "kepler"
        d = Float64(get(get(m, "orbit", Dict{String,Any}()), "apocentre", 15.0))
        d * minimum(Ns) / N_total
    else
        n = length(Ns)
        pos = [Float64.(get(m["cluster$i"], "position", [0.0, 0.0, 0.0])) for i in 1:n]
        centre = sum(Ns[i] .* pos[i] for i in 1:n) ./ N_total
        sqrt(sum(Ns[i] * sum((pos[i] .- centre) .^ 2) for i in 1:n) / N_total)
    end
    return ((r_centres + rbar_mean) / rbar_control)^1.5
end

"""
    control_merger_dict(raw::AbstractDict) -> Dict{String,Any}

Derive the control configuration from a parsed merger TOML (`raw` holds the
top-level `"merger"` table): one cluster in explicit orbit mode at the
origin and at rest, `N` the sum over the merger's clusters, `rbar` their
`N`-weighted mean, every other cluster-1 key (model, `W0`, IMF, binaries)
and the `output`, `nbody6`, `stellar`, `tidal` and `seed` entries kept
verbatim; the `orbit` table and the other clusters are dropped. Time
intervals in N-body units (`output.tcrit`, `dtadj`, `deltat`,
`stellar.dtplot`) are multiplied by [`_control_time_factor`](@ref) so the
control covers about the merger's physical span; intervals in Myr are left
as they are and match exactly.
"""
function control_merger_dict(raw::AbstractDict)
    haskey(raw, "merger") || throw(ArgumentError("control_merger_dict: missing [merger] table"))
    m = deepcopy(raw["merger"])
    n = Int(get(m, "n_clusters", 2))
    n ≥ 1 || throw(ArgumentError("control_merger_dict: n_clusters must be ≥ 1, got $n"))
    tables = [get(m, "cluster$i", nothing) for i in 1:n]
    any(isnothing, tables) && throw(
        ArgumentError(
            "control_merger_dict: every [merger.clusterN] table up to n_clusters = $n is required",
        ),
    )
    Ns = [Int(get(t, "N", 50000)) for t in tables]
    rbars = [Float64(get(t, "rbar", 2.0)) for t in tables]
    N_total = sum(Ns)
    control = deepcopy(tables[1])
    control["N"] = N_total
    control["rbar"] = sum(Ns .* rbars) / N_total
    control["position"] = [0.0, 0.0, 0.0]
    control["velocity"] = [0.0, 0.0, 0.0]
    factor = _control_time_factor(m, Ns, rbars, control["rbar"])
    out = get!(m, "output", Dict{String,Any}())
    for key in ("tcrit", "dtadj", "deltat")
        Float64(get(out, key * "_myr", 0.0)) > 0 && continue
        haskey(out, key) && (out[key] = Float64(out[key]) * factor)
    end
    # The plotting interval follows the snapshot interval (it must stay ≥
    # deltat), so its effective NB value is scaled too — also when the base
    # file leaves it at the default.
    st = get!(m, "stellar", Dict{String,Any}())
    if Float64(get(st, "dtplot_myr", 0.0)) == 0 && Float64(get(out, "deltat_myr", 0.0)) == 0
        st["dtplot"] = Float64(get(st, "dtplot", 1.0)) * factor
    end
    for i in 2:n
        delete!(m, "cluster$i")
    end
    delete!(m, "orbit")
    m["cluster1"] = control
    m["n_clusters"] = 1
    m["orbit_mode"] = "explicit"
    result = Dict{String,Any}(k => v for (k, v) in raw if k != "merger")
    result["merger"] = m
    return result
end

"""
    write_control_merger_config(src::AbstractString, dst::AbstractString) -> String

Read the merger TOML `src`, derive its control ([`control_merger_dict`](@ref)),
write it to `dst` (an existing file is backed up, never overwritten),
validate the result through [`load_merger_config`](@ref) and return `dst`.
"""
function write_control_merger_config(src::AbstractString, dst::AbstractString)
    isfile(src) || error("Merger configuration not found: $src")
    d = control_merger_dict(TOML.parsefile(src))
    _backup_existing(dst)
    _atomic_write_toml(dst, d)
    load_merger_config(dst)
    return String(dst)
end
