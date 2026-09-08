# =============================================================================
# Primordial binary populations for the merger initial conditions
# =============================================================================
# Stars are sampled from the IMF, a fraction of them are paired, the
# centres of mass of all systems are placed by the density sampler and
# virialised, and each pair is expanded into two bodies on a Keplerian
# orbit at the end. The engine reads primordial pairs as the first 2·NBIN0
# bodies of dat.10 (bodies 2i−1, 2i), so every cluster's pairs are written
# first, cluster by cluster, followed by every cluster's singles.

"""
    BinarySpec(; fraction = 0.0, pairing = "random", period = "kroupa1995",
                 a_min = 0.01, a_max = 100.0, q_min = 0.1, eccentricity = "thermal")

Primordial binary population of one cluster (`[merger.clusterN.binaries]`).

# Fields
- `fraction`: binary fraction by systems, `N_b / (N_s + N_b)`; `0 ≤ fraction < 1`
- `pairing`: `"random"` (both components drawn from the IMF; the more massive
  becomes the primary) or `"uniform_q"` (the secondary's mass is `q m₁` with
  `q` uniform in `[q_min, 1]`, replacing the drawn mass of the partner)
- `period`: `"kroupa1995"` (the Kroupa 1995 birth period distribution,
  `f(log P) ∝ (log P − 1) / (45 + (log P − 1)²)` for `1 ≤ log P/d ≤ 8.43`,
  converted to a semi-major axis with Kepler's third law) or `"loguniform"`
  (semi-major axis log-uniform between `a_min` and `a_max`)
- `a_min`, `a_max`: semi-major axis bounds [AU] for `"loguniform"`; `0 < a_min < a_max`
- `q_min`: lower mass-ratio bound for `"uniform_q"`; `0 < q_min < 1`
- `eccentricity`: `"thermal"` (`f(e) = 2e`) or `"circular"`
"""
Base.@kwdef struct BinarySpec
    fraction::Float64 = 0.0
    pairing::String = "random"
    period::String = "kroupa1995"
    a_min::Float64 = 0.01
    a_max::Float64 = 100.0
    q_min::Float64 = 0.1
    eccentricity::String = "thermal"
end

const _AU_PC = 1 / 206264.806
const _DAY_YR = 1 / 365.25
# Kroupa (1995b) birth period distribution parameters and range (log10 P in days)
const _KROUPA95_ETA = 2.5
const _KROUPA95_DELTA = 45.0
const _KROUPA95_LOGP_MIN = 1.0
const _KROUPA95_LOGP_MAX = 8.43

"""Fail-fast bounds of a [`BinarySpec`](@ref); `key` names the owning cluster table."""
function _validate_binaries(b::BinarySpec, key::String)
    (0 ≤ b.fraction < 1) ||
        error("config: $key.binaries.fraction must satisfy 0 ≤ fraction < 1; got $(b.fraction)")
    b.pairing in ("random", "uniform_q") || error(
        "config: $key.binaries.pairing must be \"random\" or \"uniform_q\"; got \"$(b.pairing)\"",
    )
    b.period in ("kroupa1995", "loguniform") || error(
        "config: $key.binaries.period must be \"kroupa1995\" or \"loguniform\"; got \"$(b.period)\"",
    )
    b.eccentricity in ("thermal", "circular") || error(
        "config: $key.binaries.eccentricity must be \"thermal\" or \"circular\"; got \"$(b.eccentricity)\"",
    )
    (0 < b.a_min < b.a_max) || error(
        "config: $key.binaries must satisfy 0 < a_min < a_max [AU]; got $(b.a_min), $(b.a_max)",
    )
    (0 < b.q_min < 1) ||
        error("config: $key.binaries.q_min must satisfy 0 < q_min < 1; got $(b.q_min)")
    return nothing
end

"""
    _sample_log_period_kroupa1995(rng) -> Float64

`log10 P` [days] from the Kroupa (1995) birth period distribution by
inverse transform: with `x = log P − 1`, `F(x) = (η/2) ln((δ + x²)/δ)`.
"""
function _sample_log_period_kroupa1995(rng::AbstractRNG)::Float64
    x_max = _KROUPA95_LOGP_MAX - _KROUPA95_LOGP_MIN
    f_max = 0.5 * _KROUPA95_ETA * log((_KROUPA95_DELTA + x_max^2) / _KROUPA95_DELTA)
    u = rand(rng) * f_max
    x = sqrt(_KROUPA95_DELTA * (exp(2u / _KROUPA95_ETA) - 1))
    return _KROUPA95_LOGP_MIN + x
end

"""
    sample_binaries(spec::BinarySpec, masses, rng) -> (; masses, primary, secondary, m1, m2, a_pc, e)

Pair `N_b = round(fraction N / (1 + fraction))` of the `N` stars. Returns the
(possibly modified, for `"uniform_q"`) mass vector, the star indices of the
primaries and secondaries, the component masses [M☉] with `m1 ≥ m2`, the
semi-major axes [pc] from the period or axis distribution, and the
eccentricities. Every star belongs to at most one pair.
"""
function sample_binaries(spec::BinarySpec, masses::AbstractVector{Float64}, rng::AbstractRNG)
    N = length(masses)
    m = collect(Float64, masses)
    N_b = round(Int, spec.fraction * N / (1 + spec.fraction))
    2N_b ≤ N || (N_b = N ÷ 2)
    perm = randperm(rng, N)
    primary = perm[1:N_b]
    secondary = perm[(N_b + 1):(2N_b)]
    for k in 1:N_b
        if spec.pairing == "uniform_q"
            q = spec.q_min + (1 - spec.q_min) * rand(rng)
            m[secondary[k]] = q * m[primary[k]]
        elseif m[secondary[k]] > m[primary[k]]
            primary[k], secondary[k] = secondary[k], primary[k]
        end
    end
    m1 = m[primary]
    m2 = m[secondary]
    a_pc = Vector{Float64}(undef, N_b)
    e = Vector{Float64}(undef, N_b)
    for k in 1:N_b
        a_au = if spec.period == "kroupa1995"
            p_yr = exp10(_sample_log_period_kroupa1995(rng)) * _DAY_YR
            cbrt((m1[k] + m2[k]) * p_yr^2)        # Kepler III in AU, yr, M☉
        else
            spec.a_min * (spec.a_max / spec.a_min)^rand(rng)
        end
        a_pc[k] = a_au * _AU_PC
        e[k] = spec.eccentricity == "thermal" ? sqrt(rand(rng)) : 0.0
    end
    return (
        masses = m,
        primary = primary,
        secondary = secondary,
        m1 = m1,
        m2 = m2,
        a_pc = a_pc,
        e = e,
    )
end

"""
    _random_rotation(rng) -> Matrix{Float64}

Uniformly distributed rotation matrix (from a random unit quaternion).
"""
function _random_rotation(rng::AbstractRNG)::Matrix{Float64}
    u1, u2, u3 = rand(rng), rand(rng), rand(rng)
    q0 = sqrt(1 - u1) * sin(2π * u2)
    q1 = sqrt(1 - u1) * cos(2π * u2)
    q2 = sqrt(u1) * sin(2π * u3)
    q3 = sqrt(u1) * cos(2π * u3)
    return [
        1-2(q2^2+q3^2) 2(q1*q2-q0*q3) 2(q1*q3+q0*q2)
        2(q1*q2+q0*q3) 1-2(q1^2+q3^2) 2(q2*q3-q0*q1)
        2(q1*q3-q0*q2) 2(q2*q3+q0*q1) 1-2(q1^2+q2^2)
    ]
end

"""
    _kepler_relative_orbit(M, a, e, rng) -> (r_vec, v_vec)

Relative position and velocity of a two-body orbit of total mass `M`,
semi-major axis `a`, and eccentricity `e` (`G = 1`) at a uniformly random
mean anomaly (Kepler's equation by Newton iteration), in a uniformly random
orientation.
"""
function _kepler_relative_orbit(M::Real, a::Real, e::Real, rng::AbstractRNG)
    mean_anomaly = 2π * rand(rng)
    E = mean_anomaly
    for _ in 1:50
        δ = (E - e * sin(E) - mean_anomaly) / (1 - e * cos(E))
        E -= δ
        abs(δ) < 1e-13 && break
    end
    r = a * (1 - e * cos(E))
    cosν = (cos(E) - e) / (1 - e * cos(E))
    sinν = sqrt(1 - e^2) * sin(E) / (1 - e * cos(E))
    vscale = sqrt(M / (a * (1 - e^2)))
    r_plane = [r * cosν, r * sinν, 0.0]
    v_plane = [-vscale * sinν, vscale * (e + cosν), 0.0]
    R = _random_rotation(rng)
    return R * r_plane, R * v_plane
end

"""
    expand_binaries(pos, vel, mass, system_ranges, cluster_binaries, kept; rng)
        -> (; pos, vel, mass, cluster_blocks, n_pairs, hard_fraction)

Turn the combined *system* arrays (code units, `G = 1`) into *body* arrays:
every pair becomes two bodies on a Keplerian relative orbit about the
system's centre of mass, written first, cluster by cluster, followed by the
singles of every cluster (the engine's primordial-pair convention). For
cluster `i`, `system_ranges[i]` are its systems in the combined arrays,
`kept[i]` their indices in that cluster's sampled system set, and
`cluster_binaries[i]` the sampled pair table (`system_binary[j]` = pair
index of sampled system `j`, 0 for a single, plus `m1`, `m2`, `a_pc`, `e`).
`cluster_blocks[i]` are the contiguous body-index blocks of cluster `i`
(its pairs block and its singles block), `n_pairs[i]` its number of pairs,
and `hard_fraction[i]` the fraction of its pairs with binding energy above
`⟨m⟩ σ²` of the cluster's systems.
"""
function expand_binaries(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    system_ranges::AbstractVector{<:AbstractVector{Int}},
    cluster_binaries::AbstractVector,
    kept::AbstractVector{<:AbstractVector{Int}};
    rng::AbstractRNG = Random.default_rng(),
)
    n_cl = length(system_ranges)
    pairs = [Int[] for _ in 1:n_cl]     # combined-array system index of each surviving pair
    pair_ids = [Int[] for _ in 1:n_cl]  # its index in the cluster's pair table
    singles = [Int[] for _ in 1:n_cl]
    for i in 1:n_cl
        cb = cluster_binaries[i]
        for (k, sys) in enumerate(system_ranges[i])
            p = cb.system_binary[kept[i][k]]
            if p > 0
                push!(pairs[i], sys)
                push!(pair_ids[i], p)
            else
                push!(singles[i], sys)
            end
        end
    end
    n_pairs = length.(pairs)
    nbin0 = sum(n_pairs)
    N_bodies = 2nbin0 + sum(length, singles)
    pos_b = zeros(Float64, 3, N_bodies)
    vel_b = zeros(Float64, 3, N_bodies)
    mass_b = zeros(Float64, N_bodies)
    cluster_blocks = Vector{Vector{UnitRange{Int}}}(undef, n_cl)
    hard_fraction = fill(NaN, n_cl)
    b = 0
    for i in 1:n_cl
        cb = cluster_binaries[i]
        # 1-D velocity dispersion of the cluster's systems for the hard/soft boundary
        idx_sys = system_ranges[i]
        M_sys = sum(mass[idx_sys])
        vc = vec(sum(vel[:, idx_sys] .* mass[idx_sys]'; dims = 2)) ./ M_sys
        σ2 = sum(mass[s] * sum((vel[:, s] .- vc) .^ 2) for s in idx_sys) / (3 * M_sys)
        m_mean = M_sys / length(idx_sys)
        n_hard = 0
        first_pair_body = b + 1
        for (sys, p) in zip(pairs[i], pair_ids[i])
            m1, m2, a, e = cb.m1[p], cb.m2[p], cb.a_pc[p], cb.e[p]
            M = m1 + m2
            r_rel, v_rel = _kepler_relative_orbit(M, a, e, rng)
            n_hard += (m1 * m2 / (2a) > m_mean * σ2) ? 1 : 0
            for (body, m, f) in ((b + 1, m1, m2 / M), (b + 2, m2, -m1 / M))
                mass_b[body] = m
                for k in 1:3
                    pos_b[k, body] = pos[k, sys] + f * r_rel[k]
                    vel_b[k, body] = vel[k, sys] + f * v_rel[k]
                end
            end
            b += 2
        end
        hard_fraction[i] = isempty(pairs[i]) ? NaN : n_hard / length(pairs[i])
        cluster_blocks[i] = [first_pair_body:b]
    end
    s = 2nbin0
    for i in 1:n_cl
        first_single = s + 1
        for sys in singles[i]
            s += 1
            mass_b[s] = mass[sys]
            pos_b[:, s] .= pos[:, sys]
            vel_b[:, s] .= vel[:, sys]
        end
        push!(cluster_blocks[i], first_single:s)
    end
    return (
        pos = pos_b,
        vel = vel_b,
        mass = mass_b,
        cluster_blocks = cluster_blocks,
        n_pairs = n_pairs,
        hard_fraction = hard_fraction,
    )
end

"""Body-index vector of a cluster from its contiguous blocks."""
_members_from_blocks(blocks::AbstractVector{<:UnitRange{Int}}) =
    reduce(vcat, (collect(r) for r in blocks); init = Int[])
