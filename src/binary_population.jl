# =============================================================================
# Binary-population diagnostics from bev.82 snapshots
# =============================================================================
#
# Hardness follows Heggie (1975): a pair of binding energy
# E_b = G m₁ m₂ / (2a) is hard when E_b > ⟨m⟩ σ², with ⟨m⟩ the mean system
# mass and σ the one-dimensional, mass-weighted velocity dispersion of the
# systems (single stars plus the centres of mass of the pairs). The
# initial-condition generator (`expand_binaries`) applies the same convention,
# so hard fractions at t = 0 are directly comparable with the generated ones.
#
# bev.82 lists KS-regularised pairs only, so every count derived here is a
# lower bound on the bound-pair population: wide pairs beyond the
# regularisation distance RMIN are invisible to these diagnostics.

"""Solar radius in pc (IAU 2015 nominal R☉ = 6.957 × 10⁸ m; 1 pc = 3.0857 × 10¹⁶ m)."""
const _RSUN_PC = 6.957e8 / 3.0857e16

"""Astronomical unit in pc (IAU 2012: 1 au = 1.495978707 × 10¹¹ m)."""
const _AU_IN_PC = 1.495978707e11 / 3.0857e16

"""
    semi_major_axis_pc(r::BinaryRecord) -> Float64

Semi-major axis in pc from the record's `log10(a / R☉)`.
"""
semi_major_axis_pc(r::BinaryRecord) = 10.0^r.log_semi_major_axis_rsun * _RSUN_PC

"""
    binding_energy(r::BinaryRecord) -> Float64

Binding energy ``G m_1 m_2 / (2a)`` of a pair in M☉ (km s⁻¹)².
"""
function binding_energy(r::BinaryRecord)
    return _G_PC_KMS2_MSUN * r.mass1 * r.mass2 / (2 * semi_major_axis_pc(r))
end

"""
    hardness_scale(snap::Snapshot, bev::BinaryEvolutionSnapshot; bound_only = true)
        -> (; m_mean, sigma_kms, n_systems)

Mean system mass ⟨m⟩ [M☉] and one-dimensional, mass-weighted velocity
dispersion σ [km s⁻¹] of the systems in `snap`: single stars plus the centres
of mass of the pairs listed in `bev`, matched to the snapshot by particle
name. A pair with a component missing from the snapshot is not reduced and
its present component counts as a single star. With `bound_only` the scale
is taken over the self-consistently bound systems ([`_bound_members`](@ref)),
which keeps escapers and kicked stellar remnants — a handful of stars at
tens of km s⁻¹ — from dominating the mass-weighted dispersion; `n_systems`
counts the systems used. `⟨m⟩ σ²` is the hard/soft energy scale used by
[`binary_hardness`](@ref).
"""
function hardness_scale(snap::Snapshot, bev::BinaryEvolutionSnapshot; bound_only::Bool = true)
    N = nparticles(snap)
    N > 0 || throw(ArgumentError("hardness_scale: empty snapshot"))
    zm = zmbar(snap.header)
    vs = vstar(snap.header)

    index_of = Dict{Int32,Int}()
    sizehint!(index_of, N)
    for i in 1:N
        index_of[snap.name[i]] = i
    end

    m_sys = Vector{Float64}(undef, N)
    v_sys = Matrix{Float64}(undef, 3, N)
    is_component = falses(N)
    pair_index = Int[]      # first component of every reduced pair (its position stands for the c.m.)
    single_index = Int[]
    n_sys = 0
    for r in bev.records
        i = get(index_of, r.name1, 0)
        j = get(index_of, r.name2, 0)
        (i == 0 || j == 0 || i == j) && continue
        (is_component[i] || is_component[j]) && continue
        mi = Float64(snap.mass[i])
        mj = Float64(snap.mass[j])
        n_sys += 1
        push!(pair_index, i)
        m_sys[n_sys] = mi + mj
        for k in 1:3
            v_sys[k, n_sys] = (mi * snap.vel[k, i] + mj * snap.vel[k, j]) / (mi + mj)
        end
        is_component[i] = true
        is_component[j] = true
    end
    n_pairs_reduced = n_sys
    for i in 1:N
        is_component[i] && continue
        n_sys += 1
        push!(single_index, i)
        m_sys[n_sys] = Float64(snap.mass[i])
        for k in 1:3
            v_sys[k, n_sys] = Float64(snap.vel[k, i])
        end
    end

    sel = collect(1:n_sys)
    if bound_only && n_sys ≥ _MIN_MEMBERS
        p_sys = Matrix{Float64}(undef, 3, n_sys)
        for s in 1:n_sys
            i = s ≤ n_pairs_reduced ? pair_index[s] : single_index[s - n_pairs_reduced]
            for k in 1:3
                p_sys[k, s] = Float64(snap.pos[k, i])
            end
        end
        sel = _bound_members(p_sys, v_sys[:, 1:n_sys], m_sys[1:n_sys])
        length(sel) ≥ _MIN_MEMBERS || (sel = collect(1:n_sys))
    end
    M = sum(m_sys[s] for s in sel)
    M > 0 || throw(ArgumentError("hardness_scale: snapshot has zero total mass"))
    vc = ntuple(k -> sum(m_sys[s] * v_sys[k, s] for s in sel) / M, 3)
    σ2 = 0.0
    for s in sel
        σ2 += m_sys[s] * sum((v_sys[k, s] - vc[k])^2 for k in 1:3)
    end
    σ2 /= 3M

    return (; m_mean = M / length(sel) * zm, sigma_kms = sqrt(σ2) * vs, n_systems = length(sel))
end

"""
    binary_hardness(bev::BinaryEvolutionSnapshot, m_mean, sigma_kms) -> Vector{Float64}

Hardness ratio ``x = E_b / (⟨m⟩ σ²)`` of every pair in `bev` for a mean
system mass `m_mean` [M☉] and one-dimensional dispersion `sigma_kms`
[km s⁻¹]; a pair is hard when `x > 1` (Heggie 1975).
"""
function binary_hardness(bev::BinaryEvolutionSnapshot, m_mean::Real, sigma_kms::Real)
    m_mean > 0 || throw(ArgumentError("binary_hardness: m_mean must be positive, got $m_mean"))
    sigma_kms > 0 ||
        throw(ArgumentError("binary_hardness: sigma_kms must be positive, got $sigma_kms"))
    scale = m_mean * sigma_kms^2
    return [binding_energy(r) / scale for r in bev.records]
end

"""
    BinaryPopulation

Time series of the regularised-binary population from a sequence of
bev.82 snapshots (see [`binary_population`](@ref)).

# Fields
- `time_myr`: epochs [Myr]
- `n_pairs`: number of regularised pairs listed at each epoch
- `n_hard`, `n_soft`: hard/soft split by the Heggie criterion (zeros when
  `classified == false`)
- `binary_fraction`: ``f_\\mathrm{b} = N_\\mathrm{b} / (N_\\mathrm{s} + N_\\mathrm{b})``,
  pairs over systems; `NaN` when the stellar count was not supplied
- `classified`: whether a hard/soft energy scale was available
"""
struct BinaryPopulation
    time_myr::Vector{Float64}
    n_pairs::Vector{Int}
    n_hard::Vector{Int}
    n_soft::Vector{Int}
    binary_fraction::Vector{Float64}
    classified::Bool
end

_per_epoch(::Nothing, n::Int, ::String) = nothing
_per_epoch(x::Real, n::Int, ::String) = fill(Float64(x), n)
function _per_epoch(x::AbstractVector, n::Int, key::String)
    length(x) == n ||
        throw(DimensionMismatch("$key has $(length(x)) entries for $n binary snapshots"))
    return Float64.(x)
end

"""
    binary_population(bevs::Vector{BinaryEvolutionSnapshot};
                      n_stars = nothing, m_mean = nothing, sigma_kms = nothing)
        -> BinaryPopulation

Reduce a time-ordered sequence of bev.82 snapshots to pair counts, the
binary fraction and the hard/soft split. `n_stars` (total number of stars,
scalar or one value per epoch) enables the binary fraction
``N_\\mathrm{b} / (N_\\star - N_\\mathrm{b})``; `m_mean` [M☉] and
`sigma_kms` [km s⁻¹] (scalars or per-epoch vectors, e.g. from
[`hardness_scale`](@ref)) enable the classification.
"""
function binary_population(
    bevs::AbstractVector{BinaryEvolutionSnapshot};
    n_stars = nothing,
    m_mean = nothing,
    sigma_kms = nothing,
)
    n = length(bevs)
    n > 0 || throw(ArgumentError("binary_population: no binary snapshots"))
    ns = _per_epoch(n_stars, n, "n_stars")
    mm = _per_epoch(m_mean, n, "m_mean")
    ss = _per_epoch(sigma_kms, n, "sigma_kms")
    classified = mm !== nothing && ss !== nothing

    t = [b.time_myr for b in bevs]
    n_pairs = [length(b.records) for b in bevs]
    n_hard = zeros(Int, n)
    f_bin = fill(NaN, n)
    for (i, b) in enumerate(bevs)
        if classified
            n_hard[i] = count(>(1.0), binary_hardness(b, mm[i], ss[i]))
        end
        if ns !== nothing
            n_systems = ns[i] - n_pairs[i]
            f_bin[i] = n_systems > 0 ? n_pairs[i] / n_systems : NaN
        end
    end
    return BinaryPopulation(t, n_pairs, n_hard, n_pairs .- n_hard, f_bin, classified)
end

"""
    binary_scales(bevs, snaps::Vector{Snapshot};
                  pair_sum_max_n = PostprocessConfig().pair_sum_max_n)
        -> (; m_mean, sigma_kms, n_stars)

Per-epoch hard/soft energy scales and stellar counts for `bevs`, each taken
from the snapshot in `snaps` closest in physical time. The scale of a
snapshot of at most `pair_sum_max_n` particles is taken over its bound
systems ([`hardness_scale`](@ref) with `bound_only`), an O(N²) selection;
above that limit it is taken over every system. Returns `nothing` when
`snaps` is empty.
"""
function binary_scales(
    bevs::AbstractVector{BinaryEvolutionSnapshot},
    snaps::AbstractVector{Snapshot};
    pair_sum_max_n::Integer = PostprocessConfig().pair_sum_max_n,
)
    isempty(snaps) && return nothing
    t_snap = [time_myr(s.header) for s in snaps]
    n = length(bevs)
    m_mean = Vector{Float64}(undef, n)
    sigma = Vector{Float64}(undef, n)
    n_stars = Vector{Int}(undef, n)
    for (i, b) in enumerate(bevs)
        j = argmin(abs.(t_snap .- b.time_myr))
        sc = hardness_scale(snaps[j], b; bound_only = nparticles(snaps[j]) ≤ pair_sum_max_n)
        m_mean[i] = sc.m_mean
        sigma[i] = sc.sigma_kms
        n_stars[i] = nparticles(snaps[j])
    end
    return (; m_mean, sigma_kms = sigma, n_stars)
end
