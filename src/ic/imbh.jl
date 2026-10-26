# =============================================================================
# Intermediate-mass black hole seeding (the kz24 route)
# =============================================================================
#
# The IMBH is the last body of dat.10, and the engine's IMBHINIT line (read
# under KZ(24) = 1 after &INSCALE and &INXTRNL0) repeats the same state, so the
# overwrite of particle N is the identity. Code units of the generator: M☉, pc,
# G = 1, velocity unit _CODE_VSTAR_KMS km s⁻¹. `IMBHSpec` is defined in
# config.jl (a field type of `ClusterSpec`).

"""Gravitational parameter G M☉ [cm³ s⁻²] of the engine's `units.f` (`GM = 6.6743D-08*1.9884D+33`)."""
const _ENGINE_GM_CGS = 6.6743e-8 * 1.9884e33

"""Parsec [cm] of the engine's `units.f` (`PC = 3.0856776D+18`)."""
const _ENGINE_PC_CGS = 3.0856776e18

"""
Speed floor of the IMBH relative to the cluster mean, in units of the central
dispersion σ_c: the engine's start-up rejects a body whose speed is below
≈ 10⁻¹⁶ N-body units (an invalid first time step aborts the block-step sort),
so a velocity within this floor of the cluster mean is replaced by a velocity
of this magnitude in a random direction (10⁻¹² of the Brownian equilibrium
kinetic energy).
"""
const _IMBH_REST_SPEED_FLOOR = 1.0e-6

"""
    _imbh_engine_vstar_kms(mass_scale_msun, rbar) -> Float64

Velocity unit `VSTAR` [km s⁻¹] as the engine forms it in `units.f`,
`1e-5 √(GM/PC) √(ZMASS ZMBAR / RBAR)`, with `mass_scale_msun` the engine's
mass scale `ZMASS · ZMBAR` [M☉] and `rbar` the length unit [pc]. Under
`KZ(22) = 2` the engine rescales the `ZMBAR` it reads (the mean body mass) to
`ZMBAR N / ZMASS` before `UNITS` (`start.F`), so for a `dat.10` whose masses
sum to one the mass scale is the total mass `ZMBAR_in · N`.
"""
_imbh_engine_vstar_kms(mass_scale_msun::Real, rbar::Real)::Float64 =
    1.0e-5 * sqrt(_ENGINE_GM_CGS / _ENGINE_PC_CGS) * sqrt(mass_scale_msun / rbar)

"""
    _imbh_centre(pos, vel, mass, members) -> (x_c, v_c)

Mass-weighted centre of mass `x_c` and mean velocity `v_c` (3-vectors) of the
bodies `members` of `pos`, `vel` (3×N) and `mass`.
"""
function _imbh_centre(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    members::AbstractVector{<:Integer},
)
    isempty(members) && throw(ArgumentError("_imbh_centre: no members"))
    x_c = zeros(3)
    v_c = zeros(3)
    M = 0.0
    for i in members
        m = mass[i]
        M += m
        for k in 1:3
            x_c[k] += m * pos[k, i]
            v_c[k] += m * vel[k, i]
        end
    end
    x_c ./= M
    v_c ./= M
    return x_c, v_c
end

"""
    _imbh_central_dispersion(pos, vel, mass, members, x_c, v_c; fraction = 0.1, n_min = 50) -> σ_c

One-dimensional velocity dispersion `σ_c = √(Σ m |v − v_c|² / (3 Σ m))` of
the innermost `max(n_min, round(fraction · length(members)))` members by
distance from `x_c` (all members when there are fewer), in the velocity unit
of `vel`.
"""
function _imbh_central_dispersion(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    members::AbstractVector{<:Integer},
    x_c::AbstractVector{Float64},
    v_c::AbstractVector{Float64};
    fraction::Real = 0.1,
    n_min::Integer = 50,
)::Float64
    isempty(members) && throw(ArgumentError("_imbh_central_dispersion: no members"))
    n_in = min(length(members), max(Int(n_min), round(Int, fraction * length(members))))
    r² = [(pos[1, i] - x_c[1])^2 + (pos[2, i] - x_c[2])^2 + (pos[3, i] - x_c[3])^2 for i in members]
    inner = view(members, partialsortperm(r², 1:n_in))
    M = 0.0
    s = 0.0
    for i in inner
        m = mass[i]
        M += m
        s += m * ((vel[1, i] - v_c[1])^2 + (vel[2, i] - v_c[2])^2 + (vel[3, i] - v_c[3])^2)
    end
    return sqrt(s / (3M))
end

"""
    _imbh_body(spec::IMBHSpec, mass_msun, m_mean, x_c, v_c, σ_c, rng) -> (pos, vel)

Position [pc] and velocity [code unit] of the IMBH of mass `mass_msun` in a
cluster of centre `x_c`, mean velocity `v_c`, central one-dimensional
dispersion `σ_c` [code unit] and mean stellar mass `m_mean` [M☉]. Position:
`x_c`, plus `spec.offset` for `placement = "offset"`. Velocity: `v_c` for
`"rest"`; `v_c + velocity_kms / _CODE_VSTAR_KMS` for `"explicit"`; for
`"equipartition"` a Gaussian draw about `v_c` of per-component variance
`σ_c² m_mean / M•`, the Brownian equilibrium ⟨v²⟩ of a massive particle in a
stellar core (Merritt 2005, Eq. 16, DOI 10.1086/429964), drawn from `rng`.
A velocity within `_IMBH_REST_SPEED_FLOOR σ_c` of `v_c` (the `"rest"` mode of
a cluster at rest, or an explicit zero) is replaced by that speed in a random
direction from `rng`: the engine's start-up rejects a body at exact rest.
"""
function _imbh_body(
    spec::IMBHSpec,
    mass_msun::Real,
    m_mean::Real,
    x_c::AbstractVector{Float64},
    v_c::AbstractVector{Float64},
    σ_c::Real,
    rng::AbstractRNG,
)
    mass_msun > 0 || throw(ArgumentError("_imbh_body: the IMBH mass must be > 0; got $mass_msun"))
    pos = spec.placement == "offset" ? x_c .+ spec.offset : copy(x_c)
    vel = if spec.velocity == "rest"
        copy(v_c)
    elseif spec.velocity == "explicit"
        v_c .+ spec.velocity_kms ./ _CODE_VSTAR_KMS
    elseif spec.velocity == "equipartition"
        v_c .+ σ_c * sqrt(m_mean / mass_msun) .* randn(rng, 3)
    else
        throw(ArgumentError("_imbh_body: unknown velocity mode \"$(spec.velocity)\""))
    end
    v_floor = _IMBH_REST_SPEED_FLOOR * σ_c
    if sqrt(sum(abs2, vel .- v_c)) < v_floor
        u = randn(rng, 3)
        u ./= sqrt(sum(abs2, u))
        vel = v_c .+ v_floor .* u
    end
    return Vector{Float64}(pos), Vector{Float64}(vel)
end

"""
    _imbh_jeans_response!(pos, vel, mass, members, v_c, imbh_mass, imbh_pos; n_shell = 0)
        -> Vector{Float64}

Isotropic Jeans rescale of the velocities of the systems `members` (before
pair expansion; code units, G = 1) for a point mass `imbh_mass` at
`imbh_pos`. With the added potential the isotropic Jeans equation gives

    ρ σ²_tot(r) = ρ σ²_K(r) + G M• ∫_r^∞ ρ(r′) r′⁻² dr′,

whose virial consequence is `2 K_new = −(W_ss + W_s•)` with
`W_s• = −G M• Σ m_i / r_i` (exact in the continuum, since
`3 ∫ P dV = G M• ∫ ρ / r dV`). This approximates the practice of
Baumgardt, Makino & Ebisuzaki (2004a, DOI 10.1086/381968).

The members are sorted by `r_i = |x_i − imbh_pos|` into shells of equal
membership, `n_s = n_shell > 0 ? n_shell : max(50, length(members) ÷ 40)`
each (the last shell takes the remainder; an `ArgumentError` for fewer than
`2 n_s` members). Per shell `k`: the mass-weighted mean radius `r_k`, the
inner radius `r_in,k` (0, else the outer radius of shell `k − 1`) and the
outer radius `r_out,k` (its largest `r`), the density
`ρ_k = Σ m / (4π/3 (r_out,k³ − r_in,k³))` and `σ²_K,k = Σ m |v − v_c|² / (3 Σ m)`.
The pressure term of the point mass,
`P_k = M• [ρ_k (1/r_k − 1/r_out,k) + Σ_{j>k} ρ_j (1/r_in,j − 1/r_out,j)]`,
is the exact integral of the piecewise-constant shell densities from the
shell's mean radius outward, which keeps the first shell finite. Every member
of shell `k` is set to `v_c + f_k (v − v_c)` with
`f_k = √((σ²_K,k + P_k/ρ_k) / σ²_K,k)`. Returns the factors `f_k`.
"""
function _imbh_jeans_response!(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    members::AbstractVector{<:Integer},
    v_c::AbstractVector{Float64},
    imbh_mass::Real,
    imbh_pos::AbstractVector{Float64};
    n_shell::Integer = 0,
)::Vector{Float64}
    n = length(members)
    n_s = n_shell > 0 ? Int(n_shell) : max(50, n ÷ 40)
    n ≥ 2n_s || throw(
        ArgumentError(
            "_imbh_jeans_response!: $n members are fewer than two shells of $n_s; " *
            "lower n_shell or seed the IMBH in a larger cluster",
        ),
    )
    imbh_mass > 0 ||
        throw(ArgumentError("_imbh_jeans_response!: the IMBH mass must be > 0; got $imbh_mass"))
    r = [
        sqrt(
            (pos[1, i] - imbh_pos[1])^2 + (pos[2, i] - imbh_pos[2])^2 + (pos[3, i] - imbh_pos[3])^2,
        ) for i in members
    ]
    order = sortperm(r)
    n_sh = n ÷ n_s
    shell_ranges = [((k - 1) * n_s + 1):(k == n_sh ? n : k * n_s) for k in 1:n_sh]
    ρ = zeros(n_sh)
    r_mean = zeros(n_sh)
    r_in = zeros(n_sh)
    r_out = zeros(n_sh)
    σ²_K = zeros(n_sh)
    for (k, rng_k) in enumerate(shell_ranges)
        M = 0.0
        mr = 0.0
        s = 0.0
        for j in rng_k
            idx = order[j]
            i = members[idx]
            m = mass[i]
            M += m
            mr += m * r[idx]
            s += m * ((vel[1, i] - v_c[1])^2 + (vel[2, i] - v_c[2])^2 + (vel[3, i] - v_c[3])^2)
        end
        r_mean[k] = mr / M
        r_in[k] = k == 1 ? 0.0 : r_out[k - 1]
        r_out[k] = r[order[last(rng_k)]]
        volume = 4π / 3 * (r_out[k]^3 - r_in[k]^3)
        (volume > 0 && r_mean[k] > 0 && s > 0) || throw(
            ArgumentError(
                "_imbh_jeans_response!: shell $k is degenerate (coincident radii or no velocity " *
                "dispersion); choose another n_shell",
            ),
        )
        ρ[k] = M / volume
        σ²_K[k] = s / (3M)
    end
    # Σ_{j>k} ρ_j (1/r_in,j − 1/r_out,j), accumulated from the outermost shell
    tail = zeros(n_sh)
    for k in (n_sh - 1):-1:1
        tail[k] = tail[k + 1] + ρ[k + 1] * (1 / r_in[k + 1] - 1 / r_out[k + 1])
    end
    f = zeros(n_sh)
    for (k, rng_k) in enumerate(shell_ranges)
        P = imbh_mass * (ρ[k] * (1 / r_mean[k] - 1 / r_out[k]) + tail[k])
        f[k] = sqrt((σ²_K[k] + P / ρ[k]) / σ²_K[k])
        for j in rng_k
            i = members[order[j]]
            for c in 1:3
                vel[c, i] = v_c[c] + f[k] * (vel[c, i] - v_c[c])
            end
        end
    end
    return f
end
