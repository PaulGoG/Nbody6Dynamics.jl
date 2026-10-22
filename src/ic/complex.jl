# =============================================================================
# Multi-cluster complexes: member masses, sizes, centres and bulk velocities
# =============================================================================
#
# `orbit_mode = "complex"` builds a configuration of many member clusters from
# a few distribution parameters instead of explicit per-cluster tables. The
# members are drawn from a mass function, sized by a size–mass relation,
# placed by one of four centre distributions, and given a bulk velocity field
# scaled to a target virial ratio of the members treated as softened point
# masses. Units: M☉, pc, km s⁻¹, G = _G_PC_KMS2_MSUN.
# =============================================================================

"""
    ComplexSpec(; n_members = 0, mass_total = 0.0, mass_function = "equal", beta = 2.0,
                  mass_lo = 100.0, mass_hi = 1.0e4, size_mass = "constant", rh = 1.0,
                  size_scatter_dex = 0.0, centres = "uniform_sphere", radius = 5.0,
                  plummer_cutoff = 5.0, fractal_dimension = 1.6, filament_length = 20.0,
                  filament_width = 1.0, virial_ratio = 0.5, rotation_fraction = 0.0,
                  softening = 0.0, member_n_min = 50, imf_partition = "member",
                  control = "none")

Parameters of a cluster complex (`[merger.complex]`, `orbit_mode = "complex"`). Masses in
M☉, lengths in pc, speeds in km s⁻¹.

| Field | Type | Default | Meaning, bounds |
|---|---|---|---|
| `n_members` | Int | 0 | number of members; ≥ 2 when > 0 |
| `mass_total` | Float64 | 0.0 | target total mass [M☉]; > 0 when used |
| `mass_function` | String | "equal" | "equal" or "power_law" |
| `beta` | Float64 | 2.0 | dN/dM ∝ M^−β; finite (power_law) |
| `mass_lo`, `mass_hi` | Float64 | 100.0, 1.0e4 | bounds [M☉]; 0 < lo < hi (power_law) |
| `size_mass` | String | "constant" | "constant", "marks_kroupa", "frost", "dobbs" |
| `rh` | Float64 | 1.0 | half-mass radius for "constant" [pc]; > 0 |
| `size_scatter_dex` | Float64 | 0.0 | log-normal scatter of r_h [dex]; ≥ 0 |
| `centres` | String | "uniform_sphere" | "uniform_sphere", "plummer", "fractal", "filament" |
| `radius` | Float64 | 5.0 | sphere radius / Plummer scale a / fractal radius [pc]; > 0 |
| `plummer_cutoff` | Float64 | 5.0 | Plummer truncation in units of a; > 1 |
| `fractal_dimension` | Float64 | 1.6 | box-fractal dimension D; 1 ≤ D ≤ 3 |
| `filament_length` | Float64 | 20.0 | [pc]; > 0 |
| `filament_width` | Float64 | 1.0 | Gaussian transverse σ [pc]; > 0 |
| `virial_ratio` | Float64 | 0.5 | target Q_cc = T/\\|W\\| of the centres; ≥ 0 |
| `rotation_fraction` | Float64 | 0.0 | T_rot/T of a solid-body rotation about +z; 0 ≤ f < 1; 0 when virial_ratio = 0 |
| `softening` | Float64 | 0.0 | Plummer softening ε of the point-mass energy [pc]; 0 = the mean member r_h; ≥ 0 |
| `member_n_min` | Int | 50 | smallest allowed member star count; ≥ 2 |
| `imf_partition` | String | "member" | "member" (each member samples its own IMF) or "global" (one star list partitioned) |
| `control` | String | "none" | "none", "monolithic" (C-E) or "cold" (C-cold); needs `imf_partition = "global"` |

# Member masses
`"equal"` needs `n_members > 0` and `mass_total > 0` and gives every member
`mass_total / n_members`. `"power_law"` takes exactly one of them: `n_members`
draws, or draws until the running sum reaches or exceeds `mass_total` (the last
draw is kept). The inverse CDF of dN/dM ∝ M^−β on `[lo, hi]` with `u` uniform in
`[0, 1)` is `M = (lo^(1−β) + u (hi^(1−β) − lo^(1−β)))^(1/(1−β))` for β ≠ 1 and
`M = lo (hi/lo)^u` for β = 1.

# Size–mass relations
`r_h = rh` (`"constant"`); `r_h = 0.10 pc (M/M☉)^0.13` (`"marks_kroupa"`, Marks &
Kroupa 2012, Eq. 7, DOI 10.1051/0004-6361/201118231); `r_h = (2.365/8) pc
(M/10⁴ M☉)^0.18` (`"frost"`, Rantala, Naab & Lahén 2024, DOI
10.1093/mnras/stae1413, with f_h = 1/8 and R₄ = 2.365 pc; 0.45 pc at 10⁵ M☉);
`r_h = 2.26 pc (M/10⁴ M☉)^0.27` (`"dobbs"`, Dobbs et al. 2022, DOI
10.1093/mnras/stac2474). Each radius is multiplied by `10^(size_scatter_dex g)`,
`g` standard normal.

# IMF partition
`"member"`: member i has `Nᵢ = round(Mᵢ / ⟨m⟩)` stars drawn from the template
IMF. `"global"`: one star list is drawn from the template IMF in batches of
`max(1000, ⌈1.1 ΣMᵢ / ⟨m⟩⌉)` until its sum reaches ΣMᵢ, shuffled, and dealt
out in order: member i takes stars until the running total reaches or exceeds
the cumulative target Σ_{j≤i} Mⱼ, the last member until it reaches ΣMᵢ, so
every member is within one star of its target, and the rest of the list is
dropped (`_complex_partition_stars`). Each member then carries its
share as an `ExplicitMassesIMF`, and its radius follows from its realised mass.

# Controls
`"monolithic"` (C-E) and `"cold"` (C-cold) replace the complex by one Plummer
cluster of the same star list, built from the same seed, whose positions and
velocities are rescaled about the centre of mass to the total energy E of the
complex's systems (primordial pairs as point masses), keeping Q = 0.5 (C-E) or
matching also the global Q = T/|W| of the complex (C-cold). The internal
energy of the pairs is not matched (the control draws its own pairs), template
pairings that change the stellar masses are refused, and the control is not
segregated.

# Centre distributions
- `"uniform_sphere"`: `r = radius u^(1/3)`, isotropic direction; velocities
  standard normal per component.
- `"plummer"`: `a = radius`, `r = a (u^(−2/3) − 1)^(−1/2)` redrawn while
  `r > plummer_cutoff a`; speed `q (r² + a²)^(−1/4)` with
  `g(q) ∝ q² (1 − q²)^(7/2)` by rejection, isotropic velocity direction
  (Aarseth, Hénon & Wielen 1974, A&A 37, 183).
- `"fractal"`: the box fractal of Goodwin & Whitworth (2004, DOI
  10.1051/0004-6361:20031529) of dimension D: each parent cube of side L spawns
  8 children at ± L/4 per coordinate with uniform noise of ± 0.05 L, inheriting
  the parent velocity plus a Gaussian component of standard deviation 2^(−g)
  in generation g (halved per generation); a child is a parent of the next generation with probability
  2^(D − 3). The construction stops at the first generation with at least
  `n` parents inside the unit sphere, of which `n` are chosen at random and
  scaled by `radius`.
- `"filament"`: `x` uniform in `[−L/2, L/2]` (`L = filament_length`), `y` and
  `z` normal with σ = `filament_width`; velocities standard normal.

# Velocity field
With `W = −Σ_{i<j} G mᵢ mⱼ / √(|rᵢ − rⱼ|² + ε²)` of the centred members and
`w = ẑ × r`, the raw velocities lose their mass-weighted projection on `w` and
become `α v + β w` with `α = √((1 − f) Q |W| / T_v)` and
`β = √(f Q |W| / T_w)`, so that `T = Q |W|` and the solid-body rotation carries
the fraction `f = rotation_fraction` of `T` (`complex_scale_velocities!`).
The spin parameter λ = |L| √|E| / (G M^(5/2)) of the centres follows Peebles
(1969, DOI 10.1086/149876).
"""
Base.@kwdef struct ComplexSpec
    n_members::Int = 0
    mass_total::Float64 = 0.0
    mass_function::String = "equal"
    beta::Float64 = 2.0
    mass_lo::Float64 = 100.0
    mass_hi::Float64 = 1.0e4
    size_mass::String = "constant"
    rh::Float64 = 1.0
    size_scatter_dex::Float64 = 0.0
    centres::String = "uniform_sphere"
    radius::Float64 = 5.0
    plummer_cutoff::Float64 = 5.0
    fractal_dimension::Float64 = 1.6
    filament_length::Float64 = 20.0
    filament_width::Float64 = 1.0
    virial_ratio::Float64 = 0.5
    rotation_fraction::Float64 = 0.0
    softening::Float64 = 0.0
    member_n_min::Int = 50
    imf_partition::String = "member"
    control::String = "none"
end

"""Mass functions of `[merger.complex]` (the `mass_function` key)."""
const _COMPLEX_MASS_FUNCTIONS = ("equal", "power_law")

"""IMF partitions of `[merger.complex]` (the `imf_partition` key)."""
const _COMPLEX_IMF_PARTITIONS = ("member", "global")

"""Monolithic controls of `[merger.complex]` (the `control` key)."""
const _COMPLEX_CONTROLS = ("none", "monolithic", "cold")

"""Size–mass relations of `[merger.complex]` (the `size_mass` key)."""
const _COMPLEX_SIZE_MASS = ("constant", "marks_kroupa", "frost", "dobbs")

"""Centre distributions of `[merger.complex]` (the `centres` key)."""
const _COMPLEX_CENTRES = ("uniform_sphere", "plummer", "fractal", "filament")

"""Generation cap of the box-fractal construction of [`complex_centres`](@ref)."""
const _FRACTAL_MAX_GENERATIONS = 20
# Realisations of the box fractal drawn before giving up (extinct or too sparse ones are redrawn)
const _FRACTAL_MAX_ATTEMPTS = 100

"""Inverse CDF of dN/dM ∝ M^−β on `[mass_lo, mass_hi]` at `u ∈ [0, 1)`, clamped to the bounds against rounding."""
function _power_law_mass(spec::ComplexSpec, u::Float64)::Float64
    lo, hi, β = spec.mass_lo, spec.mass_hi, spec.beta
    M = if β == 1
        lo * (hi / lo)^u
    else
        p = 1 - β
        (lo^p + u * (hi^p - lo^p))^(1 / p)
    end
    return clamp(M, lo, hi)
end

"""
    complex_member_masses(spec::ComplexSpec, rng::AbstractRNG) -> Vector{Float64}

Member masses [M☉] of a complex: `n_members` equal masses `mass_total / n_members`
(`"equal"`), or draws from dN/dM ∝ M^−β on `[mass_lo, mass_hi]` (`"power_law"`),
either `n_members` of them or as many as it takes for the sum to reach
`mass_total` (the last draw kept). Throws an `ArgumentError` for an unknown
mass function or a count/target combination the mass function does not take.

# Example
```julia
complex_member_masses(ComplexSpec(; n_members = 4, mass_total = 2000.0), rng)  # 4 × 500.0
```
"""
function complex_member_masses(spec::ComplexSpec, rng::AbstractRNG)::Vector{Float64}
    has_n, has_m = spec.n_members > 0, spec.mass_total > 0
    if spec.mass_function == "equal"
        (has_n && has_m) || throw(
            ArgumentError(
                "complex_member_masses: mass_function = \"equal\" needs n_members > 0 and " *
                "mass_total > 0; got $(spec.n_members), $(spec.mass_total)",
            ),
        )
        return fill(spec.mass_total / spec.n_members, spec.n_members)
    elseif spec.mass_function == "power_law"
        xor(has_n, has_m) || throw(
            ArgumentError(
                "complex_member_masses: mass_function = \"power_law\" takes exactly one of " *
                "n_members > 0 and mass_total > 0; got $(spec.n_members), $(spec.mass_total)",
            ),
        )
        has_n && return Float64[_power_law_mass(spec, rand(rng)) for _ in 1:spec.n_members]
        masses = Float64[]
        total = 0.0
        while total < spec.mass_total
            M = _power_law_mass(spec, rand(rng))
            push!(masses, M)
            total += M
        end
        return masses
    end
    throw(
        ArgumentError(
            "complex_member_masses: unknown mass_function \"$(spec.mass_function)\"; " *
            "one of " *
            join(_COMPLEX_MASS_FUNCTIONS, ", "),
        ),
    )
end

"""
    complex_member_radius(spec::ComplexSpec, M::Real) -> Float64

Half-mass radius [pc] of a member of mass `M` [M☉] from the size–mass relation
`spec.size_mass`, without scatter: `spec.rh` (`"constant"`), `0.10 M^0.13`
(`"marks_kroupa"`), `(2.365/8) (M/10⁴)^0.18` (`"frost"`), `2.26 (M/10⁴)^0.27`
(`"dobbs"`); see [`ComplexSpec`](@ref) for the sources.
"""
function complex_member_radius(spec::ComplexSpec, M::Real)::Float64
    m = Float64(M)
    s = spec.size_mass
    s == "constant" && return spec.rh
    s == "marks_kroupa" && return 0.10 * m^0.13
    s == "frost" && return (2.365 / 8) * (m / 1.0e4)^0.18
    s == "dobbs" && return 2.26 * (m / 1.0e4)^0.27
    throw(
        ArgumentError(
            "complex_member_radius: unknown size_mass \"$s\"; one of " *
            join(_COMPLEX_SIZE_MASS, ", "),
        ),
    )
end

"""
    _complex_partition_stars(imf, M::Vector{Float64}, rng::AbstractRNG) -> Vector{Vector{Float64}}

Global IMF partition of a complex with member masses `M` [M☉]: batches of
`max(1000, ⌈1.1 ΣM / ⟨m⟩⌉)` stars are drawn from the IMF specification `imf`
(`⟨m⟩ = expected_mass(imf, 1)`) and appended until the list sum reaches ΣM; the
list is shuffled and dealt out in order. Member i < n takes stars until the
running total of the list reaches or exceeds the cumulative target
`Σ_{j≤i} M[j]`; the last member takes stars until the running total reaches
ΣM. Every member is thus within one star mass of its target, the last
included. Stars beyond that point are dropped, and a member reached after the
list is exhausted keeps what it has. All draws come from `rng`.
"""
function _complex_partition_stars(imf, M::Vector{Float64}, rng::AbstractRNG)
    n = length(M)
    M_sum = sum(M)
    n_batch = max(1000, ceil(Int, 1.1 * M_sum / expected_mass(imf, 1)))
    list = Float64[]
    total = 0.0
    while total < M_sum
        batch::Vector{Float64} = sample_masses(imf, n_batch, rng)
        append!(list, batch)
        total += sum(batch)
    end
    Random.shuffle!(rng, list)
    shares = Vector{Vector{Float64}}(undef, n)
    j = 1
    running = 0.0
    target = 0.0
    for i in 1:n
        share = Float64[]
        target = i < n ? target + M[i] : M_sum
        while j ≤ length(list) && running < target
            m = list[j]
            push!(share, m)
            running += m
            j += 1
        end
        shares[i] = share
    end
    return shares
end

"""Isotropic unit vector: cos θ uniform in [−1, 1], φ uniform in [0, 2π)."""
function _isotropic_direction(rng::AbstractRNG)::NTuple{3,Float64}
    cosθ = 2 * rand(rng) - 1
    sinθ = sqrt(max(1 - cosθ^2, 0.0))
    φ = 2π * rand(rng)
    return (sinθ * cos(φ), sinθ * sin(φ), cosθ)
end

"""Speed fraction `q` of the Plummer distribution function, `g(q) ∝ q² (1 − q²)^(7/2)`, by rejection (Aarseth, Hénon & Wielen 1974)."""
function _plummer_speed_fraction(rng::AbstractRNG)::Float64
    while true
        q = rand(rng)
        y = 0.1 * rand(rng)
        y < q^2 * (1 - q^2)^3.5 && return q
    end
end

function _sphere_centres(spec::ComplexSpec, n::Int, rng::AbstractRNG)
    pos = zeros(3, n)
    for i in 1:n
        r = spec.radius * cbrt(rand(rng))
        d = _isotropic_direction(rng)
        for k in 1:3
            pos[k, i] = r * d[k]
        end
    end
    return pos, randn(rng, 3, n)
end

function _plummer_centres(spec::ComplexSpec, n::Int, rng::AbstractRNG)
    a = spec.radius
    r_max = spec.plummer_cutoff * a
    pos = zeros(3, n)
    vel = zeros(3, n)
    for i in 1:n
        r = Inf
        while !(r ≤ r_max)
            r = a / sqrt(rand(rng)^(-2 / 3) - 1)
        end
        d = _isotropic_direction(rng)
        speed = _plummer_speed_fraction(rng) * (r^2 + a^2)^(-0.25)
        e = _isotropic_direction(rng)
        for k in 1:3
            pos[k, i] = r * d[k]
            vel[k, i] = speed * e[k]
        end
    end
    return pos, vel
end

# One realisation of the box fractal; `nothing` when the branching process dies out
# or stays below `n` parents inside the unit sphere
function _fractal_attempt(spec::ComplexSpec, n::Int, rng::AbstractRNG)
    P = 2.0^(spec.fractal_dimension - 3)
    parent_pos::Matrix{Float64} = zeros(3, 1)
    parent_vel::Matrix{Float64} = zeros(3, 1)
    L = 2.0
    x = zeros(3)
    w = zeros(3)
    for g in 1:_FRACTAL_MAX_GENERATIONS
        # Goodwin & Whitworth (2004): the random component shrinks by 1/N_div = 1/2 per generation
        σ = 2.0^(-g)
        child_pos = Float64[]
        child_vel = Float64[]
        for j in axes(parent_pos, 2), octant in 0:7
            for k in 1:3
                sgn = isodd(octant >> (k - 1)) ? 1.0 : -1.0
                x[k] = parent_pos[k, j] + sgn * L / 4 + 0.1 * L * (rand(rng) - 0.5)
            end
            for k in 1:3
                w[k] = parent_vel[k, j] + σ * randn(rng)
            end
            if rand(rng) < P
                append!(child_pos, x)
                append!(child_vel, w)
            end
        end
        parent_pos = reshape(child_pos, 3, :)
        parent_vel = reshape(child_vel, 3, :)
        L /= 2
        isempty(child_pos) && break
        inside = Int[
            i for i in axes(parent_pos, 2) if
            parent_pos[1, i]^2 + parent_pos[2, i]^2 + parent_pos[3, i]^2 < 1
        ]
        if length(inside) ≥ n
            pick = inside[Random.randperm(rng, length(inside))[1:n]]
            return parent_pos[:, pick] .* spec.radius, parent_vel[:, pick]
        end
    end
    return nothing
end

# The branching process can die out (≈ 2.5 % of the roots at D = 1.6: 0.621⁸ ≈ 2.2 % in the
# first generation, the rest later, the fixed point of q = (1 − P + P q)⁸): a new
# realisation is drawn from the same stream, so the result stays a function of the seed.
function _fractal_centres(spec::ComplexSpec, n::Int, rng::AbstractRNG)
    for _ in 1:_FRACTAL_MAX_ATTEMPTS
        out = _fractal_attempt(spec, n, rng)
        out === nothing || return out
    end
    throw(
        ArgumentError(
            "complex_centres: the box fractal of dimension $(spec.fractal_dimension) did not " *
            "reach $n parents inside the unit sphere within $(_FRACTAL_MAX_GENERATIONS) " *
            "generations in $(_FRACTAL_MAX_ATTEMPTS) realisations",
        ),
    )
end

function _filament_centres(spec::ComplexSpec, n::Int, rng::AbstractRNG)
    pos = zeros(3, n)
    for i in 1:n
        pos[1, i] = spec.filament_length * (rand(rng) - 0.5)
        pos[2, i] = spec.filament_width * randn(rng)
        pos[3, i] = spec.filament_width * randn(rng)
    end
    return pos, randn(rng, 3, n)
end

"""
    complex_centres(spec::ComplexSpec, n::Int, rng::AbstractRNG)
        -> (pos::Matrix{Float64}, vel::Matrix{Float64})

Positions [pc] (3 × `n`) of the member centres from the distribution
`spec.centres` (`"uniform_sphere"`, `"plummer"`, `"fractal"`, `"filament"`;
[`ComplexSpec`](@ref)) and their raw velocities (3 × `n`, arbitrary units, to
be scaled by [`complex_scale_velocities!`](@ref)). Throws an `ArgumentError`
for an unknown distribution, `n < 1`, or a box fractal that does not reach `n`
parents inside the unit sphere within $(_FRACTAL_MAX_GENERATIONS) generations in
$(_FRACTAL_MAX_ATTEMPTS) realisations (an extinct realisation is redrawn from the same stream).
"""
function complex_centres(
    spec::ComplexSpec,
    n::Int,
    rng::AbstractRNG,
)::Tuple{Matrix{Float64},Matrix{Float64}}
    n ≥ 1 || throw(ArgumentError("complex_centres: n must be ≥ 1; got $n"))
    c = spec.centres
    c == "uniform_sphere" && return _sphere_centres(spec, n, rng)
    c == "plummer" && return _plummer_centres(spec, n, rng)
    c == "fractal" && return _fractal_centres(spec, n, rng)
    c == "filament" && return _filament_centres(spec, n, rng)
    throw(
        ArgumentError(
            "complex_centres: unknown centres \"$c\"; one of " * join(_COMPLEX_CENTRES, ", "),
        ),
    )
end

"""
    complex_scale_velocities!(vel, pos, mass; virial_ratio, rotation_fraction = 0.0,
                              softening) -> NamedTuple

Scale the bulk velocities `vel` [km s⁻¹] (3 × n) of point masses `mass` [M☉] at
`pos` [pc] (3 × n) to the virial ratio `Q = virial_ratio` of the softened
point-mass system, with the fraction `f = rotation_fraction` of the kinetic
energy in a solid-body rotation about +z. In place:

1. `pos` and `vel` are shifted to their mass-weighted means.
2. `W = −Σ_{i<j} G mᵢ mⱼ / √(|rᵢ − rⱼ|² + ε²)`, `ε = softening > 0`.
3. `Q = 0` sets `vel` to zero (`f` must be 0).
4. Otherwise, with `wᵢ = ẑ × rᵢ`, `vel` loses its mass-weighted projection on
   `w`, `vel ← vel − (Σ m v·w / Σ m |w|²) w`, and becomes `α vel + β w` with
   `α = √((1 − f) Q |W| / T_v)`, `β = √(f Q |W| / T_w)`, where
   `T_v = ½ Σ m |v|²` and `T_w = ½ Σ m |w|²`.

Returns `(virial_ratio = T/|W|, e_total = T + W, l_total = (Lx, Ly, Lz),
rotation_fraction = β² T_w / T, kinetic = T, potential = W, softening = ε)`
of the final velocities, energies in M☉ km² s⁻², `L = Σ m r × v` in
M☉ pc km s⁻¹. Throws an `ArgumentError` for `ε ≤ 0`, `Q < 0`, `f` outside
`[0, 1)`, `f > 0` with `Q = 0` or with all members on the z axis, fewer than two
members, and vanishing velocities to scale (`T_v = 0`).

# Example
```julia
out = complex_scale_velocities!(vel, pos, mass; virial_ratio = 0.5, softening = 0.5)
out.virial_ratio   # 0.5
```
"""
function complex_scale_velocities!(
    vel::AbstractMatrix{Float64},
    pos::AbstractMatrix{Float64},
    mass::AbstractVector{Float64};
    virial_ratio::Real,
    rotation_fraction::Real = 0.0,
    softening::Real,
)::NamedTuple
    n = length(mass)
    (size(pos) == (3, n) && size(vel) == (3, n)) || throw(
        DimensionMismatch(
            "complex_scale_velocities!: pos is $(size(pos)), vel is $(size(vel)) for $n masses",
        ),
    )
    n ≥ 2 || throw(ArgumentError("complex_scale_velocities!: needs ≥ 2 members; got $n"))
    ε = Float64(softening)
    ε > 0 || throw(ArgumentError("complex_scale_velocities!: softening must be > 0 [pc]; got $ε"))
    Q = Float64(virial_ratio)
    f = Float64(rotation_fraction)
    (isfinite(Q) && Q ≥ 0) ||
        throw(ArgumentError("complex_scale_velocities!: virial_ratio must be ≥ 0; got $Q"))
    (0 ≤ f < 1) || throw(
        ArgumentError("complex_scale_velocities!: rotation_fraction must be in [0, 1); got $f"),
    )

    M = sum(mass)
    for k in 1:3
        cp = 0.0
        cv = 0.0
        for i in 1:n
            cp += mass[i] * pos[k, i]
            cv += mass[i] * vel[k, i]
        end
        cp /= M
        cv /= M
        for i in 1:n
            pos[k, i] -= cp
            vel[k, i] -= cv
        end
    end

    G = _G_PC_KMS2_MSUN
    W = 0.0
    for i in 1:(n - 1), j in (i + 1):n
        d² = (pos[1, i] - pos[1, j])^2 + (pos[2, i] - pos[2, j])^2 + (pos[3, i] - pos[3, j])^2
        W -= G * mass[i] * mass[j] / sqrt(d² + ε^2)
    end

    β = 0.0
    T_w = 0.0
    if Q == 0
        f == 0 || throw(
            ArgumentError(
                "complex_scale_velocities!: rotation_fraction must be 0 with virial_ratio = 0; got $f",
            ),
        )
        fill!(vel, 0.0)
    else
        # Projection of vel on w = ẑ × r = (−y, x, 0), mass-weighted
        num = 0.0
        den = 0.0
        for i in 1:n
            num += mass[i] * (-pos[2, i] * vel[1, i] + pos[1, i] * vel[2, i])
            den += mass[i] * (pos[1, i]^2 + pos[2, i]^2)
        end
        (f > 0 && den == 0) && throw(
            ArgumentError(
                "complex_scale_velocities!: rotation about +z needs members off the z axis",
            ),
        )
        c = den > 0 ? num / den : 0.0
        for i in 1:n
            vel[1, i] += c * pos[2, i]
            vel[2, i] -= c * pos[1, i]
        end
        T_target = Q * abs(W)
        T_v = 0.0
        for i in 1:n
            T_v += 0.5 * mass[i] * (vel[1, i]^2 + vel[2, i]^2 + vel[3, i]^2)
        end
        T_w = 0.5 * den
        T_v > 0 || throw(
            ArgumentError(
                "complex_scale_velocities!: the raw velocities have no kinetic energy to scale",
            ),
        )
        α = sqrt((1 - f) * T_target / T_v)
        β = f > 0 ? sqrt(f * T_target / T_w) : 0.0
        for i in 1:n
            vel[1, i] = α * vel[1, i] - β * pos[2, i]
            vel[2, i] = α * vel[2, i] + β * pos[1, i]
            vel[3, i] = α * vel[3, i]
        end
    end

    T = 0.0
    Lx = 0.0
    Ly = 0.0
    Lz = 0.0
    for i in 1:n
        x, y, z = pos[1, i], pos[2, i], pos[3, i]
        vx, vy, vz = vel[1, i], vel[2, i], vel[3, i]
        T += 0.5 * mass[i] * (vx^2 + vy^2 + vz^2)
        Lx += mass[i] * (y * vz - z * vy)
        Ly += mass[i] * (z * vx - x * vz)
        Lz += mass[i] * (x * vy - y * vx)
    end
    return (
        virial_ratio = T / abs(W),
        e_total = T + W,
        l_total = (Lx, Ly, Lz),
        rotation_fraction = T > 0 ? β^2 * T_w / T : 0.0,
        kinetic = T,
        potential = W,
        softening = ε,
    )
end
