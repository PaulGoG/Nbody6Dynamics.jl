# =============================================================================
# Shape, segregation and mixing estimators
# =============================================================================
# Intrinsic shape from the second-moment tensor in ellipsoidal shells,
# projected ellipticities in circular apertures, the A⁺ concentration index
# between a tracer and a reference population, and two-sample measures of
# how far two populations have mixed in radius, energy and angular momentum.

using LinearAlgebra: eigen, Symmetric, dot

"""
    ShapeShell

Axis ratios of one ellipsoidal shell from [`shape_tensor`](@ref). The shell
lies between two similar ellipsoids whose semi-major axes are the Lagrangian
radii (by particle count) that seeded it; `q` and `s` are `NaN` when the shell
fails the count thresholds.

# Fields
- `r_outer`: semi-major axis of the outer boundary (the Lagrangian radius it was seeded with)
- `q`: intermediate-to-major axis ratio `b/a`
- `s`: minor-to-major axis ratio `c/a`
- `axes`: `3 × 3`, unit eigenvectors of the shape tensor as columns, ordered `(a, b, c)`
- `n`: members of the converged shell
- `iterations`: iterations performed
- `converged`: whether `|Δq|` and `|Δs|` fell below the tolerance within the iteration limit
"""
struct ShapeShell
    r_outer::Float64
    q::Float64
    s::Float64
    axes::Matrix{Float64}
    n::Int
    iterations::Int
    converged::Bool
end

"""
    shape_tensor(pos, centre; fractions = (0.1, 0.5, 0.9), n_min_shell = 3000, n_min_enclosed = 7000,
                 tol = 1.0e-3, max_iter = 50) -> Vector{ShapeShell}

Intrinsic axis ratios `q = b/a` and `s = c/a` of the particle distribution
about `centre` in ellipsoidal shells, by the iterative shape tensor of
Zemp et al. 2011 (DOI 10.1088/0067-0049/197/2/30) with fixed semi-major axis.

The radii `r_f` enclosing the `fractions` of the particle count bound the
shells: shell `k` starts as the spherical shell `r_{f_{k−1}} < r ≤ r_{f_k}`
(`r_{f_0} = 0`). At each iteration the members are the particles with
elliptical radius `r_e² = x_a² + x_b²/q² + x_c²/s²` (coordinates in the
current eigenframe) in `(r_{f_{k−1}}, r_{f_k}]`, the tensor is
`S = (1/n) Σ x xᵀ` over their positions relative to `centre`, and its
eigenvalues `a² ≥ b² ≥ c²` give the new `q = b/a`, `s = c/a` and the new
frame. The iteration starts from `q = s = 1` in the coordinate frame and
stops when both ratios change by less than `tol`, or after `max_iter`
iterations (`converged = false`).

The tensor is unweighted: the `r⁻²`-weighted (reduced) tensor gives inner
particles of a shell more weight and biases the result towards rounder
shapes (Allgood et al. 2006, DOI 10.1111/j.1365-2966.2006.10094.x).

A shell needs a few thousand members for the ratios to be meaningful: with
`n` members, a spherical sample returns `1 − E[c/a] ∝ n^{−1/2}`, so small
departures from sphericity are interpreted against a Plummer null of the
same count. Shells with fewer than `n_min_shell` members at convergence, or
fewer than `n_min_enclosed` particles inside their outer ellipsoid, are
returned with `q = s = NaN` and the remaining fields filled.
"""
function shape_tensor(
    pos::AbstractMatrix{Float64},
    centre::AbstractVector{Float64};
    fractions = (0.1, 0.5, 0.9),
    n_min_shell::Integer = 3000,
    n_min_enclosed::Integer = 7000,
    tol::Real = 1.0e-3,
    max_iter::Integer = 50,
)::Vector{ShapeShell}
    size(pos, 1) == 3 ||
        throw(ArgumentError("shape_tensor: pos must be 3 × N; got $(size(pos, 1)) rows"))
    length(centre) == 3 ||
        throw(ArgumentError("shape_tensor: centre must have 3 components; got $(length(centre))"))
    N = size(pos, 2)
    N ≥ 1 || throw(ArgumentError("shape_tensor: need at least one particle"))
    f = collect(Float64, fractions)
    isempty(f) && throw(ArgumentError("shape_tensor: fractions must not be empty"))
    all(x -> 0 < x ≤ 1, f) ||
        throw(ArgumentError("shape_tensor: fractions must lie in (0, 1]; got $(Tuple(f))"))
    all(k -> f[k] > f[k - 1], 2:length(f)) ||
        throw(ArgumentError("shape_tensor: fractions must be strictly ascending; got $(Tuple(f))"))
    n_min_shell ≥ 1 ||
        throw(ArgumentError("shape_tensor: n_min_shell must be ≥ 1; got $n_min_shell"))
    n_min_enclosed ≥ 1 ||
        throw(ArgumentError("shape_tensor: n_min_enclosed must be ≥ 1; got $n_min_enclosed"))
    tol > 0 || throw(ArgumentError("shape_tensor: tol must be > 0; got $tol"))
    max_iter ≥ 1 || throw(ArgumentError("shape_tensor: max_iter must be ≥ 1; got $max_iter"))

    X = Matrix{Float64}(undef, 3, N)
    d = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        x = pos[1, i] - centre[1]
        y = pos[2, i] - centre[2]
        z = pos[3, i] - centre[3]
        X[1, i] = x
        X[2, i] = y
        X[3, i] = z
        d[i] = sqrt(x * x + y * y + z * z)
    end
    sort!(d)
    r_f = [d[clamp(ceil(Int, fk * N), 1, N)] for fk in f]

    shells = Vector{ShapeShell}(undef, length(f))
    for k in eachindex(f)
        # The first shell includes the centre itself (r_e = 0).
        lower² = k == 1 ? -1.0 : r_f[k - 1]^2
        shells[k] = _shape_shell(
            X,
            lower²,
            r_f[k],
            Int(n_min_shell),
            Int(n_min_enclosed),
            Float64(tol),
            Int(max_iter),
        )
    end
    return shells
end

"""
    _shape_shell(X, lower², r_outer, n_min_shell, n_min_enclosed, tol, max_iter) -> ShapeShell

Fixed-semi-major-axis iteration of [`shape_tensor`](@ref) for one shell
`lower² < r_e² ≤ r_outer²`; `X` holds the positions relative to the centre.
"""
function _shape_shell(
    X::Matrix{Float64},
    lower²::Float64,
    r_outer::Float64,
    n_min_shell::Int,
    n_min_enclosed::Int,
    tol::Float64,
    max_iter::Int,
)::ShapeShell
    upper² = r_outer^2
    q = 1.0
    s = 1.0
    V = Float64[i == j for i in 1:3, j in 1:3]
    iterations = 0
    converged = false
    for it in 1:max_iter
        iterations = it
        n, S = _shape_moments(X, V, q, s, lower², upper²)
        n == 0 && break
        F = eigen(Symmetric(S))
        λa, λb, λc = F.values[3], F.values[2], F.values[1]
        # Collinear or coplanar members leave the frame undefined.
        (λa > 0 && λc > 0) || break
        q_new = sqrt(λb / λa)
        s_new = sqrt(λc / λa)
        Δ = max(abs(q_new - q), abs(s_new - s))
        q, s = q_new, s_new
        V = F.vectors[:, [3, 2, 1]]
        if Δ < tol
            converged = true
            break
        end
    end
    n_shell, n_enclosed = _shape_counts(X, V, q, s, lower², upper²)
    masked = n_shell < n_min_shell || n_enclosed < n_min_enclosed
    return ShapeShell(
        r_outer,
        masked ? NaN : q,
        masked ? NaN : s,
        V,
        n_shell,
        iterations,
        converged,
    )
end

"""
    _elliptical_radius²(X, i, V, q⁻², s⁻²) -> Float64

`x_a² + x_b²/q² + x_c²/s²` of particle `i` in the frame whose axes are the
columns of `V`.
"""
@inline function _elliptical_radius²(
    X::Matrix{Float64},
    i::Int,
    V::Matrix{Float64},
    q⁻²::Float64,
    s⁻²::Float64,
)::Float64
    @inbounds begin
        x, y, z = X[1, i], X[2, i], X[3, i]
        xa = V[1, 1] * x + V[2, 1] * y + V[3, 1] * z
        xb = V[1, 2] * x + V[2, 2] * y + V[3, 2] * z
        xc = V[1, 3] * x + V[2, 3] * y + V[3, 3] * z
    end
    return xa * xa + xb * xb * q⁻² + xc * xc * s⁻²
end

"""
    _shape_moments(X, V, q, s, lower², upper²) -> (n, S)

Member count and unweighted second-moment tensor `S = (1/n) Σ x xᵀ` of the
particles with `lower² < r_e² ≤ upper²`.
"""
function _shape_moments(
    X::Matrix{Float64},
    V::Matrix{Float64},
    q::Float64,
    s::Float64,
    lower²::Float64,
    upper²::Float64,
)
    q⁻² = 1 / q^2
    s⁻² = 1 / s^2
    n = 0
    Sxx = Sxy = Sxz = Syy = Syz = Szz = 0.0
    @inbounds for i in 1:size(X, 2)
        r_e² = _elliptical_radius²(X, i, V, q⁻², s⁻²)
        (lower² < r_e² ≤ upper²) || continue
        x, y, z = X[1, i], X[2, i], X[3, i]
        n += 1
        Sxx += x * x
        Sxy += x * y
        Sxz += x * z
        Syy += y * y
        Syz += y * z
        Szz += z * z
    end
    n == 0 && return (0, zeros(3, 3))
    S = [Sxx Sxy Sxz; Sxy Syy Syz; Sxz Syz Szz] ./ n
    return (n, S)
end

"""
    _shape_counts(X, V, q, s, lower², upper²) -> (n_shell, n_enclosed)

Particles in the shell `lower² < r_e² ≤ upper²` and inside its outer
ellipsoid `r_e² ≤ upper²`.
"""
function _shape_counts(
    X::Matrix{Float64},
    V::Matrix{Float64},
    q::Float64,
    s::Float64,
    lower²::Float64,
    upper²::Float64,
)
    q⁻² = 1 / q^2
    s⁻² = 1 / s^2
    n_shell = 0
    n_enclosed = 0
    @inbounds for i in 1:size(X, 2)
        r_e² = _elliptical_radius²(X, i, V, q⁻², s⁻²)
        r_e² ≤ upper² || continue
        n_enclosed += 1
        lower² < r_e² && (n_shell += 1)
    end
    return (n_shell, n_enclosed)
end

"""
    projected_ellipticity(pos, centre, r_h; apertures = (0.25, 0.5, 1.0)) -> Matrix{Float64}

Projected ellipticity `ε = 1 − sqrt(λ_min/λ_max)` from the `2 × 2` second-moment
matrix, about `centre`, of the particles with projected radius `R ≤ f r_h`,
for every aperture fraction `f` in `apertures` (Gavagnin et al. 2016,
DOI 10.1093/mnras/stw2064, measure ε at `R_h/4`, `R_h/2` and `R_h`).

Rows are the projections onto the `(x, y)`, `(x, z)` and `(y, z)` planes,
columns the apertures; `NaN` where fewer than 10 particles fall inside the
aperture. For an elliptical distribution with axis ratio `b/a`, `ε → 1 − b/a`
in apertures well inside the minor axis extent; the circular aperture
truncates the major axis first, so ε of an elongated system decreases
towards large apertures.
"""
function projected_ellipticity(
    pos::AbstractMatrix{Float64},
    centre::AbstractVector{Float64},
    r_h::Real;
    apertures = (0.25, 0.5, 1.0),
)::Matrix{Float64}
    size(pos, 1) == 3 ||
        throw(ArgumentError("projected_ellipticity: pos must be 3 × N; got $(size(pos, 1)) rows"))
    length(centre) == 3 || throw(
        ArgumentError(
            "projected_ellipticity: centre must have 3 components; got $(length(centre))",
        ),
    )
    (r_h > 0 && isfinite(r_h)) ||
        throw(ArgumentError("projected_ellipticity: r_h must be positive and finite; got $r_h"))
    f = collect(Float64, apertures)
    all(x -> x > 0, f) ||
        throw(ArgumentError("projected_ellipticity: apertures must be > 0; got $(Tuple(f))"))
    N = size(pos, 2)
    ε = fill(NaN, 3, length(f))
    for (p, (i, j)) in enumerate(((1, 2), (1, 3), (2, 3))), (k, fk) in enumerate(f)
        R² = (fk * r_h)^2
        n = 0
        Sxx = Sxy = Syy = 0.0
        @inbounds for m in 1:N
            dx = pos[i, m] - centre[i]
            dy = pos[j, m] - centre[j]
            dx * dx + dy * dy ≤ R² || continue
            n += 1
            Sxx += dx * dx
            Sxy += dx * dy
            Syy += dy * dy
        end
        n < 10 && continue
        a = Sxx / n
        b = Sxy / n
        c = Syy / n
        h = 0.5 * (a + c)
        δ = sqrt((0.5 * (a - c))^2 + b^2)
        λ_max = h + δ
        λ_min = max(h - δ, 0.0)
        ε[p, k] = λ_max > 0 ? 1 - sqrt(λ_min / λ_max) : NaN
    end
    return ε
end

"""
    aplus(r_tracer, r_reference, r_h; r_max = r_h)
        -> (; aplus::Float64, error::Float64, n_tracer::Int, n_reference::Int)

Concentration of a tracer population relative to a reference population:
the area between their cumulative radial distributions (Alessandrini et al.
2016, DOI 10.3847/0004-637X/833/2/252) in the `r_h`-normalised form of
Lanzoni et al. 2016 (ApJL 833, L29),

    A⁺ = ∫ [φ_T(x) − φ_R(x)] dx,   x = log10(r / r_h),

where `φ_T` and `φ_R` are the empirical cumulative distributions of the
radii `≤ r_max`, each normalised to 1 at `r_max`, integrated from the
smallest `x` of either sample to `log10(r_max / r_h)`. Positive when the
tracer is more concentrated. Dalessandro et al. 2019
(DOI 10.3847/1538-4357/ab45f7) use `A₂⁺` within `2 r_h` (`r_max = 2 r_h`).

The integrand is a step function, and integrating it exactly on the merged
grid gives every sample point `x_i` the contribution `(x_max − x_i)/n`, so
`A⁺ = ⟨x_R⟩ − ⟨x_T⟩`. The leave-one-out values over the tracer sample follow
from the same sums in O(n); `error` is their jackknife standard deviation
`sqrt((n − 1)/n Σ (A_i − Ā)²)`, `NaN` for fewer than two tracers. Both
values are `NaN` when either sample is empty after the cut.
"""
function aplus(
    r_tracer::AbstractVector{<:Real},
    r_reference::AbstractVector{<:Real},
    r_h::Real;
    r_max::Real = r_h,
)
    (r_h > 0 && isfinite(r_h)) ||
        throw(ArgumentError("aplus: r_h must be positive and finite; got $r_h"))
    (r_max > 0 && isfinite(r_max)) ||
        throw(ArgumentError("aplus: r_max must be positive and finite; got $r_max"))
    all(r -> r > 0, r_tracer) || throw(ArgumentError("aplus: tracer radii must be positive"))
    all(r -> r > 0, r_reference) || throw(ArgumentError("aplus: reference radii must be positive"))
    xT = Float64[log10(r / r_h) for r in r_tracer if r ≤ r_max]
    xR = Float64[log10(r / r_h) for r in r_reference if r ≤ r_max]
    nT = length(xT)
    nR = length(xR)
    (nT == 0 || nR == 0) && return (; aplus = NaN, error = NaN, n_tracer = nT, n_reference = nR)
    ΣT = sum(xT)
    x̄R = sum(xR) / nR
    A = x̄R - ΣT / nT
    nT < 2 && return (; aplus = A, error = NaN, n_tracer = nT, n_reference = nR)
    Ā = 0.0
    for x in xT
        Ā += x̄R - (ΣT - x) / (nT - 1)
    end
    Ā /= nT
    ss = 0.0
    for x in xT
        ss += (x̄R - (ΣT - x) / (nT - 1) - Ā)^2
    end
    σ = sqrt((nT - 1) / nT * ss)
    return (; aplus = A, error = σ, n_tracer = nT, n_reference = nR)
end

"""
    two_sample_d(x1, x2) -> (; d::Float64, floor::Float64, n1::Int, n2::Int)

Two-sample Kolmogorov–Smirnov statistic `D = sup |F₁ − F₂|` of the empirical
distribution functions, evaluated on the merged sorted grid, and its null
floor `0.87 sqrt((n1 + n2)/(n1 n2))`, the asymptotic expectation of `D` for
two samples drawn from one distribution (`sqrt(π/2) ln 2 ≈ 0.87`). The floor
is the level an effect size is compared against, not a p-value. `d` and
`floor` are `NaN` when a sample is empty.
"""
function two_sample_d(x1::AbstractVector{<:Real}, x2::AbstractVector{<:Real})
    any(isnan, x1) && throw(ArgumentError("two_sample_d: x1 contains NaN"))
    any(isnan, x2) && throw(ArgumentError("two_sample_d: x2 contains NaN"))
    n1 = length(x1)
    n2 = length(x2)
    (n1 == 0 || n2 == 0) && return (; d = NaN, floor = NaN, n1, n2)
    a = sort!(Float64.(x1))
    b = sort!(Float64.(x2))
    i = 0
    j = 0
    D = 0.0
    while i < n1 && j < n2
        v = min(a[i + 1], b[j + 1])
        while i < n1 && a[i + 1] == v
            i += 1
        end
        while j < n2 && b[j + 1] == v
            j += 1
        end
        D = max(D, abs(i / n1 - j / n2))
    end
    return (; d = D, floor = _ks_floor(n1, n2), n1, n2)
end

"""
    _ks_floor(n1, n2) -> Float64

Expected two-sample Kolmogorov–Smirnov `D` under the null; `NaN` for an empty sample.
"""
_ks_floor(n1::Int, n2::Int)::Float64 =
    (n1 == 0 || n2 == 0) ? NaN : 0.87 * sqrt((n1 + n2) / (n1 * n2))

"""
    mixing(pos, vel, mass, tags, centre, r_h; energy = nothing)
        -> (; aplus2::Float64, aplus2_error::Float64, d_radius::Float64, d_energy::Float64,
             d_angular_momentum::Float64, floor::Float64, n::NTuple{2,Int})

Degree of mixing of two populations, identified by the two distinct values
of `tags` (the first population is the smaller tag), about `centre`.

- `aplus2`, `aplus2_error`: [`aplus`](@ref) of the first population against
  the second within `2 r_h` (`A₂⁺`), as used for the radial mixing of
  multiple populations by Vesperini et al. 2013 (MNRAS 429, 1913) and
  Dalessandro et al. 2019 (DOI 10.3847/1538-4357/ab45f7)
- `d_radius`: [`two_sample_d`](@ref) of the radii within `2 r_h`
- `d_energy`: the same for the specific energies `energy` (one value per
  particle) of the particles within `2 r_h`; `NaN` without `energy`. Two
  populations can share a radial profile long before they share an energy
  distribution, so this extends the radial test to integrals of motion.
- `d_angular_momentum`: the same for the specific angular momentum
  `|r × (v − v̄)|` within `2 r_h`, with `v̄` the mass-weighted mean velocity
  of all particles within `2 r_h`
- `floor`: null floor of `D` for the two counts within `2 r_h`
- `n`: the two counts within `2 r_h`
"""
function mixing(
    pos::AbstractMatrix{Float64},
    vel::AbstractMatrix{Float64},
    mass::AbstractVector{Float64},
    tags::AbstractVector{<:Integer},
    centre::AbstractVector{Float64},
    r_h::Real;
    energy::Union{Nothing,AbstractVector{<:Real}} = nothing,
)
    size(pos, 1) == 3 || throw(ArgumentError("mixing: pos must be 3 × N; got $(size(pos, 1)) rows"))
    N = size(pos, 2)
    size(vel) == size(pos) ||
        throw(DimensionMismatch("mixing: vel is $(size(vel)), pos is $(size(pos))"))
    length(mass) == N || throw(DimensionMismatch("mixing: $(length(mass)) masses for $N particles"))
    length(tags) == N || throw(DimensionMismatch("mixing: $(length(tags)) tags for $N particles"))
    length(centre) == 3 ||
        throw(ArgumentError("mixing: centre must have 3 components; got $(length(centre))"))
    (r_h > 0 && isfinite(r_h)) ||
        throw(ArgumentError("mixing: r_h must be positive and finite; got $r_h"))
    all(m -> m ≥ 0 && isfinite(m), mass) ||
        throw(ArgumentError("mixing: masses must be finite and non-negative"))
    energy === nothing ||
        length(energy) == N ||
        throw(DimensionMismatch("mixing: $(length(energy)) energies for $N particles"))
    tag_values = sort!(unique(tags))
    length(tag_values) == 2 || throw(
        ArgumentError(
            "mixing: tags must take exactly two distinct values; got $(length(tag_values))",
        ),
    )
    t1, t2 = tag_values

    r = Vector{Float64}(undef, N)
    @inbounds for i in 1:N
        r[i] =
            sqrt((pos[1, i] - centre[1])^2 + (pos[2, i] - centre[2])^2 + (pos[3, i] - centre[3])^2)
    end
    A = aplus(r[tags .== t1], r[tags .== t2], r_h; r_max = 2 * r_h)

    r_cut = 2 * r_h
    inside = r .≤ r_cut
    sel1 = findall(i -> inside[i] && tags[i] == t1, 1:N)
    sel2 = findall(i -> inside[i] && tags[i] == t2, 1:N)
    n1 = length(sel1)
    n2 = length(sel2)

    d_radius = two_sample_d(r[sel1], r[sel2]).d
    d_energy = energy === nothing ? NaN : two_sample_d(energy[sel1], energy[sel2]).d

    M = 0.0
    v̄x = v̄y = v̄z = 0.0
    @inbounds for i in 1:N
        inside[i] || continue
        M += mass[i]
        v̄x += mass[i] * vel[1, i]
        v̄y += mass[i] * vel[2, i]
        v̄z += mass[i] * vel[3, i]
    end
    d_angular_momentum = NaN
    if M > 0
        v̄x /= M
        v̄y /= M
        v̄z /= M
        L = Vector{Float64}(undef, N)
        @inbounds for i in 1:N
            inside[i] || continue
            x = pos[1, i] - centre[1]
            y = pos[2, i] - centre[2]
            z = pos[3, i] - centre[3]
            vx = vel[1, i] - v̄x
            vy = vel[2, i] - v̄y
            vz = vel[3, i] - v̄z
            L[i] = sqrt((y * vz - z * vy)^2 + (z * vx - x * vz)^2 + (x * vy - y * vx)^2)
        end
        d_angular_momentum = two_sample_d(L[sel1], L[sel2]).d
    end

    return (;
        aplus2 = A.aplus,
        aplus2_error = A.error,
        d_radius,
        d_energy,
        d_angular_momentum,
        floor = _ks_floor(n1, n2),
        n = (n1, n2),
    )
end
