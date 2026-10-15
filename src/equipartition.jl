# Degree of energy equipartition: the Bianchini et al. (2016) equipartition-mass fit of
# the velocity-dispersion–mass relation of a star sample, its snapshot series and CSV writer.

"""
    EquipartitionFit

Result of [`equipartition_fit`](@ref): the equipartition mass `m_eq` [M☉] (Inf when the
dispersion is consistent with independence of mass), the Δχ² = 1 interval `m_eq_low`
(NaN when it reaches the lowest grid value) and `m_eq_high` (Inf when it reaches the highest),
the amplitude `sigma0` in the units of the velocities, `chi2` at the minimum, `chi2_flat` of the
mass-independent model, the number of mass bins used `n_bins` and stars in them `n_stars`, and
the bins: `bin_mass` (mean mass [M☉]), `bin_sigma` (one-dimensional dispersion), `bin_sigma_err`
(its error), `bin_n` (members).
"""
struct EquipartitionFit
    m_eq::Float64
    m_eq_low::Float64
    m_eq_high::Float64
    sigma0::Float64
    chi2::Float64
    chi2_flat::Float64
    n_bins::Int
    n_stars::Int
    bin_mass::Vector{Float64}
    bin_sigma::Vector{Float64}
    bin_sigma_err::Vector{Float64}
    bin_n::Vector{Int}
end

"""
Mass dependence `g(m; m_eq)` of the dispersion, Bianchini et al. (2016), eq. 1:
`exp(−m / (2 m_eq))` for `m ≤ m_eq`, `exp(−1/2) (m / m_eq)^{−1/2}` above.
"""
_equipartition_g(m::Float64, m_eq::Float64)::Float64 =
    m <= m_eq ? exp(-m / (2 * m_eq)) : exp(-0.5) * sqrt(m_eq / m)

"""
Weighted least-squares amplitude and χ² of `σ_b ≈ σ₀ g_b` with weights `w_b`.
"""
function _equipartition_chi2(
    σ::Vector{Float64},
    w::Vector{Float64},
    g::Vector{Float64},
)::NTuple{2,Float64}
    num = 0.0
    den = 0.0
    for b in eachindex(σ, w, g)
        num += w[b] * σ[b] * g[b]
        den += w[b] * g[b]^2
    end
    σ₀ = num / den
    χ² = 0.0
    for b in eachindex(σ, w, g)
        χ² += w[b] * (σ[b] - σ₀ * g[b])^2
    end
    return (σ₀, χ²)
end

"""
    equipartition_fit(mass::AbstractVector{<:Real}, vel::AbstractMatrix{<:Real};
                      bin_width::Real = 0.1, n_min::Integer = 50, n_grid::Integer = 240,
                      m_eq_max_factor::Real = 100.0) -> EquipartitionFit

Fit the equipartition mass of Bianchini et al. (2016), MNRAS 458, 3644,
DOI 10.1093/mnras/stw552, eq. 1, to the velocity dispersion of a star sample binned in mass.

# Arguments
- `mass`: stellar masses [M☉], length n, each positive and finite.
- `vel`: 3 × n velocities in any one unit; the sample selection (luminous members within a
  radius) is the caller's.
- `bin_width`: width of the mass bins [M☉]; the first left edge is
  `floor(minimum(mass) / bin_width) * bin_width`, and a star belongs to the bin with
  `edge_j ≤ m < edge_{j+1}`, the most massive star to the last bin.
- `n_min`: least number of members of a bin entering the fit (≥ 2).
- `n_grid`: number of points of the logarithmic `m_eq` grid (≥ 2).
- `m_eq_max_factor`: upper end of the grid in units of the mean mass of the heaviest bin (> 1);
  the lower end is half the mean mass of the lightest bin.

The dispersions are taken about the mass-weighted mean velocity of all given stars; a bin of
`n_b` members has the one-dimensional dispersion `σ_b = sqrt(Σ |vᵢ − v̄|² / (3 n_b))` with error
`σ_b / sqrt(6 n_b)`. The model `σ(m) = σ₀ g(m; m_eq)` is evaluated on the grid, `σ₀` by weighted
linear least squares at each grid value; the grid minimum of χ² gives `m_eq` (Inf at the last
grid point), and the contiguous run of grid values with χ² within 1 of the minimum gives the
interval. Fewer than three kept bins leave the fit `NaN` (`chi2_flat` needs two).

# Returns
An [`EquipartitionFit`](@ref).

Throws a `DimensionMismatch` unless `size(vel) == (3, length(mass))`, and an `ArgumentError`
for an empty sample, a non-positive or non-finite mass, `bin_width ≤ 0`, `n_min < 2`,
`n_grid < 2` or `m_eq_max_factor ≤ 1`.

# Example
```julia
g(m, m_eq) = m <= m_eq ? exp(-m / (2m_eq)) : exp(-0.5) * sqrt(m_eq / m)
mass = 0.1 .+ 1.9 .* rand(20_000)
vel = 10 .* g.(mass, 0.8)' .* randn(3, 20_000)   # km/s
fit = equipartition_fit(mass, vel)
(fit.m_eq, fit.m_eq_low, fit.m_eq_high)          # m_eq ≈ 0.8 M☉
```
"""
function equipartition_fit(
    mass::AbstractVector{<:Real},
    vel::AbstractMatrix{<:Real};
    bin_width::Real = 0.1,
    n_min::Integer = 50,
    n_grid::Integer = 240,
    m_eq_max_factor::Real = 100.0,
)::EquipartitionFit
    n = length(mass)
    size(vel) == (3, n) ||
        throw(DimensionMismatch("vel must be 3 × $n to match mass, got $(size(vel))"))
    n >= 1 || throw(ArgumentError("equipartition_fit needs at least one star"))
    bin_width > 0 || throw(ArgumentError("bin_width must be positive, got $bin_width"))
    n_min >= 2 || throw(ArgumentError("n_min must be ≥ 2, got $n_min"))
    n_grid >= 2 || throw(ArgumentError("n_grid must be ≥ 2, got $n_grid"))
    m_eq_max_factor > 1 || throw(ArgumentError("m_eq_max_factor must be > 1, got $m_eq_max_factor"))
    for m in mass
        (isfinite(m) && m > 0) ||
            throw(ArgumentError("every mass must be positive and finite, got $m"))
    end

    Δm = Float64(bin_width)
    r1, r2, r3 = axes(vel, 1)

    # Mass-weighted mean velocity of the whole sample.
    m_tot = 0.0
    v̄x = 0.0
    v̄y = 0.0
    v̄z = 0.0
    for (i, j) in zip(eachindex(mass), axes(vel, 2))
        mi = Float64(mass[i])
        m_tot += mi
        v̄x += mi * Float64(vel[r1, j])
        v̄y += mi * Float64(vel[r2, j])
        v̄z += mi * Float64(vel[r3, j])
    end
    v̄x /= m_tot
    v̄y /= m_tot
    v̄z /= m_tot

    # Mass bins: member counts, mass sums and squared-deviation sums.
    m_min = Float64(minimum(mass))
    m_max = Float64(maximum(mass))
    edge0 = floor(m_min / Δm) * Δm
    n_all = max(1, ceil(Int, (m_max - edge0) / Δm))
    bin_count = zeros(Int, n_all)
    sum_m = zeros(Float64, n_all)
    sum_dv² = zeros(Float64, n_all)
    for (i, j) in zip(eachindex(mass), axes(vel, 2))
        mi = Float64(mass[i])
        b = clamp(floor(Int, (mi - edge0) / Δm) + 1, 1, n_all)
        bin_count[b] += 1
        sum_m[b] += mi
        sum_dv²[b] +=
            (Float64(vel[r1, j]) - v̄x)^2 +
            (Float64(vel[r2, j]) - v̄y)^2 +
            (Float64(vel[r3, j]) - v̄z)^2
    end

    # Kept bins, in the order of increasing mass.
    kept = findall(c -> c >= n_min, bin_count)
    n_bins = length(kept)
    bin_n = bin_count[kept]
    bin_mass = sum_m[kept] ./ bin_n
    bin_sigma = sqrt.(sum_dv²[kept] ./ (3 .* bin_n))
    bin_sigma_err = bin_sigma ./ sqrt.(6 .* bin_n)
    n_stars = sum(bin_n; init = 0)

    w = 1 ./ bin_sigma_err .^ 2
    chi2_flat = n_bins >= 2 ? _equipartition_chi2(bin_sigma, w, ones(n_bins))[2] : NaN
    if n_bins < 3
        return EquipartitionFit(
            NaN,
            NaN,
            NaN,
            NaN,
            NaN,
            chi2_flat,
            n_bins,
            n_stars,
            bin_mass,
            bin_sigma,
            bin_sigma_err,
            bin_n,
        )
    end

    # χ² on the logarithmic m_eq grid.
    K = Int(n_grid)
    lo = bin_mass[1] / 2
    hi = Float64(m_eq_max_factor) * bin_mass[end]
    grid = exp.(range(log(lo), log(hi); length = K))
    chi2_grid = Vector{Float64}(undef, K)
    sigma0_grid = Vector{Float64}(undef, K)
    g = Vector{Float64}(undef, n_bins)
    for k in 1:K
        for b in 1:n_bins
            g[b] = _equipartition_g(bin_mass[b], grid[k])
        end
        sigma0_grid[k], chi2_grid[k] = _equipartition_chi2(bin_sigma, w, g)
    end

    k_min = argmin(chi2_grid)
    chi2 = chi2_grid[k_min]
    sigma0 = sigma0_grid[k_min]
    m_eq = k_min == K ? Inf : grid[k_min]

    # Contiguous Δχ² = 1 run around the minimum.
    k_lo = k_min
    while k_lo > 1 && chi2_grid[k_lo - 1] <= chi2 + 1
        k_lo -= 1
    end
    k_hi = k_min
    while k_hi < K && chi2_grid[k_hi + 1] <= chi2 + 1
        k_hi += 1
    end
    m_eq_low = k_lo == 1 ? NaN : grid[k_lo]
    m_eq_high = k_hi == K ? Inf : grid[k_hi]

    return EquipartitionFit(
        m_eq,
        m_eq_low,
        m_eq_high,
        sigma0,
        chi2,
        chi2_flat,
        n_bins,
        n_stars,
        bin_mass,
        bin_sigma,
        bin_sigma_err,
        bin_n,
    )
end

"""
    EquipartitionRow

One snapshot of the equipartition series: `time_nb`, `time_myr`, `n_stars` and `n_bins` used,
`m_eq`, `m_eq_low`, `m_eq_high` [M☉], `sigma0_kms` [km s⁻¹], `chi2`, `chi2_flat`, and the mass
range of the fitted stars `mass_min`, `mass_max` [M☉].
"""
struct EquipartitionRow
    time_nb::Float64
    time_myr::Float64
    n_stars::Int
    n_bins::Int
    m_eq::Float64
    m_eq_low::Float64
    m_eq_high::Float64
    sigma0_kms::Float64
    chi2::Float64
    chi2_flat::Float64
    mass_min::Float64
    mass_max::Float64
end

"""
    write_equipartition(path::AbstractString, rows::AbstractVector{EquipartitionRow}) -> String

Write the equipartition series as CSV, one line per snapshot.

# Arguments
- `path`: output file; an existing file is backed up, never overwritten.
- `rows`: the [`EquipartitionRow`](@ref)s, written in the given order.

Columns follow the fields of [`EquipartitionRow`](@ref). Non-finite values are empty cells,
except an infinite `m_eq` or `m_eq_high`, written as `Inf` (no equipartition within the grid).

# Returns
The path written.

# Example
```julia
fit = equipartition_fit(mass, vel_kms)
row = EquipartitionRow(t_nb, t_myr, fit.n_stars, fit.n_bins, fit.m_eq, fit.m_eq_low,
                       fit.m_eq_high, fit.sigma0, fit.chi2, fit.chi2_flat,
                       minimum(mass), maximum(mass))
write_equipartition("equipartition.csv", [row])
```
"""
function write_equipartition(path::AbstractString, rows::AbstractVector{EquipartitionRow})::String
    _backup_existing(path)
    open(path, "w") do io
        println(
            io,
            "time_nb,time_myr,n_stars,n_bins,m_eq,m_eq_low,m_eq_high,sigma0_kms,chi2,chi2_flat," *
            "mass_min,mass_max",
        )
        for r in rows
            cells = (
                _csv_precise(r.time_nb),
                _csv_precise(r.time_myr),
                string(r.n_stars),
                string(r.n_bins),
                isinf(r.m_eq) ? "Inf" : _csv_precise(r.m_eq),
                _csv_precise(r.m_eq_low),
                isinf(r.m_eq_high) ? "Inf" : _csv_precise(r.m_eq_high),
                _csv_precise(r.sigma0_kms),
                _csv_precise(r.chi2),
                _csv_precise(r.chi2_flat),
                _csv_precise(r.mass_min),
                _csv_precise(r.mass_max),
            )
            println(io, join(cells, ","))
        end
    end
    return String(path)
end

"""
    equipartition_series(snaps, bound_sets, units; sevs = StellarEvolutionSnapshot[],
                         settings = PostprocessConfig()) -> Vector{EquipartitionRow}

Equipartition mass of every snapshot from the luminous bound stars within the half-mass
radius: the bound set of the snapshot (`bound_sets[k].stars`) about its shrinking-sphere
density centre (centre keys of `settings`), the members inside `r_h`, and, when a
stellar-evolution snapshot lies within half the snapshot spacing of the snapshot time (nearest
in `time_nb`), those with `K* ≤ 9` (main sequence to AGB) by name; without such a snapshot all
members enter. Masses are converted with `units.zmbar` to M☉ and velocities with
`units.vstar` to km s⁻¹ before [`equipartition_fit`](@ref) with `settings.equipartition_bin_width`
and `settings.equipartition_n_min`. A snapshot whose selection holds fewer than
`2 settings.equipartition_n_min` stars gets a row of `NaN` fit values.

# Example
```julia
rows = equipartition_series(snaps, [s.bound for s in ints], extract_scaling(diag); sevs)
write_equipartition("equipartition.csv", rows)
```
"""
function equipartition_series(
    snaps::Vector{Snapshot},
    bound_sets::AbstractVector{BoundSet},
    units::UnitScaling;
    sevs::AbstractVector{StellarEvolutionSnapshot} = StellarEvolutionSnapshot[],
    settings::PostprocessConfig = PostprocessConfig(),
)::Vector{EquipartitionRow}
    length(bound_sets) == length(snaps) || throw(
        ArgumentError(
            "equipartition_series: $(length(bound_sets)) bound sets for $(length(snaps)) snapshots",
        ),
    )
    keywords = _centre_keywords(settings)
    n_min = settings.equipartition_n_min
    times = [time_nb(s.header) for s in snaps]
    spacing = length(times) ≥ 2 ? minimum(diff(times)) : Inf
    sev_times = [isempty(v.records) ? NaN : v.records[1].time_nb for v in sevs]
    rows = EquipartitionRow[]
    for (k, snap) in enumerate(snaps)
        t = times[k]
        t_myr = time_myr(snap.header)
        b = bound_sets[k].stars
        nan_row = EquipartitionRow(t, t_myr, 0, 0, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN)
        if length(b) < 2 * n_min
            push!(rows, nan_row)
            continue
        end
        pos_b = Float64.(snap.pos[:, b])
        mass_b = Float64.(snap.mass[b])
        c = shrinking_sphere_centre(pos_b, mass_b; keywords...).centre
        r_h = _lagrangian_radii(pos_b, mass_b, c, (0.5,))[1]
        inside = [
            i for i in eachindex(b) if
            (pos_b[1, i] - c[1])^2 + (pos_b[2, i] - c[2])^2 + (pos_b[3, i] - c[3])^2 ≤ r_h^2
        ]
        # Luminous stars by the nearest stellar-evolution snapshot
        j = isempty(sevs) ? 0 : argmin(abs.(sev_times .- t))
        if j > 0 && abs(sev_times[j] - t) ≤ 0.5 * spacing
            luminous = Set{Int32}(r.name for r in sevs[j].records if r.stellar_type ≤ 9)
            filter!(i -> snap.name[b[i]] in luminous, inside)
        end
        if length(inside) < 2 * n_min
            push!(rows, nan_row)
            continue
        end
        mass_solar = mass_b[inside] .* units.zmbar
        vel_kms = Float64.(snap.vel[:, b[inside]]) .* units.vstar
        fit = equipartition_fit(
            mass_solar,
            vel_kms;
            bin_width = settings.equipartition_bin_width,
            n_min = n_min,
        )
        push!(
            rows,
            EquipartitionRow(
                t,
                t_myr,
                fit.n_stars,
                fit.n_bins,
                fit.m_eq,
                fit.m_eq_low,
                fit.m_eq_high,
                fit.sigma0,
                fit.chi2,
                fit.chi2_flat,
                minimum(mass_solar),
                maximum(mass_solar),
            ),
        )
    end
    return rows
end
