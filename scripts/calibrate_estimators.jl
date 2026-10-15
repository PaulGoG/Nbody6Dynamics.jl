# Bias and variance of the structural and kinematic estimators versus N, measured on
# static models with known answers. Writes the calibration table the validity record
# of a run reads (assets/calibration/estimators.csv) and its sidecar
# (assets/calibration/estimators.toml).
#
#   julia -t 8 scripts/calibrate_estimators.jl [--quick] [--n=1000,3162,...] [--reps=8]
#                                               [--models=plummer,plummer_om15] [--out=<dir>]
#
# Every realisation is seeded from (model, N, replica) with StableRNGs, so the table
# regenerates identically on any machine and Julia version. `--quick` runs two small
# N with three replicas each (a schema and timing check, not a calibration).
include(joinpath(@__DIR__, "activate.jl"))

using Nbody6Dynamics
using StableRNGs, Statistics, LinearAlgebra, Printf, Dates, TOML

const CALIBRATION_VERSION = 1
const N_GRID = [1_000, 3_162, 10_000, 31_623, 100_000, 316_228, 1_000_000]
const REPLICAS = Dict(
    1_000 => 64,
    3_162 => 48,
    10_000 => 32,
    31_623 => 16,
    100_000 => 12,
    316_228 => 6,
    1_000_000 => 4,
)
const MODELS = [
    "plummer",            # isotropic Plummer, equal masses: every null and every known profile
    "plummer_om15",       # Osipkov–Merritt Plummer, r_a = 1.5 a: anisotropy profile
    "plummer_rot",        # Plummer with every v_φ prograde (Lynden-Bell p = 1): rotation signal
    "plummer_kroupa",     # Plummer with independent Kroupa masses: segregation nulls, mass weighting
    "plummer_kroupa_s05", # the same with Šubr degree S = 0.5: segregation signal
    "plummer_stretch",    # Plummer stretched to axis ratios q = 0.8, s = 0.6: shape signal
]
const LAGR_FRACTIONS = (0.01, 0.05, 0.1, 0.5, 0.9)
const N_TRACER = 50
const R_A = 1.5

# Plummer relations, G = M = a = 1
const R_H_PLUMMER = 1 / sqrt(2^(2 / 3) - 1)
const E_PLUMMER = -3π / 64
const K_PLUMMER = 3π / 64
plummer_density(r) = 3 / (4π) * (1 + r^2)^(-2.5)
plummer_potential(r) = -1 / sqrt(1 + r^2)
plummer_sigma_r(r) = sqrt(1 / (6 * sqrt(1 + r^2)))
plummer_lagrangian_radius(f) = 1 / sqrt(f^(-2 / 3) - 1)
om_beta(r) = r^2 / (r^2 + R_A^2)

"""Seed of one realisation; replicas are numbered below 100."""
function realisation_seed(model::AbstractString, N::Integer, replica::Integer)
    replica < 100 || throw(ArgumentError("replica index must be below 100"))
    i = findfirst(==(model), MODELS)
    i === nothing && throw(ArgumentError("unknown model $model"))
    return UInt64(100_000_000 * i + 100 * N + replica)
end

specific_energy(pos, vel) = [
    0.5 * (vel[1, i]^2 + vel[2, i]^2 + vel[3, i]^2) +
    plummer_potential(sqrt(pos[1, i]^2 + pos[2, i]^2 + pos[3, i]^2)) for i in axes(pos, 2)
]

function build_model(model::AbstractString, N::Int, rng)
    if model == "plummer_om15"
        pos, vel = sample_plummer_om(N, 1.0, R_A; rng)
    else
        pos, vel = Nbody6Dynamics.sample_plummer(N, 1.0; rng)
    end
    mass = fill(1.0 / N, N)
    if model == "plummer_rot"
        lynden_bell_flip!(vel, pos, [0.0, 0.0, 1.0], 1.0; rng)
    elseif model == "plummer_kroupa" || model == "plummer_kroupa_s05"
        mass = Nbody6Dynamics.sample_kroupa(N; m_low = 0.08, m_up = 100.0, rng)
        mass ./= sum(mass)
        model == "plummer_kroupa_s05" &&
            rank_masses_by_energy!(mass, specific_energy(pos, vel), 0.5; rng)
    elseif model == "plummer_stretch"
        stretch_positions!(pos, 0.8, 0.6)
    end
    return pos, vel, mass
end

"""Known answers and which estimators apply, per model."""
const MODEL_TRUTH = Dict(
    "plummer" => (; q = 1.0, s = 1.0, beta = r -> 0.0, isotropic = true, segregation = NaN),
    "plummer_om15" =>
        (; q = 1.0, s = 1.0, beta = om_beta, isotropic = false, segregation = NaN),
    "plummer_rot" => (; q = 1.0, s = 1.0, beta = r -> 0.0, isotropic = true, segregation = NaN),
    "plummer_kroupa" =>
        (; q = 1.0, s = 1.0, beta = r -> 0.0, isotropic = true, segregation = 1.0),
    "plummer_kroupa_s05" =>
        (; q = 1.0, s = 1.0, beta = r -> 0.0, isotropic = true, segregation = NaN),
    "plummer_stretch" =>
        (; q = 0.8, s = 0.6, beta = r -> NaN, isotropic = false, segregation = NaN),
)

"""
One realisation: a dictionary (estimator, quantity) → (value, truth, reported error,
predicted scatter). Quantities without a known answer carry `truth = NaN`.
"""
function evaluate(model::AbstractString, N::Int, replica::Int)
    rng = StableRNG(realisation_seed(model, N, replica))
    pos, vel, mass = build_model(model, N, rng)
    truth = MODEL_TRUTH[model]
    out = Dict{Tuple{String,String},NTuple{4,Float64}}()
    record!(est, q, v, t = NaN, e = NaN, p = NaN) = (out[(est, q)] = (Float64(v), t, e, p))
    origin = zeros(3)
    kroupa = startswith(model, "plummer_kroupa")
    stretched = model == "plummer_stretch"

    # Centre: the fixed point of the iteration scatters by σ_axis (a/R_f)² in a flat core
    # (design note §7.3); the expected |offset| of a 3-D isotropic error is 2√(2/π) σ_axis.
    c = shrinking_sphere_centre(pos, mass)
    chi3 = 2 * sqrt(2 / π)
    record!(
        "centre",
        "offset",
        norm(c.centre),
        0.0,
        chi3 * c.sigma,
        chi3 * c.sigma / min(1.0, c.radius^2),
    )
    record!("centre", "radius", c.radius)
    record!("centre", "n_final", c.n_final)
    if kroupa
        cm = shrinking_sphere_centre(pos, mass; weights = :mass)
        record!(
            "centre",
            "offset_mass_weighted",
            norm(cm.centre),
            0.0,
            chi3 * cm.sigma,
            chi3 * cm.sigma / min(1.0, cm.radius^2),
        )
    end

    if !stretched
        # Core radius (Casertano–Hut): √0.3 a with exact densities
        ρn = local_density(pos, mass; k = 6, weighting = :number)
        cn = core_radius(pos, ρn)
        record!("core_radius", "r_core_number", cn.r_core, sqrt(0.3))
        record!("core_radius", "centre_offset_number", norm(cn.centre), 0.0)
        ρm = local_density(pos, mass; k = 6, weighting = :mass)
        cmr = core_radius(pos, ρm)
        record!("core_radius", "r_core_mass", cmr.r_core, sqrt(0.3))
        record!("core_radius", "centre_offset_mass", norm(cmr.centre), 0.0)
        ρe = [plummer_density(norm(@view pos[:, i])) for i in 1:N]
        record!("core_radius", "r_core_exact_density", core_radius(pos, ρe).r_core, sqrt(0.3))
        record!("core_radius", "rho_max_number", maximum(ρn), NaN)

        # Lagrangian radii about the true centre and about the estimated one
        for (label, centre) in (("true_centre", origin), ("estimated_centre", c.centre))
            g = lagrangian_grid(pos, mass, centre)
            for f in LAGR_FRACTIONS
                i = findfirst(x -> isapprox(x, f), g.fractions)
                (i === nothing || g.masked[i]) && continue
                record!(
                    "lagrangian_grid",
                    "ln_r_$(f)_$(label)",
                    log(g.radii[i]),
                    log(plummer_lagrangian_radius(f)),
                    g.sigma_log[i],
                    g.sigma_log[i],
                )
            end
        end

        # Profiles about the true centre; residuals against the model at the shell radii
        p = radial_profile(pos, vel, mass, origin)
        ok = .!isnan.(p.beta)
        if any(ok)
            record!(
                "radial_profile",
                "beta_residual",
                mean(p.beta[ok] .- truth.beta.(p.r[ok])),
                0.0,
                mean(p.beta_err[ok]),
            )
            record!(
                "radial_profile",
                "beta_residual_rms",
                sqrt(mean((p.beta[ok] .- truth.beta.(p.r[ok])) .^ 2)),
                NaN,
                mean(p.beta_err[ok]),
            )
        end
        oks = .!isnan.(p.sigma_r)
        if truth.isotropic && any(oks)
            record!(
                "radial_profile",
                "sigma_r_relative",
                mean(p.sigma_r[oks] ./ plummer_sigma_r.(p.r[oks]) .- 1),
                0.0,
                mean(p.sigma_err[oks]),
            )
        end
        okd = .!isnan.(p.rho_log) .& (p.rho_log .> 0)
        any(okd) && record!(
            "radial_profile",
            "rho_relative",
            mean(p.rho_log[okd] ./ plummer_density.(p.r_log[okd]) .- 1),
            0.0,
        )
        record!("radial_profile", "r_h", p.r_h, R_H_PLUMMER)
        record!("radial_profile", "n_shells", count(ok))

        # Relaxation and crossing times from the sample
        t = relaxation_times(pos, vel, mass, origin)
        record!(
            "relaxation_times",
            "t_rh_gamma0.02",
            t.t_rh[1],
            0.138 * N * R_H_PLUMMER^1.5 / log(0.02 * N),
        )
        record!("relaxation_times", "t_cr", t.t_cr, 2 * R_H_PLUMMER / sqrt(2 * K_PLUMMER))

        # Rotation: the Plummer sample is the null (truth 0), the flipped model the signal
        rot =
            rotation_analysis(pos, vel, mass; energy = E_PLUMMER, r_h = R_H_PLUMMER, seed = replica)
        signal = model == "plummer_rot"
        record!(
            "rotation",
            "lambda_r",
            rot.lambda_r,
            signal ? NaN : 0.0,
            rot.lambda_r_floor,
            rot.lambda_r_floor,
        )
        record!(
            "rotation",
            "lambda_peebles",
            rot.lambda_peebles,
            signal ? NaN : 0.0,
            rot.lambda_peebles_floor,
            rot.lambda_peebles_floor_analytic,
        )
        record!("rotation", "lambda_peebles_floor_analytic", rot.lambda_peebles_floor_analytic)
        record!("rotation", "v_peak", rot.v_peak, signal ? NaN : 0.0)
        record!("rotation", "r_peak", rot.r_peak)
        record!("rotation", "n_shells", length(rot.profile.radius))
        okr = .!isnan.(rot.profile.v_rot_debiased)
        any(okr) && record!(
            "rotation",
            "v_rot_debiased_max",
            maximum(rot.profile.v_rot_debiased[okr]),
            signal ? NaN : 0.0,
            mean(rot.profile.v_rot_floor[okr]),
        )

        # Pairs: a Plummer sample has no binaries, so every pair found is spurious
        kT = 2 * K_PLUMMER / (3 * N)
        prs = snapshot_pairs(pos, vel, mass; kT)
        record!("snapshot_pairs", "n_pairs", length(prs), 0.0)
        record!("snapshot_pairs", "n_tidal_ok", count(q -> q.tidal_ok, prs), 0.0)
        record!(
            "snapshot_pairs",
            "n_hard_tidal_ok",
            count(q -> q.tidal_ok && q.hardness > 1, prs),
            0.0,
        )
        record!(
            "snapshot_pairs",
            "hardness_max",
            isempty(prs) ? 0.0 : maximum(q -> q.hardness, prs),
            0.0,
        )
        record!(
            "snapshot_pairs",
            "hardness_max_tidal_ok",
            any(q -> q.tidal_ok, prs) ? maximum(q -> q.hardness, filter(q -> q.tidal_ok, prs)) :
            0.0,
            0.0,
        )

        # Mixing null: two random halves of one model
        tags = rand(rng, 1:2, N)
        mx = mixing(pos, vel, mass, tags, origin, R_H_PLUMMER; energy = specific_energy(pos, vel))
        record!("mixing", "aplus2", mx.aplus2, 0.0, mx.aplus2_error)
        record!("mixing", "d_radius", mx.d_radius, NaN, NaN, mx.floor)
        record!("mixing", "d_energy", mx.d_energy, NaN, NaN, mx.floor)
        record!("mixing", "d_angular_momentum", mx.d_angular_momentum, NaN, NaN, mx.floor)
        record!("mixing", "floor", mx.floor)
    end

    if kroupa
        tracer = collect(Int, partialsortperm(mass, 1:N_TRACER; rev = true))
        ms = mass_segregation(pos, mass; tracer, n_mst = (10, 20, 50), seed = replica)
        record!("mass_segregation", "lambda_msr", ms.lambda_msr, truth.segregation, ms.lambda_err)
        record!(
            "mass_segregation",
            "lambda_median",
            ms.lambda_median,
            truth.segregation,
            ms.lambda_median_err,
        )
        record!(
            "mass_segregation",
            "ln_gamma_msr",
            log(ms.gamma_msr),
            isnan(truth.segregation) ? NaN : 0.0,
            log(ms.gamma_err),
        )
        record!("mass_segregation", "r_half_ratio", ms.r_half_ratio, truth.segregation)
        for (k, n_mst) in enumerate(ms.n_mst_used)
            record!(
                "mass_segregation",
                "lambda_msr_nmst$(n_mst)",
                ms.lambda_scan[k],
                truth.segregation,
                ms.lambda_err_scan[k],
            )
        end
        r = [norm(@view pos[:, i]) for i in 1:N]
        isT = falses(N)
        isT[tracer] .= true
        ap = aplus(r[isT], r[.!isT], R_H_PLUMMER)
        record!("aplus", "aplus_tracer50", ap.aplus, isnan(truth.segregation) ? NaN : 0.0, ap.error)
    end

    # Shape: null (truth 1) on the spherical models, signal on the stretched one
    sh = shape_tensor(pos, origin)
    for (k, s) in enumerate(sh)
        record!("shape_tensor", "q_shell$k", s.q, truth.q)
        record!("shape_tensor", "s_shell$k", s.s, truth.s)
        record!("shape_tensor", "n_shell$k", s.n)
    end
    pe = projected_ellipticity(pos, origin, R_H_PLUMMER)
    for (k, f) in enumerate((0.25, 0.5, 1.0))
        record!("projected_ellipticity", "epsilon_xy_$(f)", pe[1, k], stretched ? NaN : 0.0)
    end
    return out
end

"""Aggregate the replicas of one (model, N) into table rows."""
function aggregate(model, N, results)
    keys_all = sort!(unique(k for r in results for k in keys(r)))
    rows = NamedTuple[]
    for key in keys_all
        vals = [r[key] for r in results if haskey(r, key) && !isnan(r[key][1])]
        isempty(vals) && continue
        v = [x[1] for x in vals]
        t = [x[2] for x in vals]
        e = [x[3] for x in vals]
        p = [x[4] for x in vals]
        n = length(v)
        m = mean(v)
        sc = n ≥ 2 ? std(v) : NaN
        tr = all(isnan, t) ? NaN : mean(filter(!isnan, t))
        push!(
            rows,
            (;
                model,
                estimator = key[1],
                quantity = key[2],
                N,
                n_rep = n,
                truth = tr,
                mean = m,
                scatter = sc,
                bias = isnan(tr) ? NaN : m - tr,
                bias_err = isnan(tr) || isnan(sc) ? NaN : sc / sqrt(n),
                reported_error = all(isnan, e) ? NaN : mean(filter(!isnan, e)),
                predicted_scatter = all(isnan, p) ? NaN : mean(filter(!isnan, p)),
            ),
        )
    end
    return rows
end

fmt(x::Real) = isnan(x) ? "NaN" : @sprintf("%.6g", x)
fmt(x::Integer) = string(x)

function write_table(path, rows)
    open(path, "w") do io
        println(
            io,
            "model,estimator,quantity,N,n_rep,truth,mean,scatter,bias,bias_err,reported_error,predicted_scatter",
        )
        for r in rows
            println(
                io,
                join(
                    [
                        r.model,
                        r.estimator,
                        r.quantity,
                        string(r.N),
                        string(r.n_rep),
                        fmt(r.truth),
                        fmt(r.mean),
                        fmt(r.scatter),
                        fmt(r.bias),
                        fmt(r.bias_err),
                        fmt(r.reported_error),
                        fmt(r.predicted_scatter),
                    ],
                    ",",
                ),
            )
        end
    end
end

function parse_args(args)
    quick = "--quick" in args
    grid = quick ? [1_000, 3_162] : copy(N_GRID)
    reps = quick ? Dict(n => 3 for n in grid) : copy(REPLICAS)
    models = copy(MODELS)
    out = joinpath(@__DIR__, "..", "assets", "calibration")
    for a in args
        if startswith(a, "--n=")
            grid = parse.(Int, split(a[5:end], ","))
        elseif startswith(a, "--reps=")
            k = parse(Int, a[8:end])
            reps = Dict(n => k for n in grid)
        elseif startswith(a, "--models=")
            models = String.(split(a[10:end], ","))
            all(m -> m in MODELS, models) ||
                throw(ArgumentError("unknown model in $(models); known: $(MODELS)"))
        elseif startswith(a, "--out=")
            out = a[7:end]
        elseif a != "--quick"
            throw(ArgumentError("unknown argument $a"))
        end
    end
    for n in grid
        haskey(reps, n) || (reps[n] = 3)
    end
    return (; quick, grid, reps, models, out)
end

function main(args)
    opts = parse_args(args)
    mkpath(opts.out)
    rows = NamedTuple[]
    timings = Dict{String,Float64}()
    t_start = now()
    for model in opts.models, N in opts.grid
        n_rep = opts.reps[N]
        results = Vector{Dict{Tuple{String,String},NTuple{4,Float64}}}(undef, n_rep)
        t0 = time()
        Threads.@threads :dynamic for replica in 1:n_rep
            results[replica] = evaluate(model, N, replica)
        end
        dt = time() - t0
        timings["$(model)/N$(N)"] = dt
        append!(rows, aggregate(model, N, results))
        @printf(
            "%-20s N = %8d  replicas = %3d  %8.1f s  (%6.2f s per replica, %d threads)\n",
            model,
            N,
            n_rep,
            dt,
            dt / n_rep,
            Threads.nthreads()
        )
        flush(stdout)
    end
    table = joinpath(opts.out, "estimators.csv")
    write_table(table, rows)
    meta = Dict(
        "version" => CALIBRATION_VERSION,
        "generated" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
        "quick" => opts.quick,
        "julia" => string(VERSION),
        "package" => string(pkgversion(Nbody6Dynamics)),
        "host" => gethostname(),
        "threads" => Threads.nthreads(),
        "models" => opts.models,
        "n_grid" => opts.grid,
        "replicas" => Dict(string(n) => opts.reps[n] for n in opts.grid),
        "seed" => "StableRNG(100_000_000 i_model + 100 N + replica)",
        "wall_seconds" => round(Dates.value(now() - t_start) / 1000; digits = 1),
        "timings_seconds" => timings,
        "columns" => Dict(
            "truth" => "known answer of the static model; NaN without one",
            "mean" => "mean of the estimate over the replicas",
            "scatter" => "standard deviation over the replicas",
            "bias" => "mean − truth",
            "bias_err" => "scatter / √n_rep",
            "reported_error" => "mean of the estimator's own error estimate",
            "predicted_scatter" => "mean of the analytic noise prediction",
        ),
    )
    open(joinpath(opts.out, "estimators.toml"), "w") do io
        TOML.print(io, meta)
    end
    println("wrote $(length(rows)) rows to $table")
end

main(ARGS)
