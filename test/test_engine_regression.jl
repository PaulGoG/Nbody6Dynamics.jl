# Engine regression suite on the CPU build (roadmap §4.2, item B6): the integrator's convergence
# order and momentum drift, an isolated binary on both KS paths, a Kozai–Lidov triple on the KS
# and merger paths, the start-up step of a cold system and the Pythagorean problem through the
# chain, each against its closed form or reference trajectory. The reference values and their
# tolerances come from the full measurements (longer spans, three seeds) behind this suite.
#
# Included by runtests.jl when NBODY6_BINARY_TESTS = 1. The engine is the AVX build under
# <NBODY6_BACKEND_ROOT>/backend/Nbody6PPGPU-beijing, or the one test_engine_gated.jl built
# (ENGINE_BUILD_BASE). About 30 s of engine time on two threads.

using StableRNGs: StableRNG
using LinearAlgebra: norm

include(joinpath(@__DIR__, "engine_regression_setup.jl"))

# Reference values from the full measurements. η ladder, Hermite scheme
# alone, N = 1024 to t = 2, ETAI/ETAR = (0.01, 0.02) f for f = 1/4, 1/2, 1. The order of the
# median interval error is 2 in theory; the full ladder (t = 4, three seeds) gives 2.03 ± 0.40
# over these rungs and 1.96–2.32 seed by seed, this configuration 2.15 ± 0.64: the band is the
# measured range widened by the fit's standard error. The bounds at f = 1 are twice to three times
# the largest value of the three seeds (energy 1.4×10⁻⁷, momentum 1.7×10⁻⁷, angular momentum
# 1.7×10⁻⁷).
const B6_ORDER_BAND = (1.4, 2.8)
const B6_EPS_BASELINE_MAX = 3.0e-7
const B6_DP_BASELINE_MAX = 5.0e-7
const B6_DL_BASELINE_MAX = 5.0e-7
# Kozai–Lidov triple, direct three-body integration (dimensionless units G = m = a_in = 1).
const B6_KOZAI = (e_max = 0.9107, t_max = 9052.0, period = 18110.0, e_min = 0.0990)
# Pythagorean problem at t = 100: the binary of masses 4 and 5 (BigFloat Vern9, 10⁻¹⁸ and 10⁻²²).
const B6_PYTHAGOREAN = (a = 0.5524949, e = 0.9887100)

"""The CPU engine binary of the backend the gated tests use, or `nothing`."""
function _b6_engine()
    base =
        isdefined(Main, :ENGINE_BUILD_BASE) ? Main.ENGINE_BUILD_BASE :
        get(ENV, "NBODY6_BACKEND_ROOT", "")
    isempty(base) && return nothing
    src = joinpath(base, "backend", "Nbody6PPGPU-beijing")
    isdir(joinpath(src, "build")) || return nothing
    try
        return Nbody6Dynamics._find_binary(
            src,
            "nbody6++",
            BuildConfig(; enable_gpu = false, enable_mpi = false),
        )
    catch err
        err isa ErrorException || rethrow()
        @info "engine regression suite skipped: no CPU build of nbody6++ under $src"
        return nothing
    end
end

@testset "Engine regression suite (binary-gated)" begin
    engine = _b6_engine()
    if engine === nothing
        @info "engine regression suite skipped: no backend (set NBODY6_BACKEND_ROOT)"
        @test_skip false
    else
        ER = EngineRegression
        work = mktempdir()
        t_suite = time()

        # 1 + 5. η ladder, Hermite scheme alone (no KS regularisation): the median interval error
        # follows η² (fourth-order Hermite with Aarseth's step, Δt ∝ η^{1/2}); momentum and angular
        # momentum drift with it (Ahmad–Cohen forces are not pairwise symmetric).
        @testset "η ladder: convergence order and momentum drift" begin
            mass, pos, vel = ER.plummer_model(1024, StableRNG(20261104))
            fs = [0.25, 0.5, 1.0]
            m = map(fs) do f
                dir = joinpath(work, "eta_f$(f)")
                ER.eta_ladder_input(dir, mass, pos, vel, f; span = 2.0, mode = :hermite)
                r = ER.run_engine(engine, dir; wall = 300)
                @test r.exit_code == 0 && r.ended && !r.halted
                ER.ladder_metrics(dir)
            end
            ε = [x.median_de for x in m]
            fit = ER.fit_order(fs, ε)
            @info "B6 η ladder: median |ΔE_k|/|E₀| $(ε), p = $(round(fit.p; digits = 2)) ± $(round(fit.σp; digits = 2)); ΔP/(Mσ) $([x.dP for x in m]); ΔL/(M r_h σ) $([x.dL for x in m])"
            @test all(x -> x.nksreg == 0, m)
            @test issorted(ε)
            @test B6_ORDER_BAND[1] ≤ fit.p ≤ B6_ORDER_BAND[2]
            @test ε[end] < B6_EPS_BASELINE_MAX
            @test m[end].dP < B6_DP_BASELINE_MAX && m[end].dL < B6_DL_BASELINE_MAX
            @test m[1].dP < m[end].dP && m[1].dL < m[end].dL
        end

        # 2. Isolated binary, e = 0.9, 10³ orbits: a and e (double-precision listings) are conserved
        # to round-off on both paths, the Runge–Lenz direction to the single-precision snapshots, and
        # the perturbed KS integration keeps the Kepler phase (the Stumpff KS step is exact for an
        # unperturbed Kepler orbit).
        @testset "Kepler binary, $(leg) KS path" for leg in (:perturbed, :unperturbed)
            mass, pos, vel, sys = ER.planted_binary(StableRNG(23))
            dir = joinpath(work, "binary_$(leg)")
            ER.binary_input(dir, mass, pos, vel, sys; leg)
            r = ER.run_engine(engine, dir; wall = 300)
            @test r.exit_code == 0 && r.ended && !r.halted
            ps = ER.pair_series(dir, 1, 2)
            keep = ps.t .≤ sys.span + 1e-12
            @test ps.missing_listings == 0 && count(keep) ≥ 64
            a = -sys.m^2 ./ (2 .* ps.eb[keep])
            @test maximum(abs, a ./ sys.a .- 1) < 1e-12
            @test maximum(abs, ps.e[keep] .- sys.e) < 1e-12
            snaps = filter(
                s -> Nbody6Dynamics.time_nb(s.header) ≤ sys.span + 1e-9,
                ER.snapshot_series(dir),
            )
            ω = map(snaps) do s
                Δx, Δv = ER.relative_state(s, 1, 2)
                el = ER.orbital_elements(sys.m, sys.m, Δx, Δv)
                atan(el.evec[2], el.evec[1])
            end
            @test maximum(abs, ω) < 1e-5
            ks = ER.step_counts(dir).ks
            if leg === :perturbed
                @test 40 < ks / sys.orbits < 100
                # Phase: mean anomaly from the position angle against n t, t on the output grid. A
                # snapshot taken while a background KS pair exists holds this pair's components at
                # an earlier time (7 of 65 in the full measurement), hence the median.
                dt = sys.span / 64
                δM = map(snaps) do s
                    t = floor(Nbody6Dynamics.time_nb(s.header) / dt) * dt
                    Δx, Δv = ER.relative_state(s, 1, 2)
                    el = ER.orbital_elements(sys.m, sys.m, Δx, Δv)
                    ν = atan(Δx[2], Δx[1]) - atan(el.evec[2], el.evec[1])
                    E = 2atan(sqrt((1 - sys.e) / (1 + sys.e)) * tan(ν / 2))
                    abs(rem2pi(E - sys.e * sin(E) - 2π * t / sys.P, RoundNearest)) / 2π
                end
                @test ER.median(δM) < 1e-5
                @test count(<(1e-4), δM) ≥ 0.8 * length(δM)
            else
                @test ks < sys.orbits / 10
            end
        end

        # 4. Kozai–Lidov triple (a_out/a_in = 15, i₀ = 75°, e₀ = 0.1). KS path: the reference
        # trajectory (direct three-body integration) reaches e_max = 0.9107 at t = 9052 and 27162,
        # e_min = 0.099 between. Merger path: the inner binary of a merger is a ghost whose
        # eccentricity does not evolve (no ECCMOD without KZ(27) = 2 and stellar radii).
        @testset "Kozai–Lidov triple" begin
            mass, pos, vel, sys = ER.kozai_triple(StableRNG(31))
            dir = joinpath(work, "kozai_ks")
            interval = 2.0^-8
            ER.kozai_input(dir, mass, pos, vel, sys; path = :ks, interval)
            r = ER.run_engine(engine, dir; wall = 300)
            @test r.exit_code == 0 && r.ended && !r.halted
            ps = ER.pair_series(dir, 1, 2)
            t = ps.t ./ sys.τ
            keep = t .≤ sys.t_units + 1
            t, e = t[keep], ps.e[keep]
            h = searchsortedfirst(t, sys.t_units / 2)
            k1, k2 = argmax(e[1:(h - 1)]), h - 1 + argmax(e[h:end])
            Δt = interval / sys.τ
            @test abs(e[k1] - B6_KOZAI.e_max) < 0.01 && abs(e[k2] - B6_KOZAI.e_max) < 0.01
            @test abs((t[k2] - t[k1]) / B6_KOZAI.period - 1) < 0.03
            @test abs(minimum(e[k1:k2]) - B6_KOZAI.e_min) < 0.01
            @test abs(t[k1] - B6_KOZAI.t_max) < 2Δt
            # The pair keeps its perturber from the start, so the first interval holds as many KS
            # steps as the later ones (before engine v2026.07.26 the largest body mass was reset at
            # the first adjustment, every perturber search of that interval found nothing, the
            # first interval held ~34 KS steps against ~2300 thereafter and e(t) lagged the
            # reference by one interval).
            sc = Int[]
            for l in eachline(joinpath(dir, "out.log"))
                k = findfirst("NSTEP(I,B,R,U)", l)
                k === nothing || push!(sc, parse(Int, split(l[(last(k) + 1):end])[4]))
            end
            @test sc[2] ≥ sc[3] / 2

            mdir = joinpath(work, "kozai_merger")
            msys = merge(sys, (; span = 2.0^-5))
            ER.kozai_input(mdir, mass, pos, vel, msys; path = :merger, interval)
            r = ER.run_engine(engine, mdir; wall = 300)
            @test r.exit_code == 0 && r.ended && !r.halted
            log = read(joinpath(mdir, "out.log"), String)
            @test count("NEW MERGER", log) == 1 && count("NEW CHAIN", log) == 0
            ein = Float64[]
            for f in filter(startswith("bdat.9_"), readdir(mdir))
                parse(Float64, f[8:end]) ≤ msys.span + 1e-12 || continue
                L = redirect_stderr(devnull) do
                    read_ks_pair_listing(joinpath(mdir, f))
                end
                k = findfirst(r -> r.name_cm < 0, L.records)
                k === nothing || push!(ein, L.records[k].eccentricity)
            end
            @test length(ein) ≥ 6 && maximum(ein) - minimum(ein) < 1e-8
            # over the same span the KS path has raised e by more than 0.01
            @test maximum(e[t .≤ msys.span / sys.τ]) - minimum(ein) > 0.01
        end

        # Start-up step of a cold system (the Pythagorean triple at rest, Hermite scheme only):
        # with the full start-up derivatives (KZ(40) = 1) and with the force-only start (KZ(40) = 0,
        # bounded by the free-fall scale of the neighbour sphere since engine v2026.07.26) the
        # triple's energy holds to 10⁻⁴ until the first close approach (before, a body at rest
        # under a constant force got the largest block step and the energy was off by 6 % at
        # t = 1.3).
        @testset "Start-up step of a cold triple, KZ(40) = $kz40" for kz40 in (1, 0)
            mass, pos, vel, sys = ER.pythagorean_triple(StableRNG(17))
            dir = joinpath(work, "cold_kz40_$(kz40)")
            ER.planted_input(
                dir,
                mass,
                pos,
                vel;
                span = 2.0^-6,
                interval = 2.0^-6,
                deltat = 2.0^-9,
                n_background = length(mass) - 3,
                rmin = 1.0e-8,
                kz = Dict(15 => 0, 40 => kz40),
            )
            r = ER.run_engine(engine, dir; wall = 120)
            @test r.exit_code == 0 && r.ended && !r.halted
            E0 = -769 / 60
            dE = [
                abs(ER.triple_energy(s, sys) / E0 - 1) for
                s in ER.snapshot_series(dir) if Nbody6Dynamics.time_nb(s.header) / sys.τ ≤ 1.3
            ]
            @test length(dE) ≥ 4
            @test maximum(dE) < 1e-4
        end

        # 3. Pythagorean problem (Szebehely & Peters 1967): mass 3 escapes, masses 4 and 5 form a
        # binary with a = 0.5525, e = 0.9887 (BigFloat reference, converged to 10⁻⁸; an initial
        # displacement of 10⁻⁶ already changes the escaper). The chain forms at t ≈ 0.16 (units of
        # the problem), returns eleven times and the engine then spins in its block loop with the
        # time frozen (engine v2026.07.26; before it, the chain was terminated at its small-step limit at t ≈ 0.37 and
        # restarted, after which the triple's centre of mass moved at ~10³ N-body velocity units
        # without an energy halt, and the internal motion ended in the published class with other
        # elements, a = 0.68, e = 0.67). Engine defect, roadmap 9.41: the completion and the
        # outcome are marked broken.
        @testset "Pythagorean problem through the chain" begin
            mass, pos, vel, sys = ER.pythagorean_triple(StableRNG(17))
            dir = joinpath(work, "pythagorean")
            ER.pythagorean_input(dir, mass, pos, vel, sys)
            r = ER.run_engine(engine, dir; wall = 60)
            snaps = ER.snapshot_series(dir)
            done =
                !r.timed_out &&
                !isempty(snaps) &&
                Nbody6Dynamics.time_nb(snaps[end].header) ≥ sys.span
            @test_broken done
            if done
                o = ER.triple_outcome(snaps[end], sys)
                @test o.escaper == 3 && o.pair == (4.0, 5.0) && o.E3 > 0
                @test_broken ER.triple_cm_speed(snaps[end], sys) < 1e-2
                @test_broken ER.energy_record(dir).abs_detot < 1e-4
                @test_broken abs(o.a / B6_PYTHAGOREAN.a - 1) < 0.05 &&
                             abs(o.e - B6_PYTHAGOREAN.e) < 0.005
            end
        end
        @info "B6 engine regression suite: $(round(time() - t_suite; digits = 1)) s"
    end
end
