# Input builders, engine launcher and output readers of the engine regression suite
# (test_engine_regression.jl). Every configuration is written in N-body units (G = 1) and read
# by the engine without rescaling (KZ(22) = 2), stellar evolution off (KZ(19) = 0), escapers kept
# (KZ(23) = 0), so that the integrals of the input are the integrals the engine conserves.
module EngineRegression

using Nbody6Dynamics
using LinearAlgebra: norm, dot, cross
using Printf: @sprintf
using Random: AbstractRNG

const NB = Nbody6Dynamics

# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------

"""
    plummer_model(N, rng) -> (mass, pos, vel)

Equal-mass Plummer sphere of `N` bodies in Hénon units: total mass 1, centre of mass at rest at
the origin, virial ratio exactly 1/2 and total energy exactly −1/4 (direct pair sums), so that
the virial radius is 1. `pos` and `vel` are 3 × N.
"""
function plummer_model(N::Integer, rng::AbstractRNG)
    N ≥ 2 || throw(ArgumentError("plummer_model: N must be ≥ 2, got $N"))
    pos, vel = NB.sample_plummer(Int(N), 3π / 16; rng = rng)
    mass = fill(1.0 / N, N)
    _to_centre_of_mass!(mass, pos, vel)
    T, W = kinetic_energy(mass, vel), potential_energy(mass, pos)
    vel .*= sqrt(-0.5W / T)
    s = -2W                      # E = W/2 → −1/4 after x → s x, v → v/√s
    pos .*= s
    vel ./= sqrt(s)
    return mass, pos, vel
end

function _to_centre_of_mass!(mass, pos, vel)
    M = sum(mass)
    for k in 1:3
        pos[k, :] .-= dot(mass, pos[k, :]) / M
        vel[k, :] .-= dot(mass, vel[k, :]) / M
    end
    return nothing
end

"""Kinetic energy ½ Σ mᵢ vᵢ² of a 3 × N velocity array."""
kinetic_energy(mass, vel) = 0.5 * sum(mass[i] * sum(abs2, @view vel[:, i]) for i in eachindex(mass))

"""Potential energy −Σᵢ<ⱼ mᵢ mⱼ / rᵢⱼ (G = 1) by direct summation."""
function potential_energy(mass, pos)
    W = 0.0
    n = length(mass)
    @inbounds for i in 1:(n - 1), j in (i + 1):n
        r = sqrt((pos[1, i] - pos[1, j])^2 + (pos[2, i] - pos[2, j])^2 + (pos[3, i] - pos[3, j])^2)
        W -= mass[i] * mass[j] / r
    end
    return W
end

"""
    kepler_pair(m1, m2, a, e; phase = :pericentre) -> (Δx, Δv)

Relative position x₂ − x₁ and velocity of a bound Kepler orbit (G = 1) in the x–y plane with the
pericentre on +x and the orbital angular momentum along +z, at pericentre or apocentre.
"""
function kepler_pair(m1, m2, a, e; phase::Symbol = :pericentre)
    (a > 0 && 0 ≤ e < 1) || throw(ArgumentError("kepler_pair: needs a > 0 and 0 ≤ e < 1"))
    M = m1 + m2
    if phase === :pericentre
        r = a * (1 - e)
        return [r, 0.0, 0.0], [0.0, sqrt(M * (1 + e) / r), 0.0]
    elseif phase === :apocentre
        r = a * (1 + e)
        return [-r, 0.0, 0.0], [0.0, -sqrt(M * (1 - e) / r), 0.0]
    end
    throw(ArgumentError("kepler_pair: phase must be :pericentre or :apocentre"))
end

"""
    orbital_elements(m1, m2, Δx, Δv) -> NamedTuple

Two-body elements of the relative orbit (G = 1): semi-major axis `a`, eccentricity `e`, specific
energy `h`, specific angular momentum vector `j`, eccentricity (Runge–Lenz) vector `evec` =
Δv × j / M − Δx / |Δx|, and the period `P` (NaN when unbound).
"""
function orbital_elements(m1, m2, Δx, Δv)
    M = m1 + m2
    r = norm(Δx)
    h = 0.5 * sum(abs2, Δv) - M / r
    j = cross(Δx, Δv)
    evec = cross(Δv, j) ./ M .- Δx ./ r
    a = -M / (2h)
    P = h < 0 ? 2π * sqrt(a^3 / M) : NaN
    return (; a, e = norm(evec), h, j, evec, P)
end

"""
    background_cluster(N, rng; offset) -> (mass, pos, vel)

The equal-mass Plummer model of [`plummer_model`](@ref) translated by `offset` (N-body lengths):
the host cluster of a planted few-body system, which keeps the engine's neighbour and
time-step machinery in its normal regime while the planted system at the origin stays isolated.
"""
function background_cluster(N::Integer, rng::AbstractRNG; offset = [100.0, 0.0, 0.0])
    mass, pos, vel = plummer_model(N, rng)
    pos .+= offset
    return mass, pos, vel
end

# ---------------------------------------------------------------------------
# Engine input and launch
# ---------------------------------------------------------------------------

"""
    write_engine_input(dir, mass, pos, vel; tcrit, dtadj, deltat, etai, etar, etau, rmin, dtmin,
                       nbin0 = 0, gmin = 1e-6, gmax = 0.01, qe = 0.02, nnbopt, rs0,
                       eclose = 1.0, kz = Dict(), rbar = 1.0, zmbar = 1.0) -> String

Write `dat.10` and `merger.inp` for a raw engine run in `dir` and return the input path. The
bodies are read as given (`KZ(22) = 2`, no rescaling); the first `2 nbin0` bodies are regularised
at the start as pairs (`KZ(8) = 2`). Stellar evolution is off (`KZ(19) = 0`), escapers are kept
(`KZ(23) = 0`), the few-body options are those of `kz` over the generator's inherited block (no
chain, no hierarchies). `tcrit`, `dtadj`, `deltat` must be multiples of 2⁻ᵏ with k ≤ 9
([`engine_interval`](@ref)).
"""
function write_engine_input(
    dir::AbstractString,
    mass::Vector{Float64},
    pos::Matrix{Float64},
    vel::Matrix{Float64};
    tcrit::Real,
    dtadj::Real,
    deltat::Real,
    etai::Real,
    etar::Real,
    etau::Real = 0.1,
    rmin::Real,
    dtmin::Real,
    nbin0::Integer = 0,
    gmin::Real = 1.0e-6,
    gmax::Real = 0.01,
    qe::Real = 0.02,
    nnbopt::Integer = 100,
    rs0::Real = 0.3,
    eclose::Real = 1.0,
    kz::Dict{Int,Int} = Dict{Int,Int}(),
    rbar::Real = 1.0,
    zmbar::Real = 1.0,
)
    for (label, dt) in (("tcrit", tcrit), ("dtadj", dtadj), ("deltat", deltat))
        dt > 0 && isinteger(dt * 2^9) || throw(
            ArgumentError("write_engine_input: $label = $dt is not a positive multiple of 2⁻⁹"),
        )
    end
    mkpath(dir)
    N = length(mass)
    spec = Nbody6ParameterSpec(;
        qe = Float64(qe),
        etai = Float64(etai),
        etar = Float64(etar),
        etau = Float64(etau),
        nnbopt = Int(nnbopt),
        rs0 = Float64(rs0),
        rmin = Float64(rmin),
        dtmin = Float64(dtmin),
        eclose = Float64(eclose),
        gmin = Float64(gmin),
        gmax = Float64(gmax),
        preset = "baseline",
        escapers = "keep",
        kz = kz,
    )
    stellar = StellarSpec(; kz19 = 0, dtplot = Float64(deltat))
    path = joinpath(dir, "merger.inp")
    Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        NB.write_dat10(joinpath(dir, "dat.10"), mass, pos, vel)
        return NB.generate_merger_inp(
            path,
            N,
            Float64(rbar),
            Float64(zmbar);
            nbody6 = spec,
            stellar = stellar,
            tcrit = Float64(tcrit),
            dtadj = Float64(dtadj),
            deltat = Float64(deltat),
            nbin0 = Int(nbin0),
            mass_bounds = extrema(mass) .* (zmbar * N),
        )
    end
    return path
end

"""
    run_engine(engine, dir; threads = 2, wall = 600.0) -> NamedTuple

Run `engine < merger.inp > out.log` in `dir` with an unlimited stack, `OMP_NUM_THREADS =
threads` and unbuffered Fortran output, ending it after `wall` seconds. Returns the exit code
(−1 on time-out), the terminating signal, the wall time and whether it timed out.
"""
function run_engine(
    engine::AbstractString,
    dir::AbstractString;
    threads::Integer = 2,
    wall::Real = 600.0,
)
    isfile(engine) || throw(ArgumentError("run_engine: no engine binary at $engine"))
    script = "ulimit -s unlimited && exec $(abspath(engine)) < merger.inp > out.log 2> err.log"
    cmd = addenv(
        Cmd(`sh -c $script`; dir = dir),
        "OMP_STACKSIZE" => "4096M",
        "OMP_NUM_THREADS" => string(threads),
        "GFORTRAN_UNBUFFERED_PRECONNECTED" => "y",
    )
    t0 = time()
    p = run(cmd; wait = false)
    while process_running(p) && time() - t0 < wall
        sleep(0.2)
    end
    timed_out = process_running(p)
    if timed_out
        kill(p)
        wait(p)
    end
    return (;
        exit_code = timed_out ? -1 : p.exitcode,
        signal = p.termsignal,
        wall_s = time() - t0,
        timed_out,
    )
end

# ---------------------------------------------------------------------------
# Output readers
# ---------------------------------------------------------------------------

"""
    energy_record(dir) -> NamedTuple

Energy bookkeeping of a run from its `ADJUST` lines: times, the relative error per interval
`de` (ΔE / max(T, |E|)), the absolute change per interval `delta`, the running sum `detot`, the
initial total energy `e0`, and the summary errors `sum_abs_de` = Σ|ΔE| / |E₀| and `abs_detot`
= |ΣΔE| / |E₀|.
"""
function energy_record(dir::AbstractString)
    adj = read_diagnostics(joinpath(dir, "out.log")).adjust
    isempty(adj) && error("energy_record: no ADJUST line in $(joinpath(dir, "out.log"))")
    e0 = adj[1].e_tot
    rows = adj[2:end]
    delta = [a.de_abs for a in rows]
    return (;
        time = [a.time_nb for a in rows],
        de = [a.de_rel for a in rows],
        delta,
        detot = [a.detot for a in rows],
        e0,
        npairs = [a.npairs for a in rows],
        sum_abs_de = sum(abs, delta; init = 0.0) / abs(e0),
        abs_detot = isempty(rows) ? 0.0 : abs(rows[end].detot) / abs(e0),
    )
end

"""
    snapshot_series(dir) -> Vector{Snapshot}

The `conf.3_*` snapshots of a run in time order.
"""
snapshot_series(dir::AbstractString) = redirect_stderr(devnull) do
    read_all_conf3(dir; threaded = false)
end

"""
    momentum_integrals(s::Snapshot) -> NamedTuple

Total mass `M`, linear momentum `P`, angular momentum `L` about the origin, the mass-weighted
velocity dispersion `σ` (about the centre-of-mass velocity) and the median distance `r_h` of
the bodies from the centre of mass (the half-mass radius of an equal-mass model), in double
precision from the single-precision snapshot.
"""
function momentum_integrals(s)
    m = Float64.(s.mass)
    x = Float64.(s.pos)
    v = Float64.(s.vel)
    M = sum(m)
    P = [dot(m, @view v[k, :]) for k in 1:3]
    L = zeros(3)
    for i in eachindex(m)
        L .+= m[i] .* cross(x[:, i], v[:, i])
    end
    vcm = P ./ M
    σ = sqrt(sum(m[i] * sum(abs2, v[:, i] .- vcm) for i in eachindex(m)) / M)
    xcm = [dot(m, @view x[k, :]) / M for k in 1:3]
    r = sort([norm(x[:, i] .- xcm) for i in eachindex(m)])
    return (; t = NB.time_nb(s.header), M, P, L, σ, r_h = r[cld(length(r), 2)])
end

"""
    fit_order(η, err) -> NamedTuple

Least-squares fit log err = p log η + c; returns the order `p`, its standard error `σp` and `c`.
"""
function fit_order(η::AbstractVector, err::AbstractVector)
    length(η) == length(err) ≥ 3 ||
        throw(ArgumentError("fit_order: needs ≥ 3 points of equal length"))
    x = log.(η)
    y = log.(err)
    n = length(x)
    x̄, ȳ = sum(x) / n, sum(y) / n
    Sxx = sum(abs2, x .- x̄)
    p = sum((x .- x̄) .* (y .- ȳ)) / Sxx
    c = ȳ - p * x̄
    s² = sum(abs2, y .- (p .* x .+ c)) / (n - 2)
    return (; p, σp = sqrt(s² / Sxx), c)
end

"""Median of a non-empty real vector."""
function median(x::AbstractVector{<:Real})
    isempty(x) && throw(ArgumentError("median: empty vector"))
    y = sort(x)
    n = length(y)
    return isodd(n) ? Float64(y[cld(n, 2)]) : 0.5 * (y[n ÷ 2] + y[n ÷ 2 + 1])
end

"""
    ks_counters(dir) -> NamedTuple

The cumulative KS and few-body counters of the last `#3` line of `out.log`: KS attempts
`nkstry`, regularisations `nksreg`, hyperbolic pairs `nkshyp`, unperturbed two-body orbits
`nksper`, slow-down changes `nksmod`, and the triple, quadruple, chain and merger counts.
"""
function ks_counters(dir::AbstractString)
    v = Int[]
    for l in eachline(joinpath(dir, "out.log"))
        startswith(l, " #3 ") && (v = parse.(Int, split(l)[2:end]))
    end
    length(v) ≥ 10 || error("ks_counters: no #3 line in $(joinpath(dir, "out.log"))")
    return (;
        nkstry = v[1],
        nksreg = v[2],
        nkshyp = v[3],
        nksper = v[4],
        nksmod = v[5],
        nttry = v[6],
        ntrip = v[7],
        nquad = v[8],
        nchain = v[9],
        nmerg = v[10],
    )
end

"""
    step_counts(dir) -> NamedTuple

Irregular, block, regular and perturbed-KS step counts summed over the output intervals of
`out.log` (the `NSTEP(I,B,R,U)` field of the `TIME[NB]` lines, which counts per interval).
"""
function step_counts(dir::AbstractString)
    c = zeros(Int, 4)
    for l in eachline(joinpath(dir, "out.log"))
        k = findfirst("NSTEP(I,B,R,U)", l)
        k === nothing && continue
        v = split(l[(last(k) + 1):end])
        c .+= parse.(Int, v[1:4])
    end
    return (; irregular = c[1], block = c[2], regular = c[3], ks = c[4])
end

# ---------------------------------------------------------------------------
# Components
# ---------------------------------------------------------------------------

"""
Options of every B6 input on top of the generator's block: `KZ(40) = 1` obtains the second and
third force derivatives at the start (`FPOLY2`) and sets the first steps by the full Aarseth
criterion. Under `KZ(40) = 0` the first steps come from `TSTEP2`, whose force-only branch
(fork v2026.07.12) gives a body without relative motion a step 0.25 η |F|/|Ḟ|, unbounded by its
dynamical time: a cold few-body system moves under its initial forces until that step ends.
"""
const STARTUP_KZ = Dict(40 => 1)

"""`(etai, etar)` of the `"baseline"` preset, the centre of the η ladder."""
const ETA_BASELINE = (etai = 0.01, etar = 0.02)

"""Half-mass radius of the Plummer model in Hénon units, (3π/16)/√(2^{2/3} − 1)."""
const PLUMMER_RH = (3π / 16) / sqrt(2^(2 / 3) - 1)

"""
    eta_ladder_input(dir, mass, pos, vel, f; span = 4.0, mode = :hermite, dtadj = 0.25) -> String

Input of one rung of the η ladder: `ETAI = 0.01 f`, `ETAR = 0.02 f` (the `"baseline"` preset
scaled by `f`), `ETAU = 0.1`, energy check and snapshot every `dtadj`, run to `span`. `mode =
:hermite` sets `RMIN = 10⁻⁶`, `DTMIN = 10⁻¹²`, so that no close encounter is regularised and the
whole system follows the Hermite scheme; `mode = :preset` uses the `"baseline"` rules,
`RMIN = 12 r_h/N` and `DTMIN = 0.016 √(ETAI/0.02) √(RMIN³ N)`, and regularises close encounters by KS.
"""
function eta_ladder_input(
    dir,
    mass,
    pos,
    vel,
    f::Real;
    span = 4.0,
    mode::Symbol = :hermite,
    dtadj = 0.25,
)
    f > 0 || throw(ArgumentError("eta_ladder_input: the ladder factor must be > 0, got $f"))
    etai, etar = ETA_BASELINE.etai * f, ETA_BASELINE.etar * f
    N = length(mass)
    rmin, dtmin = if mode === :hermite
        (1.0e-6, 1.0e-12)
    elseif mode === :preset
        r = 12 * PLUMMER_RH / N
        (r, 0.016 * sqrt(etai / 0.02) * sqrt(r^3 * N))
    else
        throw(ArgumentError("eta_ladder_input: mode must be :hermite or :preset, got $(repr(mode))"))
    end
    return write_engine_input(
        dir,
        mass,
        pos,
        vel;
        tcrit = span,
        dtadj,
        deltat = dtadj,
        etai,
        etar,
        etau = 0.1,
        rmin,
        dtmin,
        nnbopt = 50,
        rs0 = 0.4,
        kz = copy(STARTUP_KZ),
    )
end

"""
    ladder_metrics(dir) -> NamedTuple

Error measures of one η-ladder run: the median over the adjustment intervals of |ΔE_k|/|E₀|
(`median_de`, insensitive to the few intervals that hold a close encounter), Σ|ΔE_k|/|E₀|,
|ΣΔE_k|/|E₀|, the largest drift of the total linear momentum over the snapshots relative to M σ
(`dP`) and of the total angular momentum relative to M r_h σ (`dL`), with σ the initial
mass-weighted velocity dispersion and r_h the initial median radius.
"""
function ladder_metrics(dir::AbstractString)
    er = energy_record(dir)
    de = abs.(er.delta) ./ abs(er.e0)
    I = momentum_integrals.(snapshot_series(dir))
    i0 = I[1]
    dP = maximum(norm(x.P .- i0.P) for x in I) / (i0.M * i0.σ)
    dL = maximum(norm(x.L .- i0.L) for x in I) / (i0.M * i0.r_h * i0.σ)
    return (;
        median_de = median(de),
        sum_abs_de = er.sum_abs_de,
        abs_detot = er.abs_detot,
        dP,
        dL,
        n_intervals = length(de),
        n_snapshots = length(I),
        nksreg = ks_counters(dir).nksreg,
    )
end

"""
    planted_binary(rng; m = 1e-3, e = 0.9, orbits = 1000, span = 0.125, n_background = 256)
        -> (mass, pos, vel, sys)

An equal-mass binary (m each, eccentricity `e`) at pericentre at the origin, whose period is
`span / orbits`, as the first two bodies, moving at a small common velocity (a body at rest
would take its first step from the force criterion alone), followed by a background Plummer
model 100 length units away. `sys` holds `m`, `a`, `e`, `P`, `span`, `orbits`.
"""
function planted_binary(
    rng::AbstractRNG;
    m = 1.0e-3,
    e = 0.9,
    orbits = 1000,
    span = 0.125,
    n_background = 256,
)
    P = span / orbits
    a = cbrt((P / 2π)^2 * 2m)
    dx, dv = kepler_pair(m, m, a, e)
    vcm = [2.0e-3, 1.0e-3, -1.0e-3]
    mb, xb, vb = background_cluster(n_background, rng)
    mass = vcat([m, m], mb)
    pos = hcat(-0.5 .* dx, 0.5 .* dx, xb)
    vel = hcat(vcm .- 0.5 .* dv, vcm .+ 0.5 .* dv, vb)
    return mass, pos, vel, (; m, a, e, P, span, orbits)
end

"""
    planted_input(dir, mass, pos, vel; span, interval, deltat = interval, nbin0, n_background, kz,
                  gmin = 1e-6) -> String

Input of a planted few-body system run to `span` with an energy check every `interval` and a
snapshot every `deltat` (snapshots at energy checks hold the KS components resolved at the
snapshot time; between checks those of a perturbed pair date from its last KS step), `RMIN` and
`DTMIN` of the background model, `TCRIT` past `span` by the 20 `DTMIN` at which `ADJUST` ends
the run early, and the start-up options of [`STARTUP_KZ`](@ref) under `kz`.
"""
function planted_input(
    dir,
    mass,
    pos,
    vel;
    span,
    interval,
    deltat = interval,
    nbin0 = 0,
    n_background,
    kz = Dict{Int,Int}(),
    gmin = 1.0e-6,
    rmin = 12 * PLUMMER_RH / n_background,
    dtmin = 0.016 * sqrt(ETA_BASELINE.etai / 0.02) * sqrt(rmin^3 * length(mass)),
)
    tcrit = span + interval * ceil(20dtmin / interval + 1)
    return write_engine_input(
        dir,
        mass,
        pos,
        vel;
        tcrit,
        dtadj = interval,
        deltat,
        etai = ETA_BASELINE.etai,
        etar = ETA_BASELINE.etar,
        rmin,
        dtmin,
        nbin0,
        gmin,
        nnbopt = 30,
        rs0 = 0.4,
        kz = merge(STARTUP_KZ, kz),
    )
end

"""
    binary_input(dir, mass, pos, vel, sys; leg) -> String

Input of the planted binary of [`planted_binary`](@ref) on one KS path: `leg = :perturbed`
(`GMIN = 10⁻³⁰`, `KZ(26) = 0`: every step a perturbed KS step) or `:unperturbed` (`GMIN = 10⁻⁶`,
`KZ(26) = 1`: two-body motion of the isolated pair). The pair is regularised at the start
(`NBIN0 = 1`, `KZ(8) = 2`), with 64 energy checks and snapshots over the span.
"""
function binary_input(dir, mass, pos, vel, sys; leg::Symbol)
    gmin, kz26 = if leg === :perturbed
        (1.0e-30, 0)
    elseif leg === :unperturbed
        (1.0e-6, 1)
    else
        throw(ArgumentError("binary_input: leg must be :perturbed or :unperturbed, got $(repr(leg))"))
    end
    return planted_input(
        dir,
        mass,
        pos,
        vel;
        span = sys.span,
        interval = sys.span / 64,
        nbin0 = 1,
        n_background = length(mass) - 2,
        gmin,
        kz = Dict(26 => kz26),
    )
end

"""
    pair_series(dir, name1, name2) -> NamedTuple

Time series of the KS pair `(name1, name2)` from the `bdat.9_<t>` listings (`KZ(8) ≥ 2`): times
from the file names, binding energy μh (`eb`), eccentricity `e`, and the c.m. name; listings
without the pair are skipped (`missing_listings` counts them).
"""
function pair_series(dir::AbstractString, name1::Integer, name2::Integer)
    files = filter(startswith("bdat.9_"), readdir(dir))
    sort!(files; by = f -> parse(Float64, f[8:end]))
    t = Float64[]
    eb = Float64[]
    e = Float64[]
    cm = Int[]
    missing_listings = 0
    for f in files
        L = redirect_stderr(devnull) do
            read_ks_pair_listing(joinpath(dir, f))
        end
        k = findfirst(r -> minmax(r.name1, r.name2) == minmax(name1, name2), L.records)
        if k === nothing
            missing_listings += 1
            continue
        end
        r = L.records[k]
        push!(t, parse(Float64, f[8:end]))
        push!(eb, r.binding_energy_nb)
        push!(e, r.eccentricity)
        push!(cm, r.name_cm)
    end
    return (; t, eb, e, cm, missing_listings, n_listings = length(files))
end

"""
    relative_state(s, name1, name2) -> (Δx, Δv)

Relative position and velocity x₂ − x₁, v₂ − v₁ of two named bodies of a snapshot, in double
precision from the single-precision snapshot.
"""
function relative_state(s, name1::Integer, name2::Integer)
    i = findfirst(==(Int32(name1)), s.name)
    j = findfirst(==(Int32(name2)), s.name)
    (i === nothing || j === nothing) &&
        throw(ArgumentError("relative_state: body $name1 or $name2 not in the snapshot"))
    return Float64.(s.pos[:, j]) .- Float64.(s.pos[:, i]),
    Float64.(s.vel[:, j]) .- Float64.(s.vel[:, i])
end

"""
    kozai_triple(rng; m = 1e-3, e_in = 0.1, ratio = 15.0, inclination = 75.0, span = 1.0,
                 t_units = 40000.0, n_background = 256) -> (mass, pos, vel, sys)

Hierarchical triple of three bodies of mass `m`: the inner binary (bodies 1, 2) at pericentre
with a = 1 and eccentricity `e_in` in the x–y plane, the outer body (3) on a circular orbit of
radius `ratio` about the inner centre of mass, its orbital plane inclined by `inclination`
degrees about the x axis (the inner line of apsides lies on the line of nodes, ω = 0).
Dimensionless units G = m = a_in = 1 map onto N-body units by the time unit τ = span/t_units and
the length unit λ = (τ² m)^{1/3}; the background Plummer model lies 100 length units away.
`sys` holds `m`, `τ`, `λ` and the elements.
"""
function kozai_triple(
    rng::AbstractRNG;
    m = 1.0e-3,
    e_in = 0.1,
    ratio = 15.0,
    inclination = 75.0,
    span = 1.0,
    t_units = 40000.0,
    n_background = 256,
)
    τ = span / t_units
    λ = cbrt(τ^2 * m)
    ι = deg2rad(inclination)
    dx, dv = kepler_pair(1.0, 1.0, 1.0, e_in)
    R = [ratio, 0.0, 0.0]
    V = [0.0, cos(ι), sin(ι)] .* sqrt(3 / ratio)
    vcm = [2.0e-3, 1.0e-3, -1.0e-3]
    xs = hcat(-0.5 .* dx .- R ./ 3, 0.5 .* dx .- R ./ 3, (2 / 3) .* R) .* λ
    vs = hcat(-0.5 .* dv .- V ./ 3, 0.5 .* dv .- V ./ 3, (2 / 3) .* V) .* (λ / τ) .+ vcm
    mb, xb, vb = background_cluster(n_background, rng)
    return vcat(fill(m, 3), mb),
    hcat(xs, xb),
    hcat(vs, vb),
    (; m, τ, λ, e_in, ratio, inclination, span, t_units)
end

"""
    kozai_input(dir, mass, pos, vel, sys; path, interval = 2.0^-8) -> String

Input of the triple of [`kozai_triple`](@ref): the inner binary is a KS pair from the start; `path
= :ks` leaves the triple to the KS and Hermite schemes (`KZ(15) = 0`: the inner pair is a perturbed
KS pair, the outer body a single star), `path = :merger` enables the few-body regularisation
(`KZ(15) = 2`, `KZ(30) = 2`, `KZ(18) = 1`), under which a stable hierarchy becomes a merger.
"""
function kozai_input(
    dir,
    mass,
    pos,
    vel,
    sys;
    path::Symbol,
    interval = 2.0^-8,
    rmin = 12 * PLUMMER_RH / (length(mass) - 3),
    dtmin = 0.016 * sqrt(ETA_BASELINE.etai / 0.02) * sqrt(rmin^3 * length(mass)),
)
    kz = if path === :ks
        Dict(15 => 0)
    elseif path === :merger
        Dict(15 => 2, 30 => 2, 18 => 1)
    else
        throw(ArgumentError("kozai_input: path must be :ks or :merger, got $(repr(path))"))
    end
    return planted_input(
        dir,
        mass,
        pos,
        vel;
        span = sys.span,
        interval,
        nbin0 = 1,
        n_background = length(mass) - 3,
        kz,
        rmin,
        dtmin,
    )
end

"""
    pythagorean_triple(rng; μ = 1e-3, span = 0.5, t_units = 80.0, n_background = 256,
                       pair_first = false) -> (mass, pos, vel, sys)

Burrau's problem (Szebehely & Peters 1967, DOI 10.1086/110355): masses 3, 4, 5 at rest at (1, 3),
(−2, −1), (1, −1) with G = 1, scaled by the mass unit μ and the time unit τ = span/t_units (length
unit λ = (τ² μ)^{1/3}), moving at a small common velocity, beside a background Plummer model 100
length units away. The bodies are named in the order of masses 3, 4, 5, or 4, 5, 3 with
`pair_first`, when the pair of masses 4 and 5 is to be regularised from the start.
`sys.masses` lists the masses by name.
"""
function pythagorean_triple(
    rng::AbstractRNG;
    μ = 1.0e-3,
    span = 0.5,
    t_units = 80.0,
    n_background = 256,
    pair_first::Bool = false,
)
    τ = span / t_units
    λ = cbrt(τ^2 * μ)
    masses = pair_first ? [4.0, 5.0, 3.0] : [3.0, 4.0, 5.0]
    place = Dict(3.0 => [1.0, 3.0, 0.0], 4.0 => [-2.0, -1.0, 0.0], 5.0 => [1.0, -1.0, 0.0])
    xs = hcat((place[m] for m in masses)...) .* λ
    vcm = [2.0e-3, 1.0e-3, -1.0e-3]
    mb, xb, vb = background_cluster(n_background, rng)
    return vcat(μ .* masses, mb),
    hcat(xs, xb),
    hcat(repeat(vcm, 1, 3), vb),
    (; μ, τ, λ, span, t_units, masses, pair_first)
end

"""
    pythagorean_input(dir, mass, pos, vel, sys; kz26 = 0) -> String

Input of the triple of [`pythagorean_triple`](@ref) with the chain enabled (`KZ(15) = 2`,
`KZ(30) = 2`) and `KZ(26) = kz26`, energy check every 2⁻⁶ and snapshot every 2⁻⁹ N-body times;
with `sys.pair_first` the pair of masses 4 and 5 is a KS pair from the start (`NBIN0 = 1`).
"""
function pythagorean_input(dir, mass, pos, vel, sys; kz26 = 0)
    kz = Dict(15 => 2, 30 => 2, 26 => kz26)
    return planted_input(
        dir,
        mass,
        pos,
        vel;
        span = sys.span,
        interval = 2.0^-6,
        deltat = 2.0^-9,
        nbin0 = sys.pair_first ? 1 : 0,
        n_background = length(mass) - 3,
        kz,
    )
end

"""
    triple_energy(s, sys) -> Float64

Energy of the three bodies of [`pythagorean_triple`](@ref) in a snapshot, in its dimensionless
units (−769/60 ≈ −12.817 at the start).
"""
function triple_energy(s, sys)
    m = sys.masses
    x = zeros(3, 3)
    v = zeros(3, 3)
    for k in 1:3
        i = findfirst(==(Int32(k)), s.name)
        i === nothing && throw(ArgumentError("triple_energy: body $k not in the snapshot"))
        x[:, k] = Float64.(s.pos[:, i]) ./ sys.λ
        v[:, k] = Float64.(s.vel[:, i]) ./ (sys.λ / sys.τ)
    end
    vc = v * m ./ sum(m)
    T = 0.5 * sum(m[k] * sum(abs2, v[:, k] .- vc) for k in 1:3)
    W = -sum(m[i] * m[j] / norm(x[:, i] .- x[:, j]) for i in 1:3 for j in (i + 1):3)
    return T + W
end

"""
    triple_cm_speed(s, sys) -> Float64

Speed of the centre of mass of the three bodies of [`pythagorean_triple`](@ref) in a snapshot
(N-body units); 2.45×10⁻³ at the start.
"""
function triple_cm_speed(s, sys)
    p = zeros(3)
    for (k, m) in enumerate(sys.masses)
        i = findfirst(==(Int32(k)), s.name)
        i === nothing && throw(ArgumentError("triple_cm_speed: body $k not in the snapshot"))
        p .+= m .* Float64.(s.vel[:, i])
    end
    return norm(p) / sum(sys.masses)
end

"""
    triple_outcome(s, sys) -> NamedTuple

The state of the Pythagorean triple in a snapshot, in its dimensionless units: the pair of
bodies with the most negative two-body energy (`pair`, by mass), its semi-major axis `a` and
eccentricity `e`, the distance `r3` of the third body from the pair's centre of mass and that
body's energy per unit reduced mass in the pair's field (`E3`; positive when it escapes).
"""
function triple_outcome(s, sys)
    m = sys.masses
    x = zeros(3, 3)
    v = zeros(3, 3)
    for k in 1:3
        i = findfirst(==(Int32(k)), s.name)
        i === nothing && throw(ArgumentError("triple_outcome: body $k not in the snapshot"))
        x[:, k] = Float64.(s.pos[:, i]) ./ sys.λ
        v[:, k] = Float64.(s.vel[:, i]) ./ (sys.λ / sys.τ)
    end
    best = (Inf, 0, 0)
    for (i, j) in ((1, 2), (1, 3), (2, 3))
        el = orbital_elements(m[i], m[j], x[:, j] - x[:, i], v[:, j] - v[:, i])
        el.h < best[1] && (best = (el.h, i, j))
    end
    _, i, j = best
    k = 6 - i - j
    el = orbital_elements(m[i], m[j], x[:, j] - x[:, i], v[:, j] - v[:, i])
    mb = m[i] + m[j]
    xc = (m[i] .* x[:, i] .+ m[j] .* x[:, j]) ./ mb
    vc = (m[i] .* v[:, i] .+ m[j] .* v[:, j]) ./ mb
    R = x[:, k] .- xc
    E3 = 0.5 * sum(abs2, v[:, k] .- vc) - (mb + m[k]) / norm(R)
    return (; pair = (m[i], m[j]), a = el.a, e = el.e, r3 = norm(R), E3, escaper = m[k])
end

end # module
