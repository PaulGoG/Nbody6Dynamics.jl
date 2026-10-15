# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Energy-jump localisation, mass budget and the Jacobi integral of a
# linearised tidal field. The numerical references are computed here, not
# taken from the implementation.

"""Snapshot of stars `mass` at `pos`, `vel` (3 × n), time `t_nb`, `tscale` Myr per NB time."""
function _eb_snapshot(mass, pos, vel, t_nb; tscale = 1.0, zmbar = 1.0)
    params = zeros(Float32, 20)
    params[1] = Float32(t_nb)
    params[3] = 1.0f0
    params[4] = Float32(zmbar)
    params[11] = Float32(tscale)
    params[12] = 1.0f0
    n = length(mass)
    hdr = SnapshotHeader(Int32(n), Int32(0), Int32(1), Int32(20), params)
    return Snapshot(
        hdr,
        Int32.(1:n),
        Float32.(mass),
        Float32.(pos),
        Float32.(vel),
        Float32[],
        Float32[],
    )
end

@testset "Jacobi integral in a linearised tidal field" begin
    # Two bodies in the rotating frame of a point-mass galaxy with Ω = 1
    # (TIDAL = 3Ω², 0, −Ω², 2Ω), integrated with the engine's force law of
    # xtrnlf.F: F = (T₁x + T₄ẏ, −T₄ẋ, T₃z) plus the mutual attraction. The
    # Coriolis term does no work and the tidal potential is static, so
    # T + W − ½Σm(T₁x² + T₃z²) is conserved while T + W is not.
    T₁, T₃, T₄ = 3.0, -1.0, 2.0
    tidal = (T₁, 0.0, T₃, T₄)
    m = [0.6, 0.4]
    function accel!(a, x, v)
        d = x[:, 1] .- x[:, 2]
        r3 = (d[1]^2 + d[2]^2 + d[3]^2)^1.5
        for i in 1:2
            s = i == 1 ? -1.0 : 1.0
            mj = m[3 - i]
            a[1, i] = s * mj * d[1] / r3 + T₁ * x[1, i] + T₄ * v[2, i]
            a[2, i] = s * mj * d[2] / r3 - T₄ * v[1, i]
            a[3, i] = s * mj * d[3] / r3 + T₃ * x[3, i]
        end
        return a
    end
    jacobi(x, v) = begin
        T = 0.5 * sum(m[i] * (v[1, i]^2 + v[2, i]^2 + v[3, i]^2) for i in 1:2)
        r = sqrt(sum((x[k, 1] - x[k, 2])^2 for k in 1:3))
        W = -m[1] * m[2] / r
        Wt = -0.5 * sum(m[i] * (T₁ * x[1, i]^2 + T₃ * x[3, i]^2) for i in 1:2)
        (T + W, T + W + Wt)
    end
    # A bound pair (a ≈ 0.32) inside the tidal radius (M/T₁)^(1/3) ≈ 0.69,
    # its centre of mass at rest at the origin, the relative velocity with
    # components off the orbital plane so every term of the force acts.
    d0 = [0.3, 0.1, 0.16]
    w0 = [0.3, 1.5, 0.4]
    x0 = hcat(0.4 .* d0, -0.6 .* d0)
    v0 = hcat(0.4 .* w0, -0.6 .* w0)
    x = copy(x0)
    v = copy(v0)
    a = similar(x)
    h = 1.0e-4
    nsteps = 20_000  # t = 2, about two orbits of the pair
    kx = ntuple(_ -> similar(x), 4)
    kv = ntuple(_ -> similar(v), 4)
    e0, j0 = jacobi(x, v)
    jmax = 0.0
    emax = 0.0
    for _ in 1:nsteps
        kx[1] .= v
        kv[1] .= accel!(a, x, v)
        kx[2] .= v .+ 0.5h .* kv[1]
        kv[2] .= accel!(a, x .+ 0.5h .* kx[1], v .+ 0.5h .* kv[1])
        kx[3] .= v .+ 0.5h .* kv[2]
        kv[3] .= accel!(a, x .+ 0.5h .* kx[2], v .+ 0.5h .* kv[2])
        kx[4] .= v .+ h .* kv[3]
        kv[4] .= accel!(a, x .+ h .* kx[3], v .+ h .* kv[3])
        x .+= h / 6 .* (kx[1] .+ 2 .* kx[2] .+ 2 .* kx[3] .+ kx[4])
        v .+= h / 6 .* (kv[1] .+ 2 .* kv[2] .+ 2 .* kv[3] .+ kv[4])
        e, j = jacobi(x, v)
        jmax = max(jmax, abs(j - j0))
        emax = max(emax, abs(e - e0))
    end
    @test jmax < 1.0e-7 * abs(j0)       # the Jacobi integral is conserved by the force law
    @test emax > 1.0e-2 * abs(j0)       # the plain energy is not: the tide does work

    # The same two states through snapshot_integrals: e_jacobi of the first
    # and last snapshot agree to the single precision of the snapshot data.
    s0 = _eb_snapshot(m, x0, v0, 0.0)
    s1 = _eb_snapshot(m, x, v, nsteps * h)
    i0 = snapshot_integrals(s0; tidal, min_members = 2)
    i1 = snapshot_integrals(s1; tidal, min_members = 2)
    wt0 = -0.5 * sum(m[i] * (T₁ * x0[1, i]^2 + T₃ * x0[3, i]^2) for i in 1:2)
    @test i0.tidal_potential ≈ wt0 rtol = 1.0e-6
    ej(s) = s.kinetic + s.potential + s.pair_internal + s.tidal_potential
    @test ej(i0) ≈ j0 rtol = 1.0e-5
    @test ej(i1) ≈ ej(i0) rtol = 1.0e-5
    @test abs((i1.kinetic + i1.potential) - (i0.kinetic + i0.potential)) > 1.0e-2 * abs(j0)
    # Without a field the tidal term vanishes and e_jacobi equals e_total.
    @test snapshot_integrals(s0; min_members = 2).tidal_potential == 0.0
    path = joinpath(mktempdir(), "snapshot_integrals.csv")
    write_snapshot_integrals(path, [i0, i1])
    header = split(readline(path), ",")
    @test header[findfirst(==("e_pairs"), header) + 1] == "w_tidal"
    @test header[findfirst(==("e_total"), header) + 1] == "e_jacobi"
    row = split(readlines(path)[2], ",")
    @test parse(Float64, row[findfirst(==("e_jacobi"), header)]) ≈ ej(i0) rtol = 1.0e-9
    @test parse(Float64, row[findfirst(==("w_tidal"), header)]) ≈ i0.tidal_potential rtol = 1.0e-9

    @test tidal_mode.(0:5) == [:none, :linearised, :linearised, :orbit, :external, :orbit]
    @test_throws ArgumentError tidal_mode(6)
    @test tidal_mode(DiagnosticsData(AdjustRecord[], Dict{String,Float64}())) == :none
    kz = zeros(Int, 50)
    kz[14] = 2
    @test tidal_mode(DiagnosticsData(AdjustRecord[], Dict{String,Float64}(), kz, tidal)) ==
          :linearised
    @test engine_energy_valid(:none) && engine_energy_valid(:orbit)
    @test !engine_energy_valid(:linearised) && !engine_energy_valid(:external)
    @test_throws ArgumentError engine_energy_valid(:tide)
end

@testset "Energy-jump localisation" begin
    # Six adjustments; the DE series has one outlier (interval 5) and the
    # correction channels move in intervals 3 (EBIN) and 5 (EMDOT).
    t = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    de = [0.0, 1.0e-6, -2.0e-6, 1.5e-6, 2.0e-3, -1.0e-6]
    ebin = [0.0, 0.0, -0.01, -0.01, -0.01, -0.01]
    emdot = [0.0, 0.0, 0.0, 0.0, 5.0e-4, 5.0e-4]
    rec(k) = AdjustRecord(
        t[k],
        t[k],
        0.5,
        de[k],
        -0.25,
        100,
        0,
        1.0,
        de[k] * 0.25,
        0.0,
        0.125,
        0.375,
        0.0,
        ebin[k],
        0.0,
        0.0,
        0.0,
        emdot[k],
        0.0,
    )
    diag = DiagnosticsData([rec(k) for k in 1:6], Dict{String,Float64}())
    ev(kind, interval; kw...) = EngineEvent(; time_nb = NaN, kind, interval, line = interval, kw...)
    events = [
        ev(:ks_rectification, 3),
        ev(:ks_rectification, 5),
        ev(:kick_full, 5; name1 = 7, velocity_kms = 120.0),
        ev(:remnant_formation, 5; name1 = 7, kstar1 = 5, kstar2 = 14),
        ev(:escape, 7; count = 2, n_total = 98),
    ]
    a = energy_intervals(diag, events)
    @test length(a.intervals) == 5
    @test [i.index for i in a.intervals] == 2:6
    @test [i.t_start for i in a.intervals] == t[1:5]
    @test [i.t_end for i in a.intervals] == t[2:6]
    @test [i.de_rel for i in a.intervals] == de[2:6]
    @test a.intervals[2].d_ebin == -0.01 && a.intervals[3].d_ebin == 0.0
    @test a.intervals[4].d_emdot == 5.0e-4 && a.intervals[1].d_emdot == 0.0
    kinds = collect(ENGINE_EVENT_KINDS)
    @test a.intervals[2].counts[findfirst(==(:ks_rectification), kinds)] == 1
    @test a.intervals[4].counts[findfirst(==(:ks_rectification), kinds)] == 1
    @test a.intervals[4].counts[findfirst(==(:kick_full), kinds)] == 1
    @test a.intervals[4].counts[findfirst(==(:remnant_formation), kinds)] == 1
    @test a.intervals[4].n_events == 3
    @test a.intervals[5].n_events == 0            # the escape lies after the last ADJUST line
    # Robust scale: x = de[2:6]; median 1e-6; |x − med| = 0, 3e-6, 5e-7, 1.999e-3, 2e-6;
    # median of those 2e-6.
    @test a.scale ≈ 1.4826 * 2.0e-6 rtol = 1.0e-12
    level = max(1.0e-4, 5 * a.scale)
    @test [i.flagged for i in a.intervals] == [abs(x) ≥ level for x in de[2:6]]
    @test count(i -> i.flagged, a.intervals) == 1 && a.intervals[4].flagged
    @test a.intervals[4].score ≈ 2.0e-3 / a.scale rtol = 1.0e-12
    # A high floor unflags the outlier; a zero floor with sigma = 0 flags all finite DE.
    @test !any(i -> i.flagged, energy_intervals(diag, events; threshold = 5.0e-3).intervals)
    @test all(
        i -> i.flagged,
        energy_intervals(diag, events; threshold = 0.0, sigma = 0.0).intervals,
    )
    # A linearised tide makes the engine's DE a drift: nothing is flagged.
    b = energy_intervals(diag, events; engine_energy_valid = false)
    @test !any(i -> i.flagged, b.intervals) && !b.engine_energy_valid
    @test occursin("not an error measure", energy_events_summary(b))
    s = energy_events_summary(a)
    @test occursin("5 adjustment intervals", s) && occursin("1 flagged", s)
    @test occursin("kick_full=1", s) && occursin("remnant_formation=1", s)
    @test_throws ArgumentError energy_intervals(diag, events; threshold = -1.0)
    @test_throws ArgumentError energy_intervals(diag, events; sigma = -1.0)
    # NaN channels (positional ADJUST lines) propagate, finite DE still flags.
    bare = [AdjustRecord(t[k], t[k], 0.5, de[k], -0.25, 100, 0, 1.0, NaN, NaN) for k in 1:6]
    c = energy_intervals(DiagnosticsData(bare, Dict{String,Float64}()), EngineEvent[])
    @test all(i -> isnan(i.d_ebin), c.intervals) && c.intervals[4].flagged
    # Fewer than two records: no interval, zero scale.
    e = energy_intervals(DiagnosticsData([rec(1)], Dict{String,Float64}()), EngineEvent[])
    @test isempty(e.intervals) && e.scale == 0.0

    dir = mktempdir()
    p1 = write_energy_intervals(joinpath(dir, "energy_intervals.csv"), a)
    lines = readlines(p1)
    header = split(lines[1], ",")
    @test length(header) == 14 + length(ENGINE_EVENT_KINDS)
    @test header[end] == "n_chain_end" && header[15] == "n_ks_rectification"
    @test length(lines) == 6
    r4 = split(lines[5], ",")
    @test r4[findfirst(==("flagged"), header)] == "true"
    @test r4[findfirst(==("n_events"), header)] == "3"
    @test parse(Float64, r4[findfirst(==("d_emdot"), header)]) == 5.0e-4
    p2 = write_engine_events(joinpath(dir, "engine_events.csv"), events)
    elines = readlines(p2)
    @test length(elines) == 6
    eh = split(elines[1], ",")
    e3 = split(elines[4], ",")
    @test e3[findfirst(==("kind"), eh)] == "kick_full"
    @test e3[findfirst(==("kstar1"), eh)] == ""          # KSTAR_ABSENT → empty cell
    @test e3[findfirst(==("time_nb"), eh)] == ""         # NaN → empty cell
    @test parse(Float64, e3[findfirst(==("velocity_kms"), eh)]) == 120.0
    e4 = split(elines[5], ",")
    @test e4[findfirst(==("kstar2"), eh)] == "14"
    # Never overwrite: a second write backs the first up.
    write_energy_intervals(joinpath(dir, "energy_intervals.csv"), a)
    @test count(f -> startswith(f, "energy_intervals"), readdir(dir)) == 2
end

@testset "Mass budget" begin
    # Five stars of 0.2 at t = 0. By t = 1 (2 Myr at T* = 2): star 5 escaped
    # (0.2), star 1 lost 0.05 to stellar evolution (ZMDOT), stars 2 + 3 merged
    # with 0.01 lost in the collision. M* = 100 M☉ per N-body mass.
    zmbar, tscale = 100.0, 2.0
    u = UnitScaling(1.0, zmbar, tscale, 1.0)
    pos0 = [1.0 2.0 3.0 4.0 5.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]
    vel0 = zeros(3, 5)
    s0 = _eb_snapshot(fill(0.2, 5), pos0, vel0, 0.0; tscale, zmbar)
    s1 = _eb_snapshot([0.15, 0.39, 0.2], pos0[:, 1:3], vel0[:, 1:3], 1.0; tscale, zmbar)
    ints = snapshot_integrals([s0, s1]; min_members = 2)
    escapers = [EscaperRecord(0.5 * tscale, 0.2 * zmbar, 1.0, 10.0, 0, 5, 0.0, 0.0)]
    counters = [
        EventCounters(
            0.0,
            zeros(Int, 13)...,
            zeros(12)...,
            zeros(Int, 9)...,
            zeros(7)...,
            ntuple(_ -> 0, 16),
        ),
        EventCounters(
            1.0 * tscale,
            zeros(Int, 13)...,
            zeros(12)...,
            zeros(Int, 9)...,
            zeros(6)...,
            0.05 * zmbar,
            ntuple(_ -> 0, 16),
        ),
    ]
    events = [
        EngineEvent(;
            time_nb = 0.25,
            kind = :escape,
            interval = 2,
            line = 10,
            count = 1,
            n_total = 4,
            mass_solar = 0.8 * zmbar,
        ),
        EngineEvent(;
            time_nb = 0.75,
            kind = :merged_star,
            interval = 2,
            line = 20,
            name1 = 2,
            name2 = 3,
            mass_solar = 0.39 * zmbar,
            dmass_solar = 0.01 * zmbar,
        ),
    ]
    rows = mass_budget(ints, escapers, counters, events, u)
    @test length(rows) == 2
    r0, r1 = rows
    @test r0.m_snapshot ≈ 1.0 rtol = 1.0e-6
    @test r0.m_escaped == 0.0 && r0.m_lost_stellar == 0.0 && r0.m_lost_collision == 0.0
    @test isnan(r0.m_engine)
    @test abs(r0.residual) < 1.0e-6 && abs(r0.residual_rel) < 1.0e-6
    @test r1.n_stars == 3
    @test r1.m_snapshot ≈ 0.74 rtol = 1.0e-6
    @test r1.m_escaped ≈ 0.2 rtol = 1.0e-12
    @test r1.m_lost_stellar ≈ 0.05 rtol = 1.0e-12
    @test r1.m_lost_collision ≈ 0.01 rtol = 1.0e-12
    @test r1.m_engine ≈ 0.8 rtol = 1.0e-12
    @test abs(r1.residual) < 2.0e-7                     # closes at single precision
    @test r1.time_myr ≈ 2.0
    # Inclusive bounds: an escape at exactly the snapshot time counts.
    late = [EscaperRecord(1.0 * tscale, 0.2 * zmbar, 1.0, 10.0, 0, 5, 0.0, 0.0)]
    @test mass_budget(ints, late, counters, events, u)[2].m_escaped ≈ 0.2
    after = [EscaperRecord(1.0 * tscale + 0.5, 0.2 * zmbar, 1.0, 10.0, 0, 5, 0.0, 0.0)]
    @test mass_budget(ints, after, counters, events, u)[2].m_escaped == 0.0
    # Without counters and events the stellar and collision terms are zero.
    r = mass_budget(ints, escapers, EventCounters[], EngineEvent[], u)[2]
    @test r.m_lost_stellar == 0.0 && r.m_lost_collision == 0.0 && isnan(r.m_engine)
    @test r.residual ≈ 0.06 atol = 1.0e-6
    @test_throws ArgumentError mass_budget(SnapshotIntegrals[], escapers, counters, events, u)
    @test_throws ArgumentError mass_budget(
        ints,
        escapers,
        counters,
        events,
        UnitScaling(1.0, 0.0, 1.0, 1.0),
    )
    s = mass_budget_summary(rows)
    @test occursin("2 snapshots", s) && occursin("escaped 20 %", s)
    @test mass_budget_summary(MassBudgetRow[]) == "mass budget: no snapshots"
    path = write_mass_budget(joinpath(mktempdir(), "mass_budget.csv"), rows)
    lines = readlines(path)
    @test split(lines[1], ",") == [
        "time_nb",
        "time_myr",
        "n_stars",
        "m_snapshot",
        "m_escaped",
        "m_lost_stellar",
        "m_lost_collision",
        "m_engine",
        "residual",
        "residual_rel",
    ]
    @test length(lines) == 3
    @test split(lines[2], ",")[8] == ""                 # NaN m_engine → empty cell
    @test parse(Float64, split(lines[3], ",")[5]) == 0.2
end

@testset "Energy budget" begin
    # Two snapshots of the same two bodies; between them an escaper of energy
    # 0.3 was removed and the engine booked a kick (ECDOT) and mass loss (EMDOT).
    zmbar, tscale = 100.0, 2.0
    u = UnitScaling(1.0, zmbar, tscale, 1.0)
    pos = [0.0 1.0; 0.0 0.0; 0.0 0.0]
    vel = [0.0 0.0; 0.5 -0.5; 0.0 0.0]
    s0 = _eb_snapshot([0.5, 0.5], pos, vel, 0.0; tscale, zmbar)
    s1 = _eb_snapshot([0.5, 0.5], pos, 0.8 .* vel, 2.0; tscale, zmbar)
    ints = snapshot_integrals([s0, s1]; min_members = 2)
    e0 = ints[1].kinetic + ints[1].potential
    e1 = ints[2].kinetic + ints[2].potential
    adj(t, ecoll, emdot, ecdot) = AdjustRecord(
        t,
        t * tscale,
        0.5,
        0.0,
        e0,
        2,
        0,
        1.0,
        0.0,
        1.5e-5 * t,
        0.1,
        0.35,
        0.0,
        0.0,
        0.0,
        0.0,
        ecoll,
        emdot,
        ecdot,
    )
    diag = DiagnosticsData(
        [adj(0.0, 0.0, 0.0, 0.0), adj(1.0, 0.0, -0.02, -0.1), adj(2.0, 0.01, -0.05, -0.1)],
        Dict{String,Float64}(),
    )
    zc = (zeros(Int, 13)...,)
    counters = [
        EventCounters(0.0, zc..., zeros(12)..., zeros(Int, 9)..., zeros(7)..., ntuple(_ -> 0, 16)),
        EventCounters(
            1.0 * tscale,
            zc...,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.2,
            0.1,
            0.0,
            0.0,
            0.0,
            0.0,
            zeros(Int, 9)...,
            zeros(7)...,
            ntuple(_ -> 0, 16),
        ),
        EventCounters(
            2.0 * tscale,
            zc...,
            0.0,
            0.0,
            0.2,
            0.3,
            0.0,
            0.0,
            0.25,
            0.1,
            0.05,
            0.0,
            0.0,
            0.0,
            zeros(Int, 9)...,
            zeros(7)...,
            ntuple(_ -> 0, 16),
        ),
    ]
    escapers = [
        EscaperRecord(0.5 * tscale, 0.1 * zmbar, 2.0, 10.0, 0, 7, 0.0, 0.0),   # m·EESC = 0.2
        EscaperRecord(1.5 * tscale, 0.1 * zmbar, 1.5, 10.0, 0, 8, 0.0, 0.0),   # m·EESC = 0.15
        EscaperRecord(2.5 * tscale, 0.1 * zmbar, 1.0, 10.0, 0, 9, 0.0, 0.0),   # after the last snapshot
    ]
    rows = energy_budget(ints, escapers, counters, diag, u)
    @test length(rows) == 2 && all(r -> r.source == :counters && r.pair_source == :snapshot, rows)
    r0, r1 = rows
    @test r0.e_snapshot ≈ e0 && r0.e_escaped == 0.0 && r0.e_channels == 0.0
    @test r0.e_total ≈ e0 && r0.residual == 0.0 && r0.residual_rel == 0.0 && r0.detot == 0.0
    @test r1.e_snapshot ≈ e1 && r1.e_pairs == 0.0
    @test r1.e_escaped == 0.25 + 0.1 + 0.05                 # ESESC + EBESC + EMESC at t = 2
    @test r1.e_channels == 0.5                              # ECOLL + EMDOT + ECDOT of the event.35 row at t = 2
    @test r1.detot == 3.0e-5                                # DETOT of the ADJUST at t = 2
    @test r1.e_total ≈ e1 + 0.4 + 0.5
    @test r1.residual ≈ r1.e_total - e0 && r1.residual_rel ≈ r1.residual / abs(e0)
    # Without counters the escaper energies come from esc.11, cumulative to the
    # snapshot time, and the corrections from the ADJUST line.
    re = energy_budget(ints, escapers, EventCounters[], diag, u)
    @test all(r -> r.source == :escapers, re)
    @test re[1].e_escaped == 0.0 && re[2].e_escaped ≈ 0.2 + 0.15
    @test re[2].e_channels == 0.01 - 0.05 - 0.1
    rn = energy_budget(ints, EscaperRecord[], EventCounters[], diag, u)
    @test all(r -> r.source == :none && r.e_escaped == 0.0, rn)
    # Pair energies from bev.82 where a snapshot time is covered: one pair of
    # 10 and 20 M☉ at a = 10³ R☉ has E_b = G m₁ m₂ / 2a; in N-body units over
    # M* V*².
    rec = BinaryRecord(
        2.0,
        Int32(1),
        Int32(2),
        Int32(1),
        Int32(2),
        Int32(0),
        Int32(0),
        Int32(0),
        0.0,
        0.1,
        0.0,
        3.0,
        10.0,
        20.0,
        zeros(18)...,
    )
    bev = BinaryEvolutionSnapshot(round(2.0 * tscale; digits = 1), 1, [rec])   # header TPHYS rounded
    rb = energy_budget(ints, escapers, counters, diag, u; bevs = [bev])
    @test rb[1].pair_source == :snapshot && rb[2].pair_source == :bev
    eb = -binding_energy(rec) / (zmbar * u.vstar^2)
    @test rb[2].e_pairs ≈ eb && rb[2].e_snapshot ≈ e1 + eb
    @test occursin("bev.82 at 1 snapshots", energy_budget_summary(rb))
    # No ADJUST record at or before a snapshot: the DETOT is NaN; without
    # counters the channels are NaN too.
    rl = energy_budget(
        ints,
        escapers,
        EventCounters[],
        DiagnosticsData([adj(3.0, 0.0, 0.0, 0.0)], Dict{String,Float64}()),
        u,
    )
    @test all(r -> isnan(r.e_channels) && isnan(r.e_total) && isnan(r.detot), rl)
    @test_throws ArgumentError energy_budget(SnapshotIntegrals[], escapers, counters, diag, u)
    @test_throws ArgumentError energy_budget(
        ints,
        escapers,
        counters,
        diag,
        u;
        time_tolerance = -1.0,
    )
    s = energy_budget_summary(rows)
    @test occursin("2 snapshots", s) &&
          occursin("from counters", s) &&
          occursin("bev.82 at 0 snapshots", s)
    @test energy_budget_summary(EnergyBudgetRow[]) == "energy budget: no snapshots"
    path = write_energy_budget(joinpath(mktempdir(), "energy_budget.csv"), rows)
    lines = readlines(path)
    @test split(lines[1], ",")[(end - 1):end] == ["source", "pair_source"] && length(lines) == 3
    @test split(lines[3], ",")[(end - 1):end] == ["counters", "snapshot"]
    @test parse(Float64, split(lines[3], ",")[5]) == 0.4
    @test parse(Float64, split(lines[3], ",")[10]) == 3.0e-5

    # Exact stamps of the event lines against the rounded Myr stamps: a
    # collision just after a snapshot is not counted into it.
    pos5 = [1.0 2.0 3.0 4.0 5.0; 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0]
    m5 = snapshot_integrals(
        [_eb_snapshot(fill(0.2, 5), pos5, zeros(3, 5), 0.0; tscale, zmbar)];
        min_members = 2,
    )
    late = [
        EngineEvent(;
            time_nb = 3.8e-6,
            kind = :merged_star,
            interval = 2,
            line = 5,
            dmass_solar = 0.01 * zmbar,
        ),
    ]
    @test mass_budget(m5, EscaperRecord[], EventCounters[], late, u)[1].m_lost_collision == 0.0
    # A Myr-stamped escaper at the snapshot time within the rounding is counted.
    esc_r = [EscaperRecord(100.0 * tscale * (1 + 4.0e-6), 0.2 * zmbar, 1.0, 10.0, 0, 5, 0.0, 0.0)]
    m100 = snapshot_integrals(
        [_eb_snapshot(fill(0.2, 5), pos5, zeros(3, 5), 100.0; tscale, zmbar)];
        min_members = 2,
    )
    @test mass_budget(m100, esc_r, EventCounters[], EngineEvent[], u)[1].m_escaped ≈ 0.2
end
