# Merger epochs: the detectors, the regimes of snapshots, progenitors and pairs, the masks of
# the remnant and progenitor-budget quantities, the CSV layouts and the pipeline products, on
# synthetic sequences of static clusters placed by hand.

using StableRNGs: StableRNG

"""Snapshot of N-body arrays at time `t` (rbar 2 pc, tscale 3 Myr), names 1:n."""
function _ep_snap(pos, vel, mass, t)
    params = zeros(Float32, 20)
    params[1] = t
    params[3] = 2.0f0
    params[4] = 1.0f0
    params[11] = 3.0f0
    params[12] = 1.0f0
    n = length(mass)
    return Snapshot(
        SnapshotHeader(Int32(n), Int32(1), Int32(1), Int32(20), params),
        Int32.(1:n),
        Float32.(mass),
        Float32.(pos),
        Float32.(vel),
        Float32[],
        Float32[],
        0,
        Int32[],
    )
end

"""Plummer sphere of `n` stars and scale `a` with the virial velocities of a cluster of mass `M`."""
function _ep_clump(rng, n, a, M)
    p, v = Nbody6Dynamics.sample_plummer(n, a; rng)
    v .*= sqrt(M)
    return p, v
end

"""Two equal clusters of 400 stars (mass ½ each) at the separations `seps` along x, at rest."""
function _ep_pair_series(seps)
    rng = StableRNG(11)
    n = 400
    p1, v1 = _ep_clump(rng, n, 0.1, 0.5)
    p2, v2 = _ep_clump(rng, n, 0.1, 0.5)
    mass = fill(1 / (2n), 2n)
    snaps = [
        _ep_snap(
            hcat(p1 .- [d / 2, 0, 0], p2 .+ [d / 2, 0, 0]),
            hcat(v1, v2),
            mass,
            Float64(k - 1),
        ) for (k, d) in enumerate(seps)
    ]
    return snaps, [1:n, (n + 1):(2n)]
end

"""
A primary of 600 stars (mass ¾) and a secondary of 200 (mass ¼): the secondary approaches
intact (snapshots 1–3), is replaced by unbound debris off the primary's centre (4–5) and
by the same debris centred on the primary (6–8).
"""
function _ep_shredded_series()
    rng = StableRNG(13)
    n_p, n_s = 600, 200
    mass = fill(1 / (n_p + n_s), n_p + n_s)
    p_p, v_p = _ep_clump(rng, n_p, 0.1, 0.75)
    p_s, v_s = _ep_clump(rng, n_s, 0.1, 0.25)
    p_d, v_d = _ep_clump(rng, n_s, 0.3, 1.0)
    v_d .*= 1.5
    states = [
        (p_s .+ [4.0, 0, 0], v_s),
        (p_s .+ [2.0, 0, 0], v_s),
        (p_s .+ [0.8, 0, 0], v_s),
        (p_d .+ [0.5, 0, 0], v_d),
        (p_d .+ [0.4, 0, 0], v_d),
        (p_d, v_d),
        (p_d, v_d),
        (p_d, v_d),
    ]
    snaps = [
        _ep_snap(hcat(p_p, p), hcat(v_p, v), mass, Float64(k - 1)) for
        (k, (p, v)) in enumerate(states)
    ]
    return snaps, [1:n_p, (n_p + 1):(n_p + n_s)]
end

"""
Clusters 1 and 2 (450 stars each) merge at snapshot 3; cluster 3, a compact satellite of
100 stars, falls in and orbits at 0.4 from the remnant centre from snapshot 3 on.
"""
function _ep_satellite_series()
    rng = StableRNG(17)
    n1, n2, n3 = 450, 450, 100
    mass = fill(1 / 1000, 1000)
    p1, v1 = _ep_clump(rng, n1, 0.1, 0.45)
    p2, v2 = _ep_clump(rng, n2, 0.1, 0.45)
    p3, v3 = _ep_clump(rng, n3, 0.05, 0.1)
    orbit = [
        ([0, 4.0, 0], [0, 0, 0.0]),
        ([0, 3.0, 0], [0, 0, 0.0]),
        ([0, 0.4, 0], [1.0, 0, 0]),
        ([0.4, 0, 0], [0, -1.0, 0]),
        ([0, -0.4, 0], [-1.0, 0, 0]),
        ([-0.4, 0, 0], [0, 1.0, 0]),
    ]
    x12 = [2.0, 0.3, 0.01, 0.01, 0.01, 0.01]
    snaps = [
        _ep_snap(
            hcat(p1 .- [x12[k], 0, 0], p2 .+ [x12[k], 0, 0], p3 .+ orbit[k][1]),
            hcat(v1, v2, v3 .+ orbit[k][2]),
            mass,
            Float64(k - 1),
        ) for k in 1:6
    ]
    return snaps, [1:n1, (n1 + 1):(n1 + n2), (n1 + n2 + 1):1000]
end

"""Rows of a CSV written by the package (comment lines dropped), split at commas."""
_ep_rows(path) = [split(l, ",") for l in readlines(path) if !startswith(l, "#")]

const _EP_SETTINGS = PostprocessConfig(; msr_n_tracer = 20, msr_n_random = 20, rotation_n_min = 100)

@testset "Merger epochs: detectors and regime labels" begin
    R_L = Nbody6Dynamics._roche_lobe_radius
    @test R_L(1.0, 1.0) ≈ 0.49 / (0.6 + log(2))
    @test R_L(2.0, 1.0) ≈ 2 * R_L(1.0, 1.0)
    @test R_L(1.0, 3.0) > R_L(1.0, 1.0) > R_L(1.0, 1 / 3)
    # q = 10⁻⁶: q^{2/3} = 10⁻⁴, q^{1/3} = 10⁻²
    @test R_L(1.0, 1e-6) ≈ 0.49 * 1e-4 / (0.6e-4 + log1p(1e-2)) rtol = 1e-12

    # Contact: d = 1.5 puts R_L = 0.568 below r = 0.9
    @test contact_index([0.5, 0.9], [0.5, 0.9], [1.0, 1.0], [1.0, 1.0], [4.0, 1.5]) == 2
    @test contact_index([0.5, 0.5, 0.5], [0.5, 0.5, 0.5], ones(3), ones(3), [4.0, 1.5, 1.0]) == 3
    @test contact_index(fill(0.1, 3), fill(0.1, 3), ones(3), ones(3), fill(4.0, 3)) == 0
    # Mass ratio 20: R_L = 0.63 for the heavier core, 0.168 for the lighter one
    @test contact_index([0.5], [0.2], [1.0], [0.05], [1.0]) == 1
    @test contact_index([0.5], [0.12], [1.0], [0.05], [1.0]) == 0
    @test contact_index([NaN, 0.9], [0.5, 0.9], [1.0, 1.0], [1.0, 1.0], [0.1, 1.5]) == 2
    @test contact_index([0.5], [0.5], [0.0], [1.0], [0.1]) == 0   # absent core
    @test_throws DimensionMismatch contact_index(ones(2), ones(3), ones(3), ones(3), ones(3))

    # Core coalescence: the trailing run below the threshold, persistence required
    @test core_index([5.0, 3.0, 1.5, 1.2, 0.8]) == 3
    @test core_index([5.0, 1.0, 3.0, 1.0, 1.0]) == 4      # a passage that separates again
    @test core_index([5.0, 1.0, 1.0, NaN, 1.0]) == 2      # an absent core is skipped
    @test core_index([1.0, 1.0, 3.0]) == 0
    @test core_index([3.0, 1.0, NaN]) == 0                # one co-present snapshot, not the last
    @test core_index([3.0, 1.0]; threshold = 0.5) == 0
    @test_throws ArgumentError core_index([1.0]; threshold = 0)

    # Dissolution: the self-bound mass stays below the fraction to the end
    @test dissolution_index([1.0, 0.6, 0.05, 0.02], 1.0) == 3
    @test dissolution_index([1.0, 0.05, 0.6, 0.02], 1.0) == 4
    @test dissolution_index([1.0, NaN, NaN], 1.0) == 2
    @test dissolution_index([1.0, 1.0], 1.0) == 0
    @test dissolution_index([1.0, 0.3], 1.0; fraction = 0.5) == 2
    @test_throws ArgumentError dissolution_index([1.0], 1.0; fraction = 1.0)
    @test_throws ArgumentError dissolution_index([1.0], 0.0)

    # Regimes and their precedence
    @test regime_label(1, 3, 6, 9) === :detached
    @test regime_label(5, 3, 6, 9) === :interacting
    @test regime_label(7, 3, 6, 9) === :coalescing
    @test regime_label(9, 3, 6, 9) === :settled
    @test regime_label(12, 3, 6, 0) === :coalescing      # window open to the end
    @test regime_label(5, 3, 6, 9; dissolved = 4) === :dissolved
    @test regime_label(7, 3, 6, 9; dissolved = 4) === :coalescing
    @test regime_label(5, 3, 0, 0; satellite = true) === :satellite
    @test regime_label(5, 0, 0, 0) === :detached
    @test Nbody6Dynamics._MERGER_REGIMES ==
          (:detached, :interacting, :coalescing, :settled, :dissolved, :satellite)
end

@testset "Merger epochs: approach to core coalescence and a fly-by" begin
    snaps, ranges = _ep_pair_series([4.0, 2.0, 0.6, 0.3, 0.02, 0.02, 0.02])
    bsets = [bound_set(s) for s in snaps]
    co = coalescence_time(snaps, ranges; bound_sets = bsets)
    ep = merger_epochs(snaps, ranges, bsets, co; settings = _EP_SETTINGS)
    @test ep isa MergerEpochs
    @test ep.contact_index == 3 && ep.t_contact == 2.0 && ep.t_contact_myr == 6.0
    @test ep.core_index == 5 && ep.t_core == 4.0 && ep.t_core_myr == 12.0
    @test ep.coalescence_index == 5 && ep.t_coal == co.time_nb
    # 5 crossing times of the merged bound set (≈ 0.1 NB) end before the next snapshot
    @test 0 < ep.t_cr_core < 0.2 && ep.settled_index == 6
    @test ep.regime ==
          [:detached, :detached, :interacting, :interacting, :coalescing, :settled, :settled]
    @test ep.tag_regime[1, :] == ep.tag_regime[2, :] == ep.regime
    @test ep.pair_regime[1, :] == ep.regime
    @test ep.pairs == [(1, 2)] &&
          ep.pair_contact == [3] &&
          ep.pair_core == [5] &&
          ep.pair_join == [5]
    @test ep.group_a == [[1]] && ep.group_b == [[2]] && ep.event_join == [5]
    @test ep.event_contact == [3] && ep.event_core == [5] && ep.event_settled == [6]
    @test ep.main_progenitor == [1] && ep.mass_ratio[1] ≈ 1.0
    @test ep.tag_contact == [3, 3] && ep.tag_core == [5, 5] && ep.tag_dissolved == [0, 0]
    @test all(m -> isapprox(m, 0.5; rtol = 1e-6), ep.self_bound_mass)
    @test all(m -> isapprox(m, 0.5; rtol = 1e-6), ep.bound_mass)
    @test all(ep.roche_fill[1, 1:2] .< 1) && all(ep.roche_fill[1, 3:end] .> 1)
    @test ep.delta_s[1, 4] > 2 && all(ep.delta_s[1, 5:end] .< 2)
    @test length(ep.l_orb) == 3 && all(isfinite, ep.l_orb)   # cores at rest: a residual only
    # The coalescence result is evaluated when not given
    same = merger_epochs(snaps, ranges, bsets; settings = _EP_SETTINGS)
    @test same.regime == ep.regime && same.delta_s == ep.delta_s
    # Thresholds from the configuration
    late = merger_epochs(
        snaps,
        ranges,
        bsets,
        co;
        settings = PostprocessConfig(; merger_core_delta = 0.1),
    )
    @test late.core_index == 0 && all(r -> r in (:detached, :interacting), late.regime)
    half = merger_epochs(
        snaps,
        ranges,
        bsets,
        co;
        settings = PostprocessConfig(; merger_contact_lagrangian = 0.5),
    )
    @test half.contact_index == 4
    quick = merger_epochs(
        snaps,
        ranges,
        bsets,
        co;
        settings = PostprocessConfig(; merger_coalescing_tcr = 0.0),
    )
    @test quick.settled_index == quick.core_index == 5 && quick.regime[5] === :settled
    @test_throws DimensionMismatch merger_epochs(snaps, ranges, bsets[1:3], co)
    @test_throws ArgumentError merger_epochs(snaps, UnitRange{Int}[], bsets)
    @test_throws ArgumentError merger_epochs(Snapshot[], ranges, BoundSet[])
    # One progenitor: one body throughout
    single = merger_epochs(snaps, [1:800], bsets)
    @test all(==(:settled), single.regime) && single.core_index == 0 && isnan(single.t_core)
    @test isempty(single.pairs) && size(single.tag_regime) == (1, 7)

    # Remnant quantities: NaN before the core coalescence, evaluated from it on
    diag = remnant_diagnostics(
        snaps,
        ranges;
        settings = _EP_SETTINGS,
        bound_sets = bsets,
        coalescence = co,
        epochs = ep,
    )
    @test diag.regime == ep.regime
    @test diag.t_contact == 2.0 && diag.t_core == 4.0 && diag.t_core_myr == 12.0
    @test diag.coalescence_time == 4.0
    for f in (
        :r_core,
        :r_core_number,
        :rho_core,
        :r_half,
        :lambda_r,
        :lambda_peebles,
        :lambda_r_floor,
        :spin_alignment,
        :lambda_msr,
        :lambda_msr_err,
        :segregation_ratio,
        :gamma_msr,
    )
        v = getfield(diag, f)
        @test all(isnan, v[1:4])
    end
    @test all(isnan, diag.spin_axis[:, 1:4])
    @test all(isfinite, diag.r_core[5:end]) && all(isfinite, diag.r_half[5:end])
    @test all(isfinite, diag.lambda_r[5:end]) && all(isfinite, diag.lambda_msr[5:end])
    @test all(diag.n_bound .== 800) && all(isfinite, diag.kt)   # bookkeeping on every snapshot
    @test !isempty(diag.profile.radius)
    # The epochs evaluated inside give the same result; mask = false evaluates every snapshot
    inner = remnant_diagnostics(snaps, ranges; settings = _EP_SETTINGS, bound_sets = bsets)
    @test inner.regime == ep.regime && isequal(inner.r_core, diag.r_core)
    raw = remnant_diagnostics(
        snaps,
        ranges;
        settings = _EP_SETTINGS,
        bound_sets = bsets,
        mask = false,
    )
    @test all(isfinite, raw.r_core) && all(isfinite, raw.r_half)
    @test raw.r_core[1] > raw.r_core[end]   # two separated clusters, then one compact body
    @test raw.r_core[5:end] == diag.r_core[5:end]

    # Progenitor budget: the mixing from the core coalescence on, r_h about the remnant centre
    rows = progenitor_budget(snaps, ranges, bsets; settings = _EP_SETTINGS, epochs = ep)
    at(k, tag) = rows[2 * (k - 1) + tag]
    @test [at(k, 1).regime for k in 1:7] == ep.regime
    for k in 1:4, tag in 1:2
        r = at(k, tag)
        @test isfinite(r.r_h) && isfinite(r.centre_offset)
        @test all(isnan, (r.aplus2, r.aplus2_err, r.d_radius, r.d_angular_momentum, r.d_floor))
    end
    for k in 5:7, tag in 1:2
        @test isfinite(at(k, tag).aplus2) && isfinite(at(k, tag).d_radius)
    end
    s = snaps[4]
    b = bsets[4].stars
    mine = [c for c in b if s.name[c] ≤ 400]
    p_i, m_i = Float64.(s.pos[:, mine]), Float64.(s.mass[mine])
    c_i = Nbody6Dynamics.shrinking_sphere_centre(p_i, m_i).centre
    @test at(4, 1).r_h ≈ Nbody6Dynamics._lagrangian_radii(p_i, m_i, c_i, (0.5,))[1]
    s = snaps[6]
    b = bsets[6].stars
    mine = [c for c in b if s.name[c] ≤ 400]
    c_g = Nbody6Dynamics.shrinking_sphere_centre(Float64.(s.pos[:, b]), Float64.(s.mass[b])).centre
    p_i, m_i = Float64.(s.pos[:, mine]), Float64.(s.mass[mine])
    @test at(6, 1).r_h ≈ Nbody6Dynamics._lagrangian_radii(p_i, m_i, c_g, (0.5,))[1]
    @test [r.regime for r in progenitor_budget(snaps, ranges, bsets)] == [r.regime for r in rows]

    # Fly-by: contact, no core coalescence, every remnant quantity masked
    fly, ranges_f = _ep_pair_series([4.0, 2.0, 0.6, 2.0, 4.0, 6.0])
    bsets_f = [bound_set(s) for s in fly]
    ep_f = merger_epochs(fly, ranges_f, bsets_f)
    @test ep_f.contact_index == 3 && ep_f.core_index == 0 && isnan(ep_f.t_core)
    @test ep_f.regime ==
          [:detached, :detached, :interacting, :interacting, :interacting, :interacting]
    @test isempty(ep_f.event_join) && ep_f.pair_core == [0] && ep_f.coalescence_index == 0
    @test all(isnan, ep_f.l_orb)
    diag_f = remnant_diagnostics(fly, ranges_f; settings = _EP_SETTINGS, bound_sets = bsets_f)
    @test all(isnan, diag_f.r_core) && all(isnan, diag_f.lambda_r) && all(isnan, diag_f.r_half)
    @test isempty(diag_f.profile.radius) && isnan(diag_f.t_rc_min)
    rows_f = progenitor_budget(fly, ranges_f, bsets_f; epochs = ep_f)
    @test all(r -> isnan(r.aplus2) && isfinite(r.r_h), rows_f)

    # Files
    mktempdir() do dir
        path = write_remnant_diagnostics(joinpath(dir, "remnant.csv"), diag)
        lines = readlines(path)
        @test occursin(
            " t_contact_nb=2.0 t_contact_myr=6.0 t_core_nb=4.0 t_core_myr=12.0",
            lines[1],
        )
        header = split(lines[2], ",")
        @test header[end] == "regime" && length(header) == 31
        @test [split(l, ",")[end] for l in lines[3:end]] == String.(ep.regime)
        @test all(l -> split(l, ",")[6] == "", lines[3:6])   # r_core masked

        path = write_progenitor_budget(joinpath(dir, "budget.csv"), rows)
        r = _ep_rows(path)
        @test r[1][end] == "regime" && length(r[1]) == 20 && length(r) == 15
        @test [x[end] for x in r[2:2:end]] == String.(ep.regime)
        @test r[2][findfirst(==("aplus2"), r[1])] == ""

        path = write_merger_epochs(joinpath(dir, "merger_epochs.csv"), ep)
        lines = readlines(path)
        @test startswith(lines[1], "# n_clusters=2 n_snapshots=7 merger_min_members=20")
        @test occursin("merger_core_delta=2.0", lines[1])
        r = _ep_rows(path)
        @test length(r) == 6 && length(r[1]) == 20
        col(name) = findfirst(==(name), r[1])
        sys = r[2]
        @test sys[1] == "system" && sys[2] == "1;2"
        @test parse(Int, sys[col("contact_index")]) == ep.contact_index
        @test parse(Float64, sys[col("t_core_nb")]) == ep.t_core
        @test parse(Int, sys[col("settled_index")]) == ep.settled_index
        @test parse(Int, sys[col("join_index")]) == ep.coalescence_index
        ev = r[3]
        @test ev[1] == "event" && ev[2] == "1" && ev[3] == "2"
        @test parse(Int, ev[col("core_index")]) == 5 && parse(Int, ev[col("main_progenitor")]) == 1
        @test parse(Float64, ev[col("mass_ratio")]) ≈ ep.mass_ratio[1]
        @test r[4][1] == "pair" && parse(Int, r[4][col("contact_index")]) == 3
        @test [x[1] for x in r[5:6]] == ["tag", "tag"]
        @test all(
            x -> parse(Int, x[col("dissolved_index")]) == 0 && x[col("t_dissolved_nb")] == "",
            r[5:6],
        )

        path = write_pair_separation(joinpath(dir, "pair_separation.csv"), ep)
        r = _ep_rows(path)
        @test r[1] ==
              ["time_nb", "time_myr", "pair", "separation", "threshold", "centre_offset", "regime"]
        @test length(r) == 8 && all(x -> x[3] == "1;2", r[2:end])
        @test [x[end] for x in r[2:end]] == String.(ep.pair_regime[1, :])
        for (k, x) in enumerate(r[2:end])
            if k < 5
                @test parse(Float64, x[4]) ≈ co.separation[1, k] rtol = 1e-10
                @test parse(Float64, x[5]) ≈ co.threshold[1, k] rtol = 1e-10
                @test x[6] == ""
            else
                @test x[4] == "" && x[5] == ""
                @test parse(Float64, x[6]) ≈ co.separation[1, k] rtol = 1e-10
            end
        end
    end
end

@testset "Merger epochs: a shredded secondary" begin
    snaps, ranges = _ep_shredded_series()
    bsets = [bound_set(s) for s in snaps]
    ep = merger_epochs(snaps, ranges, bsets; settings = _EP_SETTINGS)
    @test ep.tag_dissolved == [0, 4]
    @test isapprox(ep.self_bound_mass[2, 3], 0.25; rtol = 1e-6)
    @test all(iszero, ep.self_bound_mass[2, 4:end])
    @test all(ep.bound_mass[2, 4:end] .> 0.1)   # the debris stays bound to the system
    @test ep.contact_index == 3 && ep.core_index == 6 && ep.event_core == [6]
    @test ep.tag_regime[2, :] == [
        :detached,
        :detached,
        :interacting,
        :dissolved,
        :dissolved,
        :coalescing,
        :coalescing,
        :settled,
    ]
    @test ep.tag_regime[1, 4:5] == [:interacting, :interacting]
    @test ep.regime[4:5] == [:interacting, :interacting] && ep.regime[6] === :coalescing
    @test ep.pair_regime[1, 4:5] == [:dissolved, :dissolved]
    @test ep.main_progenitor == [1]
    @test ep.mass_ratio[1] ≈ 1 / 3 rtol = 1e-6

    tree = merger_tree(snaps, ranges; bound_sets = bsets, epochs = ep)
    e = only(tree.events)
    @test e.index == 6 && e.t_contact == 2.0 && e.t_contact_myr == 6.0
    @test e.t_core == 5.0 && e.t_core_myr == 15.0
    @test e.main_progenitor == 1
    @test e.mass_ratio ≈ 1 / 3 rtol = 1e-6
    # Epochs of another tree (no cluster present, hence no event) are refused
    co = coalescence_time(snaps, ranges; bound_sets = bsets)
    empty_co = merge(co, (; present = falses(size(co.present))))
    @test_throws ArgumentError merger_tree(
        snaps,
        ranges;
        bound_sets = bsets,
        coalescence = co,
        epochs = merger_epochs(snaps, ranges, bsets, empty_co),
    )

    rows = progenitor_budget(snaps, ranges, bsets; settings = _EP_SETTINGS, epochs = ep)
    sec = rows[2:2:end]
    @test [r.regime for r in sec] == ep.tag_regime[2, :]
    @test all(r -> isfinite(r.r_h) && isnan(r.aplus2), sec[4:5])
    @test all(r -> isfinite(r.aplus2), sec[6:8])
    @test all(
        r -> isapprox(r.m_initial, r.m_present + r.m_removed + r.m_stellar_loss; atol = 1e-12),
        rows,
    )
end

@testset "Merger epochs: three clusters with a trapped satellite" begin
    snaps, ranges = _ep_satellite_series()
    bsets = [bound_set(s) for s in snaps]
    co = coalescence_time(snaps, ranges; bound_sets = bsets)
    ep = merger_epochs(snaps, ranges, bsets, co; settings = _EP_SETTINGS)
    @test co.index === nothing && ep.coalescence_index == 0 && isnan(ep.t_coal)
    @test ep.contact_index == 2 && ep.core_index == 3 && ep.t_core == 2.0
    @test ep.regime == [:detached, :interacting, :coalescing, :settled, :settled, :settled]
    @test ep.tag_regime[3, :] ==
          [:detached, :detached, :satellite, :satellite, :satellite, :satellite]
    @test ep.tag_regime[1, :] == ep.tag_regime[2, :] == ep.regime
    @test ep.pairs == [(1, 2), (1, 3), (2, 3)]
    @test ep.pair_core == [3, 0, 0] && ep.pair_join == [3, 0, 0]
    @test ep.pair_regime[2, 3:end] == fill(:satellite, 4) && ep.pair_regime[1, 3] === :coalescing
    @test all(ep.delta_s[2:3, 3:end] .> 2)   # the satellite keeps its own core
    @test ep.tag_dissolved == [0, 0, 0]

    tree = merger_tree(snaps, ranges; bound_sets = bsets, coalescence = co, epochs = ep)
    @test length(tree.events) == 1 && tree.joined == [3, 0, 0]
    @test tree.events[1].group_a == [1] && tree.events[1].group_b == [2]
    @test tree.events[1].t_core == 2.0 && tree.events[1].mass_ratio ≈ 1.0

    rows = progenitor_budget(snaps, ranges, bsets; settings = _EP_SETTINGS, epochs = ep)
    sat = rows[3:3:end]
    @test all(r -> r.regime === :satellite && isnan(r.aplus2), sat[3:end])
    # The satellite keeps its own centre: r_h of its compact core
    @test all(r -> r.r_h < 0.1, sat)
    @test all(r -> isfinite(r.aplus2), rows[[3 * 3 + 1, 3 * 3 + 2]])

    # Validity: the epochs in [stages], the pair (1, 2) window ends at its core coalescence
    ints = snapshot_integrals(snaps)
    remn = remnant_diagnostics(
        snaps,
        ranges;
        settings = _EP_SETTINGS,
        bound_sets = [s.bound for s in ints],
        coalescence = co,
        epochs = ep,
    )
    rec = validity_record(
        Dict{Symbol,Any}(
            :snapshots => snaps,
            :integrals => ints,
            :remnant => remn,
            :merger_epochs => ep,
            :progenitor_budget => rows,
        ),
        _EP_SETTINGS;
        calibration = nothing,
    )
    st = rec.sections["stages"]
    @test isnan(st["coalescence_time_nb"]) && st["t_contact_nb"] == 1.0 && st["t_core_nb"] == 2.0
    windows = filter(e -> e.kind === :coalescence_phase, rec.epochs)
    w12 = only(filter(e -> occursin("clusters 1–2", e.detail), windows))
    @test w12.t_end == 2.0
    @test all(e -> e.t_end ≤ 5.0, windows)
    @test [s.regime for s in rec.snapshots] == ep.regime
    mktempdir() do dir
        path = write_validity(joinpath(dir, "validity.toml"), rec)
        d = Nbody6Dynamics.TOML.parsefile(path)
        @test d["stages"]["t_core_nb"] == 2.0 && d["stages"]["t_contact_nb"] == 1.0
        csv = write_validity_snapshots(joinpath(dir, "validity_snapshots.csv"), rec)
        r = _ep_rows(csv)
        @test r[1][end] == "regime" && [x[end] for x in r[2:end]] == String.(ep.regime)
        path = write_pair_separation(joinpath(dir, "pair_separation.csv"), ep)
        r = _ep_rows(path)
        @test length(r) == 1 + 6 * 3
        @test [x[3] for x in r[2:4]] == ["1;2", "1;3", "2;3"]
        @test all(x -> x[end] == "satellite" && x[4] != "", filter(x -> x[3] == "1;3", r)[3:end])
        path = write_merger_epochs(joinpath(dir, "merger_epochs.csv"), ep)
        r = _ep_rows(path)
        @test [x[1] for x in r[2:end]] == ["system", "event", "pair", "pair", "pair", "tag", "tag", "tag"]
        @test r[2][findfirst(==("join_index"), r[1])] == "0"
    end
end

@testset "Merger epochs: configuration keys" begin
    pp = PostprocessConfig()
    @test pp.merger_min_members == 20 && pp.merger_contact_lagrangian == 0.9
    @test pp.merger_core_delta == 2.0 && pp.merger_core_n == 200
    @test pp.merger_dissolved_fraction == 0.1 && pp.merger_coalescing_tcr == 5.0
    for (key, bad) in (
        (:merger_min_members, 1),
        (:merger_contact_lagrangian, 0.0),
        (:merger_contact_lagrangian, 1.0),
        (:merger_core_delta, 0.0),
        (:merger_core_n, 10),
        (:merger_dissolved_fraction, 0.0),
        (:merger_dissolved_fraction, 1.0),
        (:merger_coalescing_tcr, -1.0),
    )
        cfg_bad = PostprocessConfig(; (key => bad,)...)
        @test_throws Regex("postprocess\\.$key") Nbody6Dynamics._validate_estimators(cfg_bad)
    end
    @test Nbody6Dynamics._validate_estimators(
        PostprocessConfig(;
            merger_coalescing_tcr = 0.0,
            merger_min_members = 50,
            merger_core_n = 50,
        ),
    ) === nothing
end

@testset "Merger epochs: pipeline products" begin
    snaps, ranges = _ep_pair_series([4.0, 2.0, 0.6, 0.3, 0.02, 0.02, 0.02])
    base = mktempdir()
    run_dir = joinpath(base, "runs", "epochs")
    out = joinpath(run_dir, "output")
    mkpath(out)
    for s in snaps
        open(joinpath(out, "conf.3_$(time_nb(s.header))"), "w") do io
            _write_fortran_record(io, Int32[nparticles(s), 1, 1, 20])
            _write_fortran_record(io, s.header.params)
            for i in 1:nparticles(s)
                buf = IOBuffer()
                write(buf, s.mass[i], s.pos[:, i]..., s.vel[:, i]..., s.name[i])
                _write_fortran_record(io, take!(buf))
            end
        end
    end
    write(
        joinpath(out, "merger_summary.txt"),
        "  Cluster 1: plummer, imf=equal, N=400 (after trunc: 400), M=1 M☉\n" *
        "  Cluster 2: plummer, imf=equal, N=400 (after trunc: 400), M=1 M☉\n",
    )
    pp = PostprocessConfig(;
        parse_stdout = false,
        read_lagr = false,
        read_escapers = false,
        read_stellar_evo = false,
        read_binary_evo = false,
        mass_budget = false,
        energy_events = false,
        msr_n_tracer = 20,
        msr_n_random = 20,
        rotation_n_min = 100,
    )
    cfg = Nbody6Config(
        InstallConfig(; enabled = false),
        BuildConfig(),
        SimulationConfig(; run_test = false, runs_dir = "runs"),
        pp,
        VisualizationConfig(; enabled = false),
        MergerPipelineConfig(),
    )
    results = Nbody6Dynamics._finish_pipeline(cfg, run_dir, base, nothing, String[], time())
    ep = results[:merger_epochs]
    @test ep.core_index == 5 && ep.contact_index == 3
    @test results[:coalescence].index == 5
    @test results[:remnant].t_core == ep.t_core && results[:remnant].regime == ep.regime
    @test results[:merger_tree].events[1].t_core == ep.t_core
    @test results[:validity].sections["stages"]["t_core_nb"] == ep.t_core
    for f in (
        "merger_epochs.csv",
        "pair_separation.csv",
        "remnant_diagnostics.csv",
        "progenitor_budget.csv",
        "merger_tree.csv",
        "validity_snapshots.csv",
    )
        @test isfile(joinpath(run_dir, f))
    end
    @test _ep_rows(joinpath(run_dir, "merger_tree.csv"))[1] == [
        "index",
        "time_nb",
        "time_myr",
        "group_a",
        "group_b",
        "separation",
        "threshold",
        "t_contact_nb",
        "t_contact_myr",
        "t_core_nb",
        "t_core_myr",
        "main_progenitor",
        "mass_ratio",
    ]
    @test _ep_rows(joinpath(run_dir, "progenitor_budget.csv"))[end][end] == "settled"
end
