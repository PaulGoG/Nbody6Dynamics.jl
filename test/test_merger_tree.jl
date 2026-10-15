# Merger tree and per-progenitor budget on a synthetic three-cluster series (uses `_am_snapshot`
# of test_angular_momentum.jl).

using StableRNGs: StableRNG
using LinearAlgebra: ×

function _three_cluster_series(rng; n_per = 300, a = 0.1)
    clumps = [Nbody6Dynamics.sample_plummer(n_per, a; rng) for _ in 1:3]
    # Internal velocities of a clump of mass 1/3
    for (p, v) in clumps
        v .*= sqrt(1 / 3)
    end
    offsets = [
        ([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]),
        ([0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]),
        ([0.0, 0.0, 0.0], [0.02, 0.0, 0.0], [0.0, 3.0, 0.0]),
        ([0.0, 0.0, 0.0], [0.02, 0.0, 0.0], [0.0, 3.0, 0.0]),
        ([0.0, 0.0, 0.0], [0.02, 0.0, 0.0], [0.0, 0.02, 0.0]),
        ([0.0, 0.0, 0.0], [0.02, 0.0, 0.0], [0.0, 0.02, 0.0]),
    ]
    mass = fill(1 / (3 * n_per), 3 * n_per)
    snaps = Snapshot[]
    for (k, off) in enumerate(offsets)
        pos = hcat((clumps[i][1] .+ off[i] for i in 1:3)...)
        vel = hcat((clumps[i][2] for i in 1:3)...)
        push!(snaps, _am_snapshot(mass, pos, vel; t = Float64(k - 1)))
    end
    ranges = [1:n_per, (n_per + 1):(2 * n_per), (2 * n_per + 1):(3 * n_per)]
    return snaps, ranges
end

@testset "merger_tree on a three-cluster sequence" begin
    rng = StableRNG(3301)
    snaps, ranges = _three_cluster_series(rng)
    bsets = [bound_set(s) for s in snaps]
    tree = merger_tree(snaps, ranges; bound_sets = bsets)
    @test tree isa MergerTree
    @test tree.pairs == [(1, 2), (1, 3), (2, 3)]
    @test tree.n_distinct == [3, 3, 2, 2, 1, 1]
    @test tree.joined == [3, 5, 5]
    @test tree.dissolved == [0, 0, 0]
    @test all(tree.present)
    @test size(tree.separation) == (3, 6) && all(isfinite, tree.separation)
    @test all(tree.separation[1, 1:2] .> tree.threshold[1, 1:2])
    @test all(tree.separation[1, 3:6] .< tree.threshold[1, 3:6])
    @test length(tree.events) == 2
    e1, e2 = tree.events
    @test e1.index == 3 && e1.group_a == [1] && e1.group_b == [2] && e1.time_nb == 2.0
    @test e2.index == 5 && e2.group_a == [1, 2] && e2.group_b == [3] && e2.time_nb == 4.0
    @test e1.separation == tree.separation[1, 3] && e1.threshold == tree.threshold[1, 3]
    @test e2.separation == tree.separation[2, 5] && e2.threshold == tree.threshold[2, 5]
    # Without bound sets the same tree follows from the snapshots alone
    same = merger_tree(snaps, ranges)
    @test same.joined == tree.joined && length(same.events) == 2
    @test_throws ArgumentError merger_tree(Snapshot[], ranges)
    @test_throws ArgumentError merger_tree(snaps, ranges[1:1])

    mktempdir() do dir
        path = joinpath(dir, "merger_tree.csv")
        @test write_merger_tree(path, tree) == path
        lines = readlines(path)
        @test lines[1] == "# n_clusters=3 n_snapshots=6 dissolved=0;0;0"
        @test lines[2] == "index,time_nb,time_myr,group_a,group_b,separation,threshold"
        @test length(lines) == 4
        c = split(lines[4], ",")
        @test c[1] == "5" && c[4] == "1;2" && c[5] == "3"
    end
end

@testset "progenitor_budget mass identity and mixing" begin
    rng = StableRNG(3302)
    snaps, ranges = _three_cluster_series(rng)
    # Snapshot 4 onward: 20 stars of cluster 3 removed; snapshot 5 onward: 30 stars of
    # cluster 1 lose 20 % of their mass
    removed = 601:620
    for k in 4:6
        s = snaps[k]
        keep = [i for i in 1:nparticles(s) if !(s.name[i] in removed)]
        m = Float64.(s.mass[keep])
        if k >= 5
            for (q, i) in enumerate(keep)
                s.name[i] <= 30 && (m[q] *= 0.8)
            end
        end
        snaps[k] = _am_snapshot(
            m,
            Float64.(s.pos[:, keep]),
            Float64.(s.vel[:, keep]);
            t = Float64(k - 1),
            names = s.name[keep],
        )
    end
    bsets = [bound_set(s) for s in snaps]
    rows = progenitor_budget(snaps, ranges, bsets)
    @test length(rows) == 18
    @test [r.tag for r in rows[1:3]] == [1, 2, 3] && rows[4].time_nb == 1.0
    for r in rows
        @test r.n_initial == 300
        @test r.m_initial ≈ r.m_present + r.m_removed + r.m_stellar_loss atol = 1e-12
        @test r.n_bound <= r.n_present <= r.n_initial
    end
    @test all(r -> r.n_bound >= 0.95 * r.n_present, rows)
    at(k, tag) = rows[3 * (k - 1) + tag]
    @test at(1, 3).m_removed == 0 && at(4, 3).n_present == 280
    @test at(4, 3).m_removed ≈ 20 / 900 atol = 1e-9
    1e-9
    @test at(4, 1).m_stellar_loss == 0
    @test at(5, 1).m_stellar_loss ≈ 30 * 0.2 / 900 atol = 1e-9
    @test at(6, 2).m_stellar_loss == 0 && at(6, 2).m_removed == 0
    # Separated clusters sit at different radii from the global centre; the overlapping
    # product is mixed to the noise floor
    @test all(isfinite, (at(1, 2).r_h, at(1, 2).centre_offset, at(1, 2).d_radius, at(1, 2).d_floor))
    @test at(1, 2).d_radius > 3 * at(1, 2).d_floor
    @test at(6, 3).d_radius < 3 * at(6, 3).d_floor && at(6, 3).centre_offset < 0.1
    @test abs(at(6, 3).aplus2) < 4 * at(6, 3).aplus2_err + 0.05
    @test isnan(at(1, 1).d_energy)   # no potential column in the synthetic snapshots
    @test_throws ArgumentError progenitor_budget(Snapshot[], ranges, BoundSet[])
    @test_throws ArgumentError progenitor_budget(snaps, ranges, bsets[1:2])
    @test_throws ArgumentError progenitor_budget(snaps, ranges, bsets; min_members = 1)

    mktempdir() do dir
        path = joinpath(dir, "progenitor_budget.csv")
        @test write_progenitor_budget(path, rows) == path
        lines = readlines(path)
        @test length(lines) == 19
        @test split(lines[1], ",")[1:4] == ["time_nb", "time_myr", "tag", "n_initial"]
        @test all(l -> length(split(l, ",")) == 19, lines[2:end])
    end
end
