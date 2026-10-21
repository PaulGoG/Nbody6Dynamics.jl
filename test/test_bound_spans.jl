# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# WP10 10.2 increment 2: physically bound pairs at a snapshot (D3) and their
# spans over the snapshot sequence (D4), on synthetic snapshots in N-body
# units (G = 1).

@testset "Bound pairs and spans (level P)" begin
    nd = Nbody6Dynamics

    # A snapshot from bodies (name, mass, x, v); the first 2·npairs bodies are KS components
    function snap(t, bodies; npairs = 0)
        n = length(bodies)
        params = zeros(Float32, 20)
        params[1] = t
        params[3] = 1.0
        params[4] = 1.0
        hdr = SnapshotHeader(Int32(n), Int32(0), Int32(1), Int32(20), params)
        name = Int32[b[1] for b in bodies]
        mass = Float32[b[2] for b in bodies]
        pos = Float32.(reduce(hcat, [b[3] for b in bodies]))
        vel = Float32.(reduce(hcat, [b[4] for b in bodies]))
        return Snapshot(hdr, name, mass, pos, vel, Float32[], Float32[], npairs)
    end
    # Six distant field stars at rest
    field = [
        (100 + k, 0.01, 2.0 .* s .* (k ≤ 3 ? [1.0, 0, 0] : [0, 1.0, 0]) .+ [0, 0, 0.3k], zeros(3)) for (k, s) in enumerate((1, -1, 1, -1, 1, -1))
    ]
    # A circular pair of names n1, n2, masses m, semi-major axis a, centre c, in the x-y plane, L along +z
    function pair(n1, n2; m = 0.01, a = 1.0e-3, c = zeros(3), e = 0.0)
        M = 2m
        r = a * (1 + e)                      # placed at apocentre
        v = sqrt(M * (1 - e) / (a * (1 + e)))
        return [
            (n1, m, c .+ [-r / 2, 0, 0], [0, -v / 2, 0]),
            (n2, m, c .+ [r / 2, 0, 0], [0, v / 2, 0]),
        ]
    end
    flyby = (9, 0.1, [0, 4.0e-3, 0], [10.0, 0, 0])   # γ_apo ≈ 0.16, unbound to the pair
    apart = [(1, 0.01, [-0.25, 0, 0], [0, -1.0, 0]), (2, 0.01, [0.25, 0, 0], [0, 1.0, 0])]

    @testset "Bound pairs at a snapshot" begin
        s = snap(0.0, vcat(pair(1, 2), field))
        bp = nd.bound_pairs(s)
        p = only(filter(q -> (q.name1, q.name2) == (1, 2), bp))
        @test p.a ≈ 1.0e-3 rtol = 1e-4
        @test p.e < 1e-3 && p.gamma_apo < 1e-6 && !p.ks && !p.hierarchy
        # A close fly-by fails the perturbation test
        @test isempty(
            filter(
                q -> (q.name1, q.name2) == (1, 2),
                nd.bound_pairs(snap(0.0, vcat(pair(1, 2), [flyby], field))),
            ),
        )
        # A KS pair is a member even when its snapshot state is unbound
        ks = snap(
            0.0,
            vcat(
                [(11, 0.01, [0, 0, 1.0], [0, 0, 0]), (12, 0.01, [1.0e-3, 0, 1.0], [0, 50.0, 0])],
                field,
            );
            npairs = 1,
        )
        q = only(filter(q -> (q.name1, q.name2) == (11, 12), nd.bound_pairs(ks)))
        @test q.ks && isnan(q.a)
        # An eccentric inner pair with a stable coplanar outer body: kept as a hierarchy despite γ_apo > 0.1
        inner = pair(21, 22; e = 0.8)
        vout = sqrt(0.04 / 4.0e-3)
        outer = (23, 0.02, [0, 4.0e-3, 0], [-vout, 0, 0])
        h = only(
            filter(
                q -> (q.name1, q.name2) == (21, 22),
                nd.bound_pairs(snap(0.0, vcat(inner, [outer], field))),
            ),
        )
        @test h.hierarchy && h.outer_name == 23 && h.gamma_apo > 0.1
        # The same body inside the stability limit: excluded
        close = (23, 0.02, [0, 3.0e-3, 0], [-sqrt(0.04 / 3.0e-3), 0, 0])
        @test isempty(
            filter(
                q -> (q.name1, q.name2) == (21, 22),
                nd.bound_pairs(snap(0.0, vcat(inner, [close], field))),
            ),
        )
    end

    @testset "Spans over snapshots" begin
        far = (9, 0.1, [0, 1.5, 1.5], [0, 0, 0])
        seq(perturbed) = [
            snap(
                Float64(j - 1),
                vcat(j ≤ 6 ? pair(1, 2) : apart, [j in perturbed ? flyby : far], field),
            ) for j in 1:8
        ]
        # One excused miss: one span
        sp = filter(s -> (s.name1, s.name2) == (1, 2), nd.bound_spans(seq((4,)); nbin0 = 1))
        s = only(sp)
        @test s.primordial &&
              s.j_first == 1 &&
              s.j_last == 6 &&
              s.n_missed == 1 &&
              s.n_snapshots == 5
        @test s.left_censored && !s.right_censored && s.end_unbound
        @test s.t_end_lo == 5.0 && s.t_end_hi == 6.0 && s.t_form_hi == 0.0
        @test s.a_first ≈ 1.0e-3 rtol = 1e-4
        # Two consecutive misses: two spans
        sp2 = filter(s -> (s.name1, s.name2) == (1, 2), nd.bound_spans(seq((4, 5)); nbin0 = 1))
        @test length(sp2) == 2 && sp2[1].j_last == 3 && sp2[2].j_first == 6
        @test sp2[2].t_form_lo == 4.0 && sp2[2].t_form_hi == 5.0
        # n_miss = 0: the single miss splits the span too
        @test length(
            filter(
                s -> (s.name1, s.name2) == (1, 2),
                nd.bound_spans(seq((4,)); nbin0 = 1, n_miss = 0),
            ),
        ) == 2
        # An exchange: 2 leaves, 3 takes its place; the old span ends, a new one begins
        ex = [
            snap(
                Float64(j - 1),
                j ≤ 3 ? vcat(pair(1, 2), [(3, 0.01, [0, 1.2, -1.2], zeros(3))], field) :
                vcat([(2, 0.01, [0, 1.2, -1.2], zeros(3))], pair(1, 3), field),
            ) for j in 1:6
        ]
        se = nd.bound_spans(ex; nbin0 = 0)
        @test only(filter(s -> (s.name1, s.name2) == (1, 2), se)).j_last == 3
        s13 = only(filter(s -> (s.name1, s.name2) == (1, 3), se))
        @test s13.j_first == 4 && s13.right_censored && !s13.primordial
        @test issorted(se; by = s -> (s.t_form_hi, s.name1, s.name2))
        @test_throws ArgumentError nd.bound_spans(ex; nbin0 = -1)
    end
end
