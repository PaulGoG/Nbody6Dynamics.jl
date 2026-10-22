# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Regression tests of the post-processing corrections: merged hierarchies at
# levels P and E, epoch matching of the stellar types,
# the q-restricted fraction, the exchange end class, chain episodes, the census
# split, restart duplicates, the bev.82 pair energies, the ADJUST reader and
# the run-directory guard of run_pipeline.

@testset "Analysis corrections" begin
    nd = Nbody6Dynamics

    # A snapshot from bodies (name, mass, x, v); the first 2·npairs bodies are KS components
    function snap(t, bodies; npairs = 0, pair_name = Int32[])
        n = length(bodies)
        params = zeros(Float32, 20)
        params[1] = t
        params[3] = 1.0
        params[4] = 1.0
        params[11] = 1.0
        params[12] = 1.0
        hdr = SnapshotHeader(Int32(n), Int32(0), Int32(1), Int32(20), params)
        name = Int32[b[1] for b in bodies]
        mass = Float32[b[2] for b in bodies]
        pos = Float32.(reduce(hcat, [b[3] for b in bodies]))
        vel = Float32.(reduce(hcat, [b[4] for b in bodies]))
        return Snapshot(hdr, name, mass, pos, vel, Float32[], Float32[], npairs, Int32.(pair_name))
    end
    # A record with the given fields and the others absent (NaN, 0, false, :none)
    function filled(::Type{T}; kw...) where {T}
        values = map(fieldnames(T), fieldtypes(T)) do f, F
            haskey(kw, f) ? kw[f] :
            F <: AbstractFloat ? F(NaN) :
            F <: Bool ? false :
            F <: Integer ? zero(F) :
            F === Symbol ? :none : error("no absent value for field $f::$F")
        end
        return T(values...)
    end
    field = [
        (100 + k, 0.01, 2.0 .* s .* (k ≤ 3 ? [1.0, 0, 0] : [0, 1.0, 0]) .+ [0, 0, 0.3k], zeros(3)) for (k, s) in enumerate((1, -1, 1, -1, 1, -1))
    ]
    # Merged hierarchy in KS slot 1: star 1 (21) first inner component, star 2 (23)
    # the outer body, star 3 (22) the second inner component among the singles
    v = sqrt(0.02 / 1.0e-3)
    merged_bodies = vcat(
        [
            (21, 0.01, [-5.0e-4, 0, 0], [0, -v / 2, 0]),
            (23, 0.01, [0, 0.5, 0], zeros(3)),
            (22, 0.01, [5.0e-4, 0, 0], [0, v / 2, 0]),
        ],
        field,
    )
    merged = snap(0.0, merged_bodies; npairs = 1, pair_name = [-30])

    @testset "1. Merged hierarchies at level P" begin
        bp = nd.bound_pairs(merged)
        inner = only(filter(q -> (q.name1, q.name2) == (21, 22), bp))
        @test !inner.ks
        @test inner.a ≈ 1.0e-3 rtol = 1e-3
        @test isempty(filter(q -> (q.name1, q.name2) == (21, 23), bp))
        @test !any(q -> q.ks, bp)
        # The same slot with a positive name stays a forced KS member
        reduced = snap(0.0, merged_bodies; npairs = 1, pair_name = [30])
        @test only(filter(q -> (q.name1, q.name2) == (21, 23), nd.bound_pairs(reduced))).ks
    end

    @testset "2. Merged hierarchies at level E (population series)" begin
        comps, hier = nd._bpop_components(merged, :ks, 10, 4, 0.1)
        @test comps == [(1, 3)] && hier == [true]
        comps_pos, hier_pos = nd._bpop_components(
            snap(0.0, merged_bodies; npairs = 1, pair_name = [30]),
            :ks,
            10,
            4,
            0.1,
        )
        @test comps_pos == [(1, 2)] && hier_pos == [false]
        rows = binary_population_series([merged]; levels = (:ks,), centre = :origin)
        r = only(filter(r -> r.bin == 0 && r.class === :all, rows))
        @test r.n_binaries == 1 && r.n_hierarchies == 1
        @test r.n_singles == 7          # the outer body and the six field stars
        # A KS slot with a massless component is no binary
        ghost = copy(merged_bodies)
        ghost[2] = (23, 0.0, [0, 0.5, 0], zeros(3))
        comps_g, _ =
            nd._bpop_components(snap(0.0, ghost; npairs = 1, pair_name = [30]), :ks, 10, 4, 0.1)
        @test isempty(comps_g)
    end

    @testset "3. Epoch matching of the stellar types" begin
        sev(t) =
            StellarEvolutionSnapshot(t, 1, [StellarRecord(t, 1, 101, 1, 0.0, 1.0, 0.0, 0.0, 3.7)])
        epochs = [sev(1.0), sev(2.0)]
        @test nd._bpop_nearest_epoch(epochs, 2.0) == (2, 0.0)
        i, lag = nd._bpop_nearest_epoch(epochs, 2.5)
        @test i == 2 && lag ≈ 0.5
        i0, lag0 = nd._bpop_nearest_epoch(epochs, 0.5)
        @test i0 === nothing && isnan(lag0)
        s = snap(2.5, [(101, 0.01, zeros(3), zeros(3))])
        maps = @test_logs (:info, r"largest lag") match_mode = :any nd._bpop_pipeline_kstar(
            [s],
            epochs,
            BinaryEvolutionSnapshot[],
        )
        @test maps !== nothing && maps[1][101] == 1
    end

    @testset "4. q-restricted fraction" begin
        bodies = vcat([(1, 0.6, [0.0, 0, 0], zeros(3)), (2, 0.5, [1.0e-3, 0, 0], zeros(3))], field)
        s = snap(0.0, bodies; npairs = 1, pair_name = [40])
        fq(kstar) = only(
            filter(
                r -> r.bin == 0 && r.class === :all,
                binary_population_series([s]; levels = (:ks,), centre = :origin, kstar),
            ),
        ).f_b_q
        types(k2) = [Dict{Int,Int}(1 => 1, 2 => k2, (100 + k => 1 for k in 1:6)...)]
        @test fq(nothing) == 1.0          # without types every binary with q > q_min
        @test fq(types(1)) == 1.0         # two main-sequence components
        @test fq(types(11)) == 0.0        # a white-dwarf companion is no photometric binary
    end

    # Level-E records
    days(p_nb) = p_nb * 3.6524e8
    line = Ref(0)
    function ks(kind, t, n1, n2; p_nb = 0.3)
        line[] += 1
        return KSEvent(;
            line = line[],
            interval = 1,
            time_nb = t,
            kind = kind,
            name1 = n1,
            name2 = n2,
            name_cm = 10_000 + min(n1, n2),
            kstar1 = 0,
            kstar2 = 0,
            kstar_cm = 0,
            semi_nb = 1.0e-3,
            ecc = 0.5,
            eb_nb = -1.0e-4,
            period_d = days(p_nb),
            m1_msun = 1.0,
            m2_msun = 0.5,
        )
    end
    function eventlog(; ks = KSEvent[], fewbody = FewBodyEvent[], exchange = ExchangeEvent[])
        log = BinaryEventLog()
        append!(log.ks, ks)
        append!(log.fewbody, fewbody)
        append!(log.exchange, exchange)
        return log
    end

    @testset "7. Exchange end class" begin
        xch(t) = ExchangeEvent(;
            line = 70,
            interval = 1,
            time_nb = t,
            kind = :exchange,
            name1_initial = 11,
            name2_initial = 12,
            name1_final = 11,
            name2_final = 13,
        )
        spans_at(t_x) = binary_spans(
            eventlog(;
                ks = [ks(:ks_new, 1.0, 11, 12), ks(:ks_end, 3.0, 11, 12), ks(:ks_new, t_x, 11, 13)],
                exchange = [xch(t_x)],
            );
            nzero = 10_000,
            nbin0 = 0,
            tstar_myr = 1.0,
            t_end_nb = 10.0,
        )
        old(s) = only(filter(x -> (x.name1, x.name2) == (11, 12), s))
        @test old(spans_at(3.0)).end_class === :exchange
        # An exchange long after the close does not relabel the end
        @test old(spans_at(5.0)).end_class === :unresolved
    end

    @testset "8. Chain episodes" begin
        fb(l, t, kind, names) = FewBodyEvent(;
            line = l,
            interval = 1,
            time_nb = t,
            kind = kind,
            n_members = 3,
            names = names,
        )
        # A second chain of other members closes the first; the last runs to the end of the log
        log = eventlog(;
            ks = [ks(:ks_new, 4.0, 1, 2)],
            fewbody = [
                fb(10, 1.0, :chain_new, (5, 6, 7, 0, 0, 0)),
                fb(20, 2.0, :chain_new, (8, 9, 10, 0, 0, 0)),
            ],
        )
        chains, _ = nd._span_chains(log, 10_000)
        @test length(chains) == 2
        @test chains[1].t_end == 2.0 && chains[2].t_start == 2.0 && chains[2].t_end == 4.0
        # A closing record without a time stamp closes at the time of the preceding record
        log2 = eventlog(;
            fewbody = [
                fb(10, 1.0, :chain_new, (5, 6, 7, 0, 0, 0)),
                fb(11, 1.5, :chain_binary, (5, 6, 0, 0, 0, 0)),
                fb(12, NaN, :chain_end, (5, 6, 7, 0, 0, 0)),
            ],
        )
        chains2, _ = nd._span_chains(log2, 10_000)
        @test only(chains2).t_end == 1.5
    end

    @testset "9. Census split" begin
        P(lo, hi, elo, ehi) = filled(
            BoundSpan;
            name1 = 1,
            name2 = 2,
            t_form_lo = lo,
            t_form_hi = hi,
            t_end_lo = elo,
            t_end_hi = ehi,
        )
        # Member snapshots 0–1 and 3–4, one unexcused miss at 2: two census spans
        c = binary_census(BinarySpan[], [P(NaN, 0.0, 1.0, 2.0), P(2.0, 3.0, 4.0, 5.0)]; nbin0 = 0)
        @test length(c) == 2 && all(s -> s.levels === :P, c)
    end

    @testset "10. Restart duplicates of coll.13 / coal.24" begin
        C(t, n1, n2) = filled(CollisionRecord; time_nb = t, name1 = n1, name2 = n2)
        # Segment 2 restarted from the dump at t = 4 after the first was killed at t ≈ 6
        recs = [C(1.0, 1, 2), C(5.0, 3, 4), C(6.0, 5, 6), C(5.0, 3, 4), C(7.0, 7, 8)]
        kept = nd._span_restart_dedup(recs, [4.0])
        @test [r.time_nb for r in kept] == [1.0, 5.0, 7.0]
        @test (kept[2].name1, kept[2].name2) == (3, 4)
        # Without restart starts only identical (time, names) records are dropped
        @test length(
            nd._span_restart_dedup([C(1.0, 1, 2), C(1.0, 2, 1), C(2.0, 3, 4)], Float64[]),
        ) == 2
        info = Dict{String,Any}(
            "segments" => [
                Dict{String,Any}("index" => 1, "t_start" => 0.0),
                Dict{String,Any}("index" => 2, "t_start" => 4.0),
            ],
        )
        @test nd._span_restart_starts(info) == [4.0]
        @test isempty(nd._span_restart_starts(Dict{String,Any}()))
    end

    @testset "12. bev.82 pair energies of reduced pairs only" begin
        # Slot 1 (names 1, 2) is reduced; slot 2 (3 inner, 4 outer) is a merged
        # hierarchy whose inner binary (3, 5) bev.82 lists
        bodies = [
            (1, 0.1, [0.0, 0, 0], zeros(3)),
            (2, 0.1, [1.0e-3, 0, 0], zeros(3)),
            (3, 0.1, [1.0, 0, 0], zeros(3)),
            (4, 0.1, [1.0, 1.0, 0], zeros(3)),
            (5, 0.1, [1.001, 0, 0], zeros(3)),
            (6, 0.1, [-1.0, 0, 0], zeros(3)),
        ]
        s = snap(1.0, bodies; npairs = 2, pair_name = [50, -60])
        @test nd._budget_reduced_pairs(s) == Set([(1, 2)])
        ints = snapshot_integrals([s]; min_members = 2)
        rec(n1, n2, loga) =
            BinaryRecord(1.0, 1, 2, n1, n2, 1, 1, 1, 0.0, 0.0, 0.0, loga, 1.0, 1.0, zeros(18)...)
        r12 = rec(1, 2, 0.0)
        r35 = rec(3, 5, 1.0)
        bev = BinaryEvolutionSnapshot(1.0, 2, [r12, r35])
        u = UnitScaling(1.0, 1.0, 1.0, 1.0)
        diag = DiagnosticsData(AdjustRecord[], Dict{String,Float64}())
        row = only(
            energy_budget(
                ints,
                EscaperRecord[],
                EventCounters[],
                diag,
                u;
                bevs = [bev],
                snapshots = [s],
            ),
        )
        @test row.pair_source === :bev
        @test row.e_pairs ≈ -binding_energy(r12) rtol = 1e-12
        all_pairs =
            only(energy_budget(ints, EscaperRecord[], EventCounters[], diag, u; bevs = [bev]))
        @test all_pairs.e_pairs ≈ -(binding_energy(r12) + binding_energy(r35)) rtol = 1e-12
        # A missing M* or T* gives NaN escaper terms, not numbers in N-body units
        esc = [EscaperRecord(0.5, 1.0, 1.0, 10.0, 0, 6, 0.0, 0.0)]
        @test isnan(
            mass_budget(
                ints,
                esc,
                EventCounters[],
                EngineEvent[],
                UnitScaling(1.0, NaN, NaN, 1.0),
            )[1].m_escaped,
        )
    end

    @testset "14. ADJUST reader skips unparseable records" begin
        dir = mktempdir()
        path = joinpath(dir, "out1000")
        open(path, "w") do io
            println(io, " ADJUST:    1.0   0.5   0.50  1e-8  -0.250   500   10   1.2")
            println(io, " ADJUST:    2.0   1.0   *****  2e-8  -0.249   498   10   1.3")
            println(io, " ADJUST:    3.0   1.5   0.52  3.0D-8  -0.248   497   10   1.4")
        end
        diag =
            @test_logs (:warn, r"1 ADJUST lines skipped") match_mode = :any read_diagnostics(path)
        @test [a.time_nb for a in diag.adjust] == [1.0, 3.0]
        @test diag.adjust[2].de_rel ≈ 3.0e-8
        @test !any(a -> a.time_nb == 0 && a.n == 0, diag.adjust)
        u = extract_scaling(diag)
        @test isnan(u.zmbar) && isnan(u.tscale) && u.rbar == 1.0
    end

    @testset "17. run_pipeline refuses an existing run directory" begin
        base = mktempdir()
        run_dir = joinpath(base, "runs", "r1")
        mkpath(run_dir)
        touch(joinpath(run_dir, "RUN_INFO.toml"))
        cfg = Nbody6Config(
            InstallConfig(; enabled = false),
            BuildConfig(),
            SimulationConfig(; run_test = false, runs_dir = "runs"),
            PostprocessConfig(; enabled = false),
            VisualizationConfig(; enabled = false),
            MergerPipelineConfig(),
        )
        err = try
            run_pipeline(cfg; base_dir = base, run_id = "r1")
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("resume_pipeline", err.msg) && occursin(run_dir, err.msg)
        @test readdir(run_dir) == ["RUN_INFO.toml"]
    end
end
