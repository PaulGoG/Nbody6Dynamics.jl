# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Roadmap 6.5 (global IMF partition) and 6.6 (controls C-E and C-cold): one
# star list drawn for the whole complex and partitioned into the members, and
# the monolithic controls with the same star list and the same total energy
# (and, for C-cold, the same global virial ratio).

@testset "Complex IMF partition and controls (6.5, 6.6)" begin
    using StableRNGs
    nd = Nbody6Dynamics

    @testset "Explicit masses" begin
        imf = nd.ExplicitMassesIMF([0.5, 1.0, 2.0])
        @test nd.sample_masses(imf, 3, StableRNG(1)) == [0.5, 1.0, 2.0]
        @test_throws ArgumentError nd.sample_masses(imf, 4, StableRNG(1))
        @test nd.expected_mass(imf, 3) == 3.5
    end

    toml(extra) = """
[merger]
orbit_mode = "complex"
seed = 66

[merger.complex]
n_members = 4
mass_function = "equal"
mass_total = 800.0
size_mass = "constant"
rh = 0.5
centres = "plummer"
radius = 3.0
virial_ratio = 0.3
imf_partition = "global"
$extra

[merger.complex.member]
model = "plummer"
imf = "kroupa"
body1 = 2.0
bodyn = 0.1

[merger.output]
format = "nbody"
truncate_jacobi = false
"""
    load(t) = (p = joinpath(mktempdir(), "c.toml"); write(p, t); load_merger_config(p))
    meta(res) = Nbody6Dynamics.TOML.parsefile(joinpath(res.output_dir, "merger_ic.toml"))

    @testset "Global partition" begin
        cx = generate_merger_ic(load(toml("")); output_dir = mktempdir())
        @test length(cx.cluster_ranges) == 4
        masses = [sum(cx.mass_physical[r]) for r in cx.cluster_ranges]
        mmax = maximum(cx.mass_physical)
        @test all(200.0 .≤ masses[1:3] .< 200.0 + mmax)
        @test all(s.imf isa nd.ExplicitMassesIMF for s in cx.cluster_specs)
        # The same seed gives the same partition
        cx2 = generate_merger_ic(load(toml("")); output_dir = mktempdir())
        @test cx2.mass_physical == cx.mass_physical
        # A control needs the global partition
        @test_throws ArgumentError load(
            replace(
                toml("control = \"monolithic\""),
                "imf_partition = \"global\"" => "imf_partition = \"member\"",
            ),
        )
    end

    @testset "Monolithic and cold controls" begin
        cx = generate_merger_ic(load(toml("")); output_dir = mktempdir())
        ce = generate_merger_ic(load(toml("control = \"monolithic\"")); output_dir = mktempdir())
        cc = generate_merger_ic(load(toml("control = \"cold\"")); output_dir = mktempdir())
        for c in (ce, cc)
            @test length(c.cluster_ranges) == 1
            # The same star list
            @test sort(c.mass_physical) == sort(cx.mass_physical)
        end
        me, mc = meta(ce)["control"], meta(cc)["control"]
        @test me["kind"] == "monolithic" && mc["kind"] == "cold"
        @test me["e_realised"] ≈ me["e_target"] rtol = 1e-10
        @test me["q_realised"] ≈ 0.5 rtol = 1e-6
        @test mc["e_realised"] ≈ mc["e_target"] rtol = 1e-10
        @test mc["q_realised"] ≈ mc["q_target"] rtol = 1e-10
        @test me["e_target"] ≈ mc["e_target"]
        @test mc["q_target"] < 0.5   # the complex starts sub-virial (Q_cc = 0.3 plus internal energy)
    end
end
