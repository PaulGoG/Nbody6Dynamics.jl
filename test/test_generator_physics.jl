# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Physical conventions of the generator: two-cluster truncation at the Roche
# lobe (Eggleton 1983) instead of the Hill radius, which exceeds the inner
# Lagrange point at comparable masses; the inertial-to-rotating conversion as
# the default under a rotating-frame tide; a centre-of-mass frame after
# truncation; complex partitions that hit every member's target mass; energy
# matching of the controls on the systems; strict cluster tables.

@testset "Generator physics" begin
    nd = Nbody6Dynamics
    load(t) = (p = joinpath(mktempdir(), "g.toml"); write(p, t); load_merger_config(p))
    meta(res) = Nbody6Dynamics.TOML.parsefile(joinpath(res.output_dir, "merger_ic.toml"))
    eggleton(q) = 0.49 * q^(2 / 3) / (0.6 * q^(2 / 3) + log(1 + cbrt(q)))

    pair_toml(; output = "truncate_jacobi = true", tidal = "", cluster1 = "") = """
[merger]
n_clusters = 2
orbit_mode = "kepler"
seed = 811

[merger.cluster1]
model = "plummer"
N = 1500
rbar = 2.0
imf = "kroupa"
body1 = 2.0
bodyn = 0.1
$cluster1

[merger.cluster2]
model = "plummer"
N = 1500
rbar = 2.0
imf = "kroupa"
body1 = 2.0
bodyn = 0.1

[merger.orbit]
apocentre = 8.0
eccentricity = 0.5

[merger.output]
format = "nbody"
$output
$tidal

[merger.nbody6]
qe = 0.01
"""

    @testset "Truncation radius of a cluster pair" begin
        roche = generate_merger_ic(load(pair_toml()); output_dir = mktempdir())
        hill = generate_merger_ic(
            load(pair_toml(; output = "truncate_jacobi = true\ntruncate_rule = \"hill\"")),
            output_dir = mktempdir(),
        )
        mr, mh = meta(roche), meta(hill)
        @test mr["output"]["truncate_rule"] == "roche" && mh["output"]["truncate_rule"] == "hill"
        # Sampled masses: kept plus truncated [M☉]
        sampled(m, i) = m["cluster$i"]["mass"] + m["cluster$i"]["mass_truncated"]
        M1, M2 = sampled(mr, 1), sampled(mr, 2)
        for (i, q) in ((1, M1 / M2), (2, M2 / M1))
            # Lobe of the sampled masses at the apocentre separation
            @test mr["cluster$i"]["truncation_radius_pc"] ≈ 8.0 * eggleton(q) rtol = 0.03
            @test mh["cluster$i"]["truncation_radius_pc"] ≈ 8.0 * cbrt(q / 3) rtol = 0.03
            # The lobe lies inside the inner Lagrange point, the Hill radius does not
            @test mr["cluster$i"]["truncation_radius_pc"] <
                  4.0 <
                  mh["cluster$i"]["truncation_radius_pc"]
            @test mr["cluster$i"]["mass_truncated"] > mh["cluster$i"]["mass_truncated"] ≥ 0
        end
        untruncated = generate_merger_ic(
            load(pair_toml(; output = "truncate_jacobi = false"));
            output_dir = mktempdir(),
        )
        @test meta(untruncated)["cluster1"]["truncation_radius_pc"] == 0.0
        @test_throws ArgumentError load(
            pair_toml(; output = "truncate_jacobi = true\ntruncate_rule = \"tidal\""),
        )
        # Centre-of-mass frame after truncation
        for res in (roche, hill)
            m = res.mass_physical
            @test maximum(abs.(res.pos_physical * m)) / sum(m) <
                  1e-10 * maximum(abs.(res.pos_physical))
            @test maximum(abs.(res.vel_physical * m)) / sum(m) <
                  1e-10 * maximum(abs.(res.vel_physical))
        end
    end

    @testset "Velocity frame under a rotating-frame tide" begin
        tide(frame) =
            "\n[merger.tidal]\nkz14 = 1" *
            (isempty(frame) ? "" : "\nvelocity_frame = \"$frame\"") *
            "\n"
        gen(frame) = generate_merger_ic(
            load(pair_toml(; output = "truncate_jacobi = false", tidal = tide(frame))),
            output_dir = mktempdir(),
        )
        default, inertial, rotating = gen(""), gen("inertial"), gen("rotating")
        # The generator builds clusters and the orbit in an inertial frame: by default the velocities are
        # converted to the engine's rotating frame
        @test meta(default)["tidal"]["velocity_frame"] == "inertial"
        @test default.vel_physical == inertial.vel_physical
        @test default.vel_physical != rotating.vel_physical
    end

    complex_toml(extra; binaries = "", member_extra = "") = """
[merger]
orbit_mode = "complex"
seed = 912

[merger.complex]
n_members = 12
mass_function = "power_law"
beta = 2.0
mass_lo = 100.0
mass_hi = 400.0
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
body1 = 5.0
bodyn = 0.1
$member_extra
$binaries

[merger.output]
format = "nbody"
truncate_jacobi = false
"""

    @testset "Global partition hits every member's mass" begin
        cx = generate_merger_ic(load(complex_toml("")); output_dir = mktempdir())
        targets = meta(cx)["complex"]["realised"]["member_target_masses"]
        realised = [sum(cx.mass_physical[r]) for r in cx.cluster_ranges]
        mmax = maximum(cx.mass_physical)
        # Within one star of its target, the last member included
        @test all(abs.(realised .- targets) .≤ mmax)
        @test abs(sum(realised) - sum(targets)) ≤ mmax
    end

    @testset "Controls with primordial binaries" begin
        bins = "\n[merger.complex.member.binaries]\nfraction = 0.3\n"
        gen(extra, b) =
            generate_merger_ic(load(complex_toml(extra; binaries = b)); output_dir = mktempdir())
        rh(res) = (
            m = res.mass_physical;
            r = vec(sqrt.(sum(abs2, res.pos_physical .- res.pos_physical * m ./ sum(m); dims = 1)));
            o = sortperm(r);
            r[o[findfirst(≥(0.5 * sum(m)), cumsum(m[o]))]]
        )
        ce0 = gen("control = \"monolithic\"", "")
        ce1 = gen("control = \"monolithic\"", bins)
        m1 = meta(ce1)["control"]
        @test m1["e_realised"] ≈ m1["e_target"] rtol = 1e-10
        # The energy is matched on the systems (pairs as point masses), as the control is scaled: the
        # binaries' binding energy no longer compresses the control
        @test 0.7 < rh(ce1) / rh(ce0) < 1.4
        # A control needs the same star list: pairings that change masses are refused
        @test_throws ArgumentError load(
            complex_toml(
                "control = \"monolithic\"";
                binaries = "\n[merger.complex.member.binaries]\nfraction = 0.3\npairing = \"uniform_q\"\n",
            ),
        )
    end

    @testset "Strict cluster tables" begin
        member_error = try
            load(complex_toml(""; member_extra = "colour = 1"))
            ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("merger.complex.member", member_error) && !occursin("cluster0", member_error)
        # Values of the wrong kind are refused, not ignored
        @test_throws ArgumentError load(pair_toml(; cluster1 = "profile = \"plummer\""))
        @test_throws ArgumentError load(
            replace(pair_toml(), "[merger.cluster1]" => "[merger.cluster1]\nbinaries = 0.2"),
        )
        @test_throws ArgumentError load(pair_toml(; cluster1 = "position = [1.0, 0.0, 0.0]"))
        @test_throws ArgumentError load(
            pair_toml(; cluster1 = "\n[merger.cluster1.binaries]\nfraction = 1.5"),
        )
    end
end
