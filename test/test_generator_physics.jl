# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Physical conventions of the generator: two-cluster truncation at the Roche
# lobe (Eggleton 1983) instead of the Hill radius, which exceeds the inner
# Lagrange point at comparable masses; the inertial-to-rotating conversion as
# the default under a rotating-frame tide; a centre-of-mass frame after
# truncation; complex partitions that hit every member's target mass; energy
# matching of the controls on the systems; strict cluster tables.

using Logging: Logging

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

@testset "Generator review fixes of 2026-11" begin
    nd = Nbody6Dynamics
    cluster(; kw...) = Dict{String,Any}(string(k) => v for (k, v) in kw)
    merger(c1, c2) = Dict{String,Any}(
        "merger" => Dict{String,Any}(
            "n_clusters" => 2,
            "cluster1" => c1,
            "cluster2" => c2,
            "orbit" => Dict{String,Any}("apocentre" => 10.0, "eccentricity" => 0.5),
        ),
    )
    control_rbar(d) = control_merger_dict(d)["merger"]["cluster1"]["rbar"]

    # Equal density at eight times the mass: r₂ = 2 r₁, and the N-weighted mean of (1, 2) is 1.5
    by_mass = merger(
        cluster(; N = 100, rbar = 1.0, imf = "equal", mass_total = 100.0),
        cluster(; N = 100, density_ratio = 1.0, imf = "equal", mass_total = 800.0),
    )
    @test control_rbar(by_mass) ≈ 1.5
    # Eight times the density at eight times the mass: the same radius
    by_mass["merger"]["cluster2"]["density_ratio"] = 8.0
    @test control_rbar(by_mass) ≈ 1.0
    # The same mass function without totals: the star counts give the mass ratio,
    # r₂ = 2 r₁ and the mean is (100 + 800 × 2)/900
    by_count = merger(cluster(; N = 100, rbar = 1.0), cluster(; N = 800, density_ratio = 1.0))
    @test control_rbar(by_count) ≈ 1700 / 900
    # Equal masses of different particle mass: eight times the mass, r₂ = 2 r₁
    equal_imf(m) = Dict{String,Any}("type" => "equal", "particle_mass" => m)
    by_particle = merger(
        cluster(; N = 100, rbar = 1.0, imf = equal_imf(1.0)),
        cluster(; N = 100, density_ratio = 1.0, imf = equal_imf(8.0)),
    )
    @test control_rbar(by_particle) ≈ 1.5
    # … and the control's stars have the mean mass, (100 + 800) M☉ / 200
    @test control_merger_dict(by_particle)["merger"]["cluster1"]["imf"]["particle_mass"] ≈ 4.5
    # Equal masses against a Kroupa member: the radius the generator gives cluster 2,
    # M₂/M₁ = 100 ⟨m⟩_Kroupa / 800 M☉, and the control has the merger's total mass
    mixed = merger(
        cluster(; N = 100, rbar = 1.0, imf = "equal", mass_total = 800.0),
        cluster(; N = 100, density_ratio = 1.0),
    )
    mixed_path = joinpath(mktempdir(), "mixed.toml")
    open(io -> nd.TOML.print(io, mixed), mixed_path, "w")
    r2 = load_merger_config(mixed_path).clusters[2].rbar
    m_kroupa = 100 * nd.kroupa_mean_mass(0.08, 100.0)
    @test r2 ≈ cbrt(m_kroupa / 800)
    @test control_rbar(mixed) ≈ (100 * 1.0 + 100 * r2) / 200
    @test control_merger_dict(mixed)["merger"]["cluster1"]["mass_total"] ≈ 800 + m_kroupa
    # A Kroupa cluster 1 without a total cannot carry the equal-mass member's mass: the
    # control of its mass function would have 200 ⟨m⟩_Kroupa instead of the merger's mass
    swapped = merger(
        cluster(; N = 100, rbar = 1.0),
        cluster(; N = 100, density_ratio = 1.0, imf = "equal", mass_total = 800.0),
    )
    @test_throws r"expected mass of .* against" control_merger_dict(swapped)
    # density_ratio as the loader takes it: positive, finite, not beside rbar, not in cluster 1
    for bad in (0.0, -1.0, NaN)
        @test_throws r"density_ratio must be finite and > 0" control_merger_dict(
            merger(cluster(; N = 100, rbar = 1.0), cluster(; N = 100, density_ratio = bad)),
        )
    end
    @test_throws r"mutually exclusive" control_merger_dict(
        merger(cluster(; N = 100, rbar = 1.0), cluster(; N = 100, rbar = 1.0, density_ratio = 2.0)),
    )
    @test_throws r"cluster2\] and later only" control_merger_dict(
        merger(cluster(; N = 100, density_ratio = 2.0), cluster(; N = 100, rbar = 1.0)),
    )

    mktempdir() do dir
        write_config(name, body) = (p = joinpath(dir, name); write(p, body); p)
        # Flat table: mass bounds do not apply to equal masses
        flat = """
            [merger]
            n_clusters = 1
            orbit_mode = "single"
            [merger.cluster1]
            model = "plummer"
            N = 100
            rbar = 1.0
            imf = "equal"
            mass_total = 100.0
            EXTRA
            [merger.output]
            output_dir = "$dir"
            """
        @test load_merger_config(write_config("eq.toml", replace(flat, "EXTRA" => ""))).clusters[1].imf isa
              nd.EqualMassIMF
        @test_throws ArgumentError load_merger_config(
            write_config("eq_bodyn.toml", replace(flat, "EXTRA" => "bodyn = 0.5")),
        )
        @test_throws ArgumentError load_merger_config(
            write_config("eq_body1.toml", replace(flat, "EXTRA" => "body1 = 50.0")),
        )

        # An IMBH given by mass fraction under stellar evolution: 4 % of about 1150 M☉ is
        # below 100 M☉ and is refused at generation; 4 % of about 3400 M☉ is accepted
        imbh = """
            [merger]
            seed = 3
            n_clusters = 1
            orbit_mode = "single"
            [merger.cluster1]
            model = "plummer"
            N = NSTARS
            rbar = 1.0
            imf = "kroupa"
            [merger.cluster1.imbh]
            mass_fraction = 0.04
            [merger.stellar]
            kz19 = 3
            [merger.nbody6]
            binary_physics = "full"
            escapers = "keep"
            [merger.output]
            output_dir = "OUT"
            """
        light = joinpath(dir, "light")
        heavy = joinpath(dir, "heavy")
        mkpath(light)
        mkpath(heavy)
        cfg_light = load_merger_config(
            write_config("light.toml", replace(imbh, "NSTARS" => "2000", "OUT" => light)),
        )
        err = try
            generate_merger_ic(cfg_light)
            nothing
        catch e
            e
        end
        @test err isa ErrorException && occursin("below 100 M☉", err.msg)
        cfg_heavy = load_merger_config(
            write_config("heavy.toml", replace(imbh, "NSTARS" => "6000", "OUT" => heavy)),
        )
        # One warning on the mass fraction, the realised one: the specification's was
        # given when the configuration was loaded
        logs, result = Test.collect_test_logs(; min_level = Logging.Warn) do
            generate_merger_ic(cfg_heavy)
        end
        @test count(l -> occursin("exceeds 0.02", string(l.message)), logs) == 1
        @test result.N_total == 6001
        # The record keeps the IMBH of its cluster
        # … read back without a second validation, so without the mass-fraction warning
        reloaded =
            (@test_logs min_level = Logging.Warn load_merger_ic_result(result.output_dir)).cluster_specs[1].imbh
        @test reloaded isa IMBHSpec
        @test reloaded.mass_fraction == 0.04 && isnan(reloaded.mass)
        @test reloaded.placement == "centre" && reloaded.response == "rigid"
        # A specification given in memory meets the bounds of the loader at generation
        bad_imbh = MergerConfig(
            [
                ClusterSpec(
                    N = 200,
                    rbar = 1.0,
                    profile = PlummerProfile(),
                    imf = EqualMassIMF(1.0),
                    imbh = IMBHSpec(mass = 40.0, spin = 1.5),
                ),
            ],
            "single",
            OrbitSpec(),
            MergerOutputSpec(output_dir = joinpath(dir, "bad_imbh")),
            Nbody6ParameterSpec(),
            StellarSpec(kz19 = 0, level = "0"),
            TidalSpec(),
            1,
            200000,
            "shared",
            ComplexSpec(),
            nothing,
        )
        @test_throws r"imbh.spin must satisfy" generate_merger_ic(bad_imbh)
    end

    @testset "Half-mass radius of a sampled cluster" begin
        using StableRNGs: StableRNG
        # The target is met in the frame the cluster is handed on in: centre of mass at the
        # origin, mass-weighted half-mass radius rbar, with a mass spectrum and with pairs
        for profile in (PlummerProfile(), KingProfile(; W0 = 6.0)),
            binaries in (BinarySpec(), BinarySpec(; fraction = 0.3))

            spec = ClusterSpec(;
                N = 600,
                rbar = 1.7,
                profile = profile,
                imf = KroupaIMF(),
                binaries = binaries,
            )
            c = nd._sample_cluster(spec; rng = StableRNG(17), nmax = 600)
            centre = [sum(c.mass .* c.pos[k, :]) for k in 1:3] ./ sum(c.mass)
            @test maximum(abs, centre) < 1e-12
            @test nd.half_mass_radius(c.mass, c.pos; centre = zeros(3)) ≈ 1.7 rtol = 1e-12
            @test c.T / abs(c.W) ≈ 0.5 rtol = 1e-12
        end
    end
end
