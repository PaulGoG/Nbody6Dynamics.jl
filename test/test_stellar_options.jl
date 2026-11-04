# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.
#
# Stellar-evolution option set (roadmap 10.12): the [merger.stellar] option keys
# and their guards, the &INSSE/&INBSE/&INCOLL groups the generator writes, the
# IC record, and the comparison of a requested set with the engine's echo of it
# (READSE, input.F), whose layout the echo fixture below copies from an engine
# stdout capture of fork v2026.07.10 (Level C, no group assignments).

@testset "Stellar-evolution option set" begin
    nd = Nbody6Dynamics

    function stellar_config(stellar_lines::String; kz19::Int = 3, level::String = "C")
        text = """
        [merger]
        seed = 3
        n_clusters = 1
        orbit_mode = "single"

        [merger.cluster1]
        N = 300
        rbar = 1.0
        model = "plummer"

        [merger.stellar]
        kz19 = $kz19
        level = "$level"
        $stellar_lines
        """
        path = joinpath(mktempdir(), "stellar.toml")
        write(path, text)
        return path
    end

    # Level C echo as the engine prints it (rank 0, after the three group reads).
    level_c_echo = """

     SSE/BSE/Coll Level = C

           SSE FLAGS: EC,WD,NS,PS,MD,BH=   1   1   4   1   4   2 KMECH=   4
           MASSES: Chandra, NS0, NS1=   1.44D+00   1.80D+00   2.50D+00
           WIND PARAMS: NETA, BWIND, FLBV=   4.77D-01   0.00D+00   1.50D+00
           KICK PARAMS: DISP, ECSIG, WDSIG1,2,  WDKMAX, VFAC=   2.65D+02   3.00D+00   2.00D+00   2.00D+00   6.00D+00   0.00D+00

           BSE PARAMS: LAMBDA, ALPHA=   0.00D+00   1.00D+00 BHSPINFL, KICKTYPE=   0  -2
           K2,K3=   2.10D-01   2.10D-01 ACC1,2, BETA, XI=   3.92D+08   1.50D+00   1.25D-01   1.00D+00
           EPSNOV, EDDFAC, GAMM1=   1.00D-03   1.00D+02  -1.00D+00

           COLL PARAMS: FactorCl=   5.00D-01

    """
    # The same echo after the baseline set (psflag 2, kmech 1, eddfac 1).
    baseline_echo = replace(
        level_c_echo,
        "   1   1   4   1   4   2 KMECH=   4" => "   1   1   4   2   4   2 KMECH=   1",
        "   1.00D-03   1.00D+02  -1.00D+00" => "   1.00D-03   1.00D+00  -1.00D+00",
    )

    @testset "values written and resolved" begin
        @test all(isempty, values(nd.stellar_namelist_values(StellarSpec())))
        base = nd.stellar_namelist_values(StellarSpec(; option_set = "baseline"))
        @test base["INSSE"] == Dict(
            "ecflag" => 1,
            "nsflag" => 4,
            "psflag" => 2,
            "mdflag" => 4,
            "bhflag" => 2,
            "kmech" => 1,
            "disp" => 265.0,
            "ecsig" => 3.0,
        )
        @test base["INBSE"] == Dict("lambd1" => 0.0, "alphac" => 1.0, "eddfac" => 1.0)
        @test base["INCOLL"] == Dict("fctorcl" => 0.5)
        # explicit keys override the set; the [merger.nbody6] recoil and spin keys join INBSE
        arm = nd.stellar_namelist_values(
            StellarSpec(;
                option_set = "baseline",
                pair_instability = "strong",
                eddington_factor = 100.0,
                ns_kick_sigma = 0.0,
            ),
            Nbody6ParameterSpec(; gw_recoil = "fit_q", bh_spin = "mesa"),
        )
        @test arm["INSSE"]["psflag"] == 1 && arm["INSSE"]["disp"] == 0.0
        @test arm["INBSE"]["eddfac"] == 100.0
        @test arm["INBSE"]["kicktype"] == 0 && arm["INBSE"]["bhspinfl"] == 5
        p = nd.stellar_parameters(StellarSpec(; option_set = "baseline"))
        @test p["psflag"] == 2 && p["kmech"] == 1 && p["eddfac"] == 1.0
        @test p["nwind"] == 0.477 && p["mch"] == 1.44 && length(p) == 33
        @test nd.stellar_parameters(StellarSpec(; level = "A"))["disp"] == 190.0
        @test isempty(nd.stellar_parameters("0", Dict{String,Any}()))
        @test nd._namelist_group_line("INBSE", arm["INBSE"]) ==
              "&INBSE lambd1=0.0,alphac=1.0,kicktype=0,bhspinfl=5,eddfac=100.0 /"
        @test nd._namelist_group_line("INCOLL", Dict{String,Any}()) == "&INCOLL /"
    end

    @testset "generated input and IC record" begin
        cfg = nd.load_merger_config(stellar_config("option_set = \"baseline\""))
        dir = mktempdir()
        nd.generate_merger_ic(cfg; output_dir = dir)
        inp = read(joinpath(dir, "merger.inp"), String)
        @test occursin(
            "&INSSE ecflag=1,nsflag=4,psflag=2,mdflag=4,bhflag=2,kmech=1,disp=265.0,ecsig=3.0 /",
            inp,
        )
        @test occursin("&INBSE lambd1=0.0,alphac=1.0,eddfac=1.0 /", inp)
        @test occursin("&INCOLL fctorcl=0.5 /", inp)
        meta = nd.TOML.parsefile(joinpath(dir, "merger_ic.toml"))
        so = meta["stellar_options"]
        @test so["option_set"] == "baseline" && so["level"] == "C"
        @test so["written"]["INSSE"]["psflag"] == 2
        @test so["parameters"]["psflag"] == 2 && so["parameters"]["nwind"] == 0.477
        @test isnan(meta["stellar"]["eddington_factor"])   # the key itself was not set
        summary = read(joinpath(dir, "merger_summary.txt"), String)
        @test occursin("Stellar options: option_set=baseline; written: ecflag=1", summary)
        # the input requests exactly what the record holds
        req = nd.requested_stellar_parameters(joinpath(dir, "merger.inp"))
        @test req.level == "C" && isempty(req.unknown)
        @test req.parameters == so["parameters"]
        # no option key: empty groups, as earlier versions wrote them
        dir0 = mktempdir()
        nd.generate_merger_ic(nd.load_merger_config(stellar_config("")); output_dir = dir0)
        inp0 = read(joinpath(dir0, "merger.inp"), String)
        @test occursin("&INSSE /\n\n&INBSE /\n\n&INCOLL /\n", inp0)
        @test occursin(
            "Stellar options: option_set=level; written: none (Level C values otherwise)",
            read(joinpath(dir0, "merger_summary.txt"), String),
        )
    end

    @testset "configuration guards" begin
        load(lines; kw...) = nd.load_merger_config(stellar_config(lines; kw...))
        @test load("pair_instability = \"moderate\"").stellar.pair_instability == "moderate"
        @test load("ns_kick_sigma = 0.0").stellar.ns_kick_sigma == 0.0
        @test_throws ArgumentError load("option_set = \"science\"")
        @test_throws ArgumentError load("pair_instability = \"medium\"")
        @test_throws ArgumentError load("kick_mechanism = \"neutrino\"")
        @test_throws ArgumentError load("remnant_mass = \"delayed\"")
        @test_throws ArgumentError load("ns_kick_sigma = -100.0")
        @test_throws ArgumentError load("eddington_factor = 0.0")
        @test_throws ArgumentError load("bh_accreted_fraction = 1.5")
        @test_throws ArgumentError load("option_set = \"baseline\""; level = "B")
        @test_throws ArgumentError load("winds = \"belczynski2010\""; kz19 = 0)
        @test_throws ArgumentError load("bh_kicks = \"none\""; level = "0")
        @test_throws ArgumentError load("eddington_factor = \"level\"")
        @test_logs (:warn, r"eddington_factor = 100.0") match_mode = :any load(
            "eddington_factor = 100.0",
        )
        # without option keys the Level-0 and kz19 = 0 configurations load as before
        @test load(""; level = "0").stellar.option_set == "level"
        @test load(""; kz19 = 0).stellar.option_set == "level"
    end

    @testset "engine echo" begin
        dir = mktempdir()
        echo_file(text) = (p = joinpath(dir, "out$(rand(UInt32))"); write(p, text); p)
        inp_file(groups; level = "C") = (
            p = joinpath(dir, "in$(rand(UInt32)).inp");
            write(
                p,
                "&INNBODY6\nKSTART=1 /\n\n&ININPUT\nN=100,\nLevel='$level' /\n\n" *
                groups *
                "\n&INDATA\nZMET=0.02 /\n",
            );
            p
        )
        empty_groups = "&INSSE /\n\n&INBSE /\n\n&INCOLL /\n"
        base_groups =
            "&INSSE ecflag=1,nsflag=4,psflag=2,mdflag=4,bhflag=2,kmech=1,disp=265.0,ecsig=3.0 /\n\n" *
            "&INBSE lambd1=0.0,alphac=1.0,eddfac=1.0 /\n\n&INCOLL fctorcl=0.5 /\n"

        echo = nd.read_stellar_echo(echo_file(" N = 100\n" * level_c_echo * " ALPHAS = 2.35\n"))
        @test echo.complete && isempty(echo.invalid) && !echo.level_zero
        @test length(echo.values) == 33
        @test echo.values["psflag"] == 1 && echo.values["kicktype"] == -2
        @test echo.values["acc1"] == 3.92e8 && echo.values["gamm1"] == -1.0

        out_c = echo_file(level_c_echo)
        out_b = echo_file(baseline_echo)
        @test nd.check_stellar_echo(inp_file(empty_groups), out_c).status == "passed"
        @test nd.check_stellar_echo(inp_file(base_groups), out_b).status == "passed"
        # the baseline input against an engine that ran Level C: every departure named
        r = nd.check_stellar_echo(inp_file(base_groups), out_c)
        @test r.status == "mismatch"
        @test any(startswith("psflag:"), r.mismatches) &&
              any(startswith("kmech:"), r.mismatches) &&
              any(startswith("eddfac:"), r.mismatches) &&
              length(r.mismatches) == 3
        # a group the engine could not read
        out_inv = echo_file("  insse: invalid Namelist format\n" * level_c_echo)
        r = nd.check_stellar_echo(inp_file(base_groups), out_inv)
        @test r.status == "invalid_namelist" && r.mismatches == ["insse: invalid Namelist format"]
        # a name no group defines is a mismatch even when the echo agrees
        r = nd.check_stellar_echo(inp_file("&INSSE psflg=2 /\n\n&INBSE /\n\n&INCOLL /\n"), out_c)
        @test r.status == "mismatch" &&
              r.mismatches == ["INSSE.psflg: not a numerical member of the group"]
        # Fortran D exponents and upper-case names in a hand-written group
        @test nd.check_stellar_echo(
            inp_file("&INSSE PSFLAG=2,KMECH=1 /\n\n&INBSE EDDFAC=1.0D0 /\n\n&INCOLL /\n"),
            out_b,
        ).status == "passed"
        # pending until COLL PARAMS; not applicable for Level '0'
        partial = echo_file(first(level_c_echo, findfirst("BSE PARAMS", level_c_echo)[1] - 1))
        @test nd.check_stellar_echo(inp_file(empty_groups), partial).status == "pending"
        @test nd.check_stellar_echo(inp_file(empty_groups; level = "0"), out_c).status ==
              "not_applicable"
        @test nd.read_stellar_echo(echo_file("  Run without stellar evolution\n")).level_zero
        # an overflowed field makes the echo unreadable
        r = nd.check_stellar_echo(
            inp_file(empty_groups),
            echo_file(replace(level_c_echo, "KMECH=   4" => "KMECH=****")),
        )
        @test r.status == "mismatch" &&
              r.mismatches == ["echo unreadable (overflowed or missing fields)"]
        # a restart's echo is read from the segment's offset on
        two = echo_file(level_c_echo * " ADJUST: TIME 1.0\n" * baseline_echo)
        offset = length(codeunits(level_c_echo * " ADJUST: TIME 1.0\n"))
        @test nd.check_stellar_echo(inp_file(base_groups), two; offset = offset).status == "passed"
        @test nd.check_stellar_echo(inp_file(base_groups), two).status == "mismatch"
    end

    @testset "start-up check stops a mismatching engine" begin
        dir = mktempdir()
        inp = joinpath(dir, "merger.inp")
        write(inp, "&ININPUT\nLevel='C' /\n\n&INSSE psflag=2 /\n\n&INBSE /\n\n&INCOLL /\n")
        out = joinpath(dir, "out1000")
        write(out, "")
        # A stand-in engine: prints the Level C echo, then would run for a minute.
        script = joinpath(dir, "engine.sh")
        write(script, "cat $(joinpath(dir, "echo.txt")) >> $out\nsleep 60\n")
        write(joinpath(dir, "echo.txt"), level_c_echo)
        proc = run(detach(`bash $script`); wait = false)
        t0 = time()
        seen = Ref{Any}(nothing)
        chk = nd._start_echo_check(inp, out, proc; poll = 0.1, on_result = r -> (seen[] = r))
        wait(chk.task)
        wait(proc)
        @test chk.fired[] && chk.result[].status == "mismatch"
        @test seen[] !== nothing && seen[].status == "mismatch"   # recorded before the kill
        @test time() - t0 < 30
        # an engine that exits without an echo is recorded as absent
        write(out, "")
        proc2 = run(`true`; wait = false)
        chk2 = nd._start_echo_check(inp, out, proc2; poll = 0.1)
        wait(chk2.task)
        @test !chk2.fired[] && chk2.result[].status == "absent"
    end
end
