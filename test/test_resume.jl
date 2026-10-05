# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Restart dump selection and join" begin
    record(payload::Vector{UInt8}) = vcat(
        reinterpret(UInt8, [Int32(length(payload))]),
        payload,
        reinterpret(UInt8, [Int32(length(payload))]),
    )
    dump_bytes() = vcat(record(collect(UInt8, 1:40)), record(zeros(UInt8, 16)), record(UInt8[7]))
    marker_lines(t, name, unit) =
        "  MYDUMP    $(t)      $(name)                         $(unit) $(name)                \n" *
        "  W MYDUMP J,II,NPARTMP=         $(unit)           1     1572864\n" *
        "  NA-NS=          85         168         530\n"
    adj(t) =
        " ADJUST:  TIME    $(t)  T[Myr]   0.100E+01  Q   0.500E+00  DE   0.000E+00" *
        " DELTA   0.000E+00 DETOT   0.000E+00 E  -2.500E-01\n"
    first_line(s) = s[1:findfirst('\n', s)]
    NA_NS = "  NA-NS=          85         168         530\n"

    # The marker block parses to the dump's name and exact time
    d0 = mktempdir()
    write(joinpath(d0, "out1000"), marker_lines(1.5, "comm.2_1.5", 204))
    mk0 = Nbody6Dynamics._dump_markers(joinpath(d0, "out1000"))
    @test length(mk0) == 1 && mk0[1].file == "comm.2_1.5" && mk0[1].time_nb == 1.5

    # 1. Record walk of a dump file
    complete = Nbody6Dynamics._dump_complete
    d1 = mktempdir()
    probe(bytes) = (p = joinpath(d1, "probe"); write(p, bytes); complete(p))
    @test probe(dump_bytes())
    @test !probe(record(collect(UInt8, 1:40)))
    @test !probe(dump_bytes()[1:(end - 3)])
    @test !probe(vcat(dump_bytes(), UInt8[1, 2]))
    mismatched = vcat(
        reinterpret(UInt8, [Int32(4)]),
        zeros(UInt8, 4),
        reinterpret(UInt8, [Int32(5)]),
        record(zeros(UInt8, 16)),
    )
    @test !probe(mismatched)
    @test !probe(UInt8[])
    @test !complete(joinpath(d1, "absent"))

    # 2. Choice of the restart dump
    resume = Nbody6Dynamics._resume_dump
    d = mktempdir()
    so = joinpath(d, "out1000")
    head_a = adj("1.00000E+00") * marker_lines(1.0, "comm.2_1.0", 203) * adj("1.50000E+00")
    text_a =
        head_a *
        marker_lines(1.5, "comm.2_1.5", 204) *
        adj("2.00000E+00") *
        marker_lines(2.0, "comm.2_2.0", 205)
    write(so, text_a)
    write(joinpath(d, "comm.2_1.0"), dump_bytes())
    write(joinpath(d, "comm.2_1.5"), dump_bytes())
    write(joinpath(d, "comm.2_2.0"), dump_bytes()[1:20])
    mk_a = resume(d, so)
    @test mk_a !== nothing
    @test mk_a.file == "comm.2_1.5" && mk_a.time_nb == 1.5
    @test mk_a.line_end == sizeof(head_a * first_line(marker_lines(1.5, "comm.2_1.5", 204)))
    write(joinpath(d, "comm.2_2.0"), dump_bytes())
    @test resume(d, so).file == "comm.2_2.0"

    # 3. End of the lines printed with a dump
    block_end = Nbody6Dynamics._marker_block_end
    pos = block_end(so, mk_a)
    @test endswith(text_a[1:pos], NA_NS)
    @test startswith(text_a[(pos + 1):end], " ADJUST:  TIME    2.00000E+00")
    @test block_end(so, (file = "comm.2_2.0", time_nb = 2.0, line_end = filesize(so))) ==
          filesize(so)
    @test_throws ArgumentError block_end(so, (file = "comm.2_1.5", time_nb = 1.5, line_end = -1))

    # 2c. A name written twice is judged by its last marker only
    d = mktempdir()
    so = joinpath(d, "out1000")
    write(
        so,
        marker_lines(1.0, "comm.2_1.0", 203) *
        marker_lines(13.390625, "comm.1_13.4", 228) *
        marker_lines(13.390625, "comm.1_13.4", 201),
    )
    write(joinpath(d, "comm.2_1.0"), dump_bytes())
    write(joinpath(d, "comm.1_13.4"), dump_bytes()[1:20])
    @test resume(d, so).file == "comm.2_1.0"

    # 2d. No marker: the latest complete dump by name
    d = mktempdir()
    so = joinpath(d, "out1000")
    write(so, adj("1.00000E+00"))
    write(joinpath(d, "comm.2_1.0"), dump_bytes())
    write(joinpath(d, "comm.2_10.0"), dump_bytes())
    @test resume(d, so) == (file = "comm.2_10.0", time_nb = 10.0, line_end = -1)

    # 2e. No marker and no complete dump
    d = mktempdir()
    so = joinpath(d, "out1000")
    write(so, adj("1.00000E+00"))
    @test resume(d, so) === nothing
    write(joinpath(d, "comm.2_1.0"), dump_bytes()[1:20])
    @test resume(d, so) === nothing

    # 4. Cut offset of a series file
    cut_offset = Nbody6Dynamics._series_cut_offset
    d4 = mktempdir()
    series_path = joinpath(d4, "lagr.7")
    series_lines = [
        "## header\n",
        "TIME   a   b\n",
        "0.00000000000000000E+00 1 2\n",
        "0.10000000000000000E+01 1 2\n",
        "0.20000000000000000E+01 1 2\n",
        "0.30000000000000000E+01 1 2\n",
    ]
    write(series_path, join(series_lines))
    off_2 = sum(sizeof, series_lines[1:4])
    off_3 = sum(sizeof, series_lines[1:5])
    @test cut_offset(series_path, 1.0) == off_2
    @test cut_offset(series_path, 3.0) === nothing
    @test cut_offset(series_path, 1.99) == off_2
    @test cut_offset(series_path, 2.0 - 1.0e-12) == off_3
    @test cut_offset(joinpath(d4, "absent"), 1.0) === nothing

    # 5. Moving a tail
    move_tail = Nbody6Dynamics._move_tail
    d5 = mktempdir()
    src = joinpath(d5, "a")
    dest = joinpath(d5, "tails", "a.tail")
    write(src, "abcdef")
    @test move_tail(src, 2, dest) == 4
    @test read(src, String) == "ab" && read(dest, String) == "cdef"
    src2 = joinpath(d5, "b")
    write(src2, "xyz")
    @test move_tail(src2, 1, dest) == 2
    @test read(src2, String) == "x" && read(dest, String) == "yz"
    @test read(joinpath(d5, "tails", "a#1.tail"), String) == "cdef"
    dest3 = joinpath(d5, "tails", "c.tail")
    @test move_tail(src2, filesize(src2), dest3) == 0
    @test !isfile(dest3) && read(src2, String) == "x"

    # 6. Time unit from global.30
    d6 = mktempdir()
    write(
        joinpath(d6, "global.30"),
        "TIME[NB} TIME[Myr] TCR[Myr]\n" * "0.0 0.0 1\n" * "0.5 3.25 1\n",
    )
    @test Nbody6Dynamics._time_unit_myr(d6) == 6.5
    @test isnan(Nbody6Dynamics._time_unit_myr(mktempdir()))

    # 7. Time in the name of a time-stamped file
    stamp = Nbody6Dynamics._timestamped_time
    @test stamp("conf.3_12.5") == 12.5
    @test stamp("comm.2_10.0") == 10.0
    @test stamp("sev.83_3") == 3.0
    for name in ("dat.10", "lagr.7", "merger_ic.toml", "telemetry_2.csv", "comm.2_notatime")
        @test isnan(stamp(name))
    end

    # 8. TCRIT of an engine input
    d8 = mktempdir()
    inp = joinpath(d8, "a.inp")
    write(
        inp,
        "! TCRIT = 99 in words\n" *
        "&ININPUT\n" *
        "N=1000,DTADJ=0.5,DELTAT=1,TCRIT=40.00,QE=2.000E-04 /\n",
    )
    @test Nbody6Dynamics._tcrit_of(inp) == 40.0
    write(inp, "&ININPUT\nN=1000,TCRIT=1.5D+02,QE=2.000E-04 /\n")
    @test Nbody6Dynamics._tcrit_of(inp) == 150.0
    write(inp, "&ININPUT\nN=1000,QE=2.000E-04 /\n")
    @test_throws ErrorException Nbody6Dynamics._tcrit_of(inp)

    # 9. Join after a killed segment
    discard = Nbody6Dynamics._discard_after
    cut_names = ("out1000", "lagr.7", "global.30", "esc.11", "event.35")
    function build_killed(out)
        write(
            joinpath(out, "out1000"),
            adj("0.00000E+00") *
            adj("1.00000E+00") *
            adj("2.00000E+00") *
            marker_lines(2.0, "comm.2_2.0", 205) *
            " WD/NS/BH FORMATION T=  2.10000E+00\n" *
            adj("2.50000E+00") *
            marker_lines(2.5, "comm.2_2.5", 206) *
            adj("3.00000E+00"),
        )
        lagr_times = (
            "0.00000000000000000E+00",
            "0.10000000000000000E+01",
            "0.20000000000000000E+01",
            "0.25000000000000000E+01",
            "0.30000000000000000E+01",
        )
        write(joinpath(out, "lagr.7"), "## header\nTIME a\n" * join("$(t) 1\n" for t in lagr_times))
        times = (0.0, 1.0, 2.0, 3.0)
        write(
            joinpath(out, "global.30"),
            "TIME[NB} TIME[Myr]\n" * join("$(t) $(6.5 * t)\n" for t in times),
        )
        write(
            joinpath(out, "esc.11"),
            "         TTOT         BODY\n" * "  1.00000E+00 1\n" * "  2.50000E+00 2\n",
        )
        write(
            joinpath(out, "event.35"),
            "TIME[Myr] NDISS\n" * join("$(6.5 * t) 0\n" for t in times),
        )
        for f in ("conf.3_2", "conf.3_3", "sev.83_3", "comm.2_2.5", "comm.2_3.0")
            write(joinpath(out, f), "content of $f\n")
        end
        for f in ("coll.13", "err1000", "dat.10")
            write(joinpath(out, f), "content of $f\n")
        end
        write(joinpath(out, "comm.2_2.0"), dump_bytes())
        touch(joinpath(out, "status.36"))
        return Dict(name => read(joinpath(out, name), String) for name in cut_names)
    end

    out = mktempdir()
    orig = build_killed(out)
    untouched =
        Dict(f => read(joinpath(out, f)) for f in ("coll.13", "err1000", "dat.10", "status.36"))
    so9 = joinpath(out, "out1000")
    m = first(Nbody6Dynamics._dump_markers(so9))
    @test m.file == "comm.2_2.0"
    r = discard(out, so9, m; segment = 1, last_status = "killed")
    @test r["moved"] == sort([
        "out1000.tail",
        "lagr.7.tail",
        "global.30.tail",
        "esc.11.tail",
        "event.35.tail",
        "conf.3_3",
        "sev.83_3",
        "comm.2_2.5",
        "comm.2_3.0",
    ])
    @test r["dir"] == joinpath("discarded", "segment_1")
    @test r["t_from"] == 2.0 && r["t_to"] == 3.0
    @test r["uncut"] == ["coll.13", "err1000"]
    disc = joinpath(out, "discarded", "segment_1")
    so_text = read(so9, String)
    @test endswith(so_text, NA_NS)
    @test count("ADJUST:", so_text) == 3
    for name in cut_names
        @test read(joinpath(out, name), String) * read(joinpath(disc, name * ".tail"), String) ==
              orig[name]
    end
    lagr_kept = split(read(joinpath(out, "lagr.7"), String), '\n'; keepempty = false)
    @test startswith(last(lagr_kept), "0.20000000000000000E+01")
    @test read(joinpath(out, "esc.11"), String) == "         TTOT         BODY\n  1.00000E+00 1\n"
    @test read(joinpath(out, "event.35"), String) == "TIME[Myr] NDISS\n0.0 0\n6.5 0\n13.0 0\n"
    @test read(joinpath(out, "global.30"), String) ==
          "TIME[NB} TIME[Myr]\n0.0 0.0\n1.0 6.5\n2.0 13.0\n"
    @test isfile(joinpath(out, "conf.3_2")) && isfile(joinpath(out, "comm.2_2.0"))
    for f in ("conf.3_3", "sev.83_3", "comm.2_2.5", "comm.2_3.0")
        @test isfile(joinpath(disc, f)) && !isfile(joinpath(out, f))
    end
    for (f, bytes) in untouched
        @test read(joinpath(out, f)) == bytes
    end
    r2 = discard(out, so9, m; segment = 1, last_status = "killed")
    @test isempty(r2["moved"]) && r2["dir"] == ""

    # 10. Join after a stop
    out = mktempdir()
    so10 = joinpath(out, "out1000")
    after_block = "\n\n         COMMON SAVED AT TOFF/TIME/TTOT =  0.0\n"
    write(
        so10,
        adj("1.00000E+00") *
        "         TERMINATION BY MANUAL INTERVENTION\n" *
        marker_lines(1.390625, "comm.1_1.4", 228) *
        after_block,
    )
    write(joinpath(out, "lagr.7"), "0.0 1\n1.0 1\n")
    write(joinpath(out, "conf.3_1"), "content of conf.3_1\n")
    write(joinpath(out, "comm.1_1.4"), dump_bytes())
    snapshot(dir) = Dict(f => read(joinpath(dir, f)) for f in readdir(dir))
    before = snapshot(out)
    m10 = only(Nbody6Dynamics._dump_markers(so10))
    r = discard(out, so10, m10; segment = 2, last_status = "stopped")
    @test isempty(r["moved"]) && r["dir"] == ""
    @test snapshot(out) == before
    r = discard(out, so10, m10; segment = 2, last_status = "killed")
    @test r["moved"] == ["out1000.tail"]
    @test read(joinpath(out, "discarded", "segment_2", "out1000.tail"), String) == after_block
    @test endswith(read(so10, String), NA_NS)

    # 11. A dump without a stdout line: the stdout capture is left whole
    out = mktempdir()
    orig = build_killed(out)
    so11 = joinpath(out, "out1000")
    m11 = (file = "comm.2_2.0", time_nb = 2.0, line_end = -1)
    r = discard(out, so11, m11; segment = 1, last_status = "killed")
    @test !("out1000.tail" in r["moved"])
    @test r["moved"] == sort([
        "lagr.7.tail",
        "global.30.tail",
        "esc.11.tail",
        "event.35.tail",
        "conf.3_3",
        "sev.83_3",
        "comm.2_2.5",
        "comm.2_3.0",
    ])
    @test read(so11, String) == orig["out1000"]
    disc = joinpath(out, "discarded", "segment_1")
    for name in ("lagr.7", "global.30", "esc.11", "event.35")
        @test read(joinpath(out, name), String) * read(joinpath(disc, name * ".tail"), String) ==
              orig[name]
    end
    @test read(joinpath(out, "esc.11"), String) == "         TTOT         BODY\n  1.00000E+00 1\n"
    for f in ("conf.3_3", "sev.83_3", "comm.2_2.5", "comm.2_3.0")
        @test isfile(joinpath(disc, f)) && !isfile(joinpath(out, f))
    end
end
