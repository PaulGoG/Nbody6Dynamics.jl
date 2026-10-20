# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Fixed-width reading of FORMAT records" begin
    lit = Nbody6Dynamics.fmt_lit
    int = Nbody6Dynamics.fmt_int
    real = Nbody6Dynamics.fmt_real
    skip = Nbody6Dynamics.fmt_skip
    read_fixed = Nbody6Dynamics.read_fixed

    # NEW KSREG head (ksinit.F FORMAT 60): a ten-digit name fills its I10 field
    # and is glued to the label before it
    head = vcat(
        lit(" NEW KSREG   TIME[NB]"),
        real(1, 17),
        lit(" NM1,2,S="),
        int(3, 10),
        lit(" KW1,2,S="),
        int(3, 4),
        lit(" e,a,eb[NB]="),
        real(2, 12),
        real(1, 11),
    )
    line =
        " NEW KSREG   TIME[NB] 1.2500000000E+01 NM1,2,S=1234567890        14     20013" *
        " KW1,2,S=  14  13  -1 e,a,eb[NB]=  9.5761E-01-6.1739E-03 -1.134E-07"
    r = read_fixed(line, head)
    @test r.values ==
          [12.5, 1234567890.0, 14.0, 20013.0, 14.0, 13.0, -1.0, 0.95761, -6.1739e-3, -1.134e-7]
    @test !r.overflow

    # Overflowed and empty fields read as NaN; overflow is flagged
    r2 = read_fixed(" N=******  T=", vcat(lit(" N="), int(1, 6), lit("  T="), real(1, 9)))
    @test isnan(r2.values[1]) && isnan(r2.values[2]) && r2.overflow

    # A literal that does not match is an error naming it
    @test_throws ArgumentError read_fixed(" END KSREG", lit(" NEW KSREG"))
    @test read_fixed(" NEW  KSREG", lit(" NEW KSREG "); strict_literals = false).values == Float64[]

    # Skipped fields and a start column
    r3 = read_fixed("xxxx   42 1.5D+00", vcat(skip(4), int(1, 5), real(1, 8)))
    @test r3.values == [42.0, 1.5]
    @test read_fixed("ab   7", int(1, 4); start = 3).values == [7.0]

    # Exponents: D form, and three-digit exponents written without the letter
    pr = Nbody6Dynamics._parse_fortran_real
    @test pr("3.31324D-01") == 0.331324
    @test pr("1.08726-157") == 1.08726e-157
    @test pr("-1.5+100") == -1.5e100
    @test pr("0.58516901771347241") == 0.58516901771347241
    @test Nbody6Dynamics._fixed_int(NaN, -99) == -99
    @test Nbody6Dynamics._fixed_int(7.0) == 7
end
