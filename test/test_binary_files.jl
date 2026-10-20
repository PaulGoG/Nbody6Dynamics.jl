# Readers of the binary, hierarchy and collision files (src/io/binary_files.jl).
# Every line below is synthetic: the layouts of the engine FORMATs with
# invented values, fixed-width fields built at their exact widths.

using Nbody6Dynamics:
    read_ks_pair_listing,
    read_all_ks_pair_listings,
    read_wide_pair_listing,
    read_all_wide_pair_listings,
    read_hierarchy_listing,
    read_collision_log,
    read_coalescence_log,
    read_cmbody_log,
    read_binary_stage_log,
    read_degenerate_log,
    read_escaped_binary_log,
    read_hierarchy_log,
    read_hardest_binary_log

"""Fixed-width test line: each `text => width` right-justified, as a Fortran descriptor writes it."""
function _bf_fields(fields::Pair{String,Int}...)::String
    for (text, width) in fields
        length(text) <= width || error("test field $(repr(text)) is wider than $width")
    end
    return join(lpad(text, width) for (text, width) in fields)
end

"""Write `lines` to `dir/name` (an empty vector gives an empty file) and return the path."""
function _bf_write(dir::AbstractString, name::AbstractString, lines::Vector{String})::String
    path = joinpath(dir, name)
    write(path, isempty(lines) ? "" : join(lines, "\n") * "\n")
    return path
end

"""Header line 1 of bdat.9 (FORMAT 30) with the given NPAIRS, MODEL and TTOT fields."""
_bf_bdat_header1(npairs::String, model::String, time::String) = _bf_fields(
    npairs => 4,
    model => 4,
    "3" => 4,
    "19975" => 6,
    "451" => 4,
    "1" => 4,
    "" => 2,
    time => 7,
    "0.97" => 7,
    "12.84" => 7,
    "0.265" => 7,
    "1.4" => 8,
    "-0.0000" => 9,
    "0.0000" => 9,
)

const _BF_BDAT_HEADER2 = _bf_fields(
    (
        [
            "-0.291737",
            "-0.000041",
            "-0.185012",
            "0.000000",
            "0.000000",
            "0.000000",
            "0.000000",
            "0.000000",
            "-0.004157",
            "0.096840",
        ] .=> 11
    )...,
)

const _BF_BDAT_HEADER3 = _bf_fields(
    (
        [
            "0.00000",
            "0.00000",
            "0.18311",
            "0.37122",
            "0.00000",
            "-0.29177",
            "0.00000",
            "-0.00416",
            "-0.39654",
            "0.98650",
            "0.32411",
            "-0.00003",
            "0.09684",
        ] .=> 10
    )...,
)

const _BF_BDAT_LABEL =
    " NAME(I1)    NAME(I2)    M1[M*]                    M2[M*]                    " * "EB[NB]"

"""bwdat.19 header (FORMAT 55) with the given TTOT field (E27.16)."""
_bf_bwdat_header(time::String) =
    " WIDE PAIRS    T TPHYS N " *
    _bf_fields(time => 27, "1.35" * "0"^14 * "E+00" => 27, "19975" => 12)

const _BF_BWDAT_LABEL =
    " NAME(I1)    NAME(I2)    M(I1)[M*]                 M(I2)[M*]                 " * "EB[NB]"

@testset "Binary output files" begin
    @testset "bdat.9 KS-pair listing" begin
        dir = mktempdir()
        rec1 =
            "          21          22  0.61234000000000000       0.43210000000000000       " *
            "-2.5000000000000000E-006  0.75000000000000000        35000.500000000000        " *
            "98.500000000000000       0.93000000000000000        7.8500000000000000        " *
            "        1           0   880.50000000000000        1.2500000000000000E-004   " *
            "1.1000000000000000E-004   7.8125000000000000E-003       20021 " *
            "-0.23860000000000000                0   1.5000000000000000        " *
            "1.5000000000000000        0.0000000000000000        0.0000000000000000     "
        rec2 =
            "        2019        4556  0.18500000000000000       0.44500000000000000        " *
            "6.8500000000000000E-006   14.250000000000000        0.0000000000000000       " *
            "-8.6500000000000000       0.23400000000000000        5.2100000000000000        " *
            "        0           0                       NaN   1.0950000000000000E-003   " *
            "3.2500000000000000E-005   6.1035156250000000E-005       22019  " *
            "-1.0315000000000000                0  0.61250000000000000       " *
            "0.61250000000000000        0.0000000000000000        0.0000000000000000     "
        bad = "          31          32  0.50000000000000000"
        path = _bf_write(
            dir,
            "bdat.9_10.0",
            String[
                _bf_bdat_header1("2", "17", "10.0"),
                _BF_BDAT_HEADER2,
                _BF_BDAT_HEADER3,
                _BF_BDAT_LABEL,
                rec1,
                # A repeated header block is skipped; the first one is parsed.
                _bf_bdat_header1("2", "18", "10.0"),
                _BF_BDAT_HEADER2,
                _BF_BDAT_HEADER3,
                _BF_BDAT_LABEL,
                rec2,
                bad,
            ],
        )
        listing = @test_logs (:warn, r"1 of 3 records skipped") read_ks_pair_listing(path)
        @test listing.n_pairs == 2
        @test listing.model == 17
        @test listing.n_run == 3
        @test listing.n_total == 19975
        @test listing.n_core == 451
        @test listing.n_merge == 1
        @test listing.time_nb == 10.0
        @test listing.rscale_nb == 0.97
        @test listing.rtide_nb == 12.84
        @test listing.core_radius_nb == 0.265
        @test listing.time_myr == 1.4
        @test listing.ecrit_nb == 0.0
        @test listing.energy_nb[1] == -0.291737
        @test listing.energy_nb[10] == 0.09684
        @test listing.kinetic_energy_nb == 0.18311
        @test listing.potential_energy_nb == 0.37122
        @test listing.ebin_nb == -0.29177
        @test listing.total_energy_nb == -0.39654
        @test listing.total_mass_nb == 0.9865
        @test listing.binary_mass_nb == 0.32411
        @test listing.ecoll_nb == 0.09684
        @test !listing.overflow
        @test length(listing.records) == 2

        r1, r2 = listing.records
        @test r1.name1 == 21 && r1.name2 == 22
        @test r1.mass1_msun == 0.61234
        @test r1.mass2_msun == 0.4321
        @test r1.binding_energy_nb == -2.5e-6
        @test r1.eccentricity == 0.75
        @test r1.period_d == 35000.5
        @test r1.semi_major_axis_au == 98.5
        @test r1.distance_density_centre_pc == 0.93
        @test r1.speed_cm_kms == 7.85
        @test r1.kstar1 == 1 && r1.kstar2 == 0
        @test r1.mean_motion_nb == 880.5
        @test r1.separation_nb == 1.25e-4
        @test r1.step_ks_nb == 1.1e-4
        @test r1.step_cm_nb == 7.8125e-3
        @test r1.name_cm == 20021
        @test r1.specific_energy_cm_nb == -0.2386
        @test r1.kstar_cm == 0
        @test r1.age1_myr == 1.5
        @test r1.epoch2_myr == 0.0
        # Hyperbolic pair: PB = 0 stored as NaN, ZN printed as NaN.
        @test isnan(r2.period_d)
        @test isnan(r2.mean_motion_nb)
        @test r2.semi_major_axis_au == -8.65
        @test r2.eccentricity == 14.25
        @test r2.name_cm == 22019

        empty_path = _bf_write(dir, "bdat.9_0.5", String[])
        empty_listing = read_ks_pair_listing(empty_path)
        @test isempty(empty_listing.records)
        @test empty_listing.n_pairs == 0
        @test isnan(empty_listing.time_nb)
        @test_throws ErrorException read_ks_pair_listing(joinpath(dir, "bdat.9_absent"))
    end

    @testset "bdat.9 time series" begin
        dir = mktempdir()
        for (name, npairs, time) in (
            ("bdat.9_0.0", "0", "0.0"),
            ("bdat.9_10.0", "****", "10.0"),
            ("bdat.9_2.5", "1", "2.5"),
        )
            _bf_write(
                dir,
                name,
                String[
                    _bf_bdat_header1(npairs, "17", time),
                    _BF_BDAT_HEADER2,
                    _BF_BDAT_HEADER3,
                    _BF_BDAT_LABEL,
                ],
            )
        end
        _bf_write(dir, "bdat.9_notes", String[])
        listings = read_all_ks_pair_listings(dir)
        @test length(listings) == 3
        @test [l.time_nb for l in listings] == [0.0, 2.5, 10.0]
        @test listings[2].n_pairs == 1
        # NPAIRS overflowed to asterisks: absent and flagged.
        @test listings[3].n_pairs == 0
        @test listings[3].overflow
        @test !listings[1].overflow
        @test isempty(read_all_ks_pair_listings(joinpath(dir, "absent")))
    end

    @testset "bwdat.19 wide-pair listing" begin
        dir = mktempdir()
        rec1 =
            "       10135       10136  0.93500000000000000       0.17500000000000000       " *
            "-2.8000000000000000E-007  0.89250000000000000        3002500.0000000000        " *
            "421.50000000000000        1.2207031250000000E-004   2.4414062500000000E-004  " *
            "0.59400000000000000        5.4700000000000000                1           0"
        rec2 =
            "        8481       20047  0.21000000000000000        56.750000000000000       " *
            "-1.7300000000000000E-006  0.85000000000000000        17031800.000000000        " *
            "4985.0000000000000        3.9062500000000000E-003   7.8125000000000000E-003  " *
            "0.62000000000000000        4.4400000000000000                0           0"
        rec3 =
            "      -10313        3883  0.38500000000000000        7.9750000000000000       " *
            "-4.2500000000000000E-003   2.0000000000000000E-003  0.27250000000000000        " *
            "124.75000000000000        6.1035156250000000E-005   7.8125000000000000E-003  " *
            "0.77500000000000000        3.2350000000000000                0           1"
        bad = "       10141       10142  0.50000000000000000"
        header = _bf_bwdat_header("1." * "0"^16 * "E+01")
        path = _bf_write(
            dir,
            "bwdat.19_10.0",
            String[header, _BF_BWDAT_LABEL, rec1, header, _BF_BWDAT_LABEL, rec2, rec3, bad],
        )
        listing = @test_logs (:warn, r"1 of 4 records skipped") read_wide_pair_listing(path)
        @test listing.time_nb == 10.0
        @test listing.time_myr == 1.35
        @test listing.n_total == 19975
        @test listing.nzero == 0
        @test !listing.overflow
        @test length(listing.records) == 3

        r1, r2, r3 = listing.records
        @test r1.name1 == 10135 && r1.name2 == 10136
        @test r1.mass1_msun == 0.935
        @test r1.binding_energy_nb == -2.8e-7
        @test r1.eccentricity == 0.8925
        @test r1.period_d == 3002500.0
        @test r1.semi_major_axis_au == 421.5
        # Step columns stored under their true meaning (labels shifted by one).
        @test r1.step1_nb == 1.220703125e-4
        @test r1.step2_nb == 2.44140625e-4
        @test r1.distance_density_centre_pc == 0.594
        @test r1.speed_cm_kms == 5.47
        @test r1.kstar1 == 1 && r1.kstar2 == 0
        @test !r1.member1_is_cm && !r1.member2_is_cm
        # Without NZERO only a negative (merged-hierarchy) name is recognised.
        @test !r2.member2_is_cm
        @test r3.member1_is_cm && !r3.member2_is_cm

        with_nzero =
            @test_logs (:warn, r"records skipped") read_wide_pair_listing(path; nzero = 20000)
        @test with_nzero.nzero == 20000
        @test with_nzero.records[2].member2_is_cm
        @test !with_nzero.records[2].member1_is_cm
        @test !with_nzero.records[1].member1_is_cm

        empty_listing = read_wide_pair_listing(_bf_write(dir, "bwdat.19_0.5", String[]))
        @test isempty(empty_listing.records)
        @test isnan(empty_listing.time_nb)
        @test_throws ErrorException read_wide_pair_listing(joinpath(dir, "bwdat.19_absent"))
    end

    @testset "bwdat.19 time series" begin
        dir = mktempdir()
        for (name, time) in (
            ("bwdat.19_0.0", "0." * "0"^16 * "E+00"),
            ("bwdat.19_10.0", "1." * "0"^16 * "E+01"),
            ("bwdat.19_2.5", "2.5" * "0"^15 * "E+00"),
        )
            _bf_write(dir, name, String[_bf_bwdat_header(time), _BF_BWDAT_LABEL])
        end
        listings = read_all_wide_pair_listings(dir; nzero = 20000)
        @test [l.time_nb for l in listings] == [0.0, 2.5, 10.0]
        @test all(l -> l.nzero == 20000, listings)
        @test isempty(read_all_wide_pair_listings(joinpath(dir, "absent")))
    end

    @testset "merger.84 hierarchy listing" begin
        dir = mktempdir()
        label =
            "BINDAT: Time[NB] NAME[I1] NAME[I3] K*[I1]" *
            "K*[I3] K*[IM] ECC0 ECC1 PERI(I3)/PCR " *
            "PERI(INCM)[RSM] P0[days] P1[days] SEMI1[NB] "
        rec1 =
            label * _bf_fields(
                "5.00000E+00" => 12,
                "10311" => 12,
                "10312" => 12,
                "0" => 4,
                "1" => 4,
                "10" => 4,
                "0.003" => 7,
                "0.412" => 7,
                "2150.25000" => 12,
                "5.12500" => 12,
                "2.75000E-01" => 14,
                "2.90000E+05" => 14,
                "8.50000E-04" => 14,
            )
        rec2 =
            label * _bf_fields(
                "5.50000E+00" => 12,
                "10311" => 12,
                "10312" => 12,
                "13" => 4,
                "13" => 4,
                "13" => 4,
                "0.001" => 7,
                "0.835" => 7,
                "253.12500" => 12,
                "*"^12 => 12,
                "2.75000E-01" => 14,
                "1.76250E+05" => 14,
                "6.07500E-04" => 14,
            )
        bad = "BINDAT: Time[NB] NAME[I1] garbage"
        path = _bf_write(dir, "merger.84", String[rec1, "", rec2, bad])
        recs = @test_logs (:warn, r"1 of 3 records skipped") read_hierarchy_listing(path)
        @test length(recs) == 2
        r1, r2 = recs
        @test r1.time_nb == 5.0
        @test r1.name_primary == 10311
        @test r1.name_ghost == 10312
        @test r1.kstar_primary == 0
        @test r1.kstar_ghost == 1
        @test r1.kstar_inner_cm == 10
        @test r1.eccentricity_inner == 0.003
        @test r1.eccentricity_outer == 0.412
        @test r1.pericentre_outer_over_pcrit == 2150.25
        @test r1.pericentre_inner_over_radius == 5.125
        @test r1.period_inner_d == 0.275
        @test r1.period_outer_d == 2.9e5
        @test r1.semi_major_axis_outer_nb == 8.5e-4
        @test !r1.overflow
        # RM printed as asterisks.
        @test r2.overflow
        @test isnan(r2.pericentre_inner_over_radius)
        @test r2.period_outer_d == 1.7625e5
        @test r2.kstar_inner_cm == 13

        @test isempty(read_hierarchy_listing(_bf_write(dir, "merger.empty", String[])))
        @test_throws ErrorException read_hierarchy_listing(joinpath(dir, "absent"))
    end

    @testset "coll.13 collision log" begin
        dir = mktempdir()
        model =
            "      MODEL:    RBAR =   9.85000000000000000E-01  <M>[M*] =   " *
            "5.95000000000000000E-01  M1[M*] =   0.00000000000000000E+00  TSCALE =   " *
            "1.35000000000000000E-01  NB0 =        6650  N0 =       19950"
        label = "                  TIME[NB]    NAME(I1)    NAME(I2)      K*(I1)      K*(I2)"
        rec1 =
            "   5.8000000000000000E-006        3105        3106           1           0" *
            "           1  0.38900000000000000        8.8000000000000000E-002  " *
            "0.45000000000000000        2.6500000000000000E-002  0.36500000000000000       " *
            "0.14000000000000000        1.5800000000000000        4.8000000000000000E-002  " *
            "0.99750000000000000        15.750000000000000     "
        rec2 =
            "   7.6500000000000000E-006        7069        7070           0           0" *
            "           0  0.13800000000000000       0.10900000000000000       " *
            "0.24700000000000000        0.0000000000000000       0.17500000000000000       " *
            "0.14200000000000000        9.9200000000000000       0.18550000000000000       " *
            "-2.5000000000000000E-002  -23.050000000000000     "
        bad = "   8.0000000000000000E-006        7071"
        path = _bf_write(
            dir,
            "coll.13",
            String["", model, "", "", label, "", rec1, "", model, "", "", label, "", rec2, bad],
        )
        recs = @test_logs (:warn, r"1 of 3 records skipped") read_collision_log(path)
        @test length(recs) == 2
        r1, r2 = recs
        @test r1.time_nb == 5.8e-6
        @test r1.name1 == 3105 && r1.name2 == 3106
        @test r1.kstar1 == 1 && r1.kstar2 == 0 && r1.kstar_product == 1
        @test r1.mass1_msun == 0.389
        @test r1.mass2_msun == 0.088
        @test r1.mass_product_msun == 0.45
        @test r1.mass_lost_msun == 0.0265
        @test r1.radius1_rsun == 0.365
        @test r1.radius2_rsun == 0.14
        @test r1.distance_over_core_radius == 1.58
        @test r1.separation_rsun == 0.048
        @test r1.eccentricity == 0.9975
        @test r1.period_d == 15.75
        @test r2.eccentricity == -0.025
        @test r2.period_d == -23.05

        @test isempty(read_collision_log(_bf_write(dir, "coll.empty", String[])))
        @test_throws ErrorException read_collision_log(joinpath(dir, "absent"))
    end

    @testset "coal.24 coalescence log" begin
        dir = mktempdir()
        model =
            "      MODEL:    RBAR =   1.01250000000000000E+00  <M> =   " *
            "5.84375000000000000E-01  M1 =   8.79375000000000000E+01  TSCALE =   " *
            "2.00000000000000000E-01  NB =        3330  N0 =        9990"
        label = " WHICH1 IQCOLL TTOT           NAME(I1,I2)        K*(I1,I2,I) M(I1)[NB]"
        rec1 =
            "  ROCHE " * _bf_fields(
                "3" => 3,
                "1.6250000000E+01" => 17,
                "1851" => 10,
                "3378" => 10,
                "1" => 4,
                "0" => 4,
                "1" => 4,
                (
                    [
                        "8.75000E+01",
                        "1.50000E-01",
                        "6.62500E-07",
                        "-2.00000E-03",
                        "6.61250E-07",
                        "-2.91250E-01",
                        "8.12500E-15",
                        "2.00500E+00",
                        "8.75000E+01",
                        "1.50000E-01",
                        "8.76500E+01",
                        "0.00000E+00",
                        "2.44375E+01",
                        "0.00000E+00",
                        "2.98750E+01",
                        "0.00000E+00",
                        "9.62500E-02",
                        "5.23125E-01",
                    ] .=> 13
                )...,
            )
        rec2 =
            " HYPERB " * _bf_fields(
                "-1" => 3,
                "4.2500000000E+01" => 17,
                "5801" => 10,
                "416" => 10,
                "13" => 4,
                "14" => 4,
                "14" => 4,
                (
                    [
                        "1.40000E+00",
                        "1.30000E+00",
                        "1.25000E-08",
                        "-5.00000E-01",
                        "-2.50000E-08",
                        "1.75000E-02",
                        "*"^13,
                        "-1.00000E+00",
                        "1.40000E+00",
                        "1.30000E+00",
                        "2.70000E+00",
                        "0.00000E+00",
                        "1.25000E-05",
                        "0.00000E+00",
                        "3.00000E-05",
                        "3.12500E+01",
                        "1.50000E-01",
                        "6.25000E-01",
                    ] .=> 13
                )...,
            )
        bad = "  CHAIN   5 1.0000000000E+00"
        path = _bf_write(
            dir,
            "coal.24",
            String["", model, "", "", label, rec1, "", model, "", "", label, rec2, bad],
        )
        recs = @test_logs (:warn, r"1 of 3 records skipped") read_coalescence_log(path)
        @test length(recs) == 2
        r1, r2 = recs
        @test r1.channel === :roche
        @test r1.iqcoll == 3
        @test r1.time_nb == 16.25
        @test r1.name1 == 1851 && r1.name2 == 3378
        @test r1.kstar1 == 1 && r1.kstar2 == 0 && r1.kstar_product == 1
        # The columns labelled M(I1)[NB], M(I2)[NB] are solar masses.
        @test r1.mass1_bse_msun == 87.5
        @test r1.mass2_bse_msun == 0.15
        @test r1.separation_nb == 6.625e-7
        @test r1.eccentricity == -2.0e-3
        @test r1.semi_major_axis_nb == 6.6125e-7
        @test r1.binding_energy_nb == -0.29125
        @test r1.tidal_potential_correction_nb == 8.125e-15
        @test r1.period_d == 2.005
        @test r1.mass1_msun == 87.5
        @test r1.mass_product_msun == 87.65
        @test r1.mass_lost_msun == 0.0
        @test r1.radius_survivor_rsun == 24.4375
        @test r1.collision_distance_rsun == 29.875
        @test r1.velocity_infinity_kms == 0.0
        @test r1.distance_density_centre_nb == 0.09625
        @test r1.speed_nb == 0.523125
        @test !r1.overflow
        @test !hasproperty(r1, :radius2_rsun)
        @test r2.channel === :hyperbolic
        @test r2.iqcoll == -1
        @test r2.kstar_product == 14
        @test r2.overflow
        @test isnan(r2.tidal_potential_correction_nb)
        @test r2.velocity_infinity_kms == 31.25
        @test r2.speed_nb == 0.625

        @test isempty(read_coalescence_log(_bf_write(dir, "coal.empty", String[])))
        @test_throws ErrorException read_coalescence_log(joinpath(dir, "absent"))
    end

    @testset "cmbody.26 coalescence log" begin
        dir = mktempdir()
        model =
            "    MODEL:    RBAR =  1.0  <M> =  0.58  M1 =  87.9  TSCALE =  0.20" *
            "  NB =3330  N0 =  9990"
        label =
            "          TIME[NB]           NAME(I1)    NAME(I2)     K*(I1)      K*(I2)" *
            "     IQCOLL"
        rec1 =
            "   12.500000000000000             2041        2042          13          13" *
            "           0   1.2500000000000000        1.3000000000000000        " *
            "0.0000000000000000        1.5000000000000000E-005   0.0000000000000000        " *
            "1.7500000000000000        3.0000000000000000E-005  0.99500000000000000        " *
            "4.5000000000000000E-003   1.2500000000000000E-005"
        rec2 =
            "   13.250000000000000              511         512          10          10" *
            "           3  0.60000000000000000       0.55000000000000000        " *
            "0.0000000000000000        1.2000000000000000E-002   0.0000000000000000       " *
            "0.75000000000000000        2.5000000000000000E-002  -0.50000000000000000       " *
            "-1.6700000000000000E-002                       NaN"
        bad = "   14.000000000000000              513"
        path = _bf_write(dir, "cmbody.26", String["", model, "", "", label, rec1, label, rec2, bad])
        recs = @test_logs (:warn, r"1 of 3 records skipped") read_cmbody_log(path)
        @test length(recs) == 2
        r1, r2 = recs
        @test r1.time_nb == 12.5
        @test r1.name1 == 2041 && r1.name2 == 2042
        @test r1.kstar1 == 13 && r1.kstar2 == 13
        @test r1.iqcoll == 0
        @test r1.mass1_msun == 1.25
        @test r1.mass2_msun == 1.3
        @test r1.mass_lost_msun == 0.0
        @test r1.radius1_rsun == 1.5e-5
        @test r1.distance_over_core_radius == 1.75
        @test r1.collision_distance_rsun == 3.0e-5
        @test r1.eccentricity == 0.995
        @test r1.semi_major_axis_rsun == 4.5e-3
        @test r1.period_d == 1.25e-5
        @test r2.iqcoll == 3
        @test r2.semi_major_axis_rsun == -1.67e-2
        @test isnan(r2.period_d)

        @test isempty(read_cmbody_log(_bf_write(dir, "cmbody.empty", String[])))
        @test_throws ErrorException read_cmbody_log(joinpath(dir, "absent"))
    end

    @testset "binev.17 binary-stage log" begin
        dir = mktempdir()
        label =
            "Time[Myr]                 NAME(I1)    NAME(I2)    K*(I1)      K*(I2)      " *
            "K*(ICM)     "
        rec1 =
            "   3.2500000000000000E-006        7293        7294           1           0" *
            "           0   1.9875000000000000        8.2000000000000000E-002   " *
            "1.6050000000000000       0.14400000000000000        3.8375000000000000       " *
            "0.99975000000000000        67.625000000000000        44.750000000000000" *
            "               -2"
        rec2 =
            "   3.3125000000000000             1851        3378           1           0" *
            "          11   87.843750000000000       0.15000000000000000        " *
            "24.437500000000000       0.17100000000000000        9.7656250000000000E-002   " *
            "1.9531250000000000E-003   29.812500000000000        2.0078125000000000" *
            "                0"
        bad = "   3.5000000000000000             1853        3380           1"
        path = _bf_write(dir, "binev.17", String[label, "", rec1, label, "", rec2, bad])
        recs = @test_logs (:warn, r"1 of 3 records skipped") read_binary_stage_log(path)
        @test length(recs) == 2
        r1, r2 = recs
        @test r1.time_myr == 3.25e-6
        @test r1.name1 == 7293 && r1.name2 == 7294
        @test r1.kstar1 == 1 && r1.kstar2 == 0 && r1.kstar_cm == 0
        @test r1.mass1_msun == 1.9875
        @test r1.mass2_msun == 0.082
        @test r1.radius1_rsun == 1.605
        @test r1.radius2_rsun == 0.144
        @test r1.distance_density_centre_pc == 3.8375
        @test r1.eccentricity == 0.99975
        @test r1.semi_major_axis_rsun == 67.625
        @test r1.period_d == 44.75
        @test r1.iqcoll == -2
        @test r2.time_myr == 3.3125
        @test r2.kstar_cm == 11
        @test r2.distance_density_centre_pc == 9.765625e-2
        @test r2.iqcoll == 0

        @test isempty(read_binary_stage_log(_bf_write(dir, "binev.empty", String[])))
        @test_throws ErrorException read_binary_stage_log(joinpath(dir, "absent"))
    end

    @testset "degen.4 degenerate-binary log" begin
        dir = mktempdir()
        model =
            "      MODEL:    RBAR =   1.01250000000000000E+00  <M>[M*] =   " *
            "5.84375000000000000E-01  M1[M*] =   8.79375000000000000E+01  TSCALE =   " *
            "2.00000000000000000E-01  NB0 =        3330  N0 =        9990"
        label =
            "       ICASE                 Time[Myr]                  SEMI[AU]" *
            "                       ECC"
        rec1 =
            "           3   12.250000000000000        0.12500000000000000       " *
            "0.25000000000000000        4.5000000000000000        15.500000000000000       " *
            "0.87500000000000000        1.3750000000000000       0.62500000000000000" *
            "                13           1           0        1201        1202"
        rec2 =
            "           4   20.500000000000000        2.5000000000000000E-002   " *
            "0.0000000000000000        3.2500000000000000        1.1250000000000000       " *
            "0.43750000000000000       0.81250000000000000       0.60000000000000000" *
            "                11          10          10         733         734"
        bad = "           5   21.000000000000000"
        path = _bf_write(
            dir,
            "degen.4",
            String["", model, "", "", label, "", rec1, "", rec2, label, bad],
        )
        recs = @test_logs (:warn, r"1 of 3 records skipped") read_degenerate_log(path)
        @test length(recs) == 2
        r1, r2 = recs
        @test r1.icase == 3
        @test r1.time_myr == 12.25
        @test r1.semi_major_axis_au == 0.125
        @test r1.eccentricity == 0.25
        @test r1.pericentre_over_radius == 4.5
        @test r1.period_d == 15.5
        @test r1.cm_distance_pc == 0.875
        @test r1.mass1_msun == 1.375
        @test r1.mass2_msun == 0.625
        @test r1.kstar1 == 13 && r1.kstar2 == 1 && r1.kstar_cm == 0
        @test r1.name1 == 1201 && r1.name2 == 1202
        @test r2.icase == 4
        @test r2.eccentricity == 0.0
        @test r2.kstar_cm == 10
        @test r2.name2 == 734

        @test isempty(read_degenerate_log(_bf_write(dir, "degen.empty", String[])))
        @test_throws ErrorException read_degenerate_log(joinpath(dir, "absent"))
    end

    @testset "escbin.31 escaped-binary log" begin
        dir = mktempdir()
        header = " TTOT         BODY1,2                   NAME1,2            KW1,2  M1,2 [Sol]"
        reals1 = [
            "1.87500E-01",
            "1.68750E-01",
            "2.03125E-01",
            "1.87500E-01",
            "8.18750E-03",
            "6.50000E-03",
            "3.58750E+00",
            "3.57500E+00",
            "2.80000E+00",
            "2.80000E+00",
            "0.00000E+00",
            "0.00000E+00",
        ]
        rec1 =
            " " * _bf_fields(
                "1.45000E+01" => 13,
                "3.18750E-05" => 13,
                "2.87500E-05" => 13,
                "1025" => 10,
                "1026" => 10,
                "0" => 4,
                "1" => 4,
                (reals1 .=> 13)...,
            )
        rec2 =
            " " * _bf_fields(
                "3.20000E+01" => 13,
                "2.03125E-05" => 13,
                "1.81250E-05" => 13,
                "597" => 10,
                "*"^10 => 10,
                "0" => 4,
                "0" => 4,
                (reals1 .=> 13)...,
            )
        bad = " " * "  1.00000E+00"
        path = _bf_write(dir, "escbin.31", String[header, rec1, header, rec2, bad])
        recs = @test_logs (:warn, r"1 of 3 records skipped") read_escaped_binary_log(path)
        @test length(recs) == 2
        r1, r2 = recs
        @test r1.time_nb == 14.5
        @test r1.mass1_nb == 3.1875e-5
        @test r1.mass2_nb == 2.875e-5
        @test r1.name1 == 1025 && r1.name2 == 1026
        @test r1.kstar1 == 0 && r1.kstar2 == 1
        @test r1.mass1_msun == 0.1875
        @test r1.mass2_msun == 0.16875
        @test r1.radius1_rsun == 0.203125
        @test r1.radius2_rsun == 0.1875
        @test r1.luminosity1_lsun == 8.1875e-3
        @test r1.luminosity2_lsun == 6.5e-3
        # The columns labelled TEFF1,2 [K] are log10 T_eff.
        @test r1.log_teff1 == 3.5875
        @test r1.log_teff2 == 3.575
        @test r1.age1_myr == 2.8
        @test r1.epoch1_myr == 0.0
        @test !r1.overflow
        @test r2.time_nb == 32.0
        @test r2.name1 == 597
        @test r2.name2 == 0
        @test r2.overflow

        @test isempty(read_escaped_binary_log(_bf_write(dir, "escbin.empty", String[])))
        @test_throws ErrorException read_escaped_binary_log(joinpath(dir, "absent"))
    end

    @testset "hiarch.12 hierarchy log" begin
        dir = mktempdir()
        model =
            "      MODEL:    RBAR =   9.85000000000000000E-01  <M>[M*] =   " *
            "5.95000000000000000E-01  M1[M*] =   1.37250000000000000E+02  TSCALE =   " *
            "1.35000000000000000E-01  NB0 =        6650  N0 =       19950"
        new_label =
            "NEW HIAR  Time SEMI0 SMEI1 ECC1 " *
            "PERI0 PERI0M P1/P0 M(INCM)/M(I3) PCR/SEMI0 " *
            "M(INCM)/<M> MR INA NAME(I1) NAME(I2) NAME(I3) " *
            "K*(INCM) ECC0 ECCMIN ECCMAX K*(I1) K*(I2) RSM "
        new_line =
            new_label * _bf_fields(
                "4.65625000000E+00" => 20,
                "3.75000E-06" => 14,
                "8.96250E-04" => 14,
                "0.41" => 6,
                "1.8750E-06" => 13,
                "2.6250E-07" => 13,
                "7.8750E+02" => 13,
                "4.8750E-02" => 13,
                "1.0325E+02" => 13,
                "6.4625E-01" => 13,
                "1.5375E+00" => 13,
                "1.1125E+02" => 13,
                "10311" => 12,
                "10312" => 12,
                "3881" => 12,
                "0" => 12,
                "0.52" => 6,
                "0.94" => 6,
                "0.95" => 6,
                "0" => 12,
                "1" => 12,
                "5.78125E-09" => 14,
            )
        end_label =
            "END HIAR  Time RI/RC SEMI0 ECC0 PERI0 P0F/P0I " *
            "RC/RSCALE GAMMA(INCM) NKI NKF NPAIRS NAME(I2)"
        end_line(time, nk) =
            end_label * _bf_fields(
                time => 20,
                "1.85938E+00" => 14,
                "9.25000E-04" => 14,
                "0.39" => 6,
                "1.0150E+02" => 13,
                "8.2375E+02" => 13,
                "3.1750E-01" => 13,
                "1.0050E-01" => 13,
                nk => 4,
                "2" => 4,
                "3778" => 12,
                "3881" => 12,
            )
        bad = "NEW HIAR  Time SEMI0 SMEI1"
        path = _bf_write(
            dir,
            "hiarch.12",
            String[
                "",
                model,
                "",
                "",
                "",
                new_line,
                end_line("4.67187500000E+00", "1503"),
                "",
                model,
                "",
                end_line("4.68750000000E+00", "****"),
                bad,
            ],
        )
        recs = @test_logs (:warn, r"1 of 4 records skipped") read_hierarchy_log(path)
        @test length(recs) == 3
        rn, re, ro = recs
        @test rn.kind === :new
        @test rn.time_nb == 4.65625
        @test rn.semi_major_axis_inner_nb == 3.75e-6
        @test rn.semi_major_axis_outer_nb == 8.9625e-4
        @test rn.eccentricity_outer == 0.41
        @test rn.pericentre_inner_nb == 1.875e-6
        @test rn.pericentre_inner_min_nb == 2.625e-7
        @test rn.period_ratio == 787.5
        @test rn.mass_ratio_inner_outer == 4.875e-2
        @test rn.pcrit_over_semi_inner == 103.25
        @test rn.mass_inner_over_mean == 0.64625
        @test rn.mass_ratio_components == 1.5375
        @test rn.inclination_deg == 111.25
        @test rn.name1 == 10311 && rn.name2 == 10312 && rn.name_outer == 3881
        @test rn.kstar_inner_cm == 0
        @test rn.eccentricity_inner == 0.52
        @test rn.eccentricity_jz == 0.94
        @test rn.eccentricity_max == 0.95
        @test rn.kstar1 == 0 && rn.kstar2 == 1
        @test rn.radius_max_nb == 5.78125e-9
        @test isnan(rn.distance_over_core_radius)
        @test rn.n_pairs == 0
        @test !rn.overflow

        @test re.kind === :end
        @test re.time_nb == 4.671875
        @test re.distance_over_core_radius == 1.85938
        # SEMI0, ECC0 of END HIAR are the outer orbit.
        @test re.semi_major_axis_outer_nb == 9.25e-4
        @test re.eccentricity_outer == 0.39
        @test isnan(re.semi_major_axis_inner_nb)
        @test isnan(re.eccentricity_inner)
        @test re.pericentre_outer_over_apocentre_inner == 101.5
        @test re.period_ratio == 823.75
        @test re.core_radius_over_rscale == 0.3175
        # NKI glued to GAMMA is cut at its width.
        @test re.gamma == 0.1005
        @test re.n_inner_periods == 1503
        @test re.n_outer_periods == 2
        @test re.n_pairs == 3778
        @test re.name_outer == 3881
        @test re.kstar1 == KSTAR_ABSENT
        @test !re.overflow
        # NKI overflowed to asterisks.
        @test ro.n_inner_periods == 0
        @test ro.overflow
        @test ro.gamma == 0.1005
        @test ro.n_outer_periods == 2

        @test isempty(read_hierarchy_log(_bf_write(dir, "hiarch.empty", String[])))
        @test_throws ErrorException read_hierarchy_log(joinpath(dir, "absent"))
    end

    @testset "hbin.39 hardest-binary log" begin
        dir = mktempdir()
        label =
            " BINARY:   Time[Myr] NAME(I1) NAME(I2) " *
            "K*(ICM) NP ECC SEMI[NB] P[days] EB[NB] EM[NB] "
        rec1 =
            label * _bf_fields(
                "0." * "0"^17 * "E+00" => 26,
                "4373" => 12,
                "4374" => 12,
                "0" => 4,
                "0" => 6,
                "6.062500E-01" => 15,
                "1.525000E-05" => 15,
                "2.875000E+02" => 15,
                "-1.375000E-03" => 15,
                "-1.975000E-01" => 15,
            )
        rec2 =
            label * _bf_fields(
                "5." * "0"^17 * "E-01" => 26,
                "13247" => 12,
                "13248" => 12,
                "10" => 4,
                "*"^6 => 6,
                "1.996094E-03" => 15,
                "2.980000E-08" => 15,
                "2.187500E-01" => 15,
                "-1.125000E-02" => 15,
                "-1.965000E-01" => 15,
            )
        bad = " BINARY:   Time[Myr] NAME(I1)"
        path = _bf_write(dir, "hbin.39", String[rec1, rec2, "", bad])
        recs = @test_logs (:warn, r"1 of 3 records skipped") read_hardest_binary_log(path)
        @test length(recs) == 2
        r1, r2 = recs
        # Labelled Time[Myr], stored as the N-body time it is.
        @test r1.time_nb == 0.0
        @test r1.name1 == 4373 && r1.name2 == 4374
        @test r1.kstar_cm == 0
        @test r1.n_neighbours == 0
        @test r1.eccentricity == 0.60625
        @test r1.semi_major_axis_nb == 1.525e-5
        @test r1.period_d == 287.5
        @test r1.binding_energy_nb == -1.375e-3
        @test r1.total_energy_nb == -0.1975
        @test !r1.overflow
        @test r2.time_nb == 0.5
        @test r2.kstar_cm == 10
        @test r2.n_neighbours == 0
        @test r2.overflow
        @test r2.binding_energy_nb == -1.125e-2

        @test isempty(read_hardest_binary_log(_bf_write(dir, "hbin.empty", String[])))
        @test_throws ErrorException read_hardest_binary_log(joinpath(dir, "absent"))
    end
end
