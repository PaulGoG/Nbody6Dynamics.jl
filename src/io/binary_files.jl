# =============================================================================
# Readers of the binary, hierarchy and collision files of the engine
# =============================================================================
#
#   file           writer   layout
#   bdat.9_<t>     BINDAT   three fixed-width header lines, a label line, one
#                           list-directed line per KS pair
#   bwdat.19_<t>   BINDAT   fixed-width header, label line, list-directed pairs
#   merger.84      BINDAT   fixed width, label text inside every record
#   coll.13        MIX      MODEL and label header, list-directed records
#   coal.24        COAL     MODEL and label header, fixed-width records
#   cmbody.26      CMBODY   MODEL and label header (practically never written),
#                           list-directed records
#   binev.17       BINEV    label header, list-directed records
#   degen.4        DEGEN    MODEL and label header, list-directed records
#   escbin.31      ESCAPE   label header (usually absent), fixed-width records
#   hiarch.12      HIARCH   MODEL header, fixed-width NEW HIAR / END HIAR records
#   hbin.39        ADJUST   fixed width, label text inside every record
#
# The per-time files bdat.9_<t> and bwdat.19_<t> are rewritten when a restart
# reaches the same time; the others are opened for appending, so a restarted
# run repeats every header guarded by a FIRST flag alone. Header lines are
# therefore skipped wherever they occur. Empty lines (the `/` of a FORMAT)
# are never part of a record.

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

"""Total width of a FORMAT, in characters."""
_binfile_width(items::AbstractVector{FortranItem})::Int =
    sum(item.width for item in items; init = 0)

"""Number of value fields (`:int`, `:real`) of a FORMAT."""
_binfile_nvalues(items::AbstractVector{FortranItem})::Int =
    count(item -> item.kind === :int || item.kind === :real, items)

"""
Fixed-width record of `line` by the FORMAT `items`, from character `start`.
The last field of every record below ends in a digit, so a line shorter than
the FORMAT is truncated: `ArgumentError`.
"""
function _binfile_fixed(
    line::AbstractString,
    items::AbstractVector{FortranItem};
    start::Int = 1,
)::FortranRecord
    need = start - 1 + _binfile_width(items)
    length(line) >= need ||
        throw(ArgumentError("record of $(length(line)) characters, its FORMAT needs $need"))
    return read_fixed(line, items; start = start)
end

"""Values of a list-directed record; `ArgumentError` unless there are exactly `n`."""
function _binfile_tokens(line::AbstractString, n::Int, file::AbstractString)
    tokens = split(line)
    length(tokens) == n ||
        throw(ArgumentError("$file record with $(length(tokens)) values, expected $n"))
    return tokens
end

"""
Records of the file `path`, one per non-empty line that is not a header line
(`is_header(lstrip(line))`), built by `parser(line)`. Lines whose parser
throws `ArgumentError` are skipped and counted in one warning per file.
"""
function _binfile_read(
    parser::F,
    ::Type{T},
    path::AbstractString,
    file::AbstractString,
    is_header::H,
)::Vector{T} where {F,T,H}
    isfile(path) || error("$file file not found: $path")
    records = T[]
    n_data = 0
    n_skipped = 0
    for line in eachline(path)
        isempty(strip(line)) && continue
        is_header(lstrip(line)) && continue
        n_data += 1
        try
            push!(records, parser(line))
        catch e
            e isa ArgumentError || rethrow()
            n_skipped += 1
            @debug "Skipping malformed $file record" line exception = e
        end
    end
    n_skipped > 0 && @warn "$file: $n_skipped of $n_data records skipped (malformed)" path
    return records
end

"""Paths of the files `<prefix><t>` in `dir`, sorted by the time `t` in the name."""
function _binfile_time_files(dir::AbstractString, prefix::AbstractString)::Vector{String}
    isdir(dir) || return String[]
    stamped = Tuple{Float64,String}[]
    for f in readdir(dir)
        startswith(f, prefix) || continue
        t = tryparse(Float64, strip(chopprefix(f, prefix)))
        t === nothing || push!(stamped, (t, joinpath(dir, f)))
    end
    sort!(stamped; by = first)
    return String[last(s) for s in stamped]
end

# ---------------------------------------------------------------------------
# bdat.9_<t>: KS-pair data bank
# ---------------------------------------------------------------------------

"""Values per record of `bdat.9_<t>` (`bindat.f`, 19 labels for 23 values)."""
const _BDAT_RECORD_VALUES = 23

"""
FORMAT 30 of `bindat.f`: `(3I4,I6,2I4,2X,F7.1,2F7.2,F7.3,F8.1,2F9.4)`; the
last field, AS(7) = 0, is not read.
"""
const _BDAT_HEADER1_FORMAT = vcat(
    fmt_int(3, 4),
    fmt_int(1, 6),
    fmt_int(2, 4),
    fmt_skip(2),
    fmt_real(1, 7),
    fmt_real(2, 7),
    fmt_real(1, 7),
    fmt_real(1, 8),
    fmt_real(1, 9),
    fmt_skip(9),
)

"""FORMAT 35 of `bindat.f`: `(10F11.6)`."""
const _BDAT_HEADER2_FORMAT = fmt_real(10, 11)

"""FORMAT 40 of `bindat.f`: `(13F10.5)`."""
const _BDAT_HEADER3_FORMAT = fmt_real(13, 10)

"""
    KSPairRecord

One KS pair of a `bdat.9_<t>` listing, written list-directed by `bindat.f`
(23 values on one line), from [`read_ks_pair_listing`](@ref).

| # | Label | Field | Unit |
|---|---|---|---|
| 1–2 | NAME(I1), NAME(I2) | `name1`, `name2` | — |
| 3–4 | M1[M*], M2[M*] | `mass1_msun`, `mass2_msun` | M☉ |
| 5 | EB[NB] | `binding_energy_nb` (μH) | NB |
| 6 | ECC | `eccentricity` | — |
| 7 | P[Days] | `period_d` | d |
| 8 | SEMI[AU] | `semi_major_axis_au` | AU |
| 9 | RI[PC] | `distance_density_centre_pc` (‖x_cm − x_dens‖) | pc |
| 10 | VI[km/s] | `speed_cm_kms` (‖v_cm‖, not relative to the cluster) | km s⁻¹ |
| 11–12 | K*(I1), K*(I2) | `kstar1`, `kstar2` | — |
| 13 | ZN[NB] | `mean_motion_nb` (√((m₁ + m₂)/a³)) | NB |
| 14 | RP[NB] | `separation_nb` (R = ΣU²) | NB |
| 15 | STEP(I1)[NB] | `step_ks_nb` | NB |
| 16 | STEP(ICM)[NB] | `step_cm_nb` | NB |
| 17 | NAME(ICM) | `name_cm` | — |
| 18 | ECM[NB] | `specific_energy_cm_nb` (½v² − φ, tidal term included) | NB |
| 19 | K*(ICM) | `kstar_cm` | — |
| 20–21 | (none) | `age1_myr`, `age2_myr` (TIME·T* − EPOCH) | Myr |
| 22–23 | (none) | `epoch1_myr`, `epoch2_myr` | Myr |

Hazards. A hierarchy row (`name_cm < 0`, `kstar_cm = −10` as printed) mixes
two orbits: `binding_energy_nb`, `eccentricity` and `period_d` describe the
inner binary (the eccentricity from ‖XREL‖ only, exact at peri- and apocentre
alone), while the masses (inner c.m. body and outer body), `semi_major_axis_au`,
`mean_motion_nb` and `separation_nb` describe the outer KS pair. A ghost row
has zero masses and a NaN mean motion; when its pair is missing from the
merger table, columns 5–7 and 18 hold values of the previous BINDAT call for
the same slot [UNSURE in the dossier; not corrected]. The engine prints
P = 0 for a < 0: stored as `NaN` (correction); the mean motion of such a pair
is NaN as printed. The ages use TIME rather than TTOT and are wrong once
TOFF > 0 (not corrected).
"""
Base.@kwdef struct KSPairRecord
    name1::Int
    name2::Int
    mass1_msun::Float64
    mass2_msun::Float64
    binding_energy_nb::Float64
    eccentricity::Float64
    period_d::Float64
    semi_major_axis_au::Float64
    distance_density_centre_pc::Float64
    speed_cm_kms::Float64
    kstar1::Int
    kstar2::Int
    mean_motion_nb::Float64
    separation_nb::Float64
    step_ks_nb::Float64
    step_cm_nb::Float64
    name_cm::Int
    specific_energy_cm_nb::Float64
    kstar_cm::Int
    age1_myr::Float64
    age2_myr::Float64
    epoch1_myr::Float64
    epoch2_myr::Float64
end

"""
    KSPairListing

One `bdat.9_<t>` file (`bindat.f`, written at every OUTPUT for KZ(8) ≥ 2),
from [`read_ks_pair_listing`](@ref): the three header lines and the KS pairs.

| Header | Field | Unit |
|---|---|---|
| line 1, NPAIRS, MODEL, NRUN | `n_pairs`, `model`, `n_run` | — |
| line 1, N, NC (core members), NMERGE | `n_total`, `n_core`, `n_merge` | — |
| line 1, TTOT (F7.1) | `time_nb` | NB |
| line 1, RSCALE, RTIDE, RC | `rscale_nb`, `rtide_nb`, `core_radius_nb` | NB |
| line 1, TPHYS (F8.1) | `time_myr` | Myr |
| line 1, ECRIT = −1.5(TIDAL(1)·M²)^⅓ | `ecrit_nb` (−0 without tidal field) | NB |
| line 2, E(1:10) | `energy_nb` | NB |
| line 3, SBCOLL, BBCOLL | `sbcoll_nb`, `bbcoll_nb` | NB |
| line 3, ZKIN, POT | `kinetic_energy_nb`, `potential_energy_nb` | NB |
| line 3, EBIN0, EBIN, ESUB, EMERGE | `ebin0_nb`, `ebin_nb`, `esub_nb`, `emerge_nb` | NB |
| line 3, BE(3) | `total_energy_nb` | NB |
| line 3, ZMASS, ZMBIN (Σ c.m. masses) | `total_mass_nb`, `binary_mass_nb` | NB |
| line 3, CHCOLL, ECOLL | `chcoll_nb`, `ecoll_nb` | NB |

`overflow` is `true` when a header field was printed as asterisks (NPAIRS,
MODEL, NRUN, NC, NMERGE ≥ 10⁴, N ≥ 10⁶); such a counter reads 0, and the
number of pairs is then `length(records)`. The header time has one decimal:
the exact output time is the `<t>` of the file name. `records` holds one
[`KSPairRecord`](@ref) per KS pair.
"""
Base.@kwdef struct KSPairListing
    n_pairs::Int = 0
    model::Int = 0
    n_run::Int = 0
    n_total::Int = 0
    n_core::Int = 0
    n_merge::Int = 0
    time_nb::Float64 = NaN
    rscale_nb::Float64 = NaN
    rtide_nb::Float64 = NaN
    core_radius_nb::Float64 = NaN
    time_myr::Float64 = NaN
    ecrit_nb::Float64 = NaN
    energy_nb::NTuple{10,Float64} = ntuple(_ -> NaN, Val(10))
    sbcoll_nb::Float64 = NaN
    bbcoll_nb::Float64 = NaN
    kinetic_energy_nb::Float64 = NaN
    potential_energy_nb::Float64 = NaN
    ebin0_nb::Float64 = NaN
    ebin_nb::Float64 = NaN
    esub_nb::Float64 = NaN
    emerge_nb::Float64 = NaN
    total_energy_nb::Float64 = NaN
    total_mass_nb::Float64 = NaN
    binary_mass_nb::Float64 = NaN
    chcoll_nb::Float64 = NaN
    ecoll_nb::Float64 = NaN
    overflow::Bool = false
    records::Vector{KSPairRecord} = KSPairRecord[]
end

"""Build a [`KSPairRecord`](@ref) from one list-directed line of `bdat.9_<t>`."""
function _parse_ks_pair_record(line::AbstractString)::KSPairRecord
    t = _binfile_tokens(line, _BDAT_RECORD_VALUES, "bdat.9")
    r(k) = _parse_fortran_real(t[k])
    i(k) = parse(Int, t[k])
    period = r(7)
    return KSPairRecord(;
        name1 = i(1),
        name2 = i(2),
        mass1_msun = r(3),
        mass2_msun = r(4),
        binding_energy_nb = r(5),
        eccentricity = r(6),
        # PB = 0 is the engine's mark of a < 0.
        period_d = period == 0 ? NaN : period,
        semi_major_axis_au = r(8),
        distance_density_centre_pc = r(9),
        speed_cm_kms = r(10),
        kstar1 = i(11),
        kstar2 = i(12),
        mean_motion_nb = r(13),
        separation_nb = r(14),
        step_ks_nb = r(15),
        step_cm_nb = r(16),
        name_cm = i(17),
        specific_energy_cm_nb = r(18),
        kstar_cm = i(19),
        age1_myr = r(20),
        age2_myr = r(21),
        epoch1_myr = r(22),
        epoch2_myr = r(23),
    )
end

"""
    read_ks_pair_listing(path::AbstractString) -> KSPairListing

Read one `bdat.9_<t>` file, the KS-pair data bank written by `bindat.f` at
every OUTPUT (KZ(8) ≥ 2), into a [`KSPairListing`](@ref): the fixed-width
header lines (FORMAT 30, 35, 40), the label line, and one list-directed
[`KSPairRecord`](@ref) per KS pair. A header block repeated inside the file
(up to three non-record lines before a further label line) is skipped; only
the first is parsed. Records without 23 values or with unparseable values are
skipped and counted in one warning. A missing file raises an error, a
malformed header line an `ArgumentError`; an empty file gives an empty
listing with absent header values. The record hazards (hierarchy rows that mix
the inner and outer orbit, ghost rows, P = 0 for a < 0 stored as `NaN`) are
listed under [`KSPairRecord`](@ref).
"""
function read_ks_pair_listing(path::AbstractString)::KSPairListing
    isfile(path) || error("KS-pair listing (bdat.9) not found: $path")
    lines = filter(l -> !isempty(strip(l)), readlines(path))
    is_label(l) = startswith(lstrip(l), "NAME(I1)")
    is_record(l) = length(split(l)) == _BDAT_RECORD_VALUES

    # Header lines: up to three non-record lines directly above a label line.
    # A file cut before its label line keeps its leading lines as the header.
    header_rows = falses(length(lines))
    first_block = Int[]
    labels = findall(is_label, lines)
    starts = isempty(labels) ? [min(3, length(lines)) + 1] : labels
    for (b, k) in enumerate(starts)
        j = k - 1
        while j >= max(1, k - 3) && !is_label(lines[j]) && !is_record(lines[j])
            header_rows[j] = true
            b == 1 && pushfirst!(first_block, j)
            j -= 1
        end
    end

    formats = (_BDAT_HEADER1_FORMAT, _BDAT_HEADER2_FORMAT, _BDAT_HEADER3_FORMAT)
    header = map(1:3) do k
        k <= length(first_block) ||
            return FortranRecord(fill(NaN, _binfile_nvalues(formats[k])), false)
        try
            _binfile_fixed(lines[first_block[k]], formats[k])
        catch e
            e isa ArgumentError || rethrow()
            throw(ArgumentError("bdat.9 header line $k of $path: $(e.msg)"))
        end
    end
    v1, v2, v3 = header[1].values, header[2].values, header[3].values

    records = KSPairRecord[]
    n_data = 0
    n_skipped = 0
    for (j, line) in enumerate(lines)
        (header_rows[j] || is_label(line)) && continue
        n_data += 1
        try
            push!(records, _parse_ks_pair_record(line))
        catch e
            e isa ArgumentError || rethrow()
            n_skipped += 1
            @debug "Skipping malformed bdat.9 record" line exception = e
        end
    end
    n_skipped > 0 && @warn "bdat.9: $n_skipped of $n_data records skipped (malformed)" path

    return KSPairListing(;
        n_pairs = _fixed_int(v1[1]),
        model = _fixed_int(v1[2]),
        n_run = _fixed_int(v1[3]),
        n_total = _fixed_int(v1[4]),
        n_core = _fixed_int(v1[5]),
        n_merge = _fixed_int(v1[6]),
        time_nb = v1[7],
        rscale_nb = v1[8],
        rtide_nb = v1[9],
        core_radius_nb = v1[10],
        time_myr = v1[11],
        ecrit_nb = v1[12],
        energy_nb = ntuple(k -> v2[k], Val(10)),
        sbcoll_nb = v3[1],
        bbcoll_nb = v3[2],
        kinetic_energy_nb = v3[3],
        potential_energy_nb = v3[4],
        ebin0_nb = v3[5],
        ebin_nb = v3[6],
        esub_nb = v3[7],
        emerge_nb = v3[8],
        total_energy_nb = v3[9],
        total_mass_nb = v3[10],
        binary_mass_nb = v3[11],
        chcoll_nb = v3[12],
        ecoll_nb = v3[13],
        overflow = any(h -> h.overflow, header),
        records = records,
    )
end

"""
    read_all_ks_pair_listings(dir::AbstractString) -> Vector{KSPairListing}

Read every `bdat.9_<t>` file in `dir` with [`read_ks_pair_listing`](@ref),
in the order of the time `<t>` of the file name (the left-justified TTOT the
engine appends). Files whose suffix is not a number are ignored; a missing
directory gives an empty vector.
"""
function read_all_ks_pair_listings(dir::AbstractString)::Vector{KSPairListing}
    return KSPairListing[read_ks_pair_listing(p) for p in _binfile_time_files(dir, "bdat.9_")]
end

# ---------------------------------------------------------------------------
# bwdat.19_<t>: wide pairs
# ---------------------------------------------------------------------------

"""Values per record of `bwdat.19_<t>` (`bindat.f`)."""
const _BWDAT_RECORD_VALUES = 14

"""FORMAT 55 of `bindat.f`: `(' WIDE PAIRS    T TPHYS N ',1P,E27.16,E27.16,0P,I12)`."""
const _BWDAT_HEADER_FORMAT =
    vcat(fmt_lit(" WIDE PAIRS    T TPHYS N "), fmt_real(1, 27), fmt_real(1, 27), fmt_int(1, 12))

"""
    WidePairRecord

One wide pair of a `bwdat.19_<t>` listing, written list-directed by
`bindat.f` (14 values), from [`read_wide_pair_listing`](@ref). Every body
I ≥ IFIRST (single stars and KS c.m. bodies) is paired with its nearest
neighbour inside 0.1·RSCALE and kept when 0 < a < 0.1·RSCALE; there is no
perturbation test.

| # | Label | Field | Unit |
|---|---|---|---|
| 1–2 | NAME(I1), NAME(I2) | `name1`, `name2` (as printed) | — |
| 3–4 | M(I1)[M*], M(I2)[M*] | `mass1_msun`, `mass2_msun` | M☉ |
| 5 | EB[NB] | `binding_energy_nb` (−½m₁m₂/a) | NB |
| 6 | ECC | `eccentricity` | — |
| 7 | P[Days] | `period_d` | d |
| 8 | SEMI[AU] | `semi_major_axis_au` | AU |
| 9 | STEP(JMIN)[NB] | `step1_nb` = STEP(I), step of member 1 | NB |
| 10 | STEP(ICM)[NB] | `step2_nb` = STEP(JMIN), step of member 2 | NB |
| 11 | RI[PC] | `distance_density_centre_pc` (pair c.m. from RDENS) | pc |
| 12 | VI[km/s] | `speed_cm_kms` (pair c.m. speed, absolute) | km s⁻¹ |
| 13–14 | K*(I1), K*(I2) | `kstar1`, `kstar2` | — |

Corrections: the two step columns are labelled one position off (the first
holds STEP(I), the second STEP(JMIN)) and are stored under their true meaning.

Hazard: a member can be a KS c.m. body, whose mass is the binary mass and
whose stellar type is the c.m. type, so the pair may be the outer orbit of a
triple. The engine names a KS c.m. NZERO + NAME(J1) (> NZERO) and a merged
hierarchy c.m. with a negative name. `member1_is_cm`, `member2_is_cm` are
`true` for a negative name and, when the reader was given `nzero > 0`, for a
name above `nzero`; without `nzero` a positive c.m. name is not recognised
(the file does not carry NZERO: N in the header is the current membership).
"""
Base.@kwdef struct WidePairRecord
    name1::Int
    name2::Int
    mass1_msun::Float64
    mass2_msun::Float64
    binding_energy_nb::Float64
    eccentricity::Float64
    period_d::Float64
    semi_major_axis_au::Float64
    step1_nb::Float64
    step2_nb::Float64
    distance_density_centre_pc::Float64
    speed_cm_kms::Float64
    kstar1::Int
    kstar2::Int
    member1_is_cm::Bool
    member2_is_cm::Bool
end

"""
    WidePairListing

One `bwdat.19_<t>` file (`bindat.f`, KZ(8) ≥ 2), from
[`read_wide_pair_listing`](@ref). Header (FORMAT 55): `time_nb` = TTOT [NB],
`time_myr` = TTOT·T* [Myr], `n_total` = N; `overflow` is `true` when N was
printed as asterisks. `nzero` is the initial membership passed to the reader
(0 when not given) that sets the c.m. flags of the `records`
([`WidePairRecord`](@ref)).
"""
Base.@kwdef struct WidePairListing
    time_nb::Float64 = NaN
    time_myr::Float64 = NaN
    n_total::Int = 0
    nzero::Int = 0
    overflow::Bool = false
    records::Vector{WidePairRecord} = WidePairRecord[]
end

"""Whether a member `name` of a wide pair is a c.m. body (negative, or above a known NZERO)."""
_wide_pair_member_is_cm(name::Int, nzero::Int)::Bool = name < 0 || (nzero > 0 && name > nzero)

"""Build a [`WidePairRecord`](@ref) from one list-directed line of `bwdat.19_<t>`."""
function _parse_wide_pair_record(line::AbstractString, nzero::Int)::WidePairRecord
    t = _binfile_tokens(line, _BWDAT_RECORD_VALUES, "bwdat.19")
    r(k) = _parse_fortran_real(t[k])
    i(k) = parse(Int, t[k])
    name1 = i(1)
    name2 = i(2)
    return WidePairRecord(;
        name1 = name1,
        name2 = name2,
        mass1_msun = r(3),
        mass2_msun = r(4),
        binding_energy_nb = r(5),
        eccentricity = r(6),
        period_d = r(7),
        semi_major_axis_au = r(8),
        step1_nb = r(9),
        step2_nb = r(10),
        distance_density_centre_pc = r(11),
        speed_cm_kms = r(12),
        kstar1 = i(13),
        kstar2 = i(14),
        member1_is_cm = _wide_pair_member_is_cm(name1, nzero),
        member2_is_cm = _wide_pair_member_is_cm(name2, nzero),
    )
end

"""
    read_wide_pair_listing(path::AbstractString; nzero::Integer = 0) -> WidePairListing

Read one `bwdat.19_<t>` file, the wide pairs written by `bindat.f` at every
OUTPUT (KZ(8) ≥ 2), into a [`WidePairListing`](@ref): the ` WIDE PAIRS`
header (FORMAT 55), the label line, and one [`WidePairRecord`](@ref) per
pair. A repeated header or label line is skipped. `nzero`, the initial number
of bodies, lets the reader flag members whose name exceeds it as KS c.m.
bodies. The shifted step labels are corrected (see `WidePairRecord`).
Records without 14 values or with unparseable values are skipped and counted
in one warning; a missing file raises an error, a malformed header an
`ArgumentError`; an empty file gives an empty listing.
"""
function read_wide_pair_listing(path::AbstractString; nzero::Integer = 0)::WidePairListing
    isfile(path) || error("Wide-pair listing (bwdat.19) not found: $path")
    nz = Int(nzero)
    header = nothing
    records = WidePairRecord[]
    n_data = 0
    n_skipped = 0
    for line in eachline(path)
        s = lstrip(line)
        isempty(rstrip(s)) && continue
        if startswith(s, "WIDE PAIRS")
            if header === nothing
                header = try
                    _binfile_fixed(line, _BWDAT_HEADER_FORMAT)
                catch e
                    e isa ArgumentError || rethrow()
                    throw(ArgumentError("bwdat.19 header of $path: $(e.msg)"))
                end
            end
            continue
        end
        startswith(s, "NAME(I1)") && continue
        n_data += 1
        try
            push!(records, _parse_wide_pair_record(line, nz))
        catch e
            e isa ArgumentError || rethrow()
            n_skipped += 1
            @debug "Skipping malformed bwdat.19 record" line exception = e
        end
    end
    n_skipped > 0 && @warn "bwdat.19: $n_skipped of $n_data records skipped (malformed)" path

    header === nothing && return WidePairListing(; nzero = nz, records = records)
    v = header.values
    return WidePairListing(;
        time_nb = v[1],
        time_myr = v[2],
        n_total = _fixed_int(v[3]),
        nzero = nz,
        overflow = header.overflow,
        records = records,
    )
end

"""
    read_all_wide_pair_listings(dir::AbstractString; nzero::Integer = 0)
        -> Vector{WidePairListing}

Read every `bwdat.19_<t>` file in `dir` with [`read_wide_pair_listing`](@ref)
(`nzero` passed through), in the order of the time `<t>` of the file name.
Files whose suffix is not a number are ignored; a missing directory gives an
empty vector.
"""
function read_all_wide_pair_listings(
    dir::AbstractString;
    nzero::Integer = 0,
)::Vector{WidePairListing}
    return WidePairListing[
        read_wide_pair_listing(p; nzero = nzero) for p in _binfile_time_files(dir, "bwdat.19_")
    ]
end

# ---------------------------------------------------------------------------
# merger.84: hierarchies at each OUTPUT
# ---------------------------------------------------------------------------

"""
FORMAT 3 of `bindat.f`: three label literals, then
`1P,E12.5,0P,2I12,3I4,2F7.3,2F12.5,1P,3E14.5,0P`.
"""
const _MERGER84_FORMAT = vcat(
    fmt_lit("BINDAT: Time[NB] NAME[I1] NAME[I3] K*[I1]"),
    fmt_lit("K*[I3] K*[IM] ECC0 ECC1 PERI(I3)/PCR "),
    fmt_lit("PERI(INCM)[RSM] P0[days] P1[days] SEMI1[NB] "),
    fmt_real(1, 12),
    fmt_int(2, 12),
    fmt_int(3, 4),
    fmt_real(2, 7),
    fmt_real(2, 12),
    fmt_real(3, 14),
)

"""
    HierarchyListingRecord

One hierarchy of `merger.84` at one OUTPUT (`bindat.f`, FORMAT 3, label text
in every record), from [`read_hierarchy_listing`](@ref). Written for every
merged hierarchy with BODY(ICM) > 0 at every BINDAT call (KZ(8) ≥ 2).

| # | Label | Field | Unit |
|---|---|---|---|
| 1 | Time[NB] | `time_nb` | NB |
| 2 | NAME[I1] | `name_primary` (inner c.m. body = inner primary name) | — |
| 3 | NAME[I3] | `name_ghost` (ghost = inner secondary) | — |
| 4 | K*[I1] | `kstar_primary` (type of the inner c.m. body) | — |
| 5 | K*[I3] | `kstar_ghost` | — |
| 6 | K*[IM] | `kstar_inner_cm` (KSTARM, inner c.m. type saved at the merger) | — |
| 7 | ECC0 | `eccentricity_inner` (from ‖XREL‖, TDOT2 ignored) | — |
| 8 | ECC1 | `eccentricity_outer` | — |
| 9 | PERI(I3)/PCR | `pericentre_outer_over_pcrit` (Mardling–Aarseth PCR) | — |
| 10 | PERI(INCM)[RSM] | `pericentre_inner_over_radius` (/ max stellar radius) | — |
| 11 | P0[days] | `period_inner_d` | d |
| 12 | P1[days] | `period_outer_d` | d |
| 13 | SEMI1[NB] | `semi_major_axis_outer_nb` | NB |

Hazards: `pericentre_inner_over_radius` (F12.5) overflows for ratios above
10⁶ (compact inner members): `NaN` and `overflow = true`. The outer body's
name is not written: join to the `bdat.9_<t>` row with `name_cm < 0` of the
same time.
"""
Base.@kwdef struct HierarchyListingRecord
    time_nb::Float64
    name_primary::Int
    name_ghost::Int
    kstar_primary::Int
    kstar_ghost::Int
    kstar_inner_cm::Int
    eccentricity_inner::Float64
    eccentricity_outer::Float64
    pericentre_outer_over_pcrit::Float64
    pericentre_inner_over_radius::Float64
    period_inner_d::Float64
    period_outer_d::Float64
    semi_major_axis_outer_nb::Float64
    overflow::Bool = false
end

"""Build a [`HierarchyListingRecord`](@ref) from one line of `merger.84`."""
function _parse_hierarchy_listing_record(line::AbstractString)::HierarchyListingRecord
    rec = _binfile_fixed(line, _MERGER84_FORMAT)
    v = rec.values
    return HierarchyListingRecord(;
        time_nb = v[1],
        name_primary = _fixed_int(v[2]),
        name_ghost = _fixed_int(v[3]),
        kstar_primary = _fixed_int(v[4], KSTAR_ABSENT),
        kstar_ghost = _fixed_int(v[5], KSTAR_ABSENT),
        kstar_inner_cm = _fixed_int(v[6], KSTAR_ABSENT),
        eccentricity_inner = v[7],
        eccentricity_outer = v[8],
        pericentre_outer_over_pcrit = v[9],
        pericentre_inner_over_radius = v[10],
        period_inner_d = v[11],
        period_outer_d = v[12],
        semi_major_axis_outer_nb = v[13],
        overflow = rec.overflow,
    )
end

"""
    read_hierarchy_listing(path::AbstractString) -> Vector{HierarchyListingRecord}

Read `merger.84`, one [`HierarchyListingRecord`](@ref) per merged hierarchy
per OUTPUT, written by `bindat.f` (KZ(8) ≥ 2, appended on restart). The file
has no header: the label text is part of every fixed-width record. Malformed
lines are skipped and counted in one warning; a missing file raises an error.
"""
function read_hierarchy_listing(path::AbstractString)::Vector{HierarchyListingRecord}
    return _binfile_read(
        _parse_hierarchy_listing_record,
        HierarchyListingRecord,
        path,
        "merger.84",
        _ -> false,
    )
end

# ---------------------------------------------------------------------------
# coll.13: collisions
# ---------------------------------------------------------------------------

"""
    CollisionRecord

One collision of `coll.13`, written list-directed by `mix.f` (16 values),
from [`read_collision_log`](@ref).

| # | Label | Field | Unit |
|---|---|---|---|
| 1 | TIME[NB] | `time_nb` | NB |
| 2–3 | NAME(I1), NAME(I2) | `name1`, `name2` | — |
| 4–5 | K*(I1), K*(I2) | `kstar1`, `kstar2` (before) | — |
| 6 | K*(INEW) | `kstar_product` | — |
| 7–8 | M(I1)[M*], M(I2)[M*] | `mass1_msun`, `mass2_msun` (before) | M☉ |
| 9 | M(INEW)[M*] | `mass_product_msun` | M☉ |
| 10 | DM[M*] | `mass_lost_msun` | M☉ |
| 11–12 | RS(I1)[R*], RS(I2)[R*] | `radius1_rsun`, `radius2_rsun` (before) | R☉ |
| 13 | RI[RC] | `distance_over_core_radius` (I1 from RDENS) | — |
| 14 | R12[R*] | `separation_rsun` | R☉ |
| 15 | ECC | `eccentricity` (1 − R/a, exact at pericentre only; < 0 occurs) | — |
| 16 | P[days] | `period_d` (negative for a < 0) | d |

Hazards: I1 is the more evolved star (higher KSTAR; the heavier for
KSTAR ≤ 1), so the names are reordered with respect to the KS pair. There is
no IQCOLL column; KS collisions found in KSINT/UNPERT reach MIX with
IQCOLL = −2 whether the cause is stellar contact or a GR coalescence, so a
record does not identify a gravitational-wave merger.
"""
Base.@kwdef struct CollisionRecord
    time_nb::Float64
    name1::Int
    name2::Int
    kstar1::Int
    kstar2::Int
    kstar_product::Int
    mass1_msun::Float64
    mass2_msun::Float64
    mass_product_msun::Float64
    mass_lost_msun::Float64
    radius1_rsun::Float64
    radius2_rsun::Float64
    distance_over_core_radius::Float64
    separation_rsun::Float64
    eccentricity::Float64
    period_d::Float64
end

"""Build a [`CollisionRecord`](@ref) from one list-directed line of `coll.13`."""
function _parse_collision_record(line::AbstractString)::CollisionRecord
    t = _binfile_tokens(line, 16, "coll.13")
    r(k) = _parse_fortran_real(t[k])
    i(k) = parse(Int, t[k])
    return CollisionRecord(;
        time_nb = r(1),
        name1 = i(2),
        name2 = i(3),
        kstar1 = i(4),
        kstar2 = i(5),
        kstar_product = i(6),
        mass1_msun = r(7),
        mass2_msun = r(8),
        mass_product_msun = r(9),
        mass_lost_msun = r(10),
        radius1_rsun = r(11),
        radius2_rsun = r(12),
        distance_over_core_radius = r(13),
        separation_rsun = r(14),
        eccentricity = r(15),
        period_d = r(16),
    )
end

"""Header lines of `coll.13` and `cmbody.26`: the MODEL line and the label line."""
_collision_header_line(s::AbstractString)::Bool =
    startswith(s, "MODEL:") || startswith(s, "TIME[NB]")

"""
    read_collision_log(path::AbstractString) -> Vector{CollisionRecord}

Read `coll.13`, the collisions written by `mix.f` (KZ(19) ≥ 3, IQCOLL ≠ 3),
one [`CollisionRecord`](@ref) per collision. The header block (MODEL line and
16-label line, written on the first collision of a run started fresh) is
skipped wherever it occurs. IQCOLL = −2, the value with which every KS
collision found in KSINT/UNPERT reaches MIX, marks stellar contact and GR
coalescence alike; the file does not record it. Records without 16 values
or with unparseable values are skipped and counted in one warning; a missing
file raises an error.
"""
function read_collision_log(path::AbstractString)::Vector{CollisionRecord}
    return _binfile_read(
        _parse_collision_record,
        CollisionRecord,
        path,
        "coll.13",
        _collision_header_line,
    )
end

# ---------------------------------------------------------------------------
# coal.24: coalescences
# ---------------------------------------------------------------------------

"""FORMAT 56 of `coal.f` after its leading `A8`: `(I3,1P,E17.10,2I10,3I4,18E13.5)`."""
const _COAL_FORMAT =
    vcat(fmt_int(1, 3), fmt_real(1, 17), fmt_int(2, 10), fmt_int(3, 4), fmt_real(18, 13))

"""Channels of `coal.24`: the `WHICH1` text of `coal.f` (A8), stripped."""
const _COAL_CHANNELS =
    Dict("BINARY" => :binary, "ROCHE" => :roche, "HYPERB" => :hyperbolic, "CHAIN" => :chain)

"""
    CoalescenceRecord

One coalescence of `coal.24`, written by `coal.f` (FORMAT 56, fixed width,
26 values) when called from EXPEL, EXPEL2 or ROCHE, from
[`read_coalescence_log`](@ref).

| # | Label | Field | Unit |
|---|---|---|---|
| 1 | WHICH1 | `channel`: `:binary`, `:roche`, `:hyperbolic`, `:chain` | — |
| 2 | IQCOLL | `iqcoll` (0 binary, 3 Roche, −1 hyperbolic, 5 chain) | — |
| 3 | TTOT | `time_nb` | NB |
| 4–5 | NAME(I1,I2) | `name1`, `name2` (before the BODY0 swap) | — |
| 6–7 | K*(I1,I2,I) | `kstar1`, `kstar2` (before) | — |
| 8 | K*(I1,I2,I) | `kstar_product` | — |
| 9–10 | M(I1)[NB], M(I2)[NB] | `mass1_bse_msun`, `mass2_bse_msun` | M☉ |
| 11 | R12[NB] | `separation_nb` | NB |
| 12 | ECC | `eccentricity` (1 − RCOLL/a; negative values occur) | — |
| 13 | SEMI | `semi_major_axis_nb` | NB |
| 14 | EB | `binding_energy_nb` (μH before KSTERM) | NB |
| 15 | DP | `tidal_potential_correction_nb` (POT2 − POT1) | NB |
| 16 | P[days] | `period_d` | d |
| 17–18 | M1[M*], M2[M*] | `mass1_msun`, `mass2_msun` (BODY·M*, before) | M☉ |
| 19 | MTOT[M*] | `mass_product_msun` | M☉ |
| 20 | DM[M*] | `mass_lost_msun` | M☉ |
| 21 | RAD1[R*] | `radius_survivor_rsun` (survivor slot, possibly swapped from I2) | R☉ |
| 22 | RAD2[R*] | not stored: always 0 (zeroed before the write) | — |
| 23 | RCOLL[R*] | `collision_distance_rsun` (KS) or DMINC (chain) | R☉ |
| 24 | VINF[km/s] | `velocity_infinity_kms` (hyperbolic, else 0) | km s⁻¹ |
| 25 | RI[NB] | `distance_density_centre_nb` (‖X(I1) − RDENS‖) | NB |
| 26 | VI[NB] | `speed_nb` (‖XDOT(I1)‖) | NB |

Corrections: columns 9–10 are labelled `[NB]` but hold the masses passed by
the caller in M☉ (BSE masses after common envelope or Roche overflow); they
are stored as `_msun`. `overflow` is `true` when a field was printed as
asterisks (read as absent).
"""
Base.@kwdef struct CoalescenceRecord
    channel::Symbol
    iqcoll::Int
    time_nb::Float64
    name1::Int
    name2::Int
    kstar1::Int
    kstar2::Int
    kstar_product::Int
    mass1_bse_msun::Float64
    mass2_bse_msun::Float64
    separation_nb::Float64
    eccentricity::Float64
    semi_major_axis_nb::Float64
    binding_energy_nb::Float64
    tidal_potential_correction_nb::Float64
    period_d::Float64
    mass1_msun::Float64
    mass2_msun::Float64
    mass_product_msun::Float64
    mass_lost_msun::Float64
    radius_survivor_rsun::Float64
    collision_distance_rsun::Float64
    velocity_infinity_kms::Float64
    distance_density_centre_nb::Float64
    speed_nb::Float64
    overflow::Bool = false
end

"""Build a [`CoalescenceRecord`](@ref) from one line of `coal.24`."""
function _parse_coalescence_record(line::AbstractString)::CoalescenceRecord
    rec = _binfile_fixed(line, _COAL_FORMAT; start = 9)
    which = strip(first(line, 8))
    channel = get(_COAL_CHANNELS, which, :none)
    channel === :none && throw(ArgumentError("coal.24 channel $(repr(which)) unknown"))
    v = rec.values
    return CoalescenceRecord(;
        channel = channel,
        iqcoll = _fixed_int(v[1]),
        time_nb = v[2],
        name1 = _fixed_int(v[3]),
        name2 = _fixed_int(v[4]),
        kstar1 = _fixed_int(v[5], KSTAR_ABSENT),
        kstar2 = _fixed_int(v[6], KSTAR_ABSENT),
        kstar_product = _fixed_int(v[7], KSTAR_ABSENT),
        mass1_bse_msun = v[8],
        mass2_bse_msun = v[9],
        separation_nb = v[10],
        eccentricity = v[11],
        semi_major_axis_nb = v[12],
        binding_energy_nb = v[13],
        tidal_potential_correction_nb = v[14],
        period_d = v[15],
        mass1_msun = v[16],
        mass2_msun = v[17],
        mass_product_msun = v[18],
        mass_lost_msun = v[19],
        radius_survivor_rsun = v[20],
        # v[21] = RADIUS(I2)·SU, zeroed by the engine before the write.
        collision_distance_rsun = v[22],
        velocity_infinity_kms = v[23],
        distance_density_centre_nb = v[24],
        speed_nb = v[25],
        overflow = rec.overflow,
    )
end

"""Header lines of `coal.24`: the MODEL line and the label line."""
_coalescence_header_line(s::AbstractString)::Bool =
    startswith(s, "MODEL:") || startswith(s, "WHICH1")

"""
    read_coalescence_log(path::AbstractString) -> Vector{CoalescenceRecord}

Read `coal.24`, the coalescences written by `coal.f` (KZ(19) ≥ 3), one
fixed-width [`CoalescenceRecord`](@ref) per event. The header block (MODEL
line and label line, written with the first coalescence of the run) is
skipped wherever it occurs. The `M(I1)[NB]`, `M(I2)[NB]` columns are stored
as M☉ (their true unit) and the always-zero `RAD2` column is dropped.
Malformed records are skipped and counted in one warning; a missing file
raises an error.
"""
function read_coalescence_log(path::AbstractString)::Vector{CoalescenceRecord}
    return _binfile_read(
        _parse_coalescence_record,
        CoalescenceRecord,
        path,
        "coal.24",
        _coalescence_header_line,
    )
end

# ---------------------------------------------------------------------------
# cmbody.26: binary coalescence via CMBODY
# ---------------------------------------------------------------------------

"""
    CMBodyRecord

One coalescence of `cmbody.26`, written list-directed by `cmbody.f` (16
values) for IQCOLL = 3 or KSTAR(ICM) ≥ 10 before KSTERM (circularised,
Roche or GR-flagged pairs), from [`read_cmbody_log`](@ref).

| # | Label | Field | Unit |
|---|---|---|---|
| 1 | TIME[NB] | `time_nb` | NB |
| 2–3 | NAME(I1), NAME(I2) | `name1`, `name2` (after the swap to the heavier body) | — |
| 4–5 | K*(I1), K*(I2) | `kstar1`, `kstar2` (after MIX) | — |
| 6 | IQCOLL | `iqcoll` | — |
| 7–8 | M(I1)[M*], M(I2)[M*] | `mass1_msun`, `mass2_msun` (slot masses) | M☉ |
| 9 | DM[M*] | `mass_lost_msun` | M☉ |
| 10–11 | RS(I1)[R*], RS(I2)[R*] | `radius1_rsun`, `radius2_rsun` (after MIX) | R☉ |
| 12 | RI/RC | `distance_over_core_radius` (ghost I2 from RDENS) | — |
| 13 | R12[R*] | `collision_distance_rsun` | R☉ |
| 14 | ECC | `eccentricity` (1 − R/a) | — |
| 15 | SEMI[R*] | `semi_major_axis_rsun` | R☉ |
| 16 | P[days] | `period_d` (NaN for a < 0) | d |

Hazards: the masses are those of the slots, so they pair with swapped names
when BODY(I2) > BODY(I1); which of `kstar1`, `kstar2` holds the product type
is not settled (MIX's own I1 is the more evolved star) [UNSURE in the
dossier]. No sample run contains a record; the layout is transcribed from
the source.
"""
Base.@kwdef struct CMBodyRecord
    time_nb::Float64
    name1::Int
    name2::Int
    kstar1::Int
    kstar2::Int
    iqcoll::Int
    mass1_msun::Float64
    mass2_msun::Float64
    mass_lost_msun::Float64
    radius1_rsun::Float64
    radius2_rsun::Float64
    distance_over_core_radius::Float64
    collision_distance_rsun::Float64
    eccentricity::Float64
    semi_major_axis_rsun::Float64
    period_d::Float64
end

"""Build a [`CMBodyRecord`](@ref) from one list-directed line of `cmbody.26`."""
function _parse_cmbody_record(line::AbstractString)::CMBodyRecord
    t = _binfile_tokens(line, 16, "cmbody.26")
    r(k) = _parse_fortran_real(t[k])
    i(k) = parse(Int, t[k])
    return CMBodyRecord(;
        time_nb = r(1),
        name1 = i(2),
        name2 = i(3),
        kstar1 = i(4),
        kstar2 = i(5),
        iqcoll = i(6),
        mass1_msun = r(7),
        mass2_msun = r(8),
        mass_lost_msun = r(9),
        radius1_rsun = r(10),
        radius2_rsun = r(11),
        distance_over_core_radius = r(12),
        collision_distance_rsun = r(13),
        eccentricity = r(14),
        semi_major_axis_rsun = r(15),
        period_d = r(16),
    )
end

"""
    read_cmbody_log(path::AbstractString) -> Vector{CMBodyRecord}

Read `cmbody.26`, the binary coalescences written by `cmbody.f` together
with the stdout ` BINARY COAL    IQCOLL =` line, one [`CMBodyRecord`](@ref)
per event. The header (MODEL line and label line, written only by a process
at TTOT = 0, practically never) is skipped wherever it occurs. Records
without 16 values or with unparseable values are skipped and counted in one
warning; a missing file raises an error.
"""
function read_cmbody_log(path::AbstractString)::Vector{CMBodyRecord}
    return _binfile_read(
        _parse_cmbody_record,
        CMBodyRecord,
        path,
        "cmbody.26",
        _collision_header_line,
    )
end

# ---------------------------------------------------------------------------
# binev.17: binary-evolution stages
# ---------------------------------------------------------------------------

"""
    BinaryStageRecord

One binary-evolution stage of `binev.17`, written list-directed by `binev.f`
(15 values), from [`read_binary_stage_log`](@ref). Callers: COAL, CMBODY
(before a collision), CHAOS, SPIRAL, CHRECT, MDOT (type change), ROCHE (new
Roche stage).

| # | Label | Field | Unit |
|---|---|---|---|
| 1 | Time[Myr] | `time_myr` (TPHYS) | Myr |
| 2–3 | NAME(I1), NAME(I2) | `name1`, `name2` | — |
| 4–6 | K*(I1), K*(I2), K*(ICM) | `kstar1`, `kstar2`, `kstar_cm` | — |
| 7–8 | M(I1)[M*], M(I2)[M*] | `mass1_msun`, `mass2_msun` | M☉ |
| 9–10 | RS(I1)[R*], RS(I2)[R*] | `radius1_rsun`, `radius2_rsun` | R☉ |
| 11 | RI[PC] | `distance_density_centre_pc` (‖X(ICM) − RDENS‖) | pc |
| 12 | ECC | `eccentricity` | — |
| 13 | SEMI[R*] | `semi_major_axis_rsun` | R☉ |
| 14 | P[days] | `period_d` | d |
| 15 | IQCOLL | `iqcoll` (global indicator at the time of the call) | — |

Hazards: TPHYS is the COMMON value last set by MDOT, ROCHE, MIX or DEGEN and
may lag TTOT·T* [UNSURE in the dossier]; for a merged or ghost pair the
eccentricity is that of the inner binary, from UM with TDOT2 ignored;
`iqcoll` is not set by BINEV (−2 after a KS collision, else 0). A stage is
not written when K*(I1) + K*(I2) + K*(ICM) equals the sum of either of the
two previous calls, of any pair.
"""
Base.@kwdef struct BinaryStageRecord
    time_myr::Float64
    name1::Int
    name2::Int
    kstar1::Int
    kstar2::Int
    kstar_cm::Int
    mass1_msun::Float64
    mass2_msun::Float64
    radius1_rsun::Float64
    radius2_rsun::Float64
    distance_density_centre_pc::Float64
    eccentricity::Float64
    semi_major_axis_rsun::Float64
    period_d::Float64
    iqcoll::Int
end

"""Build a [`BinaryStageRecord`](@ref) from one list-directed line of `binev.17`."""
function _parse_binary_stage_record(line::AbstractString)::BinaryStageRecord
    t = _binfile_tokens(line, 15, "binev.17")
    r(k) = _parse_fortran_real(t[k])
    i(k) = parse(Int, t[k])
    return BinaryStageRecord(;
        time_myr = r(1),
        name1 = i(2),
        name2 = i(3),
        kstar1 = i(4),
        kstar2 = i(5),
        kstar_cm = i(6),
        mass1_msun = r(7),
        mass2_msun = r(8),
        radius1_rsun = r(9),
        radius2_rsun = r(10),
        distance_density_centre_pc = r(11),
        eccentricity = r(12),
        semi_major_axis_rsun = r(13),
        period_d = r(14),
        iqcoll = i(15),
    )
end

"""
    read_binary_stage_log(path::AbstractString) -> Vector{BinaryStageRecord}

Read `binev.17`, the binary-evolution stages written by `binev.f`
(KZ(9) ≥ 2), one [`BinaryStageRecord`](@ref) per stage. The label line
(`Time[Myr] …`, written on the first call of a run started fresh) is skipped
wherever it occurs. Records without 15 values or with unparseable values are
skipped and counted in one warning; a missing file raises an error.
"""
function read_binary_stage_log(path::AbstractString)::Vector{BinaryStageRecord}
    return _binfile_read(
        _parse_binary_stage_record,
        BinaryStageRecord,
        path,
        "binev.17",
        s -> startswith(s, "Time[Myr]"),
    )
end

# ---------------------------------------------------------------------------
# degen.4: degenerate-star binaries
# ---------------------------------------------------------------------------

"""
    DegenerateRecord

One KS pair with a degenerate member in `degen.4`, written list-directed by
`degen.f` (14 values), from [`read_degenerate_log`](@ref).

| # | Label | Field | Unit |
|---|---|---|---|
| 1 | ICASE | `icase` (3 MDOT, 4 END ROCHE, 5 CMBODY) | — |
| 2 | Time[Myr] | `time_myr` (TTOT·T*) | Myr |
| 3 | SEMI[AU] | `semi_major_axis_au` (R(IPAIR) for a ghost pair) | AU |
| 4 | ECC | `eccentricity` (0 when R > 2a, after common envelope) | — |
| 5 | PERI/RS | `pericentre_over_radius` (a(1 − e)/max R*) | — |
| 6 | P[days] | `period_d` | d |
| 7 | RI[PC] | `cm_distance_pc` (RCM·R*) | pc |
| 8–9 | M(I1)[M*], M(I2)[M*] | `mass1_msun`, `mass2_msun` | M☉ |
| 10–12 | K*(I1), K*(I2), K*(ICM) | `kstar1`, `kstar2`, `kstar_cm` | — |
| 13–14 | NAME(I1), NAME(I2) | `name1`, `name2` | — |

Hyperbolic pairs with PERI/RS > 5 are not written. No sample run contains a
record; the layout is transcribed from the source.
"""
Base.@kwdef struct DegenerateRecord
    icase::Int
    time_myr::Float64
    semi_major_axis_au::Float64
    eccentricity::Float64
    pericentre_over_radius::Float64
    period_d::Float64
    cm_distance_pc::Float64
    mass1_msun::Float64
    mass2_msun::Float64
    kstar1::Int
    kstar2::Int
    kstar_cm::Int
    name1::Int
    name2::Int
end

"""Build a [`DegenerateRecord`](@ref) from one list-directed line of `degen.4`."""
function _parse_degenerate_record(line::AbstractString)::DegenerateRecord
    t = _binfile_tokens(line, 14, "degen.4")
    r(k) = _parse_fortran_real(t[k])
    i(k) = parse(Int, t[k])
    return DegenerateRecord(;
        icase = i(1),
        time_myr = r(2),
        semi_major_axis_au = r(3),
        eccentricity = r(4),
        pericentre_over_radius = r(5),
        period_d = r(6),
        cm_distance_pc = r(7),
        mass1_msun = r(8),
        mass2_msun = r(9),
        kstar1 = i(10),
        kstar2 = i(11),
        kstar_cm = i(12),
        name1 = i(13),
        name2 = i(14),
    )
end

"""
    read_degenerate_log(path::AbstractString) -> Vector{DegenerateRecord}

Read `degen.4`, the KS pairs with a white dwarf, neutron star or black hole
written by `degen.f` (KZ(9) ≥ 3), one [`DegenerateRecord`](@ref) per call.
The header block (MODEL line and 14-label line, first call of a run started
fresh) is skipped wherever it occurs. Records without 14 values or with
unparseable values are skipped and counted in one warning; a missing file
raises an error.
"""
function read_degenerate_log(path::AbstractString)::Vector{DegenerateRecord}
    return _binfile_read(
        _parse_degenerate_record,
        DegenerateRecord,
        path,
        "degen.4",
        s -> startswith(s, "MODEL:") || startswith(s, "ICASE"),
    )
end

# ---------------------------------------------------------------------------
# escbin.31: escaping binaries
# ---------------------------------------------------------------------------

"""FORMAT 501 of `escape.F`: `(1X,1P,3E13.5,2I10,2I4,12E13.5)`."""
const _ESCBIN_FORMAT =
    vcat(fmt_skip(1), fmt_real(3, 13), fmt_int(2, 10), fmt_int(2, 4), fmt_real(12, 13))

"""
    EscapedBinaryRecord

One escaping KS binary of `escbin.31`, written by `escape.F` (FORMAT 501,
fixed width, 19 values) at ADJUST, before the stdout BINARY ESCAPE line,
from [`read_escaped_binary_log`](@ref).

| # | Label | Field | Unit |
|---|---|---|---|
| 1 | TTOT | `time_nb` | NB |
| 2–3 | BODY1,2 | `mass1_nb`, `mass2_nb` | NB |
| 4–5 | NAME1,2 | `name1`, `name2` | — |
| 6–7 | KW1,2 | `kstar1`, `kstar2` | — |
| 8–9 | M1,2 [Sol] | `mass1_msun`, `mass2_msun` (from `hrdiag`, 0 if KZ(19) = 0) | M☉ |
| 10–11 | RAD1,2[Sol] | `radius1_rsun`, `radius2_rsun` | R☉ |
| 12–13 | LUM1,2[Sol] | `luminosity1_lsun`, `luminosity2_lsun` | L☉ |
| 14–15 | TEFF1,2 [K] | `log_teff1`, `log_teff2` (log₁₀(T_eff/K)) | — |
| 16–17 | AGE1,2[Myr] | `age1_myr`, `age2_myr` (TIME·T* − EPOCH) | Myr |
| 18–19 | EPOCH1,2 [Myr] | `epoch1_myr`, `epoch2_myr` | Myr |

Correction: columns 14–15 are labelled `[K]` but hold log₁₀ T_eff; they are
stored as `log_teff1`, `log_teff2`. No orbital elements are written: a, e
and P are on the BINARY ESCAPE line of the same ADJUST, matched by names.
`overflow` is `true` when a field was printed as asterisks (read as absent).
"""
Base.@kwdef struct EscapedBinaryRecord
    time_nb::Float64
    mass1_nb::Float64
    mass2_nb::Float64
    name1::Int
    name2::Int
    kstar1::Int
    kstar2::Int
    mass1_msun::Float64
    mass2_msun::Float64
    radius1_rsun::Float64
    radius2_rsun::Float64
    luminosity1_lsun::Float64
    luminosity2_lsun::Float64
    log_teff1::Float64
    log_teff2::Float64
    age1_myr::Float64
    age2_myr::Float64
    epoch1_myr::Float64
    epoch2_myr::Float64
    overflow::Bool = false
end

"""Build an [`EscapedBinaryRecord`](@ref) from one line of `escbin.31`."""
function _parse_escaped_binary_record(line::AbstractString)::EscapedBinaryRecord
    rec = _binfile_fixed(line, _ESCBIN_FORMAT)
    v = rec.values
    return EscapedBinaryRecord(;
        time_nb = v[1],
        mass1_nb = v[2],
        mass2_nb = v[3],
        name1 = _fixed_int(v[4]),
        name2 = _fixed_int(v[5]),
        kstar1 = _fixed_int(v[6], KSTAR_ABSENT),
        kstar2 = _fixed_int(v[7], KSTAR_ABSENT),
        mass1_msun = v[8],
        mass2_msun = v[9],
        radius1_rsun = v[10],
        radius2_rsun = v[11],
        luminosity1_lsun = v[12],
        luminosity2_lsun = v[13],
        log_teff1 = v[14],
        log_teff2 = v[15],
        age1_myr = v[16],
        age2_myr = v[17],
        epoch1_myr = v[18],
        epoch2_myr = v[19],
        overflow = rec.overflow,
    )
end

"""
    read_escaped_binary_log(path::AbstractString) -> Vector{EscapedBinaryRecord}

Read `escbin.31`, the escaping KS binaries written by `escape.F` (opened for
KZ(23) = 2 or 4), one fixed-width [`EscapedBinaryRecord`](@ref) per escaper.
The label line (FORMAT 502, guarded by an uninitialised flag and absent in
the samples) is skipped wherever it occurs. The `TEFF1,2 [K]` columns are
stored as `log_teff1`, `log_teff2`. Malformed records are skipped and counted
in one warning; a missing file raises an error.
"""
function read_escaped_binary_log(path::AbstractString)::Vector{EscapedBinaryRecord}
    return _binfile_read(
        _parse_escaped_binary_record,
        EscapedBinaryRecord,
        path,
        "escbin.31",
        s -> startswith(s, "TTOT"),
    )
end

# ---------------------------------------------------------------------------
# hiarch.12: hierarchy begin and end
# ---------------------------------------------------------------------------

"""
FORMAT 10 of `hiarch.f` (NEW HIAR): four label literals, then
`1P,E20.11,2E14.5,0P,F6.2,1P,8E13.4,0P,4I12,3F6.2,2I12,1P,E14.5,0P`.
"""
const _HIARCH_NEW_FORMAT = vcat(
    fmt_lit("NEW HIAR  Time SEMI0 SMEI1 ECC1 "),
    fmt_lit("PERI0 PERI0M P1/P0 M(INCM)/M(I3) PCR/SEMI0 "),
    fmt_lit("M(INCM)/<M> MR INA NAME(I1) NAME(I2) NAME(I3) "),
    fmt_lit("K*(INCM) ECC0 ECCMIN ECCMAX K*(I1) K*(I2) RSM "),
    fmt_real(1, 20),
    fmt_real(2, 14),
    fmt_real(1, 6),
    fmt_real(8, 13),
    fmt_int(4, 12),
    fmt_real(3, 6),
    fmt_int(2, 12),
    fmt_real(1, 14),
)

"""
FORMAT 20 of `hiarch.f` (END HIAR): two label literals, then
`1P,E20.11,2E14.5,0P,F6.2,1P,4E13.4,0P,2I4,2I12`.
"""
const _HIARCH_END_FORMAT = vcat(
    fmt_lit("END HIAR  Time RI/RC SEMI0 ECC0 PERI0 P0F/P0I "),
    fmt_lit("RC/RSCALE GAMMA(INCM) NKI NKF NPAIRS NAME(I2)"),
    fmt_real(1, 20),
    fmt_real(2, 14),
    fmt_real(1, 6),
    fmt_real(4, 13),
    fmt_int(2, 4),
    fmt_int(2, 12),
)

"""
    HierarchyLogRecord

One line of `hiarch.12` (`hiarch.f`), from [`read_hierarchy_log`](@ref):
`kind = :new` for NEW HIAR (FORMAT 10, a hierarchy formed; IPAIR is the
inner binary, JCOMP the outer body) or `kind = :end` for END HIAR (FORMAT
20, termination or escape; IPAIR is the merged pair, so its elements are
those of the outer orbit). Fields a kind does not carry are `NaN`, `0` or
`KSTAR_ABSENT`.

| NEW HIAR label | Field | Unit |
|---|---|---|
| Time | `time_nb` | NB |
| SEMI0, SMEI1 | `semi_major_axis_inner_nb`, `semi_major_axis_outer_nb` | NB |
| ECC1 (F6.2) | `eccentricity_outer` | — |
| PERI0 | `pericentre_inner_nb` (a(1 − e)) | NB |
| PERI0M | `pericentre_inner_min_nb` (a(1 − e_max)) | NB |
| P1/P0 | `period_ratio` (P_out/P_in) | — |
| M(INCM)/M(I3) | `mass_ratio_inner_outer` | — |
| PCR/SEMI0 | `pcrit_over_semi_inner` (PCRIT of the last IMPACT) | — |
| M(INCM)/<M> | `mass_inner_over_mean` | — |
| MR | `mass_ratio_components` (max(m₂/m₁, m₁/m₂)) | — |
| INA | `inclination_deg` (mutual inclination) | deg |
| NAME(I1), NAME(I2), NAME(I3) | `name1`, `name2`, `name_outer` | — |
| K*(INCM) | `kstar_inner_cm` (−10 when NAME(I) < 0, as printed) | — |
| ECC0 (F6.2) | `eccentricity_inner` (Runge–Lenz vector) | — |
| ECCMIN (F6.2) | `eccentricity_jz` (√(sin²i + e²cos²i) = √(1 − j_z²), the inner eccentricity at which the dimensionless angular momentum equals its conserved z-component j_z = √(1 − e²) cos i; not a minimum, despite the label) | — |
| ECCMAX (F6.2) | `eccentricity_max` (Heggie maximum) | — |
| K*(I1), K*(I2) | `kstar1`, `kstar2` | — |
| RSM | `radius_max_nb` (max stellar radius) | NB |

| END HIAR label | Field | Unit |
|---|---|---|
| Time | `time_nb` | NB |
| RI/RC | `distance_over_core_radius` (c.m. distance / core radius) | — |
| SEMI0 | `semi_major_axis_outer_nb` | NB |
| ECC0 (F6.2) | `eccentricity_outer` | — |
| PERI0 | `pericentre_outer_over_apocentre_inner` (/ RAP) | — |
| P0F/P0I | `period_ratio` (P_out now / P_in at the merger) | — |
| RC/RSCALE | `core_radius_over_rscale` | — |
| GAMMA(INCM) | `gamma` (perturbation of the merged pair) | — |
| NKI, NKF | `n_inner_periods`, `n_outer_periods` (since the merger) | — |
| NPAIRS | `n_pairs` | — |
| NAME(I2) | `name_outer` | — |

Corrections: END HIAR prints the outer orbit under the labels SEMI0, ECC0
and PERI0; the values are stored as outer-orbit quantities. Hazards: NKI (I4)
is glued to GAMMA and overflows at 10⁴ inner periods (common): 0 and
`overflow = true`; RAP is a single saved value shared by all systems, so the
END HIAR pericentre ratio is wrong when hierarchies overlap in time [UNSURE in
the dossier; not corrected]; the F6.2 eccentricities carry two decimals.
"""
Base.@kwdef struct HierarchyLogRecord
    kind::Symbol
    time_nb::Float64
    semi_major_axis_inner_nb::Float64 = NaN
    semi_major_axis_outer_nb::Float64 = NaN
    eccentricity_inner::Float64 = NaN
    eccentricity_outer::Float64 = NaN
    eccentricity_jz::Float64 = NaN
    eccentricity_max::Float64 = NaN
    pericentre_inner_nb::Float64 = NaN
    pericentre_inner_min_nb::Float64 = NaN
    period_ratio::Float64 = NaN
    mass_ratio_inner_outer::Float64 = NaN
    pcrit_over_semi_inner::Float64 = NaN
    mass_inner_over_mean::Float64 = NaN
    mass_ratio_components::Float64 = NaN
    inclination_deg::Float64 = NaN
    radius_max_nb::Float64 = NaN
    name1::Int = 0
    name2::Int = 0
    name_outer::Int = 0
    kstar1::Int = KSTAR_ABSENT
    kstar2::Int = KSTAR_ABSENT
    kstar_inner_cm::Int = KSTAR_ABSENT
    distance_over_core_radius::Float64 = NaN
    pericentre_outer_over_apocentre_inner::Float64 = NaN
    core_radius_over_rscale::Float64 = NaN
    gamma::Float64 = NaN
    n_inner_periods::Int = 0
    n_outer_periods::Int = 0
    n_pairs::Int = 0
    overflow::Bool = false
end

"""Build a [`HierarchyLogRecord`](@ref) from one NEW HIAR or END HIAR line of `hiarch.12`."""
function _parse_hierarchy_log_record(line::AbstractString)::HierarchyLogRecord
    if startswith(line, "NEW HIAR")
        rec = _binfile_fixed(line, _HIARCH_NEW_FORMAT)
        v = rec.values
        return HierarchyLogRecord(;
            kind = :new,
            time_nb = v[1],
            semi_major_axis_inner_nb = v[2],
            semi_major_axis_outer_nb = v[3],
            eccentricity_outer = v[4],
            pericentre_inner_nb = v[5],
            pericentre_inner_min_nb = v[6],
            period_ratio = v[7],
            mass_ratio_inner_outer = v[8],
            pcrit_over_semi_inner = v[9],
            mass_inner_over_mean = v[10],
            mass_ratio_components = v[11],
            inclination_deg = v[12],
            name1 = _fixed_int(v[13]),
            name2 = _fixed_int(v[14]),
            name_outer = _fixed_int(v[15]),
            kstar_inner_cm = _fixed_int(v[16], KSTAR_ABSENT),
            eccentricity_inner = v[17],
            eccentricity_jz = v[18],
            eccentricity_max = v[19],
            kstar1 = _fixed_int(v[20], KSTAR_ABSENT),
            kstar2 = _fixed_int(v[21], KSTAR_ABSENT),
            radius_max_nb = v[22],
            overflow = rec.overflow,
        )
    elseif startswith(line, "END HIAR")
        rec = _binfile_fixed(line, _HIARCH_END_FORMAT)
        v = rec.values
        return HierarchyLogRecord(;
            kind = :end,
            time_nb = v[1],
            distance_over_core_radius = v[2],
            # SEMI0, ECC0, PERI0 of END HIAR describe the outer orbit.
            semi_major_axis_outer_nb = v[3],
            eccentricity_outer = v[4],
            pericentre_outer_over_apocentre_inner = v[5],
            period_ratio = v[6],
            core_radius_over_rscale = v[7],
            gamma = v[8],
            n_inner_periods = _fixed_int(v[9]),
            n_outer_periods = _fixed_int(v[10]),
            n_pairs = _fixed_int(v[11]),
            name_outer = _fixed_int(v[12]),
            overflow = rec.overflow,
        )
    end
    throw(ArgumentError("hiarch.12 line is neither NEW HIAR nor END HIAR"))
end

"""
    read_hierarchy_log(path::AbstractString) -> Vector{HierarchyLogRecord}

Read `hiarch.12`, the hierarchy formations (NEW HIAR) and terminations or
escapes (END HIAR) written by `hiarch.f` (KZ(18) = 1 or 3), one fixed-width
[`HierarchyLogRecord`](@ref) per line in file order. The MODEL header line
(first call of a run started fresh) is skipped wherever it occurs. END HIAR
elements are stored as outer-orbit quantities, and the NKI field glued to
GAMMA is cut at its width. Malformed lines are skipped and counted in one
warning; a missing file raises an error.
"""
function read_hierarchy_log(path::AbstractString)::Vector{HierarchyLogRecord}
    return _binfile_read(
        _parse_hierarchy_log_record,
        HierarchyLogRecord,
        path,
        "hiarch.12",
        s -> startswith(s, "MODEL:"),
    )
end

# ---------------------------------------------------------------------------
# hbin.39: hardest binary
# ---------------------------------------------------------------------------

"""
FORMAT 62 of `adjust.F`: two label literals, then
`1P,E26.17,0P,2I12,I4,I6,1P,5E15.6,0P`.
"""
const _HBIN_FORMAT = vcat(
    fmt_lit(" BINARY:   Time[Myr] NAME(I1) NAME(I2) "),
    fmt_lit("K*(ICM) NP ECC SEMI[NB] P[days] EB[NB] EM[NB] "),
    fmt_real(1, 26),
    fmt_int(2, 12),
    fmt_int(1, 4),
    fmt_int(1, 6),
    fmt_real(5, 15),
)

"""
    HardestBinaryRecord

The most bound KS pair (BODY(ICM) > 0, H < −ECLOSE) at one OUTPUT, written
by `adjust.F` to `hbin.39` (FORMAT 62, label text in every record), from
[`read_hardest_binary_log`](@ref).

| # | Label | Field | Unit |
|---|---|---|---|
| 1 | Time[Myr] | `time_nb` (TTOT) | NB |
| 2–3 | NAME(I1), NAME(I2) | `name1`, `name2` | — |
| 4 | K*(ICM) | `kstar_cm` | — |
| 5 | NP | `n_neighbours` (LIST(1,I1)) | — |
| 6 | ECC | `eccentricity` | — |
| 7 | SEMI[NB] | `semi_major_axis_nb` | NB |
| 8 | P[days] | `period_d` | d |
| 9 | EB[NB] | `binding_energy_nb` (μH) | NB |
| 10 | EM[NB] | `total_energy_nb` (E(3) = ZKIN − POT + ETIDE) | NB |

Correction: column 1 is labelled `Time[Myr]` but holds TTOT in N-body units;
it is stored as `time_nb`. `overflow` is `true` when a field was printed as
asterisks (read as absent).
"""
Base.@kwdef struct HardestBinaryRecord
    time_nb::Float64
    name1::Int
    name2::Int
    kstar_cm::Int
    n_neighbours::Int
    eccentricity::Float64
    semi_major_axis_nb::Float64
    period_d::Float64
    binding_energy_nb::Float64
    total_energy_nb::Float64
    overflow::Bool = false
end

"""Build a [`HardestBinaryRecord`](@ref) from one line of `hbin.39`."""
function _parse_hardest_binary_record(line::AbstractString)::HardestBinaryRecord
    rec = _binfile_fixed(line, _HBIN_FORMAT)
    v = rec.values
    return HardestBinaryRecord(;
        time_nb = v[1],
        name1 = _fixed_int(v[2]),
        name2 = _fixed_int(v[3]),
        kstar_cm = _fixed_int(v[4], KSTAR_ABSENT),
        n_neighbours = _fixed_int(v[5]),
        eccentricity = v[6],
        semi_major_axis_nb = v[7],
        period_d = v[8],
        binding_energy_nb = v[9],
        total_energy_nb = v[10],
        overflow = rec.overflow,
    )
end

"""
    read_hardest_binary_log(path::AbstractString) -> Vector{HardestBinaryRecord}

Read `hbin.39`, the hardest binary written by `adjust.F` (KZ(9) = 1 or 3) at
the OUTPUT cadence (DELTAT), not at every ADJUST, one fixed-width
[`HardestBinaryRecord`](@ref) per output time; the file is appended on
restart and has no header. The `Time[Myr]` column is stored as `time_nb`, its
true unit. Malformed lines are skipped and counted in one warning; a missing
file raises an error.
"""
function read_hardest_binary_log(path::AbstractString)::Vector{HardestBinaryRecord}
    return _binfile_read(
        _parse_hardest_binary_record,
        HardestBinaryRecord,
        path,
        "hbin.39",
        _ -> false,
    )
end
