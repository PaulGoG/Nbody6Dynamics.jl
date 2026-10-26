# =============================================================================
# Reader of mbh.45, the IMBH record of the engine (bhplot.f, KZ(24) = 1, 3)
# =============================================================================
#
# bhplot.f writes one label line (first token STAT, only on a run started
# fresh) and, every DTBH, one list-directed record: WHICH (A8, BINARY or
# SINGLE) followed by 20 numbers, 21 from fork v2026.07.11 on (the body name
# last). Restarted runs append to the file.

"""Numbers per record of `mbh.45` after the `WHICH` token before fork v2026.07.11 (`bhplot.f:197-199`)."""
const _MBH45_RECORD_VALUES = 20

"""Numbers per record of `mbh.45` after the `WHICH` token from fork v2026.07.11 on (the name last)."""
const _MBH45_RECORD_VALUES_NAMED = 21

"""States of the `WHICH` column of `mbh.45`."""
const _MBH45_STATES = Dict("BINARY" => :binary, "SINGLE" => :single)

"""
    MBHRow

One record of `mbh.45`, from [`read_mbh45`](@ref).

| # | Label | Field | Unit |
|---|---|---|---|
| 1 | STAT | `state`: `:binary` (the IMBH is a KS component) or `:single` | — |
| 2 | Time[Myr] | `time_myr` ((TIME + TOFF)·T*) | Myr |
| 3 | IBH | `index`: the engine's array index of the IMBH at that time | — |
| 4–6 | X(1:3)[PC] | `pos_pc` (absolute, not about the density centre) | pc |
| 7–9 | V(1:3)[km/s] | `vel_kms` (absolute) | km s⁻¹ |
| 10 | NB | `n_neighbours` (NNBT: the neighbour list, plus the KS partner) | — |
| 11–13 | XAVE(1:3)[AU] | `xave_au`, mass-weighted mean neighbour offset | au |
| 14–16 | VAVE(1:3)[km/s] | `vave_kms`, mass-weighted mean neighbour velocity | km s⁻¹ |
| 17 | DEN[M*/PC^3] | `density_msun_pc3` | M☉ pc⁻³ |
| 18 | RIJMAX[PC] | `rijmax_pc` | pc |
| 19–21 | VSIGMA(1:3)[km^2/s^2] | `vsigma_km2s2` | km² s⁻² |
| 22 | NAME | `name`: the body name of the IMBH (fork v2026.07.11 on; 0 when the record has no name column) | — |

`index` is not the name: it changes whenever the engine reorders its
particle arrays (KS initialisation and termination, escape removal).
"""
struct MBHRow
    state::Symbol
    time_myr::Float64
    index::Int
    pos_pc::NTuple{3,Float64}
    vel_kms::NTuple{3,Float64}
    n_neighbours::Int
    xave_au::NTuple{3,Float64}
    vave_kms::NTuple{3,Float64}
    density_msun_pc3::Float64
    rijmax_pc::Float64
    vsigma_km2s2::NTuple{3,Float64}
    name::Int
end

"""Integer of a list-directed token; a real-valued token is accepted when integral."""
function _mbh45_int(s::AbstractString)::Int
    v = tryparse(Int, s)
    v === nothing || return v
    x = _parse_fortran_real(s)
    isinteger(x) || throw(ArgumentError("non-integer value $(repr(s))"))
    return Int(x)
end

"""
Build an [`MBHRow`](@ref) from the tokens of one `mbh.45` record; a 22nd token
is the body name (a list-directed integer or real, rounded), 0 without it.
"""
function _parse_mbh45_record(tokens::AbstractVector{<:AbstractString})::MBHRow
    state = get(_MBH45_STATES, uppercase(tokens[1]), :none)
    state === :none && throw(ArgumentError("unknown state $(repr(tokens[1]))"))
    r(k) = _parse_fortran_real(tokens[k + 1])
    t3(k) = (r(k), r(k + 1), r(k + 2))
    name = if length(tokens) == _MBH45_RECORD_VALUES_NAMED + 1
        x = r(_MBH45_RECORD_VALUES_NAMED)
        isfinite(x) || throw(ArgumentError("non-finite name $(repr(tokens[end]))"))
        round(Int, x)
    else
        0
    end
    return MBHRow(
        state,
        r(1),
        _mbh45_int(tokens[3]),
        t3(3),
        t3(6),
        _mbh45_int(tokens[10]),
        t3(10),
        t3(13),
        r(16),
        r(17),
        t3(18),
        name,
    )
end

"""
    read_mbh45(path::AbstractString) -> Vector{MBHRow}

Read `mbh.45`, the IMBH record written by the engine's `bhplot.f` every DTBH
under `KZ(24) = 1` or `3` (fork `bhplot.f:28-40, 197-199`): a label line whose
first token is `STAT`, then one list-directed record per output, `WHICH`
(`BINARY` or `SINGLE`) followed by 20 numbers, or 21 from fork v2026.07.11 on,
the last being the body name of the IMBH (`name`, written as a list-directed
integer or real and rounded; 0 for 20-number records), read into one
[`MBHRow`](@ref) each; `name` tells the IMBHs of `KZ(24) = 3` apart. Blank
lines and label lines are skipped wherever they occur (a restarted run
appends). A record with another token count or an unreadable value raises an
`ArgumentError` naming the line number; a missing file raises an error.

Defects of the engine's writer that the reader cannot repair:

- The separation `RIJ` of `bhplot.f:133` and `bhplot.f:157` is overwritten
  by the square of each coordinate component instead of summed, so it holds
  |Δz| alone: `rijmax_pc` is the largest |Δz| over the companion and the
  neighbours, not the largest separation, and `density_msun_pc3` inherits it
  through its volume. The companion term of `VSIGMA` (`bhplot.f:182`) is
  weighted by `BODY(I)`, `I` left at the last neighbour of the list loop,
  instead of the companion's mass. The three columns are not usable as a
  local density or dispersion.
- The positions and velocities are absolute, not about the density centre;
  [`imbh_trajectory`](@ref) gives the offset from the stellar density centre
  at the snapshots.
- Nothing is written while the IMBH sits in a chain or in a merged
  hierarchy (`NAME(ICM) ≤ 0`), so gaps in the record are expected.
- `index` is the array index of the IMBH at the output, which changes as the
  engine reorders its particles; it is not the name.

# Example
```julia
rows = read_mbh45(joinpath(run_dir, "output", "mbh.45"))
binary_fraction = count(r -> r.state === :binary, rows) / length(rows)
```
"""
function read_mbh45(path::AbstractString)::Vector{MBHRow}
    isfile(path) || error("IMBH record (mbh.45) not found: $path")
    rows = MBHRow[]
    for (n, line) in enumerate(eachline(path))
        tokens = split(line)
        isempty(tokens) && continue
        uppercase(tokens[1]) == "STAT" && continue
        length(tokens) in (_MBH45_RECORD_VALUES + 1, _MBH45_RECORD_VALUES_NAMED + 1) || throw(
            ArgumentError(
                "mbh.45 line $n of $path: $(length(tokens)) tokens, expected " *
                "$(_MBH45_RECORD_VALUES + 1) or $(_MBH45_RECORD_VALUES_NAMED + 1) (WHICH and " *
                "$(_MBH45_RECORD_VALUES) values, or $(_MBH45_RECORD_VALUES_NAMED) with the name)",
            ),
        )
        row = try
            _parse_mbh45_record(tokens)
        catch e
            e isa ArgumentError || rethrow()
            throw(ArgumentError("mbh.45 line $n of $path: $(e.msg)"))
        end
        push!(rows, row)
    end
    return rows
end
