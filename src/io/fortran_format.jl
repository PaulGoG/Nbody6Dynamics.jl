# =============================================================================
# Fixed-width reading of records written through a Fortran FORMAT
# =============================================================================
#
# The engine's formatted event lines and files glue fields together whenever
# a value fills its width (I-fields of names, F-fields after labels without a
# trailing blank, a negative 1PEw.d with w = d + 7), so whitespace splitting
# misreads them. A record is read instead by walking the edit descriptors of
# its FORMAT: literal strings are checked, value fields are cut at their
# widths. `1P`/`0P` scale factors change only the printed mantissa, never a
# width, and are ignored when reading.

"""
    FortranItem

One edit descriptor of a FORMAT, for reading: a literal string (`kind =
:lit`, `text` the characters the engine prints), a skipped field of `width`
characters (`:skip`, `nX` or a field whose value is not wanted), or a value
field of `width` characters read as an integer (`:int`, `Iw`) or a real
(`:real`, `Fw.d`, `Ew.d`, `Dw.d`). Build formats with [`fmt_lit`](@ref),
[`fmt_int`](@ref), [`fmt_real`](@ref) and [`fmt_skip`](@ref).
"""
struct FortranItem
    kind::Symbol
    width::Int
    text::String
end

"""`fmt_lit(text)`: the literal `text` of a FORMAT (`'…'`), checked when reading."""
fmt_lit(text::AbstractString) = [FortranItem(:lit, length(text), String(text))]

"""`fmt_int(n, w)`: `n` integer fields of width `w` (`nIw`)."""
fmt_int(n::Integer, w::Integer) = [FortranItem(:int, w, "") for _ in 1:n]

"""`fmt_real(n, w)`: `n` real fields of width `w` (`nFw.d`, `nEw.d`, `nDw.d`)."""
fmt_real(n::Integer, w::Integer) = [FortranItem(:real, w, "") for _ in 1:n]

"""`fmt_skip(w)`: `w` characters that are not read (`wX`, or an unwanted field)."""
fmt_skip(w::Integer) = [FortranItem(:skip, w, "")]

"""
    FortranRecord

Values of one record read by [`read_fixed`](@ref), in descriptor order:
`values[k]` is the k-th value field (integers converted exactly, names up to
10 digits included), `NaN` for an empty field or one the engine overflowed
with asterisks; `overflow` is `true` when any field was overflowed.
"""
struct FortranRecord
    values::Vector{Float64}
    overflow::Bool
end

"""
    read_fixed(line, items; start = 1, strict_literals = true) -> FortranRecord

Read one record of `line` by the FORMAT `items` (a vector of
[`FortranItem`](@ref), e.g. `vcat(fmt_lit(" T ="), fmt_real(1, 12))`),
starting at character `start`. Literals must match the line character for
character (`ArgumentError` otherwise, naming the literal); with
`strict_literals = false` only their non-blank characters are compared.
A line shorter than the format ends the record early: the missing fields
read as empty. Reals accept the D exponent and the three-digit exponents
Fortran writes without a letter (`1.234-105`); integers accept a trailing
`.`-free decimal only.
"""
function read_fixed(
    line::AbstractString,
    items::AbstractVector{FortranItem};
    start::Integer = 1,
    strict_literals::Bool = true,
)::FortranRecord
    chars = collect(line)
    pos = Int(start)
    values = Float64[]
    overflow = false
    for item in items
        stop = min(pos + item.width - 1, length(chars))
        field = pos <= length(chars) ? String(chars[pos:stop]) : ""
        if item.kind === :lit
            ok =
                strict_literals ? field == item.text :
                filter(!isspace, field) == filter(!isspace, item.text)
            ok || throw(
                ArgumentError(
                    "FORMAT literal $(repr(item.text)) not found at column $pos: $(repr(field))",
                ),
            )
        elseif item.kind === :int || item.kind === :real
            s = strip(field)
            if isempty(s)
                push!(values, NaN)
            elseif all(==('*'), s)
                push!(values, NaN)
                overflow = true
            elseif item.kind === :int
                push!(values, Float64(parse(Int, s)))
            else
                push!(values, _parse_fortran_real(s))
            end
        end
        pos += item.width
    end
    return FortranRecord(values, overflow)
end

"""
Real of a Fortran output field: `D` exponents become `E`, and an exponent
written without its letter (`1.08726-157`, magnitude above 99) gets one.
"""
function _parse_fortran_real(s::AbstractString)::Float64
    t = replace(s, r"[dD]" => "E")
    if !occursin(r"[eE]", t)
        m = match(r"^([+-]?[0-9]*\.?[0-9]+)([+-][0-9]+)$", t)
        m === nothing || (t = something(m.captures[1]) * "E" * something(m.captures[2]))
    end
    return parse(Float64, t)
end

"""Integer of a value read by [`read_fixed`](@ref); `missing_value` when the field was empty or overflowed."""
_fixed_int(v::Float64, missing_value::Int = 0)::Int = isnan(v) ? missing_value : Int(v)
