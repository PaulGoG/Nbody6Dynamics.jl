# =============================================================================
# Low-level Fortran unformatted sequential binary record reader
# =============================================================================
#
# Fortran WRITE(unit) produces records bracketed by 4-byte markers:
#   [marker₁ :: Int32]  [data :: marker₁ bytes]  [marker₂ :: Int32]
# where marker₁ == marker₂ == byte length of the data payload.

"""
    read_fortran_record(io::IO) -> Vector{UInt8}

Read one Fortran unformatted sequential record and return the raw data bytes.
Throws `EOFError` at end-of-file.
"""
function read_fortran_record(io::IO)::Vector{UInt8}
    marker1 = read(io, Int32)
    data = read(io, marker1)
    marker2 = read(io, Int32)
    marker1 == marker2 || error(
        "Fortran record marker mismatch: $marker1 vs $marker2 " *
        "(file may be corrupted or use non-standard record format)",
    )
    return data
end

"""
    read_fortran_record(io::IO, ::Type{T}, n::Int) -> Vector{T}

Read a Fortran record and reinterpret as a vector of `n` elements of type `T`.
"""
function read_fortran_record(io::IO, ::Type{T}, n::Int) where {T}
    data = read_fortran_record(io)
    expected = n * sizeof(T)
    length(data) == expected || error(
        "Record size mismatch: got $(length(data)) bytes, expected $expected for $n × $(sizeof(T))-byte elements",
    )
    return reinterpret(T, data) |> collect
end

"""
    peek_record_size(io::IO) -> Int32

Read the next record's byte length without consuming the record.
Restores the stream position.
"""
function peek_record_size(io::IO)::Int32
    pos = position(io)
    marker = read(io, Int32)
    seek(io, pos)
    return marker
end
