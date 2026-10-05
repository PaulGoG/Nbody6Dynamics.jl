# =============================================================================
# Restart dump selection and the join of segments
# =============================================================================

"""
    _dump_complete(path) -> Bool

Whether `path` holds a complete Fortran unformatted sequential file: a chain
of records `[n::Int32][n bytes][n::Int32]` with matching length markers,
ending exactly at a record boundary, with at least two records. A dump cut
short by a kill ends inside a record and fails the test. `false` when the
file is absent.
"""
function _dump_complete(path::AbstractString)::Bool
    isfile(path) || return false
    nbytes = filesize(path)
    return open(path) do io
        nrec = 0
        while !eof(io)
            nbytes - position(io) < 4 && return false
            n = read(io, Int32)
            (n < 0 || nbytes - position(io) < Int(n) + 4) && return false
            skip(io, n)
            read(io, Int32) == n || return false
            nrec += 1
        end
        return position(io) == nbytes && nrec ≥ 2
    end
end

"""
    _resume_dump(out_dir, stdout_path) -> Union{Nothing,_DumpMarker}

Dump an interrupted run is restarted from: the last dump reported in the
stdout capture ([`_dump_markers`](@ref)) whose file in `out_dir` is complete
([`_dump_complete`](@ref)). The markers are searched from the last to the
first; a file name is examined once, at its last marker, because a dump
rewritten later holds what that last marker wrote. When no reported dump
qualifies, the latest dump file by name ([`_latest_dump`](@ref)) is taken if
complete, with `line_end = -1`: such a dump has no stdout line, its time is
the rounded one of its name, and the stdout capture cannot be cut at it.
`nothing` when no complete dump exists.
"""
function _resume_dump(
    out_dir::AbstractString,
    stdout_path::AbstractString,
)::Union{Nothing,_DumpMarker}
    markers = _dump_markers(stdout_path)
    seen = Set{String}()
    for mk in Iterators.reverse(markers)
        mk.file in seen && continue
        push!(seen, mk.file)
        _dump_complete(joinpath(out_dir, mk.file)) && return mk
    end
    name = _latest_dump(out_dir)
    if name !== nothing && _dump_complete(joinpath(out_dir, name))
        return (file = name, time_nb = _dump_time(name), line_end = -1)
    end
    return nothing
end

"""
    _marker_block_end(stdout_path, marker) -> Int

Byte position after the lines the engine prints with a dump: from
`marker.line_end` on, the ` W MYDUMP` and ` NA-NS=` lines that follow the
` MYDUMP` line are skipped. Returns `marker.line_end` when none follows.
Raises an `ArgumentError` for a marker without a stdout line
(`line_end < 0`).
"""
function _marker_block_end(stdout_path::AbstractString, marker::_DumpMarker)::Int
    marker.line_end ≥ 0 || throw(
        ArgumentError("dump $(marker.file) has no stdout line (line_end = $(marker.line_end))"),
    )
    return open(stdout_path) do io
        seek(io, marker.line_end)
        while !eof(io)
            pos = position(io)
            s = lstrip(readline(io; keep = true))
            (startswith(s, "W MYDUMP") || startswith(s, "NA-NS=")) || return pos
        end
        return position(io)
    end
end

"""
    _series_cut_offset(path, t_max) -> Union{Nothing,Int}

Byte offset of the first record line of a series file whose time (first
column) exceeds `t_max` by more than a relative tolerance of 1e-9; lines
that are empty or whose first token is not a number are headers and are
skipped. `nothing` when no record lies beyond `t_max` or the file is absent.
"""
function _series_cut_offset(path::AbstractString, t_max::Real)::Union{Nothing,Int}
    isfile(path) || return nothing
    t_cut = Float64(t_max) + 1.0e-9 * max(1.0, abs(Float64(t_max)))
    return open(path) do io
        while !eof(io)
            pos = position(io)
            tokens = split(readline(io; keep = true))
            isempty(tokens) && continue
            t = tryparse(Float64, first(tokens))
            t === nothing && continue
            t > t_cut && return pos
        end
        return nothing
    end
end

"""
    _move_tail(path, offset, dest) -> Int

Move the bytes of `path` from `offset` to the end into the new file `dest`
and truncate `path` to `offset`; returns the number of bytes moved. An
existing `dest` is kept as a `#k` sibling ([`_backup_existing`](@ref)), and
the tail is written before `path` is truncated. Nothing is written when
`offset` is at or beyond the end of `path`.
"""
function _move_tail(path::AbstractString, offset::Integer, dest::AbstractString)::Int
    offset ≥ 0 || throw(ArgumentError("offset must be ≥ 0; got $offset"))
    offset ≥ filesize(path) && return 0
    mkpath(dirname(dest))
    _backup_existing(dest)
    tail = open(path) do io
        seek(io, offset)
        read(io)
    end
    write(dest, tail)
    open(path, "r+") do io
        truncate(io, offset)
    end
    return length(tail)
end

"""
    _time_unit_myr(out_dir) -> Float64

Time unit T* in Myr per N-body time unit, from the first record of
`global.30` in `out_dir` with a positive N-body time (column 1 in N-body
units, column 2 in Myr). `NaN` when the file is absent or holds no such
record.
"""
function _time_unit_myr(out_dir::AbstractString)::Float64
    path = joinpath(out_dir, "global.30")
    isfile(path) || return NaN
    return open(path) do io
        for line in eachline(io)
            tokens = split(line)
            length(tokens) ≥ 2 || continue
            t_nb = tryparse(Float64, tokens[1])
            t_myr = tryparse(Float64, tokens[2])
            (t_nb === nothing || t_myr === nothing) && continue
            t_nb > 0 && return t_myr / t_nb
        end
        return NaN
    end
end

"""
    _timestamped_time(name) -> Float64

Time in the name of a time-stamped output file `<letters>.<digits>_<time>`
(`conf.3_12.5` → 12.5, `sev.83_3` → 3.0); `NaN` for any other name.
"""
function _timestamped_time(name::AbstractString)::Float64
    m = match(r"^[A-Za-z]+\.\d+_([0-9][0-9.]*(?:[eE][+-]?\d+)?)$", name)
    m === nothing && return NaN
    return something(tryparse(Float64, m.captures[1]), NaN)
end

"""
    _discard_after(out_dir, stdout_path, marker; segment, last_status) -> Dict{String,Any}

Prepare `out_dir` for a restart from the dump `marker` by moving aside what
the interrupted segment `segment` wrote after that dump, which the next
segment writes again. Moved to `discarded/segment_<segment>/` in `out_dir`:

  - the stdout capture after the dump's block of lines, as `<name>.tail`
    (only when the marker has a stdout line, `line_end ≥ 0`);
  - the records of `lagr.7`, `global.30`, `esc.11` (N-body units) and
    `event.35` (Myr, converted with [`_time_unit_myr`](@ref)) beyond the
    dump time, as `<name>.tail`;
  - the time-stamped files ([`_timestamped_time`](@ref)) with a time beyond
    the dump time, whole, except the dump itself.

A segment with `last_status == "stopped"` restarted from the last dump it
reported is left as it is: a stop writes its dump last, so nothing follows
it. Nothing is deleted: tails are written before a file is truncated, and an
earlier file of the same name in the destination is kept as a `#k` sibling.

The returned record holds `dir` (the destination relative to `out_dir`, or
`""` when nothing was moved), `moved` (sorted names under it), `t_from` (dump
time) and `t_to` (time of the last `ADJUST:` line, the dump time when none),
and `uncut`: the non-empty appended files the package cannot cut by time
(`<letters>.<digits>` outside the four series and `dat.10`, and `err1000`),
which hold the interval `(t_from, t_to]` twice after the restart.
"""
function _discard_after(
    out_dir::AbstractString,
    stdout_path::AbstractString,
    marker::_DumpMarker;
    segment::Integer,
    last_status::AbstractString,
)::Dict{String,Any}
    t_d = marker.time_nb
    t_to = _last_adjust_time(stdout_path)
    result = Dict{String,Any}(
        "dir" => "",
        "moved" => String[],
        "uncut" => String[],
        "t_from" => t_d,
        "t_to" => isnan(t_to) ? t_d : t_to,
    )
    if last_status == "stopped"
        markers = _dump_markers(stdout_path)
        if !isempty(markers) &&
           last(markers).file == marker.file &&
           last(markers).line_end == marker.line_end
            return result
        end
    end

    dest_dir = joinpath(out_dir, "discarded", "segment_$(segment)")
    tol(t) = 1.0e-9 * max(1.0, abs(t))
    moved = String[]

    # Stdout capture after the dump's block of lines
    if marker.line_end ≥ 0
        cut = _marker_block_end(stdout_path, marker)
        if cut < filesize(stdout_path)
            tail_name = basename(stdout_path) * ".tail"
            _move_tail(stdout_path, cut, joinpath(dest_dir, tail_name))
            push!(moved, tail_name)
        end
    end

    # Series files, cut by time; T* is read before global.30 is cut
    t_unit = _time_unit_myr(out_dir)
    for (name, factor) in
        (("lagr.7", 1.0), ("global.30", 1.0), ("esc.11", 1.0), ("event.35", t_unit))
        path = joinpath(out_dir, name)
        (isfile(path) && !isnan(factor)) || continue
        off = _series_cut_offset(path, t_d * factor)
        off === nothing && continue
        _move_tail(path, off, joinpath(dest_dir, name * ".tail"))
        push!(moved, name * ".tail")
    end

    # Time-stamped files written after the dump
    for f in readdir(out_dir)
        f == marker.file && continue
        src = joinpath(out_dir, f)
        isfile(src) || continue
        t = _timestamped_time(f)
        (!isnan(t) && t > t_d + tol(t_d)) || continue
        mkpath(dest_dir)
        dest = joinpath(dest_dir, f)
        _backup_existing(dest)
        mv(src, dest)
        push!(moved, f)
    end

    # Appended files without a time column
    series = ("lagr.7", "global.30", "esc.11", "event.35", "dat.10")
    uncut = sort!(
        filter(readdir(out_dir)) do f
            p = joinpath(out_dir, f)
            isfile(p) && filesize(p) > 0 && occursin(r"^[A-Za-z]+\.\d+$", f) && !(f in series)
        end,
    )
    err_path = joinpath(out_dir, "err1000")
    isfile(err_path) && filesize(err_path) > 0 && push!(uncut, "err1000")

    result["moved"] = sort(moved)
    result["uncut"] = uncut
    result["dir"] = isempty(moved) ? "" : joinpath("discarded", "segment_$(segment)")
    return result
end

"""
    _tcrit_of(inp_path) -> Float64

End time `TCRIT` of an engine input file: the first `TCRIT=` entry outside
comment lines (first non-blank character `!`), Fortran `D` exponents
accepted. Raises an error when the file has no such entry.
"""
function _tcrit_of(inp_path::AbstractString)::Float64
    tcrit = open(inp_path) do io
        for line in eachline(io)
            startswith(lstrip(line), '!') && continue
            m = match(r"\bTCRIT\s*=\s*([0-9.eEdD+-]+)", line)
            m === nothing || return _parse_fortran_float(m.captures[1])
        end
        return nothing
    end
    tcrit === nothing && error("no TCRIT entry in $inp_path")
    return tcrit
end
