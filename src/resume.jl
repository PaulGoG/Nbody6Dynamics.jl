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

A segment with `last_status` `"stopped"` or `"completed"` restarted from the
last dump it reported is left as it is: a segment that ended by itself
wrote its dump last, and the lines that close it stay. Nothing is deleted:
tails are written before a file is truncated, and an earlier file of the
same name in the destination is kept as a `#k` sibling.

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
    if last_status in ("stopped", "completed")
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

"""
    _segment_alive(segment, out_dir) -> Bool

Whether the engine of the segment record `segment` is still running on this
host: its `pid` exists and the command line of that process names `out_dir`,
from which the engine binary is launched. A record from another host cannot
be checked and counts as not alive, as does a process that ends while its
command line is read.
"""
function _segment_alive(segment::AbstractDict, out_dir::AbstractString)::Bool
    get(segment, "host", "") == gethostname() || return false
    pid = get(segment, "pid", 0)
    (pid isa Integer && pid > 0) || return false
    cmdline_path = "/proc/$(pid)/cmdline"
    isfile(cmdline_path) || return false
    return try
        occursin(abspath(out_dir), read(cmdline_path, String))
    catch e
        e isa Union{SystemError,Base.IOError} || rethrow()
        false
    end
end

"""
    _last_status(segment, stdout_path) -> String

Status of the segment record `segment`: its `status` entry when present and
not empty. A record written before that entry existed gives `completed` when
its `completed` flag is set, and otherwise the status
[`_segment_status`](@ref) reads from the stdout capture `stdout_path`, from
the record's `stdout_offset` on and with its `exit_status`.
"""
function _last_status(segment::AbstractDict, stdout_path::AbstractString)::String
    status = get(segment, "status", "")
    isempty(status) || return String(status)
    get(segment, "completed", false) === true && return "completed"
    return _segment_status(
        stdout_path;
        offset = Int(get(segment, "stdout_offset", 0)),
        exit_status = get(segment, "exit_status", nothing),
    )
end

"""
    _launch_restart(cfg, run_dir, out_dir, original_inp, dump, tcrit_extra;
                    t_start, base_dir, tcrtp0 = nothing, discarded = nothing) -> String

Launch a restart segment of the run in `run_dir` from the dump `dump` (file
name in `out_dir`) for `tcrit_extra` N-body time units. The dump is copied to
`comm.1`, the file `KSTART = 2` reads; the engine rewrites the dump it
started from under that dump's own name, so the copy is a working copy and
not a link, and it is removed when the engine exits, whatever the outcome.
The restart input `restart.inp` is written from `original_inp`
([`_write_restart_inp`](@ref)), an earlier one kept as a `#k` sibling. The
segment record carries `t_start` and, when given, the record `discarded`
of [`_discard_after`](@ref). Returns `run_dir`.
"""
function _launch_restart(
    cfg::Nbody6Config,
    run_dir::AbstractString,
    out_dir::AbstractString,
    original_inp::AbstractString,
    dump::AbstractString,
    tcrit_extra::Real;
    t_start::Real,
    base_dir::AbstractString,
    tcrtp0::Union{Nothing,Real} = nothing,
    discarded::Union{Nothing,AbstractDict} = nothing,
)::String
    target = joinpath(out_dir, "comm.1")
    # No backup: the dump stays under its own name.
    cp(joinpath(out_dir, dump), target; force = true)
    restart_inp = joinpath(out_dir, "restart.inp")
    _backup_existing(restart_inp)
    _write_restart_inp(restart_inp, original_inp, tcrit_extra; tcrtp0 = tcrtp0)
    restart =
        discarded === nothing ?
        (dump = String(dump), tcrit_extra = Float64(tcrit_extra), t_start = Float64(t_start)) :
        (
            dump = String(dump),
            tcrit_extra = Float64(tcrit_extra),
            t_start = Float64(t_start),
            discarded = discarded,
        )
    try
        _execute_simulation(
            cfg,
            run_dir,
            out_dir,
            restart_inp;
            base_dir = base_dir,
            label = "restart",
            restart = restart,
        )
    finally
        rm(target; force = true)
    end
    return run_dir
end

"""
    resume_run(run_dir; base_dir = nothing) -> String

Continue an interrupted run to the end time of its original input. The call
is idempotent: a completed run is left as it is, and a run that was
stopped at its wall budget, killed, or ended by an engine error gets one
more segment, started from the last complete restart dump.

The dump is the last one the engine reported in its stdout whose file is
complete, at the exact time of that report. If the interrupted segment
wrote anything after that dump (it was killed, not stopped), that part is
moved to `output/discarded/segment_<k>/` before the launch: the tail of the
stdout capture, the later records of `lagr.7`, `global.30`, `esc.11` and
`event.35`, and the later snapshots and dumps. Appended files the package
cannot cut by time are listed in the new segment's record (`discarded.uncut`).
Nothing is deleted.

A run halted by the engine's energy check or terminated by the start-up
watchdog is refused: the first needs a decision on the accuracy parameters,
the second would hang again. A run whose record says `running` is refused
while its engine is alive on this host; otherwise the record is closed from
what the stdout shows.

# Arguments
- `run_dir`: run directory of a pipeline run (`RUN_INFO.toml`, `config.toml`, `output/`)
- `base_dir`: directory the engine tree is resolved against; default: that of the
  run's frozen `config.toml`, which carries absolute paths

# Returns
`run_dir` as an absolute path. Raises an `ArgumentError` for a run that
cannot be resumed, with the reason.
"""
function resume_run(
    run_dir::AbstractString;
    base_dir::Union{Nothing,AbstractString} = nothing,
)::String
    run_dir = abspath(run_dir)
    out_dir = joinpath(run_dir, "output")
    info_path = joinpath(run_dir, "RUN_INFO.toml")
    isfile(info_path) ||
        throw(ArgumentError("resume: $info_path not found; the run has no segment record"))
    info = TOML.parsefile(info_path)
    segments = get(info, "segments", Any[])
    isempty(segments) && throw(ArgumentError("resume: $info_path lists no segment"))

    cfg = load_config(joinpath(run_dir, "config.toml"))
    base_dir === nothing && (base_dir = cfg.config_dir)
    stdout_path = joinpath(out_dir, cfg.postprocess.stdout_file)
    last = segments[end]
    index = Int(last["index"])
    status = _last_status(last, stdout_path)

    # A record left open: refused while the engine runs, else closed from the stdout
    if status == "running"
        if _segment_alive(last, out_dir)
            throw(
                ArgumentError(
                    "resume: segment $index of $(basename(run_dir)) is still running " *
                    "(pid $(last["pid"]) on $(last["host"]))",
                ),
            )
        end
        offset = Int(get(last, "stdout_offset", 0))
        status = _segment_status(stdout_path; offset = offset)
        t_adj = _last_adjust_time(stdout_path; offset = offset)
        _close_segment(
            run_dir,
            index;
            status = status,
            t_end = isnan(t_adj) ? nothing : t_adj,
            fields = Dict{String,Any}("reconciled" => true),
        )
        @info "Segment $index was left open; closed as $(status) from the stdout capture"
    end

    if status == "completed"
        @info "Run $(basename(run_dir)) is complete; nothing to resume"
        return run_dir
    end
    if status == "halted"
        throw(
            ArgumentError(
                "resume: $(basename(run_dir)) was halted by the engine's energy check; a " *
                "restart would need changed accuracy parameters, which is a decision and not " *
                "a resume (restart_simulation continues it explicitly)",
            ),
        )
    end
    if status == "watchdog"
        throw(
            ArgumentError(
                "resume: $(basename(run_dir)) was terminated by the start-up watchdog; it " *
                "would hang again (check the input intervals, or raise " *
                "simulation.startup_timeout)",
            ),
        )
    end

    original_inp = joinpath(out_dir, get(get(info, "run", Dict{String,Any}()), "input_file", ""))
    isfile(original_inp) || throw(
        ArgumentError(
            "resume: the original input recorded in RUN_INFO.toml (run.input_file) is " *
            "missing in $out_dir",
        ),
    )
    tcrit = _tcrit_of(original_inp)

    marker = _resume_dump(out_dir, stdout_path)
    marker === nothing && throw(
        ArgumentError(
            "resume: no complete restart dump in $out_dir; without checkpoint = true a run " *
            "can be resumed only after a stop request",
        ),
    )
    t_d = marker.time_nb
    remaining = tcrit - t_d
    remaining > 0 || throw(
        ArgumentError("resume: the dump $(marker.file) is at t = $t_d, not before TCRIT = $tcrit"),
    )

    discarded = _discard_after(out_dir, stdout_path, marker; segment = index, last_status = status)
    if !isempty(discarded["moved"])
        @info "Resume from $(marker.file) (t = $t_d): $(length(discarded["moved"])) item(s) " *
              "written after that dump moved to $(discarded["dir"])"
        if !isempty(discarded["uncut"])
            @warn "Appended files not cut at the join; they keep what the interrupted segment " *
                  "wrote after t = $(discarded["t_from"]): $(join(discarded["uncut"], ", "))"
        end
    else
        @info "Resume from $(marker.file) (t = $t_d) for $(remaining) N-body time units"
    end

    return _launch_restart(
        cfg,
        run_dir,
        out_dir,
        original_inp,
        marker.file,
        remaining;
        t_start = t_d,
        base_dir = base_dir,
        discarded = isempty(discarded["moved"]) ? nothing : discarded,
    )
end
