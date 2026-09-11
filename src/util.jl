# =============================================================================
# Shared utilities
# =============================================================================

"""
    _git_commit(dir) -> String

Short commit hash of the git repository at `dir`, with a `-dirty` suffix
when the working tree has uncommitted changes; `"unknown"` when `dir` is not
a repository or git is unavailable. Used to stamp run provenance
(tagsave-equivalent) into RUN_INFO.toml and merger_ic.toml.
"""
function _git_commit(dir::AbstractString)::String
    try
        h = strip(read(pipeline(`git -C $dir rev-parse --short HEAD`; stderr = devnull), String))
        dirty =
            !isempty(
                strip(read(pipeline(`git -C $dir status --porcelain`; stderr = devnull), String)),
            )
        return dirty ? h * "-dirty" : h
    catch
        return "unknown"
    end
end

"""
    _source_stamp(dir) -> String

Provenance identity of the source tree at `dir`: the git commit
([`_git_commit`](@ref)) when one is available, otherwise the `version` of a
`Project.toml` there suffixed with `+nogit`, otherwise `"unknown"`.

A tree deployed to a compute host by file copy carries no `.git`, so the
commit is unobtainable exactly where provenance matters most — every run of
the 2026-09-11 fleet campaign recorded `package_commit = "unknown"`. The
version fallback ties such a run to a release at least.
"""
function _source_stamp(dir::AbstractString)::String
    commit = _git_commit(dir)
    commit == "unknown" || return commit
    project = joinpath(dir, "Project.toml")
    isfile(project) || return "unknown"
    version = try
        get(TOML.parsefile(project), "version", nothing)
    catch
        nothing
    end
    return version === nothing ? "unknown" : "v$(version)+nogit"
end

"""
    _with_run_log(f, run_dir) -> result of f()

Run `f()` with the current logger teed to a plain-text, ANSI-free
`nbody6dynamics.log` inside `run_dir` (structured logging to file for
long runs). Console behaviour is unchanged; the file sink records
timestamped Info+ records and is appended to across pipeline phases.
"""
function _with_run_log(f, run_dir::AbstractString)
    mkpath(run_dir)
    open(joinpath(run_dir, "nbody6dynamics.log"), "a") do io
        file_logger = FormatLogger(io) do fio, args
            println(
                fio,
                Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
                " [",
                args.level,
                "] ",
                args.message,
            )
        end
        tee = TeeLogger(Logging.current_logger(), MinLevelLogger(file_logger, Logging.Info))
        return Logging.with_logger(f, tee)
    end
end

"""
    export_for_paper(paths, dest_dir; run_dir = nothing) -> Vector{String}

Copy finished figures to a manuscript figures directory with provenance
attached: each file is copied as `<run_id>__<name>` and a TOML sidecar
`<run_id>__<name>.provenance.toml` records the producing run, the package
and backend commits (from the run's RUN_INFO.toml when `run_dir` is given),
and the export date. `dest_dir` is caller-supplied — the package never
hardcodes paths outside its own tree. Existing files are backed up, never
overwritten. Returns the destination paths.
"""
function export_for_paper(
    paths::AbstractVector{<:AbstractString},
    dest_dir::AbstractString;
    run_dir::Union{Nothing,AbstractString} = nothing,
)
    mkpath(dest_dir)
    run_id = run_dir === nothing ? "unattributed" : basename(abspath(run_dir))
    # Provenance from the run's RUN_INFO.toml when available
    commits = Dict{String,String}()
    if run_dir !== nothing
        info_path = joinpath(run_dir, "RUN_INFO.toml")
        if isfile(info_path)
            prov = get(TOML.parsefile(info_path), "provenance", Dict{String,Any}())
            for (short, key) in (("commit", "package_commit"), ("backend", "backend_commit"))
                haskey(prov, key) && (commits[short] = String(prov[key]))
            end
        end
    end

    out = String[]
    for src in paths
        isfile(src) || error("export_for_paper: no such file: $src")
        dest = joinpath(dest_dir, run_id * "__" * basename(src))
        _backup_existing(dest)
        cp(src, dest)
        sidecar = dest * ".provenance.toml"
        _backup_existing(sidecar)
        open(sidecar, "w") do io
            TOML.print(
                io,
                Dict(
                    "run_id" => run_id,
                    "source" => abspath(src),
                    "exported_at" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
                    "package_commit" => get(commits, "commit", "unknown"),
                    "backend_commit" => get(commits, "backend", "unknown"),
                ),
            )
        end
        push!(out, dest)
        @info "Exported for paper: $(basename(dest))"
    end
    return out
end

"""
    _backup_existing(path) -> Union{Nothing,String}

Never-overwrite protection for generated results (DrWatson-`safesave` style):
if `path` exists, move it to the first free `<stem>#<k><ext>` sibling before
the caller writes the new file. Returns the backup path, or `nothing` if no
file existed.
"""
function _backup_existing(path::AbstractString)::Union{Nothing,String}
    isfile(path) || return nothing
    base, ext = splitext(path)
    k = 1
    while isfile("$(base)#$(k)$(ext)")
        k += 1
    end
    backup = "$(base)#$(k)$(ext)"
    mv(path, backup)
    @info "Backed up existing $(basename(path)) → $(basename(backup))"
    return backup
end
