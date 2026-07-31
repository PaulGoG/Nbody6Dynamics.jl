# =============================================================================
# Shared utilities
# =============================================================================

"""
    _git_commit(dir) -> String

Short commit hash of the git repository at `dir`, with a `-dirty` suffix
when the working tree has uncommitted changes; `"unknown"` when `dir` is not
a repository or git is unavailable. Used to stamp run provenance
(tagsave-equivalent) into RUN_INFO.txt and merger_ic.toml.
"""
function _git_commit(dir::AbstractString)::String
    try
        h = strip(read(`git -C $dir rev-parse --short HEAD`, String))
        dirty = !isempty(strip(read(`git -C $dir status --porcelain`, String)))
        return dirty ? h * "-dirty" : h
    catch
        return "unknown"
    end
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
