# =============================================================================
# Shared utilities
# =============================================================================

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
