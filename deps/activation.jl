# Activation of the package's environments, shared by the activate.jl scripts
# of the root, scripts/, bench/ and docs/. Manifests are not tracked, so a
# checkout that pulls a new dependency keeps a manifest resolved for the old
# project; such a manifest is resolved again before instantiation instead of
# failing at precompilation. Stale means: resolved from another version of the
# environment's own project, or of the project of a package the environment
# takes by path ([sources]).
using Pkg
using TOML

function activate_environment(dir::AbstractString)
    Pkg.activate(dir; io = devnull)
    _manifest_stale(dir) && Pkg.resolve(; io = devnull)
    Pkg.instantiate(; io = devnull)
    return nothing
end

function _manifest_stale(dir::AbstractString)::Bool
    manifest = joinpath(dir, "Manifest.toml")
    isfile(manifest) || return false
    Pkg.is_manifest_current(dir) === false && return true
    project = TOML.parsefile(joinpath(dir, "Project.toml"))
    entries = get(TOML.parsefile(manifest), "deps", Dict{String,Any}())
    for (name, source) in get(project, "sources", Dict{String,Any}())
        haskey(source, "path") || continue
        sourced = TOML.parsefile(joinpath(dir, source["path"], "Project.toml"))
        wanted = Set(keys(get(sourced, "deps", Dict{String,Any}())))
        haskey(entries, name) || return true
        recorded = get(first(entries[name]), "deps", String[])
        Set(recorded isa AbstractDict ? keys(recorded) : recorded) == wanted || return true
    end
    return false
end
