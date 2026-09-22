# =============================================================================
# Download, configure, patch, and build Nbody6PPGPU-beijing
# =============================================================================

"""
    setup_nbody6(cfg::Nbody6Config; base_dir = cfg.config_dir)

Orchestrate the full install pipeline:
  clone → configure → HDF5 patch → CUDA setup → build.

All paths are resolved relative to `base_dir` (defaults to the package root).
"""
function setup_nbody6(cfg::Nbody6Config; base_dir::AbstractString = cfg.config_dir)
    install = cfg.install
    build = cfg.build
    src_dir = _resolve_path(base_dir, install.install_dir)

    # ------------------------------------------------------------------
    # 1. Dependency check
    # ------------------------------------------------------------------
    @info "Checking build dependencies..."
    missing_deps = check_dependencies(cfg)
    if !isempty(missing_deps)
        error("Missing dependencies: $(join(missing_deps, ", ")). Install them before proceeding.")
    end
    @info "All dependencies satisfied."

    # ------------------------------------------------------------------
    # 2. Platform checks
    # ------------------------------------------------------------------
    platform = detect_platform()
    @info "Detected platform: $platform"

    if platform == :fedora && build.enable_mpi && build.enable_hdf5
        check_fedora_h5pfc()
    end

    # ------------------------------------------------------------------
    # 3. CUDA auto-detection
    # ------------------------------------------------------------------
    cuda_path = build.cuda_path
    if build.enable_gpu && isempty(cuda_path)
        cuda_path = detect_cuda_path()
        if isempty(cuda_path)
            @warn "GPU enabled but CUDA not found. Build may fail."
        else
            @info "Auto-detected CUDA: $cuda_path"
        end
    end
    # The configure script emits no architecture flag: resolve the targets
    # here and pass them to make as a CUFLAGS override (step 7).
    cuda_archs = build.enable_gpu ? resolve_cuda_arch(build) : String[]
    if !isempty(cuda_archs)
        @info "CUDA target architectures: $(join(cuda_archs, ", "))"
        # The device capability comes from nvidia-smi; nvcc knows only what
        # the toolkit can compile. Fail here, not with an nvcc fatal mid-build.
        _check_cuda_arch_support(
            cuda_archs,
            nvcc_supported_archs(cuda_path),
            nvcc_release(cuda_path),
        )
    end
    # The toolkit may reject the host compiler; find out on a trivial kernel
    # now rather than after the clone, and add the override it needs.
    nvcc_flags = String.(build.nvcc_flags)
    if build.enable_gpu && !isempty(cuda_path)
        append!(nvcc_flags, _nvcc_host_compiler_flags(cuda_path; nvcc_flags = nvcc_flags))
    end

    # ------------------------------------------------------------------
    # 4. Clone / reinstall
    # ------------------------------------------------------------------
    if install.reinstall && isdir(src_dir)
        _is_engine_checkout(src_dir) || _refuse_foreign_directory(src_dir)
        @info "Reinstall requested — removing $src_dir"
        rm(src_dir; recursive = true)
    end

    _ensure_source_tree(src_dir, install)

    # ------------------------------------------------------------------
    # 5. Configure (output suppressed — only errors shown)
    # ------------------------------------------------------------------
    # Build environment: inherit current ENV, then overlay our variables.
    # configure needs it as much as make: the engine's configure looks for
    # nvcc on PATH, and its --with-cuda fallback reuses the cached result of
    # that check, so a toolkit off PATH is found only through PATH.
    build_env = copy(ENV)
    build_env["OMP_STACKSIZE"] = "4096M"
    if build.enable_gpu && !isempty(cuda_path)
        merge!(build_env, cuda_env_vars(cuda_path))
    end

    @info "Running ./configure ..."
    configure_args = _configure_args(build, cuda_path)
    cd(src_dir) do
        _run_quiet(Cmd(`./configure $configure_args`; env = build_env); label = "configure")
    end

    # ------------------------------------------------------------------
    # 6. HDF5 Makefile patch (robust: uses separate include file)
    # ------------------------------------------------------------------
    if build.enable_hdf5
        @info "Patching for HDF5 support..."
        patch_makefile_hdf5(src_dir, build.enable_mpi)
    end

    # ------------------------------------------------------------------
    # 7. Build
    # ------------------------------------------------------------------
    makefile_dir = joinpath(src_dir, "build")
    isdir(makefile_dir) || (makefile_dir = src_dir)

    np = build.nproc > 0 ? build.nproc : nproc_available()

    cd(makefile_dir) do
        if install.clean_build
            @info "Running make clean..."
            try
                _run_quiet(`make clean`; label = "clean")
            catch
                @info "make clean skipped (no prior build artifacts)"
            end
        end
        make_args = ["-j$np"]
        if build.enable_gpu
            cuflags = _cuflags_with_arch(joinpath(makefile_dir, "Makefile"), cuda_archs, nvcc_flags)
            push!(make_args, "CUFLAGS=$cuflags")
        end
        t_build = time()
        @info "Compiling with $np processes..."
        _run_build_with_progress(Cmd(`make $make_args`; env = build_env))
        elapsed = _format_elapsed(time() - t_build)
        @info "Compilation finished ($elapsed)"
    end

    # ------------------------------------------------------------------
    # 8. Locate binary and record the build
    # ------------------------------------------------------------------
    binary = _find_binary(src_dir, cfg.simulation.binary_name, build)
    _write_build_info(
        src_dir,
        cfg,
        configure_args,
        cuda_path,
        cuda_archs,
        binary;
        nvcc_flags = nvcc_flags,
    )
    @info "Build complete. Binary: $binary"
    return binary
end

"""
Directory of the CUDA helper headers shipped with the package (`helper_cuda.h`,
`helper_string.h` from NVIDIA's cuda-samples, tag v13.0). The engine's own
copy under `extra_inc/cuda` dates from 2012 and reads `cudaDeviceProp`
fields (`clockRate`, `computeMode`) that CUDA 13.0 removed, so it no longer
compiles with a CUDA 13 toolkit; this directory precedes it in the `nvcc`
include path.
"""
const _CUDA_HELPER_DIR = joinpath(_PACKAGE_ROOT, "deps", "cuda")

"""
    _cuflags_with_arch(makefile, archs, nvcc_flags = String[];
                       helper_dir = _CUDA_HELPER_DIR) -> String

The `CUFLAGS` value for the `make CUFLAGS=...` command-line override of a
GPU build: the include path of the shipped CUDA helper headers
(`helper_dir`, ahead of the engine's stale copy), then the value
`./configure` wrote into `makefile` (`-O3`, the `CUDA_5` define, the
engine's own include path), the [`cuda_gencode_flags`](@ref) of `archs`
and finally `nvcc_flags` verbatim. Without the architecture flags `nvcc`
compiles for its default target, which the driver JIT-compiles from PTX on
every newer device.
"""
function _cuflags_with_arch(
    makefile::AbstractString,
    archs::AbstractVector{<:AbstractString},
    nvcc_flags::AbstractVector{<:AbstractString} = String[];
    helper_dir::AbstractString = _CUDA_HELPER_DIR,
)::String
    base = "-O3"
    if isfile(makefile)
        for line in eachline(makefile)
            m = match(r"^CUFLAGS\s*=\s*(.*)$", line)
            m === nothing && continue
            base = String(strip(m.captures[1]))
            break
        end
    end
    parts = String["-I $helper_dir", base]
    gencode = cuda_gencode_flags(archs)
    isempty(gencode) || push!(parts, gencode)
    append!(parts, String.(nvcc_flags))
    return join(parts, " ")
end

"""
    _write_build_info(src_dir, cfg, configure_args, cuda_path, cuda_archs, binary;
                      nvcc_flags = cfg.build.nvcc_flags) -> String

Write `BUILD_INFO.toml` next to the binary: date, host, backend commit,
configure arguments, the MPI/GPU/HDF5 switches, the binary name and, for
GPU builds, the CUDA path, the compiled architectures, the `nvcc` release,
the `nvcc` options in effect (`nvcc_flags`: the configured ones plus any
host-compiler override the build added) and the helper-header directory.
The launcher copies the file into every run directory and merges it into
`RUN_INFO.toml` as the `[build]` table, so each result records the build
that produced it. Returns the path.
"""
function _write_build_info(
    src_dir::AbstractString,
    cfg::Nbody6Config,
    configure_args::AbstractVector{<:AbstractString},
    cuda_path::AbstractString,
    cuda_archs::AbstractVector{<:AbstractString},
    binary::AbstractString;
    nvcc_flags::AbstractVector{<:AbstractString} = cfg.build.nvcc_flags,
)::String
    build = cfg.build
    d = Dict{String,Any}(
        "date" => Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"),
        "host" => gethostname(),
        "backend_commit" => _git_commit(src_dir),
        "configure_args" => String.(configure_args),
        "enable_mpi" => build.enable_mpi,
        "enable_gpu" => build.enable_gpu,
        "enable_hdf5" => build.enable_hdf5,
        "binary" => basename(binary),
    )
    if build.enable_gpu
        d["cuda_path"] = String(cuda_path)
        d["cuda_arch"] = String.(cuda_archs)
        d["nvcc_release"] = nvcc_release(cuda_path)
        d["nvcc_flags"] = String.(nvcc_flags)
        d["cuda_helper_dir"] = _CUDA_HELPER_DIR
    end
    path = joinpath(dirname(binary), "BUILD_INFO.toml")
    open(path, "w") do io
        TOML.print(io, d)
    end
    return path
end

# ---------------------------------------------------------------------------
# Robust HDF5 Makefile patching via separate include file
# ---------------------------------------------------------------------------

"""
Write HDF5 flags to `hdf5_flags.mk` in the source root and inject an
`-include` directive into `build/Makefile`.  The flags file persists
across `./configure` re-runs; only the one-line include injection needs
to be re-applied.
"""
function patch_makefile_hdf5(src_dir::AbstractString, use_mpi::Bool)
    platform = detect_platform()
    cppflags, libs = detect_hdf5_flags(platform, use_mpi)

    # --- Write persistent hdf5_flags.mk ---
    flags_file = joinpath(src_dir, "hdf5_flags.mk")
    open(flags_file, "w") do io
        println(io, "# HDF5 flags — generated by Nbody6Dynamics.jl")
        println(io, "# Re-run setup to regenerate if paths change.")
        println(io, "CPPFLAGS += -DCONFIG_HDF5 $cppflags")
        println(io, "LIBS     += $libs")
    end
    @info "HDF5 flags written to $flags_file"

    # --- Inject include into build/Makefile ---
    makefile = joinpath(src_dir, "build", "Makefile")
    if !isfile(makefile)
        @warn "build/Makefile not found — skipping include injection"
        return
    end

    include_line = "-include ../hdf5_flags.mk"
    content = read(makefile, String)
    if !occursin(include_line, content)
        open(makefile, "a") do io
            println(io)
            println(io, "# --- Nbody6Dynamics.jl HDF5 include (survives ./configure) ---")
            println(io, include_line)
        end
        @info "Injected HDF5 include into build/Makefile"
    else
        @info "HDF5 include already present in build/Makefile"
    end
end

# ---------------------------------------------------------------------------
# Binary locator
# ---------------------------------------------------------------------------

"""
    _ensure_source_tree(src_dir, install)

Leave `src_dir` holding a usable engine checkout: clone it when absent or
empty, keep it as it is when it carries a `configure` script, and repair it
when it does not. A git directory without `configure` is a truncated
checkout — an interrupted `git clone` (a killed run, a host reset) leaves
exactly that, and the build would fail later with a bare `ENOENT` on
`./configure`. The repair restores the tracked files with
`git checkout --force`, which keeps untracked work; when git cannot restore
them the directory is set aside under a new name and a fresh clone takes its
place.

Nothing is ever deleted here. `src_dir` comes from `install.install_dir`,
and a mistyped value may name a directory that holds something else
entirely: a non-empty directory that is not a git checkout is refused with
an `ArgumentError`.
"""
function _ensure_source_tree(src_dir::AbstractString, install)
    if !isdir(src_dir) || isempty(readdir(src_dir))
        _clone_source(src_dir, install)
        return nothing
    end
    if _source_tree_ready(src_dir)
        @info "Source directory already exists: $src_dir (reference left as is)"
        return nothing
    end
    _is_engine_checkout(src_dir) || _refuse_foreign_directory(src_dir)
    @warn "Engine source tree at $src_dir has no configure script: the checkout is " *
          "incomplete (an interrupted clone leaves this). Restoring it."
    ref = isempty(install.ref) ? "HEAD" : install.ref
    try
        _run_quiet(`git -C $src_dir checkout --force $ref`; label = "checkout")
    catch e
        e isa Union{ProcessFailedException,Base.IOError} || rethrow()
        @debug "restoring the checkout failed; cloning again" exception = e
    end
    _source_tree_ready(src_dir) && return nothing
    aside = _set_aside(src_dir)
    @warn "The checkout cannot be restored; kept as $aside and cloning again."
    _clone_source(src_dir, install)
    return nothing
end

"""`true` when `dir` is recognisably an engine checkout: it holds `configure` or its own `.git`."""
_is_engine_checkout(dir::AbstractString)::Bool =
    _source_tree_ready(dir) || ispath(joinpath(dir, ".git"))

function _refuse_foreign_directory(dir::AbstractString)
    throw(
        ArgumentError(
            "install.install_dir resolves to $dir, which exists, is not empty and is not an " *
            "engine checkout (no configure script, no .git); refusing to replace it. Point " *
            "install_dir at another location, or move that directory away",
        ),
    )
end

"""
    _set_aside(dir) -> String

Rename `dir` to `<dir>.incomplete-<yyyymmdd_HHMMSS>` (suffixed `-1`, `-2`, …
while that name is taken) and return the new path.
"""
function _set_aside(dir::AbstractString)::String
    base = rstrip(abspath(dir), '/') * ".incomplete-" * Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    target, n = base, 0
    while ispath(target)
        n += 1
        target = "$base-$n"
    end
    mv(dir, target)
    return target
end

"""`true` when `src_dir` holds the engine's `configure` script, which every complete checkout has (it is tracked upstream)."""
_source_tree_ready(src_dir::AbstractString)::Bool = isfile(joinpath(src_dir, "configure"))

"""Clone `install.source_url` into `src_dir` and check out `install.ref` when one is configured."""
function _clone_source(src_dir::AbstractString, install)
    @info "Cloning $(install.source_url) → $src_dir"
    run(`git clone $(install.source_url) $src_dir`)
    if !isempty(install.ref)
        @info "Checking out backend reference $(install.ref)"
        _run_quiet(`git -C $src_dir checkout --quiet $(install.ref)`; label = "checkout")
    end
    return nothing
end

"""
    _configure_args(build, cuda_path) -> Vector{String}

The engine's `configure` arguments: the configured flags, `--disable-mpi`
and `--disable-gpu` unless enabled, and, for a GPU build with a detected
toolkit, `--with-cuda=<cuda_path>` unless the flags already carry one.
"""
function _configure_args(build, cuda_path::AbstractString)::Vector{String}
    args = String.(build.configure_flags)
    build.enable_mpi || push!(args, "--disable-mpi")
    build.enable_gpu || push!(args, "--disable-gpu")
    if build.enable_gpu && !isempty(cuda_path) && !any(startswith("--with-cuda"), args)
        push!(args, "--with-cuda=$cuda_path")
    end
    return args
end

"""
Run a command with stdout/stderr redirected to a temp file.
On failure, dump the captured output as a warning and re-throw.
On success, the temp file is deleted.
"""
function _run_quiet(cmd::Cmd; label::AbstractString = "command")
    logfile = tempname()
    try
        open(logfile, "w") do f
            run(pipeline(cmd; stdout = f, stderr = f))
        end
    catch e
        if isfile(logfile)
            output = read(logfile, String)
            !isempty(output) && @warn "$label failed. Output:\n$output"
        end
        rethrow()
    finally
        isfile(logfile) && rm(logfile; force = true)
    end
end

"""
Run the `make` build command, parsing its merged stdout+stderr for compilation
progress.  Counts `.o` targets as they complete and displays a live progress
counter.  Warnings are accumulated and summarised at the end.  Errors are
shown immediately.
"""
function _run_build_with_progress(cmd::Cmd)
    # Merge stderr into stdout so gfortran warnings don't leak to the terminal
    merged_cmd = pipeline(cmd; stderr = stdout)
    proc = open(merged_cmd; read = true, write = false)
    n_compiled = 0
    n_warnings = 0
    linking = false

    for line in eachline(proc)
        if occursin(r"-o\s+\S+\.o\b", line)
            n_compiled += 1
            print(stderr, "\r\e[K  ⚙  Compiled $n_compiled objects...")
        elseif occursin(r"-o\s+\S*(nbody6|\.exe|\.avx)", line)
            linking = true
            print(stderr, "\r\e[K  🔗  Linking ($n_compiled objects compiled)...")
        elseif occursin(r"[Ww]arning:", line)
            n_warnings += 1
        elseif occursin(r"[Ee]rror:", line)
            print(stderr, "\r\e[K")
            @error "Build error: $line"
        end
    end

    # Clear progress line
    print(stderr, "\r\e[K")

    wait(proc)
    if !success(proc)
        error("Build failed with exit code $(proc.exitcode)")
    end

    if n_warnings > 0
        @info "  $n_compiled objects compiled ($n_warnings compiler warnings)"
    else
        @info "  $n_compiled objects compiled (0 warnings)"
    end
end

"""
    _find_binary(src_dir, binary_name, build::BuildConfig) -> String

Path of the compiled engine under `src_dir/build`. The build system names
the binary after its configure options (`nbody6++.avx`, `nbody6++.avx.gpu`,
`nbody6++.avx.mpi.gpu`, …) and several variants may coexist, so the
candidate is chosen by its suffix tags: `.gpu` present exactly when
`build.enable_gpu`, `.mpi` present exactly when `build.enable_mpi`. Among
several matches the most recently modified wins. Throws when no variant
matches, listing the files found.
"""
function _find_binary(
    src_dir::AbstractString,
    binary_name::AbstractString,
    build::BuildConfig,
)::String
    build_dir = joinpath(src_dir, "build")
    isdir(build_dir) || error("Build directory not found: $build_dir. Run the install phase first.")
    found = String[]
    matches = String[]
    for f in sort(readdir(build_dir))
        startswith(f, binary_name) || continue
        path = joinpath(build_dir, f)
        isfile(path) || continue
        push!(found, f)
        tags = split(chop(f; head = length(binary_name), tail = 0), '.'; keepempty = false)
        ("gpu" in tags) == build.enable_gpu || continue
        ("mpi" in tags) == build.enable_mpi || continue
        push!(matches, path)
    end
    isempty(matches) && error(
        "No $(binary_name) binary with gpu = $(build.enable_gpu), mpi = $(build.enable_mpi) " *
        "in $build_dir " *
        (
            isempty(found) ? "(no binary at all; check the build output)" :
            "(found: $(join(found, ", ")))"
        ),
    )
    return matches[argmax(mtime.(matches))]
end
