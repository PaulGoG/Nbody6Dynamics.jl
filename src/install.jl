# =============================================================================
# Download, configure, patch, and build Nbody6PPGPU-beijing
# =============================================================================

"""
    setup_nbody6(cfg::Nbody6Config; base_dir = _PROJECT_ROOT)

Orchestrate the full install pipeline:
  clone → configure → HDF5 patch → CUDA setup → build.

All paths are resolved relative to `base_dir` (defaults to the package root).
"""
function setup_nbody6(cfg::Nbody6Config; base_dir::AbstractString = _PROJECT_ROOT)
    install = cfg.install
    build = cfg.build
    src_dir = joinpath(base_dir, install.install_dir)

    # ------------------------------------------------------------------
    # 1. Dependency check
    # ------------------------------------------------------------------
    @info "Checking build dependencies..."
    missing = check_dependencies(cfg)
    if !isempty(missing)
        error("Missing dependencies: $(join(missing, ", ")). Install them before proceeding.")
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

    # ------------------------------------------------------------------
    # 4. Clone / reinstall
    # ------------------------------------------------------------------
    if install.reinstall && isdir(src_dir)
        @info "Reinstall requested — removing $src_dir"
        rm(src_dir; recursive = true)
    end

    if !isdir(src_dir)
        @info "Cloning $(install.source_url) → $src_dir"
        run(`git clone $(install.source_url) $src_dir`)
    else
        @info "Source directory already exists: $src_dir"
    end

    # ------------------------------------------------------------------
    # 5. Configure (output suppressed — only errors shown)
    # ------------------------------------------------------------------
    @info "Running ./configure ..."
    configure_args = copy(build.configure_flags)
    build.enable_mpi || push!(configure_args, "--disable-mpi")
    build.enable_gpu || push!(configure_args, "--disable-gpu")

    cd(src_dir) do
        _run_quiet(`./configure $configure_args`; label = "configure")
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

    # Build environment: inherit current ENV, then overlay our variables
    build_env = copy(ENV)
    build_env["OMP_STACKSIZE"] = "4096M"
    if build.enable_gpu && !isempty(cuda_path)
        merge!(build_env, cuda_env_vars(cuda_path))
    end

    cd(makefile_dir) do
        if install.clean_build
            @info "Running make clean..."
            try
                _run_quiet(`make clean`; label = "clean")
            catch
                @info "make clean skipped (no prior build artifacts)"
            end
        end
        t_build = time()
        @info "Compiling with $np processes..."
        _run_build_with_progress(Cmd(`make -j$np`; env = build_env))
        elapsed = _format_elapsed(time() - t_build)
        @info "Compilation finished ($elapsed)"
    end

    # ------------------------------------------------------------------
    # 8. Locate binary
    # ------------------------------------------------------------------
    binary = _find_binary(src_dir, cfg.simulation.binary_name)
    @info "Build complete. Binary: $binary"
    return binary
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
Search common locations for the compiled binary.
"""
function _find_binary(src_dir::AbstractString, binary_name::AbstractString)::String
    build_dir = joinpath(src_dir, "build")
    candidates = [
        joinpath(build_dir, binary_name),
        joinpath(build_dir, "nbody6++.gpu"),
        joinpath(build_dir, "nbody6++"),
    ]

    if isdir(build_dir)
        for f in readdir(build_dir)
            startswith(f, "nbody6++") && push!(candidates, joinpath(build_dir, f))
        end
    end

    for path in candidates
        isfile(path) && return path
    end

    error("Could not find compiled binary in $build_dir. Check build output for errors.")
end
