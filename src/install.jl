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
        if !isempty(install.ref)
            @info "Checking out backend reference $(install.ref)"
            _run_quiet(`git -C $src_dir checkout --quiet $(install.ref)`; label = "checkout")
        end
    else
        @info "Source directory already exists: $src_dir (reference left as is)"
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
        make_args = ["-j$np"]
        if !isempty(cuda_archs)
            cuflags = _cuflags_with_arch(joinpath(makefile_dir, "Makefile"), cuda_archs)
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
    _write_build_info(src_dir, cfg, configure_args, cuda_path, cuda_archs, binary)
    @info "Build complete. Binary: $binary"
    return binary
end

"""
    _cuflags_with_arch(makefile, archs) -> String

The `CUFLAGS` value `./configure` wrote into `makefile` (`-O3`, the
`CUDA_5` define, the `helper_cuda.h` include path) extended with the
[`cuda_gencode_flags`](@ref) of `archs`, for a `make CUFLAGS=...`
command-line override. Without the override `nvcc` compiles for its
default target, which the driver JIT-compiles from PTX on every newer
device.
"""
function _cuflags_with_arch(
    makefile::AbstractString,
    archs::AbstractVector{<:AbstractString},
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
    gencode = cuda_gencode_flags(archs)
    return isempty(gencode) ? base : base * " " * gencode
end

"""
    _write_build_info(src_dir, cfg, configure_args, cuda_path, cuda_archs, binary) -> String

Write `BUILD_INFO.toml` next to the binary: date, host, backend commit,
configure arguments, the MPI/GPU/HDF5 switches, the binary name and, for
GPU builds, the CUDA path, the compiled architectures and the `nvcc`
release. The launcher copies the file into every run directory and merges
it into `RUN_INFO.toml` as the `[build]` table, so each result records the
build that produced it. Returns the path.
"""
function _write_build_info(
    src_dir::AbstractString,
    cfg::Nbody6Config,
    configure_args::AbstractVector{<:AbstractString},
    cuda_path::AbstractString,
    cuda_archs::AbstractVector{<:AbstractString},
    binary::AbstractString,
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
