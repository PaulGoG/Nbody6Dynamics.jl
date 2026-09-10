# =============================================================================
# Platform detection, dependency checking, and CUDA/HDF5 tooling
# =============================================================================

"""
    detect_platform() -> Symbol

Return `:fedora`, `:ubuntu`, `:debian`, or `:unknown` based on /etc/os-release.
"""
function detect_platform()::Symbol
    release_file = "/etc/os-release"
    !isfile(release_file) && return :unknown
    content = lowercase(read(release_file, String))
    occursin("fedora", content) && return :fedora
    occursin("ubuntu", content) && return :ubuntu
    occursin("debian", content) && return :debian
    return :unknown
end

"""
    check_command(cmd::String) -> Bool

Return `true` if `cmd` is found on PATH.
"""
function check_command(cmd::String)::Bool
    try
        success(pipeline(`which $cmd`; stdout = devnull, stderr = devnull))
    catch
        false
    end
end

"""
    check_dependencies(cfg::Nbody6Config) -> Vector{String}

Return a list of missing dependencies required for building Nbody6++:
`git`, the GNU toolchain and `make` on `PATH`, the MPI wrappers when
`build.enable_mpi`, and for GPU builds `nvcc` either on `PATH` or under
`<cuda_path>/bin` of the configured or auto-detected toolkit, whose `bin`
the build exports itself.
"""
function check_dependencies(cfg::Nbody6Config)::Vector{String}
    missing_deps = String[]

    required = ["git", "gcc", "g++", "gfortran", "make"]
    if cfg.build.enable_mpi
        append!(required, ["mpicc", "mpif90", "mpirun"])
    end

    for cmd in required
        check_command(cmd) || push!(missing_deps, cmd)
    end

    if cfg.build.enable_gpu
        cuda_path = isempty(cfg.build.cuda_path) ? detect_cuda_path() : cfg.build.cuda_path
        nvcc_present =
            check_command("nvcc") ||
            (!isempty(cuda_path) && isfile(joinpath(cuda_path, "bin", "nvcc")))
        nvcc_present || push!(missing_deps, "nvcc")
    end

    return missing_deps
end

# ---------------------------------------------------------------------------
# CUDA detection
# ---------------------------------------------------------------------------

"""
    detect_cuda_path() -> String

Auto-detect CUDA installation path.  Returns empty string if not found.
"""
function detect_cuda_path()::String
    # 1. Environment variables
    for var in ("CUDA_HOME", "CUDA_PATH", "CUDA_ROOT")
        p = get(ENV, var, "")
        !isempty(p) && isdir(p) && return p
    end

    # 2. Standard locations (newest first)
    for path in [
        "/usr/local/cuda",
        "/usr/local/cuda-12",
        "/usr/local/cuda-11.8",
        "/usr/local/cuda-11",
        "/usr/lib64/cuda",
        "/opt/cuda",
    ]
        isdir(path) && return path
    end

    # 3. Locate nvcc
    try
        nvcc = strip(read(pipeline(`which nvcc`; stderr = devnull), String))
        return dirname(dirname(nvcc))  # nvcc lives in CUDA_HOME/bin/
    catch e
        e isa Union{ProcessFailedException,Base.IOError} || rethrow()
        @debug "nvcc not found on PATH" exception = e
    end

    return ""
end

"""
    cuda_env_vars(cuda_path::String) -> Dict{String,String}

Build environment variable overrides for a CUDA-enabled build/run.
"""
function cuda_env_vars(cuda_path::String)::Dict{String,String}
    env = Dict{String,String}()
    isempty(cuda_path) && return env

    bin_dir = joinpath(cuda_path, "bin")
    lib_dir = joinpath(cuda_path, "lib64")

    env["CUDA_HOME"] = cuda_path
    env["PATH"] = "$bin_dir:" * get(ENV, "PATH", "")
    env["LD_LIBRARY_PATH"] = "$lib_dir:" * get(ENV, "LD_LIBRARY_PATH", "")

    return env
end

# ---------------------------------------------------------------------------
# CUDA target architectures
# ---------------------------------------------------------------------------

"""
Largest number of CUDA devices one engine process drives: `MAX_GPU` in the
backend's `gpunb.velocity.cu`. The j-particle range is split over the
devices of `GPU_LIST` (or every visible device), one OpenMP thread each.
"""
const _MAX_GPU_PER_PROCESS = 4

"""Accepted form of a `build.cuda_arch` entry: `sm_<major><minor>` (`sm_90`, `sm_120`)."""
const _CUDA_ARCH_PATTERN = r"^sm_[1-9][0-9]{1,2}$"

"""
    cuda_arch_from_compute_cap(cap) -> String

CUDA architecture name of a compute capability as `nvidia-smi` reports it:
`"9.0"` → `"sm_90"` (Hopper: H100, H200), `"12.0"` → `"sm_120"` (consumer
Blackwell: RTX 50 series), `"8.6"` → `"sm_86"`. Throws an `ArgumentError`
for anything but `<major>.<minor>`.
"""
function cuda_arch_from_compute_cap(cap::AbstractString)::String
    m = match(r"^\s*(\d+)\.(\d+)\s*$", cap)
    m === nothing &&
        throw(ArgumentError("compute capability must read <major>.<minor>; got \"$cap\""))
    return "sm_" * m.captures[1] * m.captures[2]
end

"""
    _parse_compute_caps(output) -> Vector{String}

Distinct compute capabilities, in order of first appearance, from the
lines of `nvidia-smi --query-gpu=compute_cap --format=csv,noheader`;
blank lines and `[N/A]` entries are skipped.
"""
function _parse_compute_caps(output::AbstractString)::Vector{String}
    caps = String[]
    for line in eachline(IOBuffer(String(output)))
        cap = strip(line)
        (isempty(cap) || startswith(cap, "[")) && continue
        cap in caps || push!(caps, String(cap))
    end
    return caps
end

"""
    detect_compute_capabilities() -> Vector{String}

Distinct compute capabilities of the visible NVIDIA devices (`"9.0"`,
`"12.0"`, …) from `nvidia-smi`; empty when the tool or a device is absent.
"""
function detect_compute_capabilities()::Vector{String}
    output = try
        read(`nvidia-smi --query-gpu=compute_cap --format=csv,noheader`, String)
    catch
        return String[]
    end
    return _parse_compute_caps(output)
end

"""Numeric part of an architecture name (`"sm_120"` → 120)."""
_arch_number(arch::AbstractString)::Int = parse(Int, chop(arch; head = 3, tail = 0))

"""
    cuda_gencode_flags(archs) -> String

`nvcc` code-generation flags for the architectures `archs` (`"sm_90"`,
`"sm_120"`, …): native code for each and PTX for the highest, so the
binary also runs, JIT-compiled, on newer devices. Empty for an empty list;
throws an `ArgumentError` for a malformed name.
"""
function cuda_gencode_flags(archs::AbstractVector{<:AbstractString})::String
    isempty(archs) && return ""
    flags = String[]
    for arch in archs
        occursin(_CUDA_ARCH_PATTERN, arch) ||
            throw(ArgumentError("not a CUDA architecture name: \"$arch\""))
        cc = _arch_number(arch)
        push!(flags, "-gencode arch=compute_$(cc),code=sm_$(cc)")
    end
    cc_max = maximum(_arch_number, archs)
    push!(flags, "-gencode arch=compute_$(cc_max),code=compute_$(cc_max)")
    return join(flags, " ")
end

"""
    resolve_cuda_arch(build::BuildConfig) -> Vector{String}

Architectures to compile the GPU kernels for: `build.cuda_arch` when set,
otherwise the compute capabilities of the visible devices
([`detect_compute_capabilities`](@ref)). Empty, with a warning, when
neither is available; the build then keeps the `nvcc` default target and
the driver JIT-compiles its PTX at the first launch.
"""
function resolve_cuda_arch(build::BuildConfig)::Vector{String}
    isempty(build.cuda_arch) || return copy(build.cuda_arch)
    caps = detect_compute_capabilities()
    if isempty(caps)
        @warn "No NVIDIA device visible to nvidia-smi and build.cuda_arch is empty: the GPU " *
              "kernels keep the nvcc default target (PTX JIT at the first launch). Set " *
              "build.cuda_arch, e.g. [\"sm_90\"], for native code."
        return String[]
    end
    return cuda_arch_from_compute_cap.(caps)
end

"""
    nvcc_release(cuda_path = "") -> String

Release of the `nvcc` on `PATH`, or under `cuda_path/bin` when given
(`"12.8"`); empty when the compiler cannot be run.
"""
function nvcc_release(cuda_path::AbstractString = "")::String
    nvcc = isempty(cuda_path) ? "nvcc" : joinpath(cuda_path, "bin", "nvcc")
    output = try
        read(`$nvcc --version`, String)
    catch
        return ""
    end
    return _parse_nvcc_release(output)
end

"""Release number in the banner of `nvcc --version` (`release 12.8, V12.8.93` → `"12.8"`); empty when absent."""
function _parse_nvcc_release(output::AbstractString)::String
    m = match(r"release\s+(\d+\.\d+)", output)
    return m === nothing ? "" : String(m.captures[1])
end

"""
    nvcc_supported_archs(cuda_path = "") -> Vector{String}

Architectures the installed `nvcc` can compile for, as `sm_<cc>` names
from `nvcc --list-gpu-arch` (`compute_50 … compute_120`); empty when the
compiler cannot be run. The device side of the pair is
[`detect_compute_capabilities`](@ref) (from `nvidia-smi`, which `nvcc`
cannot replace: the compiler knows the toolkit, not the hardware); the
build checks the targets against this list before compiling.
"""
function nvcc_supported_archs(cuda_path::AbstractString = "")::Vector{String}
    nvcc = isempty(cuda_path) ? "nvcc" : joinpath(cuda_path, "bin", "nvcc")
    output = try
        read(`$nvcc --list-gpu-arch`, String)
    catch
        return String[]
    end
    return _parse_nvcc_arch_list(output)
end

"""`compute_<cc>` lines of `nvcc --list-gpu-arch` as distinct `sm_<cc>` names, in order."""
function _parse_nvcc_arch_list(output::AbstractString)::Vector{String}
    archs = String[]
    for line in eachline(IOBuffer(String(output)))
        m = match(r"^\s*compute_(\d+)\s*$", line)
        m === nothing && continue
        arch = "sm_" * m.captures[1]
        arch in archs || push!(archs, arch)
    end
    return archs
end

"""
    _check_cuda_arch_support(archs, supported, release)

Throw an `ErrorException` naming the toolkit release and the architectures
it lacks when any of `archs` is absent from `supported`
([`nvcc_supported_archs`](@ref)). A no-op when `supported` is empty (the
compiler could not be queried), leaving the failure to the build itself.
"""
function _check_cuda_arch_support(
    archs::AbstractVector{<:AbstractString},
    supported::AbstractVector{<:AbstractString},
    release::AbstractString,
)
    isempty(supported) && return nothing
    lacking = filter(a -> !(a in supported), archs)
    isempty(lacking) && return nothing
    error(
        "CUDA toolkit " *
        (isempty(release) ? "(unknown release)" : release) *
        " cannot compile for $(join(lacking, ", ")); it supports $(join(supported, ", ")). " *
        "Install a newer toolkit (sm_120 needs CUDA ≥ 12.8, sm_90 CUDA ≥ 11.8) and set " *
        "build.cuda_path, or drop the architecture from build.cuda_arch.",
    )
end

"""
    _nvcc_host_compiler_flags(cuda_path = ""; nvcc_flags = String[]) -> Vector{String}

Compile a trivial kernel with `nvcc` (under `cuda_path/bin` when given)
and `nvcc_flags` to learn whether the toolkit accepts the host compiler.
Returns the options the build must add: empty when the compilation
succeeds or `nvcc` cannot be run at all (the build reports that itself);
`["-allow-unsupported-compiler"]`, with a warning quoting the compiler's
message, when `nvcc` rejects the host compiler version and accepts it
with the override. Any other failure throws with the `nvcc` output, before
anything is cloned or built.
"""
function _nvcc_host_compiler_flags(
    cuda_path::AbstractString = "";
    nvcc_flags::AbstractVector{<:AbstractString} = String[],
    candidates::AbstractVector{<:AbstractString} = _host_compiler_candidates(),
)::Vector{String}
    nvcc = isempty(cuda_path) ? "nvcc" : joinpath(cuda_path, "bin", "nvcc")
    override = "-allow-unsupported-compiler"
    # A -ccbin in the configured flags settles the host compiler; probe as is.
    pinned = "-ccbin" in nvcc_flags || any(startswith("-ccbin="), nvcc_flags)
    mktempdir() do dir
        src = joinpath(dir, "probe.cu")
        obj = joinpath(dir, "probe.o")
        write(
            src,
            "__global__ void probe(float *x) { x[threadIdx.x] = 0.0f; }\nint main() { return 0; }\n",
        )
        attempt(extra) = _run_capture(
            `$nvcc -c $src -o $obj $(String.(nvcc_flags)) $(String.(extra))`,
            joinpath(dir, "probe.log"),
        )
        ok, output = try
            attempt(String[])
        catch e
            e isa Base.IOError || rethrow()
            @debug "nvcc host-compiler probe skipped: nvcc cannot be run" nvcc exception = e
            return String[]
        end
        ok && return String[]
        pinned && error(
            "nvcc cannot compile a trivial kernel with the configured host compiler " *
            "(build.nvcc_flags = $(String.(nvcc_flags))); output:\n" *
            _output_excerpt(output),
        )
        unsupported = _unsupported_host_compiler(output)
        tried = String[]
        if unsupported
            ok_override, _ = attempt([override])
            if ok_override
                @warn "nvcc rejects the host compiler version; building with $override. Set " *
                      "build.nvcc_flags = [\"-ccbin\", \"<older g++>\"] instead if the build misbehaves." nvcc_message =
                    _output_excerpt(output, 3, 0)
                return [override]
            end
            push!(tried, override)
        end
        for cc in candidates
            ok_cc, _ = attempt(["-ccbin", cc])
            if ok_cc
                @warn "nvcc cannot use the default host compiler; building with -ccbin $cc" nvcc_message =
                    _output_excerpt(output, 3, 0)
                return ["-ccbin", cc]
            end
            push!(tried, "-ccbin $cc")
        end
        error(
            "nvcc cannot compile a trivial kernel with the default host compiler" *
            (isempty(tried) ? "" : " nor with " * join(tried, ", ")) *
            ". Install a host compiler the toolkit supports (CUDA 13: GCC ≤ 15, e.g. Fedora's " *
            "gcc15-c++ package providing g++-15) or set build.nvcc_flags = [\"-ccbin\", " *
            "\"<path to a supported g++>\"]. nvcc output:\n" *
            _output_excerpt(output),
        )
    end
end

"""
    _host_compiler_candidates() -> Vector{String}

Host compilers to offer `nvcc` through `-ccbin` when the default one is
rejected: `CUDAHOSTCXX` when set, then the versioned GNU and Clang C++
compilers found on `PATH`, newest first.
"""
function _host_compiler_candidates()::Vector{String}
    found = String[]
    env_cc = get(ENV, "CUDAHOSTCXX", "")
    isempty(env_cc) || push!(found, env_cc)
    for name in ("g++-15", "g++-14", "g++-13", "g++-12", "clang++-20", "clang++-19", "clang++")
        check_command(name) && push!(found, name)
    end
    return found
end

"""
    _run_capture(cmd, logfile) -> (ok::Bool, output::String)

Run `cmd` with stdout and stderr merged into `logfile` and return whether
it succeeded together with the captured text. Spawn failures propagate.
"""
function _run_capture(cmd::Cmd, logfile::AbstractString)::Tuple{Bool,String}
    ok = open(logfile, "w") do f
        success(pipeline(cmd; stdout = f, stderr = f))
    end
    return ok, read(logfile, String)
end

"""
    _output_excerpt(output, head = 25, tail = 5) -> String

The first `head` and last `tail` lines of `output` with a marker for the
lines omitted between them; the whole text when it is short enough.
"""
function _output_excerpt(output::AbstractString, head::Int = 25, tail::Int = 5)::String
    lines = collect(eachline(IOBuffer(String(output))))
    length(lines) ≤ head + tail && return join(lines, "\n")
    kept = vcat(lines[1:head], ["… ($(length(lines) - head - tail) lines omitted)"])
    tail > 0 && append!(kept, lines[(end - tail + 1):end])
    return join(kept, "\n")
end

"""`true` when `nvcc` output reports an unsupported host compiler version."""
function _unsupported_host_compiler(output::AbstractString)::Bool
    return occursin(r"unsupported (GNU|clang|Microsoft Visual Studio) version"i, output)
end

# ---------------------------------------------------------------------------
# HDF5 flag detection
# ---------------------------------------------------------------------------

"""
    detect_hdf5_flags(platform::Symbol, use_mpi::Bool) -> Tuple{String,String}

Return `(cppflags, libs)` strings for HDF5 Makefile patching.
"""
function detect_hdf5_flags(platform::Symbol, use_mpi::Bool)::Tuple{String,String}
    hdf5_lib = use_mpi ? "hdf5_openmpi" : "hdf5"

    # Try pkg-config first (suppress stderr — common for pkg-config to warn)
    try
        cflags = strip(read(pipeline(`pkg-config --cflags $hdf5_lib`; stderr = devnull), String))
        libs = strip(read(pipeline(`pkg-config --libs $hdf5_lib`; stderr = devnull), String))
        libs *= " -lhdf5_fortran"
        return (cflags, libs)
    catch e
        e isa Union{ProcessFailedException,Base.IOError} || rethrow()
        @debug "pkg-config has no entry for $hdf5_lib; using the platform fallback" exception = e
    end

    # Platform-specific fallbacks
    if platform == :ubuntu || platform == :debian
        if use_mpi
            return (
                "-I/usr/include/hdf5/openmpi",
                "-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5_openmpi -lhdf5 -lhdf5_fortran",
            )
        else
            return (
                "-I/usr/include/hdf5/serial",
                "-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_fortran",
            )
        end
    elseif platform == :fedora
        lib = use_mpi ? "-lhdf5_openmpi -lhdf5 -lhdf5_fortran" : "-lhdf5 -lhdf5_fortran"
        return ("-I/usr/include", lib)
    else
        return ("-I/usr/include", "-lhdf5 -lhdf5_fortran")
    end
end

# ---------------------------------------------------------------------------
# Fedora h5pfc workaround
# ---------------------------------------------------------------------------

"""
    check_fedora_h5pfc()

Warn if Fedora's h5pfc wrapper has the known broken includedir path.
"""
function check_fedora_h5pfc()
    h5pfc = "/usr/lib64/openmpi/bin/h5pfc"
    isfile(h5pfc) || return

    content = read(h5pfc, String)
    m = match(r"includedir=\"([^\"]+)\"", content)
    isnothing(m) && return

    if occursin("/openmpi-x86_64", m.captures[1])
        @warn """
        Fedora h5pfc has a broken includedir path:
          $(m.captures[1])
        This prevents finding hdf5.mod during compilation.

        Fix with:
          sudo sed -i 's|/openmpi-x86_64||g' $h5pfc

        Or manually edit $h5pfc and remove '/openmpi-x86_64' from includedir.
        """
    end
end

# ---------------------------------------------------------------------------
# Hardware fingerprint (run provenance)
# ---------------------------------------------------------------------------

"""
    _hardware_fingerprint(; gpu_probe = false) -> Dict{String,Any}

Platform fingerprint for run metadata, using Julia's own introspection:
host, OS/kernel, CPU model and logical core count, total memory, Julia
version, Julia and BLAS thread counts. With `gpu_probe = true` an
`nvidia-smi` query records name, VRAM, driver version, and compute
capability of every visible GPU, one `;`-separated entry per device
(`"unavailable"` when the tool or a device is absent). Together with the
config and the git commits this makes every result attributable to
config + commit + hardware.
"""
function _hardware_fingerprint(; gpu_probe::Bool = false)::Dict{String,Any}
    cpu = Sys.cpu_info()
    d = Dict{String,Any}(
        "host" => gethostname(),
        "os" => "$(Sys.KERNEL) $(Sys.MACHINE)",
        "cpu_model" => isempty(cpu) ? "unknown" : cpu[1].model,
        "cpu_threads" => Sys.CPU_THREADS,
        "total_memory_gib" => round(Sys.total_memory() / 2^30; digits = 1),
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
    )
    if gpu_probe
        gpu = try
            strip(
                read(
                    `nvidia-smi --query-gpu=name,memory.total,driver_version,compute_cap --format=csv,noheader`,
                    String,
                ),
            )
        catch
            ""
        end
        d["gpu"] = isempty(gpu) ? "unavailable" : join(strip.(split(String(gpu), '\n')), "; ")
    end
    return d
end

# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

"""
    nproc_available() -> Int

Return the number of available CPU cores for parallel compilation.
"""
function nproc_available()::Int
    try
        parse(Int, strip(read(`nproc`, String)))
    catch
        Sys.CPU_THREADS
    end
end
