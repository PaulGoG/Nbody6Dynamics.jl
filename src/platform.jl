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
        success(`which $cmd`)
    catch
        false
    end
end

"""
    check_dependencies(cfg::Nbody6Config) -> Vector{String}

Return a list of missing dependencies required for building Nbody6++.
"""
function check_dependencies(cfg::Nbody6Config)::Vector{String}
    missing_deps = String[]

    required = ["git", "gcc", "g++", "gfortran", "make"]
    if cfg.build.enable_mpi
        append!(required, ["mpicc", "mpif90", "mpirun"])
    end
    if cfg.build.enable_gpu
        push!(required, "nvcc")
    end

    for cmd in required
        check_command(cmd) || push!(missing_deps, cmd)
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
        nvcc = strip(read(`which nvcc`, String))
        return dirname(dirname(nvcc))  # nvcc lives in CUDA_HOME/bin/
    catch
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
    catch
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
# Hardware fingerprint (run provenance, §6)
# ---------------------------------------------------------------------------

"""
    _hardware_fingerprint(; gpu_probe = false) -> Dict{String,Any}

Platform fingerprint for run metadata, using Julia's own introspection:
host, OS/kernel, CPU model and logical core count, total memory, Julia
version, Julia and BLAS thread counts. With `gpu_probe = true` an
`nvidia-smi` query records the GPU name, VRAM, and driver version
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
                    `nvidia-smi --query-gpu=name,memory.total,driver_version --format=csv,noheader`,
                    String,
                ),
            )
        catch
            ""
        end
        d["gpu"] = isempty(gpu) ? "unavailable" : String(gpu)
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
