#!/usr/bin/env julia
# =============================================================================
# Benchmark cells of the Nbody6++GPU backend on the device: one run of each
# (cell, N, host threads, GPU list, MPI ranks) point through the merger
# pipeline, engine only, with the run telemetry read back into one CSV row
# per run — wall time, the backend's own timing table (regular force,
# irregular force, KS, stellar evolution, the MPI barriers), the
# regular-force kernel throughput, GPU utilisation, power and memory, host
# RSS and CPU efficiency. The cells are those of bench/merger_case.jl: the
# single cell, one King cluster of N stars at rest, and the dual cell, two
# King clusters of N/2 on an eccentric orbit.
#
#   julia bench/gpu_cells.jl [--cells=single,dual] [--n=1000000] [--threads=8]
#                            [--gpus="0;0,1"] [--tcrit=2.0] [--cpu] [--mpi-ranks=1]
#                            [--startup-timeout=14400] [--label=<tag>]
#
# --gpus     ';'-separated GPU_LIST values ("0" = the first device, "0,1" = the
#            first two); "" runs no device rows.
# --cpu      adds a row on the AVX binary for every (cell, N, threads).
# --mpi-ranks=R > 1 runs the device rows on the MPI+CUDA binary with R ranks,
#            every rank applying the same GPU list (the engine has no
#            per-rank device slicing).
# --startup-timeout=S the engine's start-up watchdog in seconds (default
#            14400): the time allowed up to the first adjustment beyond t = 0,
#            which the CPU binary at 10^6 bodies exceeds within an hour.
# --label    a tag written into the run prefix and the CSV name.
#
# Binaries: NBODY6_GPU_BACKEND (default backend/Nbody6PPGPU-beijing-gpu),
# NBODY6_CPU_BACKEND (default backend/Nbody6PPGPU-beijing) and, for MPI rows,
# NBODY6_MPI_GPU_BACKEND (default backend/Nbody6PPGPU-beijing-mpi-gpu). A run
# that fails is recorded with status "failed" and the grid continues.
# Results: bench/results/gpu_cells_<label>_<timestamp>.csv; runs live under
# bench/runs/ (ignored by git).
# =============================================================================

include(joinpath(@__DIR__, "activate.jl"))

using Nbody6Dynamics, TOML, Dates, Printf
using Nbody6Dynamics: detect_compute_capabilities

include("merger_case.jl")

const BENCH = @__DIR__
const PROJ = normpath(joinpath(BENCH, ".."))

"""Value of `--name=value` in `arg`, or `nothing` when `arg` is another option."""
option(arg, name) = startswith(arg, "--$name=") ? arg[(length(name) + 4):end] : nothing
parse_list(s, T) = isempty(strip(s)) ? T[] : [parse(T, x) for x in split(s, ',')]

cells = ["single", "dual"]
N_list = [1_000_000]
threads = [8]
gpu_lists = [[0]]
tcrit = 2.0
with_cpu = false
mpi_ranks = 1
startup_timeout = 14400.0
label = ""
for a in ARGS
    if a == "--cpu"
        global with_cpu = true
    elseif (v = option(a, "cells")) !== nothing
        global cells = String.(split(v, ','))
    elseif (v = option(a, "n")) !== nothing
        global N_list = parse_list(v, Int)
    elseif (v = option(a, "threads")) !== nothing
        global threads = parse_list(v, Int)
    elseif (v = option(a, "gpus")) !== nothing
        global gpu_lists = [parse_list(l, Int) for l in split(v, ';') if !isempty(strip(l))]
    elseif (v = option(a, "tcrit")) !== nothing
        global tcrit = parse(Float64, v)
    elseif (v = option(a, "mpi-ranks")) !== nothing
        global mpi_ranks = parse(Int, v)
    elseif (v = option(a, "startup-timeout")) !== nothing
        global startup_timeout = parse(Float64, v)
    elseif (v = option(a, "label")) !== nothing
        global label = v
    else
        println(
            stderr,
            "usage: gpu_cells.jl [--cells=single,dual] [--n=N1,N2] [--threads=T1,T2] " *
            "[--gpus=\"0;0,1\"] [--tcrit=T] [--cpu] [--mpi-ranks=R] " *
            "[--startup-timeout=S] [--label=tag]",
        )
        exit(2)
    end
end
for c in cells
    c in ("single", "dual") || error("unknown cell \"$c\"; choose single or dual")
end
mpi_ranks ≥ 1 || error("--mpi-ranks must be ≥ 1")
tcrit > 0 || error("--tcrit must be > 0")
startup_timeout > 0 || error("--startup-timeout must be > 0 s")
isempty(gpu_lists) && !with_cpu && error("nothing to run: no GPU lists and no --cpu")

const BACKEND_CPU =
    abspath(get(ENV, "NBODY6_CPU_BACKEND", joinpath(PROJ, "backend", "Nbody6PPGPU-beijing")))
const BACKEND_GPU =
    abspath(get(ENV, "NBODY6_GPU_BACKEND", joinpath(PROJ, "backend", "Nbody6PPGPU-beijing-gpu")))
const BACKEND_MPI_GPU = abspath(
    get(ENV, "NBODY6_MPI_GPU_BACKEND", joinpath(PROJ, "backend", "Nbody6PPGPU-beijing-mpi-gpu")),
)

# The binaries the grid needs must exist before anything runs; the selection
# by suffix tag is the one the pipeline applies (see `_find_binary`).
with_cpu &&
    @info "CPU binary: $(Nbody6Dynamics._find_binary(BACKEND_CPU, "nbody6++", BuildConfig()))"
if !isempty(gpu_lists)
    tree = mpi_ranks > 1 ? BACKEND_MPI_GPU : BACKEND_GPU
    build = BuildConfig(; enable_gpu = true, enable_mpi = mpi_ranks > 1)
    @info "GPU binary: $(Nbody6Dynamics._find_binary(tree, "nbody6++", build))"
    caps = detect_compute_capabilities()
    isempty(caps) && error("nvidia-smi reports no device; the GPU rows cannot run")
    @info "Visible compute capabilities: $(join(caps, ", "))"
end

runs_dir = joinpath(BENCH, "runs")
results_dir = joinpath(BENCH, "results")
mkpath(runs_dir)
mkpath(results_dir)
stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
tag = isempty(label) ? "" : "$(label)_"
csv_path = joinpath(results_dir, "gpu_cells_$(tag)$(stamp).csv")

function pipeline_toml(merger_path, nthreads, prefix; gpu::Bool, gpu_list::Vector{Int}, ranks::Int)
    tree = gpu ? (ranks > 1 ? BACKEND_MPI_GPU : BACKEND_GPU) : BACKEND_CPU
    """
    [install]
    enabled = false
    install_dir = "$(tree)"

    [build]
    enable_gpu = $(gpu)
    enable_mpi = $(gpu && ranks > 1)

    [simulation]
    run_test = true
    input_file = "../../input_files/N1k_quick.inp"
    runs_dir = "$(runs_dir)"
    binary_name = "nbody6++"
    mpi_ranks = $(gpu ? ranks : 1)
    omp_threads = $nthreads
    gpu_list = $(gpu ? string(gpu_list) : "[]")
    run_id_prefix = "$prefix"
    monitor = false
    telemetry_interval = 2.0
    startup_timeout = $(startup_timeout)
    exit_grace = 600.0

    [postprocess]
    enabled = false

    [visualization]
    enabled = false

    [merger]
    enabled = true
    config_file = "$(merger_path)"
    """
end

const COLUMNS = (
    "cell,N_requested,N_total,variant,gpu_list,mpi_ranks,omp_threads,tcrit,status,wall_s," *
    "backend_total_s,backend_reg_s,backend_irr_s,backend_ks_s,backend_mdot_s,backend_barr_s," *
    "gflops_mean,gflops_peak,kernel,gpu_util_mean,gpu_util_peak,gpu_power_mean_w," *
    "gpu_mem_peak_mib,peak_rss_mib,cpu_efficiency,gpu_devices,machine,gpu,run_id"
)

rows = NamedTuple[]
open(csv_path, "w") do io
    println(io, COLUMNS)
end

"""One CSV cell: a string quoted when it holds a comma, a number as is."""
csv_field(x) = x isa AbstractString && occursin(',', x) ? "\"$x\"" : string(x)

function run_case(cell, N, nt, variant, gpu_list, ranks)
    tag_v = variant == "cpu" ? "cpu" : "gpu" * join(gpu_list, "") * (ranks > 1 ? "_r$ranks" : "")
    prefix = "cells_$(tag)$(cell)_N$(N)_t$(nt)_$(tag_v)"
    mpath = joinpath(runs_dir, "$(cell)_N$(N)_tcrit$(tcrit).toml")
    isfile(mpath) || write(mpath, cell_toml(cell, N, tcrit))
    cpath = joinpath(runs_dir, "$(prefix).toml")
    write(
        cpath,
        pipeline_toml(
            mpath,
            nt,
            prefix;
            gpu = variant == "gpu",
            gpu_list = gpu_list,
            ranks = ranks,
        ),
    )
    cfg = load_config(cpath)
    @info "Cell $cell: N = $N, threads = $nt, $variant" gpu_list ranks
    t0 = time()
    status = "completed"
    info = Dict{String,Any}()
    N_total = 0
    try
        run_pipeline(cfg; base_dir = PROJ)
        run_dir = Nbody6Dynamics._find_latest_run(cfg, PROJ)
        info = TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))
        ic = TOML.parsefile(joinpath(run_dir, "output", "merger_ic.toml"))
        N_total = Int(ic["meta"]["N_total"])
        Nbody6Dynamics._engine_completed(run_dir) === false && (status = "halted")
    catch e
        status = "failed"
        @error "Cell $cell at N = $N ($variant) failed; the grid continues" exception =
            (e, catch_backtrace())
    end
    run = get(info, "run", Dict{String,Any}())
    tel = get(info, "telemetry", Dict{String,Any}())
    hw = get(info, "hardware", Dict{String,Any}())
    bt = get(tel, "backend_timing", Dict{String,Any}())
    gf = get(tel, "force_kernel_gflops", Dict{String,Any}())
    row = (
        cell = cell,
        N_requested = N,
        N_total = N_total,
        variant = variant,
        gpu_list = variant == "cpu" ? "" : join(gpu_list, " "),
        mpi_ranks = variant == "cpu" ? 1 : ranks,
        omp_threads = nt,
        tcrit = tcrit,
        status = status,
        wall_s = get(run, "elapsed_seconds", NaN),
        backend_total_s = get(bt, "total", NaN),
        backend_reg_s = get(bt, "reg", NaN),
        backend_irr_s = get(bt, "irr", NaN),
        backend_ks_s = get(bt, "ks", NaN),
        backend_mdot_s = get(bt, "mdot", NaN),
        backend_barr_s = get(bt, "barr_p", NaN) + get(bt, "barr_r", NaN) + get(bt, "barr_i", NaN),
        gflops_mean = get(gf, "mean", NaN),
        gflops_peak = get(gf, "peak", NaN),
        kernel = get(gf, "kernel", ""),
        gpu_util_mean = get(tel, "mean_gpu_util_pct", NaN),
        gpu_util_peak = get(tel, "peak_gpu_util_pct", NaN),
        gpu_power_mean_w = get(tel, "mean_gpu_power_w", NaN),
        gpu_mem_peak_mib = get(tel, "peak_gpu_mem_used_mib", NaN),
        peak_rss_mib = get(tel, "peak_rss_mib", NaN),
        cpu_efficiency = get(tel, "cpu_efficiency", NaN),
        gpu_devices = length(get(run, "gpu_devices", String[])),
        machine = get(hw, "machine", get(hw, "host", "")),
        gpu = get(hw, "gpu", ""),
        run_id = get(run, "id", ""),
    )
    push!(rows, row)
    open(csv_path, "a") do io
        println(io, join(csv_field.(values(row)), ","))
    end
    @info @sprintf(
        "  %s: wall %.1f s, regular force %.1f s, irregular %.1f s, %.0f GFLOP/s, GPU util %.0f %% (pipeline %.1f s)",
        status,
        row.wall_s,
        row.backend_reg_s,
        row.backend_irr_s,
        row.gflops_mean,
        row.gpu_util_mean,
        time() - t0
    )
    return row
end

for cell in cells, N in N_list, nt in threads
    with_cpu && run_case(cell, N, nt, "cpu", Int[], 1)
    for gl in gpu_lists
        run_case(cell, N, nt, "gpu", gl, mpi_ranks)
    end
end

println(
    "\ncell    N_total  thr  variant  GPUs  ranks  status     wall[s]   reg[s]   irr[s]   GFLOP/s  util%",
)
for r in rows
    @printf(
        "%-7s %8d  %3d  %-7s  %-5s  %5d  %-9s %9.1f %8.1f %8.1f %9.0f  %5.1f\n",
        r.cell,
        r.N_total,
        r.omp_threads,
        r.variant,
        isempty(r.gpu_list) ? "-" : r.gpu_list,
        r.mpi_ranks,
        r.status,
        r.wall_s,
        r.backend_reg_s,
        r.backend_irr_s,
        r.gflops_mean,
        r.gpu_util_mean
    )
end
println("Results: ", csv_path)
any(r -> r.status != "completed", rows) && exit(1)
