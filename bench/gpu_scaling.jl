#!/usr/bin/env julia
# =============================================================================
# GPU-versus-CPU scaling benchmark of the Nbody6++GPU backend through the
# merger pipeline: wall time, CPU efficiency, the backend's own timing table,
# regular-force kernel throughput and GPU utilisation per (N, binary variant,
# GPU list, OMP threads), from the run telemetry. Both binaries must exist
# beforehand: the CPU one (`nbody6++.avx`) under NBODY6_CPU_BACKEND (default
# backend/Nbody6PPGPU-beijing) and the GPU one (`nbody6++.avx.gpu`) under
# NBODY6_GPU_BACKEND (default: the same tree; input_files/gpu/ builds it into
# backend/Nbody6PPGPU-beijing-gpu).
#
#   NBODY6_GPU_BACKEND=backend/Nbody6PPGPU-beijing-gpu \
#     julia bench/gpu_scaling.jl [N_list] [thread_list] [gpu_lists] [tcrit]
#   e.g. julia bench/gpu_scaling.jl 20000,50000,100000 4,8 "0;0,1" 0.25
#
# gpu_lists: ';'-separated GPU_LIST values ("0" = first device, "0,1" = the
# first two); a CPU-binary row is added for every (N, threads). Results:
# bench/results/gpu_scaling_<timestamp>.csv (one row per run) and a console
# table with the GPU speed-up over the CPU binary at equal N and threads.
# Runs live under bench/runs/ (ignored by git).
# =============================================================================

using Pkg
Pkg.activate(@__DIR__; io = devnull)
Pkg.develop(; path = joinpath(@__DIR__, ".."), io = devnull)
Pkg.instantiate(; io = devnull)

using Nbody6Dynamics, TOML, Dates, Printf

include("merger_case.jl")

const BENCH = @__DIR__
const PROJ = normpath(joinpath(BENCH, ".."))
const BACKEND_CPU =
    abspath(get(ENV, "NBODY6_CPU_BACKEND", joinpath(PROJ, "backend", "Nbody6PPGPU-beijing")))
const BACKEND_GPU = abspath(get(ENV, "NBODY6_GPU_BACKEND", BACKEND_CPU))
parse_list(s, T) = [parse(T, x) for x in split(s, ',')]
N_list = length(ARGS) ≥ 1 ? parse_list(ARGS[1], Int) : [20000, 50000]
threads = length(ARGS) ≥ 2 ? parse_list(ARGS[2], Int) : [4, 8]
gpu_lists = length(ARGS) ≥ 3 ? [parse_list(s, Int) for s in split(ARGS[3], ';')] : [[0]]
tcrit = length(ARGS) ≥ 4 ? parse(Float64, ARGS[4]) : 0.25

# Both binaries must be present; the selection by suffix tag is the one the
# pipeline applies (see `_find_binary`).
cpu_binary = Nbody6Dynamics._find_binary(BACKEND_CPU, "nbody6++", BuildConfig())
gpu_binary = Nbody6Dynamics._find_binary(BACKEND_GPU, "nbody6++", BuildConfig(; enable_gpu = true))
@info "CPU binary: $cpu_binary"
@info "GPU binary: $gpu_binary"
caps = detect_compute_capabilities()
isempty(caps) && error("nvidia-smi reports no device; the GPU rows cannot run")
@info "Visible compute capabilities: $(join(caps, ", "))"

runs_dir = joinpath(BENCH, "runs")
results_dir = joinpath(BENCH, "results")
mkpath(runs_dir)
mkpath(results_dir)
stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
csv_path = joinpath(results_dir, "gpu_scaling_$(stamp).csv")

function pipeline_toml(merger_path, nthreads, prefix; gpu::Bool, gpu_list::Vector{Int})
    """
    [install]
    enabled = false
    install_dir = "$(gpu ? BACKEND_GPU : BACKEND_CPU)"

    [build]
    enable_gpu = $(gpu)

    [simulation]
    run_test = true
    input_file = "../../input_files/N1k_quick.inp"
    runs_dir = "$(runs_dir)"
    binary_name = "nbody6++"
    mpi_ranks = 1
    omp_threads = $nthreads
    gpu_list = $(gpu ? string(gpu_list) : "[]")
    run_id_prefix = "$prefix"
    monitor = false
    telemetry_interval = 2.0
    startup_timeout = 3600.0

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
    "N_requested,N_total,variant,gpu_list,omp_threads,wall_s,cpu_user_s,cpu_system_s," *
    "cpu_efficiency,mean_cores_busy,peak_rss_mib,backend_total_s,backend_reg_s,backend_irr_s," *
    "backend_ks_s,backend_mdot_s,gflops_mean,kernel,gpu_util_mean,gpu_devices"
)

rows = NamedTuple[]
open(csv_path, "w") do io
    println(io, COLUMNS)
end

function run_case(N, nt, variant, gpu_list)
    tag = variant == "cpu" ? "cpu" : "gpu" * join(gpu_list, "")
    prefix = "gpubench_N$(N)_t$(nt)_$(tag)"
    mpath = joinpath(runs_dir, "merger_N$(N).toml")
    isfile(mpath) || write(mpath, merger_toml(N, tcrit))
    cpath = joinpath(runs_dir, "$(prefix).toml")
    write(cpath, pipeline_toml(mpath, nt, prefix; gpu = variant == "gpu", gpu_list = gpu_list))
    cfg = load_config(cpath)
    @info "Benchmark: N = $N, threads = $nt, $variant" gpu_list
    t0 = time()
    run_pipeline(cfg; base_dir = PROJ)
    run_dir = Nbody6Dynamics._find_latest_run(cfg, PROJ)
    info = TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))
    tel = info["telemetry"]
    bt = get(tel, "backend_timing", Dict{String,Any}())
    gf = get(tel, "force_kernel_gflops", Dict{String,Any}())
    ic = TOML.parsefile(joinpath(run_dir, "output", "merger_ic.toml"))
    row = (
        N_requested = N,
        N_total = ic["meta"]["N_total"],
        variant = variant,
        gpu_list = variant == "cpu" ? "" : join(gpu_list, " "),
        omp_threads = nt,
        wall_s = info["run"]["elapsed_seconds"],
        cpu_user_s = get(tel, "cpu_user_s", NaN),
        cpu_system_s = get(tel, "cpu_system_s", NaN),
        cpu_efficiency = get(tel, "cpu_efficiency", NaN),
        mean_cores_busy = get(tel, "mean_cores_busy", NaN),
        peak_rss_mib = get(tel, "peak_rss_mib", NaN),
        backend_total_s = get(bt, "total", NaN),
        backend_reg_s = get(bt, "reg", NaN),
        backend_irr_s = get(bt, "irr", NaN),
        backend_ks_s = get(bt, "ks", NaN),
        backend_mdot_s = get(bt, "mdot", NaN),
        gflops_mean = get(gf, "mean", NaN),
        kernel = get(gf, "kernel", ""),
        gpu_util_mean = get(tel, "mean_gpu_util_pct", NaN),
        gpu_devices = length(get(info["run"], "gpu_devices", String[])),
    )
    push!(rows, row)
    open(csv_path, "a") do io
        println(io, join(string.(values(row)), ","))
    end
    @info @sprintf(
        "  wall %.1f s, efficiency %.2f, regular force %.1f s of %.1f s (pipeline %.1f s)",
        row.wall_s,
        row.cpu_efficiency,
        row.backend_reg_s,
        row.backend_total_s,
        time() - t0
    )
    return row
end

for N in N_list, nt in threads
    run_case(N, nt, "cpu", Int[])
    for gl in gpu_lists
        run_case(N, nt, "gpu", gl)
    end
end

# Console table: GPU speed-up over the CPU binary at equal N and threads
println("\nN_total  threads  variant  GPUs   wall[s]  eff   reg[s]  irr[s]  Gflops  speed-up")
for N in unique(r.N_total for r in rows), nt in threads
    cpu = filter(r -> r.N_total == N && r.omp_threads == nt && r.variant == "cpu", rows)
    isempty(cpu) && continue
    t_cpu = cpu[1].wall_s
    for r in filter(r -> r.N_total == N && r.omp_threads == nt, rows)
        @printf(
            "%7d  %7d  %7s  %-5s  %7.1f  %.2f  %6.1f  %6.1f  %6.0f  %8.2f\n",
            r.N_total,
            r.omp_threads,
            r.variant,
            isempty(r.gpu_list) ? "-" : r.gpu_list,
            r.wall_s,
            r.cpu_efficiency,
            r.backend_reg_s,
            r.backend_irr_s,
            r.gflops_mean,
            t_cpu / r.wall_s
        )
    end
end
println("Results: ", csv_path)
