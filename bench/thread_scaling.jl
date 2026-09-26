#!/usr/bin/env julia
# =============================================================================
# Thread- and N-scaling benchmark of the Nbody6++GPU backend through the
# merger pipeline: wall time and CPU efficiency per (N, OMP threads) from the
# run telemetry, for the cost model of the science sweeps.
#
#   julia bench/thread_scaling.jl [N_list] [thread_list] [tcrit]
#   e.g. julia bench/thread_scaling.jl 1000,5000,20000 1,2,4,8,16,22 0.5
#
# Results: bench/results/thread_scaling_<timestamp>.csv (one row per run) and
# a console table with the fitted wall-time exponent t ∝ N^α at the best
# thread count. Runs live under bench/runs/ (ignored by git).
# =============================================================================

include(joinpath(@__DIR__, "activate.jl"))

using Nbody6Dynamics, TOML, Dates, Printf

# SIGINT raises `InterruptException` (a plain script would exit at once): the
# pipeline handler then terminates the engine, which runs in its own session.
Base.exit_on_sigint(false)

include("merger_case.jl")

const BENCH = @__DIR__
const PROJ = normpath(joinpath(BENCH, ".."))
parse_list(s, T) = [parse(T, x) for x in split(s, ',')]
N_list = length(ARGS) ≥ 1 ? parse_list(ARGS[1], Int) : [1000, 5000, 20000]
threads = length(ARGS) ≥ 2 ? parse_list(ARGS[2], Int) : [1, 2, 4, 8, 16, Sys.CPU_THREADS]
tcrit = length(ARGS) ≥ 3 ? parse(Float64, ARGS[3]) : 0.5

runs_dir = joinpath(BENCH, "runs")
results_dir = joinpath(BENCH, "results")
mkpath(runs_dir)
mkpath(results_dir)
stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
csv_path = joinpath(results_dir, "thread_scaling_$(stamp).csv")

function pipeline_toml(merger_path, nthreads, prefix)
    """
    [install]
    enabled = false
    install_dir = "backend/Nbody6PPGPU-beijing"

    [build]
    enable_gpu = false

    [simulation]
    run_test = true
    input_file = "../../input_files/N1k_quick.inp"
    runs_dir = "$(runs_dir)"
    binary_name = "nbody6++"
    mpi_ranks = 1
    omp_threads = $nthreads
    run_id_prefix = "$prefix"
    monitor = false
    telemetry_interval = 2.0
    startup_timeout = 1800.0

    [postprocess]
    enabled = false

    [visualization]
    enabled = false

    [merger]
    enabled = true
    config_file = "$(merger_path)"
    """
end

rows = NamedTuple[]
open(csv_path, "w") do io
    println(
        io,
        "N_requested,N_total,omp_threads,wall_s,cpu_user_s,cpu_system_s,cpu_efficiency,mean_cores_busy,peak_rss_mib,backend_total_s,backend_reg_s,backend_irr_s,backend_ks_s,backend_mdot_s,gflops_mean",
    )
end
for N in N_list, nt in threads
    prefix = "bench_N$(N)_t$(nt)"
    mpath = joinpath(runs_dir, "merger_N$(N).toml")
    isfile(mpath) || write(mpath, merger_toml(N, tcrit))
    cpath = joinpath(runs_dir, "$(prefix).toml")
    write(cpath, pipeline_toml(mpath, nt, prefix))
    cfg = load_config(cpath)
    @info "Benchmark: N = $N, OMP threads = $nt"
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
    )
    push!(rows, row)
    open(csv_path, "a") do io
        println(io, join(string.(values(row)), ","))
    end
    @info @sprintf(
        "  wall %.1f s, efficiency %.2f, cores busy %s (pipeline %.1f s)",
        row.wall_s,
        row.cpu_efficiency,
        string(row.mean_cores_busy),
        time() - t0
    )
end

# Console table and cost model
println("\nN_total  threads  wall[s]  eff   cores  backend[s]  reg   irr   ks    mdot")
for r in rows
    @printf(
        "%7d  %7d  %7.1f  %.2f  %5.1f  %9.1f  %5.1f %5.1f %5.1f %5.1f\n",
        r.N_total,
        r.omp_threads,
        r.wall_s,
        r.cpu_efficiency,
        r.mean_cores_busy,
        r.backend_total_s,
        r.backend_reg_s,
        r.backend_irr_s,
        r.backend_ks_s,
        r.backend_mdot_s
    )
end
# Best thread count per N (minimum wall time) and the exponent of t_wall(N)
best = Dict{Int,NamedTuple}()
for r in rows
    (!haskey(best, r.N_total) || r.wall_s < best[r.N_total].wall_s) && (best[r.N_total] = r)
end
Ns = sort(collect(keys(best)))
if length(Ns) ≥ 2
    x = log.(Float64.(Ns))
    y = log.([best[n].wall_s for n in Ns])
    α = (length(x) * sum(x .* y) - sum(x) * sum(y)) / (length(x) * sum(x .^ 2) - sum(x)^2)
    @printf(
        "\nBest thread count per N: %s\n",
        join(
            [
                "N=$(n): $(best[n].omp_threads) threads, $(round(best[n].wall_s; digits = 1)) s" for
                n in Ns
            ],
            "; ",
        )
    )
    @printf(
        "Wall time per %.2f NB time units scales as N^%.2f at the best thread count\n",
        tcrit,
        α
    )
    n_ref = Ns[end]
    @printf(
        "Estimate for 100 NB time units at N = %d: %.1f h\n",
        n_ref,
        best[n_ref].wall_s / tcrit * 100 / 3600
    )
end
println("Results: ", csv_path)
