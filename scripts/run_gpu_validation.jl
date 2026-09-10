#!/usr/bin/env julia
# =============================================================================
# Nbody6Dynamics — GPU validation of a CUDA host
# =============================================================================
#
# Usage:
#   julia scripts/run_gpu_validation.jl [--stages=suite,gpu,cpu,bench] [--dry-run]
#                                       [--n=20000,50000,100000] [--threads=4,8]
#                                       [--gpus="0;0,1"] [--tcrit=0.25]
#
# Runs the GPU-gated test suite, the GPU and CPU pipelines of input_files/gpu
# and the scaling benchmark as logged stages, and collects the host record,
# every log, the benchmark results and a summary under
# runs/gpu_validation_<host>_<timestamp>/ (see `run_gpu_validation`). The
# exit code is 1 when any stage failed.

include(joinpath(@__DIR__, "..", "activate.jl"))

using Nbody6Dynamics
using TOML

function usage()
    println(
        stderr,
        "usage: run_gpu_validation.jl [--stages=suite,gpu,cpu,bench] [--dry-run] " *
        "[--n=20000,50000,100000] [--threads=4,8] [--gpus=\"0;0,1\"] [--tcrit=0.25]",
    )
    exit(2)
end

parse_ints(s) = [parse(Int, x) for x in split(s, ',')]

function main()
    stages = [:suite, :gpu, :cpu, :bench]
    n = [20000, 50000, 100000]
    threads = [4, 8]
    gpus = nothing
    tcrit = 0.25
    dry_run = false
    for a in ARGS
        if a == "--dry-run"
            dry_run = true
        elseif startswith(a, "--stages=")
            stages = Symbol.(split(a[(length("--stages=") + 1):end], ','))
        elseif startswith(a, "--n=")
            n = parse_ints(a[(length("--n=") + 1):end])
        elseif startswith(a, "--threads=")
            threads = parse_ints(a[(length("--threads=") + 1):end])
        elseif startswith(a, "--gpus=")
            gpus = [parse_ints(s) for s in split(a[(length("--gpus=") + 1):end], ';')]
        elseif startswith(a, "--tcrit=")
            tcrit = parse(Float64, a[(length("--tcrit=") + 1):end])
        else
            usage()
        end
    end
    dir = run_gpu_validation(;
        stages = stages,
        bench_n = n,
        bench_threads = threads,
        bench_gpu_lists = gpus,
        bench_tcrit = tcrit,
        dry_run = dry_run,
    )
    results = TOML.parsefile(joinpath(dir, "VALIDATION.toml"))["results"]
    for s in stages
        r = results[String(s)]
        println(
            rpad(":" * String(s), 8),
            r["status"],
            haskey(r, "seconds") ? "  ($(r["seconds"]) s)" : "",
        )
    end
    println("Validation directory: $dir")
    exit(any(r["status"] == "failed" for r in values(results)) ? 1 : 0)
end

main()
