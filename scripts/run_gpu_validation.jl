#!/usr/bin/env julia
# =============================================================================
# Nbody6Dynamics — GPU validation of a CUDA host
# =============================================================================
#
# Usage:
#   julia scripts/run_gpu_validation.jl [--stages=suite,gpu,cpu,bench] [--dry-run]
#                                       [--n=20000,50000,100000] [--threads=4,8]
#                                       [--gpus="0;0,1"] [--tcrit=0.25]
#                                       [--max-retries=1] [--retry-stages=suite]
#                                       [--retry-signals=4,6,7,11] [--stop-on-failure]
#
# Runs the GPU-gated test suite, the GPU and CPU pipelines of input_files/gpu
# and the scaling benchmark as logged stages, and collects the host record,
# every log, the benchmark results and a summary under
# runs/gpu_validation_<host>_<timestamp>/ (see `run_gpu_validation`). The
# four probe stages `gpu_merger_600k`, `gpu_single_600k`, `cpu_merger_600k`,
# `cpu_single_600k` run `input_files/gpu/<stage>.toml`; with
# `--stop-on-failure` the stages form a gated chain. The exit code is 1 when
# any stage failed.

include(joinpath(@__DIR__, "activate.jl"))

using Nbody6Dynamics
using TOML

# SIGINT raises `InterruptException` (a plain script would exit at once): the
# handler in `_tee_run` then interrupts the running stage, which runs in its
# own session and would otherwise be orphaned together with its engine.
Base.exit_on_sigint(false)

const USAGE =
    "usage: run_gpu_validation.jl [--stages=suite,gpu,cpu,bench] [--dry-run] " *
    "[--n=N1,N2,...] [--threads=T1,T2,...] [--gpus=\"0;0,1\"] [--tcrit=T] " *
    "[--max-retries=K] [--retry-stages=suite,...] [--retry-signals=4,6,7,11] " *
    "[--stop-on-failure]"

function usage()
    println(stderr, USAGE)
    exit(2)
end

parse_ints(s) = [parse(Int, x) for x in split(s, ',')]

"""Value of `--name=value` in `arg`, or `nothing` when `arg` is another option."""
option(arg, name) = startswith(arg, "--$name=") ? arg[(length(name) + 4):end] : nothing

# Only the options given on the command line are forwarded, so the defaults
# are those of `run_gpu_validation` and stated nowhere else.
function main()
    kwargs = Dict{Symbol,Any}()
    for a in ARGS
        if a == "--dry-run"
            kwargs[:dry_run] = true
        elseif a == "--stop-on-failure"
            kwargs[:stop_on_failure] = true
        elseif (v = option(a, "stages")) !== nothing
            kwargs[:stages] = Symbol.(split(v, ','))
        elseif (v = option(a, "n")) !== nothing
            kwargs[:bench_n] = parse_ints(v)
        elseif (v = option(a, "threads")) !== nothing
            kwargs[:bench_threads] = parse_ints(v)
        elseif (v = option(a, "gpus")) !== nothing
            kwargs[:bench_gpu_lists] = [parse_ints(l) for l in split(v, ';')]
        elseif (v = option(a, "tcrit")) !== nothing
            kwargs[:bench_tcrit] = parse(Float64, v)
        elseif (v = option(a, "max-retries")) !== nothing
            kwargs[:max_retries] = parse(Int, v)
        elseif (v = option(a, "retry-stages")) !== nothing
            kwargs[:retry_stages] = Symbol.(split(v, ','))
        elseif (v = option(a, "retry-signals")) !== nothing
            kwargs[:retry_signals] = parse_ints(v)
        else
            usage()
        end
    end
    dir = run_gpu_validation(; kwargs...)
    summary = TOML.parsefile(joinpath(dir, "VALIDATION.toml"))
    results = summary["results"]
    for s in summary["stages"]
        r = results[s]
        println(
            rpad(":" * s, 18),
            r["status"],
            haskey(r, "seconds") ? "  ($(r["seconds"]) s)" : "",
            haskey(r, "retried_signals") ?
            "  (retried after signal $(join(r["retried_signals"], ", ")))" : "",
            haskey(r, "reason") ? "  ($(r["reason"]))" : "",
        )
    end
    println("Validation directory: $dir")
    exit(all(r["status"] in ("passed", "planned") for r in values(results)) ? 0 : 1)
end

main()
