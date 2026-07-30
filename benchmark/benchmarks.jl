#!/usr/bin/env julia
# =============================================================================
# Nbody6Setup Performance Benchmarks
# =============================================================================
#
# Usage:
#   julia --project=path/to/Nbody6Setup benchmark/benchmarks.jl

using BenchmarkTools
using Nbody6Setup
using Printf

const BENCHDIR = mktempdir()

# ---------------------------------------------------------------------------
# Helpers: generate synthetic data files
# ---------------------------------------------------------------------------

function _write_fortran_record(io::IO, data::AbstractVector)
    bytes = reinterpret(UInt8, collect(data))
    marker = Int32(length(bytes))
    write(io, marker)
    write(io, bytes)
    write(io, marker)
end

function generate_conf3(path::String, n::Int)
    open(path, "w") do io
        _write_fortran_record(io, Int32[n, 1, 1, 20])
        params = zeros(Float32, 20)
        params[1] = 1.0f0
        _write_fortran_record(io, params)
        for i in 1:n
            buf = IOBuffer()
            write(buf, Float32(1.0 / n))
            write(buf, Float32.(randn(3))...)
            write(buf, Float32.(0.1 .* randn(3))...)
            write(buf, Int32(i))
            _write_fortran_record(io, take!(buf))
        end
    end
end

"""Generate diagnostics file with positional ADJUST format."""
function generate_diagnostics_positional(path::String, n_lines::Int)
    open(path, "w") do io
        println(io, " PHYSICAL SCALING:  R* = 1.0  M* = 1000.0  V* = 5.0  T* = 10.0")
        for i in 1:n_lines
            t = Float64(i) / n_lines
            @printf(io,
                " ADJUST:  %10.4f  %10.2f  %7.3f  %10.2E  %8.4f  %6d  %6d  %7.3f\n",
                t, t * 100, 1.0, 1e-6 * randn(), -0.25, 10000 - i, 500 - div(i, 2), 1.0 + 0.1 * t,
            )
        end
    end
end

"""Generate diagnostics file with key-value ADJUST + TIME[NB] + RSCALE lines."""
function generate_diagnostics_keyvalue(path::String, n_epochs::Int)
    open(path, "w") do io
        println(io, " PHYSICAL SCALING:  R* = 1.0  M* = 1000.0  V* = 5.0  T* = 10.0")
        println(io, "                    <M> = 0.500  SU = 1.0  AU = 1.0")
        for i in 1:n_epochs
            t = Float64(i) / n_epochs
            t_myr = t * 100
            de = 1e-6 * randn()
            n = 10000 - i
            npairs = 500 - div(i, 2)
            rscale = 1.0 + 0.1 * t
            @printf(io, " ADJUST:  TIME  %10.4f  T[Myr]  %10.2f  Q  %7.3f  DE  %10.2E  ETOT  %8.4f\n",
                    t, t_myr, 1.0, de, -0.25)
            @printf(io, " RMIN =    0.001 RSCALE =   %7.3f\n", rscale)
            @printf(io, " TIME[NB]   %10.4f N   %6d <NB>      0 NPAIRS   %5d\n", t, n, npairs)
        end
    end
end

# ---------------------------------------------------------------------------
# Benchmarks
# ---------------------------------------------------------------------------

println("=" ^ 60)
println("  Nbody6Setup.jl — Performance Benchmarks")
println("=" ^ 60)

# --- conf.3 reader ---
for n in [1_000, 10_000, 100_000]
    fpath = joinpath(BENCHDIR, "conf3_N$(n)")
    generate_conf3(fpath, n)

    println("\n--- read_conf3  N = $n ---")
    display(@benchmark read_conf3($fpath) samples=5 evals=1)
    println()
end

# --- Diagnostics parser (positional format) ---
for n in [100, 1_000, 10_000]
    fpath = joinpath(BENCHDIR, "diag_pos_$(n)")
    generate_diagnostics_positional(fpath, n)

    println("\n--- read_diagnostics (positional)  lines = $n ---")
    display(@benchmark read_diagnostics($fpath) samples=5 evals=1)
    println()
end

# --- Diagnostics parser (key-value + TIME[NB] merging) ---
for n in [100, 1_000, 10_000]
    fpath = joinpath(BENCHDIR, "diag_kv_$(n)")
    generate_diagnostics_keyvalue(fpath, n)

    println("\n--- read_diagnostics (key-value + TIME[NB])  epochs = $n ---")
    display(@benchmark read_diagnostics($fpath) samples=5 evals=1)
    println()
end

# --- Config loading ---
cfg_path = joinpath(BENCHDIR, "bench_config.toml")
cp(joinpath(@__DIR__, "..", "config.toml"), cfg_path)
println("\n--- load_config ---")
display(@benchmark load_config($cfg_path) samples=10)
println()

# --- _format_elapsed ---
println("\n--- _format_elapsed ---")
display(@benchmark Nbody6Setup._format_elapsed(t) setup=(t=rand()*10000.0) samples=100)
println()

# --- _auto_fps ---
println("\n--- _auto_fps ---")
display(@benchmark Nbody6Setup._auto_fps(n; target_duration=12.0, min_fps=1, max_fps=30) setup=(n=rand(1:500)) samples=100)
println()

println("\nBenchmarks complete.")
