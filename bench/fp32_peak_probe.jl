# =============================================================================
# Measured FP32 peak of the device, with the peak probe of GPUDiagnostics.jl
# (the package the paper's inverse-Thomson section uses for its FP64 peaks):
# a sweep of independent fused multiply–add chains over chain count and
# launch size, never routed to matrix units. The result is the reference of
# the "fraction of peak" panel of the cross-hardware figure.
#
#   julia bench/fp32_peak_probe.jl [--backend=cuda|cpu] [--out=<dir>]
#
# Writes <out>/fp32_peak_<hostname>.toml (default out: ./runs). The script
# keeps its own shared environment, `@nb6_fp32_peak`, and installs
# GPUDiagnostics (unregistered, by URL), KernelAbstractions and CUDA into it
# on first use; nothing touches the package environments.
# =============================================================================

using Pkg, TOML, Dates

backend_name = "cuda"
out_dir = "runs"
for a in ARGS
    if startswith(a, "--backend=")
        global backend_name = a[(length("--backend=") + 1):end]
    elseif startswith(a, "--out=")
        global out_dir = a[(length("--out=") + 1):end]
    else
        println(stderr, "usage: julia bench/fp32_peak_probe.jl [--backend=cuda|cpu] [--out=<dir>]")
        exit(2)
    end
end
backend_name in ("cuda", "cpu") || error("--backend must be cuda or cpu")

# A git configuration that rewrites https to ssh makes Pkg's libgit2 clone
# ask for a key; the git CLI honours the same rewrite with the user's agent.
ENV["JULIA_PKG_USE_CLI_GIT"] = "true"
Pkg.activate("nb6_fp32_peak"; shared = true, io = devnull)
deps = Set(keys(Pkg.project().dependencies))
"GPUDiagnostics" in deps ||
    Pkg.add(url = "https://github.com/SebastianM-C/GPUDiagnostics.jl"; io = devnull)
"KernelAbstractions" in deps || Pkg.add("KernelAbstractions"; io = devnull)
backend_name == "cuda" && !("CUDA" in deps) && Pkg.add("CUDA"; io = devnull)

using GPUDiagnostics, KernelAbstractions

backend, device = if backend_name == "cuda"
    @eval using CUDA
    Base.invokelatest(CUDA.functional) || error("CUDA.jl is not functional on this host")
    dev = Base.invokelatest(CUDA.device)
    Base.invokelatest(CUDA.CUDABackend), Base.invokelatest(CUDA.name, dev)
else
    KernelAbstractions.CPU(), Sys.cpu_info()[1].model
end

record = Dict{String,Any}(
    "date" => Dates.format(now(), "yyyy-mm-dd HH:MM:SS"),
    "host" => gethostname(),
    "backend" => backend_name,
    "device" => device,
    "gpudiagnostics_version" => string(pkgversion(GPUDiagnostics)),
)
for T in (Float32, Float64)
    p = peak_flops_probe(backend, T)
    key = T === Float32 ? "fp32" : "fp64"
    record[key] = Dict{String,Any}(
        "tflops" => p.flops / 1e12,
        "chains" => p.chains,
        "n_threads" => p.n_threads,
        "workgroup" => p.workgroup,
        "trials" => p.trials,
        "sweep_chains" => collect(p.sweep_chains),
        "sweep_tflops" => collect(p.sweep_flops) ./ 1e12,
    )
    println(
        rpad(key, 6),
        round(p.flops / 1e12; digits = 3),
        " TFLOP/s  (",
        p.chains,
        " chains × ",
        p.n_threads,
        " threads, best of ",
        p.trials,
        ")",
    )
end
mkpath(out_dir)
path = joinpath(out_dir, "fp32_peak_$(gethostname()).toml")
open(path, "w") do io
    TOML.print(io, record)
end
println("written: ", path)
