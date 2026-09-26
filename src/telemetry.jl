# =============================================================================
# Runtime hardware telemetry: exact child CPU accounting plus process-tree
# and GPU sampling for benchmarking (provenance and profiling)
# =============================================================================

"""
    TelemetrySample

One sample of the simulation process tree and, when probed, the GPU(s).
`NaN` marks quantities unavailable at sampling time.

# Fields
- `elapsed_s`: wall-clock time since launch [s]
- `n_processes`: processes in the tree rooted at the launched PID
- `rss_mib`: resident set size summed over the tree [MiB]
- `hwm_mib`: resident high-water mark summed over the tree [MiB]
- `cpu_time_s`: cumulative user + system CPU time of the live tree [s]
- `cores_busy`: CPU-time rate since the previous sample [cores]
- `load_1min`: system 1-minute load average
- `gpu_util_pct`, `gpu_mem_util_pct`: mean over devices [%]
- `gpu_mem_used_mib`: sum over devices [MiB]
- `gpu_power_w`: sum over devices [W]
- `gpu_temp_c`: maximum over devices [°C]
"""
struct TelemetrySample
    elapsed_s::Float64
    n_processes::Int
    rss_mib::Float64
    hwm_mib::Float64
    cpu_time_s::Float64
    cores_busy::Float64
    load_1min::Float64
    gpu_util_pct::Float64
    gpu_mem_util_pct::Float64
    gpu_mem_used_mib::Float64
    gpu_power_w::Float64
    gpu_temp_c::Float64
end

const _NVIDIA_SMI_QUERY = "utilization.gpu,utilization.memory,memory.used,power.draw,temperature.gpu"

# ---------------------------------------------------------------------------
# Exact accounting (rusage)
# ---------------------------------------------------------------------------

"""
    _children_cpu_times() -> (user_s, system_s)

Cumulative user and system CPU time of every terminated, waited-for child of
this Julia process (`getrusage(RUSAGE_CHILDREN)`) in seconds. The difference
of two calls bracketing a launch is the exact CPU consumption of that child
tree, independent of any sampling. `(NaN, NaN)` off Linux.
"""
function _children_cpu_times()::Tuple{Float64,Float64}
    Sys.islinux() || return (NaN, NaN)
    # struct rusage on LP64 Linux: ru_utime and ru_stime (two timevals of two
    # longs each) followed by 14 longs.
    buf = zeros(Int64, 18)
    rc = ccall(:getrusage, Cint, (Cint, Ptr{Int64}), -1, buf)   # RUSAGE_CHILDREN = -1
    rc == 0 || return (NaN, NaN)
    return (buf[1] + buf[2] / 1e6, buf[3] + buf[4] / 1e6)
end

"""Kernel clock ticks per second (`sysconf(_SC_CLK_TCK)`), the unit of `/proc` CPU times."""
_clock_ticks_per_second()::Int = Sys.islinux() ? Int(ccall(:sysconf, Clong, (Cint,), 2)) : 100

# ---------------------------------------------------------------------------
# /proc readers
# ---------------------------------------------------------------------------

"""
    _parse_proc_stat(record) -> (ppid, utime_ticks, stime_ticks)

Parse one `/proc/<pid>/stat` record. The command name sits in parentheses and
may itself contain spaces or parentheses, so fields are located relative to
the last closing parenthesis: the record's fields 4 (ppid), 14 (utime), and
15 (stime) are positions 2, 12, and 13 of the remainder.
"""
function _parse_proc_stat(record::AbstractString)::Tuple{Int,Int,Int}
    close_idx = findlast(')', record)
    close_idx === nothing && throw(ArgumentError("malformed /proc stat record: $record"))
    fields = split(SubString(record, nextind(record, close_idx)))
    length(fields) ≥ 13 || throw(ArgumentError("malformed /proc stat record: $record"))
    return (parse(Int, fields[2]), parse(Int, fields[12]), parse(Int, fields[13]))
end

"""
    _read_proc_file(path) -> Union{String,Nothing}

Contents of a `/proc` file, or `nothing` when the process disappeared between
the directory scan and the read (the only error class tolerated here).
"""
function _read_proc_file(path::AbstractString)::Union{String,Nothing}
    try
        return read(path, String)
    catch e
        e isa Union{SystemError,Base.IOError} || rethrow()
        return nothing
    end
end

"""
    _read_proc_status_memory(pid) -> (rss_kib, hwm_kib)

`VmRSS` and `VmHWM` from `/proc/<pid>/status` [KiB]; zeros once the process
has exited. `VmHWM` is the kernel's own resident high-water mark, reset at
`exec`, so it is exact for the launched binary regardless of the sampling
cadence.
"""
function _read_proc_status_memory(pid::Integer)::Tuple{Int,Int}
    content = _read_proc_file("/proc/$pid/status")
    content === nothing && return (0, 0)
    rss = 0
    hwm = 0
    for line in eachline(IOBuffer(content))
        if startswith(line, "VmRSS:")
            rss = parse(Int, split(line)[2])
        elseif startswith(line, "VmHWM:")
            hwm = parse(Int, split(line)[2])
        end
    end
    return (rss, hwm)
end

"""
    _process_tree(root_pid) -> Vector{Int}

PIDs of `root_pid` and all its descendants from one scan of `/proc`. Empty
when the root has exited or `/proc` is unavailable.
"""
function _process_tree(root_pid::Integer)::Vector{Int}
    isdir("/proc") || return Int[]
    isfile("/proc/$root_pid/stat") || return Int[]
    children = Dict{Int,Vector{Int}}()
    for entry in readdir("/proc")
        all(isdigit, entry) || continue
        record = _read_proc_file("/proc/$entry/stat")
        record === nothing && continue
        ppid = _parse_proc_stat(record)[1]
        push!(get!(children, ppid, Int[]), parse(Int, entry))
    end
    tree = Int[]
    queue = [Int(root_pid)]
    while !isempty(queue)
        pid = popfirst!(queue)
        push!(tree, pid)
        append!(queue, get(children, pid, Int[]))
    end
    return tree
end

"""
    _sample_process_tree(root_pid, clk_tck) -> (; n, rss_mib, hwm_mib, cpu_time_s)

Memory and cumulative CPU time summed over the live process tree. Per-process
`utime + stime` covers all OpenMP threads of that process; MPI ranks are
separate processes and are summed through the tree.
"""
function _sample_process_tree(root_pid::Integer, clk_tck::Integer)
    pids = _process_tree(root_pid)
    rss = 0
    hwm = 0
    ticks = 0
    for pid in pids
        r, h = _read_proc_status_memory(pid)
        rss += r
        hwm += h
        record = _read_proc_file("/proc/$pid/stat")
        record === nothing && continue
        _, ut, st = _parse_proc_stat(record)
        ticks += ut + st
    end
    return (
        n = length(pids),
        rss_mib = rss / 1024,
        hwm_mib = hwm / 1024,
        cpu_time_s = ticks / clk_tck,
    )
end

"""System 1-minute load average from `/proc/loadavg`; `NaN` when unavailable."""
function _load_average_1min()::Float64
    content = _read_proc_file("/proc/loadavg")
    content === nothing && return NaN
    v = tryparse(Float64, first(split(content)))
    return v === nothing ? NaN : v
end

# ---------------------------------------------------------------------------
# GPU query
# ---------------------------------------------------------------------------

"""
    _parse_nvidia_smi(output) -> (util_pct, mem_util_pct, mem_used_mib, power_w, temp_c)

Aggregate the per-device lines of
`nvidia-smi --query-gpu=$(_NVIDIA_SMI_QUERY) --format=csv,noheader,nounits`:
mean utilizations, summed memory and power, maximum temperature. Fields the
driver reports as `[N/A]` become `NaN`; an empty output gives all `NaN`.
"""
function _parse_nvidia_smi(output::AbstractString)::NTuple{5,Float64}
    rows = NTuple{5,Float64}[]
    for line in eachline(IOBuffer(String(output)))
        isempty(strip(line)) && continue
        fields = strip.(split(line, ','))
        length(fields) == 5 || throw(ArgumentError("unexpected nvidia-smi record: $line"))
        push!(rows, ntuple(i -> something(tryparse(Float64, fields[i]), NaN), 5))
    end
    isempty(rows) && return (NaN, NaN, NaN, NaN, NaN)
    n = length(rows)
    return (
        sum(r[1] for r in rows) / n,
        sum(r[2] for r in rows) / n,
        sum(r[3] for r in rows),
        sum(r[4] for r in rows),
        maximum(r[5] for r in rows),
    )
end

"""One `nvidia-smi` query, aggregated by [`_parse_nvidia_smi`](@ref); throws when the tool fails."""
function _query_nvidia_smi()::NTuple{5,Float64}
    output =
        read(`nvidia-smi --query-gpu=$(_NVIDIA_SMI_QUERY) --format=csv,noheader,nounits`, String)
    return _parse_nvidia_smi(output)
end

# ---------------------------------------------------------------------------
# Sampler
# ---------------------------------------------------------------------------

"""
    TelemetryMonitor

Handle of the asynchronous sampler started by [`_start_telemetry`](@ref):
the root PID, launch time, sampling interval, CSV path, the collected
samples, and the task state. `gpu_probe` flips to `false` after the first
failed `nvidia-smi` call so one missing tool cannot flood the run log.
"""
mutable struct TelemetryMonitor
    const pid::Int
    const t_start::Float64
    const interval::Float64
    const csv_path::String
    const samples::Vector{TelemetrySample}
    const clk_tck::Int
    gpu_probe::Bool
    stop_requested::Bool
    task::Union{Nothing,Task}
end

"""
    _take_sample!(mon::TelemetryMonitor) -> TelemetrySample

Sample the process tree, the load average, and (when enabled) the GPU; the
CPU-core rate is the CPU-time difference to the previous sample over the
wall-clock difference, `NaN` for the first sample.
"""
function _take_sample!(mon::TelemetryMonitor)::TelemetrySample
    elapsed = time() - mon.t_start
    tree = _sample_process_tree(mon.pid, mon.clk_tck)
    cores_busy = if isempty(mon.samples)
        NaN
    else
        prev = mon.samples[end]
        Δt = elapsed - prev.elapsed_s
        Δt > 0 ? max(tree.cpu_time_s - prev.cpu_time_s, 0.0) / Δt : NaN
    end
    gpu = (NaN, NaN, NaN, NaN, NaN)
    if mon.gpu_probe
        gpu = try
            _query_nvidia_smi()
        catch e
            e isa Union{ProcessFailedException,Base.IOError,ArgumentError} || rethrow()
            @warn "GPU telemetry disabled: nvidia-smi query failed" exception = e
            mon.gpu_probe = false
            (NaN, NaN, NaN, NaN, NaN)
        end
    end
    sample = TelemetrySample(
        elapsed,
        tree.n,
        tree.rss_mib,
        tree.hwm_mib,
        tree.cpu_time_s,
        cores_busy,
        _load_average_1min(),
        gpu...,
    )
    push!(mon.samples, sample)
    return sample
end

function _write_sample(io::IO, s::TelemetrySample)
    vals = (getfield(s, k) for k in fieldnames(TelemetrySample))
    println(io, join((v isa Integer ? string(v) : @sprintf("%.3f", v) for v in vals), ","))
    return nothing
end

"""Sleep in short slices until `deadline` or until the monitor is asked to stop."""
function _sleep_until(mon::TelemetryMonitor, deadline::Float64)
    while !mon.stop_requested
        remaining = deadline - time()
        remaining ≤ 0 && return nothing
        sleep(min(remaining, 0.25))
    end
    return nothing
end

"""
    _start_telemetry(run_dir, pid, t_start; interval, gpu_probe, csv_name = "telemetry.csv",
                     interrupt_to = nothing) -> TelemetryMonitor

Start the asynchronous sampler for the process tree rooted at `pid`,
appending one row per `interval` seconds to `<run_dir>/telemetry.csv`
(header = the [`TelemetrySample`](@ref) field names). An existing CSV is
backed up, never overwritten. The task shares the
main thread and yields between samples, so it interleaves with the process
wait and the live monitor without extra threads. An `InterruptException`
landing in the task is handed to `interrupt_to`, the task waiting on the
engine ([`_forwarding_interrupts`](@ref)).
"""
function _start_telemetry(
    run_dir::AbstractString,
    pid::Integer,
    t_start::Float64;
    interval::Float64,
    gpu_probe::Bool,
    csv_name::AbstractString = "telemetry.csv",
    interrupt_to::Union{Nothing,Task} = nothing,
)::TelemetryMonitor
    interval > 0 || throw(ArgumentError("telemetry interval must be > 0 s; got $interval"))
    csv_path = joinpath(run_dir, csv_name)
    mon = TelemetryMonitor(
        Int(pid),
        t_start,
        interval,
        csv_path,
        TelemetrySample[],
        _clock_ticks_per_second(),
        gpu_probe,
        false,
        nothing,
    )
    _backup_existing(csv_path)
    io = open(csv_path, "w")
    println(io, join(string.(fieldnames(TelemetrySample)), ","))
    flush(io)
    mon.task = _forwarding_interrupts(interrupt_to) do
        try
            while !mon.stop_requested
                t_next = time() + mon.interval
                _write_sample(io, _take_sample!(mon))
                flush(io)
                _sleep_until(mon, t_next)
            end
        finally
            close(io)
        end
    end
    return mon
end

"""
    _stop_telemetry!(mon::TelemetryMonitor) -> Dict{String,Any}

Stop the sampler and return its summary table. A sampler failure is reported
as a warning, never as an error: telemetry is auxiliary to the run.
"""
function _stop_telemetry!(mon::TelemetryMonitor)::Dict{String,Any}
    mon.stop_requested = true
    if mon.task !== nothing
        try
            wait(mon.task)
        catch e
            @warn "Telemetry sampler failed; summary uses the samples collected so far" exception =
                (e, catch_backtrace())
        end
    end
    return _telemetry_summary(mon)
end

_mean(v) = sum(v) / length(v)

"""
    _telemetry_summary(mon::TelemetryMonitor) -> Dict{String,Any}

Summary statistics of the collected samples for `RUN_INFO.toml [telemetry]`:
interval, sample count, CSV name, peak RSS (the larger of the sampled RSS
and the kernel high-water mark), mean and peak CPU cores busy, peak load,
and GPU means/peaks when any GPU sample exists. Quantities without data are
omitted rather than written as `NaN`.
"""
function _telemetry_summary(mon::TelemetryMonitor)::Dict{String,Any}
    s = mon.samples
    d = Dict{String,Any}(
        "sampling_interval_s" => mon.interval,
        "samples" => length(s),
        "csv" => basename(mon.csv_path),
    )
    isempty(s) && return d
    d["peak_rss_mib"] =
        round(max(maximum(x.rss_mib for x in s), maximum(x.hwm_mib for x in s)); digits = 1)
    busy = filter(!isnan, [x.cores_busy for x in s])
    if !isempty(busy)
        d["mean_cores_busy"] = round(_mean(busy); digits = 2)
        d["peak_cores_busy"] = round(maximum(busy); digits = 2)
    end
    load = filter(!isnan, [x.load_1min for x in s])
    isempty(load) || (d["peak_load_1min"] = round(maximum(load); digits = 2))
    gpu = filter(x -> !isnan(x.gpu_util_pct), s)
    if !isempty(gpu)
        d["mean_gpu_util_pct"] = round(_mean([x.gpu_util_pct for x in gpu]); digits = 1)
        d["peak_gpu_util_pct"] = round(maximum(x.gpu_util_pct for x in gpu); digits = 1)
        mem = filter(!isnan, [x.gpu_mem_used_mib for x in gpu])
        isempty(mem) || (d["peak_gpu_mem_used_mib"] = round(maximum(mem); digits = 1))
        power = filter(!isnan, [x.gpu_power_w for x in gpu])
        isempty(power) || (d["mean_gpu_power_w"] = round(_mean(power); digits = 1))
        temp = filter(!isnan, [x.gpu_temp_c for x in gpu])
        isempty(temp) || (d["peak_gpu_temp_c"] = round(maximum(temp); digits = 1))
    end
    return d
end

"""
    _finish_telemetry(mon, cpu_before, cpu_after, elapsed, n_threads) -> Dict{String,Any}

Merge the sampler summary (`mon` may be `nothing` when sampling is disabled)
with the exact child CPU accounting: user and system seconds from the
`getrusage` deltas, `threads_total` (effective OpenMP threads × MPI ranks),
and the CPU efficiency `(user + system) / (elapsed × threads_total)`, the
fraction of the reserved CPU capacity the backend actually used.
"""
function _finish_telemetry(
    mon::Union{Nothing,TelemetryMonitor},
    cpu_before::Tuple{Float64,Float64},
    cpu_after::Tuple{Float64,Float64},
    elapsed::Float64,
    n_threads::Integer,
)::Dict{String,Any}
    d =
        mon === nothing ? Dict{String,Any}("sampling_interval_s" => 0.0, "samples" => 0) :
        _stop_telemetry!(mon)
    user = cpu_after[1] - cpu_before[1]
    system = cpu_after[2] - cpu_before[2]
    d["threads_total"] = Int(n_threads)
    if !isnan(user) && !isnan(system)
        d["cpu_user_s"] = round(user; digits = 2)
        d["cpu_system_s"] = round(system; digits = 2)
        if elapsed > 0 && n_threads > 0
            d["cpu_efficiency"] = round((user + system) / (elapsed * n_threads); digits = 3)
        end
    end
    return d
end

# ---------------------------------------------------------------------------
# Backend-reported performance (stdout timing table, force-kernel Gflops)
# ---------------------------------------------------------------------------

"""
    _timing_key(name) -> String

TOML-safe key for a column of the backend timing table: lower case, runs of
non-alphanumerics collapsed to `_`, no leading or trailing `_`
(`"Reg.GPU.S"` → `"reg_gpu_s"`, `"KS.Init.B"` → `"ks_init_b"`).
"""
function _timing_key(name::AbstractString)::String
    return strip(replace(lowercase(name), r"[^a-z0-9]+" => "_"), '_')
end

"""
    _backend_timing_table(stdout_path) -> Union{Nothing,Dict{String,Any}}

The last timing table the backend printed to stdout (`adjust.F`, once per
`DTADJ`): a header line beginning with `rank PE N Total …` followed by one
data line. Values are cumulative seconds of CPU time per code section
(regular/irregular force, KS, adjust, output, communication, …), so the last
table is the whole-run breakdown. Keys are the column names through
[`_timing_key`](@ref); the rank and PE columns are dropped, `N` becomes
`n_particles`. `nothing` when the file or the table is absent.
"""
function _backend_timing_table(stdout_path::AbstractString)::Union{Nothing,Dict{String,Any}}
    isfile(stdout_path) || return nothing
    header = nothing
    data = nothing
    lines = readlines(stdout_path)
    for i in 1:(length(lines) - 1)
        tokens = split(lines[i])
        if length(tokens) ≥ 4 && tokens[1] == "rank" && tokens[2] == "PE" && tokens[3] == "N"
            values = split(lines[i + 1])
            if length(values) == length(tokens)
                header = tokens
                data = values
            end
        end
    end
    header === nothing && return nothing
    table = Dict{String,Any}()
    for (name, value) in zip(header, data)
        name in ("rank", "PE") && continue
        key = name == "N" ? "n_particles" : _timing_key(name)
        v = tryparse(Int, value)
        table[key] = v === nothing ? something(tryparse(Float64, value), NaN) : v
    end
    return table
end

"""
    _force_kernel_gflops(stderr_path) -> Union{Nothing,Dict{String,Any}}

Regular-force kernel throughput the backend's AVX/SSE or GPU profiles print
to stderr once per `DTADJ` (`Perf.(Gflops) <value>`, reset after each
print). Returns the sample count, mean, and peak in Gflops and the kernel
labels the lines carry (`"GPU Reg.F"`, `"AVX Reg.F"`, …; several joined
with `, `), or `nothing` when the file carries no such line.
"""
function _force_kernel_gflops(stderr_path::AbstractString)::Union{Nothing,Dict{String,Any}}
    isfile(stderr_path) || return nothing
    values = Float64[]
    kernels = String[]
    for line in eachline(stderr_path)
        m = match(r"Perf\.\(Gflops\)\s+([0-9.eE+-]+)", line)
        m === nothing && continue
        v = tryparse(Float64, m.captures[1])
        v === nothing && continue
        push!(values, v)
        k = match(r"\[R\.\d+\s+(.+?)\s*\]", line)
        k === nothing && continue
        label = String(k.captures[1])
        label in kernels || push!(kernels, label)
    end
    isempty(values) && return nothing
    d = Dict{String,Any}(
        "samples" => length(values),
        "mean" => round(_mean(values); digits = 2),
        "peak" => round(maximum(values); digits = 2),
    )
    isempty(kernels) || (d["kernel"] = join(kernels, ", "))
    return d
end

"""
    _backend_performance(stdout_path, stderr_path) -> Dict{String,Any}

Backend-reported performance for the `[telemetry]` table: the
`backend_timing` sub-table from [`_backend_timing_table`](@ref) and the
`force_kernel_gflops` sub-table from [`_force_kernel_gflops`](@ref), each
present only when the corresponding output exists.
"""
function _backend_performance(stdout_path::AbstractString, stderr_path::AbstractString)
    d = Dict{String,Any}()
    timing = _backend_timing_table(stdout_path)
    timing === nothing || (d["backend_timing"] = timing)
    gflops = _force_kernel_gflops(stderr_path)
    gflops === nothing || (d["force_kernel_gflops"] = gflops)
    return d
end

# ---------------------------------------------------------------------------
# Readers for the sampler's CSVs
# ---------------------------------------------------------------------------

"""
    read_telemetry(path) -> Vector{TelemetrySample}

Read a telemetry CSV written by the sampler (header = the field names of
[`TelemetrySample`](@ref), matched by name in any order). Empty cells and
`NaN` read as `NaN`; `n_processes` is rounded to an integer. Throws an
`ArgumentError` when a field of the sample is missing from the header or a
row has the wrong number of cells.
"""
function read_telemetry(path::AbstractString)::Vector{TelemetrySample}
    lines = filter(!isempty, strip.(readlines(path)))
    isempty(lines) && return TelemetrySample[]
    header = String.(strip.(split(lines[1], ',')))
    names = string.(fieldnames(TelemetrySample))
    col = Dict(h => i for (i, h) in enumerate(header))
    missing_names = filter(n -> !haskey(col, n), names)
    isempty(missing_names) || throw(
        ArgumentError(
            "telemetry CSV $path lacks the column(s) $(join(missing_names, ", ")); " *
            "header: $(join(header, ", "))",
        ),
    )
    samples = TelemetrySample[]
    for line in lines[2:end]
        cells = strip.(split(line, ','))
        length(cells) == length(header) || throw(
            ArgumentError(
                "telemetry CSV $path: a row has $(length(cells)) cells for $(length(header)) columns",
            ),
        )
        value(n) = something(tryparse(Float64, cells[col[n]]), NaN)
        push!(
            samples,
            TelemetrySample(
                value("elapsed_s"),
                round(Int, something(tryparse(Float64, cells[col["n_processes"]]), 0.0)),
                value("rss_mib"),
                value("hwm_mib"),
                value("cpu_time_s"),
                value("cores_busy"),
                value("load_1min"),
                value("gpu_util_pct"),
                value("gpu_mem_util_pct"),
                value("gpu_mem_used_mib"),
                value("gpu_power_w"),
                value("gpu_temp_c"),
            ),
        )
    end
    return samples
end

"""
    read_run_telemetry(run_dir) -> Vector{TelemetrySample}

Every telemetry segment of a run (`telemetry.csv`, `telemetry_2.csv`, …,
one per launch) concatenated in segment order, the elapsed time of each
segment offset by the last elapsed time of the previous one, since restarts
run one after another. Empty when the run directory holds no telemetry.
"""
function read_run_telemetry(run_dir::AbstractString)::Vector{TelemetrySample}
    isdir(run_dir) || return TelemetrySample[]
    files = filter(f -> match(r"^telemetry(_\d+)?\.csv$", f) !== nothing, readdir(run_dir))
    segment_index(f) = f == "telemetry.csv" ? 1 : parse(Int, match(r"_(\d+)\.csv$", f).captures[1])
    sort!(files; by = segment_index)
    samples = TelemetrySample[]
    offset = 0.0
    for f in files
        seg = read_telemetry(joinpath(run_dir, f))
        isempty(seg) && continue
        for s in seg
            push!(
                samples,
                TelemetrySample(
                    s.elapsed_s + offset,
                    s.n_processes,
                    s.rss_mib,
                    s.hwm_mib,
                    s.cpu_time_s,
                    s.cores_busy,
                    s.load_1min,
                    s.gpu_util_pct,
                    s.gpu_mem_util_pct,
                    s.gpu_mem_used_mib,
                    s.gpu_power_w,
                    s.gpu_temp_c,
                ),
            )
        end
        offset = samples[end].elapsed_s
    end
    return samples
end
