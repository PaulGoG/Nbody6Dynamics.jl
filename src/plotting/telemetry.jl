# =============================================================================
# Run telemetry: readers for the sampler's CSVs and the time-series figure
# =============================================================================

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

"""Colours of the telemetry series (Okabe–Ito): cores busy, resident memory, GPU utilisation."""
const _TELEMETRY_COLORS = (
    cores = _OKABE_ITO[3],   # bluish green
    memory = _OKABE_ITO[2],  # orange
    gpu = _OKABE_ITO[6],     # vermillion
)

"""Wall-clock axis of a telemetry series: seconds under two minutes, minutes under two hours, hours beyond; returns `(values, label)`."""
function _telemetry_time_axis(elapsed_s::AbstractVector{<:Real})
    t_end = maximum(elapsed_s)
    if t_end < 120
        return (Float64.(elapsed_s), L"t \; [\mathrm{s}]")
    elseif t_end < 7200
        return (Float64.(elapsed_s) ./ 60, L"t \; [\mathrm{min}]")
    end
    return (Float64.(elapsed_s) ./ 3600, L"t \; [\mathrm{h}]")
end

"""Draw the finite part of `y` against `t` on `ax` and register a legend entry; returns whether anything was drawn."""
function _telemetry_series!(
    ax::Axis,
    t::AbstractVector{<:Real},
    y::AbstractVector{<:Real},
    elements::Vector{Union{LineElement,MarkerElement}},
    labels::Vector{AbstractString},
    label::AbstractString;
    color,
    linestyle = :solid,
)::Bool
    ok = .!isnan.(y)
    any(ok) || return false
    lines!(ax, t[ok], y[ok]; color = color, linestyle = linestyle)
    push!(elements, LineElement(; color = color, linestyle = linestyle))
    push!(labels, label)
    return true
end

"""Canvas of the telemetry stack: the column width by half a single-panel height per row plus a shared margin, so three panels stay a compact column figure."""
_telemetry_canvas(cfg::VisualizationConfig, nrows::Int) =
    (_figsize_px(cfg)[1], round(Int, _figsize_px(cfg)[2] * (0.5 + 0.55 * nrows)))

"""
    plot_telemetry(samples, cfg; filename = "telemetry") -> Union{Nothing,String}

Run telemetry against wall-clock time as stacked panels sharing the time
axis: cores busy, with the host's 1-minute load average dotted on a twin
axis (grey tick labels); resident memory (RSS solid, high-water mark
dashed); and, when any sample carries GPU data, GPU utilisation with the
memory utilisation dashed. The mean cores busy, the peak RSS and the mean
GPU utilisation appear as the legend entries' quantitative takeaways.
Returns the output path, or `nothing` with fewer than two samples.
"""
function plot_telemetry(
    samples::Vector{TelemetrySample},
    cfg::VisualizationConfig;
    filename::AbstractString = "telemetry",
)
    length(samples) ≥ 2 ||
        (@warn "Fewer than two telemetry samples; nothing to plot"; return nothing)
    t, tlabel = _telemetry_time_axis([s.elapsed_s for s in samples])
    ttk = _time_ticks(first(t), last(t))
    gpu = any(s -> !isnan(s.gpu_util_pct), samples)
    nrows = gpu ? 3 : 2
    fig = Figure(; size = _telemetry_canvas(cfg, nrows))
    elements = Union{LineElement,MarkerElement}[]
    labels = AbstractString[]
    finite(v) = filter(!isnan, v)

    # --- Panel 1: cores busy (left axis), host load average (twin axis) ---
    cores = [s.cores_busy for s in samples]
    load = [s.load_1min for s in samples]
    ax1 =
        Axis(fig[1, 1]; ylabel = L"\mathrm{Cores\;busy}", xticks = ttk, xticklabelsvisible = false)
    label_cores =
        isempty(finite(cores)) ? L"\mathrm{Cores\;busy}" :
        latexstring(@sprintf("\\mathrm{Cores\\;busy\\;(mean\\;%.2f)}", _mean(finite(cores))))
    _telemetry_series!(
        ax1,
        t,
        cores,
        elements,
        labels,
        label_cores;
        color = _TELEMETRY_COLORS.cores,
    )
    axes = Axis[ax1]
    if !isempty(finite(load))
        ax1b = Axis(
            fig[1, 1];
            yaxisposition = :right,
            ylabel = L"\mathrm{Load\;(1\;min)}",
            ylabelcolor = :gray40,
            yticklabelcolor = :gray40,
            ytickcolor = :gray40,
            xgridvisible = false,
            ygridvisible = false,
        )
        hidexdecorations!(ax1b)
        hidespines!(ax1b)
        _telemetry_series!(
            ax1b,
            t,
            load,
            elements,
            labels,
            L"\mathrm{Load\;(1\;min,\;host)}";
            color = :gray40,
            linestyle = :dot,
        )
        push!(axes, ax1b)
    end

    # --- Panel 2: resident memory ---
    rss = [s.rss_mib for s in samples]
    ax2 = Axis(
        fig[2, 1];
        ylabel = L"\mathrm{Memory} \; [\mathrm{MiB}]",
        xlabel = gpu ? "" : tlabel,
        xticks = ttk,
        xticklabelsvisible = !gpu,
    )
    label_rss =
        isempty(finite(rss)) ? L"\mathrm{RSS}" :
        latexstring(@sprintf("\\mathrm{RSS\\;(peak\\;%.3g\\;MiB)}", maximum(finite(rss))))
    _telemetry_series!(ax2, t, rss, elements, labels, label_rss; color = _TELEMETRY_COLORS.memory)
    _telemetry_series!(
        ax2,
        t,
        [s.hwm_mib for s in samples],
        elements,
        labels,
        L"\mathrm{High\;water\;mark}";
        color = _TELEMETRY_COLORS.memory,
        linestyle = :dash,
    )
    push!(axes, ax2)

    # --- Panel 3: GPU utilisation ---
    if gpu
        util = [s.gpu_util_pct for s in samples]
        ax3 = Axis(
            fig[3, 1];
            ylabel = L"\mathrm{GPU} \; [\%]",
            xlabel = tlabel,
            xticks = ttk,
            limits = (nothing, nothing, 0, 105),
        )
        _telemetry_series!(
            ax3,
            t,
            util,
            elements,
            labels,
            latexstring(
                @sprintf("\\mathrm{GPU\\;utilisation\\;(mean\\;%.1f\\,\\%%)}", _mean(finite(util)))
            );
            color = _TELEMETRY_COLORS.gpu,
        )
        _telemetry_series!(
            ax3,
            t,
            [s.gpu_mem_util_pct for s in samples],
            elements,
            labels,
            L"\mathrm{GPU\;memory\;utilisation}";
            color = _TELEMETRY_COLORS.gpu,
            linestyle = :dash,
        )
        push!(axes, ax3)
    end

    isempty(elements) || _top_legend!(fig, elements, labels; nbanks = length(elements) > 3 ? 2 : 1)
    linkxaxes!(axes...)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    return _save_fig(cfg, filename, fig)
end
