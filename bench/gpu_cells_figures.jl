#!/usr/bin/env julia
# =============================================================================
# Figures of the device benchmark cells: the CSV rows of bench/gpu_cells.jl,
# gathered from every host, drawn as two publication figures.
#
#   gpu_cells_scaling    the node with two devices: strong scaling at the
#                        largest N, weak scaling at a fixed N per device, and
#                        the run time against the OpenMP thread count;
#   gpu_cells_hardware   one pair of bars per card: the engine wall time with a
#                        tick at the regular-force share, and the
#                        regular-force kernel rate as a fraction of the card's
#                        FP32 peak.
#
#   julia bench/gpu_cells_figures.jl [--data=<dir>] [--cards=<toml>] [--out=<dir>]
#                                    [--synthetic]
#
# --data       searched recursively for gpu_cells_*.csv and fp32_peak_*.toml
#              (default bench/results); subdirectories named "synthetic" are
#              passed over.
# --cards      [[card]] tables with match (a substring of the CSV gpu column),
#              label, class ("datacenter" | "consumer"), fp32_peak_tflops (the
#              datasheet value) and an optional host whose fp32_peak_<host>.toml
#              under --data replaces the datasheet peak. File order is axis
#              order within each class. Without it the cards are taken from the
#              CSV and the kernel-rate panel is omitted.
# --out        output directory (default bench/results/figures).
# --synthetic  writes a fabricated data set to <out>/synthetic/ and draws from
#              it, so the layout can be checked before any run exists:
#
#   julia bench/gpu_cells_figures.jl --synthetic --out=<dir>
#
# The backend timing table is summed over the OpenMP threads; the regular and
# irregular force times drawn are therefore its shares of the wall time. The
# script runs in the scripts/ environment (the package with CairoMakie); the
# bench environment carries no Makie.
# =============================================================================

include(joinpath(@__DIR__, "..", "scripts", "activate.jl"))

using CairoMakie, Nbody6Dynamics, TOML, Printf, Statistics, Dates

const BENCH = @__DIR__
const PROJ = normpath(joinpath(BENCH, ".."))

# The palette and the edge colour live in the package's Makie extension.
const MAKIE_EXT = Base.get_extension(Nbody6Dynamics, :Nbody6DynamicsMakieExt)
MAKIE_EXT === nothing && error("the Makie extension of Nbody6Dynamics did not load")
const PALETTE = MAKIE_EXT._OKABE_ITO

const CELLS = ("single", "dual")
const CELL_COLOUR = Dict("single" => PALETTE[1], "dual" => PALETTE[2])
const CELL_LABEL = Dict("single" => "Single cell", "dual" => "Dual cell")
const FONT = 26                          # base font size of publication_theme()
const ANNOTATION_SIZE = 0.8 * FONT
const VALUE_SIZE = 0.75 * FONT
const GUIDE_WIDTH = 1.5
const GUIDE_TEXT = :grey40
const BAR_WIDTH = 0.36
const BAR_OFFSET = Dict("single" => -0.2, "dual" => 0.2)
const PREFERRED_THREADS = 8              # thread count of the hardware comparison

const COLUMNS = split(
    "cell,N_requested,N_total,variant,gpu_list,mpi_ranks,omp_threads,tcrit,status,wall_s," *
    "backend_total_s,backend_reg_s,backend_irr_s,backend_ks_s,backend_mdot_s,backend_barr_s," *
    "gflops_mean,gflops_peak,kernel,gpu_util_mean,gpu_util_peak,gpu_power_mean_w," *
    "gpu_mem_peak_mib,peak_rss_mib,cpu_efficiency,gpu_devices,machine,gpu,run_id",
    ',',
)

const USAGE = "usage: gpu_cells_figures.jl [--data=<dir>] [--cards=<toml>] [--out=<dir>] [--synthetic]"

"""Darker same-hue edge of a filled marker or bar."""
edge(colour) = MAKIE_EXT._band_edge(colour)

# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

"""Value of `--name=value` in `arg`, or `nothing` when `arg` is another option."""
option(arg, name) = startswith(arg, "--$name=") ? arg[(length(name) + 4):end] : nothing

"""
    parse_options(args) -> NamedTuple

Options `(data, cards, out, synthetic)` of the command line, paths made
absolute; an unknown option prints the usage line and exits with status 2.
"""
function parse_options(args)
    data = joinpath(PROJ, "bench", "results")
    cards = nothing
    out = joinpath(PROJ, "bench", "results", "figures")
    synthetic = false
    for a in args
        if a == "--synthetic"
            synthetic = true
        elseif (v = option(a, "data")) !== nothing
            data = abspath(expanduser(v))
        elseif (v = option(a, "cards")) !== nothing
            cards = abspath(expanduser(v))
        elseif (v = option(a, "out")) !== nothing
            out = abspath(expanduser(v))
        else
            println(stderr, USAGE)
            exit(2)
        end
    end
    return (data = data, cards = cards, out = out, synthetic = synthetic)
end

# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

"""
    split_csv(line) -> Vector{String}

Fields of one CSV line; a double-quoted field keeps its commas.
"""
function split_csv(line::AbstractString)
    fields = String[]
    buf = IOBuffer()
    quoted = false
    for c in line
        if c == '"'
            quoted = !quoted
        elseif c == ',' && !quoted
            push!(fields, String(take!(buf)))
        else
            print(buf, c)
        end
    end
    push!(fields, String(take!(buf)))
    return fields
end

"""
    find_inputs(dir) -> (csvs, peaks)

The `gpu_cells_*.csv` and `fp32_peak_*.toml` files below `dir`, sorted. Files
in a subdirectory named `synthetic` are passed over, so a rehearsal data set
written under the default output directory never mixes with measured rows.
"""
function find_inputs(dir::AbstractString)
    csvs, peaks = String[], String[]
    n_skipped = 0
    for (root, _, files) in walkdir(dir), f in files
        is_csv = startswith(f, "gpu_cells_") && endswith(f, ".csv")
        is_peak = startswith(f, "fp32_peak_") && endswith(f, ".toml")
        (is_csv || is_peak) || continue
        if "synthetic" in splitpath(relpath(root, dir))
            n_skipped += 1
            continue
        end
        push!(is_csv ? csvs : peaks, joinpath(root, f))
    end
    n_skipped > 0 && println("passed over $n_skipped file(s) in synthetic/ below $dir")
    return sort(csvs), sort(peaks)
end

"""
    share(wall, part, total) -> Float64

Share `wall × part / total` of the wall time taken by one entry of the backend
timing table, which is summed over the OpenMP threads; `NaN` without a total.
"""
share(wall, part, total) = (isnan(total) || total == 0) ? NaN : wall * part / total

"""Float of a CSV field; an empty field reads as `NaN`."""
number(s::AbstractString) = isempty(s) ? NaN : parse(Float64, s)

"""
    read_rows(paths) -> Vector{NamedTuple}

Completed device rows (`status == "completed"`, `variant == "gpu"`) of the
`gpu_cells_*.csv` files `paths`, with the device count `D` and the wall-time
shares `t_reg`, `t_irr`, `t_ks` of the backend timing table.
"""
function read_rows(paths)
    rows = NamedTuple[]
    for path in paths
        lines = readlines(path)
        isempty(lines) && continue
        header = split_csv(strip(lines[1]))
        col = Dict(name => i for (i, name) in enumerate(header))
        for name in COLUMNS
            haskey(col, name) || error("$path: column \"$name\" missing from the header")
        end
        for (k, line) in enumerate(lines[2:end])
            isempty(strip(line)) && continue
            f = split_csv(rstrip(line))
            length(f) == length(header) ||
                error("$path, line $(k + 1): $(length(f)) fields, the header has $(length(header))")
            field = name -> f[col[name]]
            (field("status") == "completed" && field("variant") == "gpu") || continue
            wall = number(field("wall_s"))
            total = number(field("backend_total_s"))
            D = parse(Int, field("gpu_devices"))
            D == 0 && (D = length(split(field("gpu_list"))))
            push!(
                rows,
                (
                    cell = field("cell"),
                    N_requested = parse(Int, field("N_requested")),
                    N_total = parse(Int, field("N_total")),
                    omp_threads = parse(Int, field("omp_threads")),
                    D = D,
                    wall_s = wall,
                    t_reg = share(wall, number(field("backend_reg_s")), total),
                    t_irr = share(wall, number(field("backend_irr_s")), total),
                    t_ks = share(wall, number(field("backend_ks_s")), total),
                    gflops_mean = number(field("gflops_mean")),
                    gpu_util_mean = number(field("gpu_util_mean")),
                    gpu_power_mean_w = number(field("gpu_power_mean_w")),
                    machine = field("machine"),
                    gpu = field("gpu"),
                    run_id = field("run_id"),
                ),
            )
        end
    end
    return rows
end

"""
    read_cards(path, peak_files) -> Vector{NamedTuple}

Cards `(match, label, class, peak, source)` of the `[[card]]` tables in
`path`, datacenter cards first, each class in file order. A card with a
`host` takes its FP32 peak from the `fp32_peak_*.toml` among `peak_files`
whose `host` equals it, unless that file names another device.
"""
function read_cards(path, peak_files)
    doc = TOML.parsefile(path)
    get(doc, "card", nothing) isa AbstractVector || error("$path: no [[card]] tables")
    measured = Dict{String,Tuple{Float64,String,String}}()
    for file in peak_files
        d = TOML.parsefile(file)
        tflops = get(get(d, "fp32", Dict()), "tflops", nothing)
        if !haskey(d, "host") || tflops === nothing
            println("$(basename(file)): no host or fp32.tflops, ignored")
            continue
        end
        measured[d["host"]] = (Float64(tflops), file, get(d, "device", ""))
    end
    cards = map(enumerate(doc["card"])) do (i, c)
        for key in ("match", "label", "class", "fp32_peak_tflops")
            haskey(c, key) || error("$path: card $i lacks the key \"$key\"")
        end
        c["class"] in ("datacenter", "consumer") ||
            error("$path: card $i: class must be one of \"datacenter\" | \"consumer\"")
        peak = Float64(c["fp32_peak_tflops"])
        source = "datasheet"
        host = get(c, "host", "")
        if haskey(measured, host)
            tflops, file, device = measured[host]
            if isempty(device) || occursin(c["match"], device)
                peak = tflops
                source = "measured $(basename(file))"
            else
                println(
                    "$(basename(file)) measured \"$device\", not \"$(c["match"])\": ",
                    "datasheet peak kept",
                )
            end
        end
        (
            match = String(c["match"]),
            label = String(c["label"]),
            class = String(c["class"]),
            peak = peak,
            source = source,
        )
    end
    return vcat(
        filter(c -> c.class == "datacenter", cards),
        filter(c -> c.class == "consumer", cards),
    )
end

"""
    infer_cards(rows) -> Vector{NamedTuple}

Cards taken from the distinct `gpu` values of `rows`: the device name before
the first comma as label, class "consumer", no FP32 peak.
"""
function infer_cards(rows)
    gpus = sort(unique(r.gpu for r in rows if !isempty(r.gpu)))
    return [
        (
            match = g,
            label = String(strip(first(split(g, ',')))),
            class = "consumer",
            peak = NaN,
            source = "none",
        ) for g in gpus
    ]
end

"""
    pick(rows, cell, N, D, threads) -> NamedTuple or nothing

The last row of `cell` at `N_requested == N` on `D` devices with `threads`
OpenMP threads, or `nothing`.
"""
function pick(rows, cell, N, D, threads)
    i = findlast(
        r -> r.cell == cell && r.N_requested == N && r.D == D && r.omp_threads == threads,
        rows,
    )
    return i === nothing ? nothing : rows[i]
end

# ---------------------------------------------------------------------------
# Formatting and drawing helpers
# ---------------------------------------------------------------------------

const SUPERSCRIPT = Dict(zip("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))

"""
    format_N(N) -> String

`N` in figure notation: `10⁶`, `5×10⁵`, `2.5×10⁵`; plain digits below 100.
"""
function format_N(N::Real)
    N < 100 && return string(round(Int, N))
    e = floor(Int, log10(N))
    m = N / exp10(e)
    if m ≥ 10 - 1e-9
        e += 1
        m /= 10
    end
    mantissa = isapprox(m, round(m); atol = 1e-9) ? string(round(Int, m)) : @sprintf("%.3g", m)
    power = "10" * map(c -> SUPERSCRIPT[c], string(e))
    return mantissa == "1" ? power : mantissa * "×" * power
end

"""Wall time label: three significant digits, an integer from 100 up."""
format_wall(y) = y ≥ 100 ? string(round(Int, y)) : @sprintf("%.3g", y)

"""Kernel-rate label: one decimal."""
format_rate(y) = @sprintf("%.1f", y)

"""Tick label of a log-axis value: integers without exponent, fractions as decimals."""
plain_decimal(v) = v ≥ 1 ? string(round(Int, v)) : string(round(v; sigdigits = 1))

"""
    decade_ticks(lo, hi) -> (values, labels)

Log-axis ticks at the decades between `lo` and `hi` with the 2× and 5×
intermediates, labelled as plain decimals.
"""
function decade_ticks(lo, hi)
    vals = Float64[]
    for e in floor(Int, log10(lo)):ceil(Int, log10(hi)), m in (1, 2, 5)
        v = m * exp10(e)
        lo * (1 - 1e-9) ≤ v ≤ hi * (1 + 1e-9) && push!(vals, v)
    end
    return (vals, plain_decimal.(vals))
end

"""
    series!(ax, x, y, colour; hollow = false, linestyle = :solid)

One series: a line with circular markers, filled with a darker edge or hollow.
"""
function series!(ax, x, y, colour; hollow::Bool = false, linestyle = :solid)
    lines!(ax, x, y; color = colour, linestyle = linestyle)
    if hollow
        scatter!(ax, x, y; color = :white, strokecolor = colour, strokewidth = 2)
    else
        scatter!(ax, x, y; color = colour, strokecolor = edge(colour), strokewidth = 1.5)
    end
    return nothing
end

"""
    series_element(colour; hollow = false, linestyle = :solid) -> Vector{LegendElement}

Legend glyph matching [`series!`](@ref): line and marker.
"""
function series_element(colour; hollow::Bool = false, linestyle = :solid)
    line = LineElement(; color = colour, linestyle = linestyle, linewidth = 3)
    marker = if hollow
        MarkerElement(;
            marker = :circle,
            markersize = 14,
            color = :white,
            strokecolor = colour,
            strokewidth = 2,
        )
    else
        MarkerElement(;
            marker = :circle,
            markersize = 14,
            color = colour,
            strokecolor = edge(colour),
            strokewidth = 1.5,
        )
    end
    return Makie.LegendElement[line, marker]
end

"""Annotation `s` in the top-left corner of `ax`."""
function corner_note!(ax, s)
    text!(
        ax,
        0.04,
        0.96;
        text = s,
        space = :relative,
        align = (:left, :top),
        fontsize = ANNOTATION_SIZE,
    )
    return nothing
end

"""Annotation `s` centred in `ax`, for a panel without data."""
function centre_note!(ax, s)
    text!(
        ax,
        0.5,
        0.5;
        text = s,
        space = :relative,
        align = (:center, :center),
        fontsize = ANNOTATION_SIZE,
    )
    return nothing
end

"""Save `fig` as `<name>.pdf` and `<name>.png` (4 px per unit) in `out`."""
function save_figure(fig, out, name)
    for ext in ("pdf", "png")
        path = joinpath(out, "$name.$ext")
        save(path, fig; px_per_unit = 4)
        println("saved: ", path)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Figure 1: scaling on the node with two devices
# ---------------------------------------------------------------------------

"""
    strong_panel!(pos, rows, cells, N, T, ids) -> Axis

Speed-up t(D = 1)/t(D) at `N` stars and `T` threads, for the regular-force
share and the wall time of each cell, against the ideal y = D; the run IDs
are appended to `ids`.
"""
function strong_panel!(pos, rows, cells, N, T, ids)
    ax = Axis(
        pos;
        xlabel = L"\mathrm{Devices}\;D",
        ylabel = "Speed-up",
        xticks = ([1, 2], ["1", "2"]),
        yticks = 1:0.5:2,
    )
    lines!(ax, [1, 2], [1, 2]; color = :grey, linestyle = :dash, linewidth = GUIDE_WIDTH)
    # Labelled at mid-line, on its upper-left side: the series lie on or
    # below the ideal, and the corner note owns the band above y = 2.
    text!(
        ax,
        1.5,
        1.5;
        text = "Ideal",
        align = (:right, :bottom),
        offset = (-8, 4),
        color = GUIDE_TEXT,
        fontsize = ANNOTATION_SIZE,
    )
    for cell in cells
        r1 = pick(rows, cell, N, 1, T)
        r2 = pick(rows, cell, N, 2, T)
        (r1 === nothing || r2 === nothing) && continue
        c = CELL_COLOUR[cell]
        series!(ax, [1, 2], [1, r1.wall_s / r2.wall_s], c; hollow = true)
        series!(ax, [1, 2], [1, r1.t_reg / r2.t_reg], c)
        push!(ids, "panel a: $cell D=1 $(r1.run_id), D=2 $(r2.run_id)")
    end
    xlims!(ax, 0.85, 2.15)
    ylims!(ax, 0.9, 2.3)
    corner_note!(ax, "Strong scaling, N = $(format_N(N)), $T threads")
    return ax
end

"""
    weak_panel!(pos, rows, cells, T, ids) -> Axis

Run time at N_w stars on one device and 2 N_w on two, `T` threads, for the
wall time and the regular-force share of each cell, with the ideal (constant)
run time dashed; N_w is the smallest N at which such a pair exists. The run
IDs are appended to `ids`.
"""
function weak_panel!(pos, rows, cells, T, ids)
    Ns = sort(unique(r.N_requested for r in rows))
    paired = (c, N) -> pick(rows, c, N, 1, T) !== nothing && pick(rows, c, 2N, 2, T) !== nothing
    i = findfirst(N -> any(c -> paired(c, N), cells), Ns)
    N_w = i === nothing ? 0 : Ns[i]
    ticks = i === nothing ? ["1", "2"] : ["1\n" * format_N(N_w), "2\n" * format_N(2N_w)]
    ax =
        Axis(pos; xlabel = L"\mathrm{Devices}\;D", ylabel = "Runtime [s]", xticks = ([1, 2], ticks))
    xlims!(ax, 0.85, 2.15)
    if i === nothing
        centre_note!(ax, "No D = 1, D = 2 pair at N and 2N")
        return ax
    end
    pairs = [(cell, pick(rows, cell, N_w, 1, T), pick(rows, cell, 2N_w, 2, T)) for cell in cells]
    filter!(p -> p[2] !== nothing && p[3] !== nothing, pairs)
    ideal = filter(isfinite, [v for (_, r1, _) in pairs for v in (r1.wall_s, r1.t_reg)])
    isempty(ideal) || hlines!(ax, ideal; color = :grey, linestyle = :dash, linewidth = GUIDE_WIDTH)
    for (cell, r1, r2) in pairs
        c = CELL_COLOUR[cell]
        series!(ax, [1, 2], [r1.wall_s, r2.wall_s], c; hollow = true)
        series!(ax, [1, 2], [r1.t_reg, r2.t_reg], c)
        push!(
            ids,
            "panel b: $cell D=1 N=$(r1.N_requested) $(r1.run_id), " *
            "D=2 N=$(r2.N_requested) $(r2.run_id)",
        )
    end
    ys = filter(
        isfinite,
        [v for (_, r1, r2) in pairs for v in (r1.wall_s, r2.wall_s, r1.t_reg, r2.t_reg)],
    )
    isempty(ys) || ylims!(ax, 0, 1.15 * maximum(ys))
    # Just below the highest guide at its left end: every series starts on
    # its own guide at D = 1 and rises from there, so this spot stays clear.
    isempty(ideal) || text!(
        ax,
        1.03,
        maximum(ideal);
        text = "Ideal",
        align = (:left, :top),
        offset = (4, -4),
        color = GUIDE_TEXT,
        fontsize = ANNOTATION_SIZE,
    )
    corner_note!(ax, "Weak scaling, N = $(format_N(N_w))·D, $T threads")
    return ax
end

"""
    thread_panel!(pos, rows, cells, N, ids) -> Axis

Wall time, regular-force share and irregular-force share at `N` stars on one
device against the OpenMP thread count, on a log₂ axis; the run IDs are
appended to `ids`.
"""
function thread_panel!(pos, rows, cells, N, ids)
    sel = filter(r -> r.N_requested == N && r.D == 1 && r.cell in cells, rows)
    ts = sort(unique(r.omp_threads for r in sel))
    if length(ts) < 2
        ax = Axis(pos; xlabel = "OpenMP threads", ylabel = "Runtime [s]")
        centre_note!(ax, "One thread count: no thread scaling")
        return ax
    end
    ax = Axis(
        pos;
        xlabel = "OpenMP threads",
        ylabel = "Runtime [s]",
        xscale = log2,
        xticks = (ts, string.(ts)),
    )
    ymax = 0.0
    for cell in cells
        present = [r for r in (pick(sel, cell, N, 1, t) for t in ts) if r !== nothing]
        isempty(present) && continue
        x = [r.omp_threads for r in present]
        c = CELL_COLOUR[cell]
        series!(ax, x, [r.wall_s for r in present], c; hollow = true)
        series!(ax, x, [r.t_reg for r in present], c)
        series!(ax, x, [r.t_irr for r in present], c; linestyle = :dash)
        ys = [v for r in present for v in (r.wall_s, r.t_reg, r.t_irr) if isfinite(v)]
        ymax = max(ymax, maximum(ys; init = 0.0))
        runs = join(["t=$(r.omp_threads) $(r.run_id)" for r in present], ", ")
        push!(ids, "panel c: $cell $runs")
    end
    xlims!(ax, ts[1] / 2^0.3, ts[end] * 2^0.3)
    # The fewest threads are the slowest: headroom keeps the corner note clear.
    ymax > 0 && ylims!(ax, 0, 1.3 * ymax)
    corner_note!(ax, "Host threads, N = $(format_N(N)), D = 1")
    return ax
end

"""
    scaling_figure(rows, out)

Strong scaling, weak scaling and host-thread scaling on the machine with
two-device rows, saved as `gpu_cells_scaling.{pdf,png}` in `out`, with the run
IDs behind each panel printed. Skipped with a notice without two-device rows.
"""
function scaling_figure(rows, out)
    machines = sort(unique(r.machine for r in rows if r.D == 2))
    if isempty(machines)
        println("no two-device rows: scaling figure skipped")
        return nothing
    end
    machine = first(machines)
    length(machines) > 1 &&
        println("two-device rows on $(join(machines, ", ")): scaling figure for $machine")
    mrows = filter(r -> r.machine == machine, rows)
    N_max = maximum(r.N_requested for r in mrows)
    threads_at = D -> Set(r.omp_threads for r in mrows if r.N_requested == N_max && r.D == D)
    common = intersect(threads_at(1), threads_at(2))
    if isempty(common)
        println("no thread count at both D = 1 and D = 2 for N = $N_max: scaling figure skipped")
        return nothing
    end
    T0 = minimum(common)
    cells = [c for c in CELLS if any(r -> r.cell == c, mrows)]
    ids = String[]

    fig = Figure(; size = (1500, 560))
    strong_panel!(fig[1, 1], mrows, cells, N_max, T0, ids)
    weak_panel!(fig[1, 2], mrows, cells, T0, ids)
    thread_panel!(fig[1, 3], mrows, cells, N_max, ids)
    Legend(
        fig[0, 1:3],
        [
            [series_element(CELL_COLOUR[c]) for c in cells],
            [
                series_element(:black),
                series_element(:black; hollow = true),
                series_element(:black; linestyle = :dash),
            ],
        ],
        [
            [CELL_LABEL[c] for c in cells],
            ["Regular force (device)", "Wall time", "Irregular force (host)"],
        ],
        ["Cell:", "Quantity:"];
        orientation = :horizontal,
        titleposition = :left,
        nbanks = 1,
        framevisible = false,
        tellheight = true,
        tellwidth = false,
        groupgap = 30,
    )
    colgap!(fig.layout, 18)
    rowgap!(fig.layout, 8)
    save_figure(fig, out, "gpu_cells_scaling")
    println("scaling figure: $machine, N_max = $N_max, $T0 threads")
    foreach(println, ids)
    return nothing
end

# ---------------------------------------------------------------------------
# Figure 2: the cards
# ---------------------------------------------------------------------------

"""
    card_selection(rows, card) -> NamedTuple or nothing

Single-device rows of `card` (its `match` in the `gpu` column) at the largest
N present and `PREFERRED_THREADS` threads if run there, otherwise the fewest:
`(N, threads, cells)` with `cells` mapping each cell to its last such row.
"""
function card_selection(rows, card)
    sel = filter(r -> occursin(card.match, r.gpu) && r.D == 1, rows)
    isempty(sel) && return nothing
    N = maximum(r.N_requested for r in sel)
    ts = Set(r.omp_threads for r in sel if r.N_requested == N)
    T = PREFERRED_THREADS in ts ? PREFERRED_THREADS : minimum(ts)
    cells = Dict{String,Any}()
    for cell in CELLS
        r = pick(sel, cell, N, 1, T)
        r === nothing || (cells[cell] = r)
    end
    return (N = N, threads = T, cells = cells)
end

"""
    bars!(ax, entries, value, fillto, label; tick = nothing)

Bars of `value(entry, row)` per card and cell, single filled on the left and
dual hollow on the right, from `fillto`, labelled above with `label(value)` in
the cell colour; `tick(entry, row)`, when given, draws a black tick across
each bar.
"""
function bars!(ax, entries, value, fillto, label; tick = nothing)
    for cell in CELLS
        xs, ys, marks = Float64[], Float64[], Point2f[]
        for (x, e) in enumerate(entries)
            r = get(e.cells, cell, nothing)
            r === nothing && continue
            y = value(e, r)
            isfinite(y) || continue
            xb = x + BAR_OFFSET[cell]
            push!(xs, xb)
            push!(ys, y)
            tick === nothing && continue
            t = tick(e, r)
            if isfinite(t) && t > fillto
                push!(marks, Point2f(xb - BAR_WIDTH / 2, t), Point2f(xb + BAR_WIDTH / 2, t))
            end
        end
        isempty(xs) && continue
        c = CELL_COLOUR[cell]
        style = if cell == "single"
            (color = c, strokecolor = edge(c), strokewidth = 1.5)
        else
            (color = (:white, 0.0), strokecolor = c, strokewidth = 2.5)
        end
        barplot!(ax, xs, ys; width = BAR_WIDTH, gap = 0, fillto = fillto, style...)
        isempty(marks) || linesegments!(ax, marks; color = :black, linewidth = 2.5)
        text!(
            ax,
            xs,
            ys;
            text = label.(ys),
            align = (:center, :bottom),
            offset = (0, 4),
            fontsize = VALUE_SIZE,
            color = c,
        )
    end
    return nothing
end

"""
    class_divider!(ax, entries, xmin, xmax)

Dashed divider between the last datacenter and the first consumer card, with
the class names at the top of the axis on either side; nothing unless both
classes are present.
"""
function class_divider!(ax, entries, xmin, xmax)
    n_dc = count(e -> e.card.class == "datacenter", entries)
    0 < n_dc < length(entries) || return nothing
    xd = n_dc + 0.5
    vlines!(ax, [xd]; color = :grey, linestyle = :dash, linewidth = GUIDE_WIDTH)
    xr = (xd - xmin) / (xmax - xmin)
    for (s, align, dx) in (("Datacenter", (:right, :top), -8), ("Consumer", (:left, :top), 8))
        text!(
            ax,
            xr,
            0.97;
            text = s,
            space = :relative,
            align = align,
            offset = (dx, 0),
            color = GUIDE_TEXT,
            fontsize = ANNOTATION_SIZE,
        )
    end
    return nothing
end

"""
    hardware_figure(rows, cards, out)

Engine wall time per card and cell with a tick at the regular-force share,
and, when FP32 peaks are known, the regular-force kernel rate as a percentage
of the peak; saved as `gpu_cells_hardware.{pdf,png}` in `out`, with the run
IDs and peaks printed. Cards without rows are skipped with a notice.
"""
function hardware_figure(rows, cards, out)
    entries = NamedTuple[]
    for card in cards
        s = card_selection(rows, card)
        if s === nothing
            println("card $(card.label): no single-device rows, skipped")
            continue
        end
        push!(entries, (card = card, N = s.N, threads = s.threads, cells = s.cells))
    end
    if isempty(entries)
        println("no card has rows: hardware figure skipped")
        return nothing
    end
    walls = [r.wall_s for e in entries for r in values(e.cells) if isfinite(r.wall_s)]
    filter!(>(0), walls)
    # The regular-force ticks sit below the bar tops; the axis floor must show them too.
    ticks = [r.t_reg for e in entries for r in values(e.cells) if isfinite(r.t_reg) && r.t_reg > 0]
    if isempty(walls)
        println("no finite wall time: hardware figure skipped")
        return nothing
    end
    n = length(entries)
    with_rate = any(e -> isfinite(e.card.peak), entries)
    xmin, xmax = 0.4, n + 0.6
    x_axis = (
        xticks = (1:n, [e.card.label for e in entries]),
        xticklabelrotation = n > 4 ? π / 6 : 0.0,
        xgridvisible = false,
    )

    fig = Figure(; size = (900, with_rate ? 950 : 600))
    y_lo = exp10(floor(log10(0.7 * minimum(vcat(walls, ticks)))))
    # Four times the tallest bar: room for its value label and, above that,
    # for the class labels beside the divider.
    y_hi = exp10(log10(maximum(walls)) + 0.6)
    ax_t = Axis(
        fig[1, 1];
        yscale = log10,
        ylabel = "Engine wall time [s]",
        yticks = decade_ticks(y_lo, y_hi),
        xticklabelsvisible = !with_rate,
        x_axis...,
    )
    bars!(ax_t, entries, (e, r) -> r.wall_s, y_lo, format_wall; tick = (e, r) -> r.t_reg)
    class_divider!(ax_t, entries, xmin, xmax)
    xlims!(ax_t, xmin, xmax)
    ylims!(ax_t, y_lo, y_hi)

    if with_rate
        rate = (e, r) -> 100 * r.gflops_mean / (e.card.peak * 1000)
        ax_b = Axis(fig[2, 1]; ylabel = "Kernel rate [% of FP32 peak]", x_axis...)
        bars!(ax_b, entries, rate, 0.0, format_rate)
        class_divider!(ax_b, entries, xmin, xmax)
        rates = [rate(e, r) for e in entries for r in values(e.cells)]
        filter!(isfinite, rates)
        ylims!(ax_b, 0, isempty(rates) ? 1.0 : 1.25 * maximum(rates))
        linkxaxes!(ax_t, ax_b)
        xlims!(ax_b, xmin, xmax)
    end

    c1, c2 = CELL_COLOUR["single"], CELL_COLOUR["dual"]
    Legend(
        fig[0, 1],
        Makie.LegendElement[
            PolyElement(; color = c1, strokecolor = edge(c1), strokewidth = 1.5),
            PolyElement(; color = (:white, 0.0), strokecolor = c2, strokewidth = 2.5),
            LineElement(; color = :black, linewidth = 2.5),
        ],
        ["Single cell", "Dual cell", "Tick: regular-force share"];
        orientation = :horizontal,
        framevisible = false,
        tellheight = true,
        tellwidth = false,
    )
    rowgap!(fig.layout, 12)
    save_figure(fig, out, "gpu_cells_hardware")
    for e in entries
        present = [cell for cell in CELLS if haskey(e.cells, cell)]
        ids = join(["$cell $(e.cells[cell].run_id)" for cell in present], ", ")
        peak = if isfinite(e.card.peak)
            @sprintf("%.1f TFLOP/s (%s)", e.card.peak, e.card.source)
        else
            "none"
        end
        println("$(e.card.label) (N = $(e.N), $(e.threads) threads): $ids; FP32 peak $peak")
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Console table
# ---------------------------------------------------------------------------

"""
    print_table(rows)

All rows, sorted by machine, cell, N, threads and device count, with the
wall-time shares of the regular force, irregular force and KS.
"""
function print_table(rows)
    sorted = sort(rows; by = r -> (r.machine, r.cell, r.N_requested, r.omp_threads, r.D))
    w = max(7, maximum(length(r.machine) for r in rows))
    println()
    print(rpad("machine", w), "  ")
    @printf(
        "%-6s %8s %2s %4s %10s %10s %10s %10s %8s %6s %9s  %s\n",
        "cell",
        "N_total",
        "D",
        "thr",
        "wall[s]",
        "reg[s]",
        "irr[s]",
        "ks[s]",
        "GFLOP/s",
        "util%",
        "power[W]",
        "run_id",
    )
    for r in sorted
        print(rpad(r.machine, w), "  ")
        @printf(
            "%-6s %8d %2d %4d %10.1f %10.1f %10.1f %10.1f %8.0f %6.1f %9.1f  %s\n",
            r.cell,
            r.N_total,
            r.D,
            r.omp_threads,
            r.wall_s,
            r.t_reg,
            r.t_irr,
            r.t_ks,
            r.gflops_mean,
            r.gpu_util_mean,
            r.gpu_power_mean_w,
            r.run_id,
        )
    end
    println()
    return nothing
end

# ---------------------------------------------------------------------------
# Synthetic rehearsal data set
# ---------------------------------------------------------------------------

# Five hosts of the rehearsal: card, device string, machine, host slowness of
# the irregular force and other work, and whether the node has two devices.
const SYNTHETIC_HOSTS = [
    (
        match = "Tesla T4",
        label = "Tesla T4",
        class = "datacenter",
        peak = 8.1,
        gpu = "Tesla T4, 15360 MiB, 575.57.08, 7.5",
        machine = "lisaisssci-EPYC7551P-T4",
        g_host = 2.5,
        two_devices = false,
    ),
    (
        match = "H200",
        label = "H200 NVL",
        class = "datacenter",
        peak = 60.0,
        gpu = "NVIDIA H200 NVL, 143771 MiB, 595.71.05, 9.0",
        machine = "issaf-0-6-H200x2",
        g_host = 1.2,
        two_devices = true,
    ),
    (
        match = "RTX 2080 Super",
        label = "RTX 2080 Super MQ",
        class = "consumer",
        peak = 6.0,
        gpu = "NVIDIA GeForce RTX 2080 Super with Max-Q Design, 8192 MiB, 610.57.04, 7.5",
        machine = "trustee-i7-2080SMQ",
        g_host = 2.0,
        two_devices = false,
    ),
    (
        match = "RTX 5070 Ti",
        label = "RTX 5070 Ti",
        class = "consumer",
        peak = 43.9,
        gpu = "NVIDIA GeForce RTX 5070 Ti, 16303 MiB, 610.57.04, 12.0",
        machine = "nicolin-i9-5070Ti",
        g_host = 1.1,
        two_devices = false,
    ),
    (
        match = "RTX 5090",
        label = "RTX 5090",
        class = "consumer",
        peak = 104.8,
        gpu = "NVIDIA GeForce RTX 5090, 32607 MiB, 610.57.04, 12.0",
        machine = "nicolin-9950X-5090",
        g_host = 1.0,
        two_devices = false,
    ),
]

"""One CSV cell: a string quoted when it holds a comma, a number as is."""
csv_field(x) = x isa AbstractString && occursin(',', x) ? "\"$x\"" : string(x)

"""
    synthetic_row(h, cell, N, threads, D) -> NamedTuple

One fabricated CSV row of host `h`, in column order: the regular force scales
with N^1.39, the card speed and 1/D, the irregular force and the remaining work
with the host speed, the irregular force also with threads^-0.7.
"""
function synthetic_row(h, cell, N, threads, D)
    f_card = (104.8 / h.peak)^0.5
    base = 82.0 * (N / 1e5)^1.39 * (cell == "single" ? 0.55 : 1.0)
    reg = 0.45 * base * f_card / D
    irr = 0.40 * base * h.g_host * (8 / threads)^0.7
    other = 0.15 * base * h.g_host
    wall = reg + irr + other
    total = wall * threads
    gflops = 0.25 * h.peak * 1000 * (D == 2 ? 1.8 : 1.0)
    util = 100 * reg / wall
    return (
        cell = cell,
        N_requested = N,
        N_total = cell == "dual" ? round(Int, 0.96 * N) : N,
        variant = "gpu",
        gpu_list = D == 2 ? "0 1" : "0",
        mpi_ranks = 1,
        omp_threads = threads,
        tcrit = 2.0,
        status = "completed",
        wall_s = wall,
        backend_total_s = total,
        backend_reg_s = reg * threads,
        backend_irr_s = irr * threads,
        backend_ks_s = 0.05 * total,
        backend_mdot_s = 0.03 * total,
        backend_barr_s = 0.0,
        gflops_mean = gflops,
        gflops_peak = 1.3 * gflops,
        kernel = "GPU Reg.F",
        gpu_util_mean = util,
        gpu_util_peak = min(100.0, 1.4 * util),
        gpu_power_mean_w = 0.6 * (h.peak > 50 ? 500 : 150),
        gpu_mem_peak_mib = 1800,
        peak_rss_mib = 15000 * N / 1e6,
        cpu_efficiency = 0.9,
        gpu_devices = D,
        machine = h.machine,
        gpu = h.gpu,
        run_id = "synthetic_$(h.machine)_$(cell)_N$(N)_t$(threads)_D$(D)",
    )
end

"""
    write_synthetic(dir)

Deterministic rehearsal data set in `dir`: `cards.toml` with the five cards and
one `gpu_cells_synthetic_<k>.csv` per host. Every card runs both cells at
N = 10⁶ on one device with 8 threads; the two-device node adds D = 2 at
N = 10⁶, N = 5×10⁵ on one device, and 16 and 32 threads at N = 10⁶.
"""
function write_synthetic(dir)
    mkpath(dir)
    open(joinpath(dir, "cards.toml"), "w") do io
        println(io, "# Cards of the synthetic rehearsal data set, in axis order.")
        for h in SYNTHETIC_HOSTS
            println(io)
            println(io, "[[card]]")
            println(io, "match = \"$(h.match)\"")
            println(io, "label = \"$(h.label)\"")
            println(io, "class = \"$(h.class)\"")
            println(io, "fp32_peak_tflops = $(h.peak)")
        end
    end
    for (k, h) in enumerate(SYNTHETIC_HOSTS)
        points = [(cell, 1_000_000, 8, 1) for cell in CELLS]
        if h.two_devices
            for cell in CELLS
                push!(
                    points,
                    (cell, 1_000_000, 8, 2),
                    (cell, 500_000, 8, 1),
                    (cell, 1_000_000, 16, 1),
                    (cell, 1_000_000, 32, 1),
                )
            end
        end
        open(joinpath(dir, "gpu_cells_synthetic_$k.csv"), "w") do io
            println(io, join(COLUMNS, ","))
            for (cell, N, threads, D) in points
                row = synthetic_row(h, cell, N, threads, D)
                @assert [String(k) for k in keys(row)] == COLUMNS
                println(io, join(csv_field.(values(row)), ","))
            end
        end
    end
    println("synthetic data set written to ", dir)
    return nothing
end

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

"""
    main(args)

Read the rows (writing the synthetic set first under `--synthetic`), print the
table and draw both figures into the output directory.
"""
function main(args)
    opts = parse_options(args)
    data, cards_path = opts.data, opts.cards
    if opts.synthetic
        data = joinpath(opts.out, "synthetic")
        write_synthetic(data)
        cards_path = joinpath(data, "cards.toml")
    end
    isdir(data) || error("data directory $data does not exist")
    csvs, peaks = find_inputs(data)
    isempty(csvs) && error("no gpu_cells_*.csv below $data")
    rows = read_rows(csvs)
    isempty(rows) &&
        error("no completed device rows in the $(length(csvs)) CSV file(s) below $data")
    mkpath(opts.out)
    print_table(rows)
    set_theme!(publication_theme())
    scaling_figure(rows, opts.out)
    cards = cards_path === nothing ? infer_cards(rows) : read_cards(cards_path, peaks)
    hardware_figure(rows, cards, opts.out)
    return nothing
end

main(ARGS)
