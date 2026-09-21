# =============================================================================
# Ensemble figures: median and central 68 % / 95 % bands over the seeds of
# every grid point, coloured by the value of one grid axis
# =============================================================================

"""Fractions of the data range reserved above and below the bands; the
annotation lives in the upper strip, which the bands never reach."""
const _ENSEMBLE_HEADROOM = 0.18
const _ENSEMBLE_FOOTROOM = 0.04

"""
    plot_sweep_ensemble(sweep_dir, cfg::VisualizationConfig;
                        quantity = :lagrangian, axis = "", fixed = Dict{String,Any}(),
                        fraction = 0.5, filename = "ensemble_<quantity>") -> String

Median (line) and central 68 % (dark band) and 95 % (light band) intervals
over the seeds of each grid point for `quantity` (`:lagrangian`,
`:energy`, `:n_stars`, `:n_pairs`; see `_run_series`), coloured by the
value of `axis` (the first grid axis by default). With several axes the
other axes are held at `fixed[key]`, or at their first value when not
given; the fixed values are annotated. `fixed` may name only other grid
axes. A sweep without grid axes yields one ensemble in a single colour.
"""
@publication function Nbody6Dynamics.plot_sweep_ensemble(
    sweep_dir::AbstractString,
    cfg::VisualizationConfig;
    quantity::Symbol = :lagrangian,
    axis::AbstractString = "",
    fixed::AbstractDict = Dict{String,Any}(),
    fraction::Real = 0.5,
    filename::AbstractString = "ensemble_" * String(quantity),
)::String
    idx = read_sweep_index(sweep_dir)
    axes = String[idx["sweep"]["axes"]...]
    ensembles = sweep_ensembles(sweep_dir, quantity; fraction = fraction)
    isempty(ensembles) && error("No completed runs with $(quantity) series in sweep $sweep_dir")

    # Hold every other axis at a fixed value (given, or the first one seen)
    held = Dict{String,Any}()
    if !isempty(axes)
        axis = _sweep_axis(idx, axis)
        for k in keys(fixed)
            String(k) in axes || throw(
                ArgumentError("fixed key \"$k\" is not a grid axis; axes: $(join(axes, ", "))"),
            )
            String(k) == axis &&
                throw(ArgumentError("fixed key \"$k\" is the comparison axis itself"))
        end
        for a in axes
            a == axis && continue
            held[a] =
                haskey(fixed, a) ? fixed[a] :
                first(
                    sort(
                        unique(e.values[findfirst(==(a) ∘ first, e.values)][2] for e in ensembles),
                    ),
                )
        end
        ensembles = filter(
            e -> all(e.values[findfirst(==(a) ∘ first, e.values)][2] == v for (a, v) in held),
            ensembles,
        )
        isempty(ensembles) &&
            error("No completed runs at the held values $(held) in sweep $sweep_dir")
    end
    axis_value(e) = isempty(axes) ? nothing : e.values[findfirst(==(axis) ∘ first, e.values)][2]
    values = isempty(axes) ? Any[nothing] : sort(unique(axis_value(e) for e in ensembles))
    colors = isempty(axes) ? Dict{Any,Any}(nothing => _OKABE_ITO[1]) : _sweep_axis_colors(values)

    fig = Figure(; size = _figsize_px(cfg))
    ax = Axis(
        fig[1, 1];
        xlabel = L"t \; [\mathrm{Myr}]",
        ylabel = _series_label(quantity, fraction),
        yscale = quantity === :energy ? log10 : identity,
    )
    # Occupancy decides between the two upper corners of the reserved
    # strip; it counts the band edges, not only the medians.
    occ_x = Float64[]
    occ_y = Float64[]
    lo_all, hi_all = Inf, 0.0
    t_lo, t_hi = Inf, -Inf
    for e in ensembles
        c = colors[axis_value(e)]
        s = e.stats
        band!(ax, s.time, s.q025, s.q975; color = (c, 0.5 * _STYLE.band_alpha))
        lines!(ax, s.time, s.q025; color = _band_edge(c), linewidth = _STYLE.envelope)
        lines!(ax, s.time, s.q975; color = _band_edge(c), linewidth = _STYLE.envelope)
        band!(ax, s.time, s.q16, s.q84; color = (c, _STYLE.band_alpha))
        lines!(ax, s.time, s.q16; color = _band_edge(c), linewidth = _STYLE.envelope)
        lines!(ax, s.time, s.q84; color = _band_edge(c), linewidth = _STYLE.envelope)
        lines!(ax, s.time, s.median; color = c, linewidth = _STYLE.data)
        for y in (s.median, s.q025, s.q975)
            append!(occ_x, s.time)
            append!(occ_y, quantity === :energy ? log10.(max.(y, eps())) : y)
        end
        lo_all = min(lo_all, minimum(s.q025))
        hi_all = max(hi_all, maximum(s.q975))
        t_lo = min(t_lo, s.time[1])
        t_hi = max(t_hi, s.time[end])
    end
    t_hi > t_lo && (ax.xticks = _time_ticks(t_lo, t_hi))
    # Reserve a data-free strip above the bands for the annotation.
    if quantity === :energy && lo_all > 0
        ax.yticks = _log_ticks(lo_all, hi_all)
        span = max(log10(hi_all / lo_all), 0.5)
        ylims!(ax, lo_all / 10^(_ENSEMBLE_FOOTROOM * span), hi_all * 10^(_ENSEMBLE_HEADROOM * span))
    elseif isfinite(lo_all) && isfinite(hi_all)
        span = max(hi_all - lo_all, eps(Float64))
        ylims!(ax, lo_all - _ENSEMBLE_FOOTROOM * span, hi_all + _ENSEMBLE_HEADROOM * span)
    end

    ns = unique(e.stats.n for e in ensembles)
    n_text = if length(ns) == 1
        "$(ns[1]) seed$(ns[1] == 1 ? "" : "s") per point"
    else
        "$(minimum(ns))–$(maximum(ns)) seeds per point"
    end
    held_text =
        isempty(held) ? "" :
        "\n" * join(
            [
                "$(_axis_short(a)) = $(_format_axis_value(v))" for
                (a, v) in sort(collect(held); by = first)
            ],
            ", ",
        )
    _annotate!(
        ax,
        n_text * held_text;
        corner = _emptiest_corner(occ_x, occ_y; corners = (:tl, :tr)),
    )

    # Legend: axis values and the meaning of the bands (typed vectors)
    band_entries = _LegendElement[]
    entries = Vector{Vector{_LegendElement}}()
    labels = Vector{Vector{AbstractString}}()
    titles = String[]
    if length(values) ≥ 2
        push!(
            entries,
            _LegendElement[
                LineElement(; color = colors[v], linewidth = _STYLE.data) for v in values
            ],
        )
        push!(labels, AbstractString[_format_axis_value(v) for v in values])
        push!(titles, _axis_short(axis) * ":")
    end
    grey = _OKABE_ITO[8]
    push!(
        entries,
        _LegendElement[
            PolyElement(;
                color = (grey, _STYLE.band_alpha),
                strokecolor = _band_edge(grey),
                strokewidth = _STYLE.envelope,
            ),
            PolyElement(;
                color = (grey, 0.5 * _STYLE.band_alpha),
                strokecolor = _band_edge(grey),
                strokewidth = _STYLE.envelope,
            ),
        ],
    )
    push!(labels, AbstractString["68 %", "95 %"])
    push!(titles, "Bands:")
    Legend(
        fig[0, :],
        entries,
        labels,
        titles;
        orientation = :horizontal,
        nbanks = length(values) > 4 ? 2 : 1,
        framevisible = false,
        titleposition = :left,
        tellheight = true,
        padding = (0, 0, 0, 0),
    )
    return _save_fig(cfg, filename, fig)
end
