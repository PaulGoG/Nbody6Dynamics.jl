# =============================================================================
# HR diagrams (Hertzsprung–Russell) of a run's stellar population
# =============================================================================
#
# Every HR figure of a run is the same layout: a legend row, one or more HR
# panels, and a census strip giving the number of stars of each class against
# time. Axes, colours, markers and legend entries are fixed by the run as a
# whole, not by the epochs a figure happens to show, so the single-epoch
# figures, the montage and the animation of one run are directly comparable.
#
# Evolved phases of massive stars last less than the output interval and hold
# a handful of stars, so a class is typically present at some epochs and
# absent at others, and compact remnants have no place on the plane at all.
# The census strip carries that history: a class enters the strip where it
# first appears, its line breaks while it is absent, and neutron stars and
# black holes are followed there although no panel can show them.

"""Colour, marker and census line style of a stellar class, the same in every figure."""
const _HRClassStyle = @NamedTuple{color::Makie.RGBAf, marker::Symbol, linestyle::Symbol}

_hr_grey(v) = Makie.RGBAf(v, v, v, 1)

# The main sequence is the backdrop, in grey; the Okabe–Ito colours go to the
# classes that come and go. Marker shape follows the family: circles for
# hydrogen-rich evolved stars, a diamond for helium stars, triangles and a
# square for remnants. Classes that are never drawn on the plane have broken
# census lines.
const _HR_CLASS_STYLE = Dict{Symbol,_HRClassStyle}(
    :pre_main_sequence => (color = _hr_grey(0.3), marker = :rect, linestyle = :solid),
    :main_sequence => (color = _hr_grey(0.62), marker = :circle, linestyle = :solid),
    :hertzsprung_gap => (color = _OKABE_ITO[2], marker = :circle, linestyle = :solid),
    :red_giant => (color = _OKABE_ITO[6], marker = :circle, linestyle = :solid),
    :core_helium_burning => (color = _OKABE_ITO[5], marker = :circle, linestyle = :solid),
    :asymptotic_giant => (color = _OKABE_ITO[4], marker = :circle, linestyle = :solid),
    :helium_star => (color = _OKABE_ITO[3], marker = :diamond, linestyle = :solid),
    :white_dwarf => (color = _OKABE_ITO[1], marker = :utriangle, linestyle = :solid),
    :neutron_star => (color = _OKABE_ITO[8], marker = :dtriangle, linestyle = :dash),
    :black_hole => (color = _OKABE_ITO[8], marker = :rect, linestyle = :dot),
    :massless_remnant => (color = _hr_grey(0.45), marker = :xcross, linestyle = :dot),
)

_hr_style(k::Int)::_HRClassStyle = _HR_CLASS_STYLE[STELLAR_CLASSES[k].key]

"""The class drawn beneath all others and left out of the census strip."""
const _HR_BACKDROP = :main_sequence
_hr_is_backdrop(k::Int) = STELLAR_CLASSES[k].key === _HR_BACKDROP

"""Marker size of the backdrop class relative to `_STYLE.marker`."""
const _HR_BACKDROP_SCALE = 0.7
"""Height of the census strip relative to a single-panel figure."""
const _HR_STRIP_FRAC = 0.3
"""Layout units one bank of the legend row takes (the group headers take one more)."""
const _HR_LEGEND_BANK_HEIGHT = 38
"""Layout units one class entry of the legend takes along the row (marker, longest label, gap)."""
const _HR_LEGEND_ENTRY_WIDTH = 270
"""Layout units the key group of the legend takes along the row."""
const _HR_LEGEND_KEY_WIDTH = 330
const _HR_MIN_SPAN_DEX = 0.2

"""`(lo, hi)` widened by 6 % on each side, or to `_HR_MIN_SPAN_DEX` about the midpoint when narrower."""
function _padded_range(lo::Real, hi::Real)
    span = hi - lo
    if span < _HR_MIN_SPAN_DEX
        mid = 0.5 * (lo + hi)
        return (mid - 0.5 * _HR_MIN_SPAN_DEX, mid + 0.5 * _HR_MIN_SPAN_DEX)
    end
    return (lo - 0.06 * span, hi + 0.06 * span)
end

"""
    _HRRun

What every HR figure of one run shares: the HR-plane population of each
epoch, the class census, the classes present at any epoch, and the axis
limits (log T_eff, log L) spanning all epochs.
"""
struct _HRRun
    populations::Vector{HRPopulation}
    census::StellarCensus
    classes::Vector{Int}
    tlims::NTuple{2,Float64}
    llims::NTuple{2,Float64}
end

function _hr_run(
    sevs::AbstractVector{StellarEvolutionSnapshot},
    bevs::AbstractVector{BinaryEvolutionSnapshot},
)::_HRRun
    isempty(sevs) && throw(ArgumentError("no stellar-evolution snapshots to plot"))
    pops = hr_populations(sevs, bevs)
    all(p -> length(p) == 0, pops) &&
        throw(ArgumentError("no star of any snapshot lies on the HR plane"))
    teff = (extrema(p.log_teff) for p in pops if length(p) > 0)
    lum = (extrema(p.log_luminosity) for p in pops if length(p) > 0)
    tlims = _padded_range(minimum(first, teff), maximum(last, teff))
    llims = _padded_range(minimum(first, lum), maximum(last, lum))
    census = stellar_census(sevs, bevs)
    return _HRRun(pops, census, classes_present(census), tlims, llims)
end

"""Classes of the run that the census strip follows: all but the backdrop."""
_hr_strip_classes(run::_HRRun) = [k for k in run.classes if !_hr_is_backdrop(k)]

"""Whether the run has a history to show: at least two epochs and one class besides the backdrop."""
_hr_has_strip(run::_HRRun) = length(run.census.time_myr) ≥ 2 && !isempty(_hr_strip_classes(run))

"""Whether an evolved star is a member of a KS pair at any of the `epochs` of the run."""
function _hr_has_evolved_members(run::_HRRun, epochs)::Bool
    return any(view(run.populations, epochs)) do p
        any(i -> p.binary_member[i] && !_hr_is_backdrop(p.class[i]), eachindex(p.class))
    end
end

function _hr_points(pop::HRPopulation, k::Int, member::Bool)::Vector{Point2f}
    return [
        Point2f(pop.log_teff[i], pop.log_luminosity[i]) for
        i in eachindex(pop.class) if pop.class[i] == k && pop.binary_member[i] == member
    ]
end

"""
    _draw_hr_population!(ax, pop, classes; scale = 1.0)

Draw the observable population `pop` on `ax`, one scatter per class of
`classes` and multiplicity, in a fixed order: the backdrop first and small,
then each evolved class at full size with a darker edge — single stars
filled, members of KS pairs open — so that one evolved star stays visible
against thousands of main-sequence stars.
"""
function _draw_hr_population!(
    ax,
    pop::Observable{HRPopulation},
    classes::Vector{Int};
    scale::Real = 1.0,
)
    drawn = [k for k in classes if STELLAR_CLASSES[k].luminous]
    for k in vcat(filter(_hr_is_backdrop, drawn), filter(!_hr_is_backdrop, drawn))
        style = _hr_style(k)
        if _hr_is_backdrop(k)
            # Multiplicity cannot be read off a saturated band; one style.
            points = @lift(vcat(_hr_points($pop, k, false), _hr_points($pop, k, true)))
            scatter!(
                ax,
                points;
                color = style.color,
                marker = style.marker,
                markersize = _HR_BACKDROP_SCALE * _STYLE.marker * scale,
                strokewidth = 0,
            )
            continue
        end
        scatter!(
            ax,
            @lift(_hr_points($pop, k, false));
            color = style.color,
            marker = style.marker,
            markersize = _STYLE.marker * scale,
            strokecolor = _band_edge(style.color),
            strokewidth = _STYLE.marker_stroke * scale,
        )
        scatter!(
            ax,
            @lift(_hr_points($pop, k, true));
            color = (style.color, 0.0),
            marker = style.marker,
            markersize = _STYLE.marker * scale,
            strokecolor = style.color,
            strokewidth = _STYLE.band_edge * scale,
        )
    end
    return nothing
end

"""
    _draw_hr_census!(ax, run, epochs)

Census strip: stars per class against time for every class of the run but
the backdrop, on a logarithmic count axis. A class absent at an epoch has no
point there, so its line starts where the class appears and breaks while it
is gone; a class present at one isolated epoch shows as a lone marker. The
observable `epochs` are marked by dashed guides.
"""
function _draw_hr_census!(ax, run::_HRRun, epochs::Observable{Vector{Float64}})
    counts = class_counts(run.census)
    vlines!(ax, epochs; color = (:grey, 0.7), linestyle = :dash, linewidth = _STYLE.guide)
    for k in _hr_strip_classes(run)
        style = _hr_style(k)
        n = [c > 0 ? Float64(c) : NaN for c in view(counts, :, k)]
        scatterlines!(
            ax,
            run.census.time_myr,
            n;
            color = style.color,
            marker = style.marker,
            linestyle = style.linestyle,
            markersize = _STYLE.marker,
            strokecolor = _band_edge(style.color),
            strokewidth = _STYLE.marker_stroke,
        )
    end
    return nothing
end

"""Smallest 1–2–5 tick value not below `n_max`, so the largest count of the census lies under a labelled tick."""
function _census_axis_top(n_max::Integer)::Float64
    decade = 10.0^floor(log10(max(n_max, 1)))
    for m in (1, 2, 5, 10)
        m * decade ≥ n_max && return m * decade
    end
    return 10 * decade
end

"""Legend groups of an HR figure: the classes of the run, then the marker and guide conventions in use."""
function _hr_legend_groups(run::_HRRun, strip::Bool, guide_label::AbstractString, drawn)
    # Without a strip nothing in the figure stands for a class that is not
    # drawn on the plane, so such a class gets no entry.
    listed = strip ? run.classes : [k for k in run.classes if STELLAR_CLASSES[k].luminous]
    class_entries = Vector{_LegendElement}[]
    class_labels = AbstractString[]
    for k in listed
        style = _hr_style(k)
        entry = _LegendElement[]
        strip &&
            !_hr_is_backdrop(k) &&
            push!(entry, LineElement(; color = style.color, linestyle = style.linestyle))
        push!(
            entry,
            MarkerElement(;
                color = style.color,
                marker = style.marker,
                markersize = _STYLE.marker,
                strokecolor = _band_edge(style.color),
                strokewidth = _hr_is_backdrop(k) ? 0 : _STYLE.marker_stroke,
            ),
        )
        push!(class_entries, entry)
        push!(class_labels, STELLAR_CLASSES[k].label)
    end
    key_entries = Vector{_LegendElement}[]
    key_labels = AbstractString[]
    # Open markers exist only on the panels, so their entry follows the epochs
    # drawn, unlike the class entries, which the strip always answers for.
    if _hr_has_evolved_members(run, drawn)
        push!(
            key_entries,
            _LegendElement[MarkerElement(;
                color = (:black, 0.0),
                marker = :circle,
                markersize = _STYLE.marker,
                strokecolor = :black,
                strokewidth = _STYLE.band_edge,
            )],
        )
        push!(key_labels, "Member of a KS pair")
    end
    if strip
        push!(
            key_entries,
            _LegendElement[LineElement(;
                color = (:grey, 0.7),
                linestyle = :dash,
                linewidth = _STYLE.guide,
            )],
        )
        push!(key_labels, guide_label)
    end
    return class_entries, class_labels, key_entries, key_labels
end

"""
    _hr_figure(run, indices, cfg; drawn = indices) -> (fig, populations)

The HR figure of `run` showing the epochs `indices`, one panel each in a
grid of at most three columns: legend row, panels on shared axes with the
time of each in a data-free band above the data, and the census strip.
Returns the figure and the observable population of every panel; an
animation steps a single panel through the run by assigning to it, and names
in `drawn` the epochs it will show.
"""
function _hr_figure(run::_HRRun, indices::Vector{Int}, cfg::VisualizationConfig; drawn = indices)
    n = length(indices)
    ncols = min(n, 3)
    nrows = cld(n, ncols)
    grid = ncols > 1
    strip = _hr_has_strip(run)

    class_entries, class_labels, key_entries, key_labels =
        _hr_legend_groups(run, strip, n > 1 ? "Panel epoch" : "Epoch shown", drawn)
    legend = length(class_entries) + length(key_entries) ≥ 2

    pw, ph = _figsize_px(cfg)
    # The legend is as wide as the canvas allows and as deep as it then needs.
    key_width = isempty(key_entries) ? 0 : _HR_LEGEND_KEY_WIDTH
    legend_cols = max(1, floor(Int, (pw - key_width) / _HR_LEGEND_ENTRY_WIDTH))
    nbanks = max(cld(length(class_entries), legend_cols), length(key_entries), 1)
    gap = grid ? _multipanel_gap(ncols; inner_ticks = false) : _TWO_PANEL_ROWGAP
    strip_height = round(Int, _HR_STRIP_FRAC * ph)
    extra =
        (legend ? (nbanks + 1) * _HR_LEGEND_BANK_HEIGHT + gap : 0) +
        (strip ? strip_height + gap + _AXIS_PROTRUSION : 0)
    size =
        grid ? _fig_multipanel(cfg, nrows, ncols; inner_ticks = false, extra_height = extra) :
        (pw, ph + extra)
    fig = Figure(; size = size)

    target_ticks = grid ? 3 : 8
    xticks = _logval_ticks(run.tlims...; target_n = target_ticks)
    yticks = _logval_ticks(run.llims...; target_n = target_ticks)
    # Data-free band above the data, where the time annotation sits
    llims = (run.llims[1], run.llims[2] + _MONTAGE_BAND_FRAC * (run.llims[2] - run.llims[1]))
    scale = grid ? max(_multipanel_scale(cfg, ncols; inner_ticks = false), 0.75) : 1.0
    populations = [Observable(run.populations[i]) for i in indices]
    for (panel, pop) in enumerate(populations)
        row, col = cld(panel, ncols), mod1(panel, ncols)
        # The last panel of an incomplete bottom row still needs its x labels.
        bottom = row == nrows || panel + ncols > n
        ax = Axis(
            fig[row, col];
            xlabel = bottom ? L"\log_{10}(T_\mathrm{eff} \, / \, \mathrm{K})" : "",
            ylabel = col == 1 ? L"\log_{10}(L \, / \, L_\odot)" : "",
            xreversed = true,   # hot → cool from left to right
            limits = (run.tlims..., llims...),
            xticks = xticks,
            yticks = yticks,
            xticklabelsvisible = bottom,
            yticklabelsvisible = col == 1,
            xgridvisible = false,
            ygridvisible = false,
        )
        time_label = @lift(latexstring("t = ", _fmt_latex_sig3($pop.time_myr), "\\;\\mathrm{Myr}"))
        _annotate!(ax, time_label; corner = :tr)
        _draw_hr_population!(ax, pop, run.classes; scale = scale)
    end

    if strip
        times = run.census.time_myr
        margin = 0.04 * (times[end] - times[1])
        n_max = maximum(view(class_counts(run.census), :, _hr_strip_classes(run)))
        n_lims = (0.7, 1.15 * _census_axis_top(n_max))
        ax = Axis(
            fig[nrows + 1, 1:ncols];
            xlabel = L"t \; [\mathrm{Myr}]",
            ylabel = L"N",
            yscale = log10,
            limits = (times[1] - margin, times[end] + margin, n_lims...),
            xticks = _time_ticks(times[1], times[end]),
            yticks = _log_ticks(n_lims...),
            xgridvisible = false,   # the only vertical guides are the epochs
        )
        epochs = lift((ps...) -> Float64[p.time_myr for p in ps], populations...)
        _draw_hr_census!(ax, run, epochs)
        rowsize!(fig.layout, nrows + 1, Fixed(strip_height))
    end

    if legend
        groups = isempty(key_entries) ? [class_entries] : [class_entries, key_entries]
        labels = isempty(key_entries) ? [class_labels] : [class_labels, key_labels]
        titles = isempty(key_entries) ? ["Stellar class:"] : ["Stellar class:", "Key:"]
        Legend(
            fig[0, 1:ncols],
            groups,
            labels,
            titles;
            orientation = :horizontal,
            nbanks = nbanks,
            framevisible = false,
            titleposition = :top,
            titlehalign = :left,
            titlesize = _STYLE.label,
            tellheight = true,
            padding = (0, 0, 0, 0),
        )
    end

    colgap!(fig.layout, gap)
    rowgap!(fig.layout, gap)
    return fig, populations
end

"""Note for a figure without a census strip: the stars of the epoch that have no place on the plane."""
function _hr_dark_note(run::_HRRun, epoch::Int)::String
    counts = class_counts(run.census)
    parts = [
        "$(counts[epoch, k]) $(lowercase(STELLAR_CLASSES[k].label))" for
        k in run.classes if !STELLAR_CLASSES[k].luminous && counts[epoch, k] > 0
    ]
    return isempty(parts) ? "" : "Not on the plane: " * join(parts, ", ")
end

"""
    plot_hr(sevs::Vector{StellarEvolutionSnapshot}, cfg::VisualizationConfig;
            epoch = length(sevs), bevs = BinaryEvolutionSnapshot[],
            filename = "hr_diagram") -> String
    plot_hr(sev::StellarEvolutionSnapshot, cfg::VisualizationConfig;
            bev = nothing, filename = "hr_diagram") -> String

Hertzsprung–Russell diagram (log T_eff against log L) of epoch `epoch` of a
run, coloured by stellar class ([`STELLAR_CLASSES`](@ref)), above a census
strip of the stars per class against time with the epoch marked. Axes and
legend are those of the whole run, so the figures of different epochs of
one run match. With `bevs`, the members of KS-regularised pairs are part of
the population and drawn as open markers.

The single-snapshot method has no run to refer to: it shows the classes of
that snapshot alone, without a strip, and notes the neutron stars and black
holes it cannot draw.

Returns the output file path.
"""
@publication function Nbody6Dynamics.plot_hr(
    sevs::Vector{StellarEvolutionSnapshot},
    cfg::VisualizationConfig;
    epoch::Int = length(sevs),
    bevs::Vector{BinaryEvolutionSnapshot} = BinaryEvolutionSnapshot[],
    filename::AbstractString = "hr_diagram",
)::String
    run = _hr_run(sevs, bevs)
    epoch in eachindex(sevs) || throw(ArgumentError("epoch = $epoch is outside 1:$(length(sevs))"))
    fig, _ = _hr_figure(run, [epoch], cfg)
    if !_hr_has_strip(run)
        note = _hr_dark_note(run, epoch)
        isempty(note) || _annotate!(content(fig[1, 1]), note; corner = :tr, dy = 0.08)
    end
    return _save_fig(cfg, filename, fig)
end

@publication function Nbody6Dynamics.plot_hr(
    sev::StellarEvolutionSnapshot,
    cfg::VisualizationConfig;
    bev::Union{Nothing,BinaryEvolutionSnapshot} = nothing,
    filename::AbstractString = "hr_diagram",
)::String
    bevs = bev === nothing ? BinaryEvolutionSnapshot[] : [bev]
    return plot_hr([sev], cfg; bevs = bevs, filename = filename)
end

"""
    plot_hr_evolution(sevs::Vector{StellarEvolutionSnapshot}, cfg::VisualizationConfig;
                      bevs = BinaryEvolutionSnapshot[], max_panels = 6,
                      epochs = nothing, filename = "hr_evolution") -> String

HR diagrams of up to `max_panels` epochs spaced evenly through the run (or
of the snapshot indices `epochs`), on shared axes under one legend, above a
census strip of the stars per class against time with the panel epochs
marked.

The legend lists every class present at any epoch of the run, not only at
the epochs drawn: evolved phases are short against the output interval, and
a class the panels miss still shows in the strip.

Returns the output file path.
"""
@publication function Nbody6Dynamics.plot_hr_evolution(
    sevs::Vector{StellarEvolutionSnapshot},
    cfg::VisualizationConfig;
    bevs::Vector{BinaryEvolutionSnapshot} = BinaryEvolutionSnapshot[],
    max_panels::Int = 6,
    epochs::Union{Nothing,Vector{Int}} = nothing,
    filename::AbstractString = "hr_evolution",
)::String
    run = _hr_run(sevs, bevs)
    max_panels ≥ 1 || throw(ArgumentError("max_panels must be at least 1; got $max_panels"))
    n = length(sevs)
    indices = if epochs !== nothing
        all(in(eachindex(sevs)), epochs) ||
            throw(ArgumentError("epochs must lie within 1:$n; got $epochs"))
        epochs
    elseif n ≤ max_panels
        collect(1:n)
    else
        unique(round.(Int, range(1, n; length = max_panels)))
    end
    fig, _ = _hr_figure(run, indices, cfg)
    return _save_fig(cfg, filename, fig)
end

"""
    animate_hr(sevs::Vector{StellarEvolutionSnapshot}, cfg::VisualizationConfig;
               bevs = BinaryEvolutionSnapshot[], filename = "hr_evolution_anim",
               fps = nothing) -> String

Animate the HR diagram through the epochs of a run, in the layout of
[`plot_hr`](@ref): the legend and axes stay fixed while the population
changes, and the epoch guide moves along the census strip.

# Arguments
- `fps`: frames per second (`nothing` = use `cfg.style.anim_fps`; `0` there
  auto-calculates for ~`cfg.style.anim_target_seconds` s, clamped to 1–8 fps).
  HR frames are information-dense, so the auto rate favours a slower pace
  than cluster animations.

Returns the output file path.
"""
@publication function Nbody6Dynamics.animate_hr(
    sevs::Vector{StellarEvolutionSnapshot},
    cfg::VisualizationConfig;
    bevs::Vector{BinaryEvolutionSnapshot} = BinaryEvolutionSnapshot[],
    filename::AbstractString = "hr_evolution_anim",
    fps::Union{Int,Nothing} = nothing,
)::String
    run = _hr_run(sevs, bevs)
    nframes = length(sevs)
    fps = something(fps, cfg.style.anim_fps)
    fps > 0 || (
        fps = _auto_fps(
            nframes;
            target_duration = cfg.style.anim_target_seconds,
            min_fps = 1,
            max_fps = 8,
        )
    )
    outpath = _anim_output_path(cfg, filename)
    @info "Animating HR diagram: $nframes frames → $outpath  ($(fps) fps, ~$(round(Int, nframes/fps)) s)"

    fig, populations = _hr_figure(run, [1], cfg; drawn = eachindex(run.populations))
    _backup_existing(outpath)
    record(fig, outpath, 1:nframes; framerate = fps, px_per_unit = cfg.style.anim_px_per_unit) do i
        populations[1][] = run.populations[i]
    end

    @info "HR animation saved: $outpath  ($nframes frames, $(fps) fps)"
    return outpath
end
