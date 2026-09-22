# =============================================================================
# Energy and virial ratio evolution plots
# =============================================================================

"""
    plot_energy(diag::DiagnosticsData, cfg::VisualizationConfig;
                filename::AbstractString = "energy")

Two-panel figure: (top) relative energy error vs time, (bottom) virial ratio.
The time axis is in Myr when `cfg.units == "physical"` (from the ADJUST
records), N-body units otherwise.
"""
@publication function Nbody6Dynamics.plot_energy(
    diag::DiagnosticsData,
    cfg::VisualizationConfig;
    filename::AbstractString = "energy",
)
    adj = diag.adjust
    isempty(adj) && (@warn "No ADJUST data to plot"; return nothing)

    physical = cfg.units == "physical"
    t = physical ? [r.time_myr for r in adj] : [r.time_nb for r in adj]
    de = [r.de_rel for r in adj]
    qvir = [r.qvir for r in adj]
    tlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]"

    fig = Figure(; size = _fig_two_panel(cfg))

    ttk = _time_ticks(first(t), last(t))

    # --- Panel 1: relative energy error (log scale) ---
    # Skip points where ΔE is exactly zero (typically the t=0 reference frame);
    # plotting them on a log axis creates a spurious spike to the axis floor.
    nz = [abs(d) > 0 for d in de]
    de_abs = abs.(de[nz])
    # Explicit limits so the peak keeps headroom on the log axis; the ticks
    # follow the padded ends, not the data extrema.
    de_lims = isempty(de_abs) ? nothing : (minimum(de_abs) * 0.8, maximum(de_abs) * 1.5)

    ax1 = Axis(
        fig[1, 1];
        ylabel = L"|\Delta E \, / \, E|",
        yscale = log10,
        xticklabelsvisible = false,
        xticks = ttk,
        yticks = de_lims === nothing ? Makie.automatic : _log_ticks(de_lims...),
    )

    if de_lims !== nothing
        lines!(ax1, t[nz], de_abs; color = _SEMANTIC_COLORS[:energy_error])
        de_hi = maximum(de_abs)
        ylims!(ax1, de_lims[1], de_lims[2])
        # Annotate the maximum error (2 significant digits)
        _annotate!(
            ax1,
            latexstring("\\max|\\Delta E/E| = " * _fmt_latex_sig(de_hi, 2));
            corner = :tr,
            color = _SEMANTIC_COLORS[:energy_error],
        )
    end

    # --- Panel 2: virial ratio ---
    # Nbody6++ reports Q = T/|W| (virial equilibrium at Q = 0.5)
    # For merger systems Q can reach 10^4–10^5 — use log scale for large Q.
    q_max = maximum(qvir)
    use_log_q = q_max > cfg.style.q_log_threshold

    # Floor only on the log axis (zero/tiny Q is invalid there); raw otherwise
    q_plot = use_log_q ? max.(qvir, cfg.style.q_floor) : qvir
    # Ticks over the padded ends of the panel, the virial guide always in range
    q_lo, q_hi = extrema(q_plot)
    q_lims =
        use_log_q ? (min(q_lo, 0.5) * 0.8, max(q_hi, 0.5) * 1.5) :
        _pad_limits(min(q_lo, 0.5), max(q_hi, 0.5); floor = 0)

    ax2 = Axis(
        fig[2, 1];
        xlabel = tlabel,
        ylabel = L"Q = T/|W|",
        xticks = ttk,
        yscale = use_log_q ? log10 : identity,
        yticks = use_log_q ? _log_ticks(q_lims...) : _nice_ticks(q_lims...),
        limits = (nothing, q_lims),
    )

    lines!(
        ax2,
        t,
        q_plot;
        color = _SEMANTIC_COLORS[:virial],
        label = L"Q = T/|W|\;\mathrm{(virial\;ratio)}",
    )
    hlines!(
        ax2,
        [0.5];
        color = :gray50,
        linestyle = :dash,
        linewidth = _STYLE.guide,
        label = L"Q = 0.5\;\mathrm{(virial\;equilibrium)}",
    )

    # Shared legend for the two-panel figure: horizontal, above the axes
    _top_legend!(fig, ax2)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, _TWO_PANEL_ROWGAP)

    _save_fig(cfg, filename, fig)
    return nothing
end

"""
    plot_particle_count(diag::DiagnosticsData, cfg::VisualizationConfig;
                        filename::AbstractString = "particle_count")

Two-panel figure: (top) bound particle count N, (bottom) KS binary pairs.
Each uses a linear y-axis since N and N_pairs evolve on different scales.
A run in which no pair is ever regularised gets the single N panel, with a
note in place of the second one. The time axis is in Myr when
`cfg.units == "physical"`.
"""
@publication function Nbody6Dynamics.plot_particle_count(
    diag::DiagnosticsData,
    cfg::VisualizationConfig;
    filename::AbstractString = "particle_count",
)
    adj = diag.adjust
    isempty(adj) && return nothing

    physical = cfg.units == "physical"
    t = physical ? [r.time_myr for r in adj] : [r.time_nb for r in adj]
    n = [r.n for r in adj]
    np = [r.npairs for r in adj]

    # A run that never regularises a pair has no second panel to draw: an empty
    # one reads as a broken figure.
    no_pairs = all(iszero, np)
    tlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]"

    fig = Figure(; size = no_pairs ? _figsize_px(cfg) : _fig_two_panel(cfg))
    ttk = _time_ticks(first(t), last(t))

    # --- Top panel: N (bound particles) ---
    # Auto-scale y to the data range with ≥5 % padding; fall back to (0, n_max)
    # only when N is genuinely pinned to zero (empty simulation).
    n_lo, n_hi = extrema(n)
    if n_hi == n_lo
        y_pad = max(1.0, 0.02 * n_hi)
        ylims_n = (n_lo - y_pad, n_hi + y_pad)
    else
        span = n_hi - n_lo
        # For very stable N (typical single-cluster), inflate the window so the
        # small variations don't look like wild oscillations.
        min_span = max(0.02 * n_hi, 1.0)
        effective_span = max(span, min_span)
        pad = 0.1 * effective_span
        ylims_n = (n_lo - pad, n_hi + pad)
    end

    ax1 = Axis(
        fig[1, 1];
        xlabel = no_pairs ? tlabel : "",
        ylabel = L"N\;\mathrm{(bound\;particles)}",
        xticklabelsvisible = no_pairs,
        xticks = ttk,
        limits = (nothing, ylims_n),
    )
    lines!(ax1, t, n; color = _SEMANTIC_COLORS[:n_particles])

    # Annotate the total particle change over the run, coloured to the series
    n0, n1 = n[1], n[end]
    pct_str = @sprintf("%+.1f", 100 * (n1 - n0) / max(n0, 1))
    _annotate!(
        ax1,
        latexstring("N: $(n0) \\rightarrow $(n1)\\;($(pct_str)\\,\\%)");
        corner = :tr,
        color = _SEMANTIC_COLORS[:n_particles],
    )

    if no_pairs
        _annotate!(
            ax1,
            L"\mathrm{No\;regularised\;pairs}";
            corner = _emptiest_corner(Float64.(t), Float64.(n); corners = (:tl, :br)),
            color = :gray30,
        )
    else
        # --- Bottom panel: N_pairs (KS binaries) ---
        np_lo, np_hi = _pad_limits(minimum(np), maximum(np); floor = 0)
        ax2 = Axis(
            fig[2, 1];
            xlabel = tlabel,
            ylabel = L"N_\mathrm{pairs}\;\mathrm{(KS\;binaries)}",
            xticks = ttk,
            yticks = _integer_ticks(np_lo, np_hi),
            limits = (nothing, (np_lo, np_hi)),
        )
        lines!(ax2, t, np; color = _SEMANTIC_COLORS[:n_pairs])

        linkxaxes!(ax1, ax2)
        rowgap!(fig.layout, _TWO_PANEL_ROWGAP)
    end

    _save_fig(cfg, filename, fig)
    return nothing
end
