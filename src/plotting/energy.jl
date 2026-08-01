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
function plot_energy(
    diag::DiagnosticsData, cfg::VisualizationConfig;
    filename::AbstractString = "energy",
)
    adj = diag.adjust
    isempty(adj) && (@warn "No ADJUST data to plot"; return nothing)

    physical = cfg.units == "physical"
    t     = physical ? [r.time_myr for r in adj] : [r.time_nb for r in adj]
    de    = [r.de_rel   for r in adj]
    qvir  = [r.qvir     for r in adj]
    tlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]"

    fig = Figure(; size = _fig_two_panel(cfg))

    ttk = _time_ticks(first(t), last(t))

    # --- Panel 1: relative energy error (log scale) ---
    # Skip points where ΔE is exactly zero (typically the t=0 reference frame);
    # plotting them on a log axis creates a spurious spike to the axis floor.
    nz = [abs(d) > 0 for d in de]
    de_abs = abs.(de[nz])

    ax1 = Axis(fig[1, 1];
        ylabel = L"|\Delta E \, / \, E|",
        yscale = log10,
        xticklabelsvisible = false,
        xticks = ttk,
        yticks = isempty(de_abs) ? Makie.automatic : _log_ticks(extrema(de_abs)...),
    )

    if !isempty(de_abs)
        lines!(ax1, t[nz], de_abs; color = _SEMANTIC_COLORS[:energy_error])
        # Explicit limits so the peak keeps headroom on the log axis
        de_lo, de_hi = extrema(de_abs)
        ylims!(ax1, de_lo * 0.8, de_hi * 1.5)
        # Annotate the maximum error (2 significant digits)
        m_str, e_str = split(@sprintf("%.1e", de_hi), 'e')
        text!(ax1, 0.96, 0.96;
            text = latexstring("\\max|\\Delta E/E| = $(m_str) \\times 10^{$(parse(Int, e_str))}"),
            space = :relative, align = (:right, :top), fontsize = 16,
            color = _SEMANTIC_COLORS[:energy_error])
    end

    # --- Panel 2: virial ratio ---
    # Nbody6++ reports Q = T/|W| (virial equilibrium at Q = 0.5)
    # For merger systems Q can reach 10^4–10^5 — use log scale for large Q.
    q_max = maximum(qvir)
    use_log_q = q_max > cfg.style.q_log_threshold

    # Floor only on the log axis (zero/tiny Q is invalid there); raw otherwise
    q_plot = use_log_q ? max.(qvir, cfg.style.q_floor) : qvir

    ax2 = Axis(fig[2, 1];
        xlabel = tlabel,
        ylabel = L"Q = T/|W|",
        xticks = ttk,
        yscale = use_log_q ? log10 : identity,
        yticks = use_log_q ? _log_ticks(extrema(q_plot)...) : Makie.automatic,
    )

    lines!(ax2, t, q_plot; color = _SEMANTIC_COLORS[:virial],
           label = L"Q = T/|W|\;\mathrm{(virial\;ratio)}")
    hlines!(ax2, [0.5]; color = :gray50, linestyle = :dash, linewidth = 1.0,
            label = L"Q = 0.5\;\mathrm{(virial\;equilibrium)}")

    # Shared legend for the two-panel figure: horizontal, above the axes
    _top_legend!(fig, ax2)
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, 12)

    _save_fig(cfg, filename, fig)
    return nothing
end

"""
    plot_particle_count(diag::DiagnosticsData, cfg::VisualizationConfig;
                        filename::AbstractString = "particle_count")

Two-panel figure: (top) bound particle count N, (bottom) KS binary pairs.
Each uses a linear y-axis since N and N_pairs evolve on different scales.
The time axis is in Myr when `cfg.units == "physical"`.
"""
function plot_particle_count(
    diag::DiagnosticsData, cfg::VisualizationConfig;
    filename::AbstractString = "particle_count",
)
    adj = diag.adjust
    isempty(adj) && return nothing

    physical = cfg.units == "physical"
    t  = physical ? [r.time_myr for r in adj] : [r.time_nb for r in adj]
    n  = [r.n       for r in adj]
    np = [r.npairs  for r in adj]

    fig = Figure(; size = _fig_two_panel(cfg))
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

    ax1 = Axis(fig[1, 1];
        ylabel = L"N\;\mathrm{(bound\;particles)}",
        xticklabelsvisible = false,
        xticks = ttk,
        limits = (nothing, ylims_n),
    )
    lines!(ax1, t, n; color = _SEMANTIC_COLORS[:n_particles])

    # Annotate the total particle change over the run, coloured to the series
    n0, n1 = n[1], n[end]
    pct_str = @sprintf("%+.1f", 100 * (n1 - n0) / max(n0, 1))
    text!(ax1, 0.96, 0.96;
        text = latexstring("N: $(n0) \\rightarrow $(n1)\\;($(pct_str)\\%)"),
        space = :relative, align = (:right, :top), fontsize = 16,
        color = _SEMANTIC_COLORS[:n_particles])

    # --- Bottom panel: N_pairs (KS binaries) ---
    np_lo, np_hi = extrema(np)
    all_zero = np_lo == 0 && np_hi == 0
    # When no binaries ever form, integer-ticks degenerates and the panel looks
    # broken — give it a small fixed range and annotate.
    np_ylims = all_zero ? (-0.5, 1.0) : nothing
    np_ytk   = all_zero ? [0.0, 1.0] : _integer_ticks(np_lo, np_hi)
    ax2 = Axis(fig[2, 1];
        xlabel = physical ? L"t \; [\mathrm{Myr}]" : L"t \; [\mathrm{NB}]",
        ylabel = L"N_\mathrm{pairs}\;\mathrm{(KS\;binaries)}",
        xticks = ttk,
        yticks = np_ytk,
        limits = (nothing, np_ylims),
    )
    lines!(ax2, t, np; color = _SEMANTIC_COLORS[:n_pairs])
    if all_zero
        text!(ax2, 0.5, 0.55;
              text = L"\mathrm{no\;KS\;binaries\;formed}",
              space = :relative, align = (:center, :center),
              color = :gray30, fontsize = 18)
    end

    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, 12)

    _save_fig(cfg, filename, fig)
    return nothing
end
