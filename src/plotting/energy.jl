# =============================================================================
# Energy and virial ratio evolution plots
# =============================================================================

"""
    plot_energy(diag::DiagnosticsData, cfg::VisualizationConfig;
                filename::AbstractString = "energy")

Two-panel figure: (top) relative energy error vs time, (bottom) virial ratio.
"""
function plot_energy(
    diag::DiagnosticsData, cfg::VisualizationConfig;
    filename::AbstractString = "energy",
)
    adj = diag.adjust
    isempty(adj) && (@warn "No ADJUST data to plot"; return nothing)

    t     = [r.time_nb  for r in adj]
    de    = [r.de_rel   for r in adj]
    qvir  = [r.qvir     for r in adj]

    fig = Figure(; size = _fig_two_panel(cfg))

    ttk = _time_ticks(first(t), last(t))

    # --- Panel 1: relative energy error (log scale) ---
    ax1 = Axis(fig[1, 1];
        ylabel = L"|\Delta \mathrm{E} \, / \, \mathrm{E}|",
        title  = L"\textbf{Energy Conservation \& Virial Equilibrium}",
        yscale = log10,
        xticklabelsvisible = false,
        yminorticksvisible = false,
        xticks = ttk,
    )

    # Skip points where ΔE is exactly zero (typically the t=0 reference frame);
    # plotting them on a log axis creates a spurious spike to the axis floor.
    nz = [abs(d) > 0 for d in de]
    lines!(ax1, t[nz], abs.(de[nz]); color = :steelblue)

    # --- Panel 2: virial ratio ---
    # Nbody6++ reports Q = T/|W| (virial equilibrium at Q = 0.5)
    # For merger systems Q can reach 10^4–10^5 — use log scale when Q > 10.
    q_max = maximum(qvir)
    use_log_q = q_max > 10.0

    ax2 = Axis(fig[2, 1];
        xlabel = L"\mathrm{t} \; \mathrm{[NB]}",
        ylabel = L"\mathrm{Q} = \mathrm{T}/|\mathrm{W}|",
        xticks = ttk,
        yscale = use_log_q ? log10 : identity,
        yminorticksvisible = !use_log_q,
    )

    lines!(ax2, t, max.(qvir, 1e-3); color = :firebrick,
           label = L"\mathrm{Q} = \mathrm{T}/|\mathrm{W}|\;\mathrm{(virial\;ratio)}")
    hlines!(ax2, [0.5]; color = :gray50, linestyle = :dash, linewidth = 1.0,
            label = L"\mathrm{Q} = 0.5\;\mathrm{(virial\;equilibrium)}")

    axislegend(ax2; position = use_log_q ? :rb : :rt, framevisible = true,
               backgroundcolor = (:white, 0.6))
    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, 12)

    save(_output_path(cfg, filename), fig; px_per_unit = cfg.dpi / 72)
    return nothing
end

"""
    plot_particle_count(diag::DiagnosticsData, cfg::VisualizationConfig;
                        filename::AbstractString = "particle_count")

Two-panel figure: (top) bound particle count N, (bottom) KS binary pairs.
Each uses a linear y-axis since N and N_pairs evolve on different scales.
"""
function plot_particle_count(
    diag::DiagnosticsData, cfg::VisualizationConfig;
    filename::AbstractString = "particle_count",
)
    adj = diag.adjust
    isempty(adj) && return nothing

    t  = [r.time_nb for r in adj]
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
        ylabel = L"\mathrm{N}\;\mathrm{(bound\;particles)}",
        title  = L"\textbf{Particle \& Binary Evolution}",
        xticklabelsvisible = false,
        xticks = ttk,
        limits = (nothing, ylims_n),
    )
    lines!(ax1, t, n; color = :steelblue)

    # --- Bottom panel: N_pairs (KS binaries) ---
    np_lo, np_hi = extrema(np)
    all_zero = np_lo == 0 && np_hi == 0
    # When no binaries ever form, integer-ticks degenerates and the panel looks
    # broken — give it a small fixed range and annotate.
    np_ylims = all_zero ? (-0.5, 1.0) : nothing
    np_ytk   = all_zero ? [0.0, 1.0] : _integer_ticks(np_lo, np_hi)
    ax2 = Axis(fig[2, 1];
        xlabel = L"\mathrm{t} \; \mathrm{[NB]}",
        ylabel = L"\mathrm{N}_\mathrm{pairs}\;\mathrm{(KS\;binaries)}",
        xticks = ttk,
        yticks = np_ytk,
        limits = (nothing, np_ylims),
    )
    lines!(ax2, t, np; color = :darkorange)
    if all_zero
        text!(ax2, 0.5, 0.55;
              text = L"\mathrm{no\;KS\;binaries\;formed}",
              space = :relative, align = (:center, :center),
              color = :gray30, fontsize = 18)
    end

    linkxaxes!(ax1, ax2)
    rowgap!(fig.layout, 12)

    save(_output_path(cfg, filename), fig; px_per_unit = cfg.dpi / 72)
    return nothing
end
