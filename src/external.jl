# =============================================================================
# External post-processing — standalone analysis of arbitrary Nbody6++ output
# =============================================================================

# ---------------------------------------------------------------------------
# Output scan result type
# ---------------------------------------------------------------------------

"""
    OutputScan

Result of scanning a directory for Nbody6++ output files.
Each field is either a `Vector{String}` of found file paths or a `String` path.
The `available` dict maps data categories to `true`/`false`.
"""
struct OutputScan
    dir::String
    conf3_files::Vector{String}
    hdf5_files::Vector{String}
    stdout_file::String          # "" if not found
    lagr_file::String            # "" if not found
    escapers_file::String        # "" if not found
    stellar_evo_files::Vector{String}
    available::Dict{Symbol,Bool}
end

# ---------------------------------------------------------------------------
# Directory scanner
# ---------------------------------------------------------------------------

# Standard Nbody6++ output file names to search for
const _STDOUT_CANDIDATES  = ["out1000", "out1", "stdout", "output.log"]
const _LAGR_CANDIDATES    = ["lagr.7"]
const _ESC_CANDIDATES     = ["esc.11"]

"""
    scan_output(dir::AbstractString) -> OutputScan

Scan `dir` for Nbody6++ output files and report what is available.

Detects:
- **conf.3 snapshots**: files matching `conf.3_*` or a bare `conf.3`
- **HDF5 snapshots**: files matching `*.h5part` or `*.hdf5`
- **Diagnostics/stdout**: `out1000`, `out1`, `stdout`, `output.log`
- **Lagrangian radii**: `lagr.7`
- **Escapers**: `esc.11`
- **Stellar evolution**: files matching `sev.83_*`

The returned [`OutputScan`](@ref) contains full paths and an `available` dict
summarising which data categories were found.

# Example
```julia
scan = scan_output("/data/cluster_sim/output")
scan.available          # Dict(:snapshots_conf3 => true, :lagr => true, ...)
println(scan)           # human-readable summary
```
"""
function scan_output(dir::AbstractString)::OutputScan
    isdir(dir) || error("Directory does not exist: $dir")
    dir = abspath(dir)
    entries = readdir(dir)

    # ── conf.3 snapshots ──
    conf3_prefix = "conf.3"
    conf3_files = sort(
        [joinpath(dir, f) for f in entries
         if startswith(f, conf3_prefix * "_") || f == conf3_prefix];
        by = f -> begin
            name = basename(f)
            suffix = replace(name, conf3_prefix * "_" => ""; count = 1)
            name == conf3_prefix ? -1.0 : something(tryparse(Float64, suffix), Inf)
        end,
    )

    # ── HDF5 snapshots ──
    hdf5_files = sort([joinpath(dir, f) for f in entries
                       if endswith(f, ".h5part") || endswith(f, ".hdf5")])

    # ── Diagnostics stdout ──
    stdout_file = ""
    for cand in _STDOUT_CANDIDATES
        if cand in entries
            stdout_file = joinpath(dir, cand)
            break
        end
    end

    # ── Lagrangian radii ──
    lagr_file = ""
    for cand in _LAGR_CANDIDATES
        if cand in entries
            lagr_file = joinpath(dir, cand)
            break
        end
    end

    # ── Escapers ──
    esc_file = ""
    for cand in _ESC_CANDIDATES
        if cand in entries
            esc_file = joinpath(dir, cand)
            break
        end
    end

    # ── Stellar evolution ──
    sev_prefix = "sev.83_"
    sev_files = sort([joinpath(dir, f) for f in entries if startswith(f, sev_prefix)])

    # ── Availability summary ──
    available = Dict{Symbol,Bool}(
        :snapshots_conf3 => !isempty(conf3_files),
        :snapshots_hdf5  => !isempty(hdf5_files),
        :diagnostics     => !isempty(stdout_file),
        :lagr            => !isempty(lagr_file),
        :escapers        => !isempty(esc_file),
        :stellar_evo     => !isempty(sev_files),
    )

    return OutputScan(dir, conf3_files, hdf5_files, stdout_file,
                      lagr_file, esc_file, sev_files, available)
end

# Pretty-print for REPL
function Base.show(io::IO, ::MIME"text/plain", s::OutputScan)
    println(io, "OutputScan: $(s.dir)")
    println(io, "─────────────────────────────────────────────")

    _section(io, "conf.3 snapshots", s.available[:snapshots_conf3],
             "$(length(s.conf3_files)) files")
    _section(io, "HDF5 snapshots",   s.available[:snapshots_hdf5],
             "$(length(s.hdf5_files)) files")
    _section(io, "Diagnostics",      s.available[:diagnostics],
             isempty(s.stdout_file) ? "" : basename(s.stdout_file))
    _section(io, "Lagrangian radii", s.available[:lagr],
             isempty(s.lagr_file) ? "" : basename(s.lagr_file))
    _section(io, "Escapers",         s.available[:escapers],
             isempty(s.escapers_file) ? "" : basename(s.escapers_file))
    _section(io, "Stellar evolution", s.available[:stellar_evo],
             "$(length(s.stellar_evo_files)) files")

    n_avail = count(values(s.available))
    n_total = length(s.available)
    println(io, "─────────────────────────────────────────────")
    print(io, "Available: $n_avail / $n_total data categories")

    # Determine which plots can be generated
    plots = String[]
    s.available[:snapshots_conf3] && push!(plots, "snapshot", "snapshot_evolution", "cluster_anim")
    s.available[:diagnostics]     && push!(plots, "energy", "particle_count")
    s.available[:lagr]            && push!(plots, "lagrangian_radii", "lagrangian_anim")
    s.available[:stellar_evo]     && push!(plots, "hr_diagram", "hr_evolution", "hr_anim")
    if !isempty(plots)
        println(io)
        print(io, "Plots available: ", join(plots, ", "))
    end
end

function _section(io::IO, label::AbstractString, found::Bool, detail::AbstractString)
    mark = found ? "✓" : "✗"
    if found
        println(io, "  $mark  $label  ($detail)")
    else
        println(io, "  $mark  $label  — not found")
    end
end

# ---------------------------------------------------------------------------
# Standalone post-processing
# ---------------------------------------------------------------------------

"""
    postprocess_external(dir::AbstractString;
                         output_dir::AbstractString = "",   # default: <dir>/../plots
                         format::AbstractString = "png",
                         dpi::Int = 300,
                         figsize::Tuple{Int,Int} = (10, 8),
                         generate_plots::Bool = true,
                         generate_animations::Bool = true) -> Dict{Symbol,Any}

Post-process Nbody6++ output from an arbitrary directory.

This is the **config-free** entry point: no `Nbody6Config` or `config.toml`
required.  The function scans `dir` for all recognised output files, reads
whatever is available, runs sanity checks, and optionally generates plots
and animations.

# Arguments
- `dir`: directory containing Nbody6++ output files (conf.3_*, out1000, lagr.7, etc.)
- `output_dir`: where to save plots (default: `<dir>/../plots/`)
- `format`: plot format — `"png"`, `"pdf"`, `"svg"`
- `dpi`: output resolution
- `figsize`: figure size in inches `(width, height)`
- `generate_plots`: set `false` to skip plot generation (returns data only)
- `generate_animations`: set `false` to skip GIF animations (faster)

# Returns
A `Dict{Symbol,Any}` with keys `:scan`, `:snapshots`, `:diagnostics`,
`:lagr`, `:escapers`, `:stellar_evo` (present only if corresponding data exists).

# Example
```julia
using Nbody6Setup

# Point at any directory with Nbody6++ output
results = postprocess_external("/scratch/sim42/output")

# Or data-only (no plots):
results = postprocess_external("/scratch/sim42/output"; generate_plots=false)

# Inspect what was found:
results[:scan]  # OutputScan with availability info
```
"""
function postprocess_external(
    dir::AbstractString;
    output_dir::AbstractString = "",
    format::AbstractString = "png",
    dpi::Int = 300,
    figsize::Tuple{Int,Int} = (10, 8),
    generate_plots::Bool = true,
    generate_animations::Bool = true,
)::Dict{Symbol,Any}
    # ── Scan ──
    scan = scan_output(dir)
    show(stdout, MIME("text/plain"), scan)
    println()

    results = Dict{Symbol,Any}(:scan => scan)

    # ── Read snapshots ──
    # HDF5/H5Part snapshots are detected but not readable: the fork's KZ(46)
    # writer uses a layout (Step#i groups, numbered datasets) the removed
    # reader never supported. conf.3 is the supported snapshot source.
    if scan.available[:snapshots_hdf5]
        @warn "Found $(length(scan.hdf5_files)) HDF5 snapshot file(s) — the " *
              "KZ(46) H5Part layout is not supported; using conf.3 snapshots instead."
    end

    if scan.available[:snapshots_conf3]
        @info "Reading conf.3 snapshots ($(length(scan.conf3_files)) files)..."
        snaps = read_all_conf3(scan.dir, "conf.3_*")
        if !isempty(snaps)
            results[:snapshots] = snaps
            _sanity_snapshots(results)
        end
    end

    # ── Read diagnostics ──
    if scan.available[:diagnostics]
        @info "Reading diagnostics from $(basename(scan.stdout_file))..."
        diag = read_diagnostics(scan.stdout_file)
        results[:diagnostics] = diag
        _sanity_diagnostics(diag)
    end

    # ── Read Lagrangian radii ──
    if scan.available[:lagr]
        @info "Reading Lagrangian radii from $(basename(scan.lagr_file))..."
        lagr = read_lagr(scan.lagr_file)
        results[:lagr] = lagr
        _sanity_lagr(lagr)
    end

    # ── Read escapers ──
    if scan.available[:escapers]
        @info "Reading escapers from $(basename(scan.escapers_file))..."
        esc_data = read_escapers(scan.escapers_file)
        if !isempty(esc_data)
            results[:escapers] = esc_data
            @info "  $(length(esc_data)) escaper records"
        end
    end

    # ── Read stellar evolution ──
    if scan.available[:stellar_evo]
        @info "Reading stellar evolution ($(length(scan.stellar_evo_files)) files)..."
        sevs = read_all_stellar_evolution(scan.dir, "sev.83_*")
        if !isempty(sevs)
            results[:stellar_evo] = sevs
            _sanity_stellar_evo(sevs)
        end
    end

    # ── Generate plots ──
    if generate_plots
        out = isempty(output_dir) ? joinpath(dirname(scan.dir), "plots") : output_dir
        vis = VisualizationConfig(;
            enabled = true, format = String(format),
            dpi = dpi, figsize = figsize, output_dir = out,
        )

        set_publication_theme!()
        @info "Generating plots → $out"

        # Static plots
        if haskey(results, :snapshots)
            snaps = results[:snapshots]::Vector{Snapshot}
            if !isempty(snaps)
                plot_snapshot(snaps[end], vis; filename = "snapshot_final")
                length(snaps) > 1 && plot_snapshot_evolution(snaps, vis)

                # Merger-specific: inter-cluster separation and virial ratio
                summary_path = joinpath(scan.dir, "merger_summary.txt")
                if isfile(summary_path)
                    ranges = parse_merger_summary(summary_path)
                    if !isempty(ranges) && length(snaps) ≥ 2
                        length(ranges) ≥ 2 && plot_cluster_separation(snaps, ranges, vis)
                        plot_cluster_virial(snaps, ranges, vis)
                    end
                end
            end
        end

        if haskey(results, :diagnostics)
            diag = results[:diagnostics]::DiagnosticsData
            if !isempty(diag.adjust)
                plot_energy(diag, vis)
                plot_particle_count(diag, vis)
            end
        end

        if haskey(results, :lagr)
            lagr = results[:lagr]::LagrangianData
            !isempty(lagr.time) && plot_lagrangian(lagr, vis)
        end

        if haskey(results, :stellar_evo)
            sevs = results[:stellar_evo]::Vector{StellarEvolutionSnapshot}
            if !isempty(sevs)
                mid = max(1, length(sevs) ÷ 2)
                for (idx, fname) in [(1, "hr_diagram_early"), (mid, "hr_diagram_mid"),
                                     (length(sevs), "hr_diagram_final")]
                    plot_hr(sevs[idx], vis; filename = fname)
                end
                length(sevs) > 1 && plot_hr_evolution(sevs, vis)
            end
        end

        # Animations
        if generate_animations
            if haskey(results, :snapshots) && length(results[:snapshots]) > 1
                animate_cluster(results[:snapshots], vis)
            end
            if haskey(results, :lagr) && length(results[:lagr].time) > 1
                animate_lagrangian(results[:lagr], vis)
            end
            if haskey(results, :stellar_evo) && length(results[:stellar_evo]) > 1
                animate_hr(results[:stellar_evo], vis)
            end
        end

        @info "All plots saved to: $out"
    end

    return results
end

# ---------------------------------------------------------------------------
# Sanity checks — warnings for unusual or suspicious data
# ---------------------------------------------------------------------------

function _sanity_snapshots(results::Dict{Symbol,Any})
    haskey(results, :snapshots) || return
    snaps = results[:snapshots]::Vector{Snapshot}
    n = length(snaps)
    @info "  $n snapshots loaded"

    if n > 0
        t_first = time_nb(snaps[1].header)
        t_last  = time_nb(snaps[end].header)
        n_first = nparticles(snaps[1])
        n_last  = nparticles(snaps[end])
        @info @sprintf("  Time range: %.4f → %.4f [NB]", t_first, t_last)
        @info @sprintf("  Particles:  %d → %d", n_first, n_last)

        # Check for time monotonicity
        times = [time_nb(s.header) for s in snaps]
        if !issorted(times)
            @warn "  Snapshot times are NOT monotonically increasing — possible ordering issue"
        end

        # Check for duplicate times
        if length(unique(times)) < length(times)
            @warn "  Duplicate snapshot times detected"
        end

        # Large particle loss
        if n_last < 0.5 * n_first
            @warn @sprintf("  >50%% particle loss: %d → %d (check for dissolution or escaper flood)",
                           n_first, n_last)
        end
    end
end

function _sanity_diagnostics(diag::DiagnosticsData)
    n = length(diag.adjust)
    @info "  $n ADJUST records"

    if n > 0
        # Check energy conservation
        de_vals = [abs(r.de_rel) for r in diag.adjust if r.de_rel != 0]
        if !isempty(de_vals)
            de_max = maximum(de_vals)
            de_med = sort(de_vals)[max(1, length(de_vals) ÷ 2)]
            @info @sprintf("  Energy error: median=%.2e, max=%.2e", de_med, de_max)
            if de_max > 1e-2
                @warn "  Very large energy error detected (max |ΔE/E| > 1e-2) — check simulation stability"
            elseif de_max > 1e-4
                @warn "  Elevated energy errors (max |ΔE/E| > 1e-4) — possibly strong encounters"
            end
        end

        # Check virial ratio
        qvir_vals = [r.qvir for r in diag.adjust if r.qvir > 0]
        if !isempty(qvir_vals)
            q_med = sort(qvir_vals)[max(1, length(qvir_vals) ÷ 2)]
            q_max = maximum(qvir_vals)
            @info @sprintf("  Virial ratio: median=%.3f, max=%.1f", q_med, q_max)
        end

        # Check for N=0 records (missing TIME[NB] lines)
        n_zero = count(r -> r.n == 0, diag.adjust)
        if n_zero > 0
            @warn "  $n_zero / $n ADJUST records have N=0 (TIME[NB] lines were sparse)"
        end
    end
end

function _sanity_lagr(lagr::LagrangianData)
    nt = length(lagr.time)
    nf = length(lagr.mass_fractions)
    @info "  $nt time steps, $nf mass fractions"

    if nt > 1
        @info @sprintf("  Time range: %.4f → %.4f [NB]", lagr.time[1], lagr.time[end])
        # Check for non-positive radii (invalid on log-scale plots)
        n_neg = count(lagr.radii .<= 0)
        if n_neg > 0
            @warn "  $n_neg non-positive radius values detected (masked as NaN in log-scale plots)"
        end
    end
end

function _sanity_stellar_evo(sevs::Vector{StellarEvolutionSnapshot})
    n = length(sevs)
    @info "  $n stellar evolution snapshots"

    if n > 0
        t_first = sevs[1].time_myr
        t_last  = sevs[end].time_myr
        @info @sprintf("  Time range: %.4f → %.4f [Myr]", t_first, t_last)

        # Count stellar types across all epochs
        all_types = Set{Int}()
        for sev in sevs
            for r in sev.records
                push!(all_types, Int(r.stellar_type))
            end
        end
        type_labels = [get(STELLAR_TYPE_LABELS, kt, "K*=$kt") for kt in sort(collect(all_types))]
        @info "  Stellar types present: $(join(type_labels, ", "))"
    end
end
