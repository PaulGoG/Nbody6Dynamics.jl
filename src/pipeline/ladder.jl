# =============================================================================
# N ladder: one configuration per particle number at fixed mass, size and orbit
# =============================================================================
#
# The rungs of a ladder differ in the number of bodies alone. Equal-mass
# bodies are given by the total mass of their cluster (flat `imf = "equal"`
# with `mass_total`, structured `imf = { type = "equal", target_mass }`), so
# the loader sets the body mass to M/N for every rung while the cluster
# masses, half-mass radii and the orbit, all in physical units, are those of
# the base file. The N-body units (M, RBAR, T*) are then the same along the
# ladder up to the sampling scatter of the realised half-mass radius, and an
# observable that coincides between rungs in t/t_cr is collisionless, one
# that coincides in t/t_rh collisional.

"""Leading dotted keys an arm of a ladder sweep may not address: the ladder's axis and its invariants."""
const _N_LADDER_ARM_REFUSED =
    ("merger.cluster", "merger.n_clusters", "merger.orbit_mode", "merger.stellar")

"""`[sweep]` keys the ladder writes itself; the `sweep` argument of [`n_ladder`](@ref) may not carry them."""
const _N_LADDER_SWEEP_OWN = ("name", "merger_config", "seeds", "grid")

"""
    _n_ladder_check_base(cfg::MergerConfig, raw::AbstractDict)

Refuse a base configuration whose rungs would differ in more than `N`: the
complex mode (its members' `N` are drawn, not set per cluster table), a
cluster whose IMF is not equal-mass or fixes the body mass instead of the
cluster mass, stellar evolution (`kz19 ≠ 0`), primordial binaries and a
seeded IMBH. Explicit `rs0`, `rmin` and `dtmin` are kept on every rung and
draw a warning, since the derived values scale with `N`.
"""
function _n_ladder_check_base(cfg::MergerConfig, raw::AbstractDict)
    cfg.orbit_mode == "complex" && _config_error(
        "merger.orbit_mode",
        "= \"complex\" is not an N ladder base: the complex draws its members' N; use a single, " *
        "kepler or explicit configuration",
    )
    m = raw["merger"]
    for (i, spec) in enumerate(cfg.clusters)
        key = "merger.cluster$i"
        spec.imf isa EqualMassIMF || _config_error(
            "$key.imf",
            "must be \"equal\" for an N ladder (equal-mass bodies, M/N each); got " *
            "\"$(imf_name(spec.imf))\", whose mass spectrum changes with the sample size",
        )
        imf_raw = get(m["cluster$i"], "imf", nothing)
        (imf_raw isa AbstractDict && !haskey(imf_raw, "target_mass")) && _config_error(
            "$key.imf.particle_mass",
            "fixes the body mass, so the cluster mass would grow with N along the ladder; give the " *
            "cluster mass as $key.imf.target_mass [M☉] instead",
        )
        (spec.binaries.fraction > 0 || spec.binaries.massive.model != "none") && _config_error(
            "$key.binaries",
            "must be absent for an N ladder: the hardness of primordial pairs, E_b/(m σ²) with " *
            "m = M/N, changes along the ladder at fixed physical periods, so the rungs would differ " *
            "in more than N",
        )
        spec.imbh === nothing || _config_error(
            "$key.imbh",
            "must be absent for an N ladder: M•/m grows ∝ N at fixed M•, so the IMBH is not the same " *
            "body on every rung",
        )
    end
    cfg.stellar.kz19 == 0 || _config_error(
        "merger.stellar.kz19",
        "must be 0 for an N ladder; got $(cfg.stellar.kz19): stellar evolution acts on the body " *
        "mass M/N, which changes along the ladder, so the rungs would differ in their physics",
    )
    explicit = [String(k) for k in (:rs0, :rmin, :dtmin) if getfield(cfg.nbody6, k) > 0]
    isempty(explicit) ||
        @warn "n_ladder: merger.nbody6 sets $(join(explicit, ", ")) explicitly; the value is kept " *
              "on every rung, while the derived one scales with N"
    return nothing
end

"""
    _n_ladder_split(n_total, n_base) -> Vector{Int}

The cluster sizes of a rung of `n_total` bodies in the ratio of the base
sizes `n_base`; an `ArgumentError` unless the split is exact, so the body
mass of every cluster stays `M_i / N_i` with `N_i / N_j` fixed.
"""
function _n_ladder_split(n_total::Integer, n_base::AbstractVector{Int})::Vector{Int}
    total = sum(n_base)
    unit = total ÷ gcd(n_base)
    n_total % unit == 0 || throw(
        ArgumentError(
            "n_ladder: n_total = $n_total does not split in the base ratio " *
            "$(join(n_base .÷ gcd(n_base), ':')) of the cluster sizes $(n_base); give a " *
            "multiple of $unit",
        ),
    )
    return Int[n_total * n ÷ total for n in n_base]
end

"""Write `table` as TOML to `path` (sorted keys) after the comment lines `header`, through a temporary file and a rename."""
function _write_toml_with_header(
    path::AbstractString,
    header::AbstractVector{<:AbstractString},
    table,
)
    tmp = path * ".tmp"
    open(tmp, "w") do io
        for line in header
            println(io, "# ", line)
        end
        println(io)
        TOML.print(io, table; sorted = true)
    end
    mv(tmp, path; force = true)
    return String(path)
end

"""
    n_ladder(base_merger, output_dir; n_total, seeds, sweep, name = "") -> Vector{NamedTuple}

Write an N ladder of the merger TOML `base_merger`: per rung `k` a merger
TOML `<name>_n<N>.toml`, the base with the cluster sizes replaced and
`merger.seed = seeds[k][1]`, and a sweep TOML `<name>_n<N>_sweep.toml`, the
seed ensemble `seeds[k]` of that merger TOML, both in `output_dir`. The rungs
differ in `N` only: the cluster masses, half-mass radii, profiles and the
orbit are the base's, in physical units.

# Arguments
- `base_merger`: merger TOML whose clusters have equal-mass bodies given by
  the cluster mass (flat `imf = "equal"` with `mass_total`, or a structured
  `imf` table with `type = "equal"` and `target_mass`). Refused: a body mass
  given as `particle_mass`, any other IMF, stellar evolution
  (`merger.stellar.kz19 ≠ 0`), primordial binaries, a seeded IMBH and the
  complex mode, each with the reason. Explicit `rs0`, `rmin`, `dtmin` warn.
- `output_dir`: directory of the written files (created; an existing file is
  moved to a `#k` backup, never overwritten).
- `n_total`: total body count of every rung, distinct; with several clusters
  it is split in the ratio of the base's cluster sizes and must split exactly
  (a multiple of `sum(N_i) / gcd(N_i)`).
- `seeds`: one nonempty list of distinct seeds per rung (the ensembles may
  differ in size).
- `sweep`: the further `[sweep]` keys of every rung's sweep TOML, written as
  given, so relative paths resolve against `output_dir`; `pipeline_config`
  is required, `name`, `merger_config`, `seeds` and `grid` are the ladder's
  own and refused. Arms (`"arms" => Dict(name => Dict(key => value))`) may
  not address `merger.cluster…`, `merger.n_clusters`, `merger.orbit_mode` or
  `merger.stellar…`.
- `name`: file and sweep name prefix; default the base file name without
  `.toml` and a trailing `_merger`.

Every rung's merger TOML is loaded back through [`load_merger_config`](@ref)
before anything is written (a rung above `merger.virial_max_n` fails there),
and every written sweep TOML through [`load_sweep_config`](@ref).

# Returns
One named tuple per rung, in the order of `n_total`: `n_total`, `n` (the
cluster sizes), `seeds`, `merger_config` and `sweep_config` (absolute paths).

# Example
```julia
n_ladder(
    "v3_b_king.toml",
    "v3_b";
    n_total = [25_000, 50_000, 100_000],
    seeds = [[3211, 3212], [3221], [3231]],
    sweep = Dict("pipeline_config" => "v3_b_pipeline.toml", "runs_dir" => "runs",
                 "devices" => [0, 1], "engines_per_device" => 4),
)
```
"""
function n_ladder(
    base_merger::AbstractString,
    output_dir::AbstractString;
    n_total::AbstractVector{<:Integer},
    seeds::AbstractVector{<:AbstractVector{<:Integer}},
    sweep::AbstractDict,
    name::AbstractString = "",
)
    isfile(base_merger) ||
        throw(ArgumentError("n_ladder: base configuration not found: $base_merger"))
    stem =
        isempty(name) ? replace(splitext(basename(base_merger))[1], r"_merger$" => "") :
        String(name)
    occursin(_SWEEP_NAME_PATTERN, stem) || throw(
        ArgumentError("n_ladder: name must match $(_SWEEP_NAME_PATTERN.pattern); got \"$stem\""),
    )
    isempty(n_total) && throw(ArgumentError("n_ladder: n_total is empty"))
    allunique(n_total) || throw(ArgumentError("n_ladder: n_total must be distinct; got $n_total"))
    length(seeds) == length(n_total) || throw(
        ArgumentError(
            "n_ladder: seeds holds $(length(seeds)) lists for $(length(n_total)) rungs; give one per rung",
        ),
    )
    for (k, s) in enumerate(seeds)
        (!isempty(s) && allunique(s)) || throw(
            ArgumentError(
                "n_ladder: seeds[$k] must be a nonempty list of distinct integers; got $s",
            ),
        )
    end

    # The [sweep] keys: the ladder's own are refused, the pipeline config must exist,
    # and no arm may undo the ladder's invariants
    for k in _N_LADDER_SWEEP_OWN
        haskey(sweep, k) && throw(
            ArgumentError("n_ladder: sweep.$k is set by the ladder for every rung; remove it"),
        )
    end
    haskey(sweep, "pipeline_config") ||
        throw(ArgumentError("n_ladder: sweep.pipeline_config is required"))
    pipeline =
        _resolve_path(output_dir, _typed(sweep["pipeline_config"], String, "sweep.pipeline_config"))
    isfile(pipeline) || throw(
        ArgumentError(
            "n_ladder: sweep.pipeline_config not found: $pipeline (relative paths resolve against output_dir)",
        ),
    )
    arms = get(sweep, "arms", Dict{String,Any}())
    arms isa AbstractDict ||
        throw(ArgumentError("n_ladder: sweep.arms must be a table of arm tables"))
    for (arm, table) in arms, key in (table isa AbstractDict ? keys(table) : ())
        any(p -> startswith(String(key), p), _N_LADDER_ARM_REFUSED) && throw(
            ArgumentError(
                "n_ladder: sweep.arms.$arm.\"$key\" addresses a key the ladder holds fixed or varies " *
                "(cluster tables, n_clusters, orbit_mode, stellar)",
            ),
        )
    end

    cfg = load_merger_config(base_merger)
    raw = TOML.parsefile(base_merger)
    _n_ladder_check_base(cfg, raw)
    n_base = Int[c.N for c in cfg.clusters]
    splits = [_n_ladder_split(n, n_base) for n in n_total]

    rungs = map(eachindex(n_total)) do k
        d = deepcopy(raw)
        m = d["merger"]
        m["seed"] = Int(first(seeds[k]))
        for (i, n) in enumerate(splits[k])
            m["cluster$i"]["N"] = n
        end
        d
    end
    # Every rung through the loader before a file is written
    mktempdir() do dir
        for (k, d) in enumerate(rungs)
            path = joinpath(dir, "rung$k.toml")
            open(io -> TOML.print(io, d; sorted = true), path, "w")
            try
                load_merger_config(path)
            catch err
                err isa ArgumentError || rethrow()
                throw(ArgumentError("n_ladder: rung n_total = $(n_total[k]): $(err.msg)"))
            end
        end
    end

    mkpath(output_dir)
    K = length(n_total)
    origin = "Written by Nbody6Dynamics.n_ladder from $(basename(base_merger)); regenerate instead of editing."
    return map(eachindex(n_total)) do k
        N = Int(n_total[k])
        sizes = join(["cluster$i = $n" for (i, n) in enumerate(splits[k])], ", ")
        merger_file = "$(stem)_n$(N).toml"
        merger_path = abspath(joinpath(output_dir, merger_file))
        _backup_existing(merger_path)
        _write_toml_with_header(
            merger_path,
            [
                "N ladder $stem, rung $k of $K: N = $N ($sizes).",
                "The sweep $(stem)_n$(N)_sweep.toml replaces the seed per point.",
                origin,
            ],
            rungs[k],
        )
        table = Dict{String,Any}(String(key) => v for (key, v) in sweep)
        table["name"] = "$(stem)_n$(N)"
        table["merger_config"] = merger_file
        table["seeds"] = Int[s for s in seeds[k]]
        sweep_path = abspath(joinpath(output_dir, "$(stem)_n$(N)_sweep.toml"))
        _backup_existing(sweep_path)
        _write_toml_with_header(
            sweep_path,
            ["N ladder $stem, rung $k of $K: $(length(seeds[k])) seeds of $merger_file.", origin],
            Dict{String,Any}("sweep" => table),
        )
        load_sweep_config(sweep_path)
        (
            n_total = N,
            n = splits[k],
            seeds = Int[s for s in seeds[k]],
            merger_config = merger_path,
            sweep_config = sweep_path,
        )
    end
end
