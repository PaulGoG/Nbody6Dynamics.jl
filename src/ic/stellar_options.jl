# Stellar-evolution option set of the engine: the SSE, BSE and collision parameters that
# READSE (input.F) sets per Level and then reads from the &INSSE, &INBSE and &INCOLL groups,
# their resolution from `[merger.stellar]`, and the comparison with the engine's echo of them.

"""Members of the three namelist groups of `READSE`, in the order the generator writes them: the engine's declaration order, except `kicktype` before `bhspinfl`, as earlier versions wrote them (generated inputs stay byte-identical)."""
const _SE_GROUPS = (
    INSSE = (
        "ecflag",
        "wdflag",
        "nsflag",
        "psflag",
        "mdflag",
        "bhflag",
        "kmech",
        "mch",
        "mxns0",
        "mxns1",
        "nwind",
        "bwind",
        "flbv",
        "disp",
        "ecsig",
        "wdsig1",
        "wdsig2",
        "wdkmax",
        "vvfac",
    ),
    INBSE = (
        "lambd1",
        "alphac",
        "kicktype",
        "bhspinfl",
        "xk2",
        "xk3",
        "acc1",
        "acc2",
        "xbeta",
        "xxi",
        "epsnov",
        "eddfac",
        "gamm1",
    ),
    INCOLL = ("fctorcl",),
)

"""Order of the values in the engine's echo (`input.F` formats 30, 40, 50)."""
const _SE_ECHO_ORDER = (
    "ecflag",
    "wdflag",
    "nsflag",
    "psflag",
    "mdflag",
    "bhflag",
    "kmech",
    "mch",
    "mxns0",
    "mxns1",
    "nwind",
    "bwind",
    "flbv",
    "disp",
    "ecsig",
    "wdsig1",
    "wdsig2",
    "wdkmax",
    "vvfac",
    "lambd1",
    "alphac",
    "bhspinfl",
    "kicktype",
    "xk2",
    "xk3",
    "acc1",
    "acc2",
    "xbeta",
    "xxi",
    "epsnov",
    "eddfac",
    "gamm1",
    "fctorcl",
)

const _SEValue = Union{Int,Float64}

"""Parameters `READSE` sets for every Level other than `'0'`."""
const _SE_LEVEL_INDEPENDENT = Dict{String,_SEValue}(
    "mch" => 1.44,
    "mxns0" => 1.8,
    "mxns1" => 2.5,
    "xk2" => 0.21,
    "xk3" => 0.21,
    "acc1" => 3.920659e8,
    "acc2" => 1.5,
    "fctorcl" => 0.5,
)

"""
The parameters `READSE` (`input.F`) assigns per Level before it reads the three groups,
Kamlah et al. (2022) levels A, B and C; the Level-independent ones are in
`_SE_LEVEL_INDEPENDENT`. `mch` is a single-precision literal in the engine (1.4399999…);
the comparison with the echo works at the echo's three significant digits.
"""
const _SE_LEVEL_TABLE = Dict{String,Dict{String,_SEValue}}(
    "A" => Dict{String,_SEValue}(
        "ecflag" => 0,
        "wdflag" => 1,
        "nsflag" => 1,
        "psflag" => 0,
        "mdflag" => 1,
        "bhflag" => 0,
        "nwind" => 0.5,
        "bwind" => 0.0,
        "flbv" => 1.5,
        "kmech" => 1,
        "disp" => 190.0,
        "ecsig" => 20.0,
        "wdsig1" => 2.0,
        "wdsig2" => 2.0,
        "wdkmax" => 6.0,
        "vvfac" => 0.0,
        "lambd1" => 0.5,
        "alphac" => 3.0,
        "bhspinfl" => 0,
        "kicktype" => -2,
        "xbeta" => 0.125,
        "xxi" => 1.0,
        "epsnov" => 1.0e-3,
        "eddfac" => 100.0,
        "gamm1" => -1.0,
    ),
    "B" => Dict{String,_SEValue}(
        "ecflag" => 1,
        "wdflag" => 1,
        "nsflag" => 3,
        "psflag" => 0,
        "mdflag" => 3,
        "bhflag" => 2,
        "nwind" => 0.5,
        "bwind" => 0.0,
        "flbv" => 1.5,
        "kmech" => 1,
        "disp" => 265.0,
        "ecsig" => 3.0,
        "wdsig1" => 2.0,
        "wdsig2" => 2.0,
        "wdkmax" => 6.0,
        "vvfac" => 0.0,
        "lambd1" => 0.5,
        "alphac" => 3.0,
        "bhspinfl" => 0,
        "kicktype" => -2,
        "xbeta" => 0.125,
        "xxi" => 1.0,
        "epsnov" => 1.0e-3,
        "eddfac" => 100.0,
        "gamm1" => -1.0,
    ),
    "C" => Dict{String,_SEValue}(
        "ecflag" => 1,
        "wdflag" => 1,
        "nsflag" => 4,
        "psflag" => 1,
        "mdflag" => 4,
        "bhflag" => 2,
        "nwind" => 0.477,
        "bwind" => 0.0,
        "flbv" => 1.5,
        "kmech" => 4,
        "disp" => 265.0,
        "ecsig" => 3.0,
        "wdsig1" => 2.0,
        "wdsig2" => 2.0,
        "wdkmax" => 6.0,
        "vvfac" => 0.0,
        "lambd1" => 0.0,
        "alphac" => 1.0,
        "bhspinfl" => 0,
        "kicktype" => -2,
        "xbeta" => 0.125,
        "xxi" => 1.0,
        "epsnov" => 1.0e-3,
        "eddfac" => 100.0,
        "gamm1" => -1.0,
    ),
)

"""
The parameters `option_set = "baseline"` writes on top of Level C: the Level C values of
`ecflag`, `nsflag`, `mdflag`, `bhflag`, `disp`, `ecsig`, `lambd1`, `alphac` and `fctorcl`,
written explicitly, and the departures from Level C adopted for the campaign: the moderate
pulsational pair-instability fit of Leung et al. (2019) (`psflag = 2`), Eddington-limited
accretion onto compact objects (`eddfac = 1`), and the momentum-conserving kick (`kmech = 1`,
which Level C's `kmech = 4` already reduces to in `kick.F`).
"""
const _SE_BASELINE_SET = (
    "ecflag" => 1,
    "nsflag" => 4,
    "psflag" => 2,
    "mdflag" => 4,
    "bhflag" => 2,
    "kmech" => 1,
    "disp" => 265.0,
    "ecsig" => 3.0,
    "lambd1" => 0.0,
    "alphac" => 1.0,
    "eddfac" => 1.0,
    "fctorcl" => 0.5,
)

"""Values of `[merger.stellar] option_set`."""
const _SE_OPTION_SETS = ("level", "baseline")

"""`psflag` per `pair_instability` value (`hrdiag.f`: 1 Belczynski et al. 2016, 2 and 3 the moderate and weak fits of Leung et al. 2019)."""
const _PAIR_INSTABILITY_PSFLAG = Dict("none" => 0, "strong" => 1, "moderate" => 2, "weak" => 3)
"""`nsflag` per `remnant_mass` value (`hrdiag.f`)."""
const _REMNANT_MASS_NSFLAG = Dict(
    "hurley2000" => 0,
    "belczynski2002" => 1,
    "belczynski2008" => 2,
    "fryer2012_rapid" => 3,
    "fryer2012_delayed" => 4,
    "eldridge2004" => 5,
)
"""`mdflag` per `winds` value (`mlwind.f`)."""
const _WINDS_MDFLAG = Dict(
    "hurley2000" => 1,
    "hurley2000_lbv" => 2,
    "belczynski2010" => 3,
    "belczynski2010_no_bistability" => 4,
)
"""`ecflag` per `electron_capture` value."""
const _ELECTRON_CAPTURE_ECFLAG = Dict("off" => 0, "on" => 1)
"""`bhflag` per `bh_kicks` value (`kick.F`)."""
const _BH_KICKS_BHFLAG = Dict("none" => 0, "full" => 1, "fallback" => 2)
"""`kmech` per `kick_mechanism` value (`kick.F`; 4 is a disabled branch and is not offered)."""
const _KICK_MECHANISM_KMECH = Dict("momentum" => 1, "convection" => 2, "collapse" => 3)

"""Enumerated keys of `[merger.stellar]` that write one engine parameter, with their value maps."""
const _SE_ENUM_KEYS = (
    (:pair_instability, "psflag", _PAIR_INSTABILITY_PSFLAG),
    (:remnant_mass, "nsflag", _REMNANT_MASS_NSFLAG),
    (:winds, "mdflag", _WINDS_MDFLAG),
    (:electron_capture, "ecflag", _ELECTRON_CAPTURE_ECFLAG),
    (:bh_kicks, "bhflag", _BH_KICKS_BHFLAG),
    (:kick_mechanism, "kmech", _KICK_MECHANISM_KMECH),
)

"""Real-valued keys of `[merger.stellar]` that write one engine parameter (`NaN` = the Level value)."""
const _SE_REAL_KEYS = (
    (:ns_kick_sigma, "disp"),
    (:ecsn_kick_sigma, "ecsig"),
    (:eddington_factor, "eddfac"),
    (:bh_accreted_fraction, "fctorcl"),
)

"""`true` when `s` sets any stellar-evolution option, i.e. writes into the three groups."""
function _sets_se_options(s::StellarSpec)::Bool
    s.option_set == "level" || return true
    any(k -> getfield(s, k[1]) != "level", _SE_ENUM_KEYS) && return true
    return any(k -> !isnan(getfield(s, k[1])::Float64), _SE_REAL_KEYS)
end

"""
    stellar_namelist_values(stellar, nbody6 = Nbody6ParameterSpec())
        -> Dict{String,Dict{String,Union{Int,Float64}}}

The values the generator writes into the `&INSSE`, `&INBSE` and `&INCOLL` groups of
`merger.inp`, keyed by group name (`"INSSE"`, `"INBSE"`, `"INCOLL"`) and the engine's
parameter names. Precedence: `stellar.option_set` (`"baseline"` writes
`_SE_BASELINE_SET`), then the explicit keys of `[merger.stellar]`, then the recoil and
spin keys of `[merger.nbody6]` (`gw_recoil` → `kicktype`, `bh_spin` → `bhspinfl`, written
only when they depart from the Level value as before). Every parameter not listed keeps
its Level value in the engine.

# Example
```julia
v = stellar_namelist_values(StellarSpec(; option_set = "baseline"))
v["INSSE"]["psflag"] == 2 && v["INBSE"]["eddfac"] == 1.0
```
"""
function stellar_namelist_values(
    s::StellarSpec,
    p::Nbody6ParameterSpec = Nbody6ParameterSpec(),
)::Dict{String,Dict{String,_SEValue}}
    set = Dict{String,_SEValue}()
    s.option_set == "baseline" && for (k, v) in _SE_BASELINE_SET
        set[k] = v
    end
    for (field, name, map) in _SE_ENUM_KEYS
        v = getfield(s, field)::String
        v == "level" || (set[name] = map[v])
    end
    for (field, name) in _SE_REAL_KEYS
        v = getfield(s, field)::Float64
        isnan(v) || (set[name] = v)
    end
    p.gw_recoil == "off" || (set["kicktype"] = _GW_RECOIL_KICKTYPE[p.gw_recoil])
    p.bh_spin == "zero" || (set["bhspinfl"] = _BH_SPIN_FLAG[p.bh_spin])
    groups = Dict{String,Dict{String,_SEValue}}()
    for g in keys(_SE_GROUPS)
        members = getfield(_SE_GROUPS, g)
        groups[String(g)] = Dict{String,_SEValue}(k => v for (k, v) in set if k in members)
    end
    return groups
end

"""
    stellar_parameters(level, written) -> Dict{String,Union{Int,Float64}}

Every parameter of `READSE` as the engine holds it after reading the groups: the Level
table of `level` (`"A"`, `"B"`, `"C"`), overwritten by `written` (a flat
`Dict(name => value)` or the per-group form of [`stellar_namelist_values`](@ref)). Level
`"0"` reads no group and sets nothing: an empty `Dict` is returned.
"""
function stellar_parameters(level::AbstractString, written::AbstractDict)
    level == "0" && return Dict{String,_SEValue}()
    haskey(_SE_LEVEL_TABLE, level) || throw(
        ArgumentError(
            "stellar_parameters: level must be one of \"A\", \"B\", \"C\", \"0\"; got $level",
        ),
    )
    params = merge(Dict{String,_SEValue}(), _SE_LEVEL_INDEPENDENT, _SE_LEVEL_TABLE[level])
    for (k, v) in written
        if v isa AbstractDict
            for (kk, vv) in v
                params[String(kk)] = vv
            end
        else
            params[String(k)] = v
        end
    end
    return params
end

stellar_parameters(s::StellarSpec, p::Nbody6ParameterSpec = Nbody6ParameterSpec()) =
    stellar_parameters(s.level, stellar_namelist_values(s, p))

"""A real namelist value with round-trip precision (`repr`: `265.0`, `0.001`, `3.920659e8`)."""
_namelist_real(x::Real) = repr(Float64(x))

"""
    _namelist_group_line(name, values) -> String

One namelist group as `merger.inp` carries it: `&NAME a=1,b=2.0 /`, members in the order of
`_SE_GROUPS`; `&NAME /` when `values` is empty.
"""
function _namelist_group_line(name::AbstractString, values::AbstractDict)::String
    members = getfield(_SE_GROUPS, Symbol(name))
    items = String[]
    for k in members
        haskey(values, k) || continue
        v = values[k]
        push!(items, k * "=" * (v isa Integer ? string(v) : _namelist_real(v)))
    end
    return isempty(items) ? "&$name /" : "&$name " * join(items, ",") * " /"
end

# -----------------------------------------------------------------------------
# The engine's echo
# -----------------------------------------------------------------------------

"""Numbers of a Fortran output field (`D` exponents accepted); `nothing` when `s` is none."""
function _fortran_number(s::AbstractString)
    t = replace(s, r"[dD]" => "e")
    v = tryparse(Int, t)
    v === nothing || return v
    return tryparse(Float64, t)
end

"""
    _parse_namelist_text(text) -> Dict{String,Union{Int,Float64}}

The `name=value` assignments of the body of a namelist group (names lower-cased, Fortran
`D` exponents read). Assignments that are not numbers are kept out and returned as
`"<name>"` with value `NaN`, so the caller can report them.
"""
function _parse_namelist_text(text::AbstractString)::Dict{String,_SEValue}
    out = Dict{String,_SEValue}()
    for m in eachmatch(r"([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([^,\s/]+)", text)
        name = lowercase(String(something(m.captures[1])))
        v = _fortran_number(String(something(m.captures[2])))
        out[name] = v === nothing ? NaN : v
    end
    return out
end

"""
    requested_stellar_parameters(input_path) -> NamedTuple

What an engine input requests of `READSE`: `(; level, written, unknown, parameters)` with
`level` the `Level` of the `&ININPUT` group (`""` when absent), `written` the assignments
of the `&INSSE`, `&INBSE` and `&INCOLL` groups (flat, engine names), `unknown` the assigned
names that are no member of their group or whose value is not a number (the engine's
namelist read fails on them), and `parameters` the resulting full set
([`stellar_parameters`](@ref); empty for Level `'0'` or without a Level).
"""
function requested_stellar_parameters(input_path::AbstractString)
    text = read(input_path, String)
    ml = match(r"Level\s*=\s*'([^']*)'"i, text)
    level = ml === nothing ? "" : String(something(ml.captures[1]))
    written = Dict{String,_SEValue}()
    unknown = String[]
    for g in keys(_SE_GROUPS)
        m = match(Regex("&" * String(g) * "\\b([^/]*)/", "i"), text)
        m === nothing && continue
        members = getfield(_SE_GROUPS, g)
        for (k, v) in _parse_namelist_text(String(something(m.captures[1])))
            if !(k in members) || (v isa Float64 && isnan(v))
                push!(unknown, "$(g).$(k)")
            else
                written[k] = v
            end
        end
    end
    parameters =
        haskey(_SE_LEVEL_TABLE, level) ? stellar_parameters(level, written) :
        Dict{String,_SEValue}()
    return (; level, written, unknown = sort!(unknown), parameters)
end

"""
    read_stellar_echo(stdout_path; offset = 0) -> NamedTuple

The engine's echo of the SSE/BSE/collision parameters (`READSE`, `input.F`) in the stdout
capture from byte `offset` on: `(; complete, values, invalid, level_zero)`. `complete` is
`true` once the `COLL PARAMS` line is present; `values` maps the 33 echoed parameters to
their printed values (three significant digits for the reals; empty until complete or
when a field is unreadable); `invalid` lists the groups the engine reported as an
`invalid Namelist format` (`"insse"`, `"inbse"`, `"incoll"`); `level_zero` is `true` when
the engine printed `Run without stellar evolution` (Level `'0'`: no group read, no echo).
"""
function read_stellar_echo(stdout_path::AbstractString; offset::Integer = 0)
    values = Dict{String,_SEValue}()
    invalid = String[]
    level_zero = false
    complete = false
    isfile(stdout_path) || return (; complete, values, invalid, level_zero)
    lines = open(stdout_path) do io
        offset > 0 && seek(io, offset)
        readlines(io)
    end
    numbers = _SEValue[]
    in_block = false
    readable = true
    for line in lines
        if !in_block
            mi = match(r"^\s*(insse|inbse|incoll):\s*invalid Namelist format", line)
            mi === nothing || push!(invalid, String(something(mi.captures[1])))
            occursin("Run without stellar evolution", line) && (level_zero = true)
            occursin("SSE FLAGS:", line) || continue
            in_block = true
        end
        segments = split(line, '=')
        for seg in segments[2:end]
            for tok in split(seg)
                v = _fortran_number(tok)
                if v === nothing
                    # A field of asterisks is an overflowed value; any other token ends the field list.
                    occursin(r"^\*+$", tok) && (readable = false)
                    break
                end
                push!(numbers, v)
            end
        end
        if occursin("COLL PARAMS:", line)
            complete = true
            break
        end
    end
    if complete && readable && length(numbers) == length(_SE_ECHO_ORDER)
        for (k, v) in zip(_SE_ECHO_ORDER, numbers)
            values[k] = v
        end
    end
    return (; complete, values, invalid, level_zero)
end

"""Echo `e` agrees with the requested `r`: integers exactly, reals at the three significant digits of the engine's `1P D11.2` field."""
_se_agrees(e::Integer, r::Integer) = e == r
_se_agrees(e::Real, r::Real) = isapprox(Float64(e), Float64(r); rtol = 5.01e-3, atol = 0.0)

"""
    check_stellar_echo(input_path, stdout_path; offset = 0) -> NamedTuple

Compare the SSE/BSE/collision parameters an engine input requests
([`requested_stellar_parameters`](@ref)) with the engine's echo of them in the stdout
capture from byte `offset` on ([`read_stellar_echo`](@ref)). Returns
`(; status, mismatches)`:

- `"passed"`: every echoed parameter equals the requested value (integers exactly, reals
  to the echo's three significant digits);
- `"mismatch"`: an echoed value differs, the input assigns a name no group defines, or
  the echo is unreadable; `mismatches` names each case;
- `"invalid_namelist"`: the engine reported a group as an invalid namelist and continued
  with the Level values;
- `"pending"`: the echo is not complete yet;
- `"not_applicable"`: Level `'0'` or no `Level` in the input (no group is read).

`READSE` reports a failed group read only by a line on stdout and goes on with the Level
values, so this comparison is the only evidence that a run used the option set it was
generated with.
"""
function check_stellar_echo(
    input_path::AbstractString,
    stdout_path::AbstractString;
    offset::Integer = 0,
)
    req = requested_stellar_parameters(input_path)
    haskey(_SE_LEVEL_TABLE, req.level) ||
        return (; status = "not_applicable", mismatches = String[])
    echo = read_stellar_echo(stdout_path; offset)
    isempty(echo.invalid) || return (;
        status = "invalid_namelist",
        mismatches = String["$(g): invalid Namelist format" for g in echo.invalid],
    )
    echo.complete || return (; status = "pending", mismatches = String[])
    mismatches = String["$(u): not a numerical member of the group" for u in req.unknown]
    if isempty(echo.values)
        push!(mismatches, "echo unreadable (overflowed or missing fields)")
    else
        for k in _SE_ECHO_ORDER
            r = req.parameters[k]
            e = echo.values[k]
            _se_agrees(e, r) || push!(mismatches, "$k: requested $r, engine $e")
        end
    end
    return (; status = isempty(mismatches) ? "passed" : "mismatch", mismatches)
end

"""
    _start_echo_check(input_path, stdout_path, process; offset = 0, interrupt_to = nothing,
                      poll = 0.5, on_result = nothing) -> (; stop, result, fired, task)

Asynchronous check of the engine's echo of the SSE/BSE/collision parameters
([`check_stellar_echo`](@ref)). It polls the stdout capture from byte `offset` on until
the echo is complete or the engine has exited; a `"mismatch"` or `"invalid_namelist"`
kills the engine's process group at once (SIGKILL rather than a stop request, so the engine
integrates nothing with the wrong set; the start-up dump at t = 0 may exist already, and a
resume from it is checked the same way, since `MODIFY` reads the groups again) and sets
`fired`. `result[]` holds the last outcome (`status = "pending"` while the engine has not
echoed; `"absent"` when it exited without an echo). `on_result`, when given, is called
with the outcome once the echo is judged (`passed`, `mismatch`, `invalid_namelist`), before
a mismatching engine is killed; a failure of it is logged and ignored. Setting `stop` ends
the task quietly.
"""
function _start_echo_check(
    input_path::AbstractString,
    stdout_path::AbstractString,
    process::Base.Process;
    offset::Integer = 0,
    interrupt_to::Union{Nothing,Task} = nothing,
    poll::Real = 0.5,
    on_result::Union{Nothing,Function} = nothing,
)
    stop = Ref(false)
    report(r) =
        on_result === nothing || try
            on_result(r)
        catch err
            @warn "Stellar-option echo: the outcome could not be recorded: $(sprint(showerror, err))"
        end
    fired = Ref(false)
    result = Ref{Any}((; status = "pending", mismatches = String[]))
    task = errormonitor(
        _forwarding_interrupts(interrupt_to) do
            while !stop[]
                running = process_running(process)
                r = check_stellar_echo(input_path, stdout_path; offset)
                result[] = r
                r.status in ("passed", "mismatch", "invalid_namelist") && report(r)
                if r.status in ("mismatch", "invalid_namelist")
                    fired[] = true
                    # SIGKILL of the group, not a stop request: the engine must not
                    # integrate, nor dump, with the wrong option set.
                    running && try
                        _kill_group(getpid(process), Base.SIGKILL)
                    catch
                        process_running(process) && rethrow()
                    end
                    @error "Stellar-evolution options: the engine does not run the requested " *
                           "set; killed at start-up\n  " *
                           join(r.mismatches, "\n  ")
                    return nothing
                end
                r.status == "pending" || return nothing
                if !running
                    result[] = (; status = "absent", mismatches = String[])
                    return nothing
                end
                sleep(poll)
            end
        end,
    )
    return (; stop, result, fired, task)
end
