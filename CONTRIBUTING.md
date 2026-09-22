# Contributing

I develop this package in the open and take issues and pull requests. This
page says what a change needs so that it can be merged without a second
round.

## Setting up

```bash
git clone https://github.com/PaulGoG/Nbody6Dynamics.jl.git
cd Nbody6Dynamics.jl
julia activate.jl            # resolves and precompiles the package environment
julia scripts/activate.jl    # the entry scripts' environment (adds the figure backend)
```

The test suite runs with

```bash
julia -e 'include("activate.jl"); Pkg.test()'
```

and needs no engine binary. The engine-dependent tests (build, launch,
restart, tidal field, telemetry) run when `NBODY6_BINARY_TESTS=1` is set; they
clone and build Nbody6++GPU in a temporary directory, or use an existing build
under `NBODY6_BACKEND_ROOT`. Julia 1.13 or later is required; a Linux host with
`gfortran` and `make` is needed for the engine.

## What a pull request needs

- The suite green, including the static checks it carries (Aqua, JET,
  ExplicitImports). Add tests for what the change does; remove tests for what
  it no longer does.
- Formatting per `.JuliaFormatter.toml` (`using JuliaFormatter; format(".")`).
- Docstrings for every public function or type, and the manual or the input
  file reference updated when behaviour or configuration keys change.
- An entry under `[Unreleased]` in `CHANGELOG.md`.
- Physics first: a change to a diagnostic, a sampler or a generated engine
  parameter says in its description what the formulation is and where it comes
  from (a reference with a DOI where one exists), and its tests check a
  limiting case or a published value.
- No run outputs, engine sources, manifests or figures in the commit; `runs/`,
  `backend/` and `Manifest.toml` are ignored for that reason.

## Reporting a problem

Open an issue with the configuration you ran, the engine commit
(`install.ref`, or `backend_commit` in `RUN_INFO.toml`), the Julia version,
the `[hardware]` table of `RUN_INFO.toml` when a run is involved, and the log
(`nbody6dynamics.log` in the run directory). A halted integration is often the
engine's own energy check; the input-file reference explains the tolerance.
