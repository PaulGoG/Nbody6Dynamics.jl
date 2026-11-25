# Implementation notes

The package's corrections to the engine source, its probe of the CUDA host compiler and the path of a signal from a terminal or a scheduler to the engine are described here, apart from the [User Manual](manual.md), because none of it is needed to configure or run a simulation. It serves for reading a run record closely, for bringing up a new host and for changing the code.

## Build

### HDF5 build patching

With `enable_hdf5 = true` the build adds the HDF5 compiler and linker flags and `-DCONFIG_HDF5` to the engine Makefile. This engine version reads no `CONFIG_HDF5` symbol: its H5Part writer is compiled only under `H5OUTPUT`, which the package does not define, so the setting changes nothing in what the binary writes. The package reads no `.h5part` file in any case; `scan_output` reports them.

### Source amendments

The package's engine fork carries the corrections below as commits on its `main` branch (the default `install.source_url` and `install.ref`), and a build from it finds them in place. For a build from an upstream checkout, the install phase applies the same corrections before `configure` from the package's `deps/engine/amendments.toml`. Each is an exact text substitution in one file:

| Identifier | File | Correction |
|---|---|---|
| `escape-summary-format` | `src/Main/escape.F` | The escaper summary writes its name list in a repeatable format group (`1000(I10)` for `1000I10`). Without it a step that removes more than 1000 escapers ends the run with a format error; up to 1000 the output is the same |
| `tidal-energy-total` | `src/Main/energy.F` | The whole-system call of `xtrnlv.F` is restored for a linearised external field (`CALL XTRNLV(1,N)` in a branch for `KZ(14) ≤ 2`; the call is commented out in this engine version). The potential energy of the field then enters `ETIDE`, the total energy and the energy check, and the external virial term enters `Q`, as in the original NBODY6. The orbit modes (`KZ(14) = 3, 5`) keep their code path: on the Milky-Way potential the restored whole-system call corrupts the energy bookkeeping at the first adjustment |
| `tidal-energy-total-mpi` | `src/Main/energy_mpi.F` | The same restoration in the MPI variant of the energy evaluation |

The step is idempotent (a substitution found in place counts as applied), and the identifiers in effect are recorded as `source_amendments` in `BUILD_INFO.toml`, hence in the `[build]` table of every run. An amendment whose text does not occur exactly once in the tree, as with another engine revision, is skipped with a warning and listed under `source_amendments_skipped`; the build then proceeds on the unamended source. On an upstream checkout `backend_commit` carries `-dirty`: the substitutions are local modifications, and `configure` rewrites two tracked files (`config.log`, `include/params.h`). The fork tracks neither and ignores the build outputs, so a build from its tag records a clean commit, and `-dirty` there means a hand-edited tree. Further corrections are made as commits of the fork, not as amendments; the fork's `NOTICE.md` lists all of them by tag (among them the hang of an output interval without a short exact binary expansion, the end-of-run dumps skipped under `NCOMM > 1`, and the thread limit of the AVX build). An existing checkout under `install_dir` is reused as it is (whatever its origin); `reinstall = true` replaces it with a fresh clone of `source_url` at `ref`.

### Host-compiler probe

The build probes `nvcc` on a trivial kernel first: when the toolkit rejects the host compiler it tries `-allow-unsupported-compiler`, then `-ccbin` with `CUDAHOSTCXX` and the versioned `g++-15` … `g++-12`, `clang++` found on `PATH`, and records what it settled on; a configured `-ccbin` is used as is. When the host's glibc (≥ 2.42) declares `rsqrt`/`rsqrtf` with an exception specification the CUDA headers lack, the probe retries the host-compiler choice with `-U_GNU_SOURCE -D_DEFAULT_SOURCE` and keeps it. When no combination works the error quotes the `nvcc` output of every attempt. The engine's `configure` finds `nvcc` only on `PATH` (its `--with-cuda` fallback reuses the cached result of the `PATH` check), which is why the exported `PATH` matters.

## Signals and interrupts

The engine runs in its own session, so neither the hangup nor the interrupt of a terminal above the pipeline reaches it. A stop request gets to it through the launch script and through the driver.

**Launch script.** The script stays the parent of the engine: it starts the engine (or `mpirun`) in a subshell that ignores SIGTERM, writes the engine's process id to `output/engine.pid` while the engine runs, and waits. A SIGTERM sent to the script becomes the `STOP` file of a stop request ([Wall budget and stop requests](@ref)); the script exits with the engine's exit status and re-raises on itself the signal that ended the engine, so the segment record shows the engine's end.

**SIGTERM.** A termination signal is handled in two places. The launch script turns it into the `STOP` file while the engine it started ignores it, so a job-wide signal stops the engine with a dump instead of killing it. The driver's exit hook writes `STOP` for every engine the process has running, waits up to `stop_margin` seconds for them to exit, terminates those that have not (SIGTERM to the script, SIGKILL of the process group 15 s later), and closes their segment records (`stop_requested = "signal"`). The engine runs in its own session, so without the hook it would outlive the driver.

The hook relies on Julia's handling of SIGTERM, which serves the signal wherever the main thread happens to be: with the driver waiting on the engine, its state for all but moments of a run, the stop is orderly; a signal that interrupts it in the first seconds after start-up, while code is still being compiled, can block the exit instead. The engine then stops all the same when the signal also reached the launch script, as a scheduler's does; a run whose driver blocked is left with a `running` record and is continued from the dump like one that crashed.

**SIGINT.** An interrupt of a pipeline process is forwarded to the engine as a stop request. The entry scripts set `Base.exit_on_sigint(false)` for that, because a Julia script exits on SIGINT at once by default; an interrupt that lands in one of the pipeline's helper tasks (telemetry sampler, watchdogs) is handed to the task waiting on the engine. The validation driver forwards SIGINT to its running stage the same way.

**Watchdogs.** The start-up watchdog and the completion monitor terminate the engine before they log the fact, so a pipeline whose console has gone (its launcher died) still ends the engine it decided to end, and writes its run record and file log.

## The interval hang of the upstream engine

An upstream build given an output interval without an exact binary expansion prints the t = 0 adjustment and never advances ([Run never advances past t = 0](@ref)). The cause is the upstream `string_left.f`, which counts the decimal digits of `DELTAT`, `DTADJ` and `DTPLOT` by multiplying by ten until the value is an integer, with a default-kind `int`. A value such as `0.6302` never becomes an integer in binary arithmetic; past 2³¹ the conversion overflows and the loop never ends, inside the very first output.
