# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "GPU validation driver" begin
    # The dependency check accepts an nvcc that sits under the toolkit rather
    # than on PATH (the build exports <cuda_path>/bin itself).
    cuda_dir = mktempdir()
    mkpath(joinpath(cuda_dir, "bin"))
    touch(joinpath(cuda_dir, "bin", "nvcc"))
    dep_cfg(body) = (p = joinpath(cuda_dir, "cfg.toml"); write(p, body); load_config(p))
    @test !(
        "nvcc" in
        check_dependencies(dep_cfg("[build]\nenable_gpu = true\ncuda_path = \"$cuda_dir\"\n"))
    )
    @test !("nvcc" in check_dependencies(dep_cfg("[build]\nenable_gpu = false\n")))
    if !Nbody6Dynamics.check_command("nvcc")
        absent = joinpath(cuda_dir, "absent")
        @test "nvcc" in
              check_dependencies(dep_cfg("[build]\nenable_gpu = true\ncuda_path = \"$absent\"\n"))
    end

    # nvcc host-compiler probe: classification of the compiler message
    @test Nbody6Dynamics._unsupported_host_compiler(
        "#error -- unsupported GNU version! gcc versions later than 14 are not supported!",
    )
    @test Nbody6Dynamics._unsupported_host_compiler("unsupported clang version")
    @test !Nbody6Dynamics._unsupported_host_compiler("probe.cu(1): error: expected a \";\"")
    # No nvcc under the path: the probe stands aside and the build reports it
    @test Nbody6Dynamics._nvcc_host_compiler_flags("/nonexistent/cuda") == String[]

    # Logged stage runner: merged output, exit code, escape sequences stripped
    vdir = mktempdir()
    log = joinpath(vdir, "stage.log")
    julia = Base.julia_cmd()
    @test Nbody6Dynamics._tee_run(
        `$julia -e 'print("\e[31mred\e[0m\n"); println(stderr, "err line")'`,
        log,
    ) == (exitcode = 0, signal = 0)
    lines = readlines(log)
    @test "red" in lines && "err line" in lines
    @test Nbody6Dynamics._tee_run(`$julia -e 'exit(3)'`, log) == (exitcode = 3, signal = 0)
    # A genuine exit code ≥ 128 is not a signal.
    @test Nbody6Dynamics._tee_run(`sh -c 'exit 139'`, log) == (exitcode = 139, signal = 0)
    # A process killed by a signal reports exitcode 0 through libuv, so
    # the runner reports the signal separately or a crash is recorded as
    # a success.
    @test Nbody6Dynamics._tee_run(`sh -c 'kill -s SEGV $$'`, log) == (exitcode = 0, signal = 11)
    # The child runs in its own session (detach), so the hangup of a terminal
    # closing above the driver never reaches it: it is its own session leader.
    Sys.islinux() && @test Nbody6Dynamics._tee_run(
        `sh -c '[ "$(ps -o sid= -p $$ | tr -d " ")" = "$$" ]'`,
        log,
    ) == (exitcode = 0, signal = 0)

    # An operator interrupt is forwarded to a detached stage (SIGINT first).
    child = run(detach(`sleep 30`); wait = false)
    Nbody6Dynamics._interrupt_child(child; grace = 5.0)
    @test !process_running(child) && Nbody6Dynamics._exit_status(child) == -2
    # Every entry script that starts a stage or an engine makes SIGINT an
    # exception, so that the forwarding above and the pipeline's handler run
    # at all: a plain script exits at once and orphans its detached children.
    root = normpath(joinpath(@__DIR__, ".."))
    for script in (
        "scripts/run_gpu_validation.jl",
        "scripts/run_setup.jl",
        "scripts/run_sweep.jl",
        "scripts/run_verif_suite.jl",
        "bench/gpu_scaling.jl",
        "bench/gpu_cells.jl",
        "bench/thread_scaling.jl",
    )
        @test occursin("Base.exit_on_sigint(false)", read(joinpath(root, script), String))
    end

    # A crash signal is retried up to max_retries, keeping every output
    sdir = mktempdir()
    slog = joinpath(sdir, "stage.log")
    segv = `sh -c 'kill -s SEGV $$'`
    result, retried =
        @test_logs (:warn, r"killed by signal 11") Nbody6Dynamics._run_stage_with_retry(
            segv,
            slog,
            :suite,
        )
    @test result.signal == 11 && retried == [11]
    @test isfile(joinpath(sdir, "stage.attempt1.signal11.log"))
    result, retried =
        @test_logs (:warn, r"retry 1 of 2") (:warn, r"retry 2 of 2") Nbody6Dynamics._run_stage_with_retry(
            segv,
            joinpath(sdir, "twice.log"),
            :suite;
            max_retries = 2,
        )
    @test retried == [11, 11]
    @test isfile(joinpath(sdir, "twice.attempt2.signal11.log"))
    # A stage that crashes once and then succeeds is reported as passed,
    # with the signal it met on record
    marker = joinpath(sdir, "once")
    flaky = `sh -c "if [ -e $marker ]; then exit 0; else touch $marker; kill -s SEGV \$\$; fi"`
    result, retried =
        @test_logs (:warn, r"killed by signal 11") Nbody6Dynamics._run_stage_with_retry(
            flaky,
            joinpath(sdir, "flaky.log"),
            :suite,
        )
    @test result == (exitcode = 0, signal = 0) && retried == [11]
    # Not retried: a kill from outside, retries switched off, an ordinary failure
    result, retried = Nbody6Dynamics._run_stage_with_retry(`sh -c 'kill -s KILL $$'`, slog, :suite)
    @test result.signal == 9 && isempty(retried)
    result, retried = Nbody6Dynamics._run_stage_with_retry(segv, slog, :gpu; max_retries = 0)
    @test result.signal == 11 && isempty(retried)
    result, retried = Nbody6Dynamics._run_stage_with_retry(`$julia -e 'exit(3)'`, slog, :gpu)
    @test result.exitcode == 3 && isempty(retried)
    # The policy is validated at the public interface
    @test_throws ArgumentError run_gpu_validation(; max_retries = -1, dry_run = true)
    @test_throws ArgumentError run_gpu_validation(; retry_stages = [:plots], dry_run = true)
    @test_throws ArgumentError run_gpu_validation(; retry_signals = [0], dry_run = true)
    @test Nbody6Dynamics._tool_banner(`$julia -e 'println("banner line"); println("second")'`) ==
          "banner line"
    @test Nbody6Dynamics._tool_banner(`/nonexistent/tool --version`) == "unavailable"

    # Probe against a stand-in nvcc: a script that mimics the toolkit's
    # verdicts (rejects the default compiler, accepts an override or a
    # -ccbin), so the search logic runs without a CUDA toolkit.
    function fake_nvcc(dir, body)
        mkpath(joinpath(dir, "bin"))
        path = joinpath(dir, "bin", "nvcc")
        write(path, "#!/bin/sh\n" * body)
        chmod(path, 0o755)
        return dir
    end
    accepts_override = fake_nvcc(
        mktempdir(),
        "case \" \$* \" in *' -allow-unsupported-compiler '*) exit 0;; esac\n" *
        "echo '#error -- unsupported GNU version! gcc versions later than 15 are not supported!' >&2\nexit 1\n",
    )
    flags =
        @test_logs (:warn, r"building with -allow-unsupported-compiler") Nbody6Dynamics._nvcc_host_compiler_flags(
            accepts_override;
            candidates = String[],
        )
    @test flags == ["-allow-unsupported-compiler"]
    needs_ccbin = fake_nvcc(
        mktempdir(),
        "case \" \$* \" in *' -ccbin /usr/bin/true '*) exit 0;; esac\n" *
        "echo 'unsupported GNU version! gcc versions later than 15 are not supported!' >&2\n" *
        "for i in \$(seq 1 60); do echo \"type_traits(\$i): error: identifier char8_t is undefined\" >&2; done\nexit 1\n",
    )
    flags =
        @test_logs (:warn, r"building with -ccbin /usr/bin/true") Nbody6Dynamics._nvcc_host_compiler_flags(
            needs_ccbin;
            candidates = ["/nonexistent/g++-99", "/usr/bin/true"],
        )
    @test flags == ["-ccbin", "/usr/bin/true"]
    hopeless = try
        Nbody6Dynamics._nvcc_host_compiler_flags(needs_ccbin; candidates = String[])
    catch err
        err
    end
    @test hopeless isa ErrorException
    @test occursin("default host compiler; -allow-unsupported-compiler", hopeless.msg)
    @test occursin("gcc15-c++", hopeless.msg) && occursin("lines omitted", hopeless.msg)
    @test occursin("--- default host compiler ---", hopeless.msg) &&
          occursin("--- -allow-unsupported-compiler ---", hopeless.msg)
    # Every attempt's own output is reported, so a rejected -ccbin candidate
    # can be diagnosed from the message alone.
    echoing = fake_nvcc(
        mktempdir(),
        "echo \"nvcc args: \$*\" >&2\n" *
        "echo 'unsupported GNU version! gcc versions later than 15 are not supported!' >&2\nexit 1\n",
    )
    per_attempt = try
        Nbody6Dynamics._nvcc_host_compiler_flags(echoing; candidates = ["/nonexistent/g++-99"])
    catch err
        err
    end
    @test per_attempt isa ErrorException
    @test occursin(
        r"--- -ccbin /nonexistent/g\+\+-99 ---\nnvcc args: [^\n]* -ccbin /nonexistent/g\+\+-99\n",
        per_attempt.msg,
    )
    @test occursin(
        r"--- -allow-unsupported-compiler ---\nnvcc args: [^\n]* -allow-unsupported-compiler\n",
        per_attempt.msg,
    )
    # glibc ≥ 2.42 declares rsqrt/rsqrtf with an exception specification
    # the CUDA headers lack: the probe retries with the feature-macro
    # override, alone or on top of the host-compiler choice.
    glibc_line =
        "/usr/include/bits/mathcalls.h(206): error: exception specification is incompatible " *
        "with that of previous function \"rsqrt\" (declared at line 629 of crt/math_functions.h)"
    @test Nbody6Dynamics._glibc_c2y_conflict(glibc_line)
    @test !Nbody6Dynamics._glibc_c2y_conflict("unsupported GNU version! gcc versions later than 15")
    glibc_msg = "echo '$glibc_line' >&2\nexit 1\n"
    glibc_only = fake_nvcc(
        mktempdir(),
        "case \" \$* \" in *' -U_GNU_SOURCE -D_DEFAULT_SOURCE '*) exit 0;; esac\n" * glibc_msg,
    )
    flags =
        @test_logs (:warn, r"glibc declares rsqrt and rsqrtf") Nbody6Dynamics._nvcc_host_compiler_flags(
            glibc_only;
            candidates = String[],
        )
    @test flags == ["-U_GNU_SOURCE", "-D_DEFAULT_SOURCE"]
    # A configured -ccbin stays in nvcc_flags; only the override is added
    flags =
        @test_logs (:warn, r"glibc declares rsqrt and rsqrtf") Nbody6Dynamics._nvcc_host_compiler_flags(
            glibc_only;
            nvcc_flags = ["-ccbin", "/usr/bin/true"],
            candidates = String[],
        )
    @test flags == ["-U_GNU_SOURCE", "-D_DEFAULT_SOURCE"]
    # Host compiler rejected AND the glibc conflict (Fedora 44 with CUDA
    # 13.1): the first accepted -ccbin candidate, with the override
    both = fake_nvcc(
        mktempdir(),
        "case \" \$* \" in\n" *
        "  *' -ccbin /usr/bin/true -U_GNU_SOURCE -D_DEFAULT_SOURCE '*) exit 0;;\n" *
        "  *' -ccbin /usr/bin/true '*) " *
        glibc_msg *
        ";;\n" *
        "esac\n" *
        "echo 'unsupported GNU version! gcc versions later than 15 are not supported!' >&2\nexit 1\n",
    )
    flags = @test_logs (:warn, r"building with -ccbin /usr/bin/true") (
        :warn,
        r"glibc declares rsqrt and rsqrtf",
    ) Nbody6Dynamics._nvcc_host_compiler_flags(
        both;
        candidates = ["/nonexistent/g++-99", "/usr/bin/true"],
    )
    @test flags == ["-ccbin", "/usr/bin/true", "-U_GNU_SOURCE", "-D_DEFAULT_SOURCE"]
    # The override does not help: every attempt is reported, with the
    # header-patch advice
    glibc_stuck = fake_nvcc(mktempdir(), glibc_msg)
    stuck = try
        Nbody6Dynamics._nvcc_host_compiler_flags(glibc_stuck; candidates = String[])
    catch err
        err
    end
    @test stuck isa ErrorException
    @test occursin("default host compiler with -U_GNU_SOURCE -D_DEFAULT_SOURCE", stuck.msg)
    @test occursin("math_functions.h", stuck.msg) && occursin("noexcept(true)", stuck.msg)
    stuck_pinned = try
        Nbody6Dynamics._nvcc_host_compiler_flags(
            glibc_stuck;
            nvcc_flags = ["-ccbin", "/usr/bin/false"],
            candidates = String[],
        )
    catch err
        err
    end
    @test stuck_pinned isa ErrorException && occursin("configured host compiler", stuck_pinned.msg)
    @test occursin("math_functions.h", stuck_pinned.msg)
    # Host-record verdict of the probe
    @test Nbody6Dynamics._nvcc_probe_record("/nonexistent/cuda", "") ==
          (String[], "skipped: nvcc unavailable")
    flags, verdict = Nbody6Dynamics._nvcc_probe_record(echoing, "13.1")
    @test flags == String[] && occursin("--- default host compiler ---", verdict)
    flags, verdict =
        @test_logs (:warn, r"building with -allow-unsupported-compiler") Nbody6Dynamics._nvcc_probe_record(
            accepts_override,
            "13.1",
        )
    @test flags == ["-allow-unsupported-compiler"] && verdict == "passed"
    # A configured -ccbin is final: no search, the failure is reported as is
    pinned_fail = try
        Nbody6Dynamics._nvcc_host_compiler_flags(
            needs_ccbin;
            nvcc_flags = ["-ccbin", "/usr/bin/false"],
            candidates = ["/usr/bin/true"],
        )
    catch err
        err
    end
    @test pinned_fail isa ErrorException && occursin("configured host compiler", pinned_fail.msg)
    @test Nbody6Dynamics._host_compiler_candidates() isa Vector{String}
    @test Nbody6Dynamics._output_excerpt("a\nb\nc", 1, 1) == "a\n… (1 lines omitted)\nc"
    @test Nbody6Dynamics._output_excerpt("a\nb", 5, 5) == "a\nb"
    @test Nbody6Dynamics._run_capture(
        `$(Base.julia_cmd()) -e 'println(stderr, "e"); print("o")'`,
        joinpath(vdir, "cap.log"),
    ) == (true, "e\no")

    # Dry run: host record and planned commands, nothing executed
    @test_throws ArgumentError run_gpu_validation(;
        base_dir = vdir,
        stages = [:nope],
        dry_run = true,
    )
    @test_throws ArgumentError run_gpu_validation(;
        base_dir = vdir,
        stages = Symbol[],
        dry_run = true,
    )
    @test_throws ArgumentError run_gpu_validation(;
        base_dir = vdir,
        bench_tcrit = 0.0,
        dry_run = true,
    )
    out = run_gpu_validation(;
        base_dir = vdir,
        stages = [:suite, :gpu, :bench],
        bench_n = [1000],
        bench_threads = [2],
        bench_gpu_lists = [[0], [0, 1]],
        bench_tcrit = 0.5,
        dry_run = true,
    )
    @test startswith(basename(out), "gpu_validation_") && dirname(out) == joinpath(vdir, "runs")
    host = Nbody6Dynamics.TOML.parsefile(joinpath(out, "HOST_INFO.toml"))
    for key in (
        "host",
        "gpu",
        "compute_capabilities",
        "nvcc_release",
        "gcc",
        "gfortran",
        "package_commit",
        "glibc",
        "host_compilers",
        "nvcc_host_flags",
        "nvcc_probe",
    )
        @test haskey(host, key)
    end
    @test host["nvcc_probe"] isa String && host["nvcc_host_flags"] isa Vector
    summary = Nbody6Dynamics.TOML.parsefile(joinpath(out, "VALIDATION.toml"))
    @test summary["dry_run"] && summary["stages"] == ["suite", "gpu", "bench"]
    @test summary["results"]["suite"]["status"] == "planned"
    @test occursin("NBODY6_GPU_TESTS=1", summary["results"]["suite"]["command"])
    bench_cmd = summary["results"]["bench"]["command"]
    @test occursin("gpu_scaling.jl 1000 2 0;0,1 0.5", bench_cmd)
    @test occursin("NBODY6_GPU_BACKEND=", bench_cmd) &&
          occursin("Nbody6PPGPU-beijing-gpu", bench_cmd)
    @test !haskey(summary["results"], "cpu") && isempty(filter(endswith(".log"), readdir(out)))
    # The pipeline stages are told the run ID whose artefacts decide their verdict.
    @test haskey(summary, "run_ids") &&
          occursin("--run-id=gpu_validation_", summary["results"]["gpu"]["command"])
    # A stage whose prerequisites are missing is skipped with the reason, no process spawned
    bare = mktempdir()
    skipped = run_gpu_validation(; base_dir = bare, stages = [:bench])
    skipped_summary = Nbody6Dynamics.TOML.parsefile(joinpath(skipped, "VALIDATION.toml"))
    @test skipped_summary["results"]["bench"]["status"] == "skipped"
    @test occursin("Nbody6PPGPU-beijing-gpu", skipped_summary["results"]["bench"]["reason"])
    @test !isfile(joinpath(skipped, "bench.log")) && haskey(skipped_summary, "finished")
    @test Nbody6Dynamics._stage_prerequisite(:suite, bare) === nothing
    # The probe stages are opt-in pipelines with their own run IDs
    probes = run_gpu_validation(;
        base_dir = vdir,
        stages = [:gpu_merger_600k, :cpu_single_600k],
        dry_run = true,
    )
    probe_summary = Nbody6Dynamics.TOML.parsefile(joinpath(probes, "VALIDATION.toml"))
    @test probe_summary["results"]["gpu_merger_600k"]["status"] == "planned"
    @test probe_summary["results"]["cpu_single_600k"]["status"] == "planned"
    probe_cmd = probe_summary["results"]["gpu_merger_600k"]["command"]
    @test occursin("gpu_merger_600k.toml", probe_cmd) &&
          occursin("--run-id=gpu_validation_", probe_cmd)
    @test haskey(probe_summary["run_ids"], "gpu_merger_600k") &&
          haskey(probe_summary["run_ids"], "cpu_single_600k")
    @test probe_summary["stop_on_failure"] == false
    # Each probe needs the tree its binary is built in
    @test occursin(
        "Nbody6PPGPU-beijing-gpu",
        Nbody6Dynamics._stage_prerequisite(:gpu_merger_600k, bare),
    )
    cpu_unmet = Nbody6Dynamics._stage_prerequisite(:cpu_single_600k, bare)
    @test occursin("Nbody6PPGPU-beijing", cpu_unmet) && !occursin("-gpu", cpu_unmet)
    # Gated chain: the :cpu stage fails (no scripts/run_setup.jl under bare),
    # so the probe after it is skipped for that reason and never spawned
    chain = run_gpu_validation(;
        base_dir = bare,
        stages = [:cpu, :cpu_single_600k],
        stop_on_failure = true,
    )
    chain_summary = Nbody6Dynamics.TOML.parsefile(joinpath(chain, "VALIDATION.toml"))
    chain_results = chain_summary["results"]
    @test chain_results["cpu"]["status"] == "failed"
    @test chain_results["cpu_single_600k"]["status"] == "skipped"
    @test occursin("after failed stage :cpu", chain_results["cpu_single_600k"]["reason"])
    @test !isfile(joinpath(chain, "cpu_single_600k.log"))
    @test chain_summary["stop_on_failure"] == true
    # Without the gate the probe is judged on its own prerequisite
    open_chain = run_gpu_validation(; base_dir = bare, stages = [:cpu, :cpu_single_600k])
    open_results = Nbody6Dynamics.TOML.parsefile(joinpath(open_chain, "VALIDATION.toml"))["results"]
    @test open_results["cpu"]["status"] == "failed"
    @test open_results["cpu_single_600k"]["status"] == "skipped"
    open_reason = open_results["cpu_single_600k"]["reason"]
    @test occursin("Nbody6PPGPU-beijing", open_reason) &&
          !occursin("after failed stage", open_reason)
    mkpath(joinpath(bare, "backend", "Nbody6PPGPU-beijing", "build"))
    @test Nbody6Dynamics._stage_prerequisite(:cpu_single_600k, bare) === nothing
    # A host whose compiler probe did not pass runs nothing when a stage needs nvcc
    if !Nbody6Dynamics.check_command("nvcc")
        bare2 = mktempdir()
        gated =
            run_gpu_validation(; base_dir = bare2, stages = [:gpu, :cpu], stop_on_failure = true)
        gated_host = Nbody6Dynamics.TOML.parsefile(joinpath(gated, "HOST_INFO.toml"))
        # A toolkit off PATH may still pass the probe; the gate is then open.
        if gated_host["nvcc_probe"] != "passed"
            gated_results =
                Nbody6Dynamics.TOML.parsefile(joinpath(gated, "VALIDATION.toml"))["results"]
            for s in ("gpu", "cpu")
                @test gated_results[s]["status"] == "skipped"
                @test occursin("host-compiler probe", gated_results[s]["reason"])
                @test !isfile(joinpath(gated, "$s.log"))
            end
        end
    end
    # Benchmark artefact collection on an empty bench tree is a no-op
    @test Nbody6Dynamics._collect_bench_artefacts(joinpath(vdir, "bench"), out, 0.0) == String[]
end

# =====================================================================

@testset "Validation stage completeness check" begin
    mktempdir() do base
        runs = joinpath(base, "runs")
        mkpath(runs)
        run_id = "merger_cpu_20260911_194002_b4b0"

        # The suite and benchmark stages are judged by exit code alone.
        @test Nbody6Dynamics._stage_incomplete(:suite, base, run_id) === nothing
        @test Nbody6Dynamics._stage_incomplete(:bench, base, run_id) === nothing

        # A pipeline stage that produced nothing at all.
        reason = Nbody6Dynamics._stage_incomplete(:cpu, base, run_id)
        @test reason !== nothing
        @test occursin("no run directory", reason)

        # A run directory whose summary stops at the engine phase.
        run_dir = joinpath(runs, run_id)
        mkpath(run_dir)
        open(joinpath(run_dir, "RUN_INFO.toml"), "w") do io
            Nbody6Dynamics.TOML.print(io, Dict("run" => Dict("id" => run_id)))
        end
        reason = Nbody6Dynamics._stage_incomplete(:cpu, base, run_id)
        @test reason !== nothing
        @test occursin("completed", reason)

        # Once the pipeline stamps completion the stage is accepted.
        @test Nbody6Dynamics._stamp_pipeline_completion(run_dir, ["simulation"], 1.0)
        @test Nbody6Dynamics._stage_incomplete(:cpu, base, run_id) === nothing

        # Only the run directory the stage was assigned counts.
        @test Nbody6Dynamics._stage_incomplete(:cpu, base, "absent_id") !== nothing
        @test Nbody6Dynamics._stage_incomplete(:cpu_single_600k, base, "absent_id") !== nothing

        # The pipeline also completes on partial output: an engine that
        # ended without END RUN fails the stage despite the marker.
        info = joinpath(run_dir, "RUN_INFO.toml")
        write_summary(completed) = open(info, "w") do io
            Nbody6Dynamics.TOML.print(
                io,
                Dict(
                    "run" => Dict("id" => basename(run_dir)),
                    "segments" => [Dict("index" => 1, "completed" => completed)],
                ),
            )
        end
        write_summary(false)
        @test Nbody6Dynamics._engine_completed(run_dir) === false
        @test Nbody6Dynamics._stamp_pipeline_completion(run_dir, ["simulation"], 1.0)
        @test Nbody6Dynamics.TOML.parsefile(info)["pipeline"]["engine_completed"] === false
        reason = Nbody6Dynamics._stage_incomplete(:cpu, base, run_id)
        @test reason !== nothing
        @test occursin("END RUN", reason)

        write_summary(true)
        @test Nbody6Dynamics._engine_completed(run_dir) === true
        @test Nbody6Dynamics._stamp_pipeline_completion(run_dir, ["simulation"], 1.0)
        @test Nbody6Dynamics._stage_incomplete(:cpu, base, run_id) === nothing

        # No segment (post-processing only): nothing to judge.
        @test Nbody6Dynamics._engine_completed(mktempdir()) === nothing
    end
end

# =====================================================================
