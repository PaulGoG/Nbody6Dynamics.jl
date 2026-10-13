# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Cluster jobs (Slurm layer)" begin
    # Stand-in for the engine, copied from the "resume_run" testset of
    # test_resume.jl. Integer time, one step per 0.1 s. Fresh start: t = 0 …
    # TCRIT. Restart (KSTART=2 in the input): t0 from the payload of comm.1,
    # end = t0 + TCRIT. STOP is honoured with a dump comm.1_<t>.0; END RUN is
    # followed by a dump comm.1_<end>.0.
    fake_engine = raw"""
    #!/bin/bash
    input=$(cat)
    tcrit=$(printf '%s' "$input" | grep -o 'TCRIT=[0-9.]*' | head -n 1 | cut -d= -f2)
    tcrit=${tcrit%.*}
    dump() {  # dump <file> <t>: two Fortran records, the first holding the time
        printf '\x04\x00\x00\x00%4d\x04\x00\x00\x00\x04\x00\x00\x00abcd\x04\x00\x00\x00' "$2" > "$1"
    }
    if printf '%s' "$input" | grep -q 'KSTART=2'; then
        t0=$(dd if=comm.1 bs=1 skip=4 count=4 2>/dev/null | tr -d ' ')
        dump "comm.1_$t0.0" "$t0"
        printf '  MYDUMP    %d.0000000000000000      %d.0                         201 comm.1_%d.0                \n' "$t0" "$t0" "$t0"
        printf '  W MYDUMP J,II,NPARTMP=         201\n  NA-NS=          85\n'
        t=$((t0 + 1))
        end=$((t0 + tcrit))
    else
        printf '## header\nTIME a\n' > lagr.7
        t=0
        end=$tcrit
    fi
    while [ "$t" -le "$end" ]; do
        printf ' ADJUST:  TIME    %d.00000E+00  T[Myr]   0.100E+01  Q   0.500E+00  DE   0.000E+00 DELTA   0.000E+00 DETOT   0.000E+00 E  -2.500E-01\n' "$t"
        printf '%d.0 1\n' "$t" >> lagr.7
        : > "conf.3_$t"
        # like the engine: the end time is tested in the adjustment, before a
        # stop request is looked at
        [ "$t" -eq "$end" ] && break
        if [ -f STOP ]; then
            printf '\n         TERMINATION BY MANUAL INTERVENTION\n'
            dump "comm.1_$t.0" "$t"
            printf '  MYDUMP    %d.0000000000000000      %d.0                         203 comm.1_%d.0                \n' "$t" "$t" "$t"
            printf '  W MYDUMP J,II,NPARTMP=         203\n  NA-NS=          85\n'
            printf '\n\n         COMMON SAVED AT TOFF/TIME/TTOT =  0.00000000E+00  TCOMP =  0.0\n'
            exit 0
        fi
        sleep 0.1
        t=$((t + 1))
    done
    printf '\n         END RUN    TIME[Myr] =    1.00  TOFF/TIME/TTOT=      0.00000000      %d.00000000      %d.00000000  CPUTOT =    0.0  ERRTOT = 0.00000D+00  DETOT = 0.00000D+00\n' "$end" "$end"
    dump "comm.1_$end.0" "$end"
    printf '  MYDUMP    %d.0000000000000000      %d.0                         204 comm.1_%d.0                \n' "$end" "$end" "$end"
    printf '  W MYDUMP J,II,NPARTMP=         204\n  NA-NS=          85\n'
    """
    base = mktempdir()
    engine_path = joinpath(base, "engine", "build", "nbody6++")
    mkpath(dirname(engine_path))
    # `#!/bin/bash` at column 1, whether or not the literal was dedented.
    write(engine_path, replace(fake_engine, r"^ {4}"m => ""))
    chmod(engine_path, 0o755)
    # TCRIT = 40 is 4 s of steps: a wall budget of 3 s with a stop margin of
    # 1 s stops the first segment part-way.
    write(
        joinpath(base, "in.inp"),
        """
        &INNBODY6
        KSTART=1,TCOMP=1E+08,TCRTP0=3600 /

        &ININPUT
        N=100,NFIX=1,NCRIT=10,NRAND=1,NNBOPT=10,NRUN=1,NCOMM=1,
        DTADJ=1,DELTAT=1,TCRIT=40.0,QE=1.0E-02
        Level='C' /
        """,
    )
    function cluster_cfg(name; chain = "resubmit", max_segments = 0)
        path = joinpath(base, name)
        write(
            path,
            """
            [install]
            enabled = false
            install_dir = "engine"

            [build]
            enable_mpi = false
            enable_gpu = false
            enable_hdf5 = false

            [simulation]
            input_file = "in.inp"
            runs_dir = "runs"
            binary_name = "nbody6++"
            omp_threads = 1
            monitor = false
            telemetry_interval = 0.0
            startup_timeout = 0.0
            exit_grace = 0.0
            wall_budget = 3.0
            stop_margin = 1.0

            [postprocess]
            enabled = false

            [visualization]
            enabled = false

            [cluster]
            enabled = true
            time_limit = "0:10:00"
            overhead = 0.0
            chain = "$chain"
            max_segments = $max_segments
            """,
        )
        return load_config(path)
    end
    # Stand-in for sbatch: appends its arguments to sbatch.log, counts the
    # submissions and prints the count as the job id (`--parsable` form).
    function stand_in_sbatch()
        dir = mktempdir()
        path = joinpath(dir, "sbatch")
        write(
            path,
            """
            #!/bin/bash
            dir=\$(dirname "\$0")
            printf '%s\\n' "\$*" >> "\$dir/sbatch.log"
            n=\$(cat "\$dir/count" 2>/dev/null || echo 0)
            n=\$((n + 1))
            printf '%d\\n' "\$n" > "\$dir/count"
            printf '%d\\n' "\$n"
            """,
        )
        chmod(path, 0o755)
        return dir
    end
    on_path(dir) = "PATH" => dir * ":" * get(ENV, "PATH", "")
    submissions(dir) = isfile(joinpath(dir, "count")) ? read(joinpath(dir, "count"), String) : ""
    ledger(path, run_id) = filter(r -> r.run_id == run_id, Nbody6Dynamics.read_job_ledger(path))
    run_info(run_dir) = Nbody6Dynamics.TOML.parsefile(joinpath(run_dir, "RUN_INFO.toml"))
    runs_ledger = joinpath(base, "runs", "cluster_jobs.csv")

    @testset "Config" begin
        cl = ClusterConfig()
        @test cl.enabled == false
        @test cl.time_limit == "36:00:00"
        @test cl.chain == "resubmit"
        @test cl.overhead == 1800.0
    end

    @testset "Script rendering" begin
        full = (;
            enabled = true,
            partition = "gpu",
            account = "proj42",
            qos = "long",
            job_name = "merger25k",
            time_limit = "1-12:00:00",
            cpus_per_task = 16,
            mem = "64G",
            gres = "gpu:1",
            exclusive = true,
            sbatch_options = ["--hint=nomultithread", "--mail-type=END"],
            prologue = ["module load cuda/13.1", "export JULIA_DEPOT_PATH=/scratch/depot:"],
            julia = "/opt/julia/bin/julia",
        )
        sim = SimulationConfig(; wall_budget = 120_000.0, stop_margin = 600.0)
        with_cluster(cl; simulation = sim) = Nbody6Config(
            InstallConfig(),
            BuildConfig(),
            simulation,
            PostprocessConfig(),
            VisualizationConfig(),
            MergerPipelineConfig(),
            cl,
            CostConfig(),
            base,
        )
        run_dir = joinpath(base, "runs", "render")
        text = job_script(with_cluster(ClusterConfig(; full...)), run_dir)
        lines = split(text, '\n')
        @test lines[1] == "#!/bin/bash"
        for expected in (
            "#SBATCH --job-name=merger25k",
            "#SBATCH --partition=gpu",
            "#SBATCH --account=proj42",
            "#SBATCH --qos=long",
            "#SBATCH --mem=64G",
            "#SBATCH --gres=gpu:1",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=1",
            "#SBATCH --cpus-per-task=16",
            "#SBATCH --exclusive",
            "#SBATCH --time=1-12:00:00",
            "#SBATCH --signal=B:TERM@600",
            "#SBATCH --output=$(joinpath(run_dir, "cluster", "slurm_%j.out"))",
            "#SBATCH --hint=nomultithread",
            "#SBATCH --mail-type=END",
            "export JULIA_NUM_THREADS=\"\${SLURM_CPUS_PER_TASK:-16}\"",
        )
            @test expected in lines
        end
        # Every #SBATCH line precedes the first command, or sbatch ignores it.
        first_command = findfirst(==("set -o pipefail"), lines)
        @test findlast(startswith("#SBATCH"), lines) < first_command
        i_module = findfirst(==("module load cuda/13.1"), lines)
        i_depot = findfirst(==("export JULIA_DEPOT_PATH=/scratch/depot:"), lines)
        i_exec = findfirst(startswith("exec "), lines)
        @test first_command < i_module < i_depot < i_exec
        segment_script = joinpath(Nbody6Dynamics._PACKAGE_ROOT, "scripts", "run_cluster_segment.jl")
        @test lines[i_exec] ==
              "exec '/opt/julia/bin/julia' --startup-file=no '$(segment_script)' '$(run_dir)'"
        @test isfile(segment_script)
        script_path = joinpath(base, "render.sbatch")
        write(script_path, text)
        @test success(`sh -n $script_path`)

        # Explicit Julia threads; a stop margin below 1 s is clamped
        threaded = job_script(
            with_cluster(
                ClusterConfig(; full..., julia_threads = 4);
                simulation = SimulationConfig(; wall_budget = 120_000.0, stop_margin = 0.0),
            ),
            run_dir,
        )
        @test occursin("\nexport JULIA_NUM_THREADS=4\n", threaded)
        @test occursin("\n#SBATCH --signal=B:TERM@1\n", threaded)

        # Omitted options are absent; the job name defaults to the run id
        minimal = job_script(with_cluster(ClusterConfig(; enabled = true)), run_dir)
        for absent in ("--partition", "--account", "--qos", "--mem", "--gres", "--exclusive")
            @test !occursin("#SBATCH $absent", minimal)
        end
        @test occursin("\n#SBATCH --job-name=render\n", minimal)
        @test occursin("\nexport JULIA_NUM_THREADS=\"\${SLURM_CPUS_PER_TASK:-8}\"\n", minimal)
        @test occursin("exec $(Nbody6Dynamics._sh_quote(joinpath(Sys.BINDIR, "julia"))) ", minimal)
    end

    @testset "Dry run, dependency chain of three" begin
        cfg = cluster_cfg("chain.toml"; chain = "dependency", max_segments = 3)
        sub = submit_run(cfg; base_dir = base, run_id = "chain", dry_run = true)
        @test isdir(sub.run_dir) && basename(sub.run_dir) == "chain"
        frozen = load_config(joinpath(sub.run_dir, "cluster", "config.toml"))
        @test isabspath(frozen.simulation.runs_dir)
        @test isabspath(frozen.simulation.input_file)
        @test isabspath(frozen.install.install_dir)
        @test frozen.cluster.chain == "dependency" && frozen.cluster.max_segments == 3
        @test sub.script == joinpath(sub.run_dir, "cluster", "job.sbatch")
        @test isfile(sub.script) && filemode(sub.script) & 0o111 != 0
        @test length(sub.commands) == 3 && isempty(sub.job_ids)
        @test !occursin("--dependency", sub.commands[1])
        @test occursin("--dependency=afterok:<job 1>", sub.commands[2])
        @test occursin("--dependency=afterok:<job 2>", sub.commands[3])
        @test occursin("--kill-on-invalid-dep=yes", sub.commands[2])
        @test !isfile(joinpath(sub.run_dir, "cluster", "jobs.csv"))
        @test !isfile(runs_ledger)

        # A run in progress is not submitted again
        write(joinpath(sub.run_dir, "RUN_INFO.toml"), "[run]\nid = \"chain\"\n")
        @test_throws ArgumentError submit_run(
            cfg;
            base_dir = base,
            run_id = "chain",
            dry_run = true,
        )
        # Nor is a configuration without [cluster] enabled
        local_cfg = Nbody6Config(
            InstallConfig(),
            BuildConfig(),
            SimulationConfig(),
            PostprocessConfig(),
            VisualizationConfig(),
            MergerPipelineConfig(),
        )
        @test_throws ArgumentError submit_run(local_cfg; base_dir = base, dry_run = true)
    end

    @testset "Stand-in sbatch" begin
        sb = stand_in_sbatch()
        cfg = cluster_cfg("dep.toml"; chain = "dependency", max_segments = 3)
        sub = withenv(on_path(sb)) do
            submit_run(cfg; base_dir = base, run_id = "dep")
        end
        @test sub.job_ids == ["1", "2", "3"]
        @test length(sub.commands) == 3
        log = read(joinpath(sb, "sbatch.log"), String)
        @test occursin("--dependency=afterok:1", log)
        @test occursin("--dependency=afterok:2", log)
        @test count("--kill-on-invalid-dep=yes", log) == 2
        for path in (joinpath(sub.run_dir, "cluster", "jobs.csv"), runs_ledger)
            rows = ledger(path, "dep")
            @test length(rows) == 3
            @test all(r -> r.event == "submitted", rows)
            @test [r.job_id for r in rows] == ["1", "2", "3"]
            @test [r.depends_on for r in rows] == ["", "1", "2"]
            @test all(r -> isempty(r.segment), rows)
            @test rows[2].note == sub.commands[2]
        end
    end

    @testset "Segments through cluster_segment" begin
        sb = stand_in_sbatch()
        cfg = cluster_cfg("resubmit.toml")
        sub = submit_run(cfg; base_dir = base, run_id = "segs", dry_run = true)
        run_dir = sub.run_dir
        run_ledger = joinpath(run_dir, "cluster", "jobs.csv")

        # First job: the pipeline starts, stops at the wall budget, and the
        # job submits its successor.
        first_result = withenv(on_path(sb), "SLURM_JOB_ID" => "77") do
            cluster_segment(run_dir)
        end
        @test isempty(first_result)
        info = run_info(run_dir)
        @test length(info["segments"]) == 1
        @test info["segments"][1]["status"] == "stopped"
        @test info["run"]["status"] == "stopped"
        rows = ledger(run_ledger, "segs")
        @test [r.event for r in rows] == ["started", "ended", "submitted"]
        @test rows[1].job_id == "77" && rows[1].segment == "0" && rows[1].host == gethostname()
        @test rows[2].job_id == "77" && rows[2].status == "stopped" && rows[2].segment == "1"
        @test rows[3].job_id == "1" && rows[3].depends_on == "77"
        @test startswith(rows[3].note, "resubmitted by job 77: ")
        @test occursin("--dependency=afterany:77", rows[3].note)
        @test occursin("--dependency=afterany:77", read(joinpath(sb, "sbatch.log"), String))
        @test ledger(runs_ledger, "segs") == rows

        # Changes between jobs go to <run_dir>/config.toml: a budget the rest
        # of the run fits into.
        run_config = joinpath(run_dir, "config.toml")
        c = load_config(run_config)
        sim = SimulationConfig(;
            (k => getfield(c.simulation, k) for k in fieldnames(SimulationConfig))...,
            wall_budget = 60.0,
        )
        save_config(
            Nbody6Config(
                c.install,
                c.build,
                sim,
                c.postprocess,
                c.visualization,
                c.merger,
                c.cluster,
                c.cost,
                c.config_dir,
            ),
            run_config,
        )

        # Second job: the run reaches its end time and the pipeline finishes.
        withenv(on_path(sb), "SLURM_JOB_ID" => "78") do
            cluster_segment(run_dir)
        end
        info = run_info(run_dir)
        @test info["run"]["status"] == "completed"
        @test length(info["segments"]) == 2
        @test Nbody6Dynamics._pipeline_completed(run_dir)
        rows = ledger(run_ledger, "segs")
        @test rows[end].event == "ended" && rows[end].job_id == "78"
        @test rows[end].status == "completed" && rows[end].segment == "2"
        @test count(r -> r.event == "submitted", rows) == 1

        # Third job: nothing left to do, nothing submitted.
        third_result = withenv(on_path(sb), "SLURM_JOB_ID" => "79") do
            cluster_segment(run_dir)
        end
        @test third_result isa Dict{Symbol,Any} && isempty(third_result)
        rows = ledger(run_ledger, "segs")
        @test rows[end].event == "ended" && rows[end].job_id == "79"
        @test rows[end].status == "completed"
        @test count(r -> r.event == "submitted", rows) == 1
        @test submissions(sb) == "1\n"

        # max_segments = 1: the first stopped segment ends the chain.
        single_cfg = cluster_cfg("single.toml"; max_segments = 1)
        single = submit_run(single_cfg; base_dir = base, run_id = "single", dry_run = true)
        withenv(on_path(sb), "SLURM_JOB_ID" => "90") do
            cluster_segment(single.run_dir)
        end
        @test Nbody6Dynamics._run_status(single.run_dir) == "stopped"
        single_rows = ledger(joinpath(single.run_dir, "cluster", "jobs.csv"), "single")
        @test [r.event for r in single_rows] == ["started", "ended"]
        @test submissions(sb) == "1\n"
    end

    @testset "read_job_ledger" begin
        dir = mktempdir()
        @test isempty(Nbody6Dynamics.read_job_ledger(joinpath(dir, "absent.csv")))
        @test Nbody6Dynamics._csv_field("plain") == "plain"
        @test Nbody6Dynamics._csv_field("a,b") == "\"a,b\""
        @test Nbody6Dynamics._csv_field("say \"x\"") == "\"say \"\"x\"\"\""
        note = "sbatch --comment=\"a, b\" job.sbatch\nsecond line"
        Nbody6Dynamics._record_job(
            joinpath(dir, "r1"),
            dir;
            event = "submitted",
            job_id = "5",
            note = note,
        )
        Nbody6Dynamics._record_job(joinpath(dir, "r1"), dir; event = "started", job_id = "5")
        rows = Nbody6Dynamics.read_job_ledger(joinpath(dir, "cluster_jobs.csv"))
        @test length(rows) == 2
        @test fieldnames(typeof(rows[1])) ==
              (:time, :run_id, :event, :job_id, :depends_on, :segment, :status, :host, :note)
        @test all(r -> all(v -> v isa String, values(r)), rows)
        @test rows[1].note == note
        @test rows[1].run_id == "r1" && rows[1].job_id == "5" && rows[1].depends_on == ""
        @test rows[2].event == "started" && rows[2].note == ""
        @test occursin(r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}$", rows[1].time)
        @test Nbody6Dynamics.read_job_ledger(joinpath(dir, "r1", "cluster", "jobs.csv")) == rows
        bad = joinpath(dir, "bad.csv")
        write(bad, "time,run,event\n1,2,3\n")
        @test_throws ArgumentError Nbody6Dynamics.read_job_ledger(bad)
    end
end
