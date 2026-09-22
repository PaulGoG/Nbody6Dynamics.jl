# Included by runtests.jl inside its top-level test set; the helpers and
# constants of that file are in scope.

@testset "Static QA — Aqua" begin
    using Aqua
    # Method ambiguities are checked for this package only — recursing
    # into the Makie/SciML dependency tree reports upstream noise.
    Aqua.test_all(Nbody6Dynamics; ambiguities = false, persistent_tasks = false)
    # Aqua's ambiguity check spawns a process that `require`s each module
    # it is given; an extension is not loadable that way, so the package
    # is checked here and the pair in the figure-extension testset.
    Aqua.test_ambiguities(Nbody6Dynamics)
end

@testset "Static QA — ExplicitImports" begin
    using ExplicitImports
    # Pragmatic subset: the package uses plain `using` for its small,
    # stable dependency surface (full explicit-import migration is
    # tracked in the roadmap). These checks catch the real hazards:
    # stale explicit imports, self-qualified names, and accesses of
    # non-owning modules. The extension imports the core internals it
    # needs explicitly, so its list is held to the same standard.
    @test check_no_stale_explicit_imports(Nbody6Dynamics) === nothing
    @test check_no_self_qualified_accesses(Nbody6Dynamics) === nothing
    makie_ext = Base.get_extension(Nbody6Dynamics, :Nbody6DynamicsMakieExt)
    @test check_no_stale_explicit_imports(makie_ext) === nothing
    @test check_no_self_qualified_accesses(makie_ext) === nothing
end

@testset "Static QA — JET" begin
    using JET
    # Reports scoped to this package's own frames — Base/dependency
    # internals (e.g. @sync's sync_end, tuple broadcasting) produce
    # known false positives outside our control.
    JET.test_package(Nbody6Dynamics; target_modules = (Nbody6Dynamics,))
    # The figure layer is not covered here: `report_package` takes a
    # package and an extension is not one, and `report_file` on the
    # extension analyses it against the package's own project, where the
    # Makie trigger packages are weak dependencies and cannot be
    # resolved. Aqua's ambiguity check above and the figure smoke tests
    # are what the extension has instead.
end
