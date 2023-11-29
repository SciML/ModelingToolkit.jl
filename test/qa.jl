using ModelingToolkit, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(ModelingToolkit)
    Aqua.test_ambiguities(ModelingToolkit, recursive = false)
    Aqua.test_deps_compat(ModelingToolkit)
    Aqua.test_piracies(ModelingToolkit,
        treat_as_own = [])
    Aqua.test_project_extras(ModelingToolkit)
    Aqua.test_stale_deps(ModelingToolkit)
    Aqua.test_unbound_args(ModelingToolkit)
    Aqua.test_undefined_exports(ModelingToolkit)
end
