using ModelingToolkitStandardLibrary, Aqua
@testset "Aqua" begin
    Aqua.find_persistent_tasks_deps(ModelingToolkitStandardLibrary)
    Aqua.test_ambiguities(ModelingToolkitStandardLibrary, recursive = false)
    Aqua.test_deps_compat(ModelingToolkitStandardLibrary)
    Aqua.test_piracies(ModelingToolkitStandardLibrary)
    Aqua.test_project_extras(ModelingToolkitStandardLibrary)
    Aqua.test_stale_deps(ModelingToolkitStandardLibrary)
    Aqua.test_unbound_args(ModelingToolkitStandardLibrary)
    Aqua.test_undefined_exports(ModelingToolkitStandardLibrary)
end
