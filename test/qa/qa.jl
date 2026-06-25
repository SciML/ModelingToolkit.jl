using SciMLTesting, ModelingToolkit, Test
using JET

run_qa(
    ModelingToolkit;
    jet_kwargs = (; target_defined_modules = true),
    # JET's `report_package` re-parses the package toplevel, and the
    # `@recompile_invalidations begin using StaticArrays; using DiffEqBase; ... end`
    # blocks make it `require` ModelingToolkit's own (transitive, not directly listed)
    # dependencies, which surface as JET toplevel errors. Tracked in
    # https://github.com/SciML/ModelingToolkit.jl/issues/4670; marked broken so the lane
    # records `Broken` and auto-flags once JET is clean.
    jet_broken = true,
    # Pre-existing, tracked Aqua findings (https://github.com/SciML/ModelingToolkit.jl/issues/4670):
    #  * undefined_exports: `ModelingToolkit.Variable` is reexported but not defined in
    #    the namespace.
    #  * piracies: ModelingToolkit deliberately defines SciMLBase problem constructors
    #    (e.g. `SCCNonlinearProblem(::System, op)`) and a handful of MTKBase/Symbolics
    #    methods on upstream types -- intentional, long-standing "piracy".
    aqua_broken = (:undefined_exports, :piracies),
    explicit_imports = true,
    ei_kwargs = (;
        # Names imported/accessed through a reexport hub rather than their defining
        # owner (e.g. `unwrap`/`Term` reexported by `Symbolics` from `SymbolicUtils`,
        # `value`/`var_from_nested_derivative` reexported by `ModelingToolkitBase`, the
        # `BipartiteGraphs`/`Graphs` graph helpers reexported across the hub).
        all_explicit_imports_via_owners = (;
            ignore = (
                Symbol("@add_kwonly"), :Operator, :Term, :_isone, :_iszero,
                :getname, :maketerm, :metadata, :scalarize, :unwrap,
                :IncrementalCycleTracker, :add_edge_checked!, :schedule,
                :topological_sort, :value, :var_from_nested_derivative,
            ),
        ),
        all_qualified_accesses_via_owners = (;
            ignore = (:NullParameters, :unwrap, :value),
        ),
    ),
    # All four are pre-existing, tracked findings
    # (https://github.com/SciML/ModelingToolkit.jl/issues/4670):
    #  * no_implicit_imports / no_stale_explicit_imports: the `StructuralTransformations`
    #    submodule relies on a large block of implicit imports (whole packages such as
    #    `BipartiteGraphs`/`Graphs`) and carries an equally large block of stale explicit
    #    imports; the top-level module is additionally unanalyzable (a dynamic `include`
    #    in `src/precompile.jl`).
    #  * the public-API checks are dominated by names that are not (yet) declared public
    #    in their defining packages (Symbolics, SymbolicUtils, SciMLBase, StateSelection,
    #    ModelingToolkitBase, ModelingToolkitTearing, ...); the fix is upstream `public`
    #    declarations (the SciML make-public effort).
    # Marked broken so a clean lane records `Broken` and auto-flags an Unexpected Pass
    # once each underlying finding is resolved.
    ei_broken = (
        :no_implicit_imports, :no_stale_explicit_imports,
        :all_qualified_accesses_are_public, :all_explicit_imports_are_public,
    ),
)
