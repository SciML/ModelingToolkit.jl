using SciMLTesting, ModelingToolkitBase, Test

run_qa(
    ModelingToolkitBase;
    # The bespoke JET type-stability suite lives in `jet.jl` (run separately in the same
    # session). Because `using JET` there registers JET process-wide, `run_qa` would
    # otherwise default `jet = true` and run a second, hard `JET.report_package` typo
    # check on top -- one that never ran under the previous `Aqua.test_all` QA. Keep the
    # JET coverage to `jet.jl` alone; `run_qa` here is Aqua + ExplicitImports only.
    jet = false,
    # Pre-existing, tracked Aqua findings (https://github.com/SciML/ModelingToolkit.jl/issues/4670):
    #  * ambiguities / unbound_args: method ambiguities and unbound type parameters
    #    (e.g. `_remake_buffer` with an unbound `P`).
    #  * undefined_exports: `ModelingToolkitBase.Variable` is reexported but not defined.
    #  * stale_deps: `SimpleNonlinearSolve` is declared but unused.
    #  * deps_compat: the `Optimization` test extra lacks a `[compat]` entry.
    #  * piracies: ModelingToolkitBase deliberately defines `search_variables!`/`toexpr`
    #    methods on upstream types.
    # Marked broken so a clean lane records `Broken` and the placeholder prompts a fix.
    aqua_broken = (
        :ambiguities, :unbound_args, :undefined_exports,
        :stale_deps, :deps_compat, :piracies,
    ),
    explicit_imports = true,
    ei_kwargs = (;
        # Names imported/accessed through a reexport hub rather than their defining
        # owner (e.g. `unwrap`/`Term`/`scalarize` reexported by `Symbolics` from
        # `SymbolicUtils`, `NullParameters`/`@add_kwonly` reexported by `DiffEqBase`
        # from `SciMLBase`, `endpoints` reexported by `DomainSets` from `IntervalSets`).
        all_explicit_imports_via_owners = (;
            ignore = (
                Symbol("@add_kwonly"), :Operator, :Term, :_isone, :_iszero,
                :getname, :maketerm, :metadata, :scalarize, :unwrap,
            ),
        ),
        all_qualified_accesses_via_owners = (;
            ignore = (
                :AbstractParameterizedFunction, :DISCRETE_INPLACE_DEFAULT, :FnType,
                :NullParameters, :Operator, :_iszero, :allowedkeywords, :anyeltypedual,
                :endpoints, :get_updated_symbolic_problem, :metadata, :promote_u0,
                :scalarize, :search_variables!, :shape, :unwrap, :wrapfun_iip,
            ),
        ),
    ),
    # Pre-existing, tracked findings (https://github.com/SciML/ModelingToolkit.jl/issues/4670):
    #  * no_implicit_imports / no_stale_explicit_imports: ModelingToolkitBase relies on a
    #    large block of implicit imports (whole packages such as `BipartiteGraphs`,
    #    `Graphs`, `DataStructures`, `DiffEqBase`, ...) and carries an equally large block
    #    of stale explicit imports.
    #  * the public-API checks are dominated by names that are not (yet) declared public
    #    in their defining packages (Symbolics, SymbolicUtils, SciMLBase, ...); the fix is
    #    upstream `public` declarations (the SciML make-public effort).
    # Marked broken so a clean lane records `Broken` and auto-flags an Unexpected Pass
    # once each underlying finding is resolved.
    ei_broken = (
        :no_implicit_imports, :no_stale_explicit_imports,
        :all_qualified_accesses_are_public, :all_explicit_imports_are_public,
    ),
)
