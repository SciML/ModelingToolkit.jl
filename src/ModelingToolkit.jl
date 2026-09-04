"""
$(DocStringExtensions.README)
"""
module ModelingToolkit

import PrecompileTools
using PrecompileTools: @recompile_invalidations
@recompile_invalidations begin
    import StaticArrays
    import Symbolics
    # ONLY here for the invalidations
    import REPL
    import BlockArrays: BlockVector, undef_blocks
end

import SymbolicUtils
import SymbolicUtils as SU
import SymbolicUtils: iscall, arguments, operation,
    issym, BSImpl, Operator,
    @rule, Rewriters, substitute, BasicSymbolic,
    symtype, _iszero, unwrap
import TermInterface: maketerm, metadata
import SymbolicUtils.Code
import SymbolicUtils.Code: Assignment, AtIndex, Let, MakeArray, SetArray, toexpr
import DocStringExtensions
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
@recompile_invalidations begin
    import DiffEqBase, SciMLBase, ForwardDiff
end
import Graphs
using Graphs: SimpleGraph, bfs_parents, rem_edge!, topological_sort_by_dfs
import OrderedCollections
using OrderedCollections: OrderedSet

import SymbolicIndexingInterface
using SymbolicIndexingInterface: NotSymbolic, ParameterIndexingProxy, ProblemState,
    current_time, getp, getsym, getu, hasname, is_parameter, is_time_dependent,
    parameter_index, parameter_values, setp, setp_oop, setsym, setu, state_values,
    symbolic_type, variable_symbols
using SymbolicIndexingInterface: getname
import LinearAlgebra, SparseArrays
using LinearAlgebra: I, UniformScaling, UpperTriangular, cond, issuccess, lu, mul!, svd
using SparseArrays: SparseMatrixCSC, blockdiag, findnz, nnz, sparse, spzeros
import InteractiveUtils
import DataStructures
using DataStructures: Queue
import Base.Threads
import Setfield, ConstructionBase
using Setfield: @set!, @set
import Libdl
import Combinatorics
using SciMLBase: TimeDomain, Clock, SolverStepClock, ContinuousClock, OverrideInit, NoInit,
    LinearProblem, NonlinearFunction, NonlinearLeastSquaresProblem, NonlinearProblem,
    ODEFunction, ODEProblem, SCCNonlinearProblem, SplitFunction, SplitODEProblem, remake
import Moshi
import SCCNonlinearSolve
import Graphs: add_edge!
import CommonSolve
using CommonSolve: solve

import RuntimeGeneratedFunctions

using Symbolics: VartypeT, SymbolicT
using Symbolics: value, @derivatives, get_variables, symbolic_linear_solve, wrap,
    COMMON_ZERO, NAMESPACE_SEPARATOR, fixpoint_sub
const NAMESPACE_SEPARATOR_SYMBOL = Symbol(NAMESPACE_SEPARATOR)
import Symbolics: rename, islinear,
    tosymbol, build_function
const derivative = Symbolics.derivative
import ModelingToolkitBase as MTKBase
import SimpleNonlinearSolve

using UnPack: @unpack
RuntimeGeneratedFunctions.init(@__MODULE__)

import DifferentiationInterface as DI
using ADTypes: AutoForwardDiff
import SciMLPublic: @public
import PreallocationTools
import PreallocationTools: DiffCache
import FillArrays
import BipartiteGraphs
using BipartiteGraphs: BipartiteEdge, BipartiteGraph, DiCMOBiGraph,
    MatchedCondensationGraph, Unassigned, invview, maximal_matching, ndsts, nsrcs,
    𝑑neighbors, 𝑑vertices, 𝑠neighbors
import SciMLStructures

@recompile_invalidations begin
    import StateSelection
    import StateSelection: CLIL
    import ModelingToolkitTearing as MTKTearing
    using ModelingToolkitTearing: TearingState, SystemStructure

    MTKBase.complete(dg::StateSelection.DiffGraph) = BipartiteGraphs.complete(dg)
end

macro import_mtkbase()
    allnames = names(MTKBase; all = true)
    banned_names = Set{Symbol}([:eval, :include, :Variable, :__init__])
    using_expr = Expr(:using, Expr(:(:), Expr(:., :ModelingToolkitBase)))
    inner_using_expr = using_expr.args[1]

    public_expr = :(@public)
    inner_public_expr = Expr(:tuple)
    push!(public_expr.args, inner_public_expr)
    export_expr = Expr(:export)

    for name in allnames
        name in banned_names && continue
        startswith(string(name), '#') && continue
        push!(inner_using_expr.args, Expr(:., name))
        # Base.ispublic was added in Julia 1.11
        @static if VERSION >= v"1.11"
            if Base.ispublic(MTKBase, name) && !Base.isexported(MTKBase, name)
                push!(inner_public_expr.args, name)
            end
        end
        Base.isexported(MTKBase, name) && push!(export_expr.args, name)
    end

    return quote
        $using_expr
        $export_expr
        $(esc(public_expr))
    end
end

@import_mtkbase

using ModelingToolkitBase: COMMON_NOTHING, COMMON_MISSING, COMMON_TRUE
using ModelingToolkitBase: build_function_wrapper, BuildFunctionWrapperOptions,
    GeneratedFunctionOptions

# `@recompile_invalidations` evaluates its body through `Core.eval` at run time, which
# hides these files from static analyzers and resolves the relative paths against the
# caller's source directory rather than `src/`. Keep the includes at real top level.
include("linearization.jl")
include("systems/analysis_points.jl")
include("systems/solver_nlprob.jl")

include("problems/docs.jl")
include("systems/codegen.jl")
include("systems/codegen_compat.jl")
include("problems/semilinearodeproblem.jl")
include("problems/sccnonlinearproblem.jl")

include("discretedomain.jl")
include("systems/systemstructure.jl")
include("initialization.jl")
include("systems/systems.jl")
include("systems/clock_inference.jl")
include("systems/if_lifting.jl")
include("systems/substitute_component.jl")

include("systems/alias_elimination.jl")
include("structural_transformation/StructuralTransformations.jl")

import .StructuralTransformations
using .StructuralTransformations: tearing, dae_index_lowering, dummy_derivative,
    sorted_incidence_matrix, pantelides_reassemble, find_solvables!,
    tearing_substitution, but_ordered_incidence, lowest_order_variable_mask,
    highest_order_variable_mask
export tearing, dae_index_lowering, dummy_derivative,
    sorted_incidence_matrix, pantelides_reassemble, find_solvables!,
    tearing_substitution, but_ordered_incidence, lowest_order_variable_mask,
    highest_order_variable_mask
export StructuralTransformations

export SemilinearODEFunction, SemilinearODEProblem
export analyze_initialization_jacobian
export alias_elimination
export linearize, linearization_function,
    LinearizationProblem, LinearizationOpPoint, linearization_ap_transform
export solve
export map_variables_to_equations, substitute_component

export TearingState

export Clock, SolverStepClock, TimeDomain
export get_sensitivity_function, get_comp_sensitivity_function,
    get_looptransfer_function, get_sensitivity, get_comp_sensitivity, get_looptransfer
export isolate_subsystem

"""
    FMIComponent(version; fmu, type, name, kwargs...)

Construct a ModelingToolkit component from a Functional Mock-up Unit (FMU).

# Arguments

- `version`: a `Val` selecting the FMI major version, for example `Val(2)`.

# Keyword Arguments

- `fmu`: the loaded FMU object.
- `type`: the FMI interface type, such as `:ME` for model exchange.
- `name`: the component name used by `@named`.
- `kwargs...`: additional FMI extension options.

# Returns

An `AbstractSystem` component representing the FMU interface.

# Examples

```julia
using ModelingToolkit

@named model = ModelingToolkit.FMIComponent(Val(2); fmu, type = :ME)
```

!!! note
    `FMIComponent` is implemented by the FMI extension. Load the FMI dependency before
    constructing FMU-backed components.
"""
function FMIComponent end

@public linearize_symbolic, reorder_unknowns
@public similarity_transform
@public precompile_ode_problem, precompile_dae_problem

include("precompile.jl")

function __init__()
    SU.hashcons(StructuralTransformations.NOTHING_EQ.lhs, true)
    SU.hashcons(StructuralTransformations.NOTHING_EQ.rhs, true)
    SU.hashcons(unwrap(ODE_GAMMA[1]), true)
    SU.hashcons(unwrap(ODE_GAMMA[2]), true)
    SU.hashcons(unwrap(ODE_GAMMA[3]), true)
    SU.hashcons(unwrap(ODE_C), true)
    return SU.hashcons(SCC_EXPLICITFUN_CACHE_OUT, true)
end

end # module
