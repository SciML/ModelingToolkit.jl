"""
$(DocStringExtensions.README)
"""
module ModelingToolkit
using PrecompileTools, Reexport
@recompile_invalidations begin
    using StaticArrays
    using Symbolics
    # ONLY here for the invalidations
    import REPL
end

import SymbolicUtils
import SymbolicUtils as SU
import SymbolicUtils: iscall, arguments, operation, maketerm, promote_symtype,
                      isadd, ismul, ispow, issym, FnType, isconst, BSImpl,
                      @rule, Rewriters, substitute, metadata, BasicSymbolic,
                      symtype
using SymbolicUtils.Code
import SymbolicUtils.Code: toexpr
import SymbolicUtils.Rewriters: Chain, Postwalk, Prewalk, Fixpoint
using DocStringExtensions
@recompile_invalidations begin
    using DiffEqBase, SciMLBase, ForwardDiff
end
using Graphs
import OrderedCollections

using SymbolicIndexingInterface
using LinearAlgebra, SparseArrays
using InteractiveUtils
using DataStructures
@static if pkgversion(DataStructures) >= v"0.19"
    import DataStructures: IntDisjointSet
else
    import DataStructures: IntDisjointSets
    const IntDisjointSet = IntDisjointSets
end
using Base.Threads
using Setfield
import Libdl
using DocStringExtensions
using Base: RefValue
using Combinatorics
using SciMLBase: StandardODEProblem, StandardNonlinearProblem, handle_varmap, TimeDomain,
                 PeriodicClock, Clock, SolverStepClock, ContinuousClock, OverrideInit,
                 NoInit
import Moshi
using Moshi.Data: @data
import SCCNonlinearSolve
using Reexport
import Graphs: SimpleDiGraph, add_edge!, incidence_matrix
using OffsetArrays: Origin
import CommonSolve

using RuntimeGeneratedFunctions
using RuntimeGeneratedFunctions: drop_expr

using Symbolics: degree, VartypeT, SymbolicT
using Symbolics: parse_vars, value, @derivatives, get_variables,
                 exprs_occur_in, symbolic_linear_solve, unwrap, wrap,
                 VariableSource, getname, variable, COMMON_ZERO,
                 NAMESPACE_SEPARATOR, setdefaultval, Arr,
                 hasnode, fixpoint_sub, CallAndWrap, SArgsT, SSym, STerm
const NAMESPACE_SEPARATOR_SYMBOL = Symbol(NAMESPACE_SEPARATOR)
import Symbolics: rename, get_variables!, _solve, hessian_sparsity,
                  jacobian_sparsity, isaffine, islinear, _iszero, _isone,
                  tosymbol, lower_varname, diff2term, var_from_nested_derivative,
                  BuildTargets, JuliaTarget, StanTarget, CTarget, MATLABTarget,
                  ParallelForm, SerialForm, MultithreadedForm, build_function,
                  rhss, lhss, gradient,
                  jacobian, hessian, derivative, sparsejacobian, sparsehessian,
                  scalarize, hasderiv
import ModelingToolkitBase as MTKBase

import DiffEqBase: @add_kwonly
@reexport using Symbolics
@reexport using UnPack
@reexport using ModelingToolkitBase
RuntimeGeneratedFunctions.init(@__MODULE__)

import DifferentiationInterface as DI
using ADTypes: AutoForwardDiff
import SciMLPublic: @public
import PreallocationTools
import PreallocationTools: DiffCache
import FillArrays
using BipartiteGraphs
import BlockArrays: BlockArray, BlockedArray, Block, blocksize, blocksizes, blockpush!,
                    undef_blocks, blocks

@recompile_invalidations begin
    import StateSelection
    import StateSelection: CLIL
    import ModelingToolkitTearing as MTKTearing
    using ModelingToolkitTearing: TearingState, SystemStructure

    ModelingToolkitBase.complete(dg::StateSelection.DiffGraph) = BipartiteGraphs.complete(dg)
end

macro import_mtkbase()
    allnames = names(MTKBase; all = true)
    banned_names = Set{Symbol}([:eval, :include, :Variable])
    using_expr = Expr(:using, Expr(:(:), Expr(:., :ModelingToolkitBase)))
    inner_using_expr = using_expr.args[1]

    public_expr = :(@public)
    inner_public_expr = Expr(:tuple)
    push!(public_expr.args, inner_public_expr)

    for name in allnames
        name in banned_names && continue
        startswith(string(name), '#') && continue
        push!(inner_using_expr.args, Expr(:., name))
        if Base.ispublic(MTKBase, name) && !Base.isexported(MTKBase, name)
            push!(inner_public_expr.args, name)
        end
    end

    quote
        $using_expr
        $(esc(public_expr))
    end
end

@import_mtkbase

using ModelingToolkitBase: COMMON_SENTINEL, COMMON_NOTHING, COMMON_MISSING,
                           COMMON_TRUE, COMMON_FALSE, COMMON_INF

@recompile_invalidations begin
    include("linearization.jl")
    include("systems/analysis_points.jl")
    include("systems/solver_nlprob.jl")

    include("problems/docs.jl")
    include("systems/codegen.jl")
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
end

@reexport using .StructuralTransformations

export SemilinearODEFunction, SemilinearODEProblem
export alias_elimination
export linearize, linearization_function,
       LinearizationProblem, linearization_ap_transform
export solve
export map_variables_to_equations, substitute_component

export TearingState

export Clock, SolverStepClock, TimeDomain
export get_sensitivity_function, get_comp_sensitivity_function,
       get_looptransfer_function, get_sensitivity, get_comp_sensitivity, get_looptransfer

function FMIComponent end

@public linearize_symbolic, reorder_unknowns
@public similarity_transform

include(pkgdir(ModelingToolkitBase, "src", "precompile.jl"))
end # module
