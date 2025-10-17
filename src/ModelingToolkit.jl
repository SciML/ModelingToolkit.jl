""""""
module ModelingToolkit
using PrecompileTools, Reexport
@recompile_invalidations begin
    using StaticArrays
    using JumpProcesses
    import REPL
end
import SymbolicUtils
import SymbolicUtils as SU
import SymbolicUtils: iscall, arguments, operation, maketerm, promote_symtype,
                      isadd, ismul, ispow, issym, FnType, isconst, BSImpl,
                      @rule, Rewriters, substitute, metadata, BasicSymbolic
using SpecialFunctions, NaNMath
@recompile_invalidations begin
using DiffEqBase, SciMLBase, ForwardDiff
end
using Graphs
import ExprTools: splitdef, combinedef
import OrderedCollections
using SymbolicIndexingInterface
using LinearAlgebra, SparseArrays
using InteractiveUtils
using DataStructures
@static if pkgversion(DataStructures) >= v"0.19"
end
using Base.Threads
using Latexify, ArrayInterface
using Setfield, ConstructionBase
using SciMLBase: StandardODEProblem, StandardNonlinearProblem, handle_varmap, TimeDomain,
                 PeriodicClock, Clock, SolverStepClock, ContinuousClock, OverrideInit,
                 NoInit
using Distributed
using MLStyle
import Moshi
using Moshi.Data: @data
import SCCNonlinearSolve
import Graphs: SimpleDiGraph, add_edge!, incidence_matrix
import BlockArrays: BlockArray, BlockedArray, Block, blocksize, blocksizes, blockpush!,
                    undef_blocks, blocks
using OffsetArrays: Origin
using RuntimeGeneratedFunctions
using RuntimeGeneratedFunctions: drop_expr
using Symbolics: degree, VartypeT, SymbolicT
using Symbolics: parse_vars, value, @derivatives, get_variables,
                 exprs_occur_in, symbolic_linear_solve, unwrap, wrap,
                 VariableSource, getname, variable,
                 NAMESPACE_SEPARATOR, setdefaultval, Arr,
                 hasnode, fixpoint_sub, CallAndWrap, SArgsT, SSym, STerm
const NAMESPACE_SEPARATOR_SYMBOL = Symbol(NAMESPACE_SEPARATOR)
import Symbolics: rename, get_variables!, _solve, hessian_sparsity,
                  jacobian_sparsity, isaffine, islinear, _iszero, _isone,
                  scalarize, hasderiv
import DiffEqBase: @add_kwonly
export independent_variables, unknowns, observables, parameters, full_parameters,
       continuous_events, discrete_events
@reexport using Symbolics
@reexport using UnPack
RuntimeGeneratedFunctions.init(@__MODULE__)
import DynamicQuantities
const DQ = DynamicQuantities
for fun in [:toexpr]
    @eval begin
        function $fun(eq::Equation; kw...)
            Expr(:call, :(==), $fun(eq.lhs; kw...), $fun(eq.rhs; kw...))
            if ineq.relational_op == Symbolics.leq
                Expr(:call, :(<=), $fun(ineq.lhs; kw...), $fun(ineq.rhs; kw...))
            end
        end
    end
end
abstract type AbstractSystem end
abstract type IntermediateDeprecationSystem <: AbstractSystem end
function complete end
include("bipartite_graph.jl")
using .BipartiteGraphs
include("variables.jl")
include("parameters.jl")
include("independent_variables.jl")
const SymmapT = Dict{SymbolicT, SymbolicT}
include("utils.jl")
include("systems/index_cache.jl")
include("systems/abstractsystem.jl")
include("systems/connectors.jl")
include("systems/callbacks.jl")
include("systems/system.jl")
include("systems/sparsematrixclil.jl")
include("clock.jl")
include("systems/systemstructure.jl")
include("systems/systems.jl")
include("systems/alias_elimination.jl")
function linear_subsys_adjmat!(state::TransformationState; kwargs...)
    graph = state.structure.graph
    if state.structure.solvable_graph === nothing
        state.structure.solvable_graph = BipartiteGraph(nsrcs(graph), ndsts(graph))
    end
    linear_equations = Int[]
    eadj = Vector{Int}[]
    cadj = Vector{Int}[]
    mm = SparseMatrixCLIL(nsrcs(graph),
        ndsts(graph),
        linear_equations, eadj, cadj)
end
const t_nounits = let
    only(@independent_variables t)
end
const D_nounits = Differential(t_nounits)
export System, Pre, complete
PrecompileTools.@compile_workload begin
    @variables x(ModelingToolkit.t_nounits) y(ModelingToolkit.t_nounits)
    sys = System([ModelingToolkit.D_nounits(x) ~ x * y, y ~ 3x + 4 * ModelingToolkit.D_nounits(y)], ModelingToolkit.t_nounits, [x, y], Num[]; name = :sys)
    TearingState(sys)
end
end
