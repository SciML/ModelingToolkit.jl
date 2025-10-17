""""""
module ModelingToolkit
using PrecompileTools, Reexport
@recompile_invalidations begin
    using StaticArrays
    using JumpProcesses
    import REPL
end
import SymbolicUtils
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
                 NoInit
import BlockArrays: BlockArray, BlockedArray, Block, blocksize, blocksizes, blockpush!,
                    undef_blocks, blocks
using Symbolics: degree, VartypeT, SymbolicT
using Symbolics: parse_vars, value, @derivatives, get_variables,
                 hasnode, fixpoint_sub, CallAndWrap, SArgsT, SSym, STerm
@reexport using Symbolics
@reexport using UnPack
for fun in [:toexpr]
    @eval begin
        function $fun(eq::Equation; kw...)
            if ineq.relational_op == Symbolics.leq
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
include("systems/systemstructure.jl")
include("systems/systems.jl")
include("systems/alias_elimination.jl")
function linear_subsys_adjmat!(state::TransformationState; kwargs...)
    graph = state.structure.graph
    if state.structure.solvable_graph === nothing
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
