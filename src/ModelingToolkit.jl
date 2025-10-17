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
abstract type AbstractSystem end
abstract type IntermediateDeprecationSystem <: AbstractSystem end
function complete end
include("bipartite_graph.jl")
using .BipartiteGraphs
function default_toterm(x::SymbolicT)
end
import SymbolicUtils: symtype, term, hasmetadata, issym
@enum VariableType VARIABLE PARAMETER BROWNIAN
struct MTKVariableTypeCtx end
getvariabletype(x, def = VARIABLE) = safe_getmetadata(MTKVariableTypeCtx, unwrap(x), def)::Union{typeof(def), VariableType}
isparameter(x::Union{Num, Symbolics.Arr, Symbolics.CallAndWrap}) = isparameter(unwrap(x))
function isparameter(x::SymbolicT)
    varT = getvariabletype(x, nothing)
    return varT === PARAMETER
end
isparameter(x) = false
function toparam(s)
    if s isa Symbolics.Arr
        Symbolics.wrap(toparam(Symbolics.unwrap(s)))
    elseif s isa AbstractArray
        map(toparam, s)
    else
        setmetadata(s, MTKVariableTypeCtx, PARAMETER)
    end
end
toparam(s::Num) = wrap(toparam(value(s)))
tovar(s::SymbolicT) = setmetadata(s, MTKVariableTypeCtx, VARIABLE)
tovar(s::Union{Num, Symbolics.Arr}) = wrap(tovar(unwrap(s)))
macro parameters(xs...)
    Symbolics.parse_vars(:parameters,
        Real,
        xs,
        toparam)
end
macro independent_variables(ts...)
    Symbolics.parse_vars(:independent_variables,
        Real,
        ts,
        toiv)
end
toiv(s::SymbolicT) = GlobalScope(setmetadata(s, MTKVariableTypeCtx, PARAMETER))
toiv(s::Symbolics.Arr) = wrap(toiv(value(s)))
toiv(s::Num) = Num(toiv(value(s)))
const SymmapT = Dict{SymbolicT, SymbolicT}
include("systems/abstractsystem.jl")
include("systems/callbacks.jl")
include("systems/system.jl")
include("systems/systemstructure.jl")
include("systems/systems.jl")
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
