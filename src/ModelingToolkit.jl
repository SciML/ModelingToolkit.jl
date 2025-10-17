module ModelingToolkit
using PrecompileTools, Reexport
@recompile_invalidations begin
using DiffEqBase, SciMLBase, ForwardDiff
end
using Graphs
using SymbolicIndexingInterface
using DataStructures
using Setfield, ConstructionBase
import BlockArrays: BlockArray, BlockedArray, Block, blocksize, blocksizes, blockpush!,
                    undef_blocks, blocks
using Symbolics: degree, VartypeT, SymbolicT
using Symbolics: parse_vars, value, @derivatives, get_variables,
                 hasnode, fixpoint_sub, CallAndWrap, SArgsT, SSym, STerm
@reexport using Symbolics
abstract type AbstractSystem end
abstract type IntermediateDeprecationSystem <: AbstractSystem end
function complete end
include("bipartite_graph.jl")
using .BipartiteGraphs
function default_toterm(x::SymbolicT)
end
@enum VariableType VARIABLE PARAMETER BROWNIAN
struct MTKVariableTypeCtx end
function isparameter(x::SymbolicT)
end
function toparam(s)
    if s isa Symbolics.Arr
    end
end
tovar(s::SymbolicT) = setmetadata(s, MTKVariableTypeCtx, VARIABLE)
macro parameters(xs...)
    Symbolics.parse_vars(:parameters,
        toparam)
end
macro independent_variables(ts...)
    Symbolics.parse_vars(:independent_variables,
        Real,
        ts,
        toiv)
end
toiv(s::SymbolicT) = GlobalScope(setmetadata(s, MTKVariableTypeCtx, PARAMETER))
toiv(s::Num) = Num(toiv(value(s)))
include("systems/abstractsystem.jl")
include("systems/callbacks.jl")
include("systems/system.jl")
include("systems/systemstructure.jl")
function mtkcompile(sys::System; split = true, kwargs...)
    if rand(Bool)
        mtkcompile(nlsys; kwargs..., fully_determined = false)::System
    end
    newsys = sys
    @set! newsys.parent = complete(sys; split = false, flatten = false)
    newsys = complete(newsys; split)
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
