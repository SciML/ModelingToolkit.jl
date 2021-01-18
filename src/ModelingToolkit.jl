module ModelingToolkit

using DiffEqBase, Distributed
using StaticArrays, LinearAlgebra, SparseArrays, LabelledArrays
using Latexify, Unitful, ArrayInterface
using MacroTools
using UnPack: @unpack
using DiffEqJump
using DataStructures: OrderedDict, OrderedSet
using SpecialFunctions, NaNMath
using RuntimeGeneratedFunctions
using Base.Threads
import MacroTools: splitdef, combinedef, postwalk, striplines
import Libdl
using DocStringExtensions
using Base: RefValue
import IfElse

import Distributions

RuntimeGeneratedFunctions.init(@__MODULE__)

using RecursiveArrayTools

import SymbolicUtils
import SymbolicUtils: Term, Add, Mul, Pow, Sym, to_symbolic, FnType,
                      @rule, Rewriters, substitute, similarterm,
                      promote_symtype

import SymbolicUtils.Rewriters: Chain, Postwalk, Prewalk, Fixpoint

using LinearAlgebra: LU, BlasInt

import LightGraphs: SimpleDiGraph, add_edge!

import TreeViews

using Requires

export Num, Variable
"""
$(TYPEDEF)

Wrap anything in a type that is a subtype of Real
"""
struct Num <: Real
    val
end

const show_numwrap = Ref(false)

Num(x::Num) = x # ideally this should never be called
(n::Num)(args...) = Num(value(n)(map(value,args)...))
value(x) = x
value(x::Num) = x.val


using SymbolicUtils: to_symbolic
SymbolicUtils.to_symbolic(n::Num) = value(n)
SymbolicUtils.@number_methods(Num,
                              Num(f(value(a))),
                              Num(f(value(a), value(b))))

SymbolicUtils.simplify(n::Num; kw...) = Num(SymbolicUtils.simplify(value(n); kw...))

SymbolicUtils.symtype(n::Num) = symtype(n.val)

function Base.iszero(x::Num)
    _x = SymbolicUtils.to_mpoly(value(x))[1]
    return (_x isa Number || _x isa SymbolicUtils.MPoly) && iszero(_x)
end

import SymbolicUtils: <ₑ, Symbolic, Term, operation, arguments

Base.show(io::IO, n::Num) = show_numwrap[] ? print(io, :(Num($(value(n))))) : Base.show(io, value(n))

Base.promote_rule(::Type{<:Number}, ::Type{<:Num}) = Num
Base.promote_rule(::Type{<:Symbolic{<:Number}}, ::Type{<:Num}) = Num
function Base.getproperty(t::Union{Add, Mul, Pow, Term}, f::Symbol)
    if f === :op
        Base.depwarn("`x.op` is deprecated, use `operation(x)` instead", :getproperty, force=true)
        operation(t)
    elseif f === :args
        Base.depwarn("`x.args` is deprecated, use `arguments(x)` instead", :getproperty, force=true)
        arguments(t)
    else
        getfield(t, f)
    end
end
<ₑ(s::Num, x) = value(s) <ₑ value(x)
<ₑ(s, x::Num) = value(s) <ₑ value(x)
<ₑ(s::Num, x::Num) = value(s) <ₑ value(x)

for T in (Integer, Rational)
    @eval Base.:(^)(n::Num, i::$T) = Num(value(n)^i)
end

macro num_method(f, expr, Ts=nothing)
    if Ts === nothing
        Ts = [Any]
    else
        @assert Ts.head == :tuple
        # e.g. a tuple or vector
        Ts = Ts.args
    end

    ms = [quote
              $f(a::$T, b::$Num) = $expr
              $f(a::$Num, b::$T) = $expr
          end for T in Ts]
    quote
        $f(a::$Num, b::$Num) = $expr
        $(ms...)
    end |> esc
end

"""
    tosymbolic(a::Union{Sym,Num}) -> Sym{Real}
    tosymbolic(a::T) -> T
"""
tosymbolic(a::Num) = tosymbolic(value(a))
tosymbolic(a::Sym) = tovar(a)
tosymbolic(a) = a
@num_method Base.isless isless(tosymbolic(a), tosymbolic(b)) (Real,)
@num_method Base.:(<) (tosymbolic(a) < tosymbolic(b)) (Real,)
@num_method Base.:(<=) (tosymbolic(a) <= tosymbolic(b)) (Real,)
@num_method Base.:(>) (tosymbolic(a) > tosymbolic(b)) (Real,)
@num_method Base.:(>=) (tosymbolic(a) >= tosymbolic(b)) (Real,)
@num_method Base.isequal isequal(tosymbolic(a), tosymbolic(b)) (AbstractFloat, Number, Symbolic)
@num_method Base.:(==) tosymbolic(a) == tosymbolic(b) (AbstractFloat,Number)

Base.hash(x::Num, h::UInt) = hash(value(x), h)

Base.convert(::Type{Num}, x::Symbolic{<:Number}) = Num(x)
Base.convert(::Type{Num}, x::Number) = Num(x)
Base.convert(::Type{Num}, x::Num) = x

Base.convert(::Type{<:Array{Num}}, x::AbstractArray) = map(Num, x)
Base.convert(::Type{<:Array{Num}}, x::AbstractArray{Num}) = x
Base.convert(::Type{Sym}, x::Num) = value(x) isa Sym ? value(x) : error("cannot convert $x to Sym")

LinearAlgebra.lu(x::Array{Num}; check=true, kw...) = sym_lu(x; check=check)

"""
$(TYPEDEF)

TODO
"""
abstract type AbstractSystem end
abstract type AbstractODESystem <: AbstractSystem end

"""
$(TYPEDSIGNATURES)

Get the set of independent variables for the given system.
"""
function independent_variables end

"""
$(TYPEDSIGNATURES)

Get the set of states for the given system.
"""
function states end

"""
$(TYPEDSIGNATURES)

Get the set of parameters variables for the given system.
"""
function parameters end

include("bipartite_graph.jl")

include("variables.jl")
include("context_dsl.jl")
include("differentials.jl")

include("equations.jl")
include("utils.jl")
include("linearity.jl")
include("solve.jl")
include("direct.jl")
include("domains.jl")
include("register_function.jl")

include("systems/abstractsystem.jl")
include("systems/systemstructure.jl")

include("systems/diffeqs/odesystem.jl")
include("systems/diffeqs/sdesystem.jl")
include("systems/diffeqs/abstractodesystem.jl")
include("systems/diffeqs/first_order_transform.jl")
include("systems/diffeqs/modelingtoolkitize.jl")
include("systems/diffeqs/validation.jl")
include("systems/diffeqs/basic_transformations.jl")

include("systems/jumps/jumpsystem.jl")

include("systems/nonlinear/nonlinearsystem.jl")

include("systems/optimization/optimizationsystem.jl")

include("systems/control/controlsystem.jl")

include("systems/pde/pdesystem.jl")

include("systems/reaction/reactionsystem.jl")
include("systems/dependency_graphs.jl")

include("systems/reduction.jl")

include("latexify_recipes.jl")
include("build_function.jl")
include("extra_functions.jl")

export ODESystem, ODEFunction, ODEFunctionExpr, ODEProblemExpr
export SDESystem, SDEFunction, SDEFunctionExpr, SDESystemExpr
export SystemStructure
export JumpSystem
export ODEProblem, SDEProblem
export NonlinearProblem, NonlinearProblemExpr
export OptimizationProblem, OptimizationProblemExpr
export SteadyStateProblem, SteadyStateProblemExpr
export JumpProblem, DiscreteProblem
export NonlinearSystem, OptimizationSystem
export ControlSystem
export ode_order_lowering, liouville_transform
export runge_kutta_discretize
export PDESystem
export Reaction, ReactionSystem, ismassaction, oderatelaw, jumpratelaw
export Differential, expand_derivatives, @derivatives
export IntervalDomain, ProductDomain, ⊗, CircleDomain
export Equation, ConstrainedEquation
export Term, Sym
export independent_variable, states, parameters, equations, controls, pins, observed

export calculate_jacobian, generate_jacobian, generate_function
export calculate_tgrad, generate_tgrad
export calculate_gradient, generate_gradient
export calculate_factorized_W, generate_factorized_W
export calculate_hessian, generate_hessian
export calculate_massmatrix, generate_diffusion_function
export stochastic_integral_transform

export BipartiteGraph, equation_dependencies, variable_dependencies
export eqeq_dependencies, varvar_dependencies
export asgraph, asdigraph

export toexpr, get_variables
export simplify, substitute
export build_function
export @register
export modelingtoolkitize
export @variables, @parameters

const HAS_DAGGER = Ref{Bool}(false)
function __init__()
    @require Dagger="d58978e5-989f-55fb-8d15-ea34adc7bf54" include("dagger.jl")
end

end # module
