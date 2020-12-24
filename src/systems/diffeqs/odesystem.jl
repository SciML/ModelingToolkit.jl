"""
$(TYPEDEF)

A system of ordinary differential equations.

# Fields
$(FIELDS)

# Example

```
using ModelingToolkit

@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

de = ODESystem(eqs,t,[x,y,z],[σ,ρ,β])
```
"""
struct ODESystem <: AbstractODESystem
    """The ODEs defining the system."""
    eqs::Vector{Equation}
    """Independent variable."""
    iv::Sym
    """Dependent (state) variables."""
    states::Vector
    """Parameter variables."""
    ps::Vector
    pins::Vector{Num}
    observed::Vector{Equation}
    """
    Time-derivative matrix. Note: this field will not be defined until
    [`calculate_tgrad`](@ref) is called on the system.
    """
    tgrad::RefValue{Vector{Num}}
    """
    Jacobian matrix. Note: this field will not be defined until
    [`calculate_jacobian`](@ref) is called on the system.
    """
    jac::RefValue{Any}
    """
    `Wfact` matrix. Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact::RefValue{Matrix{Num}}
    """
    `Wfact_t` matrix. Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact_t::RefValue{Matrix{Num}}
    """
    Name: the name of the system
    """
    name::Symbol
    """
    systems: The internal systems
    """
    systems::Vector{ODESystem}
end

function ODESystem(deqs::AbstractVector{<:Equation}, iv, dvs, ps;
                   pins = Num[],
                   observed = Num[],
                   systems = ODESystem[],
                   name=gensym(:ODESystem))
    iv′ = value(iv)
    dvs′ = value.(dvs)
    ps′ = value.(ps)
    tgrad = RefValue(Vector{Num}(undef, 0))
    jac = RefValue{Any}(Matrix{Num}(undef, 0, 0))
    Wfact   = RefValue(Matrix{Num}(undef, 0, 0))
    Wfact_t = RefValue(Matrix{Num}(undef, 0, 0))
    ODESystem(deqs, iv′, dvs′, ps′, pins, observed, tgrad, jac, Wfact, Wfact_t, name, systems)
end

var_from_nested_derivative(x, i=0) = (missing, missing)
var_from_nested_derivative(x::Term,i=0) = x.op isa Differential ? var_from_nested_derivative(x.args[1],i+1) : (x,i)
var_from_nested_derivative(x::Sym,i=0) = (x,i)

iv_from_nested_derivative(x::Term) = x.op isa Differential ? iv_from_nested_derivative(x.args[1]) : x.args[1]
iv_from_nested_derivative(x::Sym) = x
iv_from_nested_derivative(x) = missing

vars(exprs::Term) = vars([exprs])
vars(exprs) = foldl(vars!, exprs; init = Set())
function vars!(vars, O)
    isa(O, Sym) && return push!(vars, O)
    !isa(O, Term) && return vars

    O.op isa Sym && push!(vars, O)
    for arg ∈ O.args
        vars!(vars, arg)
    end

    return vars
end

function ODESystem(eqs, iv=nothing; kwargs...)
    # NOTE: this assumes that the order of algebric equations doesn't matter
    diffvars = OrderedSet()
    allstates = OrderedSet()
    ps = OrderedSet()
    # reorder equations such that it is in the form of `diffeq, algeeq`
    diffeq = Equation[]
    algeeq = Equation[]
    # initial loop for finding `iv`
    if iv === nothing
        for eq in eqs
            if !(eq.lhs isa Number) # assume eq.lhs is either Differential or Number
                iv = iv_from_nested_derivative(eq.lhs)
                break
            end
        end
    end
    iv === nothing && throw(ArgumentError("Please pass in independent variables."))
    for eq in eqs
        for var in vars(eq.rhs for eq ∈ eqs)
            if isparameter(var) || isparameter(var.op)
                isequal(var, iv) || push!(ps, var)
            else
                push!(allstates, var)
            end
        end
        if !(eq.lhs isa Symbolic)
            push!(algeeq, eq)
        else
            diffvar = first(var_from_nested_derivative(eq.lhs))
            isequal(iv, iv_from_nested_derivative(eq.lhs)) || throw(ArgumentError("An ODESystem can only have one independent variable."))
            diffvar in diffvars && throw(ArgumentError("The differential variable $diffvar is not unique in the system of equations."))
            push!(diffvars, diffvar)
            push!(diffeq, eq)
        end
    end
    algevars = setdiff(allstates, diffvars)
    # the orders here are very important!
    return ODESystem(append!(diffeq, algeeq), iv, vcat(collect(diffvars), collect(algevars)), ps; kwargs...)
end

Base.:(==)(sys1::ODESystem, sys2::ODESystem) =
    _eq_unordered(sys1.eqs, sys2.eqs) && isequal(sys1.iv, sys2.iv) &&
    _eq_unordered(sys1.states, sys2.states) && _eq_unordered(sys1.ps, sys2.ps)
# NOTE: equality does not check cached Jacobian

function rename(sys::ODESystem,name)
    ODESystem(sys.eqs, sys.iv, sys.states, sys.ps, sys.pins, sys.observed, sys.tgrad, sys.jac, sys.Wfact, sys.Wfact_t, name, sys.systems)
end

isdiffeq(eq) = eq.lhs isa Term && operation(eq.lhs) isa Differential
