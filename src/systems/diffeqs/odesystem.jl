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
D = Differential(t)

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
    systems: The internal systems. These are required to have unique names.
    """
    systems::Vector{ODESystem}
    """
    default_u0: The default initial conditions to use when initial conditions
    are not supplied in `ODEProblem`.
    """
    default_u0::Dict
    """
    default_p: The default parameters to use when parameters are not supplied
    in `ODEProblem`.
    """
    default_p::Dict
    """
    structure: structural information of the system
    """
    structure::Any
end

function ODESystem(
                   deqs::AbstractVector{<:Equation}, iv, dvs, ps;
                   observed = Num[],
                   systems = ODESystem[],
                   name=gensym(:ODESystem),
                   default_u0=Dict(),
                   default_p=Dict(),
                  )
    iv′ = value(iv)
    dvs′ = value.(dvs)
    ps′ = value.(ps)

    default_u0 isa Dict || (default_u0 = Dict(default_u0))
    default_p isa Dict || (default_p = Dict(default_p))
    default_u0 = Dict(value(k) => value(default_u0[k]) for k in keys(default_u0))
    default_p = Dict(value(k) => value(default_p[k]) for k in keys(default_p))

    tgrad = RefValue(Vector{Num}(undef, 0))
    jac = RefValue{Any}(Matrix{Num}(undef, 0, 0))
    Wfact   = RefValue(Matrix{Num}(undef, 0, 0))
    Wfact_t = RefValue(Matrix{Num}(undef, 0, 0))
    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end
    ODESystem(deqs, iv′, dvs′, ps′, observed, tgrad, jac, Wfact, Wfact_t, name, systems, default_u0, default_p, nothing)
end

var_from_nested_derivative(x, i=0) = (missing, missing)
var_from_nested_derivative(x::Term,i=0) = operation(x) isa Differential ? var_from_nested_derivative(arguments(x)[1],i+1) : (x,i)
var_from_nested_derivative(x::Sym,i=0) = (x,i)

iv_from_nested_derivative(x::Term) = operation(x) isa Differential ? iv_from_nested_derivative(arguments(x)[1]) : arguments(x)[1]
iv_from_nested_derivative(x::Sym) = x
iv_from_nested_derivative(x) = missing

vars(x::Sym) = [x]
vars(exprs::Symbolic) = vars([exprs])
vars(exprs) = foldl(vars!, exprs; init = Set())
vars!(vars, eq::Equation) = (vars!(vars, eq.lhs); vars!(vars, eq.rhs); vars)
function vars!(vars, O)
    isa(O, Sym) && return push!(vars, O)
    !istree(O) && return vars

    operation(O) isa Differential && return push!(vars, O)

    operation(O) isa Sym && push!(vars, O)
    for arg in arguments(O)
        vars!(vars, arg)
    end

    return vars
end

find_derivatives!(vars, expr::Equation, f=identity) = (find_derivatives!(vars, expr.lhs, f); find_derivatives!(vars, expr.rhs, f); vars)
function find_derivatives!(vars, expr, f)
    !istree(O) && return vars
    operation(O) isa Differential && push!(vars, f(O))
    for arg in arguments(O)
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
    iv = value(iv)
    iv === nothing && throw(ArgumentError("Please pass in independent variables."))
    for eq in eqs
        collect_vars!(allstates, ps, eq.lhs, iv)
        collect_vars!(allstates, ps, eq.rhs, iv)
        if isdiffeq(eq)
            diffvar, _ = var_from_nested_derivative(eq.lhs)
            isequal(iv, iv_from_nested_derivative(eq.lhs)) || throw(ArgumentError("An ODESystem can only have one independent variable."))
            diffvar in diffvars && throw(ArgumentError("The differential variable $diffvar is not unique in the system of equations."))
            push!(diffvars, diffvar)
            push!(diffeq, eq)
        else
            push!(algeeq, eq)
        end
    end
    algevars = setdiff(allstates, diffvars)
    # the orders here are very important!
    return ODESystem(append!(diffeq, algeeq), iv, vcat(collect(diffvars), collect(algevars)), ps; kwargs...)
end

function collect_vars!(states, parameters, expr, iv)
    if expr isa Sym
        collect_var!(states, parameters, expr, iv)
    else
        for var in vars(expr)
            if istree(var) && operation(var) isa Differential
                var, _ = var_from_nested_derivative(var)
            end
            collect_var!(states, parameters, var, iv)
        end
    end
    return nothing
end

function collect_var!(states, parameters, var, iv)
    isequal(var, iv) && return nothing
    if isparameter(var) || (istree(var) && isparameter(operation(var)))
        push!(parameters, var)
    else
        push!(states, var)
    end
    return nothing
end

# NOTE: equality does not check cached Jacobian
function Base.:(==)(sys1::ODESystem, sys2::ODESystem)
    iv1 = independent_variable(sys1)
    iv2 = independent_variable(sys2)
    isequal(iv1, iv2) &&
    _eq_unordered(get_eqs(sys1), get_eqs(sys2)) &&
    _eq_unordered(get_states(sys1), get_states(sys2)) &&
    _eq_unordered(get_ps(sys1), get_ps(sys2)) &&
    all(s1 == s2 for (s1, s2) in zip(get_systems(sys1), get_systems(sys2)))
end

function flatten(sys::ODESystem)
    systems = get_systems(sys)
    if isempty(systems)
        return sys
    else
        return ODESystem(
                         equations(sys),
                         independent_variable(sys),
                         states(sys),
                         parameters(sys),
                         observed=observed(sys),
                         default_u0=default_u0(sys),
                         default_p=default_p(sys),
                        )
    end
end

ODESystem(eq::Equation, args...; kwargs...) = ODESystem([eq], args...; kwargs...)
