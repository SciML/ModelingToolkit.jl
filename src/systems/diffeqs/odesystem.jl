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
    iv::Variable
    """Dependent (state) variables."""
    states::Vector{Variable}
    """Parameter variables."""
    ps::Vector{Variable}
    pins::Vector{Variable}
    observed::Vector{Equation}
    """
    Time-derivative matrix. Note: this field will not be defined until
    [`calculate_tgrad`](@ref) is called on the system.
    """
    tgrad::RefValue{Vector{Expression}}
    """
    Jacobian matrix. Note: this field will not be defined until
    [`calculate_jacobian`](@ref) is called on the system.
    """
    jac::RefValue{Any}
    """
    `Wfact` matrix. Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact::RefValue{Matrix{Expression}}
    """
    `Wfact_t` matrix. Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact_t::RefValue{Matrix{Expression}}
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
                   pins = Variable[],
                   observed = Operation[],
                   systems = ODESystem[],
                   name=gensym(:ODESystem))
    iv′ = convert(Variable,iv)
    dvs′ = convert.(Variable,dvs)
    ps′ = convert.(Variable,ps)
    tgrad = RefValue(Vector{Expression}(undef, 0))
    jac = RefValue{Any}(Matrix{Expression}(undef, 0, 0))
    Wfact   = RefValue(Matrix{Expression}(undef, 0, 0))
    Wfact_t = RefValue(Matrix{Expression}(undef, 0, 0))
    ODESystem(deqs, iv′, dvs′, ps′, pins, observed, tgrad, jac, Wfact, Wfact_t, name, systems)
end

var_from_nested_derivative(x::Constant) = (missing, missing)
var_from_nested_derivative(x,i=0) = x.op isa Differential ? var_from_nested_derivative(x.args[1],i+1) : (x.op,i)

iv_from_nested_derivative(x) = x.op isa Differential ? iv_from_nested_derivative(x.args[1]) : x.args[1].op
iv_from_nested_derivative(x::Constant) = missing

function ODESystem(eqs, iv=nothing; kwargs...)
    # NOTE: this assumes that the order of algebric equations doesn't matter
    diffvars = OrderedSet{Variable}()
    allstates = OrderedSet{Variable}()
    ps = OrderedSet{Variable}()
    # reorder equations such that it is in the form of `diffeq, algeeq`
    diffeq = Equation[]
    algeeq = Equation[]
    # initial loop for finding `iv`
    if iv === nothing
        for eq in eqs
            if !(eq.lhs isa Constant) # assume eq.lhs is either Differential or Constant
                iv = iv_from_nested_derivative(eq.lhs)
                break
            end
        end
    else
        iv = convert(Variable, iv)
    end
    iv === nothing && throw(ArgumentError("Please pass in independent variables."))
    for eq in eqs
        for var in vars(eq.rhs for eq ∈ eqs)
            var isa Variable || continue
            if isparameter(var)
                isequal(var, iv) || push!(ps, var)
            else
                push!(allstates, var)
            end
        end
        if eq.lhs isa Constant
            push!(algeeq, eq)
        else
            diffvar = first(var_from_nested_derivative(eq.lhs))
            iv == iv_from_nested_derivative(eq.lhs) || throw(ArgumentError("An ODESystem can only have one independent variable."))
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
