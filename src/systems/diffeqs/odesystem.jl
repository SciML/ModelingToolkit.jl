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

de = ODESystem(eqs,[x,y,z],[σ,ρ,β])
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
    """
    Time-derivative matrix. Note: this field will not be defined until
    [`calculate_tgrad`](@ref) is called on the system.
    """
    tgrad::RefValue{Vector{Expression}}
    """
    Jacobian matrix. Note: this field will not be defined until
    [`calculate_jacobian`](@ref) is called on the system.
    """
    jac::RefValue{Matrix{Expression}}
    """
    Wfact matrix. Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact::RefValue{Matrix{Expression}}
    """
    Wfact_t matrix. Note: this field will not be defined until
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
                   systems = ODESystem[],
                   name=gensym(:ODESystem))
    iv′ = convert(Variable,iv)
    dvs′ = convert.(Variable,dvs)
    ps′ = convert.(Variable,ps)
    tgrad = RefValue(Vector{Expression}(undef, 0))
    jac = RefValue(Matrix{Expression}(undef, 0, 0))
    Wfact   = RefValue(Matrix{Expression}(undef, 0, 0))
    Wfact_t = RefValue(Matrix{Expression}(undef, 0, 0))
    ODESystem(deqs, iv′, dvs′, ps′, tgrad, jac, Wfact, Wfact_t, name, systems)
end

var_from_nested_derivative(x) = var_from_nested_derivative(x,0)
var_from_nested_derivative(x,i) = x.op isa Differential ? var_from_nested_derivative(x.args[1],i+1) : (x.op,i)
iv_from_nested_derivative(x) = x.op isa Differential ? iv_from_nested_derivative(x.args[1]) : x.args[1].op

function ODESystem(eqs; kwargs...)
    ivs = unique(iv_from_nested_derivative(eq.lhs) for eq ∈ eqs)
    length(ivs) == 1 || throw(ArgumentError("one independent variable currently supported"))
    iv = first(ivs)

    dvs = unique(var_from_nested_derivative(eq.lhs)[1] for eq ∈ eqs)
    ps = filter(vars(eq.rhs for eq ∈ eqs)) do x
        isparameter(x) & !isequal(x, iv)
    end |> collect
    ODESystem(eqs, iv, dvs, ps; kwargs...)
end

Base.:(==)(sys1::ODESystem, sys2::ODESystem) =
    _eq_unordered(sys1.eqs, sys2.eqs) && isequal(sys1.iv, sys2.iv) &&
    _eq_unordered(sys1.states, sys2.states) && _eq_unordered(sys1.ps, sys2.ps)
# NOTE: equality does not check cached Jacobian

function rename(sys::ODESystem,name)
    ODESystem(sys.eqs, sys.iv, sys.states, sys.ps, sys.tgrad, sys.jac, sys.Wfact, sys.Wfact_t, name, sys.systems)
end
