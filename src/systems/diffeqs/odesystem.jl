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
var_from_nested_derivative(x::Constant) = (missing, missing)
var_from_nested_derivative(x,i) = x.op isa Differential ? var_from_nested_derivative(x.args[1],i+1) : (x.op,i)

function extract_eqs_states_ps(eqs::AbstractArray{<:Equation}, iv)
    # NOTE: this assumes that the order of algebric equations doesn't matter
    diffvars = OrderedSet{Variable}()
    allstates = OrderedSet{Variable}()
    ps = OrderedSet{Variable}()
    # reorder equations such that it is in the form of `diffeq, algeeq`
    diffeq = Equation[]
    algeeq = Equation[]
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
            diffvar in diffvars && throw(ArgumentError("The differential variable $diffvar is not unique in the system of equations."))
            push!(diffvars, diffvar)
            push!(diffeq, eq)
        end
    end
    algevars = setdiff(allstates, diffvars)
    # the orders here are very important!
    return append!(diffeq, algeeq), vcat(collect(diffvars), collect(algevars)), ps
end

iv_from_nested_derivative(x) = x.op isa Differential ? iv_from_nested_derivative(x.args[1]) : x.args[1].op
iv_from_nested_derivative(x::Constant) = missing

function ODESystem(eqs; kwargs...)
    ivs = unique(skipmissing(iv_from_nested_derivative(eq.lhs) for eq ∈ eqs))
    length(ivs) == 1 || throw(ArgumentError("An ODESystem can only have one independent variable."))
    iv = first(ivs)
    eqs, dvs, ps = extract_eqs_states_ps(eqs, iv)
    return ODESystem(eqs, iv, dvs, ps; kwargs...)
end

Base.:(==)(sys1::ODESystem, sys2::ODESystem) =
    _eq_unordered(sys1.eqs, sys2.eqs) && isequal(sys1.iv, sys2.iv) &&
    _eq_unordered(sys1.states, sys2.states) && _eq_unordered(sys1.ps, sys2.ps)
# NOTE: equality does not check cached Jacobian

function rename(sys::ODESystem,name)
    ODESystem(sys.eqs, sys.iv, sys.states, sys.ps, sys.tgrad, sys.jac, sys.Wfact, sys.Wfact_t, name, sys.systems)
end
