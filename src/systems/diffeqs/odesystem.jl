export ODESystem, ODEFunction

isintermediate(eq::Equation) = !(isa(eq.lhs, Operation) && isa(eq.lhs.op, Differential))

function flatten_differential(O::Operation)
    @assert is_derivative(O) "invalid differential: $O"
    is_derivative(O.args[1]) || return (O.args[1], O.op.x, 1)
    (x, t, order) = flatten_differential(O.args[1])
    isequal(t, O.op.x) || throw(ArgumentError("non-matching differentials on lhs: $t, $(O.op.x)"))
    return (x, t, order + 1)
end

struct ODEExpr  # dⁿx/dtⁿ = rhs
    x::Variable
    n::Int
    rhs::Expression
end
function Base.convert(::Type{ODEExpr},eq::Equation)
    isintermediate(eq) && throw(ArgumentError("intermediate equation received"))
    (x, t, n) = flatten_differential(eq.lhs)
    (isa(t, Operation) && isa(t.op, Variable) && isempty(t.args)) ||
        throw(ArgumentError("invalid independent variable $t"))
    (isa(x, Operation) && isa(x.op, Variable) && length(x.args) == 1 && isequal(first(x.args), t)) ||
        throw(ArgumentError("invalid dependent variable $x"))
    return t.op, ODEExpr(x.op, n, eq.rhs)
end
Base.:(==)(a::ODEExpr, b::ODEExpr) = isequal((a.x, a.n, a.rhs), (b.x, b.n, b.rhs))

"""
$(TYPEDEF)

A system of ordinary differential equations.

# Fields
* `eqs` - The ODEs defining the system.

# Examples

```
using ModelingToolkit

@parameters t σ ρ β
@variables x(t) y(t) z(t)
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

de = ODESystem(eqs)
```
"""
struct ODESystem <: AbstractODESystem
    """The ODEs defining the system."""
    eqs::Vector{ODEExpr}
    """Independent variable."""
    iv::Variable
    """Dependent (state) variables."""
    dvs::Vector{Variable}
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
end

function ODESystem(eqs)
    reformatted = convert.(ODEExpr,eqs)

    ivs = unique(r[1] for r ∈ reformatted)
    length(ivs) == 1 || throw(ArgumentError("one independent variable currently supported"))
    iv = first(ivs)

    deqs = [r[2] for r ∈ reformatted]

    dvs = unique(deq.x for deq ∈ deqs)
    ps = filter(vars(deq.rhs for deq ∈ deqs)) do x
        x.known & !isequal(x, iv)
    end |> collect

    ODESystem(deqs, iv, dvs, ps)
end

function ODESystem(deqs::AbstractVector{ODEExpr}, iv, dvs, ps)
    tgrad = RefValue(Vector{Expression}(undef, 0))
    jac = RefValue(Matrix{Expression}(undef, 0, 0))
    Wfact   = RefValue(Matrix{Expression}(undef, 0, 0))
    Wfact_t = RefValue(Matrix{Expression}(undef, 0, 0))
    ODESystem(deqs, iv, dvs, ps, tgrad, jac, Wfact, Wfact_t)
end

function ODESystem(deqs::AbstractVector{<:Equation}, iv, dvs, ps)
    _dvs = [deq.op for deq ∈ dvs]
    _iv = iv.op
    _ps = [p.op for p ∈ ps]
    ODESystem(getindex.(convert.(ODEExpr,deqs),2), _iv, _dvs, _ps)
end

function _eq_unordered(a, b)
    length(a) === length(b) || return false
    n = length(a)
    idxs = Set(1:n)
    for x ∈ a
        idx = findfirst(isequal(x), b)
        idx === nothing && return false
        idx ∈ idxs      || return false
        delete!(idxs, idx)
    end
    return true
end
Base.:(==)(sys1::ODESystem, sys2::ODESystem) =
    _eq_unordered(sys1.eqs, sys2.eqs) && isequal(sys1.iv, sys2.iv) &&
    _eq_unordered(sys1.dvs, sys2.dvs) && _eq_unordered(sys1.ps, sys2.ps)
# NOTE: equality does not check cached Jacobian

independent_variables(sys::ODESystem) = Set{Variable}([sys.iv])
dependent_variables(sys::ODESystem) = Set{Variable}(sys.dvs)
parameters(sys::ODESystem) = Set{Variable}(sys.ps)
