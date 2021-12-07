using Symbolics: Operator, Num, Term, value, recursive_hasoperator
using Symbolics: Differential, Difference
using Symbolics: hasderiv, hasdiff

abstract type TimeDomain end
abstract type AbstractDiscrete <: TimeDomain end

Base.Broadcast.broadcastable(d::TimeDomain) = Ref(d)


struct Inferred <: TimeDomain end
struct InferredDiscrete <: AbstractDiscrete end
struct Continuous <: TimeDomain end

const UnknownDomain = Union{Nothing, Inferred, InferredDiscrete}
const InferredDomain = Union{Inferred, InferredDiscrete}

Symbolics.option_to_metadata_type(::Val{:timedomain}) = TimeDomain

"""
    is_continuous_domain(x::Sym)

Determine if variable `x` is a continuous-time variable.
"""
is_continuous_domain(x::Sym) = getmetadata(x, TimeDomain, false) isa Continuous
# is_continuous_domain(x::Sym) = isvarkind(Continuous, x)

"""
    is_discrete_domain(x::Sym)

Determine if variable `x` is a discrete-time variable.
"""
is_discrete_domain(x::Sym) = getmetadata(x, TimeDomain, false) isa Discrete
# is_discrete_domain(x::Sym) = isvarkind(Discrete, x)


has_continuous_domain(x::Sym) = is_continuous_domain(x)
has_discrete_domain(x::Sym) = is_discrete_domain(x)

get_time_domain(x) = getmetadata(x, TimeDomain, nothing)
get_time_domain(x::Num) = get_time_domain(value(x))

"""
    has_time_domain(x::Sym)

Determine if variable `x` has a time-domain attributed to it.
"""
function has_time_domain(x::Union{Sym, Term})
    # getmetadata(x, Continuous, nothing) !== nothing || 
    # getmetadata(x, Discrete,   nothing) !== nothing
    getmetadata(x, TimeDomain, nothing) !== nothing
end
has_time_domain(x::Num) = has_time_domain(value(x))
has_time_domain(x) = false




for op in [Differential, Difference]
    @eval input_timedomain(::$op, arg=nothing) = Continuous()
    @eval output_timedomain(::$op, arg=nothing) = Continuous()
end


"""
    sampletime(op::Operator)

Return the sample time for operators that are clocked. Returns `nothing` if the operator does not have a sample time.
"""
sampletime(op::Operator) = hasfield(typeof(op), :dt) ? op.dt : nothing

"""
    has_discrete_domain(x)

true if `x` contains discrete signals (`x` may or may not contain continuous-domain signals). `x` may be an expression or equation.
See also [`is_discrete_domain`](@ref)
"""
has_discrete_domain(x) = hasshift(x) || hassample(x) || hashold(x)

"""
    has_continuous_domain(x)

true if `x` contains continuous signals (`x` may or may not contain discrete-domain signals). `x` may be an expression or equation.
See also [`is_continuous_domain`](@ref)
"""
has_continuous_domain(x) = hasderiv(x) || hasdiff(x) || hassample(x) || hashold(x)


"""
    is_hybrid_domain(x)

true if `x` contains both discrete and continuous-domain signals. `x` may be an expression or equation.
"""
is_hybrid_domain(x) = has_discrete_domain(x) && has_continuous_domain(x)

"""
    is_discrete_domain(x)

true if `x` contains only discrete-domain signals.
See also [`has_discrete_domain`](@ref)
"""
is_discrete_domain(x) = has_discrete_domain(x) && !has_continuous_domain(x)

"""
    is_continuous_domain(x)

true if `x` contains only continuous-domain signals.
See also [`has_continuous_domain`](@ref)
"""
is_continuous_domain(x) = !has_discrete_domain(x) && has_continuous_domain(x)