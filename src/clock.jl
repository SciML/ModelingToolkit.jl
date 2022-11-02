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

"""
    is_discrete_domain(x::Sym)
Determine if variable `x` is a discrete-time variable.
"""
is_discrete_domain(x::Sym) = getmetadata(x, TimeDomain, false) isa Discrete

# is_discrete_domain(x::Sym) = isvarkind(Discrete, x)

has_continuous_domain(x::Sym) = is_continuous_domain(x)
has_discrete_domain(x::Sym) = is_discrete_domain(x)

function get_time_domain(x)
    if istree(x) && operation(x) isa Operator
        output_timedomain(x)
    else
        getmetadata(x, TimeDomain, nothing)
    end
end
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
    @eval input_timedomain(::$op, arg = nothing) = Continuous()
    @eval output_timedomain(::$op, arg = nothing) = Continuous()
end

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

struct ClockInferenceException <: Exception
    msg::Any
end

function Base.showerror(io::IO, cie::ClockInferenceException)
    print(io, "ClockInferenceException: ", cie.msg)
end

abstract type AbstractClock <: AbstractDiscrete end

"""
Clock <: AbstractClock
Clock(t; dt)
The default periodic clock with independent variables `t` and tick interval `dt`.
If `dt` is left unspecified, it will be inferred (if possible).
"""
struct Clock <: AbstractClock
    "Independent variable"
    t::Any
    "Period"
    dt::Any
    Clock(t, dt = nothing) = new(value(t), dt)
end

sampletime(c) = isdefined(c, :dt) ? c.dt : nothing
