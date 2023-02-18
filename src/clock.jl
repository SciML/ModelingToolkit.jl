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
    is_continuous_domain(x)

true if `x` contains only continuous-domain signals.
See also [`has_continuous_domain`](@ref)
"""
function is_continuous_domain(x)
    issym(x) && return getmetadata(x, TimeDomain, false) isa Continuous
    !has_discrete_domain(x) && has_continuous_domain(x)
end

function get_time_domain(x)
    if istree(x) && operation(x) isa Operator
        output_timedomain(x)
    else
        getmetadata(x, TimeDomain, nothing)
    end
end
get_time_domain(x::Num) = get_time_domain(value(x))

"""
    has_time_domain(x)

Determine if variable `x` has a time-domain attributed to it.
"""
function has_time_domain(x::Symbolic)
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
function has_discrete_domain(x)
    issym(x) && return is_discrete_domain(x)
    hasshift(x) || hassample(x) || hashold(x)
end

"""
    has_continuous_domain(x)

true if `x` contains continuous signals (`x` may or may not contain discrete-domain signals). `x` may be an expression or equation.
See also [`is_continuous_domain`](@ref)
"""
function has_continuous_domain(x)
    issym(x) && return is_continuous_domain(x)
    hasderiv(x) || hasdiff(x) || hassample(x) || hashold(x)
end

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
function is_discrete_domain(x)
    issym(x) && return getmetadata(x, TimeDomain, false) isa Discrete
    !has_discrete_domain(x) && has_continuous_domain(x)
end

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
Base.:(==)(c1::Clock, c2::Clock) = isequal(c1.t, c2.t) && c1.dt == c2.dt
