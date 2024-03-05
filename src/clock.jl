abstract type TimeDomain end
abstract type AbstractDiscrete <: TimeDomain end

Base.Broadcast.broadcastable(d::TimeDomain) = Ref(d)

struct Inferred <: TimeDomain end
struct InferredDiscrete <: AbstractDiscrete end
struct Continuous <: TimeDomain end

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

for op in [Differential]
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
    if hasmetadata(x, TimeDomain) || issym(x)
        return getmetadata(x, TimeDomain, false) isa AbstractDiscrete
    end
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
    Clock([t]; dt)

The default periodic clock with independent variables `t` and tick interval `dt`.
If `dt` is left unspecified, it will be inferred (if possible).
"""
struct Clock <: AbstractClock
    "Independent variable"
    t::Union{Nothing, Symbolic}
    "Period"
    dt::Union{Nothing, Float64}
    Clock(t::Union{Num, Symbolic}, dt = nothing) = new(value(t), dt)
    Clock(t::Nothing, dt = nothing) = new(t, dt)
end
Clock(dt::Real) = Clock(nothing, dt)
Clock() = Clock(nothing, nothing)

sampletime(c) = isdefined(c, :dt) ? c.dt : nothing
Base.hash(c::Clock, seed::UInt) = hash(c.dt, seed ⊻ 0x953d7a9a18874b90)
function Base.:(==)(c1::Clock, c2::Clock)
    ((c1.t === nothing || c2.t === nothing) || isequal(c1.t, c2.t)) && c1.dt == c2.dt
end

is_concrete_time_domain(x) = x isa Union{AbstractClock, Continuous}

"""
    SolverStepClock <: AbstractClock
    SolverStepClock()
    SolverStepClock(t)

A clock that ticks at each solver step (sometimes referred to as "continuous sample time"). This clock **does generally not have equidistant tick intervals**, instead, the tick interval depends on the adaptive step-size slection of the continuous solver, as well as any continuous event handling. If adaptivity of the solver is turned off and there are no continuous events, the tick interval will be given by the fixed solver time step `dt`. 

Due to possibly non-equidistant tick intervals, this clock should typically not be used with discrete-time systems that assume a fixed sample time, such as PID controllers and digital filters.
"""
struct SolverStepClock <: AbstractClock
    "Independent variable"
    t::Union{Nothing, Symbolic}
    "Period"
    SolverStepClock(t::Union{Num, Symbolic}) = new(value(t))
end
SolverStepClock() = SolverStepClock(nothing)

Base.hash(c::SolverStepClock, seed::UInt) = seed ⊻ 0x953d7b9a18874b91
function Base.:(==)(c1::SolverStepClock, c2::SolverStepClock)
    ((c1.t === nothing || c2.t === nothing) || isequal(c1.t, c2.t))
end
