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
    Clock()
    Clock(dt)
    Clock(t, dt)

The default periodic clock with independent variables `t` and tick interval `dt`.
If `dt` is left unspecified, it will be inferred (if possible).

See also [`RealtimeClock`](@ref).
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
    RealtimeClock()
    RealtimeClock(dt)
    RealtimeClock(t, dt)

Similar to [`Clock`](@ref), but with with the additional property that the simulation is run no faster than real-time, i.e., a simulation with `tspan = (0, 10)` will take at least 10 seconds of wallclock time to run.

To achieve real-time execution, the wall-clock duration of each integration step is measured, and a call to `Libc.systemsleep` is made to ensure that the next integration step is not started before `dt` seconds of wall-clock time has elapsed. This leads to a simulation that takes at least as long as the length of the time span, but may take longer if the computational load is high enough that one integration step takes more than `dt` seconds.
"""
struct RealtimeClock <: AbstractClock
    clock::Clock
    RealtimeClock(args...) = new(Clock(args...))
end

sampletime(c) = sampletime(c.clock)
Base.hash(c::RealtimeClock, seed::UInt) = hash(c.dt, seed ⊻ 0x9d3d7a9a18874b90)
Base.:(==)(c1::RealtimeClock, c2::RealtimeClock) = c1.clock == c2.clock
function Base.getproperty(c::RealtimeClock, s::Symbol)
    s ∈ fieldnames(typeof(c)) && return getfield(c, s)
    getproperty(getfield(c, :clock), s)
end
