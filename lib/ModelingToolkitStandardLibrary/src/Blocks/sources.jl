using DiffEqBase
import ChainRulesCore
using PreallocationTools

# Define and register smooth functions
# These are "smooth" aka differentiable and avoid Gibbs effect
# These follow: `offset` + `smooth_wave` * `smooth_step` with zero output for `t < start_time`
function smooth_cos(x, δ, f, amplitude, ϕ, offset, start_time)
    offset +
    amplitude * cos(2 * π * f * (x - start_time) + ϕ) *
    smooth_step(x, δ, one(x), zero(x), start_time)
end

function smooth_damped_sin(x, δ, f, amplitude, damping, ϕ, offset, start_time)
    offset +
    exp((start_time - x) * damping) * amplitude * sin(2 * π * f * (x - start_time) + ϕ) *
    smooth_step(x, δ, one(x), zero(x), start_time)
end

function smooth_ramp(x, δ, height, duration, offset, start_time)
    offset +
    height / (duration) *
    (smooth_xH(x, δ, start_time) - smooth_xH(x, δ, start_time + duration))
end

function smooth_sin(x, δ, f, amplitude, ϕ, offset, start_time)
    offset +
    amplitude * sin(2 * pi * f * (x - start_time) + ϕ) *
    smooth_step(x, δ, one(x), zero(x), start_time)
end

function smooth_square(x, δ, f, amplitude, offset, start_time)
    offset +
    amplitude * 2atan(sin(2π * (x - start_time) * f) / δ) / π *
    smooth_step(x, δ, one(x), zero(x), start_time)
end

function smooth_step(x, δ, height, offset, start_time)
    offset + height * (atan((x - start_time) / δ) / π + 0.5)
end

function smooth_triangular(x, δ, f, amplitude, offset, start_time)
    offset +
    amplitude * (1 - 2acos((1 - δ)sin(2π * (x - start_time) * f)) / π) *
    smooth_step(x, δ, one(x), zero(x), start_time)
end

function smooth_xH(x, δ, tₒ)
    0.5 * (x - tₒ) * (1 + ((x - tₒ) / sqrt((x - tₒ)^2 + δ^2)))
end

function square(x, f, amplitude, offset, start_time)
    offset +
    (x > start_time) * (amplitude *
     (4 * floor(f * (x - start_time)) - 2 * floor(2 * (x - start_time) * f) + 1))
end

function triangular(x, f, amplitude, offset, start_time)
    p = 1 / f # period
    offset +
    (x > start_time) *
    (4 * amplitude * f * abs(abs((x - p / 4 - start_time) % p) - p / 2) - amplitude)
end

"""
    Constant(; name, k = 0.0)

Generate constant signal.

# Parameters:

  - `k`: Constant output value

# Connectors:

  - `output`
"""
@mtkmodel Constant begin
    @components begin
        output = RealOutput()
    end
    @parameters begin
        k = 0.0, [description = "Constant output value of block"]
    end
    @equations begin
        output.u ~ k
    end
end

"""
    TimeVaryingFunction(f; name)

Outputs ``f(t)``.

The input variable `t` can be changed by passing a different variable as the keyword argument `t`.

# Connectors:
- `output`
"""
@mtkmodel TimeVaryingFunction begin
    @structural_parameters begin
        f
    end
    @components begin
        output = RealOutput()
    end
    @equations begin
        output.u ~ f(t)
    end
end
TimeVaryingFunction.f(f; name) = TimeVaryingFunction(; f, name)

"""
    Sine(; name, frequency, amplitude = 1, phase = 0, offset = 0, start_time = 0,
    smooth = false)

Generate sine signal.

# Parameters:

  - `frequency`: [Hz] Frequency of sine wave
  - `amplitude`: Amplitude of sine wave
  - `phase`: [rad] Phase of sine wave
  - `offset`: Offset of output signal
  - `start_time`: [s] Output `y = offset` for `t < start_time`
  - `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function Sine(; name,
        frequency,
        amplitude = 1,
        phase = 0,
        offset = 0,
        start_time = 0,
        smooth = false)
    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time amplitude=amplitude frequency=frequency phase=phase
    equation = if smooth == false
        offset + ifelse(t < start_time, 0,
            amplitude * sin(2 * pi * frequency * (t - start_time) + phase))
    else
        smooth === true && (smooth = 1e-5)
        smooth_sin(t, smooth, frequency, amplitude, phase, offset, start_time)
    end

    eqs = [
        output.u ~ equation
    ]

    compose(System(eqs, t, [], pars; name = name), [output])
end

"""
    Cosine(; name, frequency, amplitude = 1, phase = 0, offset = 0, start_time = 0,
    smooth = false)

Generate cosine signal.

# Parameters:
- `frequency`: [Hz] Frequency of cosine wave
- `amplitude`: Amplitude of cosine wave
- `phase`: [rad] Phase of cosine wave
- `offset`: Offset of output signal
- `start_time`: [s] Output `y = offset` for `t < start_time`
- `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
             It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:
- `output`
"""
@component function Cosine(; name,
        frequency,
        amplitude = 1,
        phase = 0,
        offset = 0,
        start_time = 0,
        smooth = false)
    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time amplitude=amplitude frequency=frequency phase=phase
    equation = if smooth == false
        offset + ifelse(t < start_time, zero(t),
            amplitude * cos(2 * pi * frequency * (t - start_time) + phase))
    else
        smooth === true && (smooth = 1e-5)
        smooth_cos(t, smooth, frequency, amplitude, phase, offset, start_time)
    end
    eqs = [
        output.u ~ equation
    ]

    compose(System(eqs, t, [], pars; name = name), [output])
end

"""
    ContinuousClock(; name, offset = 0, start_time = 0)

Generate current time signal.

# Parameters:

  - `offset`: Offset of output signal
  - `start_time`: [s] Output `y = offset` for `t < start_time`

# Connectors:

  - `output`
"""
@component function ContinuousClock(; name, offset = 0, start_time = 0)
    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time
    eqs = [
        output.u ~ offset + ifelse(t < start_time, zero(t), t - start_time)
    ]

    compose(System(eqs, t, [], pars; name = name), [output])
end

"""
Ramp(; name, height = 1, duration = 1, offset = 0, start_time = 0, smooth = false)

Generate ramp signal.

# Parameters:

  - `height`: Height of ramp
  - `duration`: [s] Duration of ramp (= 0.0 gives a Step)
  - `offset`: Offset of output signal
  - `start_time`: [s] Output `y = offset` for `t < start_time`
  - `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function Ramp(; name,
        height = 1.0,
        duration = 1.0,
        offset = 0.0,
        start_time = 0.0,
        smooth = false)
    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time height=height duration=duration
    equation = if smooth == false
        offset + ifelse(t < start_time, zero(height),
            ifelse(t < (start_time + duration), (t - start_time) * height / duration,
                height))
    else
        smooth === true && (smooth = 1e-5)
        smooth_ramp(t, smooth, height, duration, offset, start_time)
    end

    eqs = [
        output.u ~ equation
    ]

    compose(System(eqs, t, [], pars; name = name), [output])
end

"""
    Square(; name, frequency = 1.0, amplitude = 1.0, offset = 0.0, start_time = 0.0,
    smooth = false)
Generate smooth square signal.

# Parameters:

  - `frequency`: [Hz] Frequency of square wave
  - `amplitude`: Amplitude of square wave
  - `offset`: Offset of output signal
  - `start_time`: [s] Output `y = offset` for `t < start_time`
  - `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function Square(; name, frequency = 1.0, amplitude = 1.0,
        offset = 0.0, start_time = 0.0, smooth = false)
    @named output = RealOutput()
    pars = @parameters begin
        frequency = frequency
        amplitude = amplitude
        offset = offset
        start_time = start_time
    end

    equation = if smooth == false
        square(t, frequency, amplitude, offset, start_time)
    else
        smooth === true && (smooth = 1e-5)
        smooth_square(t, smooth, frequency, amplitude, offset, start_time)
    end

    eqs = [
        output.u ~ equation
    ]

    compose(System(eqs, t, [], pars; name = name), [output])
end

"""
    Step(;name, height=1, offset=0, start_time=0, duration=Inf, smooth=true)

Generate step signal.

# Parameters:

  - `height`: Height of step
  - `offset`: Offset of output signal
  - `start_time`: [s] Output `y = offset` for `t < start_time` and thereafter `offset+height`.
  - `duration`: [s] If `duration < Inf` is supplied, the output will revert to `offset` after `duration` seconds.
  - `smooth`:  If `true`, returns a smooth wave. Defaults to `true`
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function Step(;
        name, height = 1.0, offset = 0.0, start_time = 0.0, duration = Inf,
        smooth = 1e-5)
    @named output = RealOutput()
    duration_numeric = duration
    pars = @parameters offset=offset start_time=start_time height=height duration=duration
    equation = if smooth == false # use comparison in case smooth is a float
        offset +
        ifelse((start_time <= t) & (t < start_time + duration), height, zero(height))
    else
        smooth === true && (smooth = 1e-5)
        if duration_numeric == Inf
            smooth_step(t, smooth, height, offset, start_time)
        else
            smooth_step(t, smooth, height, offset, start_time) -
            smooth_step(t, smooth, height, zero(start_time), start_time + duration)
        end
    end

    eqs = [
        output.u ~ equation
    ]

    compose(System(eqs, t, [], pars; name = name), [output])
end

"""
    ExpSine(; name, frequency, amplitude = 1, damping = 0.1, phase = 0, offset = 0, start_time = 0, smooth = false)

Exponentially damped sine signal.

# Parameters:

  - `frequency`: [Hz] Frequency of sine wave
  - `amplitude`: Amplitude of sine wave
  - `damping`: [1/s] Damping coefficient of sine wave
  - `phase`: [rad] Phase of sine wave
  - `offset`: Offset of output signal
  - `start_time`: [s] Output `y = offset` for `t < start_time`
  - `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function ExpSine(; name,
        frequency,
        amplitude = 1.0,
        damping = 0.1,
        phase = 0.0,
        offset = 0.0,
        start_time = 0.0,
        smooth = false)
    @named output = RealOutput()
    pars = @parameters offset=offset start_time=start_time amplitude=amplitude frequency=frequency phase=phase damping=damping

    equation = if smooth == false
        offset + ifelse(t < start_time, zero(amplitude),
            amplitude * exp(-damping * (t - start_time)) *
            sin(2 * pi * frequency * (t - start_time) + phase))
    else
        smooth === true && (smooth = 1e-5)
        smooth_damped_sin(t, smooth, frequency, amplitude, damping, phase, offset,
            start_time)
    end

    eqs = [
        output.u ~ equation
    ]

    compose(System(eqs, t, [], pars; name = name), [output])
end

"""
    Triangular(; name, amplitude = 1.0, frequency = 1.0, offset = 0.0,
    start_time = 0.0, smooth = false)

Generate smooth triangular signal for frequencies less than or equal to 25 Hz

# Parameters:

  - `frequency`: [Hz] Frequency of square wave
  - `amplitude`: Amplitude of square wave
  - `offset`: Offset of output signal.
  - `start_time`: [s] Output `y = offset` for `t < start_time`
  - `smooth`:  If `true`, returns a smooth wave. Defaults to `false`
    It uses a default smoothing factor of `δ=1e-5`, but this can be changed by supplying `smooth=δ`.

# Connectors:

  - `output`
"""
@component function Triangular(; name, amplitude = 1.0, frequency = 1.0,
        offset = 0.0, start_time = 0.0, smooth = false)
    @named output = RealOutput()
    pars = @parameters begin
        amplitude = amplitude
        frequency = frequency
        offset = offset
        start_time = start_time
    end

    equation = if smooth == false
        triangular(t, frequency, amplitude, offset, start_time)
    else
        smooth === true && (smooth = 1e-5)
        smooth_triangular(t, smooth, frequency, amplitude, offset, start_time)
    end

    eqs = [
        output.u ~ equation
    ]

    compose(System(eqs, t, [], pars; name = name), [output])
end

# TODO:
# - Exponentials    Generate a rising and falling exponential signal
# - Pulse   Generate pulse signal of type Real
# - SawTooth    Generate saw tooth signal
# - Trapezoid   Generate trapezoidal signal of type Real

# SampledData Parameter struct ----------------

struct Parameter{T <: Real}
    data::Vector{T}
    ref::T
    circular_buffer::Bool
end

Parameter(data::Vector{T}, ref::T) where {T <: Real} = Parameter(data, ref, true)
Parameter(x::Parameter) = x
function Parameter(x::T; tofloat = true) where {T <: Real}
    if tofloat
        x = float(x)
        P = typeof(x)
    else
        P = T
    end

    return Parameter(P[], x)
end

function Base.isequal(x::Parameter, y::Parameter)
    b0 = length(x.data) == length(y.data)
    if b0
        b1 = all(x.data .== y.data)
        b2 = x.ref == y.ref
        return b1 & b2
    else
        return false
    end
end

Base.:*(x::Number, y::Parameter) = x * y.ref
Base.:*(y::Parameter, x::Number) = Base.:*(x, y)
Base.:*(x::Parameter, y::Parameter) = x.ref * y.ref

Base.:/(x::Number, y::Parameter) = x / y.ref
Base.:/(y::Parameter, x::Number) = y.ref / x
Base.:/(x::Parameter, y::Parameter) = x.ref / y.ref

Base.:+(x::Number, y::Parameter) = x + y.ref
Base.:+(y::Parameter, x::Number) = Base.:+(x, y)
Base.:+(x::Parameter, y::Parameter) = x.ref + y.ref

Base.:-(y::Parameter) = -y.ref
Base.:-(x::Number, y::Parameter) = x - y.ref
Base.:-(y::Parameter, x::Number) = y.ref - x
Base.:-(x::Parameter, y::Parameter) = x.ref - y.ref

Base.:^(x::Number, y::Parameter) = Base.:^(x, y.ref)
Base.:^(y::Parameter, x::Number) = Base.:^(y.ref, x)
Base.:^(x::Parameter, y::Parameter) = Base.:^(x.ref, y.ref)

Base.isless(x::Parameter, y::Number) = Base.isless(x.ref, y)
Base.isless(y::Number, x::Parameter) = Base.isless(y, x.ref)

Base.copy(x::Parameter{T}) where {T} = Parameter{T}(copy(x.data), x.ref)

ifelse(c::Bool, x::Parameter, y::Parameter) = ifelse(c, x.ref, y.ref)
ifelse(c::Bool, x::Parameter, y::Number) = ifelse(c, x.ref, y)
ifelse(c::Bool, x::Number, y::Parameter) = ifelse(c, x, y.ref)

Base.max(x::Number, y::Parameter) = max(x, y.ref)
Base.max(x::Parameter, y::Number) = max(x.ref, y)
Base.max(x::Parameter, y::Parameter) = max(x.ref, y.ref)

Base.min(x::Number, y::Parameter) = min(x, y.ref)
Base.min(x::Parameter, y::Number) = min(x.ref, y)
Base.min(x::Parameter, y::Parameter) = min(x.ref, y.ref)

function Base.show(io::IO, m::MIME"text/plain", p::Parameter)
    if !isempty(p.data)
        print(io, p.data)
    else
        print(io, p.ref)
    end
end

get_sample_time(memory::Parameter) = memory.ref
Symbolics.@register_symbolic get_sample_time(memory::Parameter)

Base.convert(::Type{T}, x::Parameter{T}) where {T <: Real} = x.ref
function Base.convert(::Type{<:Parameter{T}}, x::Number) where {T <: Real}
    Parameter{T}(T[], x, true)
end

# SampledData utilities ----------------

function linear_interpolation(x1::Real, x2::Real, t1::Real, t2::Real, t)
    if t1 != t2
        slope = (x2 - x1) / (t2 - t1)
        intercept = x1 - slope * t1

        return slope * t + intercept
    else
        @assert x1==x2 "x1 ($x1) and x2 ($x2) should be equal if t1 == t2"

        return x2
    end
end

function first_order_backwards_difference(t, memory)
    Δt = get_sample_time(memory)
    x1 = get_sampled_data(t, memory)
    x0 = get_sampled_data(t - Δt, memory)

    return (x1 - x0) / Δt
end

function first_order_backwards_difference(t, buffer, Δt, circular_buffer)
    x1 = get_sampled_data(t, buffer, Δt, circular_buffer)
    x0 = get_sampled_data(t - Δt, buffer, Δt, circular_buffer)

    return (x1 - x0) / Δt
end

function get_sampled_data(t,
        buffer::Vector{T},
        dt::T,
        circular_buffer = true) where {T <: Real}
    if t < 0
        t = zero(t)
    end

    if isempty(buffer)
        if T <: AbstractFloat
            return T(NaN)
        else
            return zero(T)
        end
    end

    i1 = floor(Int, t / dt) + 1 #expensive
    i2 = i1 + 1

    t1 = (i1 - 1) * dt
    x1 = @inbounds buffer[i1]

    if t == t1
        return x1
    else
        n = length(buffer)

        if circular_buffer
            i1 = (i1 - 1) % n + 1
            i2 = (i2 - 1) % n + 1
        else
            if i2 > n
                i2 = n
                i1 = i2 - 1
            end
        end

        t2 = (i2 - 1) * dt
        x2 = @inbounds buffer[i2]
        return linear_interpolation(x1, x2, t1, t2, t)
    end
end
function get_sampled_data(t, buffer)
    get_sampled_data(t, buffer.data, buffer.ref, buffer.circular_buffer)
end
Symbolics.@register_symbolic Parameter(data::Vector, ref, circular_buffer::Bool)
Symbolics.@register_symbolic get_sampled_data(t, buffer::Parameter)
Symbolics.@register_symbolic get_sampled_data(t, buffer::Vector, dt, circular_buffer) false

function Symbolics.derivative(::typeof(get_sampled_data), args::NTuple{2, Any}, ::Val{1})
    t = @inbounds args[1]
    buffer = @inbounds args[2]
    first_order_backwards_difference(t, buffer)
end
function ChainRulesCore.frule((_, ẋ, _), ::typeof(get_sampled_data), t, buffer)
    first_order_backwards_difference(t, buffer) * ẋ
end

function Symbolics.derivative(::typeof(get_sampled_data), args::NTuple{4, Any}, ::Val{1})
    t = @inbounds args[1]
    buffer = @inbounds args[2]
    sample_time = @inbounds args[3]
    circular_buffer = @inbounds args[4]
    first_order_backwards_difference(t, buffer, sample_time, circular_buffer)
end
function ChainRulesCore.frule((_, ẋ, _),
        ::typeof(get_sampled_data),
        t,
        buffer,
        sample_time,
        circular_buffer)
    first_order_backwards_difference(t, buffer, sample_time, circular_buffer) * ẋ
end

# SampledData component ----------------

module SampledDataType
@enum Option vector_based struct_based
end

"""
    SampledData(; name, buffer, sample_time, circular_buffer=true)

data input component.

# Parameters:
  - `buffer::Vector{Real}`:  holds the data sampled at `sample_time`
  - `sample_time::Real`
  - `circular_buffer::Bool = true`: how to handle `t > length(buffer)*sample_time`.  If true data is considered circular, otherwise last data point is held.

# Connectors:
  - `output`
"""
@component function SampledData(::Val{SampledDataType.vector_based};
        name,
        buffer,
        sample_time,
        circular_buffer = true)
    T = eltype(buffer)
    pars = @parameters begin
        buffer::Vector{T} = buffer #::Vector{Real}
        sample_time::T = sample_time #::Real
        circular_buffer::Bool = circular_buffer #::Bool
    end
    @parameters p::Parameter{T} = Parameter(buffer, sample_time, circular_buffer)
    vars = []
    systems = @named begin
        output = RealOutput()
    end
    eqs = [
        output.u ~ get_sampled_data(t, p)
    ]
    return System(eqs, t, vars, [pars; p]; name, systems)
end

"""
    SampledData(; name, buffer)

data input component.

# Parameters:
  - `buffer`: a `Parameter` type which holds the data and sample time

# Connectors:
  - `output`
"""
@component function SampledData(
        ::Val{SampledDataType.struct_based}; name, buffer::Parameter)
    pars = @parameters begin
        buffer::typeof(buffer) = buffer #::Parameter
    end
    vars = []
    systems = @named begin
        output = RealOutput()
    end
    eqs = [
        output.u ~ get_sampled_data(t, buffer)
    ]
    return System(eqs, t, vars, pars; name, systems)
end

SampledData(x::SampledDataType.Option; kwargs...) = SampledData(Val(x); kwargs...)

# struct_based
function SampledData(T::Type, circular_buffer = true; name)
    SampledData(SampledDataType.struct_based;
        name,
        buffer = Parameter(T[], zero(T), circular_buffer))
end

# vector_based
function SampledData(sample_time::T, circular_buffer = true; name) where {T <: Real}
    SampledData(SampledDataType.vector_based;
        name,
        buffer = T[],
        sample_time,
        circular_buffer)
end
function SampledData(buffer::Vector{<:Real},
        sample_time::Real,
        circular_buffer = true;
        name)
    SampledData(SampledDataType.vector_based; name, buffer, sample_time, circular_buffer)
end
function SampledData(; name, buffer, sample_time, circular_buffer)
    SampledData(SampledDataType.vector_based; name, buffer, sample_time, circular_buffer)
end

"""
    Interpolation(interp_type, u, x, args...; name)

Represent function interpolation symbolically as a block component.
By default interpolation types from [`DataInterpolations.jl`](https://github.com/SciML/DataInterpolations.jl) are supported,
but in general any callable type that builds the interpolation object via `itp = interpolation_type(u, x, args...)` and calls
the interpolation with `itp(t)` should work. This does not need to represent an interpolation, it can be any type that satisfies
the interface, such as lookup tables.
# Arguments:
  - `interp_type`: the type of the interpolation. For `DataInterpolations`,
these would be any of [the available interpolations](https://github.com/SciML/DataInterpolations.jl?tab=readme-ov-file#available-interpolations),
such as `LinearInterpolation`, `ConstantInterpolation` or `CubicSpline`.
  - `u`: the data used for interpolation. For `DataInterpolations` this will be an `AbstractVector`
  - `x`: the values that each data points correspond to, usually the times corresponding to each value in `u`.
  - `args`: any other arguments needed to build the interpolation
# Keyword arguments:
  - `name`: the name of the component

# Parameters:
  - `interpolator`: the symbolic representation of the interpolation object, callable as `interpolator(t)`

# Connectors:
  - `input`: a [`RealInput`](@ref) connector corresponding to the input variable
  - `output`: a [`RealOutput`](@ref) connector corresponding to the interpolated value
"""
function Interpolation(interp_type, u, x, args...; name)
    itp = interp_type(u, x, args...)
    Interpolation(; itp, name)
end

@deprecate Interpolation(itp; name) Interpolation(; itp, name)

function Interpolation(; itp, name)
    @parameters (interpolator::typeof(itp))(..) = itp
    @named input = RealInput()
    @named output = RealOutput()

    eqs = [output.u ~ interpolator(input.u)]

    System(
        eqs, t, [], [interpolator]; name, systems = [input, output])
end

"""
    CachedInterpolation

This callable struct caches the calls to an interpolation object via PreallocationTools.
"""
struct CachedInterpolation{T, I, U, X, C}
    interpolation_type::I
    prev_u::U
    prev_x::X
    cache::C

    function CachedInterpolation(interpolation_type, u, x, args)
        # we need to copy the inputs to avoid aliasing
        prev_u = DiffCache(copy(u))
        # Interpolation points can be a range, but we want to be able
        # to update the cache if needed (and setindex! is not defined on ranges)
        # with a view from MTKParameters, so we collect to get a vector
        prev_x = DiffCache(collect(copy(x)))
        cache = GeneralLazyBufferCache() do (u, x)
            interpolation_type(get_tmp(prev_u, u), get_tmp(prev_x, x), args...)
        end
        T = typeof(cache[(get_tmp(prev_u, u), get_tmp(prev_x, x))])
        I = typeof(interpolation_type)
        U = typeof(prev_u)
        X = typeof(prev_x)
        C = typeof(cache)

        new{T, I, U, X, C}(interpolation_type, prev_u, prev_x, cache)
    end
end

function (f::CachedInterpolation{T})(u, x, args) where {T}
    (; prev_u, prev_x, cache, interpolation_type) = f

    interp = @inbounds if (u, x) ≠ (get_tmp(prev_u, u), get_tmp(prev_x, x))
        get_tmp(prev_u, u) .= u
        get_tmp(prev_x, x) .= x
        cache.bufs[(u, x)] = interpolation_type(
            get_tmp(prev_u, u), get_tmp(prev_x, x), args...)
    else
        cache[(u, x)]
    end

    return interp
end

Base.nameof(::CachedInterpolation) = :CachedInterpolation

@register_symbolic (f::CachedInterpolation)(u::AbstractArray, x::AbstractArray, args::Tuple)

"""
    ParametrizedInterpolation(interp_type, u, x, args...; name, t = ModelingToolkit.t_nounits)

Represent function interpolation symbolically as a block component, with the interpolation data represented parametrically.
By default interpolation types from [`DataInterpolations.jl`](https://github.com/SciML/DataInterpolations.jl) are supported,
but in general any callable type that builds the interpolation object via `itp = interpolation_type(u, x, args...)` and calls
the interpolation with `itp(t)` should work. This does not need to represent an interpolation, it can be any type that satisfies
the interface, such as lookup tables.
# Arguments:
  - `interp_type`: the type of the interpolation. For `DataInterpolations`,
these would be any of [the available interpolations](https://github.com/SciML/DataInterpolations.jl?tab=readme-ov-file#available-interpolations),
such as `LinearInterpolation`, `ConstantInterpolation` or `CubicSpline`.
  - `u`: the data used for interpolation. For `DataInterpolations` this will be an `AbstractVector`
  - `x`: the values that each data points correspond to, usually the times corresponding to each value in `u`.
  - `args`: any other arguments beeded to build the interpolation
# Keyword arguments:
  - `name`: the name of the component

# Parameters:
  - `data`: the symbolic representation of the data passed at construction time via `u`.
  - `ts`: the symbolic representation of times corresponding to the data passed at construction time via `x`.

# Connectors:
  - `input`: a [`RealInput`](@ref) connector corresponding to the independent variable
  - `output`: a [`RealOutput`](@ref) connector corresponding to the interpolated value
"""
function ParametrizedInterpolation(
        interp_type::T, u::AbstractVector, x::AbstractVector, args...;
        name) where {T}
    build_interpolation = CachedInterpolation(interp_type, u, x, args)

    @parameters data[1:length(x)] = u
    @parameters ts[1:length(x)] = x
    @parameters interpolation_type::T=interp_type [tunable = false]
    @parameters (interpolator::interp_type)(..)::eltype(u)

    @named input = RealInput()
    @named output = RealOutput()

    eqs = [output.u ~ interpolator(input.u)]

    System(eqs, ModelingToolkit.t_nounits, [],
        [data, ts, interpolation_type, interpolator];
        parameter_dependencies = [
            interpolator ~ build_interpolation(data, ts, args)
        ],
        systems = [input, output],
        name)
end

function ParametrizedInterpolation(; interp_type, u::AbstractVector, x::AbstractVector, name)
    ParametrizedInterpolation(interp_type, u, x; name)
end
