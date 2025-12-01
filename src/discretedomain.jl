using Symbolics: Operator, Num, Term, value, recursive_hasoperator

"""
    $(TYPEDSIGNATURES)

Trait to be implemented for operators which determines whether application of the operator
generates a semantically different variable or not. For example, `Differential` and `Shift`
are not transparent but `Sample` and `Hold` are. Defaults to `false` if not implemented.
"""
is_transparent_operator(x) = is_transparent_operator(typeof(x))
is_transparent_operator(::Type) = false

"""
    $(TYPEDSIGNATURES)

Trait to be implemented for operators which determines whether the operator is applied to
a time-varying quantity and results in a time-varying quantity. For example, `Initial` and
`Pre` are not time-varying since while they are applied to variables, the application
results in a non-discrete-time parameter. `Differential`, `Shift`, `Sample` and `Hold` are
all time-varying operators. All time-varying operators must implement `input_timedomain` and
`output_timedomain`.
"""
is_timevarying_operator(x) = is_timevarying_operator(typeof(x))
is_timevarying_operator(::Type{<:Symbolics.Operator}) = true
is_timevarying_operator(::Type{Initial}) = false
is_timevarying_operator(::Type{Pre}) = false
is_timevarying_operator(::Type) = false

MTKBase.ShiftIndex() = MTKBase.ShiftIndex(Inferred())

"""
    function SampleTime()

`SampleTime()` can be used in the equations of a hybrid system to represent time sampled
at the inferred clock for that equation.
"""
struct SampleTime <: Operator
    SampleTime() = SymbolicUtils.term(SampleTime, type = Real)
end
SymbolicUtils.promote_symtype(::Type{SampleTime}, ::Type{T}) where {T} = Real
SymbolicUtils.promote_shape(::Type{SampleTime}, @nospecialize(x::SU.ShapeT)) = x
Base.nameof(::SampleTime) = :SampleTime
SymbolicUtils.isbinop(::SampleTime) = false

function MTKBase.validate_operator(op::SampleTime, args, iv; context = nothing) end

# Sample

"""
$(TYPEDEF)

Represents a sample operator. A discrete-time signal is created by sampling a continuous-time signal.

# Constructors
`Sample(clock::Union{TimeDomain, InferredTimeDomain} = InferredDiscrete())`
`Sample(dt::Real)`

`Sample(x::Num)`, with a single argument, is shorthand for `Sample()(x)`.

# Fields
$(FIELDS)

# Examples

```jldoctest
julia> using Symbolics

julia> t = ModelingToolkit.t_nounits

julia> Δ = Sample(0.01)
(::Sample) (generic function with 2 methods)
```
"""
struct Sample <: Operator
    clock::Any
    Sample(clock::Union{TimeDomain, InferredTimeDomain} = InferredDiscrete()) = new(clock)
end

is_transparent_operator(::Type{Sample}) = true

function Sample(arg::Real)
    arg = unwrap(arg)
    if symbolic_type(arg) == NotSymbolic()
        Sample(Clock(arg))
    else
        Sample()(arg)
    end
end
(D::Sample)(x) = STerm(D, SArgsT((x,)); type = symtype(x), shape = SU.shape(x))
(D::Sample)(x::Num) = Num(D(value(x)))
SymbolicUtils.promote_symtype(::Sample, ::Type{T}) where {T} = T
SymbolicUtils.promote_shape(::Sample, @nospecialize(x::SU.ShapeT)) = x
Base.nameof(::Sample) = :Sample
SymbolicUtils.isbinop(::Sample) = false

Base.show(io::IO, D::Sample) = print(io, "Sample(", D.clock, ")")

Base.:(==)(D1::Sample, D2::Sample) = isequal(D1.clock, D2.clock)
Base.hash(D::Sample, u::UInt) = hash(D.clock, xor(u, 0x055640d6d952f101))

function MTKBase.validate_operator(op::Sample, args, iv; context = nothing)
    arg = unwrap(only(args))
    if !is_variable_floatingpoint(arg)
        throw(ContinuousOperatorDiscreteArgumentError(op, arg, context))
    end
    if isparameter(arg)
        throw(ArgumentError("""
        Expected argument of $op to be an unknown, found $arg which is a parameter.
        """))
    end
end

"""
    hassample(O)

Returns true if the expression or equation `O` contains [`Sample`](@ref) terms.
"""
hassample(O) = recursive_hasoperator(Sample, unwrap(O))

# Hold

"""
$(TYPEDEF)

Represents a hold operator. A continuous-time signal is produced by holding a discrete-time signal `x` with zero-order hold.

```
cont_x = Hold()(disc_x)
```
"""
struct Hold <: Operator
end

is_transparent_operator(::Type{Hold}) = true

(D::Hold)(x) = STerm(D, SArgsT((x,)); type = symtype(x), shape = SU.shape(x))
(D::Hold)(x::Number) = x
(D::Hold)(x::Num) = Num(D(value(x)))
SymbolicUtils.promote_symtype(::Hold, ::Type{T}) where {T} = T
SymbolicUtils.promote_shape(::Hold, @nospecialize(x::SU.ShapeT)) = x
Base.nameof(::Hold) = :Hold
SymbolicUtils.isbinop(::Hold) = false

Hold(x) = Hold()(x)

function MTKBase.validate_operator(op::Hold, args, iv; context = nothing)
    # TODO: maybe validate `VariableTimeDomain`?
    return nothing
end

"""
    hashold(O)

Returns true if the expression or equation `O` contains [`Hold`](@ref) terms.
"""
hashold(O) = recursive_hasoperator(Hold, unwrap(O))

const InputTimeDomainElT = Union{TimeDomain, InferredTimeDomain, IntegerSequence}

"""
    input_timedomain(op::Operator)

Return the time-domain type (`ContinuousClock()` or `InferredDiscrete()`) that `op` operates on.
Should return a tuple containing the time domain type for each argument to the operator.
"""
function input_timedomain(s::Shift, arg = nothing)
    if has_time_domain(arg)
        return InputTimeDomainElT[get_time_domain(arg)]
    end
    InputTimeDomainElT[InferredDiscrete()]
end

"""
    output_timedomain(op::Operator)

Return the time-domain type (`ContinuousClock()` or `InferredDiscrete()`) that `op` results in.
"""
function output_timedomain(s::Shift, arg = nothing)
    if has_time_domain(t, arg)
        return get_time_domain(t, arg)
    end
    InferredDiscrete()
end

input_timedomain(::Sample, _ = nothing) = InputTimeDomainElT[ContinuousClock()]
output_timedomain(s::Sample, _ = nothing) = s.clock

function input_timedomain(::Hold, arg = nothing)
    if has_time_domain(arg)
        return InputTimeDomainElT[get_time_domain(arg)]
    end
    InputTimeDomainElT[InferredDiscrete()] # the Hold accepts any discrete
end
output_timedomain(::Hold, _ = nothing) = ContinuousClock()

sampletime(op::Sample, _ = nothing) = sampletime(op.clock)
sampletime(op::ShiftIndex, _ = nothing) = sampletime(op.clock)

function output_timedomain(x)
    if isoperator(x, Operator)
        args = arguments(x)
        return output_timedomain(operation(x), if length(args) == 1
            args[]
        else
            args
        end)
    else
        throw(ArgumentError("$x of type $(typeof(x)) is not an operator expression"))
    end
end

function input_timedomain(x)
    if isoperator(x, Operator)
        args = arguments(x)
        return input_timedomain(operation(x), if length(args) == 1
            args[]
        else
            args
        end)
    else
        throw(ArgumentError("$x of type $(typeof(x)) is not an operator expression"))
    end
end

function ZeroCrossing(expr; name = gensym(), up = true, down = true, kwargs...)
    return SymbolicContinuousCallback(
        [expr ~ 0], up ? ImperativeAffect(Returns(nothing)) : nothing;
        affect_neg = down ? ImperativeAffect(Returns(nothing)) : nothing,
        kwargs..., zero_crossing_id = name)
end

function SciMLBase.Clocks.EventClock(cb::SymbolicContinuousCallback)
    return SciMLBase.Clocks.EventClock(cb.zero_crossing_id)
end

MTKBase.distribute_shift_into_operator(::Sample) = false
MTKBase.distribute_shift_into_operator(::Hold) = false
