using Symbolics: Operator, Num, Term, value, recursive_hasoperator
""""""
is_transparent_operator(x) = is_transparent_operator(typeof(x))
is_transparent_operator(::Type) = false
""""""
is_timevarying_operator(x) = is_timevarying_operator(typeof(x))
is_timevarying_operator(::Type{<:Symbolics.Operator}) = true
is_timevarying_operator(::Type) = false
""""""
struct SampleTime <: Operator
    SampleTime() = SymbolicUtils.term(SampleTime, type = Real)
end
SymbolicUtils.promote_symtype(::Type{SampleTime}, ::Type{T}) where {T} = Real
SymbolicUtils.promote_shape(::Type{SampleTime}, @nospecialize(x::SU.ShapeT)) = x
Base.nameof(::SampleTime) = :SampleTime
SymbolicUtils.isbinop(::SampleTime) = false
function validate_operator(op::SampleTime, args, iv; context = nothing) end
""""""
struct Shift <: Operator
    """Fixed Shift"""
    t::Union{Nothing, SymbolicT}
    steps::Int
    Shift(t, steps = 1) = new(value(t), steps)
end
Shift(steps::Int) = new(nothing, steps)
normalize_to_differential(s::Shift) = Differential(s.t)^s.steps
Base.nameof(::Shift) = :Shift
SymbolicUtils.isbinop(::Shift) = false
function (D::Shift)(x::Equation, allow_zero = false)
    D(x.lhs, allow_zero) ~ D(x.rhs, allow_zero)
end
function (D::Shift)(x, allow_zero = false)
    !allow_zero && D.steps == 0 && return x
    term(D, x; type = symtype(x), shape = SU.shape(x))
end
function (D::Shift)(x::Union{Num, Symbolics.Arr}, allow_zero = false)
    !allow_zero && D.steps == 0 && return x
    vt = value(x)
    if iscall(vt)
        op = operation(vt)
        if op isa Sample
            error("Cannot shift a `Sample`. Create a variable to represent the sampled value and shift that instead")
        elseif op isa Shift
            if D.t === nothing || isequal(D.t, op.t)
                arg = arguments(vt)[1]
                newsteps = D.steps + op.steps
                return wrap(newsteps == 0 ? arg : Shift(D.t, newsteps)(arg))
            end
        end
    end
    wrap(D(vt, allow_zero))
end
SymbolicUtils.promote_symtype(::Shift, ::Type{T}) where {T} = T
SymbolicUtils.promote_shape(::Shift, @nospecialize(x::SU.ShapeT)) = x
Base.show(io::IO, D::Shift) = print(io, "Shift(", D.t, ", ", D.steps, ")")
Base.:(==)(D1::Shift, D2::Shift) = isequal(D1.t, D2.t) && isequal(D1.steps, D2.steps)
Base.hash(D::Shift, u::UInt) = hash(D.steps, hash(D.t, xor(u, 0x055640d6d952f101)))
Base.:^(D::Shift, n::Integer) = Shift(D.t, D.steps * n)
Base.literal_pow(f::typeof(^), D::Shift, ::Val{n}) where {n} = Shift(D.t, D.steps * n)
function validate_operator(op::Shift, args, iv; context = nothing)
    isequal(op.t, iv) || throw(OperatorIndepvarMismatchError(op, iv, context))
    op.steps <= 0 || error("""
    Only non-positive shifts are allowed. Found shift of $(op.steps) in $context.
    """)
end
hasshift(eq::Equation) = hasshift(eq.lhs) || hasshift(eq.rhs)
""""""
hasshift(O) = recursive_hasoperator(Shift, O)
""""""
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
function validate_operator(op::Sample, args, iv; context = nothing)
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
""""""
hassample(O) = recursive_hasoperator(Sample, unwrap(O))
""""""
struct Hold <: Operator
end
is_transparent_operator(::Type{Hold}) = true
(D::Hold)(x) = STerm(D, SArgsT((x,)); type = symtype(x), shape = SU.shape(x))
(D::Hold)(x::Num) = Num(D(value(x)))
SymbolicUtils.promote_symtype(::Hold, ::Type{T}) where {T} = T
SymbolicUtils.promote_shape(::Hold, @nospecialize(x::SU.ShapeT)) = x
Base.nameof(::Hold) = :Hold
SymbolicUtils.isbinop(::Hold) = false
Hold(x) = Hold()(x)
function validate_operator(op::Hold, args, iv; context = nothing)
    return nothing
end
""""""
hashold(O) = recursive_hasoperator(Hold, unwrap(O))
""""""
struct ShiftIndex
    clock::Union{InferredTimeDomain, TimeDomain, IntegerSequence}
    steps::Int
    function ShiftIndex(
            clock::Union{TimeDomain, InferredTimeDomain, IntegerSequence} = Inferred(), steps::Int = 0)
        new(clock, steps)
    end
    ShiftIndex(dt::Real, steps::Int = 0) = new(Clock(dt), steps)
    ShiftIndex(::Num, steps::Int) = new(IntegerSequence(), steps)
end
function (xn::Num)(k::ShiftIndex)
    @unpack clock, steps = k
    x = unwrap(xn)
    vars = Set{SymbolicT}()
    SU.search_variables!(vars, x)
    if length(vars) != 1
        error("Cannot shift a multivariate expression $x. Either create a new unknown and shift this, or shift the individual variables in the expression.")
    end
    var = only(vars)
    if !iscall(var)
        throw(ArgumentError("Cannot shift time-independent variable $var"))
    end
    if operation(var) == getindex
        var = first(arguments(var))
    end
    if length(arguments(var)) != 1
        error("Cannot shift an expression with multiple independent variables $x.")
    end
    xn = setmetadata(xn, VariableTimeDomain, k.clock)
    if steps == 0
        return xn
    end
    Shift(t, steps)(xn)
end
function (xn::Symbolics.Arr)(k::ShiftIndex)
    @unpack clock, steps = k
    x = value(xn)
    vars = ModelingToolkit.vars(x)
    if length(vars) != 1
        error("Cannot shift a multivariate expression $x. Either create a new unknown and shift this, or shift the individual variables in the expression.")
    end
    var = only(vars)
    if !iscall(var)
        throw(ArgumentError("Cannot shift time-independent variable $var"))
    end
    if length(arguments(var)) != 1
        error("Cannot shift an expression with multiple independent variables $x.")
    end
    xn = wrap(setmetadata(unwrap(xn), VariableTimeDomain, k.clock))
    if steps == 0
        return xn
    end
    Shift(t, steps)(xn)
end
Base.:+(k::ShiftIndex, i::Int) = ShiftIndex(k.clock, k.steps + i)
Base.:-(k::ShiftIndex, i::Int) = k + (-i)
const InputTimeDomainElT = Union{TimeDomain, InferredTimeDomain}
""""""
function input_timedomain(s::Shift, arg = nothing)
    if has_time_domain(arg)
        return InputTimeDomainElT[get_time_domain(arg)]
    end
    InputTimeDomainElT[InferredDiscrete()]
end
""""""
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
    InputTimeDomainElT[InferredDiscrete()]
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
