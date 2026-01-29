using Symbolics: Operator, Num, Term, value, recursive_hasoperator

# Shift

"""
$(TYPEDEF)

Represents a shift operator.

# Fields
$(FIELDS)

# Examples

```jldoctest
julia> using Symbolics

julia> Δ = Shift(t)
(::Shift) (generic function with 2 methods)
```
"""
struct Shift <: Operator
    """Fixed Shift"""
    t::Union{Nothing, SymbolicT}
    steps::Int
    Shift(t, steps = 1) = new(unwrap(t), steps)
end
Shift(steps::Int) = new(nothing, steps)
normalize_to_differential(s::Shift) = Differential(s.t)^s.steps
Base.nameof(::Shift) = :Shift
SymbolicUtils.isbinop(::Shift) = false

function (D::Shift)(x::Equation, allow_zero = false)
    return D(x.lhs, allow_zero) ~ D(x.rhs, allow_zero)
end
function (D::Shift)(x, allow_zero = false)
    !allow_zero && D.steps == 0 && return x
    return term(D, x; type = symtype(x), shape = SU.shape(x))
end
function (D::Shift)(x::Union{Num, Symbolics.Arr}, allow_zero = false)
    !allow_zero && D.steps == 0 && return x
    vt = value(x)
    if iscall(vt)
        op = operation(vt)
        if op isa Shift
            if D.t === nothing || isequal(D.t, op.t)
                arg = arguments(vt)[1]
                newsteps = D.steps + op.steps
                return wrap(newsteps == 0 ? arg : Shift(D.t, newsteps)(arg))
            end
        end
    end
    return wrap(D(vt, allow_zero))
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
    return op.steps <= 0 || error(
        """
        Only non-positive shifts are allowed. Found shift of $(op.steps) in $context.
        """
    )
end

hasshift(eq::Equation) = hasshift(eq.lhs) || hasshift(eq.rhs)

"""
    hasshift(O)

Returns true if the expression or equation `O` contains [`Shift`](@ref) terms.
"""
hasshift(O) = recursive_hasoperator(Shift, O)

# ShiftIndex

struct IntegerSequence end

"""
    ShiftIndex

The `ShiftIndex` operator allows you to index a signal and obtain a shifted discrete-time signal. If the signal is continuous-time, the signal is sampled before shifting.

# Examples

```
julia> t = ModelingToolkitBase.t_nounits;

julia> @variables x(t);

julia> k = ShiftIndex(t, 0.1);

julia> x(k)      # no shift
x(t)

julia> x(k+1)    # shift
Shift(1)(x(t))
```
"""
struct ShiftIndex
    clock::Any
    steps::Int
    ShiftIndex(clock, steps::Int = 0) = new(clock, steps)
    ShiftIndex(dt::Real, steps::Int = 0) = new(Clock(dt), steps)
    ShiftIndex(::Num, steps::Int) = new(IntegerSequence(), steps)
end

function (xn::Num)(k::ShiftIndex)
    @unpack clock, steps = k
    x = unwrap(xn)
    # Verify that the independent variables of k and x match and that the expression doesn't have multiple variables
    vars = Set{SymbolicT}()
    SU.search_variables!(vars, x)
    if length(vars) != 1
        error(lazy"Cannot shift a multivariate expression $x. Either create a new unknown and shift this, or shift the individual variables in the expression.")
    end
    var = only(vars)
    if operation(var) === getindex
        var = arguments(var)[1]
    end
    if !iscall(var)
        throw(ArgumentError(lazy"Cannot shift time-independent variable $var"))
    end
    if length(arguments(var)) != 1
        error(lazy"Cannot shift an expression with multiple independent variables $x.")
    end
    t = only(arguments(var))

    # d, _ = propagate_time_domain(xn)
    # if d != clock # this is only required if the variable has another clock
    #     xn = Sample(t, clock)(xn)
    # end
    # QUESTION: should we return a variable with time domain set to k.clock?
    xn = setmetadata(xn, VariableTimeDomain, k.clock)
    if steps == 0
        return xn # x(k) needs no shift operator if the step of k is 0
    end
    return Shift(t, steps)(xn) # a shift of k steps
end

function (xn::Symbolics.Arr)(k::ShiftIndex)
    @unpack clock, steps = k
    x = unwrap(xn)
    # Verify that the independent variables of k and x match and that the expression doesn't have multiple variables
    vars = Set{SymbolicT}()
    SU.search_variables!(vars, x)
    if length(vars) != 1
        error(lazy"Cannot shift a multivariate expression $x. Either create a new unknown and shift this, or shift the individual variables in the expression.")
    end
    var = only(vars)
    if !iscall(var)
        throw(ArgumentError(lazy"Cannot shift time-independent variable $var"))
    end
    if length(arguments(var)) != 1
        error(lazy"Cannot shift an expression with multiple independent variables $x.")
    end
    t = only(arguments(var))

    # d, _ = propagate_time_domain(xn)
    # if d != clock # this is only required if the variable has another clock
    #     xn = Sample(t, clock)(xn)
    # end
    # QUESTION: should we return a variable with time domain set to k.clock?
    xn = wrap(setmetadata(unwrap(xn), VariableTimeDomain, k.clock))
    if steps == 0
        return xn # x(k) needs no shift operator if the step of k is 0
    end
    return Shift(t, steps)(xn) # a shift of k steps
end

Base.:+(k::ShiftIndex, i::Int) = ShiftIndex(k.clock, k.steps + i)
Base.:-(k::ShiftIndex, i::Int) = k + (-i)

# SampleTime

"""
    function SampleTime()

`SampleTime()` can be used in the equations of a hybrid system to represent time sampled
at the inferred clock for that equation.
"""
struct SampleTime <: Operator
    SampleTime() = SymbolicUtils.term(SampleTime, type = Real)
end
SymbolicUtils.promote_symtype(::Type{SampleTime}) = Real
SymbolicUtils.promote_shape(::Type{SampleTime}) = SU.ShapeVecT()
Base.nameof(::SampleTime) = :SampleTime
SymbolicUtils.isbinop(::SampleTime) = false

function validate_operator(op::SampleTime, args, iv; context = nothing) end

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
end

is_transparent_operator(::Type{Sample}) = true

Sample(x::Num) = Sample()(unwrap(x))
Sample(arg::Real) = Sample(Clock(arg))
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
    return if isparameter(arg)
        throw(
            ArgumentError(
                """
                Expected argument of $op to be an unknown, found $arg which is a parameter.
                """
            )
        )
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
struct Hold <: Operator end

(D::Hold)(x) = STerm(D, SArgsT((x,)); type = symtype(x), shape = SU.shape(x))
(D::Hold)(x::Number) = x
(D::Hold)(x::Num) = Num(D(value(x)))
SymbolicUtils.promote_symtype(::Hold, ::Type{T}) where {T} = T
SymbolicUtils.promote_shape(::Hold, @nospecialize(x::SU.ShapeT)) = x
Base.nameof(::Hold) = :Hold
SymbolicUtils.isbinop(::Hold) = false

Hold(x) = Hold()(x)

function validate_operator(op::Hold, args, iv; context = nothing)
    # TODO: maybe validate `VariableTimeDomain`?
    return nothing
end

"""
    hashold(O)

Returns true if the expression or equation `O` contains [`Hold`](@ref) terms.
"""
hashold(O) = recursive_hasoperator(Hold, unwrap(O))

function ZeroCrossing(expr; name = gensym(), up = true, down = true, kwargs...)
    return SymbolicContinuousCallback(
        [expr ~ 0], up ? ImperativeAffect(Returns(nothing)) : nothing;
        affect_neg = down ? ImperativeAffect(Returns(nothing)) : nothing,
        kwargs..., zero_crossing_id = name
    )
end

function SciMLBase.Clocks.EventClock(cb::SymbolicContinuousCallback)
    return SciMLBase.Clocks.EventClock(cb.zero_crossing_id)
end

distribute_shift_into_operator(::Sample) = false
distribute_shift_into_operator(::Hold) = false
