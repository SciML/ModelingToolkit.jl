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

julia> @variables t;

julia> Δ = Shift(t)
(::Shift) (generic function with 2 methods)
```
"""
struct Shift <: Operator
    """Fixed Shift"""
    t
    steps::Int
    Shift(t, steps=1) = new(value(t), steps)
end
function (D::Shift)(x, allow_zero = false)
    !allow_zero && D.steps == 0 && return x
    Term{symtype(x)}(D, [x])
end
function (D::Shift)(x::Num, allow_zero = false)
    !allow_zero && D.steps == 0 && return x
    vt = value(x)
    if vt isa Term
        op = operation(vt)
        if op isa Shift
            if isequal(D.t, op.t)
                arg = arguments(vt)[1]
                newsteps = D.steps + op.steps
                return Num(newsteps == 0 ? arg : Shift(D.t, newsteps)(arg))
            end
        end
    end
    Num(D(vt, allow_zero))
end
SymbolicUtils.promote_symtype(::Shift, t) = t

Base.show(io::IO, D::Shift) = print(io, "Shift(", D.t, ", ", D.steps, ")")

Base.:(==)(D1::Shift, D2::Shift) = isequal(D1.t, D2.t) && isequal(D1.steps, D2.steps)
Base.hash(D::Shift, u::UInt) = hash(D.steps, hash(D.t, xor(u, 0x055640d6d952f101)))

Base.:^(D::Shift, n::Integer) = Shift(D.t, D.steps*n)
Base.literal_pow(f::typeof(^), D::Shift, ::Val{n}) where n = Shift(D.t, D.steps*n)

hasshift(eq::Equation) = hasshift(eq.lhs) || hasshift(eq.rhs)


"""
    hasshift(O)

Returns true if the expression or equation `O` contains [`Shift`](@ref) terms.
"""
hasshift(O) = recursive_hasoperator(Shift, O)


# Sample

"""
$(TYPEDEF)

Represents a sample operator. A discrete-time signal is created by sampling a continuous-time signal.

# Fields
$(FIELDS)

# Examples

```jldoctest
julia> using Symbolics

julia> @variables t;

julia> Δ = Sample(t; clock=Clock(0.01))
(::Sample) (generic function with 2 methods)
```
"""
struct Sample <: Operator
    t # TODO: remove the need to keep the iv, it can be inferred from the varibale that is sampled
    clock
    Sample(t, clock::TimeDomain = InferredDiscrete()) = new(value(t), clock)
    Sample(t, dt::Real) = new(value(t), Clock(t, dt))
end
Sample(clock::AbstractClock) = Sample(clock.t, clock)
(D::Sample)(x) = Term{symtype(x)}(D, [x])
(D::Sample)(x::Num) = Num(D(value(x)))
SymbolicUtils.promote_symtype(::Sample, x) = x

Base.show(io::IO, D::Sample) = print(io, "Sample(", D.t, "; clock=", D.clock, ")")

Base.:(==)(D1::Sample, D2::Sample) = isequal(D1.t, D2.t) && isequal(D1.clock, D2.clock)
Base.hash(D::Sample, u::UInt) = hash(D.clock, hash(D.t, xor(u, 0x055640d6d952f101)))

"""
    hassample(O)

Returns true if the expression or equation `O` contains [`Sample`](@ref) terms.
"""
hassample(O) = recursive_hasoperator(Sample, O)


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
(D::Hold)(x) = Term{symtype(x)}(D, [x])
(D::Hold)(x::Num) = Num(D(value(x)))
SymbolicUtils.promote_symtype(::Hold, x) = x

Hold(x) = Hold()(x)

"""
    hashold(O)

Returns true if the expression or equation `O` contains [`Hold`](@ref) terms.
"""
hashold(O) = recursive_hasoperator(Hold, O)



# SampledTime

"""
    SampledTime

The `SampledTime` operator allows you to index a signal and obtain a shifted discrete-time signal. If the signal is continuous-time, the signal is sampled before shifting.

# Examples 
```
julia> @variables t x(t);

julia> k = SampledTime(t, 0.1);

julia> x(k)                                         # no shift, only sample
Sample(t; clock=Clock(t, 0.1))(x(t))

julia> x(k+1)                                       # sample and shift
Shift(t, 1)(Sample(t; clock=Clock(t, 0.1))(x(t)))

julia> @variables t xd(t) [timedomain=k.clock];

julia> xd(k)                                        # Discrete variables are not sampled
xd(t)

julia> xd(k-1)                                      # only shift, no sample
Shift(t, -1)(xd(t))
```
"""
struct SampledTime
    t
    clock::TimeDomain
    steps::Int
    SampledTime(t, clock=Inferred(), steps=0) = new(value(t), clock, steps)
    SampledTime(t, dt::Real, steps=0) = new(value(t), Clock(t, dt), steps)
end
SampledTime(d::AbstractClock) = SampledTime(d.t, d)


function (xn::Num)(k::SampledTime)
    @unpack t, clock, steps = k
    x = value(xn)
    t = k.t
    # Verify that the independent variables of k and x match and that the expression doesn't have multiple variables
    vars = Symbolics.get_variables(x)
    length(vars) == 1 ||
        error("Cannot sample a multivariate expression $x. Create a new state and sample this instead.")
    args = Symbolics.arguments(vars[]) # args should be one element vector with the t in x(t)
    length(args) == 1 ||
        error("Cannot sample an expression with multiple independent variables $x.")
    isequal(args[], t) ||
        error("Independent variable of $xn is not the same as that of the SampledTime $(k.t)")

    # d, _ = propagate_time_domain(xn)
    # if d != clock # this is only required if the variable has another clock
    #     xn = Sample(t, clock)(xn) 
    # end
    if steps == 0
        return xn # x(k) needs no shift operator if the step of k is 0
    end
    Shift(t, steps)(xn) # a shift of k steps
end

Base.:+(k::SampledTime, i::Int) = SampledTime(k.t, k.clock, k.steps + i)
Base.:-(k::SampledTime, i::Int) = k + (-i)




"""
    input_timedomain(op::Operator)

Return the time-domain type (`Continuous()` or `Discrete()`) that `op` operates on. 
"""
function input_timedomain(s::Shift, arg=nothing)
    if has_time_domain(arg)
        return get_time_domain(arg)
    end
    InferredDiscrete()
end

"""
    output_timedomain(op::Operator)

Return the time-domain type (`Continuous()` or `Discrete()`) that `op` results in. 
"""
function output_timedomain(s::Shift, arg=nothing)
    if has_time_domain(arg)
        return get_time_domain(arg)
    end
    InferredDiscrete()
end

input_timedomain(::Sample, arg=nothing) = Continuous()
output_timedomain(s::Sample, arg=nothing) = s.clock

function input_timedomain(h::Hold, arg=nothing)
    if has_time_domain(arg)
        return get_time_domain(arg)
    end
    InferredDiscrete() # the Hold accepts any discrete
end
output_timedomain(::Hold, arg=nothing) = Continuous()

sampletime(op::Sample, arg=nothing) = sampletime(op.clock)
sampletime(op::SampledTime, arg=nothing) = sampletime(op.clock)

changes_domain(op) = isoperator(op, Union{Sample, Hold})

function output_timedomain(x)
    if isoperator(x, Operator)
        return output_timedomain(operation(x), arguments(x)[])
    else
        throw(ArgumentError("$x of type $(typeof(x)) is not an operator expression"))
    end
end

function input_timedomain(x)
    if isoperator(x, Operator)
        return input_timedomain(operation(x), arguments(x)[])
    else
        throw(ArgumentError("$x of type $(typeof(x)) is not an operator expression"))
    end
end