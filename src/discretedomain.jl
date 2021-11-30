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
(D::Shift)(t) = Term{symtype(t)}(D, [t])
function (D::Shift)(x::Num)
    vt = value(x)
    if vt isa Term
        op = operation(vt)
        if op isa Shift
            if isequal(D.t, op.t)
                arg = arguments(vt)[1]
                return Num(Shift(D.t, D.steps + op.steps)(arg)) # Add the steps together
            end
        end
    end
    Num(D(vt))
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
    t
    clock
    Sample(t, clock::TimeDomain) = new(value(t), clock)
    Sample(t, dt::Real) = new(value(t), Clock(t, dt))
end
(D::Sample)(t) = Term{symtype(t)}(D, [t])
(D::Sample)(t::Num) = Num(D(value(t)))
SymbolicUtils.promote_symtype(::Sample, t) = t

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

julia> k = SampledTime(t, clock=0.1);

julia> x(k)
Sample(t; dt=0.1)(x(t))

julia> x(k+1)
Shift(t, 1; dt=0.1)(Sample(t; dt=0.1)(x(t)))
```
"""
struct SampledTime
    t
    clock
    steps::Int
    SampledTime(t, clock=Inferred(), steps=0) = new(value(t), clock, steps)
    SampledTime(t, dt::Real, steps=0) = new(value(t), Clock(t, dt), steps)
end


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

    d = propagate_time_domain(xn)
    if d != clock # this is only required if the variable has another clock
        xn = Sample(t, clock)(xn) 
    end
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
input_timedomain(s::Shift) = InferredDiscrete()

"""
    output_timedomain(op::Operator)

Return the time-domain type (`Continuous()` or `Discrete()`) that `op` results in. 
"""
output_timedomain(s::Shift) = InferredDiscrete()

input_timedomain(::Sample) = Inferred() # TODO: change name to inferred, because this operator can be used on discrete variables as well
output_timedomain(s::Sample) = s.clock

input_timedomain(h::Hold) = InferredDiscrete() # the Hold accepts any discrete
output_timedomain(::Hold) = Continuous()

sampletime(op::Sample) = sampletime(op.clock)
sampletime(op::SampledTime) = sampletime(op.clock)