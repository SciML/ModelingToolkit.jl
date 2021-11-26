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

julia> Δ = Shift(t; dt=0.01)
(::Shift) (generic function with 2 methods)
```
"""
struct Shift <: Operator
    """Fixed Shift"""
    t
    dt
    steps::Int
    Shift(t, steps=1; dt) = new(value(t), dt, steps)
end
(D::Shift)(t) = Term{symtype(t)}(D, [t])
function (D::Shift)(x::Num)
    vt = value(x)
    if vt isa Term
        op = operation(vt)
        if op isa Shift
            op.dt == D.dt || error("Multi rate shifts are currently not supported.")
            if isequal(D.t, op.t)
                arg = arguments(vt)[1]
                return Num(Shift(D.t, D.steps + op.steps, dt=D.dt)(arg)) # Add the steps together
            end
        end
    end
    Num(D(vt))
end
SymbolicUtils.promote_symtype(::Shift, t) = t

Base.show(io::IO, D::Shift) = print(io, "Shift(", D.t, ", ", D.steps, "; dt=", D.dt, ")")

Base.:(==)(D1::Shift, D2::Shift) = isequal(D1.t, D2.t) && isequal(D1.steps, D2.steps) && isequal(D1.dt, D2.dt)
Base.hash(D::Shift, u::UInt) = hash(D.steps, hash(D.dt, hash(D.t, xor(u, 0x055640d6d952f101))))

Base.:^(D::Shift, n::Integer) = Shift(D.t, D.steps*n; dt=D.dt)
Base.literal_pow(f::typeof(^), D::Shift, ::Val{n}) where n = Shift(D.t, D.steps*n; dt=D.dt)

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

julia> Δ = Sample(t; dt=0.01)
(::Sample) (generic function with 2 methods)
```
"""
struct Sample <: Operator
    t
    dt
    Sample(t; dt) = new(value(t), dt)
end
(D::Sample)(t) = Term{symtype(t)}(D, [t])
(D::Sample)(t::Num) = Num(D(value(t)))
SymbolicUtils.promote_symtype(::Sample, t) = t

Base.show(io::IO, D::Sample) = print(io, "Sample(", D.t, "; dt=", D.dt, ")")

Base.:(==)(D1::Sample, D2::Sample) = isequal(D1.t, D2.t) && isequal(D1.dt, D2.dt)
Base.hash(D::Sample, u::UInt) = hash(D.dt, hash(D.t, xor(u, 0x055640d6d952f101)))

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

julia> k = SampledTime(t, dt=0.1);

julia> x(k)
Sample(t; dt=0.1)(x(t))

julia> x(k+1)
Shift(t, 1; dt=0.1)(Sample(t; dt=0.1)(x(t)))
```
"""
struct SampledTime
    t
    dt
    steps::Int
    SampledTime(t, steps=0; dt) = new(value(t), dt, steps)
end


function (xn::Num)(k::SampledTime)
    @unpack t, dt, steps = k
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

    sample = Sample(t; dt=dt)
    if steps == 0
        return sample(xn) # x(k) needs no shift operator if the step of k is 0
    end
    z = Shift(t, steps; dt=dt) # a shift of k steps
    z(sample(xn))
end

Base.:+(k::SampledTime, i::Int) = SampledTime(k.t, k.steps + i; dt = k.dt)

