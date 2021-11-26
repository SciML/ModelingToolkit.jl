using Symbolics: Operator, Num, Term, value

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
    Shift(t; dt) = new(value(t), dt)
end
(D::Shift)(t) = Term{symtype(t)}(D, [t])
(D::Shift)(t::Num) = Num(D(value(t)))
SymbolicUtils.promote_symtype(::Shift, t) = t

Base.show(io::IO, D::Shift) = print(io, "Shift(", D.t, "; dt=", D.dt, ")")

Base.:(==)(D1::Shift, D2::Shift) = isequal(D1.t, D2.t) && isequal(D1.dt, D2.dt)
Base.hash(D::Shift, u::UInt) = hash(D.dt, hash(D.t, xor(u, 0x055640d6d952f101)))

Base.:^(D::Shift, n::Integer) = Symbolics._repeat_apply(D, n)

hasshift(eq::Equation) = hasshift(eq.lhs) || hasshift(eq.rhs)


"""
    hasshift(O)

Returns true if the expression or equation `O` contains [`Shift`](@ref) terms.
"""
function hasshift(O)
    istree(O) || return false
    if operation(O) isa Shift
        return true
    else
        if O isa Union{Add, Mul}
            any(hasshift, keys(O.dict))
        elseif O isa Pow
            hasshift(O.base) || hasshift(O.exp)
        elseif O isa SymbolicUtils.Div
            hasshift(O.num) || hasshift(O.den)
        else
            any(hasshift, arguments(O))
        end
    end
end


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


hassample(eq::Equation) = hassample(eq.lhs) || hassample(eq.rhs)
hassample(h::Sample) = true

"""
    hassample(O)

Returns true if the expression or equation `O` contains [`Sample`](@ref) terms.
"""
function hassample(O)
    istree(O) || return false
    if operation(O) isa Sample
        return true
    else
        if O isa Union{Add, Mul}
            any(hassample, keys(O.dict))
        elseif O isa Pow
            hassample(O.base) || hassample(O.exp)
        elseif O isa SymbolicUtils.Div
            hassample(O.num) || hassample(O.den)
        else
            any(hassample, arguments(O))
        end
    end
end


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

hashold(eq::Equation) = hashold(eq.lhs) || hashold(eq.rhs)
hashold(h::Hold) = true

"""
    hashold(O)

Returns true if the expression or equation `O` contains [`Hold`](@ref) terms.
"""
function hashold(O)
    istree(O) || return false
    if operation(O) isa Hold
        return true
    else
        if O isa Union{Add, Mul}
            any(hashold, keys(O.dict))
        elseif O isa Pow
            hashold(O.base) || hashold(O.exp)
        elseif O isa SymbolicUtils.Div
            hashold(O.num) || hashold(O.den)
        else
            any(hashold, arguments(O))
        end
    end
end