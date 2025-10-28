struct LogicTable{N} <: AbstractArray{Logic, N}
    logic::Array{Logic, N}
    function LogicTable(l::Array{Logic})
        any(i -> i != 9, size(l)) &&
            throw(ArgumentError("Incorrect number of logic values are passed. A variable of type
`LogicTable` must have nine entries corresponding to 9 logic levels"))
        new{ndims(l)}(l)
    end
end

Base.size(l::LogicTable) = size(l.logic)

Base.axes(l::LogicTable) = axes(l.logic)

getindex(s::LogicTable, l::Logic) = getindex(s.logic, get_logic_level(l))
function Base.getindex(s::LogicTable, i1::Logic, i2::Logic)
    getindex(s.logic, get_logic_level(i1), get_logic_level(i2))
end
function Base.getindex(s::LogicTable, i1::Logic, i2::Logic,
        I::Logic...)
    getindex(s.logic, get_logic_level(i1), get_logic_level(i2), get_logic_level(I...)...)
end

getindex(s::LogicTable, l::Int) = getindex(s.logic, l)
function getindex(s::LogicTable, i1::Int, i2::Int, I::Int...)
    getindex(s.logic, i1, i2, I...)
end

function Base.setindex!(A::LogicTable, x::Logic, i1::Int)
    setindex!(A.logic, x, i1)
end
function Base.setindex!(A::LogicTable, x::Logic, i1::Int, i2::Int, I::Int...)
    setindex!(A.logic, x, i1, i2, I...)
end

get_logic_level(l::LogicTable) = Int.(l.logic)

# AND gate
const AndTable = LogicTable([
                             # U X F0 F1 Z W L H DC
                             U U F0 U U U F0 U U        # U
                             U X F0 X X X F0 X X        # X
                             F0 F0 F0 F0 F0 F0 F0 F0 F0 # F0
                             U X F0 F1 X X F0 F1 X      # F1
                             U X F0 X X X F0 X X        # Z
                             U X F0 X X X F0 X X        # W
                             F0 F0 F0 F0 F0 F0 F0 F0 F0 # L
                             U X F0 F1 X X F0 F1 X      # H
                             U X F0 X X X F0 X X])       # DC

function _and2(a::Logic, b::Logic)
    AndTable[a, b]
end
_and2(a::Number, b::Logic) = _and2(convert(Logic, a), b)
_and2(a::Logic, b::Number) = _and2(a, convert(Logic, b))
_and2(a::Number, b::Number) = _and2(convert(Logic, a), convert(Logic, b))

function _and(x...)
    y = x[1]
    for i in 2:lastindex(x)
        y = _and2(y, x[i])
    end
    return y
end

@register_symbolic _and(a, b)

# NOT gate
const NotTable = LogicTable([U, X, F1, F0, X, X, F1, F0, X])

_not(x::Logic) = NotTable[x]
_not(x::Number) = _not(convert(Logic, x))

@register_symbolic _not(x)

# OR gate
const OrTable = LogicTable([
                            # U X F0 F1 Z W L H DC
                            U U U F1 U U U F1 U        # U
                            U X X F1 X X X F1 X        # X
                            U X F0 F1 X X F0 F1 X      # F0
                            F1 F1 F1 F1 F1 F1 F1 F1 F1 # F1
                            U X X F1 X X X F1 X        # Z
                            U X X F1 X X X F1 X        # W
                            U X F0 F1 X X F0 F1 X      # L
                            F1 F1 F1 F1 F1 F1 F1 F1 F1 # H
                            U X X F1 X X X F1 X])       # DC

function _or2(a::Logic, b::Logic)
    OrTable[a, b]
end
_or2(a::Number, b::Logic) = _or2(convert(Logic, a), b)
_or2(a::Logic, b::Number) = _or2(a, convert(Logic, b))
_or2(a::Number, b::Number) = _or2(convert(Logic, a), convert(Logic, b))

function _or(x...)
    y = x[1]
    for i in 2:lastindex(x)
        y = _or2(y, x[i])
    end
    return y
end

@register_symbolic _or(a, b)

# XOR gate
const XorTable = LogicTable([
                             # U X F0 F1 Z W L H DC
                             U U U U U U U U U     # U
                             U X X X X X X X X     # X
                             U X F0 F1 X X F0 F1 X # F0
                             U X F1 F0 X X F1 F0 X # F1
                             U X X X X X X X X     # Z
                             U X X X X X X X X     # W
                             U X F0 F1 X X F0 F1 X # L
                             U X F1 F0 X X F1 F0 X # H
                             U X X X X X X X X])    # DC

function _xor2(a::Logic, b::Logic)
    XorTable[a, b]
end
_xor2(a::Number, b::Logic) = _xor2(convert(Logic, a), b)
_xor2(a::Logic, b::Number) = _xor2(a, convert(Logic, b))
_xor2(a::Number, b::Number) = _xor2(convert(Logic, a), convert(Logic, b))

function _xor(x...)
    y = x[1]
    for i in 2:lastindex(x)
        y = _xor2(y, x[i])
    end
    return y
end

@register_symbolic _xor(a, b)
