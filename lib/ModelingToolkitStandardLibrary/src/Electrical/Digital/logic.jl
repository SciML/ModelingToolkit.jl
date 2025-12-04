@enum Logic Uninitialized=1 ForcingUnknown ForcingZero ForcingOne HighImpedance WeakUnknown WeakZero WeakOne DontCare

const U = Uninitialized
const X = ForcingUnknown
const F0 = ForcingZero
const F1 = ForcingOne
const Z = HighImpedance
const W = WeakUnknown
const L = WeakZero
const H = WeakOne
const DC = DontCare

function Base.show(io::IO, ::MIME"text/plain", l::Logic)
    if Int(l) == 1
        print(io, "U")
    elseif Int(l) == 2
        print(io, "X")
    elseif Int(l) == 3
        print(io, "F0")
    elseif Int(l) == 4
        print(io, "F1")
    elseif Int(l) == 5
        print(io, "Z")
    elseif Int(l) == 6
        print(io, "W")
    elseif Int(l) == 7
        print(io, "L")
    elseif Int(l) == 8
        print(io, "H")
    elseif Int(l) == 9
        print(io, "DC")
    else
        print(io, "Invalid logic level: $l")
    end
end

Base.zero(::Logic) = F0
Base.zero(::Type{Logic}) = F0
Base.one(::Logic) = F1
Base.one(::Type{Logic}) = F1

# Helpers to convert 1 and 0 to their `Logic` counterparts
function Base.convert(l::Type{Logic}, i::Number)
    if i == zero(i)
        zero(l)
    elseif i == one(i)
        one(l)
    else
        throw("$i isn't a valid `Logic` value")
    end
end
Base.convert(l::Type{Logic}, i::Logic) = i

get_logic_level(l::Logic) = Int(l)
