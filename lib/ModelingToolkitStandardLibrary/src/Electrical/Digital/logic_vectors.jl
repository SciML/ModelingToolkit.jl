import Base: size, axes, getindex, setindex!

const LogicOrNumber = Union{Logic, Number}

struct StdULogicVector{N} <: AbstractArray{Logic, N}
    logic::Array{Logic, N}
    function StdULogicVector(l::Array)
        new{ndims(l)}(Array{Logic}(convert.(Logic, l)))
    end
end

struct StdLogicVector{N} <: AbstractArray{Logic, N}
    logic::Array{Logic, N}
    function StdLogicVector(l::Array)
        new{ndims(l)}(Array{Logic}(convert.(Logic, l)))
    end
end

const LogicVector = Union{StdULogicVector, StdLogicVector}

size(l::LogicVector) = size(l.logic)

axes(l::LogicVector) = axes(l.logic)

getindex(s::LogicVector, i::Int) = getindex(s.logic, i)
function Base.getindex(s::LogicVector, i1::Int, i2::Int,
        I::Int...)
    getindex(s.logic, i1, i2, I...)
end

setindex!(A::LogicVector, x::Logic, i1::Int) = setindex!(A.logic, x, i1)
function Base.setindex!(A::LogicVector, x::Logic, i1::Int, i2::Int, I::Int...)
    setindex!(A.logic, x, i1, i2, I...)
end

get_logic_level(s::LogicVector) = Int.(s.logic)

# predefined vectors
const std_ulogic = StdULogicVector([U, X, F0, F1, Z, W, L, H, DC])
const UX01 = StdULogicVector([U, X, F0, F1])
const UX01Z = StdULogicVector([U, X, F0, F1, Z])
const X01 = StdULogicVector([X, F0, F1])
const X01Z = StdULogicVector([X, F0, F1, Z])
