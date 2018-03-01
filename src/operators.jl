abstract type AbstractOperator end

struct Differential
    x::Variable
    order::Int
end
Differential(x) = Differential(x,1)

function Derivative end
Base.:*(D::Differential,x::Expression) = Operation(Derivative,Expression[x,D.x])

export Differential,Derivative,AbstractOperator
