struct Differential{V<:AbstractVariable}
    x::V
end
function Derivative end
Base.:*(D::Differential,x::Expression) = Operation(Derivative,Expression[x,D.x])

export Differential,Derivative
