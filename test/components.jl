using ModelingToolkit
using Test

struct Lorenz <: AbstractComponent
    x::Variable
    y::Variable
    z::Variable
    σ::Variable
    ρ::Variable
    β::Variable
    eqs::Vector{Term}
end
function generate_lorenz_eqs(t,x,y,z,σ,ρ,β)
    D = Differential(t)
    [@term D(x) ~ σ*(y-x)
     @term D(y) ~ x*(ρ-z)-y
     @term D(z) ~ x*y - β*z]
end
function Lorenz(t)
    @DVar x(t) y(t) z(t)
    @Param σ ρ β
    Lorenz(x, y, z, σ, ρ, β, generate_lorenz_eqs(t, x, y, z, σ, ρ, β))
end

@IVar t
lz1 = Lorenz(t)
lz2 = Lorenz(t)
Term[
    @term lz1.x ~ lz2.x
    @term lz1
    @term lz2
]
