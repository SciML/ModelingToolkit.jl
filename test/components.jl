using ModelingToolkit
using Test

struct Lorenz <: AbstractComponent
    x::Variable
    y::Variable
    z::Variable
    σ::Variable
    ρ::Variable
    β::Variable
    eqs::Vector{Equation}
end
function generate_lorenz_eqs(t,x,y,z,σ,ρ,β)
    D = Differential(t)
    [D(x) ~ σ*(y-x)
     D(y) ~ x*(ρ-z)-y
     D(z) ~ x*y - β*z]
end
function Lorenz(t)
    @Unknown x(t) y(t) z(t)
    @Param σ ρ β
    Lorenz(x, y, z, σ, ρ, β, generate_lorenz_eqs(t, x, y, z, σ, ρ, β))
end

@IVar t
lz1 = Lorenz(t)
lz2 = Lorenz(t)
lz1.x ~ lz2.x
