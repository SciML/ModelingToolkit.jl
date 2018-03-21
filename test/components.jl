using ModelingToolkit
using Base.Test

struct Lorenz <: AbstractComponent
  x::Variable
  y::Variable
  z::Variable
  σ::Variable
  ρ::Variable
  β::Variable
  eqs::Vector{Expression}
end
function generate_lorenz_eqs(t,x,y,z,σ,ρ,β)
  D = Differential(t)
  [D*x ~ σ*(y-x)
   D*y ~ x*(ρ-z)-y
   D*z ~ x*y - β*z]
end
Lorenz(t) = Lorenz(first(@DVar(x(t))),first(@DVar(y(t))),first(@DVar(z(t))),first(@Param(σ)),first(@Param(ρ)),first(@Param(β)),generate_lorenz_eqs(t,x,y,z,σ,ρ,β))

@IVar t
lz1 = Lorenz(t)
lz2 = Lorenz(t)
Expression[lz1.x ~ lz2.x
           lz1
           lz2]
