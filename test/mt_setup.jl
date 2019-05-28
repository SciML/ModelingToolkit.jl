using ModelingToolkit
using Test

f = @mt_setup begin
    @parameters t σ ρ β
    @variables x(t) y(t) z(t)
    @derivatives D'~t

    eqs = [D(x) ~ σ*(y-x),
           D(y) ~ x*(ρ-z)-y,
           D(z) ~ x*y - β*z]
    de = ODESystem(eqs)
    ODEFunction(de, [x,y,z], [σ,ρ,β])
end

du = zeros(3)
u  = collect(1:3)
p  = collect(4:6)
f(du, u, p, 0.1)
