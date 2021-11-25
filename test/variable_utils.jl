using ModelingToolkit, Test
using ModelingToolkit: value
using SymbolicUtils: <ₑ
@parameters α β δ
expr = (((1 / β - 1) + δ) / α) ^ (1 / (α - 1))
ref = sort([β, δ, α], lt = <ₑ)
sol = sort(Num.(ModelingToolkit.get_variables(expr)), lt = <ₑ)
@test all(x->x isa Num, sol[i] == ref[i] for i in 1:3)
@test all(simplify∘value, sol[i] == ref[i] for i in 1:3)

@parameters γ
s = α => γ
expr = (((1 / β - 1) + δ) / α) ^ (1 / (α - 1))
sol = ModelingToolkit.substitute(expr, s)
new = (((1 / β - 1) + δ) / γ) ^ (1 / (γ - 1))
@test iszero(sol - new)


@variables t
function UnitDelay(dt; name)
    @variables u(t)=0.0 [input=true] y(t)=0.0 [output=true]
    Dₜ = Difference(t; dt=dt, update=true)
    eqs = [
        Dₜ(y) ~ u
    ]
    DiscreteSystem(eqs, t, name=name)
end

dt = 0.1
@named int = UnitDelay(dt)
ModelingToolkit.collect_difference_variables(int) == Set(Any[@nonamespace(int.y)])