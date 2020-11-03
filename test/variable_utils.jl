using ModelingToolkit, Test
using ModelingToolkit: value
@parameters α β δ
expr = (((1 / β - 1) + δ) / α) ^ (1 / (α - 1))
ref = [β, δ, α]
sol = Num.(ModelingToolkit.get_variables(expr))
@test all(simplify, sol[i] == ref[i] for i in 1:3)

s = value(α) => value(γ)
@parameters γ
expr = (((1 / β - 1) + δ) / α) ^ (1 / (α - 1))
sol = ModelingToolkit.substitute(expr, Dict(s))
new = (((1 / β - 1) + δ) / γ) ^ (1 / (γ - 1))
@test iszero(sol - new)
