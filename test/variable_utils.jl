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



# Continuous
using ModelingToolkit: isdifferential, isdifference, vars, difference_vars, collect_difference_variables, collect_differential_variables, collect_ivs
@variables t u(t) y(t)
D = Differential(t)
eq = D(y) ~ u
v = vars(eq)
@test v == Set([D(y), u])

ov = collect_differential_variables(eq)
@test ov == Set(Any[y])

aov = ModelingToolkit.collect_applied_operators(eq, Differential)
@test aov == Set(Any[D(y)])

ts = collect_ivs([eq])
@test ts == Set([t])




# Discrete
z = Difference(t; dt=1, update=true)
eq = z(y) ~ u
v = difference_vars(eq)
@test v == Set([z(y), u])

ov = collect_difference_variables(eq)
@test ov == Set(Any[y])

aov = ModelingToolkit.collect_applied_operators(eq, Difference)
@test aov == Set(Any[z(y)])


ts = collect_ivs([eq], Difference)
@test ts == Set([t])


@variables t
function UnitDelay(dt; name)
    @variables u(t)=0.0 [input=true] y(t)=0.0 [output=true]
    z = Difference(t; dt=dt, update=true)
    eqs = [
        z(y) ~ u
    ]
    DiscreteSystem(eqs, t, name=name)
end

dt = 0.1
@named int = UnitDelay(dt)
collect_difference_variables(int) == Set(Any[@nonamespace(int.y)])
