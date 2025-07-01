using ModelingToolkit, SciMLBase, Serialization, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables x(t)

@named sys = System([D(x) ~ -0.5 * x], t, defaults = Dict(x => 1.0))
sys = complete(sys)
for prob in [
    eval(ModelingToolkit.ODEProblem{false}(sys, nothing, nothing)),
    eval(ModelingToolkit.ODEProblem{false}(sys, nothing, nothing; expression = Val{true}))
]
    _fn = tempname()

    open(_fn, "w") do f
        serialize(f, prob)
    end

    _cmd = "using ModelingToolkit, Serialization; deserialize(\"$_fn\")"

    run(`$(Base.julia_cmd()) -e $(_cmd)`)
end

include("common/rc_model.jl")
@unpack capacitor = rc_model
io = IOBuffer()
write(io, expand_connections(rc_model))
str = String(take!(io))

sys = include_string(@__MODULE__, str)
rc2 = expand_connections(rc_model)
@test isapprox(sys, rc2)
@test issetequal(equations(sys), equations(rc2))
@test issetequal(unknowns(sys), unknowns(rc2))
@test issetequal(parameters(sys), parameters(rc2))

# check answer
ss = mtkcompile(rc_model)
all_obs = observables(ss)
prob = ODEProblem(ss, [capacitor.v => 0.0], (0, 0.1))
sol = solve(prob, ImplicitEuler())

## Check System with Observables ----------
ss_exp = ModelingToolkit.toexpr(ss)
ss_ = complete(eval(ss_exp))
prob_ = ODEProblem(ss_, [capacitor.v => 0.0], (0, 0.1))
sol_ = solve(prob_, ImplicitEuler())
@test sol[all_obs] == sol_[all_obs]

## Check ODEProblemExpr with Observables -----------

# build the observable function expression
# ODEProblemExpr with observedfun_exp included
probexpr = ODEProblem{true}(ss, [capacitor.v => 0.0], (0, 0.1); expr = Val{true});
prob_obs = eval(probexpr)
sol_obs = solve(prob_obs, ImplicitEuler())
@show all_obs
@test sol_obs[all_obs] == sol[all_obs]
