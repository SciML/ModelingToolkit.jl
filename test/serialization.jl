using ModelingToolkit, SciMLBase, Serialization, OrdinaryDiffEq

@parameters t
@variables x(t)
D = Differential(t)

@named sys = ODESystem([D(x) ~ -0.5 * x], defaults = Dict(x => 1.0))
for prob in [
    eval(ModelingToolkit.ODEProblem{false}(sys, nothing, nothing,
        SciMLBase.NullParameters())),
    eval(ModelingToolkit.ODEProblemExpr{false}(sys, nothing, nothing,
        SciMLBase.NullParameters())),
]
    _fn = tempname()

    open(_fn, "w") do f
        serialize(f, prob)
    end

    _cmd = "using ModelingToolkit, Serialization; deserialize(\"$_fn\")"

    run(`$(Base.julia_cmd()) -e $(_cmd)`)
end

include("../examples/rc_model.jl")
io = IOBuffer()
write(io, rc_model)
str = String(take!(io))
sys = include_string(@__MODULE__, str)
@test sys == flatten(rc_model) # this actually kind of works, but the variables would have different identities.

# check answer
ss = structural_simplify(rc_model)
all_obs = [o.lhs for o in observed(ss)]
prob = ODEProblem(ss, [], (0, 0.1))
sol = solve(prob, ImplicitEuler())

## Check ODESystem with Observables ----------
ss_exp = ModelingToolkit.toexpr(ss)
ss_ = eval(ss_exp)
prob_ = ODEProblem(ss_, [], (0, 0.1))
sol_ = solve(prob_, ImplicitEuler())
@test sol[all_obs] == sol_[all_obs]

## Check ODEProblemExpr with Observables -----------

# build the observable function expression
obs_exps = []
for var in all_obs
    f = ModelingToolkit.build_explicit_observed_function(ss, var; expression = true)
    sym = ModelingToolkit.getname(var) |> string
    ex = :(if name == Symbol($sym)
        return $f(u0, p, t)
    end)
    push!(obs_exps, ex)
end
# observedfun expression for ODEFunctionExpr
observedfun_exp = :(function (var, u0, p, t)
    name = ModelingToolkit.getname(var)
    $(obs_exps...)
end)

# ODEProblemExpr with observedfun_exp included
probexpr = ODEProblemExpr{true}(ss, [], (0, 0.1); observedfun_exp);
prob_obs = eval(probexpr)
sol_obs = solve(prob_obs, ImplicitEuler())

@test sol_obs[all_obs] == sol[all_obs]
