"""
    $(TYPEDSIGNATURES)

Generates an DiscreteProblem from an DiscreteSystem.
"""
function DiffEqBase.DiscreteProblem(sys::DiscreteSystem,u0map,tspan,parammap=DiffEqBase.NullParameters();kwargs...)
    dvs = states(sys)
    ps = parameters(sys)
    eqs = equations(sys)
    # defs = defaults(sys)
    t = get_iv(sys)
    u0 = varmap_to_vars(u0map,dvs)
    rhss = [eq.rhs for eq in eqs]
    u = dvs
    p = varmap_to_vars(parammap,ps)
    eval_module = @__MODULE__
    eval_expression = true
    f_gen = build_function(rhss, dvs, ps, t; expression=Val{eval_expression}, expression_module=eval_module)
    f_oop,f_iip = (@RuntimeGeneratedFunction(eval_module, ex) for ex in f_gen)
    f(u,p,t) = f_oop(u,p,t)
    DiscreteProblem(f,u0,tspan,p)
end
