module MTKPyomoDynamicOptExt
using ModelingToolkit
using Pyomo
using DiffEqBase
using UnPack
using NaNMath
const MTK = ModelingToolkit

struct PyomoDynamicOptModel
    model::ConcreteModel
    U::PyomoVar
    V::PyomoVar
    tₛ::Union{Int, PyomoVar}
    is_free_final::Bool
    dU::PyomoVar
    model_sym::Union{Num, Symbolics.BasicSymbolic}
    t_sym::Union{Num, Symbolics.BasicSymbolic}
    idx_sym::Union{Num, Symbolics.BasicSymbolic}

    function PyomoDynamicOptModel(model, U, V, tₛ, is_free_final)
        @variables MODEL_SYM::Symbolics.symstruct(PyomoDynamicOptModel) IDX_SYM::Int T_SYM
        model.dU = dae.DerivativeVar(U, wrt = model.t, initialize = 0)
        new(model, U, V, tₛ, is_free_final, PyomoVar(model.dU), MODEL_SYM, T_SYM, IDX_SYM)
    end
end

struct PyomoDynamicOptProblem{uType, tType, isinplace, P, F, K} <:
       AbstractDynamicOptProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    wrapped_model::PyomoDynamicOptModel
    kwargs::K

    function PyomoDynamicOptProblem(f, u0, tspan, p, model, kwargs...)
        new{typeof(u0), typeof(tspan), SciMLBase.isinplace(f, 5),
            typeof(p), typeof(f), typeof(kwargs)}(f, u0, tspan, p, model, kwargs)
    end
end

function MTK.PyomoDynamicOptProblem(sys::ODESystem, u0map, tspan, pmap;
        dt = nothing, steps = nothing,
        guesses = Dict(), kwargs...)
    prob = MTK.process_DynamicOptProblem(PyomoDynamicOptProblem, PyomoDynamicOptModel, sys, u0map, tspan, pmap; dt, steps, guesses, kwargs...)
    conc_model = prob.wrapped_model.model
    MTK.add_equational_constraints!(prob.wrapped_model, sys, pmap, tspan)
    prob
end

MTK.generate_internal_model(m::Type{PyomoDynamicOptModel}) = ConcreteModel(pyomo.ConcreteModel())

function MTK.generate_time_variable!(m::ConcreteModel, tspan, tsteps)
    m.steps = length(tsteps)
    m.t = dae.ContinuousSet(initialize = tsteps, bounds = tspan)
end

function MTK.generate_state_variable!(m::ConcreteModel, u0, ns, ts) 
    m.u_idxs = pyomo.RangeSet(1, ns)
    m.U = pyomo.Var(m.u_idxs, m.t, initialize = 0)
    PyomoVar(m.U)
end

function MTK.generate_input_variable!(m::ConcreteModel, c0, nc, ts)
    m.v_idxs = pyomo.RangeSet(1, nc)
    m.V = pyomo.Var(m.v_idxs, m.t, initialize = 0)
    PyomoVar(m.V)
end

function MTK.generate_timescale!(m::ConcreteModel, guess, is_free_t)
    m.tₛ = is_free_t ? PyomoVar(pyomo.Var(initialize = guess, bounds = (0, Inf))) : 1
end

function MTK.add_constraint!(pmodel::PyomoDynamicOptModel, cons)
    @unpack model, model_sym, idx_sym, t_sym = pmodel
    @show model.dU
    expr = if cons isa Equation
        cons.lhs - cons.rhs == 0
    elseif cons.relational_op === Symbolics.geq
        cons.lhs - cons.rhs ≥ 0
    else
        cons.lhs - cons.rhs ≤ 0
    end
    constraint_f = Symbolics.build_function(expr, model_sym, idx_sym, t_sym, expression = Val{false})
    @show typeof(constraint_f)
    @show typeof(Pyomo.pyfunc(constraint_f))
    cons_sym = gensym()
    setproperty!(model, cons_sym, pyomo.Constraint(model.u_idxs, model.t, rule = Pyomo.pyfunc(constraint_f)))
end

function MTK.set_objective!(m::PyomoDynamicOptModel, expr)
    m.model.obj = pyomo.Objective(expr = expr)
end

function MTK.add_initial_constraints!(model::PyomoDynamicOptModel, u0, u0_idxs, ts)
    for i in u0_idxs
        model.U[i, 0].fix(u0[i])
    end
end

function MTK.lowered_integral(m::PyomoDynamicOptModel, arg, lo, hi)
    @unpack model, model_sym, t_sym = m
    arg_f = Symbolics.build_function(arg, model_sym, t_sym)
    dae.Integral(model.t, wrt = model.t, rule=arg_f)
end

MTK.process_integral_bounds(model, integral_span, tspan) = integral_span

function MTK.lowered_derivative(m::PyomoDynamicOptModel, i)
    mdU = Symbolics.symbolic_getproperty(m.model_sym, :dU).val
    Symbolics.unwrap(mdU[i, m.t_sym])
end

function MTK.lowered_var(m::PyomoDynamicOptModel, uv, i, t)
    X = Symbolics.symbolic_getproperty(m.model_sym, uv).val
    var = t isa Union{Num, Symbolics.Symbolic} ? X[i, m.t_sym] : X[i, t]
    Symbolics.unwrap(var)
end

struct PyomoCollocation <: AbstractCollocation
    solver::Union{String, Symbol}
    derivative_method::Pyomo.DiscretizationMethod
end

MTK.PyomoCollocation(solver, derivative_method = LagrangeRadau(5)) = PyomoCollocation(solver, derivative_method)

function MTK.prepare_and_optimize!(prob::PyomoDynamicOptProblem, collocation; verbose, kwargs...)
    m = prob.wrapped_model.model
    dm = collocation.derivative_method
    discretizer = TransformationFactory(dm)
    ncp = Pyomo.is_finite_difference(dm) ? 1 : dm.np
    discretizer.apply_to(m, wrt = m.t, nfe = m.steps, scheme = Pyomo.scheme_string(dm))
    solver = SolverFactory(string(collocation.solver))
    solver.solve(m, tee = true)
    Main.xx[] = solver
end

MTK.get_U_values(m::PyomoDynamicOptModel) = [pyomo.value(m.model.U[i]) for i in m.model.U.index_set()]
MTK.get_V_values(m::PyomoDynamicOptModel) = [pyomo.value(m.model.V[i]) for i in m.model.V.index_set()]
MTK.get_t_values(m::PyomoDynamicOptModel) = Pyomo.get_results(m.model, :t)

function MTK.successful_solve(m::PyomoDynamicOptModel)
    ss = m.solver.status
    tc = m.solver.termination_condition
    if ss == opt.SolverStatus.ok && (tc == opt.TerminationStatus.optimal || tc == opt.TerminationStatus.locallyOptimal)
        return true
    else
        return false
    end
end
end
