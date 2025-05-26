module MTKPyomoDynamicOptExt
using ModelingToolkit
using Pyomo
using DiffEqBase
using UnPack
using NaNMath
using Setfield
const MTK = ModelingToolkit

struct PyomoDynamicOptModel
    model::ConcreteModel
    U::PyomoVar
    V::PyomoVar
    tₛ::PyomoVar
    is_free_final::Bool
    solver_model::Union{Nothing, ConcreteModel}
    dU::PyomoVar
    model_sym::Union{Num, Symbolics.BasicSymbolic}
    t_sym::Union{Num, Symbolics.BasicSymbolic}
    uidx_sym::Union{Num, Symbolics.BasicSymbolic}
    vidx_sym::Union{Num, Symbolics.BasicSymbolic}

    function PyomoDynamicOptModel(model, U, V, tₛ, is_free_final)
        @variables MODEL_SYM::Symbolics.symstruct(ConcreteModel) U_IDX_SYM::Int V_IDX_SYM::Int T_SYM
        model.dU = dae.DerivativeVar(U, wrt = model.t, initialize = 0)
        new(model, U, V, tₛ, is_free_final, nothing, PyomoVar(model.dU), MODEL_SYM, T_SYM, U_IDX_SYM, V_IDX_SYM)
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

pysym_getproperty(s, name::Symbol) = Symbolics.wrap(SymbolicUtils.term(_getproperty, s, Val{name}(), type = Symbolics.Struct{PyomoVar}))
_getproperty(s, name::Val{fieldname}) where fieldname = getproperty(s, fieldname)

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
    m.tₛ = is_free_t ? PyomoVar(pyomo.Var(initialize = guess, bounds = (0, Inf))) : PyomoVar(Pyomo.Py(1))
end

function MTK.add_constraint!(pmodel::PyomoDynamicOptModel, cons; n_idxs = 1)
    @unpack model, model_sym, t_sym = pmodel
    expr = if cons isa Equation
        cons.lhs - cons.rhs == 0
    elseif cons.relational_op === Symbolics.geq
        cons.lhs - cons.rhs ≥ 0
    else
        cons.lhs - cons.rhs ≤ 0
    end
    f_expr = Symbolics.build_function(expr, model_sym, t_sym)
    cons_sym = Symbol("cons", hash(cons)) 
    constraint_f = eval(:(cons_sym = $f_expr))
    setproperty!(model, cons_sym, pyomo.Constraint(model.t, rule = Pyomo.pyfunc(constraint_f)))
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
    mdU = Symbolics.value(pysym_getproperty(m.model_sym, :dU))
    Symbolics.unwrap(mdU[i, m.t_sym])
end

function MTK.lowered_var(m::PyomoDynamicOptModel, uv, i, t)
    X = Symbolics.value(pysym_getproperty(m.model_sym, uv))
    var = t isa Union{Num, Symbolics.Symbolic} ? X[i, m.t_sym] : X[i, t]
    Symbolics.unwrap(var)
end

struct PyomoCollocation <: AbstractCollocation
    solver::Union{String, Symbol}
    derivative_method::Pyomo.DiscretizationMethod
end

MTK.PyomoCollocation(solver, derivative_method = LagrangeRadau(5)) = PyomoCollocation(solver, derivative_method)

function MTK.prepare_and_optimize!(prob::PyomoDynamicOptProblem, collocation; verbose, kwargs...)
    solver_m = prob.wrapped_model.model.clone()
    dm = collocation.derivative_method
    discretizer = TransformationFactory(dm)
    ncp = Pyomo.is_finite_difference(dm) ? 1 : dm.np
    discretizer.apply_to(solver_m, wrt = solver_m.t, nfe = solver_m.steps, scheme = Pyomo.scheme_string(dm))
    solver = SolverFactory(string(collocation.solver))
    solver.solve(solver_m, tee = true)
    solver_m
end

MTK.get_U_values(m::ConcreteModel) = [[pyomo.value(m.U[i, t]) for i in m.u_idxs] for t in m.t]
MTK.get_V_values(m::ConcreteModel) = [[pyomo.value(m.V[i, t]) for i in m.v_idxs] for t in m.t]
MTK.get_t_values(m::ConcreteModel) = [t for t in m.t]

function MTK.successful_solve(m::ConcreteModel)
    ss = m.solver.status
    tc = m.solver.termination_condition
    if ss == opt.SolverStatus.ok && (tc == opt.TerminationStatus.optimal || tc == opt.TerminationStatus.locallyOptimal)
        return true
    else
        return false
    end
end
end
