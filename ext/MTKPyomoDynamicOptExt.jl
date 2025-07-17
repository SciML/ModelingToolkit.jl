module MTKPyomoDynamicOptExt
using ModelingToolkit
using Pyomo
using DiffEqBase
using UnPack
using NaNMath
using Setfield
const MTK = ModelingToolkit

const SPECIAL_FUNCTIONS_DICT = Dict([acos => Pyomo.py_acos,
    acosh => Pyomo.py_acosh,
    asin => Pyomo.py_asin,
    tan => Pyomo.py_tan,
    atanh => Pyomo.py_atanh,
    cos => Pyomo.py_cos,
    log => Pyomo.py_log,
    sin => Pyomo.py_sin,
    sqrt => Pyomo.py_sqrt,
    exp => Pyomo.py_exp])

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
    dummy_sym::Union{Num, Symbolics.BasicSymbolic}

    function PyomoDynamicOptModel(model, U, V, tₛ, is_free_final)
        @variables MODEL_SYM::Symbolics.symstruct(ConcreteModel) T_SYM DUMMY_SYM
        model.dU = dae.DerivativeVar(U, wrt = model.t, initialize = 0)
        new(model, U, V, tₛ, is_free_final, nothing,
            PyomoVar(model.dU), MODEL_SYM, T_SYM, DUMMY_SYM)
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

function pysym_getproperty(s::Union{Num, Symbolics.Symbolic}, name::Symbol)
    Symbolics.wrap(SymbolicUtils.term(
        _getproperty, Symbolics.unwrap(s), Val{name}(), type = Symbolics.Struct{PyomoVar}))
end
_getproperty(s, name::Val{fieldname}) where {fieldname} = getproperty(s, fieldname)

function MTK.PyomoDynamicOptProblem(sys::System, op, tspan;
        dt = nothing, steps = nothing,
        guesses = Dict(), kwargs...)
    prob,
    pmap = MTK.process_DynamicOptProblem(PyomoDynamicOptProblem, PyomoDynamicOptModel,
        sys, op, tspan; dt, steps, guesses, kwargs...)
    conc_model = prob.wrapped_model.model
    MTK.add_equational_constraints!(prob.wrapped_model, sys, pmap, tspan)
    prob
end

function MTK.generate_internal_model(m::Type{PyomoDynamicOptModel})
    ConcreteModel(pyomo.ConcreteModel())
end

function MTK.generate_time_variable!(m::ConcreteModel, tspan, tsteps)
    m.steps = length(tsteps)
    m.t = dae.ContinuousSet(initialize = tsteps, bounds = tspan)
    m.time = pyomo.Var(m.t)
end

function MTK.generate_state_variable!(m::ConcreteModel, u0, ns, ts)
    m.u_idxs = pyomo.RangeSet(1, ns)
    init_f = Pyomo.pyfunc((m, i, t) -> (u0[Pyomo.pyconvert(Int, i)]))
    m.U = pyomo.Var(m.u_idxs, m.t, initialize = init_f)
    PyomoVar(m.U)
end

function MTK.generate_input_variable!(m::ConcreteModel, c0, nc, ts)
    m.v_idxs = pyomo.RangeSet(1, nc)
    init_f = Pyomo.pyfunc((m, i, t) -> (c0[Pyomo.pyconvert(Int, i)]))
    m.V = pyomo.Var(m.v_idxs, m.t, initialize = init_f)
    PyomoVar(m.V)
end

function MTK.generate_timescale!(m::ConcreteModel, guess, is_free_t)
    m.tₛ = is_free_t ? pyomo.Var(initialize = guess, bounds = (0, Inf)) : Pyomo.Py(1)
    PyomoVar(m.tₛ)
end

function MTK.add_constraint!(pmodel::PyomoDynamicOptModel, cons; n_idxs = 1)
    @unpack model, model_sym, t_sym, dummy_sym = pmodel
    expr = if cons isa Equation
        cons.lhs - cons.rhs == 0
    elseif cons.relational_op === Symbolics.geq
        cons.lhs - cons.rhs ≥ 0
    else
        cons.lhs - cons.rhs ≤ 0
    end
    expr = Symbolics.substitute(
        Symbolics.unwrap(expr), SPECIAL_FUNCTIONS_DICT, fold = false)

    cons_sym = Symbol("cons", hash(cons))
    if occursin(Symbolics.unwrap(t_sym), expr)
        f = eval(Symbolics.build_function(expr, model_sym, t_sym))
        setproperty!(model, cons_sym, pyomo.Constraint(model.t, rule = Pyomo.pyfunc(f)))
    else
        f = eval(Symbolics.build_function(expr, model_sym, dummy_sym))
        setproperty!(model, cons_sym, pyomo.Constraint(rule = Pyomo.pyfunc(f)))
    end
end

function MTK.set_objective!(pmodel::PyomoDynamicOptModel, expr)
    @unpack model, model_sym, t_sym, dummy_sym = pmodel
    expr = Symbolics.substitute(expr, SPECIAL_FUNCTIONS_DICT, fold = false)
    if occursin(Symbolics.unwrap(t_sym), expr)
        f = eval(Symbolics.build_function(expr, model_sym, t_sym))
        model.obj = pyomo.Objective(model.t, rule = Pyomo.pyfunc(f))
    else
        f = eval(Symbolics.build_function(expr, model_sym, dummy_sym))
        model.obj = pyomo.Objective(rule = Pyomo.pyfunc(f))
    end
end

function MTK.add_initial_constraints!(model::PyomoDynamicOptModel, u0, u0_idxs, ts)
    for i in u0_idxs
        model.U[i, 0].fix(u0[i])
    end
end

function MTK.lowered_integral(m::PyomoDynamicOptModel, arg, lo, hi)
    @unpack model, model_sym, t_sym, dummy_sym = m
    total = 0
    dt = Pyomo.pyconvert(Float64, (model.t.at(-1) - model.t.at(1))/(model.steps - 1))
    f = Symbolics.build_function(arg, model_sym, t_sym, expression = false)
    for (i, t) in enumerate(model.t)
        if Bool(lo < t) && Bool(t < hi)
            t_p = model.t.at(i-1)
            Δt = min(t - lo, t - t_p)
            total += 0.5*Δt*(f(model, t) + f(model, t_p))
        elseif Bool(t >= hi) && Bool(t - dt < hi)
            t_p = model.t.at(i-1)
            Δt = hi - t + dt
            total += 0.5*Δt*(f(model, t) + f(model, t_p))
        end
    end
    PyomoVar(model.tₛ * total)
end

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

function MTK.PyomoCollocation(solver, derivative_method = LagrangeRadau(5))
    PyomoCollocation(solver, derivative_method)
end

function MTK.prepare_and_optimize!(
        prob::PyomoDynamicOptProblem, collocation; verbose, kwargs...)
    solver_m = prob.wrapped_model.model.clone()
    dm = collocation.derivative_method
    discretizer = TransformationFactory(dm)
    if MTK.is_free_final(prob.wrapped_model) && !Pyomo.is_finite_difference(dm)
        error("The Lagrange-Radau and Lagrange-Legendre collocations currently cannot be used for free final problems.")
    end
    ncp = Pyomo.is_finite_difference(dm) ? 1 : dm.np
    discretizer.apply_to(solver_m, wrt = solver_m.t, nfe = solver_m.steps - 1,
        scheme = Pyomo.scheme_string(dm))

    solver = SolverFactory(string(collocation.solver))
    results = solver.solve(solver_m, tee = true)
    PyomoOutput(results, solver_m)
end

struct PyomoOutput
    result::Pyomo.Py
    model::Pyomo.Py
end

function MTK.get_U_values(output::PyomoOutput)
    m = output.model
    [[Pyomo.pyconvert(Float64, pyomo.value(m.U[i, t])) for i in m.u_idxs] for t in m.t]
end
function MTK.get_V_values(output::PyomoOutput)
    m = output.model
    [[Pyomo.pyconvert(Float64, pyomo.value(m.V[i, t])) for i in m.v_idxs] for t in m.t]
end
function MTK.get_t_values(output::PyomoOutput)
    m = output.model
    Pyomo.pyconvert(Float64, pyomo.value(m.tₛ)) * [Pyomo.pyconvert(Float64, t) for t in m.t]
end

function MTK.objective_value(output::PyomoOutput)
    Pyomo.pyconvert(Float64, pyomo.value(output.model.obj))
end

function MTK.successful_solve(output::PyomoOutput)
    r = output.result
    ss = r.solver.status
    tc = r.solver.termination_condition
    if Bool(ss == opt.SolverStatus.ok) && (Bool(tc == opt.TerminationCondition.optimal) ||
        Bool(tc == opt.TerminationCondition.locallyOptimal))
        return true
    else
        return false
    end
end
end
