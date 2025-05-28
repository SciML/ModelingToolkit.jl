module MTKPyomoDynamicOptExt
using ModelingToolkit
using Pyomo
using DiffEqBase
using UnPack
using NaNMath
using Setfield
const MTK = ModelingToolkit

SPECIAL_FUNCTIONS_DICT = Dict([acos => Pyomo.py_acos,
                               log1p => Pyomo.py_log1p,
                               acosh => Pyomo.py_acosh,
                               log2 => Pyomo.py_log2,
                               asin => Pyomo.py_asin,
                               tan => Pyomo.py_tan,
                               atanh => Pyomo.py_atanh,
                               cos => Pyomo.py_cos,
                               log => Pyomo.py_log,
                               sin => Pyomo.py_sin,
                               log10 => Pyomo.py_log10,
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
        new(model, U, V, tₛ, is_free_final, nothing, PyomoVar(model.dU), MODEL_SYM, T_SYM, DUMMY_SYM)
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

pysym_getproperty(s::Union{Num, Symbolics.Symbolic}, name::Symbol) = Symbolics.wrap(SymbolicUtils.term(_getproperty, Symbolics.unwrap(s), Val{name}(), type = Symbolics.Struct{PyomoVar}))
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
    expr = Symbolics.substitute(expr, SPECIAL_FUNCTIONS_DICT)

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
    @show expr
    expr = Symbolics.substitute(expr, SPECIAL_FUNCTIONS_DICT)
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
    @unpack model, model_sym, t_sym = m
    arg_f = eval(Symbolics.build_function(arg, model_sym, t_sym))
    int_sym = Symbol(:int, hash(arg))
    setproperty!(model, int_sym, dae.Integral(model.t, wrt = model.t, rule=Pyomo.pyfunc(arg_f)))
    PyomoVar(model.tₛ * model.int_sym)
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
    discretizer.apply_to(solver_m, wrt = solver_m.t, nfe = solver_m.steps - 1, scheme = Pyomo.scheme_string(dm))
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

function MTK.successful_solve(output::PyomoOutput)
    r = output.result
    ss = r.solver.status
    Main.xx[] = ss
    tc = r.solver.termination_condition
    if Bool(ss == opt.SolverStatus.ok) && (Bool(tc == opt.TerminationCondition.optimal) || Bool(tc == opt.TerminationCondition.locallyOptimal))
        return true
    else
        return false
    end
end
end
