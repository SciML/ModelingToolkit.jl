module MTKPyomoDynamicOptExt
using ModelingToolkit
using PythonCall
using DiffEqBase
using UnPack
using NaNMath
const MTK = ModelingToolkit

# import pyomo
const pyomo = PythonCall.pynew()
PythonCall.pycopy!(pyomo, pyimport("pyomo.environ"))

struct PyomoDAEVar
    v::Py
end
(v::PyomoDAEVar)(t) = v.v[:, t]
getindex(v::PyomoDAEVar, i::Union{Num, Symbolic}, t::Union{Num, Symbolic}) = wrap(Term{symeltype(v)}(getindex, [v, unwrap(i), unwrap(t)]))
getindex(v::PyomoDAEVar, i::Int) = wrap(Term{symeltype(v)}(getindex, [v, unwrap(i), Colon()]))

for ff in [acos, log1p, acosh, log2, asin, tan, atanh, cos, log, sin, log10, sqrt]
    f = nameof(ff)
    @eval NaNMath.$f(x::PyomoDAEVar) = Base.$f(x)
end

const SymbolicConcreteModel = Symbolics.symstruct(ConcreteModel)

struct PyomoModel
    model::ConcreteModel
    U::PyomoDAEVar
    V::PyomoDAEVar
    tₛ::Union{Int, Py}
    is_free_final::Bool
    model_sym::SymbolicConcreteModel
    t_sym::Union{Num, BasicSymbolic}
    idx_sym::Union{Num, BasicSymbolic}

    function PyomoModel(model, U, V, tₛ, is_free_final)
        @variables MODEL_SYM::SymbolicConcreteModel IDX_SYM::Int T_SYM
        PyomoModel(model, U, V, tₛ, is_free_final, MODEL_SYM, T_SYM, INDEX_SYM)
    end
end

struct PyomoDynamicOptProblem{uType, tType, isinplace, P, F, K} <:
       AbstractDynamicOptProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    wrapped_model::ConcreteModel
    kwargs::K

    function PyomoDynamicOptProblem(f, u0, tspan, p, model, kwargs...)
        new{typeof(u0), typeof(tspan), SciMLBase.isinplace(f, 5),
            typeof(p), typeof(f), typeof(kwargs)}(f, u0, tspan, p, model, kwargs)
    end
end

function MTK.PyomoDynamicOptProblem(sys::ODESystem, u0map, tspan, pmap;
        dt = nothing, steps = nothing,
        guesses = Dict(), kwargs...)
    prob = MTK.process_DynamicOptProblem(PyomoDynamicOptProblem, PyomoModel, sys, u0map, tspan, pmap; dt, steps, guesses, kwargs...)
    prob.wrapped_model.model.dU = pyomo.DerivativeVar(prob.wrapped_model.model.U, wrt = model.t)
    MTK.add_equational_constraints!(prob.wrapped_model, sys, pmap, tspan)
    prob
end

MTK.generate_internal_model(m::Type{PyomoModel}) = pyomo.ConcreteModel()

function MTK.generate_time_variable!(m::ConcreteModel, tspan, tsteps)
    m.steps = length(tsteps)
    m.t = pyomo.ContinuousSet(initialize = collect(tsteps), bounds = tspan)
end

function MTK.generate_state_variable!(m::ConcreteModel, u0, ns, ts) 
    m.u_idxs = pyomo.RangeSet(1, ns)
    PyomoDAEVar(pyomo.Var(m.u_idxs, m.t))
end

function MTK.generate_input_variable!(m::ConcreteModel, u0, nc, ts)
    m.v_idxs = pyomo.RangeSet(1, nc)
    PyomoDAEVar(pyomo.Var(m.v_idxs, m.t))
end

function MTK.generate_timescale(m::ConcreteModel, guess, is_free_t)
    m.tₛ = is_free_t ? pyomo.Var(initialize = guess, bounds = (0, Inf)) : 1
end

function MTK.add_constraint!(pmodel::PyomoModel, cons)
    @unpack model, model_sym, idx_sym, t_sym = pmodel
    expr = if cons isa Equation
        cons.lhs - cons.rhs == 0
    elseif cons.relational_op === Symbolics.geq
        cons.lhs - cons.rhs ≥ 0
    else
        cons.lhs - cons.rhs ≤ 0
    end
    constraint_f = Symbolics.build_function(expr, model_sym, idx_sym, t_sym)
    pyomo.Constraint(rule = constraint_f)
end

function MTK.set_objective!(m::PyomoModel, expr) = pyomo.Objective(expr = expr)

function add_initial_constraints!(model::PyomoModel, u0, u0_idxs)
    for i in u0_idxs
        model.U[i, 0].fix(u0[i])
    end
end

function fixed_t_map(m::PyomoModel, sys, exprs)
    stidxmap = Dict([v => i for (i, v) in x_ops])
    ctidxmap = Dict([v => i for (i, v) in c_ops])
    mU = Symbolics.symbolic_getproperty(model_sym, :U)
    mV = Symbolics.symbolic_getproperty(model_sym, :V)
    fixed_t_map = Dict()
    for expr in exprs
        vars = MTK.vars(exprs)
        for st in vars
            MTK.iscall(st) || continue
            x = MTK.operation(st)
            t = only(MTK.arguments(st))
            MTK.symbolic_type(t) === MTK.NotSymbolic() || continue
            m.model.t.add(t)
            fixed_t_map[x(t)] = haskey(stidxmap, x) ? mU[stidxmap[x], t] : mV[ctidxmap[x], t]
        end
    end
    fixed_t_map
end

function MTK.free_t_map(model, x_ops, c_ops)
    mU = Symbolics.symbolic_getproperty(model_sym, :U)
    mV = Symbolics.symbolic_getproperty(model_sym, :V)
    Dict([[x(tₛ) => mU[i, end] for (i, x) in enumerate(x_ops)];
         [c(tₛ) => mV[i, end] for (i, c) in enumerate(c_ops)]])
end

function MTK.whole_t_map(model)
    mU = Symbolics.symbolic_getproperty(model_sym, :U)
    mV = Symbolics.symbolic_getproperty(model_sym, :V)
    Dict([[v => mU[i, t_sym] for (i, v) in enumerate(sts)];
         [v => mV[i, t_sym] for (i, v) in enumerate(cts)]])
end

function MTK.lowered_integral(m::PyomoModel, arg)
    @unpack model, model_sym, t_sym = m
    arg_f = Symbolics.build_function(arg, model_sym, t_sym)
    Integral(model.t, wrt = model.t, rule=arg_f)
end

MTK.process_integral_bounds(model, integral_span, tspan) = integral_span

function MTK.lowered_derivative(m::PyomoModel, i)
    mdU = Symbolics.symbolic_getproperty(model_sym, :dU)
    mdU[i, t_sym]
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
    ncp = is_finite_difference(dm) ? 1 : dm.np
    discretizer.apply_to(model, wrt = m.t, nfe = m.steps, ncp = ncp, scheme = scheme_string(dm))
    solver = SolverFactory(string(collocation.solver))
    solver.solve(m)
end

function MTK.get_U_values(m::PyomoModel)
    [pyomo.value(model.U[i]) for i in model.U]
end

function MTK.get_V_values(m::PyomoModel)
    [pyomo.value(model.V[i]) for i in model.V]
end

function MTK.get_t_values(m::PyomoModel)
    [pyomo.value(model.t[i]) for i in model.t]
end

function MTK.successful_solve(m::PyomoModel)
    ss = m.solver.status
    tc = m.solver.termination_condition
    if ss == opt.SolverStatus.ok && (tc == opt.TerminationStatus.optimal || tc == opt.TerminationStatus.locallyOptimal)
        return true
    else
        return false
    end
end
end
