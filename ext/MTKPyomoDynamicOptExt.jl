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
getindex(v::PyomoDAEVar, i::Union{Num, Symbolic}, t::Union{Num, Symbolic}) = wrap(Term{symeltype(A)}(getindex, [A, unwrap(i), unwrap(t)]))

for ff in [acos, log1p, acosh, log2, asin, tan, atanh, cos, log, sin, log10, sqrt]
    f = nameof(ff)
    @eval NaNMath.$f(x::PyomoDAEVar) = Base.$f(x)
end

const SymbolicConcreteModel = Symbolics.symstruct(ConcreteModel)

struct PyomoModel
    model::ConcreteModel
    U
    V
    tₛ::Union{Int}
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

"""
    PyomoDynamicOptProblem(sys::ODESystem, u0, tspan, p; dt, steps)

Convert an ODESystem representing an optimal control system into a Pyomo model
for solving using optimization. Must provide either `dt`, the timestep between collocation 
points (which, along with the timespan, determines the number of points), or directly 
provide the number of points as `steps`.
"""
function MTK.PyomoDynamicOptProblem(sys::ODESystem, u0map, tspan, pmap;
        dt = nothing,
        steps = nothing,
        guesses = Dict(), kwargs...)
    prob = MTK.process_DynamicOptProblem(PyomoDynamicOptProblem, PyomoModel, sys, u0map, tspan, pmap; dt, steps, guesses, kwargs...)
    MTK.add_equational_constraints!(prob.wrapped_model, sys, pmap, tspan)
    prob
end

MTK.generate_internal_model(m::Type{PyomoModel}) = pyomo.ConcreteModel()
function MTK.generate_time_variable!(m::ConcreteModel, tspan, tsteps)
    m.t = pyomo.ContinuousSet(initialize = collect(tsteps), bounds = tspan)
end

function MTK.generate_state_variable!(m::ConcreteModel, u0, ns, ts) 
    m.u_idxs = pyomo.RangeSet(1, ns)
    pyomo.Var(m.u_idxs, m.t)
end

function MTK.generate_input_variable!(m::ConcreteModel, u0, nc, ts)
    m.v_idxs = pyomo.RangeSet(1, nc)
    pyomo.Var(m.v_idxs, m.t)
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

function substitute_fixed_t_vars!(model::PyomoModel, sys, exprs)
    stidxmap = Dict([v => i for (i, v) in enumerate(unknowns(sys))])
    ctidxmap = Dict([v => i for (i, v) in enumerate(MTK.unbound_inputs(sys))])
    iv = MTK.get_iv(sys)

    for cons in jconstraints
        consvars = MTK.vars(cons)
        for st in consvars
            MTK.iscall(st) || continue
            x = MTK.operation(st)
            t = only(MTK.arguments(st))
            MTK.symbolic_type(t) === MTK.NotSymbolic() || continue
            if haskey(stidxmap, x(iv))
                idx = stidxmap[x(iv)]
                cv = :U
            else
                idx = ctidxmap[x(iv)]
                cv = :V
            end
            model.t.add(t)
            cons = Symbolics.substitute(cons, Dict(x(t) => model.cv[idx, t]))
        end
    end
end

function MTK.substitute_model_vars(pmodel::PyomoModel, sys, pmap, exprs, tspan)
    @unwrap model, model_sym, idx_sym, t_sym = pmodel
    x_ops = [MTK.operation(MTK.unwrap(st)) for st in unknowns(sys)]
    c_ops = [MTK.operation(MTK.unwrap(ct)) for ct in MTK.unbound_inputs(sys)]
    mU = Symbolics.symbolic_getproperty(model_sym, :U)
    mV = Symbolics.symbolic_getproperty(model_sym, :V)

    (ti, tf) = tspan
    if MTK.symbolic_type(tf) === MTK.ScalarSymbolic()
        _tf = model.tₛ + ti
        exprs = map(c -> Symbolics.fast_substitute(c, Dict(tf => _tf)), exprs)
        free_t_map = Dict([[x(tₛ) => mU[i, end] for (i, x) in enumerate(x_ops)];
                           [c(tₛ) => mV[i, end] for (i, c) in enumerate(c_ops)]])
        exprs = map(c -> Symbolics.fixpoint_sub(c, free_t_map), exprs)
    end

    whole_interval_map = Dict([[v => mU[i, t_sym] for (i, v) in enumerate(sts)];
                               [v => mV[i, t_sym] for (i, v) in enumerate(cts)]])
    exprs = map(c -> Symbolics.fixpoint_sub(c, whole_interval_map), exprs)
    exprs
end

function MTK.substitute_integral!(model::PyomoModel, exprs, tspan)
    intmap = Dict()
    for int in MTK.collect_applied_operators(exprs, Symbolics.Integral)
        op = MTK.operation(int)
        arg = only(arguments(MTK.value(int)))
        lo, hi = MTK.value.((op.domain.domain.left, op.domain.domain.right))
        if MTK.is_free_final(model) && isequal((lo, hi), tspan)
            (lo, hi) = (0, 1)
        elseif MTK.is_free_final(model)
            error("Free final time problems cannot handle partial timespans.")
        end
        intmap[int] = model.tₛ * InfiniteOpt.∫(arg, model.model[:t], lo, hi)
    end
    exprs = map(c -> Symbolics.substitute(c, intmap), exprs)
end

function MTK.substitute_differentials(model::PyomoModel, sys, eqs)
    pmodel = prob.model
    @unpack model, model_sym, t_sym, idx_sym = pmodel
    model.dU = pyomo.DerivativeVar(model.U, wrt = model.t)

    mdU = Symbolics.symbolic_getproperty(model_sym, :dU)
    mU = Symbolics.symbolic_getproperty(model_sym, :U)
    mtₛ = Symbolics.symbolic_getproperty(model_sym, :tₛ)
    diffsubmap = Dict([D(mU[i, t_sym]) => mdU[i, t_sym] for i in 1:length(unknowns(sys))])

    diff_eqs = substitute_model_vars(model, sys, pmap, diff_equations(sys))
    diff_eqs = map(e -> Symbolics.substitute(e, diffsubmap), diff_eqs)
    [mtₛ * eq.rhs - eq.lhs == 0 for eq in diff_eqs]
end

struct PyomoCollocation <: AbstractCollocation
    solver::Any
    derivative_method
end
MTK.PyomoCollocation(solver, derivative_method = 1) = PyomoCollocation(solver, derivative_method)

function MTK.prepare_and_optimize!()
end
function MTK.get_U_values()
end
function MTK.get_V_values()
end
function MTK.get_t_values()
end
function MTK.successful_solve()
end

end
