module MTKCasADiDynamicOptExt
using ModelingToolkit
using CasADi
using DiffEqDevTools, DiffEqBase
using DataInterpolations
const MTK = MOdelingToolkit

struct CasADiDynamicOptProblem{uType, tType, isinplace, P, F, K} <:
       AbstractDynamicOptProblem{uType, tType, isinplace}
    f::F
    u0::uType
    tspan::tType
    p::P
    model::Opti
    kwargs::K

    function CasADiDynamicOptProblem(f, u0, tspan, p, model, kwargs...)
        new{typeof(u0), typeof(tspan), SciMLBase.isinplace(f, 5),
            typeof(p), typeof(f), typeof(kwargs)}(f, u0, tspan, p, model, kwargs)
    end
end

struct CasADiModel
    opti::Opti
    U::MX
    V::MX
end

struct TimedMX
end

function MTK.CasADiDynamicOptProblem(sys::ODESystem, u0map, tspan, pmap;
        dt = nothing,
        steps = nothing, 
        guesses = Dict(), kwargs...) 
    MTK.warn_overdetermined(sys, u0map)
    _u0map = has_alg_eqs(sys) ? u0map : merge(Dict(u0map), Dict(guesses))
    f, u0, p = MTK.process_SciMLProblem(ODEInputFunction, sys, _u0map, pmap;
        t = tspan !== nothing ? tspan[1] : tspan, kwargs...)

    pmap = Dict{Any, Any}(pmap)
    steps, is_free_t = MTK.process_tspan(tspan, dt, steps)
    model = init_model()
end

function init_model(sys, tspan, steps, u0map, pmap, u0; is_free_t)
    ctrls = MTK.unbound_inputs(sys)
    states = unknowns(sys)
    model = CasADi.Opti()
    
    U = CasADi.variable!(model, length(states), steps)
    V = CasADi.variable!(model, length(ctrls), steps)
end

function add_initial_constraints!()
    
end

function add_user_constraints!(model::CasADiModel, sys, pmap; is_free_t = false)
    
end

function add_cost_function!(model)

end

function add_solve_constraints!(prob, tableau; is_free_t)
    A = tableau.A
    α = tableau.α
    c = tableau.c
    model = prob.model
    f = prob.f
    p = prob.p

    opti = model.opti
    t = model[:t]
    tsteps = supports(t)
    tmax = tsteps[end]
    pop!(tsteps)
    tₛ = is_free_t ? model[:tf] : 1
    dt = tsteps[2] - tsteps[1]

    U = model.U
    V = model.V
    nᵤ = length(U)
    nᵥ = length(V)

    if is_explicit(tableau)
        K = Any[]
        for τ in tsteps
            for (i, h) in enumerate(c)
                ΔU = sum([A[i, j] * K[j] for j in 1:(i - 1)], init = zeros(nᵤ))
                Uₙ = [U[i](τ) + ΔU[i] * dt for i in 1:nᵤ]
                Vₙ = [V[i](τ) for i in 1:nᵥ]
                Kₙ = tₛ * f(Uₙ, Vₙ, p, τ + h * dt) # scale the time
                push!(K, Kₙ)
            end
            ΔU = dt * sum([α[i] * K[i] for i in 1:length(α)])
            subject_to!(model, U[n](τ) + ΔU[n]==U[n](τ + dt))
            empty!(K)
        end
    else
    end
end

function DiffEqBase.solve(prob::CasADiDynamicOptProblem, solver::Union{String, Symbol}, ode_solver::Union{String, Symbol}; silent = false)
    model = prob.model
    opti = model.opti

    solver!(opti, solver)
    sol = solve(opti)
    DynamicOptSolution(model, sol, input_sol)
end

end
