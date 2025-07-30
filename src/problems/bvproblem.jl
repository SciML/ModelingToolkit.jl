"""
    BVProblem(sys::System, op, tspan; kwargs...)

Construct a `BVProblem` from a ModelingToolkit `System` for boundary value problems.

## Additional Keyword Arguments

Beyond the arguments listed below, this constructor accepts all keyword arguments
supported by the DifferentialEquations.jl `solve` function. For a complete list
and detailed descriptions, see the [ModelingToolkit.jl problem building documentation](https://docs.sciml.ai/ModelingToolkit/stable/API/problems/) 
and the [DifferentialEquations.jl solve documentation](https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/).

## Arguments
- `sys::System`: The ModelingToolkit system to convert (with boundary conditions)
- `op`: Operating point/initial conditions
- `tspan`: Time span for the problem

## Keywords
- `check_compatibility=true`: Whether to check system compatibility
- `cse=true`: Whether to perform common subexpression elimination
- `checkbounds=false`: Whether to enable bounds checking
- `eval_expression=false`: Whether to evaluate expressions
- `eval_module=@__MODULE__`: Module for expression evaluation
- `expression=Val{false}`: Expression evaluation mode
- `guesses=Dict()`: Initial guesses for boundary value problem
- `callback=nothing`: Callback functions (note: BVP solvers do not support callbacks)
- `kwargs...`: Additional keyword arguments passed to the solver
"""
@fallback_iip_specialize function SciMLBase.BVProblem{iip, spec}(
        sys::System, op, tspan;
        check_compatibility = true, cse = true,
        checkbounds = false, eval_expression = false, eval_module = @__MODULE__,
        expression = Val{false}, guesses = Dict(), callback = nothing,
        kwargs...) where {iip, spec}
    check_complete(sys, BVProblem)
    check_compatibility && check_compatible_system(BVProblem, sys)
    isnothing(callback) || error("BVP solvers do not support callbacks.")

    dvs = unknowns(sys)
    op = to_varmap(op, dvs)
    # Systems without algebraic equations should use both fixed values + guesses
    # for initialization.
    _op = has_alg_eqs(sys) ? op : merge(Dict(op), Dict(guesses))

    fode, u0,
    p = process_SciMLProblem(
        ODEFunction{iip, spec}, sys, _op; guesses,
        t = tspan !== nothing ? tspan[1] : tspan, check_compatibility = false, cse,
        checkbounds, time_dependent_init = false, expression, kwargs...)

    stidxmap = Dict([v => i for (i, v) in enumerate(dvs)])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(dvs)) :
              [stidxmap[k] for (k, v) in op if haskey(stidxmap, k)]
    fbc = generate_boundary_conditions(
        sys, u0, u0_idxs, tspan[1]; expression = Val{false},
        wrap_gfw = Val{true}, cse, checkbounds)

    if (length(constraints(sys)) + length(op) > length(dvs))
        @warn "The BVProblem is overdetermined. The total number of conditions (# constraints + # fixed initial values given by op) exceeds the total number of states. The BVP solvers will default to doing a nonlinear least-squares optimization."
    end

    kwargs = process_kwargs(sys; expression, kwargs...)
    args = (; fode, fbc, u0, tspan, p)

    return maybe_codegen_scimlproblem(expression, BVProblem{iip}, args; kwargs...)
end

function check_compatible_system(T::Type{BVProblem}, sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    check_is_continuous(sys, T)

    if !isempty(discrete_events(sys)) || !isempty(continuous_events(sys))
        throw(SystemCompatibilityError("BVP solvers do not support events."))
    end
end
