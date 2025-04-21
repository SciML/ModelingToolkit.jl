"""
```julia
SciMLBase.BVProblem{iip}(sys::AbstractSystem, u0map, tspan,
                         parammap = DiffEqBase.NullParameters();
                         constraints = nothing, guesses = nothing,
                         version = nothing, tgrad = false,
                         jac = true, sparse = true,
                         simplify = false,
                         kwargs...) where {iip}
```

Create a boundary value problem from the [`System`](@ref). 

`u0map` is used to specify fixed initial values for the states. Every variable 
must have either an initial guess supplied using `guesses` or a fixed initial 
value specified using `u0map`.

Boundary value conditions are supplied to Systems in the form of a list of constraints.
These equations  should specify values that state variables should take at specific points,
as in `x(0.5) ~ 1`). More general constraints that  should hold over the entire solution,
such as `x(t)^2 + y(t)^2`, should be  specified as one of the equations used to build the
`System`.

If a `System` without `constraints` is specified, it will be treated as an initial value problem. 

```julia
    @parameters g t_c = 0.5
    @variables x(..) y(t) λ(t)
    eqs = [D(D(x(t))) ~ λ * x(t)
           D(D(y)) ~ λ * y - g
           x(t)^2 + y^2 ~ 1]
    cstr = [x(0.5) ~ 1]
    @mtkbuild pend = System(eqs, t; constraints = cstrs)

    tspan = (0.0, 1.5)
    u0map = [x(t) => 0.6, y => 0.8]
    parammap = [g => 1]
    guesses = [λ => 1]

    bvp = SciMLBase.BVProblem{true, SciMLBase.AutoSpecialize}(pend, u0map, tspan, parammap; guesses, check_length = false)
```

If the `System` has algebraic equations, like `x(t)^2 + y(t)^2`, the resulting 
`BVProblem` must be solved using BVDAE solvers, such as Ascher.
"""
@fallback_iip_specialize function SciMLBase.BVProblem{iip, spec}(
        sys::System, u0map, tspan, parammap = SciMLBase.NullParameters();
        check_compatibility = true, cse = true, checkbounds = false, eval_expression = false,
        eval_module = @__MODULE__, guesses = Dict(), kwargs...) where {iip, spec}
    check_complete(sys, BVProblem)
    check_compatibility && check_compatible_system(BVProblem, sys)

    # ODESystems without algebraic equations should use both fixed values + guesses
    # for initialization.
    _u0map = has_alg_eqs(sys) ? u0map : merge(Dict(u0map), Dict(guesses))
    fode, u0, p = process_SciMLProblem(
        ODEFunction{iip, spec}, sys, _u0map, parammap; guesses,
        t = tspan !== nothing ? tspan[1] : tspan, check_compatibility = false, cse, checkbounds,
        time_dependent_init = false, kwargs...)

    dvs = unknowns(sys)
    stidxmap = Dict([v => i for (i, v) in enumerate(dvs)])
    u0_idxs = has_alg_eqs(sys) ? collect(1:length(dvs)) : [stidxmap[k] for (k, v) in u0map]
    fbc = generate_boundary_conditions(
        sys, u0, u0_idxs, tspan; expression = Val{false}, cse, checkbounds)
    kwargs = process_kwargs(sys; kwargs...)
    # Call `remake` so it runs initialization if it is trivial
    return remake(BVProblem{iip}(fode, fbc, u0, tspan[1], p; kwargs...))
end

function check_compatible_system(T::Union{Type{BVPFunction}, Type{BVProblem}}, sys::System)
    check_time_dependent(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
    check_is_continuous(sys, T)
end
