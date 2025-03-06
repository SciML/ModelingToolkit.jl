"""
$(TYPEDSIGNATURES)

Generates the Liouville transformed set of ODEs, which is the original
ODE system with a new variable `trJ` appended, corresponding to the
-tr(Jacobian). This variable is used for properties like uncertainty
propagation from a given initial distribution density.

For example, if ``u'=p*u`` and `p` follows a probability distribution
``f(p)``, then the probability density of a future value with a given
choice of ``p`` is computed by setting the initial `trJ = f(p)`, and
the final value of `trJ` is the probability of ``u(t)``.

Example:

```julia
using ModelingToolkit, OrdinaryDiffEq

@independent_variables t
@parameters α β γ δ
@variables x(t) y(t)
D = Differential(t)
eqs = [D(x) ~ α*x - β*x*y, D(y) ~ -δ*y + γ*x*y]
@named sys = ODESystem(eqs, t)

sys2 = liouville_transform(sys)
sys2 = complete(sys2)
u0 = [x => 1.0, y => 1.0, sys2.trJ => 1.0]
prob = ODEProblem(sys2, u0, tspan, p)
sol = solve(prob, Tsit5())
```

Where `sol[3,:]` is the evolution of `trJ` over time.

Sources:

Probabilistic Robustness Analysis of F-16 Controller Performance: An
Optimal Transport Approach

Abhishek Halder, Kooktae Lee, and Raktim Bhattacharya
https://abhishekhalder.bitbucket.io/F16ACC2013Final.pdf
"""
function liouville_transform(sys::AbstractODESystem; kwargs...)
    t = get_iv(sys)
    @variables trJ
    D = ModelingToolkit.Differential(t)
    neweq = D(trJ) ~ trJ * -tr(calculate_jacobian(sys))
    neweqs = [equations(sys); neweq]
    vars = [unknowns(sys); trJ]
    ODESystem(neweqs, t, vars, parameters(sys); checks = false, name = nameof(sys), kwargs...)
end

"""
    change_independent_variable(sys::AbstractODESystem, iv, eqs = []; dummies = false, simplify = true, fold = false, kwargs...)

Transform the independent variable (e.g. ``t``) of the ODE system `sys` to a dependent variable `iv` (e.g. ``u(t)``).
An equation in `sys` must define the rate of change of the new independent variable (e.g. ``du(t)/dt``).
This or other additional equations can also be specified through `eqs`.

The transformation is well-defined when the mapping between the new and old independent variables are one-to-one.
This is satisfied if one is a strictly increasing function of the other (e.g. ``du(t)/dt > 0`` or ``du(t)/dt < 0``).

Keyword arguments
=================
If `dummies`, derivatives of the new independent variable with respect to the old one are expressed through dummy equations; otherwise they are explicitly inserted into the equations.
If `simplify`, these dummy expressions are simplified and often give a tidier transformation.
If `fold`, internal substitutions will evaluate numerical expressions.
Additional keyword arguments `kwargs...` are forwarded to the constructor that rebuilds `sys`.

Usage before structural simplification
======================================
The variable change must take place before structural simplification.
Subsequently, consider passing `allow_symbolic = true` to `structural_simplify(sys)` to reduce the number of unknowns, with the understanding that the transformation is well-defined.

Usage with non-autonomous systems
=================================
If `sys` is autonomous (i.e. ``t`` appears explicitly in its equations), it is often desirable to also pass an algebraic equation relating the new and old independent variables (e.g. ``t = f(u(t))``).
Otherwise the transformed system will be underdetermined and cannot be structurally simplified without additional changes.

Usage with hierarchical systems
===============================
It is recommended that `iv` is a non-namespaced variable in `sys`.
This means it can belong to the top-level system or be a variable in a subsystem declared with `GlobalScope`.

Example
=======
Consider a free fall with constant horizontal velocity.
Physics naturally describes position as a function of time.
By changing the independent variable, it can be reformulated for vertical position as a function of horizontal position:
```julia
julia> @variables x(t) y(t);

julia> @named M = ODESystem([D(D(y)) ~ -9.81, D(x) ~ 10.0], t);

julia> M′ = change_independent_variable(complete(M), x);

julia> unknowns(M′)
1-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 y(x)
```
"""
function change_independent_variable(sys::AbstractODESystem, iv, eqs = []; dummies = false, simplify = true, fold = false, kwargs...)
    iv2_of_iv1 = unwrap(iv) # e.g. u(t)
    iv1 = get_iv(sys) # e.g. t

    if !iscomplete(sys)
        error("System $(nameof(sys)) is incomplete. Complete it first!")
    elseif is_dde(sys)
        error("System $(nameof(sys)) contains delay differential equations (DDEs). This is currently not supported!")
    elseif isscheduled(sys)
        error("System $(nameof(sys)) is structurally simplified. Change independent variable before structural simplification!")
    elseif !iscall(iv2_of_iv1) || !isequal(only(arguments(iv2_of_iv1)), iv1)
        error("Variable $iv is not a function of the independent variable $iv1 of the system $(nameof(sys))!")
    end

    iv1name = nameof(iv1) # e.g. :t
    iv2name = nameof(operation(iv2_of_iv1)) # e.g. :u
    iv2, = @independent_variables $iv2name # e.g. u
    iv1_of_iv2, = @variables $iv1name(iv2) # inverse in case sys is autonomous; e.g. t(u)
    iv1_of_iv2 = GlobalScope(iv1_of_iv2) # do not namespace old independent variable as new dependent variable
    D1 = Differential(iv1) # e.g. d/d(t)
    D2 = Differential(iv2_of_iv1) # e.g. d/d(u(t))

    # 1) Utility that performs the chain rule on an expression, e.g. (d/dt)(f(t)) -> (d/dt)(f(u(t))) -> df(u(t))/du(t) * du(t)/dt
    function chain_rule(ex)
        vars = get_variables(ex)
        for var_of_iv1 in vars # loop through e.g. f(t)
            if iscall(var_of_iv1) && isequal(only(arguments(var_of_iv1)), iv1) && !isequal(var_of_iv1, iv2_of_iv1) # handle e.g. f(t) -> f(u(t)), but not u(t) -> u(u(t))
                var_of_iv2 = substitute(var_of_iv1, iv1 => iv2_of_iv1) # e.g. f(t) -> f(u(t))
                ex = substitute(ex, var_of_iv1 => var_of_iv2)
            end
        end
        ex = expand_derivatives(ex, simplify) # expand chain rule, e.g. (d/dt)(f(u(t)))) -> df(u(t))/du(t) * du(t)/dt
        return ex
    end

    # 2) Find e.g. du/dt in equations, then calculate e.g. d²u/dt², ...
    eqs = [eqs; get_eqs(sys)] # all equations (system-defined + user-provided) we may use
    idxs = findall(eq -> isequal(eq.lhs, D1(iv2_of_iv1)), eqs)
    if length(idxs) != 1
        error("Exactly one equation for $D1($iv2_of_iv1) was not specified!")
    end
    div2_of_iv1_eq = popat!(eqs, only(idxs)) # get and remove e.g. du/dt = ... (may be added back later as a dummy)
    div2_of_iv1 = chain_rule(div2_of_iv1_eq.rhs)
    if isequal(div2_of_iv1, 0) # e.g. du/dt ~ 0
        error("Independent variable transformation $(div2_of_iv1_eq) is singular!")
    end
    ddiv2_of_iv1 = chain_rule(D1(div2_of_iv1)) # TODO: implement higher orders (order >= 3) derivatives with a loop

    # 3) If requested, insert extra dummy equations for e.g. du/dt, d²u/dt², ...
    #    Otherwise, replace all these derivatives by their explicit expressions
    if dummies
        div2name = Symbol(iv2name, :_, iv1name) # e.g. :u_t # TODO: customize
        ddiv2name = Symbol(div2name, iv1name) # e.g. :u_tt # TODO: customize
        div2, ddiv2 = @variables $div2name(iv2) $ddiv2name(iv2) # e.g. u_t(u), u_tt(u)
        div2, ddiv2 = GlobalScope.([div2, ddiv2]) # do not namespace dummies in new system
        eqs = [[div2 ~ div2_of_iv1, ddiv2 ~ ddiv2_of_iv1]; eqs] # add dummy equations
        @set! sys.unknowns = [get_unknowns(sys); [div2, ddiv2]] # add dummy variables
    else
        div2 = div2_of_iv1
        ddiv2 = ddiv2_of_iv1
    end
    @set! sys.eqs = eqs # add extra equations we derived before starting transformation process

    # 4) Transform everything from old to new independent variable, e.g. t -> u.
    #    Substitution order matters! Must begin with highest order to get D(D(u(t))) -> u_tt(u).
    #    If we had started with the lowest order, we would get D(D(u(t))) -> D(u_t(u)) -> 0!
    iv1_to_iv2_subs = [ # a vector ensures substitution order
        D1(D1(iv2_of_iv1)) => ddiv2 # order 2, e.g. D(D(u(t))) -> u_tt(u)
        D1(iv2_of_iv1) => div2 # order 1, e.g. D(u(t)) -> u_t(u)
        iv2_of_iv1 => iv2 # order 0, e.g. u(t) -> u
        iv1 => iv1_of_iv2 # in case sys was autonomous, e.g. t -> t(u)
    ]
    function transform(ex)
        ex = chain_rule(ex)
        for sub in iv1_to_iv2_subs
            ex = substitute(ex, sub; fold)
        end
        return ex
    end
    function transform(sys::AbstractODESystem)
        eqs = map(transform, get_eqs(sys))
        unknowns = map(transform, get_unknowns(sys))
        unknowns = filter(var -> !isequal(var, iv2), unknowns) # remove e.g. u
        ps = map(transform, get_ps(sys))
        ps = filter(!isinitial, ps) # remove Initial(...) # # TODO: shouldn't have to touch this
        observed = map(transform, get_observed(sys))
        initialization_eqs = map(transform, get_initialization_eqs(sys))
        parameter_dependencies = map(transform, get_parameter_dependencies(sys))
        defaults = Dict(transform(var) => transform(val) for (var, val) in get_defaults(sys))
        guesses = Dict(transform(var) => transform(val) for (var, val) in get_guesses(sys))
        assertions = Dict(transform(condition) => msg for (condition, msg) in get_assertions(sys))
        systems = get_systems(sys) # save before reconstructing system
        wascomplete = iscomplete(sys) # save before reconstructing system
        sys = typeof(sys)( # recreate system with transformed fields
            eqs, iv2, unknowns, ps; observed, initialization_eqs, parameter_dependencies, defaults, guesses,
            assertions, name = nameof(sys), description = description(sys), kwargs...
        )
        systems = map(transform, systems) # recurse through subsystems
        sys = compose(sys, systems) # rebuild hierarchical system
        if wascomplete
            wasflat = isempty(systems)
            sys = complete(sys; flatten = wasflat) # complete output if input was complete
        end
        return sys
    end
    return transform(sys)
end
