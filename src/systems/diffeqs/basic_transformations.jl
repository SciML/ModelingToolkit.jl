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
    ODESystem(
        neweqs, t, vars, parameters(sys);
        checks = false, name = nameof(sys), kwargs...
    )
end

"""
    change_independent_variable(
        sys::AbstractODESystem, iv, eqs = [];
        add_old_diff = false, simplify = true, fold = false
    )

Transform the independent variable (e.g. ``t``) of the ODE system `sys` to a dependent variable `iv` (e.g. ``u(t)``).
The transformation is well-defined when the mapping between the new and old independent variables are one-to-one.
This is satisfied if one is a strictly increasing function of the other (e.g. ``du(t)/dt > 0`` or ``du(t)/dt < 0``).

Any extra equations `eqs` involving the new and old independent variables will be taken into account in the transformation.

# Keyword arguments

- `add_old_diff`: Whether to add a differential equation for the old independent variable in terms of the new one using the inverse function rule ``dt/du = 1/(du/dt)``.
- `simplify`: Whether expanded derivative expressions are simplified. This can give a tidier transformation.
- `fold`: Whether internal substitutions will evaluate numerical expressions.

# Usage before structural simplification

The variable change must take place before structural simplification.
In following calls to `structural_simplify`, consider passing `allow_symbolic = true` to avoid undesired constraint equations between between dummy variables.

# Usage with non-autonomous systems

If `sys` is non-autonomous (i.e. ``t`` appears explicitly in its equations), consider passing an algebraic equation relating the new and old independent variables (e.g. ``t = f(u(t))``).
Otherwise the transformed system can be underdetermined.
If an algebraic relation is not known, consider using `add_old_diff` instead.

# Usage with hierarchical systems

It is recommended that `iv` is a non-namespaced variable in `sys`.
This means it can belong to the top-level system or be a variable in a subsystem declared with `GlobalScope`.

# Example

Consider a free fall with constant horizontal velocity.
Physics naturally describes position as a function of time.
By changing the independent variable, it can be reformulated for vertical position as a function of horizontal position:
```julia
julia> @variables x(t) y(t);

julia> @named M = ODESystem([D(D(y)) ~ -9.81, D(D(x)) ~ 0.0], t);

julia> M = change_independent_variable(M, x);

julia> M = structural_simplify(M; allow_symbolic = true);

julia> unknowns(M)
3-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 xˍt(x)
 y(x)
 yˍx(x)
```
"""
function change_independent_variable(
        sys::AbstractODESystem, iv, eqs = [];
        add_old_diff = false, simplify = true, fold = false
)
    iv2_of_iv1 = unwrap(iv) # e.g. u(t)
    iv1 = get_iv(sys) # e.g. t

    if is_dde(sys)
        error("System $(nameof(sys)) contains delay differential equations (DDEs). This is currently not supported!")
    elseif isscheduled(sys)
        error("System $(nameof(sys)) is structurally simplified. Change independent variable before structural simplification!")
    elseif !iscall(iv2_of_iv1) || !isequal(only(arguments(iv2_of_iv1)), iv1)
        error("Variable $iv is not a function of the independent variable $iv1 of the system $(nameof(sys))!")
    end

    # Set up intermediate and final variables for the transformation
    iv1name = nameof(iv1) # e.g. :t
    iv2name = nameof(operation(iv2_of_iv1)) # e.g. :u
    D1 = Differential(iv1) # e.g. d/d(t)

    # construct new terms, e.g:
    #   iv2 -> u
    #   iv1_of_iv2 -> t(u), (inverse, global because iv1 has no namespacing in sys)
    #   div2_of_iv1 -> uˍt(t)
    iv2_unit = getmetadata(iv2_of_iv1, VariableUnit, nothing)
    if isnothing(iv2_unit)
        iv2, = @independent_variables $iv2name
        iv1_of_iv2, = GlobalScope.(@variables $iv1name(iv2))
        div2_of_iv1 = GlobalScope(default_toterm(D1(iv2_of_iv1)))
    else
        iv2, = @independent_variables $iv2name [unit = iv2_unit]
        iv1_of_iv2, = GlobalScope.(@variables $iv1name(iv2) [unit = get_unit(iv1)])
        div2_of_iv1 = GlobalScope(diff2term_with_unit(D1(iv2_of_iv1), iv1))
    end

    div2_of_iv2 = substitute(div2_of_iv1, iv1 => iv2) # e.g. uˍt(u)
    div2_of_iv2_of_iv1 = substitute(div2_of_iv2, iv2 => iv2_of_iv1) # e.g. uˍt(u(t))

    # If requested, add a differential equation for the old independent variable as a function of the old one
    if add_old_diff
        eqs = [eqs; Differential(iv2)(iv1_of_iv2) ~ 1 / div2_of_iv2] # e.g. dt(u)/du ~ 1 / uˍt(u) (https://en.wikipedia.org/wiki/Inverse_function_rule)
    end
    @set! sys.eqs = [get_eqs(sys); eqs] # add extra equations we derived
    @set! sys.unknowns = [get_unknowns(sys); [iv1, div2_of_iv1]] # add new variables, will be transformed to e.g. t(u) and uˍt(u)

    # Create a utility that performs the chain rule on an expression, followed by insertion of the new independent variable:
    # e.g. (d/dt)(f(t)) -> (d/dt)(f(u(t))) -> df(u(t))/du(t) * du(t)/dt -> df(u)/du * uˍt(u)
    function transform(ex)
        # 1) Replace the argument of every function; e.g. f(t) -> f(u(t))
        for var in vars(ex; op = Nothing) # loop over all variables in expression (op = Nothing prevents interpreting "D(f(t))" as one big variable)
            is_function_of_iv1 = iscall(var) && isequal(only(arguments(var)), iv1) # of the form f(t)?
            if is_function_of_iv1 && !isequal(var, iv2_of_iv1) # prevent e.g. u(t) -> u(u(t))
                var_of_iv1 = var # e.g. f(t)
                var_of_iv2_of_iv1 = substitute(var_of_iv1, iv1 => iv2_of_iv1) # e.g. f(u(t))
                ex = substitute(ex, var_of_iv1 => var_of_iv2_of_iv1; fold)
            end
        end
        # 2) Repeatedly expand chain rule until nothing changes anymore
        orgex = nothing
        while !isequal(ex, orgex)
            orgex = ex # save original
            ex = expand_derivatives(ex, simplify) # expand chain rule, e.g. (d/dt)(f(u(t)))) -> df(u(t))/du(t) * du(t)/dt
            ex = substitute(ex, D1(iv2_of_iv1) => div2_of_iv2_of_iv1; fold) # e.g. du(t)/dt -> uˍt(u(t))
        end
        # 3) Set new independent variable
        ex = substitute(ex, iv2_of_iv1 => iv2; fold) # set e.g. u(t) -> u everywhere
        ex = substitute(ex, iv1 => iv1_of_iv2; fold) # set e.g. t -> t(u) everywhere
        return ex
    end

    # Use the utility function to transform everything in the system!
    function transform(sys::AbstractODESystem)
        eqs = map(transform, get_eqs(sys))
        unknowns = map(transform, get_unknowns(sys))
        unknowns = filter(var -> !isequal(var, iv2), unknowns) # remove e.g. u
        ps = map(transform, get_ps(sys))
        ps = filter(!isinitial, ps) # remove Initial(...) # TODO: shouldn't have to touch this
        observed = map(transform, get_observed(sys))
        initialization_eqs = map(transform, get_initialization_eqs(sys))
        parameter_dependencies = map(transform, get_parameter_dependencies(sys))
        defaults = Dict(transform(var) => transform(val)
        for (var, val) in get_defaults(sys))
        guesses = Dict(transform(var) => transform(val) for (var, val) in get_guesses(sys))
        assertions = Dict(transform(ass) => msg for (ass, msg) in get_assertions(sys))
        systems = get_systems(sys) # save before reconstructing system
        wascomplete = iscomplete(sys) # save before reconstructing system
        sys = typeof(sys)( # recreate system with transformed fields
            eqs, iv2, unknowns, ps; observed, initialization_eqs,
            parameter_dependencies, defaults, guesses,
            assertions, name = nameof(sys), description = description(sys)
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
