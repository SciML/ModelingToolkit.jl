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
    function change_independent_variable(sys::AbstractODESystem, iv, eqs = []; dummies = false, simplify = true, verbose = false, kwargs...)

Transform the independent variable (e.g. ``t``) of the ODE system `sys` to a dependent variable `iv` (e.g. ``f(t)``).
An equation in `sys` must define the rate of change of the new independent variable (e.g. ``df(t)/dt``).
Alternatively, `eq` can specify such an equation.

The transformation is well-defined when the mapping between the new and old independent variables are one-to-one.
This is satisfied if one is a strictly increasing function of the other (e.g. ``df(t)/dt > 0`` or ``df(t)/dt < 0``).

Keyword arguments
=================
If `dummies`, derivatives of the new independent variable are expressed through dummy equations; otherwise they are explicitly inserted into the equations.
If `simplify`, these dummy expressions are simplified and often give a tidier transformation.
If `verbose`, the function prints intermediate transformations of equations to aid debugging.
Any additional keyword arguments `kwargs...` are forwarded to the constructor that rebuilds the system.

Usage before structural simplification
======================================
The variable change must take place before structural simplification.
Subsequently, consider passing `allow_symbolic = true` to `structural_simplify(sys)` to reduce the number of unknowns, with the understanding that the transformation is well-defined.

Example
=======
Consider a free fall with constant horizontal velocity.
The laws of physics naturally describes position as a function of time.
By changing the independent variable, it can be reformulated for vertical position as a function of horizontal distance:
```julia
julia> @variables x(t) y(t);

julia> @named M = ODESystem([D(D(y)) ~ -9.81, D(x) ~ 10.0], t);

julia> M′ = change_independent_variable(complete(M), x);

julia> unknowns(M′)
1-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 y(x)
```
"""
function change_independent_variable(sys::AbstractODESystem, iv, eqs = []; dummies = false, simplify = true, verbose = false, kwargs...)
    if !iscomplete(sys)
        error("Cannot change independent variable of incomplete system $(nameof(sys))")
    elseif isscheduled(sys)
        error("Cannot change independent variable of structurally simplified system $(nameof(sys))")
    elseif !isempty(get_systems(sys))
        error("Cannot change independent variable of hierarchical system $(nameof(sys)). Flatten it first.") # TODO: implement
    end

    iv = unwrap(iv)
    iv1 = get_iv(sys) # e.g. t
    if !iscall(iv) || !isequal(only(arguments(iv)), iv1)
        error("New independent variable $iv is not a function of the independent variable $iv1 of the system $(nameof(sys))")
    elseif !isautonomous(sys) && isempty(findall(eq -> isequal(eq.lhs, iv1), eqs))
        error("System $(nameof(sys)) is autonomous in $iv1. An equation of the form $iv1 ~ F($iv) must be provided.")
    end

    iv2func = iv # e.g. a(t)
    iv2name = nameof(operation(iv))
    iv2, = @independent_variables $iv2name # e.g. a
    D1 = Differential(iv1)

    iv1name = nameof(iv1) # e.g. t
    iv1func, = @variables $iv1name(iv2) # e.g. t(a)

    eqs = [get_eqs(sys); eqs] # copies system equations to avoid modifying original system

    # 1) Find and compute all necessary expressions for e.g. df/dt, d²f/dt², ...
    # 1.1) Find the 1st order derivative of the new independent variable (e.g. da(t)/dt = ...), ...
    div2_div1_idxs = findall(eq -> isequal(eq.lhs, D1(iv2func)), eqs) # index of e.g. da/dt = ...
    if length(div2_div1_idxs) != 1
        error("Exactly one equation for $D1($iv2func) was not specified.")
    end
    div2_div1_eq = popat!(eqs, only(div2_div1_idxs)) # get and remove e.g. df/dt = ... (may be added back later)
    div2_div1 = div2_div1_eq.rhs
    if isequal(div2_div1, 0)
        error("Cannot change independent variable from $iv1 to $iv2 with singular transformation $div2_div1_eq.")
    end
    # 1.2) ... then compute the 2nd order derivative of the new independent variable
    div1_div2 = 1 / div2_div1 # TODO: URL reference for clarity
    ddiv2_ddiv1 = expand_derivatives(-Differential(iv2func)(div1_div2) / div1_div2^3, simplify) # e.g. https://math.stackexchange.com/questions/249253/second-derivative-of-the-inverse-function # TODO: higher orders # TODO: pass simplify here
    # 1.3) # TODO: handle higher orders (3+) derivatives ...

    # 2) If requested, insert extra dummy equations for e.g. df/dt, d²f/dt², ...
    #    Otherwise, replace all these derivatives by their explicit expressions
    if dummies
        div2name = Symbol(iv2name, :_t) # TODO: not always t
        div2, = @variables $div2name(iv2) # e.g. a_t(a)
        ddiv2name = Symbol(iv2name, :_tt) # TODO: not always t
        ddiv2, = @variables $ddiv2name(iv2) # e.g. a_tt(a)
        eqs = [eqs; [div2 ~ div2_div1, ddiv2 ~ ddiv2_ddiv1]] # add dummy equations
        derivsubs = [D1(D1(iv2func)) => ddiv2, D1(iv2func) => div2] # order is crucial!
    else
        derivsubs = [D1(D1(iv2func)) => ddiv2_ddiv1, D1(iv2func) => div2_div1] # order is crucial!
    end
    derivsubs = [derivsubs; [iv2func => iv2, iv1 => iv1func]]

    if verbose
        # Explain what we just did
        println("Order 1 (found):    $div2_div1_eq")
        println("Order 2 (computed): $(D1(div2_div1_eq.lhs) ~ ddiv2_ddiv1)")
        println()
        println("Substitutions will be made in this order:")
        for (n, sub) in enumerate(derivsubs)
            println("$n: $(sub[1]) => $(sub[2])")
        end
        println()
    end

    # 3) Define a transformation function that performs the change of variable on any expression/equation
    function transform(ex)
        verbose && println("Step 0: ", ex)

        # Step 1: substitute f(t₁) => f(t₂(t₁)) in all variables in the expression
        vars = Symbolics.get_variables(ex)
        for var1 in vars
            if Symbolics.iscall(var1) && !isequal(var1, iv2func) # && isequal(only(arguments(var1)), iv1) # skip e.g. constants
                name = nameof(operation(var1))
                var2, = @variables $name(iv2func)
                ex = substitute(ex, var1 => var2; fold = false)
            end
        end
        verbose && println("Step 1: ", ex)

        # Step 2: expand out all chain rule derivatives
        ex = expand_derivatives(ex) # expand out with chain rule to get d(iv2)/d(iv1)
        verbose && println("Step 2: ", ex)

        # Step 3: substitute d²f/dt², df/dt, ... (to dummy variables or explicit expressions, depending on dummies)
        for sub in derivsubs
            ex = substitute(ex, sub)
        end
        verbose && println("Step 3: ", ex)
        verbose && println()

        return ex
    end

    # 4) Transform all fields
    eqs = map(transform, eqs)
    observed = map(transform, get_observed(sys))
    initialization_eqs = map(transform, get_initialization_eqs(sys))
    parameter_dependencies = map(transform, get_parameter_dependencies(sys))
    defaults = Dict(transform(var) => transform(val) for (var, val) in get_defaults(sys))
    guesses = Dict(transform(var) => transform(val) for (var, val) in get_guesses(sys))
    assertions = Dict(transform(condition) => msg for (condition, msg) in get_assertions(sys))
    # TODO: handle subsystems

    # 5) Recreate system with transformed fields
    return typeof(sys)(
        eqs, iv2;
        observed, initialization_eqs, parameter_dependencies, defaults, guesses, assertions,
        name = nameof(sys), description = description(sys), kwargs...
    ) |> complete # original system had to be complete
end
