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
    function change_independent_variable(sys::AbstractODESystem, iv, eq = nothing; dummies = false, simplify = true, verbose = false, kwargs...)

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
function change_independent_variable(sys::AbstractODESystem, iv, eq = nothing; dummies = false, simplify = true, verbose = false, kwargs...)
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
    end

    iv2func = iv # e.g. a(t)
    iv2name = nameof(operation(iv))
    iv2, = @independent_variables $iv2name # e.g. a
    D1 = Differential(iv1)
    D2 = Differential(iv2)

    div2name = Symbol(iv2name, :_t)
    div2, = @variables $div2name(iv2) # e.g. a_t(a)

    ddiv2name = Symbol(iv2name, :_tt)
    ddiv2, = @variables $ddiv2name(iv2) # e.g. a_tt(a)

    eqs = ModelingToolkit.get_eqs(sys) |> copy # don't modify original system
    !isnothing(eq) && push!(eqs, eq)
    vars = []
    div2_div1 = nothing
    for (i, eq) in enumerate(eqs)
        verbose && println("1. ", eq)

        # 1) Substitute f(t₁) => f(t₂(t₁)) in all variables
        vars = Symbolics.get_variables(eq)
        for var1 in vars
            if Symbolics.iscall(var1) && !isequal(var1, iv2func) # && isequal(only(arguments(var1)), iv1) # skip e.g. constants
                name = nameof(operation(var1))
                var2, = @variables $name(iv2func)
                eq = substitute(eq, var1 => var2; fold = false)
            end
        end
        verbose && println("2. ", eq)

        # 2) Substitute in dummy variables for dⁿt₂/dt₁ⁿ
        eq = expand_derivatives(eq) # expand out with chain rule to get d(iv2)/d(iv1)
        verbose && println("3. ", eq)
        eq = substitute(eq, D1(D1(iv2func)) => ddiv2) # order 2 # TODO: more orders
        eq = substitute(eq, D1(iv2func) => div2) # order 1; e.g. D(a(t)) => a_t(t)
        eq = substitute(eq, iv2func => iv2) # order 0; make iv2 independent
        verbose && println("4. ", eq)
        verbose && println()

        eqs[i] = eq

        if isequal(eq.lhs, div2)
            div2_div1 = eq.rhs
        end
    end

    verbose && println("Found $div2 = $div2_div1")
    if isnothing(div2_div1)
        error("No equation for $D1($iv2func) was specified.")
    elseif isequal(div2_div1, 0)
        error("Cannot change independent variable from $iv1 to $iv2 with singular transformation $(div2 ~ div2_div1).")
    end

    # 3) Add equations for dummy variables
    div1_div2 = 1 / div2_div1
    ddiv1_ddiv2 = expand_derivatives(-D2(div1_div2) / div1_div2^3)
    if simplify
        ddiv1_ddiv2 = Symbolics.simplify(ddiv1_ddiv2)
    end
    push!(eqs, ddiv2 ~ ddiv1_ddiv2) # e.g. https://math.stackexchange.com/questions/249253/second-derivative-of-the-inverse-function # TODO: higher orders

    # 4) If requested, instead remove and insert dummy equations
    if !dummies
        dummyidxs = findall(eq -> isequal(eq.lhs, div2) || isequal(eq.lhs, ddiv2), eqs)
        dummyeqs = splice!(eqs, dummyidxs) # return and remove dummy equations
        dummysubs = Dict(eq.lhs => eq.rhs for eq in dummyeqs)
        eqs = substitute.(eqs, Ref(dummysubs)) # don't iterate over dummysubs
    end

    # 5) Recreate system with new equations
    sys2 = typeof(sys)(eqs, iv2; name = nameof(sys), description = description(sys), kwargs...)
    return sys2
end
