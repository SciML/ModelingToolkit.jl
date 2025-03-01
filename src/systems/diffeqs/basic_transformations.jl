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
using ModelingToolkit, OrdinaryDiffEq, Test

@independent_variables t
@parameters α β γ δ
@variables x(t) y(t)
D = Differential(t)

eqs = [D(x) ~ α*x - β*x*y,
       D(y) ~ -δ*y + γ*x*y]

sys = ODESystem(eqs)
sys2 = liouville_transform(sys)
@variables trJ

u0 = [x => 1.0,
      y => 1.0,
      trJ => 1.0]

prob = ODEProblem(complete(sys2),u0,tspan,p)
sol = solve(prob,Tsit5())
```

Where `sol[3,:]` is the evolution of `trJ` over time.

Sources:

Probabilistic Robustness Analysis of F-16 Controller Performance: An
Optimal Transport Approach

Abhishek Halder, Kooktae Lee, and Raktim Bhattacharya
https://abhishekhalder.bitbucket.io/F16ACC2013Final.pdf
"""
function liouville_transform(sys::AbstractODESystem)
    t = get_iv(sys)
    @variables trJ
    D = ModelingToolkit.Differential(t)
    neweq = D(trJ) ~ trJ * -tr(calculate_jacobian(sys))
    neweqs = [equations(sys); neweq]
    vars = [unknowns(sys); trJ]
    ODESystem(neweqs, t, vars, parameters(sys), checks = false)
end

# TODO: handle case when new iv is a variable already in the system
function change_independent_variable(sys::AbstractODESystem, iv, iv1_of_iv2, iv2_of_iv1; kwargs...)
    iv1 = ModelingToolkit.get_iv(sys) # old independent variable
    iv2 = iv # new independent variable

    name2 = nameof(iv2)
    iv2func, = @variables $name2(iv1)

    eqs = ModelingToolkit.get_eqs(sys) |> copy # don't modify original system
    vars = []
    for (i, eq) in enumerate(eqs)
        vars = Symbolics.get_variables(eq)
        for var1 in vars
            if Symbolics.iscall(var1) # skip e.g. constants
                name = nameof(operation(var1))
                var2, = @variables $name(iv2func)
                eq = substitute(eq, var1 => var2; fold = false)
            end
        end
        eq = expand_derivatives(eq) # expand out with chain rule to get d(iv2)/d(iv1)
        eq = substitute(eq, Differential(iv1)(iv2func) => Differential(iv1)(iv2_of_iv1)) # substitute in d(iv2)/d(iv1)
        eq = expand_derivatives(eq)
        eq = substitute(eq, iv2func => iv2) # make iv2 independent
        eq = substitute(eq, iv1 => iv1_of_iv2) # substitute any remaining old ivars
        eqs[i] = eq
    end

    sys2 = typeof(sys)(eqs, iv2; name = nameof(sys), description = description(sys), kwargs...)
    return sys2
end

function change_independent_variable(sys::AbstractODESystem, iv; kwargs...)
    iv1 = get_iv(sys) # e.g. t
    iv2func = iv # e.g. a(t)
    iv2name = nameof(operation(unwrap(iv))) # TODO: handle namespacing?
    iv2, = @independent_variables $iv2name

    # TODO: find iv2func in system and replace it with some dummy variable
    # TODO: not just to 1st, but to all orders

    div2name = Symbol(iv2name, :_t)
    div2, = @variables $div2name(iv2)

    eqs = ModelingToolkit.get_eqs(sys) |> copy # don't modify original system
    vars = []
    for (i, eq) in enumerate(eqs)
        vars = Symbolics.get_variables(eq)
        for var1 in vars
            if Symbolics.iscall(var1) && !isequal(var1, iv2func) # && isequal(only(arguments(var1)), iv1) # skip e.g. constants
                name = nameof(operation(var1))
                var2, = @variables $name(iv2func)
                eq = substitute(eq, var1 => var2; fold = false)
            end
        end
        eq = expand_derivatives(eq) # expand out with chain rule to get d(iv2)/d(iv1)
        #eq = substitute(eq, Differential(iv1)(iv2func) => Differential(iv1)(iv2_of_iv1)) # substitute in d(iv2)/d(iv1)
        #eq = expand_derivatives(eq)
        eq = substitute(eq, Differential(iv1)(Differential(iv1)(iv2func)) => div2) # e.g. D(a(t)) => a_t(t) # TODO: more orders
        eq = substitute(eq, Differential(iv1)(iv2func) => div2) # e.g. D(a(t)) => a_t(t) # TODO: more orders
        eq = substitute(eq, iv2func => iv2) # make iv2 independent
        #eq = substitute(eq, iv1 => iv1_of_iv2) # substitute any remaining old ivars
        eqs[i] = eq
        println(eq)
    end

    sys2 = typeof(sys)(eqs, iv2; name = nameof(sys), description = description(sys), kwargs...)
    return sys2
end
