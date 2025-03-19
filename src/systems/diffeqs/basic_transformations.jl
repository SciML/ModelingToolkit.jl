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

"""
$(TYPEDSIGNATURES)

Transforms an ODE system by applying a set of variable transformations and their inverses (back-transformations). 
Returns the new set of differential equations and the updated list of variables after substitution and simplification. This is useful for redefining an ODE system in terms of new dependent variables while preserving the dynamics.

For example, transforming `x` to `z = exp(x)` redefines the system in terms of `z`, adjusting the equations accordingly.

Example:

```julia
using ModelingToolkit

@parameters t α
@variables x(t) z(t)
D = Differential(t)
eqs = [D(x) ~ α*x]
sys = ODESystem(eqs, t, [x], [α])
transformations = [z => exp(x)]
back_transformations = [x => log(z)]
new_eqs, new_vars = transform_eqs(sys, transformations, back_transformations)
# new_eqs: [D(z) ~ α*z], new_vars: [z]
"""


function transform_eqs(sys::ODESystem, transformations, back_transformations)
    all_dept_vars = unknowns(sys)
    diff_eqs = equations(sys)
    diff_eqs_map = Dict(eq.lhs => eq.rhs for eq in diff_eqs)
    eq_transformations = [k ~ v for (k, v) in transformations]

    new_vars = [] # transformed variables
    mutable_vars = [] # mutable variables
    for tr in transformations
        push!(new_vars, tr.first)
        change_var = get_variables(tr.second)

        # change_var can also be a parameter like α. hence the check. 
        for every_change_var in change_var
            for all_var in all_dept_vars
                if isequal(every_change_var,all_var)
                    push!(mutable_vars, all_var)
                end
            end
        end
    end

    current_vars = setdiff(vcat(all_dept_vars, new_vars), mutable_vars)

    # Differentiating on both sides of the transformation
    dzdt = Symbolics.derivative(last.(transformations), t)
    dz = Symbolics.derivative(first.(transformations), t)
    new_eqs = Equation[]
    for (dzi, dzdti) in zip(dz, dzdt)
        updated_dzdti = substitute(dzdti, diff_eqs_map)
        exp = dzi ~ updated_dzdti
        updated_exp = substitute(exp, Dict(back_transformations))
        push!(new_eqs, Symbolics.simplify(updated_exp))
    end

    return new_eqs, current_vars
end

### `substitute_defaults`


"""
$(TYPEDSIGNATURES)

Substitutes default values into transformation expressions for an ODE system, returning a dictionary of new default values for transformed variables. Optionally evaluates at a specific time `t_val`. Useful for initializing transformed systems with consistent initial conditions.

For example, if `x => 1.0` is a default and `z = exp(x)` is a transformation, it computes `z => exp(1.0)`.

Example:

```julia
using ModelingToolkit

@parameters t α
@variables x(t) z(t)
eqs = [Differential(t)(x) ~ α*x]
sys = ODESystem(eqs, t, [x], [α]; defaults=Dict(x => 1.0))
transformations = [z => exp(x)]
new_defaults = substitute_defaults(sys, transformations)
# new_defaults: Dict(z => 2.718281828459045)
"""

function substitute_defaults(sys::ODESystem, transformations, t_val=missing)
    existing_defaults = ModelingToolkit.get_defaults(sys)
    new_defs = Dict()
    for tr in transformations
        ex = substitute(last(tr), existing_defaults)
        if !ismissing(t_val)
            ex = substitute(ex, t => t_val)
        end
        vars = get_variables(ex)
        if length(vars) == 0
            new_defs[first(tr)] = ex
        end
    end

    return new_defs
end
