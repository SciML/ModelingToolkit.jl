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
@named sys = System(eqs, t)

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
function liouville_transform(sys::System; kwargs...)
    t = get_iv(sys)
    @variables trJ
    D = Differential(t)
    neweq = D(trJ) ~ trJ * -tr(calculate_jacobian(sys))
    neweqs = [equations(sys); neweq]
    vars = [unknowns(sys); trJ]
    System(
        neweqs, t, vars, parameters(sys);
        checks = false, name = nameof(sys), kwargs...
    )
end

"""
$(TYPEDSIGNATURES)

Generates the set of ODEs after change of variables.


Example:

```julia
using ModelingToolkit, OrdinaryDiffEq, Test

# Change of variables: z = log(x)
# (this implies that x = exp(z) is automatically non-negative)

@independent_variables t
@parameters α
@variables x(t)
D = Differential(t)
eqs = [D(x) ~ α*x]

tspan = (0., 1.)
def = [x => 1.0, α => -0.5]

@mtkcompile sys = System(eqs, t;defaults=def)
prob = ODEProblem(sys, [], tspan)
sol = solve(prob, Tsit5())

@variables z(t)
forward_subs  = [log(x) => z]
backward_subs = [x => exp(z)]
new_sys = change_of_variables(sys, t, forward_subs, backward_subs)
@test equations(new_sys)[1] == (D(z) ~ α)

new_prob = ODEProblem(new_sys, [], tspan)
new_sol = solve(new_prob, Tsit5())

@test isapprox(new_sol[x][end], sol[x][end], atol=1e-4)
```

"""
function change_of_variables(
        sys::System, iv, forward_subs, backward_subs;
        simplify = true, t0 = missing, isSDE = false
)
    t = iv

    old_vars = first.(backward_subs)
    new_vars = last.(forward_subs)

    # use: f = Y(t, X)
    # use: dY = (∂f/∂t + μ∂f/∂x + (1/2)*σ^2*∂2f/∂x2)dt + σ∂f/∂xdW
    old_eqs = equations(sys)
    neqs = get_noise_eqs(sys)
    brownvars = brownians(sys)

    if neqs === nothing && length(brownvars) === 0
        neqs = ones(1, length(old_eqs))
    elseif neqs !== nothing
        isSDE = true
        neqs = [neqs[i, :] for i in 1:size(neqs, 1)]

        brownvars = map([Symbol(:B, :_, i) for i in 1:length(neqs[1])]) do name
            unwrap(only(@brownians $name))
        end
    else
        isSDE = true
        neqs = Vector{Any}[]
        for (i, eq) in enumerate(old_eqs)
            neq = Any[]
            right = eq.rhs
            for Bv in brownvars
                lin_exp = linear_expansion(right, Bv)
                right = lin_exp[2]
                push!(neq, lin_exp[1])
            end
            push!(neqs, neq)
            old_eqs[i] = eq.lhs ~ right
        end
    end

    # df/dt = ∂f/∂x dx/dt + ∂f/∂t
    dfdt = Symbolics.derivative(first.(forward_subs), t)
    ∂f∂x = [Symbolics.derivative(first(f_sub), old_var)
            for (f_sub, old_var) in zip(forward_subs, old_vars)]
    ∂2f∂x2 = Symbolics.derivative.(∂f∂x, old_vars)
    new_eqs = Equation[]

    for (new_var, ex, first, second) in zip(new_vars, dfdt, ∂f∂x, ∂2f∂x2)
        for (eqs, neq) in zip(old_eqs, neqs)
            if occursin(value(eqs.lhs), value(ex))
                ex = substitute(ex, eqs.lhs => eqs.rhs)
                if isSDE
                    for (noise, B) in zip(neq, brownvars)
                        ex = ex + 1/2 * noise^2 * second + noise*first*B
                    end
                end
            end
        end
        ex = substitute(ex, Dict(forward_subs))
        ex = substitute(ex, Dict(backward_subs))
        if simplify
            ex = Symbolics.simplify(ex, expand = true)
        end
        push!(new_eqs, Differential(t)(new_var) ~ ex)
    end

    defs = get_defaults(sys)
    new_defs = Dict()
    for f_sub in forward_subs
        ex = substitute(first(f_sub), defs)
        if !ismissing(t0)
            ex = substitute(ex, t => t0)
        end
        new_defs[last(f_sub)] = ex
    end
    for para in parameters(sys)
        if haskey(defs, para)
            new_defs[para] = defs[para]
        end
    end

    @named new_sys = System(
        vcat(new_eqs, first.(backward_subs) .~ last.(backward_subs)), t;
        defaults = new_defs,
        observed = observed(sys)
    )
    if simplify
        return mtkcompile(new_sys)
    end
    return new_sys
end

"""
Generates the system of ODEs to find solution to FDEs.

Example:

```julia
@independent_variables t
@variables x(t)
D = Differential(t)
tspan = (0., 1.)

α = 0.5
eqs = (9*gamma(1 + α)/4) - (3*t^(4 - α/2)*gamma(5 + α/2)/gamma(5 - α/2))
eqs += (gamma(9)*t^(8 - α)/gamma(9 - α)) + (3/2*t^(α/2)-t^4)^3 - x^(3/2)
sys = fractional_to_ordinary(eqs, x, α, 10^-7, 1)

prob = ODEProblem(sys, [], tspan)
sol = solve(prob, radau5(), abstol = 1e-10, reltol = 1e-10)
```
"""
function fractional_to_ordinary(
        eqs, variables, alphas, epsilon, T;
        initials = 0, additional_eqs = [], iv = only(@independent_variables t), matrix=false
)
    D = Differential(iv)
    i = 0
    all_eqs = Equation[]
    all_def = Pair[]
    
    function fto_helper(sub_eq, sub_var, α; initial=0)
        alpha_0 = α

        if (α > 1)
            coeff = 1/(α - 1)
            m = 2
            while (α - m > 0)
                coeff /= α - m
                m += 1
            end
            alpha_0 = α - m + 1
        end

        δ = (gamma(alpha_0+1) * epsilon)^(1/alpha_0)
        a = pi/2*(1-(1-alpha_0)/((2-alpha_0) * log(epsilon^-1)))
        h = 2*pi*a / log(1 + (2/epsilon * (cos(a))^(alpha_0 - 1)))

        x_sub = (gamma(2-alpha_0) * epsilon)^(1/(1-alpha_0))
        x_sup = -log(gamma(1-alpha_0) * epsilon)
        M = floor(Int, log(x_sub / T) / h)
        N = ceil(Int, log(x_sup / δ) / h)

        function c_i(index)
            h * sin(pi * alpha_0) / pi * exp((1-alpha_0)*h*index)
        end

        function γ_i(index)
            exp(h * index)
        end

        new_eqs = Equation[]
        def = Pair[]

        if matrix
            new_z = Symbol(:ʐ, :_, i)
            i += 1
            γs = diagm([γ_i(index) for index in M:N-1])
            cs = [c_i(index) for index in M:N-1]

            if (α < 1)
                new_z = only(@variables $new_z(iv)[1:N-M])
                new_eq = D(new_z) ~ -γs*new_z .+ sub_eq
                rhs = dot(cs, new_z) + initial
                push!(def, new_z=>zeros(N-M))
            else
                new_z = only(@variables $new_z(iv)[1:N-M, 1:m])
                new_eq = D(new_z) ~ -γs*new_z + hcat(fill(sub_eq, N-M, 1), collect(new_z[:, 1:m-1]*diagm(1:m-1)))
                rhs = coeff*sum(cs[i]*new_z[i, m] for i in 1:N-M)
                for (index, value) in enumerate(initial)
                    rhs += value * iv^(index - 1) / gamma(index)
                end
                push!(def, new_z=>zeros(N-M, m))
            end
            push!(new_eqs, new_eq)
        else
            if (α < 1)  
                rhs = initial
                for index in range(M, N-1; step=1)
                    new_z = Symbol(:ʐ, :_, i)
                    i += 1
                    new_z = ModelingToolkit.unwrap(only(@variables $new_z(iv)))
                    new_eq = D(new_z) ~ sub_eq - γ_i(index)*new_z
                    push!(new_eqs, new_eq)
                    push!(def, new_z=>0)
                    rhs += c_i(index)*new_z
                end
            else
                rhs = 0
                for (index, value) in enumerate(initial)
                    rhs += value * iv^(index - 1) / gamma(index)
                end
                for index in range(M, N-1; step=1)
                    new_z = Symbol(:ʐ, :_, i)
                    i += 1
                    γ = γ_i(index)
                    base = sub_eq
                    for k in range(1, m; step=1)
                        new_z = Symbol(:ʐ, :_, index-M, :_, k)
                        new_z = ModelingToolkit.unwrap(only(@variables $new_z(iv)))
                        new_eq = D(new_z) ~ base - γ*new_z
                        base = k * new_z
                        push!(new_eqs, new_eq)
                        push!(def, new_z=>0)
                    end
                    rhs += coeff*c_i(index)*new_z
                end
            end
        end
        push!(new_eqs, sub_var ~ rhs)
        return (new_eqs, def)
    end

    for (eq, cur_var, alpha, init) in zip(eqs, variables, alphas, initials)
        (new_eqs, def) = fto_helper(eq, cur_var, alpha; initial=init)
        append!(all_eqs, new_eqs)
        append!(all_def, def)
    end
    append!(all_eqs, additional_eqs)
    @named sys = System(all_eqs, iv; defaults=all_def)
    return mtkcompile(sys)
end

"""
Generates the system of ODEs to find solution to FDEs.

Example:

```julia
@independent_variables t
@variables x_0(t)
D = Differential(t)
tspan = (0., 5000.)

function expect(t)
    return sqrt(2) * sin(t + pi/4)
end

sys = linear_fractional_to_ordinary([3, 2.5, 2, 1, .5, 0], [1, 1, 1, 4, 1, 4], 6*cos(t), 10^-5, 5000; initials=[1, 1, -1])
prob = ODEProblem(sys, [], tspan)
sol = solve(prob, radau5(), abstol = 1e-5, reltol = 1e-5)
```
"""
function linear_fractional_to_ordinary(
        degrees, coeffs, rhs, epsilon, T;
        initials = 0, symbol = :x, iv = only(@independent_variables t), matrix=false
)
    previous = Symbol(symbol, :_, 0)
    previous = ModelingToolkit.unwrap(only(@variables $previous(iv)))
    @variables x_0(iv)
    D = Differential(iv)
    i = 0
    all_eqs = Equation[]
    all_def = Pair[]

    function fto_helper(sub_eq, α)
        δ = (gamma(α+1) * epsilon)^(1/α)
        a = pi/2*(1-(1-α)/((2-α) * log(epsilon^-1)))
        h = 2*pi*a / log(1 + (2/epsilon * (cos(a))^(α - 1)))

        x_sub = (gamma(2-α) * epsilon)^(1/(1-α))
        x_sup = -log(gamma(1-α) * epsilon)
        M = floor(Int, log(x_sub / T) / h)
        N = ceil(Int, log(x_sup / δ) / h)

        function c_i(index)
            h * sin(pi * α) / pi * exp((1-α)*h*index)
        end

        function γ_i(index)
            exp(h * index)
        end

        new_eqs = Equation[]
        def = Pair[]
        if matrix
            new_z = Symbol(:ʐ, :_, i)
            i += 1
            γs = diagm([γ_i(index) for index in M:N-1])
            cs = [c_i(index) for index in M:N-1]

            new_z = only(@variables $new_z(iv)[1:N-M])
            new_eq = D(new_z) ~ -γs*new_z .+ sub_eq
            sum = dot(cs, new_z)
            push!(def, new_z=>zeros(N-M))
            push!(new_eqs, new_eq)
        else
            sum = 0
            for index in range(M, N-1; step=1)
                new_z = Symbol(:ʐ, :_, i)
                i += 1
                new_z = ModelingToolkit.unwrap(only(@variables $new_z(iv)))
                new_eq = D(new_z) ~ sub_eq - γ_i(index)*new_z
                push!(new_eqs, new_eq)
                push!(def, new_z=>0)
                sum += c_i(index)*new_z
            end
        end
        return (new_eqs, def, sum)
    end

    for i in range(1, ceil(Int, degrees[1]); step=1)
        new_x = Symbol(symbol, :_, i)
        new_x = ModelingToolkit.unwrap(only(@variables $new_x(iv)))
        push!(all_eqs, D(previous) ~ new_x)
        push!(all_def, previous => initials[i])
        previous = new_x
    end

    new_rhs = -rhs
    for (degree, coeff) in zip(degrees, coeffs)
        rounded = ceil(Int, degree)
        new_x = Symbol(symbol, :_, rounded)
        new_x = ModelingToolkit.unwrap(only(@variables $new_x(iv)))
        if isinteger(degree)
            new_rhs += coeff * new_x
        else
            (new_eqs, def, sum) = fto_helper(new_x, rounded - degree)
            append!(all_eqs, new_eqs)
            append!(all_def, def)
            new_rhs += coeff * sum
        end
    end
    push!(all_eqs, 0 ~ new_rhs)
    @named sys = System(all_eqs, iv; defaults=all_def)
    return mtkcompile(sys)
end

"""
    change_independent_variable(
        sys::System, iv, eqs = [];
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
In following calls to `mtkcompile`, consider passing `allow_symbolic = true` to avoid undesired constraint equations between between dummy variables.

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

julia> @named M = System([D(D(y)) ~ -9.81, D(D(x)) ~ 0.0], t);

julia> M = change_independent_variable(M, x);

julia> M = mtkcompile(M; allow_symbolic = true);

julia> unknowns(M)
3-element Vector{SymbolicUtils.BasicSymbolic{Real}}:
 xˍt(x)
 y(x)
 yˍx(x)
```
"""
function change_independent_variable(
        sys::System, iv, eqs = [];
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

    # A utility function that returns whether var (e.g. f(t)) is a function of iv (e.g. t)
    function is_function_of(var, iv)
        # Peel off outer calls to find the argument of the function of
        if iscall(var) && operation(var) === getindex # handle array variables
            var = arguments(var)[1] # (f(t))[1] -> f(t)
        end
        if iscall(var)
            var = only(arguments(var)) # e.g. f(t) -> t
            return isequal(var, iv)
        end
        return false
    end

    # Create a utility that performs the chain rule on an expression, followed by insertion of the new independent variable:
    # e.g. (d/dt)(f(t)) -> (d/dt)(f(u(t))) -> df(u(t))/du(t) * du(t)/dt -> df(u)/du * uˍt(u)
    function transform(ex::T) where {T}
        # 1) Replace the argument of every function; e.g. f(t) -> f(u(t))
        for var in vars(ex; op = Nothing) # loop over all variables in expression (op = Nothing prevents interpreting "D(f(t))" as one big variable)
            if is_function_of(var, iv1) && !isequal(var, iv2_of_iv1) # of the form f(t)? but prevent e.g. u(t) -> u(u(t))
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
        return ex::T
    end

    # overload to specifically handle equations, which can be an equation or a connection
    function transform(eq::Equation, systems_map)
        if eq.rhs isa Connection
            eq = connect((systems_map[nameof(s)] for s in eq.rhs.systems)...)
        else
            eq = transform(eq)
        end
        return eq::Equation
    end

    # Use the utility function to transform everything in the system!
    function transform(sys::System)
        systems = map(transform, get_systems(sys)) # recurse through subsystems
        # transform equations and connections
        systems_map = Dict(get_name(s) => s for s in systems)
        eqs = map(eq -> transform(eq, systems_map)::Equation, get_eqs(sys))
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
        connector_type = get_connector_type(sys)
        assertions = Dict(transform(ass) => msg for (ass, msg) in get_assertions(sys))
        wascomplete = iscomplete(sys) # save before reconstructing system
        wassplit = is_split(sys)
        wasflat = isempty(systems)
        sys = typeof(sys)( # recreate system with transformed fields
            eqs, iv2, unknowns, ps; observed, initialization_eqs,
            defaults, guesses, connector_type,
            assertions, name = nameof(sys), description = description(sys)
        )
        sys = compose(sys, systems) # rebuild hierarchical system
        if wascomplete
            sys = complete(sys; split = wassplit, flatten = wasflat) # complete output if input was complete
            @set! sys.parameter_dependencies = parameter_dependencies
        end
        return sys
    end
    return transform(sys)
end

"""
$(TYPEDSIGNATURES)

Choose correction_factor=-1//2 (1//2) to convert Ito -> Stratonovich (Stratonovich->Ito).
"""
function stochastic_integral_transform(sys::System, correction_factor)
    if !isempty(get_systems(sys))
        throw(ArgumentError("The system must be flattened."))
    end
    if get_noise_eqs(sys) === nothing
        throw(ArgumentError("""
        `$stochastic_integral_transform` expects a system with noise_eqs. If your \
        noise is specified using brownian variables, consider calling \
        `mtkcompile`.
        """))
    end
    name = nameof(sys)
    noise_eqs = get_noise_eqs(sys)
    eqs = equations(sys)
    dvs = unknowns(sys)
    ps = parameters(sys)
    # use the general interface
    if noise_eqs isa Vector
        _eqs = reduce(vcat, [eqs[i].lhs ~ noise_eqs[i] for i in eachindex(dvs)])
        de = System(_eqs, get_iv(sys), dvs, ps, name = name, checks = false)

        jac = calculate_jacobian(de, sparse = false, simplify = false)
        ∇σσ′ = simplify.(jac * noise_eqs)
    else
        dimunknowns, m = size(noise_eqs)
        _eqs = reduce(vcat, [eqs[i].lhs ~ noise_eqs[i] for i in eachindex(dvs)])
        de = System(_eqs, get_iv(sys), dvs, ps, name = name, checks = false)

        jac = calculate_jacobian(de, sparse = false, simplify = false)
        ∇σσ′ = simplify.(jac * noise_eqs[:, 1])
        for k in 2:m
            __eqs = reduce(vcat,
                [eqs[i].lhs ~ noise_eqs[Int(i + (k - 1) * dimunknowns)]
                 for i in eachindex(dvs)])
            de = System(__eqs, get_iv(sys), dvs, dvs, name = name, checks = false)

            jac = calculate_jacobian(de, sparse = false, simplify = false)
            ∇σσ′ = ∇σσ′ + simplify.(jac * noise_eqs[:, k])
        end
    end
    deqs = reduce(vcat,
        [eqs[i].lhs ~ eqs[i].rhs + correction_factor * ∇σσ′[i] for i in eachindex(dvs)])

    # reduce(vcat, [1]) == 1 for some reason
    if deqs isa Equation
        deqs = [deqs]
    end
    return @set sys.eqs = deqs
end

"""
$(TYPEDSIGNATURES)

Measure transformation method that allows for a reduction in the variance of an estimator `Exp(g(X_t))`.
Input:  Original SDE system and symbolic function `u(t,x)` with scalar output that
        defines the adjustable parameters `d` in the Girsanov transformation. Optional: initial
        condition for `θ0`.
Output: Modified SDE System with additional component `θ_t` and initial value `θ0`, as well as
        the weight `θ_t/θ0` as observed equation, such that the estimator `Exp(g(X_t)θ_t/θ0)`
        has a smaller variance.

Reference:
Kloeden, P. E., Platen, E., & Schurz, H. (2012). Numerical solution of SDE through computer
experiments. Springer Science & Business Media.

# Example

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters α β
@variables x(t) y(t) z(t)

eqs = [D(x) ~ α*x]
noiseeqs = [β*x]

@named de = System(eqs,t,[x],[α,β]; noise_eqs = noiseeqs)

# define u (user choice)
u = x
θ0 = 0.1
g(x) = x[1]^2
demod = ModelingToolkit.Girsanov_transform(de, u; θ0=0.1)

u0modmap = [
    x => x0
]

parammap = [
    α => 1.5,
    β => 1.0
]

probmod = SDEProblem(complete(demod),u0modmap,(0.0,1.0),parammap)
ensemble_probmod = EnsembleProblem(probmod;
          output_func = (sol,i) -> (g(sol[x,end])*sol[demod.weight,end],false),
          )

simmod = solve(ensemble_probmod,EM(),dt=dt,trajectories=numtraj)
```

"""
function Girsanov_transform(sys::System, u; θ0 = 1.0)
    name = nameof(sys)

    # register new variable θ corresponding to 1D correction process θ(t)
    t = get_iv(sys)
    D = Differential(t)
    @variables θ(t), weight(t)

    # determine the adjustable parameters `d` given `u`
    # gradient of u with respect to unknowns
    grad = Symbolics.gradient(u, unknowns(sys))

    noiseeqs = copy(get_noise_eqs(sys))
    if noiseeqs isa Vector
        d = simplify.(-(noiseeqs .* grad) / u)
        drift_correction = noiseeqs .* d
    else
        d = simplify.(-noiseeqs * grad / u)
        drift_correction = noiseeqs * d
    end

    eqs = equations(sys)
    dvs = unknowns(sys)
    # transformation adds additional unknowns θ: newX = (X,θ)
    # drift function for unknowns is modified
    # θ has zero drift
    deqs = reduce(
        vcat, [eqs[i].lhs ~ eqs[i].rhs - drift_correction[i] for i in eachindex(dvs)])
    if deqs isa Equation
        deqs = [deqs]
    end
    deqsθ = D(θ) ~ 0
    push!(deqs, deqsθ)

    # diffusion matrix is of size d x m (d unknowns, m noise), with diagonal noise represented as a d-dimensional vector
    # for diagonal noise processes with m>1, the noise process will become non-diagonal; extra unknown component but no new noise process.
    # new diffusion matrix is of size d+1 x M
    # diffusion for state is unchanged

    noiseqsθ = θ * d

    if noiseeqs isa Vector
        m = size(noiseeqs)
        if m == 1
            push!(noiseeqs, noiseqsθ)
        else
            noiseeqs = [Array(Diagonal(wrap.(noiseeqs))); noiseqsθ']
        end
    else
        noiseeqs = [Array(noiseeqs); noiseqsθ']
    end

    unknown_vars = [dvs; θ]

    # return modified SDE System
    @set! sys.eqs = deqs
    @set! sys.noise_eqs = noiseeqs
    @set! sys.unknowns = unknown_vars
    get_defaults(sys)[θ] = θ0
    obs = observed(sys)
    @set! sys.observed = [weight ~ θ / θ0; obs]
    if get_parent(sys) !== nothing
        @set! sys.parent.unknowns = [get_unknowns(get_parent(sys)); [θ, weight]]
    end
    return sys
end

"""
    $(TYPEDSIGNATURES)

Add accumulation variables for `vars`. For every unknown `x` in `vars`, add
`D(accumulation_x) ~ x` as an equation.
"""
function add_accumulations(sys::System, vars = unknowns(sys))
    avars = [rename(v, Symbol(:accumulation_, getname(v))) for v in vars]
    return add_accumulations(sys, avars .=> vars)
end

"""
    $(TYPEDSIGNATURES)

Add accumulation variables for `vars`. `vars` is a vector of pairs in the form
of

```julia
[cumulative_var1 => x + y, cumulative_var2 => x^2]
```
Then, cumulative variables `cumulative_var1` and `cumulative_var2` that computes
the cumulative `x + y` and `x^2` would be added to `sys`.

All accumulation variables have a default of zero.
"""
function add_accumulations(sys::System, vars::Vector{<:Pair})
    eqs = get_eqs(sys)
    avars = map(first, vars)
    ints = intersect(avars, unknowns(sys))
    if !isempty(ints)
        error("$ints already exist in the system!")
    end
    D = Differential(get_iv(sys))
    @set! sys.eqs = [eqs; Equation[D(a) ~ v[2] for (a, v) in zip(avars, vars)]]
    @set! sys.unknowns = [get_unknowns(sys); avars]
    @set! sys.defaults = merge(get_defaults(sys), Dict(a => 0.0 for a in avars))
    return sys
end

"""
    $(TYPEDSIGNATURES)

Given a system with noise in the form of noise equation (`get_noise_eqs(sys) !== nothing`)
return an equivalent system which represents the noise using brownian variables.

# Keyword Arguments

- `names`: The name(s) to use for the brownian variables. If this is a `Symbol`, variables
  with the given name and successive numeric `_i` suffixes will be used. If a `Vector`,
  this must have appropriate length for the noise equations of the system. The
  corresponding number of brownian variables are created with the given names.
"""
function noise_to_brownians(sys::System; names::Union{Symbol, Vector{Symbol}} = :α)
    neqs = get_noise_eqs(sys)
    if neqs === nothing
        throw(ArgumentError("Expected a system with `noise_eqs`."))
    end
    if !isempty(get_systems(sys))
        throw(ArgumentError("The system must be flattened."))
    end
    # vector means diagonal noise
    nbrownians = ndims(neqs) == 1 ? length(neqs) : size(neqs, 2)
    if names isa Symbol
        names = [Symbol(names, :_, i) for i in 1:nbrownians]
    end
    if length(names) != nbrownians
        throw(ArgumentError("""
        The system has $nbrownians brownian variables. Received $(length(names)) names \
        for the brownian variables. Provide $nbrownians names or a single `Symbol` to use \
        an array variable of the appropriately length.
        """))
    end
    brownvars = map(names) do name
        unwrap(only(@brownians $name))
    end

    terms = if ndims(neqs) == 1
        neqs .* brownvars
    else
        neqs * brownvars
    end

    eqs = map(get_eqs(sys), terms) do eq, term
        eq.lhs ~ eq.rhs + term
    end

    @set! sys.eqs = eqs
    @set! sys.brownians = brownvars
    @set! sys.noise_eqs = nothing

    return sys
end

"""
    $(TYPEDSIGNATURES)

Function which takes a system `sys` and an independent variable `t` and changes the
independent variable of `sys` to `t`. This is different from
[`change_independent_variable`](@ref) since this function only does a symbolic substitution
of the independent variable. `sys` must not be a reduced system (`observed(sys)` must be
empty). If `sys` is time-independent, this can be used to turn it into a time-dependent
system.

# Keyword arguments

- `name`: The name of the returned system.
"""
function convert_system_indepvar(sys::System, t; name = nameof(sys))
    isempty(observed(sys)) ||
        throw(ArgumentError("""
        `convert_system_indepvar` cannot handle reduced model (i.e. observed(sys) is non-\
        empty).
        """))
    t = value(t)
    varmap = Dict()
    sts = unknowns(sys)
    newsts = similar(sts, Any)
    for (i, s) in enumerate(sts)
        if iscall(s)
            args = arguments(s)
            length(args) == 1 ||
                throw(InvalidSystemException("Illegal unknown: $s. The unknown can have at most one argument like `x(t)`."))
            arg = args[1]
            if isequal(arg, t)
                newsts[i] = s
                continue
            end
            ns = maketerm(typeof(s), operation(s), Any[t],
                SymbolicUtils.metadata(s))
            newsts[i] = ns
            varmap[s] = ns
        else
            ns = variable(getname(s); T = FnType)(t)
            newsts[i] = ns
            varmap[s] = ns
        end
    end
    sub = Base.Fix2(substitute, varmap)
    if is_time_dependent(sys)
        iv = only(independent_variables(sys))
        sub.x[iv] = t # otherwise the Differentials aren't fixed
    end
    neweqs = map(sub, equations(sys))
    defs = Dict(sub(k) => sub(v) for (k, v) in defaults(sys))
    neqs = get_noise_eqs(sys)
    if neqs !== nothing
        neqs = map(sub, neqs)
    end
    cstrs = map(sub, get_constraints(sys))
    costs = Vector{Union{Real, BasicSymbolic}}(map(sub, get_costs(sys)))
    @set! sys.eqs = neweqs
    @set! sys.iv = t
    @set! sys.unknowns = newsts
    @set! sys.defaults = defs
    @set! sys.name = name
    @set! sys.noise_eqs = neqs
    @set! sys.constraints = cstrs
    @set! sys.costs = costs

    var_to_name = Dict(k => get(varmap, v, v) for (k, v) in get_var_to_name(sys))
    @set! sys.var_to_name = var_to_name
    return sys
end

"""
    $(TYPEDSIGNATURES)

Shorthand for `respecialize(sys, []; all = true)`
"""
respecialize(sys::AbstractSystem) = respecialize(sys, []; all = true)

"""
    $(TYPEDSIGNATURES)

Specialize nonnumeric parameters in `sys` by changing their symtype to a concrete type.
`mapping` is an iterable, where each element can be a parameter or a pair mapping a parameter
to a value. If the element is a parameter, it must have a default. Each specified parameter
is updated to have the symtype of the value associated with it (either in `mapping` or in
the defaults). This operation can only be performed on nonnumeric, non-array parameters. The
defaults of respecialized parameters are set to the associated values.

This operation can only be performed on `complete`d systems.

# Keyword arguments

- `all`: Specialize all nonnumeric parameters in the system. This will error if any such
  parameter does not have a default.
"""
function respecialize(sys::AbstractSystem, mapping; all = false)
    if !iscomplete(sys)
        error("""
        This operation can only be performed on completed systems. Use `complete(sys)` or
        `mtkcompile(sys)`.
        """)
    end
    if !is_split(sys)
        error("""
        This operation can only be performed on split systems. Use `complete(sys)` or
        `mtkcompile(sys)` with the `split = true` keyword argument.
        """)
    end

    new_ps = copy(get_ps(sys))
    @set! sys.ps = new_ps

    extras = []
    if all
        for x in filter(!is_variable_numeric, get_ps(sys))
            if any(y -> isequal(x, y) || y isa Pair && isequal(x, y[1]), mapping) ||
               symbolic_type(x) === ArraySymbolic() ||
               iscall(x) && operation(x) === getindex
                continue
            end
            push!(extras, x)
        end
    end
    ps_to_specialize = Iterators.flatten((extras, mapping))

    defs = copy(defaults(sys))
    @set! sys.defaults = defs
    final_defs = copy(defs)
    evaluate_varmap!(final_defs, ps_to_specialize)

    subrules = Dict()

    for element in ps_to_specialize
        if element isa Pair
            k, v = element
        else
            k = element
            v = get(final_defs, k, nothing)
            @assert v !== nothing """
            Parameter $k needs an associated value to be respecialized.
            """
            @assert symbolic_type(v) == NotSymbolic() && !is_array_of_symbolics(v) """
            Parameter $k needs an associated value to be respecialized. Found symbolic \
            default $v.
            """
        end

        k = unwrap(k)
        T = typeof(v)

        @assert !is_variable_numeric(k) """
        Numeric types cannot be respecialized - tried to respecialize $k.
        """
        @assert symbolic_type(k) !== ArraySymbolic() """
        Cannot respecialize array symbolics - tried to respecialize $k.
        """
        @assert !iscall(k) || operation(k) !== getindex """
        Cannot respecialized scalarized array variables - tried to respecialize $k.
        """
        idx = findfirst(isequal(k), get_ps(sys))
        @assert idx !== nothing """
        Parameter $k does not exist in the system.
        """

        if iscall(k)
            op = operation(k)::BasicSymbolic
            @assert !iscall(op)
            op = SymbolicUtils.Sym{SymbolicUtils.FnType{Tuple{Any}, T}}(nameof(op))
            args = arguments(k)
            new_p = op(args...)
        else
            new_p = SymbolicUtils.Sym{T}(getname(k))
        end

        get_ps(sys)[idx] = new_p
        defaults(sys)[new_p] = v
        subrules[unwrap(k)] = unwrap(new_p)
    end

    substituter = Base.Fix2(fast_substitute, subrules)
    @set! sys.eqs = map(substituter, get_eqs(sys))
    @set! sys.observed = map(substituter, get_observed(sys))
    @set! sys.initialization_eqs = map(substituter, get_initialization_eqs(sys))
    if get_noise_eqs(sys) !== nothing
        @set! sys.noise_eqs = map(substituter, get_noise_eqs(sys))
    end
    @set! sys.assertions = Dict([substituter(k) => v for (k, v) in assertions(sys)])
    @set! sys.parameter_dependencies = map(substituter, get_parameter_dependencies(sys))
    @set! sys.defaults = Dict([substituter(k) => substituter(v) for (k, v) in defaults(sys)])
    @set! sys.guesses = Dict([k => substituter(v) for (k, v) in guesses(sys)])
    @set! sys.continuous_events = map(get_continuous_events(sys)) do cev
        SymbolicContinuousCallback(
            map(substituter, cev.conditions), substituter(cev.affect),
            substituter(cev.affect_neg), substituter(cev.initialize),
            substituter(cev.finalize), cev.rootfind,
            cev.reinitializealg, cev.zero_crossing_id)
    end
    @set! sys.discrete_events = map(get_discrete_events(sys)) do dev
        SymbolicDiscreteCallback(map(substituter, dev.conditions), substituter(dev.affect),
            substituter(dev.initialize), substituter(dev.finalize), dev.reinitializealg)
    end
    if get_schedule(sys) !== nothing
        sched = get_schedule(sys)
        @set! sys.schedule = Schedule(
            sched.var_sccs, AnyDict(k => substituter(v) for (k, v) in sched.dummy_sub))
    end
    @set! sys.constraints = map(substituter, get_constraints(sys))
    @set! sys.tstops = map(substituter, get_tstops(sys))
    @set! sys.costs = Vector{Union{Real, BasicSymbolic}}(map(substituter, get_costs(sys)))
    sys = complete(sys; split = is_split(sys))
    return sys
end
