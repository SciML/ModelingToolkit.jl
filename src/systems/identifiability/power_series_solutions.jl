using LinearAlgebra, Symbolics, SymbolicUtils
using Symbolics:value
using Primes
using Nemo

# struct ODEIdProblem
#     state_equations
#     output_equations
#     states
#     inputs
#     parameters
#     polynomial_ring
# end

function switch_ring(var, ring)
    ind = findfirst(vv -> vv == var, gens(parent(var)))
    var_s = symbols(parent(var))[ind]
    ind = findfirst(vv -> (vv == var_s), symbols(ring))
    if ind == nothing
        throw(Base.KeyError("Variable $var_s is not found in $ring"))
    end
    return gens(ring)[ind]
end

function get_numerator(f)
    if f isa Generic.Frac
        return numerator(f)
    elseif f isa MPolyElem
        return f
    end
end

function get_denominator(f)
    if f isa Generic.Frac
        return denominator(f)
    elseif f isa MPolyElem
        return denominator(f)# one(parent(f))
    end
end

function get_height_and_coeff(f)
    if isequal(length(f), 0)
        return (0, 1)
    end
    max_coef = 1
    for c in coeffs(f)
        max_coef = max(max_coef, 2 * height_bits(c))
    end
    return (total_degree(f), max_coef) 
end

function numerator_height_coeff(f)
    return get_height_and_coeff(get_numerator(f))
end

function denominator_height_coeff(f)
    return get_height_and_coeff(get_denominator(f)) 
end


""" function PreprocessODE(ODE, outputs, x, )
Accepts ODE as array [D(x)~f(x,u,y)]
"""
function PreprocessODE(ODE, outputs, x, y, u, θ, t)
    D = Differential(t) 
    @parameters ẋ[1:length(x)]

    input_symbols = vcat(ẋ, x, u, y, θ)
    generators = string.(input_symbols)
    R, gens_ = Nemo.PolynomialRing(Nemo.QQ, generators)
    
    state_eqs = [substitute(eqn.lhs, D.(x) .=> ẋ) ~ substitute(eqn.rhs, D.(x) .=> ẋ) for eqn in ODE];
    out_eqs = [substitute(value(eqn.lhs - eqn.rhs), input_symbols .=> gens_) for eqn in outputs]

    state_eqn_dict = Dict(substitute(value(eqn.lhs), input_symbols .=> gens_) => substitute(value(eqn.rhs), input_symbols .=> gens_) for eqn in state_eqs)
    state_eqs = [k - v for (k, v) in state_eqn_dict]

    states = [substitute(value(each),  input_symbols .=> gens_) for each in x]
    params = [substitute(value(each),  input_symbols .=> gens_) for each in θ]
    inputs = [substitute(value(each),  input_symbols .=> gens_) for each in u]
    outputs = [substitute(value(each),  input_symbols .=> gens_) for each in y] 

    return (state_eqs, out_eqs, states, params, inputs, outputs, state_eqn_dict)
end

function Initialize(state_eqs, out_eqs, states, inputs, outputs, params, proba)
    # Proposition 3.3 in https://doi.org/10.1006/jsco.2002.0532
    d, h = 1, 1
    for f in vcat(state_eqs, out_eqs)
        numer_df, numer_hf = numerator_height_coeff(f)
        denom_df, denom_hf = denominator_height_coeff(f)
        d = max(d, max(numer_df, denom_df))
        h = max(h, max(numer_hf, denom_hf))
    end
    p_per_func = 1 - (1 - proba) / (length(states) + length(params)) # TODO: add capability for checking different functions
    mu = ceil(1 / (1 - sqrt(p_per_func)))
    solution = Dict()
    n = length(states)
    m = length(outputs)
    r = length(inputs)
    ℓ = length(params)
    D = 4 * (n + ℓ)^2 * (n + m) * d
    Dprime = D * (2 * log(n + ℓ + r + 1) + log(mu * D)) + 4 * (n + ℓ)^2 * ((n + m) * h + log(2 * n * D))
    prime_number = Primes.nextprime(Int(ceil(2 * mu * Dprime)))
    F = Nemo.GF(prime_number)
    prec = n + ℓ # max precision
    # params_vals = Dict(p => F(rand(1:prime_number)) for p in params)
    # inputs = Dict(u => [F(rand(1:prime_number)) for i in 1:prec] for u in inputs)
    # initial_conditions = Dict(x => F(rand(1:prime_number)) for x in states)
    # return params_vals, inputs, initial_conditions, prime_number, prec
    return prime_number, prec
end

function ReduceODEModP(state_eqs, out_eqs, states, inputs, outputs, params, prime_number)
    equations = vcat(state_eqs, out_eqs)
    GaloisF = Nemo.GF(prime_number)
    original_ring = parent(equations[1])
    new_ring, new_gens = Nemo.PolynomialRing(GaloisF, string.(gens(original_ring)))

    new_states = map(xx -> switch_ring(xx, new_ring), states)
    new_inputs = map(u -> switch_ring(u, new_ring), inputs)
    new_outputs = map(y -> switch_ring(y, new_ring), outputs)
    new_params = map(y -> switch_ring(y, new_ring), params)

    new_state_eqs = []
    new_out_eqs = []
    for i in 1:length(state_eqs)
        numer = get_numerator(equations[i])
        denom = get_denominator(equations[i])
        if isequal(denom, 0)
            throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
        end
        push!(new_state_eqs, change_base_ring(GaloisF, numer) // change_base_ring(GaloisF, denom)) # TODO: fix error on base ring change
    end

    for i in 1:length(out_eqs)
        numer = get_numerator(equations[i])
        denom = get_denominator(equations[i])
        if isequal(denom, 0)
            throw(Base.ArgumentError("Prime $p divides the denominator of $poly"))
        end
        push!(new_out_eqs, change_base_ring(GaloisF, numer) // change_base_ring(GaloisF, denom)) # TODO: (possibly same) fix error on base ring change
    end
    return new_state_eqs, new_out_eqs, new_states, new_inputs, new_outputs, new_params
end

# ModelingToolkit.@parameters θ[1:4]
# ModelingToolkit.@variables t, x[1:2], u[1:1], y[1:1]
# ModelingToolkit.@parameters ẋ[1:length(x)]
# D = Differential(t)
# using ModelingToolkit
# eqs = [D(x[1]) ~ x[1]^2 * θ[1] + θ[2] * x[1] * x[2] + u[1], D(x[2]) ~ θ[3] * x[1]^2 + θ[4] * x[1] * x[2]];
# outputs = [y[1] ~ x[1]];

# state_eqs, out_eqs, states, params, inputs, outputs = PreprocessODE(eqs, outputs, x, y, u, θ, t)

# # params_vals, input_values, initial_conditions, prime_number, ν = Initialize(state_eqs, out_eqs, states, inputs, outputs, params, 0.99)
# prime_number, ν = Initialize(state_eqs, out_eqs, states, inputs, outputs, params, 0.99)

# state_eqs, out_eqs, states, inputs, outputs, params = ReduceODEModP(state_eqs, out_eqs, states, inputs, outputs, params, prime_number)


function evaluate_poly(poly::MPolyElem, eval_dict)
    R = parent(first(values(eval_dict)))
    point = [get(eval_dict, v, zero(R)) for v in gens(parent(poly))]
    return evaluate(poly, point)
end

function evaluate_poly(poly::Generic.Frac{<: P}, eval_dict) where P <: MPolyElem
    numer, denom = get_numerator(poly), get_denominator(poly)
    return evaluate_poly(numer, eval_dict) // evaluate_poly(denom, eval_dict)
end

function PowerSeriesSolutionODE(
        state_eqs, out_eqs, states, params, inputs, outputs
    )

    # initialize get the prime and precision
    prime_number, ν = Initialize(state_eqs, out_eqs, states, inputs, outputs, params, proba)
    F = Nemo.GF(prime_number)

    # # get the input ring before reduction
    # poly_ring = parent(state_eqs[1])

    # reduce the ODE modulo prime and convert polynomials accordingly
    state_eqs, out_eqs, states, inputs, outputs, params = ReduceODEModP(state_eqs, out_eqs, states, inputs, outputs, params, prime_number) 
    poly_ring = parent(state_eqs[1])
    
    # random parameters, input p.s., and initial conditions
    params_vals = Dict(p => F(rand(1:prime_number)) for p in params)
    inputs_vals = Dict(u => [F(rand(1:prime_number)) for i in 1:ν] for u in inputs)
    initial_conditions = Dict(x => F(rand(1:prime_number)) for x in states)

    # create power series ring over rationals
    power_series_ring, τ = PowerSeriesRing(base_ring(poly_ring), ν, "τ"; model=:capped_absolute)

    # we need to switch the ring here since we will specify the parameters
    derivatives = [gen(base_ring(poly_ring), i) for i in 1:length(states)]
    new_ring, new_gens = Nemo.PolynomialRing(base_ring(poly_ring), string.(vcat(derivatives, states, inputs)))

    # 1. switch ring on each symbol
    # 2. for each state equation:
    #   2.1. evaluate numerator and denominator at parameters (substitute parameter values into the system)
    #   2.2. replace symbols (states, inputs) with the new-ring symbols in state and output equations
    
    # evaluation = Dict(k => v for (k, v) in params_vals) # key value pairs
    evaluation = Dict(k => new_ring(v) for (k, v) in params_vals)
    for v in vcat(derivatives, states, inputs)
        evaluation[v] = switch_ring(v, new_ring)
    end

    equations = [] # Array{Any}()
    for i in 1:length(state_eqs)
        num, den = map(p -> evaluate_poly(p, evaluation), [get_numerator(state_eqs[i]), get_denominator(state_eqs[i])])
        push!(equations, derivatives[i] * den - num)
    end






    


    








    MS_n_by_n = AbstractAlgebra.MatrixSpace(poly_ring, n, n)
    MS_n_by_1 = AbstractAlgebra.MatrixSpace(poly_ring, n, 1)
    Const_Space_n_by_1 = AbstractAlgebra.MatrixSpace(base_ring(poly_ring), n, 1)
    
    # compute Jacobians
    P = MS_n_by_1(eqs)
    ∂P∂ẋ = MS_n_by_n([derivative(p, deriv) for p in eqs, deriv in derivatives])
    ∂P∂x = MS_n_by_n([derivative(p, state) for p in eqs, state in states])
    ∂P∂θ = MS_n_by_n([derivative(p, param) for p in eqs, param in parameters])

    ν_current = 1

    # begin power series computation
    while ν_current < ν
        ν_new = min(ν, 2 * ν_current)

        for i in 1:length(states)
            set_precision!(solution[states[i]], ν_new)
            set_precision!(solution[derivatives[i]], ν_new)
        end
	
	# get a point to the right (power series) precision
        point = [solution[v] for v in gens(poly_ring)]
        set_precision!.(point, ν_new)

	# evaluate jacobians
        P_at_point = map(p -> evaluate(p, point), P) # P
	    ∂P∂x_at_point = map(p -> evaluate(p, point), ∂P∂x) # ∂P∂ẋ_at_point
        ∂P∂ẋ_at_point = map(p -> evaluate(p, point), ∂P∂ẋ) # ∂P∂x_at_point
        ∂P∂θ_at_point = map(p -> evaluate(p, point), ∂P∂θ) # ∂P∂x_at_point

	####
	#### Stuck Here
	####
	#### We're solving: ∂P∂ẋ_at_point * Ė + ∂P∂x_at_point * E + ∇P = 0
	#### where E is the error (see the original paper). Next, linear ode:
	#### Ė = -(∂P∂ẋ⁻¹ * ∂P∂x_at_point) * E - ∂P∂ẋ⁻¹ * ∇P
	
	#### get inverse of ∂P∂ẋ_at_point:
	    ∂P∂ẋ⁻¹ = InversePowerSeriesMatrix(∂P∂ẋ_at_point) # TO BE IMPLEMENTED

	#### Resolve Ė = -(∂P∂ẋ⁻¹ * ∂P∂x_at_point) * E - ∂P∂ẋ⁻¹ * ∇P
	#### A = -(∂P∂ẋ⁻¹ * ∂P∂x_at_point), B = ∂P∂ẋ⁻¹ * ∇P 
        A = - ∂P∂ẋ⁻¹ * ∂P∂x_at_point
        B = - ∂P∂ẋ⁻¹ * P_at_point
        InitialCondition = zero(Const_Space_n_by_1)

	#### This function will use method of variation of constants to resolve the 
	#### linear ode Ė = A * E + B
	    E = LinearSolution(A, B, InitialCondition) # TO BE IMPLEMENTED

	#### update solution dict via E (to be implemented)

	#### return solution dictionary
    end	
    return solution #### this will be a tuple
end

function InversePowerSeriesMatrix(J)
	#### returns inverse of J
end

function LinearSolution(A, B, IC)
	#### solve linear ODE Ė = A * E + B, with initial condition IC
	#### solve via variation of constant
    
	#### 1. Find Homogeneous Solution
    E_hom = HomogeneousResolution(A, B)
	
    #### 2. Find Particular Solution 
    E_part = ConstantVariation(E_hom, B)

    return E_part
end


