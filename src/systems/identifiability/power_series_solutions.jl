using LinearAlgebra, Symbolics, SymbolicUtils
using Symbolics:value, polygamma
using Primes
using Nemo

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
    elseif f isa fmpq_mpoly
        return f
    end
end

function get_denominator(f)
    if f isa Generic.Frac
        return denominator(f)
    elseif f isa fmpq_mpoly
        return one(parent(f))
    end
end

function get_height_and_coeff(f)
    if length(f)=0
        return (0, 1)
    end
    max_ = 1
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
    ModelingToolkit.@parameters ẋ[1:length(x)]

    input_symbols = vcat(ẋ,x,u,y,θ,t)
    generators = string.(input_symbols)
    R, gens_ = Nemo.PolynomialRing(Nemo.QQ, generators)
    
    state_eqs = [substitute(eqn.lhs, D.(x).=>ẋ) ~ substitute(eqn.rhs, D.(x).=>ẋ) for eqn in ODE];
    out_eqs = [substitute(value(eqn.lhs - eqn.rhs), input_symbols.=>gens_) for eqn in outputs]

    state_eqn_dict = Dict(substitute(value(eqn.lhs), input_symbols .=> gens_) => substitute(value(eqn.rhs), input_symbols .=> gens_) for eqn in state_eqs)
    state_eqs = [k - v for (k, v) in state_eqn_dict]

    states = [substitute(value(each),  input_symbols.=>gens_) for each in x]
    params = [substitute(value(each),  input_symbols.=>gens_) for each in θ]
    inputs = [substitute(value(each),  input_symbols.=>gens_) for each in u]
    outputs = [substitute(value(each),  input_symbols.=>gens_) for each in y] 

    return (state_eqs, out_eqs, states, params, inputs, outputs, state_eqn_dict)
end

function Initialize(x_eqs, y_eqs, states, inputs, outputs, params, proba)
    # Proposition 3.3 in https://doi.org/10.1006/jsco.2002.0532
    d, h = 1, 1
    for f in vcat(x_eqs, y_eqs)
        numer_df, numer_hf = numerator_height_coeff(f)
        denom_df, denom_hf = denominator_height_coeff(f)
        d = max(d, max(numer_df, denom_df))
        h = max(h, max(numer_hf, denom_hf))
    end
    p_per_func = 1 - (1 - proba) / (length(states)+length(params)) # TODO: add capability for checking different functions
    mu = ceil(1 / (1 - sqrt(p_per_func)))
    solution = Dict()
    n = length(states)
    m = length(outputs)
    r = length(inputs)
    ℓ = length(params)
    D = 4 * (n + ℓ)^2 * (n + m) * d
    Dprime = D * (2 * log(n + ℓ + r + 1) + log(mu * D)) + 4 * (n + ℓ)^2 * ((n + m) * h + log(2 * n * D))
    prime = Primes.nextprime(Int(ceil(2 * mu * Dprime)))
    F = Nemo.GF(prime)

    # TODO: reduce the x_eqs, y_eqs modulo prime: change polynomials to F-based

    prec = n + ℓ # max precision

    params_vals = Dict(p => F(rand(1:prime)) for p in params)
    inputs = Dict(u => [F(rand(1:prime)) for i in 1:prec] for u in inputs)
    initial_conditions = Dict(x => F(rand(1:prime)) for x in states)
    return params_vals, inputs, initial_conditions
end


#### this file contains code that computes power series solution to ODE
#### The local identifiability algorithm will compute power series solution with θ (parameters)
#### specialized to random values
#### u is specialized to a random power series of degree ν

#### The original ode: ẋ=F(x, θ, u), F is A rational function
#### Let P = numerator(ẋ-F(x, θ, u)).
#### Compute two Jacobians: ∂P∂ẋ, ∂P∂x
#### Specialize them to random (integer) values of θ and u (power series with int. coefficients): 
####	P(θ̂, û), ∂P∂ẋ_at_point, ∂P∂x_at_point (see code below)

function PowerSeriesSolution(
        x_eqs, y_eqs, states, inputs, outputs, parameters, proba, prec
    )

    params_vals, inputs, initial_conditions = Initialize(x_eqs, y_eqs, states, inputs, outputs, params, proba)
    n = length(eqs)

    poly_ring = parent(x_eqs[1]) # equations must be polynomials in a polynomial ring over rationals
    power_series_ring, τ = PowerSeriesRing(base_ring(poly_ring), prec, "τ"; model=:capped_absolute)
    
    # we need to switch the ring here since we will specify the parameters
    derivatives = [gen(poly_ring, i) for i in 1:length(states)]
    new_ring, new_gens = Nemo.PolynomialRing(base_ring(poly_ring), string.(vcat(derivatives, states, inputs)))

    # 1. switch ring on each symbol
    # 2. for each state equation:
    #   2.1. evaluate numerator and denominator at parameters (substitute parameter values into the system)
    #   2.2. replace symbols (states, inputs) with the new-ring symbols in state and output equations
    
    evaluation = Dict(k => new_ring(v) for (k, v) in param_values) # key value pairs
    


    








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


