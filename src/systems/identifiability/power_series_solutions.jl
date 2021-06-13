using AbstractAlgebra, LinearAlgebra, Symbolics, SymbolicUtils
using Symbolics:value
using Primes

function Initialize(eqs, states, inputs, outputs, parameters)
    # Proposition 3.3 in https://doi.org/10.1006/jsco.2002.0532
    solution = Dict()
    n = length(states)
    m = lengh(outputs)
    r = length(inputs)
    ℓ = length(parameters)
    D = 4 * (n + ell)^2 * (n + m) * d



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
        eqs, derivatives, initial_conditions, states, inputs, parameters, ν
    )

    solution = Initialize()
    n = length(eqs)
    poly_ring = parent(eqs[1]) # equations must be polynomials in a polynomial ring over rationals
    power_series_ring, τ = PowerSeriesRing(base_ring(poly_ring), ν, "τ"; model=:capped_absolute)
    
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


