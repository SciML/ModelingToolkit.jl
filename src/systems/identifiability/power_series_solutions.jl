using AbstractAlgebra, LinearAlgebra, Symbolics, SymbolicUtils
using Symbolics:value

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
        eqs, states, derivatives, initial_conditions, inputs, prec
    )
    n = length(eqs)
    poly_ring = parent(eqs[1])
    power_series_ring, τ = PowerSeriesRing(base_ring(poly_ring), prec, "τ"; model=:capped_absolute)
    
    MS_n_by_n = AbstractAlgebra.MatrixSpace(poly_ring, n, n)
    MS_n_by_1 = AbstractAlgebra.MatrixSpace(poly_ring, n, 1)
    Const_Space_n_by_1 = AbstractAlgebra.MatrixSpace(base_ring(poly_ring), n, 1)
    
    # compute Jacobians
    P = MS_n_by_1(eqs)
    ∂P∂ẋ = MS_n_by_n([derivative(p, deriv) for p in eqs, deriv in derivatives])
    ∂P∂x = MS_n_by_n([derivative(p, state) for p in eqs, state in states])


    solution = Dict()
    for (u, coeffs) in inputs
        solution[u] = sum([coeffs[i] * τ^(i - 1) for i in 1:length(coeffs)])
    end

    for state in states
        solution[state] = power_series_ring(initial_conditions[state])
    end

    for dxdt in derivatives
        solution[dxdt] = power_series_ring(0)
        set_precision!(solution[dxdt], 1)
    end

    cur_prec = 1

    # begin power series computation
    while cur_prec < prec
        new_prec = min(prec, 2 * cur_prec)

        for i in 1:length(states)
            set_precision!(solution[states[i]], new_prec)
            set_precision!(solution[derivatives[i]], new_prec)
        end
	
	# get a point to the right precision
        eval_point = [solution[v] for v in gens(poly_ring)]
        set_precision!.(eval_point, 2 * cur_prec)
	
	# evaluate jacobians
        eqs_eval = map(p -> evaluate(p, eval_point), P) # ∇P
	∂P∂x_at_point = map(p -> evaluate(p, eval_point), ∂P∂x) # ∂P∂ẋ_at_point
        ∂P∂ẋ_at_point = map(p -> evaluate(p, eval_point), ∂P∂ẋ) # ∂P∂x_at_point

	####
	#### Stuck Here
	####
	#### We're solving: ∂P∂ẋ_at_point * Ė + ∂P∂x_at_point * E + ∇P = 0
	#### where E is the error (see the original paper). Next, linear ode:
	#### Ė = -(∂P∂ẋ⁻¹ * ∂P∂x_at_point) * E - ∂P∂ẋ⁻¹ * ∇P
	
	#### get inverse of ∂P∂ẋ_at_point:
	∂P∂ẋ⁻¹ = PowerSeriesInverseMatrix(∂P∂ẋ_at_point) # TO BE IMPLEMENTED

	#### Resolve Ė = -(∂P∂ẋ⁻¹ * ∂P∂x_at_point) * E - ∂P∂ẋ⁻¹ * ∇P
	#### A = -(∂P∂ẋ⁻¹ * ∂P∂x_at_point), B = ∂P∂ẋ⁻¹ * ∇P 
	A = - ∂P∂ẋ⁻¹ * ∂P∂x_at_point
	B = - ∂P∂ẋ⁻¹ * eqs_eval
	InitialCondition = zero(Const_Space_n_by_1)

	#### This function will use method of variation of constants to resolve the 
	#### linear ode Ė = A * E + B
	E = LinearSolution(A, B, InitialCondition) # TO BE IMPLEMENTED

	#### update solution dict via E (to be implemented)

	#### 	
    end	
    return solution #### this will be a tuple
end
