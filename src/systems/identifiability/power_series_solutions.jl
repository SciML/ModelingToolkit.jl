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
        return one(parent(f))
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
    
    state_eqn_dict = Dict([x[i] => substitute(ODE[i].rhs, D.(x) .=> ẋ) for i in 1:length(ODE)])
    state_eqn_dict = Dict(substitute(value(k), input_symbols .=> gens_) => substitute(value(v), input_symbols .=> gens_) for (k, v) in state_eqn_dict)
    
    out_eqn_dict = Dict([y[i] => substitute(value(outputs[i].lhs), input_symbols .=> gens_) - substitute(value(outputs[i].rhs), input_symbols .=> gens_) for i in 1:length(outputs)])
    out_eqn_dict = Dict(substitute(value(k), input_symbols .=> gens_) => substitute(value(v), input_symbols .=> gens_) for (k, v) in out_eqn_dict)
    

    states = [substitute(value(each),  input_symbols .=> gens_) for each in x]
    params = [substitute(value(each),  input_symbols .=> gens_) for each in θ]
    inputs = [substitute(value(each),  input_symbols .=> gens_) for each in u]
    outputs = [substitute(value(each),  input_symbols .=> gens_) for each in y] 

    return (state_eqn_dict, out_eqn_dict, states, params, inputs, outputs)
end

function Initialize(state_eqs, out_eqs, states, inputs, outputs, params, proba)
    # Proposition 3.3 in https://doi.org/10.1006/jsco.2002.0532
    d, h = 1, 1
    for f in vcat([x for x in values(state_eqs)], [x for x in values(out_eqs)])
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
    GaloisF = Nemo.GF(prime_number)
    original_ring = parent(first(values(state_eqs)))
    new_ring, new_gens = Nemo.PolynomialRing(GaloisF, string.(gens(original_ring)))

    new_state_eqs = Dict(switch_ring(k, new_ring) => map_coefficients(x -> divexact(GaloisF(numerator(x)), GaloisF(denominator(x))), v) for (k, v) in state_eqs)
    new_out_eqs = Dict(switch_ring(k, new_ring) => map_coefficients(x -> divexact(GaloisF(numerator(x)), GaloisF(denominator(x))), v) for (k, v) in out_eqs)
    
    new_states = [switch_ring(x, new_ring) for x in states]
    new_outputs = [x for x in keys(new_out_eqs)]
    new_inputs = map(u -> switch_ring(u, new_ring), inputs)
    new_params = map(y -> switch_ring(y, new_ring), params)

    return new_state_eqs, new_out_eqs, new_states, new_inputs, new_outputs, new_params
end

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

    # reduce the ODE modulo prime and convert polynomials accordingly
    state_eqs, out_eqs, states, inputs, outputs, params = ReduceODEModP(state_eqs, out_eqs, states, inputs, outputs, params, prime_number) 
    poly_ring = parent(first(values(state_eqs)))
    
    # random parameters, input p.s., and initial conditions
    params_vals = Dict(p => F(rand(1:prime_number)) for p in params)
    inputs_vals = Dict(u => [F(rand(1:prime_number)) for i in 1:ν] for u in inputs)
    initial_conditions = Dict(x => F(rand(1:prime_number)) for x in states)

    # create power series ring over rationals
    power_series_ring, τ = PowerSeriesRing(base_ring(poly_ring), ν, "τ"; model=:capped_absolute)

    # we need to switch the ring here since we will specify the parameters
    derivatives_dict = Dict(states[i] => gen(poly_ring, i) for i in 1:length(states))
    derivatives = [x for x in values(derivatives_dict)]
    new_ring, new_gens = Nemo.PolynomialRing(base_ring(poly_ring), string.(vcat(derivatives, states, inputs)))
    derivatives_dict = Dict(k => switch_ring(v, new_ring) for (k, v) in derivatives_dict)
    
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
    for i in 1:length(states)
        eqn = state_eqs[states[i]]
        num, den = map(p -> evaluate_poly(p, evaluation), [get_numerator(eqn), get_denominator(eqn)])
        push!(equations, derivatives_dict[states[i]] * den - num)
    end

    new_inputs = Dict(switch_ring(k, new_ring) => v for (k, v) in input_vals)
    new_init = Dict(switch_ring(k, new_ring) => v for (k, v) in initial_conditions)



    # return #solution #### this will be a tuple
end


