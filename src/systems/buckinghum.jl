using ModelingToolkit
using DynamicQuantities; 
using Unitful
using LinearAlgebra


"""
    no_of_fundamental_dims(mp::Vector{Any})
    Finds number of fundamental dimensions
"""
function no_of_fundamental_dims(mp)
    fundamental_dimensions = 0
    for val in mp
        if val == 1
            fundamental_dimensions += 1
        end
    end

    return fundamental_dimensions
end


"""
    get_dims_of_vars(Vector{Any}, Number, Vector{Any})

    Gets units of each variable in an array. Returns a matrix where each row corresponds
    to units represented in binary values
    Example : kg*m^3 - [3, 1, 0, 0, 0, 0]
    In the returned value each row corresponds to length, mass, time, current, luminosity 
    and temperature respectively. 
"""
function get_dims_of_vars(dims_vars, total_vars, dim_map)
    # For every single variable it contains row of 1s and 0s mentioning which unit is present
    dims_of_all_vars = zeros(total_vars, 6)

    for (ind, dim) in enumerate(dims_vars)
        temp_dims = 0
        temp_dims_arr = zeros(6)
        if ulength(dim) != 0
            dim_map[1] = 1
            temp_dims_arr[1] = ulength(dim)
            temp_dims += 1
        end 
        if umass(dim) != 0
            dim_map[2] = 1
            temp_dims_arr[2] = umass(dim)
            temp_dims += 1 
        end 
        if utime(dim) != 0
            dim_map[3] = 1
            temp_dims_arr[3] = utime(dim)
            temp_dims +=1
        end 
        if ucurrent(dim) != 0
            dim_map[4] = 1
            temp_dims_arr[4] = ucurrent(dim)
            temp_dims +=1
        end 
        if utemperature(dim) != 0
            dim_map[5] = 1
            temp_dims_arr[5] = utemperature(dim)

            temp_dims +=1 
        end 
        if uluminosity(dim) != 0
            dim_map[6] = 1
            temp_dims_arr[6] = uluminosity(dim)
            temp_dims +=1 
        end 
        dims_of_all_vars[ind,:] = temp_dims_arr'
    end

    return dims_of_all_vars
end


"""
    # find_pi_term_exponents(Matrix{Float64}, Number, Number)

    Finds PI terms and returns the exponents and indices of repeating variables and 
    non repeating index in the form of a dictionary 
"""
function find_pi_term_exponents(dims_of_all_vars, total_vars, fundamental_dimensions)
    pi_terms_data = []

    # We are considering the repeating variables as starting variables for now. 
    repeating_idx = []  
    for i in 1:fundamental_dimensions
        push!(repeating_idx, i)
    end
    non_repeating_idx = []    # V1
    for i in fundamental_dimensions+1:total_vars
        push!(non_repeating_idx, i)
    end

    for idx in non_repeating_idx
        # Form system of equations for exponents
        A = dims_of_all_vars[repeating_idx, 1:fundamental_dimensions]'  # Transpose for solving
        b = -dims_of_all_vars[idx, 1:fundamental_dimensions]
        
        exponents = A \ b  # Linear solve

        pi_term = Dict("var_idx" => idx, "exponents" => exponents, "repeating_idx" => repeating_idx)
        push!(pi_terms_data, pi_term)
    end
    
    return pi_terms_data
end


"""
    retrieve_pi_terms(Dict{Any}, Vector{Num})

    Gets the actual PI term provided non repeating index, repeating indices and the original variables
"""
function retrieve_pi_terms(pi_term,  var_names)
    repeating_idx = pi_term["repeating_idx"]
    exp_arr = pi_term["exponents"]
    final_pi_term = 1

    @. exp_arr = round(exp_arr, digits=2)

    for (ind, val) in enumerate(repeating_idx)
        final_pi_term *= var_names[val]^exp_arr[ind]
    end

    # multiplying with non repeating variable for the pi term
    final_pi_term *= var_names[pi_term["var_idx"]]

    return final_pi_term
end


"""
    buckinghumFun(Vector{DynamicQuantities.Quantity}, Vector{Num})

    Takes an array of DynamicQuantities.Quantity type and variable names separately. Gets the buckinghum 
    PI terms and returns the array of PI terms
"""
function buckinghumFun(vars_quants, var_names)

    # vars_quants : [m^3/kg, s^-1, ms^-1]
    # var_names  : [ρ, μ, u, v]

    # Number of variables
    total_vars = length(vars_quants)

    # Required for counting fundamental dimensions
    # says whicheever units are in picture. In this case : length, mass, time
    # [1, 1, 1, 0, 0, 0]
    dim_map = zeros(6)
    

    # [kg^-3* m, kg*s, s^-1]
    dims_vars = []
    for u_arr in vars_quants
        push!(dims_vars, u_arr.dimensions)
    end

    dims_of_all_vars = get_dims_of_vars(dims_vars, total_vars, dim_map)

    fundamental_dimensions = no_of_fundamental_dims(dim_map)

    pi_terms = find_pi_term_exponents(dims_of_all_vars, total_vars, fundamental_dimensions)

    pis_arr = []
    for val in pi_terms
        push!(pis_arr, retrieve_pi_terms(val, var_names))
    end

    return pis_arr
end


"""
    transform_sys(Vector{Any}, system::ODESystem)

    Substitutes the PI terms in the equations to form new set of equations and returns the
    new ODESystem
"""

function transform_sys(pi_eqs, sys::ODESystem, pis_vars)
    original_equations = equations(sys)
    orginal_eqs_mp = Dict{Any, Any}(eq.lhs => eq.rhs for eq in original_equations)
    transformed_eqs = []
    for eq in pi_eqs
        pi_term = eq
        for each_var in get_variables(pi_term)
            println(each_var)
            if haskey(orginal_eqs_mp, D(each_var)) != true
                orginal_eqs_mp[D(each_var)] = 0
            end
        end
    
        der_eq = expand_derivatives(D(pi_term))
        sub_eq = substitute(der_eq, orginal_eqs_mp)
        push!(transformed_eqs, sub_eq)
    end

    equations_for_system = Equation[]
    dependent_vars = Any[unknowns(sys)...]
    for (i, eq) in pairs(transformed_eqs)
        push!(dependent_vars, pis_vars[i])
        push!(equations_for_system, D(pis_vars[i]) ~ eq)
    end

    @named new_sys = ODESystem(equations_for_system, ModelingToolkit.get_iv(sys),dependent_vars, [α]; defaults=defaults(sys))

    return new_sys
end