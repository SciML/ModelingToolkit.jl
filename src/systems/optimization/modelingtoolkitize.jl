"""
$(TYPEDSIGNATURES)

Generate `OptimizationSystem`, dependent variables, and parameters from an `OptimizationProblem`.
"""
function modelingtoolkitize(prob::DiffEqBase.OptimizationProblem;
        u_names = nothing, p_names = nothing, kwargs...)
    num_cons = isnothing(prob.lcons) ? 0 : length(prob.lcons)
    if prob.p isa Tuple || prob.p isa NamedTuple
        p = [x for x in prob.p]
    else
        p = prob.p
    end
    has_p = !(p isa Union{DiffEqBase.NullParameters, Nothing})
    if u_names !== nothing
        varnames_length_check(prob.u0, u_names; is_unknowns = true)
        _vars = [variable(name) for name in u_names]
    elseif SciMLBase.has_sys(prob.f)
        varnames = getname.(variable_symbols(prob.f.sys))
        varidxs = variable_index.((prob.f.sys,), varnames)
        invpermute!(varnames, varidxs)
        _vars = [variable(name) for name in varnames]
        if prob.f.sys isa OptimizationSystem
            for (i, sym) in enumerate(variable_symbols(prob.f.sys))
                if hasbounds(sym)
                    _vars[i] = Symbolics.setmetadata(
                        _vars[i], VariableBounds, getbounds(sym))
                end
            end
        end
    else
        _vars = [variable(:x, i) for i in eachindex(prob.u0)]
    end
    _vars = reshape(_vars, size(prob.u0))
    vars = ArrayInterface.restructure(prob.u0, _vars)
    if prob.ub !== nothing # lb is also !== nothing
        vars = map(vars, prob.lb, prob.ub) do sym, lb, ub
            if iszero(lb) && iszero(ub) || isinf(lb) && lb < 0 && isinf(ub) && ub > 0
                sym
            else
                Symbolics.setmetadata(sym, VariableBounds, (lb, ub))
            end
        end
    end
    params = if has_p
        if p_names === nothing && SciMLBase.has_sys(prob.f)
            p_names = Dict(parameter_index(prob.f.sys, sym) => sym
            for sym in parameter_symbols(prob.f.sys))
        end
        if p isa MTKParameters
            old_to_new = Dict()
            for sym in parameter_symbols(prob)
                idx = parameter_index(prob, sym)
                old_to_new[unwrap(sym)] = unwrap(p_names[idx])
            end
            order = reorder_parameters(prob.f.sys, parameters(prob.f.sys))
            for arr in order
                for i in eachindex(arr)
                    arr[i] = old_to_new[arr[i]]
                end
            end
            _params = order
        else
            _params = define_params(p, p_names)
        end
        p isa Number ? _params[1] :
        (p isa Tuple || p isa NamedTuple || p isa AbstractDict || p isa MTKParameters ?
         _params :
         ArrayInterface.restructure(p, _params))
    else
        []
    end

    if p isa MTKParameters
        eqs = prob.f(vars, params...)
    else
        eqs = prob.f(vars, params)
    end

    if DiffEqBase.isinplace(prob) && !isnothing(prob.f.cons)
        lhs = Array{Num}(undef, num_cons)
        if p isa MTKParameters
            prob.f.cons(lhs, vars, params...)
        else
            prob.f.cons(lhs, vars, params)
        end
        cons = Union{Equation, Inequality}[]

        if !isnothing(prob.lcons)
            for i in 1:num_cons
                if !isinf(prob.lcons[i])
                    if prob.lcons[i] != prob.ucons[i]
                        push!(cons, prob.lcons[i] ≲ lhs[i])
                    else
                        push!(cons, lhs[i] ~ prob.ucons[i])
                    end
                end
            end
        end

        if !isnothing(prob.ucons)
            for i in 1:num_cons
                if !isinf(prob.ucons[i]) && prob.lcons[i] != prob.ucons[i]
                    push!(cons, lhs[i] ≲ prob.ucons[i])
                end
            end
        end
        if (isnothing(prob.lcons) || all(isinf, prob.lcons)) &&
           (isnothing(prob.ucons) || all(isinf, prob.ucons))
            throw(ArgumentError("Constraints passed have no proper bounds defined.
            Ensure you pass equal bounds (the scalar that the constraint should evaluate to) for equality constraints
            or pass the lower and upper bounds for inequality constraints."))
        end
    elseif !isnothing(prob.f.cons)
        cons = p isa MTKParameters ? prob.f.cons(vars, params...) :
               prob.f.cons(vars, params)
    else
        cons = []
    end
    params = values(params)
    params = if params isa Number || (params isa Array && ndims(params) == 0)
        [params[1]]
    elseif p isa MTKParameters
        reduce(vcat, params)
    else
        vec(collect(params))
    end

    sts = vec(collect(vars))
    default_u0 = Dict(sts .=> vec(collect(prob.u0)))
    default_p = if has_p
        if prob.p isa AbstractDict
            Dict(v => prob.p[k] for (k, v) in pairs(_params))
        elseif prob.p isa MTKParameters
            Dict(params .=> reduce(vcat, prob.p))
        else
            Dict(params .=> vec(collect(prob.p)))
        end
    else
        Dict()
    end
    de = OptimizationSystem(eqs, sts, params;
        name = gensym(:MTKizedOpt),
        constraints = cons,
        defaults = merge(default_u0, default_p),
        kwargs...)
    de
end
