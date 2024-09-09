"""
$(TYPEDSIGNATURES)

Takes a Nth order ODESystem and returns a new ODESystem written in first order
form by defining new variables which represent the N-1 derivatives.
"""
function ode_order_lowering(sys::ODESystem)
    iv = get_iv(sys)
    eqs_lowered, new_vars = ode_order_lowering(equations(sys), iv, unknowns(sys))
    @set! sys.eqs = eqs_lowered
    @set! sys.unknowns = new_vars
    return sys
end

function dae_order_lowering(sys::ODESystem)
    iv = get_iv(sys)
    eqs_lowered, new_vars = dae_order_lowering(equations(sys), iv, unknowns(sys))
    @set! sys.eqs = eqs_lowered
    @set! sys.unknowns = new_vars
    return sys
end

function ode_order_lowering(eqs, iv, unknown_vars)
    var_order = OrderedDict{Any, Int}()
    D = Differential(iv)
    diff_eqs = Equation[]
    diff_vars = []
    alge_eqs = Equation[]

    for (i, eq) in enumerate(eqs)
        if !isdiffeq(eq)
            push!(alge_eqs, eq)
        else
            var, maxorder = var_from_nested_derivative(eq.lhs)
            maxorder > get(var_order, var, 1) && (var_order[var] = maxorder)
            var′ = lower_varname(var, iv, maxorder - 1)
            rhs′ = diff2term_with_unit(eq.rhs, iv)
            push!(diff_vars, var′)
            push!(diff_eqs, D(var′) ~ rhs′)
        end
    end

    for (var, order) in var_order
        for o in (order - 1):-1:1
            lvar = lower_varname(var, iv, o - 1)
            rvar = lower_varname(var, iv, o)
            push!(diff_vars, lvar)

            rhs = rvar
            eq = Differential(iv)(lvar) ~ rhs
            push!(diff_eqs, eq)
        end
    end

    # we want to order the equations and variables to be `(diff, alge)`
    return (vcat(diff_eqs, alge_eqs), vcat(diff_vars, setdiff(unknown_vars, diff_vars)))
end

function dae_order_lowering(eqs, iv, unknown_vars)
    var_order = OrderedDict{Any, Int}()
    D = Differential(iv)
    diff_eqs = Equation[]
    diff_vars = OrderedSet()
    alge_eqs = Equation[]
    vars = Set()
    subs = Dict()

    for (i, eq) in enumerate(eqs)
        vars!(vars, eq)
        n_diffvars = 0
        for vv in vars
            isdifferential(vv) || continue
            var, maxorder = var_from_nested_derivative(vv)
            isparameter(var) && continue
            n_diffvars += 1
            order = get(var_order, var, nothing)
            seen = order !== nothing
            if !seen
                order = 1
            end
            maxorder > order && (var_order[var] = maxorder)
            var′ = lower_varname(var, iv, maxorder - 1)
            subs[vv] = D(var′)
            if !seen
                push!(diff_vars, var′)
            end
        end
        n_diffvars == 0 && push!(alge_eqs, eq)
        empty!(vars)
    end

    for (var, order) in var_order
        for o in (order - 1):-1:1
            lvar = lower_varname(var, iv, o - 1)
            rvar = lower_varname(var, iv, o)
            push!(diff_vars, lvar)

            rhs = rvar
            eq = Differential(iv)(lvar) ~ rhs
            push!(diff_eqs, eq)
        end
    end

    return ([diff_eqs; substitute.(eqs, (subs,))],
        vcat(collect(diff_vars), setdiff(unknown_vars, diff_vars)))
end
