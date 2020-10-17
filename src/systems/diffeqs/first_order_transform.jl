function lower_varname(var::Term, idv, order)
    order == 0 && return var
    name = Symbol(nameof(var.op), :ˍ, string(idv)^order)
    #name = Symbol(var.name, :ˍ, string(idv.name)^order)
    return Sym{symtype(var.op)}(name)(var.args[1])
end

function lower_varname(t::Term, iv)
    var, order = var_from_nested_derivative(t)
    lower_varname(var, iv, order)
end
lower_varname(t::Sym, iv) = t

function flatten_differential(O::Term)
    @assert is_derivative(O) "invalid differential: $O"
    is_derivative(O.args[1]) || return (O.args[1], O.op.x, 1)
    (x, t, order) = flatten_differential(O.args[1])
    isequal(t, O.op.x) || throw(ArgumentError("non-matching differentials on lhs: $t, $(O.op.x)"))
    return (x, t, order + 1)
end

"""
$(TYPEDSIGNATURES)

Takes a Nth order ODESystem and returns a new ODESystem written in first order
form by defining new variables which represent the N-1 derivatives.
"""
function ode_order_lowering(sys::ODESystem)
    eqs_lowered, new_vars = ode_order_lowering(equations(sys), sys.iv, states(sys))
    return ODESystem(eqs_lowered, sys.iv, new_vars, sys.ps)
end

function ode_order_lowering(eqs, iv, states)
    var_order = OrderedDict{Any,Int}()
    D = Differential(iv)
    diff_eqs = Equation[]
    diff_vars = []
    alge_eqs = Equation[]
    alge_vars = []

    for (i, (eq, ss)) ∈ enumerate(zip(eqs, states))
        if _iszero(eq.lhs)
            push!(alge_vars, ss)
            push!(alge_eqs, eq)
        else
            var, maxorder = var_from_nested_derivative(eq.lhs)
            # only save to the dict when we need to lower the order to save memory
            maxorder > get(var_order, var, 1) && (var_order[var] = maxorder)
            var′ = lower_varname(var, iv, maxorder - 1)
            rhs′ = rename_lower_order(eq.rhs)
            push!(diff_vars, var′)
            push!(diff_eqs, D(var′) ~ rhs′)
        end
    end

    for (var, order) ∈ var_order
        for o in (order-1):-1:1
            lvar = lower_varname(var, iv, o-1)
            rvar = lower_varname(var, iv, o)
            push!(diff_vars, lvar)

            rhs = rvar
            eq = Differential(iv)(lvar) ~ rhs
            push!(diff_eqs, eq)
        end
    end

    # we want to order the equations and variables to be `(diff, alge)`
    return (vcat(diff_eqs, alge_eqs), vcat(diff_vars, alge_vars))
end

function rename_lower_order(O)
    isa(O, Term) || return O
    if is_derivative(O)
        (x, t, order) = flatten_differential(O)
        return lower_varname(x, t, order)
    end
    return Term(O.op, rename_lower_order.(O.args))
end
