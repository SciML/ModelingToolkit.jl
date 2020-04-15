function lower_varname(var::Variable, idv, order)
    order == 0 && return var
    name = Symbol(var.name, :ˍ, string(idv.name)^order)
    return Variable{vartype(var)}(name)
end

function flatten_differential(O::Operation)
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
    var_order = Dict{Variable,Int}()
    vars = Variable[]
    new_eqs = Equation[]
    new_vars = Variable[]
    D = Differential(iv())

    for (i, (eq, ss)) ∈ enumerate(zip(eqs, states))
        if isequal(eq.lhs, Constant(0))
            push!(new_vars, ss)
            push!(new_eqs, eq)
        else
            var, maxorder = var_from_nested_derivative(eq.lhs)
            if maxorder > get(var_order, var, 0)
                var_order[var] = maxorder
                any(isequal(var), vars) || push!(vars, var)
            end
            var′ = lower_varname(var, iv, maxorder - 1)
            rhs′ = rename_lower_order(eq.rhs)
            push!(new_vars, var′)
            push!(new_eqs, D(var′(iv())) ~ rhs′)
        end
    end

    for var ∈ vars
        order = var_order[var]
        for o in (order-1):-1:1
            lvar = lower_varname(var, iv, o-1)
            rvar = lower_varname(var, iv, o)
            push!(new_vars, lvar)

            rhs = rvar(iv())
            eq = Differential(iv())(lvar(iv())) ~ rhs
            push!(new_eqs, eq)
        end
    end

    return (new_eqs, new_vars)
end

function rename_lower_order(O::Expression)
    isa(O, Operation) || return O
    if is_derivative(O)
        (x, t, order) = flatten_differential(O)
        return lower_varname(x.op, t.op, order)(x.args...)
    end
    return Operation(O.op, rename_lower_order.(O.args))
end
