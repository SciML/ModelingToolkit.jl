export ode_order_lowering


function lower_varname(var::Variable, idv, order)
    order == 0 && return var
    name = Symbol(var.name, :_, string(idv.name)^order)
    return Variable(name, var.dependents; known = var.known)
end

function ode_order_lowering(sys::DiffEqSystem)
    eqs_lowered, vars_lowered = ode_order_lowering(sys.eqs, sys.iv)
    DiffEqSystem(eqs_lowered, sys.iv, vars_lowered, sys.ps)
end
function ode_order_lowering(eqs, iv)
    var_order = Dict{Variable,Int}()
    vars = Variable[]
    new_vars = Variable[]
    new_eqs = similar(eqs, DiffEq)

    for (i, eq) ∈ enumerate(eqs)
        var, maxorder = eq.x, eq.n
        if maxorder > get(var_order, var, 0)
            var_order[var] = maxorder
            any(isequal(var), vars) || push!(vars, var)
        end
        var′ = lower_varname(eq.x, eq.t, eq.n - 1)
        push!(new_vars, var′)
        rhs′ = rename(eq.rhs)
        new_eqs[i] = DiffEq(var′, iv, 1, rhs′)
    end

    for var ∈ vars
        order = var_order[var]
        for o in (order-1):-1:1
            lvar = lower_varname(var, iv, o-1)
            push!(new_vars, lvar)
            rhs = lower_varname(var, iv, o)
            eq = DiffEq(lvar, iv, 1, rhs)
            push!(new_eqs, eq)
        end
    end

    return new_eqs, new_vars
end

function rename(O::Expression)
    isa(O, Operation) || return O
    if is_derivative(O)
        (x, t, order) = flatten_differential(O)
        return lower_varname(x, t, order)
    end
    return Operation(O.op, rename.(O.args))
end
