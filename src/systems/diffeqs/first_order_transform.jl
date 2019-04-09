export ode_order_lowering


function lower_varname(var::Variable, idv, order)
    order == 0 && return var
    name = Symbol(var.name, :_, string(idv.name)^order)
    return Variable(name; known = var.known)
end

function ode_order_lowering(sys::ODESystem)
    (eqs_lowered, new_vars) = ode_order_lowering(sys.eqs, sys.iv)
    ODESystem(eqs_lowered, sys.iv, [sys.dvs; new_vars], sys.ps)
end
function ode_order_lowering(eqs, iv)
    var_order = Dict{Variable,Int}()
    vars = Variable[]
    new_eqs = similar(eqs, DiffEq)
    new_vars = Variable[]

    for (i, eq) ∈ enumerate(eqs)
        var, maxorder = eq.x, eq.n
        if maxorder > get(var_order, var, 0)
            var_order[var] = maxorder
            any(isequal(var), vars) || push!(vars, var)
        end
        var′ = lower_varname(eq.x, iv, eq.n - 1)
        rhs′ = rename(eq.rhs)
        new_eqs[i] = DiffEq(var′, 1, rhs′)
    end

    for var ∈ vars
        order = var_order[var]
        for o in (order-1):-1:1
            lvar = lower_varname(var, iv, o-1)
            rvar = lower_varname(var, iv, o)
            push!(new_vars, rvar)

            rhs = rvar(iv())
            eq = DiffEq(lvar, 1, rhs)
            push!(new_eqs, eq)
        end
    end

    return (new_eqs, new_vars)
end

function rename(O::Expression)
    isa(O, Operation) || return O
    if is_derivative(O)
        (x, t, order) = flatten_differential(O)
        return lower_varname(x.op, t.op, order)(x.args...)
    end
    return Operation(O.op, rename.(O.args))
end
