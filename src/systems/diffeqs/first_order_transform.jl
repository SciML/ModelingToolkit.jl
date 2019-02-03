function lower_varname(D::Differential, x; lower=false)
    order = lower ? D.order-1 : D.order
    return lower_varname(x, D.x, order)
end
function lower_varname(var::Variable, idv, order::Int)
    sym = var.name
    name = order == 0 ? sym : Symbol(sym, :_, string(idv.name)^order)
    return Variable(name, var.subtype, var.dependents)
end

function ode_order_lowering(sys::DiffEqSystem)
    eqs_lowered = ode_order_lowering(sys.eqs, sys.iv)
    DiffEqSystem(eqs_lowered, sys.iv)
end
function ode_order_lowering(eqs, iv)
    D = Differential(iv, 1)
    var_order = Dict{Variable,Int}()
    vars = Variable[]
    new_eqs = similar(eqs, DiffEq)

    for (i, eq) ∈ enumerate(eqs)
        var, maxorder = eq.var, eq.D.order
        maxorder == 1 && continue # fast pass
        if maxorder > get(var_order, var, 0)
            var_order[var] = maxorder
            var ∈ vars || push!(vars, var)
        end
        var′ = lower_varname(eq.D, eq.var, lower = true)
        rhs′ = rename(eq.rhs)
        new_eqs[i] = DiffEq(D, var′, rhs′)
    end

    for var ∈ vars
        order = var_order[var]
        for o in (order-1):-1:1
            lvar = lower_varname(var, iv, o-1)
            rhs = lower_varname(var, iv, o)
            eq = DiffEq(D, lvar, rhs)
            push!(new_eqs, eq)
        end
    end

    return new_eqs
end

function rename(O::Expression)
    isa(O, Operation) || return O
    isa(O.op, Differential) && return lower_varname(O.op, O.args[1])
    return Operation(O.op, rename.(O.args))
end

export ode_order_lowering
