extract_idv(eq::DiffEq) = eq.D.x

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
        var, maxorder = extract_var_order(eq)
        maxorder == 1 && continue # fast pass
        if maxorder > get(var_order, var, 0)
            var_order[var] = maxorder
            var ∈ vars || push!(vars, var)
        end
        eq = lhs_renaming(eq, D)
        eq = rhs_renaming(eq)
        new_eqs[i] = eq
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

lhs_renaming(eq::DiffEq, D) = DiffEq(D, lower_varname(eq.D, eq.var, lower=true), eq.rhs)
rhs_renaming(eq::DiffEq) = DiffEq(eq.D, eq.var, _rec_renaming(eq.rhs))

function _rec_renaming(rhs)
    isa(rhs, Operation) || return rhs
    isa(rhs.op, Differential) && return lower_varname(rhs.op, rhs.args[1])
    return Operation(rhs.op, _rec_renaming.(rhs.args))
end

extract_var_order(eq::DiffEq) = (eq.var, eq.D.order)

export ode_order_lowering
