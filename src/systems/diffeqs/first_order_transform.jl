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

function ode_order_lowering(sys::DiffEqSystem; kwargs...)
    eqs = sys.eqs
    eqs_lowered = ode_order_lowering(eqs; kwargs...)
    DiffEqSystem(eqs_lowered, sys.iv)
end
ode_order_lowering(eqs) = ode_order_lowering!(deepcopy(eqs))
function ode_order_lowering!(eqs)
    idv = extract_idv(eqs[1])
    D   = Differential(idv, 1)
    var_order = Dict{Variable,Int}()
    vars = Variable[]
    dv_name = eqs[1].var.subtype

    for eq in eqs
        var, maxorder = extract_var_order(eq)
        maxorder == 1 && continue # fast pass
        if maxorder > get(var_order, var, 0)
            var_order[var] = maxorder
            var ∈ vars || push!(vars, var)
        end
        lhs_renaming!(eq, D)
        rhs_renaming!(eq)
    end

    for var ∈ vars
        order = var_order[var]
        for o in (order-1):-1:1
            lvar = lower_varname(var, idv, o-1)
            rhs = lower_varname(var, idv, o)
            eq = DiffEq(D, lvar, rhs)
            push!(eqs, eq)
        end
    end

    return eqs
end

function lhs_renaming!(eq::DiffEq, D)
    eq.var = lower_varname(eq.D, eq.var, lower=true)
    eq.D = D
    return eq
end
rhs_renaming!(eq::DiffEq) = _rec_renaming!(eq.rhs)

function _rec_renaming!(rhs)
    isa(rhs, Operation) && isa(rhs.op, Differential) &&
        return lower_varname(rhs.op, rhs.args[1])
    if rhs isa Operation
        args = rhs.args
        for i in eachindex(args)
            args[i] = _rec_renaming!(args[i])
        end
    end
    rhs
end

extract_var_order(eq::DiffEq) = (eq.var, eq.D.order)

export ode_order_lowering
