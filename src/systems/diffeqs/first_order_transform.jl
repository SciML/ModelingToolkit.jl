extract_idv(eq::Equation) = eq.lhs.op.x

function lower_varname(O::Operation, naming_scheme; lower=false)
    @assert isa(O.op, Differential)

    D, x = O.op, O.args[1]
    order = lower ? D.order-1 : D.order

    sym = x.name
    name = order == 0 ? sym : Symbol(sym, naming_scheme, string(D.x.name)^order)

    Variable(name, x.subtype, x.dependents)
end
function lower_varname(var::Variable, idv, order::Int, naming_scheme)
    sym = var.name
    name = order == 0 ? sym : Symbol(sym, naming_scheme, string(idv.name)^order)
    return Variable(name, var.subtype, var.dependents)
end

function ode_order_lowering(sys::DiffEqSystem; kwargs...)
    eqs = sys.eqs
    ivs = sys.ivs
    eqs_lowered = ode_order_lowering(eqs; kwargs...)
    DiffEqSystem(eqs_lowered, ivs)
end
ode_order_lowering(eqs; naming_scheme = "_") = ode_order_lowering!(deepcopy(eqs), naming_scheme)
function ode_order_lowering!(eqs, naming_scheme)
    idv = extract_idv(eqs[1])
    D   = Differential(idv, 1)
    var_order = Dict{Variable,Int}()
    dv_name = eqs[1].lhs.args[1].subtype

    for eq in eqs
        var, maxorder = extract_var_order(eq)
        maxorder == 1 && continue # fast pass
        if maxorder > get(var_order, var, 0)
            var_order[var] = maxorder
        end
        lhs_renaming!(eq, D, naming_scheme)
        rhs_renaming!(eq, naming_scheme)
    end

    for var ∈ keys(var_order)
        order = var_order[var]
        for o in (order-1):-1:1
            lhs = D(lower_varname(var, idv, o-1, naming_scheme))
            rhs = lower_varname(var, idv, o, naming_scheme)
            eq = Equation(lhs, rhs)
            push!(eqs, eq)
        end
    end

    return eqs
end

function lhs_renaming!(eq, D, naming_scheme)
    eq.lhs = D(lower_varname(eq.lhs, naming_scheme, lower=true))
    return eq
end
rhs_renaming!(eq, naming_scheme) = _rec_renaming!(eq.rhs, naming_scheme)

function _rec_renaming!(rhs, naming_scheme)
    isa(rhs, Operation) && isa(rhs.op, Differential) && return lower_varname(rhs, naming_scheme)
    if rhs isa Operation
        args = rhs.args
        for i in eachindex(args)
            args[i] = _rec_renaming!(args[i], naming_scheme)
        end
    end
    rhs
end

function extract_var_order(eq)
    # We assume that the differential with the highest order is always going to be in the LHS
    dv = eq.lhs
    var = dv.args[1]
    order = dv.op.order
    return (var, order)
end

export ode_order_lowering
