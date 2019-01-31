extract_idv(eq::DiffEq) = eq.D.x

function lower_varname(D::Differential, x, naming_scheme; lower=false)
    order = lower ? D.order-1 : D.order
    return lower_varname(x, D.x, order, naming_scheme)
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
    vars = Variable[]
    dv_name = eqs[1].var.subtype

    for eq in eqs
        var, maxorder = extract_var_order(eq)
        maxorder == 1 && continue # fast pass
        if maxorder > get(var_order, var, 0)
            var_order[var] = maxorder
            var ∈ vars || push!(vars, var)
        end
        lhs_renaming!(eq, naming_scheme)
        rhs_renaming!(eq, naming_scheme)
    end

    for var ∈ vars
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

function lhs_renaming!(eq::DiffEq, naming_scheme)
    eq.var = lower_varname(eq.D, eq.var, naming_scheme, lower=true)
    return eq
end
rhs_renaming!(eq::DiffEq, naming_scheme) = _rec_renaming!(eq.rhs, naming_scheme)

function _rec_renaming!(rhs, naming_scheme)
    isa(rhs, Operation) && isa(rhs.op, Differential) &&
        return lower_varname(rhs.op, rhs.args[1], naming_scheme)
    if rhs isa Operation
        args = rhs.args
        for i in eachindex(args)
            args[i] = _rec_renaming!(args[i], naming_scheme)
        end
    end
    rhs
end

extract_var_order(eq::DiffEq) = (eq.var, eq.D.order)

export ode_order_lowering
