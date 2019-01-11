function lower_varname(var::Variable, naming_scheme; lower=false)
    D = var.diff
    D == nothing && return var
    order = lower ? D.order-1 : D.order
    lower_varname(var.name, D.x, order, var.subtype, naming_scheme)
end
function lower_varname(sym::Symbol, idv, order::Int, subtype::Symbol, naming_scheme)
    order == 0 && return Variable(sym, subtype=subtype)
    name = Symbol(String(sym)*naming_scheme*String(idv.name)^order)
    Variable(name, subtype=subtype)
end

function ode_order_lowering(sys::DiffEqSystem; kwargs...)
    eqs = sys.eqs
    ivs = sys.ivs
    eqs_lowered = ode_order_lowering(eqs; kwargs...)
    DiffEqSystem(eqs_lowered, ivs)
end
ode_order_lowering(eqs; naming_scheme = "_") = ode_order_lowering!(deepcopy(eqs), naming_scheme)
function ode_order_lowering!(eqs, naming_scheme)
    ind = findfirst(x->!(isintermediate(x)), eqs)
    idv = extract_idv(eqs[ind])
    D   = Differential(idv, 1)
    sym_order = Dict{Symbol, Int}()
    dv_name = eqs[1].args[1].subtype
    for eq in eqs
        isintermediate(eq) && continue
        sym, maxorder = extract_symbol_order(eq)
        maxorder == 1 && continue # fast pass
        if maxorder > get(sym_order, sym, 0)
            sym_order[sym] = maxorder
        end
        eq = lhs_renaming!(eq, D, naming_scheme)
        eq = rhs_renaming!(eq, naming_scheme)
    end
    for sym in keys(sym_order)
        order = sym_order[sym]
        for o in (order-1):-1:1
            lhs = D(lower_varname(sym, idv, o-1, dv_name, naming_scheme))
            rhs = lower_varname(sym, idv, o, dv_name, naming_scheme)
            eq = Operation(==, [lhs, rhs])
            push!(eqs, eq)
        end
    end
    eqs
end

function lhs_renaming!(eq, D, naming_scheme)
    eq.args[1] = D(lower_varname(eq.args[1], naming_scheme, lower=true))
    return eq
end
function rhs_renaming!(eq, naming_scheme)
    rhs = eq.args[2]
    _rec_renaming!(rhs, naming_scheme)
end

function _rec_renaming!(rhs, naming_scheme)
    rhs isa Variable && rhs.diff != nothing && return lower_varname(rhs, naming_scheme)
    if rhs isa Operation
        args = rhs.args
        for i in eachindex(args)
            args[i] = _rec_renaming!(args[i], naming_scheme)
        end
    end
    rhs
end

function extract_symbol_order(eq)
    # We assume that the differential with the highest order is always going to be in the LHS
    dv = eq.args[1]
    sym = dv.name
    order = dv.diff.order
    sym, order
end

export ode_order_lowering
