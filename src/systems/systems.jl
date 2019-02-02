export generate_jacobian, generate_function


abstract type AbstractSystem end

function system_eqs    end
function system_vars   end
function system_params end

function generate_jacobian(sys::AbstractSystem)
    vs, ps = system_vars(sys), system_params(sys)
    var_exprs = [:($(vs[i].name) = u[$i]) for i in eachindex(vs)]
    param_exprs = [:($(ps[i].name) = p[$i]) for i in eachindex(ps)]
    jac = calculate_jacobian(sys)
    jac_exprs = [:(J[$i,$j] = $(convert(Expr, jac[i,j]))) for i in 1:size(jac,1), j in 1:size(jac,2)]
    exprs = vcat(var_exprs, param_exprs, vec(jac_exprs))
    block = expr_arr_to_block(exprs)
    :((J,u,p,t) -> $(block))
end

function generate_function(sys::AbstractSystem; version::FunctionVersion = ArrayFunction)
    sys_eqs = system_eqs(sys)
    vs, ps = system_vars(sys), system_params(sys)

    var_pairs   = [(u.name, :(u[$i])) for (i, u) ∈ enumerate(vs)]
    param_pairs = [(p.name, :(p[$i])) for (i, p) ∈ enumerate(ps)]
    (ls, rs) = collect(zip(var_pairs..., param_pairs...))

    var_eqs = Expr(:(=), build_expr(:tuple, ls), build_expr(:tuple, rs))
    sys_exprs = build_expr(:tuple, [convert(Expr, eq.rhs) for eq ∈ sys_eqs])
    let_expr = Expr(:let, var_eqs, sys_exprs)

    if version === ArrayFunction
        :((du,u,p,t) -> du .= $let_expr)
    elseif version === SArrayFunction
        :((u,p,t) -> begin
            du = $let_expr
            T = StaticArrays.similar_type(typeof(u), eltype(du))
            T(du)
        end)
    end
end
