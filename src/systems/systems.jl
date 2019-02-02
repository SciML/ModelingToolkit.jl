export generate_jacobian


abstract type AbstractSystem end


function system_vars end
function system_params end

function generate_jacobian(sys::AbstractSystem, simplify = true)
    vs, ps = system_vars(sys), system_params(sys)
    var_exprs = [:($(vs[i].name) = u[$i]) for i in eachindex(vs)]
    param_exprs = [:($(ps[i].name) = p[$i]) for i in eachindex(ps)]
    jac = calculate_jacobian(sys, simplify)
    jac_exprs = [:(J[$i,$j] = $(convert(Expr, jac[i,j]))) for i in 1:size(jac,1), j in 1:size(jac,2)]
    exprs = vcat(var_exprs, param_exprs, vec(jac_exprs))
    block = expr_arr_to_block(exprs)
    :((J,u,p,t) -> $(block))
end
