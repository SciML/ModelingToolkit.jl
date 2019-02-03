export generate_jacobian, generate_function


abstract type AbstractSystem end

function system_eqs    end
function system_vars   end
function system_params end

function generate_jacobian(sys::AbstractSystem; version = ArrayFunction)
    vs, ps = system_vars(sys), system_params(sys)
    jac = calculate_jacobian(sys)
    return build_function(jac, vs, ps, (:t,); version = version)
end

function generate_function(sys::AbstractSystem; version::FunctionVersion = ArrayFunction)
    sys_eqs = system_eqs(sys)
    vs, ps = system_vars(sys), system_params(sys)
    return build_function([eq.rhs for eq ∈ sys_eqs], vs, ps, (:t,); version = version)
end

function build_function(rhss, vs, ps, args; version::FunctionVersion)
    var_pairs   = [(u.name, :(u[$i])) for (i, u) ∈ enumerate(vs)]
    param_pairs = [(p.name, :(p[$i])) for (i, p) ∈ enumerate(ps)]
    (ls, rs) = collect(zip(var_pairs..., param_pairs...))

    var_eqs = Expr(:(=), build_expr(:tuple, ls), build_expr(:tuple, rs))

    if version === ArrayFunction
        X = gensym()
        sys_exprs = [:($X[$i] = $(convert(Expr, rhs))) for (i, rhs) ∈ enumerate(rhss)]
        let_expr = Expr(:let, var_eqs, build_expr(:block, sys_exprs))
        :(($X,u,p,$(args...)) -> $let_expr)
    elseif version === SArrayFunction
        sys_expr = build_expr(:tuple, [convert(Expr, rhs) for rhs ∈ rhss])
        let_expr = Expr(:let, var_eqs, sys_expr)
        :((u,p,$(args...)) -> begin
            X = $let_expr
            T = StaticArrays.similar_type(typeof(u), eltype(X))
            T(X)
        end)
    end
end
