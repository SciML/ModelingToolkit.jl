struct NonlinearSystem <: AbstractSystem
    eqs::Vector{Equation}
    vs::Vector{Variable}
    ps::Vector{Variable}
end

function NonlinearSystem(eqs)
    vs, ps = extract_elements(eqs, [_subtype(:Unknown), _subtype(:Parameter)])
    NonlinearSystem(eqs, vs, ps)
end

function generate_nlsys_function(sys::NonlinearSystem)
    var_exprs = [:($(sys.vs[i].name) = u[$i]) for i in 1:length(sys.vs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]
    sys_eqs, calc_eqs = partition(eq -> isequal(eq.lhs, Constant(0)), sys.eqs)
    calc_exprs = [:($(eq.lhs.name) = $(eq.rhs)) for eq in calc_eqs if isa(eq.lhs, Variable)]
    sys_exprs = [:($(Symbol("resid[$i]")) = $(sys_eqs[i].rhs)) for i in eachindex(sys_eqs)]

    exprs = vcat(var_exprs,param_exprs,calc_exprs,sys_exprs)
    block = expr_arr_to_block(exprs)
    :((du,u,p)->$(block))
end

function calculate_jacobian(sys::NonlinearSystem,simplify=true)
    sys_eqs, calc_eqs = partition(eq -> isequal(eq.lhs, Constant(0)), sys.eqs)
    rhs = [eq.rhs for eq in sys_eqs]

    for calc_eq ∈ calc_eqs
        find_replace!.(rhs, calc_eq.lhs, calc_eq.rhs)
    end

    sys_exprs = calculate_jacobian(rhs,sys.vs)
    sys_exprs = Expression[expand_derivatives(expr) for expr in sys_exprs]
    sys_exprs
end

function generate_nlsys_jacobian(sys::NonlinearSystem,simplify=true)
    var_exprs = [:($(sys.vs[i].name) = u[$i]) for i in 1:length(sys.vs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]
    jac = calculate_jacobian(sys,simplify)
    jac_exprs = [:(J[$i,$j] = $(convert(Expr, jac[i,j]))) for i in 1:size(jac,1), j in 1:size(jac,2)]
    exprs = vcat(var_exprs,param_exprs,vec(jac_exprs))
    block = expr_arr_to_block(exprs)
    :((J,u,p,t)->$(block))
end

export NonlinearSystem
export generate_nlsys_function
