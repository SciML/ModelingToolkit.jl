struct NonlinearSystem <: AbstractSystem
    eqs::Vector{Operation}
    vs::Vector{Variable}
    ps::Vector{Variable}
    v_name::Vector{Symbol}
    p_name::Symbol
end

function NonlinearSystem(eqs, vs, ps;
                         v_name = :Variable,
                         dv_name = :DependentVariable,
                         p_name = :Parameter)
    NonlinearSystem(eqs, vs, ps, [v_name,dv_name], p_name)
end

function NonlinearSystem(eqs;
                         v_name = :Variable,
                         dv_name = :DependentVariable,
                         p_name = :Parameter)
    # Allow the use of :DependentVariable to make it seamless with DE use
    dvs, vs, ps = extract_elements(eqs, (dv_name, v_name, p_name))
    vs = [dvs;vs]
    NonlinearSystem(eqs, vs, ps, [v_name,dv_name], p_name)
end

function generate_nlsys_function(sys::NonlinearSystem)
    var_exprs = [:($(sys.vs[i].name) = u[$i]) for i in 1:length(sys.vs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]
    sys_idxs = map(eq->isequal(eq.args[1],Constant(0)),sys.eqs)
    sys_eqs = sys.eqs[sys_idxs]
    calc_eqs = sys.eqs[.!(sys_idxs)]
    calc_exprs = [:($(Symbol("$(eq.args[1].name)")) = $(eq.args[2])) for eq in calc_eqs]
    sys_exprs = [:($(Symbol("resid[$i]")) = $(sys_eqs[i].args[2])) for i in eachindex(sys_eqs)]

    exprs = vcat(var_exprs,param_exprs,calc_exprs,sys_exprs)
    block = expr_arr_to_block(exprs)
    :((du,u,p)->$(block))
end

function generate_nlsys_jacobian(sys::NonlinearSystem,simplify=true)
    var_exprs = [:($(sys.vs[i].name) = u[$i]) for i in 1:length(sys.vs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]

    sys_idxs = map(eq->isequal(eq.args[1],Constant(0)),sys.eqs)
    sys_eqs = sys.eqs[sys_idxs]
    calc_eqs = sys.eqs[.!(sys_idxs)]
    sys_exprs = [:($(Symbol("resid[$i]")) = $(sys_eqs[i].args[2])) for i in eachindex(sys_eqs)]
    rhs = [eq.args[2] for eq in sys_eqs]

    for i in 1:length(calc_eqs)
        find_replace!.(rhs,calc_eqs[i].args[1],calc_eqs[i].args[2])
    end

    sys_exprs = calculate_jacobian(rhs,sys.vs)
    sys_exprs = Expression[expand_derivatives(expr) for expr in sys_exprs]
    if simplify
        sys_exprs = Expression[simplify_constants(expr) for expr in sys_exprs]
    end
    sys_exprs
end

export NonlinearSystem
export generate_nlsys_function
