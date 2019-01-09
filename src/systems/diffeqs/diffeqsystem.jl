mutable struct DiffEqSystem <: AbstractSystem
    eqs::Vector{Operation}
    ivs::Vector{Variable}
    dvs::Vector{Variable}
    vs::Vector{Variable}
    ps::Vector{Variable}
    iv_name::Symbol
    dv_name::Symbol
    p_name::Symbol
    jac::Matrix{Expression}
end

function DiffEqSystem(eqs, ivs, dvs, vs, ps)
    iv_name = ivs[1].subtype
    dv_name = dvs[1].subtype
    p_name = isempty(ps) ? :Parameter : ps[1].subtype
    DiffEqSystem(eqs, ivs, dvs, vs, ps, iv_name, dv_name, p_name, Matrix{Expression}(undef,0,0))
end

function DiffEqSystem(eqs; iv_name = :IndependentVariable,
                           dv_name = :DependentVariable,
                           v_name = :Variable,
                           p_name = :Parameter)
    targetmap =  Dict(iv_name => iv_name, dv_name => dv_name, v_name => v_name,
                       p_name => p_name)
    ivs, dvs, vs, ps = extract_elements(eqs, targetmap)
    DiffEqSystem(eqs, ivs, dvs, vs, ps, iv_name, dv_name, p_name, Matrix{Expression}(0,0))
end

function DiffEqSystem(eqs, ivs;
                      dv_name = :DependentVariable,
                      v_name = :Variable,
                      p_name = :Parameter)
    targetmap =  Dict(dv_name => dv_name, v_name => v_name, p_name => p_name)
    dvs, vs, ps = extract_elements(eqs, targetmap)
    DiffEqSystem(eqs, ivs, dvs, vs, ps, ivs[1].subtype, dv_name, p_name, Matrix{Expression}(undef,0,0))
end

function generate_ode_function(sys::DiffEqSystem;version = ArrayFunction)
    var_exprs = [:($(sys.dvs[i].name) = u[$i]) for i in 1:length(sys.dvs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]
    sys_exprs = build_equals_expr.(sys.eqs)
    if version == ArrayFunction
      dvar_exprs = [:(du[$i] = $(Symbol("$(sys.dvs[i].name)_$(sys.ivs[1].name)"))) for i in 1:length(sys.dvs)]
      exprs = vcat(var_exprs,param_exprs,sys_exprs,dvar_exprs)
      block = expr_arr_to_block(exprs)
      :((du,u,p,t)->$(toexpr(block)))
    elseif version == SArrayFunction
      dvar_exprs = [:($(Symbol("$(sys.dvs[i].name)_$(sys.ivs[1].name)"))) for i in 1:length(sys.dvs)]
      svector_expr = quote
        E = eltype(tuple($(dvar_exprs...)))
        T = StaticArrays.similar_type(typeof(u), E)
        T($(dvar_exprs...))
      end
      exprs = vcat(var_exprs,param_exprs,sys_exprs,svector_expr)
      block = expr_arr_to_block(exprs)
      :((u,p,t)->$(toexpr(block)))
    end
end

isintermediate(eq) = eq.args[1].diff == nothing

function build_equals_expr(eq)
    @assert typeof(eq.args[1]) <: Variable
    if !(isintermediate(eq))
        # Differential statement
        :($(Symbol("$(eq.args[1].name)_$(eq.args[1].diff.x.name)")) = $(eq.args[2]))
    else
        # Intermediate calculation
        :($(Symbol("$(eq.args[1].name)")) = $(eq.args[2]))
    end
end

function calculate_jacobian(sys::DiffEqSystem,simplify=true)
    diff_idxs = map(eq->eq.args[1].diff !=nothing,sys.eqs)
    diff_exprs = sys.eqs[diff_idxs]
    rhs = [eq.args[2] for eq in diff_exprs]
    # Handle intermediate calculations by substitution
    calcs = sys.eqs[.!(diff_idxs)]
    for i in 1:length(calcs)
        find_replace!.(rhs,calcs[i].args[1],calcs[i].args[2])
    end
    sys_exprs = calculate_jacobian(rhs,sys.dvs)
    sys_exprs = Expression[expand_derivatives(expr) for expr in sys_exprs]
    if simplify
        sys_exprs = Expression[simplify_constants(expr) for expr in sys_exprs]
    end
    sys_exprs
end

function generate_ode_jacobian(sys::DiffEqSystem,simplify=true)
    var_exprs = [:($(sys.dvs[i].name) = u[$i]) for i in 1:length(sys.dvs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]
    diff_idxs = map(eq->eq.args[1].diff !=nothing,sys.eqs)
    diff_exprs = sys.eqs[diff_idxs]
    jac = calculate_jacobian(sys,simplify)
    sys.jac = jac
    jac_exprs = [:(J[$i,$j] = $(convert(Expr, jac[i,j]))) for i in 1:size(jac,1), j in 1:size(jac,2)]
    exprs = vcat(var_exprs,param_exprs,vec(jac_exprs))
    block = expr_arr_to_block(exprs)
    :((J,u,p,t)->$(block))
end

function generate_ode_iW(sys::DiffEqSystem,simplify=true)
    var_exprs = [:($(sys.dvs[i].name) = u[$i]) for i in 1:length(sys.dvs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]
    diff_idxs = map(eq->eq.args[1].diff !=nothing,sys.eqs)
    diff_exprs = sys.eqs[diff_idxs]
    jac = sys.jac

    gam = Variable(:gam)

    W = LinearAlgebra.I - gam*jac
    W = SMatrix{size(W,1),size(W,2)}(W)
    iW = inv(W)

    if simplify
        iW = simplify_constants.(iW)
    end

    W = inv(LinearAlgebra.I/gam - jac)
    W = SMatrix{size(W,1),size(W,2)}(W)
    iW_t = inv(W)
    if simplify
        iW_t = simplify_constants.(iW_t)
    end

    iW_exprs = [:(iW[$i,$j] = $(convert(Expr, iW[i,j]))) for i in 1:size(iW,1), j in 1:size(iW,2)]
    exprs = vcat(var_exprs,param_exprs,vec(iW_exprs))
    block = expr_arr_to_block(exprs)

    iW_t_exprs = [:(iW[$i,$j] = $(convert(Expr, iW_t[i,j]))) for i in 1:size(iW_t,1), j in 1:size(iW_t,2)]
    exprs = vcat(var_exprs,param_exprs,vec(iW_t_exprs))
    block2 = expr_arr_to_block(exprs)
    :((iW,u,p,gam,t)->$(block)),:((iW,u,p,gam,t)->$(block2))
end

function DiffEqBase.ODEFunction(sys::DiffEqSystem;version = ArrayFunction,kwargs...)
    expr = generate_ode_function(sys;version=version,kwargs...)
    if version == ArrayFunction
      ODEFunction{true}(eval(expr))
    elseif version == SArrayFunction
      ODEFunction{false}(eval(expr))
    end
end

export DiffEqSystem, ODEFunction
export generate_ode_function
