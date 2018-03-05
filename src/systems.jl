struct DiffEqSystem <: AbstractSystem
    eqs::Vector{Operation}
    ivs::Vector{Variable}
    dvs::Vector{Variable}
    vs::Vector{Variable}
    ps::Vector{Variable}
end
function DiffEqSystem(eqs)
    ivs, dvs, vs, ps = extract_elements(eqs, (:IndependentVariable, :DependentVariable, :Variable, :Parameter))
    DiffEqSystem(eqs, ivs, dvs, vs, ps)
end
function DiffEqSystem(eqs, ivs)
    dvs, vs, ps = extract_elements(eqs, (:DependentVariable, :Variable, :Parameter))
    DiffEqSystem(eqs, ivs, dvs, vs, ps)
end

function generate_ode_function(sys::DiffEqSystem)
    var_exprs = [:($(sys.dvs[i].name) = u[$i]) for i in 1:length(sys.dvs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]
    sys_exprs = map(eq->:($(Symbol("$(eq.args[1].name)_$(eq.args[1].diff.x.name)")) = $(eq.args[2])),sys.eqs)
    dvar_exprs = [:(du[$i] = $(Symbol("$(sys.dvs[i].name)_$(sys.ivs[1].name)"))) for i in 1:length(sys.dvs)]
    exprs = vcat(var_exprs,param_exprs,sys_exprs,dvar_exprs)
    block = expr_arr_to_block(exprs)
    :((du,u,p,t)->$(block))
end

function generate_ode_jacobian(sys::DiffEqSystem,simplify=true)
    var_exprs = [:($(sys.dvs[i].name) = u[$i]) for i in 1:length(sys.dvs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]
    rhs = [eq.args[2] for eq in sys.eqs]
    sys_exprs = calculate_jacobian(rhs,sys.dvs)
    sys_exprs = expand_derivatives.(sys_exprs)
    if simplify
        sys_exprs = simplify_constants.(sys_exprs)
    end
    sys_exprs
end

function DiffEqBase.DiffEqFunction(sys::DiffEqSystem)
    expr = generate_ode_function(sys)
    DiffEqFunction{true}(eval(expr))
end

struct NonlinearSystem <: AbstractSystem
    eqs::Vector{Operation}
    vs::Vector{Variable}
    ps::Vector{Variable}
end

function NonlinearSystem(eqs)
    # Allow the use of :DependentVariable to make it seamless with DE use
    dvs, vs, ps = extract_elements(eqs, (:DependentVariable, :Variable, :Parameter))
    vs = [dvs;vs]
    NonlinearSystem(eqs, vs, ps)
end

function generate_nlsys_function(sys::NonlinearSystem)
    var_exprs = [:($(sys.vs[i].name) = u[$i]) for i in 1:length(sys.vs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]
    sys_exprs = [:($(Symbol("resid[$i]")) = $(sys.eqs[i].args[2])) for i in eachindex(sys.eqs)]
    exprs = vcat(var_exprs,param_exprs,sys_exprs)
    block = expr_arr_to_block(exprs)
    :((du,u,p)->$(block))
end

function generate_nlsys_jacobian(sys::NonlinearSystem,simplify=true)
    var_exprs = [:($(sys.vs[i].name) = u[$i]) for i in 1:length(sys.vs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]
    rhs = [eq.args[2] for eq in sys.eqs]
    sys_exprs = calculate_jacobian(rhs,sys.vs)
    sys_exprs = expand_derivatives.(sys_exprs)
    if simplify
        sys_exprs = simplify_constants.(sys_exprs)
    end
    sys_exprs
end

export DiffEqSystem, NonlinearSystem, DiffEqFunction
export generate_ode_function, generate_nlsys_function
