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
    sys_exprs = build_equals_expr.(sys.eqs)
    dvar_exprs = [:(du[$i] = $(Symbol("$(sys.dvs[i].name)_$(sys.ivs[1].name)"))) for i in 1:length(sys.dvs)]
    exprs = vcat(var_exprs,param_exprs,sys_exprs,dvar_exprs)
    block = expr_arr_to_block(exprs)
    :((du,u,p,t)->$(block))
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

function generate_ode_jacobian(sys::DiffEqSystem,simplify=true)
    var_exprs = [:($(sys.dvs[i].name) = u[$i]) for i in 1:length(sys.dvs)]
    param_exprs = [:($(sys.ps[i].name) = p[$i]) for i in 1:length(sys.ps)]
    diff_idxs = map(eq->eq.args[1].diff !=nothing,sys.eqs)
    diff_exprs = sys.eqs[diff_idxs]
    rhs = [eq.args[2] for eq in diff_exprs]
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

function DiffEqBase.DiffEqFunction(sys::DiffEqSystem)
    expr = generate_ode_function(sys)
    DiffEqFunction{true}(eval(expr))
end

export DiffEqSystem, DiffEqFunction
export generate_ode_function
