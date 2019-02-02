struct NonlinearSystem <: AbstractSystem
    eqs::Vector{Equation}
    vs::Vector{Variable}
    ps::Vector{Variable}
end

function NonlinearSystem(eqs)
    vs, ps = extract_elements(eqs, [_subtype(:Unknown), _subtype(:Parameter)])
    NonlinearSystem(eqs, vs, ps)
end

iscalc(eq) = isequal(eq.lhs, Constant(0))

function generate_nlsys_function(sys::NonlinearSystem)
    sys_eqs, calc_eqs = partition(iscalc, sys.eqs)

    var_pairs   = [(u.name, :(u[$i])) for (i, u) ∈ enumerate(sys.vs)]
    param_pairs = [(p.name, :(p[$i])) for (i, p) ∈ enumerate(sys.ps)]
    calc_pairs  = [(eq.lhs.name, convert(Expr, eq.rhs)) for eq ∈ calc_eqs if isa(eq.lhs, Variable)]
    (ls, rs) = collect(zip(var_pairs..., param_pairs..., calc_pairs...))

    var_eqs = Expr(:(=), build_expr(:tuple, ls), build_expr(:tuple, rs))
    sys_exprs = build_expr(:tuple, [convert(Expr, eq.rhs) for eq ∈ sys_eqs])
    let_expr = Expr(:let, var_eqs, sys_exprs)

    :((du,u,p) -> du .= $let_expr)
end

function calculate_jacobian(sys::NonlinearSystem,simplify=true)
    sys_eqs, calc_eqs = partition(iscalc, sys.eqs)
    rhs = [eq.rhs for eq in sys_eqs]

    for calc_eq ∈ calc_eqs
        find_replace!.(rhs, calc_eq.lhs, calc_eq.rhs)
    end

    sys_exprs = calculate_jacobian(rhs,sys.vs)
    sys_exprs = Expression[expand_derivatives(expr) for expr in sys_exprs]
    sys_exprs
end

system_vars(sys::NonlinearSystem) = sys.vs
system_params(sys::NonlinearSystem) = sys.ps

export NonlinearSystem
export generate_nlsys_function
