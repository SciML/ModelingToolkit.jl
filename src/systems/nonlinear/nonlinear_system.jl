export NonlinearSystem


struct NonlinearSystem <: AbstractSystem
    eqs::Vector{Equation}
    vs::Vector{Variable}
    ps::Vector{Variable}
end

function NonlinearSystem(eqs)
    vs, ps = extract_elements(eqs, [_subtype(:Unknown), _subtype(:Parameter)])
    NonlinearSystem(eqs, vs, ps)
end


function calculate_jacobian(sys::NonlinearSystem, simplify=true)
    sys_eqs, calc_eqs = system_eqs(sys), filter(iscalc, sys.eqs)
    rhs = [eq.rhs for eq in sys_eqs]

    for calc_eq ∈ calc_eqs
        find_replace!.(rhs, calc_eq.lhs, calc_eq.rhs)
    end

    sys_exprs = calculate_jacobian(rhs,sys.vs)
    sys_exprs = Expression[expand_derivatives(expr) for expr in sys_exprs]
    sys_exprs
end

iscalc(eq) = !isequal(eq.lhs, Constant(0))

system_eqs(sys::NonlinearSystem) = filter(!iscalc, sys.eqs)
system_extras(sys::NonlinearSystem) = filter(eq -> isa(eq.lhs, Variable), sys.eqs)
system_vars(sys::NonlinearSystem) = sys.vs
system_params(sys::NonlinearSystem) = sys.ps
