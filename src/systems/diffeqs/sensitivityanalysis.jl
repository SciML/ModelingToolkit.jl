function calculate_paramderivs(de::DiffEqSystem)
    eqns = []
    for i in de.eqs
        push!(eqns,i.args[2])
    end
    derivs = calculate_jacobian(eqns, de.ps)
    derivs
end

function calculate_dvarderivs(de::DiffEqSystem)
    dvars = []
    for i in de.dvars
        push!(dvars,i)
    end
    derivs = calculate_jacobian(dvars, de.ps)
    derivs
end