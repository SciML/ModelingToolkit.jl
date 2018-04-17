calculate_paramderivs(de::DiffEqSystem) = calculate_jacobian(de,de.ps)
function calculate_dvarderivs(de::DiffEqSystem)
    params = de.ps
    dvars = de.dvs
    new_params = [DependentVariable(i.name,dependents=de.ivs) for i in de.ps]
    calculate_jacobian(dvars,new_params)
end

function get_sensitivity(de::DiffEqSystem)
    jac = calculate_jacobian(de)
    param_jac = calculate_paramderivs(de)
    dvar_jac = calculate_dvarderivs(de)
    sensitivity = Operation[]
    for i in 1:length(de.ps)
        push!(sensitivity,jac*dvar_jac[i]+param_jac[i])
    end
    sensitivity
end 