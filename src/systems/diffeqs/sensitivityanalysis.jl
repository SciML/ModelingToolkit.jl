calculate_paramderivs(de::DiffEqSystem) = calculate_jacobian(de,de.ps)
calculate_dvarderivs(de::DiffEqSystem) = calculate_jacobian(de.dvs,de.ps)

function get_sensitivity(de::DiffEqSystem)
    jac = calculate_jacobian(de)
    param_jac = calculate_paramderivs(de)
    dvar_jac = calculate_dvarderivs(de)
    sensitivity = []
    for i in 1:length(de.ps)
        push!(sensitivity,jac*dvar_jac[i]+param_jac[i])
    end
    sensitivity
end 