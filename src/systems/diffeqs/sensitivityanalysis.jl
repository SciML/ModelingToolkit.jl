calculate_paramderivs(de::DiffEqSystem) = calculate_jacobian(de,de.ps)
calculate_dvarderivs(de::DiffEqSystem) = calculate_jacobian(de.dvs,de.ps)