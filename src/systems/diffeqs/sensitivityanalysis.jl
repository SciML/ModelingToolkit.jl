function calculate_paramderivs(de::DiffEqSystem)
    derivs = []
    for i in 1:length(de.ps)
        deriv_param = []
        for j in 1:length(de.eqs)
            @Deriv D'~de.ps[i]
            deriv_eq = de.eqs[j].args[2]
            deriv_eq = D*deriv_eq
            push!(deriv_param,expand_derivatives(deriv_eq))
        end
        push!(derivs,deriv_param)
    end
    derivs
end

function calculate_dvarderivs(de::DiffEqSystem)
    derivs = []
    for i in 1:length(de.ps)
        deriv_param = []
        @Deriv D'~de.ps[i]
        for j in 1:length(de.dvs)
            l = DependentVariable(de.dvs[j].name,dependents=[de.ps[i]])
            push!(deriv_param,D*l)
        end
        push!(derivs,deriv_param)
    end
    derivs
end