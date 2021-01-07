"""
$(TYPEDSIGNATURES)

Generates the Liouville transformed set of ODEs, which is the original
ODE system with a new variable `trJ` appended, corresponding to the
-tr(Jacobian). This variable is used for properties like uncertainty
propagation and should be used with a zero initial condition.
"""
function liouville_transform(sys)
      t = independent_variable(sys)
      @variables trJ
      D = ModelingToolkit.Differential(t)
      neweq = D(trJ) ~ -tr(calculate_jacobian(sys))
      neweqs = [equations(sys);neweq]
      vars = [states(sys);trJ]
      ODESystem(neweqs,t,vars,parameters(sys))
end
