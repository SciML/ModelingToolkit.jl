"""
$(TYPEDSIGNATURES)

Generates the Liouville transformed set of ODEs, which is the original
ODE system with a new variable `trJ` appended, corresponding to the
-tr(Jacobian). This variable is used for properties like uncertainty
propagation and should be used with a zero initial condition.



Sources:

Probabilistic Robustness Analysis of F-16 Controller Performance: An
Optimal Transport Approach

Abhishek Halder, Kooktae Lee, and Raktim Bhattacharya
https://abhishekhalder.bitbucket.io/F16ACC2013Final.pdf


"""
function liouville_transform(sys)
      t = independent_variable(sys)
      @variables trJ
      D = ModelingToolkit.Differential(t)
      neweq = D(trJ) ~ trJ*-tr(calculate_jacobian(sys))
      neweqs = [equations(sys);neweq]
      vars = [states(sys);trJ]
      ODESystem(neweqs,t,vars,parameters(sys))
end
