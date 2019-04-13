using ModelingToolkit
using Test

# Define some variables
@parameters t σ ρ β, α(..)
@variables x(t) y(t) z(t)
@derivatives D'~t

eq1 = D(x) ~ σ*(y-x)
eq2 = D(x) ~ α(t-2)*(y-x)

find_parameter_calls(eq1.rhs,Variable[])
r = find_parameter_calls(eq2.rhs,Variable[])
!isempty(r) && r[1].name == :α
