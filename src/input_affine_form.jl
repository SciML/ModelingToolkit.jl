"""
    input_affine_form(eqs, inputs)

Extract control-affine (input-affine) form from symbolic equations.

Given a system of equations of the form `ẋ = F(x, u)` where `x` is the state
and `u` are the inputs, this function extracts the control-affine representation:
`ẋ = f(x) + g(x)*u`

where:
- `f(x)` is the drift term (dynamics without inputs)
- `g(x)` is the input matrix

# Arguments
- `eqs`: Vector of symbolic equations representing the system dynamics
- `inputs`: Vector of input/control variables

# Returns
- `f`: The drift vector f(x)
- `g`: The input matrix g(x)

# Example
```julia
using ModelingToolkit, Symbolics

@variables x1 x2 u1 u2
state  = [x1, x2]
inputs = [u1, u2]

# Define system equations: ẋ = F(x, u)
eqs = [
    -x1 + 2*x2 + u1,
    x1*x2 - x2 + u1 + 2*u2
]

# Extract control-affine form
f, g = input_affine_form(eqs, inputs)
```

# Notes
The function assumes that the equations are affine in the inputs. If the equations
are nonlinear in the inputs, the result may not be meaningful.
"""
function input_affine_form(eqs, inputs)
    # Extract the input matrix g(x) by taking coefficients of each input
    g = [Symbolics.coeff(Symbolics.simplify(eq, expand=true), u) for eq in eqs, u in inputs]
    g = Symbolics.simplify.(g, expand=true)
    
    # Extract the drift term f(x) by substituting inputs = 0
    f = Symbolics.substitute.(eqs, Ref(Dict(inputs .=> 0)))
    f = Symbolics.simplify.(f, expand=true)
    
    return f, g
end