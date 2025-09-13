using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using OrdinaryDiffEq
using Plots
using Test
using StaticArrays



# -----------------------------------------
# ----- example ---------------------------
# -----------------------------------------

vars = @variables begin
    x(t)=1, [input=true]

    # states
    y(t) = 0
end

eqs = [  
    # equations
    D(y) ~ x

]

@mtkcompile sys = System(eqs, t, vars, []) inputs=[x]
sys, set_input!, finalize! = ModelingToolkit.setup_inputs(sys);
prob = ODEProblem(sys, [], (0, 4))

# indeterminate form -----------------------

integrator = init(prob, Tsit5())

set_input!(integrator, sys.x, 1.0)
step!(integrator, 1.0, true)

set_input!(integrator, sys.x, 2.0)
step!(integrator, 1.0, true)

set_input!(integrator, sys.x, 3.0)
step!(integrator, 1.0, true)

set_input!(integrator, sys.x, 4.0)
step!(integrator, 1.0, true)

finalize!(integrator)

@test integrator.sol(0.0; idxs=sys.x) == 1.0
@test integrator.sol(1.0; idxs=sys.x) == 2.0
@test integrator.sol(2.0; idxs=sys.x) == 3.0
@test integrator.sol(3.0; idxs=sys.x) == 4.0
@test integrator.sol(4.0; idxs=sys.y) ≈ 10.0

# determinate form -----------------------
input = ModelingToolkit.Input(sys.x, SA[1,2,3,4], SA[0,1,2,3])
sol = solve(prob, input, Tsit5(); input_funs = (set_input!, finalize!));


@test sol(0.0; idxs=sys.x) == 1.0
@test sol(1.0; idxs=sys.x) == 2.0
@test sol(2.0; idxs=sys.x) == 3.0
@test sol(3.0; idxs=sys.x) == 4.0
@test sol(4.0; idxs=sys.y) ≈ 10.0