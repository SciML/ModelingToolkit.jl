using ModelingToolkit, SparseArrays, Test, GalacticOptim, Optim

@variables x y
@parameters a b
loss = (a - x)^2 + b * (y - x^2)^2
sys1 = OptimizationSystem(loss,[x,y],[a,b],name=:sys1)
sys2 = OptimizationSystem(loss,[x,y],[a,b],name=:sys2)

@variables z
@parameters β
loss2 = sys1.x - sys2.y + z*β
combinedsys = OptimizationSystem(loss2,[z],[β],systems=[sys1,sys2],name=:combinedsys)

equations(combinedsys)
states(combinedsys)
parameters(combinedsys)

calculate_gradient(combinedsys)
calculate_hessian(combinedsys)
generate_function(combinedsys)
generate_gradient(combinedsys)
generate_hessian(combinedsys)
ModelingToolkit.hessian_sparsity(combinedsys)

u0 = [
    sys1.x=>1.0
    sys1.y=>2.0
    sys2.x=>3.0
    sys2.y=>4.0
    z=>5.0
]
p = [
    sys1.a => 6.0
    sys1.b => 7.0
    sys2.a => 8.0
    sys2.b => 9.0
    β => 10.0
]

prob = OptimizationProblem(combinedsys,u0,p,grad=true)
sol = solve(prob,NelderMead())
@test sol.minimum < -1e5

prob2 = remake(prob,u0=sol.minimizer)
sol = solve(prob,BFGS(initial_stepnorm=0.0001),allow_f_increases=true)
@test sol.minimum < -1e8
sol = solve(prob2,BFGS(initial_stepnorm=0.0001),allow_f_increases=true)
@test sol.minimum < -1e9
