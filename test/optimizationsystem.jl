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

rosenbrock(x, p) =  (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
_p  = [1.0, 100.0]

f = OptimizationFunction(rosenbrock,ModelingToolkit.AutoModelingToolkit(),x0,_p,grad=true,hess=true)
prob = OptimizationProblem(f,x0,_p)
sol = solve(prob,Optim.Newton())

# issue #819
@testset "Combined system name collisions" begin
    sys2 = OptimizationSystem(loss,[x,y],[a,b],name=:sys1)
    @test_throws ArgumentError OptimizationSystem(loss2,[z],[β],systems=[sys1,sys2])
end
