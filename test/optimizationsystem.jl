using ModelingToolkit, SparseArrays, Test, GalacticOptim, GalacticOptimJL

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

f = OptimizationFunction(rosenbrock,GalacticOptim.AutoModelingToolkit())
prob = OptimizationProblem(f,x0,_p)
sol = solve(prob,Newton())

# issue #819
@testset "Combined system name collisions" begin
    sys2 = OptimizationSystem(loss, [x, y], [a, b], name = :sys1)
    @test_throws ArgumentError OptimizationSystem(loss2, [z], [β], systems = [sys1, sys2])
end

# observed variable handling
@variables OBS
@named sys2 = OptimizationSystem(loss, [x,y], [a,b]; observed=[OBS ~ x+y])
OBS2 = OBS
@test isequal(OBS2, @nonamespace sys2.OBS)
@unpack OBS = sys2
@test isequal(OBS2,OBS)
