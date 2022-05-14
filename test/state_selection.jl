using ModelingToolkit

@variables t
sts = @variables x1(t) x2(t) x3(t) x4(t)
params = @parameters u1(t) u2(t) u3(t) u4(t)
D = Differential(t)
eqs = [
    x1 + x2 + u1 ~ 0
    x1 + x2 + x3 + u2 ~ 0
    x1 + D(x3) + x4 + u3 ~ 0
    2*D(D(x1)) + D(D(x2)) + D(D(x3)) + D(x4) + u4 ~ 0
]
@named sys = ODESystem(eqs, t)

let dd = dummy_derivative(sys)
    has_dx1 = has_dx2 = false
    for eq in equations(dd)
        vars = ModelingToolkit.vars(eq)
        has_dx1 |= D(x1) in vars || D(D(x1)) in vars
        has_dx2 |= D(x2) in vars || D(D(x2)) in vars
    end
    @test has_dx1 ⊻ has_dx2 # only one of x1 and x2 can be a dummy derivative
    @test length(states(dd)) == length(equations(dd)) == 9
    @test length(states(structural_simplify(dd))) <= 6
end

let pss = partial_state_selection(sys)
    @test length(equations(pss)) == 1
    @test length(states(pss)) == 2
    @test length(equations(ode_order_lowering(pss))) == 2
end

@parameters t σ ρ β
@variables x(t) y(t) z(t) a(t) u(t) F(t)
D = Differential(t)

eqs = [
       D(x) ~ σ*(y-x)
       D(y) ~ x*(ρ-z)-y + β
       0 ~ z - x + y
       0 ~ a + z
       u ~ z + a
    ]

lorenz1 = ODESystem(eqs,t,name=:lorenz1)
let al1 = alias_elimination(lorenz1)
    let lss = partial_state_selection(al1)
        @test length(equations(lss)) == 2
    end
end
