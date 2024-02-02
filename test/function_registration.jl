# Each test here builds an ODEFunction including some user-registered
# Operations. The test simply checks that calling the ODEFunction
# appropriately calls the registered functions, whether the call is
# qualified (with a module name) or not.

# TEST: Function registration in a module.
# ------------------------------------------------
module MyModule
using ModelingToolkit, DiffEqBase, LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as Dt
@parameters x
@variables u(t)

function do_something(a)
    a + 10
end
@register_symbolic do_something(a)

eq = Dt(u) ~ do_something(x) + MyModule.do_something(x)
@named sys = ODESystem([eq], t, [u], [x])
sys = complete(sys)
fun = ODEFunction(sys)

u0 = 5.0
@test fun([0.5], [u0], 0.0) == [do_something(u0) * 2]
end

# TEST: Function registration in a nested module.
# ------------------------------------------------
module MyModule2
module MyNestedModule
using ModelingToolkit, DiffEqBase, LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as Dt
@parameters x
@variables u(t)

function do_something_2(a)
    a + 20
end
@register_symbolic do_something_2(a)

eq = Dt(u) ~ do_something_2(x) + MyNestedModule.do_something_2(x)
@named sys = ODESystem([eq], t, [u], [x])
sys = complete(sys)
fun = ODEFunction(sys)

u0 = 3.0
@test fun([0.5], [u0], 0.0) == [do_something_2(u0) * 2]
end
end

# TEST: Function registration outside any modules.
# ------------------------------------------------
using ModelingToolkit, DiffEqBase, LinearAlgebra, Test
using ModelingToolkit: t_nounits as t, D_nounits as Dt
@parameters x
@variables u(t)

function do_something_3(a)
    a + 30
end
@register_symbolic do_something_3(a)

eq = Dt(u) ~ do_something_3(x) + (@__MODULE__).do_something_3(x)
@named sys = ODESystem([eq], t, [u], [x])
sys = complete(sys)
fun = ODEFunction(sys)

u0 = 7.0
@test fun([0.5], [u0], 0.0) == [do_something_3(u0) * 2]

# TEST: Function registration works with derivatives.
# ---------------------------------------------------
foo(x, y) = sin(x) * cos(y)
@variables x(t) y(t) z(t);
D = Dt
@register_symbolic foo(x, y)

using ModelingToolkit: value, arguments, operation
expr = value(foo(x, y))
@test operation(expr) === foo
@test arguments(expr)[1] === value(x)
@test arguments(expr)[2] === value(y)
ModelingToolkit.derivative(::typeof(foo), (x, y), ::Val{1}) = cos(x) * cos(y) # derivative w.r.t. the first argument
ModelingToolkit.derivative(::typeof(foo), (x, y), ::Val{2}) = -sin(x) * sin(y) # derivative w.r.t. the second argument
@test isequal(expand_derivatives(D(foo(x, y))), expand_derivatives(D(sin(x) * cos(y))))

# TEST: Function registration run from inside a function.
# -------------------------------------------------------
# This tests that we can get around the world age issue by falling back to
# GeneralizedGenerated instead of function expressions.
# Might be useful in cases where someone wants to define functions that build
# up and use ODEFunctions given some parameters.
function do_something_4(a)
    a + 30
end
@register_symbolic do_something_4(a)
function build_ode()
    @parameters x
    @variables u(t)
    eq = Dt(u) ~ do_something_4(x) + (@__MODULE__).do_something_4(x)
    @named sys = ODESystem([eq], t, [u], [x])
    sys = complete(sys)
    fun = ODEFunction(sys, eval_expression = false)
end
function run_test()
    fun = build_ode()
    u0 = 10.0
    @test fun([0.5], [u0], 0.0) == [do_something_4(u0) * 2]
end
run_test()

using ModelingToolkit: arguments
@variables a
@register_symbolic foo(x, y, z)
@test 1 * foo(a, a, a) * Num(1) isa Num
@test !any(x -> x isa Num, arguments(value(1 * foo(a, a, a) * Num(1))))
