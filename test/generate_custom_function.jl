using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@variables x(t) y(t)[1:3]
@parameters p1=1.0 p2[1:3]=[1.0, 2.0, 3.0] p3::Int=1 p4::Bool=false

sys = complete(ODESystem(Equation[], t, [x; y], [p1, p2, p3, p4]; name = :sys))
u0 = [1.0, 2.0, 3.0, 4.0]
p = ModelingToolkit.MTKParameters(sys, [])

fn1 = generate_custom_function(sys, x + y[1] + p1 + p2[1] + p3 * t; expression = Val(false))
@test fn1(u0, p, 0.0) == 5.0

fn2 = generate_custom_function(
    sys, x + y[1] + p1 + p2[1] + p3 * t, [x], [p1, p2, p3]; expression = Val(false))
@test fn1(u0, p, 0.0) == 5.0

fn3_oop, fn3_iip = generate_custom_function(
    sys, [x + y[2], y[3] + p2[2], p1 + p3, 3t]; expression = Val(false))

buffer = zeros(4)
fn3_iip(buffer, u0, p, 1.0)
@test buffer == [4.0, 6.0, 2.0, 3.0]
@test fn3_oop(u0, p, 1.0) == [4.0, 6.0, 2.0, 3.0]

fn4 = generate_custom_function(sys, ifelse(p4, p1, p2[2]); expression = Val(false))
@test fn4(u0, p, 1.0) == 2.0
fn5 = generate_custom_function(sys, ifelse(!p4, p1, p2[2]); expression = Val(false))
@test fn5(u0, p, 1.0) == 1.0

@variables x y[1:3]
sys = complete(NonlinearSystem(Equation[], [x; y], [p1, p2, p3, p4]; name = :sys))

fn1 = generate_custom_function(sys, x + y[1] + p1 + p2[1] + p3; expression = Val(false))
@test fn1(u0, p) == 6.0

fn2 = generate_custom_function(
    sys, x + y[1] + p1 + p2[1] + p3, [x], [p1, p2, p3]; expression = Val(false))
@test fn1(u0, p) == 6.0

fn3_oop, fn3_iip = generate_custom_function(
    sys, [x + y[2], y[3] + p2[2], p1 + p3]; expression = Val(false))

buffer = zeros(3)
fn3_iip(buffer, u0, p)
@test buffer == [4.0, 6.0, 2.0]
@test fn3_oop(u0, p, 1.0) == [4.0, 6.0, 2.0]

fn4 = generate_custom_function(sys, ifelse(p4, p1, p2[2]); expression = Val(false))
@test fn4(u0, p, 1.0) == 2.0
fn5 = generate_custom_function(sys, ifelse(!p4, p1, p2[2]); expression = Val(false))
@test fn5(u0, p, 1.0) == 1.0
