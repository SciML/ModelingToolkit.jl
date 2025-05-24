using ModelingToolkit
using Test
MT = ModelingToolkit

@independent_variables t
@variables x
struct MyNLS <: MT.AbstractSystem
    name::Any
    systems::Any
end
tmp = independent_variables(MyNLS("sys", []))
@test tmp == []

struct MyTDS <: MT.AbstractSystem
    iv::Any
    name::Any
    systems::Any
end
iv = independent_variables(MyTDS(t, "sys", []))
@test all(isequal.(iv, [t]))

struct MyMVS <: MT.AbstractSystem
    ivs::Any
    name::Any
    systems::Any
end
ivs = independent_variables(MyMVS([t, x], "sys", []))
@test all(isequal.(ivs, [t, x]))
