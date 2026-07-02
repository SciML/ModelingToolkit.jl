using ModelingToolkitBase, FunctionProperties, Test
using ModelingToolkitBase: t_nounits as t, D_nounits as D

# `getindex(::MTKParameters, ::Int)` selects a parameter buffer by index. Its branch is on the
# index (constant at every real generated-code call site), never on values, so `hasbranching`
# should not report it. The MTKFunctionPropertiesExt extension marks it a leaf via `is_leaf_sig`.

@parameters a = 2.0
@variables x(t) = 1.0
@named sys = System([D(x) ~ a * x], t)
sys = mtkcompile(sys)
prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0))
p = prob.p

@test p isa ModelingToolkitBase.MTKParameters

# The extension's signature-level exemption is active for the parameter-container indexing...
@test FunctionProperties.is_leaf_sig(Tuple{typeof(getindex), typeof(p), Int})
# ...and does not over-reach to ordinary containers.
@test !FunctionProperties.is_leaf_sig(Tuple{typeof(getindex), Vector{Float64}, Int})

# A function that only indexes into the parameter container therefore has no value-dependent branch.
index_param(pp) = @inbounds pp[1]
@test !hasbranching(index_param, p)
