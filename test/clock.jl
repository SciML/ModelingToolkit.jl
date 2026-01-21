using ModelingToolkit, Test, Setfield, OrdinaryDiffEq, DiffEqCallbacks
using OrderedCollections
using ModelingToolkit: ContinuousClock
using ModelingToolkit: t_nounits as t, D_nounits as D
using Symbolics, SymbolicUtils
using Symbolics: SymbolicT, VartypeT
import ModelingToolkitTearing as MTKTearing

function infer_clocks(sys)
    ts = TearingState(sys)
    ci = MTKTearing.ClockInference(ts)
    return MTKTearing.infer_clocks!(ci), Dict(ci.ts.fullvars .=> ci.var_domain)
end

@info "Testing hybrid system"
dt = 0.1
@variables x(t) y(t) u(t) yd(t) ud(t) r(t)
@parameters kp
# u(n + 1) := f(u(n))

eqs = [
    yd ~ Sample(dt)(y)
    ud ~ kp * (r - yd)
    r ~ 1.0

    # plant (time continuous part)
    u ~ Hold(ud)
    D(x) ~ -x + u
    y ~ x
]
@named sys = System(eqs, t)
# compute equation and variables' time domains
#TODO: test linearize

#=
 Differential(t)(x(t)) ~ u(t) - x(t)
 0 ~ Sample(Clock(t, 0.1))(y(t)) - yd(t)
 0 ~ kp*(r(t) - yd(t)) - ud(t)
 0 ~ Hold()(ud(t)) - u(t)
 0 ~ x(t) - y(t)

====
By inference:

 Differential(t)(x(t)) ~ u(t) - x(t)
 0 ~ Hold()(ud(t)) - u(t) # Hold()(ud(t)) is constant except in an event
 0 ~ x(t) - y(t)

 0 ~ Sample(Clock(t, 0.1))(y(t)) - yd(t)
 0 ~ kp*(r(t) - yd(t)) - ud(t)

====

 Differential(t)(x(t)) ~ u(t) - x(t)
 0 ~ Hold()(ud(t)) - u(t)
 0 ~ x(t) - y(t)

 yd(t) := Sample(Clock(t, 0.1))(y(t))
 ud(t) := kp*(r(t) - yd(t))
=#

#=
     D(x) ~ Shift(x, 0, dt) + 1 # this should never meet with continuous variables
=>   (Shift(x, 0, dt) - Shift(x, -1, dt))/dt ~ Shift(x, 0, dt) + 1
=>   Shift(x, 0, dt) - Shift(x, -1, dt) ~ Shift(x, 0, dt) * dt + dt
=>   Shift(x, 0, dt) - Shift(x, 0, dt) * dt ~ Shift(x, -1, dt) + dt
=>   (1 - dt) * Shift(x, 0, dt) ~ Shift(x, -1, dt) + dt
=>   Shift(x, 0, dt) := (Shift(x, -1, dt) + dt) / (1 - dt) # Discrete system
=#

ci, varmap = infer_clocks(sys)
eqmap = ci.eq_domain
tss, inputs, continuous_id = MTKTearing.split_system(deepcopy(ci))
sss = ModelingToolkit._mtkcompile!(
    deepcopy(tss[continuous_id]), inputs = OrderedSet{SymbolicT}(inputs[continuous_id])
)
@test equations(sss) == [D(x) ~ u - x]
sss = ModelingToolkit._mtkcompile!(
    deepcopy(tss[1]), inputs = OrderedSet{SymbolicT}(inputs[1])
)
@test isempty(equations(sss))
d = Clock(dt)
k = ShiftIndex(d)
@test issetequal(
    observed(sss),
    [
        yd ~ Sample(dt)(y); r ~ 1.0;
        ud ~ kp * (r - yd)
    ]
)

canonical_eqs = map(eqs) do eq
    if iscall(eq.lhs) && operation(eq.lhs) isa Differential
        return eq
    else
        return 0 ~ eq.rhs - eq.lhs
    end
end
eqs_idxs = findfirst.(isequal.(canonical_eqs), (equations(ci.ts),))
d = Clock(dt)
# Note that TearingState reorders the equations
@test eqmap[eqs_idxs[1]] == d
@test eqmap[eqs_idxs[2]] == d
@test eqmap[eqs_idxs[3]] == d
@test eqmap[eqs_idxs[4]] == ContinuousClock()
@test eqmap[eqs_idxs[5]] == ContinuousClock()
@test eqmap[eqs_idxs[6]] == ContinuousClock()

@test varmap[yd] == d
@test varmap[ud] == d
@test varmap[r] == d
@test varmap[x] == ContinuousClock()
@test varmap[y] == ContinuousClock()
@test varmap[u] == ContinuousClock()

@info "Testing shift normalization"
dt = 0.1
@variables x(t) y(t) u(t) yd(t) ud(t)
@parameters kp
d = Clock(dt)
k = ShiftIndex(d)

eqs = [
    yd ~ Sample(dt)(y)
    ud ~ kp * yd + ud(k - 2)

    # plant (time continuous part)
    u ~ Hold(ud)
    D(x) ~ -x + u
    y ~ x
]
@named sys = System(eqs, t)
@test_throws ModelingToolkit.HybridSystemNotSupportedException ss = mtkcompile(sys);

@testset "Clock inference uses and splits initialization equations" begin
    @variables x(t) y(t) z(t)
    k = ShiftIndex()
    clk = Clock(0.1)
    eqs = [D(x) ~ x, y ~ Sample(clk)(x), z ~ z(k - 1) + 1]
    initialization_eqs = [y ~ z, x ~ 1]
    @named sys = System(eqs, t; initialization_eqs)
    ts = TearingState(sys)
    ci = MTKTearing.ClockInference(ts)
    @test length(ci.init_eq_domain) == 2
    MTKTearing.infer_clocks!(ci)
    canonical_eqs = map(eqs) do eq
        if iscall(eq.lhs) && operation(eq.lhs) isa Differential
            return eq
        else
            return 0 ~ eq.rhs - eq.lhs
        end
    end
    eqs_idxs = findfirst.(isequal.(canonical_eqs), (equations(ci.ts),))
    @test ci.eq_domain[eqs_idxs[1]] == ContinuousClock()
    @test ci.eq_domain[eqs_idxs[2]] == clk
    @test ci.eq_domain[eqs_idxs[3]] == clk
    varmap = Dict(ci.ts.fullvars .=> ci.var_domain)
    @test varmap[x] == ContinuousClock()
    @test varmap[y] == clk
    @test varmap[z] == clk
end

struct ZeroArgOp <: Symbolics.Operator end
(o::ZeroArgOp)() = SymbolicUtils.Term{VartypeT}(o, Any[]; type = Bool, shape = [])
SymbolicUtils.promote_symtype(::ZeroArgOp, T) = Union{Bool, T}
SymbolicUtils.isbinop(::ZeroArgOp) = false
Base.nameof(::ZeroArgOp) = :ZeroArgOp
MTKTearing.input_timedomain(::ZeroArgOp, _ = nothing) = MTKTearing.InputTimeDomainElT[]
MTKTearing.output_timedomain(::ZeroArgOp, _ = nothing) = Clock(0.1)
ModelingToolkit.validate_operator(::ZeroArgOp, args, iv; context = nothing) = nothing
SciMLBase.is_discrete_time_domain(::ZeroArgOp) = true

@testset "Zero-argument clock operators" begin
    @variables x(t) y(t)
    clk = Clock(0.1)
    eqs = [
        D(x) ~ x
        y ~ ZeroArgOp()()
    ]
    @named sys = System(eqs, t)
    @test issetequal(unknowns(sys), [x, y])
    ts = TearingState(sys)
    @test issetequal(ts.fullvars, [D(x), x, y, ZeroArgOp()()])
    ci, clkmap = infer_clocks(sys)
    @test clkmap[ZeroArgOp()()] == clk
end
