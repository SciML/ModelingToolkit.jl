### Preparations ###

# Fetch packages.
using ModelingToolkit, Test
using ModelingToolkit: t_nounits as t, D_nounits as D
import ModelingToolkit: get_ps, get_unknowns, get_observed, get_eqs, get_continuous_events,
                        get_discrete_events, namespace_equations

# Creates helper functions.
function all_sets_equal(args...)
    for arg in args[2:end]
        issetequal(args[1], arg) || return false
    end
    return true
end
function sym_issubset(set1, set2)
    for sym1 in set1
        any(isequal(sym1, sym2) for sym2 in set2) || return false
    end
    return true
end

### Basic Tests ###

# Checks `toplevel = false` argument for various accessors (currently only for `ODESystem`s).
# Compares to `toplevel = true` version, and `get_` functions.
# Checks  accessors for parameters, unknowns, equations, observables, and events.
let
    # Prepares model components.
    @parameters p_top p_mid1 p_mid2 p_bot d
    @variables X_top(t) X_mid1(t) X_mid2(t) X_bot(t) Y(t) O(t)

    # Creates the systems (individual and hierarchical).
    eqs_top = [
        D(X_top) ~ p_top - d * X_top,
        D(Y) ~ log(X_top) - Y^2 + 3.0,
        O ~ (p_top + d) * X_top + Y
    ]
    eqs_mid1 = [
        D(X_mid1) ~ p_mid1 - d * X_mid1^2,
        D(Y) ~ D(X_mid1) - Y^3,
        O ~ (p_mid1 + d) * X_mid1 + Y
    ]
    eqs_mid2 = [
        D(X_mid2) ~ p_mid2 - d * X_mid2,
        X_mid2^3 ~ log(X_mid2 + Y) - Y^2 + 3.0,
        O ~ (p_mid2 + d) * X_mid2 + Y
    ]
    eqs_bot = [
        D(X_bot) ~ p_bot - d * X_bot,
        D(Y) ~ -Y^3,
        O ~ (p_bot + d) * X_bot + Y
    ]
    cevs = [[t ~ 1.0] => [Y ~ Y + 2.0]]
    devs = [(t == 2.0) => [Y ~ Y + 2.0]]
    @named sys_bot = ODESystem(
        eqs_bot, t; systems = [], continuous_events = cevs, discrete_events = devs)
    @named sys_mid2 = ODESystem(
        eqs_mid2, t; systems = [], continuous_events = cevs, discrete_events = devs)
    @named sys_mid1 = ODESystem(
        eqs_mid1, t; systems = [sys_bot], continuous_events = cevs, discrete_events = devs)
    @named sys_top = ODESystem(eqs_top, t; systems = [sys_mid1, sys_mid2],
        continuous_events = cevs, discrete_events = devs)
    sys_bot_comp = complete(sys_bot)
    sys_mid2_comp = complete(sys_mid2)
    sys_mid1_comp = complete(sys_mid1)
    sys_top_comp = complete(sys_top)
    sys_bot_ss = structural_simplify(sys_bot)
    sys_mid2_ss = structural_simplify(sys_mid2)
    sys_mid1_ss = structural_simplify(sys_mid1)
    sys_top_ss = structural_simplify(sys_top)

    # Checks `parameters` for `toplevel = false`.
    @test all_sets_equal(parameters.([sys_bot, sys_bot_comp, sys_bot_ss])..., [d, p_bot])
    @test all_sets_equal(parameters.([sys_mid1, sys_mid1_comp, sys_mid1_ss])...,
        [d, p_mid1, sys_bot.d, sys_bot.p_bot])
    @test all_sets_equal(
        parameters.([sys_mid2, sys_mid2_comp, sys_mid2_ss])..., [d, p_mid2])
    @test all_sets_equal(parameters.([sys_top, sys_top_comp, sys_top_ss])...,
        [d, p_top, sys_mid1.d, sys_mid1.p_mid1, sys_mid1.sys_bot.d,
            sys_mid1.sys_bot.p_bot, sys_mid2.d, sys_mid2.p_mid2])

    # Checks `parameters`` for `toplevel = true`. Compares to known parameters and also checks that
    # these are subset of what `get_ps` returns.
    @test all_sets_equal(
        parameters.([sys_bot, sys_bot_comp, sys_bot_ss]; toplevel = true)..., [d, p_bot])
    @test_broken all_sets_equal(
        parameters.([sys_mid1, sys_mid1_comp, sys_mid1_ss]; toplevel = true)...,
        [d, p_mid1])
    @test all_sets_equal(
        parameters.([sys_mid2, sys_mid2_comp, sys_mid2_ss]; toplevel = true)...,
        [d, p_mid2])
    @test_broken all_sets_equal(
        parameters.([sys_top, sys_top_comp, sys_top_ss]; toplevel = true)..., [d, p_top])
    @test all(sym_issubset(parameters(sys; toplevel = true), get_ps(sys))
    for sys in [sys_bot, sys_mid2, sys_mid1, sys_top])

    # Checks `unknowns` for `toplevel = false`. O(t) is eliminated by `structural_simplify` and
    # must be considered separately.
    @test all_sets_equal(unknowns.([sys_bot, sys_bot_comp])..., [O, Y, X_bot])
    @test all_sets_equal(unknowns.([sys_bot_ss])..., [Y, X_bot])
    @test all_sets_equal(unknowns.([sys_mid1, sys_mid1_comp])...,
        [O, Y, X_mid1, sys_bot.Y, sys_bot.O, sys_bot.X_bot])
    @test all_sets_equal(unknowns.([sys_mid1_ss])..., [Y, X_mid1, sys_bot.Y, sys_bot.X_bot])
    @test all_sets_equal(unknowns.([sys_mid2, sys_mid2_comp])..., [O, Y, X_mid2])
    @test all_sets_equal(unknowns.([sys_mid2_ss])..., [Y, X_mid2])
    @test all_sets_equal(unknowns.([sys_top, sys_top_comp])...,
        [O, Y, X_top, sys_mid1.O, sys_mid1.Y, sys_mid1.X_mid1,
            sys_mid1.sys_bot.O, sys_mid1.sys_bot.Y, sys_mid1.sys_bot.X_bot,
            sys_mid2.O, sys_mid2.Y, sys_mid2.X_mid2])
    @test all_sets_equal(unknowns.([sys_top_ss])...,
        [Y, X_top, sys_mid1.Y, sys_mid1.X_mid1, sys_mid1.sys_bot.Y,
            sys_mid1.sys_bot.X_bot, sys_mid2.Y, sys_mid2.X_mid2])

    # Checks `unknowns` for `toplevel = true`. Note that O is not eliminated here (as we go back
    # to original parent system). Also checks that outputs are subsets of what `get_ps` returns..
    @test all_sets_equal(
        unknowns.([sys_bot, sys_bot_comp, sys_bot_ss]; toplevel = true)..., [O, Y, X_bot])
    @test all_sets_equal(
        unknowns.([sys_mid1, sys_mid1_comp]; toplevel = true)..., [O, Y, X_mid1])
    @test all_sets_equal(
        unknowns.([sys_mid2, sys_mid2_comp]; toplevel = true)..., [O, Y, X_mid2])
    @test all_sets_equal(
        unknowns.([sys_top, sys_top_comp]; toplevel = true)..., [O, Y, X_top])
    @test all(sym_issubset(unknowns(sys; toplevel = true), get_unknowns(sys))
    for sys in [sys_bot, sys_mid1, sys_mid2, sys_top])

    # Checks `equations` for `toplevel = false`. Do not check ss equations as these might potentially
    # be structurally simplified to new equations.
    @test all_sets_equal(equations.([sys_bot, sys_bot_comp])..., eqs_bot)
    @test all_sets_equal(
        equations.([sys_mid1, sys_mid1_comp])..., [eqs_mid1; namespace_equations(sys_bot)])
    @test all_sets_equal(equations.([sys_mid2, sys_mid2_comp])..., eqs_mid2)
    @test all_sets_equal(equations.([sys_top, sys_top_comp])...,
        [eqs_top; namespace_equations(sys_mid1); namespace_equations(sys_mid2)])

    # Checks `equations` for `toplevel = true`. Do not check ss equations  directly as these
    # might potentially be structurally simplified to new equations.
    @test all_sets_equal(equations.([sys_bot, sys_bot_comp]; toplevel = true)..., eqs_bot)
    @test all_sets_equal(
        equations.([sys_mid1, sys_mid1_comp]; toplevel = true)..., eqs_mid1)
    @test all_sets_equal(
        equations.([sys_mid2, sys_mid2_comp]; toplevel = true)..., eqs_mid2)
    @test all_sets_equal(equations.([sys_top, sys_top_comp]; toplevel = true)..., eqs_top)
    @test all(sym_issubset(equations(sys; toplevel = true), get_eqs(sys))
    for sys in [sys_bot, sys_mid2, sys_mid1, sys_top])

    # Checks `observed for `toplevel = false`. For non-ss systems, there are no observables.
    @test all_sets_equal(
        observed.([sys_bot, sys_bot_comp, sys_mid1, sys_mid1_comp,
            sys_mid2, sys_mid2_comp, sys_top, sys_top_comp])...,
        [])
    @test issetequal(observed(sys_bot_ss), eqs_bot[3:3])
    @test issetequal(
        observed(sys_mid1_ss), [eqs_mid1[3:3]; namespace_equations(sys_bot)[3:3]])
    @test issetequal(observed(sys_mid2_ss), eqs_mid2[3:3])
    @test issetequal(observed(sys_top_ss),
        [eqs_top[3:3]; namespace_equations(sys_mid1)[3:3:6];
         namespace_equations(sys_mid2)[3:3]])

    # Checks `observed` for `toplevel = true`. Again, for non-ss systems, there are no observables.
    # Also checks that `toplevel = true` yields subset of `get_observed`.
    @test all_sets_equal(
        observed.(
            [sys_bot, sys_bot_comp, sys_mid1, sys_mid1_comp,
                sys_mid2, sys_mid2_comp, sys_top, sys_top_comp];
            toplevel = true)...,
        [])
    @test_broken issetequal(observed(sys_bot_ss; toplevel = true), eqs_bot[3:3])
    @test_broken issetequal(observed(sys_mid1_ss; toplevel = true), eqs_mid1[3:3])
    @test_broken issetequal(observed(sys_mid2_ss; toplevel = true), eqs_mid2[3:3])
    @test_broken issetequal(observed(sys_top_ss; toplevel = true), eqs_top[3:3])
    @test all(sym_issubset(observed(sys; toplevel = true), get_observed(sys))
    for sys in [sys_bot, sys_mid2, sys_mid1, sys_top])

    # Checks `continuous_events` and  `discrete_events` for `toplevel = true` (more straightforward
    # as I stored the same singe event in all systems). Don't check for `toplevel = false` as
    # technically not needed for these tests and name spacing the events is a mess.
    mtk_cev = ModelingToolkit.SymbolicContinuousCallback.(cevs)[1]
    mtk_dev = ModelingToolkit.SymbolicDiscreteCallback.(devs)[1]
    @test all_sets_equal(
        continuous_events.(
            [sys_bot, sys_bot_comp, sys_bot_ss, sys_mid1, sys_mid1_comp, sys_mid1_ss,
                sys_mid2, sys_mid2_comp, sys_mid2_ss, sys_top, sys_top_comp, sys_top_ss];
            toplevel = true)...,
        [mtk_cev])
    @test all_sets_equal(
        discrete_events.(
            [sys_bot, sys_bot_comp, sys_bot_ss, sys_mid1, sys_mid1_comp, sys_mid1_ss,
                sys_mid2, sys_mid2_comp, sys_mid2_ss, sys_top, sys_top_comp, sys_top_ss];
            toplevel = true)...,
        [mtk_dev])
    @test all(sym_issubset(
                  continuous_events(sys; toplevel = true), get_continuous_events(sys))
    for sys in [sys_bot, sys_mid2, sys_mid1, sys_top])
    @test all(sym_issubset(discrete_events(sys; toplevel = true), get_discrete_events(sys))
    for sys in [sys_bot, sys_mid2, sys_mid1, sys_top])
end
