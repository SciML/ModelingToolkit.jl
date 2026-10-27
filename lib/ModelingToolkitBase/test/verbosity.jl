using ModelingToolkitBase, Test, Logging
using ModelingToolkitBase: topsort_equations, t_nounits as t, D_nounits as D, unwrap,
    _override_toggle, _route_problem_verbose, _toggle_enabled
using SciMLLogging: SciMLLogging, Silent, InfoLevel, WarnLevel,
    None, Minimal, Standard, Detailed, All
using SymbolicIndexingInterface: remake_buffer
using SciMLBase: SciMLBase, BVProblem, remake
using DiffEqBase: DEVerbosity
using JumpProcesses: Direct
using NonlinearSolve: NonlinearVerbosity
using SteadyStateDiffEq: DynamicSS
using OrdinaryDiffEqRosenbrock: Rodas5P

# The warnings of `f()` other than deprecation warnings.
function nondeprecation_warnings(f)
    logs, _ = Test.collect_test_logs(f; min_level = Logging.Warn)
    return filter(log -> log.group !== :depwarn, logs)
end

@testset "MTKVerbosity construction" begin
    @test MTKVerbosity() isa MTKVerbosity{true}
    @test MTKVerbosity(None()) isa MTKVerbosity{false}
    for preset in (None(), Minimal(), Standard(), Detailed(), All())
        @test MTKVerbosity(preset) isa MTKVerbosity
    end
    @test MTKVerbosity() == MTKVerbosity(Standard())
    verb = MTKVerbosity(observed_equation_cycle = InfoLevel)
    @test verb.observed_equation_cycle == InfoLevel
    @test verb.state_priority_tie == WarnLevel
    verb = MTKVerbosity(compilation = Silent)
    @test verb.if_lifting_condition_grammar == Silent
    # groups set all their toggles
    verb = MTKVerbosity(initialization = Silent)
    @test verb.singular_initialization == Silent
    @test verb.overdetermined_initialization == Silent
    @test verb.underdetermined_initialization == Silent
    @test verb.scc_initialization_unavailable == Silent
    @test verb.state_priority_tie == WarnLevel
    verb = MTKVerbosity(problem_construction = WarnLevel)
    @test verb.cyclic_dependency == WarnLevel
    verb = MTKVerbosity(analysis = Silent)
    @test verb.empty_operating_point == Silent
    @test verb.initialization_analysis == Silent
    # sub-specifier defaults and overrides
    @test MTKVerbosity().initialization_verbosity isa Minimal
    @test MTKVerbosity(None()).initialization_verbosity isa None
    @test MTKVerbosity(initialization_verbosity = None()).initialization_verbosity isa None
end

@testset "_override_toggle" begin
    verb = MTKVerbosity()
    @test _override_toggle(verb, nothing) === verb
    silenced = _override_toggle(verb, false, :singular_initialization => WarnLevel)
    @test silenced.singular_initialization == Silent
    @test silenced.overdetermined_initialization == WarnLevel
    # `true` re-enables inside a `{false}` specifier and only the named toggle
    enabled = _override_toggle(MTKVerbosity(None()), true, :cyclic_dependency => WarnLevel)
    @test enabled isa MTKVerbosity{true}
    @test enabled.cyclic_dependency == WarnLevel
    @test enabled.singular_initialization == Silent
end

@testset "observed_equation_cycle toggle" begin
    @variables x(t) y(t) z(t) k(t)
    @named sys = System(Equation[], t, [x, y, z, k], [])
    cyclic_eqs = [
        x ~ y + z
        z ~ 2
        y ~ 2z + x
    ]
    vars = unwrap.([x, y, z, k])
    @test_logs @test_throws ArgumentError topsort_equations(sys, cyclic_eqs, vars)
    @test_logs (:info, r"cycle") match_mode = :any @test_throws ArgumentError topsort_equations(
        sys, cyclic_eqs, vars; verbose = MTKVerbosity(observed_equation_cycle = InfoLevel)
    )
end

@testset "mtkcompile accepts verbose" begin
    @variables x y
    @named sys = System([0 ~ x^2 + y, 0 ~ y - 2x], [x, y], [])
    for verbose in (Standard(), true, false, MTKVerbosity(), MTKVerbosity(None()))
        @test mtkcompile(sys; verbose) isa System
    end
end

@testset "problem constructor verbose type-routing" begin
    @variables x(t) y(t)
    @parameters τ = 1.0
    @mtkcompile sys = System([D(x) ~ (y - x) / τ, y ~ 2x], t)
    op = [x => 1.0]
    tspan = (0.0, 1.0)
    # MTKVerbosity forwards its solver sub-specifier for the problem type
    prob = ODEProblem(sys, op, tspan; verbose = MTKVerbosity())
    @test prob.kwargs[:verbose] isa Standard
    prob = ODEProblem(sys, op, tspan; verbose = MTKVerbosity(de_verbosity = None()))
    @test prob.kwargs[:verbose] isa None
    # Bool and preset are forwarded to prob.kwargs (solver-facing) as before
    prob = ODEProblem(sys, op, tspan; verbose = false)
    @test prob.kwargs[:verbose] === false
    prob = ODEProblem(sys, op, tspan; verbose = Minimal())
    @test prob.kwargs[:verbose] isa Minimal
    # no verbose kwarg: nothing stored
    prob = ODEProblem(sys, op, tspan)
    @test !haskey(prob.kwargs, :verbose)
end

@testset "initialization_verbosity sub-specifier forwarding" begin
    @variables a(t) b(t)
    @mtkcompile sys = System([D(a) ~ -a + b, b ~ 2a], t)
    tspan = (0.0, 1.0)
    getinit(prob) = prob.f.initialization_data.initializeprob.kwargs[:verbose]
    prob = ODEProblem(sys, [a => 1.0], tspan)
    @test getinit(prob) isa Minimal
    prob = ODEProblem(sys, [a => 1.0], tspan; verbose = false)
    @test getinit(prob) isa None
    prob = ODEProblem(
        sys, [a => 1.0], tspan; verbose = MTKVerbosity(initialization_verbosity = None())
    )
    @test getinit(prob) isa None

    # the inner problem keeps the outer toggles, and gets `initialization_verbosity` as
    # its solver verbosity
    initverb = ModelingToolkitBase._initialization_problem_verbosity
    inner = initverb(
        MTKVerbosity(problem_construction = Silent, initialization_verbosity = None())
    )
    @test inner.parameter_not_in_system == Silent
    @test inner.nonlinear_verbosity isa None
    @test inner.linear_verbosity isa None
    nlverb = NonlinearVerbosity()
    inner = initverb(MTKVerbosity(initialization_verbosity = nlverb))
    @test inner.parameter_not_in_system == MTKVerbosity().parameter_not_in_system
    @test inner.nonlinear_verbosity === nlverb
    @test inner.linear_verbosity === MTKVerbosity().linear_verbosity

    # a linear initialization problem does not get a `NonlinearVerbosity`
    prob = ODEProblem(
        sys, [b => 2.0], tspan; guesses = [a => 1.0],
        verbose = MTKVerbosity(initialization_verbosity = nlverb)
    )
    # `LinearInitializationProblem` wraps the `LinearProblem` in an `SCCNonlinearProblem`
    linprob = only(prob.f.initialization_data.initializeprob.probs)
    @test linprob isa SciMLBase.LinearProblem
    @test linprob.kwargs[:verbose] === MTKVerbosity().linear_verbosity
end

@testset "cyclic_dependency toggle" begin
    @variables x(t) y(t)
    @mtkcompile sys = System([D(x) ~ -x, D(y) ~ -y], t)
    op = Dict(unwrap(x) => 2unwrap(y), unwrap(y) => 2unwrap(x))
    tspan = (0.0, 1.0)
    make(kws...) = ODEProblem(sys, op, tspan; kws...)
    @test_logs min_level = Logging.Warn make()
    @test_logs (:warn, r"Cycles in unknowns") match_mode = :any make(
        :verbose => MTKVerbosity(Detailed())
    )
    @test_logs (:warn, r"Cycles in unknowns") match_mode = :any make(
        :verbose => MTKVerbosity(cyclic_dependency = WarnLevel)
    )
    @test_logs (:warn, r"Cycles in unknowns") match_mode = :any make(
        :warn_cyclic_dependency => true
    )
    @test_deprecated r"warn_cyclic_dependency" make(:warn_cyclic_dependency => true)
    # with the toggle on, MTK passes `warn_cyclic_dependency` to the inner initialization
    # problem itself, which is not a deprecated use
    logs, _ = Test.collect_test_logs(min_level = Logging.Warn) do
        ODEProblem(sys, [x => 1.0, y => 1.0], tspan; verbose = MTKVerbosity(Detailed()))
    end
    @test !any(log -> log.group === :depwarn, logs)
end

@testset "determinacy toggles and deprecated override" begin
    @variables x(t)
    @mtkcompile sys = System([D(x) ~ -x], t; initialization_eqs = [x ~ 1, x^2 ~ 1])
    tspan = (0.0, 1.0)
    guesses = [x => 1.0]
    @test_logs (:warn, r"overdetermined") match_mode = :any ODEProblem(
        sys, [], tspan; guesses
    )
    @test_logs min_level = Logging.Warn ODEProblem(
        sys, [], tspan; guesses, verbose = MTKVerbosity(overdetermined_initialization = Silent)
    )
    @test isempty(
        nondeprecation_warnings() do
            ODEProblem(sys, [], tspan; guesses, warn_initialize_determined = false)
        end
    )
    @test_deprecated r"warn_initialize_determined" ODEProblem(
        sys, [], tspan; guesses, warn_initialize_determined = false
    )
    # explicit deprecated `true` overrides even a `None()` specifier
    @test_logs (:warn, r"overdetermined") match_mode = :any ODEProblem(
        sys, [], tspan; guesses, verbose = MTKVerbosity(None()),
        warn_initialize_determined = true
    )
end

@testset "_problem_verbosity applies the deprecated keywords" begin
    problem_verbosity = ModelingToolkitBase._problem_verbosity
    @test problem_verbosity((;)) == MTKVerbosity()
    @test problem_verbosity((; verbose = false)) == MTKVerbosity(None())
    verb = problem_verbosity(
        (;
            verbose = MTKVerbosity(None()), warn_initialize_determined = true,
            warn_cyclic_dependency = true,
        )
    )
    @test verb.overdetermined_initialization == WarnLevel
    @test verb.singular_initialization == WarnLevel
    @test verb.cyclic_dependency == WarnLevel
    @test verb.parameter_not_in_system == Silent
end

@testset "get_initialization_problem_type deprecated keyword" begin
    @variables x(t)
    @mtkcompile sys = System([D(x) ~ -x], t; initialization_eqs = [x ~ 1, x^2 ~ 1])
    prob = ODEProblem(sys, [], (0.0, 1.0); guesses = [x => 1.0], verbose = false)
    isys = prob.f.initialization_data.initializeprob.f.sys
    probtype = ModelingToolkitBase.get_initialization_problem_type
    @test_logs (:warn, r"overdetermined") match_mode = :any probtype(sys, isys)
    @test_logs min_level = Logging.Warn probtype(sys, isys; verbose = false)
    @test_deprecated r"warn_initialize_determined" probtype(
        sys, isys; warn_initialize_determined = false
    )
    # the deprecated keyword still overrides the toggles
    @test_logs (:warn, r"overdetermined") match_mode = :any probtype(
        sys, isys; verbose = false, warn_initialize_determined = true
    )
end

@testset "no_unbound_inputs toggle" begin
    @variables x(t)
    @named sys = System([D(x) ~ -x], t)
    @test_logs (:warn, r"No unbound inputs") match_mode = :any generate_control_function(sys)
    @test_logs min_level = Logging.Warn generate_control_function(sys; verbose = None())
    @test_logs min_level = Logging.Warn generate_control_function(
        sys; verbose = MTKVerbosity(no_unbound_inputs = Silent)
    )
end

@testset "analysis_point_causality toggle" begin
    @variables vin(t) [input = true] vout(t) [output = true]
    # reversed causality: the first argument should be an output
    @test_logs (:warn, r"was not a output|was not an output") match_mode = :any connect(
        vin, :ap, vout
    )
    @test_logs min_level = Logging.Warn connect(vin, :ap, vout; verbose = false)
    @test_logs min_level = Logging.Warn connect(
        vin, :ap, vout; verbose = MTKVerbosity(analysis_point_causality = Silent)
    )
end

@testset "system_construction group" begin
    verb = MTKVerbosity(system_construction = Silent)
    @test verb.independent_variable_not_parameter == Silent
    @test verb.constraint_variable_without_argument == Silent
    @test verb.unbalanced_connector == Silent
    @test verb.affect_default_independent_variable == Silent
    @test verb.substitute_skips_events == Silent
    @test verb.imperative_affect_specification == WarnLevel
end

@testset "System constructor stores verbose" begin
    @variables τ
    Dτ = Differential(τ)
    @variables u(τ)
    @test_logs (:warn, r"@independent_variables") match_mode = :any System(
        [Dτ(u) ~ -u], τ; name = :s
    )
    quiet = @test_logs min_level = Logging.Warn System(
        [Dτ(u) ~ -u], τ; name = :s, verbose = false
    )
    @test ModelingToolkitBase._system_verbosity(quiet) == MTKVerbosity(None())
    # a system rebuilt from the metadata keeps the stored value
    @test_logs min_level = Logging.Warn System(
        [Dτ(u) ~ -u], τ; name = :s, metadata = ModelingToolkitBase.get_metadata(quiet)
    )
    @variables y(t)
    @named loud = System([D(y) ~ -y], t)
    @test ModelingToolkitBase._system_verbosity(loud) == MTKVerbosity()
end

@testset "mtkcompile keeps the verbosity stored on the system" begin
    @variables y(t)
    @named quiet = System([D(y) ~ -y], t; verbose = false)
    @test ModelingToolkitBase._system_verbosity(mtkcompile(quiet)) == MTKVerbosity(None())
    # an explicit `verbose` wins
    csys = mtkcompile(quiet; verbose = true)
    @test ModelingToolkitBase._system_verbosity(csys) == MTKVerbosity()
    @named loud = System([D(y) ~ -y], t)
    @test ModelingToolkitBase._system_verbosity(mtkcompile(loud)) == MTKVerbosity()
    # the compiled system stores the verbosity once
    n_stored(sys) = count(
        entry -> first(entry) === ModelingToolkitBase.MTKVerbosityCtx,
        ModelingToolkitBase.get_metadata(sys)
    )
    @test n_stored(mtkcompile(quiet)) == 1
    @test n_stored(mtkcompile(loud)) == 1
end

@testset "SciMLLogging is exported for the MTKVerbosity levels" begin
    @test Base.isexported(ModelingToolkitBase, :SciMLLogging)
    @test ModelingToolkitBase.SciMLLogging.Silent isa SciMLLogging.MessageLevel
end

@testset "constraint_variable_without_argument toggle" begin
    @variables xf(..)
    z = unwrap(xf())
    @test_logs (:warn, r"has no argument") match_mode = :any ModelingToolkitBase.validate_vars_and_find_ps!(
        Set([z]), Set(), Set([z]), unwrap(t)
    )
    @test_logs min_level = Logging.Warn ModelingToolkitBase.validate_vars_and_find_ps!(
        Set([z]), Set(), Set([z]), unwrap(t), MTKVerbosity(None())
    )
end

# `balanced = true` gives a balanced instance of the same connector type
@connector function UnbalancedPin(; name, balanced = false)
    @variables v(t) w(t) i(t) [connect = Flow]
    return System(Equation[], t, balanced ? [v, i] : [v, w, i], []; name)
end

@connector function QuietUnbalancedPin(; name)
    @variables v(t) w(t) i(t) [connect = Flow]
    return System(
        Equation[], t, [v, w, i], []; name,
        verbose = MTKVerbosity(unbalanced_connector = Silent)
    )
end

@testset "unbalanced_connector toggle" begin
    # creating a connector does not warn
    @test_logs min_level = Logging.Warn UnbalancedPin(; name = :p)
    @variables x(t)
    # a balanced model: the flow variable of an unconnected pin is zero, and the
    # other variables of a pin get an equation here
    function n_unbalanced_warnings(pins...; kws...)
        eqs = [D(x) ~ -x]
        for pin in pins, u in unknowns(pin)
            ModelingToolkitBase.get_connection_type(u) === ModelingToolkitBase.Flow && continue
            push!(eqs, getproperty(pin, ModelingToolkitBase.getname(u)) ~ 1.0)
        end
        model = System(eqs, t; name = :model, systems = collect(pins))
        logs, _ = Test.collect_test_logs(min_level = Logging.Warn) do
            mtkcompile(model; kws...)
        end
        return count(log -> occursin("flow variables, yet", string(log.message)), logs)
    end
    n_unbalanced_warnings(pin::Function; kws...) =
        n_unbalanced_warnings(pin(; name = :p), pin(; name = :q); kws...)
    # `mtkcompile` warns once for each connector type
    @test n_unbalanced_warnings(UnbalancedPin) == 1
    @test n_unbalanced_warnings(UnbalancedPin; verbose = false) == 0
    @test n_unbalanced_warnings(
        UnbalancedPin; verbose = MTKVerbosity(unbalanced_connector = Silent)
    ) == 0
    # the author of a connector can silence it
    @test n_unbalanced_warnings(QuietUnbalancedPin) == 0
    # a balanced instance does not hide an unbalanced instance of the same type
    bad = UnbalancedPin(; name = :bad)
    good = UnbalancedPin(; name = :good, balanced = true)
    @test n_unbalanced_warnings(bad, good) == 1
    @test n_unbalanced_warnings(good, bad) == 1
end

@testset "affect_default_independent_variable toggle" begin
    @variables x(t)
    @parameters p
    @named loud = System([D(x) ~ -x], t)
    @named quiet = System([D(x) ~ -x], t; verbose = false)
    aff = ModelingToolkitBase.SymbolicAffect([x ~ p])
    @test_logs (:warn, r"Defaulting to t_nounits") match_mode = :any ModelingToolkitBase.AffectSystem(
        aff; parent_sys = loud
    )
    @test_logs min_level = Logging.Warn ModelingToolkitBase.AffectSystem(
        aff; parent_sys = quiet
    )
end

@testset "substitute_skips_events toggle" begin
    @variables x(t)
    @parameters p q
    @named sys = System(
        [D(x) ~ -p * x], t; continuous_events = [[x ~ 0.5] => [x ~ Pre(x) + 1]]
    )
    @test_logs (:warn, r"has events") match_mode = :any substitute(sys, [p => q])
    @test_logs min_level = Logging.Warn substitute(sys, [p => q]; verbose = false)
end

@testset "imperative_affect_specification toggle" begin
    @variables temp(t)
    @parameters furnace_on::Bool = false
    cb = ModelingToolkitBase.SymbolicContinuousCallback(
        [temp ~ 0.7],
        ModelingToolkitBase.ImperativeAffect(
            modified = (; furnace_on), observed = (; furnace_on)
        ) do x, o, c, i
            return (; furnace_on = false)
        end
    )
    @named sys = System([D(temp) ~ furnace_on - temp], t; continuous_events = [cb])
    sys = mtkcompile(sys)
    op = [temp => 0.0, furnace_on => true]
    @test_logs (:warn, r"both observed and modified") match_mode = :any ODEProblem(
        sys, op, (0.0, 1.0)
    )
    @test_logs min_level = Logging.Warn ODEProblem(
        sys, op, (0.0, 1.0); verbose = MTKVerbosity(imperative_affect_specification = Silent)
    )
    # a preset silences MTK and is still forwarded to the solver
    prob = @test_logs min_level = Logging.Warn ODEProblem(sys, op, (0.0, 1.0); verbose = None())
    @test prob.kwargs[:verbose] isa None
end

@testset "parameter_not_in_system toggle" begin
    @variables x(t)
    @parameters p q
    @named sys = System([D(x) ~ -p * x], t)
    loud = mtkcompile(sys)
    quiet = mtkcompile(sys; verbose = false)
    @test ModelingToolkitBase._system_verbosity(quiet) == MTKVerbosity(None())
    for (csys, logs) in ((loud, true), (quiet, false))
        prob = ODEProblem(csys, [x => 1.0, p => 1.0], (0.0, 1.0))
        if logs
            @test_logs (:warn, r"not a \(non-dependent\) parameter") match_mode = :any remake_buffer(
                csys, prob.p, [q], [2.0]
            )
        else
            @test_logs min_level = Logging.Warn remake_buffer(csys, prob.p, [q], [2.0])
        end
    end
end

@testset "initialization system stores the problem verbosity" begin
    @variables x(t) y(t)
    @named sys = System([D(x) ~ -x, y ~ 2x], t; guesses = [x => 1.0])
    sys = mtkcompile(sys)
    initsys(prob) = prob.f.initialization_data.initializeprob.f.sys
    verb = MTKVerbosity(None())
    prob = ODEProblem(sys, [y => 4.0], (0.0, 1.0); verbose = verb)
    @test ModelingToolkitBase._system_verbosity(initsys(prob)) == verb
    prob = ODEProblem(sys, [y => 4.0], (0.0, 1.0))
    @test ModelingToolkitBase._system_verbosity(initsys(prob)) == MTKVerbosity()
end

@testset "remake rebuilds initialization with the system verbosity" begin
    @variables x(t) y(t)
    # a nonlinear initialization without SCCs gives a plain `NonlinearProblem`
    @named sys = System([D(x) ~ -x, y ~ x^3 + x], t; guesses = [x => 1.0])
    sys = mtkcompile(sys; verbose = false)
    prob = ODEProblem(sys, [y => 2.0], (0.0, 1.0); verbose = false, use_scc = false)
    # a symbolic `u0` makes `remake` rebuild the initialization problem
    iprob = remake(prob; u0 = [y => 10.0]).f.initialization_data.initializeprob
    @test iprob isa NonlinearProblem
    @test iprob.kwargs[:verbose] == None()
    @test ModelingToolkitBase._system_verbosity(iprob.f.sys) == MTKVerbosity(None())
end

@testset "overdetermined_constraints toggle" begin
    @variables x(..)
    @parameters k = 1.0
    @mtkcompile sys = System([D(x(t)) ~ -k * x(t)], t; constraints = [x(0.5) ~ 0.5])
    # 1 constraint + 1 fixed initial value > 1 unknown
    op = [x(t) => 1.0]
    @test_logs (:warn, r"BVProblem is overdetermined") match_mode = :any BVProblem(
        sys, op, (0.0, 1.0)
    )
    @test_logs min_level = Logging.Warn BVProblem(
        sys, op, (0.0, 1.0); verbose = MTKVerbosity(overdetermined_constraints = Silent)
    )
end

@testset "dynamic_opt_time_grid toggle" begin
    process_tspan = ModelingToolkitBase.process_tspan
    quiet = MTKVerbosity(dynamic_opt_time_grid = Silent)
    # fixed time span: `steps` is ignored
    @test_logs (:warn, r"number of steps") match_mode = :any process_tspan(
        (0.0, 1.0), 0.1, 10
    )
    @test_logs min_level = Logging.Warn process_tspan((0.0, 1.0), 0.1, 10, quiet)
    # free final time: `dt` is ignored
    @parameters tf
    @test_logs (:warn, r"Specified dt") match_mode = :any process_tspan((0.0, tf), 0.1, 10)
    @test_logs min_level = Logging.Warn process_tspan((0.0, tf), 0.1, 10, quiet)
end

@testset "solver sub-specifiers" begin
    field = ModelingToolkitBase._solver_verbosity_field
    @test field(SciMLBase.ODEProblem) === :de_verbosity
    @test field(SciMLBase.SDEProblem) === :de_verbosity
    @test field(SciMLBase.DiscreteProblem) === :de_verbosity
    @test field(JumpProblem) === :de_verbosity
    @test field(BVProblem) === :bvp_verbosity
    @test field(SciMLBase.NonlinearProblem) === :nonlinear_verbosity
    @test field(SciMLBase.SteadyStateProblem) === :nonlinear_verbosity
    @test field(SciMLBase.IntervalNonlinearProblem) === :nonlinear_verbosity
    @test field(SciMLBase.LinearProblem) === :linear_verbosity
    @test field(SciMLBase.OptimizationProblem) === :optimization_verbosity
    @test field(Int) === nothing
    # the dynamic optimization `solve` methods take their own `verbose::Bool`
    @test field(SciMLBase.AbstractDynamicOptProblem) === nothing
    @test ModelingToolkitBase._solver_kwargs(
        SciMLBase.AbstractDynamicOptProblem, (; verbose = MTKVerbosity(), a = 1)
    ) == (; a = 1)
    # presets set every solver sub-specifier to the same preset
    @test MTKVerbosity().de_verbosity isa Standard
    @test MTKVerbosity(None()).de_verbosity isa None
    @test MTKVerbosity(Detailed()).nonlinear_verbosity isa Detailed

    @variables x(t)
    @parameters k = 1.0
    @mtkcompile sys = System([D(x) ~ -k * x], t)

    # a `SteadyStateProblem` gets `nonlinear_verbosity`; `DynamicSS` does not give it to
    # its ODE solve
    @mtkcompile sssys = System([D(x) ~ k * (1 - x)], t)
    ssprob = SteadyStateProblem(
        sssys, [x => 0.5];
        verbose = MTKVerbosity(nonlinear_verbosity = NonlinearVerbosity())
    )
    @test ssprob.kwargs[:verbose] isa NonlinearVerbosity
    @test SciMLBase.successful_retcode(solve(ssprob, DynamicSS(Rodas5P())))

    # a nested solver specifier gets to the solver unchanged
    de = DEVerbosity(linear_verbosity = None())
    prob = ODEProblem(sys, [x => 1.0], (0.0, 1.0); verbose = MTKVerbosity(de_verbosity = de))
    @test prob.kwargs[:verbose] === de

    @variables a b
    @parameters c = 4.0
    @mtkcompile nsys = System([0 ~ a^2 - c, 0 ~ b - a])
    prob = NonlinearProblem(
        nsys, [a => 1.0, b => 1.0]; verbose = MTKVerbosity(nonlinear_verbosity = Detailed())
    )
    @test prob.kwargs[:verbose] isa Detailed

    @variables y(..)
    @mtkcompile bsys = System([D(y(t)) ~ -k * y(t)], t; constraints = [y(0.5) ~ 0.5])
    prob = BVProblem(
        bsys, [y(t) => 1.0], (0.0, 1.0); verbose = MTKVerbosity(bvp_verbosity = None())
    )
    @test prob.kwargs[:verbose] isa None

    @parameters β = 0.1
    @variables S(t) = 10
    jsys = complete(JumpSystem([SymbolicMassActionJump(β, [S => 1], [S => -1])], t, [S], [β]; name = :jsys))
    jprob = JumpProblem(
        jsys, [S => 10], (0.0, 1.0); aggregator = Direct(),
        verbose = MTKVerbosity(de_verbosity = None())
    )
    # `solve` gets the keywords of the inner problem
    @test jprob.prob.kwargs[:verbose] isa None
end
