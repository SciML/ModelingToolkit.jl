using ModelingToolkitBase, Test, Logging
using ModelingToolkitBase: topsort_equations, t_nounits as t, D_nounits as D, unwrap,
    _override_toggle, _route_problem_verbose, _toggle_enabled
using SciMLLogging: SciMLLogging, Silent, InfoLevel, WarnLevel,
    None, Minimal, Standard, Detailed, All
using SymbolicIndexingInterface: remake_buffer
using SciMLBase: BVProblem

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
    # MTKVerbosity is consumed and never lands in prob.kwargs
    prob = ODEProblem(sys, op, tspan; verbose = MTKVerbosity())
    @test !haskey(prob.kwargs, :verbose)
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
end

@testset "cyclic_dependency toggle" begin
    @variables x(t) y(t)
    @mtkcompile sys = System([D(x) ~ -x, D(y) ~ -y], t)
    op = Dict(unwrap(x) => 2unwrap(y), unwrap(y) => 2unwrap(x))
    tspan = (0.0, 1.0)
    make(kws...) = try
        ODEProblem(sys, op, tspan; kws...)
    catch
        nothing # the cyclic ICs may error after the warning; only the logs matter
    end
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
    @test_logs min_level = Logging.Warn ODEProblem(
        sys, [], tspan; guesses, warn_initialize_determined = false
    )
    # explicit deprecated `true` overrides even a `None()` specifier
    @test_logs (:warn, r"overdetermined") match_mode = :any ODEProblem(
        sys, [], tspan; guesses, verbose = MTKVerbosity(None()),
        warn_initialize_determined = true
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

@connector function UnbalancedPin(; name)
    @variables v(t) w(t) i(t) [connect = Flow]
    return System(Equation[], t, [v, w, i], []; name)
end

@connector function QuietUnbalancedPin(; name)
    @variables v(t) w(t) i(t) [connect = Flow]
    return System(
        Equation[], t, [v, w, i], []; name,
        verbose = MTKVerbosity(unbalanced_connector = Silent)
    )
end

@testset "unbalanced_connector toggle" begin
    @test_logs (:warn, r"flow variables, yet") match_mode = :any UnbalancedPin(; name = :p)
    @test_logs min_level = Logging.Warn QuietUnbalancedPin(; name = :p)
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
