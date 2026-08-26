using ModelingToolkit, Test, Logging
using ModelingToolkit: t_nounits as t, D_nounits as D, unwrap, topsort_equations, IfLifting
using SciMLLogging: SciMLLogging, MessageLevel, Silent, InfoLevel, WarnLevel,
    None, Minimal, Standard, Detailed, All
const MTKBase = ModelingToolkit.MTKBase

@testset "MTKVerbosity construction" begin
    @test MTKVerbosity() isa MTKVerbosity{true}
    @test MTKVerbosity(None()) isa MTKVerbosity{false}
    for preset in (None(), Minimal(), Standard(), Detailed(), All())
        @test MTKVerbosity(preset) isa MTKVerbosity
    end
    # `Standard` is the no-argument default
    @test MTKVerbosity() == MTKVerbosity(Standard())
    # per-toggle kwarg overrides
    verb = MTKVerbosity(observed_equation_cycle = InfoLevel)
    @test verb.observed_equation_cycle == InfoLevel
    @test verb.state_priority_tie == WarnLevel
    # group kwarg sets every toggle in the group
    verb = MTKVerbosity(compilation = Silent)
    @test verb.state_priority_tie == Silent
    @test verb.underconstrained_variables == Silent
    @test verb.if_lifting_condition_grammar == Silent
    @test verb.observed_equation_cycle == Silent
end

@testset "_process_verbose_param" begin
    process = MTKBase._process_verbose_param
    verb = MTKVerbosity(state_priority_tie = Silent)
    @test process(verb) === verb
    @test process(Detailed()) == MTKVerbosity(Detailed())
    @test process(true) == MTKVerbosity()
    @test process(false) isa MTKVerbosity{false}
end

@testset "state_priority_tie toggle" begin
    function tied_priority_sys()
        @variables x(t) y(t)
        return System(
            [D(x) ~ -x, x ~ y], t;
            state_priorities = [x => 105, y => 105], name = :sys
        )
    end
    @test_logs (:warn, r"state_priority") match_mode = :any mtkcompile(tied_priority_sys())
    # silenced via preset, Bool, and per-toggle override
    @test_logs min_level = Logging.Warn mtkcompile(tied_priority_sys(); verbose = None())
    @test_logs min_level = Logging.Warn mtkcompile(tied_priority_sys(); verbose = false)
    @test_logs min_level = Logging.Warn mtkcompile(
        tied_priority_sys(); verbose = MTKVerbosity(state_priority_tie = Silent)
    )
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
    # silent by default, but always throws
    @test_logs @test_throws ArgumentError topsort_equations(sys, cyclic_eqs, vars)
    # the cycle is reported when the toggle is raised; the error is still thrown
    @test_logs (:info, r"cycle") match_mode = :any @test_throws ArgumentError topsort_equations(
        sys, cyclic_eqs, vars; verbose = MTKVerbosity(observed_equation_cycle = InfoLevel)
    )
    @test_logs (:info, r"cycle") match_mode = :any @test_throws ArgumentError topsort_equations(
        sys, cyclic_eqs, vars; verbose = MTKVerbosity(All())
    )
end

@testset "if_lifting_condition_grammar toggle" begin
    # A non-call constant that is not `Bool`/`Int` falls outside the limited conditional
    # grammar `CondRewriter` supports and triggers the diagnostic.
    bad_expr = MTKBase.COMMON_INF
    dep = MTKBase.COMMON_TRUE
    cw = ModelingToolkit.CondRewriter(unwrap(t))
    @test_logs (:warn, r"limited conditional grammar") match_mode = :any cw(bad_expr, dep)
    cw_silent = ModelingToolkit.CondRewriter(unwrap(t), MTKVerbosity(None()))
    @test_logs min_level = Logging.Warn cw_silent(bad_expr, dep)

    # end-to-end: IfLifting as an additional pass receives the verbosity specifier
    @variables x(t)
    @named abs_sys = System([D(x) ~ abs(x)], t)
    @test_nowarn mtkcompile(abs_sys; additional_passes = [IfLifting], verbose = None())
end

@testset "verbose threads through the full pipeline" begin
    @variables x(t) y(t)
    @parameters τ
    @named sys = System([D(x) ~ (y - x) / τ, y ~ 2x], t)
    for verbose in (
        Standard(), Detailed(), All(), true, false, MTKVerbosity(),
        MTKVerbosity(None()),
    )
        @test mtkcompile(sys; verbose) isa System
    end
end

@testset "empty_operating_point toggle" begin
    @variables x(t) = 0.0 u(t) [input = true] y(t) [output = true]
    @named sys = System([D(x) ~ -x + u, y ~ x], t)
    @test_logs (:warn, r"empty operating point") match_mode = :any linearization_function(
        sys, [u], [y]
    )
    @test_logs min_level = Logging.Warn linearization_function(
        sys, [u], [y]; verbose = SciMLLogging.None()
    )
    @test_logs min_level = Logging.Warn linearization_function(
        sys, [u], [y];
        verbose = MTKVerbosity(
            empty_operating_point = SciMLLogging.Silent,
            initialization = SciMLLogging.Silent
        )
    )
    # deprecated flags still work (the algebraic-only initialization here is
    # underdetermined, so both flags are needed for full silence)
    @test_logs min_level = Logging.Warn linearization_function(
        sys, [u], [y]; warn_empty_op = false, warn_initialize_determined = false
    )
end

@testset "initialization_analysis toggle" begin
    @variables a(t)
    @mtkcompile sys = System([D(a) ~ -a], t)
    prob = ODEProblem(sys, [a => 1.0], (0.0, 1.0))
    @test_logs (:info, r"No initialization problem|no unknowns to solve") match_mode = :any analyze_initialization_jacobian(prob)
    @test_logs min_level = Logging.Info analyze_initialization_jacobian(prob; verbose = false)
    @test_logs min_level = Logging.Info analyze_initialization_jacobian(
        prob; verbose = MTKVerbosity(initialization_analysis = SciMLLogging.Silent)
    )
end

@testset "scc_initialization_unavailable toggle" begin
    @variables x(t)
    @named sys = System(
        [D(x) ~ -x], t; initialization_eqs = [x^2 ~ 4], guesses = [x => 1.0]
    )
    csys = mtkcompile(sys; split = false)
    tspan = (0.0, 1.0)
    @test_logs (:warn, r"SCCNonlinearProblem") match_mode = :any ODEProblem(csys, [], tspan)
    @test_logs min_level = Logging.Warn ODEProblem(
        csys, [], tspan;
        verbose = MTKVerbosity(scc_initialization_unavailable = SciMLLogging.Silent)
    )
    # the message promises `use_scc = false` also disables it
    @test_logs min_level = Logging.Warn ODEProblem(csys, [], tspan; use_scc = false)
end
