using ModelingToolkit, Symbolics, Test
using ModelingToolkit: get_connector_type, get_defaults, get_gui_metadata,
                       get_systems, get_ps, getdefault, getname, readable_code,
                       scalarize, symtype, VariableDescription, RegularConnector,
                       get_unit
using URIs: URI
using Distributions
using DynamicQuantities, OrdinaryDiffEq
using ModelingToolkit: t, D

ENV["MTK_ICONS_DIR"] = "$(@__DIR__)/icons"

# Mock module used to test if the `@mtkmodel` macro works with fully-qualified names as well.
module MyMockModule
using ModelingToolkit, DynamicQuantities
using ModelingToolkit: t, D

export Pin
@connector Pin begin
    v(t), [unit = u"V"]                    # Potential at the pin [V]
    i(t), [connect = Flow, unit = u"A"]    # Current flowing into the pin [A]
    @icon "pin.png"
end

ground_logo = read(abspath(ENV["MTK_ICONS_DIR"], "ground.svg"), String)
@mtkmodel Ground begin
    @components begin
        g = Pin()
    end
    @icon ground_logo
    @equations begin
        g.v ~ 0
    end
end
end

using .MyMockModule

@connector RealInput begin
    u(t), [input = true, unit = u"V"]
end
@connector RealOutput begin
    u(t), [output = true, unit = u"V"]
end
@mtkmodel Constant begin
    @components begin
        output = RealOutput()
    end
    @parameters begin
        k, [description = "Constant output value of block", unit = u"V"]
    end
    @equations begin
        output.u ~ k
    end
end

@named p = Pin(; v = π * u"V")

@test getdefault(p.v) ≈ π
@test Pin.isconnector == true

@mtkmodel OnePort begin
    @components begin
        p = Pin()
        n = Pin()
    end
    @variables begin
        v(t), [unit = u"V"]
        i(t), [unit = u"A"]
    end
    @icon "oneport.png"
    @equations begin
        v ~ p.v - n.v
        0 ~ p.i + n.i
        i ~ p.i
    end
end

@test OnePort.isconnector == false

@mtkmodel Resistor begin
    @extend v, i = oneport = OnePort()
    @parameters begin
        R, [unit = u"Ω"]
    end
    @icon """<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="80" height="30">
<path d="M10 15
l15 0
l2.5 -5
l5 10
l5 -10
l5 10
l5 -10
l5 10
l2.5 -5
l15 0" stroke="black" stroke-width="1" stroke-linejoin="bevel" fill="none"></path>
</svg>
"""
    @equations begin
        v ~ i * R
    end
end

@mtkmodel Capacitor begin
    @parameters begin
        C, [unit = u"F"]
    end
    @extend OnePort(; v = 0.0u"V")
    @icon "https://upload.wikimedia.org/wikipedia/commons/7/78/Capacitor_symbol.svg"
    @equations begin
        D(v) ~ i / C
    end
end

@named capacitor = Capacitor(C = 10u"F", v = 10.0u"V")
@test getdefault(capacitor.v) == 10.0

@mtkmodel Voltage begin
    @extend v, i = oneport = OnePort()
    @components begin
        V = RealInput()
    end
    @equations begin
        v ~ V.u
    end
end

@mtkmodel RC begin
    @description "An RC circuit."
    @structural_parameters begin
        R_val = 10u"Ω"
        C_val = 10u"F"
        k_val = 10u"V"
    end
    @components begin
        resistor = Resistor(; R = R_val)
        capacitor = Capacitor(; C = C_val)
        source = Voltage()
        constant = Constant(; k = k_val)
        ground = MyMockModule.Ground()
    end
    @equations begin
        connect(constant.output, source.V)
        connect(source.p, resistor.p)
        connect(resistor.n, capacitor.p)
        connect(capacitor.n, source.n, ground.g)
    end
end

C_val = 20u"F"
R_val = 20u"Ω"
res__R = 100u"Ω"
@mtkbuild rc = RC(; C_val, R_val, resistor.R = res__R)
prob = ODEProblem(rc, [], (0, 1e9))
sol = solve(prob)
defs = ModelingToolkit.defaults(rc)
@test sol[rc.capacitor.v, end] ≈ defs[rc.constant.k]
resistor = getproperty(rc, :resistor; namespace = false)
@test ModelingToolkit.description(rc) == "An RC circuit."
@test getname(rc.resistor) === getname(resistor)
@test getname(rc.resistor.R) === getname(resistor.R)
@test getname(rc.resistor.v) === getname(resistor.v)
# Test that `resistor.R` overrides `R_val` in the argument.
@test getdefault(rc.resistor.R) * get_unit(rc.resistor.R) == res__R != R_val
# Test that `C_val` passed via argument is set as default of C.
@test getdefault(rc.capacitor.C) * get_unit(rc.capacitor.C) == C_val
# Test that `k`'s default value is unchanged.
@test getdefault(rc.constant.k) * get_unit(rc.constant.k) ==
      eval(RC.structure[:kwargs][:k_val][:value])
@test getdefault(rc.capacitor.v) == 0.0

@test get_gui_metadata(rc.resistor).layout == Resistor.structure[:icon] ==
      read(joinpath(ENV["MTK_ICONS_DIR"], "resistor.svg"), String)
@test get_gui_metadata(rc.ground).layout ==
      read(abspath(ENV["MTK_ICONS_DIR"], "ground.svg"), String)
@test get_gui_metadata(rc.capacitor).layout ==
      URI("https://upload.wikimedia.org/wikipedia/commons/7/78/Capacitor_symbol.svg")
@test OnePort.structure[:icon] ==
      URI("file:///" * abspath(ENV["MTK_ICONS_DIR"], "oneport.png"))
@test ModelingToolkit.get_gui_metadata(rc.resistor.p).layout == Pin.structure[:icon] ==
      URI("file:///" * abspath(ENV["MTK_ICONS_DIR"], "pin.png"))

@test length(equations(rc)) == 1

@testset "Constants" begin
    @mtkmodel PiModel begin
        @constants begin
            _p::Irrational = π, [description = "Value of Pi.", unit = u"V"]
        end
        @parameters begin
            p = _p, [description = "Assign constant `_p` value."]
            e, [unit = u"V"]
        end
        @equations begin
            # This validates units; indirectly verifies that metadata was correctly passed.
            e ~ _p
        end
    end

    @named pi_model = PiModel()

    @test typeof(ModelingToolkit.getdefault(pi_model.p)) <:
          SymbolicUtils.BasicSymbolic{Irrational}
    @test getdefault(getdefault(pi_model.p)) == π
end

@testset "Parameters and Structural parameters in various modes" begin
    @mtkmodel MockModel begin
        @parameters begin
            a
            a2[1:2]
            b(t)
            b2(t)[1:2]
            cval
            jval
            kval
            c(t) = cval + jval
            d = 2
            d2[1:2] = 2
            e, [description = "e"]
            e2[1:2], [description = "e2"]
            f = 3, [description = "f"]
            h(t), [description = "h(t)"]
            h2(t)[1:2], [description = "h2(t)"]
            i(t) = 4, [description = "i(t)"]
            j(t) = jval, [description = "j(t)"]
            k = kval, [description = "k"]
            l(t)[1:2, 1:3] = 2, [description = "l is more than 1D"]
            n # test defaults with Number input
            n2 # test defaults with Function input
        end
        @structural_parameters begin
            m = 1
            func
        end
        begin
            g() = 5
        end
        @defaults begin
            n => 1.0
            n2 => g()
        end
    end

    kval = 5
    @named model = MockModel(; b2 = [1, 3], kval, cval = 1, func = identity)

    @test lastindex(parameters(model)) == 31

    @test all(getdescription.([model.e2...]) .== "e2")
    @test all(getdescription.([model.h2...]) .== "h2(t)")

    @test hasmetadata(model.e, VariableDescription)
    @test hasmetadata(model.f, VariableDescription)
    @test hasmetadata(model.h, VariableDescription)
    @test hasmetadata(model.i, VariableDescription)
    @test hasmetadata(model.j, VariableDescription)
    @test hasmetadata(model.k, VariableDescription)
    @test all(collect(hasmetadata.(model.l, ModelingToolkit.VariableDescription)))

    @test all(lastindex.([model.a2, model.b2, model.d2, model.e2, model.h2]) .== 2)
    @test size(model.l) == (2, 3)
    @test MockModel.structure[:parameters][:l][:size] == (2, 3)

    model = complete(model)
    @test getdefault(model.cval) == 1
    @test isequal(getdefault(model.c), model.cval + model.jval)
    @test getdefault(model.d) == 2
    @test_throws ErrorException getdefault(model.e)
    @test getdefault(model.f) == 3
    @test getdefault(model.i) == 4
    @test all(getdefault.(scalarize(model.b2)) .== [1, 3])
    @test all(getdefault.(scalarize(model.l)) .== 2)
    @test isequal(getdefault(model.j), model.jval)
    @test isequal(getdefault(model.k), model.kval)
    @test get_defaults(model)[model.n] == 1.0
    @test get_defaults(model)[model.n2] == 5

    @test MockModel.structure[:defaults] == Dict(:n => 1.0, :n2 => "g()")
end

@testset "Arrays using vanilla-@variable syntax" begin
    @mtkmodel TupleInArrayDef begin
        @structural_parameters begin
            N
            M
        end
        @parameters begin
            (l(t)[1:2, 1:3] = 1), [description = "l is more than 1D"]
            (l2(t)[1:N, 1:M] = 2),
            [description = "l is more than 1D, with arbitrary length"]
            (l3(t)[1:3] = 3), [description = "l2 is 1D"]
            (l4(t)[1:N] = 4), [description = "l2 is 1D, with arbitrary length"]
            (l5(t)[1:3]::Int = 5), [description = "l3 is 1D and has a type"]
            (l6(t)[1:N]::Int = 6),
            [description = "l3 is 1D and has a type, with arbitrary length"]
        end
    end

    N, M = 4, 5
    @named arr = TupleInArrayDef(; N, M)
    @test getdefault(arr.l) == 1
    @test getdefault(arr.l2) == 2
    @test getdefault(arr.l3) == 3
    @test getdefault(arr.l4) == 4
    @test getdefault(arr.l5) == 5
    @test getdefault(arr.l6) == 6

    @test size(arr.l2) == (N, M)
    @test size(arr.l4) == (N,)
    @test size(arr.l6) == (N,)
end

@testset "Type annotation" begin
    @mtkmodel TypeModel begin
        @structural_parameters begin
            flag::Bool = true
        end
        @parameters begin
            par0::Bool = true
            par1::Int = 1
            par2(t)::Int,
            [description = "Enforced `par4` to be an Int by setting the type to the keyword-arg."]
            par3(t)::BigFloat = 1.0
            par4(t)::Float64 = 1 # converts 1 to 1.0 of Float64 type
            par5[1:3]::BigFloat
            par6(t)[1:3]::BigFloat
            par7(t)[1:3, 1:3]::BigFloat = 1.0, [description = "with description"]
        end
    end

    @named type_model = TypeModel()

    @test symtype(type_model.par1) == Int
    @test symtype(type_model.par2) == Int
    @test symtype(type_model.par3) == BigFloat
    @test symtype(type_model.par4) == Float64
    @test symtype(type_model.par5[1]) == BigFloat
    @test symtype(type_model.par6[1]) == BigFloat
    @test symtype(type_model.par7[1, 1]) == BigFloat

    @test_throws TypeError TypeModel(; name = :throws, flag = 1)
    @test_throws TypeError TypeModel(; name = :throws, par0 = 1)
    @test_throws TypeError TypeModel(; name = :throws, par1 = 1.5)
    @test_throws TypeError TypeModel(; name = :throws, par2 = 1.5)
    @test_throws TypeError TypeModel(; name = :throws, par3 = true)
    @test_throws TypeError TypeModel(; name = :throws, par4 = true)
    # par7 should be an AbstractArray of BigFloat.
    @test_throws MethodError TypeModel(; name = :throws, par7 = rand(Int, 3, 3))

    # Test that array types are correctly added.
    @named type_model2 = TypeModel(; par5 = rand(BigFloat, 3))
    @test symtype(type_model2.par5[1]) == BigFloat

    @named type_model3 = TypeModel(; par7 = rand(BigFloat, 3, 3))
    @test symtype(type_model3.par7[1, 1]) == BigFloat

    # Ensure that instances of models with conditional arrays with types can be created.
    @mtkmodel TypeCondition begin
        @structural_parameters begin
            flag
        end
        if flag
            @parameters begin
                k_if(t)[1:3, 1:3]::Float64, [description = "when true"]
            end
        else
            @parameters begin
                k_else[1:3]::Float64, [description = "when false"]
            end
        end
    end

    @named type_condition1 = TypeCondition(; flag = true, k_if = rand(Float64, 3, 3))
    @test symtype(type_condition1.k_if[1, 2]) == Float64

    @named type_condition2 = TypeCondition(; flag = false, k_else = rand(Float64, 3))
    @test symtype(type_condition2.k_else[1]) == Float64
end

@testset "Defaults of subcomponents MTKModel" begin
    @mtkmodel A begin
        @parameters begin
            p
        end
        @components begin
            b = B(i = p, j = 1 / p, k = 1)
        end
    end

    @mtkmodel B begin
        @parameters begin
            i
            j
            k
        end
    end

    @named a = A(p = 10)
    params = get_ps(a)
    @test isequal(getdefault(a.b.i), params[1])
    @test isequal(getdefault(a.b.j), 1 / params[1])
    @test getdefault(a.b.k) == 1

    @named a = A(p = 10, b.i = 20, b.j = 30, b.k = 40)
    @test getdefault(a.b.i) == 20
    @test getdefault(a.b.j) == 30
    @test getdefault(a.b.k) == 40
end

@testset "Metadata in variables" begin
    metadata = Dict(:description => "Variable to test metadata in the Model.structure",
        :input => true, :bounds => :((-1, 1)), :connection_type => :Flow,
        :tunable => false, :disturbance => true, :dist => :(Normal(1, 1)))

    @connector MockMeta begin
        m(t),
        [description = "Variable to test metadata in the Model.structure",
            input = true, bounds = (-1, 1), connect = Flow,
            tunable = false, disturbance = true, dist = Normal(1, 1)]
    end

    for (k, v) in metadata
        @test MockMeta.structure[:variables][:m][k] == v
    end
end

@testset "Connector with parameters, equations..." begin
    @connector A begin
        @extend (e,) = extended_e = E()
        @icon "pin.png"
        @parameters begin
            p
        end
        @variables begin
            v(t)
        end
        @components begin
            cc = C()
        end
        @equations begin
            e ~ 0
        end
    end

    @connector C begin
        c(t)
    end

    @connector E begin
        e(t)
    end

    @named aa = A()
    @test get_connector_type(aa) == RegularConnector()

    @test A.isconnector == true

    @test A.structure[:parameters] == Dict(:p => Dict(:type => Real))
    @test A.structure[:extend] == [[:e], :extended_e, :E]
    @test A.structure[:equations] == ["e ~ 0"]
    @test A.structure[:kwargs] == Dict{Symbol, Dict}(
        :p => Dict{Symbol, Union{Nothing, DataType}}(:value => nothing, :type => Real),
        :v => Dict{Symbol, Union{Nothing, DataType}}(:value => nothing, :type => Real))
    @test A.structure[:components] == [[:cc, :C]]
end

using ModelingToolkit: D_nounits
@testset "Event handling in MTKModel" begin
    @mtkmodel M begin
        @variables begin
            x(t)
            y(t)
            z(t)
        end
        @equations begin
            x ~ -D_nounits(x)
            D_nounits(y) ~ 0
            D_nounits(z) ~ 0
        end
        @continuous_events begin
            [x ~ 1.5] => [x ~ 5, y ~ 1]
        end
        @discrete_events begin
            (t == 1.5) => [x ~ x + 5, z ~ 2]
        end
    end

    @mtkbuild model = M()
    u0 = [model.x => 10, model.y => 0, model.z => 0]

    prob = ODEProblem(model, u0, (0, 5.0))
    sol = solve(prob, Tsit5(), tstops = [1.5])

    @test isequal(sol[model.y][end], 1.0)
    @test isequal(sol[model.z][end], 2.0)
end

# Ensure that modules consisting MTKModels with component arrays and icons of
# `Expr` type and `unit` metadata can be precompiled.
module PrecompilationTest
push!(LOAD_PATH, joinpath(@__DIR__, "precompile_test"))
using Unitful, Test, ModelParsingPrecompile, ModelingToolkit
using ModelingToolkit: getdefault, scalarize
@testset "Precompile packages with MTKModels" begin
    using ModelParsingPrecompile: ModelWithComponentArray

    @named model_with_component_array = ModelWithComponentArray()

    @test eval(ModelWithComponentArray.structure[:parameters][:r][:unit]) ==
          eval(u"Ω")
    @test lastindex(parameters(model_with_component_array)) == 3

    # Test the constant `k`. Manually k's value should be kept in sync here
    # and the ModelParsingPrecompile.
    @test all(getdefault.(getdefault.(scalarize(model_with_component_array.r))) .== 1)

    pop!(LOAD_PATH)
end
end

@testset "Conditional statements inside the blocks" begin
    @mtkmodel C begin end

    # Conditional statements inside @components, @equations
    # Conditional default value of parameters and variables
    @mtkmodel InsideTheBlock begin
        @structural_parameters begin
            flag = 1
        end
        @parameters begin
            eq = flag == 1 ? 1 : 0
            if flag == 1
                if_parameter = 100
            elseif flag == 2
                elseif_parameter = 101
            else
                else_parameter = 102
            end
        end
        @components begin
            default_sys = C()
            if flag == 1
                if_sys = C()
            elseif flag == 2
                elseif_sys = C()
            else
                else_sys = C()
            end
        end
        @equations begin
            eq ~ 0
            if flag == 1
                eq ~ 1
            elseif flag == 2
                eq ~ 2
            else
                eq ~ 3
            end
            flag == 1 ? eq ~ 4 : eq ~ 5
        end
    end

    @named if_in_sys = InsideTheBlock()
    if_in_sys = complete(if_in_sys; flatten = false)
    @named elseif_in_sys = InsideTheBlock(flag = 2)
    elseif_in_sys = complete(elseif_in_sys; flatten = false)
    @named else_in_sys = InsideTheBlock(flag = 3)
    else_in_sys = complete(else_in_sys; flatten = false)

    @test sort(getname.(parameters(if_in_sys))) == [:eq, :if_parameter]
    @test sort(getname.(parameters(elseif_in_sys))) == [:elseif_parameter, :eq]
    @test sort(getname.(parameters(else_in_sys))) == [:else_parameter, :eq]

    @test getdefault(if_in_sys.if_parameter) == 100
    @test getdefault(elseif_in_sys.elseif_parameter) == 101
    @test getdefault(else_in_sys.else_parameter) == 102

    @test nameof.(get_systems(if_in_sys)) == [:if_sys, :default_sys]
    @test nameof.(get_systems(elseif_in_sys)) == [:elseif_sys, :default_sys]
    @test nameof.(get_systems(else_in_sys)) == [:else_sys, :default_sys]

    @test all([
        if_in_sys.eq ~ 0,
        if_in_sys.eq ~ 1,
        if_in_sys.eq ~ 4
    ] .∈ [equations(if_in_sys)])
    @test all([
        elseif_in_sys.eq ~ 0,
        elseif_in_sys.eq ~ 2,
        elseif_in_sys.eq ~ 5
    ] .∈ [equations(elseif_in_sys)])
    @test all([
        else_in_sys.eq ~ 0,
        else_in_sys.eq ~ 3,
        else_in_sys.eq ~ 5
    ] .∈ [equations(else_in_sys)])

    @test getdefault(if_in_sys.eq) == 1
    @test getdefault(elseif_in_sys.eq) == 0
end

@testset "Conditional statements outside the blocks" begin
    @mtkmodel C begin end

    # Branching statement outside the begin blocks
    @mtkmodel OutsideTheBlock begin
        @structural_parameters begin
            condition = 0
        end

        @parameters begin
            default_parameter
        end
        @components begin
            default_sys = C()
        end
        @equations begin
            default_parameter ~ 0
        end

        if condition == 1
            @parameters begin
                if_parameter = 100
            end
            @equations begin
                if_parameter ~ 0
            end
            @components begin
                if_sys = C()
            end
        elseif condition == 2
            @parameters begin
                elseif_parameter = 101
            end
            @equations begin
                elseif_parameter ~ 0
            end
            @components begin
                elseif_sys = C()
            end
        else
            @parameters begin
                else_parameter = 102
            end
            @equations begin
                else_parameter ~ 0
            end
            @components begin
                else_sys = C()
            end
        end
    end

    @named if_out_sys = OutsideTheBlock(condition = 1)
    if_out_sys = complete(if_out_sys; flatten = false)
    @named elseif_out_sys = OutsideTheBlock(condition = 2)
    elseif_out_sys = complete(elseif_out_sys; flatten = false)
    @named else_out_sys = OutsideTheBlock(condition = 10)
    else_out_sys = complete(else_out_sys; flatten = false)
    @named ternary_out_sys = OutsideTheBlock(condition = 4)
    else_out_sys = complete(else_out_sys; flatten = false)

    @test getname.(parameters(if_out_sys)) == [:if_parameter, :default_parameter]
    @test getname.(parameters(elseif_out_sys)) == [:elseif_parameter, :default_parameter]
    @test getname.(parameters(else_out_sys)) == [:else_parameter, :default_parameter]

    @test getdefault(if_out_sys.if_parameter) == 100
    @test getdefault(elseif_out_sys.elseif_parameter) == 101
    @test getdefault(else_out_sys.else_parameter) == 102

    @test nameof.(get_systems(if_out_sys)) == [:if_sys, :default_sys]
    @test nameof.(get_systems(elseif_out_sys)) == [:elseif_sys, :default_sys]
    @test nameof.(get_systems(else_out_sys)) == [:else_sys, :default_sys]

    @test Equation[if_out_sys.if_parameter ~ 0
                   if_out_sys.default_parameter ~ 0] == equations(if_out_sys)
    @test Equation[elseif_out_sys.elseif_parameter ~ 0
                   elseif_out_sys.default_parameter ~ 0] == equations(elseif_out_sys)
    @test Equation[else_out_sys.else_parameter ~ 0
                   else_out_sys.default_parameter ~ 0] == equations(else_out_sys)

    @mtkmodel TernaryBranchingOutsideTheBlock begin
        @structural_parameters begin
            condition = true
        end
        condition ? begin
            @parameters begin
                ternary_parameter_true
            end
            @equations begin
                ternary_parameter_true ~ 0
            end
            @components begin
                ternary_sys_true = C()
            end
        end : begin
            @parameters begin
                ternary_parameter_false
            end
            @equations begin
                ternary_parameter_false ~ 0
            end
            @components begin
                ternary_sys_false = C()
            end
        end
    end

    @named ternary_true = TernaryBranchingOutsideTheBlock()
    ternary_true = complete(ternary_true; flatten = false)

    @named ternary_false = TernaryBranchingOutsideTheBlock(condition = false)
    ternary_false = complete(ternary_false; flatten = false)

    @test getname.(parameters(ternary_true)) == [:ternary_parameter_true]
    @test getname.(parameters(ternary_false)) == [:ternary_parameter_false]

    @test nameof.(get_systems(ternary_true)) == [:ternary_sys_true]
    @test nameof.(get_systems(ternary_false)) == [:ternary_sys_false]

    @test Equation[ternary_true.ternary_parameter_true ~ 0] == equations(ternary_true)
    @test Equation[ternary_false.ternary_parameter_false ~ 0] == equations(ternary_false)
end

_b = Ref{Any}()
@mtkmodel MyModel begin
    @variables begin
        x___(t) = 0
    end
    begin
        _b[] = x___
    end
end
@named m = MyModel()
@variables x___(t)
@test isequal(x___, _b[])

@testset "Component array" begin
    @mtkmodel SubComponent begin
        @parameters begin
            sc
        end
    end

    @mtkmodel Component begin
        @structural_parameters begin
            N = 2
        end
        @components begin
            comprehension = [SubComponent(sc = i) for i in 1:N]
            written_out_for = for i in 1:N
                sc = i + 1
                SubComponent(; sc)
            end
            single_sub_component = SubComponent()
        end
    end

    @named component = Component()
    component = complete(component; flatten = false)

    @test nameof.(ModelingToolkit.get_systems(component)) == [
        :comprehension_1,
        :comprehension_2,
        :written_out_for_1,
        :written_out_for_2,
        :single_sub_component
    ]

    @test getdefault(component.comprehension_1.sc) == 1
    @test getdefault(component.comprehension_2.sc) == 2
    @test getdefault(component.written_out_for_1.sc) == 2
    @test getdefault(component.written_out_for_2.sc) == 3

    @mtkmodel ConditionalComponent begin
        @structural_parameters begin
            N = 2
        end
        @components begin
            if N == 2
                if_comprehension = [SubComponent(sc = i) for i in 1:N]
            elseif N == 3
                elseif_comprehension = [SubComponent(sc = i) for i in 1:N]
            else
                else_comprehension = [SubComponent(sc = i) for i in 1:N]
            end
        end
    end

    @named if_component = ConditionalComponent()
    @test nameof.(get_systems(if_component)) == [:if_comprehension_1, :if_comprehension_2]

    @named elseif_component = ConditionalComponent(; N = 3)
    @test nameof.(get_systems(elseif_component)) ==
          [:elseif_comprehension_1, :elseif_comprehension_2, :elseif_comprehension_3]

    @named else_component = ConditionalComponent(; N = 4)
    @test nameof.(get_systems(else_component)) ==
          [:else_comprehension_1, :else_comprehension_2,
        :else_comprehension_3, :else_comprehension_4]
end

@testset "Parent module of Models" begin
    @test parentmodule(MyMockModule.Ground) == MyMockModule
end

@testset "Guesses with expression" begin
    @mtkmodel GuessModel begin
        @variables begin
            k(t)
            l(t) = 10, [guess = k, unit = u"A"]
            i(t), [guess = k, unit = u"A"]
            j(t), [guess = k + l / i]
        end
    end

    @named guess_model = GuessModel()

    j_guess = getguess(guess_model.j)
    @test typeof(j_guess) == Num
    @test readable_code(j_guess) == "l(t) / i(t) + k(t)"

    i_guess = getguess(guess_model.i)
    @test typeof(i_guess) == Num
    @test readable_code(i_guess) == "k(t)"

    l_guess = getguess(guess_model.l)
    @test typeof(l_guess) == Num
    @test readable_code(l_guess) == "k(t)"
end

@testset "Argument order" begin
    @mtkmodel OrderModel begin
        @structural_parameters begin
            b = 1 # reverse alphabetical order to test that the order is preserved
            a = b
        end
        @parameters begin
            c = a
            d = b
        end
    end
    @named ordermodel = OrderModel()
    ordermodel = complete(ordermodel)
    defs = ModelingToolkit.defaults(ordermodel)
    @test defs[ordermodel.c] == 1
    @test defs[ordermodel.d] == 1

    @test_nowarn @named ordermodel = OrderModel(a = 2)
    ordermodel = complete(ordermodel)
    defs = ModelingToolkit.defaults(ordermodel)
    @test defs[ordermodel.c] == 2
    @test defs[ordermodel.d] == 1
end

@testset "Vector defaults" begin
    @mtkmodel VectorDefaultWithMetadata begin
        @parameters begin
            n[1:3] = [1, 2, 3], [description = "Vector defaults"]
        end
    end

    @named vec = VectorDefaultWithMetadata()
    for i in 1:3
        @test getdefault(vec.n[i]) == i
    end

    @mtkmodel VectorConditionalDefault begin
        @structural_parameters begin
            flag = true
        end
        @parameters begin
            n[1:3] = if flag
                [2, 2, 2]
            else
                [1, 1, 1]
            end
        end
    end

    @named vec_true = VectorConditionalDefault()
    for i in 1:3
        @test getdefault(vec_true.n[i]) == 2
    end
    @named vec_false = VectorConditionalDefault(flag = false)
    for i in 1:3
        @test getdefault(vec_false.n[i]) == 1
    end
end

@testset "Duplicate names" begin
    mod = @__MODULE__
    @test_throws ErrorException ModelingToolkit._model_macro(mod, :ATest,
        :(begin
            @variables begin
                a(t)
                a(t)
            end
        end),
        false)
    @test_throws ErrorException ModelingToolkit._model_macro(mod, :ATest,
        :(begin
            @variables begin
                a(t)
            end
            @parameters begin
                a
            end
        end),
        false)
end

@mtkmodel BaseSys begin
    @parameters begin
        p1
        p2
    end
    @variables begin
        v1(t)
    end
end

@testset "Arguments of base system" begin
    @mtkmodel MainSys begin
        @extend BaseSys(p1 = 1)
    end

    @test names(MainSys) == [:p2, :p1, :v1]
    @named main_sys = MainSys(p1 = 11, p2 = 12, v1 = 13)
    @test getdefault(main_sys.p1) == 11
    @test getdefault(main_sys.p2) == 12
    @test getdefault(main_sys.v1) == 13
end

@mtkmodel InnerModel begin
    @parameters begin
        p
    end
end

@mtkmodel MidModel begin
    @components begin
        inmodel = InnerModel()
    end
end

@mtkmodel MidModelB begin
    @parameters begin
        b
    end
    @components begin
        inmodel_b = InnerModel()
    end
end

@mtkmodel OuterModel begin
    @extend MidModel()
    @equations begin
        inmodel.p ~ 0
    end
end

# The base system is fetched from the module while extending implicitly. This
# way of defining fails when defined inside the `@testset`. So, it is moved out.
@testset "Test unpacking of components in implicit extend" begin
    @named out = OuterModel()
    @test OuterModel.structure[:extend][1] == [:inmodel]
end

@mtkmodel MultipleExtend begin
    @extend MidModel()
    @extend MidModelB()
end

@testset "Multiple extend statements" begin
    @named multiple_extend = MultipleExtend()
    @test collect(nameof.(multiple_extend.systems)) == [:inmodel_b, :inmodel]
    @test MultipleExtend.structure[:extend][1] == [:inmodel, :b, :inmodel_b]
    @test tosymbol.(parameters(multiple_extend)) == [:b, :inmodel_b₊p, :inmodel₊p]
end

struct CustomStruct end
@testset "Nonnumeric parameters" begin
    @mtkmodel MyModel begin
        @parameters begin
            p::CustomStruct
        end
    end
    @named sys = MyModel(p = CustomStruct())
    @test ModelingToolkit.defaults(sys)[@nonamespace sys.p] == CustomStruct()
end

@testset "Variables are not callable symbolics" begin
    @mtkmodel Example begin
        @variables begin
            x(t)
            y(t)
        end
        @equations begin
            x ~ y
        end
    end
    @named ex = Example()
    vars = Symbolics.get_variables(only(equations(ex)))
    @test length(vars) == 2
    for u in Symbolics.unwrap.(unknowns(ex))
        @test !Symbolics.hasmetadata(u, Symbolics.CallWithParent)
        @test any(isequal(u), vars)
    end
end
