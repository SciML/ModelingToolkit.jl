PrecompileTools.@compile_workload begin
        fold1 = Val{false}()
        using SymbolicUtils
        using SymbolicUtils: shape
        using Symbolics
        @syms x y f(t) q[1:5]
        SymbolicUtils.Sym{SymReal}(:a; type = Real, shape = SymbolicUtils.ShapeVecT())
        x + y
        x * y
        x / y
        x ^ y
        x ^ 5
        6 ^ x
        x - y
        x * x * q[1]
        -y
        2y
        z = 2
        dict = SymbolicUtils.ACDict{VartypeT}()
        dict[x] = 1
        dict[y] = 1
        type::typeof(DataType) = rand() < 0.5 ? Real : Float64
        nt = (; type, shape, unsafe = true)
        Base.pairs(nt)
        BSImpl.AddMul{VartypeT}(1, dict, SymbolicUtils.AddMulVariant.MUL; type, shape = SymbolicUtils.ShapeVecT(), unsafe = true)
        *(y, z)
        *(z, y)
        SymbolicUtils.symtype(y)
        f(x)
        (5x / 5)
        expand((x + y) ^ 2)
        simplify(x ^ (1//2) + (sin(x) ^ 2 + cos(x) ^ 2) + 2(x + y) - x - y)
        ex = x + 2y + sin(x)
        rules1 = Dict(x => y)
        rules2 = Dict(x => 1)
        Dx = Differential(x)
        Differential(y)(ex)
        uex = unwrap(ex)
        Symbolics.executediff(Dx, uex)
        # Running `fold = Val(true)` invalidates the precompiled statements
        # for `fold = Val(false)` and itself doesn't precompile anyway.
        # substitute(ex, rules1)
        substitute(ex, rules1; fold = fold1)
        substitute(ex, rules2; fold = fold1)
        @variables foo
        f(foo)
        @variables x y f(::Real) q[1:5]
        x + y
        x * y
        x / y
        x ^ y
        x ^ 5
        # 6 ^ x
        x - y
        -y
        2y
        symtype(y)
        z = 2
        *(y, z)
        *(z, y)
        f(x)
        (5x / 5)
        [x, y]
        [x, f, f]
        promote_type(Int, Num)
        promote_type(Real, Num)
        promote_type(Float64, Num)
        # expand((x + y) ^ 2)
        # simplify(x ^ (1//2) + (sin(x) ^ 2 + cos(x) ^ 2) + 2(x + y) - x - y)
        ex = x + 2y + sin(x)
        rules1 = Dict(x => y)
        # rules2 = Dict(x => 1)
        # Running `fold = Val(true)` invalidates the precompiled statements
        # for `fold = Val(false)` and itself doesn't precompile anyway.
        # substitute(ex, rules1)
        substitute(ex, rules1; fold = fold1)
        Symbolics.linear_expansion(ex, y)
        # substitute(ex, rules2; fold = fold1)
        # substitute(ex, rules2)
        # substitute(ex, rules1; fold = fold2)
        # substitute(ex, rules2; fold = fold2)
        q[1]
        q'q
     using ModelingToolkitBase
    @parameters g
    @variables x(ModelingToolkitBase.t_nounits)
    @variables y(ModelingToolkitBase.t_nounits) [state_priority = 10]
    @variables λ(ModelingToolkitBase.t_nounits)
    eqs = [
        ModelingToolkitBase.D_nounits(ModelingToolkitBase.D_nounits(x)) ~ λ * x
        ModelingToolkitBase.D_nounits(ModelingToolkitBase.D_nounits(y)) ~ λ * y - g
        x^2 + y^2 ~ 1
    ]
    dvs = Num[x, y, λ]
    ps = Num[g]
    ics = Dict{SymbolicT, SymbolicT}()
    ics[y] = -1.0
    ics[ModelingToolkitBase.D_nounits(x)] = 0.5
    isequal(ModelingToolkitBase.D_nounits.x, ModelingToolkitBase.t_nounits)
    sys = System(eqs, ModelingToolkitBase.t_nounits, dvs, ps; initial_conditions = ics, guesses = ics, name = :sys)
    complete(sys)
    @static if @isdefined(ModelingToolkit)
        TearingState(sys)
    end
    mtkcompile(sys)
    @syms p[1:2]
    ndims(p)
    size(p)
    axes(p)
    length(p)
    v = [p]
    isempty(v)
    # mtkcompile(sys)
end

precompile(Tuple{typeof(SymbolicUtils.isequal_somescalar), Float64, Float64})
precompile(Tuple{typeof(Base.:(var"==")), ModelingToolkitBase.Initial, ModelingToolkitBase.Initial})
