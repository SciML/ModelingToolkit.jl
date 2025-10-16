"""
    lin_fun, simplified_sys = linearization_function(sys::AbstractSystem, inputs, outputs; simplify = false, initialize = true, initialization_solver_alg = nothing, kwargs...)

Return a function that linearizes the system `sys`. The function [`linearize`](@ref) provides a higher-level and easier to use interface.

`lin_fun` is a function `(variables, p, t) -> (; f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u)`, i.e., it returns a NamedTuple with the Jacobians of `f,g,h` for the nonlinear `sys` (technically for `simplified_sys`) on the form

```math
\\begin{aligned}
ẋ &= f(x, z, u) \\\\
0 &= g(x, z, u) \\\\
y &= h(x, z, u)
\\end{aligned}
```

where `x` are differential unknown variables, `z` algebraic variables, `u` inputs and `y` outputs. To obtain a linear statespace representation, see [`linearize`](@ref). The input argument `variables` is a vector defining the operating point, corresponding to `unknowns(simplified_sys)` and `p` is a vector corresponding to the parameters of `simplified_sys`. Note: all variables in `inputs` have been converted to parameters in `simplified_sys`.

The `simplified_sys` has undergone [`mtkcompile`](@ref) and had any occurring input or output variables replaced with the variables provided in arguments `inputs` and `outputs`. The unknowns of this system also indicate the order of the unknowns that holds for the linearized matrices.

# Arguments:

  - `sys`: A [`System`](@ref) of ODEs. This function will automatically apply simplification passes on `sys` and return the resulting `simplified_sys`.
  - `inputs`: A vector of variables that indicate the inputs of the linearized input-output model.
  - `outputs`: A vector of variables that indicate the outputs of the linearized input-output model.
  - `simplify`: Apply simplification in tearing.
  - `initialize`: If true, a check is performed to ensure that the operating point is consistent (satisfies algebraic equations). If the op is not consistent, initialization is performed.
  - `initialization_solver_alg`: A NonlinearSolve algorithm to use for solving for a feasible set of state and algebraic variables that satisfies the specified operating point.
  - `autodiff`: An `ADType` supported by DifferentiationInterface.jl to use for calculating the necessary jacobians. Defaults to using `AutoForwardDiff()`
  - `kwargs`: Are passed on to `find_solvables!`

See also [`linearize`](@ref) which provides a higher-level interface.
"""
function linearization_function(sys::AbstractSystem, inputs,
        outputs; simplify = false,
        initialize = true,
        initializealg = nothing,
        initialization_abstol = 1e-5,
        initialization_reltol = 1e-3,
        op = Dict(),
        p = DiffEqBase.NullParameters(),
        zero_dummy_der = false,
        initialization_solver_alg = nothing,
        autodiff = AutoForwardDiff(),
        eval_expression = false, eval_module = @__MODULE__,
        warn_initialize_determined = true,
        guesses = Dict(),
        warn_empty_op = true,
        t = 0.0,
        kwargs...)
    op = Dict(op)
    if isempty(op) && warn_empty_op
        @warn "An empty operating point was passed to `linearization_function`. An operating point containing the variables that will be changed in `linearize` should be provided. Disable this warning by passing `warn_empty_op = false`."
    end
    inputs isa AbstractVector || (inputs = [inputs])
    outputs isa AbstractVector || (outputs = [outputs])
    inputs = mapreduce(vcat, inputs; init = []) do var
        symbolic_type(var) == ArraySymbolic() ? collect(var) : [var]
    end
    outputs = mapreduce(vcat, outputs; init = []) do var
        symbolic_type(var) == ArraySymbolic() ? collect(var) : [var]
    end
    ssys = mtkcompile(sys; inputs, outputs, simplify, kwargs...)
    diff_idxs, alge_idxs = eq_idxs(ssys)
    if zero_dummy_der
        dummyder = setdiff(unknowns(ssys), unknowns(sys))
        defs = Dict(x => 0.0 for x in dummyder)
        @set! ssys.defaults = merge(defs, defaults(ssys))
        op = merge(defs, op)
    end
    sys = ssys

    if initializealg === nothing
        initializealg = initialize ? OverrideInit() : NoInit()
    end

    prob = ODEProblem{true, SciMLBase.FullSpecialize}(
        sys, merge(op, anydict(p)), (t, t); allow_incomplete = true,
        algebraic_only = true, guesses)
    u0 = state_values(prob)

    ps = parameters(sys)
    h = build_explicit_observed_function(sys, outputs; eval_expression, eval_module)

    initialization_kwargs = (;
        abstol = initialization_abstol, reltol = initialization_reltol,
        nlsolve_alg = initialization_solver_alg)

    p = parameter_values(prob)
    t0 = current_time(prob)
    inputvals = [prob.ps[i] for i in inputs]

    hp_fun = let fun = h, setter = setp_oop(sys, inputs)
        function hpf(du, input, u, p, t)
            p = setter(p, input)
            fun(du, u, p, t)
            return du
        end
    end
    if u0 === nothing
        uf_jac = h_jac = pf_jac = nothing
        T = p isa MTKParameters ? eltype(p.tunable) : eltype(p)
        hp_jac = PreparedJacobian{true}(
            hp_fun, zeros(T, size(outputs)), autodiff, inputvals,
            DI.Constant(prob.u0), DI.Constant(p), DI.Constant(t0))
    else
        uf_fun = let fun = prob.f
            function uff(du, u, p, t)
                SciMLBase.UJacobianWrapper(fun, t, p)(du, u)
            end
        end
        uf_jac = PreparedJacobian{true}(
            uf_fun, similar(prob.u0), autodiff, prob.u0, DI.Constant(p), DI.Constant(t0))
        # observed function is a `GeneratedFunctionWrapper` with iip component
        h_jac = PreparedJacobian{true}(h, similar(prob.u0, size(outputs)), autodiff,
            prob.u0, DI.Constant(p), DI.Constant(t0))
        pf_fun = let fun = prob.f, setter = setp_oop(sys, inputs)
            function pff(du, input, u, p, t)
                p = setter(p, input)
                SciMLBase.ParamJacobianWrapper(fun, t, u)(du, p)
            end
        end
        pf_jac = PreparedJacobian{true}(pf_fun, similar(prob.u0), autodiff, inputvals,
            DI.Constant(prob.u0), DI.Constant(p), DI.Constant(t0))
        hp_jac = PreparedJacobian{true}(
            hp_fun, similar(prob.u0, size(outputs)), autodiff, inputvals,
            DI.Constant(prob.u0), DI.Constant(p), DI.Constant(t0))
    end

    lin_fun = LinearizationFunction(
        diff_idxs, alge_idxs, inputs, length(unknowns(sys)),
        prob, h, u0 === nothing ? nothing : similar(u0), uf_jac, h_jac, pf_jac,
        hp_jac, initializealg, initialization_kwargs)
    return lin_fun, sys
end

"""
Return the set of indexes of differential equations and algebraic equations in the simplified system.
"""
function eq_idxs(sys::AbstractSystem)
    eqs = equations(sys)
    alge_idxs = findall(!isdiffeq, eqs)
    diff_idxs = setdiff(1:length(eqs), alge_idxs)

    diff_idxs, alge_idxs
end

"""
    $(TYPEDEF)

Callable struct which stores a function and its prepared `DI.jacobian`. Calling with the
appropriate arguments for DI returns the jacobian.

# Fields

$(TYPEDFIELDS)
"""
struct PreparedJacobian{iip, P, F, B, A}
    """
    The preparation object.
    """
    prep::P
    """
    The function whose jacobian is calculated.
    """
    f::F
    """
    Buffer for in-place functions.
    """
    buf::B
    """
    ADType to use for differentiation.
    """
    autodiff::A
end

function PreparedJacobian{true}(f, buf, autodiff, args...)
    prep = DI.prepare_jacobian(f, buf, autodiff, args...; strict = Val(false))
    return PreparedJacobian{true, typeof(prep), typeof(f), typeof(buf), typeof(autodiff)}(
        prep, f, buf, autodiff)
end

function PreparedJacobian{false}(f, autodiff, args...)
    prep = DI.prepare_jacobian(f, autodiff, args...; strict = Val(false))
    return PreparedJacobian{true, typeof(prep), typeof(f), Nothing, typeof(autodiff)}(
        prep, f, nothing)
end

function (pj::PreparedJacobian{true})(args...)
    DI.jacobian(pj.f, pj.buf, pj.prep, pj.autodiff, args...)
end

function (pj::PreparedJacobian{false})(args...)
    DI.jacobian(pj.f, pj.prep, pj.autodiff, args...)
end

"""
    $(TYPEDEF)

A callable struct which linearizes a system.

# Fields

$(TYPEDFIELDS)
"""
struct LinearizationFunction{
    DI <: AbstractVector{Int}, AI <: AbstractVector{Int}, I, P <: ODEProblem,
    H, C, J1, J2, J3, J4, IA <: SciMLBase.DAEInitializationAlgorithm, IK}
    """
    The indexes of differential equations in the linearized system.
    """
    diff_idxs::DI
    """
    The indexes of algebraic equations in the linearized system.
    """
    alge_idxs::AI
    """
    The indexes of parameters in the linearized system which represent
    input variables.
    """
    inputs::I
    """
    The number of unknowns in the linearized system.
    """
    num_states::Int
    """
    The `ODEProblem` of the linearized system.
    """
    prob::P
    """
    A function which takes `(u, p, t)` and returns the outputs of the linearized system.
    """
    h::H
    """
    Any required cache buffers.
    """
    caches::C
    """
    `PreparedJacobian` for calculating jacobian of `prob.f` w.r.t. `u`
    """
    uf_jac::J1
    """
    `PreparedJacobian` for calculating jacobian of `h` w.r.t. `u`
    """
    h_jac::J2
    """
    `PreparedJacobian` for calculating jacobian of `prob.f` w.r.t. `p`
    """
    pf_jac::J3
    """
    `PreparedJacobian` for calculating jacobian of `h` w.r.t. `p`
    """
    hp_jac::J4
    """
    The initialization algorithm to use.
    """
    initializealg::IA
    """
    Keyword arguments to be passed to `SciMLBase.get_initial_values`.
    """
    initialize_kwargs::IK
end

SymbolicIndexingInterface.symbolic_container(f::LinearizationFunction) = f.prob
SymbolicIndexingInterface.state_values(f::LinearizationFunction) = state_values(f.prob)
function SymbolicIndexingInterface.parameter_values(f::LinearizationFunction)
    parameter_values(f.prob)
end
SymbolicIndexingInterface.current_time(f::LinearizationFunction) = current_time(f.prob)

function Base.show(io::IO, mime::MIME"text/plain", lf::LinearizationFunction)
    printstyled(io, "LinearizationFunction"; bold = true, color = :blue)
    println(io, " which wraps:")
    show(io, mime, lf.prob)
end

"""
    $(TYPEDSIGNATURES)

Linearize the wrapped system at the point given by `(unknowns, p, t)`.
"""
function (linfun::LinearizationFunction)(u, p, t)
    if eltype(p) <: Pair
        p = todict(p)
        newps = copy(parameter_values(linfun.prob))
        for (k, v) in p
            if is_parameter(linfun, k)
                v = fixpoint_sub(v, p)
                setp(linfun, k)(newps, v)
            end
        end
        p = newps
    end

    fun = linfun.prob.f
    input_vals = [linfun.prob.ps[i] for i in linfun.inputs]
    if u !== nothing # Handle systems without unknowns
        linfun.num_states == length(u) ||
            error("Number of unknown variables ($(linfun.num_states)) does not match the number of input unknowns ($(length(u)))")
        integ_cache = (linfun.caches,)
        integ = MockIntegrator{true}(u, p, t, fun, integ_cache, nothing)
        u, p,
        success = SciMLBase.get_initial_values(
            linfun.prob, integ, fun, linfun.initializealg, Val(true);
            linfun.initialize_kwargs...)
        if !success
            error("Initialization algorithm $(linfun.initializealg) failed with `unknowns = $u` and `p = $p`.")
        end
        fg_xz = linfun.uf_jac(u, DI.Constant(p), DI.Constant(t))
        h_xz = linfun.h_jac(u, DI.Constant(p), DI.Constant(t))
        fg_u = linfun.pf_jac(input_vals,
            DI.Constant(u), DI.Constant(p), DI.Constant(t))
    else
        linfun.num_states == 0 ||
            error("Number of unknown variables (0) does not match the number of input unknowns ($(length(u)))")
        fg_xz = zeros(0, 0)
        h_xz = fg_u = zeros(0, length(linfun.inputs))
    end
    h_u = linfun.hp_jac(input_vals,
        DI.Constant(u), DI.Constant(p), DI.Constant(t))
    (f_x = fg_xz[linfun.diff_idxs, linfun.diff_idxs],
        f_z = fg_xz[linfun.diff_idxs, linfun.alge_idxs],
        g_x = fg_xz[linfun.alge_idxs, linfun.diff_idxs],
        g_z = fg_xz[linfun.alge_idxs, linfun.alge_idxs],
        f_u = fg_u[linfun.diff_idxs, :],
        g_u = fg_u[linfun.alge_idxs, :],
        h_x = h_xz[:, linfun.diff_idxs],
        h_z = h_xz[:, linfun.alge_idxs],
        h_u = h_u,
        x = u,
        p,
        t)
end

"""
    $(TYPEDEF)

Mock `DEIntegrator` to allow using `CheckInit` without having to create a new integrator
(and consequently depend on `OrdinaryDiffEq`).

# Fields

$(TYPEDFIELDS)
"""
struct MockIntegrator{iip, U, P, T, F, C, O} <: SciMLBase.DEIntegrator{Nothing, iip, U, T}
    """
    The state vector.
    """
    u::U
    """
    The parameter object.
    """
    p::P
    """
    The current time.
    """
    t::T
    """
    The wrapped `SciMLFunction`.
    """
    f::F
    """
    The integrator cache.
    """
    cache::C
    """
    Integrator "options" for `CheckInit`.
    """
    opts::O
end

function MockIntegrator{iip}(
        u::U, p::P, t::T, f::F, cache::C, opts::O) where {iip, U, P, T, F, C, O}
    return MockIntegrator{iip, U, P, T, F, C, O}(u, p, t, f, cache, opts)
end

SymbolicIndexingInterface.state_values(integ::MockIntegrator) = integ.u
SymbolicIndexingInterface.parameter_values(integ::MockIntegrator) = integ.p
SymbolicIndexingInterface.current_time(integ::MockIntegrator) = integ.t
SciMLBase.get_tmp_cache(integ::MockIntegrator) = integ.cache

"""
    $(TYPEDEF)

A struct representing a linearization operation to be performed. Can be symbolically
indexed to efficiently update the operating point for multiple linearizations in a loop.
The value of the independent variable can be set by mutating the `.t` field of this struct.
"""
mutable struct LinearizationProblem{F <: LinearizationFunction, T}
    """
    The wrapped `LinearizationFunction`
    """
    const f::F
    t::T
end

function Base.show(io::IO, mime::MIME"text/plain", prob::LinearizationProblem)
    printstyled(io, "LinearizationProblem"; bold = true, color = :blue)
    println(io, " at time ", prob.t, " which wraps:")
    show(io, mime, prob.f.prob)
end

"""
    $(TYPEDSIGNATURES)

Construct a `LinearizationProblem` for linearizing the system `sys` with the given
`inputs` and `outputs`.

# Keyword arguments

- `t`: The value of the independent variable

All other keyword arguments are forwarded to `linearization_function`.
"""
function LinearizationProblem(sys::AbstractSystem, inputs, outputs; t = 0.0, kwargs...)
    linfun, _ = linearization_function(sys, inputs, outputs; kwargs...)
    return LinearizationProblem(linfun, t)
end

SymbolicIndexingInterface.symbolic_container(p::LinearizationProblem) = p.f
SymbolicIndexingInterface.state_values(p::LinearizationProblem) = state_values(p.f)
SymbolicIndexingInterface.parameter_values(p::LinearizationProblem) = parameter_values(p.f)
SymbolicIndexingInterface.current_time(p::LinearizationProblem) = p.t

function Base.getindex(prob::LinearizationProblem, idx)
    getu(prob, idx)(prob)
end

function Base.setindex!(prob::LinearizationProblem, val, idx)
    setu(prob, idx)(prob, val)
end

function Base.getproperty(prob::LinearizationProblem, x::Symbol)
    if x == :ps
        return ParameterIndexingProxy(prob)
    end
    return getfield(prob, x)
end

function CommonSolve.solve(prob::LinearizationProblem; allow_input_derivatives = false)
    u0 = state_values(prob)
    p = parameter_values(prob)
    t = current_time(prob)
    linres = prob.f(u0, p, t)
    f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u, x, p, t = linres

    nx, nu = size(f_u)
    nz = size(f_z, 2)
    ny = size(h_x, 1)

    D = h_u

    if isempty(g_z)
        A = f_x
        B = f_u
        C = h_x
        @assert iszero(g_x)
        @assert iszero(g_z)
        @assert iszero(g_u)
    else
        gz = lu(g_z; check = false)
        issuccess(gz) ||
            error("g_z not invertible, this indicates that the DAE is of index > 1.")
        gzgx = -(gz \ g_x)
        A = [f_x f_z
             gzgx*f_x gzgx*f_z]
        B = [f_u
             gzgx * f_u] # The cited paper has zeros in the bottom block, see derivation in https://github.com/SciML/ModelingToolkit.jl/pull/1691 for the correct formula

        C = [h_x h_z]
        Bs = -(gz \ g_u) # This equation differ from the cited paper, the paper is likely wrong since their equaiton leads to a dimension mismatch.
        if !iszero(Bs)
            if !allow_input_derivatives
                der_inds = findall(vec(any(!=(0), Bs, dims = 1)))
                error("Input derivatives appeared in expressions (-g_z\\g_u != 0), the following inputs appeared differentiated: $(inputs(prob.f.prob.f.sys)[der_inds]). Call `linearize` with keyword argument `allow_input_derivatives = true` to allow this and have the returned `B` matrix be of double width ($(2nu)), where the last $nu inputs are the derivatives of the first $nu inputs.")
            end
            B = [B [zeros(nx, nu); Bs]]
            D = [D zeros(ny, nu)]
        end
    end

    (; A, B, C, D), (; x, p, t)
end

"""
    (; A, B, C, D), simplified_sys = linearize_symbolic(sys::AbstractSystem, inputs, outputs; simplify = false, allow_input_derivatives = false, kwargs...)

Similar to [`linearize`](@ref), but returns symbolic matrices `A,B,C,D` rather than numeric. While `linearize` uses ForwardDiff to perform the linearization, this function uses `Symbolics.jacobian`.

See [`linearize`](@ref) for a description of the arguments.

# Extended help
The named tuple returned as the first argument additionally contains the jacobians `f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u` of
```math
\\begin{aligned}
ẋ &= f(x, z, u) \\\\
0 &= g(x, z, u) \\\\
y &= h(x, z, u)
\\end{aligned}
```
where `x` are differential unknown variables, `z` algebraic variables, `u` inputs and `y` outputs.
"""
function linearize_symbolic(sys::AbstractSystem, inputs,
        outputs; simplify = false, allow_input_derivatives = false,
        eval_expression = false, eval_module = @__MODULE__, split = true,
        kwargs...)
    sys = mtkcompile(sys; inputs, outputs, simplify, split, kwargs...)
    diff_idxs, alge_idxs = eq_idxs(sys)
    sts = unknowns(sys)
    t = get_iv(sys)
    ps = parameters(sys; initial_parameters = true)
    p = reorder_parameters(sys, ps)

    fun_expr = generate_rhs(sys; expression = Val{true})[1]
    fun = eval_or_rgf(fun_expr; eval_expression, eval_module)
    
    h = build_explicit_observed_function(sys, outputs; eval_expression, eval_module)
    if split
        dx = fun(sts, p, t)
        y = h(sts, p, t)
    else
        dx = fun(sts, p..., t)
        y = h(sts, p..., t)
    end

    fg_xz = Symbolics.jacobian(dx, sts)
    fg_u = Symbolics.jacobian(dx, inputs)
    h_xz = Symbolics.jacobian(y, sts)
    h_u = Symbolics.jacobian(y, inputs)
    f_x = fg_xz[diff_idxs, diff_idxs]
    f_z = fg_xz[diff_idxs, alge_idxs]
    g_x = fg_xz[alge_idxs, diff_idxs]
    g_z = fg_xz[alge_idxs, alge_idxs]
    f_u = fg_u[diff_idxs, :]
    g_u = fg_u[alge_idxs, :]
    h_x = h_xz[:, diff_idxs]
    h_z = h_xz[:, alge_idxs]

    nx, nu = size(f_u)
    nz = size(f_z, 2)
    ny = size(h_x, 1)

    D = h_u

    if isempty(g_z) # ODE
        A = f_x
        B = f_u
        C = h_x
    else
        gz = lu(g_z; check = false)
        issuccess(gz) ||
            error("g_z not invertible, this indicates that the DAE is of index > 1.")
        gzgx = -(gz \ g_x)
        A = [f_x f_z
             gzgx*f_x gzgx*f_z]
        B = [f_u
             gzgx * f_u] # The cited paper has zeros in the bottom block, see derivation in https://github.com/SciML/ModelingToolkit.jl/pull/1691 for the correct formula

        C = [h_x h_z]
        Bs = -(gz \ g_u) # This equation differ from the cited paper, the paper is likely wrong since their equaiton leads to a dimension mismatch.
        if !iszero(Bs)
            if !allow_input_derivatives
                der_inds = findall(vec(any(!iszero, Bs, dims = 1)))
                error("Input derivatives appeared in expressions (-g_z\\g_u != 0), the following inputs appeared differentiated: $(ModelingToolkit.inputs(sys)[der_inds]). Call `linearize_symbolic` with keyword argument `allow_input_derivatives = true` to allow this and have the returned `B` matrix be of double width ($(2nu)), where the last $nu inputs are the derivatives of the first $nu inputs.")
            end
            B = [B [zeros(nx, nu); Bs]]
            D = [D zeros(ny, nu)]
        end
    end

    (; A, B, C, D, f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u), sys
end

struct IONotFoundError <: Exception
    variant::String
    sysname::Symbol
    not_found::Any
end

function Base.showerror(io::IO, err::IONotFoundError)
    println(io, "The following $(err.variant) provided to `mtkcompile` were not found in the system:")
    maybe_namespace_issue = false
    for var in err.not_found
        println(io, "  ", var)
        if hasname(var) && startswith(string(getname(var)), string(err.sysname))
            maybe_namespace_issue = true
        end
    end
    if maybe_namespace_issue
        println(io, """
        Some of the missing variables are namespaced with the name of the system \
        `$(err.sysname)` passed to `mtkcompile`. This may be indicative of a namespacing \
        issue. `mtkcompile` requires that the $(err.variant) provided are not namespaced \
        with the name of the root system. This issue can occur when using `getproperty` \
        to access the variables passed as $(err.variant). For example:

        ```julia
        @named sys = MyModel()
        inputs = [sys.input_var]
        mtkcompile(sys; inputs)
        ```

        Here, `mtkcompile` expects the input to be named `input_var`, but since `sys`
        performs namespacing, it will be named `sys$(NAMESPACE_SEPARATOR)input_var`. To \
        fix this issue, namespacing can be temporarily disabled:

        ```julia
        @named sys = MyModel()
        sys_nns = toggle_namespacing(sys, false)
        inputs = [sys_nns.input_var]
        mtkcompile(sys; inputs)
        ```
        """)
    end
end

"""
Modify the variable metadata of system variables to indicate which ones are inputs, outputs, and disturbances. Needed for `inputs`, `outputs`, `disturbances`, `unbound_inputs`, `unbound_outputs` to return the proper subsets.
"""
function markio!(state, orig_inputs, inputs, outputs, disturbances; check = true)
    fullvars = get_fullvars(state)
    inputset = Dict{Any, Bool}(i => false for i in inputs)
    outputset = Dict{Any, Bool}(o => false for o in outputs)
    disturbanceset = Dict{Any, Bool}(d => false for d in disturbances)
    for (i, v) in enumerate(fullvars)
        if v in keys(inputset)
            if v in keys(outputset)
                v = setio(v, true, true)
                outputset[v] = true
            else
                v = setio(v, true, false)
            end
            inputset[v] = true
            fullvars[i] = v
        elseif v in keys(outputset)
            v = setio(v, false, true)
            outputset[v] = true
            fullvars[i] = v
        else
            if isinput(v)
                push!(orig_inputs, v)
            end
            v = setio(v, false, false)
            fullvars[i] = v
        end

        if v in keys(disturbanceset)
            v = setio(v, true, false)
            v = setdisturbance(v, true)
            disturbanceset[v] = true
            fullvars[i] = v
        end
    end
    if check
        ikeys = keys(filter(!last, inputset))
        if !isempty(ikeys)
            throw(IONotFoundError("inputs", nameof(state.sys), ikeys))
        end
        dkeys = keys(filter(!last, disturbanceset))
        if !isempty(dkeys)
            throw(IONotFoundError("disturbance inputs", nameof(state.sys), ikeys))
        end
        okeys = keys(filter(!last, outputset))
        if !isempty(okeys)
            throw(IONotFoundError("outputs", nameof(state.sys), okeys))
        end
    end
    state, orig_inputs
end

"""
    (; A, B, C, D), simplified_sys, extras = linearize(sys, inputs, outputs;    t=0.0, op = Dict(), allow_input_derivatives = false, zero_dummy_der=false, kwargs...)
    (; A, B, C, D), extras                 = linearize(simplified_sys, lin_fun; t=0.0, op = Dict(), allow_input_derivatives = false, zero_dummy_der=false)

Linearize `sys` between `inputs` and `outputs`, both vectors of variables. Return a NamedTuple with the matrices of a linear statespace representation
on the form

```math
\\begin{aligned}
ẋ &= Ax + Bu\\\\
y &= Cx + Du
\\end{aligned}
```

The first signature automatically calls [`linearization_function`](@ref) internally,
while the second signature expects the outputs of [`linearization_function`](@ref) as input.

`op` denotes the operating point around which to linearize. If none is provided,
the default values of `sys` are used.

If `allow_input_derivatives = false`, an error will be thrown if input derivatives (``u̇``) appear as inputs in the linearized equations. If input derivatives are allowed, the returned `B` matrix will be of double width, corresponding to the input `[u; u̇]`.

`zero_dummy_der` can be set to automatically set the operating point to zero for all dummy derivatives.

The return value `extras` is a NamedTuple `(; x, p, t)` containing the result of the initialization problem that was solved to determine the operating point.

See also [`linearization_function`](@ref) which provides a lower-level interface, [`linearize_symbolic`](@ref) and [`ModelingToolkit.reorder_unknowns`](@ref).

See extended help for an example.

The implementation and notation follows that of
["Linear Analysis Approach for Modelica Models", Allain et al. 2009](https://ep.liu.se/ecp/043/075/ecp09430097.pdf)

# Extended help

This example builds the following feedback interconnection and linearizes it from the input of `F` to the output of `P`.

```

  r ┌─────┐       ┌─────┐     ┌─────┐
───►│     ├──────►│     │  u  │     │
    │  F  │       │  C  ├────►│  P  │ y
    └─────┘     ┌►│     │     │     ├─┬─►
                │ └─────┘     └─────┘ │
                │                     │
                └─────────────────────┘
```

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
function plant(; name)
    @variables x(t) = 1
    @variables u(t)=0 y(t)=0
    eqs = [D(x) ~ -x + u
           y ~ x]
    System(eqs, t; name = name)
end

function ref_filt(; name)
    @variables x(t)=0 y(t)=0
    @variables u(t)=0 [input = true]
    eqs = [D(x) ~ -2 * x + u
           y ~ x]
    System(eqs, t, name = name)
end

function controller(kp; name)
    @variables y(t)=0 r(t)=0 u(t)=0
    @parameters kp = kp
    eqs = [
        u ~ kp * (r - y),
    ]
    System(eqs, t; name = name)
end

@named f = ref_filt()
@named c = controller(1)
@named p = plant()

connections = [f.y ~ c.r # filtered reference to controller reference
               c.u ~ p.u # controller output to plant input
               p.y ~ c.y]

@named cl = System(connections, t, systems = [f, c, p])

lsys0, ssys = linearize(cl, [f.u], [p.x])
desired_order = [f.x, p.x]
lsys = ModelingToolkit.reorder_unknowns(lsys0, unknowns(ssys), desired_order)

@assert lsys.A == [-2 0; 1 -2]
@assert lsys.B == [1; 0;;]
@assert lsys.C == [0 1]
@assert lsys.D[] == 0

## Symbolic linearization
lsys_sym, _ = ModelingToolkit.linearize_symbolic(cl, [f.u], [p.x])

@assert substitute(lsys_sym.A, ModelingToolkit.defaults(cl)) == lsys.A
```
"""
function linearize(sys, lin_fun::LinearizationFunction; t = 0.0,
        op = Dict(), allow_input_derivatives = false,
        p = DiffEqBase.NullParameters())
    prob = LinearizationProblem(lin_fun, t)
    op = anydict(op)
    evaluate_varmap!(op, keys(op))
    for (k, v) in op
        v === nothing && continue
        if symbolic_type(v) != NotSymbolic() || is_array_of_symbolics(v)
            v = getu(prob, v)(prob)
        end
        if is_parameter(prob, Initial(k))
            setu(prob, Initial(k))(prob, v)
        else
            setu(prob, k)(prob, v)
        end
    end
    p = anydict(p)
    for (k, v) in p
        setu(prob, k)(prob, v)
    end
    return solve(prob; allow_input_derivatives)
end

function linearize(sys, inputs, outputs; op = Dict(), t = 0.0,
        allow_input_derivatives = false,
        zero_dummy_der = false,
        kwargs...)
    lin_fun,
    ssys = linearization_function(sys,
        inputs,
        outputs;
        zero_dummy_der,
        op, t,
        kwargs...)
    mats, extras = linearize(ssys, lin_fun; op, t, allow_input_derivatives)
    mats, ssys, extras
end

"""
    (; Ã, B̃, C̃, D̃) = similarity_transform(sys, T; unitary=false)

Perform a similarity transform `T : Tx̃ = x` on linear system represented by matrices in NamedTuple `sys` such that

```
Ã = T⁻¹AT
B̃ = T⁻¹ B
C̃ = CT
D̃ = D
```

If `unitary=true`, `T` is assumed unitary and the matrix adjoint is used instead of the inverse.
"""
function similarity_transform(sys::NamedTuple, T; unitary = false)
    if unitary
        A = T'sys.A * T
        B = T'sys.B
    else
        Tf = lu(T)
        A = Tf \ sys.A * T
        B = Tf \ sys.B
    end
    C = sys.C * T
    D = sys.D
    (; A, B, C, D)
end

"""
    reorder_unknowns(sys::NamedTuple, old, new)

Permute the state representation of `sys` obtained from [`linearize`](@ref) so that the state unknown is changed from `old` to `new`
Example:

```
lsys, ssys = linearize(pid, [reference.u, measurement.u], [ctr_output.u])
desired_order = [int.x, der.x] # Unknowns that are present in unknowns(ssys)
lsys = ModelingToolkit.reorder_unknowns(lsys, unknowns(ssys), desired_order)
```

See also [`ModelingToolkit.similarity_transform`](@ref)
"""
function reorder_unknowns(sys::NamedTuple, old, new)
    nx = length(old)
    length(new) == nx || error("old and new must have the same length")
    perm = [findfirst(isequal(n), old) for n in new]
    issorted(perm) && return sys # shortcut return, no reordering
    P = zeros(Int, nx, nx)
    for i in 1:nx # Build permutation matrix
        P[perm[i], i] = 1
    end
    similarity_transform(sys, P; unitary = true)
end
