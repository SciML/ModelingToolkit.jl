"""
    $(TYPEDEF)

Wraps an `ODESolution` and a time point (or vector of time points) `t`. When passed as
`op` to [`linearize`](@ref), an operating point is constructed from the values of
differential state variables and parameters of `sol` evaluated at `t`. Algebraic
variables are not set and will be determined by the initialization algorithm.

When `t` is an `AbstractVector`, [`linearize`](@ref) calls [`linearization_function`](@ref)
once and evaluates the linearization at each time point, returning vectors of matrices and
extras.

The `op` keyword argument provides additional operating-point values that are merged into
the solution-derived operating point at every time point (taking precedence). This is how
values for variables that are not present in `sol` are supplied — in particular the
parameters created by `loop_openings`, e.g.
`LinearizationOpPoint(sol, t; op = Dict(opened_signal => 0))`.

# Fields

$(TYPEDFIELDS)
"""
struct LinearizationOpPoint{S <: SciMLBase.AbstractODESolution, T, D <: AbstractDict}
    """
    The solution to extract operating point values from.
    """
    sol::S
    """
    The time point (or vector of time points) at which to evaluate the solution.
    """
    t::T
    """
    Additional operating-point values merged into the solution-derived operating point at
    each time point (takes precedence over solution-derived values).
    """
    op::D
end

function LinearizationOpPoint(sol::SciMLBase.AbstractODESolution, t; op = Dict{SymbolicT, SymbolicT}())
    return LinearizationOpPoint(sol, t, op)
end

function _build_op_from_solution(op::LinearizationOpPoint)
    sol_sys = MTKBase.indp_to_system(op.sol)
    eqs = equations(sol_sys)
    sts = unknowns(sol_sys)
    u = op.sol(op.t)
    result = Dict{SymbolicT, SymbolicT}()
    for (i, eq) in enumerate(eqs)
        isdiffeq(eq) || continue
        result[sts[i]] = u[i]
    end
    for p in parameters(sol_sys)
        result[p] = getp(op.sol, p)(op.sol)
    end
    for (k, v) in op.op
        result[unwrap(k)] = v
    end
    return result
end

function _build_op_from_solution(op::LinearizationOpPoint{S, <:AbstractVector}) where {S <: SciMLBase.AbstractODESolution}
    sol_sys = MTKBase.indp_to_system(op.sol)
    eqs = equations(sol_sys)
    sts = unknowns(sol_sys)
    # Find differential equation indices and extract parameters once — both are
    # time-independent, so we only do this work once regardless of how many time
    # points are requested.
    diff_idxs = findall(isdiffeq, eqs)
    param_vals = Dict{SymbolicT, SymbolicT}()
    for p in parameters(sol_sys)
        param_vals[p] = getp(op.sol, p)(op.sol)
    end
    # Interpolate once per time point to build the per-point operating-point dict.
    extra_op = Dict{SymbolicT, SymbolicT}(unwrap(k) => v for (k, v) in op.op)
    return map(op.t) do ti
        u = op.sol(ti)
        result = copy(param_vals)
        for i in diff_idxs
            result[sts[i]] = u[i]
        end
        merge!(result, extra_op)
        result
    end
end

function _linearization_wrap_odeproblem_f(@nospecialize(prob::ODEProblem), ::Type{T}) where {T}
    f = SciMLBase.Void{Any}(prob.f.f)
    u0 = T.(prob.u0)
    t = T(prob.tspan[1])
    f = SciMLBase.wrapfun_iip(f, (u0, u0, prob.p, t))
    odef = remake(prob.f; f = f)
    return remake(prob; f = odef)
end

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
  - `ignore_system_initial_conditions`: Whether to ignore `initial_conditions(sys)` and only use `op`.
  - `kwargs`: Are passed on to `find_solvables!`

See also [`linearize`](@ref) which provides a higher-level interface.
"""
function linearization_function(
        sys::AbstractSystem, inputs,
        outputs; simplify = false,
        initialize = true,
        initializealg = nothing,
        initialization_abstol = 1.0e-5,
        initialization_reltol = 1.0e-3,
        op = Dict{SymbolicT, SymbolicT}(),
        p = DiffEqBase.NullParameters(),
        zero_dummy_der = false,
        initialization_solver_alg = nothing,
        autodiff = AutoForwardDiff(),
        eval_expression = false, eval_module = @__MODULE__,
        warn_initialize_determined = true,
        guesses = Dict{SymbolicT, SymbolicT}(),
        warn_empty_op = true,
        missing_guess_value = MTKBase.default_missing_guess_value(),
        t = 0.0,
        ignore_system_initial_conditions = false,
        loop_opening_params = SymbolicT[],
        kwargs...
    )
    op = Dict(op)
    if isempty(op) && warn_empty_op
        @warn "An empty operating point was passed to `linearization_function`. An operating point containing the variables that will be changed in `linearize` should be provided. Disable this warning by passing `warn_empty_op = false`."
    end
    inputs isa AbstractVector || (inputs = [inputs])
    outputs isa AbstractVector || (outputs = [outputs])
    ssys = mtkcompile(sys; inputs, outputs, simplify, kwargs...)
    if ignore_system_initial_conditions
        ics = copy(initial_conditions(ssys))
        filter!(Base.Fix2(SU.hasmetadata, MTKBase.AnalysisVariable) ∘ first, ics)
        @set! ssys.initial_conditions = ics
    end
    diff_idxs, alge_idxs = eq_idxs(ssys)
    if zero_dummy_der
        dummyder = setdiff(unknowns(ssys), unknowns(sys))
        ics = initial_conditions(ssys)
        for x in dummyder
            ics[x] = Symbolics.COMMON_ZERO
        end
    end

    _inputs = SymbolicT[]
    _outputs = SymbolicT[]
    for x in inputs
        if SU.is_array_shape(SU.shape(x))
            append!(_inputs, vec(collect(x)::Array{SymbolicT})::Vector{SymbolicT})
        else
            push!(_inputs, x)
        end
    end
    for x in outputs
        if SU.is_array_shape(SU.shape(x))
            append!(_outputs, vec(collect(x)::Array{SymbolicT})::Vector{SymbolicT})
        else
            push!(_outputs, x)
        end
    end
    inputs = _inputs
    outputs = _outputs
    sys = ssys

    if initializealg === nothing
        initializealg = initialize ? OverrideInit() : NoInit()
    end

    prob = ODEProblem{true}(
        sys, merge(op, anydict(p)), (t, t); allow_incomplete = true,
        algebraic_only = true, guesses, missing_guess_value
    )
    initial_idxs_for_unknowns = ParameterIndex{SciMLStructures.Initials, Int}[]
    for v in unknowns(sys)
        push!(initial_idxs_for_unknowns, parameter_index(sys, Initial(v))::eltype(initial_idxs_for_unknowns))
    end
    u0 = state_values(prob)

    ps = parameters(sys)
    h = build_explicit_observed_function(
        sys, outputs,
        GeneratedFunctionOptions(; expression = Val{false}, eval_expression, eval_module)
    )
    h = SciMLBase.Void{Any}(h)

    initialization_kwargs = (;
        abstol = initialization_abstol, reltol = initialization_reltol,
        nlsolve_alg = initialization_solver_alg,
    )

    p = parameter_values(prob)
    t0 = current_time(prob)
    inputvals = [prob.ps[i] for i in inputs]

    if u0 === nothing
        T = typeof(t0)
    else
        T = promote_type(eltype(u0), typeof(t0))
    end
    prob = _linearization_wrap_odeproblem_f(prob, T)
    ct0 = DI.Constant(T(t0))
    u0T = if u0 === nothing
        u0
    else
        T.(u0)
    end
    h = SciMLBase.wrapfun_iip(h, (u0T, u0T, p, T(t0)))
    hp_fun = HPFun(h, setp_oop(sys, inputs))

    cu0T = DI.Constant(u0T)
    cp = DI.Constant(p)

    if u0 === nothing
        uf_jac = h_jac = pf_jac = nothing
        Tp = promote_type(p isa MTKParameters ? eltype(p.tunable) : eltype(p), typeof(t0))
        hp_jac = PreparedJacobian{true}(
            hp_fun, zeros(Tp, size(outputs)), autodiff, inputvals,
            cu0T, cp, DI.Constant(t0)
        )
    else
        uf_fun = UFFun(prob.f)

        uf_jac = PreparedJacobian{true}(
            uf_fun, similar(prob.u0, T), autodiff, u0T, cp, ct0
        )
        # observed function is a `GeneratedFunctionWrapper` with iip component
        h_jac = PreparedJacobian{true}(
            h, similar(prob.u0, T, size(outputs)), autodiff,
            u0T, cp, ct0
        )

        pf_fun = PFFun(prob.f, setp_oop(sys, inputs))
        pf_jac = PreparedJacobian{true}(
            pf_fun, similar(prob.u0, T), autodiff, inputvals,
            cu0T, cp, ct0
        )
        hp_jac = PreparedJacobian{true}(
            hp_fun, similar(prob.u0, T, size(outputs)), autodiff, inputvals,
            cu0T, cp, ct0
        )
    end

    input_getter = getsym(prob, inputs)

    lin_fun = LinearizationFunction(
        diff_idxs, alge_idxs, input_getter, length(inputs), length(unknowns(sys)),
        prob, h, u0 === nothing ? nothing : similar(u0, T), uf_jac, h_jac, pf_jac,
        hp_jac, initializealg, initialization_kwargs, initial_idxs_for_unknowns,
        collect(SymbolicT, loop_opening_params)
    )
    return lin_fun, sys
end

struct HPFun{F, S}
    fn::F
    setter::S
end

function (hpf::HPFun)(du, input, u, p, t)
    p = hpf.setter(p, input)
    hpf.fn(du, u, p, t)
    return du
end

struct UFFun{F}
    fn::F
end

function (uff::UFFun)(du, u, p, t)
    return SciMLBase.UJacobianWrapper(uff.fn, t, p)(du, u)
end

struct PFFun{F, S}
    fn::F
    setter::S
end

function (pff::PFFun)(du, input, u, p, t)
    p = pff.setter(p, input)
    return SciMLBase.ParamJacobianWrapper(pff.fn, t, u)(du, p)
end

"""
Return the set of indexes of differential equations and algebraic equations in the simplified system.
"""
function eq_idxs(sys::AbstractSystem)
    eqs = equations(sys)
    alge_idxs = findall(!isdiffeq, eqs)
    diff_idxs = setdiff(1:length(eqs), alge_idxs)

    return diff_idxs, alge_idxs
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
        prep, f, buf, autodiff
    )
end

function PreparedJacobian{false}(f, autodiff, args...)
    prep = DI.prepare_jacobian(f, autodiff, args...; strict = Val(false))
    return PreparedJacobian{false, typeof(prep), typeof(f), Nothing, typeof(autodiff)}(
        prep, f, nothing, autodiff
    )
end

function (pj::PreparedJacobian{true})(args...)
    return DI.jacobian(pj.f, pj.buf, pj.prep, pj.autodiff, args...)
end

function (pj::PreparedJacobian{false})(args...)
    return DI.jacobian(pj.f, pj.prep, pj.autodiff, args...)
end

"""
    $(TYPEDEF)

A callable struct which linearizes a system.

# Fields

$(TYPEDFIELDS)
"""
mutable struct LinearizationFunction{
        I, P <: ODEProblem,
        H, C, J1, J2, J3, J4, IA <: SciMLBase.DAEInitializationAlgorithm, IK,
    }
    """
    The indexes of differential equations in the linearized system.
    """
    const diff_idxs::Vector{Int}
    """
    The indexes of algebraic equations in the linearized system.
    """
    const alge_idxs::Vector{Int}
    """
    Getter function for parameters in the linearized system which represent
    input variables.
    """
    const inputs_getter::I
    """
    Number of input variables.
    """
    const num_inputs::Int
    """
    The number of unknowns in the linearized system.
    """
    const num_states::Int
    """
    The `ODEProblem` of the linearized system.
    """
    prob::P
    """
    A function which takes `(u, p, t)` and returns the outputs of the linearized system.
    """
    const h::H
    """
    Any required cache buffers.
    """
    const caches::C
    """
    `PreparedJacobian` for calculating jacobian of `prob.f` w.r.t. `u`
    """
    const uf_jac::J1
    """
    `PreparedJacobian` for calculating jacobian of `h` w.r.t. `u`
    """
    const h_jac::J2
    """
    `PreparedJacobian` for calculating jacobian of `prob.f` w.r.t. `p`
    """
    const pf_jac::J3
    """
    `PreparedJacobian` for calculating jacobian of `h` w.r.t. `p`
    """
    const hp_jac::J4
    """
    The initialization algorithm to use.
    """
    const initializealg::IA
    """
    Keyword arguments to be passed to `SciMLBase.get_initial_values`.
    """
    const initialize_kwargs::IK
    """
    Index of `Initial(x)` for every `x` in unknowns.
    """
    const initial_idxs_for_unknowns::Vector{ParameterIndex{SciMLStructures.Initials, Int}}
    """
    Variables turned into parameters by `loop_openings`. Their operating-point values are
    not implied by the rest of the system, so they must be provided explicitly in the `op`
    passed to `linearize`; otherwise an error is thrown.
    """
    const loop_opening_params::Vector{SymbolicT}
end

SymbolicIndexingInterface.symbolic_container(f::LinearizationFunction) = f.prob
SymbolicIndexingInterface.state_values(f::LinearizationFunction) = state_values(f.prob)
function SymbolicIndexingInterface.parameter_values(f::LinearizationFunction)
    return parameter_values(f.prob)
end
SymbolicIndexingInterface.current_time(f::LinearizationFunction) = current_time(f.prob)
function SymbolicIndexingInterface.set_state!(f::LinearizationFunction, val, idx)
    f.prob.u0[idx] = val
    return f.prob.p[f.initial_idxs_for_unknowns[idx]] = val
end

function Base.show(io::IO, mime::MIME"text/plain", lf::LinearizationFunction)
    printstyled(io, "LinearizationFunction"; bold = true, color = :blue)
    println(io, " which wraps:")
    return show(io, mime, lf.prob)
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
                setp(linfun, k)(newps, value(v))
            end
        end
        p = newps
    end

    fun = linfun.prob.f
    input_vals = linfun.inputs_getter(linfun.prob)
    if u !== nothing # Handle systems without unknowns
        linfun.num_states == length(u) ||
            error("Number of unknown variables ($(linfun.num_states)) does not match the number of input unknowns ($(length(u)))")
        integ_cache = (linfun.caches,)
        integ = MockIntegrator{true}(u, p, t, fun, integ_cache, nothing)
        u, p,
            success = SciMLBase.get_initial_values(
            linfun.prob, integ, fun, linfun.initializealg, Val(true);
            linfun.initialize_kwargs...
        )
        u = u::typeof(linfun.prob.u0)
        p = p::typeof(linfun.prob.p)
        success = success::Bool
        if !success
            error("Initialization algorithm $(linfun.initializealg) failed with `unknowns = $u` and `p = $p`.")
        end
        fg_xz = linfun.uf_jac(u, DI.Constant(p), DI.Constant(t))
        h_xz = linfun.h_jac(u, DI.Constant(p), DI.Constant(t))
        fg_u = linfun.pf_jac(
            input_vals,
            DI.Constant(u), DI.Constant(p), DI.Constant(t)
        )
    else
        linfun.num_states == 0 ||
            error("Number of unknown variables (0) does not match the expected number of unknowns ($(linfun.num_states))")
        fg_xz = zeros(0, 0)
        h_xz = fg_u = zeros(0, length(linfun.num_inputs))
    end
    h_u = linfun.hp_jac(
        input_vals,
        DI.Constant(u), DI.Constant(p), DI.Constant(t)
    )
    return (
        f_x = fg_xz[linfun.diff_idxs, linfun.diff_idxs],
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
        t,
    )
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
        u::U, p::P, t::T, f::F, cache::C, opts::O
    ) where {iip, U, P, T, F, C, O}
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
    return show(io, mime, prob.f.prob)
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
function SymbolicIndexingInterface.set_state!(p::LinearizationProblem, val, idx)
    return SymbolicIndexingInterface.set_state!(p.f, val, idx)
end


function Base.getindex(prob::LinearizationProblem, idx)
    return getu(prob, idx)(prob)
end

function Base.setindex!(prob::LinearizationProblem, val, idx)
    return setu(prob, idx)(prob, val)
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
        A = [
            f_x f_z
            gzgx * f_x gzgx * f_z
        ]
        B = [
            f_u
            gzgx * f_u
        ] # The cited paper has zeros in the bottom block, see derivation in https://github.com/SciML/ModelingToolkit.jl/pull/1691 for the correct formula

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

    return (; A, B, C, D), (; x, p, t)
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
function linearize_symbolic(
        sys::AbstractSystem, inputs,
        outputs; simplify = false, allow_input_derivatives = false,
        eval_expression = false, eval_module = @__MODULE__, split = true,
        kwargs...
    )
    # We cannot use `inline_linear_sccs` since it prevents symbolic AD
    reassemble_alg = MTKTearing.DefaultReassembleAlgorithm(; inline_linear_sccs = false)
    sys = mtkcompile(sys; inputs, outputs, simplify, split, reassemble_alg, kwargs...)
    check_symbolic_ad_allowed(sys)
    diff_idxs, alge_idxs = eq_idxs(sys)
    sts = unknowns(sys)
    t = get_iv(sys)
    ps = parameters(sys; initial_parameters = true)
    p = Tuple(reorder_parameters(sys, ps))

    fun_result = generate_rhs(sys, GeneratedFunctionOptions(; expression = Val{true}))
    fun_expr = fun_result isa Tuple ? fun_result[1] : fun_result
    fun = eval_or_rgf(fun_expr; eval_expression, eval_module)

    h = build_explicit_observed_function(
        sys, outputs,
        GeneratedFunctionOptions(; expression = Val{false}, eval_expression, eval_module)
    )
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
        gz = lu(Num.(g_z); check = false)
        issuccess(gz) ||
            error("g_z not invertible, this indicates that the DAE is of index > 1.")
        gzgx = -(gz \ Num.(g_x))
        A = [
            f_x f_z
            gzgx * f_x gzgx * f_z
        ]
        B = [
            f_u
            gzgx * f_u
        ] # The cited paper has zeros in the bottom block, see derivation in https://github.com/SciML/ModelingToolkit.jl/pull/1691 for the correct formula

        C = [h_x h_z]
        Bs = -(gz \ Num.(g_u)) # This equation differ from the cited paper, the paper is likely wrong since their equaiton leads to a dimension mismatch.
        if !iszero(Bs)
            if !allow_input_derivatives
                der_inds = findall(vec(any(!iszero, Bs, dims = 1)))
                error("Input derivatives appeared in expressions (-g_z\\g_u != 0), the following inputs appeared differentiated: $(ModelingToolkit.inputs(sys)[der_inds]). Call `linearize_symbolic` with keyword argument `allow_input_derivatives = true` to allow this and have the returned `B` matrix be of double width ($(2nu)), where the last $nu inputs are the derivatives of the first $nu inputs.")
            end
            B = [B [zeros(nx, nu); Bs]]
            D = [D zeros(ny, nu)]
        end
    end

    return (; A, B, C, D, f_x, f_z, g_x, g_z, f_u, g_u, h_x, h_z, h_u), sys
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
    return if maybe_namespace_issue
        println(
            io, """
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
            """
        )
    end
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
function linearize(
        sys, lin_fun::LinearizationFunction; t = 0.0,
        op = Dict(), allow_input_derivatives = false,
        p = DiffEqBase.NullParameters()
    )
    if op isa LinearizationOpPoint && op.t isa AbstractVector
        ops = _build_op_from_solution(op)
        results = map(zip(ops, op.t)) do (op_i, ti)
            linearize(sys, lin_fun; t = ti, op = op_i, allow_input_derivatives, p)
        end
        return first.(results), last.(results)
    end
    if op isa LinearizationOpPoint
        t = op.t
        op = _build_op_from_solution(op)
    end
    prob = LinearizationProblem(lin_fun, t)
    op = as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(op), COMMON_NOTHING)
    evaluate_varmap!(op, keys(op))
    _check_loop_opening_op(lin_fun.loop_opening_params, op)
    for (k, v) in op
        isequal(v, COMMON_NOTHING) && continue
        v = _resolve_op_value(prob, v)
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

# Resolve an operating-point value to something that can be passed to a setter.
# A numeric value stored in a `Dict{SymbolicT, SymbolicT}` is wrapped as a symbolic
# constant; unwrap it to a plain number/array. Routing a constant through `getu` would
# build a fresh `RuntimeGeneratedFunction` (with the value baked into the generated code)
# on every call, which dominates runtime when linearizing along a trajectory. Only
# genuinely symbolic values (expressions referencing the problem state) need `getu`.
function _resolve_op_value(prob, v)
    if SU.isconst(v)
        return SU.unwrap_const(v)
    elseif symbolic_type(v) != NotSymbolic() || is_array_of_symbolics(v)
        return getu(prob, v)(prob)
    else
        return v
    end
end

# Variables turned into parameters by `loop_openings` have no operating-point value implied
# by the rest of the system, so they must be supplied explicitly in `op`. Error (rather than
# silently using a stale/default value) if any of them is missing.
function _check_loop_opening_op(loop_opening_params, op)
    isempty(loop_opening_params) && return nothing
    missing_params = SymbolicT[]
    for p in loop_opening_params
        v = get(op, p, COMMON_NOTHING)
        isequal(v, COMMON_NOTHING) && push!(missing_params, p)
    end
    isempty(missing_params) && return nothing
    params_str = join(string.(missing_params), ", ")
    error(
        """
        The operating point does not provide values for the loop-opening parameter(s): \
        $(params_str). When `loop_openings` is used, the opened signals become \
        parameters whose operating-point values are not implied by the rest of the \
        system, so they must be provided explicitly in `op` (e.g. set to zero). When \
        linearizing along a trajectory with `LinearizationOpPoint`, pass them via its \
        `op` keyword argument: `LinearizationOpPoint(sol, t; op = Dict(signal => value))`.
        """
    )
end

function __linearize_multiple_op_barrier(ssys, lin_fun; ops, ts, allow_input_derivatives)
    T = eltype(lin_fun.prob.u0)
    results = @NamedTuple{A::Matrix{T}, B::Matrix{T}, C::Matrix{T}, D::Matrix{T}}[]
    xpts = @NamedTuple{x::typeof(lin_fun.prob.u0), p::typeof(lin_fun.prob.p), t::typeof(lin_fun.prob.tspan[1])}[]
    isempty(ops) && return results, xpts

    # Build the linearization problem once and reuse it across all time points, mutating
    # only the operating point and `.t`. This avoids reconstructing the problem and
    # rebuilding the symbolic setters on every iteration. The set of operating-point keys
    # is identical across time points (only the values change), so resolve the first op to
    # obtain the keys and build the setters once.
    prob = LinearizationProblem(lin_fun, ts[1])
    op1 = as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(ops[1]), COMMON_NOTHING)
    evaluate_varmap!(op1, keys(op1))
    # The op keys are identical across time points, so checking the first one suffices.
    _check_loop_opening_op(lin_fun.loop_opening_params, op1)
    op_keys = collect(keys(op1))
    setters = map(op_keys) do k
        is_parameter(prob, Initial(k)) ? setu(prob, Initial(k)) : setu(prob, k)
    end

    for (op, t) in zip(ops, ts)
        op = as_atomic_dict_with_defaults(Dict{SymbolicT, SymbolicT}(op), COMMON_NOTHING)
        evaluate_varmap!(op, keys(op))
        for (setter, k) in zip(setters, op_keys)
            v = get(op, k, COMMON_NOTHING)
            isequal(v, COMMON_NOTHING) && continue
            setter(prob, _resolve_op_value(prob, v))
        end
        prob.t = t
        res, xpt = solve(prob; allow_input_derivatives)::Tuple{eltype(results), eltype(xpts)}
        push!(results, res)
        push!(xpts, xpt)
    end
    return results, xpts
end

function linearize(
        sys, inputs, outputs; op = Dict(), t = 0.0,
        allow_input_derivatives = false,
        zero_dummy_der = false,
        kwargs...
    )
    if op isa LinearizationOpPoint && op.t isa AbstractVector
        ops = _build_op_from_solution(op)
        ts = op.t
        # Build the linearization function once using the first operating point, then
        # reuse it for all subsequent time points — this avoids redundant `mtkcompile`
        # and Jacobian preparation work.
        lin_fun, ssys = linearization_function(
            sys, inputs, outputs;
            zero_dummy_der, op = ops[1], t = ts[1],
            ignore_system_initial_conditions = true, kwargs...
        )
        ress, ops = __linearize_multiple_op_barrier(ssys, lin_fun; ops, ts, allow_input_derivatives)
        return ress, ssys, ops
    end
    ignore_system_ics = false
    if op isa LinearizationOpPoint
        t = op.t
        op = _build_op_from_solution(op)
        ignore_system_ics = true
    end
    lin_fun,
        ssys = linearization_function(
        sys,
        inputs,
        outputs;
        zero_dummy_der,
        op, t,
        ignore_system_initial_conditions = ignore_system_ics,
        kwargs...
    )
    mats, extras = linearize(ssys, lin_fun; op, t, allow_input_derivatives)
    return mats, ssys, extras
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
    return (; A, B, C, D)
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
    return similarity_transform(sys, P; unitary = true)
end
