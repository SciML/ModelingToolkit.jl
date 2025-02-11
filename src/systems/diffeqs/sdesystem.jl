"""
$(TYPEDEF)

A system of stochastic differential equations.

# Fields
$(FIELDS)

# Example

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(x) ~ σ*(y-x),
       D(y) ~ x*(ρ-z)-y,
       D(z) ~ x*y - β*z]

noiseeqs = [0.1*x,
            0.1*y,
            0.1*z]

@named de = SDESystem(eqs,noiseeqs,t,[x,y,z],[σ,ρ,β]; tspan = (0, 1000.0))
```
"""
struct SDESystem <: AbstractODESystem
    """
    A tag for the system. If two systems have the same tag, then they are
    structurally identical.
    """
    tag::UInt
    """The expressions defining the drift term."""
    eqs::Vector{Equation}
    """The expressions defining the diffusion term."""
    noiseeqs::AbstractArray
    """Independent variable."""
    iv::BasicSymbolic{Real}
    """Dependent variables. Must not contain the independent variable."""
    unknowns::Vector
    """Parameter variables. Must not contain the independent variable."""
    ps::Vector
    """Time span."""
    tspan::Union{NTuple{2, Any}, Nothing}
    """Array variables."""
    var_to_name::Any
    """Control parameters (some subset of `ps`)."""
    ctrls::Vector
    """Observed variables."""
    observed::Vector{Equation}
    """
    Time-derivative matrix. Note: this field will not be defined until
    [`calculate_tgrad`](@ref) is called on the system.
    """
    tgrad::RefValue
    """
    Jacobian matrix. Note: this field will not be defined until
    [`calculate_jacobian`](@ref) is called on the system.
    """
    jac::RefValue
    """
    Control Jacobian matrix. Note: this field will not be defined until
    [`calculate_control_jacobian`](@ref) is called on the system.
    """
    ctrl_jac::RefValue{Any}
    """
    Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact::RefValue
    """
    Note: this field will not be defined until
    [`generate_factorized_W`](@ref) is called on the system.
    """
    Wfact_t::RefValue
    """
    The name of the system.
    """
    name::Symbol
    """
    A description of the system.
    """
    description::String
    """
    The internal systems. These are required to have unique names.
    """
    systems::Vector{SDESystem}
    """
    The default values to use when initial conditions and/or
    parameters are not supplied in `ODEProblem`.
    """
    defaults::Dict
    """
    The guesses to use as the initial conditions for the
    initialization system.
    """
    guesses::Dict
    """
    The system for performing the initialization.
    """
    initializesystem::Union{Nothing, NonlinearSystem}
    """
    Extra equations to be enforced during the initialization sequence.
    """
    initialization_eqs::Vector{Equation}
    """
    Type of the system.
    """
    connector_type::Any
    """
    A `Vector{SymbolicContinuousCallback}` that model events.
    The integrator will use root finding to guarantee that it steps at each zero crossing.
    """
    continuous_events::Vector{SymbolicContinuousCallback}
    """
    A `Vector{SymbolicDiscreteCallback}` that models events. Symbolic
    analog to `SciMLBase.DiscreteCallback` that executes an affect when a given condition is
    true at the end of an integration step.
    """
    discrete_events::Vector{SymbolicDiscreteCallback}
    """
    Topologically sorted parameter dependency equations, where all symbols are parameters and
    the LHS is a single parameter.
    """
    parameter_dependencies::Vector{Equation}
    """
    Mapping of conditions which should be true throughout the solution process to corresponding error
    messages. These will be added to the equations when calling `debug_system`.
    """
    assertions::Dict{BasicSymbolic, String}
    """
    Metadata for the system, to be used by downstream packages.
    """
    metadata::Any
    """
    Metadata for MTK GUI.
    """
    gui_metadata::Union{Nothing, GUIMetadata}
    """
    If a model `sys` is complete, then `sys.x` no longer performs namespacing.
    """
    complete::Bool
    """
    Cached data for fast symbolic indexing.
    """
    index_cache::Union{Nothing, IndexCache}
    """
    The hierarchical parent system before simplification.
    """
    parent::Any
    """
    Signal for whether the noise equations should be treated as a scalar process. This should only
    be `true` when `noiseeqs isa Vector`. 
    """
    is_scalar_noise::Bool
    """
    A boolean indicating if the given `ODESystem` represents a system of DDEs.
    """
    is_dde::Bool
    isscheduled::Bool

    function SDESystem(tag, deqs, neqs, iv, dvs, ps, tspan, var_to_name, ctrls, observed,
            tgrad, jac, ctrl_jac, Wfact, Wfact_t, name, description, systems, defaults,
            guesses, initializesystem, initialization_eqs, connector_type,
            cevents, devents, parameter_dependencies, assertions = Dict{
                BasicSymbolic, Nothing},
            metadata = nothing, gui_metadata = nothing,
            complete = false, index_cache = nothing, parent = nothing, is_scalar_noise = false,
            is_dde = false,
            isscheduled = false;
            checks::Union{Bool, Int} = true)
        if checks == true || (checks & CheckComponents) > 0
            check_independent_variables([iv])
            check_variables(dvs, iv)
            check_parameters(ps, iv)
            check_equations(deqs, iv)
            check_equations(neqs, dvs)
            if size(neqs, 1) != length(deqs)
                throw(ArgumentError("Noise equations ill-formed. Number of rows must match number of drift equations. size(neqs,1) = $(size(neqs,1)) != length(deqs) = $(length(deqs))"))
            end
            check_equations(equations(cevents), iv)
            if is_scalar_noise && neqs isa AbstractMatrix
                throw(ArgumentError("Noise equations ill-formed. Received a matrix of noise equations of size $(size(neqs)), but `is_scalar_noise` was set to `true`. Scalar noise is only compatible with an `AbstractVector` of noise equations."))
            end
        end
        if checks == true || (checks & CheckUnits) > 0
            u = __get_unit_type(dvs, ps, iv)
            check_units(u, deqs, neqs)
        end
        new(tag, deqs, neqs, iv, dvs, ps, tspan, var_to_name, ctrls, observed, tgrad, jac,
            ctrl_jac, Wfact, Wfact_t, name, description, systems,
            defaults, guesses, initializesystem, initialization_eqs, connector_type, cevents,
            devents, parameter_dependencies, assertions, metadata, gui_metadata, complete,
            index_cache, parent, is_scalar_noise, is_dde, isscheduled)
    end
end

function SDESystem(deqs::AbstractVector{<:Equation}, neqs::AbstractArray, iv, dvs, ps;
        controls = Num[],
        observed = Num[],
        systems = SDESystem[],
        tspan = nothing,
        default_u0 = Dict(),
        default_p = Dict(),
        defaults = _merge(Dict(default_u0), Dict(default_p)),
        guesses = Dict(),
        initializesystem = nothing,
        initialization_eqs = Equation[],
        name = nothing,
        description = "",
        connector_type = nothing,
        checks = true,
        continuous_events = nothing,
        discrete_events = nothing,
        parameter_dependencies = Equation[],
        assertions = Dict{BasicSymbolic, String}(),
        metadata = nothing,
        gui_metadata = nothing,
        complete = false,
        index_cache = nothing,
        parent = nothing,
        is_scalar_noise = false,
        is_dde = nothing)
    name === nothing &&
        throw(ArgumentError("The `name` keyword must be provided. Please consider using the `@named` macro"))
    iv′ = value(iv)
    dvs′ = value.(dvs)
    ps′ = value.(ps)
    ctrl′ = value.(controls)
    parameter_dependencies, ps′ = process_parameter_dependencies(
        parameter_dependencies, ps′)

    sysnames = nameof.(systems)
    if length(unique(sysnames)) != length(sysnames)
        throw(ArgumentError("System names must be unique."))
    end
    if !(isempty(default_u0) && isempty(default_p))
        Base.depwarn(
            "`default_u0` and `default_p` are deprecated. Use `defaults` instead.",
            :SDESystem, force = true)
    end

    defaults = Dict{Any, Any}(todict(defaults))
    guesses = Dict{Any, Any}(todict(guesses))
    var_to_name = Dict()
    process_variables!(var_to_name, defaults, guesses, dvs′)
    process_variables!(var_to_name, defaults, guesses, ps′)
    process_variables!(
        var_to_name, defaults, guesses, [eq.lhs for eq in parameter_dependencies])
    process_variables!(
        var_to_name, defaults, guesses, [eq.rhs for eq in parameter_dependencies])
    defaults = Dict{Any, Any}(value(k) => value(v)
    for (k, v) in pairs(defaults) if v !== nothing)
    guesses = Dict{Any, Any}(value(k) => value(v)
    for (k, v) in pairs(guesses) if v !== nothing)

    isempty(observed) || collect_var_to_name!(var_to_name, (eq.lhs for eq in observed))

    tgrad = RefValue(EMPTY_TGRAD)
    jac = RefValue{Any}(EMPTY_JAC)
    ctrl_jac = RefValue{Any}(EMPTY_JAC)
    Wfact = RefValue(EMPTY_JAC)
    Wfact_t = RefValue(EMPTY_JAC)
    cont_callbacks = SymbolicContinuousCallbacks(continuous_events)
    disc_callbacks = SymbolicDiscreteCallbacks(discrete_events)
    if is_dde === nothing
        is_dde = _check_if_dde(deqs, iv′, systems)
    end
    assertions = Dict{BasicSymbolic, Any}(unwrap(k) => v for (k, v) in assertions)
    SDESystem(Threads.atomic_add!(SYSTEM_COUNT, UInt(1)),
        deqs, neqs, iv′, dvs′, ps′, tspan, var_to_name, ctrl′, observed, tgrad, jac,
        ctrl_jac, Wfact, Wfact_t, name, description, systems, defaults, guesses,
        initializesystem, initialization_eqs, connector_type,
        cont_callbacks, disc_callbacks, parameter_dependencies, assertions, metadata, gui_metadata,
        complete, index_cache, parent, is_scalar_noise, is_dde; checks = checks)
end

function SDESystem(sys::ODESystem, neqs; kwargs...)
    SDESystem(equations(sys), neqs, get_iv(sys), unknowns(sys), parameters(sys); kwargs...)
end

function SDESystem(eqs::Vector{Equation}, noiseeqs::AbstractArray, iv; kwargs...)
    diffvars, allunknowns, ps, eqs = process_equations(eqs, iv)

    for eq in get(kwargs, :parameter_dependencies, Equation[])
        collect_vars!(allunknowns, ps, eq, iv)
    end

    for ssys in get(kwargs, :systems, ODESystem[])
        collect_scoped_vars!(allunknowns, ps, ssys, iv)
    end

    for v in allunknowns
        isdelay(v, iv) || continue
        collect_vars!(allunknowns, ps, arguments(v)[1], iv)
    end

    new_ps = OrderedSet()
    for p in ps
        if iscall(p) && operation(p) === getindex
            par = arguments(p)[begin]
            if Symbolics.shape(Symbolics.unwrap(par)) !== Symbolics.Unknown() &&
               all(par[i] in ps for i in eachindex(par))
                push!(new_ps, par)
            else
                push!(new_ps, p)
            end
        else
            push!(new_ps, p)
        end
    end

    # validate noise equations
    noisedvs = OrderedSet()
    noiseps = OrderedSet()
    collect_vars!(noisedvs, noiseps, noiseeqs, iv)
    for dv in noisedvs
        dv ∈ allunknowns ||
            throw(ArgumentError("Variable $dv in noise equations is not an unknown of the system."))
    end
    algevars = setdiff(allunknowns, diffvars)

    return SDESystem(eqs, noiseeqs, iv, Iterators.flatten((diffvars, algevars)),
        [ps; collect(noiseps)]; kwargs...)
end

function SDESystem(eq::Equation, noiseeqs::AbstractArray, args...; kwargs...)
    SDESystem([eq], noiseeqs, args...; kwargs...)
end
function SDESystem(eq::Equation, noiseeq, args...; kwargs...)
    SDESystem([eq], [noiseeq], args...; kwargs...)
end

function Base.:(==)(sys1::SDESystem, sys2::SDESystem)
    sys1 === sys2 && return true
    iv1 = get_iv(sys1)
    iv2 = get_iv(sys2)
    isequal(iv1, iv2) &&
        isequal(nameof(sys1), nameof(sys2)) &&
        _eq_unordered(get_eqs(sys1), get_eqs(sys2)) &&
        _eq_unordered(get_noiseeqs(sys1), get_noiseeqs(sys2)) &&
        isequal(get_is_scalar_noise(sys1), get_is_scalar_noise(sys2)) &&
        _eq_unordered(get_unknowns(sys1), get_unknowns(sys2)) &&
        _eq_unordered(get_ps(sys1), get_ps(sys2)) &&
        _eq_unordered(continuous_events(sys1), continuous_events(sys2)) &&
        _eq_unordered(discrete_events(sys1), discrete_events(sys2)) &&
        all(s1 == s2 for (s1, s2) in zip(get_systems(sys1), get_systems(sys2)))
end

"""
    function ODESystem(sys::SDESystem)

Convert an `SDESystem` to the equivalent `ODESystem` using `@brownian` variables instead
of noise equations. The returned system will not be `iscomplete` and will not have an
index cache, regardless of `iscomplete(sys)`.
"""
function ODESystem(sys::SDESystem)
    neqs = get_noiseeqs(sys)
    eqs = equations(sys)
    is_scalar_noise = get_is_scalar_noise(sys)
    nbrownian = if is_scalar_noise
        length(neqs)
    else
        size(neqs, 2)
    end
    brownvars = map(1:nbrownian) do i
        name = gensym(Symbol(:brown_, i))
        only(@brownian $name)
    end
    if is_scalar_noise
        brownterms = reduce(+, neqs .* brownvars; init = 0)
        neweqs = map(eqs) do eq
            eq.lhs ~ eq.rhs + brownterms
        end
    else
        if neqs isa AbstractVector
            neqs = reshape(neqs, (length(neqs), 1))
        end
        brownterms = neqs * brownvars
        neweqs = map(eqs, brownterms) do eq, brown
            eq.lhs ~ eq.rhs + brown
        end
    end
    newsys = ODESystem(neweqs, get_iv(sys), unknowns(sys), parameters(sys);
        parameter_dependencies = parameter_dependencies(sys), defaults = defaults(sys),
        continuous_events = continuous_events(sys), discrete_events = discrete_events(sys),
        assertions = assertions(sys),
        name = nameof(sys), description = description(sys), metadata = get_metadata(sys))
    @set newsys.parent = sys
end

function __num_isdiag_noise(mat)
    for i in axes(mat, 1)
        nnz = 0
        for j in axes(mat, 2)
            if !isequal(mat[i, j], 0)
                nnz += 1
            end
        end
        if nnz > 1
            return (false)
        end
    end
    true
end
function __get_num_diag_noise(mat)
    map(axes(mat, 1)) do i
        for j in axes(mat, 2)
            mij = mat[i, j]
            if !isequal(mij, 0)
                return mij
            end
        end
        0
    end
end

function generate_diffusion_function(sys::SDESystem, dvs = unknowns(sys),
        ps = parameters(sys); isdde = false, kwargs...)
    eqs = get_noiseeqs(sys)
    p = reorder_parameters(sys, ps)
    return build_function_wrapper(sys, eqs, dvs, p..., get_iv(sys); kwargs...)
end

"""
$(TYPEDSIGNATURES)

Choose correction_factor=-1//2 (1//2) to convert Ito -> Stratonovich (Stratonovich->Ito).
"""
function stochastic_integral_transform(sys::SDESystem, correction_factor)
    name = nameof(sys)
    # use the general interface
    if typeof(get_noiseeqs(sys)) <: Vector
        eqs = vcat([equations(sys)[i].lhs ~ get_noiseeqs(sys)[i]
                    for i in eachindex(unknowns(sys))]...)
        de = ODESystem(eqs, get_iv(sys), unknowns(sys), parameters(sys), name = name,
            checks = false)

        jac = calculate_jacobian(de, sparse = false, simplify = false)
        ∇σσ′ = simplify.(jac * get_noiseeqs(sys))

        deqs = vcat([equations(sys)[i].lhs ~ equations(sys)[i].rhs +
                                             correction_factor * ∇σσ′[i]
                     for i in eachindex(unknowns(sys))]...)
    else
        dimunknowns, m = size(get_noiseeqs(sys))
        eqs = vcat([equations(sys)[i].lhs ~ get_noiseeqs(sys)[i]
                    for i in eachindex(unknowns(sys))]...)
        de = ODESystem(eqs, get_iv(sys), unknowns(sys), parameters(sys), name = name,
            checks = false)

        jac = calculate_jacobian(de, sparse = false, simplify = false)
        ∇σσ′ = simplify.(jac * get_noiseeqs(sys)[:, 1])
        for k in 2:m
            eqs = vcat([equations(sys)[i].lhs ~ get_noiseeqs(sys)[Int(i +
                                                                      (k - 1) *
                                                                      dimunknowns)]
                        for i in eachindex(unknowns(sys))]...)
            de = ODESystem(eqs, get_iv(sys), unknowns(sys), parameters(sys), name = name,
                checks = false)

            jac = calculate_jacobian(de, sparse = false, simplify = false)
            ∇σσ′ = ∇σσ′ + simplify.(jac * get_noiseeqs(sys)[:, k])
        end

        deqs = vcat([equations(sys)[i].lhs ~ equations(sys)[i].rhs +
                                             correction_factor * ∇σσ′[i]
                     for i in eachindex(unknowns(sys))]...)
    end

    SDESystem(deqs, get_noiseeqs(sys), get_iv(sys), unknowns(sys), parameters(sys),
        name = name, description = description(sys),
        parameter_dependencies = parameter_dependencies(sys), checks = false)
end

"""
$(TYPEDSIGNATURES)

Measure transformation method that allows for a reduction in the variance of an estimator `Exp(g(X_t))`.
Input:  Original SDE system and symbolic function `u(t,x)` with scalar output that
        defines the adjustable parameters `d` in the Girsanov transformation. Optional: initial
        condition for `θ0`.
Output: Modified SDESystem with additional component `θ_t` and initial value `θ0`, as well as
        the weight `θ_t/θ0` as observed equation, such that the estimator `Exp(g(X_t)θ_t/θ0)`
        has a smaller variance.

Reference:
Kloeden, P. E., Platen, E., & Schurz, H. (2012). Numerical solution of SDE through computer
experiments. Springer Science & Business Media.

# Example

```julia
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D

@parameters α β
@variables x(t) y(t) z(t)

eqs = [D(x) ~ α*x]
noiseeqs = [β*x]

@named de = SDESystem(eqs,noiseeqs,t,[x],[α,β])

# define u (user choice)
u = x
θ0 = 0.1
g(x) = x[1]^2
demod = ModelingToolkit.Girsanov_transform(de, u; θ0=0.1)

u0modmap = [
    x => x0
]

parammap = [
    α => 1.5,
    β => 1.0
]

probmod = SDEProblem(complete(demod),u0modmap,(0.0,1.0),parammap)
ensemble_probmod = EnsembleProblem(probmod;
          output_func = (sol,i) -> (g(sol[x,end])*sol[demod.weight,end],false),
          )

simmod = solve(ensemble_probmod,EM(),dt=dt,trajectories=numtraj)
```

"""
function Girsanov_transform(sys::SDESystem, u; θ0 = 1.0)
    name = nameof(sys)

    # register new variable θ corresponding to 1D correction process θ(t)
    t = get_iv(sys)
    D = Differential(t)
    @variables θ(t), weight(t)

    # determine the adjustable parameters `d` given `u`
    # gradient of u with respect to unknowns
    grad = Symbolics.gradient(u, unknowns(sys))

    noiseeqs = get_noiseeqs(sys)
    if noiseeqs isa Vector
        d = simplify.(-(noiseeqs .* grad) / u)
        drift_correction = noiseeqs .* d
    else
        d = simplify.(-noiseeqs * grad / u)
        drift_correction = noiseeqs * d
    end

    # transformation adds additional unknowns θ: newX = (X,θ)
    # drift function for unknowns is modified
    # θ has zero drift
    deqs = vcat([equations(sys)[i].lhs ~ equations(sys)[i].rhs - drift_correction[i]
                 for i in eachindex(unknowns(sys))]...)
    deqsθ = D(θ) ~ 0
    push!(deqs, deqsθ)

    # diffusion matrix is of size d x m (d unknowns, m noise), with diagonal noise represented as a d-dimensional vector
    # for diagonal noise processes with m>1, the noise process will become non-diagonal; extra unknown component but no new noise process.
    # new diffusion matrix is of size d+1 x M
    # diffusion for state is unchanged

    noiseqsθ = θ * d

    if noiseeqs isa Vector
        m = size(noiseeqs)
        if m == 1
            push!(noiseeqs, noiseqsθ)
        else
            noiseeqs = [Array(Diagonal(noiseeqs)); noiseqsθ']
        end
    else
        noiseeqs = [Array(noiseeqs); noiseqsθ']
    end

    unknown_vars = [unknowns(sys); θ]

    # return modified SDE System
    SDESystem(deqs, noiseeqs, get_iv(sys), unknown_vars, parameters(sys);
        defaults = Dict(θ => θ0), observed = [weight ~ θ / θ0],
        name = name, description = description(sys),
        parameter_dependencies = parameter_dependencies(sys),
        checks = false)
end

function DiffEqBase.SDEFunction{iip, specialize}(sys::SDESystem, dvs = unknowns(sys),
        ps = parameters(sys),
        u0 = nothing;
        version = nothing, tgrad = false, sparse = false,
        jac = false, Wfact = false, eval_expression = false,
        eval_module = @__MODULE__,
        checkbounds = false, initialization_data = nothing,
        kwargs...) where {iip, specialize}
    if !iscomplete(sys)
        error("A completed `SDESystem` is required. Call `complete` or `structural_simplify` on the system before creating an `SDEFunction`")
    end
    dvs = scalarize.(dvs)

    f_gen = generate_function(sys, dvs, ps; expression = Val{true}, kwargs...)
    f_oop, f_iip = eval_or_rgf.(f_gen; eval_expression, eval_module)
    g_gen = generate_diffusion_function(sys, dvs, ps; expression = Val{true},
        kwargs...)
    g_oop, g_iip = eval_or_rgf.(g_gen; eval_expression, eval_module)

    f(u, p, t) = f_oop(u, p, t)
    f(du, u, p, t) = f_iip(du, u, p, t)
    g(u, p, t) = g_oop(u, p, t)
    g(du, u, p, t) = g_iip(du, u, p, t)

    if tgrad
        tgrad_gen = generate_tgrad(sys, dvs, ps; expression = Val{true},
            kwargs...)
        tgrad_oop, tgrad_iip = eval_or_rgf.(tgrad_gen; eval_expression, eval_module)

        _tgrad(u, p, t) = tgrad_oop(u, p, t)
        _tgrad(J, u, p, t) = tgrad_iip(J, u, p, t)
    else
        _tgrad = nothing
    end

    if jac
        jac_gen = generate_jacobian(sys, dvs, ps; expression = Val{true},
            sparse = sparse, kwargs...)
        jac_oop, jac_iip = eval_or_rgf.(jac_gen; eval_expression, eval_module)

        _jac(u, p, t) = jac_oop(u, p, t)
        _jac(J, u, p, t) = jac_iip(J, u, p, t)
    else
        _jac = nothing
    end

    if Wfact
        tmp_Wfact, tmp_Wfact_t = generate_factorized_W(sys, dvs, ps, true;
            expression = Val{true}, kwargs...)
        Wfact_oop, Wfact_iip = eval_or_rgf.(tmp_Wfact; eval_expression, eval_module)
        Wfact_oop_t, Wfact_iip_t = eval_or_rgf.(tmp_Wfact_t; eval_expression, eval_module)

        _Wfact(u, p, dtgamma, t) = Wfact_oop(u, p, dtgamma, t)
        _Wfact(W, u, p, dtgamma, t) = Wfact_iip(W, u, p, dtgamma, t)
        _Wfact_t(u, p, dtgamma, t) = Wfact_oop_t(u, p, dtgamma, t)
        _Wfact_t(W, u, p, dtgamma, t) = Wfact_iip_t(W, u, p, dtgamma, t)
    else
        _Wfact, _Wfact_t = nothing, nothing
    end

    M = calculate_massmatrix(sys)
    _M = (u0 === nothing || M == I) ? M : ArrayInterface.restructure(u0 .* u0', M)

    observedfun = ObservedFunctionCache(
        sys; eval_expression, eval_module, checkbounds = get(kwargs, :checkbounds, false))

    SDEFunction{iip, specialize}(f, g;
        sys = sys,
        jac = _jac === nothing ? nothing : _jac,
        tgrad = _tgrad === nothing ? nothing : _tgrad,
        Wfact = _Wfact === nothing ? nothing : _Wfact,
        Wfact_t = _Wfact_t === nothing ? nothing : _Wfact_t,
        mass_matrix = _M, initialization_data,
        observed = observedfun)
end

"""
```julia
DiffEqBase.SDEFunction{iip}(sys::SDESystem, dvs = sys.unknowns, ps = sys.ps;
                            version = nothing, tgrad = false, sparse = false,
                            jac = false, Wfact = false, kwargs...) where {iip}
```

Create an `SDEFunction` from the [`SDESystem`](@ref). The arguments `dvs` and `ps`
are used to set the order of the dependent variable and parameter vectors,
respectively.
"""
function DiffEqBase.SDEFunction(sys::SDESystem, args...; kwargs...)
    SDEFunction{true}(sys, args...; kwargs...)
end

function DiffEqBase.SDEFunction{true}(sys::SDESystem, args...;
        kwargs...)
    SDEFunction{true, SciMLBase.AutoSpecialize}(sys, args...; kwargs...)
end

function DiffEqBase.SDEFunction{false}(sys::SDESystem, args...;
        kwargs...)
    SDEFunction{false, SciMLBase.FullSpecialize}(sys, args...; kwargs...)
end

"""
```julia
DiffEqBase.SDEFunctionExpr{iip}(sys::AbstractODESystem, dvs = unknowns(sys),
                                ps = parameters(sys);
                                version = nothing, tgrad = false,
                                jac = false, Wfact = false,
                                skipzeros = true, fillzeros = true,
                                sparse = false,
                                kwargs...) where {iip}
```

Create a Julia expression for an `SDEFunction` from the [`SDESystem`](@ref).
The arguments `dvs` and `ps` are used to set the order of the dependent
variable and parameter vectors, respectively.
"""
struct SDEFunctionExpr{iip} end

function SDEFunctionExpr{iip}(sys::SDESystem, dvs = unknowns(sys),
        ps = parameters(sys), u0 = nothing;
        version = nothing, tgrad = false,
        jac = false, Wfact = false,
        sparse = false, linenumbers = false,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed `SDESystem` is required. Call `complete` or `structural_simplify` on the system before creating an `SDEFunctionExpr`")
    end
    idx = iip ? 2 : 1
    f = generate_function(sys, dvs, ps; expression = Val{true}, kwargs...)[idx]
    g = generate_diffusion_function(sys, dvs, ps; expression = Val{true}, kwargs...)[idx]
    if tgrad
        _tgrad = generate_tgrad(sys, dvs, ps; expression = Val{true}, kwargs...)[idx]
    else
        _tgrad = :nothing
    end

    if jac
        _jac = generate_jacobian(sys, dvs, ps; sparse = sparse, expression = Val{true},
            kwargs...)[idx]
    else
        _jac = :nothing
    end

    if Wfact
        tmp_Wfact, tmp_Wfact_t = generate_factorized_W(
            sys, dvs, ps; expression = Val{true},
            kwargs...)
        _Wfact = tmp_Wfact[idx]
        _Wfact_t = tmp_Wfact_t[idx]
    else
        _Wfact, _Wfact_t = :nothing, :nothing
    end

    M = calculate_massmatrix(sys)

    _M = (u0 === nothing || M == I) ? M : ArrayInterface.restructure(u0 .* u0', M)

    ex = quote
        f = $f
        g = $g
        tgrad = $_tgrad
        jac = $_jac
        Wfact = $_Wfact
        Wfact_t = $_Wfact_t
        M = $_M
        SDEFunction{$iip}(f, g,
            jac = jac,
            tgrad = tgrad,
            Wfact = Wfact,
            Wfact_t = Wfact_t,
            mass_matrix = M)
    end
    !linenumbers ? Base.remove_linenums!(ex) : ex
end

function SDEFunctionExpr(sys::SDESystem, args...; kwargs...)
    SDEFunctionExpr{true}(sys, args...; kwargs...)
end

function DiffEqBase.SDEProblem{iip, specialize}(
        sys::SDESystem, u0map = [], tspan = get_tspan(sys),
        parammap = DiffEqBase.NullParameters();
        sparsenoise = nothing, check_length = true,
        callback = nothing, kwargs...) where {iip, specialize}
    if !iscomplete(sys)
        error("A completed `SDESystem` is required. Call `complete` or `structural_simplify` on the system before creating an `SDEProblem`")
    end

    f, u0, p = process_SciMLProblem(
        SDEFunction{iip, specialize}, sys, u0map, parammap; check_length,
        t = tspan === nothing ? nothing : tspan[1], kwargs...)
    cbs = process_events(sys; callback, kwargs...)
    sparsenoise === nothing && (sparsenoise = get(kwargs, :sparse, false))

    noiseeqs = get_noiseeqs(sys)
    is_scalar_noise = get_is_scalar_noise(sys)
    if noiseeqs isa AbstractVector
        noise_rate_prototype = nothing
        if is_scalar_noise
            noise = WienerProcess(0.0, 0.0, 0.0)
        else
            noise = nothing
        end
    elseif sparsenoise
        I, J, V = findnz(SparseArrays.sparse(noiseeqs))
        noise_rate_prototype = SparseArrays.sparse(I, J, zero(eltype(u0)))
        noise = nothing
    else
        noise_rate_prototype = zeros(eltype(u0), size(noiseeqs))
        noise = nothing
    end

    kwargs = filter_kwargs(kwargs)

    SDEProblem{iip}(f, u0, tspan, p; callback = cbs, noise,
        noise_rate_prototype = noise_rate_prototype, kwargs...)
end

function DiffEqBase.SDEProblem(sys::ODESystem, args...; kwargs...)
    if any(ModelingToolkit.isbrownian, unknowns(sys))
        error("SDESystem constructed by defining Brownian variables with @brownian must be simplified by calling `structural_simplify` before a SDEProblem can be constructed.")
    else
        error("Cannot construct SDEProblem from a normal ODESystem.")
    end
end

"""
```julia
DiffEqBase.SDEProblem{iip}(sys::SDESystem, u0map, tspan, p = parammap;
                           version = nothing, tgrad = false,
                           jac = false, Wfact = false,
                           checkbounds = false, sparse = false,
                           sparsenoise = sparse,
                           skipzeros = true, fillzeros = true,
                           linenumbers = true, parallel = SerialForm(),
                           kwargs...)
```

Generates an SDEProblem from an SDESystem and allows for automatically
symbolically calculating numerical enhancements.
"""
function DiffEqBase.SDEProblem(sys::SDESystem, args...; kwargs...)
    SDEProblem{true}(sys, args...; kwargs...)
end

function DiffEqBase.SDEProblem(sys::SDESystem,
        u0map::StaticArray,
        args...;
        kwargs...)
    SDEProblem{false, SciMLBase.FullSpecialize}(sys, u0map, args...; kwargs...)
end

function DiffEqBase.SDEProblem{true}(sys::SDESystem, args...; kwargs...)
    SDEProblem{true, SciMLBase.AutoSpecialize}(sys, args...; kwargs...)
end

function DiffEqBase.SDEProblem{false}(sys::SDESystem, args...; kwargs...)
    SDEProblem{false, SciMLBase.FullSpecialize}(sys, args...; kwargs...)
end

"""
```julia
DiffEqBase.SDEProblemExpr{iip}(sys::AbstractODESystem, u0map, tspan,
                               parammap = DiffEqBase.NullParameters();
                               version = nothing, tgrad = false,
                               jac = false, Wfact = false,
                               checkbounds = false, sparse = false,
                               linenumbers = true, parallel = SerialForm(),
                               kwargs...) where {iip}
```

Generates a Julia expression for constructing an ODEProblem from an
ODESystem and allows for automatically symbolically calculating
numerical enhancements.
"""
struct SDEProblemExpr{iip} end

function SDEProblemExpr{iip}(sys::SDESystem, u0map, tspan,
        parammap = DiffEqBase.NullParameters();
        sparsenoise = nothing, check_length = true,
        kwargs...) where {iip}
    if !iscomplete(sys)
        error("A completed `SDESystem` is required. Call `complete` or `structural_simplify` on the system before creating an `SDEProblemExpr`")
    end
    f, u0, p = process_SciMLProblem(
        SDEFunctionExpr{iip}, sys, u0map, parammap; check_length,
        kwargs...)
    linenumbers = get(kwargs, :linenumbers, true)
    sparsenoise === nothing && (sparsenoise = get(kwargs, :sparse, false))

    noiseeqs = get_noiseeqs(sys)
    is_scalar_noise = get_is_scalar_noise(sys)
    if noiseeqs isa AbstractVector
        noise_rate_prototype = nothing
        if is_scalar_noise
            noise = WienerProcess(0.0, 0.0, 0.0)
        else
            noise = nothing
        end
    elseif sparsenoise
        I, J, V = findnz(SparseArrays.sparse(noiseeqs))
        noise_rate_prototype = SparseArrays.sparse(I, J, zero(eltype(u0)))
        noise = nothing
    else
        T = u0 === nothing ? Float64 : eltype(u0)
        noise_rate_prototype = zeros(T, size(get_noiseeqs(sys)))
        noise = nothing
    end
    ex = quote
        f = $f
        u0 = $u0
        tspan = $tspan
        p = $p
        noise_rate_prototype = $noise_rate_prototype
        noise = $noise
        SDEProblem(
            f, u0, tspan, p; noise_rate_prototype = noise_rate_prototype, noise = noise,
            $(kwargs...))
    end
    !linenumbers ? Base.remove_linenums!(ex) : ex
end

function SDEProblemExpr(sys::SDESystem, args...; kwargs...)
    SDEProblemExpr{true}(sys, args...; kwargs...)
end
