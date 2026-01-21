const REPEATED_SIMPLIFICATION_MESSAGE = "Structural simplification cannot be applied to a completed system. Double simplification is not allowed."

struct RepeatedStructuralSimplificationError <: Exception end

function Base.showerror(io::IO, e::RepeatedStructuralSimplificationError)
    return print(io, REPEATED_SIMPLIFICATION_MESSAGE)
end

function canonicalize_io(iovars, type::String)
    iobuffer = OrderedSet{SymbolicT}()
    arrsyms = AtomicArrayDict{OrderedSet{SymbolicT}}()
    for var in iovars
        sh = SU.shape(var)
        if SU.is_array_shape(sh)
            if sh isa SU.ShapeVecT
                union!(iobuffer, vec(collect(var)::Array{SymbolicT})::Vector{SymbolicT})
                continue
            end
            throw(
                ArgumentError(
                    """
                    All $(type)s must have known shape. Found $var with unknown shape.
                    """
                )
            )
        end
        arr, isarr = split_indexed_var(var)
        if isarr
            tmp = get!(OrderedSet{SymbolicT}, arrsyms, arr)
            push!(tmp, var)
        end
        push!(iobuffer, var)
    end

    for (k, v) in arrsyms
        if !symbolic_has_known_size(k)
            throw(
                ArgumentError(
                    """
                    All $(type)s must have known shape. Found $k with unknown shape.
                    """
                )
            )
        end
        if type != "output" && length(k) != length(v)
            throw(
                ArgumentError(
                    """
                    Part of an array variable cannot be made an $type. The entire array must be \
                    an $type. Found $k which has $(length(v)) elements out of $(length(k)) in \
                    the $(type)s. Either pass all scalarized elements in sorted order as $(type)s \
                    or simply pass $k as an $type.
                    """
                )
            )
        end
        if type != "output" && !isequal(vec(collect(k)::Array{SymbolicT})::Vector{SymbolicT}, collect(v))
            throw(
                ArgumentError(
                    """
                    Elements of scalarized array variables must be in sorted order in $(type)s. \
                    Either pass all scalarized elements in sorted order as $(type)s \
                    or simply pass $k as an $type.
                    """
                )
            )
        end
    end

    return iobuffer
end

"""
$(SIGNATURES)

Compile the given system into a form that ModelingToolkitBase can generate code for. Also
performs order reduction for ODEs and handles simple discrete/implicit-discrete systems.

# Keyword Arguments

+ `fully_determined=true` controls whether or not an error will be thrown if the number of equations don't match the number of inputs, outputs, and equations.
+ `inputs`, `outputs` and `disturbance_inputs` are passed as keyword arguments.` All inputs` get converted to parameters and are allowed to be unconnected, allowing models where `n_unknowns = n_equations - n_inputs`.
"""
function mtkcompile(
        sys::System; additional_passes = (),
        inputs = SymbolicT[], outputs = SymbolicT[],
        disturbance_inputs = SymbolicT[],
        split = true, kwargs...
    )
    isscheduled(sys) && throw(RepeatedStructuralSimplificationError())
    # Canonicalize types of arguments to prevent repeated compilation of inner methods
    inputs = canonicalize_io(unwrap_vars(inputs), "input")
    outputs = canonicalize_io(unwrap_vars(outputs), "output")
    disturbance_inputs = canonicalize_io(unwrap_vars(disturbance_inputs), "disturbance input")
    newsys = _mtkcompile(
        sys;
        inputs, outputs, disturbance_inputs, additional_passes,
        kwargs...
    )
    for pass in additional_passes
        newsys = pass(newsys)
    end
    @set! newsys.parent = complete(sys; split = false, flatten = false)
    newsys = complete(newsys; split)
    return newsys
end

function scalarized_vars(vars)
    scal = SymbolicT[]
    for var in vars
        if !SU.is_array_shape(SU.shape(var))
            push!(scal, var)
            continue
        end
        for i in SU.stable_eachindex(var)
            push!(scal, var[i])
        end
    end
    return scal
end

function _mtkcompile(sys::AbstractSystem; kwargs...)
    # TODO: convert noise_eqs to brownians for simplification
    if has_noise_eqs(sys) && get_noise_eqs(sys) !== nothing
        sys = noise_to_brownians(sys; names = :αₘₜₖ)
    end
    if !isempty(jumps(sys))
        return sys
    end
    if isempty(equations(sys)) && !is_time_dependent(sys) && !_iszero(cost(sys))
        return simplify_optimization_system(sys; kwargs...)::System
    end
    if !isempty(brownians(sys))
        return simplify_sde_system(sys; kwargs...)
    end
    return __mtkcompile(sys; kwargs...)
end

function __mtkcompile(
        sys::AbstractSystem;
        inputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        outputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        disturbance_inputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        fully_determined = true,
        kwargs...
    )
    sys = expand_connections(sys)
    sys = discrete_unknowns_to_parameters(sys)
    sys = discover_globalscoped(sys)
    flat_dvs = scalarized_vars(unknowns(sys))
    original_vars = Set{SymbolicT}(flat_dvs)
    eqs = flatten_equations(equations(sys))
    all_dvs = Set{SymbolicT}()
    for eq in eqs
        SU.search_variables!(all_dvs, eq; is_atomic = OperatorIsAtomic{Union{Initial, Pre}}())
    end
    _all_dvs = Set{SymbolicT}()
    for v in all_dvs
        if Symbolics.isarraysymbolic(v)
            for i in SU.stable_eachindex(v)
                push!(_all_dvs, v[i])
            end
        else
            push!(_all_dvs, v)
        end
    end
    all_dvs = _all_dvs
    filter!(all_dvs) do v
        v in original_vars || split_indexed_var(v)[1] in original_vars
    end

    new_binds = copy(parent(bindings(sys)))
    new_ics = copy(initial_conditions(sys))
    for unused in setdiff(original_vars, all_dvs)
        arr, isarr = split_indexed_var(unused)
        arr in all_dvs && continue
        delete!(new_binds, arr)
        delete!(new_ics, arr)
    end

    setdiff!(all_dvs, inputs, disturbance_inputs)
    if fully_determined === nothing
        fully_determined = false
    end
    if fully_determined && length(eqs) > length(all_dvs)
        throw(
            ExtraEquationsSystemException(
                """
                The system is unbalanced. There are $(length(eqs)) equations and \
                $(length(all_dvs)) unknowns.
                """
            )
        )
    elseif fully_determined && length(eqs) < length(all_dvs)
        throw(
            ExtraVariablesSystemException(
                """
                The system is unbalanced. There are $(length(eqs)) equations and \
                $(length(all_dvs)) unknowns. This may also be a high-index DAE, which \
                ModelingToolkitBase.jl cannot handle. Consider using ModelingToolkit.jl to \
                simplify this system.
                """
            )
        )
    end

    flat_dvs = collect(all_dvs)
    has_derivatives = any(hasderiv, eqs)
    has_shifts = any(hasshift, eqs)
    if has_derivatives && has_shifts
        throw(
            HybridSystemNotSupportedException(
                """
                ModelingToolkitBase.jl cannot simplify systems with both `Shift` and \
                `Differential` operators.
                """
            )
        )
    end
    # Nonlinear system
    if !has_derivatives && !has_shifts
        obseqs = Equation[]
        get_trivial_observed_equations!(Equation[], eqs, obseqs, all_dvs, nothing)
        add_array_observed!(obseqs)
        obseqs = topsort_equations(obseqs, [eq.lhs for eq in obseqs])
        map!(eq -> Symbolics.COMMON_ZERO ~ (eq.rhs - eq.lhs), eqs, eqs)
        observables = Set{SymbolicT}()
        for eq in obseqs
            push!(observables, eq.lhs)
        end
        setdiff!(flat_dvs, observables)
        @set! sys.eqs = eqs
        @set! sys.unknowns = flat_dvs
        @set! sys.observed = obseqs
        return sys
    end
    iv = get_iv(sys)::SymbolicT
    total_sub = Dict{SymbolicT, SymbolicT}()
    subst = SU.Substituter{false}(total_sub, SU.default_substitute_filter)
    if has_derivatives
        D = Differential(iv)

        diffeq_idxs = isdiffeq.(eqs)
        diffeqs = eqs[diffeq_idxs]
        alg_eqs = eqs[.!diffeq_idxs]
        for i in eachindex(diffeqs)
            eq = diffeqs[i]
            var, order = Moshi.Match.@match eq.lhs begin
                BSImpl.Term(; f, args) && if f isa Differential end => (args[1], f.order::Int)
            end

            @assert order >= 1
            # Simple order reduction
            cur = var
            for i in 1:(order - 1)
                lhs = D(cur)
                rhs = default_toterm(lhs)
                push!(diffeqs, lhs ~ rhs)
                cur = rhs
            end

            diffeqs[i] = D(cur) ~ eq.rhs
        end

        obseqs = Equation[]
    else
        # The "most differentiated" variable in `x(k) ~ x(k - 1) + x(k - 2)` is `x(k)`.
        # To find how many times it is "differentiated", find the lowest shift.
        lowest_shift = Dict{SymbolicT, Int}()
        varsbuf = Set{SymbolicT}()
        for eq in eqs
            SU.search_variables!(varsbuf, eq.lhs; is_atomic = OperatorIsAtomic{Shift}())
            SU.search_variables!(varsbuf, eq.rhs; is_atomic = OperatorIsAtomic{Shift}())
        end
        for v in varsbuf
            Moshi.Match.@match v begin
                BSImpl.Term(; f, args) && if f isa Shift end => begin
                    if f.steps > 0
                        throw(
                            ArgumentError(
                                """
                                Positive shifts are disallowed in unsimplified equations. Found $v.
                                """
                            )
                        )
                    end
                    var = args[1]
                    lowest_shift[var] = min(get(lowest_shift, var, 0), f.steps)
                end
                _ => nothing
            end
        end

        # "differential" equations are ones with shifted variables on the LHS
        diffeq_idxs = falses(length(eqs))
        for i in eachindex(eqs)
            eq = eqs[i]
            if eq.lhs in all_dvs && !haskey(lowest_shift, eq.lhs)
                lowest_shift[eq.lhs] = 0
            end
            diffeq_idxs[i] = get(lowest_shift, eqs[i].lhs, typemax(Int)) <= 0
        end
        # They actually become observed.
        obseqs = eqs[diffeq_idxs]
        alg_eqs = eqs[.!diffeq_idxs]
        diffeqs = Equation[]
        for (var, order) in lowest_shift
            order = -order
            @assert order >= 0
            # A variable shifted back `order` times requires `order` elements of
            # history.
            for i in 1:order
                lhs = Shift(iv, 1)(default_toterm(Shift(iv, -i)(var)))
                rhs = default_toterm(Shift(iv, -i + 1)(var))
                push!(diffeqs, lhs ~ rhs)
                total_sub[Shift(iv, -i)(var)] = default_toterm(Shift(iv, -i)(var))
            end
        end

        _obseqs = topsort_equations(obseqs, collect(all_dvs); check = false)
        _algeqs = setdiff!(obseqs, _obseqs)
        for i in eachindex(_algeqs)
            _algeqs[i] = Symbolics.COMMON_ZERO ~ _algeqs[i].rhs - _algeqs[i].lhs
        end
        obseqs = _obseqs
        append!(alg_eqs, _algeqs)
    end

    # Substitute derivatives used in RHS of equations
    for eq in diffeqs
        total_sub[eq.lhs] = eq.rhs
    end
    for i in eachindex(diffeqs)
        eq = diffeqs[i]
        diffeqs[i] = eq.lhs ~ fixpoint_sub(eq.rhs, total_sub)
    end
    diffvars = SymbolicT[]
    # Store fixpoint subbed mapping
    for eq in diffeqs
        total_sub[eq.lhs] = eq.rhs
        push!(
            diffvars, Moshi.Match.@match eq.lhs begin
                BSImpl.Term(; args) => args[1]
            end
        )
    end
    get_trivial_observed_equations!(diffeqs, alg_eqs, obseqs, all_dvs, iv)
    add_array_observed!(obseqs)
    obseqs = topsort_equations(obseqs, [eq.lhs for eq in obseqs])
    for i in eachindex(alg_eqs)
        eq = alg_eqs[i]
        alg_eqs[i] = 0 ~ subst(eq.rhs - eq.lhs)
    end
    for i in eachindex(obseqs)
        eq = obseqs[i]
        obseqs[i] = eq.lhs ~ subst(eq.rhs)
    end
    alg_vars = setdiff!(flat_dvs, diffvars, [eq.lhs for eq in obseqs], inputs, disturbance_inputs)

    new_eqs = [diffeqs; alg_eqs]
    new_dvs = [diffvars; alg_vars]
    new_ps = [get_ps(sys); collect(inputs)]

    for eq in new_eqs
        if SU.query(eq.rhs) do v
                Moshi.Match.@match v begin
                    BSImpl.Term(; f) && if f isa Union{Differential, Shift} end => true
                    _ => false
                end
            end
            throw(
                ArgumentError(
                    """
                    ModelingToolkitBase.jl is unable to simplify such systems. Encountered \
                    derivative in RHS of equation $eq. Please consider using ModelingToolkit.jl \
                    for such systems.
                    """
                )
            )
        end
    end

    dummy_sub = Dict{SymbolicT, SymbolicT}()
    for eq in diffeqs
        dummy_sub[eq.lhs] = eq.rhs
    end
    var_sccs = [collect(eachindex(new_eqs))]
    schedule = Schedule(var_sccs, dummy_sub)

    @set! sys.eqs = new_eqs
    @set! sys.observed = obseqs
    @set! sys.unknowns = new_dvs
    @set! sys.ps = new_ps
    @set! sys.inputs = inputs
    @set! sys.outputs = outputs
    @set! sys.schedule = schedule
    @set! sys.isscheduled = true
    @set! sys.bindings = new_binds
    @set! sys.initial_conditions = new_ics
    return sys
end

"""
    $TYPEDSIGNATURES

For explicit algebraic equations in `algeqs`, find ones where the RHS is a function of
differential variables or other observed variables. These equations are removed from
`algeqs` and appended to `obseqs`. The process runs iteratively until a fixpoint is
reached.
"""
function get_trivial_observed_equations!(
        diffeqs::Vector{Equation}, algeqs::Vector{Equation},
        obseqs::Vector{Equation}, all_dvs::Set{SymbolicT},
        @nospecialize(iv::Union{SymbolicT, Nothing})
    )
    # Maximum number of times to loop over all algebraic equations
    maxiters = 100
    # Whether it's worth doing another loop, or we already reached a fixpoint
    active = true

    current_observed = Set{SymbolicT}()
    for eq in obseqs
        push!(current_observed, eq.lhs)
    end
    diffvars = Set{SymbolicT}()
    for eq in diffeqs
        push!(
            diffvars, Moshi.Match.@match eq.lhs begin
                BSImpl.Term(; f, args) && if f isa Union{Shift, Differential} end => args[1]
            end
        )
    end
    # Incidence information
    vars_in_each_algeq = Set{SymbolicT}[]
    sizehint!(vars_in_each_algeq, length(algeqs))
    for eq in algeqs
        buffer = Set{SymbolicT}()
        SU.search_variables!(buffer, eq.rhs)
        # We only care for variables
        intersect!(buffer, all_dvs)
        # If `eq.lhs` is only dependent on differential or other observed variables,
        # we can tear it. So we don't care about those either.
        setdiff!(buffer, diffvars)
        setdiff!(buffer, current_observed)
        if iv isa SymbolicT
            delete!(buffer, iv)
        end
        push!(vars_in_each_algeq, buffer)
    end
    # Algebraic equations that we still consider for elimination
    active_alg_eqs = trues(length(algeqs))
    # The number of equations we're considering for elimination
    candidate_eqs_count = length(algeqs)
    # Algebraic equations that we still consider algebraic
    alg_eqs_mask = trues(length(algeqs))
    # Observed variables added by this process
    new_observed_variables = Set{SymbolicT}()
    while active && maxiters > 0 && candidate_eqs_count > 0
        # We've reached a fixpoint unless the inner loop adds an observed equation
        active = false
        for i in eachindex(algeqs)
            # Ignore if we're not considering this for elimination or it is already eliminated
            active_alg_eqs[i] || continue
            alg_eqs_mask[i] || continue
            eq = algeqs[i]
            candidate_var = eq.lhs
            # LHS must be an unknown and must not be another observed
            if !(candidate_var in all_dvs) || candidate_var in new_observed_variables
                active_alg_eqs[i] = false
                candidate_eqs_count -= 1
                continue
            end
            # Remove newly added observed variables
            vars_in_algeq = vars_in_each_algeq[i]
            setdiff!(vars_in_algeq, new_observed_variables)
            # If the incidence is empty, it is a function of observed and diffvars
            isempty(vars_in_algeq) || continue

            # We added an observed equation, so we haven't reached a fixpoint yet
            active = true
            push!(new_observed_variables, candidate_var)
            push!(obseqs, eq)
            # This is no longer considered for elimination
            active_alg_eqs[i] = false
            candidate_eqs_count -= 1
            # And is no longer algebraic
            alg_eqs_mask[i] = false
        end
        # Safeguard against infinite loops, because `while true` is potentially dangerous
        maxiters -= 1
    end

    return keepat!(algeqs, alg_eqs_mask)
end

function offset_array(origin, arr)
    if all(isone, origin)
        return arr
    end
    return Origin(origin)(arr)
end

@register_array_symbolic offset_array(origin::Any, arr::AbstractArray) begin
    size = size(arr)
    eltype = eltype(arr)
    ndims = ndims(arr)
end

function add_array_observed!(obseqs::Vector{Equation})
    array_obsvars = Set{SymbolicT}()
    for eq in obseqs
        arr, isarr = split_indexed_var(eq.lhs)
        isarr && push!(array_obsvars, arr)
    end
    for var in array_obsvars
        firstind = first(SU.stable_eachindex(var))::SU.StableIndex{Int}
        firstind = Tuple(firstind.idxs)
        scal = SymbolicT[]
        for i in SU.stable_eachindex(var)
            push!(scal, var[i])
        end
        push!(obseqs, var ~ offset_array(firstind, reshape(scal, size(var))))
    end
    return
end

function simplify_sde_system(sys::AbstractSystem; kwargs...)
    brown_vars = brownians(sys)
    @set! sys.brownians = SymbolicT[]
    sys = __mtkcompile(sys; kwargs...)

    new_eqs = copy(equations(sys))
    Is = Int[]
    Js = Int[]
    vals = SymbolicT[]
    for (i, eq) in enumerate(new_eqs)
        resid = eq.rhs
        for (j, bvar) in enumerate(brown_vars)
            coeff, resid, islin = Symbolics.linear_expansion(resid, bvar)
            if !islin
                throw(
                    ArgumentError(
                        """
                        Expected brownian variables to appear linearly in equations. Brownian $bvar \
                        appears non-linearly in equation $eq.
                        """
                    )
                )
            end
            _iszero(coeff) && continue

            push!(Is, i)
            push!(Js, j)
            push!(vals, coeff)
        end
        new_eqs[i] = eq.lhs ~ resid
    end

    g = Matrix(sparse(Is, Js, vals, length(new_eqs), length(brown_vars)))
    @set! sys.eqs = new_eqs
    # Fix for https://github.com/SciML/ModelingToolkit.jl/issues/2490
    if size(g, 2) == 1
        # If there's only one brownian variable referenced across all the equations,
        # we get a Nx1 matrix of noise equations, which is a special case known as scalar noise
        noise_eqs = reshape(g[:, 1], (:, 1))
        is_scalar_noise = true
    elseif __num_isdiag_noise(g)
        # If each column of the noise matrix has either 0 or 1 non-zero entry, then this is "diagonal noise".
        # In this case, the solver just takes a vector column of equations and it interprets that to
        # mean that each noise process is independent
        noise_eqs = __get_num_diag_noise(g)
        is_scalar_noise = false
    else
        noise_eqs = g
        is_scalar_noise = false
    end

    dummy_sub = Dict{SymbolicT, SymbolicT}()
    for eq in new_eqs
        isdiffeq(eq) || continue
        dummy_sub[eq.lhs] = eq.rhs
    end
    var_sccs = [collect(eachindex(new_eqs))]
    schedule = Schedule(var_sccs, dummy_sub)

    @set! sys.eqs = new_eqs
    @set! sys.noise_eqs = noise_eqs
    @set! sys.schedule = schedule
    return sys
end

function simplify_optimization_system(sys::System; split = true, kwargs...)
    sys = flatten(sys)
    cons = constraints(sys)
    econs = Equation[]
    icons = Inequality[]
    for e in cons
        if e isa Equation
            push!(econs, e)
        elseif e isa Inequality
            push!(icons, e)
        end
    end
    irreducible_subs = Dict{SymbolicT, SymbolicT}()
    dvs = SymbolicT[]
    for var in unknowns(sys)
        sh = SU.shape(var)::SU.ShapeVecT
        if isempty(sh)
            push!(dvs, var)
        else
            append!(dvs, vec(collect(var)::Array{SymbolicT})::Vector{SymbolicT})
        end
    end
    for i in eachindex(dvs)
        var = dvs[i]
        if hasbounds(var)
            irreducible_subs[var] = irrvar = setirreducible(var, true)::SymbolicT
            dvs[i] = irrvar
        end
    end
    subst = SU.Substituter{false}(irreducible_subs, SU.default_substitute_filter)
    for i in eachindex(econs)
        econs[i] = subst(econs[i])
    end
    nlsys = System(econs, dvs, parameters(sys); name = :___tmp_nlsystem)
    snlsys = mtkcompile(nlsys; kwargs..., fully_determined = false)::System
    obs = observed(snlsys)
    seqs = equations(snlsys)
    trueobs = observed(unhack_system(snlsys))
    subs = Dict{SymbolicT, SymbolicT}()
    for eq in trueobs
        subs[eq.lhs] = eq.rhs
    end
    cons_simplified = Union{Equation, Inequality}[]
    for eq in seqs
        push!(cons_simplified, fixpoint_sub(eq, subs))
    end
    for eq in icons
        push!(cons_simplified, fixpoint_sub(eq, subs))
    end
    setdiff!(dvs, keys(subs))
    newsts = dvs
    @set! sys.constraints = cons_simplified
    newobs = copy(observed(sys))
    append!(newobs, obs)
    @set! sys.observed = newobs
    newcosts = copy(get_costs(sys))
    for i in eachindex(newcosts)
        newcosts[i] = fixpoint_sub(newcosts[i], subs)
    end
    @set! sys.costs = newcosts
    @set! sys.unknowns = newsts
    return sys
end

function __num_isdiag_noise(mat)
    for i in axes(mat, 1)
        nnz = 0
        for j in axes(mat, 2)
            nnz += !_iszero(mat[i, j])
        end
        if nnz > 1
            return (false)
        end
    end
    return true
end

function __get_num_diag_noise(mat::Matrix{SymbolicT})
    result = fill(Symbolics.COMMON_ZERO, size(mat, 1))
    for i in axes(mat, 1)
        for j in axes(mat, 2)
            mij = mat[i, j]
            _iszero(mij) && continue
            result[i] = mij
            break
        end
    end
    return result
end

"""
    $TYPEDSIGNATURES

Given observed equations `eqs` and a list of variables `unknowns`, construct the incidence
graph for the equations. Also construct a `Vector{Int}` mapping indices of `eqs` to the
index in `unknowns` of the observed variable on the LHS of each equation. Return the
constructed incidence graph and index mapping.
"""
function observed2graph(eqs::Vector{Equation}, unknowns::Vector{SymbolicT})::Tuple{BipartiteGraph{Int, Nothing}, Vector{Int}}
    graph = BipartiteGraph(length(eqs), length(unknowns))
    v2j = Dict{SymbolicT, Int}(unknowns .=> 1:length(unknowns))

    # `assigns: eq -> var`, `eq` defines `var`
    assigns = similar(eqs, Int)
    vars = Set{SymbolicT}()
    for (i, eq) in enumerate(eqs)
        lhs_j = get(v2j, eq.lhs, nothing)
        lhs_j === nothing &&
            throw(ArgumentError("The lhs $(eq.lhs) of $eq, doesn't appear in unknowns."))
        assigns[i] = lhs_j
        empty!(vars)
        SU.search_variables!(vars, eq.rhs; is_atomic = OperatorIsAtomic{SU.Operator}())
        for v in vars
            j = get(v2j, v, nothing)
            if j isa Int
                add_edge!(graph, i, j)
            end
        end
    end

    return graph, assigns
end

"""
    $(TYPEDSIGNATURES)

Use Kahn's algorithm to topologically sort observed equations.

Example:
```julia
julia> t = ModelingToolkit.t_nounits

julia> @variables x(t) y(t) z(t) k(t)
(x(t), y(t), z(t), k(t))

julia> eqs = [
           x ~ y + z
           z ~ 2
           y ~ 2z + k
       ];

julia> ModelingToolkit.topsort_equations(eqs, [x, y, z, k])
3-element Vector{Equation}:
 Equation(z(t), 2)
 Equation(y(t), k(t) + 2z(t))
 Equation(x(t), y(t) + z(t))
```
"""
function topsort_equations(eqs::Vector{Equation}, unknowns::Vector{SymbolicT}; check = true)
    graph, assigns = observed2graph(eqs, unknowns)
    neqs = length(eqs)
    degrees = zeros(Int, neqs)

    for 𝑠eq in 1:length(eqs)
        var = assigns[𝑠eq]
        for 𝑑eq in 𝑑neighbors(graph, var)
            # 𝑠eq => 𝑑eq
            degrees[𝑑eq] += 1
        end
    end

    q = Queue{Int}(neqs)
    for (i, d) in enumerate(degrees)
        @static if pkgversion(DataStructures) >= v"0.19"
            d == 0 && push!(q, i)
        else
            d == 0 && enqueue!(q, i)
        end
    end

    idx = 0
    ordered_eqs = similar(eqs, 0)
    sizehint!(ordered_eqs, neqs)
    while !isempty(q)
        @static if pkgversion(DataStructures) >= v"0.19"
            𝑠eq = popfirst!(q)
        else
            𝑠eq = dequeue!(q)
        end
        idx += 1
        push!(ordered_eqs, eqs[𝑠eq])
        var = assigns[𝑠eq]
        for 𝑑eq in 𝑑neighbors(graph, var)
            degree = degrees[𝑑eq] = degrees[𝑑eq] - 1
            @static if pkgversion(DataStructures) >= v"0.19"
                degree == 0 && push!(q, 𝑑eq)
            else
                degree == 0 && enqueue!(q, 𝑑eq)
            end
        end
    end

    (check && idx != neqs) && throw(ArgumentError("The equations have at least one cycle."))

    return ordered_eqs
end
