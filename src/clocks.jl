#=
# Clock inferece functionality
The semantics of clocked variables roughly follows
https://specification.modelica.org/v3.4/Ch16.html

In particular
> The base clock partitions are identified as clocked or as continuous-time partitions according to the following properties:
> A variable u in sample(u) and a variable y in y = hold(ud) is in a continuous-time partition.
> Correspondingly, variables u and y in y = sample(uc), y = subSample(u), y = superSample(u), y = shiftSample(u), y = backSample(u), y = previous(u), are in a clocked partition. Equations in a clocked when clause are also in a clocked partition. Other partitions where none of the variables in the partition are associated with any of the operators above have an unspecified partition kind and are considered continuous-time partitions.
>Since sample(u) returns the left limit of u, and the left limit of u is a known value, all inputs to a base-clock partition are treated as known during sorting.


=#
using Symbolics: value
using SymbolicUtils.Rewriters: Postwalk, Prewalk, PassThrough, Empty
import SymbolicUtils.Rewriters

struct ClockInferenceException <: Exception
    msg
end

abstract type AbstractClock <: AbstractDiscrete end


"""
    Clock <: AbstractClock
    Clock(t; dt)

The default periodic clock with independent variables `t` and tick interval `dt`.
If `dt` is left unspecified, it will be inferred (if possible).
"""
struct Clock <: AbstractClock
    "Independent variable"
    t
    "Period"
    dt
    Clock(t, dt=nothing) = new(value(t), dt)
end

sampletime(c::AbstractClock) = c.dt
#=
TODO: A clock may need a unique id since otherwise Clock(t, 0.1) === Clock(t, 0.1)
=#

# The abortable prewal lets you do the transformation `f(t) + t => f(t) + g(t)`
struct AbortablePrewalk{C, F, CO}
    rw::C
    similarterm::F
    cond::CO
end

function AbortablePrewalk(rw; similarterm=similarterm, cond)
    AbortablePrewalk{typeof(rw), typeof(similarterm), typeof(cond)}(rw, similarterm, cond)
end

function (p::AbortablePrewalk{C, F})(x) where {C, F}
    p.cond(x) && return x
    if istree(x)
        x = p.rw(x)
        p.cond(x) && return x
        if istree(x)
            x = p.similarterm(x, operation(x), map(PassThrough(p), Rewriters.unsorted_arguments(x)))
            p.cond(x) && return x
        end
        return x
    else
        return p.rw(x)
    end
end

struct ClockPartitioning
    eqs::Vector{Equation}
    "A vector of TimeDomain of the same length as ccs"
    clocks::Vector{TimeDomain}
    "A vector of TimeDomain of the same length as the number of equations"
    eqmap::Vector{TimeDomain}
    "A vector of all variables in the inferred equation set"
    vars
    "A dict that maps variables to TimeDomain"
    varmap::Dict
end
struct ClockInferenceResult
    clockpart::ClockPartitioning
    "A Vector{Vector{Int}}. incidences[i] points to the unknown variables, as well as variables x in Differential(x), and Shift(x), which lexically appear in eqs[i] except as first argument of base-clock conversion operators: Sample() and Hold(). Example, if j ∈ incidences[eq_ind], then allvars[j] appears in eqs[eq_ind]"
    incidences
    "A vector of vectors with equation indices that form each connected component"
    ccs::Vector{Vector{Int}}
    bg::BipartiteGraph
end

function Base.getproperty(cir::ClockInferenceResult, s::Symbol)
    s ∈ fieldnames(typeof(cir)) && return getfield(cir, s)
    getproperty(getfield(cir, :clockpart), s)
end



# The rules are, `nothing` loses to everything, Inferred loses to everything else, then InferredDiscrete and so on
function merge_inferred(x, y)
    x === nothing && return y
    y === nothing && return x
    x isa Inferred && return y
    y isa Inferred && return x
    x isa InferredDiscrete && !(y isa Continuous) && return y
    y isa InferredDiscrete && !(x isa Continuous) && return x
    x == y || throw(ClockInferenceException("Cannot merge $x and $y"))
    x
end


# nothing cases
merge_domains(::Nothing, ::Nothing, args...) = nothing
merge_domains(x::Any, ::Nothing, args...) = x
merge_domains(::Nothing, x::Any, args...) = x


"""
    merge_domains(above::A, below::B, ex)

Takes in two time domains and returns the most specific one. The specificity rules are
1.  A concrete domain like `Continuous` or `Clock`
2. InferredDiscrete
3. Inferred
4. nothing

If both arguments are concrete domains, they must be equal, otherwise a ClockInferenceException is thrown.
"""
function merge_domains(above::A, below::B, ex) where {A <: TimeDomain, B <: TimeDomain}
    emsg = "In expression $ex, the domain was inferred to $above from the root of the tree but $below from the leaves."
    above == below && return above
    if above isa InferredDomain || below isa InferredDomain
        return merge_inferred(above, below)
    end
    if above isa Continuous || below isa Continuous
        above === below || throw(ClockInferenceException(emsg))
        return Continuous()
    end

    throw(ClockInferenceException(emsg))
end


"""
    preprocess_hybrid_equations(eqs, states)

1. Perform clock inference using [`clock_inference`](@ref).
2. Normalize `Shift` operations so that the largest positive shift becomes 1, and is the only term on the LHS of the equation. This uses [`normalize_shifts`](@ref).
3. For each clock partition, expand shift equations to companion form where the only remaining shifts are 1 and 0. This introduces `n` new states and equations where `n` = max_shift - min_shift - 1. 
4. Rewrite Shift(1) to the continuous-time ` [`DiscreteUpdate`](@ref) operator and strip away hybrid and discrete operators from equations (domains are known from the clock inference). Uses [`strip_operator`](@ref)
"""
function preprocess_hybrid_equations(eqs::AbstractVector{Equation}, original_vars)
    cires = clock_inference(eqs)
    @unpack eqmap, varmap = cires
    if any(d isa UnknownDomain for d in eqmap)
        display(eqs .=> eqmap)
        throw(ClockInferenceException("Clock inference failed to infer the time domain of all equations, the inference result was printed above."))
    end
    eqs = map(zip(eqs, eqmap)) do (eq, domain)
        if domain isa Continuous
            eq = strip_operator(eq, Hold)
        elseif domain isa AbstractDiscrete
            eq = normalize_shifts(eq, varmap)
            @assert !hassample(eq)
        end
        eq
    end

    clocks = unique(eqmap)
    sort!(clocks, by = c -> c === Continuous() ? 0.0 : sampletime(c)) # Continuous partition will be first, then ascending sampletime
    # QUESTION: we now process each clock partition independently. We possibly want to perform additional simplification on clock partitions independently?
    # TODO: If we want to handle discrete and continuous states separately, we definitely want to separate them
    new_eqs_and_vars = map(clocks) do clock
        cp_inds  = findall(==(clock), eqmap) 
        cp_eqs   = eqs[cp_inds] # get all equations that belong to this particular clock partition
        clock isa Continuous && return (cp_eqs, []) # only work on discrete eqs
        alg_inds = [!isoperator(eq.lhs, Operator) for eq in cp_eqs] # only expand shifts for non-algeraic equations
        new_eqs, new_vars = shift_order_lowering(cp_eqs[.!alg_inds], clock)
        for nv in new_vars
            varmap[nv] = clock # add the new variables to the varmap
        end
        @assert !any(hassample, new_eqs)
        new_eqs  = discrete2continuous_operators.(new_eqs, clock)
        # new_eqs = reduce(vcat, new_eqs)
        @assert !any(hassample, new_eqs)
        new_eqs  = [new_eqs; cp_eqs[alg_inds]]
        
        new_eqs, new_vars
    end

    new_eqs, new_vars = first.(new_eqs_and_vars), last.(new_eqs_and_vars)
    eqmap = reduce(vcat, [fill(c, length(eqs)) for (c,eqs) in zip(clocks, new_eqs)]) # recreate shuffled eqmap
    new_eqs = reduce(vcat, new_eqs)
    new_vars = reduce(vcat, new_vars)
    # new_eqs, eqmap, removed_vars = substitute_algebraic_eqs(new_eqs, eqmap) # broken
    removed_vars = []
    # sort equations so that all operator equations come first
    perm = sortperm(new_eqs, by = eq->!isoperator(eq.lhs, Operator))
    new_eqs = new_eqs[perm]
    eqmap = eqmap[perm]


    allvars = [setdiff(original_vars, removed_vars); new_vars]
    part = ClockPartitioning(new_eqs, clocks, eqmap, allvars, varmap)
end

function discrete_var2param(part::ClockPartitioning)
    @unpack eqs, clocks, eqmap, vars, varmap = part
    vars = map(enumerate(vars)) do (i, var)
        varmap[var] isa AbstractDiscrete || return var # we only need to change variables in the continuous partition
        param = toparam(var)
        eqs .= substitute.(eqs, Ref(Dict(var => param)))
        vars[i] = param
    end
    ClockPartitioning(eqs, clocks, eqmap, vars, varmap)
end

"""
    normalize_shifts(eq, varmap)

Normalize `Shift` operations so that the largest positive shift becomes 1, and is the only term on the LHS of the equation. This function also removes `Sample` operators from the equation.
Example:
```julia
julia> eq = Shift(t, 1)(u) + Shift(t, 3)(u) ~ 0
Shift(t, 1)(u(t)) + Shift(t, 3)(u(t)) ~ 0

julia> normalize_shifts(eq)
Shift(t, 1)(u(t)) ~ -Shift(t, -1)(u(t))
```
Equations have a "floating" shift index, meaning that a shift equation is equivalent to another shift equation where all shifts are related by a constant offset. Equations may be shift-normalized independently, but different variables in the same equation are *not* normalized independently.
"""
function normalize_shifts(eq0::Equation, varmap)
    eq = insert_zero_shifts(eq0, varmap)
    is_discrete_algebraic(eq) && !hassample(eq) && return eq0
    ops = collect(collect_applied_operators(eq, Shift)) 

    maxshift, maxind = findmax(op->op.f.steps, ops)
    minshift = minimum(op->op.f.steps, ops)
    maxshift == minshift == 0 && return strip_operator(eq, Sample) # This implies an algebraic equation without delay
    newops = map(ops) do op
        s = op.f
        Shift(s.t, s.steps-maxshift+1)(op.arguments[1]) # arguments[1] is the shifted variable
    end
    subs = Dict(ops .=> newops)
    eq = substitute(eq, subs) # eq now has normalized shifts
    highest_order_term = ops[maxind].arguments[1]
    @variables __placeholder__
    eq = substitute(eq, Dict(newops[maxind]=>__placeholder__)) # solve_for can not solve for general terms, so we replace temporarily
    rhs = solve_for(eq, __placeholder__)
    strip_operator(newops[maxind] ~ rhs, Sample)
end


"""
    insert_zero_shifts(eq)

Inserts `Shift(t, 0)` terms, i.e., replace `u(t)` with `Shift(t, 0)(u(t))`.
Does not insert shifts inside other operators, e.g., `Shift(t, 1)(u(t))` does not become
`Shift(t, 1)(Shift(t, 0)(u(t)))`
"""
function insert_zero_shifts(eq, varmap)
    discrete_vars = []
    local t
    for k in keys(varmap)
        if varmap[k] isa AbstractDiscrete
            push!(discrete_vars, k)
            t = varmap[k].t
        end
    end
    discrete_vars = Set(discrete_vars)
    predicate = var -> var ∈ discrete_vars
    r1 = @rule ~var::predicate => Shift(t,0)(~var, true)
    r = AbortablePrewalk(PassThrough(r1), cond=isoperator(Operator)) # We should not insert a shift inside another operator
    r(eq.lhs) ~ r(eq.rhs)
end

"""
    new_eqs, new_vars = shift_order_lowering(eqs::Vector{Equation}, clock)

Takes a vector of equations, all belonging to the same discrete clock partition, and introduces new variables for terms that are shifted by more than -1. For example
```
u(k+1) ~ u(k) + u(k-1) + u(k-2)
```
becomes
```
u(k+1)  ~  u(k) + u_delay[1](k) + u_delay[2](k)
u_delay[1](k+1) ~ u(k)
u_delay[2](k+1) ~ u_delay[1](k)
```
This corresponds to [`ode_order_lowering`](@ref) for continuous systems.

This function assumes that shifts are already normalized using [`normalize_shifts`](@ref).
"""
function shift_order_lowering(eqs::Vector{Equation}, clock::AbstractDiscrete)
    t = clock.t
    shift(x) = Shift(t)(x)
    vars = collect(union(collect_operator_variables.(eqs, Ref(Shift))...))
    # eqs = insert_zero_shifts(eqs, vars)
    ops  = collect(union(collect_applied_operators.(eqs, Ref(Shift))...))
    
    new_eqs = Equation[]
    new_vars = []
    # for each variable that appears shifted, find the maximum (negative) shift, introduce new states and substitute shifted occurances by new variables
    for var in vars       
        varops = filter(ops) do op # extract shifted var terms 
            isequal(arguments(op)[1], var)
        end
        maximum(op->op.f.steps, varops) <= 1 || throw(ArgumentError("Shift equations must be normalized to have shifts <= 1, see normalize_shifts(eq)"))
        n_new_states = -minimum(op->op.f.steps, varops) # number of new states required is maximum negative shift
        n_new_states < 1 && continue
        varname = Symbol(string(value(var).f)*"_delay")
        new_eq_vars = @variables $varname[1:n_new_states](t) [timedomain=clock]
        new_eqvars = collect(new_eq_vars[]) # @variables returns a vector wrapper
        new_vareqs = map(1:n_new_states-1) do i
            shift(new_eqvars[i+1]) ~ new_eqvars[i]
        end
        pushfirst!(new_vareqs, shift(new_eqvars[1]) ~ var)
        sort!(varops, by=op->operation(op).steps, rev=true) # highest shift first
        startind = findfirst(op->operation(op).steps <= -1, varops) # we should not replace shifts 1 and 0

        subs = Dict(varops[startind:end] .=> new_eqvars)
        eqs = substitute.(eqs, Ref(subs))
        
        append!(new_eqs, new_vareqs)
        append!(new_vars, new_eqvars)
    end

    [eqs; new_eqs], new_vars
end

"""
    strip_operator(eq, operator)

Removes operators from equation `eq`, example
```julia
julia> eq = u ~ Hold(ud) + 1
u(t) ~ 1 + Hold()(ud(t))

julia> strip_operator(eq, Hold)
u(t) ~ 1 + ud(t)
```
"""
function strip_operator(eq::Equation, operator)
    ops = collect(collect_applied_operators(eq, operator))
    args = map(ops) do op
        op.arguments[1]
    end
    subs = Dict(ops .=> args)
    eq = substitute(eq, subs)
end

strip_operator(ex, operator) = strip_operator(ex~0, operator).lhs


"""
    is_discrete_algebraic(eq)

Returns true if the equation is a discrete algebraic equation. This function assumes that the equation has been shift normalized using [`normalize_shifts`](@ref).
"""
function is_discrete_algebraic(eq)
    t = eq.lhs
    t isa Term && operation(t) isa Shift && operation(t).steps == 0
end


"""
    discrete2continuous_operators(eq, clock)

Rewrites `eq` by removing all shift operators, and replaces `Shift(1)` on the eq.lhs with the continuous-time [`DiscreteUpdate`](@ref) operator which lowers to a [`DiscreteCallback`](@ref) when an [`ODEProblem`](@ref) is created.
"""
function discrete2continuous_operators(eq, clock)
    t = clock.t
    dt = sampletime(clock)
    term = value(eq.lhs)
    var = arguments(term)[1]
    @assert !hassample(eq)
    if is_discrete_algebraic(eq)
        deq = DiscreteUpdate(t; dt=dt)(var) ~ strip_operator(eq.rhs, Shift)
    else
        deq = DiscreteUpdate(t; dt=dt)(var) ~ eq.rhs
    end
    # [deq; Differential(t)(var) ~ 0] # get segfault if simulating with D(x) ~ 0 equaiton added
    deq
end

function substitute_algebraic_eqs(eqs, eqmap)
    n = typemax(Int)
    removed_states = []
    while length(eqs) < n
        n = length(eqs)
        for (i, eq) in enumerate(eqs)
            isoperator(eq.lhs, Operator) && continue
            eqmap[i] === Continuous() && continue # may be possible to omit this
            for (j, other_eq) in pairs(eqs)
                other_eq === eq && continue
                eqmap[j] === Continuous() && continue
                eqs[j] = substitute(other_eq, Dict(eq.lhs => eq.rhs))
            end
            if eqmap[i] isa AbstractDiscrete
                push!(removed_states, eq.lhs)
                eqs[i] = 0 ~ 0
            end
        end
        keep_inds = map(eq->!(isequal(eq.rhs, 0) && isequal(eq.lhs, 0)), eqs)
        eqs = eqs[keep_inds] # remove empty equations
        eqmap = eqmap[keep_inds]
    end
    eqs, eqmap, removed_states
end


"""
    hybrid_simplify(sys::ODESystem)

Simplification of hybrid systems (containing both discrete and continuous states).
The equations of the returned system will not contain any discrete operators, and shift equations have been rewritten into [`DiscreteUpdate`](@ref) equations. This simplification pass may introduce new variables and equations for equations with relative shifts larger than 1, corresponding to [`ode_order_lowering`](@ref) for continuous systems.
"""
function hybrid_simplify(sys::ODESystem; param=false)
    part = preprocess_hybrid_equations(equations(sys), states(sys))
    # @set! sys.eqs = new_eqs # this only modifies the equations of the outer system, not inner systems
    # @set! sys.states = [states(sys); new_vars]
    # TODO: propagate all other fields to ODESystem
    if param
        part = discrete_var2param(part)
        new_params = filter(isparameter, part.vars)
        new_states = filter(!isparameter, part.vars)
        return ODESystem(part.eqs, get_iv(sys), new_states, [new_params; parameters(sys)], name=sys.name, clock_partitioning=part)
    else
        return ODESystem(part.eqs, get_iv(sys), part.vars, parameters(sys), name=sys.name, clock_partitioning=part)
    end
end

function get_eq_domain(x::Num, domain_above=Inferred())
    v = value(x)
    if v isa Term && operation(v) isa Sym
        return merge_domains(domain_above, get_time_domain(v), x)
    else
        return get_eq_domain(v, domain_above)
    end
end

"""
    domain = get_eq_domain(e, domain_above = Inferred())

Determine the time-domain (clock inference) of an expression or equation `e`.
"""
function get_eq_domain(eq, domain_above=Inferred())
    ovars = collect_applied_operators(eq, Operator)
    for op in ovars
        while isoperator(op, Operator) # handle nested operators
            odom = output_timedomain(op)
            domain_above = merge_domains(domain_above, odom, eq)
            changes_domain(op) && break 
            op = arguments(op)[1] # unwrap operator in case there are wrapped operators with a more specific output
        end
        if op isa Num
            if has_time_domain(op)
                domain_above = merge_domains(domain_above, get_time_domain(op), op)
            end
        end
    end
    domain_above
end

# incidence(e) = the unknown variables, as well as variables x in der(x), pre(x), and previous(x), which lexically appear in e except as first argument of base-clock conversion operators: sample() and hold().
"""
    clock_inference(eqs::Vector{Equation}, domain_above::TimeDomain = Inferred())

Performs clock-inference for a collection of equations, like `equations(sys)`
`domain_above` is the knowledge at the starting point. Returns a [`ClockInferenceResult`](@ref) containing `clockpart::ClockPartitioning`.
"""
function clock_inference(eqs::Vector{Equation}, domain_above::TimeDomain = Inferred())
    varmap = Dict{Any, Any}()
    eq2dom = Vector{Any}(undef, length(eqs)) .= domain_above

    innervars = vars(eqs; op=Nothing) |> collect
    filter!(istree, innervars) # only keep variables of time
    varinds = Dict(innervars .=> 1:length(innervars))

    # Build graph
    incidences = [Int[] for _ in eachindex(eqs)]
    for (eqind, eq) in pairs(eqs)
        eqvars = vars(eq; op=Operator)
        for var in eqvars
            changes_domain(var) && continue # Sample and Hold change the time domain
            while isoperator(var, Operator) # handle nested operators
                var = arguments(var)[1]
            end
            changes_domain(var) && continue # Sample and Hold change the time domain
            istree(var) || continue # skip if not function of time
            push!(incidences[eqind], varinds[var])
        end
    end

    # Infer variable and equation time domains
    allvars = vars(eqs; op=Operator) |> collect
    for var in allvars
        if has_time_domain(var)
            varmap[var] = get_time_domain(var)
        elseif isoperator(var, Operator)
            inner_var = arguments(var)[1]
            while isoperator(inner_var, Operator) # handle nested operators
                var = inner_var
                inner_var = arguments(var)[1]
            end

            idom = input_timedomain(var)
            varmap[inner_var] = idom
        end
    end
    for (eqind, eq) in pairs(eqs)
        eq2dom[eqind] = get_eq_domain(eq, eq2dom[eqind])
    end

    # Some equations and variables will still be uninferred, we now propagate
    # inferred domains to all other equations and variables
    bg = BipartiteGraph(incidences)
    im = incidence_matrix(bg)
    g2 = SimpleGraph([0I im; im' 0I]) # expand BipartiteGraph to full graph (I can't figure out how to use connected_components on BipartiteGraph)
    ccs = connected_components(g2)
    ccs = map(ccs) do cc
        filter!(<=(length(eqs)), cc) # "undo" graph expansion, keep only equation part
    end

    clocks = map(ccs) do cc
        dom = domain_above
        for eqind in cc
            dom = merge_domains(dom, eq2dom[eqind], eqs[eqind])
        end
        eq2dom[cc] .= dom
        for eqind in cc
            for varind in incidences[eqind]
                var = innervars[varind]
                varmap[var] = merge_domains(get(varmap, var, Inferred()), dom, var)
            end
        end
        dom
    end
    part = ClockPartitioning(eqs, clocks, eq2dom, allvars, varmap)
    ClockInferenceResult(part, incidences, ccs, bg)
end


function get_clocked_partitions(sys,  part = sys.clock_partitioning)
    t = independent_variable(sys)
    @unpack eqs, clocks, eqmap, vars, varmap = part
    sort!(clocks, by = c->c === Continuous() ? 0.0 : sampletime(c)) # cont., then ascending sampletime

    if clocks[1] === Continuous()
        continds = eqmap .== Continuous()
        cont = ODESystem(eqs[continds], t; name=sys.name) # TODO: propagate other properties
        clocks = clocks[2:end]
    else
        error("This function should not be called for completely discrete systems.")
    end
    disc = map(clocks) do clock
        inds = eqmap .== clock
        DiscreteSystem(eqs[inds], t; name=sys.name) # TODO: propagate other properties
    end
    cont, disc
end


# struct ClockPartitioning
#     eqs::Vector{Equation}
#     "A vector of TimeDomain of the same length as ccs"
#     clocks::Vector{TimeDomain}
#     "A vector of TimeDomain of the same length as the number of equations"
#     eqmap::Vector{TimeDomain}
#     "A vector of all variables in the inferred equation set"
#     vars
#     "A dict that maps variables to TimeDomain"
#     varmap::Dict
# end