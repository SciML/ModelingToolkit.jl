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


#=
TODO: A clock may need a unique id since otherwise Clock(t, 0.1) === Clock(t, 0.1)
=#

function add_varmap!(varmap, e, ::Nothing)
    return get!(varmap, e, nothing) # throw away the nothing value
end
function add_varmap!(varmap, e, domain::TimeDomain)
    if  e isa Sym
        return domain # This discards variables without independet variable from the varmap
    end
    if e ∈ keys(varmap)
        if varmap[e] isa UnknownDomain
            varmap[e] = domain
        elseif domain isa UnknownDomain
            return varmap[e]
        else
            domain == varmap[e] || throw(ClockInferenceException("Variable $e was associate with more than one domain, $domain and $(varmap[e])"))
        end
    else
        varmap[e] = domain
    end
    domain
end

function propagate_time_domain(eq::Equation, domain_above::Union{TimeDomain, Nothing} = nothing, varmap = Dict{Any, Any}())
    last_hash = hash(varmap)
    ex = eq.lhs - eq.rhs
    for i = 1:10 # max 10 step fixed-point iteration
        # The fix-point iteration is performed in order to make 
        d, varmap = propagate_time_domain(ex, domain_above, varmap)
        domain_above = merge_domains(domain_above, d, eq)
        # d, varmap = propagate_time_domain(eq.rhs, domain_above, varmap)
        # domain_above = merge_domains(domain_above, d, eq)
        h = hash(varmap)
        h == last_hash && break
        h = last_hash 
    end
    domain_above, varmap
end


# Nums are just unwrapped
propagate_time_domain(e::Num, domain_above::Union{TimeDomain, Nothing} = nothing, args...) = propagate_time_domain(value(e), domain_above, args...)


"""
    domain, varmap = propagate_time_domain(e, domain_above::Union{TimeDomain, Nothing} = nothing, varmap = Dict{Any, Any}())

Determine the time-domain (clock inference) of an expression or equation `e`. `domain` is the time domain for the compound expression `e`, while `varmap` is a dict that maps variables to time domains / clocks.

# Details
At the expression root, the domain is typically unknown (nothing), we thus start from the root and recurse down until qe encounter an expression with known domain. When we find a domain, we continue to propagate all the way down to the leaves, checking consistency on the way if the domain is specified also further down the tree. We then propagate the domain back up to the root. At a middle point, there might be one domain coming from above, and one from below, if they differ it's an error.

In a DiscreteSystem, the root starts out as [`Inferred()`](@ref). If it remains inferred after inference, it's set to the default `Clock(t, 1)` which corresponds to the default time step for discrete systems in DifferentialEquations.
"""
function propagate_time_domain(e, domain_above::Union{TimeDomain, Nothing} = nothing, varmap = Dict{Any, Any}())
    e isa Symbolics.Symbolic || return (domain_above, varmap)
    if istree(e)
        f = operation(e)
        if f isa Operator # operators are the only things that can change the domain
            domain_below = output_timedomain(f)
            return_domain = merge_domains(domain_above, domain_below, e) # this is the domain we return upwards, but on the way down we use the operator input time domain
            downward_domain = input_timedomain(f) # we expect this domain further down
            domain_below, varmap = propagate_time_domain(only(arguments(e)), downward_domain, varmap) # the received domain from below only matters if the operator has inferred output domain
            if return_domain isa InferredDomain
                # This is where the Inferred is inferred
                return_domain = merge_domains(return_domain, domain_below, e)
            end
            return return_domain, varmap
        end
        if !has_time_domain(e) && !(f isa Sym) # vars created with @variables x(t) [timedomain=d] will be terms with metadata
            return propagate_time_domain(arguments(e), domain_above, varmap)
        end
    end
    
    # If we reach here, e is a Sym or Term with metadata
    if has_time_domain(e)
        domain_below = get_time_domain(e)
        domain_above = merge_domains(domain_above, domain_below, e)
    end

    stored_domain = add_varmap!(varmap, e, domain_above)
    domain_above = merge_domains(domain_above ,stored_domain, e)
    # QUESTION: if the domain_above is nothing here, should we set continuous, do everything one more time now that the domain of the root is known, or error?

    # NOTE: the metadata is immutable so we can not set it, we are thus currently only able to classify the tome domain expressions, but not set the domain or variables
    # setmetadata!(e, TimeDomain, domain_above)
    return domain_above, varmap
end


function propagate_time_domain(v::AbstractVector, domain_above::Union{TimeDomain, Nothing} = nothing, varmap = Dict{Any, Any}())
    for arg in v
        domain_below, varmap = propagate_time_domain(arg, domain_above, varmap)
        domain_above = merge_domains(domain_above, domain_below, v)
    end
    return domain_above, varmap
end


function equation_and_variable_time_domains(eqs::Vector{Equation}, domain_above::Union{TimeDomain, Nothing} = nothing)
    varmap = Dict{Any, Any}()
    eq2dom = Vector{Any}(undef, length(eqs)) .= domain_above
    last_hash = UInt(0)
    while true # fix-point iteration: iterate until the variable map no longer changes. This allows domain information to propagate between equations.
        for (j, eq) in pairs(eqs)
            d, varmap = propagate_time_domain(eq, eq2dom[j], varmap)
            # domain_above = merge_domains(domain_above, d, eq) # we should not merge the domain here since different eqs can have different domains
            eq2dom[j] = d
        end
        h = hash(varmap)
        h == last_hash && break
        last_hash = h
    end
    eq2dom, varmap
end

function propagate_time_domain(sys::DiscreteSystem)
    d = propagate_time_domain(equations(sys), nothing)
    d === nothing && return Clock(get_iv(sys), 1)
    d
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
    preprocess_hybrid_equations(eqs)

1. Perform clock inference using [`equation_and_variable_time_domains`](@ref).
2. Normalize `Shift` operations so that the largest positive shift becomes 0, and is the only term on the LHS of the equation. This uses [`normalize_shifts`](@ref).
3. For each clock partition, expand shift equations to companion form where the largest shift is -1. This introduces max_shift-1 new states and equations. 
4. Strip away hybrid and discrete operators from equations (domains are known from the clock inference). Uses [`strip_operator`](@ref)

"""
function preprocess_hybrid_equations(eqs::AbstractVector{Equation})
    eqmap, varmap = equation_and_variable_time_domains(eqs)

    eqs = map(zip(eqs, eqmap)) do (eq, domain)
        if domain isa Continuous
            eq = strip_operator(eq, Hold)
        elseif domain isa AbstractDiscrete
            eq = strip_operator(eq, Sample)
            if hasshift(eq)
                eq = normalize_shifts(eq)
            end
        end
        eq
    end

    clocks = unique(eqmap)
    new_eqs_and_vars = map(clocks) do clock
        cp_inds = findall(==(clock), eqmap) 
        cp_eqs = eqs[cp_inds] # get all equations that belong to this particular clock partition
        clock isa Continuous && return cp_eqs # only work on discrete eqs
        expand_shifts(cp_eqs, clock)
    end
    new_eqs, new_vars = first.(new_eqs_and_vars), last.(new_eqs_and_vars)
    reduce(vcat, new_eqs), reduce(vcat, new_vars)
end

"""
    normalize_shifts(eq)

Normalize `Shift` operations so that the largest positive shift becomes 1, and is the only term on the LHS of the equation.
Example:
```julia
julia> eq = Shift(t, 1)(u) + Shift(t, 3)(u) ~ 0
Shift(t, 1)(u(t)) + Shift(t, 3)(u(t)) ~ 0

julia> normalize_shifts(eq)
Shift(t, 1)(u(t)) ~ -Shift(t, -1)(u(t))
```
"""
function normalize_shifts(eq::Equation)
    eq = insert_zero_shifts(eq)
    ops = collect(collect_applied_operators(eq, Shift)) # this leaves out u(k) (no shift)

    maxshift, maxind = findmax(op->op.f.steps, ops)
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
    newops[maxind] ~ rhs
end

is_function_of(t) = function(var)
    var isa Term                       &&
        operation(var) isa Sym         &&
        length(arguments(var)) == 1    &&
        isequal(arguments(var)[1], t)
end


"""
    insert_zero_shifts(eq)

Inserts `Shift(t, 0)` terms, i.e., replace `u(t)` with `Shift(t, 0)(u(t))`.
"""
function insert_zero_shifts(eq)
    vars = collect(collect_operator_variables(eq, Shift))
    ops = collect(collect_applied_operators(eq, Shift))
    t = operation(ops[1]).t
    r1 = @rule ~var::is_function_of(t) => Shift(t,0)(~var, true)
    r = AbortablePrewalk(PassThrough(r1), cond=isoperator(Shift)) # We should not insert a shift inside another shift
    r(eq.lhs) ~ r(eq.rhs)
end

"""
    expand_shifts(eqs::Vector{Equation}, clock)

Takes a vector of equations, all belonging to the same discrete clock partition, and introduces new variables for terms that are shifted by more than -1. For example
```
u(k+1) ~ u(k) + u(k-1) + u(k-2)
```
becomes
```
u(k+1)  ~  u(k) + u1(k) + u2(k)
u1(k+1) ~  u(k)
u2(k+1) ~ u1(k)
```
"""
function expand_shifts(eqs::Vector{Equation}, clock::AbstractDiscrete)
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
        new_eq_vars = @variables $varname[1:n_new_states](t)
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
"""
function strip_operator(eq::Equation, operator)
    ops = collect(collect_applied_operators(eq, operator))
    args = map(ops) do op
        op.arguments[1]
    end
    subs = Dict(ops .=> args)
    eq = substitute(eq, subs)
end