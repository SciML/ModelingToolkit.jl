using Symbolics: value

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


function propagate_time_domain(eq::Equation, domain_above::Union{TimeDomain, Nothing} = nothing)
    d = propagate_time_domain(eq.lhs, domain_above)
    domain_above = merge_domains(domain_above, d, eq)
    d = propagate_time_domain(eq.rhs, domain_above)
    merge_domains(domain_above, d, eq)
end


# Nums are just unwrapped
propagate_time_domain(e::Num, domain_above::Union{TimeDomain, Nothing} = nothing) = propagate_time_domain(value(e), domain_above)


"""
    propagate_time_domain(e, domain_above::Union{TimeDomain, Nothing} = nothing)

Determine the time-domain (clock inference) of an expression or equation `e`.

# Details
At the expression root, the domain is typicaööy unknown (nothing), we thus start from the root and recurse down until qe encounter an expression with known domain. When we find a domain, we continue to propagate all the way down to the leaves, checking consistency on the way if the domain is specified also further down the tree. We then propagate the domain back up to the root. At a middle point, there might be one domain coming from above, and one from below, if they differ it's an error.

In a DiscreteSystem, the root starts out as [`Inferred()`](@ref). If it remains inferred after inference, it's set to the default `Clock(t, 1)` which corresponds to the default time step for discrete systems in DifferentialEquations.
"""
function propagate_time_domain(e, domain_above::Union{TimeDomain, Nothing} = nothing)
    if istree(e)
        f = operation(e)
        if f isa Operator # operators are the only things that can change the domain
            domain_below = output_timedomain(f)
            return_domain = merge_domains(domain_above, domain_below, e) # this is the domain we return upwards, but on the way down we use the operator input time domain
            downward_domain = input_timedomain(f) # we expect this domain further down
            domain_below = propagate_time_domain(only(arguments(e)), downward_domain) # the received domain from below only matters if the operator has inferred output domain
            if return_domain isa Inferred
                # This is where the Inferred is inferred
                return_domain = merge_domains(return_domain, domain_below, e)
            end
            return return_domain
        end
        if !has_time_domain(e) # vars created with @variables x(t) [timedomain=d] will be terms with metadata
            return propagate_time_domain(arguments(e), domain_above)
        end
    end
    
    # If we reach here, e is a Sym or Term with metadata
    if has_time_domain(e)
        domain_below = get_time_domain(e)
        return merge_domains(domain_above, domain_below, e)
    else
        # if the domain_above is nothing here, should we set continuous, do everything one more time now that the domain of the root is known, or error?
        if domain_above === nothing
            return nothing
        else
            # the metadata is immutable so we can not set it, we are thus currently only able to classify the tome domain expressions, but not set the domain or variables
            # setmetadata!(e, TimeDomain, domain_above)
            return domain_above
        end
    end
    
end


function propagate_time_domain(v::AbstractVector, domain_above::Union{TimeDomain, Nothing} = nothing)
    for arg in v
        domain_below = propagate_time_domain(arg, domain_above)
        domain_above = merge_domains(domain_above, domain_below, v)
    end
    return domain_above
end

function propagate_time_domain(sys::DiscreteSystem)
    d = propagate_time_domain(equations(sys), nothing)
    d === nothing && return Clock(get_iv(sys), 1)
    d
end


# The rules are, `nothing` loses to everything, DiscreteInferred loses to everything else
function merge_inferred(x, y)
    x === nothing && return y
    y === nothing && return x
    x isa Inferred && return y
    x == y || throw(ClockInferenceException("Cannot merge $x and $y"))
    x
end


# nothing cases
merge_domains(::Nothing, ::Nothing, args...) = nothing
merge_domains(x::Any, ::Nothing, args...) = x
merge_domains(::Nothing, x::Any, args...) = x

function merge_domains(above::A, below::B, ex) where {A <: TimeDomain, B <: TimeDomain}
    emsg = "TimeDomain inference error. In expression $ex, the domain was inferred to $above from the root of the tree but $below from the leaves."
    above == below && return above
    if above isa Continuous || below isa Continuous
        above === below || throw(ClockInferenceException(emsg))
        return Continuous()
    end
    if above isa Inferred || below isa Inferred
        return merge_inferred(above, below)
    end

    throw(ClockInferenceException(emsg))
end

