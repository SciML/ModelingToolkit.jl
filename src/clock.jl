@data InferredClock begin
    Inferred
    InferredDiscrete(Int)
end
const InferredTimeDomain = InferredClock.Type
using .InferredClock: Inferred, InferredDiscrete
function InferredClock.InferredDiscrete()
    return InferredDiscrete(0)
end
Base.Broadcast.broadcastable(x::InferredTimeDomain) = Ref(x)
struct VariableTimeDomain end
Symbolics.option_to_metadata_type(::Val{:timedomain}) = VariableTimeDomain
is_concrete_time_domain(::TimeDomain) = true
is_concrete_time_domain(_) = false
""""""
function is_continuous_domain(x)
    issym(x) && return getmetadata(x, VariableTimeDomain, false) == ContinuousClock()
    !has_discrete_domain(x) && has_continuous_domain(x)
end
get_time_domain(_, x) = get_time_domain(x)
function get_time_domain(x)
    if iscall(x) && operation(x) isa Operator
        output_timedomain(x)
    else
        getmetadata(x, VariableTimeDomain, nothing)
    end
end
get_time_domain(x::Num) = get_time_domain(value(x))
has_time_domain(_, x) = has_time_domain(x)
""""""
function has_time_domain(x::SymbolicT)
    getmetadata(x, VariableTimeDomain, nothing) !== nothing
end
has_time_domain(x::Num) = has_time_domain(value(x))
has_time_domain(x) = false
input_timedomain(::Differential, arg = nothing) = InputTimeDomainElT[ContinuousClock()]
output_timedomain(::Differential, arg = nothing) = ContinuousClock()
""""""
function has_discrete_domain(x)
    issym(x) && return is_discrete_domain(x)
    hasshift(x) || hassample(x) || hashold(x)
end
""""""
function has_continuous_domain(x)
    issym(x) && return is_continuous_domain(x)
    hasderiv(x) || hassample(x) || hashold(x)
end
""""""
is_hybrid_domain(x) = has_discrete_domain(x) && has_continuous_domain(x)
""""""
function is_discrete_domain(x)
    if hasmetadata(x, VariableTimeDomain) || issym(x)
        return is_discrete_time_domain(getmetadata(x, VariableTimeDomain, false))
    end
    !has_discrete_domain(x) && has_continuous_domain(x)
end
sampletime(c) = Moshi.Match.@match c begin
    x::SciMLBase.AbstractClock => nothing
    PeriodicClock(dt) => dt
    _ => nothing
end
struct ClockInferenceException <: Exception
    msg::Any
end
function Base.showerror(io::IO, cie::ClockInferenceException)
    print(io, "ClockInferenceException: ", cie.msg)
end
struct IntegerSequence end
