using Symbolics: Differential, Difference
using Symbolics: hasderiv, hasdiff

abstract type TimeDomain end

struct Continuous <: TimeDomain end
struct Discrete <: TimeDomain end


for op in [Differential, Difference]
    @eval input_timedomain(::$op) = Continuous()
    @eval output_timedomain(::$op) = Continuous()
end


"""
    input_timedomain(op::Operator)

Return the time-domain type (`Continuous()` or `Discrete()`) that `op` operates on. 
"""
input_timedomain(::Shift) = Discrete()

"""
    output_timedomain(op::Operator)

Return the time-domain type (`Continuous()` or `Discrete()`) that `op` results in. 
"""
output_timedomain(::Shift) = Discrete()

input_timedomain(::Sample) = Continuous()
output_timedomain(::Sample) = Discrete()

input_timedomain(::Hold) = Discrete()
output_timedomain(::Hold) = Continuous()


"""
    sampletime(op::Operator)

Return the sample time for operators that are clocked. Returns `nothing` if the operator does not have a sample time.
"""
sampletime(op::Operator) = hasfield(typeof(op), :dt) ? op.dt : nothing

"""
    has_discrete_domain(x)

true if `x` contains discrete signals (`x` may or may not contain continuous-domain signals). `x` may be an expression or equation.
See also [`is_discrete_domain`](@ref)
"""
has_discrete_domain(x) = hasshift(x) || hassample(x) || hashold(x)

"""
    has_continuous_domain(x)

true if `x` contains continuous signals (`x` may or may not contain discrete-domain signals). `x` may be an expression or equation.
See also [`is_continuous_domain`](@ref)
"""
has_continuous_domain(x) = hasderiv(x) || hasdiff(x) || hassample(x) || hashold(x)


"""
    transitions_timedomain(x)

true if `x` contains both discrete and continuous-domain signals. `x` may be an expression or equation.
"""
transitions_timedomain(x) = has_discrete_domain(x) && has_continuous_domain(x)

"""
    is_discrete_domain(x)

true if `x` contains only discrete-domain signals.
See also [`has_discrete_domain`](@ref)
"""
is_discrete_domain(x) = has_discrete_domain(x) && !has_continuous_domain(x)

"""
    is_continuous_domain(x)

true if `x` contains only continuous-domain signals.
See also [`has_continuous_domain`](@ref)
"""
is_continuous_domain(x) = !has_discrete_domain(x) && has_continuous_domain(x)