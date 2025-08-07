module ModelingToolkitUnitfulExt

using ModelingToolkit, Symbolics, SciMLBase, Unitful, RecursiveArrayTools
using ModelingToolkit: ValidationError, Connection, instream, JumpType, VariableUnit,
                         get_systems, Conditional, Comparison, Integral, Differential
using JumpProcesses: MassActionJump, ConstantRateJump, VariableRateJump
using Symbolics: Symbolic, value, issym, isadd, ismul, ispow, iscall, operation, arguments, getmetadata

using Unitful
using SciMLBase

# Import necessary types and functions from ModelingToolkit
import ModelingToolkit: ValidationError, _get_unittype, get_unit, screen_unit, 
                        equivalent, _is_dimension_error, convert_units, check_units

const MT = ModelingToolkit

# Add Unitful-specific unit type detection
function MT._get_unittype(u::Unitful.Unitlike)
    return Val(:Unitful)
end

# Base operations for mixing Symbolic and Unitful
Base.:*(x::Union{MT.Num, Symbolic}, y::Unitful.AbstractQuantity) = x * y
Base.:/(x::Union{MT.Num, Symbolic}, y::Unitful.AbstractQuantity) = x / y

# Unitful-specific get_unit method
function MT.get_unit(x::Unitful.Quantity)
    return screen_unit(Unitful.unit(x))
end

# Unitful-specific screen_unit method
function MT.screen_unit(result::Unitful.Unitlike)
    result isa Unitful.ScalarUnits ||
        throw(ValidationError("Non-scalar units such as $result are not supported. Use a scalar unit instead."))
    result == Unitful.u"°" &&
        throw(ValidationError("Degrees are not supported. Use radians instead."))
    return result
end

# Unitful-specific equivalence check
function MT.equivalent(x::Unitful.Unitlike, y::Unitful.Unitlike)
    return isequal(1 * x, 1 * y)
end

# Mixed equivalence checks
MT.equivalent(x::Unitful.Unitlike, y) = isequal(1 * x, y)
MT.equivalent(x, y::Unitful.Unitlike) = isequal(x, 1 * y)

# The safe_get_unit function stays in the main package and already handles DQ.DimensionError
# We just need to make sure it can handle Unitful.DimensionError too
# This will be handled by the main function's MethodError catch

# Unitful-specific dimension error detection for model parsing
MT._is_dimension_error(e::Unitful.DimensionError) = true

# Unitful-specific convert_units methods for model parsing
function MT.convert_units(varunits::Unitful.FreeUnits, value)
    Unitful.ustrip(varunits, value)
end

MT.convert_units(::Unitful.FreeUnits, value::MT.NoValue) = MT.NO_VALUE

function MT.convert_units(varunits::Unitful.FreeUnits, value::AbstractArray{T}) where {T}
    Unitful.ustrip.(varunits, value)
end

MT.convert_units(::Unitful.FreeUnits, value::MT.Num) = value

# Unitful-specific check_units method
function MT.check_units(::Val{:Unitful}, eqs...)
    # Use the main package's validate function
    MT.validate(eqs...) ||
        throw(ValidationError("Some equations had invalid units. See warnings for details."))
end

# Define Unitful time variables (moved from main module)
const t_unitful = let
    MT.only(MT.@independent_variables t [unit = Unitful.u"s"])
end
const D_unitful = MT.Differential(t_unitful)

Base.:*(x::Union{MT.Num, Symbolic}, y::Unitful.AbstractQuantity) = x * y
Base.:/(x::Union{MT.Num, Symbolic}, y::Unitful.AbstractQuantity) = x / y

"""
Throw exception on invalid unit types, otherwise return argument.
"""
function screen_unit(result)
    result isa Unitful.Unitlike ||
        throw(ValidationError("Unit must be a subtype of Unitful.Unitlike, not $(typeof(result))."))
    result isa Unitful.ScalarUnits ||
        throw(ValidationError("Non-scalar units such as $result are not supported. Use a scalar unit instead."))
    result == Unitful.u"°" &&
        throw(ValidationError("Degrees are not supported. Use radians instead."))
    result
end

"""
Test unit equivalence.
"""
equivalent(x, y) = isequal(1 * x, 1 * y)
const unitless = Unitful.unit(1)

"""
Find the unit of a symbolic item.
"""
get_unit(x::Real) = unitless
get_unit(x::Unitful.Quantity) = screen_unit(Unitful.unit(x))
get_unit(x::AbstractArray) = map(get_unit, x)
get_unit(x::MT.Num) = get_unit(value(x))
function get_unit(x::Union{Symbolics.ArrayOp, Symbolics.Arr, Symbolics.CallWithMetadata})
    get_literal_unit(x)
end
get_unit(op::Differential, args) = get_unit(args[1]) / get_unit(op.x)
get_unit(op::typeof(getindex), args) = get_unit(args[1])
get_unit(x::SciMLBase.NullParameters) = unitless
get_unit(op::typeof(instream), args) = get_unit(args[1])

get_literal_unit(x) = screen_unit(getmetadata(x, VariableUnit, unitless))

function get_unit(op, args) # Fallback
    result = op(1 .* get_unit.(args)...)
    try
        Unitful.unit(result)
    catch
        throw(ValidationError("Unable to get unit for operation $op with arguments $args."))
    end
end

function get_unit(op::Integral, args)
    unit = 1
    if op.domain.variables isa Vector
        for u in op.domain.variables
            unit *= get_unit(u)
        end
    else
        unit *= get_unit(op.domain.variables)
    end
    return get_unit(args[1]) * unit
end

function get_unit(op::Conditional, args)
    terms = get_unit.(args)
    terms[1] == unitless ||
        throw(ValidationError(", in $op, [$(terms[1])] is not dimensionless."))
    equivalent(terms[2], terms[3]) ||
        throw(ValidationError(", in $op, units [$(terms[2])] and [$(terms[3])] do not match."))
    return terms[2]
end

function get_unit(op::typeof(Symbolics._mapreduce), args)
    if args[2] == +
        get_unit(args[3])
    else
        throw(ValidationError("Unsupported array operation $op"))
    end
end

function get_unit(op::Comparison, args)
    terms = get_unit.(args)
    equivalent(terms[1], terms[2]) ||
        throw(ValidationError(", in comparison $op, units [$(terms[1])] and [$(terms[2])] do not match."))
    return unitless
end

function get_unit(x::Symbolic)
    if issym(x)
        get_literal_unit(x)
    elseif isadd(x)
        terms = get_unit.(arguments(x))
        firstunit = terms[1]
        for other in terms[2:end]
            termlist = join(map(repr, terms), ", ")
            equivalent(other, firstunit) ||
                throw(ValidationError(", in sum $x, units [$termlist] do not match."))
        end
        return firstunit
    elseif ispow(x)
        pargs = arguments(x)
        base, expon = get_unit.(pargs)
        @assert expon isa Unitful.DimensionlessUnits
        if base == unitless
            unitless
        else
            pargs[2] isa Number ? base^pargs[2] : (1 * base)^pargs[2]
        end
    elseif iscall(x)
        op = operation(x)
        if issym(op) || (iscall(op) && iscall(operation(op))) # Dependent variables, not function calls
            return screen_unit(getmetadata(x, VariableUnit, unitless)) # Like x(t) or x[i]
        elseif iscall(op) && !iscall(operation(op))
            gp = getmetadata(x, Symbolics.GetindexParent, nothing) # Like x[1](t)
            return screen_unit(getmetadata(gp, VariableUnit, unitless))
        end  # Actual function calls:
        args = arguments(x)
        return get_unit(op, args)
    else # This function should only be reached by Terms, for which `iscall` is true
        throw(ArgumentError("Unsupported value $x."))
    end
end

end # module UnitfulUnitCheck
