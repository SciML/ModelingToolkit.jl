module ModelingToolkitUnitfulExt

using ModelingToolkit
using Unitful
using Symbolics: Symbolic, value
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
function MT.equivalent(x::Unitful.Unitlike, y)
    if y isa MT.DQ.AbstractQuantity
        # Handle dimensionless case
        if Unitful.dimension(x) == Unitful.NoDims && MT.DQ.is_unitless(y)
            return true
        end
        # For mixed unit systems, we can't reliably compare
        # This would require a full dimensional analysis system
        return false
    else
        try
            return isequal(1 * x, y)
        catch
            return false
        end
    end
end

function MT.equivalent(x, y::Unitful.Unitlike)
    if x isa MT.DQ.AbstractQuantity
        # Handle dimensionless case
        if Unitful.dimension(y) == Unitful.NoDims && MT.DQ.is_unitless(x)
            return true
        end
        # For mixed unit systems, we can't reliably compare
        return false
    else
        try
            return isequal(x, 1 * y)
        catch
            return false
        end
    end
end

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

# Extension loaded - all Unitful-specific functionality is now available

end # module