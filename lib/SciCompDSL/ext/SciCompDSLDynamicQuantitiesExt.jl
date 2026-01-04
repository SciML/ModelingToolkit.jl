module SciCompDSLDynamicQuantitiesExt

import DynamicQuantities
const DQ = DynamicQuantities

using ModelingToolkitBase, Symbolics
using ModelingToolkitBase: VariableUnit, setdefault
using SciCompDSL
using SciCompDSL: convert_units, NoValue, NO_VALUE

import ModelingToolkitBase as MTK

function SciCompDSL.convert_units(varunits::DynamicQuantities.Quantity, value)
    return DynamicQuantities.ustrip(
        DynamicQuantities.uconvert(
            DynamicQuantities.SymbolicUnits.as_quantity(varunits), value
        )
    )
end

SciCompDSL.convert_units(::DynamicQuantities.Quantity, value::NoValue) = NO_VALUE

function SciCompDSL.convert_units(
        varunits::DynamicQuantities.Quantity, value::AbstractArray{T}
    ) where {T}
    return DynamicQuantities.ustrip.(
        DynamicQuantities.uconvert.(
            DynamicQuantities.SymbolicUnits.as_quantity(varunits), value
        )
    )
end

SciCompDSL.convert_units(::DynamicQuantities.Quantity, value::AbstractArray{Num}) = value

SciCompDSL.convert_units(::DynamicQuantities.Quantity, value::Num) = value

function SciCompDSL.__generate_variable_with_unit(metadata_with_exprs, name, vv, def)
    unit = metadata_with_exprs[VariableUnit]
    return quote
        $name = if $name === $(NO_VALUE)
            $setdefault($vv, $def)
        else
            try
                $setdefault($vv, $convert_units($unit, $name))
            catch e
                if isa(e, $(DynamicQuantities.DimensionError))
                    error("Unable to convert units for \'" * string(:($$vv)) * "\'")
                elseif isa(e, MethodError)
                    error(
                        "No or invalid units provided for \'" * string(:($$vv)) *
                            "\'"
                    )
                else
                    rethrow(e)
                end
            end
        end
    end
end

end
