# Unit checking utilities
check_units(_...) = true
__get_unit_type(_...) = nothing
get_unit(_...) = nothing
validate(_...) = nothing
_validate(_...) = true
struct ValidationError <: Exception
    message::String
end

struct PleaseImportDynamicQuantities end

"""
    t

Default unitful independent variable used with DynamicQuantities.

Import DynamicQuantities before using this binding. For unitless models, use
[`t_nounits`](@ref) instead.

# Examples

```julia
using DynamicQuantities, ModelingToolkitBase

t = ModelingToolkitBase.t
```
"""
global t::Union{PleaseImportDynamicQuantities, Num} = PleaseImportDynamicQuantities()

function Base.show(io::IO, ::PleaseImportDynamicQuantities)
    return __import_dynamic_quantities()
end

function __import_dynamic_quantities(_...)
    error(
        """
        Please import DynamicQuantities.jl to use this `t` and `D`.
        """
    )
end

"""
    D

Default unitful differential operator with respect to [`t`](@ref).

Import DynamicQuantities before using this binding. For unitless models, use
[`D_nounits`](@ref) instead.

# Examples

```julia
using DynamicQuantities, ModelingToolkitBase

D = ModelingToolkitBase.D
```
"""
global D::Union{typeof(__import_dynamic_quantities), Differential} = __import_dynamic_quantities
