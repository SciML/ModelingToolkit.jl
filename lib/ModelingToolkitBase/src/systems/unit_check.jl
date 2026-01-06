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
global t::Union{PleaseImportDynamicQuantities, Num} = PleaseImportDynamicQuantities()

function Base.show(io::IO, ::PleaseImportDynamicQuantities)
    return __import_dynamic_quantities()
end

function __import_dynamic_quantities(_...)
    error(
        """
        Please import DynamicQuantites.jl to use this `t` and `D`.
        """
    )
end
global D::Union{typeof(__import_dynamic_quantities), Differential} = __import_dynamic_quantities
