module ModelingToolkitStandardLibrary
import Symbolics: unwrap

"""
  @symcheck J > 0 || throw(ArgumentError("Expected `J` to be positive"))
  
Omits the check expression if the argument `J` is symbolic.
"""
macro symcheck(ex)
    ex.args[1].head === :call ||
        error("Expected an expression on the form sym > val || error()")
    sym = ex.args[1].args[2]
    quote
        _issymbolic(x) = !(unwrap(x) isa Real)
        _issymbolic($(esc(sym))) || ($(esc(ex)))
    end
end

include("Blocks/Blocks.jl")
include("Mechanical/Mechanical.jl")
include("Thermal/Thermal.jl")
include("Electrical/Electrical.jl")
include("Magnetic/Magnetic.jl")
include("Hydraulic/Hydraulic.jl")

end
