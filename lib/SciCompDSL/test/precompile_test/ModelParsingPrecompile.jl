module ModelParsingPrecompile

using ModelingToolkitBase, DynamicQuantities
using ModelingToolkitBase: t
using SciCompDSL

@mtkmodel ModelWithComponentArray begin
    @constants begin
        k = 1, [description = "Default val of R"]
    end
    @parameters begin
        r(t)[1:3] = [k, k, k], [description = "Parameter array", unit = u"Ω"]
    end
end

end
