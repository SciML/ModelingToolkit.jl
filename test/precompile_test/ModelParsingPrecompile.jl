module ModelParsingPrecompile

using ModelingToolkit, Unitful
using ModelingToolkit: t

@mtkmodel ModelWithComponentArray begin
    @constants begin
        k = 1, [description = "Default val of R"]
    end
    @parameters begin
        r(t)[1:3] = k, [description = "Parameter array", unit = u"Ω"]
    end
end

end
