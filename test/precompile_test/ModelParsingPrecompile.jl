module ModelParsingPrecompile

using ModelingToolkit, Unitful

@mtkmodel ModelWithComponentArray begin
    @parameters begin
        R(t)[1:3] = 1, [description = "Parameter array", unit = u"Ω"]
    end
end

end
