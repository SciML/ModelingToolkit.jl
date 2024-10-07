module MTKLabelledArraysExt

using ModelingToolkit, LabelledArrays
using ModelingToolkit: _defvar, toparam, variable, varnames_length_check
function ModelingToolkit.define_vars(u::Union{SLArray, LArray}, t)
    [ModelingToolkit._defvar(x)(t) for x in LabelledArrays.symnames(typeof(u))]
end

function ModelingToolkit.define_params(p::Union{SLArray, LArray}, names = nothing)
    if names === nothing
        [toparam(variable(x)) for x in LabelledArrays.symnames(typeof(p))]
    else
        varnames_length_check(p, names)
        [toparam(variable(names[i])) for i in eachindex(p)]
    end
end

end
