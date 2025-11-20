module MTKLabelledArraysExt

using ModelingToolkitBase, LabelledArrays
using ModelingToolkitBase: _defvar, toparam, variable, varnames_length_check
function ModelingToolkitBase.define_vars(u::Union{SLArray, LArray}, t)
    [ModelingToolkitBase._defvar(x)(t) for x in LabelledArrays.symnames(typeof(u))]
end

function ModelingToolkitBase.define_params(p::Union{SLArray, LArray}, t, names = nothing)
    if names === nothing
        [toparam(variable(x)) for x in LabelledArrays.symnames(typeof(p))]
    else
        varnames_length_check(p, names)
        [toparam(variable(names[i])) for i in eachindex(p)]
    end
end

end
