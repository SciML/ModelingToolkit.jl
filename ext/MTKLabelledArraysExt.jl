module MTKLabelledArraysExt

using ModelingToolkit, LabelledArrays

function ModelingToolkit.define_vars(u::Union{SLArray, LArray}, t)
    [_defvar(x)(t) for x in LabelledArrays.symnames(typeof(u))]
end

end
