module MTKFunctionPropertiesExt

import ModelingToolkitBase as MTKBase
import FunctionProperties

# Integer indexing into an `MTKParameters` buffer selects one of the parameter portions by index.
# The selecting branch — `getindex(::MTKParameters, ::Int)` is a `@generated` `if idx == k` chain —
# is value-independent: every real call site in generated right-hand-side code passes a literal
# index that constant-folds the branch away, but `hasbranching`'s recursion only sees the widened
# `Int`, so it reports the branch. Mark this plumbing as a leaf so `hasbranching` reflects the
# user's RHS (e.g. a registered branching activation) rather than parameter-container indexing.
FunctionProperties.is_leaf_sig(
    ::Type{<:Tuple{typeof(Base.getindex), <:MTKBase.MTKParameters, Vararg}}
) = true

end
