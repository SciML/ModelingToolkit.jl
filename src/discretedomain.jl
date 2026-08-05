using Symbolics: Num, value, recursive_hasoperator
using SymbolicUtils: Operator, Term

MTKBase.ShiftIndex() = MTKBase.ShiftIndex(MTKTearing.Inferred())
MTKBase.Sample() = MTKBase.Sample(MTKTearing.InferredDiscrete())
