using SciCompDSL
using Base.Test

@DVar x y z
@IVar t

struct __MyType__ end
Base.:~(::__MyType__,::Number) = 2

x ~ 2x + y
