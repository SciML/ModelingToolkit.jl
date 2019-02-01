using ModelingToolkit
using Test

@Param t()
@Unknown x(t) y(t) z(t)

struct __MyType__ end
Base.:~(::__MyType__,::Number) = 2
Base.:~(::__MyType__,::Any) = 2
x ~ 2x + y
x ~ 2
