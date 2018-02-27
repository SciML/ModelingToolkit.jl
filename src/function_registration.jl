Base.:-(x::Expression,y::Expression) = Operation(-,Expression[x,y])
Base.:*(x::Expression,y::Expression) = Operation(*,Expression[x,y])
Base.:(==)(x::Expression,y::Expression) = Statement(x,y)
