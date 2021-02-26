"""
named(ex):

Used for naming an AbstractODESystem from an expression of the form `var = value`.
A keyword argument (name = :var) is added as the first keyword argument.

example:

`@named sys = ODESystem(eqs, args...; kwargs...)`

"""
macro named(ex)
    if !(ex isa Expr && ex.head == :(=))
      throw(ArgumentError("expression should be of the form `var = value`"))
    end
    pushfirst!(ex.args[2].args[2].args, :(name = :foo))
    :($(esc(ex.args[1])) = $(esc(ex.args[2])))
end
