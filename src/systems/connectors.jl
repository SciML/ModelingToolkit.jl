struct Connection
    systems::Any
end
"""
    $(TYPEDSIGNATURES)

Given a hierarchical system with [`connect`](@ref) equations, expand the connection
equations and return the new system. `tol` is the tolerance for handling the singularities
in stream connection equations that happen when a flow variable approaches zero.
"""
function expand_connections(sys::AbstractSystem; tol = 1e-10)
    return flatten(sys)
end
