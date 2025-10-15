struct Connection
    systems::Any
end
""""""
function expand_connections(sys::AbstractSystem; tol = 1e-10)
    return flatten(sys)
end
