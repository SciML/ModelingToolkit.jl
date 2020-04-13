"""
calculate_tgrad(sys::AbstractSystem)

Calculate the time gradient of a system.

Returns a vector of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_tgrad end

"""
calculate_grad(sys::AbstractSystem)

Calculate the gradient of a scalar system.

Returns a vector of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_grad end

"""
calculate_jacobian(sys::AbstractSystem)

Calculate the jacobian matrix of a system.

Returns a matrix of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_jacobian end

"""
calculate_factorized_W(sys::AbstractSystem)

Calculate the factorized W-matrix of a system.

Returns a matrix of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_factorized_W end

"""
calculate_hessian(sys::AbstractSystem)

Calculate the hessian matrix of a scalar system.

Returns a matrix of [`Expression`](@ref) instances. The result from the first
call will be cached in the system object.
"""
function calculate_hessian end

"""
generate_tgrad(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)

Generates a function for the time gradient of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_tgrad end

"""
generate_grad(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)

Generates a function for the gradient of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_grad end

"""
generate_jacobian(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; sparse = false, kwargs...)

Generates a function for the jacobian matrix matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_jacobian end

"""
generate_factorized_W(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; sparse = false, kwargs...)

Generates a function for the factorized W-matrix matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_factorized_W end

"""
generate_hessian(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; sparse = false, kwargs...)

Generates a function for the hessian matrix matrix of a system. Extra arguments control
the arguments to the internal [`build_function`](@ref) call.
"""
function generate_hessian end

"""
generate_function(sys::AbstractSystem, dvs = states(sys), ps = parameters(sys), expression = Val{true}; kwargs...)

Generate a function to evaluate the system's equations.
"""
function generate_function end

function Base.getproperty(sys::AbstractSystem, name::Symbol)
    if name ∈ fieldnames(typeof(sys))
        return getfield(sys,name)
    elseif !isempty(sys.systems)
        i = findfirst(x->x.name==name,sys.systems)
        if i !== nothing
            return rename(sys.systems[i],renamespace(sys.name,name))
        end
    end
    i = findfirst(x->x.name==name,sys.states)
    if i !== nothing
        x = rename(sys.states[i],renamespace(sys.name,name))
        if :iv ∈ fieldnames(typeof(sys))
            return x(getfield(sys,:iv)())
        else
            return x()
        end
    end
    if :ps ∈ fieldnames(typeof(sys))
        i = findfirst(x->x.name==name,sys.ps)
        if i !== nothing
            return rename(sys.ps[i],renamespace(sys.name,name))()
        end
    end
    throw(error("Variable $name does not exist"))
end

renamespace(namespace,name) = Symbol(namespace,:₊,name)

function namespace_variables(sys::AbstractSystem)
    [rename(x,renamespace(sys.name,x.name)) for x in states(sys)]
end

function namespace_parameters(sys::AbstractSystem)
    [rename(x,renamespace(sys.name,x.name)) for x in parameters(sys)]
end

namespace_equations(sys::AbstractSystem) = namespace_equation.(equations(sys),sys.name,sys.iv.name)

function namespace_equation(eq::Equation,name,ivname)
    _lhs = namespace_operation(eq.lhs,name,ivname)
    _rhs = namespace_operation(eq.rhs,name,ivname)
    _lhs ~ _rhs
end

function namespace_operation(O::Operation,name,ivname)
    if O.op isa Variable && O.op.name != ivname
        Operation(rename(O.op,renamespace(name,O.op.name)),namespace_operation.(O.args,name,ivname))
    else
        Operation(O.op,namespace_operation.(O.args,name,ivname))
    end
end
namespace_operation(O::Constant,name,ivname) = O

independent_variable(sys::AbstractSystem) = sys.iv
states(sys::AbstractSystem) = isempty(sys.systems) ? sys.states : [sys.states;reduce(vcat,namespace_variables.(sys.systems))]
parameters(sys::AbstractSystem) = isempty(sys.systems) ? sys.ps : [sys.ps;reduce(vcat,namespace_parameters.(sys.systems))]

function equations(sys::AbstractSystem)
    isempty(sys.systems) ? sys.eqs : [sys.eqs;reduce(vcat,namespace_equations.(sys.systems))]
end

function states(sys::AbstractSystem,name::Symbol)
    x = sys.states[findfirst(x->x.name==name,sys.states)]
    rename(x,renamespace(sys.name,x.name))(sys.iv())
end

function parameters(sys::AbstractSystem,name::Symbol)
    x = sys.ps[findfirst(x->x.name==name,sys.ps)]
    rename(x,renamespace(sys.name,x.name))()
end

function states(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(Symbol,[Symbol(:₊,x.name) for x in args[1:end-1]])
    newname = renamespace(extra_names,name)
    rename(x,renamespace(sys.name,newname))(sys.iv())
end

function parameters(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(Symbol,[Symbol(:₊,x.name) for x in args[1:end-1]])
    newname = renamespace(extra_names,name)
    rename(x,renamespace(sys.name,newname))()
end

struct AbstractSysToExpr
    sys::AbstractSystem
    states::Vector{Variable}
end
AbstractSysToExpr(sys) = AbstractSysToExpr(sys,states(sys))
function (f::AbstractSysToExpr)(O::Operation)
    any(isequal(O), f.states) && return O.op.name  # variables
    if isa(O.op, Variable)
        isempty(O.args) && return O.op.name  # 0-ary parameters
        return build_expr(:call, Any[O.op.name; f.(O.args)])
    end
    return build_expr(:call, Any[O.op; f.(O.args)])
end
(f::AbstractSysToExpr)(x) = convert(Expr, x)
