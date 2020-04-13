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

renamespace(namespace,name) = Symbol(string(namespace)*"₊"*string(name))

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
    Variable(Symbol(string(sys.name)*"₊"*string(x.name)))(sys.iv())
end

function parameters(sys::AbstractSystem,name::Symbol)
    x = sys.ps[findfirst(x->x.name==name,sys.ps)]
    Variable(Symbol(string(sys.name)*"₊"*string(x.name)))(sys.iv())
end

function states(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(*,["₊$(x.name)" for x in args[1:end-1]])
    Variable(Symbol(string(sys.name)*extra_names*"₊"*string(name)))(sys.iv())
end

function parameters(sys::AbstractSystem,args...)
    name = last(args)
    extra_names = reduce(*,["₊$(x.name)" for x in args[1:end-1]])
    Variable(Symbol(string(sys.name)*extra_names*"₊"*string(name)))(sys.iv())
end
