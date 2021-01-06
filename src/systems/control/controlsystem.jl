abstract type AbstractControlSystem <: AbstractSystem end

function namespace_controls(sys::AbstractSystem)
    [rename(x,renamespace(sys.name,x.name)) for x in controls(sys)]
end

function controls(sys::AbstractControlSystem,args...)
    name = last(args)
    extra_names = reduce(Symbol,[Symbol(:₊,x.name) for x in args[1:end-1]])
    newname = renamespace(extra_names,name)
    rename(x,renamespace(sys.name,newname))(sys.iv())
end

function controls(sys::AbstractControlSystem,name::Symbol)
    x = sys.controls[findfirst(x->x.name==name,sys.ps)]
    rename(x,renamespace(sys.name,x.name))()
end

controls(sys::AbstractControlSystem) = isempty(sys.systems) ? sys.controls : [sys.controls;reduce(vcat,namespace_controls.(sys.systems))]

"""
$(TYPEDEF)

A system describing an optimal control problem. This contains a loss function
and ordinary differential equations with control variables that describe the
dynamics.

# Fields
$(FIELDS)

# Example

```
using ModelingToolkit

@variables t x(t) v(t) u(t)
@derivatives D'~t

loss = (4-x)^2 + 2v^2 + u^2
eqs = [
    D(x) ~ v
    D(v) ~ u^3
]

sys = ControlSystem(loss,eqs,t,[x,v],[u],[])
```
"""
struct ControlSystem <: AbstractControlSystem
    """The Loss function"""
    loss::Term
    """The ODEs defining the system."""
    eqs::Vector{Equation}
    """Independent variable."""
    iv::Sym
    """Dependent (state) variables."""
    states::Vector
    """Control variables."""
    controls::Vector
    """Parameter variables."""
    ps::Vector
    pins::Vector
    observed::Vector{Equation}
    """
    Name: the name of the system
    """
    name::Symbol
    """
    systems: The internal systems
    """
    systems::Vector{ControlSystem}
end

function ControlSystem(loss, deqs::AbstractVector{<:Equation}, iv, dvs, controls, ps;
                   pins = [],
                   observed = [],
                   systems = ODESystem[],
                   name=gensym(:ControlSystem))
    iv′ = value(iv)
    dvs′ = value.(dvs)
    controls′ = value.(controls)
    ps′ = value.(ps)
    ControlSystem(value(loss), deqs, iv′, dvs′, controls′,
                  ps′, pins, observed, name, systems)
end

struct ControlToExpr
    sys::AbstractControlSystem
    states::Vector
    controls::Vector
end
ControlToExpr(@nospecialize(sys)) = ControlToExpr(sys,states(sys),controls(sys))
function (f::ControlToExpr)(O::Term)
    res = if isa(operation(O), Sym)
        # normal variables and control variables
        (any(isequal(O), f.states) || any(isequal(O), f.controls)) && return tosymbol(O)
        build_expr(:call, Any[operation(O).name; f.(arguments(O))])
    else
        build_expr(:call, Any[Symbol(operation(O)); f.(arguments(O))])
    end
end
(f::ControlToExpr)(x::Sym) = x.name
(f::ControlToExpr)(x) = x

function constructRadauIIA5(T::Type = Float64)
  sq6 = sqrt(6)
  A = [11//45-7sq6/360 37//225-169sq6/1800 -2//225+sq6/75
       37//225+169sq6/1800 11//45+7sq6/360 -2//225-sq6/75
       4//9-sq6/36 4//9+sq6/36 1//9]
  c = [2//5-sq6/10;2/5+sq6/10;1]
  α = [4//9-sq6/36;4//9+sq6/36;1//9]
  A = map(T,A)
  α = map(T,α)
  c = map(T,c)
  return(DiffEqBase.ImplicitRKTableau(A,c,α,5))
end


"""
```julia
runge_kutta_discretize(sys::ControlSystem,dt,tspan;
                       tab = ModelingToolkit.constructRadauIIA5())
```

Transforms a nonlinear optimal control problem into a constrained
`OptimizationProblem` according to a Runge-Kutta tableau that describes
a collocation method. Requires a fixed `dt` over a given `timespan`.
Defaults to using the 5th order RadauIIA tableau, and altnerative tableaus
can be specified using the SciML tableau style.
"""
function runge_kutta_discretize(sys::ControlSystem,dt,tspan;
                       tab = ModelingToolkit.constructRadauIIA5())
    n = length(tspan[1]:dt:tspan[2]) - 1
    m = length(tab.α)

    f = @RuntimeGeneratedFunction(build_function([x.rhs for x in equations(sys)],sys.states,sys.controls,sys.ps,sys.iv,conv = ModelingToolkit.ControlToExpr(sys))[1])
    L = @RuntimeGeneratedFunction(build_function(sys.loss,sys.states,sys.controls,sys.ps,sys.iv,conv = ModelingToolkit.ControlToExpr(sys)))

    var(n, i...) = var(nameof(n), i...)
    var(n::Symbol, i...) = Sym{FnType{Tuple{symtype(sys.iv)}, Real}}(nameof(Variable(n, i...)))
    # Expand out all of the variables in time and by stages
    timed_vars = [[var(operation(x),i)(sys.iv) for i in 1:n+1] for x in states(sys)]
    k_vars = [[var(Symbol(:ᵏ,nameof(operation(x))),i,j)(sys.iv) for i in 1:m, j in 1:n] for x in states(sys)]
    states_timeseries = [getindex.(timed_vars,j) for j in 1:n+1]
    k_timeseries = [[Num.(getindex.(k_vars,i,j)) for i in 1:m] for j in 1:n]
    control_timeseries = [[[var(operation(x),i,j)(sys.iv) for x in controls(sys)] for i in 1:m] for j in 1:n]
    ps = parameters(sys)
    iv = sys.iv

    # Calculate all of the update and stage equations
    mult = [tab.A * k_timeseries[i] for i in 1:n]
    tmps = [[states_timeseries[i] .+ mult[i][j] for j in 1:m] for i in 1:n]

    bs = [states_timeseries[i] .+ dt .* reduce(+, tab.α .* k_timeseries[i],dims=1)[1] for i in 1:n]
    updates = reduce(vcat,[states_timeseries[i+1] .~ bs[i] for i in 1:n])

    df = [[dt .* Base.invokelatest(f,tmps[j][i],control_timeseries[j][i],ps,iv) for i in 1:m] for j in 1:n]
    stages = reduce(vcat,[k_timeseries[i][j] .~ df[i][j] for i in 1:n for j in 1:m])

    # Enforce equalities in the controls
    control_equality = reduce(vcat,[control_timeseries[i][end] .~ control_timeseries[i+1][1] for i in 1:n-1])

    # Create the loss function
    losses = [Base.invokelatest(L,states_timeseries[i],control_timeseries[i][1],ps,iv) for i in 1:n]
    losses = vcat(losses,[Base.invokelatest(L,states_timeseries[n+1],control_timeseries[n][end],ps,iv)])

    # Calculate final pieces
    equalities = vcat(stages,updates,control_equality)
    opt_states = vcat(reduce(vcat,reduce(vcat,states_timeseries)),reduce(vcat,reduce(vcat,k_timeseries)),reduce(vcat,reduce(vcat,control_timeseries)))

    OptimizationSystem(reduce(+,losses, init=0),opt_states,ps,equality_constraints = equalities)
end
