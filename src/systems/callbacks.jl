#################################### system operations #####################################
get_continuous_events(sys::AbstractSystem) = Equation[]
get_continuous_events(sys::AbstractODESystem) = getfield(sys, :continuous_events)
has_continuous_events(sys::AbstractSystem) = isdefined(sys, :continuous_events)

has_discrete_events(sys::AbstractSystem) = isdefined(sys, :discrete_events)
function get_discrete_events(sys::AbstractSystem)
    has_discrete_events(sys) || return SymbolicDiscreteCallback[]
    getfield(sys, :discrete_events)
end

#################################### continuous events #####################################

const NULL_AFFECT = Equation[]
struct SymbolicContinuousCallback
    eqs::Vector{Equation}
    affect::Vector{Equation}
    function SymbolicContinuousCallback(eqs::Vector{Equation}, affect = NULL_AFFECT)
        new(eqs, affect)
    end # Default affect to nothing
end

function Base.:(==)(e1::SymbolicContinuousCallback, e2::SymbolicContinuousCallback)
    isequal(e1.eqs, e2.eqs) && isequal(e1.affect, e2.affect)
end
Base.isempty(cb::SymbolicContinuousCallback) = isempty(cb.eqs)
function Base.hash(cb::SymbolicContinuousCallback, s::UInt)
    s = foldr(hash, cb.eqs, init = s)
    foldr(hash, cb.affect, init = s)
end

to_equation_vector(eq::Equation) = [eq]
to_equation_vector(eqs::Vector{Equation}) = eqs
function to_equation_vector(eqs::Vector{Any})
    isempty(eqs) || error("This should never happen")
    Equation[]
end

function SymbolicContinuousCallback(args...)
    SymbolicContinuousCallback(to_equation_vector.(args)...)
end # wrap eq in vector
SymbolicContinuousCallback(p::Pair) = SymbolicContinuousCallback(p[1], p[2])
SymbolicContinuousCallback(cb::SymbolicContinuousCallback) = cb # passthrough

SymbolicContinuousCallbacks(cb::SymbolicContinuousCallback) = [cb]
SymbolicContinuousCallbacks(cbs::Vector{<:SymbolicContinuousCallback}) = cbs
SymbolicContinuousCallbacks(cbs::Vector) = SymbolicContinuousCallback.(cbs)
function SymbolicContinuousCallbacks(ve::Vector{Equation})
    SymbolicContinuousCallbacks(SymbolicContinuousCallback(ve))
end
function SymbolicContinuousCallbacks(others)
    SymbolicContinuousCallbacks(SymbolicContinuousCallback(others))
end
SymbolicContinuousCallbacks(::Nothing) = SymbolicContinuousCallbacks(Equation[])

equations(cb::SymbolicContinuousCallback) = cb.eqs
function equations(cbs::Vector{<:SymbolicContinuousCallback})
    reduce(vcat, [equations(cb) for cb in cbs])
end
affect_equations(cb::SymbolicContinuousCallback) = cb.affect
function affect_equations(cbs::Vector{SymbolicContinuousCallback})
    reduce(vcat, [affect_equations(cb) for cb in cbs])
end
function namespace_equation(cb::SymbolicContinuousCallback, s)::SymbolicContinuousCallback
    SymbolicContinuousCallback(namespace_equation.(equations(cb), (s,)),
                               namespace_equation.(affect_equations(cb), (s,)))
end

function continuous_events(sys::AbstractSystem)
    obs = get_continuous_events(sys)
    filter(!isempty, obs)
    systems = get_systems(sys)
    cbs = [obs;
           reduce(vcat,
                  (map(o -> namespace_equation(o, s), continuous_events(s))
                   for s in systems),
                  init = SymbolicContinuousCallback[])]
    filter(!isempty, cbs)
end

#################################### continuous events #####################################

struct SymbolicDiscreteCallback
    condition
    affects::Vector{Equation}
    function SymbolicDiscreteCallback(condition, affects = NULL_AFFECT)
        c = value(scalarize(condition))
        a = scalarize(affects)
        new(c, a)
    end # Default affect to nothing
end

SymbolicDiscreteCallback(p::Pair) = SymbolicDiscreteCallback(p[1], p[2])
SymbolicDiscreteCallback(cb::SymbolicDiscreteCallback) = cb # passthrough

function Base.show(io::IO, db::SymbolicDiscreteCallback)
    println(io, "condition: ", db.condition)
    println(io, "affects:")
    for affect in db.affects
        println(io, "  ", affect)
    end
end

function Base.:(==)(e1::SymbolicDiscreteCallback, e2::SymbolicDiscreteCallback)
    isequal(e1.condition, e2.condition) && isequal(e1.affects, e2.affects)
end
function Base.hash(cb::SymbolicDiscreteCallback, s::UInt)
    s = foldr(hash, cb.condition, init = s)
    foldr(hash, cb.affects, init = s)
end

condition(cb::SymbolicDiscreteCallback) = cb.condition
function conditions(cbs::Vector{<:SymbolicDiscreteCallback})
    reduce(vcat, condition(cb) for cb in cbs)
end

affect_equations(cb::SymbolicDiscreteCallback) = cb.affects
function affect_equations(cbs::Vector{SymbolicDiscreteCallback})
    reduce(vcat, affect_equations(cb) for cb in cbs)
end
function namespace_equation(cb::SymbolicDiscreteCallback, s)::SymbolicDiscreteCallback
    SymbolicDiscreteCallback(namespace_expr(condition(cb), s),
                             namespace_equation.(affect_equations(cb), Ref(s)))
end

SymbolicDiscreteCallbacks(cb::Pair) = SymbolicDiscreteCallback[SymbolicDiscreteCallback(cb)]
SymbolicDiscreteCallbacks(cbs::Vector) = SymbolicDiscreteCallback.(cbs)
SymbolicDiscreteCallbacks(cb::SymbolicDiscreteCallback) = [cb]
SymbolicDiscreteCallbacks(cbs::Vector{<:SymbolicDiscreteCallback}) = cbs
SymbolicDiscreteCallbacks(::Nothing) = SymbolicDiscreteCallback[]

function discrete_events(sys::AbstractSystem)
    obs = get_discrete_events(sys)
    systems = get_systems(sys)
    cbs = [obs;
           reduce(vcat,
                  (map(o -> namespace_equation(o, s), discrete_events(s)) for s in systems),
                  init = SymbolicDiscreteCallback[])]
    cbs
end

################################# compilation functions ####################################

# handles ensuring that affect! functions work with integrator arguments
function add_integrator_header(out=:u)
    integrator = gensym(:MTKIntegrator)

    expr -> Func([DestructuredArgs(expr.args, integrator, inds = [:u, :p, :t])], [],
                 expr.body),
    expr -> Func([DestructuredArgs(expr.args, integrator, inds = [out, :u, :p, :t])], [],
                 expr.body)
end

function condition_header()
    integrator = gensym(:MTKIntegrator)
    expr -> Func([expr.args[1], expr.args[2],
                  DestructuredArgs(expr.args[3:end], integrator, inds = [:p])], [], expr.body)
end

"""
    compile_condition(cb::SymbolicDiscreteCallback, sys, dvs, ps; expression, kwargs...)

Returns a function `condition(u,p,t)` returning the `condition(cb)`.

Notes
- `expression = Val{true}`, causes the generated function to be returned as an expression.
  If  set to `Val{false}` a `RuntimeGeneratedFunction` will be returned.
- `kwargs` are passed through to `Symbolics.build_function`.
"""
function compile_condition(cb::SymbolicDiscreteCallback, sys, dvs, ps;
                           expression = Val{true}, kwargs...)

    u = map(x -> time_varying_as_func(value(x), sys), dvs)
    p = map(x -> time_varying_as_func(value(x), sys), ps)
    t = get_iv(sys)
    condit = condition(cb)
    build_function(condit, u, t, p; expression, wrap_code = condition_header(), kwargs...)
end

function compile_affect(cb::SymbolicContinuousCallback, args...; kwargs...)
    compile_affect(affect_equations(cb), args...; kwargs...)
end

"""
    compile_affect(eqs::Vector{Equation}, sys, dvs, ps; expression, outputidxs, kwargs...)
    compile_affect(cb::SymbolicContinuousCallback, args...; kwargs...)

Returns a function that takes an integrator as argument and modifies the state with the
affect. The generated function has the signature `affect!(integrator)`.

Notes
- `expression = Val{true}`, causes the generated function to be returned as an expression.
  If  set to `Val{false}` a `RuntimeGeneratedFunction` will be returned.
- `outputidxs`, a vector of indices of the states that correspond to outputs.
- `kwargs` are passed through to `Symbolics.build_function`.
"""
function compile_affect(eqs::Vector{Equation}, sys, dvs, ps; outputidxs = nothing,
                        expression = Val{true},
                        kwargs...)
    if isempty(eqs)
        if expression == Val{true}
            return :((args...) -> ())
        else
            return (args...) -> () # We don't do anything in the callback, we're just after the event
        end
    else
        rhss = map(x -> x.rhs, eqs)
        outvar = :u
        if outputidxs === nothing
            lhss = map(x -> x.lhs, eqs)
            update_vars = collect(Iterators.flatten(map(ModelingToolkit.vars, lhss))) # these are the ones we're chaning
            length(update_vars) == length(unique(update_vars)) == length(eqs) ||
                error("affected variables not unique, each state can only be affected by one equation for a single `root_eqs => affects` pair.")
            alleq = all(isequal(isparameter(first(update_vars))),
                        Iterators.map(isparameter, update_vars))
            if !isparameter(first(lhss)) && alleq
                stateind(sym) = findfirst(isequal(sym), dvs)
                update_inds = stateind.(update_vars)
            elseif isparameter(first(lhss)) && alleq
                psind(sym) = findfirst(isequal(sym), ps)
                update_inds = psind.(update_vars)
                outvar = :p
            else
                error("Error, building an affect function for a callback that wants to modify both parameters and states. This is not currently allowed in one individual callback.")
            end
        else
            update_inds = outputidxs
        end

        u = map(x -> time_varying_as_func(value(x), sys), dvs)
        p = map(x -> time_varying_as_func(value(x), sys), ps)
        t = get_iv(sys)
        rf_oop, rf_ip = build_function(rhss, u, p, t; expression = expression,
                                       wrap_code = add_integrator_header(outvar),
                                       outputidxs = update_inds,
                                       kwargs...)
        rf_ip
    end
end

function generate_rootfinding_callback(sys::AbstractODESystem, dvs = states(sys),
                                       ps = parameters(sys); kwargs...)
    cbs = continuous_events(sys)
    isempty(cbs) && return nothing
    generate_rootfinding_callback(cbs, sys, dvs, ps; kwargs...)
end

function generate_rootfinding_callback(cbs, sys::AbstractODESystem, dvs = states(sys),
                                       ps = parameters(sys); kwargs...)
    eqs = map(cb -> cb.eqs, cbs)
    num_eqs = length.(eqs)
    (isempty(eqs) || sum(num_eqs) == 0) && return nothing
    # fuse equations to create VectorContinuousCallback
    eqs = reduce(vcat, eqs)
    # rewrite all equations as 0 ~ interesting stuff
    eqs = map(eqs) do eq
        isequal(eq.lhs, 0) && return eq
        0 ~ eq.lhs - eq.rhs
    end

    rhss = map(x -> x.rhs, eqs)
    root_eq_vars = unique(collect(Iterators.flatten(map(ModelingToolkit.vars, rhss))))

    u = map(x -> time_varying_as_func(value(x), sys), dvs)
    p = map(x -> time_varying_as_func(value(x), sys), ps)
    t = get_iv(sys)
    rf_oop, rf_ip = build_function(rhss, u, p, t; expression = Val{false}, kwargs...)

    affect_functions = map(cbs) do cb # Keep affect function separate
        eq_aff = affect_equations(cb)
        affect = compile_affect(eq_aff, sys, dvs, ps; expression = Val{false}, kwargs...)
    end

    if length(eqs) == 1
        cond = function (u, t, integ)
            if DiffEqBase.isinplace(integ.sol.prob)
                tmp, = DiffEqBase.get_tmp_cache(integ)
                rf_ip(tmp, u, integ.p, t)
                tmp[1]
            else
                rf_oop(u, integ.p, t)
            end
        end
        ContinuousCallback(cond, affect_functions[])
    else
        cond = function (out, u, t, integ)
            rf_ip(out, u, integ.p, t)
        end

        # since there may be different number of conditions and affects,
        # we build a map that translates the condition eq. number to the affect number
        eq_ind2affect = reduce(vcat,
                               [fill(i, num_eqs[i]) for i in eachindex(affect_functions)])
        @assert length(eq_ind2affect) == length(eqs)
        @assert maximum(eq_ind2affect) == length(affect_functions)

        affect = let affect_functions = affect_functions, eq_ind2affect = eq_ind2affect
            function (integ, eq_ind) # eq_ind refers to the equation index that triggered the event, each event has num_eqs[i] equations
                affect_functions[eq_ind2affect[eq_ind]](integ)
            end
        end
        VectorContinuousCallback(cond, affect, length(eqs))
    end
end

function generate_discrete_callbacks(sys::AbstractSystem, dvs = states(sys),
                                    ps = parameters(sys); kwargs...)
    has_discrete_events(sys) || return nothing
    symcbs = discrete_events(sys)
    isempty(symcbs) && return nothing

    dbs = map(symcbs) do cb
        c = compile_condition(cb, sys, dvs, ps; expression=Val{false}, kwargs...)
        as = compile_affect(affect_equations(cb), sys, dvs, ps; expression = Val{false},
                            kwargs...)
        DiscreteCallback(c, as)
    end

    dbs
end

merge_cb(::Nothing, ::Nothing) = nothing
merge_cb(::Nothing, x) = merge_cb(x, nothing)
merge_cb(x, ::Nothing) = x
merge_cb(x, y) = CallbackSet(x, y)

function process_events(sys; callback = nothing, has_difference = false, kwargs...)
    if has_continuous_events(sys)
        contin_cb = generate_rootfinding_callback(sys; kwargs...)
    else
        contin_cb = nothing
    end
    if has_discrete_events(sys)
        discrete_cb = generate_discrete_callbacks(sys; kwargs...)
    else
        discrete_cb = nothing
    end
    difference_cb = has_difference ? generate_difference_cb(sys; kwargs...) : nothing

    cb = CallbackSet(contin_cb, difference_cb, callback, discrete_cb...)
    # cb = merge_cb(contin_cb, discrete_cb...)
    # cb = merge_cb(cb, difference_cb)
    # cb = merge_cb(cb, callback)
    # @show typeof(cbs),cbs
    cb
end