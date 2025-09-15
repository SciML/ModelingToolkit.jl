_nameof(s) = nameof(s)
_nameof(s::Union{Int, Symbol}) = s
abstract type StateMachineOperator end
Base.broadcastable(x::StateMachineOperator) = Ref(x)
Symbolics.hide_lhs(_::StateMachineOperator) = true
struct InitialState <: StateMachineOperator
    s::Any
end
Base.show(io::IO, s::InitialState) = print(io, "initial_state(", _nameof(s.s), ")")
initial_state(s) = Equation(InitialState(nothing), InitialState(s))

Base.@kwdef struct Transition{A, B, C} <: StateMachineOperator
    from::A = nothing
    to::B = nothing
    cond::C = nothing
    immediate::Bool = true
    reset::Bool = true
    synchronize::Bool = false
    priority::Int = 1
    function Transition(from, to, cond, immediate, reset, synchronize, priority)
        cond = unwrap(cond)
        new{typeof(from), typeof(to), typeof(cond)}(from, to, cond, immediate,
            reset, synchronize,
            priority)
    end
end
function Base.:(==)(transition1::Transition, transition2::Transition)
    transition1.from == transition2.from &&
        transition1.to == transition2.to &&
        isequal(transition1.cond, transition2.cond) &&
        transition1.immediate == transition2.immediate &&
        transition1.reset == transition2.reset &&
        transition1.synchronize == transition2.synchronize &&
        transition1.priority == transition2.priority
end

"""
    transition(from, to, cond; immediate::Bool = true, reset::Bool = true, synchronize::Bool = false, priority::Int = 1)

Create a transition from state `from` to state `to` that is enabled when transitioncondition `cond` evaluates to `true`.

# Arguments:
- `from`: The source state of the transition.
- `to`: The target state of the transition.
- `cond`: A transition condition that evaluates to a Bool, such as `ticksInState() >= 2`.
- `immediate`: If `true`, the transition will fire at the same tick as it becomes true, if `false`, the actions of the state are evaluated first, and the transition fires during the next tick.
- `reset`: If true, the destination state `to` is reset to its initial condition when the transition fires.
- `synchronize`: If true, the transition will only fire if all sub-state machines in the source state are in their final (terminal) state. A final state is one that has no outgoing transitions.
- `priority`: If a state has more than one outgoing transition, all outgoing transitions must have a unique priority. The transitions are evaluated in priority order, i.e., the transition with priority 1 is evaluated first.
"""
function transition(from, to, cond;
        immediate::Bool = true, reset::Bool = true, synchronize::Bool = false,
        priority::Int = 1)
    Equation(
        Transition(), Transition(; from, to, cond, immediate, reset,
            synchronize, priority))
end
function Base.show(io::IO, s::Transition)
    print(io, _nameof(s.from), " → ", _nameof(s.to), " if (", s.cond, ") [")
    print(io, "immediate: ", Int(s.immediate), ", ")
    print(io, "reset: ", Int(s.reset), ", ")
    print(io, "sync: ", Int(s.synchronize), ", ")
    print(io, "prio: ", s.priority, "]")
end

function activeState end
function entry end
function ticksInState end
function timeInState end

for (s, T) in [(:timeInState, :Real),
    (:ticksInState, :Integer),
    (:entry, :Bool),
    (:activeState, :Bool)]
    seed = hash(s)
    @eval begin
        $s(x) = wrap(term($s, x))
        SymbolicUtils.promote_symtype(::typeof($s), _...) = $T
        function SymbolicUtils.show_call(io, ::typeof($s), args)
            if isempty(args)
                print(io, $s, "()")
            else
                arg = only(args)
                print(io, $s, "(", arg isa Number ? arg : nameof(arg), ")")
            end
        end
    end
    if s != :activeState
        @eval $s() = wrap(term($s))
    end
end

@doc """
    timeInState()
    timeInState(state)

Get the time (in seconds) spent in a state in a finite state machine.

When used to query the time spent in the enclosing state, the method without arguments is used, i.e.,
```julia
@mtkmodel FSM begin
    ...
    @equations begin
        var(k+1) ~ timeInState() >= 2 ? 0.0 : var(k)
    end
end
```

If used to query the residence time of another state, the state is passed as an argument.

This operator can be used in both equations and transition conditions.

See also [`ticksInState`](@ref) and [`entry`](@ref)
""" timeInState

@doc """
    ticksInState()
    ticksInState(state)

Get the number of ticks spent in a state in a finite state machine.

When used to query the number of ticks spent in the enclosing state, the method without arguments is used, i.e.,
```julia
@mtkmodel FSM begin
    ...
    @equations begin
        var(k+1) ~ ticksInState() >= 2 ? 0.0 : var(k)
    end
end
```

If used to query the number of ticks in another state, the state is passed as an argument.

This operator can be used in both equations and transition conditions.

See also [`timeInState`](@ref) and [`entry`](@ref)
""" ticksInState

@doc """
    entry()
    entry(state)

When used in a finite-state machine, this operator returns true at the first tick when the state is active, and false otherwise.

When used to query the entry of the enclosing state, the method without arguments is used, when used to query the entry of another state, the state is passed as an argument.

This can be used to perform a unique action when entering a state.
"""
entry

@doc """
    activeState(state)

When used in a finite state machine, this operator returns `true` if the queried state is active and false otherwise. 
""" activeState

function vars!(vars, O::Transition; op = Differential)
    vars!(vars, O.from)
    vars!(vars, O.to)
    vars!(vars, O.cond; op)
    return vars
end
function vars!(vars, O::InitialState; op = Differential)
    vars!(vars, O.s; op)
    return vars
end
function vars!(vars, O::StateMachineOperator; op = Differential)
    error("Unhandled state machine operator")
end

function namespace_expr(
        O::Transition, sys, n = nameof(sys); ivs = independent_variables(sys))
    return Transition(
        O.from === nothing ? O.from : renamespace(sys, O.from),
        O.to === nothing ? O.to : renamespace(sys, O.to),
        O.cond === nothing ? O.cond : namespace_expr(O.cond, sys),
        O.immediate, O.reset, O.synchronize, O.priority
    )
end

function namespace_expr(
        O::InitialState, sys, n = nameof(sys); ivs = independent_variables(sys))
    return InitialState(O.s === nothing ? O.s : renamespace(sys, O.s))
end

function namespace_expr(O::StateMachineOperator, sys, n = nameof(sys); kwargs...)
    error("Unhandled state machine operator")
end
