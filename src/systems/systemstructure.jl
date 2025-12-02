export mtkcompile!

"""
    $(TYPEDSIGNATURES)

Descend through the system hierarchy and look for statemachines. Remove equations from
the inner statemachine systems. Return the new `sys` and an array of top-level
statemachines.
"""
function extract_top_level_statemachines(sys::System)
    eqs = get_eqs(sys)
    predicate = Base.Fix2(isa, MTKTearing.StateMachineOperator) ∘ SU.unwrap_const
    if !isempty(eqs) && all(predicate, eqs)
        # top-level statemachine
        with_removed = @set sys.systems = map(remove_child_equations, get_systems(sys))
        return with_removed, [sys]
    elseif !isempty(eqs) && any(predicate, eqs)
        # error: can't mix
        error("Mixing statemachine equations and standard equations in a top-level statemachine is not allowed.")
    else
        # descend
        subsystems = get_systems(sys)
        newsubsystems = System[]
        statemachines = System[]
        for subsys in subsystems
            newsubsys, sub_statemachines = extract_top_level_statemachines(subsys)
            push!(newsubsystems, newsubsys)
            append!(statemachines, sub_statemachines)
        end
        @set! sys.systems = newsubsystems
        return sys, statemachines
    end
end

"""
    $(TYPEDSIGNATURES)

Return `sys` with all equations (including those in subsystems) removed.
"""
function remove_child_equations(sys::System)
    @set! sys.eqs = Equation[]
    @set! sys.systems = map(remove_child_equations, get_systems(sys))
    return sys
end

function make_eqs_zero_equals!(ts::TearingState)
    neweqs = map(enumerate(get_eqs(ts.sys))) do kvp
        i, eq = kvp
        isalgeq = true
        for j in 𝑠neighbors(ts.structure.graph, i)
            isalgeq &= invview(ts.structure.var_to_diff)[j] === nothing
        end
        if isalgeq
            return 0 ~ eq.rhs - eq.lhs
        else
            return eq
        end
    end
    copyto!(get_eqs(ts.sys), neweqs)
end


"""
Turn input variables into parameters of the system.
"""
function inputs_to_parameters!(state::TearingState, inputsyms::OrderedSet{SymbolicT}, outputsyms::OrderedSet{SymbolicT})
    @unpack structure, fullvars, sys = state
    @unpack var_to_diff, graph, solvable_graph = structure
    @assert solvable_graph === nothing

    var_reidx = zeros(Int, length(fullvars))
    nvar = 0
    new_fullvars = SymbolicT[]
    for (i, v) in enumerate(fullvars)
        if v in inputsyms
            if var_to_diff[i] !== nothing
                error("Input $(fullvars[i]) is differentiated!")
            end
            var_reidx[i] = -1
        else
            nvar += 1
            var_reidx[i] = nvar
            push!(new_fullvars, v)
        end
    end
    ninputs = length(inputsyms)
    @set! sys.inputs = inputsyms
    @set! sys.outputs = outputsyms
    if ninputs == 0
        state.sys = sys
        return state
    end

    nvars = ndsts(graph) - ninputs
    new_graph = BipartiteGraph(nsrcs(graph), nvars, Val(false))

    for ie in 1:nsrcs(graph)
        for iv in 𝑠neighbors(graph, ie)
            iv = var_reidx[iv]
            iv > 0 || continue
            add_edge!(new_graph, ie, iv)
        end
    end

    new_var_to_diff = StateSelection.DiffGraph(nvars, true)
    for (i, v) in enumerate(var_to_diff)
        new_i = var_reidx[i]
        (new_i < 1 || v === nothing) && continue
        new_v = var_reidx[v]
        @assert new_v > 0
        new_var_to_diff[new_i] = new_v
    end
    @set! structure.var_to_diff = complete(new_var_to_diff)
    @set! structure.graph = complete(new_graph)

    @set! sys.unknowns = setdiff(unknowns(sys), inputsyms)
    ps = copy(parameters(sys))
    append!(ps, inputsyms)
    @set! sys.ps = ps
    @set! state.sys = sys
    @set! state.fullvars = Vector{SymbolicT}(new_fullvars)
    @set! state.structure = structure
    return state
end

function mtkcompile!(state::TearingState;
        check_consistency = true, fully_determined = true,
        inputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        outputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        disturbance_inputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        kwargs...)
    if !is_time_dependent(state.sys)
        return _mtkcompile!(state; check_consistency,
            inputs, outputs, disturbance_inputs,
            fully_determined, kwargs...)
    end
    # split_system returns one or two systems and the inputs for each
    # mod clock inference to be binary
    # if it's continuous keep going, if not then error unless given trait impl in additional passes
    ci = MTKTearing.ClockInference(state)
    ci = MTKTearing.infer_clocks!(ci)
    tss, clocked_inputs, continuous_id, id_to_clock = MTKTearing.split_system(ci)
    if !isempty(tss) && continuous_id == 0
        # do a trait check here - handle fully discrete system
        additional_passes = get(kwargs, :additional_passes, nothing)
        if !isnothing(additional_passes) && any(discrete_compile_pass, additional_passes)
            # take the first discrete compilation pass given for now
            discrete_pass_idx = findfirst(discrete_compile_pass, additional_passes)
            discrete_compile = additional_passes[discrete_pass_idx]
            deleteat!(additional_passes, discrete_pass_idx)
            return discrete_compile(tss, clocked_inputs, ci)
        end
        throw(HybridSystemNotSupportedException("""
        Discrete systems with multiple clocks are not supported with the standard \
        MTK compiler.
        """))
    end
    if length(tss) > 1
        make_eqs_zero_equals!(tss[continuous_id])
        # simplify as normal
        sys = _mtkcompile!(tss[continuous_id]; simplify,
            inputs, outputs, disturbance_inputs,
            discrete_inputs = OrderedSet{SymbolicT}(clocked_inputs[continuous_id]),
            check_consistency, fully_determined,
            kwargs...)
        additional_passes = get(kwargs, :additional_passes, nothing)
        if !isnothing(additional_passes) && any(discrete_compile_pass, additional_passes)
            discrete_pass_idx = findfirst(discrete_compile_pass, additional_passes)
            discrete_compile = additional_passes[discrete_pass_idx]
            deleteat!(additional_passes, discrete_pass_idx)
            # in the case of a hybrid system, the discrete_compile pass should take the currents of sys.discrete_subsystems
            # and modifies discrete_subsystems to bea tuple of the io and anything else, while adding or manipulating the rest of sys as needed
            return discrete_compile(
                sys, tss[[i for i in eachindex(tss) if i != continuous_id]],
                clocked_inputs, ci, id_to_clock)
        end
        throw(HybridSystemNotSupportedException("""
        Hybrid continuous-discrete systems are currently not supported with \
        the standard MTK compiler. This system requires JuliaSimCompiler.jl, \
        see https://help.juliahub.com/juliasimcompiler/stable/
        """))
    end
    if get_is_discrete(state.sys) ||
       continuous_id == 1 && any(Base.Fix2(isoperator, Shift), state.fullvars)
        state.structure.only_discrete = true
        state = MTKTearing.shift_discrete_system(state)
        sys = state.sys
        @set! sys.is_discrete = true
        state.sys = sys
    end

    sys = _mtkcompile!(state; check_consistency,
        inputs, outputs, disturbance_inputs,
        fully_determined, kwargs...)
    return sys
end

function _mtkcompile!(state::TearingState;
        check_consistency = true, fully_determined = true,
        dummy_derivative = true,
        discrete_inputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        inputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        outputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        disturbance_inputs::OrderedSet{SymbolicT} = OrderedSet{SymbolicT}(),
        kwargs...)
    if fully_determined isa Bool
        check_consistency &= fully_determined
    else
        check_consistency = true
    end
    orig_inputs = Set{SymbolicT}()
    validate_io!(state, orig_inputs, inputs, discrete_inputs, outputs, disturbance_inputs)
    # ModelingToolkit.markio!(state, orig_inputs, inputs, outputs, disturbance_inputs)
    union!(inputs, disturbance_inputs)
    state = ModelingToolkit.inputs_to_parameters!(state, discrete_inputs, OrderedSet{SymbolicT}())
    state = ModelingToolkit.inputs_to_parameters!(state, inputs, outputs)
    StateSelection.trivial_tearing!(state)
    sys, mm = ModelingToolkit.alias_elimination!(state; fully_determined, kwargs...)
    if check_consistency
        fully_determined = StateSelection.check_consistency(
            state, orig_inputs; nothrow = fully_determined === nothing)
    end
    # This phrasing avoids making the `kwcall` dynamic dispatch due to the type of a
    # keyword (`mm`) being non-concrete
    if mm isa CLIL.SparseMatrixCLIL{BigInt, Int}
        sys = _mtkcompile_worker!(state, sys, mm; fully_determined, dummy_derivative, kwargs...)
    else
        sys =_mtkcompile_worker!(state, sys, mm; fully_determined, dummy_derivative, kwargs...)
    end
    fullunknowns = [observables(sys); unknowns(sys)]
    @set! sys.observed = MTKBase.topsort_equations(observed(sys), fullunknowns)

    MTKBase.invalidate_cache!(sys)
end

function _mtkcompile_worker!(state::TearingState, sys::System, mm::CLIL.SparseMatrixCLIL{T, Int};
                             fully_determined::Bool, dummy_derivative::Bool,
                             kwargs...) where {T}
    if fully_determined && dummy_derivative
        sys = ModelingToolkit.dummy_derivative(
            sys, state; mm, kwargs...)
    elseif fully_determined
        var_eq_matching = StateSelection.pantelides!(state; finalize = false, kwargs...)
        sys = pantelides_reassemble(state, var_eq_matching)
        state = TearingState(sys)
        sys, mm::CLIL.SparseMatrixCLIL{T, Int} = ModelingToolkit.alias_elimination!(state; fully_determined, kwargs...)
        sys = ModelingToolkit.dummy_derivative(
            sys, state; mm, fully_determined, kwargs...)
    else
        sys = ModelingToolkit.tearing(
            sys, state; mm, fully_determined, kwargs...)
    end
    return sys
end

function validate_io!(state::TearingState, orig_inputs::Set{SymbolicT}, inputs::OrderedSet{SymbolicT},
                      discrete_inputs::OrderedSet{SymbolicT}, outputs::OrderedSet{SymbolicT},
                      disturbance_inputs::OrderedSet{SymbolicT})
    for v in state.fullvars
        isinput(v) && push!(orig_inputs, v)
    end
    fullvars_set = OrderedSet{SymbolicT}(state.fullvars)
    missings = OrderedSet{SymbolicT}()
    union!(missings, inputs)
    setdiff!(missings, fullvars_set)
    isempty(missings) || throw(IONotFoundError("inputs", nameof(state.sys), missings))

    union!(missings, discrete_inputs)
    setdiff!(missings, fullvars_set)
    isempty(missings) || throw(IONotFoundError("discrete inputs", nameof(state.sys), missings))

    union!(missings, outputs)
    setdiff!(missings, fullvars_set)
    isempty(missings) || throw(IONotFoundError("outputs", nameof(state.sys), missings))

    union!(missings, disturbance_inputs)
    setdiff!(missings, fullvars_set)
    isempty(missings) || throw(IONotFoundError("disturbance inputs", nameof(state.sys), missings))

    return nothing
end

struct DifferentiatedVariableNotUnknownError <: Exception
    differentiated::Any
    undifferentiated::Any
end

function Base.showerror(io::IO, err::DifferentiatedVariableNotUnknownError)
    undiff = err.undifferentiated
    diff = err.differentiated
    print(io,
        "Variable $undiff occurs differentiated as $diff but is not an unknown of the system.")
    scope = getmetadata(undiff, SymScope, LocalScope())
    depth = expected_scope_depth(scope)
    if depth > 0
        print(io,
            "\nVariable $undiff expects $depth more levels in the hierarchy to be an unknown.")
    end
end
