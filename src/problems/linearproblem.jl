struct LinearFunction{iip, I} <: SciMLBase.AbstractSciMLFunction{iip}
    interface::I
    A::AbstractMatrix
    b::AbstractVector
end

function LinearFunction{iip}(
        sys::System; expression = Val{false}, check_compatibility = true,
        sparse = false, eval_expression = false, eval_module = @__MODULE__,
        checkbounds = false, cse = true, kwargs...) where {iip}
    check_complete(sys, LinearProblem)
    check_compatibility && check_compatible_system(LinearProblem, sys)

    A, b = calculate_A_b(sys; sparse)
    update_A = generate_update_A(sys, A; expression, wrap_gfw = Val{true}, eval_expression,
        eval_module, checkbounds, cse, kwargs...)
    update_b = generate_update_b(sys, b; expression, wrap_gfw = Val{true}, eval_expression,
        eval_module, checkbounds, cse, kwargs...)
    observedfun = ObservedFunctionCache(
        sys; steady_state = false, expression, eval_expression, eval_module, checkbounds,
        cse)

    if expression == Val{true}
        symbolic_interface = quote
            update_A = $update_A
            update_b = $update_b
            sys = $sys
            observedfun = $observedfun
            $(SciMLBase.SymbolicLinearInterface)(
                update_A, update_b, sys, observedfun, nothing)
        end
    else
        symbolic_interface = SciMLBase.SymbolicLinearInterface(
            update_A, update_b, sys, observedfun, nothing)
    end

    return LinearFunction{iip, typeof(symbolic_interface)}(symbolic_interface, A, b)
end

function SciMLBase.LinearProblem(sys::System, op; kwargs...)
    SciMLBase.LinearProblem{true}(sys, op; kwargs...)
end

function SciMLBase.LinearProblem(sys::System, op::StaticArray; kwargs...)
    SciMLBase.LinearProblem{false}(sys, op; kwargs...)
end

function SciMLBase.LinearProblem{iip}(
        sys::System, op; check_length = true, expression = Val{false},
        check_compatibility = true, sparse = false, eval_expression = false,
        eval_module = @__MODULE__, u0_constructor = identity, u0_eltype = nothing,
        kwargs...) where {iip}
    check_complete(sys, LinearProblem)
    check_compatibility && check_compatible_system(LinearProblem, sys)

    f, u0,
    p = process_SciMLProblem(
        LinearFunction{iip}, sys, op; check_length, expression,
        build_initializeprob = false, symbolic_u0 = true, u0_constructor, u0_eltype,
        kwargs...)

    if any(x -> symbolic_type(x) != NotSymbolic(), u0)
        u0 = nothing
    end

    u0Type = typeof(op)
    floatT = if u0 === nothing
        calculate_float_type(op, u0Type)
    else
        eltype(u0)
    end
    u0_eltype = something(u0_eltype, floatT)

    u0_constructor = get_p_constructor(u0_constructor, u0Type, u0_eltype)
    symbolic_interface = f.interface
    A,
    b = get_A_b_from_LinearFunction(
        sys, f, p; eval_expression, eval_module, expression, u0_constructor, sparse)

    kwargs = (; u0, process_kwargs(sys; kwargs...)..., f = symbolic_interface)
    args = (; A, b, p)

    return maybe_codegen_scimlproblem(expression, LinearProblem{iip}, args; kwargs...)
end

function get_A_b_from_LinearFunction(
        sys::System, f::LinearFunction, p; eval_expression = false,
        eval_module = @__MODULE__, expression = Val{false}, u0_constructor = identity,
        u0_eltype = float, sparse = false)
    @unpack A, b, interface = f
    if expression == Val{true}
        get_A = build_explicit_observed_function(
            sys, A; param_only = true, eval_expression, eval_module)
        get_b = build_explicit_observed_function(
            sys, b; param_only = true, eval_expression, eval_module)
        A = u0_constructor(u0_eltype.(get_A(p)))
        b = u0_constructor(u0_eltype.(get_b(p)))
    else
        A = u0_constructor(u0_eltype.(interface.update_A!(p)))
        b = u0_constructor(u0_eltype.(interface.update_b!(p)))
    end
    if sparse
        A = SparseArrays.sparse(A)
    end

    return A, b
end

# For remake
function SciMLBase.get_new_A_b(
        sys::AbstractSystem, f::SciMLBase.SymbolicLinearInterface, p, A, b; kw...)
    if ArrayInterface.ismutable(A)
        f.update_A!(A, p)
        f.update_b!(b, p)
    else
        # The generated function has both IIP and OOP variants
        A = StaticArraysCore.similar_type(A)(f.update_A!(p))
        b = StaticArraysCore.similar_type(b)(f.update_b!(p))
    end
    return A, b
end

function check_compatible_system(T::Type{LinearProblem}, sys::System)
    check_time_independent(sys, T)
    check_affine(sys, T)
    check_not_dde(sys)
    check_no_cost(sys, T)
    check_no_constraints(sys, T)
    check_no_jumps(sys, T)
    check_no_noise(sys, T)
end
