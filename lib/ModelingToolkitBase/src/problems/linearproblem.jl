struct LinearFunction{iip, I} <: SciMLBase.AbstractSciMLFunction{iip}
    interface::I
    A::Union{Matrix{SymbolicT}, SparseMatrixCSC{SymbolicT, Int}}
    b::Vector{SymbolicT}
end

function LinearFunction{iip}(
        sys::System; expression = Val{false}, check_compatibility = true,
        sparse = false, eval_expression = false, eval_module = @__MODULE__,
        checkbounds = false, cse = true, kwargs...
    ) where {iip}
    check_complete(sys, LinearProblem)
    check_compatibility && check_compatible_system(LinearProblem, sys)

    A, b = calculate_A_b(sys; sparse)
    update_A = generate_update_A(
        sys, A; expression, wrap_gfw = Val{true}, eval_expression,
        eval_module, checkbounds, cse, kwargs...
    )
    update_b = generate_update_b(
        sys, b; expression, wrap_gfw = Val{true}, eval_expression,
        eval_module, checkbounds, cse, kwargs...
    )
    observedfun = ObservedFunctionCache(
        sys; steady_state = false, expression, eval_expression, eval_module, checkbounds,
        cse
    )

    if expression == Val{true}
        symbolic_interface = quote
            update_A = $update_A
            update_b = $update_b
            sys = $sys
            observedfun = $observedfun
            $(SciMLBase.SymbolicLinearInterface)(
                update_A, update_b, sys, observedfun, nothing
            )
        end
    else
        symbolic_interface = SciMLBase.SymbolicLinearInterface(
            update_A, update_b, sys, observedfun, nothing
        )
    end

    return LinearFunction{iip, typeof(symbolic_interface)}(symbolic_interface, A, b)
end

function SciMLBase.LinearProblem(sys::System, op; kwargs...)
    return SciMLBase.LinearProblem{true}(sys, op; kwargs...)
end

function SciMLBase.LinearProblem(sys::System, op::StaticArray; kwargs...)
    return SciMLBase.LinearProblem{false}(sys, op; kwargs...)
end

function SciMLBase.LinearProblem{iip}(
        sys::System, op; check_length = true, expression = Val{false},
        check_compatibility = true, sparse = false, eval_expression = false,
        eval_module = @__MODULE__, u0_constructor = identity, u0_eltype = nothing,
        kwargs...
    ) where {iip}
    check_complete(sys, LinearProblem)
    check_compatibility && check_compatible_system(LinearProblem, sys)

    u0Type = typeof(op)
    f, u0, p, op = process_SciMLProblem(
        LinearFunction{iip}, sys, op; check_length, expression,
        build_initializeprob = false, symbolic_u0 = true, u0_constructor, u0_eltype,
        return_operating_point = true, sparse, kwargs...
    )

    if any(x -> symbolic_type(x) != NotSymbolic() || x === nothing, u0)
        u0 = nothing
    end

    floatT = if u0 === nothing
        calculate_float_type(op, u0Type)
    else
        eltype(u0)
    end
    u0_eltype = something(u0_eltype, floatT)

    u0_constructor = get_p_constructor(u0_constructor, u0Type, u0_eltype)
    symbolic_interface = f.interface
    A, b = get_A_b_from_LinearFunction(
        sys, f, op; eval_expression, eval_module, expression, u0_constructor
    )

    kwargs = (; u0, process_kwargs(sys; kwargs...)..., f = symbolic_interface)
    args = (; A, b, p)

    return maybe_codegen_scimlproblem(expression, LinearProblem{iip}, args; kwargs...)
end

function get_A_b_from_LinearFunction(
        sys::System, f::LinearFunction, op; kws...
    )
    return get_A_b_from_LinearFunction(sys, f, Symbolics.FixpointSubstituter{true}(op); kws...)
end

function get_A_b_from_LinearFunction(
        sys::System, f::LinearFunction, subber::Symbolics.FixpointSubstituter{true}; eval_expression = false,
        eval_module = @__MODULE__, expression = Val{false}, u0_constructor = identity,
        u0_eltype = float
    )
    @unpack A, b, interface = f
    if A isa Matrix{SymbolicT}
        _A = similar(A, Any)
        for i in eachindex(A)
            _A[i] = unwrap_const(subber(A[i]))
        end
        A = u0_eltype.(_A)
        _A = u0_constructor(A)
        if ArrayInterface.ismutable(_A)
            A = similar(_A, size(A))
            copyto!(A, _A)
        else
            A = StaticArraysCore.similar_type(_A, StaticArraysCore.Size(size(A)))(_A)
        end
    else
        I, J, V = findnz(A)
        _V = similar(V, Any)
        for i in eachindex(V)
            _V[i] = unwrap_const(subber(V[i]))
        end
        V = u0_constructor(u0_eltype.(_V))
        A = SparseArrays.sparse(I, J, V, size(A)...)
    end
    _b = similar(b, Any)
    for i in eachindex(b)
        _b[i] = unwrap_const(subber(b[i]))
    end
    b = u0_constructor(u0_eltype.(_b))

    return A, b
end

# For remake
function SciMLBase.get_new_A_b(
        sys::AbstractSystem, f::SciMLBase.SymbolicLinearInterface, p, A, b; kw...
    )
    if ArrayInterface.ismutable(A)
        T = eltype(SciMLStructures.canonicalize(SciMLStructures.Tunable(), p)[1])
        if eltype(A) !== T
            A = similar(A, T)
            b = similar(b, T)
        end
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
    return check_no_noise(sys, T)
end
