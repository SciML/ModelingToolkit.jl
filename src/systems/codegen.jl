function MTKBase.torn_system_jacobian_sparsity(sys::System)
    state = get_tearing_state(sys)
    state isa TearingState || return nothing
    @unpack structure = state
    @unpack graph, var_to_diff = structure

    neqs = nsrcs(graph)
    nsts = ndsts(graph)
    states_idxs = findall(!Base.Fix1(StateSelection.isdervar, structure), 1:nsts)
    var2idx = StructuralTransformations.uneven_invmap(nsts, states_idxs)
    I = Int[]
    J = Int[]
    for ieq in 1:neqs
        for ivar in 𝑠neighbors(graph, ieq)
            nivar = get(var2idx, ivar, 0)
            nivar == 0 && continue
            push!(I, ieq)
            push!(J, nivar)
        end
    end
    return sparse(I, J, true, neqs, neqs)
end

"""
    $(TYPEDSIGNATURES)

Represent the equations of the system `sys` in semiquadratic form. Returns 3 matrices
referred to as `A`, `B` and `C`. Hereon, `x` refers to the state vector.

`A` contains coefficients for all linear terms in the equations. `A * x` is the linear
part of the RHS of the equations. `A` is `nothing` if none of the equations have
linear terms, in which case the corresponding term in the mathematical expression below
should be ignored.

`B` is a vector of matrices where the `i`th matrix contains coefficients for the quadratic
terms in the `i`th equation's RHS. The quadratic terms in the `i`th equation are obtained
as `transpose(x) * B[i] * x`. Each matrix `B[i]` will be `nothing` if the `i`th equation
does not have quadratic term. in which case the corresponding term in the mathematical
expression below should be ignored. In case all `B[i]` are `nothing`, `B` will be
`nothing`.

`C` is a vector of all non-linear and non-quadratic terms in the equations. `C` is
`nothing` if none of the equations have nonlinear and non-quadratic terms, in which case
the corresponding term in the mathematical expression below should be ignored.

Mathematically, the right hand side of the `i`th equation is

```math
\\mathtt{row}_i(\\mathbf{A})x + x^T(\\mathbf{B}_i)x + \\mathbf{C}_i
```

Note that any of `A`, `B` or `C` can be `nothing` if the coefficients/values are all zeros.

## Keyword arguments

- `sparse`: Return sparse matrices for `A`, `B` and `C`.
"""
function calculate_semiquadratic_form(sys::System; sparse = false)
    rhss = [eq.rhs for eq in full_equations(sys)]
    dvs = unknowns(sys)
    A, B, x2, C = semiquadratic_form(rhss, dvs)
    if nnz(B) == 0
        B = nothing
        B2 = nothing
    else
        B2 = if sparse
            Any[spzeros(Num, length(dvs), length(dvs)) for _ in 1:length(rhss)]
        else
            Any[zeros(Num, length(dvs), length(dvs)) for _ in 1:length(rhss)]
        end
        idxmap = vec([CartesianIndex(i, j) for j in 1:length(dvs) for i in 1:j])
        for (i, j, val) in zip(findnz(B)...)
            B2[i][idxmap[j]] = val
        end
        for i in eachindex(B2)
            if all(_iszero, B2[i])
                B2[i] = nothing
            end
        end
        B2 = map(Broadcast.BroadcastFunction(unwrap), B2)
    end
    if nnz(A) == 0
        A = nothing
    else
        if !sparse
            A = collect(A)
        end
        A = unwrap.(A)
    end
    if all(_iszero, C)
        C = nothing
    else
        C = unwrap.(C)
    end

    return A, B2, C
end

const DIFFCACHE_PARAM_NAME = :__mtk_diffcache

"""
    $(TYPEDSIGNATURES)

Return a symbolic variable representing a `PreallocationTools.DiffCache` with
floating-point type `T`.
"""
function get_diffcache_param(::Type{T}) where {T}
    toconstant(Symbolics.variable(
        DIFFCACHE_PARAM_NAME; T = DiffCache{Vector{T}, Vector{T}}))
end

const LINEAR_MATRIX_PARAM_NAME = :linear_Aₘₜₖ

"""
    $(TYPEDSIGNATURES)

Return a symbolic variable representing the `A` matrix returned from
[`calculate_semiquadratic_form`](@ref).
"""
function get_linear_matrix_param(size::NTuple{2, Int})
    m, n = size
    unwrap(only(@constants $LINEAR_MATRIX_PARAM_NAME[1:m, 1:n]))
end

"""
    $(TYPEDSIGNATURES)

Return the name of the `i`th matrix in `B` returned from
[`calculate_semiquadratic_form`](@ref).
"""
function get_quadratic_form_name(i::Int)
    return Symbol(:quadratic_Bₘₜₖ_, i)
end

"""
    $(TYPEDSIGNATURES)

Return a symbolic variable representing the `i`th matrix in `B` returned from
[`calculate_semiquadratic_form`](@ref).
"""
function get_quadratic_form_param(sz::NTuple{2, Int}, i::Int)
    m, n = sz
    name = get_quadratic_form_name(i)
    unwrap(only(@constants $name[1:m, 1:n]))
end

"""
    $(TYPEDSIGNATURES)

Return the parameter in `sys` corresponding to the one returned from
[`get_linear_matrix_param`](@ref), or `nothing` if `A === nothing`.
"""
function get_linear_matrix_param_from_sys(sys::System, A)
    A === nothing && return nothing
    return unwrap(getproperty(sys, LINEAR_MATRIX_PARAM_NAME))
end

"""
    $(TYPEDSIGNATURES)

Return the list of parameters in `sys` corresponding to the ones returned from
[`get_quadratic_form_param`](@ref) for each non-`nothing` matrix in `B`, or `nothing` if
`B === nothing`. If `B[i] === nothing`, the returned list will have `nothing` as the `i`th
entry.
"""
function get_quadratic_form_params_from_sys(sys::System, B)
    B === nothing && return nothing
    return map(eachindex(B)) do i
        B[i] === nothing && return nothing
        return unwrap(getproperty(sys, get_quadratic_form_name(i)))
    end
end

"""
    $(TYPEDSIGNATURES)

Return the parameter in `sys` corresponding to the one returned from
[`get_diffcache_param`](@ref).
"""
function get_diffcache_param_from_sys(sys::System, B)
    B === nothing && return nothing
    return unwrap(getproperty(sys, DIFFCACHE_PARAM_NAME))
end

"""
    $(TYPEDSIGNATURES)

Generate `f1` and `f2` for [`SemilinearODEFunction`](@ref) (internally represented as a
`SplitFunction`). `A`, `B`, `C` are the matrices returned from
[`calculate_semiquadratic_form`](@ref). This expects that the system has the necessary
extra parmameters added by [`add_semiquadratic_parameters`](@ref).

## Keyword Arguments

$SEMILINEAR_A_B_C_KWARGS
$(MTKBase.GENERATE_X_KWARGS)

All other keyword arguments are forwarded to [`build_function_wrapper`](@ref).
$SEMILINEAR_A_B_C_CONSTRAINT

$(MTKBase.EXPERIMENTAL_WARNING)
"""
function generate_semiquadratic_functions(sys::System, A, B, C; stiff_linear = true,
        stiff_quadratic = false, stiff_nonlinear = false, expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...)
    if A === nothing && B === nothing
        throw(ArgumentError("Cannot generate split form for the system - it has no linear or quadratic part."))
    end

    if (stiff_linear || A === nothing) && (stiff_quadratic || B === nothing) &&
       stiff_nonlinear
        throw(ArgumentError("All of `A`, `B` and `C` cannot be stiff at the same time."))
    end
    if (!stiff_linear || A === nothing) && (!stiff_quadratic || B === nothing) &&
       !stiff_nonlinear
        throw(ArgumentError("All of `A`, `B` and `C` cannot be non-stiff at the same time."))
    end
    linear_matrix_param = get_linear_matrix_param_from_sys(sys, A)
    quadratic_forms = get_quadratic_form_params_from_sys(sys, B)
    diffcache_par = get_diffcache_param_from_sys(sys, B)
    eqs = equations(sys)
    dvs = unknowns(sys)
    ps = reorder_parameters(sys)
    iv = get_iv(sys)
    # Codegen is a bit manual, and we're manually creating an efficient IIP function.
    # Since we explicitly provide Symbolics.DEFAULT_OUTSYM, the `u` is actually the second
    # argument.
    iip_x = generated_argument_name(2)
    oop_x = generated_argument_name(1)

    shape = SU.Unknown(-1)
    ## iip
    f1_iip_ir = Assignment[]
    f2_iip_ir = Assignment[]
    # C
    if C !== nothing
        C_ir = stiff_nonlinear ? f1_iip_ir : f2_iip_ir
        push!(C_ir, Assignment(:__tmp_C, SetArray(false, Symbolics.DEFAULT_OUTSYM, C)))
    end
    # B
    if B !== nothing
        B_ir = stiff_quadratic ? f1_iip_ir : f2_iip_ir
        B_vals = map(eachindex(eqs)) do i
            B[i] === nothing && return nothing
            tmp_buf = term(
                PreallocationTools.get_tmp, diffcache_par, Symbolics.DEFAULT_OUTSYM;
                type = Vector{Real}, shape)
            tmp_buf = term(view, tmp_buf, 1:length(dvs); type = Vector{Real}, shape)

            result = term(*, term(transpose, iip_x; type = Matrix{Real}, shape),
                          :__tmp_B_1; type = Vector{Real}, shape)
            # if both write to the same buffer, don't overwrite
            if stiff_quadratic == stiff_nonlinear && C !== nothing
                result = term(+, result, term(getindex, Symbolics.DEFAULT_OUTSYM, i;
                                              type = Real, shape);
                              type = Real, shape)
            end
            intermediates = [
                Assignment(:__tmp_B_buffer, tmp_buf),
                Assignment(:__tmp_B_1,
                    term(mul!, :__tmp_B_buffer,
                        term(UpperTriangular, quadratic_forms[i]; type = Matrix{Real}, shape),
                        iip_x; type = Any, shape))
            ]
            return AtIndex(i, Let(intermediates, result))
        end
        filter!(x -> x !== nothing, B_vals)
        push!(B_ir, Assignment(:__tmp_B, SetArray(false, Symbolics.DEFAULT_OUTSYM, B_vals)))
    end
    # A
    if A !== nothing
        A_ir = stiff_linear ? f1_iip_ir : f2_iip_ir
        retain_old = stiff_linear == stiff_quadratic && B !== nothing ||
                     stiff_linear == stiff_nonlinear && C !== nothing
        push!(A_ir,
            Assignment(:__tmp_A,
                term(mul!, Symbolics.DEFAULT_OUTSYM,
                    linear_matrix_param, iip_x, true, retain_old; type = Any, shape)))
    end
    ## oop
    f1_terms = []
    f2_terms = []
    if A !== nothing
        push!(stiff_linear ? f1_terms : f2_terms,
              term(*, linear_matrix_param, oop_x; type = Vector{Real}, shape))
    end
    if B !== nothing
        B_elems = map(eachindex(eqs)) do i
            B[i] === nothing && return 0
            term(
                *, term(transpose, oop_x; type = Matrix{Real}, shape),
                term(UpperTriangular, quadratic_forms[i]; type = Matrix{Real}, shape),
                oop_x; type = Matrix{Real}, shape)
        end
        push!(stiff_quadratic ? f1_terms : f2_terms, MakeArray(B_elems, oop_x))
    end
    if C !== nothing
        push!(stiff_nonlinear ? f1_terms : f2_terms, MakeArray(C, oop_x))
    end
    f1_expr = if length(f1_terms) == 1
        only(f1_terms)
    else
        term(+, f1_terms...; type = Vector{Real}, shape)
    end
    f2_expr = if length(f2_terms) == 1
        only(f2_terms)
    else
        term(+, f2_terms...; type = Vector{Real}, shape)
    end

    f1_iip = build_function_wrapper(
        sys, nothing, Symbolics.DEFAULT_OUTSYM, dvs, ps..., iv; p_start = 3,
        extra_assignments = f1_iip_ir, expression = Val{true}, kwargs...)
    f2_iip = build_function_wrapper(
        sys, nothing, Symbolics.DEFAULT_OUTSYM, dvs, ps..., iv; p_start = 3,
        extra_assignments = f2_iip_ir, expression = Val{true}, kwargs...)
    f1_oop = build_function_wrapper(
        sys, f1_expr, dvs, ps..., iv; expression = Val{true}, kwargs...)
    if f1_oop isa NTuple{2, Expr}
        f1_oop = f1_oop[1]
    end
    f2_oop = build_function_wrapper(
        sys, f2_expr, dvs, ps..., iv; expression = Val{true}, kwargs...)
    if f2_oop isa NTuple{2, Expr}
        f2_oop = f2_oop[1]
    end

    f1 = maybe_compile_function(expression, wrap_gfw, (2, 3, is_split(sys)),
        (f1_oop, f1_iip); eval_expression, eval_module)
    f2 = maybe_compile_function(expression, wrap_gfw, (2, 3, is_split(sys)),
        (f2_oop, f2_iip); eval_expression, eval_module)

    return f1, f2
end

"""
    $(TYPEDSIGNATURES)

Generate the jacobian of `f1` for [`SemilinearODEFunction`](@ref) (internally represented as a
`SplitFunction`). `A`, `B`, `C` are the matrices returned from
[`calculate_semiquadratic_form`](@ref). `Cjac` is the jacobian of `C` with respect to the
unknowns of the system, or `nothing` if `C === nothing`. This expects that the system has the
necessary extra parmameters added by [`add_semiquadratic_parameters`](@ref).

## Keyword Arguments

$SEMILINEAR_A_B_C_KWARGS
$GENERATE_X_KWARGS

All other keyword arguments are forwarded to [`build_function_wrapper`](@ref).
$SEMILINEAR_A_B_C_CONSTRAINT

$EXPERIMENTAL_WARNING
"""
function generate_semiquadratic_jacobian(
        sys::System, A, B, C, Cjac; sparse = false, stiff_linear = true, stiff_quadratic = false,
        stiff_nonlinear = false, expression = Val{true}, wrap_gfw = Val{false},
        eval_expression = false, eval_module = @__MODULE__, kwargs...)
    if sparse
        error("Sparse analytical jacobians for split ODEs is not implemented.")
    end
    if A === nothing && B === nothing
        throw(ArgumentError("Cannot generate split form for the system - it has no linear or quadratic part."))
    end

    if (stiff_linear || A === nothing) && (stiff_quadratic || B === nothing) &&
       stiff_nonlinear
        throw(ArgumentError("All of `A`, `B` and `C` cannot be stiff at the same time."))
    end
    if (!stiff_linear || A === nothing) && (!stiff_quadratic || B === nothing) &&
       !stiff_nonlinear
        throw(ArgumentError("All of `A`, `B` and `C` cannot be non-stiff at the same time."))
    end
    linear_matrix_param = get_linear_matrix_param_from_sys(sys, A)
    quadratic_forms = get_quadratic_form_params_from_sys(sys, B)
    diffcache_par = get_diffcache_param_from_sys(sys, B)
    eqs = equations(sys)
    dvs = unknowns(sys)
    ps = reorder_parameters(sys)
    iv = get_iv(sys)
    M = length(eqs)
    N = length(dvs)
    # Codegen is a bit manual, and we're manually creating an efficient IIP function.
    # Since we explicitly provide Symbolics.DEFAULT_OUTSYM, the `u` is actually the second
    # argument.
    iip_x = generated_argument_name(2)
    oop_x = generated_argument_name(1)

    # iip
    iip_ir = Assignment[]
    if A !== nothing && stiff_linear
        push!(iip_ir,
            Assignment(
                :__A_jac, term(copyto!, Symbolics.DEFAULT_OUTSYM, linear_matrix_param)))
    end
    if B !== nothing && stiff_quadratic
        cachebuf_name = :__B_cache
        cachelen = M * N
        push!(iip_ir,
            Assignment(cachebuf_name,
                term(PreallocationTools.get_tmp, diffcache_par, Symbolics.DEFAULT_OUTSYM)))
        push!(iip_ir,
            Assignment(
                cachebuf_name, term(reshape, term(view, cachebuf_name, 1:cachelen), M, N)))
        for (i, quadpar) in enumerate(quadratic_forms)
            B[i] === nothing && continue
            coeffvar = Symbol(:__B_matrix_, i)
            # B + B'
            push!(iip_ir, Assignment(coeffvar, term(UpperTriangular, quadpar)))
            push!(iip_ir, Assignment(:__tmp_B_1, term(copyto!, cachebuf_name, coeffvar)))
            # mul! with scalar `B` does addition
            push!(iip_ir,
                Assignment(:__tmp_B_2,
                    term(mul!, cachebuf_name, true, term(transpose, coeffvar), true, true)))
            # view the row of the jacobian this will write to
            target_name = Symbol(:__jac_row_, i)
            push!(
                iip_ir, Assignment(target_name, term(view, Symbolics.DEFAULT_OUTSYM, i, :)))
            # (B + B') * x, written directly to the jacobian. Retain the value in the jacobian if we've already
            # written to it.
            retain_old = A !== nothing && stiff_linear
            push!(iip_ir,
                Assignment(:__tmp_B_3,
                    term(mul!, target_name, cachebuf_name, iip_x, true, retain_old)))
        end
    end
    if C !== nothing && stiff_nonlinear
        @assert Cjac !== nothing
        @assert size(Cjac) == (M, N)
        if A !== nothing && stiff_linear || B !== nothing && stiff_quadratic
            idxs = map(eachindex(Cjac)) do idx
                _iszero(Cjac[idx]) && return nothing
                AtIndex(
                    idx, term(+, term(getindex, Symbolics.DEFAULT_OUTSYM, idx), Cjac[idx]))
            end
            filter!(x -> x !== nothing, idxs)
            push!(iip_ir,
                Assignment(
                    :__tmp_C, SetArray(false, Symbolics.DEFAULT_OUTSYM, idxs, false)))
        end
    end

    # oop
    terms = []
    if A !== nothing && stiff_linear
        push!(terms, linear_matrix_param)
    end
    if B !== nothing && stiff_quadratic
        B_terms = map(eachindex(quadratic_forms)) do i
            B[i] === nothing && return term(FillArrays.Falses, 1, N)
            var = quadratic_forms[i]
            var = term(UpperTriangular, var)
            return term(*, term(+, var, term(transpose, var)), oop_x)
        end
        push!(terms, term(transpose, term(hcat, B_terms...)))
    end
    if C !== nothing && stiff_nonlinear
        push!(terms, MakeArray(Cjac, oop_x))
    end
    oop_expr = length(terms) == 1 ? only(terms) : term(+, terms...)

    j_iip = build_function_wrapper(
        sys, nothing, Symbolics.DEFAULT_OUTSYM, dvs, ps..., iv; p_start = 3,
        extra_assignments = iip_ir, expression = Val{true}, kwargs...)
    j_oop,
    _ = build_function_wrapper(
        sys, oop_expr, dvs, ps..., iv; expression = Val{true}, kwargs...)
    return maybe_compile_function(expression, wrap_gfw, (2, 3, is_split(sys)),
        (j_oop, j_iip); eval_expression, eval_module)
end

"""
    $(TYPEDSIGNATURES)

Return the sparsity pattern of the  jacobian of `f1` for [`SemilinearODEFunction`](@ref)
(internally represented as a `SplitFunction`). `A`, `B`, `C` are the matrices returned from
[`calculate_semiquadratic_form`](@ref). `Cjac` is the jacobian of `C` with respect to the
unknowns of the system, or `nothing` if `C === nothing`. This expects that the system has the
necessary extra parmameters added by [`add_semiquadratic_parameters`](@ref).

## Keyword Arguments

$SEMILINEAR_A_B_C_KWARGS
$GENERATE_X_KWARGS
- `mm`: The mass matrix of `sys`.

$SEMILINEAR_A_B_C_CONSTRAINT

$EXPERIMENTAL_WARNING
"""
function get_semiquadratic_W_sparsity(
        sys::System, A, B, C, Cjac; stiff_linear = true, stiff_quadratic = false,
        stiff_nonlinear = false, mm = calculate_massmatrix(sys))
    eqs = equations(sys)
    dvs = unknowns(sys)
    M = length(eqs)
    N = length(dvs)
    jac = spzeros(Num, M, N)
    if stiff_linear && A !== nothing
        tmp = wrap.(A)
        jac .+= tmp
    end
    if stiff_quadratic && B !== nothing
        for (i, mat) in enumerate(B)
            mat === nothing && continue
            jac[i, :] .+= (mat + transpose(mat)) * dvs
        end
    end
    if stiff_nonlinear && C !== nothing
        jac .+= Cjac
    end
    M_sparsity = mm isa UniformScaling ? sparse(I, M, N) :
                 SparseMatrixCSC{Bool, Int64}((!iszero).(mm))
    return (!_iszero).(jac) .| M_sparsity
end

