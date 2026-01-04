struct SemilinearODEFunction{iip, spec} end
struct SemilinearODEProblem{iip, spec} end

const SEMILINEAR_EXTRA_BODY = """
This is a special form of an ODE which uses a `SplitFunction` internally. The equations are
separated into linear, quadratic and general terms and phrased as matrix operations. See
[`calculate_semiquadratic_form`](@ref) for information on how the equations are split. This
formulation allows leveraging split ODE solvers such as `KenCarp4` and is useful for systems
where the stiff and non-stiff terms can be separated out in such a manner. Typically the linear
part of the equations is the stiff part, but the keywords `stiff_linear`, `stiff_quadratic` and `stiff_nonlinear` can
be used to control which parts are considered as stiff.
"""

const SEMILINEAR_A_B_C_KWARGS = """
- `stiff_linear`: Whether the linear part of the equations should be part of the stiff function
  in the split form. Has no effect if the equations have no linear part.
- `stiff_quadratic`: Whether the quadratic part of the equations should be part of the stiff
  function in the split form. Has no effect if the equations have no quadratic part.
- `stiff_nonlinear`: Whether the non-linear non-quadratic part of the equations should be part of
  the stiff function in the split form. Has no effect if the equations have no such
  non-linear non-quadratic part.
"""

const SEMILINEAR_A_B_C_CONSTRAINT = """
Note that all three of `stiff_linear`, `stiff_quadratic`, `stiff_nonlinear` cannot be identical, and at least
two of `A`, `B`, `C` returned from [`calculate_semiquadratic_form`](@ref) must be
non-`nothing`. In other words, both of the functions in the split form must be non-empty.
"""

for (mod, prob, func, istd, kws) in [
        (SciMLBase, :SCCNonlinearProblem, NonlinearFunction, false, (; init = false)),
        (
            ModelingToolkit,
            :SemilinearODEProblem,
            :SemilinearODEFunction,
            true,
            (;
                extra_body = SEMILINEAR_EXTRA_BODY, extra_kwargs = SEMILINEAR_A_B_C_KWARGS,
                extra_kwargs_desc = SEMILINEAR_A_B_C_CONSTRAINT,
            ),
        ),
    ]
    kwexpr = Expr(:parameters)
    for (k, v) in pairs(kws)
        push!(kwexpr.args, Expr(:kw, k, v))
    end
    @eval @doc MTKBase.problem_docstring($kwexpr, $mod.$prob, $func, $istd) $mod.$prob
end

for (mod, func, istd, optionals, kws) in [
        (
            ModelingToolkit,
            :SemilinearODEFunction,
            true,
            [:jac],
            (;
                extra_body = SEMILINEAR_EXTRA_BODY, extra_kwargs = SEMILINEAR_A_B_C_KWARGS,
                extra_kwargs_desc = SEMILINEAR_A_B_C_CONSTRAINT,
            ),
        ),
    ]
    kwexpr = Expr(:parameters)
    for (k, v) in pairs(kws)
        push!(kwexpr.args, Expr(:kw, k, v))
    end
    @eval @doc MTKBase.function_docstring($kwexpr, $mod.$func, $istd, $optionals) $mod.$func
end
