abstract type BuildTargets end
struct JuliaTarget <: BuildTargets end
struct StanTarget <: BuildTargets end
struct CTarget <: BuildTargets end
struct MATLABTarget <: BuildTargets end

"""
`build_function`

Generates a numerically-usable function from a ModelingToolkit `Expression`.
If the `Expression` is an `Operation`, the generated function is a function
with a scalar output, otherwise if it's an `AbstractArray{Operation}` the output
is two functions, one for out-of-place AbstractArray output and a second which
is a mutating function. The outputted functions match the given argument order,
i.e. f(u,p,args...) for the out-of-place and scalar functions and
`f!(du,u,p,args..)` for the in-place version.

```julia
build_function(ex vs, ps = (), args = (),
               conv = simplified_expr, expression = Val{true};
               checkbounds = false, constructor=nothing,
               linenumbers = true, target = JuliaTarget())
```

Arguments:

- `ex`: The `Expression` to compile
- `vs`: The variables of the expression
- `ps`: The parameters of the expression
- `args`: Extra arguments to the function
- `conv`: The conversion function of the Operation to Expr. By default this uses
  the `simplified_expr` function utilized in `convert(Expr,x)`.
- `expression`: Whether to generate code or whether to generate the compiled form.
  By default, `expression = Val{true}`, which means that the code for the
  function is returned. If `Val{false}`, then the returned value is a compiled
  Julia function which utilizes GeneralizedGenerated.jl in order to world-age
  free.

Keyword Arguments:

- `checkbounds`: For whether to enable bounds checking inside of the generated
  function. Defaults to false, meaning that `@inbounds` is applied.
- `constructor`: Allows for an arbitrary constructor function to be passed in
  for handling expressions of "weird" types. Defaults to nothing.
- `linenumbers`: Determines whether the generated function expression retains
  the line numbers. Defaults to true.
- `target`: The output target of the compilation process. Possible options are:
    - `JuliaTarget`: Generates a Julia function
    - `CTarget`: Generates a C function
    - `StanTarget`: Generates a function for compiling with the Stan probabilistic
      programming language
    - `MATLABTarget`: Generates an anonymous function for use in MATLAB and Octave
      environments
"""
function build_function(args...;target = JuliaTarget(),kwargs...)
  _build_function(target,args...;kwargs...)
end

# Scalar output
function _build_function(target::JuliaTarget, op::Operation, vs, ps = (), args = (),
                         conv = simplified_expr, expression = Val{true};
                         checkbounds = false, constructor=nothing,
                         linenumbers = true)
    _vs = convert.(Variable,vs)
    _ps = convert.(Variable,ps)
    var_pairs   = [(u.name, :(u[$i])) for (i, u) ∈ enumerate(_vs)]
    param_pairs = [(p.name, :(p[$i])) for (i, p) ∈ enumerate(_ps)]
    (ls, rs) = zip(var_pairs..., param_pairs...)
    var_eqs = Expr(:(=), build_expr(:tuple, ls), build_expr(:tuple, rs))

    fname = gensym(:ModelingToolkitFunction)
    out_expr = conv(op)
    let_expr = Expr(:let, var_eqs, out_expr)
    bounds_block = checkbounds ? let_expr : :(@inbounds begin $let_expr end)

    fargs = ps == () ? :(u,$(args...)) : :(u,p,$(args...))

    oop_ex = :(
        ($(fargs.args...),) -> begin
            $bounds_block
        end
    )

    if !linenumbers
        oop_ex = striplines(oop_ex)
    end

    if expression == Val{true}
        return oop_ex
    else
        return GeneralizedGenerated.mk_function(@__MODULE__,oop_ex)
    end
end

function _build_function(target::JuliaTarget, rhss, vs, ps = (), args = (),
                         conv = simplified_expr, expression = Val{true};
                         checkbounds = false, constructor=nothing,
                         linenumbers = true, multithread=false)
    _vs = convert.(Variable,vs)
    _ps = convert.(Variable,ps)
    var_pairs   = [(u.name, :(u[$i])) for (i, u) ∈ enumerate(_vs)]
    param_pairs = [(p.name, :(p[$i])) for (i, p) ∈ enumerate(_ps)]
    (ls, rs) = zip(var_pairs..., param_pairs...)

    var_eqs = Expr(:(=), build_expr(:tuple, ls), build_expr(:tuple, rs))

    fname = gensym(:ModelingToolkitFunction)

    X = gensym(:MTIIPVar)
    if rhss isa SparseMatrixCSC
        ip_sys_exprs = [:($X.nzval[$i] = $(conv(rhs))) for (i, rhs) ∈ enumerate(rhss.nzval)]
    else
        ip_sys_exprs = [:($X[$i] = $(conv(rhs))) for (i, rhs) ∈ enumerate(rhss)]
    end

    ip_let_expr = Expr(:let, var_eqs, build_expr(:block, ip_sys_exprs))

    if multithread
        lens = Int(ceil(length(ip_let_expr.args[2].args)/Threads.nthreads()))
        threaded_exprs = vcat([quote
           Threads.@spawn begin
              $(ip_let_expr.args[2].args[((i-1)*lens+1):i*lens]...)
           end
        end for i in 1:Threads.nthreads()-1],
           quote
              Threads.@spawn begin
                 $(ip_let_expr.args[2].args[((Threads.nthreads()-1)*lens+1):end]...)
              end
           end)
        ip_let_expr.args[2] =  ModelingToolkit.build_expr(:block, threaded_exprs)
    end

    tuple_sys_expr = build_expr(:tuple, [conv(rhs) for rhs ∈ rhss])

    if rhss isa Matrix
        arr_sys_expr = build_expr(:vcat, [build_expr(:row,[conv(rhs) for rhs ∈ rhss[i,:]]) for i in 1:size(rhss,1)])
        # : x because ??? what to do in the general case?
        _constructor = constructor === nothing ? :(u isa ModelingToolkit.StaticArrays.StaticArray ? ModelingToolkit.StaticArrays.SMatrix{$(size(rhss)...)} :  x->(out = similar(typeof(u),$(size(rhss)...)); out .= x)) : constructor
    elseif typeof(rhss) <: Array && !(typeof(rhss) <: Vector)
        vector_form = build_expr(:vect, [conv(rhs) for rhs ∈ rhss])
        arr_sys_expr = :(reshape($vector_form,$(size(rhss)...)))
        _constructor = constructor === nothing ? :(u isa ModelingToolkit.StaticArrays.StaticArray ? ModelingToolkit.StaticArrays.SArray{$(size(rhss)...)} :  x->(out = similar(typeof(u),$(size(rhss)...)); out .= x)) : constructor
    elseif rhss isa SparseMatrixCSC
        vector_form = build_expr(:vect, [conv(rhs) for rhs ∈ nonzeros(rhss)])
        arr_sys_expr = :(SparseMatrixCSC{eltype(u),Int}($(size(rhss)...), $(rhss.colptr), $(rhss.rowval), $vector_form))
        # Static and sparse? Probably not a combo that will actually be hit, but give a default anyways
        _constructor = constructor === nothing ? :(u isa ModelingToolkit.StaticArrays.StaticArray ? ModelingToolkit.StaticArrays.SMatrix{$(size(rhss)...)} : x->x) : constructor
    else # Vector
        arr_sys_expr = build_expr(:vect, [conv(rhs) for rhs ∈ rhss])
        # Handle vector constructor separately using `typeof(u)` to support things like LabelledArrays
        _constructor = constructor === nothing ? :(u isa ModelingToolkit.StaticArrays.StaticArray ? ModelingToolkit.StaticArrays.similar_type(typeof(u), eltype(X)) : x->convert(typeof(u),x)) : constructor
    end

    let_expr = Expr(:let, var_eqs, tuple_sys_expr)
    arr_let_expr = Expr(:let, var_eqs, arr_sys_expr)
    bounds_block = checkbounds ? let_expr : :(@inbounds begin $let_expr end)
    arr_bounds_block = checkbounds ? arr_let_expr : :(@inbounds begin $arr_let_expr end)
    ip_bounds_block = checkbounds ? ip_let_expr : :(@inbounds begin $ip_let_expr end)

    fargs = ps == () ? :(u,$(args...)) : :(u,p,$(args...))

    oop_ex = :(
        ($(fargs.args...),) -> begin
            # If u is a weird non-StaticArray type and we want a sparse matrix, just do the optimized sparse anyways
            if $(fargs.args[1]) isa Array || (!(typeof($(fargs.args[1])) <: StaticArray) && $(rhss isa SparseMatrixCSC))
                return $arr_bounds_block
            else
                X = $bounds_block
                construct = $_constructor
                return construct(X)
            end
        end
    )

    iip_ex = :(
        ($X,$(fargs.args...)) -> begin
            $ip_bounds_block
            nothing
        end
    )

    if !linenumbers
        oop_ex = striplines(oop_ex)
        iip_ex = striplines(iip_ex)
    end

    if expression == Val{true}
        return oop_ex, iip_ex
    else
        return GeneralizedGenerated.mk_function(@__MODULE__,oop_ex), GeneralizedGenerated.mk_function(@__MODULE__,iip_ex)
    end
end

get_varnumber(varop::Operation,vars::Vector{Operation}) =  findfirst(x->isequal(x,varop),vars)
get_varnumber(varop::Operation,vars::Vector{<:Variable})  =  findfirst(x->isequal(x,varop.op),vars)

function numbered_expr(O::Equation,args...;kwargs...)
  :($(numbered_expr(O.lhs,args...;kwargs...)) = $(numbered_expr(O.rhs,args...;kwargs...)))
end

function numbered_expr(O::Operation,vars,parameters;
                       derivname=:du,
                       varname=:u,paramname=:p)
  if isa(O.op, ModelingToolkit.Differential)
    varop = O.args[1]
    i = get_varnumber(varop,vars)
    return :($derivname[$i])
  elseif isa(O.op, ModelingToolkit.Variable)
    i = get_varnumber(O,vars)
    if i == nothing
      i = get_varnumber(O,parameters)
      return :($paramname[$i])
    else
      return :($varname[$i])
    end
  end
  return Expr(:call, Symbol(O.op),
         [numbered_expr(x,vars,parameters;derivname=derivname,
                        varname=varname,paramname=paramname) for x in O.args]...)
end

function numbered_expr(de::ModelingToolkit.Equation,vars::Vector{<:Variable},parameters;
                       derivname=:du,varname=:u,paramname=:p)
    i = findfirst(x->isequal(x.name,var_from_nested_derivative(de.lhs)[1].name),vars)
    :($derivname[$i] = $(numbered_expr(de.rhs,vars,parameters;
                                     derivname=derivname,
                                     varname=varname,paramname=paramname)))
end
function numbered_expr(de::ModelingToolkit.Equation,vars::Vector{Operation},parameters;
                       derivname=:du,varname=:u,paramname=:p)
    i = findfirst(x->isequal(x.op.name,var_from_nested_derivative(de.lhs)[1].name),vars)
    :($derivname[$i] = $(numbered_expr(de.rhs,vars,parameters;
                                     derivname=derivname,
                                     varname=varname,paramname=paramname)))
end
numbered_expr(c::ModelingToolkit.Constant,args...;kwargs...) = c.value

function _build_function(target::StanTarget, eqs, vs, ps, iv,
                         conv = simplified_expr, expression = Val{true};
                         fname = :diffeqf, derivname=:internal_var___du,
                         varname=:internal_var___u,paramname=:internal_var___p)
    differential_equation = string(join([numbered_expr(eq,vs,ps,derivname=derivname,
                                   varname=varname,paramname=paramname) for
                                   (i, eq) ∈ enumerate(eqs)],";\n  "),";")
    """
    real[] $fname(real $iv,real[] $varname,real[] $paramname,real[] x_r,int[] x_i) {
      real $derivname[$(length(eqs))];
      $differential_equation
      return $derivname;
    }
    """
end

function _build_function(target::CTarget, eqs, vs, ps, iv,
                         conv = simplified_expr, expression = Val{true};
                         fname = :diffeqf, derivname=:internal_var___du,
                         varname=:internal_var___u,paramname=:internal_var___p)
    differential_equation = string(join([numbered_expr(eq,vs,ps,derivname=derivname,
                                  varname=varname,paramname=paramname) for
                                  (i, eq) ∈ enumerate(eqs)],";\n  "),";")
    """
    void $fname(double* $derivname, double* $varname, double* $paramname, $iv) {
      $differential_equation
    }
    """
end

function _build_function(target::MATLABTarget, eqs, vs, ps, iv,
                         conv = simplified_expr, expression = Val{true};
                         fname = :diffeqf, derivname=:internal_var___du,
                         varname=:internal_var___u,paramname=:internal_var___p)
    matstr = join([numbered_expr(eq.rhs,vs,ps,derivname=derivname,
                                  varname=varname,paramname=paramname) for
                                  (i, eq) ∈ enumerate(eqs)],"; ")

    matstr = replace(matstr,"["=>"(")
    matstr = replace(matstr,"]"=>")")
    matstr = "$fname = @(t,$varname) ["*matstr*"];"
    matstr
end
