using SymbolicUtils.Code

abstract type BuildTargets end
struct JuliaTarget <: BuildTargets end
struct StanTarget <: BuildTargets end
struct CTarget <: BuildTargets end
struct MATLABTarget <: BuildTargets end

abstract type ParallelForm end
struct SerialForm <: ParallelForm end
struct MultithreadedForm <: ParallelForm end
struct DistributedForm <: ParallelForm end
struct DaggerForm <: ParallelForm end

"""
`build_function`

Generates a numerically-usable function from a ModelingToolkit `Num`.

```julia
build_function(ex, args...;
               expression = Val{true},
               target = JuliaTarget(),
               kwargs...)
```

Arguments:

- `ex`: The `Num` to compile
- `args`: The arguments of the function
- `expression`: Whether to generate code or whether to generate the compiled form.
  By default, `expression = Val{true}`, which means that the code for the
  function is returned. If `Val{false}`, then the returned value is compiled.

Keyword Arguments:

- `target`: The output target of the compilation process. Possible options are:
    - `JuliaTarget`: Generates a Julia function
    - `CTarget`: Generates a C function
    - `StanTarget`: Generates a function for compiling with the Stan probabilistic
      programming language
    - `MATLABTarget`: Generates an anonymous function for use in MATLAB and Octave
      environments
- `fname`: Used by some targets for the name of the function in the target space.

Note that not all build targets support the full compilation interface. Check the
individual target documentation for details.
"""
function build_function(args...;target = JuliaTarget(),kwargs...)
  _build_function(target,args...;kwargs...)
end


function add_integrator_header(ex, fargs, iip; X=gensym(:MTIIPVar))
  integrator = gensym(:MTKIntegrator)
  if iip
    wrappedex = :(
        $integrator -> begin
        ($X,$(fargs.args...)) = (($integrator).u,($integrator).u,($integrator).p,($integrator).t)
        $ex
        nothing
      end
    )
  else
    wrappedex = :(
      $integrator -> begin
      ($(fargs.args...),) = (($integrator).u,($integrator).p,($integrator).t)
        $ex
      end
    )
  end
  wrappedex
end

function unflatten_long_ops(op, N=4)
    rule1 = @rule((+)((~~x)) => length(~~x) > N ?
                 +(+((~~x)[1:N]...) + (+)((~~x)[N+1:end]...)) : nothing)
    rule2 = @rule((*)((~~x)) => length(~~x) > N ?
                 *(*((~~x)[1:N]...) * (*)((~~x)[N+1:end]...)) : nothing)

    op = value(op)
    Rewriters.Fixpoint(Rewriters.Postwalk(Rewriters.Chain([rule1, rule2])))(op)
end

function observed_let(eqs)
    process -> ex -> begin
        isempty(eqs) && return ex

        assignments = map(eq -> :($(process(eq.lhs)) = $(process(eq.rhs))), eqs)
        letexpr = :(let $(assignments...)
                    end)
        # avoid a superfluous `begin ... end` block
        letexpr.args[2] = ex
        return letexpr
    end
end

# Scalar output

destructure_arg(arg::Union{AbstractArray, Tuple}) = DestructuredArgs(map(value, arg))
destructure_arg(arg) = value(arg)

function _build_function(target::JuliaTarget, op, args...;
                         conv = toexpr,
                         expression = Val{true},
                         checkbounds = false,
                         linenumbers = true)

    dargs = map(destructure_arg, args)
    expr = toexpr(Func(dargs, [], value(op)))

    if expression == Val{true}
        expr
    else
        _build_and_inject_function(@__MODULE__, expr)
    end
end

function _build_and_inject_function(mod::Module, ex)
    @RuntimeGeneratedFunction(ex)
end

# Detect heterogeneous element types of "arrays of matrices/sparce matrices"
function is_array_matrix(F)
    return isa(F, AbstractVector) && all(x->isa(x, AbstractArray), F)
end
function is_array_sparse_matrix(F)
    return isa(F, AbstractVector) && all(x->isa(x, AbstractSparseMatrix), F)
end
# Detect heterogeneous element types of "arrays of arrays of matrices/sparce matrices"
function is_array_array_matrix(F)
    return isa(F, AbstractVector) && all(x->isa(x, AbstractArray{<:AbstractMatrix}), F)
end
function is_array_array_sparse_matrix(F)
    return isa(F, AbstractVector) && all(x->isa(x, AbstractArray{<:AbstractSparseMatrix}), F)
end

toexpr(n::Num, st) = toexpr(value(n), st)

function fill_array_with_zero!(x::AbstractArray)
    if eltype(x) <: AbstractArray
        foreach(fill_array_with_zero!, x)
    else
        fill!(x, false)
    end
    return x
end

"""
Build function target: JuliaTarget

```julia
function _build_function(target::JuliaTarget, rhss, args...;
                         conv = toexpr, expression = Val{true},
                         checkbounds = false,
                         linenumbers = false, multithread=nothing,
                         headerfun = addheader, outputidxs=nothing,
                         convert_oop = true, force_SA = false,
                         skipzeros = outputidxs===nothing,
                         fillzeros = skipzeros && !(typeof(rhss)<:SparseMatrixCSC),
                         parallel=SerialForm(), kwargs...)
```

Generates a Julia function which can then be utilized for further evaluations.
If expression=Val{false}, the return is a Julia function which utilizes
RuntimeGeneratedFunctions.jl in order to be free of world-age issues.

If the `Num` is an `Operation`, the generated function is a function
with a scalar output, otherwise if it's an `AbstractArray{Operation}`, the output
is two functions, one for out-of-place AbstractArray output and a second which
is a mutating function. The outputted functions match the given argument order,
i.e., f(u,p,args...) for the out-of-place and scalar functions and
`f!(du,u,p,args..)` for the in-place version.

Special Keyword Argumnets:

- `parallel`: The kind of parallelism to use in the generated function. Defaults
  to `SerialForm()`, i.e. no parallelism. Note that the parallel forms are not
  exported and thus need to be chosen like `ModelingToolkit.SerialForm()`.
  The choices are:
  - `SerialForm()`: Serial execution.
  - `MultithreadedForm()`: Multithreaded execution with a static split, evenly
    splitting the number of expressions per thread.
  - `DistributedForm()`: Multiprocessing using Julia's Distributed with a static
    schedule, evenly splitting the number of expressions per process.
  - `DaggerForm()`: Multithreading and multiprocessing using Julia's Dagger.jl
    for dynamic scheduling and load balancing.
- `conv`: The conversion function of the Operation to Expr. By default this uses
  the `toexpr` function.
- `checkbounds`: For whether to enable bounds checking inside of the generated
  function. Defaults to false, meaning that `@inbounds` is applied.
- `linenumbers`: Determines whether the generated function expression retains
  the line numbers. Defaults to true.
- `convert_oop`: Determines whether the OOP version should try to convert
  the output to match the type of the first input. This is useful for
  cases like LabelledArrays or other array types that carry extra
  information. Defaults to true.
- `force_SA`: Forces the output of the OOP version to be a StaticArray.
  Defaults to `false`, and outputs a static array when the first argument
  is a static array.
- `skipzeros`: Whether to skip filling zeros in the in-place version if the
  filling function is 0.
- `fillzeros`: Whether to perform `fill(out,0)` before the calculations to ensure
  safety with `skipzeros`.
  """
  function _build_function(target::JuliaTarget, rhss::AbstractArray, args...;
                           conv = toexpr, expression = Val{true},
                           checkbounds = false,
                           linenumbers = false, multithread=nothing,
                           outputidxs=nothing,
                           skipzeros = false,
                           fillzeros = skipzeros && !(typeof(rhss)<:SparseMatrixCSC),
                           parallel=SerialForm(), kwargs...)

    dargs = map(destructure_arg, args)
    i = findfirst(x->x isa DestructuredArgs, dargs)
    similarto = i === nothing ? Array : dargs[i].name
    array_expr = _make_array(rhss, similarto)
    oop_expr = Func(dargs, [], array_expr)

    out = Sym{Any}(gensym("out"))
    if rhss isa SparseMatrixCSC
        I,J, _ = findnz(rhss)
        outputidxs = CartesianIndex.(I, J)
    elseif rhss isa SparseVector
        I,_ = findnz(rhss)
        outputidxs = I
    elseif isnothing(outputidxs)
        outputidxs = collect(eachindex(rhss))
    end

    if skipzeros
        ii = findall(i->!_iszero(rhss[i]), outputidxs)
        array = AtIndex.(outputidxs[ii], rhss[ii])
    else
        array = AtIndex.(outputidxs, rhss)
    end
    ip_expr = Func([out, dargs...], [], SetArray(false, out, array))

    if expression == Val{true}
        return toexpr(oop_expr), toexpr(ip_expr)
    else
        return _build_and_inject_function(@__MODULE__, toexpr(oop_expr)),
               _build_and_inject_function(@__MODULE__, toexpr(ip_expr))
    end
end

function _make_array(rhss::AbstractSparseArray, similarto)
    MakeSparseArray(map(x->_make_array(x, similarto), rhss))
end

function _make_array(rhss::AbstractArray, similarto)
    MakeArray(map(x->_make_array(x, similarto), rhss), similarto)
end

_make_array(x, similarto) = x

function vars_to_pairs(name,vs::Union{Tuple, AbstractArray}, symsdict=Dict())
    vs_names = tosymbol.(vs)
    for (v,k) in zip(vs_names, vs)
        symsdict[k] = Sym{symtype(k)}(v)
    end
    exs = [:($name[$i]) for (i, u) ∈ enumerate(vs)]
    vs_names,exs
end
function vars_to_pairs(name,vs, symsdict)
    symsdict[vs] = Sym{symtype(vs)}(tosymbol(vs))
    [tosymbol(vs)], [name]
end

get_varnumber(varop, vars::Vector) =  findfirst(x->isequal(x,varop),vars)

function numbered_expr(O::Symbolic,args...;varordering = args[1],offset = 0,
                       lhsname=gensym("du"),rhsnames=[gensym("MTK") for i in 1:length(args)])
    O = value(O)
    if O isa Sym || isa(operation(O), Sym)
        for j in 1:length(args)
            i = get_varnumber(O,args[j])
            if i !== nothing
                return :($(rhsnames[j])[$(i+offset)])
            end
        end
    end
  return Expr(:call, O isa Sym ? tosymbol(O, escape=false) : Symbol(operation(O)),
         [numbered_expr(x,args...;offset=offset,lhsname=lhsname,
                        rhsnames=rhsnames,varordering=varordering) for x in arguments(O)]...)
end

function numbered_expr(de::ModelingToolkit.Equation,args...;varordering = args[1],
                       lhsname=gensym("du"),rhsnames=[gensym("MTK") for i in 1:length(args)],offset=0)

    varordering = value.(args[1])
    var = var_from_nested_derivative(de.lhs)[1]
    i = findfirst(x->isequal(tosymbol(x isa Sym ? x : operation(x), escape=false), tosymbol(var, escape=false)),varordering)
    :($lhsname[$(i+offset)] = $(numbered_expr(de.rhs,args...;offset=offset,
                                              varordering = varordering,
                                              lhsname = lhsname,
                                              rhsnames = rhsnames)))
end
numbered_expr(c,args...;kwargs...) = c
numbered_expr(c::Num,args...;kwargs...) = error("Num found")

"""
Build function target: CTarget

```julia
function _build_function(target::CTarget, eqs::Array{<:Equation}, args...;
                         conv = toexpr, expression = Val{true},
                         fname = :diffeqf,
                         lhsname=:du,rhsnames=[Symbol("RHS\$i") for i in 1:length(args)],
                         libpath=tempname(),compiler=:gcc)
```

This builds an in-place C function. Only works on arrays of equations. If
`expression == Val{false}`, then this builds a function in C, compiles it,
and returns a lambda to that compiled function. These special keyword arguments
control the compilation:

- libpath: the path to store the binary. Defaults to a temporary path.
- compiler: which C compiler to use. Defaults to :gcc, which is currently the
  only available option.
"""
function _build_function(target::CTarget, eqs::Array{<:Equation}, args...;
                         conv = toexpr, expression = Val{true},
                         fname = :diffeqf,
                         lhsname=:du,rhsnames=[Symbol("RHS$i") for i in 1:length(args)],
                         libpath=tempname(),compiler=:gcc)

    differential_equation = string(join([numbered_expr(eq,args...,lhsname=lhsname,
                                  rhsnames=rhsnames,offset=-1) for
                                  (i, eq) ∈ enumerate(eqs)],";\n  "),";")

    argstrs = join(vcat("double* $(lhsname)",[typeof(args[i])<:Array ? "double* $(rhsnames[i])" : "double $(rhsnames[i])" for i in 1:length(args)]),", ")
    ex = """
    void $fname($(argstrs...)) {
      $differential_equation
    }
    """

    if expression == Val{true}
        return ex
    else
        @assert compiler == :gcc
        ex = build_function(eqs,args...;target=ModelingToolkit.CTarget())
        open(`gcc -fPIC -O3 -msse3 -xc -shared -o $(libpath * "." * Libdl.dlext) -`, "w") do f
            print(f, ex)
        end
        @RuntimeGeneratedFunction(:((du::Array{Float64},u::Array{Float64},p::Array{Float64},t::Float64) -> ccall(("diffeqf", $libpath), Cvoid, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64), du, u, p, t)))
    end
end

"""
Build function target: StanTarget

```julia
function _build_function(target::StanTarget, eqs::Array{<:Equation}, vs, ps, iv;
                         conv = toexpr, expression = Val{true},
                         fname = :diffeqf, lhsname=:internal_var___du,
                         rhsnames=[:internal_var___u,:internal_var___p,:internal_var___t])
```

This builds an in-place Stan function compatible with the Stan differential equation solvers.
Unlike other build targets, this one requestions (vs, ps, iv) as the function arguments.
Only allowed on arrays of equations.
"""
function _build_function(target::StanTarget, eqs::Array{<:Equation}, vs, ps, iv;
                         conv = toexpr, expression = Val{true},
                         fname = :diffeqf, lhsname=:internal_var___du,
                         rhsnames=[:internal_var___u,:internal_var___p,:internal_var___t])
    @assert expression == Val{true}
    differential_equation = string(join([numbered_expr(eq,vs,ps,lhsname=lhsname,
                                   rhsnames=rhsnames) for
                                   (i, eq) ∈ enumerate(eqs)],";\n  "),";")
    """
    real[] $fname(real $(conv(iv)),real[] $(rhsnames[1]),real[] $(rhsnames[2]),real[] x_r,int[] x_i) {
      real $lhsname[$(length(eqs))];
      $differential_equation
      return $lhsname;
    }
    """
end

"""
Build function target: MATLABTarget

```julia
function _build_function(target::MATLABTarget, eqs::Array{<:Equation}, args...;
                         conv = toexpr, expression = Val{true},
                         lhsname=:internal_var___du,
                         rhsnames=[:internal_var___u,:internal_var___p,:internal_var___t])
```

This builds an out of place anonymous function @(t,rhsnames[1]) to be used in MATLAB.
Compatible with the MATLAB differential equation solvers. Only allowed on arrays
of equations.
"""
function _build_function(target::MATLABTarget, eqs::Array{<:Equation}, args...;
                         conv = toexpr, expression = Val{true},
                         fname = :diffeqf, lhsname=:internal_var___du,
                         rhsnames=[:internal_var___u,:internal_var___p,:internal_var___t])
    @assert expression == Val{true}
    matstr = join([numbered_expr(eq.rhs,args...,lhsname=lhsname,
                                  rhsnames=rhsnames) for
                                  (i, eq) ∈ enumerate(eqs)],"; ")

    matstr = replace(matstr,"["=>"(")
    matstr = replace(matstr,"]"=>")")
    matstr = "$fname = @(t,$(rhsnames[1])) ["*matstr*"];"
    matstr
end
