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

Generates a numerically-usable function from a ModelingToolkit `Expression`.

```julia
build_function(ex, args...;
               expression = Val{true},
               target = JuliaTarget(),
			   kwargs...)
```

Arguments:

- `ex`: The `Expression` to compile
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

function addheader(ex, fargs, iip; X=gensym(:MTIIPVar))
  if iip
    wrappedex = :(
        ($X,$(fargs.args...)) -> begin
        $ex
        nothing
      end
    )
  else
    wrappedex = :(
      ($(fargs.args...),) -> begin
        $ex
      end
    )
  end
  wrappedex
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

# Scalar output
function _build_function(target::JuliaTarget, op, args...;
                         conv = toexpr, expression = Val{true},
                         checkbounds = false,
                         linenumbers = true, headerfun=addheader)

    argnames = [gensym(:MTKArg) for i in 1:length(args)]
    arg_pairs = map(vars_to_pairs,zip(argnames,args))
    ls = reduce(vcat,first.(arg_pairs))
    rs = reduce(vcat,last.(arg_pairs))
    var_eqs = Expr(:(=), ModelingToolkit.build_expr(:tuple, ls), ModelingToolkit.build_expr(:tuple, rs))

    fname = gensym(:ModelingToolkitFunction)
    out_expr = conv(op)
    let_expr = Expr(:let, var_eqs, Expr(:block, out_expr))
    bounds_block = checkbounds ? let_expr : :(@inbounds begin $let_expr end)

    fargs = Expr(:tuple,argnames...)
    oop_ex = headerfun(bounds_block, fargs, false)

    if !linenumbers
        oop_ex = striplines(oop_ex)
    end

    if expression == Val{true}
        return ModelingToolkit.inject_registered_module_functions(oop_ex)
    else
        _build_and_inject_function(@__MODULE__, oop_ex)
    end
end

function _build_and_inject_function(mod::Module, ex)
    # Generate the function, which will process the expression
    runtimefn = GeneralizedGenerated.mk_function(mod, ex)

    # Extract the processed expression of the function body
    params = typeof(runtimefn).parameters
    fn_expr = GeneralizedGenerated.NGG.from_type(params[3])

    # Inject our externally registered module functions
    new_expr = ModelingToolkit.inject_registered_module_functions(fn_expr)

    # Reconstruct the RuntimeFn's Body
    new_body = GeneralizedGenerated.NGG.to_type(new_expr)
    return GeneralizedGenerated.RuntimeFn{params[1:2]..., new_body, params[4]}()
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
GeneralizedGenerated.jl in order to be free of world-age issues.

If the `Expression` is an `Operation`, the generated function is a function
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
  -	`SerialForm()`: Serial execution.
  - `MultithreadedForm()`: Multithreaded execution with a static split, evenly
    splitting the number of expressions per thread.
  - `DistributedForm()`: Multiprocessing using Julia's Distributed with a static
    schedule, evenly splitting the number of expressions per process.
  - `DaggerForm()`: Multithreading and multiprocessing using Julia's Dagger.jl
    for dynamic scheduling and load balancing.
- `conv`: The conversion function of the Operation to Expr. By default this uses
  the `toexpr` function utilized in `convert(Expr,x)`.
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
                         headerfun = addheader, outputidxs=nothing,
						 convert_oop = true, force_SA = false,
                         skipzeros = outputidxs===nothing,
						 fillzeros = skipzeros && !(typeof(rhss)<:SparseMatrixCSC),
						 parallel=SerialForm(), kwargs...)
    conv = conv ∘ rm_calls_with_iv
	if multithread isa Bool
		@warn("multithraded is deprecated for the parallel argument. See the documentation.")
		parallel = multithread ? MultithreadedForm() : SerialForm()
	end

    argnames = [gensym(:MTKArg) for i in 1:length(args)]
    arg_pairs = map(vars_to_pairs,zip(argnames,args))
    ls = reduce(vcat,first.(arg_pairs))
    rs = reduce(vcat,last.(arg_pairs))
    var_eqs = Expr(:(=), ModelingToolkit.build_expr(:tuple, ls), ModelingToolkit.build_expr(:tuple, rs))

    fname = gensym(:ModelingToolkitFunction)
    fargs = Expr(:tuple,argnames...)


    oidx = isnothing(outputidxs) ? (i -> i) : (i -> outputidxs[i])
    X = gensym(:MTIIPVar)

	rhs_length = rhss isa SparseMatrixCSC ? length(rhss.nzval) : length(rhss)

	if parallel isa DistributedForm
		numworks = Distributed.nworkers()
		reducevars = [Variable(gensym(:MTReduceVar))() for i in 1:numworks]
		lens = Int(ceil(rhs_length/numworks))
		finalsize = rhs_length - (numworks-1)*lens
		_rhss = vcat(reduce(vcat,[[getindex(reducevars[i],j) for j in 1:lens] for i in 1:numworks-1],init=Expr[]),
						 [getindex(reducevars[end],j) for j in 1:finalsize])
    elseif parallel isa DaggerForm
		computevars = [Variable(gensym(:MTComputeVar))() for i in axes(rhss,1)]
        reducevar = Variable(gensym(:MTReduceVar))()
        _rhss = [getindex(reducevar,i) for i in axes(rhss,1)]
	elseif rhss isa SparseMatrixCSC
		_rhss = rhss.nzval
	else
		_rhss = rhss
	end

    ip_sys_exprs = Expr[]
    # we cannot reliably fill the array with the presence of index translation
    if is_array_array_sparse_matrix(rhss) # Array of arrays of sparse matrices
        for (i, rhsel) ∈ enumerate(_rhss)
            for (j, rhsel2) ∈ enumerate(rhsel)
                for (k, rhs) ∈ enumerate(rhsel2.nzval)
                    rhs′ = conv(rhs)
                    (skipzeros && rhs′ isa Number && iszero(rhs′)) && continue
                    push!(ip_sys_exprs, :($X[$i][$j].nzval[$k] = $rhs′))
                end
            end
        end
    elseif is_array_array_matrix(rhss) # Array of arrays of arrays
        for (i, rhsel) ∈ enumerate(_rhss)
            for (j, rhsel2) ∈ enumerate(rhsel)
                for (k, rhs) ∈ enumerate(rhsel2)
                    rhs′ = conv(rhs)
                    (skipzeros && rhs′ isa Number && iszero(rhs′)) && continue
                    push!(ip_sys_exprs, :($X[$i][$j][$k] = $rhs′))
                end
            end
        end
    elseif is_array_sparse_matrix(rhss) # Array of sparse matrices
        for (i, rhsel) ∈ enumerate(_rhss)
            for (j, rhs) ∈ enumerate(rhsel.nzval)
                rhs′ = conv(rhs)
                (skipzeros && rhs′ isa Number && iszero(rhs′)) && continue
                push!(ip_sys_exprs, :($X[$i].nzval[$j] = $rhs′))
            end
        end
    elseif is_array_matrix(rhss) # Array of arrays
        for (i, rhsel) ∈ enumerate(_rhss)
            for (j, rhs) ∈ enumerate(rhsel)
                rhs′ = conv(rhs)
                (skipzeros && rhs′ isa Number && iszero(rhs′)) && continue
                push!(ip_sys_exprs, :($X[$i][$j] = $rhs′))
            end
        end
    elseif rhss isa SparseMatrixCSC
        for (i, rhs) ∈ enumerate(_rhss)
            rhs′ = conv(rhs)
            (skipzeros && rhs′ isa Number && iszero(rhs′)) && continue
            push!(ip_sys_exprs, :($X.nzval[$i] = $rhs′))
        end
    else
        for (i, rhs) ∈ enumerate(_rhss)
            rhs′ = conv(rhs)
            (skipzeros && rhs′ isa Number && iszero(rhs′)) && continue
            push!(ip_sys_exprs, :($X[$(oidx(i))] = $rhs′))
        end
    end

    ip_let_expr = Expr(:let, var_eqs, build_expr(:block, ip_sys_exprs))

    if parallel isa MultithreadedForm
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
		ip_let_expr = :(@sync begin $ip_let_expr end)
    elseif parallel isa DistributedForm
		numworks = Distributed.nworkers()
		lens = Int(ceil(length(ip_let_expr.args[2].args)/numworks))
		spawnvars = [gensym(:MTSpawnVar) for i in 1:numworks]
		rhss_flat = rhss isa SparseMatrixCSC ? rhss.nzval : rhss
		spawnvectors = vcat(
					   [build_expr(:vect, [conv(rhs) for rhs ∈ rhss_flat[((i-1)*lens+1):i*lens]]) for i in 1:numworks-1],
					   build_expr(:vect, [conv(rhs) for rhs ∈ rhss_flat[((numworks-1)*lens+1):end]]))

        spawn_exprs = [quote
           $(spawnvars[i]) = ModelingToolkit.Distributed.remotecall($(i+1)) do
              $(spawnvectors[i])
           end
        end for i in 1:numworks]
        spawn_exprs = ModelingToolkit.build_expr(:block, spawn_exprs)
		resunpack_exprs = [:($(Symbol(reducevars[iter])) = fetch($(spawnvars[iter]))) for iter in 1:numworks]

		ip_let_expr.args[2] = quote
			@sync begin
				$spawn_exprs
				$(resunpack_exprs...)
				$(ip_let_expr.args[2])
			end
		end
    elseif parallel isa DaggerForm
        @assert HAS_DAGGER[] "Dagger.jl is not loaded; please do `using Dagger`"
        dagwrap(x) = x
        dagwrap(ex::Expr) = dagwrap(ex, Val(ex.head))
        dagwrap(ex::Expr, ::Val) = ex
        dagwrap(ex::Expr, ::Val{:call}) = :(Dagger.delayed($(ex.args[1]))($(dagwrap.(ex.args[2:end])...)))
        new_rhss = dagwrap.(conv.(rhss))
        delayed_exprs = build_expr(:block, [:($(Symbol(computevars[i])) = Dagger.delayed(identity)($(new_rhss[i]))) for i in axes(computevars,1)])
        # TODO: treereduce?
        reduce_expr = quote
            $(Symbol(reducevar)) = collect(Dagger.delayed(vcat)($(computevars...)))
        end
        ip_let_expr.args[2] = quote
			@sync begin
	            $delayed_exprs
	            $reduce_expr
	            $(ip_let_expr.args[2])
			end
        end
	end

	if fillzeros && outputidxs === nothing
        ip_let_expr = quote
			$fill_array_with_zero!($X)
			$ip_let_expr
		end
    end


    tuple_sys_expr = build_expr(:tuple, [conv(rhs) for rhs ∈ rhss])

    if rhss isa Matrix
        arr_sys_expr = build_expr(:vcat, [build_expr(:row,[conv(rhs) for rhs ∈ rhss[i,:]]) for i in 1:size(rhss,1)])
    elseif typeof(rhss) <: Array && !(typeof(rhss) <: Vector)
        vector_form = build_expr(:vect, [conv(rhs) for rhs ∈ rhss])
        arr_sys_expr = :(reshape($vector_form,$(size(rhss)...)))
    elseif rhss isa SparseMatrixCSC
        vector_form = build_expr(:vect, [conv(rhs) for rhs ∈ nonzeros(rhss)])
        arr_sys_expr = :(SparseMatrixCSC{eltype($(first(argnames))),Int}($(size(rhss)...), $(rhss.colptr), $(rhss.rowval), $vector_form))
    else # Vector
        arr_sys_expr = build_expr(:vect, [conv(rhs) for rhs ∈ rhss])
    end

	xname = gensym(:MTK)

	arr_sys_expr = (typeof(rhss) <: Vector || typeof(rhss) <: Matrix) && !(eltype(rhss) <: AbstractArray) ? quote
		if $force_SA || typeof($(fargs.args[1])) <: Union{ModelingToolkit.StaticArrays.SArray,ModelingToolkit.LabelledArrays.SLArray}
			$xname = ModelingToolkit.StaticArrays.@SArray $arr_sys_expr
			if $convert_oop && !(typeof($(fargs.args[1])) <: Number) && $(typeof(rhss) <: Vector) # Only try converting if it should match `u`
				return similar_type($(fargs.args[1]),eltype($xname))($xname)
			else
				return $xname
			end
		else
			$xname = $arr_sys_expr
			if $convert_oop && $(typeof(rhss) <: Vector)
				if !(typeof($(fargs.args[1])) <: Array) && !(typeof($(fargs.args[1])) <: Number) && eltype($(fargs.args[1])) <: eltype($xname)
					# Last condition: avoid known error because this doesn't change eltypes!
					return convert(typeof($(fargs.args[1])),$xname)
				elseif typeof($(fargs.args[1])) <: ModelingToolkit.LabelledArrays.LArray
					# LArray just needs to add the names back!
					return ModelingToolkit.LabelledArrays.LArray{ModelingToolkit.LabelledArrays.symnames(typeof($(fargs.args[1])))}($xname)
				else
					return $xname
				end
			else
				return $xname
			end
		end
	end : arr_sys_expr

    let_expr = Expr(:let, var_eqs, tuple_sys_expr)
    arr_let_expr = Expr(:let, var_eqs, arr_sys_expr)
    bounds_block = checkbounds ? let_expr : :(@inbounds begin $let_expr end)
    oop_bounds_block = checkbounds ? arr_let_expr : :(@inbounds begin $arr_let_expr end)
    ip_bounds_block = checkbounds ? ip_let_expr : :(@inbounds begin $ip_let_expr end)

    oop_ex = headerfun(oop_bounds_block, fargs, false)
    iip_ex = headerfun(ip_bounds_block, fargs, true; X=X)

    if !linenumbers
        oop_ex = striplines(oop_ex)
        iip_ex = striplines(iip_ex)
    end

    if expression == Val{true}
        return ModelingToolkit.inject_registered_module_functions(oop_ex), ModelingToolkit.inject_registered_module_functions(iip_ex)
    else
        return _build_and_inject_function(@__MODULE__, oop_ex), _build_and_inject_function(@__MODULE__, iip_ex)
    end
end

vars_to_pairs(args) = vars_to_pairs(args[1],args[2])
function vars_to_pairs(name,vs::AbstractArray)
    vs_names = [term_to_symbol(value(u)) for u ∈ vs]
	exs = [:($name[$i]) for (i, u) ∈ enumerate(vs)]
	vs_names,exs
end

function term_to_symbol(t::Term)
    if operation(t) isa Sym
        s = nameof(operation(t))
    else
        error("really?")
    end
end

term_to_symbol(s::Sym) = nameof(s)

function vars_to_pairs(name,vs)
    [term_to_symbol(value(vs))], [name]
end

function rm_calls_with_iv(expr)
    Rewriters.Prewalk(Rewriters.Chain([@rule((~f::(x->x isa Sym))(~t::(x->x isa Sym)) => Sym{symtype((~f)(~t))}((term_to_symbol(~f))))]))(value(expr))
end

get_varnumber(varop::Operation,vars::Vector{Operation}) =  findfirst(x->isequal(x,varop),vars)
get_varnumber(varop::Operation,vars::Vector{<:Variable})  =  findfirst(x->isequal(x,varop.op),vars)

function numbered_expr(O::Operation,args...;varordering = args[1],offset = 0,
                       lhsname=gensym("du"),rhsnames=[gensym("MTK") for i in 1:length(args)])
  if isa(O.op, ModelingToolkit.Variable)
	for j in 1:length(args)
		i = get_varnumber(O,args[j])
		if i !== nothing
			return :($(rhsnames[j])[$(i+offset)])
		end
	end
  end
  return Expr(:call, Symbol(O.op),
         [numbered_expr(x,args...;offset=offset,lhsname=lhsname,
                        rhsnames=rhsnames,varordering=varordering) for x in O.args]...)
end

function numbered_expr(de::ModelingToolkit.Equation,args...;varordering = args[1],
                       lhsname=gensym("du"),rhsnames=[gensym("MTK") for i in 1:length(args)],offset=0)
    i = findfirst(x->isequal(x isa Sym ? term_to_symbol(x) : term_to_symbol(x.op),term_to_symbol(var_from_nested_derivative(de.lhs)[1])),varordering)
    :($lhsname[$(i+offset)] = $(numbered_expr(de.rhs,args...;offset=offset,
											  varordering = varordering,
											  lhsname = lhsname,
											  rhsnames = rhsnames)))
end
numbered_expr(c::ModelingToolkit.Constant,args...;kwargs...) = c.value

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
		eval(:((du::Array{Float64},u::Array{Float64},p::Array{Float64},t::Float64) -> ccall(("diffeqf", $libpath), Cvoid, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64), du, u, p, t)))
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
    real[] $fname(real $iv,real[] $(rhsnames[1]),real[] $(rhsnames[2]),real[] x_r,int[] x_i) {
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
