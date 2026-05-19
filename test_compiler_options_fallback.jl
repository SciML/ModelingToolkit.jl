# Test that module-level compiler options are inherited by functions eval'd in those modules,
# including RuntimeGeneratedFunctions.
#
# Run with: julia --project test_compiler_options_fallback.jl

using RuntimeGeneratedFunctions
using InteractiveUtils: code_native, code_typed

# Create modules with different compiler options
module EvalOpt0
    Base.Experimental.@compiler_options optimize=0
    using RuntimeGeneratedFunctions
    RuntimeGeneratedFunctions.init(@__MODULE__)
end

module EvalOpt1
    Base.Experimental.@compiler_options optimize=1
    using RuntimeGeneratedFunctions
    RuntimeGeneratedFunctions.init(@__MODULE__)
end

module EvalOpt0NoInfer
    Base.Experimental.@compiler_options optimize=0 infer=false
    using RuntimeGeneratedFunctions
    RuntimeGeneratedFunctions.init(@__MODULE__)
end

module EvalDefault
    using RuntimeGeneratedFunctions
    RuntimeGeneratedFunctions.init(@__MODULE__)
end

# A non-trivial expression to compile
test_expr = :(function (x)
    y = sin(x) + cos(x)
    return y * y + exp(x)
end)

# Test eval path
println("=== Testing eval path ===")
f_default = EvalDefault.eval(test_expr)
f_opt0 = EvalOpt0.eval(test_expr)
f_opt1 = EvalOpt1.eval(test_expr)
f_opt0_noinfer = EvalOpt0NoInfer.eval(test_expr)

x = 1.5
println("Default:         f($x) = $(f_default(x))")
println("Opt0:            f($x) = $(f_opt0(x))")
println("Opt1:            f($x) = $(f_opt1(x))")
println("Opt0+no infer:   f($x) = $(f_opt0_noinfer(x))")
println("Results match: ", f_default(x) == f_opt0(x) == f_opt1(x) == f_opt0_noinfer(x))

# Test RGF path
println("\n=== Testing RuntimeGeneratedFunction path ===")
rgf_default = RuntimeGeneratedFunction(EvalDefault, EvalDefault, test_expr)
rgf_opt0 = RuntimeGeneratedFunction(EvalOpt0, EvalOpt0, test_expr)
rgf_opt1 = RuntimeGeneratedFunction(EvalOpt1, EvalOpt1, test_expr)
rgf_opt0_noinfer = RuntimeGeneratedFunction(EvalOpt0NoInfer, EvalOpt0NoInfer, test_expr)

println("RGF Default:         f($x) = $(rgf_default(x))")
println("RGF Opt0:            f($x) = $(rgf_opt0(x))")
println("RGF Opt1:            f($x) = $(rgf_opt1(x))")
println("RGF Opt0+no infer:   f($x) = $(rgf_opt0_noinfer(x))")
println("Results match: ", rgf_default(x) == rgf_opt0(x) == rgf_opt1(x) == rgf_opt0_noinfer(x))

# Check that the compiler options actually take effect by inspecting code_native.
# At opt level 0, the native code should contain more instructions (no inlining/folding).
println("\n=== Checking native code differences ===")
native_default = sprint(code_native, rgf_default, Tuple{Float64})
native_opt0 = sprint(code_native, rgf_opt0, Tuple{Float64})
native_opt1 = sprint(code_native, rgf_opt1, Tuple{Float64})

# Count call instructions — opt0 should have more calls (nothing inlined)
count_calls(s) = count(r"callq?\s"i, s)
println("Call instructions - Default: $(count_calls(native_default)), Opt0: $(count_calls(native_opt0)), Opt1: $(count_calls(native_opt1))")
println("Raw instruction count - Default: $(count('\n', native_default)), Opt0: $(count('\n', native_opt0)), Opt1: $(count('\n', native_opt1))")
if count_calls(native_opt0) > count_calls(native_opt1)
    println("Opt0 has more calls than Opt1 -> module-level options ARE taking effect")
elseif count_calls(native_opt0) == count_calls(native_opt1)
    println("Same call count — check raw instruction counts above for differences")
else
    println("WARNING: Opt0 has fewer calls than Opt1 — unexpected")
end

# Check inference works/doesn't work
println("\n=== Checking infer=false effect ===")
ci_default = only(code_typed(rgf_default, Tuple{Float64}))
ci_noinfer = only(code_typed(rgf_opt0_noinfer, Tuple{Float64}))
println("Default return type:       $(ci_default[2])")
println("Opt0+no infer return type: $(ci_noinfer[2])")
if ci_noinfer[2] === Any
    println("infer=false IS taking effect (return type is Any)")
else
    println("infer=false may NOT be taking effect (return type is concrete)")
end

println("\nDone.")
