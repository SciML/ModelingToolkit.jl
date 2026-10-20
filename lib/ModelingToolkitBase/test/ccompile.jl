using ModelingToolkitBase, Test, Libdl
using Symbolics: CTarget
using ModelingToolkitBase: t_nounits as t, D_nounits as D

@parameters a
@variables x y
eqs = [
    D(x) ~ a * x - x * y,
    D(y) ~ -3y + x * y,
]
ccode = build_function([x.rhs for x in eqs], [x, y], [a], t, target = CTarget())
@test ccode isa String
libpath = tempname() * "." * Libdl.dlext
open(`gcc -fPIC -O3 -xc -shared -o $libpath -`, "w") do io
    print(io, ccode)
end
lib = Libdl.dlopen(libpath)
fptr = Libdl.dlsym(lib, :diffeqf)
f(du, u, p, t) = ccall(
    fptr, Cvoid, (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64), du, u, p, t
)
f2 = eval(build_function([x.rhs for x in eqs], [x, y], [a], t)[2])
du = rand(2);
du2 = rand(2);
u = rand(2)
p = rand(1)
_t = rand()
f(du, u, p, _t)
f2(du2, u, p, _t)
@test du ≈ du2
Libdl.dlclose(lib)
