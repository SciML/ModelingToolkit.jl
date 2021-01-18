using Test
using Latexify
using ModelingToolkit

### Tips for generating latex tests:
### Latexify has an unexported macro:
###
### Latexify.@generate_test
###
### which generates a test using a given latexify function.
### For example:
###
### Latexify.@generate_test latexify([1, 2, 3], [4, 5, 6]; env=:mdtable)
###
### This puts a ready-made test in your clipboard which you can paste into the
### test file.
###
### Just be sure to remove all such macros before you commit a change since it
### will cause issues with Travis.

@parameters t σ ρ β
@variables x(t) y(t) z(t)
D = Differential(t)

eqs = [D(x) ~ σ*(y-x)*D(x-y)/D(z),
       0 ~ σ*x*(ρ-z)/10-y,
       D(z) ~ x*y^(2//3) - β*z]


# Latexify.@generate_test latexify(eqs)
@test latexify(eqs) == replace(
raw"\begin{align}
\frac{dx(t)}{dt} =& \sigma \frac{d\left(x\left( t \right) -1 \cdot y\left( t \right)\right)}{dt} \left( y\left( t \right) -1 x\left( t \right) \right) \left( \frac{dz(t)}{dt} \right)^{-1} \\
0 =& -1 y\left( t \right) + 0.1 \sigma x\left( t \right) \left( \rho -1 z\left( t \right) \right) \\
\frac{dz(t)}{dt} =& x\left( t \right) \left( y\left( t \right) \right)^{\frac{2}{3}} -1 \beta z\left( t \right)
\end{align}
", "\r\n"=>"\n")


@variables u[1:3](t)
@parameters p[1:3]
eqs = [D(u[1]) ~ p[3]*(u[2]-u[1]),
       0 ~ p[2]*p[3]*u[1]*(p[1]-u[1])/10-u[2],
       D(u[3]) ~ u[1]*u[2]^(2//3) - p[3]*u[3]]

@test latexify(eqs) == replace(
raw"\begin{align}
\frac{du{_1}(t)}{dt} =& p{_3} \left( \mathrm{u{_2}}\left( t \right) -1 \mathrm{u{_1}}\left( t \right) \right) \\
0 =& -1 \mathrm{u{_2}}\left( t \right) + 0.1 p{_2} p{_3} \mathrm{u{_1}}\left( t \right) \left( p{_1} -1 \mathrm{u{_1}}\left( t \right) \right) \\
\frac{du{_3}(t)}{dt} =& \mathrm{u{_1}}\left( t \right) \left( \mathrm{u{_2}}\left( t \right) \right)^{\frac{2}{3}} -1 p{_3} \mathrm{u{_3}}\left( t \right)
\end{align}
", "\r\n"=>"\n")

eqs = [D(u[1]) ~ p[3]*(u[2]-u[1]),
       D(u[2]) ~ p[2]*p[3]*u[1]*(p[1]-u[1])/10-u[2],
       D(u[3]) ~ u[1]*u[2]^(2//3) - p[3]*u[3]]

@test latexify(eqs) == replace(
raw"\begin{align}
\frac{du{_1}(t)}{dt} =& p{_3} \left( \mathrm{u{_2}}\left( t \right) -1 \mathrm{u{_1}}\left( t \right) \right) \\
\frac{du{_2}(t)}{dt} =& -1 \mathrm{u{_2}}\left( t \right) + 0.1 p{_2} p{_3} \mathrm{u{_1}}\left( t \right) \left( p{_1} -1 \mathrm{u{_1}}\left( t \right) \right) \\
\frac{du{_3}(t)}{dt} =& \mathrm{u{_1}}\left( t \right) \left( \mathrm{u{_2}}\left( t \right) \right)^{\frac{2}{3}} -1 p{_3} \mathrm{u{_3}}\left( t \right)
\end{align}
", "\r\n"=>"\n")


@parameters t
@variables x(t)
D = Differential(t)
eqs = [D(x) ~ (1+cos(t))/(1+2*x)]

@test latexify(eqs) == replace(
raw"\begin{align}
\frac{dx(t)}{dt} =& \left( 1 + \cos\left( t \right) \right) \left( 1 + 2 x\left( t \right) \right)^{-1}
\end{align}
", "\r\n"=>"\n")
