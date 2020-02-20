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
@derivatives D'~t

eqs = [D(x) ~ σ*(y-x)*D(x-y)/D(z),
       0 ~ σ*x*(ρ-z)/10-y,
       D(z) ~ x*y^(2//3) - β*z]


@test latexify(eqs) ==
raw"\begin{align}
\mathrm{derivative}\left( \mathrm{x}\left( t \right), t \right) =& \frac{\sigma \cdot \left( \mathrm{y}\left( t \right) - \mathrm{x}\left( t \right) \right) \cdot \mathrm{derivative}\left( \mathrm{x}\left( t \right) - \mathrm{y}\left( t \right), t \right)}{\mathrm{derivative}\left( \mathrm{z}\left( t \right), t \right)} \\
0 =& \frac{\sigma \cdot \mathrm{x}\left( t \right) \cdot \left( \rho - \mathrm{z}\left( t \right) \right)}{10} - \mathrm{y}\left( t \right) \\
\mathrm{derivative}\left( \mathrm{z}\left( t \right), t \right) =& \mathrm{x}\left( t \right) \cdot \left( \mathrm{y}\left( t \right) \right)^{\frac{2}{3}} - \beta \cdot \mathrm{z}\left( t \right)
\end{align}
"

@parameters t p[1:3]
@variables u[1:3](t)
@derivatives D'~t

eqs = [D(u[1]) ~ p[3]*(u[2]-u[1]),
       0 ~ p[2]*p[3]*u[1]*(p[1]-u[1])/10-u[2],
       D(u[3]) ~ u[1]*u[2]^(2//3) - p[3]*u[3]]

@test latexify(eqs) ==
raw"\begin{align}
\mathrm{derivative}\left( \mathrm{u_1}\left( t \right), t \right) =& p_3 \cdot \left( \mathrm{u_2}\left( t \right) - \mathrm{u_1}\left( t \right) \right) \\
0 =& \frac{p_2 \cdot p_3 \cdot \mathrm{u_1}\left( t \right) \cdot \left( p_1 - \mathrm{u_1}\left( t \right) \right)}{10} - \mathrm{u_2}\left( t \right) \\
\mathrm{derivative}\left( \mathrm{u_3}\left( t \right), t \right) =& \mathrm{u_1}\left( t \right) \cdot \left( \mathrm{u_2}\left( t \right) \right)^{\frac{2}{3}} - p_3 \cdot \mathrm{u_3}\left( t \right)
\end{align}
"

eqs = [D(u[1]) ~ p[3]*(u[2]-u[1]),
       D(u[2]) ~ p[2]*p[3]*u[1]*(p[1]-u[1])/10-u[2],
       D(u[3]) ~ u[1]*u[2]^(2//3) - p[3]*u[3]]

sys = ODESystem(eqs)

@test latexify(sys.eqs) ==
raw"\begin{align}
\frac{du_{1}}{dt} =& p_3 \cdot \left( \mathrm{u_2}\left( t \right) - \mathrm{u_1}\left( t \right) \right) \\
\frac{du_{2}}{dt} =& \frac{p_2 \cdot p_3 \cdot \mathrm{u_1}\left( t \right) \cdot \left( p_1 - \mathrm{u_1}\left( t \right) \right)}{10} - \mathrm{u_2}\left( t \right) \\
\frac{du_{3}}{dt} =& \mathrm{u_1}\left( t \right) \cdot \left( \mathrm{u_2}\left( t \right) \right)^{\frac{2}{3}} - p_3 \cdot \mathrm{u_3}\left( t \right)
\end{align}
"

@test latexify(sys) ==
raw"$\begin{align}
\frac{du_{1}}{dt} =& p_3 \cdot \left( \mathrm{u_2}\left( t \right) - \mathrm{u_1}\left( t \right) \right) \\
\frac{du_{2}}{dt} =& \frac{p_2 \cdot p_3 \cdot \mathrm{u_1}\left( t \right) \cdot \left( p_1 - \mathrm{u_1}\left( t \right) \right)}{10} - \mathrm{u_2}\left( t \right) \\
\frac{du_{3}}{dt} =& \mathrm{u_1}\left( t \right) \cdot \left( \mathrm{u_2}\left( t \right) \right)^{\frac{2}{3}} - p_3 \cdot \mathrm{u_3}\left( t \right)
\end{align}
$"
