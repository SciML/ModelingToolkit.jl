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
       D(z) ~ x*y - β*z]


@test latexify(eqs) == 
raw"\begin{align}
\frac{dx(t)}{dt} =& \frac{\sigma \cdot \left( \mathrm{y}\left( t \right) - \mathrm{x}\left( t \right) \right) \cdot \frac{d\left(\mathrm{x}\left( t \right) - \mathrm{y}\left( t \right)\right)}{dt}}{\frac{dz(t)}{dt}} \\
0 =& \frac{\sigma \cdot \mathrm{x}\left( t \right) \cdot \left( \rho - \mathrm{z}\left( t \right) \right)}{10} - \mathrm{y}\left( t \right) \\
\frac{dz(t)}{dt} =& \mathrm{x}\left( t \right) \cdot \mathrm{y}\left( t \right) - \beta \cdot \mathrm{z}\left( t \right)
\end{align}
"

