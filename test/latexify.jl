using Test
using Latexify
using ModelingToolkit
using ReferenceTests
using ModelingToolkit: t_nounits as t, D_nounits as D

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

@parameters σ ρ β
@variables x(t) y(t) z(t)

eqs = [D(x) ~ σ * (y - x) * D(x - y) / D(z),
    0 ~ σ * x * (ρ - z) / 10 - y,
    D(z) ~ x * y^(2 // 3) - β * z]

# Latexify.@generate_test latexify(eqs)
@test_reference "latexify/10.tex" latexify(eqs)

@variables u(t)[1:3]
@parameters p[1:3]
eqs = [D(u[1]) ~ p[3] * (u[2] - u[1]),
    0 ~ p[2] * p[3] * u[1] * (p[1] - u[1]) / 10 - u[2],
    D(u[3]) ~ u[1] * u[2]^(2 // 3) - p[3] * u[3]]

@test_reference "latexify/20.tex" latexify(eqs)

eqs = [D(u[1]) ~ p[3] * (u[2] - u[1]),
    D(u[2]) ~ p[2] * p[3] * u[1] * (p[1] - u[1]) / 10 - u[2],
    D(u[3]) ~ u[1] * u[2]^(2 // 3) - p[3] * u[3]]

@test_reference "latexify/30.tex" latexify(eqs)
@variables x(t)
eqs = [D(x) ~ (1 + cos(t)) / (1 + 2 * x)]

@test_reference "latexify/40.tex" latexify(eqs)
