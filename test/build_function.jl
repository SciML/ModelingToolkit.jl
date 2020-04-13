using ModelingToolkit, Test
@variables a b c1 c2 c3 d e g

# Multiple argument matrix
h = [a + b + c1 + c2; c3 + d + e + g] # uses the same number of arguments as our application
h_julia(a, b, c, d, e, g) = [a[1] + b[1] + c[1] + c[2]; c[3] + d[1] + e[1] + g[1]]
function h_julia!(out, a, b, c, d, e, g)
    out .= [a[1] + b[1] + c[1] + c[2]; c[3] + d[1] + e[1] + g[1]]
end

h_str = ModelingToolkit.build_function(h, [a], [b], ([c1, c2, c3], [d], [e], [g]))
h_oop = eval(h_str[1])
h_ip! = eval(h_str[2])
inputs = ([1], [2], [3, 4, 5], [6], [7], [8])

@test h_oop(inputs...) == h_julia(inputs...)
out_1 = Array{Int64}(undef, 2)
out_2 = similar(out_1)
h_ip!(out_1, inputs...)
h_julia!(out_2, inputs...)
@test out_1 == out_2

# Multiple input matrix, some unused arguments
h_skip = [a + b + c1; c2 + c3 + g] # skip d, e
h_julia_skip(a, b, c, d, e, g) = [a[1] + b[1] + c[1]; c[2] + c[3] + g[1]]
function h_julia_skip!(out, a, b, c, d, e, g)
    out .= [a[1] + b[1] + c[1]; c[2] + c[3] + g[1]]
end

h_str_skip = ModelingToolkit.build_function(h_skip, [a], [b], ([c1, c2, c3], [], [], [g]))
h_oop_skip = eval(h_str[1])
h_ip!_skip = eval(h_str[2])
inputs_skip = ([1], [2], [3, 4, 5], [], [], [8])

@test h_oop_skip(inputs_skip...) == h_julia_skip(inputs_skip...)
out_1_skip = Array{Int64}(undef, 2)
out_2_skip = similar(out_1_skip)
h_ip!_skip(out_1_skip, inputs_skip...)
h_julia!_skip(out_2_skip, inputs_skip...)
@test out_1_skip == out_2_skip

# Same as above, except test ability to call with non-matrix arguments (i.e., for `nt`)
inputs_skip_2 = ([1], [2], [3, 4, 5], [], (a = 1, b = 2), [8])
@test h_oop_skip_2(inputs_skip_2...) == h_julia_skip_2(inputs_skip_2...)
out_1_skip_2 = Array{Int64}(undef, 2)
out_2_skip_2 = similar(out_1_skip_2)
h_ip!_skip_2(out_1_skip_2, inputs_skip_2...)
h_julia!_skip_2(out_2_skip_2, inputs_skip_2...)
@test out_1_skip_2 == out_2_skip_2

# Multiple input scalar
h_scalar = a + b + c1 + c2 + c3 + d + e + g
h_julia_scalar(a, b, c, d, e, g) = a[1] + c[1]; c[2] + c[3] + d[1] + e[1] + g[1]
h_str_scalar = ModelingToolkit.build_function(h_scalar, [a], [b], ([c1, c2, c3], [d], [e], [g]))
h_oop_scalar = eval(h_str[1])
@test h_oop_scalar(inputs...) == h_julia_scalar(inputs...)
