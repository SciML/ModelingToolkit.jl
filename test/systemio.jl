using ModelingToolkit, Test

include("../examples/rc_model.jl")

# I can't get buffer to work with read
# io = IOBuffer() 

open("rc.jl", "w") do io 
    write(io, rc_model)
end

open("rc.jl", "r") do io 
    ser = read(io, ODESystem)
end

@tset ser == rc_model
