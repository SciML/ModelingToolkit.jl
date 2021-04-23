using ModelingToolkit, Test

include("../examples/rc_model.jl")

# I can't get buffer to work with read
# io = IOBuffer() 

open("rc.jl", "w") do io 
    write(io, rc_model)
end

ser = open("rc.jl", "r") do io 
    read(io, AbstractSystem)
end

@test isequal(ser, rc_model)
