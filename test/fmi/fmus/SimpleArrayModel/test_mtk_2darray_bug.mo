model SimpleArrayModel
    "Simple model with a 2D array variable and its derivative"
    
    parameter Integer n = 3 "First dimension size";
    parameter Integer m = 2 "Second dimension size";
    
    Real x[n,m] "Two-dimensional array variable";
    Real dx[n,m] "Derivative of the array variable";
    
equation
    // Define the derivative relationship
    der(x) = dx;
    
    // Simple example dynamics for the derivative
    for i in 1:n loop
        for j in 1:m loop
            dx[i,j] = -x[i,j] + sin(time);
        end for;
    end for;
    
initial equation
    // Initial conditions
    x = zeros(n,m);
    
end SimpleArrayModel;