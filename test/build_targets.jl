using ModelingToolkit, Test
@parameters t a
@variables x(t) y(t)
@derivatives D'~t
eqs = [D(x) ~ a*x - x*y,
       D(y) ~ -3y + x*y]
@test ModelingToolkit.build_function(eqs,[x,y],[a],t,target = ModelingToolkit.StanTarget()) ==
    """
    real[] diffeqf(real t,real[] internal_var___u,real[] internal_var___p,real[] x_r,int[] x_i) {
      real internal_var___du[2];
      internal_var___du[1] = internal_var___p[1] * internal_var___u[1] - internal_var___u[1] * internal_var___u[2];
      internal_var___du[2] = -3 * internal_var___u[2] + internal_var___u[1] * internal_var___u[2];
      return internal_var___du;
    }
    """

@test ModelingToolkit.build_function(eqs,[x,y],[a],t,target = ModelingToolkit.CTarget()) ==
  """
  void diffeqf(double* internal_var___du, double* internal_var___u, double* internal_var___p, t) {
    internal_var___du[1] = internal_var___p[1] * internal_var___u[1] - internal_var___u[1] * internal_var___u[2];
    internal_var___du[2] = -3 * internal_var___u[2] + internal_var___u[1] * internal_var___u[2];
  }
  """

@test ModelingToolkit.build_function(eqs,[x,y],[a],t,target = ModelingToolkit.MATLABTarget()) ==
  """
  diffeqf = @(t,internal_var___u) [internal_var___p(1) * internal_var___u(1) - internal_var___u(1) * internal_var___u(2); -3 * internal_var___u(2) + internal_var___u(1) * internal_var___u(2)];"""

sys = ODESystem(eqs,t,[x,y],[a])

@test ModelingToolkit.build_function(eqs,[x,y],[a],t,target = ModelingToolkit.CTarget()) ==
      ModelingToolkit.build_function(sys.eqs,[x,y],[a],t,target = ModelingToolkit.CTarget())

@test ModelingToolkit.build_function(eqs,[x,y],[a],t,target = ModelingToolkit.StanTarget()) ==
      ModelingToolkit.build_function(sys.eqs,[x,y],[a],t,target = ModelingToolkit.StanTarget())

@test ModelingToolkit.build_function(eqs,[x,y],[a],t,target = ModelingToolkit.MATLABTarget()) ==
      ModelingToolkit.build_function(sys.eqs,[x,y],[a],t,target = ModelingToolkit.MATLABTarget())
