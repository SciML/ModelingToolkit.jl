# ModelingToolkit.jl Documentation


ModelingToolkit.jl is an intermediate representation (IR) of computational graphs
for scientific computing problems. Its purpose is to be a common target for
modeling DSLs in order to allow for a common platform for model inspection and
transformation. It uses a tagged variable IR in order to allow specification of
complex models and allow for transformations of models. It has ways to plug into
its function registration and derivative system so that way it can interact
nicely with user-defined routines. Together, this is an abstract form of a
scientific model that is easy for humans to generate but also easy for programs
to manipulate.
