# [Components and Connectors](@id mtkmodel_connector)

## MTK Model

MTK represents components and connectors with `Model`.

```@docs
ModelingToolkit.Model
```

## Components

Components are models from various domains. These models contain states and their
equations.

### [Defining components with `@mtkmodel`](@id mtkmodel)

`@mtkmodel` is a convenience macro to define components. It returns
`ModelingToolkit.Model`, which includes a constructor that returns the ODESystem, a
`structure` dictionary with metadata, and flag `isconnector` which is set to `false`.

## What can an MTK-Model definition have?

`@mtkmodel` definition contains begin blocks of

  - `@components`: for listing sub-components of the system
  - `@equations`: for the list of equations
  - `@extend`: for extending a base system and unpacking its states
  - `@parameters`: for specifying the symbolic parameters
  - `@structural_parameters`: for specifying non-symbolic parameters
  - `@variables`: for specifing the states

Let's explore these in more detail with the following example:

```@example mtkmodel-example
using ModelingToolkit

@mtkmodel ModelA begin
    @parameters begin
        k1
        k2
    end
end

@mtkmodel ModelB begin
    @parameters begin
        p1 = 1.0, [description = "Parameter of ModelB"]
        p2 = 1.0, [description = "Parameter of ModelB"]
    end
end

@mtkmodel ModelC begin
    @structural_parameters begin
        f = sin
    end
    begin
        v_var = 1.0
    end
    @variables begin
        v(t) = v_var
    end
    @extend p1, p2 = model_b = ModelB(; p1)
    @components begin
        model_a = ModelA(; k1)
    end
    @equations begin
        model_a.k1 ~ f(v)
    end
end
```

### `@parameters` and `@variables` begin block

  - Parameters and variables are declared with respective begin blocks.
  - Variables must be functions of an independent variable.
  - Optionally, initial guess and metadata can be specified for these parameters and variables. See `ModelB` in the above example.
  - Along with creating parameters and variables, keyword arguments of same name with default value `nothing` are created.
  - Whenever a parameter or variable has initial value, for example `v(t) = 0.0`, a symbolic variable named `v` with initial value 0.0 and a keyword argument `v`, with default value `nothing` are created. <br> This way, users can optionally pass new value of `v` while creating a component.

```julia
julia > @named model_c = ModelC(; v = 2.0);

julia > ModelingToolkit.getdefault(model_c.v)
2.0
```

### `@structural_parameters` begin block

  - This block is for non symbolic input arguements. These are for inputs that usually are not meant to be part of components; but influence how they are defined. One can list inputs like boolean flags, functions etc... here.
  - Whenever default values are specified, unlike parameters/variables, they are reflected in the keyword argument list.

#### `@extend` begin block

To extend a partial system,

  - List the variables to unpack. If there is a single variable, explicitly specify it as a tuple.
  - Give a name to the base system
  - List the kwargs of the base system that should be listed as kwargs of the main component.
  - Note that in above example, `p1` is promoted as an argument of `ModelC`. Users can set the value of `p1` as

```julia
julia> @named model_c = ModelC(; p1 = 2.0)

```

However, as `p2` isn't listed in the model definition, its initial guess can't
specified while creating an instance of `ModelC`.

### `@components` begin block

  - Declare the subcomponents within `@components` begin block.
  - The arguments in these subcomponents are promoted as keyword arguments as `subcomponent_name__argname` with `nothing` as default value.
  - Whenever components are created with `@named` macro, these can be accessed with `.` operator as `subcomponent_name.argname`
  - In the above example, `k1` of `model_a` can be set in following ways:

```julia
julia> @named model_c1 = ModelC(; model_a.k1 = 1);

```

And as `k2` isn't listed in the sub-component definition of `ModelC`, its default value can't be modified by users.

### `@equations` begin block

  - List all the equations here

#### A begin block

  - Any other Julia operations can be included with dedicated begin blocks.

## Connectors

Connectors are special models that can be used to connect different components together.
MTK provides 3 distinct connectors:

  - `DomainConnector`: A connector which has only one state which is of `Flow` type,
    specified by `[connect = Flow]`.
  - `StreamConnector`: A connector which has atleast one stream variable, specified by
    `[connect = Stream]`. A `StreamConnector` must have exactly one flow variable.
  - `RegularConnector`: Connectors that don't fall under above categories.

### [Defining connectors with `@connector`](@id connector)

`@connector` returns `ModelingToolkit.Model`. It includes a constructor that returns
a connector ODESystem, a `structure` dictionary with metadata, and flag `isconnector`
which is set to `true`.

A simple connector can be defined with syntax similar to following example:

```@example connector
using ModelingToolkit

@connector Pin begin
    v(t) = 0.0, [description = "Voltage"]
    i(t), [connect = Flow]
end
```

Variables (as functions of independent variable) are listed out in the definition. These variables can optionally have initial values and metadata like `description`, `connect` and so on. For more details on setting metadata, check out [Symbolic Metadata](@id symbolic_metadata).

Similar to `@mtkmodel`, `@connector` accepts begin blocks of `@components`, `@equations`, `@extend`, `@parameters`, `@structural_parameters`, `@variables`. These keywords mean the same as described above for `@mtkmodel`.
For example, the following `HydraulicFluid` connector is defined with parameters, variables and equations.

```@example connector
@connector HydraulicFluid begin
    @parameters begin
        ρ = 997
        β = 2.09e9
        μ = 0.0010016
        n = 1
        let_gas = 1
        ρ_gas = 0.0073955
        p_gas = -1000
    end
    @variables begin
        dm(t) = 0.0, [connect = Flow]
    end
    @equations begin
        dm ~ 0
    end
end
```

!!! note
    
    For more examples of usage, checkout [ModelingToolkitStandardLibrary.jl](https://github.com/SciML/ModelingToolkitStandardLibrary.jl/)

## More on `Model.structure`

`structure` stores metadata that describes composition of a model. It includes:

  - `:components`: List of sub-components in the form of [[name, sub_component_name],...].
  - `:extend`: The list of extended states, name given to the base system, and name of the base system.
  - `:structural_parameters`: Dictionary of structural parameters mapped to their default values.
  - `:parameters`: Dictionary of symbolic parameters mapped to their metadata.
  - `:variables`: Dictionary of symbolic variables mapped to their metadata.
  - `:kwargs`: Dictionary of keyword arguments mapped to their default values.
  - `:independent_variable`: Independent variable, which is added while generating the Model.
  - `:equations`: List of equations (represented as strings).

For example, the structure of `ModelC` is:

```julia
julia> ModelC.structure
Dict{Symbol, Any} with 6 entries:
  :components           => [[:model_a, :ModelA]]
  :variables            => Dict{Symbol, Dict{Symbol, Any}}(:v=>Dict(:default=>:v_var))
  :kwargs               => Dict{Symbol, Any}(:f=>:sin, :v=>:v_var, :p1=>nothing, :model_a__k1=>nothing)
  :independent_variable => t
  :extend               => Any[[:p1, :p2], :model_b, :ModelB]
  :equations            => ["model_a.k1 ~ f(v)"]
```
