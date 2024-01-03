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

### What can an MTK-Model definition have?

`@mtkmodel` definition contains begin blocks of

  - `@components`: for listing sub-components of the system
  - `@equations`: for the list of equations
  - `@extend`: for extending a base system and unpacking its states
  - `@icon` : for embedding the model icon
  - `@parameters`: for specifying the symbolic parameters
  - `@structural_parameters`: for specifying non-symbolic parameters
  - `@variables`: for specifying the states

Let's explore these in more detail with the following example:

```@example mtkmodel-example
using ModelingToolkit

@mtkmodel ModelA begin
    @parameters begin
        k
        k_array[1:2]
    end
end

@mtkmodel ModelB begin
    @parameters begin
        p1 = 1.0, [description = "Parameter of ModelB"]
        p2 = 1.0, [description = "Parameter of ModelB"]
    end
end

@mtkmodel ModelC begin
    @icon "https://github.com/SciML/SciMLDocs/blob/main/docs/src/assets/logo.png"
    @structural_parameters begin
        f = sin
    end
    begin
        v_var = 1.0
    end
    @variables begin
        v(t) = v_var
        v_array(t)[1:2, 1:3]
    end
    @extend ModelB(; p1)
    @components begin
        model_a = ModelA(; k_array)
    end
    @equations begin
        model_a.k ~ f(v)
    end
end
```

#### `@icon`

An icon can be embedded in 3 ways:

  - URI
  - Path to a valid image-file.<br>
    It can be an absolute path. Or, a path relative to an icon directory; which is
    `DEPOT_PATH[1]/mtk_icons` by default and can be changed by setting
    `ENV["MTK_ICONS_DIR"]`.<br>
    Internally, it is saved in the _File URI_ scheme.

```julia
@mtkmodel WithPathtoIcon begin
    @icon "/home/user/.julia/dev/mtk_icons/icon.png"
    # Rest of the model definition
end
```

  - Inlined SVG.

```julia
@mtkmodel WithInlinedSVGIcon begin
    @icon """<svg height="100" width="100">
    <circle cx="50" cy="50" r="40" stroke="green" fill="none" stroke-width="3"/>
    </svg>
    """
    # Rest of the model definition
end
```

#### `@structural_parameters` begin block

  - This block is for non symbolic input arguments. These are for inputs that usually are not meant to be part of components; but influence how they are defined. One can list inputs like boolean flags, functions etc... here.
  - Whenever default values are specified, unlike parameters/variables, they are reflected in the keyword argument list.

#### `@parameters` and `@variables` begin block

  - Parameters and variables are declared with respective begin blocks.
  - Variables must be functions of an independent variable.
  - Optionally, initial guess and metadata can be specified for these parameters and variables. See `ModelB` in the above example.
  - Along with creating parameters and variables, keyword arguments of same name with default value `nothing` are created.
  - Whenever a parameter or variable has initial value, for example `v(t) = 0.0`, a symbolic variable named `v` with initial value 0.0 and a keyword argument `v`, with default value `nothing` are created. <br> This way, users can optionally pass new value of `v` while creating a component.

```julia
julia> @mtkbuild model_c1 = ModelC(; v = 2.0);

julia> ModelingToolkit.getdefault(model_c1.v)
2.0
```

#### `@extend` begin block

  - Partial systems can be extended in a higher system as `@extend PartialSystem(; kwargs)`.
  - Keyword arguments pf partial system in the `@extend` definition are added as the keyword arguments of the base system.
  - Note that in above example, `p1` is promoted as an argument of `ModelC`. Users can set the value of `p1`. However, as `p2` isn't listed in the model definition, its initial guess can't be specified while creating an instance of `ModelC`.

```julia
julia> @mtkbuild model_c2 = ModelC(; p1 = 2.0)

```

#### `@components` begin block

  - Declare the subcomponents within `@components` begin block.
  - The arguments in these subcomponents are promoted as keyword arguments as `subcomponent_name__argname` with `nothing` as default value.
  - Whenever components are created with `@named` macro, these can be accessed with `.` operator as `subcomponent_name.argname`
  - In the above example, as `k` of `model_a` isn't listed while defining the sub-component in `ModelC`, its default value can't be modified by users. While `k_array` can be set as:

```@example mtkmodel-example
using ModelingToolkit: getdefault

@mtkbuild model_c3 = ModelC(; model_a.k_array = [1.0, 2.0])

getdefault(model_c3.model_a.k_array[1])
# 1.0
getdefault(model_c3.model_a.k_array[2])
# 2.0

@mtkbuild model_c4 = ModelC(model_a.k_array = 3.0)

getdefault(model_c4.model_a.k_array[1])
# 3.0
getdefault(model_c4.model_a.k_array[2])
# 3.0
```

#### `@equations` begin block

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

Variables (as functions of independent variable) are listed out in the definition. These variables can optionally have initial values and metadata like `description`, `connect` and so on. For more details on setting metadata, check out [Symbolic Metadata](@ref symbolic_metadata).

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
  - `:parameters`: Dictionary of symbolic parameters mapped to their metadata. For
    parameter arrays, length is added to the metadata as `:size`.
  - `:variables`: Dictionary of symbolic variables mapped to their metadata. For
    variable arrays, length is added to the metadata as `:size`.
  - `:kwargs`: Dictionary of keyword arguments mapped to their default values.
  - `:independent_variable`: Independent variable, which is added while generating the Model.
  - `:equations`: List of equations (represented as strings).

For example, the structure of `ModelC` is:

```julia
julia> ModelC.structure
Dict{Symbol, Any} with 7 entries:
  :components           => [[:model_a, :ModelA]]
  :variables            => Dict{Symbol, Dict{Symbol, Any}}(:v=>Dict(:default=>:v_var), :v_array=>Dict(:size=>(2, 3)))
  :icon                 => URI("https://github.com/SciML/SciMLDocs/blob/main/docs/src/assets/logo.png")
  :kwargs               => Dict{Symbol, Any}(:f=>:sin, :v=>:v_var, :v_array=>nothing, :model_a__k_array=>nothing, :p1=>nothing)
  :independent_variable => t
  :extend               => Any[[:p2, :p1], Symbol("#mtkmodel__anonymous__ModelB"), :ModelB]
  :equations            => ["model_a.k ~ f(v)"]
```

### Using conditional statements

#### Conditional elements of the system

Both `@mtkmodel` and `@connector` support conditionally defining parameters,
variables, equations, and components.

The if-elseif-else statements can be used inside `@equations`, `@parameters`,
`@variables`, `@components`.

```@example branches-in-components
using ModelingToolkit

@mtkmodel C begin end

@mtkmodel BranchInsideTheBlock begin
    @structural_parameters begin
        flag = true
    end
    @parameters begin
        if flag
            a1
        else
            a2
        end
    end
    @components begin
        if flag
            sys1 = C()
        else
            sys2 = C()
        end
    end
end
```

Alternatively, the `@equations`, `@parameters`, `@variables`, `@components` can be
used inside the if-elseif-else statements.

```@example branches-in-components
@mtkmodel BranchOutsideTheBlock begin
    @structural_parameters begin
        flag = true
    end
    if flag
        @parameters begin
            a1
        end
        @components begin
            sys1 = C()
        end
        @equations begin
            a1 ~ 0
        end
    else
        @parameters begin
            a2
        end
        @equations begin
            a2 ~ 0
        end
    end
end
```

The conditional parts are reflected in the `structure`. For `BranchOutsideTheBlock`, the metadata is:

```julia
julia> BranchOutsideTheBlock.structure
Dict{Symbol, Any} with 5 entries:
  :components           => Any[(:if, :flag, [[:sys1, :C]], Any[])]
  :kwargs               => Dict{Symbol, Any}(:flag=>true)
  :independent_variable => t
  :parameters           => Dict{Symbol, Dict{Symbol, Any}}(:a1=>Dict(:condition=>(:if, :flag, Dict{Symbol, Any}(:kwargs => Dict{Any, Any}(:a1 => nothing), :parameters => Any[Dict{Symbol, Dict{Symbol, Any}}(:a1 => Dict())]), Dict{Symbol, Any}(:kwargs => Dict{Any, Any}(:a2 => nothing), :parameters => Any[Dict{Symbol, Dict{Symbol, Any}}(:a2 => Dict())]))
  :equations            => Any[(:if, :flag, ["a1 ~ 0"], ["a2 ~ 0"])]
```

Conditional entries are entered in the format of `(branch, condition, [case when it is true], [case when it is false])`;
where `branch` is either `:if` or `:elseif`.<br>
The `[case when it is false]` is either an empty vector or `nothing` when only if branch is
present; it is a vector or dictionary whenever else branch is present; it is a conditional tuple
whenever elseif branches are present.

For the conditional components and equations these condition tuples are added
directly, while for parameters and variables these are added as `:condition` metadata.

#### Conditional initial guess of symbolic variables

Using ternary operator or if-elseif-else statement, conditional initial guesses can be assigned to parameters and variables.

```@example branches-in-components
@mtkmodel DefaultValues begin
    @structural_parameters begin
        flag = true
    end
    @parameters begin
        p = flag ? 1 : 2
    end
end
```
