# AstrodynamicalModels.jl
_Common models within astrodynamics!_

## Overview

This package extends `ModelingToolkit` to represent
common astrodynamical models. All available models
are shown on the [Docstrings](docstrings.md) page.
Consult the **Models** pages for more detail about
each model in this package!

## Usage

If you're familiar with [`ModelingToolkit.jl`](https://mtk.sciml.ai/dev/),
then you'll be able to use this package! Some 
`AstrodynamicalModels`-specific usage instructions are provided 
here. Please don't be shy about making [Discourse](https://discourse.julialang.org)
posts, or filing [issues](https://github.com/cadojo/AstrodynamicalModels.jl) on
GitHub!

### Installation & Setup

This package can be installed just like any other 
[registered](https://juliahub.com) Julia package.

```julia
# To install wherever Julia code runs...
import Pkg
Pkg.add("AstrodynamicalModels") # or ]add AstrodynamicalModels in Julia's REPL
```

To load the package, simply enter `using AstrodynamicalModels`.

```@repl main
using AstrodynamicalModels
```

### Retrieving a Model

Each model within this package is implemented with a function –
each function returns some `AbstractSystem` from `ModelingToolkit.jl`.
Typically, this will be an `ODESystem`. If you're worried about
overhead from calling each function every time you need a particular
model, don't! Each function is implemented with 
[`@memoize`](https://github.com/JuliaCollections/Memoize.jl), so all 
results are cached the first time you call a model's function 
_with a particular function signature_. 

```@repl main
R2BPModel = R2BP() # Restricted Two-body Problem dynamics

CR3BPModel = CR3BP() # Circular Restricted Three-body Problem dynamics

CR3BPModelWithSTM = CR3BP(; stm=true) # Optionally include state transition matrix dynamics
```

### Using a Model

To actually _use_ each model, you probably _also_ want to load 
`ModelingToolkit` (and any other [SciML](https://sciml.ai) 
packages of your choice).

```@repl main
using ModelingToolkit
```

Now you can use any method defined for `ModelingToolkit.AbstractSystem`
instances. Once again, the [`ModelingToolkit` Documentation](https://mtk.sciml.ai)
are the best place to learn how to interact with `AbstractSystem` instances!
Some quick examples are shown below. 

#### Check the Equations of Motion
```@repl main
eqs = equations(R2BP())
```

#### List the States and Parameters
```@repl main
x = states(R2BP())
p = parameters(R2BP())
```

#### Calculate the Jacobian
```@repl main
J = calculate_jacobian(R2BP())
```

#### Generate Code to Replicate the Model
```@repl main
print(build_function(R2BP()))
```

#### Generate Code which Implements the Dynamics
```@repl main
print(R2BPFunction())
```

#### Generate C/C++ and MATLAB Code
```@repl main
print(build_function([eq.rhs for eq in equations(R2BP())], states(R2BP()), parameters(R2BP()); target=Symbolics.CTarget()))
print(build_function([eq.rhs for eq in equations(R2BP())], states(R2BP()), parameters(R2BP()); target=Symbolics.MATLABTarget()))
```

If you're interested in learning a bit about each astrodynamical model, or you'd like
more specific examples which show how to use each model, consult the __Models__ 
pages! 



