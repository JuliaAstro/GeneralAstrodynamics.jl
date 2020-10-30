# `NBody` Data Structures

As with [`TwoBody` Data Structures](@ref), the `NBody` submodule includes data structures for storing multibody orbital states. The `Body` structure holds position, velocity, and mass information for a single body. A `MultibodySystem` contains an array of `Body` structures, and is used to completely describe a multibody orbital state.

```@docs
Body
MultibodySystem
```