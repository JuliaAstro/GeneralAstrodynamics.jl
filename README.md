[![Tests](https://github.com/cadojo/AstrodynamicalModels.jl/workflows/UnitTests/badge.svg)](https://github.com/cadojo/AstrodynamicalModels.jl/actions?query=workflow%3AUnitTests)
[![Docs](https://github.com/cadojo/AstrodynamicalModels.jl/workflows/Documentation/badge.svg)](https://cadojo.github.io/AstrodynamicalModels.jl)

# AstrodynamicalModels.jl

_An extension of
[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) which provides
common astrodynamics models._

## Installation

Choose one of the two lines below!

```julia
Pkg.add("AstrodynamicalModels")  # in Julia code
```

```julia
pkg> add AstrodynamicalModels    # in Julia's REPL
```

## Currently Implemented

Note – for all non-entry models below, you can optionally append state transition matrix
dynamics.

- Restricted Two-body Problem
- Circular Restricted Three-body Problem
- N-body Problem
- Planar Entry
- Attitude Kinematics & Dynamics

## Future Additions

- Aspherical Restricted Two-body Problem
- Solar radiation pressure dynamics
- _Others? Let me know about, or submit a PR with, your desired astrodynamics
  models!_
