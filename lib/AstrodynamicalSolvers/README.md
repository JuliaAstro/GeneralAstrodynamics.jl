# `AstrodynamicalSolvers.jl`

_Common solvers relating to orbital mechanics and astrodynamics._

## Installation

```julia
pkg> add AstrodynamicalSolvers
```

## Overview

This package contains most of the solvers you (I) encountered in your (my) first year of graduate astrodynamics coursework.
Right now, the primary functionality is iterative solvers for periodic orbits within CR3BP dynamics.
## Getting Stated

Please refer to the [documentation](https://cadojo.github.io/AstrodynamicalSolvers.jl)
for more detailed instructions, and usage examples.

```julia
julia> μ = 0.012150584395829193
0.012150584395829193

julia> u, T = halo(μ, 1) # lyapunov (planar) orbit

([0.8567678285004178, 0.0, 0.0, 0.0, -0.14693135696819282, 0.0], 2.7536820160579087)

julia> u, T = halo(μ, 2; amplitude=0.005) # halo (non-planar) orbit
([1.180859455641048, 0.0, -0.006335144846688764, 0.0, -0.15608881601817765, 0.0], 3.415202902714686)
```
