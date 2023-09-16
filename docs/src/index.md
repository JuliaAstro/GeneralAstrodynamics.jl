# `AstrodynamicalSolvers.jl`

_Common solvers within orbital mechanics and astrodynamics._

```@docs
AstrodynamicalSolvers
```


## Installation

```julia
pkg> add AstrodynamicalSolvers
```

## Getting Started

This package contains differential correctors, and helpful wrapper functions, for 
finding periodic orbits within Circular Restricted Three Body Problem dynamics.

```jldoctest usage
julia> using AstrodynamicalSolvers

julia> μ = 0.012150584395829193
0.012150584395829193

julia> u, T = halo(μ, 1) # lyapunov (planar) orbit

([0.8567678285004178, 0.0, 0.0, 0.0, -0.14693135696819282, 0.0], 2.7536820160579087)

julia> 

julia> u, T = halo(μ, 2; amplitude=0.005) # halo (non-planar) orbit
([1.180859455641048, 0.0, -0.006335144846688764, 0.0, -0.15608881601817765, 0.0], 3.415202902714686)
```