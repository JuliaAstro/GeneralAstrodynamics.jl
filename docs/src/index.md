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

julia> r = [0.0, 11681.0, 0.0];     # km

julia> v = [5.134, 4.226, 2.787];   # km/s

julia> μ = 398600.4354360959;       # km³ s⁻²

julia> e, a, i, Ω, ω, ν = cartesian_to_keplerian(r, v, μ) 
(0.723452708202361, 24509.265399338536, 2.6442542356744734, 1.5707963267948966, 4.712449617676915, 1.5707356895026707)

julia> T = orbital_period(a, μ) 
38186.19850882009

julia> rₙ, vₙ = kepler(r, v, μ, 2.5T)
([36872.96963754877, -2574.241491333036, 20016.549742861007], [-0.3141726028666592, -1.6044679459972122, -0.17054909314167882])
```