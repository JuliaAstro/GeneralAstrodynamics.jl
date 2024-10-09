# `AstrodynamicalCalculations.jl`

_Common calculations within orbital mechanics and astrodynamics._

```@docs
AstrodynamicalCalculations
```


## Installation

```julia
pkg> add AstrodynamicalCalculations
```

## Getting Started

One common task within first-year orbital mechanics courses is converting to, and from,
Keplerian parameters (orbital elements). Kepler's prediction algorithm is another staple!
These simple calculations are provided by the `R2BPCalculations` submodule within 
`AstrodynamicalCalculations`.

!!! tip
    The size of all of the vectors in these calculations are known at compile time. 
    You can _drastically_ improve performance by using `StaticArrays` types, i.e. `SVector`.

```julia
julia> using AstrodynamicalCalculations

julia> r = [0.0, 11681.0, 0.0];     # km

julia> v = [5.134, 4.226, 2.787];   # km/s

julia> μ = 398600.4354360959;       # km³ s⁻²

julia> e, a, i, Ω, ω, ν = cartesian_to_keplerian(r, v, μ) 
(0.723452708202361, 24509.265399338536, 2.6442542356744734, 1.5707963267948966, 4.712449617676915, 1.5707356895026707)

julia> dt = orbital_period(a, μ)
38186.19850882009

julia> rₙ, vₙ = kepler(r..., v..., μ, dt)
(x = 0.0, y = 11681.0, z = 0.0, ẋ = 5.134, ẏ = 4.226, ż = 2.787)
```