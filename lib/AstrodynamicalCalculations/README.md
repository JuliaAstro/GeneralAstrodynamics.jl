[![Tests](https://github.com/cadojo/AstrodynamicalCalculations.jl/workflows/Tests/badge.svg)](https://github.com/cadojo/AstrodynamicalCalculations.jl/actions?query=workflow%3ATests)
[![Docs](https://github.com/cadojo/AstrodynamicalCalculations.jl/workflows/Documentation/badge.svg)](https://cadojo.github.io/AstrodynamicalCalculations.jl)

# `AstrodynamicalCalculations.jl`

_Common calculations relating to orbital mechanics and astrodynamics._

## Installation

```julia
pkg> add AstrodynamicalCalculations
```

## Overview

Almost all of the equations used in your (my) first-year graduate astrodynamics coursework,
in Julia! This package includes equations which are valid within the Restricted Two Body
Problem (R2BP), and the Circular Restricted Three Body Problem (CR3BP). Each function's
docstring clearly labels the models in which the underlying equations are valid. Iterative
solvers, including Kepler's prediction algorithm, and the universal Lambert solver, are
also provided.

Two newer Lambert solver implementations (Lancaster / Blanchard, and
Oldenhuis) are also provided, but they are currently not working. Something seems to have
not transferred well while porting the code from MATLAB to Julia. I hope to work on this
more in 2024.

## Getting Stated

Please refer to the [documentation](https://cadojo.github.io/AstrodynamicalCalculations.jl)
for more detailed instructions, and usage examples.

```julia
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
