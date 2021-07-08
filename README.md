[![Tests](https://github.com/cadojo/GeneralAstrodynamics.jl/workflows/Tests/badge.svg)](https://github.com/cadojo/GeneralAstrodynamics.jl/actions?query=workflow%3ATests)
[![Docs](https://github.com/cadojo/GeneralAstrodynamics.jl/workflows/Documentation/badge.svg)](https://cadojo.github.io/GeneralAstrodynamics.jl/dev)

# GeneralAstrodynamics.jl
_Common astrodynamics calculations, with units!_

## Features

### Restricted Two-body Problem (R2BP)
* Structures for Cartesian and Keplerian states, and R2BP systems
* Functions which implement common R2BP equations
* Kepler and Lambert solvers
* Orbit propagation and plotting

### Circular Restricted Three-body Problem (CR3BP)
* Structures for dimensioned and normalized Cartesian states, and dimensioned and normalized CR3BP systems
* Functions which implement common CR3BP equations
* Analytical and iterative (numerical) Halo orbit solvers
* Unstable and stable Halo orbit manifold computation
* Orbit propagation and plotting
* Zero-velocity curve computation and plotting

### N-body Problem (NBP)
* This was implemented in a previous package version, and is currently being refactored

## Usage

Check out the [Getting Started](https://cadojo.github.io/GeneralAstrodynamics.jl/dev/) documentation for code examples, and more detail about using this package. 
