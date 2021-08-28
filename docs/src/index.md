# GeneralAstrodynamics.jl
_Common astrodynamics calculations in Julia, with units!_

!!! note
    This package is fairly new, and documentation is even newer! 
    Thanks for your patience as we get these docs up and running. 

## Introduction

Welcome! This package contains types, functions, and abstractions
for common astrodynamics models, including the _Restricted Two-body Problem_, 
the _Circular Restricted Three-body Problem_, and the _N-body Problem_, \
as well as some abstractions for interpolating ephemeris data. 
If you're not familiar with those models by name, that's okay. We'll walk
through common concepts in astrodynamics in this documentation, alongside
`GeneralAstrodynamics` usage and syntax!

You might imagine we want to describe how a spacecraft moves in space.
As with all physical systems, reality is really complicated – we can 
make simplifying assumptions to write _models_, or _equations of motion_
which _approximately_ describe how a spacecraft moves in space. 
Some models are better than others in a variety of circumstances. For 
example, if your spacecraft is only $200$ km above the Earth's surface,
then the gravity due to Earth will _far_ outweigh the gravity of any other
celestial body in our universe. On the other hand, if you're 
around halfway between the Earth and the Moon, then a model which only
includes the gravitational affect of _one_ of the two bodies may not be 
accurate enough to be useful. 

Several models are provided by 
[`AstrodynamicalModels.jl`](https://github.com/cadojo/AstrodynamicalModels.jl).
Each are incorporated into this package. **This documentation will soon get a lot better!**

## Restricted Two-body Problem
_A massless spacecraft orbiting __one__ point mass._


## Circular Restricted Three-body Problem
_A massless spacecraft orbiting __two__ point masses which orbit their common center of mass._


## N-body Problem
_A collection of __N__ point masses which are each affected by the gravitational pull of the others._
