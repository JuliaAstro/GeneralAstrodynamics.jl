# Orbits.jl

## Overview

Welcome to the _new_ iteration of `UnitfulAstrodynamics`!
This package contains types, functions, and abstractions
for common astrodynamics problems, including the 
_Restricted Two-body Problem_, the 
_Circular Restricted Three-body Problem_, and the 
_N-body Problem_, as well as some abstractions 
for interpolating ephemeris data.

These new docs are entirely new! The old version was 
completely deleted. Thanks for your patience as we 
get this page up and running!

## Examples!

```julia
using Orbits

using Plots
using Unitful, UnitfulAngles

# Find a numerically periodic halo orbit
halo_orbit, halo_period = halo(SunEarth; Az=200_000u"km", L=2)

# Propagate the orbit in time
trajectory = propagate(halo_orbit, halo_period)

# And plot!
plotpositions(trajectory; title="Try Plots.jl `kwargs`")

# Make a simple Restricted Two-body orbit about Earth
orbit = Orbit(randn(3), randn(3), Earth) # defaults to km, km/s

# Specify units with...
orbit = Orbit(randn(3), randn(3), Earth; lengthunits=u"km", timeunits=u"s")

# Or with...
orbit = Orbit(randn(3) * u"km", randn(3) * u"km/s", Earth)

# Or, specify an orbit with Keplerian elements! (e, a, i, Ω, ω, ν)
orbit = KeplerianOrbit(0.4, 10_000u"km", 0u"rad", 0u"rad", 0u"rad", 0u"rad", Earth)

# Propagate for one orbital period
trajectory = propagate(orbit, period(orbit))

# And plot! 
plotpositions(trajectory)
```

__There are many more features. More documentation to come!__