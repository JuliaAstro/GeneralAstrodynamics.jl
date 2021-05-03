"""
Provides functions to interact with Cartesian ASCII ephemeris data.

This is BYOE (bring your own ephemeris) module! Check out
JPL Horizons and 
[wrapper scripts](https://github.com/cadojo/Halo-Explorations/blob/main/scripts/ephemeris/fetch_ephemeris.sh)
to learn how to fetch ephemeris files.
"""
module Ephemeris

# Module Exports
export loadascii, Interpolator, interpolator

# Module Dependencies
using StaticArrays
using Interpolations
using CSV, DataFrames
using ..OrbitsBase

# Module source code
include("LoadASCII.jl")

end