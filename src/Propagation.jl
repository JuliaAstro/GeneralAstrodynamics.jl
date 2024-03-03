"""
Wrappers around SciML differential equation solvers for fast and convenient 
orbit propagation.

# Extended Help

## Exports
$(EXPORTS)

## Imports
$(IMPORTS)
"""
module Propagation

using DocStringExtensions

@template (
    FUNCTIONS,
    METHODS,
    MACROS,
) = """
    $(SIGNATURES)

    $(DOCSTRING)
    """

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)

                               $(DOCSTRING)
                               """

using AstrodynamicalModels
using OrdinaryDiffEq

"""
Numerically integrate the orbit forward (or backward) in time, and return a new 
`AstrodynamicalOrbit` instance with identical parameters to the provided orbit.
"""
function propagate(orbit::AstrodynamicalModels.AstrodynamicalOrbit, Δt; algorithm=Vern7(), kwargs...)



end


end