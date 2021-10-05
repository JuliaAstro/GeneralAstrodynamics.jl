"""
Provides astrodynamical models as `ModelingToolkit.ODESystems`. 
Check out the `ModelingToolkit` docs to learn how to use these 
systems for orbit propagation with `DifferentialEquations`, or
see `GeneralAstrodynamics` for some convenient orbit propagation 
wrappers.

# Extended help

## Exports

$(EXPORTS)

## Imports

$(IMPORTS)
"""
module AstrodynamicalModels

# Export every model!
export R2BP, CR3BP, NBP

# Export every `ODEFunction`!
export R2BPFunction, CR3BPFunction, NBPFunction

using Symbolics
using StaticArrays
using LinearAlgebra
using ModelingToolkit

import Memoize: @memoize

using DocStringExtensions
@template (FUNCTIONS, METHODS, MACROS) =
    """
    $(SIGNATURES)

    $(DOCSTRING)
    """

@template (TYPES, CONSTANTS) =
    """
    $(TYPEDEF)

    $(DOCSTRING)
    """

include("R2BP.jl")
include("CR3BP.jl")
include("NBP.jl")

end # module
