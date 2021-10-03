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

# NOTE right at the top. DO NOT CHANGE THE ORDER
# OF THE VARIABLES IN THE MODELS BELOW. 
# Downstream users do not have any way to 
# specify individual states when constructing 
# an `ODEProblem` from each `ODESystem`.

# Export every model!
export R2BP, CR3BP

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

end # module
