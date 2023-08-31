"""
Common calculations within orbital mechanics and astrodynamics.

# Extended Help

## Imports

$(IMPORTS)

## Exports

$(EXPORTS)

"""
module AstrodynamicalCalculations

using Reexport
using DocStringExtensions

@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(SIGNATURES)
                                         $(DOCSTRING)
                                         """

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)
                               $(DOCSTRING)
                               """

include("R2BPCalculations.jl")
@reexport using .R2BPCalculations

include("CR3BPCalculations.jl")
@reexport using .CR3BPCalculations

end # module AstrodynamicalCalculations
