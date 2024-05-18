"""
Provides astrodynamical solvers, including Lyapunov and halo orbit correctors.

# Extended help

## License
$(LICENSE)

## Exports
$(EXPORTS)

## Imports
$(IMPORTS)
"""
module AstrodynamicalSolvers

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


include("CR3BSolvers.jl")
include("Propagation.jl")

@reexport using .CR3BSolvers
@reexport using .Propagation

end