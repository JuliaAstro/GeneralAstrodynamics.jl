"""
Provides astrodynamical models as `ModelingToolkit.ODESystems`.
Check out the `ModelingToolkit` docs to learn how to use these
systems for orbit propagation with `DifferentialEquations`, or
see `GeneralAstrodynamics` for some convenient orbit propagation
wrappers.

# Extended help

## License
$(LICENSE)

## Exports
$(EXPORTS)

## Imports
$(IMPORTS)
"""
module AstrodynamicalModels

# Export every model!
export R2BP, CR3BP, NBP, PlanarEntry, Attitude

# Export every `ODEFunction`!
export R2BPFunction, CR3BPFunction, NBPFunction, PlanarEntryFunction, AttitudeFunction

using Symbolics
using LinearAlgebra
using ModelingToolkit

using DocStringExtensions
@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(SIGNATURES)

                                         $(DOCSTRING)
                                         """

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)

                               $(DOCSTRING)
                               """

include("R2BP.jl")
include("CR3BP.jl")
include("NBP.jl")
include("Entry.jl")
include("Attitude.jl")

end # module
