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

include("R2BPCalculations.jl")
@reexport using .R2BPCalculations

include("CR3BPCalculations.jl")
@reexport using .CR3BPCalculations

end # module AstrodynamicalCalculations
