"""
A superpackage for handling common astrodynamics 
problems. See the __Extended Help__ section
for more information!

# Extended Help

## License
$(LICENSE)

## Imports
$(IMPORTS)

## Exports
$(EXPORTS)
"""
module GeneralAstrodynamics

using Reexport 
using DocStringExtensions

include("AstrodynamicalFrames/AstrodynamicalFrames.jl")
@reexport using .AstrodynamicalFrames

include("AstrodynamicalStates/AstrodynamicalStates.jl")
@reexport using .AstrodynamicalStates

include("AstrodynamicalCalculations/AstrodynamicalCalculations.jl")
@reexport using .AstrodynamicalCalculations

end # module
