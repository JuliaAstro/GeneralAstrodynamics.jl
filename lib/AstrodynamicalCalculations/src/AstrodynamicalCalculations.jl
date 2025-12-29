"""
Common calculations within orbital mechanics and astrodynamics.

# Extended Help

## Documentation

See the project documentation for more information: https://JuliaAstro.org/AstrodynamicalCalculations.jl.

## License

$(LICENSE)
"""
module AstrodynamicalCalculations

using Reexport: @reexport
using DocStringExtensions: EXPORTS, IMPORTS, LICENSE

include("R2BPCalculations.jl")
@reexport using .R2BPCalculations

include("CR3BPCalculations.jl")
@reexport using .CR3BPCalculations

end # module AstrodynamicalCalculations
