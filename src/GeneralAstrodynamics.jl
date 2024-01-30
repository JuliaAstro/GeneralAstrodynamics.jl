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

@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(SIGNATURES)
                                         $(DOCSTRING)
                                         $(METHODLIST)
                                         """

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)
                               $(DOCSTRING)
                               """

@reexport using AstrodynamicalCalculations
@reexport using AstrodynamicalModels
@reexport using AstrodynamicalSolvers
@reexport using SPICE, SPICEKernels, SPICEBodies

end # module
