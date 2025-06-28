"""
A superpackage for handling common astrodynamics problems.
See the **Extended Help** section for more information!

# Extended Help

## License
$(LICENSE)

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
@reexport using SPICE: furnsh
@reexport using EphemerisSources

end # module
