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

"""
Load several small SPICE kernels to support simple operations. The kernels
loaded into memory sum to less than 11MB.

# Extended Help

Calls `SPICE.furnsh` on the following kernels.

- `de432s.bsp` (approximately 10MB)
- `latest_leapseconds_lsk.tls` (approximately 5KB)
- `gm_de440.tpc` (approximately 12KB)
- `pck00011.tpc` (approximately 130KB)
"""
function load_generic_spice_kernels!()
    return furnsh(
        de432s(),                   # position and velocity data for nearby planets
        latest_leapseconds_lsk(),   # timekeeping, parsing epochs
        gm_de440(),                 # mass parameters for major solar system bodies
        pck00011(),                 # physical properties of major solar system bodies
    )
end

@reexport using AstrodynamicalCalculations
@reexport using AstrodynamicalModels
@reexport using AstrodynamicalSolvers
@reexport using SPICE: furnsh
@reexport using SPICEKernels
@reexport using SPICEBodies

end # module
