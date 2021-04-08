#
# Run all unit tests in UnitulAstrodynamics.jl
#

include("TwoBody/test_twobody.jl")
include("ThreeBody/test_threebody.jl")
include("NBody/test_nbody.jl")
include("Propagators/test_propagators.jl")
include("OrbitPlots/test_plots.jl")
