#
# Check Halo.jl
#

using UnitfulAstrodynamics

μ = nondimensionalize(Earth.μ, Moon.μ)
halo(μ; max_iter=100)



