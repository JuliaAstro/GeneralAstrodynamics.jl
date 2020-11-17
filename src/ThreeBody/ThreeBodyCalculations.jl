#
# ThreeBodyCalculations.jl
# Calculations for the Circular Restricted
# Three Body Problem.
#

"""
Returns `vec` non-dimensionalized (`vec` divided by `scalar`).
If two scalar arguments are provided, then the first is assumed
to be the `length_scalar`, and the  second is assumed to be the
`time_scalar`.
"""
nondimensionalize(vec, scalar) = @. upreferred(vec / scalar)
nondimensionalize(vec, length_scalar, time_scalar) = nondimensionalize(vec, length_scalar / time_scalar)

"""
Returns the position w.r.t. body 1 (or 2).
"""
position(r, xᵢ=0) = @views @. √( (r[1]-xᵢ)^2 + r[2]^2 + r[3]^2 )

"""
Returns the potential energy `U`.
"""
potential_energy(r, μ, x₁, x₂) = @views @. (r[1]^2 + r[2]^2) + (2(1-μ)/position(r,x₁)) + (2μ/position(r,x₂))

"""
Returns the Jacobi Constant `C`.
"""
jacobi_constant(r, v, μ, x₁, x₂) = potential_energy(r, μ, x₁, x₂) - (v⋅v)