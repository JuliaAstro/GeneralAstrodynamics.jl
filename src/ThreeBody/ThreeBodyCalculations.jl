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
nondimensionalize(vec::T, scalar) where T<:AbstractVector = @. upreferred(vec / scalar)
nondimensionalize(vec::T, length_scalar, time_scalar) where T<:AbstractVector = nondimensionalize(vec, length_scalar / time_scalar)
nondimensionalize(time::T, scalar) where T<:Number = time / scalar

"""
Returns `vec` re-dimensionalized (`vec` multiplied by `scalar`).
If two scalar arguments are provided, then the first is assumed
to be the `length_scalar`, and the second is assumed to be the
`time_scalar`.
"""
redimensionalize(vec::T, scalar) where T<:AbstractVector = @. upreferred(vec * scalar)
redimensionalize(vec::T, length_scalar, time_scalar) where T<:AbstractVector = redimensionalize(vec, length_scalar / time_scalar)
redimensionalize(time::T, scalar) where T<:Number = time * scalar

"""
Returns the position w.r.t. body 1 (or 2).
"""
position(r, xᵢ=0) = √( (r[1]-xᵢ)^2 + r[2]^2 + r[3]^2 )

"""
Returns the potential energy `U`.
"""
potential_energy(r, μ, x₁, x₂) = (r[1]^2 + r[2]^2) + (2(1-μ)/position(r,x₁)) + (2μ/position(r,x₂))

"""
Returns the Jacobi Constant `C`.
"""
jacobi_constant(r, v, μ, x₁, x₂) = 2*potential_energy(r, μ, x₁, x₂) - (v⋅v)

"""
Returns the position and velocity vectors in the inertial reference frame.
"""
function inertial(vecₛ, t, ω=1.0)

    θ = ω*t
    ᴵTₛ = [
        cos(θ) sin(θ) 0
       -sin(θ) cos(θ) 0
        0      0      1
    ]

    return  ᴵTₛ * vecₛ

end

"""
Returns the position and velocity vectors in the synodic (rotating) reference frame.
"""
synodic(rᵢ, vᵢ, a, Tₛ) =  nondimensionalize(rᵢ, a), nondimensionalize(vᵢ, a, Tₛ)