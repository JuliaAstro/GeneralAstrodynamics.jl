#
# ThreeBodyCalculations.jl
# Calculations for the Circular Restricted
# Three Body Problem.
#

"""
Normalizes a CR3BP orbit in the rotating reference frame.
"""
function LinearAlgebra.normalize(r::AbstractVector{<:Real}, v::AbstractVector{<:Real}, t::Real, a::Real, μs::Tuple{<:Real, <:Real}, u"km", u"s")

    rₙ = r ./ a
    Tₛ = period(a, sum(μs))
    vₙ = v ./ (a / Tₛ)
    tₙ = t / Tₛ
    μₙ = min(μs) / sum(μs)

    DU = a * lengthunit
    DT = Tₛ * timeunit

    return rₙ, vₙ, tₙ, μₙ, DU, DT

end

"""
Redimensionalizes a CR3BP orbit in the rotating reference frame.
"""
function redimension(rₙ::AbstractVector{<:Real}, vₙ::AbstractVector{<:Real}, tₙ::Real, μₙ::Real, DU::Number, TU::Number)

    r = rₙ .* DU
    v = vₙ .* DU ./ TU
    t = tₙ * TU

    sum_μs = DU^3 / ((TU / 2π)^2)
    μ₂ = μₙ * sum_μs
    μ₁ = sum_μ - μ₂

    return r, v, t, DU, (μ₁,μ₂)

end

"""
Returns dimensional form of (`Unitful`) scalar posiion.
"""
redimensionalize(rᵢ::R, a::A) where {
        R<:Real, A<:Unitful.Length
    } = redimensionalize_length(rᵢ, a)

"""
Returns dimensional form of (`Unitful`) velocity vector.
"""
redimensionalize(vᵢ::U, a::A, Tₛ::T) where {
        U<:Real, A<:Unitful.Length, T<:Unitful.Time
    } = redimensionalize_velocity(vᵢ, a, Tₛ)

"""
Returns dimensional form of (`Unitful`) time duration.
"""
redimensionalize(t::T1, Tₛ::T2) where {
        T1<:Real, T2<:Unitful.Time
    } = redimensionalize_time(t, Tₛ)

"""
Returns dimensional form of (`Unitful`) time duration.
"""
redimensionalize(t::T1, a::A, μ₁::U1, μ₂::U2) where {
        T1<:Real, A<:Unitful.Length, U1<:MassParameter, U2<:MassParameter
    } = redimensionalize(t, time_scale_factor(a, μ₁, μ₂))

"""
Returns the spacecraft's nondimensional position w.r.t. body 1 (or 2).
"""
nondimensional_radius(r, xᵢ=0) = √( (r[1]-xᵢ)^2 + r[2]^2 + r[3]^2 )

#=
"""
Returns synodic distance to primary body.
"""
distance_to_primary(orbit::CR3BPOrbit) = nondimensional_radius((position_vector ∘ synodic ∘ normalize)(orbit), -normalized_mass_parameter(orbit.system))

"""
Returns synodic distance to secondary body.
"""
distance_to_secondary(orbit::CR3BPOrbit) = nondimensional_radius((position_vector ∘ synodic ∘ normalize)(orbit), 1-normalized_mass_parameter(orbit.system))
=#

"""
Returns the potential energy `U` in the Synodic frame with Normalized units.
"""
potential_energy(r, μ) = (r[1]^2 + r[2]^2)/2 + ((1-μ)/nondimensional_radius(r,-μ)) + (μ/nondimensional_radius(r,1-μ))

#=
"""
Returns the potential energy `U` in the Synodic frame with Normalized units.
"""
function potential_energy(orbit::CircularRestrictedThreeBodyOrbit) 
    orb = (normalize ∘ synodic)(orbit)
    return potential_energy(position_vector(orb), normalized_mass_parameter(orb.system))
end
=#

"""
Returns the Jacobi Constant, `C` in the Synodic frame with Normalized units.
"""
jacobi_constant(r, v, μ) = 2*potential_energy(r, μ) - (v⋅v)

#=
"""
Returns the Jacobi Constant, `C` in the Synodic frame with Normalized units.
"""
function jacobi_constant(orbit::CircularRestrictedThreeBodyOrbit)
    orb = (normalize ∘ synodic)(orbit)
    return jacobi_constant(position_vector(orb), velocity_vector(orb), normalized_mass_parameter(orb.system))
end
=#

"""
Given the Synodic frame vector, returns the vector in the barycenteric-inertial reference frame.
"""
function inertial(vecₛ::AbstractVector, t, ω::Unitful.AbstractQuantity=1.0u"rad"/unit(t))

    θ = ω*t
    ᴵTₛ = @SMatrix [
        cos(θ) sin(θ) 0
       -sin(θ) cos(θ) 0
        0      0      1
    ]

    return  ᴵTₛ * vecₛ

end

#=
"""
Given an `InertialCartesianState`, returns the state in the synodic (rotating) reference frame.
"""
function synodic(state::NormalizedCartesianState{F, Inertial}, ω::Unitful.AbstractQuantity=1.0u"rad"/AstrodynamicsCore.timeunit(state)) where F

    t = epoch(state)
    θ = ω*t
    ˢTᵢ = inv(SMatrix{3,3}([
        cos(θ) sin(θ) 0
       -sin(θ) cos(θ) 0
        0      0      1
    ]))

    return CartesianState(ˢTᵢ * state.r, ˢTᵢ * state.v, state.t, Synodic; lengthunit = AstrodynamicsCore.lengthunit(state), timeunit = AstrodynamicsCore.timeunit(state))
end
=#

#=
"""
Given an `InertialCartesianState`, returns the state in the synodic (rotating) reference frame.
"""
function synodic(state::SynodicCartesianSTMState{F}, ω::Unitful.AbstractQuantity=1.0u"rad"/AstrodynamicsCore.timeunit(state)) where F
    return state
end

"""
Given a `SynodicCartesianState`, returns the state in the inertial reference frame.
"""
function inertial(state::NormalizedCartesianState{F, Synodic}, ω::Unitful.AbstractQuantity=(1.0u"rad")/AstrodynamicsCore.timeunit(state)) where F

    t = epoch(state)
    θ = ω*t
    ⁱTₛ = @MMatrix [
        cos(θ) sin(θ) 0
       -sin(θ) cos(θ) 0
        0      0      1
    ]

    return CartesianState(ⁱTₛ * state.r, ⁱTₛ * state.v, state.t, Inertial; lengthunit = AstrodynamicsCore.lengthunit(state), timeunit = AstrodynamicsCore.timeunit(state))
end

"""
Given a `SynodicCartesianState`, returns the state in the inertial reference frame.
"""
function inertial(state::NormalizedSynodicCartesianSTMState{F}, ω::Unitful.AbstractQuantity=(1.0u"rad")/AstrodynamicsCore.timeunit(state)) where F

    t = epoch(state.cart)
    θ = ω*t
    ⁱTₛ = @MMatrix [
        cos(θ) sin(θ) 0
       -sin(θ) cos(θ) 0
        0      0      1
    ]

    return CartesianState(ⁱTₛ * state.cart.r, ⁱTₛ * state.cart.v, state.cart.t, Inertial; lengthunit = AstrodynamicsCore.lengthunit(state.cart), timeunit = AstrodynamicsCore.timeunit(state.cart))
end

"""
If given an inertial `CartesianState`, no operation needed.
"""
function inertial(state::InertialCartesianState, ω::Unitful.AbstractQuantity=1.0u"rad"/AstrodynamicsCore.timeunit(state))
    return state
end

"""
If given a Synodic `CartesianState`, no operation needed.
"""
function synodic(state::SynodicCartesianState, ω::Unitful.AbstractQuantity=1.0u"rad"/AstrodynamicsCore.timeunit(state))
    return state
end

"""
Given a CR3BP orbit, returns the orbit in the inertial reference frame.
"""
function inertial(orb::CircularRestrictedThreeBodyOrbit)
    if !(orb.state isa NormalizedCartesianState) # TODO not type stable
        state = inertial(normalize(orb).state)
        return redimensionalize(CircularRestrictedThreeBodyOrbit(state, orb.system))
    else
        state = inertial(orb.state)
        return CircularRestrictedThreeBodyOrbit(state, orb.system)
    end
end

"""
Given a CR3BP orbit, returns the orbit in the Synodic (rotating) reference frame.
"""
function synodic(orb::CircularRestrictedThreeBodyOrbit)
    if !(orb.state isa NormalizedCartesianState) # TODO not type stable
        state = synodic(normalize(orb).state)
        return redimensionalize(CircularRestrictedThreeBodyOrbit(state, orb.system))
    else
        state = synodic(orb.state)
        return CircularRestrictedThreeBodyOrbit(state, orb.system)
    end
end

"""
Position of primary body.
"""
function primary_synodic_position(orb::NormalizedSynodicCR3BPOrbit)
    μ = normalized_mass_parameter(orb)
    SVector{3}(
        -μ, 
        zero(μ),
        zero(μ)
    )
end

"""
Position of secondary body.
"""
function secondary_synodic_position(orb::NormalizedSynodicCR3BPOrbit)
    μ = normalized_mass_parameter(orb)
    SVector{3}(
        one(μ) - μ, 
        zero(μ),
        zero(μ)
    )
end

"""
Translates the state from one inertial frame to another.
"""
function transform(orb::CircularRestrictedThreeBodyOrbit, inertial_position, inertial_velocity = [0, 0, 0]u"km/s")
    orbit = (inertial ∘ redimensionalize)(orb)

    r = position_vector(orbit) - inertial_position
    v = velocity_vector(orbit) # - inertial_velocity

    return CircularRestrictedThreeBodyOrbit(r, v, orb.system, epoch(orbit); frame = Inertial)

end

"""
Translates the state to the primary body center of mass (inertial frame).
"""
function transform_to_primary(orb::CircularRestrictedThreeBodyOrbit)
    return transform(orb, inertial(
            (primary_synodic_position ∘ synodic ∘ normalize)(orb) .* 
                normalized_length_unit(orb.system),  
            (epoch ∘ redimensionalize)(orb)))
end

"""
Translates the state to the secondary body center of mass (inertial frame).
"""
function transform_to_secondary(orb::CircularRestrictedThreeBodyOrbit)
    return transform(orb, inertial(
            (secondary_synodic_position ∘ synodic ∘ normalize)(orb) .* 
                normalized_length_unit(orb.system),  
            (epoch ∘ redimensionalize)(orb)))
end
=#

"""
Returns the position and velocity vectors in the synodic (rotating) reference frame.
"""
LinearAlgebra.normalize(rᵢ, vᵢ, a, Tₛ) =  nondimensionalize(rᵢ, a), nondimensionalize(vᵢ, a, Tₛ)

#=
"""
Normalize a `CartesianState`, given a `CircularRestrictedThreeBodySystem`.
"""
function LinearAlgebra.normalize(cart::CartesianState{F, LU, TU}, sys::CircularRestrictedThreeBodySystem) where {F, LU, TU}
    if LU == NormalizedLengthUnit && TU == NormalizedTimeUnit
        return cart 
    else
       r = nondimensionalize(position_vector(cart), normalized_length_unit(sys))
       v = nondimensionalize(velocity_vector(cart), normalized_length_unit(sys),  normalized_time_unit(sys))
       t = nondimensionalize(epoch(cart), normalized_time_unit(sys))
       return CartesianState(r, v, t, Synodic; lengthunit = NormalizedLengthUnit(), timeunit = NormalizedTimeUnit())
    end
end

"""
Normalize a `CartesianState`, given a `CircularRestrictedThreeBodySystem`.
"""
function LinearAlgebra.normalize(orbit::CircularRestrictedThreeBodyOrbit)
    cart = orbit.state
    sys  = orbit.system
    LU = typeof(lengthunit(cart))
    TU = typeof(timeunit(cart))
    if LU == NormalizedLengthUnit && TU == NormalizedTimeUnit
        return orbit 
    else
        r = nondimensionalize(position_vector(cart), normalized_length_unit(sys))
        v = nondimensionalize(velocity_vector(cart), normalized_length_unit(sys),  normalized_time_unit(sys))
        t = nondimensionalize(epoch(cart), normalized_time_unit(sys))
       return CircularRestrictedThreeBodyOrbit(
           CartesianState(r, v,t, Synodic; 
                          lengthunit = NormalizedLengthUnit(), 
                          timeunit = NormalizedTimeUnit()), 
           sys
        )
    end
end

"""
Un-normalizes (or re-dimensionalizes) a `NormalizedCartesianState` back into a `CartesianState`.
"""
function redimensionalize(cart::CartesianState{F, LU, TU}, sys::CircularRestrictedThreeBodySystem) where {F, LU, TU}
    if LU != NormalizedLengthUnit && TU != NormalizedTimeUnit
        return cart
    else
        r = redimensionalize_length(position_vector(cart), normalized_length_unit(sys))
        v = redimensionalize_velocity(velocity_vector(cart), normalized_length_unit(sys),  normalized_time_unit(sys))
        t = redimensionalize_time(epoch(cart), normalized_time_unit(sys))
       return CartesianState(r, v,t, coordinateframe(cart))
    end
end

"""
Un-normalizes (or re-dimensionalizes) a `NormalizedCartesianState` back into a `CartesianState`.
"""
function redimensionalize(orb::CircularRestrictedThreeBodyOrbit)
    cart = orb.state
    sys  = orb.system
    LU = typeof(lengthunit(cart))
    TU = typeof(timeunit(cart))
    if LU != NormalizedLengthUnit && TU != NormalizedTimeUnit
        return orb
    else
        r = redimensionalize_length(position_vector(cart), normalized_length_unit(sys))
        v = redimensionalize_velocity(velocity_vector(cart), normalized_length_unit(sys),  normalized_time_unit(sys))
        t = redimensionalize_time(epoch(cart), normalized_time_unit(sys))
       return CircularRestrictedThreeBodyOrbit(
           CartesianState(r, v,t, coordinateframe(cart)),
           sys
       )
    end
end
=#

"""
Returns the lagrange points for a CR3BP system.

__Arguments:__ 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `L`: Langrange points requested, must be in range [1,5]

__Outputs:__
- Tuple of Lagrange points
- Throws `ArgumentError` if L is out of range [1,5]

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/)
"""
function lagrange(μ::Real, L=1:5)
    
    if !all(L[i] ∈ (1,2,3,4,5) for i ∈ 1:length(L))
        throw(ArgumentError("Requested lagrange points must be in range [1,5]"))
    end

	expressions = @SVector [
		x -> x - (1-μ)/(x+μ)^2 + μ/(x+μ-1)^2,
		x -> x - (1-μ)/(x+μ)^2 - μ/(x+μ-1)^2,
		x -> x + (1-μ)/(x+μ)^2 + μ/(x+μ+1)^2
	]
	
	return  (map(f->[find_zero(f, (-3,3)), 0, 0], expressions)..., 
			[(1/2) - μ, √(3)/2, 0], [(1/2) - μ, -√(3)/2, 0])[L]
	
end

#=
"""
Returns the lagrange points for a CR3BP system.

__Arguments:__ 
- `sys`: Non-dimensional mass parameter for the CR3BP system.
- `L`: Langrange points requested, must be in range [1,5]

__Outputs:__
- Tuple of Lagrange points
- Throws `ArgumentError` if L is out of range [1,5]

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/)
"""
function lagrange(sys::CircularRestrictedThreeBodySystem, L=1:5)
    
    return lagrange(normalized_mass_parameter(sys), L)
end

"""
Returns non-dimensional acceleration for CR3BP state.
"""
function accel!(aₛ, rₛ, vₛ, μ)

    x₁ = -μ
    x₂ = 1-μ
    r₁ = nondimensional_radius(rₛ, x₁)
    r₂ = nondimensional_radius(rₛ, x₂)

    aₛ[1] =  2vₛ[2] + rₛ[1] - (1-μ)*(rₛ[1] - x₁) / r₁^3 - μ*(rₛ[1] - x₂) / r₂^3
    aₛ[2] = -2vₛ[1] + rₛ[2] - ((1-μ) / r₁^3 + (μ / r₂^3)) * rₛ[2]
    aₛ[3] = -((1-μ) / r₁^3 + (μ / r₂^3)) * rₛ[3]
    return nothing
end

"""
Returns non-dimensional acceleration for CR3BP state.
"""
function accel(rₛ, vₛ, μ)
    aₛ = similar(vₛ)
    accel!(aₛ, rₛ, vₛ, μ)
    return aₛ
end

=#
"""
Returns an analytical solution for a Halo orbit about `L`.

__Arguments:__ 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `Az`: Desired non-dimensional Z-amplitude for Halo orbit.
- `ϕ`: Desired Halo orbit phase.
- `steps`: Number of non-dimensional timepoints in returned state.
- `L`: Lagrange point to orbit (L1 or L2).
- `hemisphere`: Specifies northern or southern Halo orbit.

__Outputs:__
- Synodic position vector `r::Array{<:AbstractFloat}`
- Synodic velocity vector `v::Array{<:Abstractfloat}`.
- Halo orbit period `Τ`.
- Throws `ArgumentError` if L is not `1` or `2`.

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function analyticalhalo(μ; Az=0.00, ϕ=0.0, steps=1,
                       L=1, hemisphere=:northern)

    if L == 1
        point = first(lagrange(μ, 1))
        γ = abs(one(μ) - μ - point)
        n = collect(1:4) .* one(μ)
        c = @. (μ + (-one(1))^n * (one(μ)-μ)γ^(n+1)) / (γ^3 * (one(μ) - γ^(n+1)))
    elseif L == 2
        point = first(lagrange(μ, 2))
        γ = abs(point - one(μ) + μ)
        n = collect(1:4) .* one(μ)
        c = @. ((-one(μ))^n * μ + (-one(μ))^n * (one(μ)-μ)γ^(n+1)) / (γ^3 * (one(μ) + γ^(n+1)))
    else
        throw(ArgumentError("Only Halo orbits about L1 or L2 are supported."))
    end

    ωₚ  = √((2 - c[2] + √((9c[2]^2 - 8c[2])))/2)
    k   = (ωₚ^2 + 1 + 2c[2]) / (2ωₚ)

    d₁  = (3ωₚ^2 / k) * (k*(6ωₚ^2 - 1) - 2ωₚ)
    d₂  = (8ωₚ^2 / k) * (k*(11ωₚ^2 - 1) - 2ωₚ)
    a₂₁ = (3c[3] * (k^2 - 2)) / (4(1 + 2c[2]))
    a₂₂ = (3c[3]) / (4(1 + 2c[2]))
    a₂₃ = (-3c[3]ωₚ / (4k*d₁)) * (3k^3 * ωₚ - 6k*(k-ωₚ) + 4)
    a₂₄ = (-3c[3]ωₚ / (4k*d₁)) * (2 + 3k*ωₚ)
    b₂₁ = (-3c[3]ωₚ / (2d₁)) * (3k*ωₚ - 4)
    b₂₂ = -3c[3]*ωₚ / d₁
    d₂₁ = -c[3] / (2ωₚ^2)
    a₃₁ = (-9ωₚ / (4d₂)) * (4c[3] * (k*a₂₃-b₂₁) + k*c[4]*(4+k^2)) + 
          ((9ωₚ^2 + 1 - c[2]) / (2d₂)) * (3c[3]*(2a₂₃-k*b₂₁) + c[4]*(2+3k^2))
    a₃₂ = (-9ωₚ / (4d₂)) * (4c[3] * (3k*a₂₄-b₂₂) + k*c[4]) -
          (3 / (2d₂)) * (9ωₚ^2 + 1 - c[2]) * (c[3]*(k*b₂₂+d₂₁-2a₂₄) - c[4])
    b₃₁ = (3 / (8d₂)) * 8ωₚ * (3c[3] * (k*b₂₁ - 2a₂₃) - c[4]*(2+3k^2)) + 
          (3/(8d₂)) * (9ωₚ^2 + 1 + 2c[2]) * (4c[3]*(k*a₂₃-b₂₁) + k*c[4]*(4+k^2))
    b₃₂ = (9ωₚ/d₂)*(c[3]*(k*b₂₂+d₂₁-2a₂₄)-c[4]) + 
          (3(9ωₚ^2 + 1 + 2c[2]) / (8d₂) * (4c[3]*(k*a₂₄-b₂₂)+k*c[4]))
    d₃₁ = (3 / (64ωₚ^2)) * (4c[3]*a₂₄ + c[4])
    d₃₂ = (3 / (64 + ωₚ^2)) * (4c[3]*(a₂₃ - d₂₁) + c[4]*(4+k^2))

    s₁  = (1 / (2ωₚ*(ωₚ*(1+k^2) - 2k))) * 
          (3c[3]/2 * (2a₂₁*(k^2 - 2) - a₂₃*(k^2 + 2) - 2k*b₂₁) - 
            (3c[4]/8) * (3k^4 - 8k^2 + 8))
    s₂  = (1 / (2ωₚ*(ωₚ*(1+k^2) - 2k))) * 
          (3c[3]/2 * (2a₂₂*(k^2-2) + a₂₄*(k^2 + 2) + 2k*b₂₂ + 5d₂₁) + 
            (3c[4]/8) * (12 - k^2))
    l₁  = (-3c[3] / 2) * (2a₂₁ + a₂₃ + 5d₂₁) - (3c[4]/8)*(12 - k^2) + 2ωₚ^2 * s₁
    l₂  = (3c[3]/2) * (a₂₄ - 2a₂₂) + (9c[4]/8) + 2ωₚ^2 * s₂
    Δ   = ωₚ^2 - c[2]

    Aᵧ  = Az / γ
    Aₓ  = √((-l₂*Aᵧ^2 - Δ) / l₁)

    ν   = 1 + s₁*Aₓ^2 + s₂*Aᵧ^2
    Τ   = 2π / (ωₚ*ν)
    τ   = ν .* (steps > 1 ? range(0, stop=Τ, length=steps) : range(0, stop=Τ, length=1000))

    if hemisphere == :northern
        m = 1.0
    elseif hemisphere == :southern
        m = 3.0
    else
        throw(ArgumentError("`hemisphere` must be `:northern` or `:southern`."))
    end

    δₘ = 2 - m
    τ₁ = @. ωₚ*τ + ϕ

    x = @. γ * (a₂₁*Aₓ^2 + a₂₂*Aᵧ^2 - Aₓ*cos(τ₁) + (a₂₃*Aₓ^2 - 
                    a₂₄*Aᵧ^2)*cos(2τ₁) + (a₃₁*Aₓ^3 - a₃₂*Aₓ*Aᵧ^2)*cos(3τ₁)) + 1 - μ - (L == 1 ? γ : -γ)
    y = @. γ * (k*Aₓ*sin(τ₁) + (b₂₁*Aₓ^2 - b₂₂*Aᵧ^2)*sin(2τ₁) + 
                    (b₃₁*Aₓ^3 - b₃₂*Aₓ*Aᵧ^2)*sin(3τ₁))
    z = @. γ * (δₘ*Aᵧ*cos(τ₁) + δₘ*d₂₁*Aₓ*Aᵧ*(cos(2τ₁)-3) + 
                    δₘ*(d₃₂*Aᵧ*Aₓ^2 - d₃₁*Aᵧ^3)*cos(3τ₁))

    ẋ = @. γ * (ωₚ*ν*Aₓ*sin(τ₁) - 2ωₚ*ν*(a₂₃*Aₓ^2 - a₂₄*Aᵧ^2)*sin(2τ₁) - 
                    3ωₚ*ν*(a₃₁*Aₓ^3 - a₃₂*Aₓ*Aᵧ^2)*sin(3τ₁))
    ẏ = @. γ * (ωₚ*ν*k*Aₓ*cos(τ₁) + 2ωₚ*ν*(b₂₁*Aₓ^2 - b₂₂*Aᵧ^2)*cos(2τ₁) + 
                    3ωₚ*ν*(b₃₁*Aₓ^3 - b₃₂*Aₓ*Aᵧ^2)*cos(3τ₁))
    ż = @. γ * (-ωₚ*ν*δₘ*Aᵧ*sin(τ₁) - 2ωₚ*ν*δₘ*d₂₁*Aₓ*Aᵧ*sin(2τ₁) - 
                    3ωₚ*ν*δₘ*(d₃₂*Aᵧ*Aₓ^2 - d₃₁*Aᵧ^2)*sin(3τ₁))

    return hcat(x, y, z)[1:steps, :], hcat(ẋ, ẏ, ż)[1:steps, :], Τ

end

#=
"""
Returns a `Vector` of `Matrix` values.
Each `Matrix` contains a 3-column nondimensional
position trajectory in the Synodic frame which
represents a Zero Velocity Curve.
"""
function zerovelocity_curves(orbit::CircularRestrictedThreeBodyOrbit;
                             nondimensional_range = range(-2; stop=2, length=1000))

    x = nondimensional_range
    y = nondimensional_range

    rᵢ = orbit |> synodic |> normalize |> position_vector 
    vᵢ = orbit |> synodic |> normalize |> velocity_vector 
    μ  = normalized_mass_parameter(orbit.system)
    LU = string(normalized_length_unit(orbit.system))

    z = [
        2 * potential_energy([xs, ys, 0.0], μ) - jacobi_constant(rᵢ, vᵢ, μ)
        for xs ∈ x, ys ∈ y
    ]

    zvc_contour = contours(x,y,z,[0.0])

    curves = Vector{Matrix{Float64}}()

    for zvc_level ∈ Contour.levels(zvc_contour)

        x = Float64[]
        y = Float64[]
        for zvc_line ∈ Contour.lines(zvc_level)
            xs, ys = coordinates(zvc_line)
            if length(x) == length(y) == 0
                x = xs
                y = ys
            else
                x = vcat(x, xs)
                y = vcat(y, ys)
            end
        end

        if length(curves) == 0
            curves = [hcat(x,y)]
        else
            curves = vcat(curves, [hcat(x, y)])
        end
    end

    return curves
                         
end
=#