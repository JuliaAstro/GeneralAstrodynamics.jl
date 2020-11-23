#
# ThreeBodyCalculations.jl
# Calculations for the Circular Restricted
# Three Body Problem.
#

"""
Returns time scale factor, `Tₛ`.
"""
time_scale_factor(a, μ₁, μ₂) = orbital_period(a, μ₁+μ₂)

"""
Returns nondimensional length unit, `DU`.
"""
nondimensionalize_length(rᵢ, a) = upreferred(rᵢ .÷ a)

"""
Returns nondimensional velocity unit, `DU/DT`.
"""
nondimensionalize_velocity(vᵢ, a, Tₛ) = vᵢ ./ (a ÷ Tₛ)

"""
Returns nondimensional time unit, `DT`.
"""
nondimensionalize_time(t, a, μ₁, μ₂) = t ÷ time_scale_factor(a, μ₁, μ₂)

"""
Returns nondimensional mass parameter, `μ`.
"""
nondimensionalize_mass_parameter(μ₁, μ₂) = min(μ₁,μ₂) ÷ (μ₁+μ₂)

"""
Returns nondimensional form of (`Unitful`) position or velocity vectors, 
time duration, or mass parameter for the Circular Restricted Three-body
problem.
"""
nondimensionalize(rᵢ::R, a::A) where {
        R<:Length, A<:Length
    } = nondimensionalize_length(rᵢ, a)
nondimensionalize(rᵢ::R, a::A) where {
        U<:Length, R<:AbstractVector{U}, A<:Length
    } = nondimensionalize_length(rᵢ, a)
nondimensionalize(vᵢ::V, a::A) where {
        V<:Velocity, A<:Length
    } = nondimensionalize_length(vᵢ, a)
nondimensionalize(vᵢ::V, a::A, Tₛ::T) where {
        U<:Velocity, V<:AbstractVector{U}, A<:Length, T<:Time
    } = nondimensionalize_velocity(vᵢ, a, Tₛ)
nondimensionalize(vᵢ::V, a::A, μ₁::U1, μ₂::U2) where {
        U<:Velocity, V<:AbstractVector{U}, A<:Length, U1<:MassParameter, U2<:MassParameter
    } = nondimensionalize_velocity(vᵢ, a, time_scale_factor(a, μ₁, μ₂))
nondimensionalize(t::T1, Tₛ::T2) where {
        T1<:Time, T2<:Time
    } = t ÷ Tₛ
nondimensionalize(t::T1, a::A, μ₁::U1, μ₂::U2) where {
        T1<:Time, A<:Length, U1<:MassParameter, U2<:MassParameter
    } = nondimensionalize(t, time_scale_factor(a, μ₁, μ₂))
nondimensionalize(μ₁::U1, μ₂::U2) where {
        U1<:MassParameter, U2<:MassParameter
    } = min(μ₁, μ₂) ÷ (μ₁+μ₂)
function nondimensionalize(r₃::R, v₃::V, Δt::T, μ₁::U1, μ₂::U2, a::A) where {
        RT<:Length, R<:AbstractVector{RT},
        VT<:Velocity, V<:AbstractVector{VT},
        T<:Time, U1<:MassParameter, U2<:MassParameter,
        A<:Length
    }

    Tₛ = time_scale_factor(a, μ₁, μ₁)
    return nondimensionalize(r₃, a), 
           nondimensionalize(v₃, a, Tₛ),
           nondimensionalize(Δt, Tₛ),
           nondimensionalize(μ₁, μ₂)

end

"""
Returns dimensional length units.
"""
redimensionalize_length(rᵢ, a) = upreferred(rᵢ .* a)

"""
Returns dimensional velocity units.
"""
redimensionalize_velocity(vᵢ, a, Tₛ) = vᵢ .* (a ÷ Tₛ)

"""
Returns dimensional time unit.
"""
redimensionalize_time(t, a, μ₁, μ₂) = t * time_scale_factor(a, μ₁, μ₂)

"""
Returns dimensional (inertial) form of (`Unitful`) position or velocity vectors, 
time duration, or mass parameter for the Circular Restricted Three-body
problem.
"""
redimensionalize(rᵢ::R, a::A) where {
        R<:Length, A<:Length
    } = redimensionalize_length(rᵢ, a)
redimensionalize(rᵢ::R, a::A) where {
        U<:Length, R<:AbstractVector{U}, A<:Length
    } = redimensionalize_length(rᵢ, a)
redimensionalize(vᵢ::V, a::A) where {
        V<:Velocity, A<:Length
    } = redimensionalize_length(vᵢ, a)
redimensionalize(vᵢ::V, a::A, Tₛ::T) where {
        U<:Velocity, V<:AbstractVector{U}, A<:Length, T<:Time
    } = redimensionalize_velocity(vᵢ, a, Tₛ)
redimensionalize(vᵢ::V, a::A, μ₁::U1, μ₂::U2) where {
        U<:Velocity, V<:AbstractVector{U}, A<:Length, U1<:MassParameter, U2<:MassParameter
    } = redimensionalize_velocity(vᵢ, a, time_scale_factor(a, μ₁, μ₂))
redimensionalize(t::T1, Tₛ::T2) where {
        T1<:Time, T2<:Time
    } = t ÷ Tₛ
redimensionalize(t::T1, a::A, μ₁::U1, μ₂::U2) where {
        T1<:Time, A<:Length, U1<:MassParameter, U2<:MassParameter
    } = redimensionalize(t, time_scale_factor(a, μ₁, μ₂))

"""
Returns the spacecraft's nondimensional position w.r.t. body 1 (or 2).
"""
nondimensional_radius(r, xᵢ=0) = √( (r[1]-xᵢ)^2 + r[2]^2 + r[3]^2 )

"""
Returns the potential energy `U`.
"""
potential_energy(r, μ, x₁, x₂) = (r[1]^2 + r[2]^2) + (2(1-μ)/nondimensional_radius(r,x₁)) + (2μ/nondimensional_radius(r,x₂))

"""
Returns the Jacobi Constant `C`.
"""
jacobi_constant(r, v, μ, x₁, x₂) = 2*potential_energy(r, μ, x₁, x₂) - (v⋅v)

"""
Given the Synodic frame vector, returns the vector in the inertial reference frame.
"""
function inertial(vecₛ, t, ω=1.0u"rad"÷unit(t))

    θ = ω*t
    ᴵTₛ = @SMatrix [
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