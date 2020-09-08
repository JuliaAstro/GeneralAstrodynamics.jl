#= Orbits.jl

Define structures & functions for orbits

=#

## Dependencies 
using StaticArrays
using LinearAlgebra: ×, ⋅, norm
using Unitful, UnitfulAstro, UnitfulAngles

## Data Structures

# Orbits ~orbit~ arround a central body
@enum CelestialBody earth sun

# Abstract type for all Orbits
abstract type AbstractOrbit end

# Cartesian parameters
struct CartesianOrbit <: AbstractOrbit

    r̲::SVector{3,Unitful.Length}
    v̲::SVector{3,Unitful.Velocity}
    body::CelestialBody

end

# Keplarian Parameters
struct CanonicalOrbit <: AbstractOrbit

    e::Unitful.DimensionlessQuantity
    a::Unitful.Length
    i::Unitful.Quantity
    Ω::Unitful.Quantity
    ω::Unitful.Quantity
    ν::Unitful.Quantity
    body::CelestialBody

end

## Functions

# Constructors

function CartesianOrbit(r̲::AbstractArray{Unitful.Length}, 
                        v̲::AbstractArray{Unitful.Velocity}, 
                        body::CelestialBody)

    return CartesianOrbit(
            SVector{3,Unitful.Quantity}(r̲),
            SVector{3,Unitful.Quantity}(v̲),
            body)
end

function CanonicalOrbit(e::Unitful.DimensionlessQuantity, 
                        a::Unitful.Length, 
                        i::Unitful.DimensionlessQuantity, 
                        Ω::Unitful.DimensionlessQuantity, 
                        ω::Unitful.DimensionlessQuantity, 
                        ν::Unitful.DimensionlessQuantity, 
                        body::CelestialBody)

    return CanonicalOrbit(
            e, a, i, 
            Ω, ω, ν, body)

end

# Calculation Functions

function semimajor_axis(r, v, μ)
   
    return inv( (2 / r) - (v^2 / μ) )

end

function semimajor_axis(orbit::CanonicalOrbit)

    return orbit.a

end

function semimajor_axis(orbit::CartesianOrbit; 
                        μ=(SVector(
                            UnitfulAstro.GMearth, 
                            UnitfulAstro.GMsun)[Int(orbit.body)+1]))

    return semimajor_axis(
                norm(orbit.r̲),
                norm(orbit.v̲),
                μ)

end

function specific_angular_momentum_vector(r̲, v̲)

    return r̲ × v̲

end

function specific_angular_momentum_vector(orbit::CartesianOrbit)

    return specific_angular_momentum_vector(orbit.r̲, orbit.v̲)

end

function specific_angular_momentum(r̲, v̲)

    return norm(specific_angular_momentum_vector(r̲, v̲))

end

function specific_angular_momentum(orbit::CartesianOrbit)

    return specific_angular_momentum(orbit.r̲, orbit.v̲)

end

function specific_angular_momentum(
            orbit::AbstractOrbit, 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return apoapsis_radius(orbit) * apoapsis_velocity(orbit)

end

function specific_energy(a, μ)

    return ( -μ / (2 * a) )

end

function specific_energy(r, v, μ)

    return (v^2 / 2) - (μ / r)

end

function specific_energy(
            orbit::AbstractOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return specific_energy(semimajor_axis(orbit), μ)

end

function eccentricity(r̲, v̲, μ)

    return (1 / μ) * 
            ((v̲ × specific_angular_momentum_vector(r̲, v̲)) - 
                μ * r̲ / norm(r̲))

end

function eccentricity_vector(
            orbit::CartesianOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return eccentricity(orbit.r̲, orbit.v̲, μ)

end

function eccentricity(
            orbit::CartesianOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return eccentricity(orbit.r̲, orbit.v̲, μ)

end

function eccentricity(orbit::CanonicalOrbit)

    return orbit.e

end

function semi_parameter(a, e)

    return a * (1 - e^2)

end

function semi_parameter(orbit::AbstractOrbit)

    return semimajor_axis(orbit) * 
            (1 - eccentricity(orbit)^2)

end

function semi_parameter(orbit::CanonicalOrbit)

    return orbit.a * (1 - orbit.e^2)

end

function periapsis_radius(a, e)

    return a * (1 - e)

end

function periapsis_radius(orbit::AbstractOrbit)

    return periapsis_radius(semimajor_axis(orbit), 
                            norm(eccentricity(orbit)))

end

function apoapsis_radius(a, e)

    return a * (1 + e)

end

function apoapsis_radius(orbit::AbstractOrbit)

    return apoapsis_radius(semimajor_axis(orbit), 
                           norm(eccentricity(orbit)))

end

function instantaneous_velocity(r, a, μ)

    return √( (2 * μ / r) - (μ / a))

end

function periapsis_velocity(
            orbit::AbstractOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return instantaneous_velocity(
                periapsis_radius(orbit), 
                semimajor_axis(orbit),
                μ)

end

function apoapsis_velocity(
            orbit::AbstractOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return instantaneous_velocity(
                apoapsis_radius(orbit), 
                semimajor_axis(orbit),
                μ)

end

function orbital_period(
            orbit::AbstractOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return 2 * π * √(semimajor_axis(orbit)^3 / μ)

end

function true_anomoly(r, h, e, μ)

    return acos( (h^2 - μ * r) / (μ * r * e) )

end

function true_anomoly(
            orbit::CartesianOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return true_anomoly(
            norm(orbit.r̲),
            specific_angular_momentum(orbit),
            eccentricity(orbit),
            μ)

end

function true_anomoly(orbit::CanonicalOrbit)

    return orbit.ν

end

function instantaneous_radius(p, e, ν)

    return p / (1 + e * cos(ν))

end

function instantaneous_radius(orbit::AbstractOrbit)

    return instantaneous_radius(semi_parameter(orbit),
                                eccentricity(orbit),
                                true_anomoly(orbit))

end

function mean_motion(a, μ)

    return √(μ / a^3)

end

function mean_motion_vector(
            orbit::CartesianOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1],
            î=SVector{3, Float64}([1, 0, 0]), 
            ĵ=SVector{3, Float64}([0, 1, 0]), 
            k̂=SVector{3, Float64}([0, 0, 1]))

    return k̂ × specific_angular_momentum_vector(orbit)

end

function mean_motion(
            orbit::AbstractOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return mean_motion(
            semimajor_axis(orbit),
            μ)

end

function eccentric_anomoly(e, ν)

    return atan( u"rad", (√(1 - e^2) * sin(ν)) / (e + cos(ν)) )

end

function eccentric_anomoly(orbit::AbstractOrbit)

    return eccentric_anomoly(
            norm(eccentricity(orbit)),
            true_anomoly(orbit))

end

function time_since_periapsis(n, e, E)

    return (E - e * sin(E)) / (n)

end

function time_since_periapsis(orbit::AbstractOrbit)

    return time_since_periapsis(
            mean_motion(orbit),
            norm(eccentricity(orbit)),
            eccentric_anomoly(orbit))

end

function inclination(orbit::CartesianOrbit)

    h_vec = specific_angular_momentum_vector(orbit)
    return acos(u"rad", h_vec[3] / norm(h_vec))

end

function inclination(orbit::CanonicalOrbit)

    return orbit.i

end

# Transformation Functions

function CartesianToCanonical(
            orbit::CartesianOrbit; 
            î=SVector{3, Float64}([1, 0, 0]), 
            ĵ=SVector{3, Float64}([0, 1, 0]), 
            k̂=SVector{3, Float64}([0, 0, 1]))

    h_vec = specific_angular_momentum_vector(orbit)
    n_vec = mean_motion_vector(orbit)
    e_vec = eccentricity_vector(orbit)

    Ω = ustrip(n_vec ⋅ ĵ) > 0 ? 
            acos(u"rad", (n_vec ⋅ î) / norm(n_vec) ) :
            2 * π * u"rad" - 
                acos(u"rad", (n_vec ⋅ î) / norm(n_vec) )

    ω = ustrip(e_vec ⋅ k̂) > 0 ?
            acos(u"rad", (n_vec ⋅ e_vec) / (norm(n_vec) * norm(e_vec)) ) :
            2 * π * u"rad" - 
                acos(u"rad", (n_vec ⋅ e_vec) / (norm(n_vec) * norm(e_vec)) )

    ν = ustrip(orbit.r̲ ⋅  e_vec) > 0 ? 
            acos(u"rad", (e_vec ⋅ orbit.r̲) / (norm(e_vec) * norm(orbit.r̲)) ) :
            2 * π * u"rad" - 
                acos(u"rad", (e_vec ⋅ orbit.r̲) / (norm(e_vec) * norm(orbit.r̲)) )

    return CanonicalOrbit( 
            norm(e_vec),
            upreferred(semimajor_axis(orbit)),
            acos(u"rad", (h_vec ⋅ k̂) / norm(h_vec) ),
            Ω,
            ω,
            ν,
            orbit.body)

end

function CanonicalToCartesian(orbit::CanonicalOrbit)



end

# Export data structures & constructors
export CelestialBody, earth, sun
export CartesianOrbit, CanonicalOrbit

# Export functions
export  semimajor_axis, 
        eccentricity, 
        eccentricity_vector,
        time_since_periapsis, 
        semi_parameter, 
        mean_motion,
        mean_motion_vector, 
        periapsis_radius, 
        apoapsis_radius,
        specific_angular_momentum_vector,
        specific_angular_momentum, 
        specific_energy,
        instantaneous_velocity, 
        periapsis_vecolity, 
        apoapsis_velocity,
        orbital_period, 
        inclination,
        true_anomoly, 
        eccentric_anomoly,
        CartesianToCanonical, 
        CanonicalToCartesian