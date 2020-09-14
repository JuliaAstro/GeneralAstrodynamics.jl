""" 
    TwoBody.jl

Provides structures & functions for two-body orbits.
"""

### Dependencies 

using Base: isapprox, isequal
using Pkg, Logging
using StaticArrays
using LinearAlgebra: ×, ⋅, norm
using Unitful, UnitfulAstro, UnitfulAngles
using DifferentialEquations
using ComponentArrays

### Export data structures & constructors
export CelestialBody, earth, sun
export CartesianOrbit, CanonicalOrbit

#### Export functions
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
        instantaneous_radius,
        instantaneous_velocity, 
        periapsis_velocity, 
        apoapsis_velocity,
        orbital_period, 
        inclination,
        true_anomoly, 
        eccentric_anomoly,
        CartesianToCanonical, 
        CanonicalToCartesian,
        isapprox,
        isequal,
        propagate,
        parse_propagation

### Data Structures

"""
    CelestialBody

Enum representing large bodies near Earth.
Currently only Earth and the Sun are supported.
"""
@enum CelestialBody earth sun

"""
    AbstractOrbit

Abstract type for all orbits.
"""
abstract type AbstractOrbit end

"""
    CartesianOrbit

Cartesian representation for orbital state.
"""
struct CartesianOrbit <: AbstractOrbit

    r̅::SVector{3, Unitful.Length}
    v̅::SVector{3, Unitful.Velocity}
    body::CelestialBody

end

"""
    CanonicalOrbit

Keplarian representation for orbital state.
"""
struct CanonicalOrbit <: AbstractOrbit

    e
    a::Unitful.Length
    i::Unitful.Quantity
    Ω::Unitful.Quantity
    ω::Unitful.Quantity
    ν::Unitful.Quantity
    body::CelestialBody

end

### Functions

## Constructors

"""
    CartesianOrbit

Constructor for Cartesian orbital representation.
"""
function CartesianOrbit(r̅::AbstractArray{Unitful.Length}, 
                        v̅::AbstractArray{Unitful.Velocity}, 
                        body::CelestialBody)

    return CartesianOrbit(
            SVector{3,Unitful.Quantity{Float64}}(r̅),
            SVector{3,Unitful.Quantity{Float64}}(v̅),
            body)
end

"""
    CanonicalOrbit

Constructor for Keplarian orbital representation.
"""
function CanonicalOrbit(e, 
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

"""
    CanonicalOrbit

Returns a Keplarian representation of an orbital state.
"""
function CanonicalOrbit(orbit::CartesianOrbit)

    î=SVector{3, Float64}([1, 0, 0]) 
    ĵ=SVector{3, Float64}([0, 1, 0]) 
    k̂=SVector{3, Float64}([0, 0, 1])

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

    ν = ustrip(orbit.r̅ ⋅  e_vec) > 0 ? 
            acos(u"rad", (e_vec ⋅ orbit.r̅) / (norm(e_vec) * norm(orbit.r̅)) ) :
            2 * π * u"rad" - 
                acos(u"rad", (e_vec ⋅ orbit.r̅) / (norm(e_vec) * norm(orbit.r̅)) )

    return CanonicalOrbit( 
            norm(e_vec),
            upreferred(semimajor_axis(orbit)),
            acos(u"rad", (h_vec ⋅ k̂) / norm(h_vec) ),
            Ω,
            ω,
            ν,
            orbit.body)

end
CanonicalOrbit(orbit::CanonicalOrbit) = orbit

"""
    CartesianOrbit

Returns a Cartesian representation of an orbital state.
"""
function CartesianOrbit(
            orbit::CanonicalOrbit,
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    # Find semilatus parameter
    p = semi_parameter(orbit.a, orbit.e)

    # Find scalar radius
    r = instantaneous_radius(p, orbit.e, orbit.ν)

    # Set perifocal axes
    P̂=SVector{3, Float64}([1, 0, 0])
    Q̂=SVector{3, Float64}([0, 1, 0]) 
    Ŵ=SVector{3, Float64}([0, 0, 1])

    # Find state in Perifocal frame
    r̅_perifocal = (r * cos(orbit.ν) * P̂ + r * sin(orbit.ν) * Q̂)
    v̅_perifocal = √(μ/p) * (-sin(orbit.ν) * P̂ + (orbit.e + cos(orbit.ν) * Q̂))

    # Set up Perifocal ⟶ Cartesian conversion # TODO - move this to a function
    R_3Ω =  SMatrix{3,3,Float64}(
            [cos(orbit.Ω)           sin(orbit.Ω)            0.;
            -sin(orbit.Ω)           cos(orbit.Ω)            0.;
             0.                     0.                      1.])
    R_1i = SMatrix{3,3,Float64}(
            [1.                     0.                      0.;
             0.                     cos(orbit.i)            sin(orbit.i);
             0.                    -sin(orbit.i)            cos(orbit.i)])
    R_3ω = SMatrix{3,3,Float64}(
            [cos(orbit.ω)           sin(orbit.ω)            0.
            -sin(orbit.ω)           cos(orbit.ω)            0.
             0.                     0.                      1.])

    ᴵTₚ = (R_3ω * R_1i * R_3Ω)' 

    return CartesianOrbit(ᴵTₚ * r̅_perifocal, ᴵTₚ * v̅_perifocal, orbit.body)

end
CartesianOrbit(orbit::CartesianOrbit) = orbit

## Calculation (Helper) Functions

"""
    semimajor_axis

Returns semimajor axis parameter, a.
"""
function semimajor_axis(r, v, μ)
   
    return inv( (2 / r) - (v^2 / μ) )

end

"""
    semimajor_axis

Returns semimajor axis parameter, a, for a Keplarian representation.
"""
function semimajor_axis(orbit::CanonicalOrbit)

    return orbit.a

end

"""
    semimajor_axis

Returns semimajor axis parameter, a, for a Cartesian representation.
"""
function semimajor_axis(orbit::CartesianOrbit; 
                        μ=(SVector(
                            UnitfulAstro.GMearth, 
                            UnitfulAstro.GMsun)[Int(orbit.body)+1]))

    return semimajor_axis(
                norm(orbit.r̅),
                norm(orbit.v̅),
                μ)

end

"""
    specific_angular_momentum_vector

Returns specific angular momentum vector, h̅.
"""
function specific_angular_momentum_vector(r̅, v̅)

    return r̅ × v̅

end

"""
    specific_angular_momentum_vector

Returns specific angular momentum vector, h̅, for a Cartesian representation.
"""
function specific_angular_momentum_vector(orbit::CartesianOrbit)

    return specific_angular_momentum_vector(orbit.r̅, orbit.v̅)

end

"""
    specific_angular_momentum

Returns scalar specific angular momentum vector, h.
"""
function specific_angular_momentum(r̅, v̅)

    return norm(specific_angular_momentum_vector(r̅, v̅))

end

"""
    specific_angular_momentum

Returns scalar specific angular momentum, h, for a Cartesian representation.
"""
function specific_angular_momentum(orbit::CartesianOrbit)

    return specific_angular_momentum(orbit.r̅, orbit.v̅)

end

"""
    specific_angular_momentum

Returns scalar specific angular momentum, h, for any orbital representation.
"""
function specific_angular_momentum(
            orbit::AbstractOrbit, 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return apoapsis_radius(orbit) * apoapsis_velocity(orbit)

end

"""
    specific_energy

Returns specific orbital energy, ϵ.
"""
function specific_energy(a, μ)

    return ( -μ / (2 * a) )

end

"""
    specific_energy

Returns specific orbital energy, ϵ.
"""
function specific_energy(r, v, μ)

    return (v^2 / 2) - (μ / r)

end

"""
    specific_energy

Returns specific orbital energy, ϵ, for any orbital representation.
"""
function specific_energy(
            orbit::AbstractOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return specific_energy(semimajor_axis(orbit), μ)

end

"""
    eccentricity_vector

Returns orbital eccentricity_vector, e̅.
"""
function eccentricity_vector(
            orbit::CartesianOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return eccentricity(orbit.r̅, orbit.v̅, μ)

end

"""
    eccentricity

Returns orbital eccentricity, e.
"""
function eccentricity(r̅, v̅, μ)

    return (1 / μ) * 
            ((v̅ × specific_angular_momentum_vector(r̅, v̅)) - 
                μ * r̅ / norm(r̅))

end

"""
    eccentricity

Returns orbital eccentricity, e.
"""
function eccentricity(
            orbit::CartesianOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return eccentricity(orbit.r̅, orbit.v̅, μ)

end

"""
    eccentricity

Returns orbital eccentricity, e.
"""
function eccentricity(orbit::CanonicalOrbit)

    return orbit.e

end

"""
    semi_parameter

Returns semilatus parameter, p.
"""
function semi_parameter(a, e)

    return a * (1 - e^2)

end

"""
    semi_parameter

Returns semilatus parameter, p, for any orbital representation.
"""
function semi_parameter(orbit::AbstractOrbit)

    return semimajor_axis(orbit) * 
            (1 - eccentricity(orbit)^2)

end

"""
    semi_parameter

Returns semilatus parameter, p, for a Keplarian representation.
"""
function semi_parameter(orbit::CanonicalOrbit)

    return orbit.a * (1 - orbit.e^2)

end

"""
    instantaneous_radius

Returns instantaneous radius, r.
"""
function instantaneous_radius(p, e, ν)

    return p / (1 + e * cos(ν))

end

"""
    instantaneous_radius

Returns instantaneous radius, r, for any orbital representation.
"""
function instantaneous_radius(orbit::AbstractOrbit)

    return instantaneous_radius(semi_parameter(orbit),
                                eccentricity(orbit),
                                true_anomoly(orbit))

end

"""
    instantaneous_velocity

Returns instantaneous velocity, v, for any orbital representation.
"""
function instantaneous_velocity(r, a, μ)

    return √( (2 * μ / r) - (μ / a))

end

"""
    periapsis_radius

Returns periapsis radius, r̅_p.
"""
function periapsis_radius(a, e)

    return a * (1 - e)

end

"""
    periapsis_radius

Returns periapsis radius, r_p, for any orbital representation.
"""
function periapsis_radius(orbit::AbstractOrbit)

    return periapsis_radius(semimajor_axis(orbit), 
                            norm(eccentricity(orbit)))

end

"""
    apoapsis_radius

Returns periapsis radius, r_a.
"""
function apoapsis_radius(a, e)

    return a * (1 + e)

end

"""
    apoapsis_radius

Returns periapsis radius, r_a, for any orbital representation.
"""
function apoapsis_radius(orbit::AbstractOrbit)

    return apoapsis_radius(semimajor_axis(orbit), 
                           eccentricity(orbit))

end

"""
    periapsis_velocity

Returns periapsis velocity, v_p, for any orbital representation.
"""
function periapsis_velocity(
            orbit::AbstractOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return instantaneous_velocity(
                periapsis_radius(orbit), 
                semimajor_axis(orbit),
                μ)

end

"""
    apoapsis_velocity

Returns apoapsis velocity, v_a, for any orbital representation.
"""
function apoapsis_velocity(
            orbit::AbstractOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return instantaneous_velocity(
                apoapsis_radius(orbit), 
                semimajor_axis(orbit),
                μ)

end

"""
    orbital_period

Returns orbital period, Τ, for any orbital representation.
"""
function orbital_period(
            orbit::AbstractOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    P = 2 * π * √(semimajor_axis(orbit)^3 / μ)

end

"""
    true_anomoly

Returns true anomoly, ν.
"""
function true_anomoly(r, h, e, μ)

    return acos( (h^2 - μ * r) / (μ * r * e) )

end

"""
    true_anomoly

Returns true anomoly, ν, for a Cartesian representation.
"""
function true_anomoly(
            orbit::CartesianOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

    return true_anomoly(
            norm(orbit.r̅),
            specific_angular_momentum(orbit),
            eccentricity(orbit),
            μ)

end

"""
    true_anomoly

Returns true anomoly, ν, for a Keplarian representation.
"""
function true_anomoly(orbit::CanonicalOrbit)

    return orbit.ν

end

"""
    mean_motion

Returns mean motion, n.
"""
function mean_motion(a, μ)

    return √(μ / a^3)

end

"""
    mean_motion

Returns mean motion, n, for any orbital representation.
"""
function mean_motion(
    orbit::AbstractOrbit; 
    μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1])

return mean_motion(
    semimajor_axis(orbit),
    μ)

end

"""
    mean_motion_vector

Returns mean motion vector, n̄, for a Cartesian representation.
"""
function mean_motion_vector(
            orbit::CartesianOrbit; 
            μ=SVector(UnitfulAstro.GMearth, UnitfulAstro.GMsun)[Int(orbit.body)+1],
            î=SVector{3, Float64}([1, 0, 0]), 
            ĵ=SVector{3, Float64}([0, 1, 0]), 
            k̂=SVector{3, Float64}([0, 0, 1]))

    return k̂ × specific_angular_momentum_vector(orbit)

end

"""
    eccentric_anomoly

Returns eccentric anomoly, E.
"""
function eccentric_anomoly(e, ν)

    return atan( u"rad", (√(1 - e^2) * sin(ν)) / (e + cos(ν)) )

end

"""
    eccentric_anomoly

Returns eccentric anomoly, E, for any orbital representation.
"""
function eccentric_anomoly(orbit::AbstractOrbit)

    return eccentric_anomoly(
            norm(eccentricity(orbit)),
            true_anomoly(orbit))

end

"""
    time_since_periapsis

Returns time since periapsis, t.
"""
function time_since_periapsis(n, e, E)

    return (E - e * sin(E)) / (n)

end

"""
    time_since_periapsis

Returns time since periapsis, t, for any orbital representation.
"""
function time_since_periapsis(orbit::AbstractOrbit)

    return time_since_periapsis(
            mean_motion(orbit),
            norm(eccentricity(orbit)),
            eccentric_anomoly(orbit))

end

"""
    inclination

Returns orbital inclination, i, for a Cartesian representation.
"""
function inclination(orbit::CartesianOrbit)

    h_vec = specific_angular_momentum_vector(orbit)
    return acos(u"rad", h_vec[3] / norm(h_vec))

end

"""
    inclination

Returns orbital inclination, i, for a Keplarian representation.
"""
function inclination(orbit::CanonicalOrbit)

    return orbit.i

end

function Base.isapprox(c1::CartesianOrbit, c2::CartesianOrbit; tolerance=1e-8)

    return all(ustrip(c1.r̅ - c2.r̅) .< tolerance) &&
           all(ustrip(c1.r̅ - c2.r̅) .< tolerance) &&
           (c1.body == c2.body)

end

function Base.isapprox(c1::CanonicalOrbit, c2::CanonicalOrbit; tolerance=1e-8)

    return ustrip(c1.e - c2.e) .< tolerance &&
           ustrip(c1.a - c2.a) .< tolerance &&
           ustrip(c1.i - c2.i) .< tolerance &&
           ustrip(c1.Ω - c2.Ω) .< tolerance &&
           ustrip(c1.ω - c2.ω) .< tolerance &&
           ustrip(c1.ν - c2.ν) .< tolerance &&
           c1.body == c2.body

end

function Base.isequal(c1::CartesianOrbit, c2::CartesianOrbit)

    return all(c1.r̅ .== c2.r̅) &&
           all(c1.v̅ .== c2.v̅) &&
           c1.body == c2.body

end

function Base.isequal(c1::CanonicalOrbit, c2::CanonicalOrbit)

    return c1.e == c2.e &&
           c1.a == c2.a &&
           c1.i == c2.i &&
           c1.Ω == c2.Ω &&
           c1.ω == c2.ω &&
           c1.ν == c2.ν &&
           c1.body == c2.body

end

function propagate(orbit::AbstractOrbit;
                   μ=SVector(1*UnitfulAstro.GMearth, 1*UnitfulAstro.GMsun)[Int(orbit.body)+1],
                   tspan=(0.0u"s", upreferred(orbital_period(orbit))),
                   reltol=1e-8,
                   saveat=0.25u"s",
                   alg=Tsit5())

    ### Referencing:
    # [1] https://diffeq.sciml.ai/v4.0/tutorials/ode_example.html
    # [2] https://github.com/SciML/DifferentialEquations.jl/issues/393#issuecomment-658210231
    
    # Note: Code from [2] was copied and modified below.
    # I previously couldn't figure out how to use numerical solvers
    # with mixed units (the state vector has length and velocity units).
    # Apparently I was not the only one - one solution (as shown by [2]),
    # is to use ComponentArrays. This also allows me to define a 
    # dynamical tic function that is vectorized.
    ## Set up problem

    # Ensure Cartesian representation
    cart = CartesianOrbit(orbit)
	r₀ = Array(ustrip.(uconvert.(u"m",cart.r̅)))
    v₀ = Array(ustrip.(uconvert.(u"m/s", cart.v̅)))
    
    # Define the problem
    problem = ODEProblem(
                orbit_tic, 
                ComponentArray((r̅=r₀, v̅=v₀)), 
                ustrip.(uconvert.(u"s",tspan)), 
                ComponentArray((μ=ustrip.(uconvert(u"m^3 / s^2", μ)))))

    # Solve the problem! 
    solution = solve(problem, alg,
                     saveat=ustrip(uconvert(u"s",saveat)))
    
    return solution

end

# Note the citation [2] above - this function was copied and
# modified from [2]. I originally had a 6 state function,
# but I had trouble with mixing units. The solution shown
# in [2] allows the use of mixed units through the ComponentArrays
# package.
function orbit_tic(du, u, p, t)
    du.r̅ =  u.v̅
    du.v̅ = -p.μ .* u.r̅ / norm(u.r̅,2)^3
end

function parse_propagation(sol)

    # Parse solution for state variables
    r̅ = map(a->a.r̅, sol.u)
    v̅ = map(a->a.v̅, sol.u)

    rx = map(r->r[1], r̅)
    ry = map(r->r[2], r̅)
    rz = map(r->r[3], r̅)

    vx = map(v->v[1], v̅)
    vy = map(v->v[2], v̅)
    vz = map(v->v[3], v̅)
    
    println(size(norm.(rx + ry + rz).^3))
    ax = -sol.prob.p.μ .* rx ./ (norm.(rx + ry + rz, Ref(2)).^3)
    ay = -sol.prob.p.μ .* ry ./ (norm.(rx + ry + rz, Ref(2)).^3)
    az = -sol.prob.p.μ .* rz ./ (norm.(rx + ry + rz, Ref(2)).^3)

    return hcat(rx,ry,rz), 
           hcat(vx,vy,vz),
           hcat(ax,ay,az)

end