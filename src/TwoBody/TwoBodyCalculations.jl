#
#   TwoBodyCalculations.jl
#
#   Includes simple calculations relevant to the Two Body Problem.

# Constructors

"""
    CartesianState(r̅, v̅, body)

Constructor for Cartesian orbital representation.
"""
function CartesianState(r̅::AbstractArray{Unitful.Length}, 
                        v̅::AbstractArray{Unitful.Velocity}, 
                        body::Body)

    return CartesianState(
            SVector{3,Unitful.Quantity{Float64}}(r̅),
            SVector{3,Unitful.Quantity{Float64}}(v̅),
            body)
end

"""
    KeplerianState(e, a, i, Ω, ω, ν, body)

Constructor for Keplarian orbital representation.
"""
function KeplerianState(e::Unitful.DimensionlessQuantity, 
                        a::Unitful.Length, 
                        i::Unitful.DimensionlessQuantity, 
                        Ω::Unitful.DimensionlessQuantity, 
                        ω::Unitful.DimensionlessQuantity, 
                        ν::Unitful.DimensionlessQuantity, 
                        body::Body)

    return KeplerianState(
            e, a, i, 
            Ω, ω, ν, body)

end

"""
    KeplerianState(orbit::CartesianState)

Returns a Keplarian representation of an orbital state.
"""
function KeplerianState(orbit::CartesianState)

    î = SVector{3, Float64}([1, 0, 0]) 
    ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3, Float64}([0, 0, 1])

    h_vec = specific_angular_momentum_vector(orbit)
    n_vec = mean_motion_vector(orbit)
    e_vec = eccentricity_vector(orbit)

    acos_chop(unit, arg) = arg ≈  1.0 ? acos(unit,  1.0) : 
                                 arg ≈ -1.0 ? acos(unit, -1.0) : 
                                 acos(unit, arg)

    Ω = ustrip(n_vec ⋅ ĵ) > 0 ? 
            acos_chop(u"rad", (n_vec ⋅ î) / norm(n_vec)) :
            2π * u"rad" - acos_chop(u"rad", (n_vec ⋅ î) / norm(n_vec))

    ω = ustrip(e_vec ⋅ k̂) > 0 ?
            acos_chop(u"rad", (n_vec ⋅ e_vec) / (norm(n_vec) * norm(e_vec))) :
            2π * u"rad" - acos_chop(u"rad", (n_vec ⋅ e_vec) / (norm(n_vec) * norm(e_vec)))

    ν = ustrip(orbit.r̅ ⋅  orbit.v̅) > 0 ? 
            acos_chop(u"rad", (e_vec ⋅ orbit.r̅) / (norm(e_vec) * norm(orbit.r̅))) :
            2π * u"rad" - acos_chop(u"rad", (e_vec ⋅ orbit.r̅) / (norm(e_vec) * norm(orbit.r̅)))

    return KeplerianState( 
            norm(e_vec) * u"rad",
            upreferred(semimajor_axis(orbit)),
            acos_chop(u"rad", (h_vec ⋅ k̂) / norm(h_vec)),
            Ω,
            ω,
            ν,
            orbit.body)

end
KeplerianState(orbit::KeplerianState) = orbit

"""
    CartesianState(orbit::KeplerianState)

Returns a Cartesian representation of an orbital state.
"""
function CartesianState(orbit::KeplerianState)

    # Find semilatus parameter
    p = semi_parameter(orbit)

    # Find scalar radius
    r = instantaneous_radius(p, orbit.e, orbit.ν)

    # Set perifocal axes
    P̂=SVector{3, Float64}([1, 0, 0])
    Q̂=SVector{3, Float64}([0, 1, 0]) 
    Ŵ=SVector{3, Float64}([0, 0, 1])

    # Find state in Perifocal frame
    r̅_perifocal = (r * cos(orbit.ν) .* P̂ .+ r * sin(orbit.ν) .* Q̂)
    v̅_perifocal = √(orbit.body.μ/p) * ((-sin(orbit.ν) * P̂) .+ ((orbit.e + cos(orbit.ν)) .* Q̂))

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

    return CartesianState(ᴵTₚ * r̅_perifocal, ᴵTₚ * v̅_perifocal, orbit.body)

end
CartesianState(orbit::CartesianState) = orbit

# Calculation (Helper) Functions

"""
    semimajor_axis(r, v, μ)

Returns semimajor axis parameter, a.
"""
function semimajor_axis(r, v, μ)
   
    return inv( (2 / r) - (v^2 / μ) )

end

"""
    semimajor_axis(orbit::KeplerianState)

Returns semimajor axis parameter, a, for a Keplarian representation.
"""
function semimajor_axis(orbit::KeplerianState)

    return orbit.a

end

"""
    semimajor_axis(orbit::CartesianState)

Returns semimajor axis parameter, a, for a Cartesian representation.
"""
function semimajor_axis(orbit::CartesianState)

    return semimajor_axis(
                norm(orbit.r̅),
                norm(orbit.v̅),
                orbit.body.μ)

end

"""
    specific_angular_momentum_vector(r̅, v̅)

Returns specific angular momentum vector, h̅.
"""
function specific_angular_momentum_vector(r̅, v̅)

    return r̅ × v̅

end

"""
    specific_angular_momentum_vector(orbit::CartesianState)

Returns specific angular momentum vector, h̅, for a Cartesian representation.
"""
function specific_angular_momentum_vector(orbit::CartesianState)

    return specific_angular_momentum_vector(orbit.r̅, orbit.v̅)

end

"""
    specific_angular_momentum(r̅, v̅)

Returns scalar specific angular momentum vector, h.
"""
function specific_angular_momentum(r̅, v̅)

    return norm(specific_angular_momentum_vector(r̅, v̅))

end

"""
    specific_angular_momentum(orbit::CartesianState)

Returns scalar specific angular momentum, h, for a Cartesian representation.
"""
function specific_angular_momentum(orbit::CartesianState)

    return specific_angular_momentum(orbit.r̅, orbit.v̅)

end

"""
    specific_angular_momentum(orbit::TwoBodyOrbit)

Returns scalar specific angular momentum, h, for any orbital representation.
"""
function specific_angular_momentum(orbit::TwoBodyOrbit)

    return apoapsis_radius(orbit) * apoapsis_velocity(orbit)

end

"""
    specific_energy(a, μ)

Returns specific orbital energy, ϵ.
"""
function specific_energy(a, μ)

    return ( -μ / (2 * a) )

end

"""
    specific_energy(r, v, μ)

Returns specific orbital energy, ϵ.
"""
function specific_energy(r, v, μ)

    return (v^2 / 2) - (μ / r)

end

"""
    specific_energy(orbit::TwoBodyOrbit)

Returns specific orbital energy, ϵ, for any orbital representation.
"""
function specific_energy(orbit::TwoBodyOrbit)

    return specific_energy(semimajor_axis(orbit), orbit.body.μ)

end

"""
    eccentricity_vector(r̅, v̅, μ)

Returns orbital eccentricity, e.
"""
function eccentricity_vector(r̅, v̅, μ)

    return (1 / μ) * ((v̅ × specific_angular_momentum_vector(r̅, v̅)) - μ * r̅ / norm(r̅))

end

"""
    eccentricity_vector(orbit::CartesianState)

Returns orbital eccentricity_vector, e̅.
"""
function eccentricity_vector(orbit::CartesianState)

    return eccentricity_vector(orbit.r̅, orbit.v̅, orbit.body.μ)

end

"""
    eccentricity(orbit::CartesianState)

Returns orbital eccentricity, e.
"""
function eccentricity(orbit::CartesianState)

    return norm(eccentricity_vector(orbit.r̅, orbit.v̅, orbit.body.μ))

end

"""
    eccentricity(orbit::KeplerianState)

Returns orbital eccentricity, e.
"""
function eccentricity(orbit::KeplerianState)

    return orbit.e

end

"""
    semi_parameter(a, e)

Returns semilatus parameter, p.
"""
function semi_parameter(a, e)

    return a * (1 - e^2)

end

"""
    semi_parameter(orbit::TwoBodyOrbit

Returns semilatus parameter, p, for any orbital representation.
"""
function semi_parameter(orbit::TwoBodyOrbit)

    return semimajor_axis(orbit) * 
            (1 - eccentricity(orbit)^2)

end

"""
    semi_parameter(orbit::KeplerianState)

Returns semilatus parameter, p, for a Keplarian representation.
"""
function semi_parameter(orbit::KeplerianState)

    return orbit.a * (1 - orbit.e^2)

end

"""
    instantaneous_radius(p, e, ν)

Returns instantaneous radius, r.
"""
function instantaneous_radius(p, e, ν)

    return p / (1 + e * cos(ν))

end

"""
    instantaneous_radius(orbit::TwoBodyOrbit)

Returns instantaneous radius, r, for any orbital representation.
"""
function instantaneous_radius(orbit::TwoBodyOrbit)

    return instantaneous_radius(semi_parameter(orbit),
                                eccentricity(orbit),
                                true_anomoly(orbit))

end

"""
    instantaneous_velocity(r, a, μ)

Returns instantaneous velocity, v, for any orbital representation.
"""
function instantaneous_velocity(r, a, μ)

    return √( (2 * μ / r) - (μ / a))

end

"""
    periapsis_radius(a, e)

Returns periapsis radius, r̅_p.
"""
function periapsis_radius(a, e)

    return a * (1 - e)

end

"""
    periapsis_radius(orbit::TwoBodyOrbit)

Returns periapsis radius, r_p, for any orbital representation.
"""
function periapsis_radius(orbit::TwoBodyOrbit)

    return periapsis_radius(semimajor_axis(orbit), 
                            eccentricity(orbit))

end

"""
    apoapsis_radius(a, e)

Returns periapsis radius, r_a.
"""
function apoapsis_radius(a, e)

    return a * (1 + e)

end

"""
    apoapsis_radius(orbit::TwoBodyOrbit)

Returns periapsis radius, r_a, for any orbital representation.
"""
function apoapsis_radius(orbit::TwoBodyOrbit)

    return apoapsis_radius(semimajor_axis(orbit), 
                           eccentricity(orbit))

end

"""
    periapsis_velocity(orbit::TwoBodyOrbit)

Returns periapsis velocity, v_p, for any orbital representation.
"""
function periapsis_velocity(orbit::TwoBodyOrbit)

    return instantaneous_velocity(
                periapsis_radius(orbit), 
                semimajor_axis(orbit),
                orbit.body.μ)

end

"""
    apoapsis_velocity(orbit::TwoBodyOrbit)

Returns apoapsis velocity, v_a, for any orbital representation.
"""
function apoapsis_velocity(orbit::TwoBodyOrbit)

    return instantaneous_velocity(
                apoapsis_radius(orbit), 
                semimajor_axis(orbit),
                orbit.body.μ)

end

"""
    orbital_period(orbit::TwoBodyOrbit)

Returns orbital period, Τ, for any orbital representation.
"""
function orbital_period(orbit::TwoBodyOrbit)

    P = 2π * u"rad" * √(semimajor_axis(orbit)^3 / orbit.body.μ)

end

"""
    true_anomoly(r, h, e, μ)

Returns true anomoly, ν.
"""
function true_anomoly(r, h, e, μ)

    return acos( (h^2 - μ * r) / (μ * r * e) )

end

"""
    true_anomoly(orbit::CartesianState)

Returns true anomoly, ν, for a Cartesian representation.
"""
function true_anomoly(orbit::CartesianState)

    return true_anomoly(
            norm(orbit.r̅),
            specific_angular_momentum(orbit),
            eccentricity(orbit),
            orbit.body.μ)

end

"""
    true_anomoly(orbit::KeplerianState)

Returns true anomoly, ν, for a Keplarian representation.
"""
function true_anomoly(orbit::KeplerianState)

    return orbit.ν

end

"""
    mean_motion(a, μ)

Returns mean motion, n.
"""
function mean_motion(a, μ)

    return √(μ / a^3)

end

"""
    mean_motion(orbit::TwoBodyOrbit)

Returns mean motion, n, for any orbital representation.
"""
function mean_motion(orbit::TwoBodyOrbit)

return mean_motion(
    semimajor_axis(orbit),
    orbit.body.μ)

end

"""
    mean_motion_vector(orbit::CartesianState)

Returns mean motion vector, n̄, for a Cartesian representation.
"""
function mean_motion_vector(orbit::CartesianState)

#   î = SVector{3, Float64}([1, 0, 0]) 
#   ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3, Float64}([0, 0, 1])

    return k̂ × specific_angular_momentum_vector(orbit)

end

"""
    eccentric_anomoly(e, ν)

Returns eccentric anomoly, E.
"""
function eccentric_anomoly(e, ν)

    return atan( u"rad", (√(1 - e^2) * sin(ν)) / (e + cos(ν)) )

end

"""
    eccentric_anomoly(orbit::TwoBodyOrbit)

Returns eccentric anomoly, E, for any orbital representation.
"""
function eccentric_anomoly(orbit::TwoBodyOrbit)

    return eccentric_anomoly(
            eccentricity(orbit),
            true_anomoly(orbit))

end

"""
    time_since_periapsis(n, e, E)

Returns time since periapsis, t.
"""
function time_since_periapsis(n, e, E)

    return (E - e * sin(E)) / (n)

end

"""
    time_since_periapsis(orbit::TwoBodyOrbit)

Returns time since periapsis, t, for any orbital representation.
"""
function time_since_periapsis(orbit::TwoBodyOrbit)

    return time_since_periapsis(
            mean_motion(orbit),
            eccentricity(orbit),
            eccentric_anomoly(orbit))

end

"""
    inclination(orbit::CartesianState)

Returns orbital inclination, i, for a Cartesian representation.
"""
function inclination(orbit::CartesianState)

    h_vec = specific_angular_momentum_vector(orbit)
    return acos(u"rad", h_vec[3] / norm(h_vec))

end

"""
    inclination(orbit::KeplerianState)

Returns orbital inclination, i, for a Keplarian representation.
"""
function inclination(orbit::KeplerianState)

    return orbit.i

end

function Base.isapprox(c1::CartesianState, c2::CartesianState; atol=1e-8)

    return all(ustrip.(c1.r̅ - c2.r̅) .< atol) &&
           all(ustrip.(c1.r̅ - c2.r̅) .< atol) &&
           (c1.body == c2.body)

end

function Base.isapprox(c1::KeplerianState, c2::KeplerianState; atol=1e-8)

    return ustrip(upreferred(c1.e - c2.e)) .< atol &&
           ustrip(upreferred(c1.a - c2.a)) .< atol &&
           ustrip(upreferred(mod(c1.i, 180u"°") - mod(c2.i, 180u"°"))) .< atol &&
           ustrip(upreferred(mod(c1.Ω, 360u"°") - mod(c2.Ω, 360u"°"))) .< atol &&
           ustrip(upreferred(mod(c1.ω, 360u"°") - mod(c2.ω, 360u"°"))) .< atol &&
           ustrip(upreferred(mod(c1.ν, 360u"°") - mod(c2.ν, 360u"°"))) .< atol &&
           c1.body == c2.body

end

function Base.isequal(c1::CartesianState, c2::CartesianState)

    return all(c1.r̅ .== c2.r̅) &&
           all(c1.v̅ .== c2.v̅) &&
           c1.body == c2.body

end

function Base.isequal(c1::KeplerianState, c2::KeplerianState)

    return c1.e == c2.e &&
           c1.a == c2.a &&
           mod(c1.i, 180u"°") == mod(c2.i, 180u"°") &&
           mod(c1.Ω, 360u"°") == mod(c2.Ω, 360u"°") &&
           mod(c1.ω, 360u"°") == mod(c2.ω, 360u"°") &&
           mod(c1.ν, 360u"°") == mod(c2.ν, 360u"°") &&
           c1.body == c2.body

end