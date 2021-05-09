#
#   TwoBodyCalculations.jl
#
#   Includes simple calculations relevant to the Two Body Problem.

"""
Returns the conic section, as specified by eccentricity `e`.
"""
function conic(e::T) where T<:Real

    if e ≈ zero(T)
        return Circular
    elseif e ≈ one(T)
        return Parabolic
    elseif zero(T) < e && e < one(T)
        return Elliptical
    elseif e > one(T)
        return Hyperbolic
    else
        return Invalid
    end

end
conic(orbit::T) where T<:RestrictedTwoBodyOrbit = conic(eccentricity(orbit))

"""
Returns a Keplarian representation of a Cartesian orbital state.
Algorithm taught in ENAE601.
"""
function keplerian(rᵢ, vᵢ, μ)

    safe_acos(unit, num) = isapprox(num, one(num)) ? 
                            acos(unit, one(num)) : 
                                isapprox(num, -one(num)) ? 
                                    acos(unit, -one(num)) : 
                                        acos(unit, num)

    î = SVector{3, Float64}([1, 0, 0]) 
    ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3, Float64}([0, 0, 1])

    h̅ = specific_angular_momentum_vector(rᵢ, vᵢ)
    a = semimajor_axis(norm(rᵢ), norm(vᵢ), μ)
    n̅ = k̂ × specific_angular_momentum_vector(rᵢ, vᵢ)
    e̅ = eccentricity_vector(rᵢ, vᵢ, μ)
    e = norm(e̅) |> upreferred

    i = safe_acos(u"rad", (h̅ ⋅ k̂) / norm(h̅))

    Ω = ustrip(n̅ ⋅ ĵ) > 0 ? 
            safe_acos(u"rad", (n̅ ⋅ î) / norm(n̅)) :
            2π * u"rad" - safe_acos(u"rad", (n̅ ⋅ î) / norm(n̅))

    ω = ustrip(e̅ ⋅ k̂) > 0 ?
            safe_acos(u"rad", (n̅ ⋅ e̅) / (norm(n̅) * e)) :
            2π * u"rad" - safe_acos(u"rad", (n̅ ⋅ e̅) / (norm(n̅) * e))

    ν = ustrip(rᵢ ⋅  vᵢ) > 0 ? 
            safe_acos(u"rad", (e̅ ⋅ rᵢ) / (e * norm(rᵢ))) :
            2π * u"rad" - safe_acos(u"rad", (e̅ ⋅ rᵢ) / (e * norm(rᵢ)))

    return e, uconvert(u"km", a), uconvert(u"°", i), 
           uconvert(u"°", Ω), uconvert(u"°", ω), 
           uconvert(u"°", ν)

end
keplerian(rᵢ, vᵢ, body::RestrictedTwoBodySystem) = keplerian(rᵢ, vᵢ, mass_parameter(body))
keplerian(orbit::KeplerianOrbit) = eccentricity(orbit.state), semimajor_axis(orbit.state), inclination(orbit.state), RAAN(orbit.state), argument_of_periapsis(orbit.state), true_anomoly(orbit.state)
keplerian(orbit::CartesianOrbit) = keplerian(position_vector(orbit.state), velocity_vector(orbit.state), orbit.system)

"""
Returns a Cartesian representation of a Keplerian two-body orbital state
in an inertial frame, centered at the center of mass of the central body.
Algorithm taught in ENAE601.
"""
function cartesian(e, a, i, Ω, ω, ν, μ)
    rᵢ, vᵢ = cartesian(i, Ω, ω, perifocal(a, e, ν, μ)...)
    return  uconvert.(u"km",    rᵢ), 
            uconvert.(u"km/s",  vᵢ)
end
cartesian(e, a, i, Ω, ω, ν, body::RestrictedTwoBodySystem) = cartesian(e, a, i, Ω, ω, ν, mass_parameter(body))
cartesian(orbit::CartesianOrbit) = position_vector(orbit.state), velocity_vector(orbit.state)
cartesian(orbit::KeplerianOrbit) = cartesian(eccentricity(orbit.state), semimajor_axis(orbit.state), inclination(orbit.state), RAAN(orbit.state), argument_of_periapsis(orbit.state), true_anomoly(orbit.state), orbit.body)

"""
Returns a Cartesian (inertial) representation of the provied Perifocal state.
"""
function cartesian(i, Ω, ω, rₚ, vₚ)

    # Set up Perifocal ⟶ Cartesian conversion
    R_3Ω =  SMatrix{3,3}(
            [cos(Ω)          -sin(Ω)            0.;
             sin(Ω)           cos(Ω)            0.;
             0.               0.                1.])
    R_1i = SMatrix{3,3}(
            [1.               0.                0.;
             0.               cos(i)           -sin(i);
             0.               sin(i)            cos(i)])
    R_3ω = SMatrix{3,3}(
            [cos(ω)          -sin(ω)            0.
             sin(ω)           cos(ω)            0.
             0.               0.                1.])

    ᴵTₚ = R_3Ω * R_1i * R_3ω

    return ᴵTₚ * rₚ, ᴵTₚ * vₚ

end

"""
Returns position and velocity vectors in the Perifocal frame.
"""
function perifocal(a, e, ν, μ)

        p = semi_parameter(a, e)
        r = scalar_position(p, e, ν)
        
        P̂=SVector{3, Float64}([1, 0, 0])
        Q̂=SVector{3, Float64}([0, 1, 0])
        Ŵ=SVector{3, Float64}([0, 0, 1])
        
        rₚ = (r * cos(ν) .* P̂ .+ r * sin(ν) .* Q̂)
        vₚ = √(μ/p) * ((-sin(ν) * P̂) .+ ((e + cos(ν)) .* Q̂))
        
        return rₚ, vₚ

end

function perifocal(i, Ω, ω, rᵢ, vᵢ)

    # Set up Cartesian ⟶ Perifocal conversion
    R_3Ω =  SMatrix{3,3}(
            [cos(Ω)          -sin(Ω)            0.;
             sin(Ω)           cos(Ω)            0.;
             0.               0.                1.])
    R_1i = SMatrix{3,3}(
            [1.               0.                0.;
             0.               cos(i)           -sin(i);
             0.               sin(i)            cos(i)])
    R_3ω = SMatrix{3,3}(
            [cos(ω)          -sin(ω)            0.
             sin(ω)           cos(ω)            0.
             0.               0.                1.])

    ᵖTᵢ = inv(R_3Ω * R_1i * R_3ω)

    return ᵖTᵢ*rᵢ, ᵖTᵢ*vᵢ

end

function perifocal(orbit::T) where T <: RestrictedTwoBodyOrbit
    return perifocal(
        inclination(orbit),
        RAAN(orbit),
        argument_of_periapsis(orbit),
        position_vector(orbit),
        velocity_vector(orbit)
    )
end

"""
Returns semimajor axis parameter, a.
"""
semimajor_axis(r, v, μ) = inv( (2 / r) - (v^2 / μ) )
semimajor_axis(orbit::CartesianOrbit) = semimajor_axis(scalar_position(orbit), scalar_velocity(orbit), mass_parameter(orbit.system)) 
semimajor_axis(orbit::KeplerianOrbit) = orbit.state.a * lengthunit(orbit.state)

"""
Returns specific angular momentum vector, h̅.
"""
specific_angular_momentum_vector(rᵢ, vᵢ) = rᵢ × vᵢ
specific_angular_momentum_vector(orbit::CartesianOrbit) = specific_angular_momentum_vector(position_vector(orbit.state), velocity_vector(orbit.state))
specific_angular_momentum_vector(orbit::KeplerianOrbit) = specific_angular_momentum_vector(CartesianOrbit(orbit))

"""
Returns scalar specific angular momentum vector, h.
"""
specific_angular_momentum(rᵢ, vᵢ) = norm(specific_angular_momentum_vector(rᵢ, vᵢ))
specific_angular_momentum(orbit::CartesianOrbit) = specific_angular_momentum(position_vector(orbit.state), velocity_vector(orbit.state))
specific_angular_momentum(orbit::KeplerianOrbit) = specific_angular_momentum(cartesian(orbit)...)

"""
Returns specific orbital energy, ϵ.
"""
specific_energy(a, μ) = ( -μ / (2 * a) )
specific_energy(r, v, μ) = (v^2 / 2) - (μ / r)
specific_energy(orbit::CartesianOrbit) = specific_energy(scalar_position(orbit.state), scalar_velocity(orbit.state), mass_parameter(orbit.system))
specific_energy(orbit::KeplerianOrbit) = specific_energy(semimajor_axis(orbit.state), mass_parameter(orbit.system))

"""
Returns C3 value.
"""
C3(r, v, μ) = v^2 - 2μ/r
C3(orbit::RestrictedTwoBodyOrbit) = C3(scalar_position(orbit), scalar_velocity(orbit), mass_parameter(orbit.system))

"""
Returns v∞.
"""
v_infinity(r, v, μ) = √C3(r, v, μ)
v_infinity(orbit::RestrictedTwoBodyOrbit) = v_infinity(scalar_position(orbit), scalar_velocity(orbit), mass_parameter(orbit.system))

"""
Returns potential energy for an orbit about a `RestrictedTwoBodySystem`.
"""
specific_potential_energy(r, μ) = (μ/r)
specific_potential_energy(r, μ, R, J₂, ϕ) = (μ/r) * (1 - J₂ * (R/r)^2 * ((3/2) * (sin(ϕ))^2 - (1/2)))
specific_potential_energy(orbit::CartesianOrbit) = specific_potential_energy(scalar_position(orbit.state), mass_parameter(orbit.system))
specific_potential_energy(orbit::KeplerianOrbit) = specific_potential_energy(CartesianOrbit(orbit))

"""
Returns orbital eccentricity vector e̅.
"""
function eccentricity_vector(rᵢ, vᵢ, μ)

    return map(x-> abs(x) < eps(typeof(x)) ? zero(x) : x, (1 / μ) * ((vᵢ × specific_angular_momentum_vector(rᵢ, vᵢ)) - μ * rᵢ / norm(rᵢ)))

end
eccentricity_vector(orbit::CartesianOrbit) = eccentricity_vector(position_vector(orbit.state), velocity_vector(orbit.state), mass_parameter(orbit.system))
eccentricity_vector(orbit::KeplerianOrbit) = eccentricity_vector(CartesianOrbit(orbit))

"""
Returns orbital eccentricity, e.
"""
eccentricity(rᵢ, vᵢ, μ) = norm(eccentricity_vector(rᵢ, vᵢ, μ)) |> upreferred
eccentricity(orbit::CartesianOrbit) = eccentricity(position_vector(orbit.state), velocity_vector(orbit.state), mass_parameter(orbit.system))
eccentricity(orbit::KeplerianOrbit) = orbit.state.e

"""
Returns semilatus parameter, p.
"""
semi_parameter(a, e) = a * (1 - e^2)
semi_parameter(orbit::RestrictedTwoBodyOrbit) = semi_parameter(semimajor_axis(orbit), eccentricity(orbit))

"""
Returns scalar_position, r.
"""
scalar_position(p, e, ν) = upreferred(p / (1 + e * cos(ν)))
scalar_position(orbit::KeplerianOrbit) = scalar_position(semi_parameter(orbit), eccentricity(orbit), true_anomoly(orbit))
scalar_position(orbit::CartesianOrbit) = norm(position_vector(orbit))

"""
Returns position vector, rᵢ.
"""
position_vector(orbit::CartesianOrbit) = position_vector(orbit.state)
position_vector(orbit::KeplerianOrbit) = position_vector(CartesianOrbit(orbit))

"""
Returns instantaneous velocity, v, for any orbital representation.
"""
scalar_velocity(r, a, μ) =  upreferred(√( (2 * μ / r) - (μ / a)))
scalar_velocity(orbit::KeplerianOrbit) = scalar_velocity(scalar_position(orbit), semimajor_axis(orbit), mass_parameter(orbit.system))
scalar_velocity(orbit::CartesianOrbit) = norm(velocity_vector(orbit))

"""
Returns velocity vector, vᵢ.
"""
velocity_vector(orbit::CartesianOrbit) = velocity_vector(orbit.state)
velocity_vector(orbit::KeplerianOrbit) = velocity_vector(CartesianOrbit(orbit))

"""
Returns periapsis scalar_position, rₚ.
"""
periapsis_radius(a, e) = a * (1 - e)
periapsis_radius(orbit::T) where T<:RestrictedTwoBodyOrbit = periapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

"""
Returns apoapsis scalar_position, rₐ.
"""
apoapsis_radius(a, e) = a * (1 + e)
apoapsis_radius(orbit::T) where T<:RestrictedTwoBodyOrbit = apoapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

"""
Returns periapsis velocity, vₚ, for any orbital representation.
"""
periapsis_velocity(orbit::T) where T<:RestrictedTwoBodyOrbit = scalar_velocity(periapsis_radius(orbit), semimajor_axis(orbit), mass_parameter(orbit.system))


"""
Returns apoapsis velocity, v_a, for any orbital representation.
"""
apoapsis_velocity(orbit::T) where T<:RestrictedTwoBodyOrbit = scalar_velocity(apoapsis_radius(orbit), semimajor_axis(orbit), mass_parameter(orbit.system))

"""
Returns mass `m`.
"""
mass(body::RestrictedTwoBodySystem) = mass_parameter(body) / G

"""
Returns mass parameter `μ`.
"""
mass_parameter(body::RestrictedTwoBodySystem) = body.μ * massparameterunit(body)

"""
Returns the orbital period.
"""
period(a, μ) = 2π * √(upreferred(a^3 / μ))

"""
Returns the orbital period.
"""
period(orbit::RestrictedTwoBodyOrbit{Parabolic}) = Inf * timeunit(orbit.state)

"""
Returns the orbital period.
"""
period(orbit::RestrictedTwoBodyOrbit{Hyperbolic}) = NaN * timeunit(orbit.state)

"""
Returns the orbital period.
"""
period(orbit::RestrictedTwoBodyOrbit{Invalid}) = NaN * timeunit(orbit.state)

"""
Returns the orbital period.
"""
period(orbit::RestrictedTwoBodyOrbit{C}) where C<:Union{Elliptical, Circular} = period(semimajor_axis(orbit), mass_parameter(orbit.system))

"""
Returns true anomoly, ν.
"""
function true_anomoly(r, h, e, μ)
    val = (h^2 - μ * r) / (μ * r * e)
    acos(u"rad", isapprox(val, one(val)) ? one(val) : val)
end

"""
Returns true anomoly, ν.
"""
true_anomoly(orbit::KeplerianOrbit) = orbit.state.ν * angularunit(orbit.state)

"""
Returns true anomoly, ν.
"""
true_anomoly(orbit::CartesianOrbit) = true_anomoly(scalar_position(orbit), specific_angular_momentum(orbit), eccentricity(orbit), mass_parameter(orbit.system))

"""
Returns mean motion, n.
"""
mean_motion(a, μ) = √(μ / a^3)
mean_motion(orbit::T) where T<:RestrictedTwoBodyOrbit = mean_motion(semimajor_axis(orbit), mass_parameter(orbit.system))

"""
Returns mean motion vector, n̄.
"""
function mean_motion_vector(orbit::T) where T<:RestrictedTwoBodyOrbit
#   î = SVector{3, Float64}([1, 0, 0]) 
#   ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3, Float64}([0, 0, 1])
    return k̂ × specific_angular_momentum_vector(orbit)
end

"""
Returns eccentric anomoly, E, parabolic anomoly, B, or hyperbolic 
anomoly, H. 
"""
function eccentric_anomoly(orbit::T) where T <: RestrictedTwoBodyOrbit
    e = eccentricity(orbit)
    ν = true_anomoly(orbit)
    return acos(u"rad", (e + cos(ν) / (1 + e * cos(ν))))
end

"""
Returns time since periapsis, t.
"""
time_since_periapsis(n, e, E) = (E - e * sin(E)) / (n)

"""
Returns time since periapsis, t.
"""
time_since_periapsis(orbit::T) where T <: RestrictedTwoBodyOrbit = time_since_periapsis(mean_motion(orbit), eccentricity(orbit), eccentric_anomoly(orbit))

"""
Returns orbital inclination, i.
"""
inclination(orbit::KeplerianOrbit) =  orbit.state.i * angularunit(orbit.state)

"""
Returns orbital inclination, i.
"""
inclination(orbit::CartesianOrbit) = inclination(KeplerianOrbit(orbit))

"""
Returns Right Ascension of the Ascending Node, Ω.
"""
RAAN(orbit::KeplerianOrbit) = orbit.state.Ω * angularunit(orbit.state)

"""
Returns Right Ascension of the Ascending Node, Ω.
"""
RAAN(orbit::CartesianOrbit) = RAAN(KeplerianOrbit(orbit))

"""
Returns the Argument of Periapsis, ω.
"""
argument_of_periapsis(orbit::KeplerianOrbit) = orbit.state.ω * angularunit(orbit.state)

"""
Returns the Argument of Periapsis, ω.
"""
argument_of_periapsis(orbit::CartesianOrbit) = argument_of_periapsis(KeplerianOrbit(orbit))

"""
Returns true if all elements in each system are within `atol` of the other.
"""
function Base.isapprox(c1::RestrictedTwoBodyOrbit, c2::RestrictedTwoBodyOrbit; atol=1e-6)

    return all(ustrip.(abs.(position_vector(c1) - position_vector(c2))) .< atol) &&
           all(ustrip.(abs.(velocity_vector(c1) - velocity_vector(c2))) .< atol) &&
           ustrip(upreferred(abs(eccentricity(c1) - eccentricity(c2)))) < atol &&
           ustrip(upreferred(abs(semimajor_axis(c1) - semimajor_axis(c2)))) < atol &&
           ustrip(upreferred(abs(mod(inclination(c1), 180u"°") - mod(inclination(c2), 180u"°")))) < atol &&
           ustrip(upreferred(abs(mod(RAAN(c1), 360u"°") - mod(RAAN(c2), 360u"°")))) < atol &&
           ustrip(upreferred(abs(mod(argument_of_periapsis(c1), 360u"°") - mod(argument_of_periapsis(c2), 360u"°")))) < atol &&
           ustrip(upreferred(abs(mod(true_anomoly(c1), 360u"°") - mod(true_anomoly(c2), 360u"°")))) < atol &&
           isapprox(c1.system, c2.system; atol=atol)

end

"""
Returns true if all elements of each system are identically equal.
"""
function Base.isequal(c1::RestrictedTwoBodyOrbit, c2::RestrictedTwoBodyOrbit)

    return all(position_vector(c1) .== position_vector(c2)) &&
           all(velocity_vector(c1) .== velocity_vector(c2)) &&
           eccentricity(c1) == eccentricity(c2) &&
           semimajor_axis(c1) == semimajor_axis(c2) &&
           mod(inclination(c1), 180u"°") == mod(inclination(c2), 180u"°") &&
           mod(RAAN(c1), 360u"°") == mod(RAAN(c2), 360u"°") &&
           mod(argument_of_periapsis(c1), 360u"°") == mod(argument_of_periapsis(c2), 360u"°") &&
           mod(true_anomoly(c1), 360u"°") == mod(true_anomoly(c2), 360u"°") &&
           c1.system == c2.system
           
end

"""
Returns true if all elements are within `atol` of the other.
"""
function Base.isapprox(b1::RestrictedTwoBodySystem, b2::RestrictedTwoBodySystem; atol=1e-6)
    return ustrip(upreferred(abs(mass_parameter(b1) - mass_parameter(b2)))) < atol
end

"""
Returns true if all elements are identically equal.
"""
function Base.isequal(b1::RestrictedTwoBodySystem, b2::RestrictedTwoBodySystem)
    return mass_parameter(b1) == mass_parameter(b2)
end

"""
Sphere of influence.
"""
SOI(a, m, M) = a * (m / M)^(2/5)

"""
Sphere of activity.
"""
SOA(a, m, M) = a * (m / 3M)^(1/3)