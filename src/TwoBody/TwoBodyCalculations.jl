#
#   TwoBodyCalculations.jl
#
#   Includes simple calculations relevant to the Two Body Problem.

"""
Returns the conic section, as specified by eccentricity `e`.
"""
function conic(e::T) where T<:Number

    if e ≈ 0
        return Circular
    elseif e ≈ 1
        return Parabolic
    elseif 0 < e && e < 1
        return Elliptical
    elseif e > 1
        return Hyperbolic
    else
        return Invalid
    end

end
conic(orbit::T) where T<:TwoBodySystem = conic(eccentricity(orbit))

"""
Returns a Keplarian representation of a Cartesian orbital state.
Algorithm taught in ENAE601.
"""
function keplerian(rᵢ, vᵢ, μ)

    safe_acos(unit, num) = isapprox(num, 1) ? acos(one(num)) * unit : acos(num) * unit

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
keplerian(rᵢ, vᵢ, body::CelestialBody) = keplerian(rᵢ, vᵢ, body.μ)
keplerian(orbit::KeplerianState) = orbit.e, orbit.a, orbit.i, orbit.Ω, orbit.ω, orbit.ν
keplerian(orbit::TwoBodyState) = keplerian(orbit.r, orbit.v, orbit.body)

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
cartesian(e, a, i, Ω, ω, ν, body::CelestialBody) = cartesian(e, a, i, Ω, ω, ν, body.μ)
cartesian(orbit::TwoBodyState) = orbit.r, orbit.v
cartesian(orbit::KeplerianState) = cartesian(orbit.e, orbit.a, orbit.i, orbit.Ω, orbit.ω, orbit.ν, orbit.body)

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
        r = radius(p, e, ν)
        
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

function perifocal(orbit::T) where T <: TwoBodySystem
    return perifocal(
        inclination(orbit),
        RAAN(orbit),
        argument_of_periapsis(orbit),
        radius_vector(orbit),
        velocity_vector(orbit)
    )
end

"""
Returns semimajor axis parameter, a.
"""
semimajor_axis(r, v, μ) = inv( (2 / r) - (v^2 / μ) )
semimajor_axis(orbit::TwoBodyState) = semimajor_axis(radius(orbit), velocity(orbit), orbit.body.μ)
semimajor_axis(orbit::KeplerianState) = orbit.a

"""
Returns specific angular momentum vector, h̅.
"""
specific_angular_momentum_vector(rᵢ, vᵢ) = rᵢ × vᵢ
specific_angular_momentum_vector(orbit::TwoBodyState) = specific_angular_momentum_vector(orbit.r, orbit.v)
specific_angular_momentum_vector(orbit::KeplerianState) = specific_angular_momentum_vector(TwoBodyState(orbit))

"""
Returns scalar specific angular momentum vector, h.
"""
specific_angular_momentum(rᵢ, vᵢ) = norm(specific_angular_momentum_vector(rᵢ, vᵢ))
specific_angular_momentum(orbit::TwoBodyState) = specific_angular_momentum(orbit.r, orbit.v)
specific_angular_momentum(orbit::KeplerianState) = specific_angular_momentum(cartesian(orbit)...)

"""
Returns specific orbital energy, ϵ.
"""
specific_energy(a, μ) = ( -μ / (2 * a) )
specific_energy(r, v, μ) = (v^2 / 2) - (μ / r)
specific_energy(orbit::TwoBodyState) = specific_energy(orbit.r, orbit.v, orbit.body.μ)
specific_energy(orbit::KeplerianState) = specific_energy(orbit.a, orbit.body.μ)

"""
Returns potential energy for an orbit about a `CelestialBody`.
"""
specific_potential_energy(r, μ) = (μ/r)
specific_potential_energy(r, μ, R, J₂, ϕ) = (μ/r) * (1 - J₂ * (R/r)^2 * ((3/2) * (sin(ϕ))^2 - (1/2)))
specific_potential_energy(orbit::TwoBodyState) = specific_potential_energy(orbit.r, orbit.body.μ)
specific_potential_energy(orbit::KeplerianState) = specific_potential_energy(TwoBodyState(orbit))

"""
Returns orbital eccentricity vector e̅.
"""
function eccentricity_vector(rᵢ, vᵢ, μ)

    return map(x-> abs(x) < eps(typeof(x)) ? zero(x) : x, (1 / μ) * ((vᵢ × specific_angular_momentum_vector(rᵢ, vᵢ)) - μ * rᵢ / norm(rᵢ)))

end
eccentricity_vector(orbit::TwoBodyState) = eccentricity_vector(orbit.r, orbit.v, orbit.body.μ)
eccentricity_vector(orbit::KeplerianState) = eccentricity_vector(TwoBodyState(orbit))

"""
Returns orbital eccentricity, e.
"""
eccentricity(rᵢ, vᵢ, μ) = norm(eccentricity_vector(rᵢ, vᵢ, μ)) |> upreferred
eccentricity(orbit::TwoBodyState) = eccentricity(orbit.r, orbit.v, orbit.body.μ)
eccentricity(orbit::KeplerianState) = orbit.e

"""
Returns semilatus parameter, p.
"""
semi_parameter(a, e) = a * (1 - e^2)
semi_parameter(orbit::KeplerianState) = semi_parameter(orbit.a, orbit.e)
semi_parameter(orbit::TwoBodyState) = semi_parameter(KeplerianState(orbit))

"""
Returns radius, r.
"""
radius(p, e, ν) = upreferred(p / (1 + e * cos(ν)))
radius(orbit::KeplerianState) = radius(semi_parameter(orbit), orbit.e, orbit.ν)
radius(orbit::TwoBodyState) = norm(orbit.r)
radius(body::CelestialBody) = body.R

"""
Returns radius vector, rᵢ.
"""
radius_vector(orbit::TwoBodyState) = orbit.r
radius_vector(orbit::KeplerianState) = radius_vector(TwoBodyState(orbit))

"""
Returns instantaneous velocity, v, for any orbital representation.
"""
velocity(r, a, μ) =  upreferred(√( (2 * μ / r) - (μ / a)))
velocity(orbit::KeplerianState) = velocity(radius(orbit), orbit.a, orbit.body.μ)
velocity(orbit::TwoBodyState) = norm(orbit.v)

"""
Returns velocity vector, vᵢ.
"""
velocity_vector(orbit::TwoBodyState) = orbit.v
velocity_vector(orbit::KeplerianState) = velocity_vector(TwoBodyState(orbit))

"""
Returns periapsis radius, rₚ.
"""
periapsis_radius(a, e) = a * (1 - e)
periapsis_radius(orbit::T) where T<:TwoBodySystem = periapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

"""
Returns apoapsis radius, rₐ.
"""
apoapsis_radius(a, e) = a * (1 + e)
apoapsis_radius(orbit::T) where T<:TwoBodySystem = apoapsis_radius(semimajor_axis(orbit), eccentricity(orbit))

"""
Returns periapsis velocity, vₚ, for any orbital representation.
"""
periapsis_velocity(orbit::T) where T<:TwoBodySystem = velocity(periapsis_radius(orbit), semimajor_axis(orbit), orbit.body.μ)


"""
Returns apoapsis velocity, v_a, for any orbital representation.
"""
apoapsis_velocity(orbit::T) where T<:TwoBodySystem = velocity(apoapsis_radius(orbit), semimajor_axis(orbit), orbit.body.μ)

"""
Returns mass `m`.
"""
mass(body::CelestialBody) = body.μ / G

"""
Returns mass parameter `μ`.
"""
mass_parameter(body::CelestialBody) = body.μ

"""
Returns the orbital period.
"""
period(a, μ) = 2π * √(upreferred(a^3 / μ))
period(orbit::T) where T<:TwoBodySystem = period(semimajor_axis(orbit), orbit.body.μ)

"""
Returns true anomoly, ν.
"""
function true_anomoly(r, h, e, μ)
    val = (h^2 - μ * r) / (μ * r * e)
    acos(u"rad", isapprox(val, one(val)) ? one(val) : val)
end
true_anomoly(orbit::KeplerianState) = orbit.ν
true_anomoly(orbit::TwoBodyState) = true_anomoly(radius(orbit), specific_angular_momentum(orbit), eccentricity(orbit), orbit.body.μ)

"""
Returns mean motion, n.
"""
mean_motion(a, μ) = √(μ / a^3)
mean_motion(orbit::T) where T<:TwoBodySystem = mean_motion(semimajor_axis(orbit), orbit.μ)

"""
Returns mean motion vector, n̄.
"""
function mean_motion_vector(orbit::T) where T<:TwoBodySystem
#   î = SVector{3, Float64}([1, 0, 0]) 
#   ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3, Float64}([0, 0, 1])
    return k̂ × specific_angular_momentum_vector(orbit)
end

"""
Returns eccentric anomoly, E, parabolic anomoly, B, or hyperbolic 
anomoly, H. 
"""
function eccentric_anomoly(orbit::T) where T <: TwoBodySystem
    e = eccentricity(orbit)
    ν = true_anomoly(orbit)
    return acos(u"rad", (e + cos(ν) / (1 + e * cos(ν))))
end

"""
Returns time since periapsis, t.
"""
time_since_periapsis(n, e, E) = (E - e * sin(E)) / (n)
time_since_periapsis(orbit::T) where T <: TwoBodySystem = time_since_periapsis(mean_motion(orbit), eccentricity(orbit), eccentric_anomoly(orbit))

"""
Returns orbital inclination, i.
"""
inclination(orbit::KeplerianState) =  orbit.i
inclination(orbit::TwoBodyState) = inclination(KeplerianState(orbit))

"""
Returns Right Ascension of the Ascending Node, Ω.
"""
RAAN(orbit::KeplerianState) = orbit.Ω
RAAN(orbit::TwoBodyState) = RAAN(KeplerianState(orbit))

"""
Returns the Argument of Periapsis, ω.
"""
argument_of_periapsis(orbit::KeplerianState) = orbit.ω
argument_of_periapsis(orbit::TwoBodyState) = argument_of_periapsis(KeplerianState(orbit))

"""
Returns true if all elements in each system are within `atol` of the other.
"""
function Base.isapprox(c1::TwoBodySystem, c2::TwoBodySystem; atol=1e-6)

    return all(ustrip.(abs.(radius_vector(c1) - radius_vector(c2))) .< atol) &&
           all(ustrip.(abs.(velocity_vector(c1) - velocity_vector(c2))) .< atol) &&
           ustrip(upreferred(abs(eccentricity(c1) - eccentricity(c2)))) < atol &&
           ustrip(upreferred(abs(semimajor_axis(c1) - semimajor_axis(c2)))) < atol &&
           ustrip(upreferred(abs(mod(inclination(c1), 180u"°") - mod(inclination(c2), 180u"°")))) < atol &&
           ustrip(upreferred(abs(mod(RAAN(c1), 360u"°") - mod(RAAN(c2), 360u"°")))) < atol &&
           ustrip(upreferred(abs(mod(argument_of_periapsis(c1), 360u"°") - mod(argument_of_periapsis(c2), 360u"°")))) < atol &&
           ustrip(upreferred(abs(mod(true_anomoly(c1), 360u"°") - mod(true_anomoly(c2), 360u"°")))) < atol &&
           isapprox(c1.body, c2.body; atol=atol)

end

"""
Returns true if all elements of each system are identically equal.
"""
function Base.isequal(c1::TwoBodySystem, c2::TwoBodySystem)

    return all(radius_vector(c1) .== radius_vector(c2)) &&
           all(velocity_vector(c1) .== velocity_vector(c2)) &&
           eccentricity(c1) == eccentricity(c2) &&
           semimajor_axis(c1) == semimajor_axis(c2) &&
           mod(inclination(c1), 180u"°") == mod(inclination(c2), 180u"°") &&
           mod(RAAN(c1), 360u"°") == mod(RAAN(c2), 360u"°") &&
           mod(argument_of_periapsis(c1), 360u"°") == mod(argument_of_periapsis(c2), 360u"°") &&
           mod(true_anomoly(c1), 360u"°") == mod(true_anomoly(c2), 360u"°") &&
           c1.body == c2.body

end

"""
Returns true if all elements are within `atol` of the other.
"""
function Base.isapprox(b1::CelestialBody, b2::CelestialBody; atol=1e-6)
    return ustrip(upreferred(abs(b1.R - b2.R))) < atol &&
           ustrip(upreferred(abs(b1.μ - b2.μ))) < atol
end

"""
Returns true if all elements are identically equal.
"""
function Base.isequal(b1::CelestialBody, b2::CelestialBody)
    return b1.R == b2.R && b1.μ == b2.μ
end