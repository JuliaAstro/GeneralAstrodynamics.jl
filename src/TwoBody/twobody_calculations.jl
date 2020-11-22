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
conic(orbit::Orbit) = conic(eccentricity(orbit))

"""
Construct `Orbit` from Cartesian elements (in the inertial frame).
"""
function Orbit(rᵢ, vᵢ, body)

    e, a, i, Ω, ω, ν = keplerian(rᵢ, vᵢ, body)
    rₚ, vₚ = perifocal(a, e, ν, body.μ)

    return Orbit{
            conic(eccentricity(rᵢ, vᵢ, body.μ)),
            Float64
        }(SVector{3}(Float64.(rᵢ)), 
          SVector{3}(Float64.(vᵢ)), 
          SVector{3}(Float64.(rₚ)),
          SVector{3}(Float64.(vₚ)),
          map(Float64, [e,a,i,Ω,ω,ν])...,
          body)

end

"""
Construct `Orbit` from Keplerian elements.
"""
function Orbit(e, a, i, Ω, ω, ν, body)

    rᵢ, vᵢ = cartesian(e, a, i, Ω, ω, ν, body)
    rₚ, vₚ = perifocal(a, e, ν, body.μ)
    
    return Orbit{
            conic(e),
            Float64
        }(SVector{3}(Float64.(rᵢ)), 
          SVector{3}(Float64.(vᵢ)), 
          SVector{3}(Float64.(rₚ)),
          SVector{3}(Float64.(vₚ)),
          map(Float64, [e,a,i,Ω,ω,ν])...,
          body)

end

"""
Returns a Keplarian representation of a Cartesian orbital state.
Algorithm taught in ENAE601.
"""
function keplerian(rᵢ, vᵢ, μ)

    î = SVector{3, Float64}([1, 0, 0]) 
    ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3, Float64}([0, 0, 1])

    h̅ = specific_angular_momentum_vector(rᵢ, vᵢ)
    a = semimajor_axis(norm(rᵢ), norm(vᵢ), μ)
    n̅ = k̂ × specific_angular_momentum_vector(rᵢ, vᵢ)
    e̅ = eccentricity_vector(rᵢ, vᵢ, μ)
    e = norm(e̅)
    i = acos(u"rad", (h̅ ⋅ k̂) / norm(h̅))

    Ω = ustrip(n̅ ⋅ ĵ) > 0 ? 
            acos(u"rad", (n̅ ⋅ î) / norm(n̅)) :
            2π * u"rad" - acos(u"rad", (n̅ ⋅ î) / norm(n̅))

    ω = ustrip(e̅ ⋅ k̂) > 0 ?
            acos(u"rad", (n̅ ⋅ e̅) / (norm(n̅) * e)) :
            2π * u"rad" - acos(u"rad", (n̅ ⋅ e̅) / (norm(n̅) * e))

    ν = ustrip(rᵢ ⋅  vᵢ) > 0 ? 
            acos(u"rad", (e̅ ⋅ rᵢ) / (e * norm(rᵢ))) :
            2π * u"rad" - acos(u"rad", (e̅ ⋅ rᵢ) / (e * norm(rᵢ)))

    return e, uconvert(u"km", a), uconvert(u"°", i), 
           uconvert(u"°", Ω), uconvert(u"°", ω), 
           uconvert(u"°", ν)

end
keplerian(rᵢ, vᵢ, body::CelestialBody) = keplerian(rᵢ, vᵢ, body.μ)
keplerian(orbit::Orbit) = orbit.e, orbit.a, orbit.i, 
                                 orbit.Ω, orbit.ω, orbit.ν

"""
Returns a Cartesian representation of a Keplerian two-body orbital state
in an inertial frame, centered at the center of mass of the central body.
Algorithm taught in ENAE601.
"""
function cartesian(e, a, i, Ω, ω, ν, μ)

    rᵢ, vᵢ = inertial(i, Ω, ω, perifocal(a, e, ν, μ)...)
    return  uconvert.(u"km",    rᵢ), 
            uconvert.(u"km/s",  vᵢ)

end
cartesian(e, a, i, Ω, ω, ν, body::CelestialBody) = cartesian(e, a, i, Ω, ω, ν, body.μ)
cartesian(orbit::Orbit) = orbit.rᵢ, orbit.vᵢ

"""
Transforms 3-vector from Perifocal frame to Cartesian space (x,y,z).
"""
function inertial(i, Ω, ω, rₚ, vₚ)

    # Set up Perifocal ⟶ Cartesian conversion
    R_3Ω =  SMatrix{3,3,Float64}(
            [cos(Ω)           sin(Ω)            0.;
            -sin(Ω)           cos(Ω)            0.;
             0.               0.                1.])
    R_1i = SMatrix{3,3,Float64}(
            [1.               0.                0.;
             0.               cos(i)            sin(i);
             0.              -sin(i)            cos(i)])
    R_3ω = SMatrix{3,3,Float64}(
            [cos(ω)           sin(ω)            0.
            -sin(ω)           cos(ω)            0.
             0.               0.                1.])

    ᴵTₚ = transpose(R_3ω * R_1i * R_3Ω)

    return ᴵTₚ * rₚ, ᴵTₚ * vₚ

end
inertial(orbit::Orbit) = orbit.rᵢ, orbit.vᵢ

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
perifocal(orbit::Orbit) = orbit.rₚ, orbit.vₚ

"""
Returns semimajor axis parameter, a.
"""
function semimajor_axis(r, v, μ)
   
    return inv( (2 / r) - (v^2 / μ) )

end
semimajor_axis(orbit::Orbit) = orbit.a

"""
Returns specific angular momentum vector, h̅.
"""
function specific_angular_momentum_vector(rᵢ, vᵢ)

    return rᵢ × vᵢ

end
specific_angular_momentum_vector(orbit::Orbit) = 
    specific_angular_momentum_vector(orbit.rᵢ, orbit.vᵢ)

"""
Returns scalar specific angular momentum vector, h.
"""
specific_angular_momentum(rᵢ, vᵢ) = norm(specific_angular_momentum_vector(rᵢ, vᵢ))
specific_angular_momentum(orbit::Orbit) = specific_angular_momentum(orbit.rᵢ, orbit.vᵢ)

"""
Returns specific orbital energy, ϵ.
"""
specific_energy(a, μ) = ( -μ / (2 * a) )
specific_energy(r, v, μ) = (v^2 / 2) - (μ / r)
specific_energy(orbit::Orbit) = specific_energy(orbit.a, orbit.body.μ)

"""
Returns orbital eccentricity vector e̅.
"""
function eccentricity_vector(rᵢ, vᵢ, μ)

    return map(x-> abs(x) < eps(typeof(x)) ? 0.0 : x, 
               (1 / μ) * ((vᵢ × specific_angular_momentum_vector(rᵢ, vᵢ)) - μ * rᵢ / norm(rᵢ)))

end
eccentricity_vector(orbit::Orbit) = eccentricity_vector(orbit.rᵢ, orbit.vᵢ, orbit.body.μ)

"""
Returns orbital eccentricity, e.
"""
eccentricity(rᵢ, vᵢ, μ) = norm(eccentricity_vector(rᵢ, vᵢ, μ))
eccentricity(orbit::Orbit) = orbit.e

"""
Returns semilatus parameter, p.
"""
semi_parameter(a, e) = a * (1 - e^2)
semi_parameter(orbit::Orbit) = semi_parameter(orbit.a, orbit.e)

"""
Returns radius, r.
"""
radius(p, e, ν) = upreferred(p / (1 + e * cos(ν)))
radius(orbit::Orbit) = radius(semi_parameter(orbit), orbit.e, orbit.ν)
radius(body::CelestialBody) = body.R

"""
Returns instantaneous velocity, v, for any orbital representation.
"""
velocity(r, a, μ) =  upreferred(√( (2 * μ / r) - (μ / a)))
velocity(orbit::Orbit) = velocity(radius(orbit), orbit.a, orbit.body.μ)

"""
Returns periapsis radius, rₚ.
"""
periapsis_radius(a, e) = a * (1 - e)
periapsis_radius(orbit::Orbit) = periapsis_radius(orbit.a, orbit.e)

"""
Returns periapsis radius, rₐ.
"""
apoapsis_radius(a, e) = a * (1 + e)
apoapsis_radius(orbit::Orbit) = apoapsis_radius(orbit.a, orbit.e)

"""
Returns periapsis velocity, vₚ, for any orbital representation.
"""
function periapsis_velocity(orbit::Orbit)

    return velocity(periapsis_radius(orbit), orbit.a, orbit.body.μ)

end

"""
Returns apoapsis velocity, v_a, for any orbital representation.
"""
function apoapsis_velocity(orbit::Orbit)

    return velocity(
                apoapsis_radius(orbit), 
                semimajor_axis(orbit),
                orbit.body.μ)

end

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
orbital_period(a, μ) = upreferred(2π * √(a^3 / μ))
orbital_period(orbit::Orbit) = orbital_period(orbit.a, orbit.body.μ)

"""
Returns true anomoly, ν.
"""
true_anomoly(r, h, e, μ) = acos( (h^2 - μ * r) / (μ * r * e) )
true_anomoly(orbit::Orbit) = orbit.ν

"""
Returns mean motion, n.
"""
mean_motion(a, μ) = √(μ / a^3)
mean_motion(orbit::Orbit) = mean_motion(orbit.a, orbit.μ)

"""
Returns mean motion vector, n̄.
"""
function mean_motion_vector(orbit::Orbit)

#   î = SVector{3, Float64}([1, 0, 0]) 
#   ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3, Float64}([0, 0, 1])

    return k̂ × specific_angular_momentum_vector(orbit)

end

"""
Returns eccentric anomoly, E, parabolic anomoly, B, or hyperbolic 
anomoly, H. 
"""
function conic_anomoly(orbit::Orbit{Elliptical})

    e = eccentricity(orbit)
    ν = true_anomoly(orbit)

    return acos(u"rad", (e + cos(ν) / (1 + e * cos(ν))))

end
function conic_anomoly(orbit::Orbit{Hyperbolic})

    e = eccentricity(orbit)
    ν = true_anomoly(orbit)

    return acosh(u"rad", (e + cos(ν) / (1 + e * cos(ν))))
    
end

"""
Returns time since periapsis, t.
"""
time_since_periapsis(n, e, E) = (E - e * sin(E)) / (n)
time_since_periapsis(orbit::Orbit) = 
    time_since_periapsis(
        mean_motion(orbit),
        orbit.e,
        eccentric_anomoly(orbit))

"""
Returns orbital inclination, i.
"""
inclination(orbit::Orbit) =  orbit.i

"""
Returns true if all elements in each system are within `atol` of the other.
"""
function Base.isapprox(c1::Orbit, c2::Orbit; atol=1e-8)

    return all(ustrip.(c1.rᵢ - c2.rᵢ) .< atol) &&
           all(ustrip.(c1.rᵢ - c2.rᵢ) .< atol) &&
           ustrip(upreferred(c1.e - c2.e)) .< atol &&
           ustrip(upreferred(c1.a - c2.a)) .< atol &&
           ustrip(upreferred(mod(c1.i, 180u"°") - mod(c2.i, 180u"°"))) .< atol &&
           ustrip(upreferred(mod(c1.Ω, 360u"°") - mod(c2.Ω, 360u"°"))) .< atol &&
           ustrip(upreferred(mod(c1.ω, 360u"°") - mod(c2.ω, 360u"°"))) .< atol &&
           ustrip(upreferred(mod(c1.ν, 360u"°") - mod(c2.ν, 360u"°"))) .< atol &&
           c1.body == c2.body

end

"""
Returns true if all elements of each system are identically equal.
"""
function Base.isequal(c1::Orbit, c2::Orbit)

    return all(c1.rᵢ .== c2.rᵢ) &&
           all(c1.vᵢ .== c2.vᵢ) &&
           c1.e == c2.e &&
           c1.a == c2.a &&
           mod(c1.i, 180u"°") == mod(c2.i, 180u"°") &&
           mod(c1.Ω, 360u"°") == mod(c2.Ω, 360u"°") &&
           mod(c1.ω, 360u"°") == mod(c2.ω, 360u"°") &&
           mod(c1.ν, 360u"°") == mod(c2.ν, 360u"°") &&
           c1.body == c2.body

end