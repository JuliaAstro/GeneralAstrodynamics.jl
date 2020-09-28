#
#   TwoBodyCalculations.jl
#
#   Includes simple calculations relevant to the Two Body Problem.


"""
    conic(e::T) where T<:Number
    conic(orbit::TwoBodyOrbit)

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
conic(orbit::TwoBodyOrbit) = conic(eccentricity(orbit))

"""
    TwoBodyOrbit(r̅, v̅, body)
    
Construct `TwoBodyOrbit` from Cartesian elements.
"""
function TwoBodyOrbit(r̅, v̅, body)


    return TwoBodyOrbit{conic(eccentricity(r̅, v̅, body.μ))}(
                SVector{3}(Float64.(r̅)), 
                SVector{3}(Float64.(v̅)), 
                orbital_elements(r̅, v̅, body)...,
                body)

end

"""
    TwoBodyOrbit(e, a, i, Ω, ω, ν, body)

Construct `TwoBodyOrbit` from Keplerian elements.
"""
function TwoBodyOrbit(e, a, i, Ω, ω, ν, body)

    r̅, v̅ = cartesian(e, a, i, Ω, ω, ν, body)
    return TwoBodyOrbit{conic(e)}(
                SVector{3}(Float64.(r̅)), 
                SVector{3}(Float64.(v̅)),
                map(Float64, [e,a,i,Ω,ω,ν])...,
                body)

end

"""
    orbital_elements(r̅, v̅, body::CelestialBody)

Returns a Keplarian representation of a Cartesian orbital state.
Algorithm taught in ENAE601.
"""
function orbital_elements(r̅, v̅, μ)

    î = SVector{3, Float64}([1, 0, 0]) 
    ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3, Float64}([0, 0, 1])

    acos_chop(un, val) = val ≈  1.0 ? acos(un,  1.0) : 
                         val ≈ -1.0 ? acos(un, -1.0) : 
                         acos(un, val)

    h̅ = specific_angular_momentum_vector(r̅, v̅)
    a = semimajor_axis(norm(r̅), norm(v̅), μ)
    n̅ = k̂ × specific_angular_momentum_vector(r̅, v̅)
    e̅ = eccentricity_vector(r̅, v̅, μ)
    e = norm(e̅)
    i = acos_chop(u"rad", (h̅ ⋅ k̂) / norm(h̅))

    Ω = ustrip(n̅ ⋅ ĵ) > 0 ? 
            acos_chop(u"rad", (n̅ ⋅ î) / norm(n̅)) :
            2π * u"rad" - acos_chop(u"rad", (n̅ ⋅ î) / norm(n̅))

    ω = ustrip(e̅ ⋅ k̂) > 0 ?
            acos_chop(u"rad", (n̅ ⋅ e̅) / (norm(n̅) * e)) :
            2π * u"rad" - acos_chop(u"rad", (n̅ ⋅ e̅) / (norm(n̅) * e))

    ν = ustrip(r̅ ⋅  v̅) > 0 ? 
            acos_chop(u"rad", (e̅ ⋅ r̅) / (e * norm(r̅))) :
            2π * u"rad" - acos_chop(u"rad", (e̅ ⋅ r̅) / (e * norm(r̅)))

    return e, uconvert(u"km", a), uconvert(u"°", i), 
           uconvert(u"°", Ω), uconvert(u"°", ω), 
           uconvert(u"°", ν)

end
orbital_elements(r̅, v̅, body::CelestialBody) = orbital_elements(r̅, v̅, body.μ)


"""
    cartesian(e, a, i, Ω, ω, ν, μ)

Returns a Cartesian representation of a Keplerian two-body orbital state.
Algorithm taught in ENAE601.
"""
function cartesian(e, a, i, Ω, ω, ν, μ)

    # Find semilatus parameter
    p = semi_parameter(a, e)

    # Find scalar radius
    r = radius(p, e, ν)

    # Set perifocal axes
    P̂=SVector{3, Float64}([1, 0, 0])
    Q̂=SVector{3, Float64}([0, 1, 0]) 
    Ŵ=SVector{3, Float64}([0, 0, 1])

    # Find state in Perifocal frame
    r̅ₚ = (r * cos(ν) .* P̂ .+ r * sin(ν) .* Q̂)
    v̅ₚ = √(μ/p) * ((-sin(ν) * P̂) .+ ((e + cos(ν)) .* Q̂))

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

    ᴵTₚ = (R_3ω * R_1i * R_3Ω)' 

    return uconvert.(u"km", ᴵTₚ * r̅ₚ), uconvert.(u"km/s", ᴵTₚ * v̅ₚ)

end
cartesian(e, a, i, Ω, ω, ν, body::CelestialBody) = cartesian(e, a, i, Ω, ω, ν, body.μ)


"""
    semimajor_axis(r, v, μ)
    semimajor_axis(orbit::TwoBodyOrbit)

Returns semimajor axis parameter, a.
"""
function semimajor_axis(r, v, μ)
   
    return inv( (2 / r) - (v^2 / μ) )

end
semimajor_axis(orbit::TwoBodyOrbit) = orbit.a

"""
    specific_angular_momentum_vector(r̅, v̅)
    specific_angular_momentum_vector(orbit::TwoBodyOrbit)

Returns specific angular momentum vector, h̅.
"""
function specific_angular_momentum_vector(r̅, v̅)

    return r̅ × v̅

end
specific_angular_momentum_vector(orbit::TwoBodyOrbit) = 
    specific_angular_momentum_vector(orbit.r̅, orbit.v̅)

"""
    specific_angular_momentum(r̅, v̅)
    specific_angular_momentum(orbit::TwoBodyOrbit)

Returns scalar specific angular momentum vector, h.
"""
specific_angular_momentum(r̅, v̅) = norm(specific_angular_momentum_vector(r̅, v̅))
specific_angular_momentum(orbit::TwoBodyOrbit) = specific_angular_momentum(orbit.r̅, orbit.v̅)

"""
    specific_energy(a, μ)
    specific_energy(r, v, μ)
    specific_energy(orbit::TwoBodyOrbit)

Returns specific orbital energy, ϵ.
"""
specific_energy(a, μ) = ( -μ / (2 * a) )
specific_energy(r, v, μ) = (v^2 / 2) - (μ / r)
specific_energy(orbit::TwoBodyOrbit) = specific_energy(orbit.a, orbit.body.μ)

"""
    eccentricity_vector(r̅, v̅, μ)
    eccentricity_vector(orbit::TwoBodyOrbit)

Returns orbital eccentricity vector e̅.
"""
function eccentricity_vector(r̅, v̅, μ)

    return (1 / μ) * ((v̅ × specific_angular_momentum_vector(r̅, v̅)) - μ * r̅ / norm(r̅))

end
eccentricity_vector(orbit::TwoBodyOrbit) = eccentricity_vector(orbit.r̅, orbit.v̅, orbit.body.μ)

"""
    eccentricity(r̅, v̅, μ)
    eccentricity(orbit::TwoBodyOrbit)

Returns orbital eccentricity, e.
"""
eccentricity(r̅, v̅, μ) = norm(eccentricity_vector(r̅, v̅, μ))
eccentricity(orbit::TwoBodyOrbit) = orbit.e

"""
    semi_parameter(a, e)
    semi_parameter(orbit::TwoBodyOrbit)

Returns semilatus parameter, p.
"""
semi_parameter(a, e) = a * (1 - e^2)
semi_parameter(orbit::TwoBodyOrbit) = semi_parameter(orbit.a, orbit.e)

"""
    radius(p, e, ν)
    radius(orbit::TwoBodyOrbit)

Returns instantaneous radius, r.
"""
radius(p, e, ν) = p / (1 + e * cos(ν))
radius(orbit::TwoBodyOrbit) = radius(semi_parameter(orbit), orbit.e, orbit.ν)

"""
    velocity(r, a, μ)
    velocity(orbit::TwoBodyOrbit)

Returns instantaneous velocity, v, for any orbital representation.
"""
velocity(r, a, μ) =  √( (2 * μ / r) - (μ / a))
velocity(orbit::TwoBodyOrbit) = velocity(radius(orbit), orbit.a, orbit.body.μ)

"""
    periapsis_radius(a, e)
    periapsis_radius(orbit::TwoBodyOrbit)

Returns periapsis radius, r̅_p.
"""
periapsis_radius(a, e) = a * (1 - e)
periapsis_radius(orbit::TwoBodyOrbit) = periapsis_radius(orbit.a, orbit.e)

"""
    apoapsis_radius(a, e)
    apoapsis_radius(orbit::TwoBodyOrbit)

Returns periapsis radius, r_a.
"""
apoapsis_radius(a, e) = a * (1 + e)
apoapsis_radius(orbit::TwoBodyOrbit) = apoapsis_radius(orbit.a, orbit.e)

"""
    periapsis_velocity(orbit::T) where T<:TwoBodyOrbit

Returns periapsis velocity, v_p, for any orbital representation.
"""
function periapsis_velocity(orbit::TwoBodyOrbit)

    return velocity(periapsis_radius(orbit), orbit.a, orbit.body.μ)

end

"""
    apoapsis_velocity(orbit::T) where T<:TwoBodyOrbit

Returns apoapsis velocity, v_a, for any orbital representation.
"""
function apoapsis_velocity(orbit::TwoBodyOrbit)

    return velocity(
                apoapsis_radius(orbit), 
                semimajor_axis(orbit),
                orbit.body.μ)

end

"""
    orbital_period(a, μ)
    orbital_period(orbit::TwoBodyOrbit)

Returns the orbital period.
"""
orbital_period(a, μ) = 2π * √(a^3 / μ)
orbital_period(orbit::TwoBodyOrbit) = orbital_period(orbit.a, orbit.body.μ)

"""
    true_anomoly(r, h, e, μ)

Returns true anomoly, ν.
"""
true_anomoly(r, h, e, μ) = acos( (h^2 - μ * r) / (μ * r * e) )
true_anomoly(orbit::TwoBodyOrbit) = orbit.ν

"""
    mean_motion(a, μ)
    mean_motion(orbit::TwoBodyOrbit)

Returns mean motion, n.
"""
mean_motion(a, μ) = √(μ / a^3)
mean_motion(orbit::TwoBodyOrbit) = mean_motion(orbit.a, orbit.μ)

"""
    mean_motion̅tor(orbit::TwoBodyOrbit)

Returns mean motion vector, n̄.
"""
function mean_motion̅tor(orbit::TwoBodyOrbit)

#   î = SVector{3, Float64}([1, 0, 0]) 
#   ĵ = SVector{3, Float64}([0, 1, 0]) 
    k̂ = SVector{3, Float64}([0, 0, 1])

    return k̂ × specific_angular_momentum_vector(orbit)

end

"""
    conic_anomoly(orbit::TwoBodyOrbit{Elliptical})
    conic_anomoly(orbit::TwoBodyOrbit{Hyperbolic})

Returns eccentric anomoly, E, parabolic anomoly, B, or hyperbolic 
anomoly, H. 
"""
function conic_anomoly(orbit::TwoBodyOrbit{Elliptical})

    e = eccentricity(orbit)
    ν = true_anomoly(orbit)

    return acos(u"rad", (e + cos(ν) / (1 + e * cos(ν))))

end
function conic_anomoly(orbit::TwoBodyOrbit{Hyperbolic})

    e = eccentricity(orbit)
    ν = true_anomoly(orbit)

    return acosh(u"rad", (e + cos(ν) / (1 + e * cos(ν))))
    
end

"""
    time_since_periapsis(n, e, E)
    time_since_periapsis(orbit::TwoBodyOrbit)

Returns time since periapsis, t.
"""
time_since_periapsis(n, e, E) = (E - e * sin(E)) / (n)
time_since_periapsis(orbit::TwoBodyOrbit) = 
    time_since_periapsis(
        mean_motion(orbit),
        orbit.e,
        eccentric_anomoly(orbit))

"""
    inclination(orbit::TwoBodyOrbit)

Returns orbital inclination, i.
"""
inclination(orbit::TwoBodyOrbit) =  orbit.i

"""
    isapprox(::TwoBodyOrbit, ::TwoBodyOrbit; atol=1e-8)

Returns true if all elements in each system are within `atol` of the other.
"""
function Base.isapprox(c1::TwoBodyOrbit, c2::TwoBodyOrbit; atol=1e-8)

    return all(ustrip.(c1.r̅ - c2.r̅) .< atol) &&
           all(ustrip.(c1.r̅ - c2.r̅) .< atol) &&
           ustrip(upreferred(c1.e - c2.e)) .< atol &&
           ustrip(upreferred(c1.a - c2.a)) .< atol &&
           ustrip(upreferred(mod(c1.i, 180u"°") - mod(c2.i, 180u"°"))) .< atol &&
           ustrip(upreferred(mod(c1.Ω, 360u"°") - mod(c2.Ω, 360u"°"))) .< atol &&
           ustrip(upreferred(mod(c1.ω, 360u"°") - mod(c2.ω, 360u"°"))) .< atol &&
           ustrip(upreferred(mod(c1.ν, 360u"°") - mod(c2.ν, 360u"°"))) .< atol &&
           c1.body == c2.body

end

"""
    isequal(::TwoBodyOrbit, ::TwoBodyOrbit)

Returns true if all elements of each system are identically equal.
"""
function Base.isequal(c1::TwoBodyOrbit, c2::TwoBodyOrbit)

    return all(c1.r̅ .== c2.r̅) &&
           all(c1.v̅ .== c2.v̅) &&
           c1.e == c2.e &&
           c1.a == c2.a &&
           mod(c1.i, 180u"°") == mod(c2.i, 180u"°") &&
           mod(c1.Ω, 360u"°") == mod(c2.Ω, 360u"°") &&
           mod(c1.ω, 360u"°") == mod(c2.ω, 360u"°") &&
           mod(c1.ν, 360u"°") == mod(c2.ν, 360u"°") &&
           c1.body == c2.body

end