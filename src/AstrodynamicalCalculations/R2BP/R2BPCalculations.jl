#
# Provides Restricted Two-body Problem equations
#

"""
Returns the conic section, as specified by eccentricity `e`.
"""
function conic(e::T) where T<:Real

    !isnan(e) || throw(ArgumentError("Provided eccentricity is not a number (NaN)."))

    if e ≈ zero(T)
        return Circular
    elseif e ≈ one(T)
        return Parabolic
    elseif zero(T) < e && e < one(T)
        return Elliptical
    elseif e > one(T)
        return Hyperbolic
    else
        throw(ArgumentError("Provided eccentricity, $e, is not valid."))
    end

end

"""
Returns a Keplarian representation of a Cartesian orbital state.
Algorithm taught in ENAE601.
"""
function keplerian(rᵢ, vᵢ, μ)

    safe_acos(num) = isapprox(num, one(num)) ? 
                            acos(one(num)) : 
                                isapprox(num, -one(num)) ? 
                                    acos(-one(num)) : 
                                        acos(num)

    î = SVector{3}(1, 0, 0) 
    ĵ = SVector{3}(0, 1, 0) 
    k̂ = SVector{3}(0, 0, 1)

    h̅ = specific_angular_momentum_vector(rᵢ, vᵢ)
    a = semimajor_axis(norm(rᵢ), norm(vᵢ), μ)
    n̅ = k̂ × specific_angular_momentum_vector(rᵢ, vᵢ)
    e̅ = eccentricity_vector(rᵢ, vᵢ, μ)
    e = norm(e̅)

    i = safe_acos((h̅ ⋅ k̂) / norm(h̅))

    Ω = ustrip(n̅ ⋅ ĵ) > 0 ? 
            safe_acos((n̅ ⋅ î) / norm(n̅)) :
            2π - safe_acos((n̅ ⋅ î) / norm(n̅))

    ω = ustrip(e̅ ⋅ k̂) > 0 ?
            safe_acos((n̅ ⋅ e̅) / (norm(n̅) * e)) :
            2π - safe_acos((n̅ ⋅ e̅) / (norm(n̅) * e))

    ν = ustrip(rᵢ ⋅  vᵢ) > 0 ? 
            safe_acos((e̅ ⋅ rᵢ) / (e * norm(rᵢ))) :
            2π - safe_acos((e̅ ⋅ rᵢ) / (e * norm(rᵢ)))

    return upreferred(e), upreferred(a), upreferred(i), 
           upreferred(Ω), upreferred(ω), 
           upreferred(ν)

end

"""
Returns a Cartesian representation of a Keplerian two-body orbital state
in an inertial frame, centered at the center of mass of the central body.
Algorithm taught in ENAE601.
"""
function cartesian(e, a, i, Ω, ω, ν, μ)
    rᵢ, vᵢ = spatial(i, Ω, ω, perifocal(a, e, ν, μ)...)
    return  upreferred.(rᵢ), 
            upreffered.(vᵢ)
end

"""
Returns a spatial representation of the provied Perifocal state.
"""
function spatial(i, Ω, ω, rₚ, vₚ)

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
        r = distance(p, e, ν)
        
        P̂=SVector{3, Float64}(1, 0, 0)
        Q̂=SVector{3, Float64}(0, 1, 0)
        # Ŵ=SVector{3, Float64}(0, 0, 1)
        
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



"""
Returns semimajor axis parameter, a.
"""
semimajor_axis(r, v, μ) = inv( (2 / r) - (v^2 / μ) )

"""
Returns specific angular momentum vector, h̅.
"""
specific_angular_momentum_vector(rᵢ, vᵢ) = rᵢ × vᵢ

"""
Returns scalar specific angular momentum vector, h.
"""
specific_angular_momentum(rᵢ, vᵢ) = norm(specific_angular_momentum_vector(rᵢ, vᵢ))

"""
Returns specific orbital energy, ϵ.
"""
specific_energy(a, μ) = ( -μ / (2 * a) )
specific_energy(r, v, μ) = (v^2 / 2) - (μ / r)

"""
Returns C3 value.
"""
C3(r, v, μ) = v^2 - 2μ/r

"""
Returns v∞.
"""
v_infinity(r, v, μ) = √C3(r, v, μ)

"""
Returns potential energy for an orbit about a `RestrictedTwoBodySystem`.
"""
specific_potential_energy(r, μ) = (μ/r)
specific_potential_energy(r, μ, R, J₂, ϕ) = (μ/r) * (1 - J₂ * (R/r)^2 * ((3/2) * (sin(ϕ))^2 - (1/2)))

"""
Returns orbital eccentricity vector e̅.
"""
function eccentricity_vector(rᵢ, vᵢ, μ)

    return map(x-> abs(x) < eps(typeof(x)) ? zero(x) : x, (1 / μ) * ((vᵢ × specific_angular_momentum_vector(rᵢ, vᵢ)) - μ * rᵢ / norm(rᵢ)))

end

"""
Returns orbital eccentricity, e.
"""
eccentricity(rᵢ, vᵢ, μ) = norm(eccentricity_vector(rᵢ, vᵢ, μ)) |> upreferred

"""
Returns semilatus parameter, p.
"""
semi_parameter(a, e) = a * (1 - e^2)

"""
Returns distance, r.
"""
distance(p, e, ν) = upreferred(p / (1 + e * cos(ν)))

"""
Returns instantaneous velocity, v, for any orbital representation.
"""
speed(r, a, μ) =  upreferred(√( (2 * μ / r) - (μ / a)))

"""
Returns the instantaneous velocity, `v`, for any orbital representation.
"""
speed(p, e, ν, a, μ) = speed(distance(p,e,ν), a, μ)

"""
Returns periapsis distance, rₚ.
"""
periapsis_radius(a, e) = a * (1 - e)

"""
Returns apoapsis distance, rₐ.
"""
apoapsis_radius(a, e) = a * (1 + e)

"""
Returns the orbital period.
"""
period(a, μ) = 2π * √(upreferred(a^3 / μ))

"""
Returns true anomoly, ν.
"""
function true_anomoly(r, h, e, μ)
    val = (h^2 - μ * r) / (μ * r * e)
    acos(isapprox(val, one(val)) ? one(val) : val)
end

"""
Returns mean motion, n.
"""
mean_motion(a, μ) = √(μ / a^3)

"""
Returns time since periapsis, t.
"""
time_since_periapsis(n, e, E) = (E - e * sin(E)) / (n)

"""
Sphere of influence.
"""
SOI(a, m, M) = a * (m / M)^(2/5)

"""
Sphere of activity.
"""
SOA(a, m, M) = a * (m / 3M)^(1/3)

"""
Computes a Hohmann transfer, and returns the departure and 
arrival velocity vectors. 
"""
function hohmann(r₁, r₂, μ)
		
    vₐ = √((2μ/r₁) - (2μ/(r₁+r₂)))
    vₚ = √((2μ/r₂) - (2μ/(r₁+r₂)))
    
    return vₐ, vₚ
    
end
