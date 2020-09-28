#
#   Kepler.jl
#
#   Solves Kepler's problem for TwoBody orbits.
#

"""
    kepler(orbit::T, Δtᵢ::N = orbital_period(orbit)) where {T<:TwoBodyOrbit, N<:Number}

Solves Kepler's Problem for `orbit` and `Δtᵢ`.
"""
function kepler(orbit::T, Δtᵢ::N = orbital_period(orbit)) where {T<:TwoBodyOrbit, N<:Number}

    conic_section = conic(orbit)

    # Guess χ₀
    if conic_section ∈ (Circular, Elliptical)

        Δt = mod(Δtᵢ, orbital_period(orbit))
        χ₀ = √(orbit.body.μ) * Δt / orbit.a

    elseif conic_section == Hyperbolic
        
        Δt = Δtᵢ
        χ₀ = sign(Δt) * √(-orbit.a) * log(ℯ, (-2 * orbit.body.μ / orbit.a * Δt) / 
                (orbit.r̅ ⋅ orbit.v̅₀ + (sign(Δt) * √(-orbit.body.μ * orbit.a) * (1 - norm(orbit.r̅) / orbit.a))))

    elseif conic_section == Parabolic

        Δt = Δtᵢ
        χ₀ = √(semi_parameter(orbit)) * tan(orbit.ν / 2)

    else

        @warn "Kepler's problem failed to converge."
        return InvalidTwoBodyOrbit(orbit.body)

    end

    # Iteratively solve for χ
    # TODO: Compare loop vs. recursion performance here.
    # There shouldn't be too large of a difference, since this tends
    # to converge with only a few iterations.
    χₙ, r, ψ, C₂, C₃ = χ(χ₀, Δt, orbit.r̅, orbit.v̅, orbit.a, orbit.body.μ)

    # Convert to a TwoBodyOrbit
    f = 1 - χₙ^2 / norm(orbit.r̅) * C₂
    ḟ = √(orbit.body.μ) / (norm(orbit.r̅) * r) * χₙ * (ψ * C₃ - 1)
    g = Δt - (χₙ^3 / √(orbit.body.μ)) * C₃
    ġ = 1 - (χₙ^2 / r) * C₂

    return TwoBodyOrbit(f * orbit.r̅ + g * orbit.v̅, ḟ * orbit.r̅ + ġ * orbit.v̅, orbit.body)

end

function χ(χₙ, Δt, r̅₀, v̅₀, a, μ; iter=1, tol=1e-14, max_iter=100)
    
    r₀ = norm(r̅₀)
    ψ = χₙ^2 / a

    if ψ > 1e-6
        C₂ = (1 - cos(√(ψ))) /  ψ
        C₃ = (√(ψ) - sin(√(ψ))) / √(ψ^3)
    elseif ψ < -1e-6
        C₂ = (1 - cosh(√(-ψ))) / ψ
        C₃ = (sinh(√(-ψ)) - √(-ψ)) / √((-ψ)^3)
    else
        C₂ = 1.0 / 2.0
        C₃ = 1.0 / 6.0
    end

    r = χₙ^2 * C₂ + (r̅₀ ⋅ v̅₀) * χₙ / √(μ) * (1 - ψ*C₃) + r₀ * (1 - ψ * C₂)
    χₙ₊₁ = χₙ + ((√(μ) * Δt - χₙ^3 * C₃ - (r̅₀ ⋅ v̅₀) / √(μ) * χₙ^2 * C₂ - r₀ * χₙ * (1 - ψ * C₃)) / r)

    if iter > max_iter
        return NaN, NaN, NaN, NaN, NaN
    elseif abs(χₙ₊₁ - χₙ) < oneunit(χₙ) * tol
        return  χₙ, r, ψ, C₂, C₃
    else
        return χ(χₙ₊₁, Δt, r̅₀, v̅₀, a, μ; iter=iter+1, tol=tol, max_iter=max_iter)
    end
    
end