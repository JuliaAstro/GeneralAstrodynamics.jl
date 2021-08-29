#
#   Kepler.jl
#
#   Solves Kepler's problem for TwoBody orbits.
#

"""
Solves Kepler's Problem for `orbit` and `Δtᵢ`.
"""
function kepler(r::AbstractVector, v::AbstractVector, μ::Number, Δt::Number = period(semimajor_axis(r,v,μ), μ); tol=1e-6, max_iter=100) 

    e, a, i, Ω, ω, ν = keplerian(r, v, μ)
    T = period(a, μ)
    conic_section = conic(e)

    # Guess χ₀
    if conic_section == Hyperbolic
        χ₀ = sign(Δt) * √(-a) * log(ℯ, (-2 * μ / a * Δt) / (r ⋅ v + (sign(Δt) * √(-μ * a) * (1 - norm(r) / a))))
    elseif conic_section == Parabolic
        χ₀ = √(a) * tan(ν / 2)
    else
        Δt = mod(Δt, T)
        χ₀ = √(μ) * Δt / a
    end

    # Iteratively solve for χ
    # TODO: Compare loop vs. recursion performance here.
    # There shouldn't be too large of a difference, since this tends
    # to converge with only a few iterations.
    χₙ, rₙ, ψ, C₂, C₃ = χₖ(χ₀, Δt, r, v, a, μ, tol=tol, max_iter=max_iter)

    # Convert to a Orbit
    f = 1 - χₙ^2 / norm(r) * C₂
    ḟ = √(μ) / (norm(r) * rₙ) * χₙ * (ψ * C₃ - 1)
    g = Δt - (χₙ^3 / √(μ)) * C₃
    ġ = 1 - (χₙ^2 / rₙ) * C₂

    return ((f * r) .+ (g * v), (ḟ * r) .+ (ġ * v))
end


function χₖ(χₙ, Δt, rᵢ₀, vᵢ₀, a, μ; iter=1, tol=1e-14, max_iter=100)
    
    r₀ = norm(rᵢ₀)
    ψ = upreferred(χₙ^2 / a)

    if ψ > tol
        C₂ = (1 - cos(√(ψ))) /  ψ
        C₃ = (√(ψ) - sin(√(ψ))) / √(ψ^3)
    elseif ψ < -tol
        C₂ = (1 - cosh(√(-ψ))) / ψ
        C₃ = (sinh(√(-ψ)) - √(-ψ)) / √((-ψ)^3)
    else
        C₂ = 1.0 / 2.0
        C₃ = 1.0 / 6.0
    end

    r = χₙ^2 * C₂ + (rᵢ₀ ⋅ vᵢ₀) * χₙ / √(μ) * (1 - ψ*C₃) + r₀ * (1 - ψ * C₂)
    χₙ₊₁ = χₙ + ((√(μ) * Δt - χₙ^3 * C₃ - (rᵢ₀ ⋅ vᵢ₀) / √(μ) * χₙ^2 * C₂ - r₀ * χₙ * (1 - ψ * C₃)) / r)

    if iter > max_iter
        @error "Failed to converge!"
        return χₙ, r, ψ, C₂, C₃
    elseif abs(χₙ₊₁ - χₙ) < oneunit(χₙ) * tol
        return  χₙ, r, ψ, C₂, C₃
    else
        return χₖ(χₙ₊₁, Δt, rᵢ₀, vᵢ₀, a, μ; iter=iter+1, tol=tol, max_iter=max_iter)
    end
    
end