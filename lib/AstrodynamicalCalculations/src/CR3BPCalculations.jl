"""
Common calculations within Circular Restricted Three Body Problem dynamics.

# Extended Help

## Imports

$(IMPORTS)

## Exports

$(EXPORTS)

"""
module CR3BPCalculations

using StaticArrays
using LinearAlgebra
using DocStringExtensions
using Roots

@template (
    FUNCTIONS,
    METHODS,
    MACROS,
) = """
    $(SIGNATURES)

    !!! warning "CR3BP Dynamics"
        This computation is valid for Circular Restricted Three Body Problem dynamics.

    $(DOCSTRING)
    """
export distance_scaling,
    time_scaling,
    reduced_mass,
    velocity_scaling,
    nondimensional,
    redimensioned,
    nondimensional_radii,
    distance_to_primary,
    distance_to_secondary,
    potential_energy,
    jacobi_constant,
    primary_synodic_position,
    secondary_synodic_position,
    lagrange_point,
    zero_velocity_curves,
    richardson_halo,
    richardson_ic,
    perturbation,
    convergent_direction,
    divergent_direction,
    diverge,
    diverge!,
    converge,
    converge!,
    perturb,
    perturb!

"""
The length scale factor used to nondimensionalize CR3BP states.
"""
Base.@pure distance_scaling(a) = a
Base.@pure distance_scaling(a, μ₁, μ₂) = a

"""
The time scale factor used to nondimensionalize CR3BP states.
"""
time_scaling(a, μ₁, μ₂) = 2π * √(a^3 / (μ₁ + μ₂))

"""
The velocity scale factor used to nondimensionalize CR3BP states.
"""
velocity_scaling(a, μ₁, μ₂) = distance_scaling(a) / time_scaling(a, μ₁, μ₂)

"""
Return the reduced mass.
"""
reduced_mass(μ₁, μ₂) = min(μ₁, μ₂) / (μ₁ + μ₂)

"""
Normalizes a CR3BP orbit in the rotating reference frame.
"""
function nondimensional(x, y, z, ẋ, ẏ, ż, t, a, μ₁, μ₂)
    r = SVector(x, y, z)
    v = SVector(ẋ, ẏ, ż)

    x, y, z = r / a
    μₜ = μ₁ + μ₂
    Tₛ = 2π * √(a^3 / μₜ)
    ẋ, ẏ, ż = v / (a / Tₛ)
    t = t / Tₛ
    μ = min(μ₁, μ₂) / μₜ

    return (; x, y, z, ẋ, ẏ, ż, t, μ, L = a, T = Tₛ)
end

"""
Redimensionalizes a CR3BP orbit in the rotating reference frame.
"""
function redimensioned(x, y, z, ẋ, ẏ, ż, t, μ, a, T)
    rₙ = SVector(x, y, z)
    vₙ = SVector(ẋ, ẏ, ż)

    x, y, z = rₙ * a
    ẋ, ẏ, ż = vₙ * a / T
    t = t * T

    sum_μs = a^3 / ((T / 2π)^2)
    μ₂ = μ * sum_μs
    μ₁ = sum_μ - μ₂

    return (; x, y, z, ẋ, ẏ, ż, t, μ₁, μ₂, L = a, T)

end

"""
Returns the spacecraft's nondimensional position w.r.t. body 1 (or 2).
"""
nondimensional_radii(x, y, z, xᵢ) = sqrt((x - xᵢ)^2 + y^2 + z^2)

"""
Returns synodic distance to primary body.
"""
distance_to_primary(r, μ) = r + μ

"""
Returns synodic distance to secondary body.
"""
distance_to_secondary(r, μ) = r - (one(μ) - μ)

"""
Returns the potential energy `U` in the Synodic frame with Normalized units.
"""
function potential_energy(x, y, z, μ)
    return (x^2 + y^2) / 2 +
           ((1 - μ) / nondimensional_radii(x, y, z, -μ)) +
           (μ / nondimensional_radii(x, y, z, 1 - μ))
end

"""
Returns the Jacobi Constant, `C` in the Synodic frame with Normalized units.
"""
function jacobi_constant(x, y, z, ẋ, ẏ, ż, μ)
    v = SVector(ẋ, ẏ, ż)
    return 2 * potential_energy(x, y, z, μ) - (v ⋅ v)
end

"""
Given the Synodic frame vector, returns the vector in the barycentric-inertial reference frame.
"""
function synodic_to_inertial(x, y, z, t, ω)
    rₛ = SVector(x, y, z)

    θ = ω * t
    O = zero(θ)
    l = one(θ)

    ᴵTₛ = SMatrix{3,3}(cos(θ), -sin(θ), O, sin(θ), cos(θ), O, O, O, l)

    x, y, z = ᴵTₛ * rₛ
    return (; x, y, z)
end

"""
Given a barycentric-inertial Cartesian state, returns the state in the synodic (rotating) reference frame.
"""
function inertial_to_synodic(x, y, z, t, ω)

    rᵢ = SVector(x, y, z)
    θ = ω * t
    O = zero(θ)
    l = one(θ)

    ᴵTₛ = SMatrix{3,3}(cos(θ), -sin(θ), O, sin(θ), cos(θ), O, O, O, l)

    ˢTᵢ = inv(ᴵTₛ)

    x, y, z = ˢTᵢ * rᵢ

    return (; x, y, z)
end

"""
Position of primary body.
"""
function primary_synodic_position(μ)
    return SVector{3}(-μ, zero(μ), zero(μ))
end

"""
Position of primary body.
"""
function secondary_synodic_position(μ)
    return SVector{3}(one(μ) - μ, zero(μ), zero(μ))
end

"""
Returns the lagrange points for a CR3BP system.

__Arguments:__ 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `L`: Langrange points requested, must be in range [1,5]

__Outputs:__
- Tuple of Lagrange points
- Throws `ArgumentError` if L is out of range [1,5]

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/)
"""
function lagrange_point(μ::Real, L::Int)
    if L < 1 || L > 5
        error(
            "Lagrange point index must be a value between 1 and 5 (inclusive). You provided: $L.",
        )
    end

    expressions = SVector{3}(
        x -> x - (1 - μ) / (x + μ)^2 + μ / (x + μ - 1)^2,
        x -> x - (1 - μ) / (x + μ)^2 - μ / (x + μ - 1)^2,
        x -> x + (1 - μ) / (x + μ)^2 + μ / (x + μ + 1)^2,
    )

    return (
        map(f -> [find_zero(f, (-3, 3)), 0, 0], expressions)...,
        [(1 / 2) - μ, √(3) / 2, 0],
        [(1 / 2) - μ, -√(3) / 2, 0],
    )[L]
end

"""
Returns an analytical solution for a Halo orbit about `L`.

# Extended Help

## Arguments 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `L`: Lagrange point to orbit (L1 or L2).
- `Z`: Desired non-dimensional Z-amplitude for Halo orbit.
- `hemisphere`: Specifies northern or southern Halo orbit.
- `ϕ`: Desired Halo orbit phase.
- `steps`: Number of non-dimensional timepoints in returned state.

## Outputs
- Near-periodic initial condition `u`
- Halo orbit period `T`.
- Throws `ArgumentError` if L is not `1` or `2`.

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function richardson_ic(μ, L::Int; Z = 0.0, hemisphere = :northern, ϕ = 0.0)
    if L == 1
        point = first(lagrange_point(μ, 1))
        γ = abs(one(μ) - μ - point)
        n = collect(1:4) .* one(μ)
        c = @. (μ + (-one(1))^n * (one(μ) - μ)γ^(n + 1)) / (γ^3 * (one(μ) - γ^(n + 1)))
    elseif L == 2
        point = first(lagrange_point(μ, 2))
        γ = abs(point - one(μ) + μ)
        n = collect(1:4) .* one(μ)
        c = @. ((-one(μ))^n * μ + (-one(μ))^n * (one(μ) - μ)γ^(n + 1)) /
           (γ^3 * (one(μ) + γ^(n + 1)))
    else
        throw(ArgumentError("Only Halo orbits about L1 or L2 are supported."))
    end

    ωₚ = √((2 - c[2] + √((9c[2]^2 - 8c[2]))) / 2)
    k = (ωₚ^2 + 1 + 2c[2]) / (2ωₚ)

    d₁ = (3ωₚ^2 / k) * (k * (6ωₚ^2 - 1) - 2ωₚ)
    d₂ = (8ωₚ^2 / k) * (k * (11ωₚ^2 - 1) - 2ωₚ)
    a₂₁ = (3c[3] * (k^2 - 2)) / (4(1 + 2c[2]))
    a₂₂ = (3c[3]) / (4(1 + 2c[2]))
    a₂₃ = (-3c[3]ωₚ / (4k * d₁)) * (3k^3 * ωₚ - 6k * (k - ωₚ) + 4)
    a₂₄ = (-3c[3]ωₚ / (4k * d₁)) * (2 + 3k * ωₚ)
    b₂₁ = (-3c[3]ωₚ / (2d₁)) * (3k * ωₚ - 4)
    b₂₂ = -3c[3] * ωₚ / d₁
    d₂₁ = -c[3] / (2ωₚ^2)
    a₃₁ =
        (-9ωₚ / (4d₂)) * (4c[3] * (k * a₂₃ - b₂₁) + k * c[4] * (4 + k^2)) +
        ((9ωₚ^2 + 1 - c[2]) / (2d₂)) * (3c[3] * (2a₂₃ - k * b₂₁) + c[4] * (2 + 3k^2))
    a₃₂ =
        (-9ωₚ / (4d₂)) * (4c[3] * (3k * a₂₄ - b₂₂) + k * c[4]) -
        (3 / (2d₂)) * (9ωₚ^2 + 1 - c[2]) * (c[3] * (k * b₂₂ + d₂₁ - 2a₂₄) - c[4])
    b₃₁ =
        (3 / (8d₂)) * 8ωₚ * (3c[3] * (k * b₂₁ - 2a₂₃) - c[4] * (2 + 3k^2)) +
        (3 / (8d₂)) * (9ωₚ^2 + 1 + 2c[2]) * (4c[3] * (k * a₂₃ - b₂₁) + k * c[4] * (4 + k^2))
    b₃₂ =
        (9ωₚ / d₂) * (c[3] * (k * b₂₂ + d₂₁ - 2a₂₄) - c[4]) +
        (3(9ωₚ^2 + 1 + 2c[2]) / (8d₂) * (4c[3] * (k * a₂₄ - b₂₂) + k * c[4]))
    d₃₁ = (3 / (64ωₚ^2)) * (4c[3] * a₂₄ + c[4])
    d₃₂ = (3 / (64 + ωₚ^2)) * (4c[3] * (a₂₃ - d₂₁) + c[4] * (4 + k^2))

    s₁ =
        (1 / (2ωₚ * (ωₚ * (1 + k^2) - 2k))) * (
            3c[3] / 2 * (2a₂₁ * (k^2 - 2) - a₂₃ * (k^2 + 2) - 2k * b₂₁) -
            (3c[4] / 8) * (3k^4 - 8k^2 + 8)
        )
    s₂ =
        (1 / (2ωₚ * (ωₚ * (1 + k^2) - 2k))) * (
            3c[3] / 2 * (2a₂₂ * (k^2 - 2) + a₂₄ * (k^2 + 2) + 2k * b₂₂ + 5d₂₁) +
            (3c[4] / 8) * (12 - k^2)
        )
    l₁ = (-3c[3] / 2) * (2a₂₁ + a₂₃ + 5d₂₁) - (3c[4] / 8) * (12 - k^2) + 2ωₚ^2 * s₁
    l₂ = (3c[3] / 2) * (a₂₄ - 2a₂₂) + (9c[4] / 8) + 2ωₚ^2 * s₂
    Δ = ωₚ^2 - c[2]

    Aᵧ = Z / γ
    Aₓ = √((-l₂ * Aᵧ^2 - Δ) / l₁)

    ν = 1 + s₁ * Aₓ^2 + s₂ * Aᵧ^2
    T = 2π / (ωₚ * ν)
    τ = 0

    if hemisphere == :northern
        m = 1
    elseif hemisphere == :southern
        m = 3
    else
        throw(ArgumentError("`hemisphere` must be `:northern` or `:southern`."))
    end

    δₘ = 2 - m
    τ₁ = @. ωₚ * τ + ϕ

    X = @. γ * (
        a₂₁ * Aₓ^2 + a₂₂ * Aᵧ^2 - Aₓ * cos(τ₁) +
        (a₂₃ * Aₓ^2 - a₂₄ * Aᵧ^2) * cos(2τ₁) +
        (a₃₁ * Aₓ^3 - a₃₂ * Aₓ * Aᵧ^2) * cos(3τ₁)
    ) + 1 - μ - (L == 1 ? γ : -γ)
    Y = @. γ * (
        k * Aₓ * sin(τ₁) +
        (b₂₁ * Aₓ^2 - b₂₂ * Aᵧ^2) * sin(2τ₁) +
        (b₃₁ * Aₓ^3 - b₃₂ * Aₓ * Aᵧ^2) * sin(3τ₁)
    )
    Z = @. γ * (
        δₘ * Aᵧ * cos(τ₁) +
        δₘ * d₂₁ * Aₓ * Aᵧ * (cos(2τ₁) - 3) +
        δₘ * (d₃₂ * Aᵧ * Aₓ^2 - d₃₁ * Aᵧ^3) * cos(3τ₁)
    )

    Ẋ = @. γ * (
        ωₚ * ν * Aₓ * sin(τ₁) - 2ωₚ * ν * (a₂₃ * Aₓ^2 - a₂₄ * Aᵧ^2) * sin(2τ₁) -
        3ωₚ * ν * (a₃₁ * Aₓ^3 - a₃₂ * Aₓ * Aᵧ^2) * sin(3τ₁)
    )
    Ẏ = @. γ * (
        ωₚ * ν * k * Aₓ * cos(τ₁) +
        2ωₚ * ν * (b₂₁ * Aₓ^2 - b₂₂ * Aᵧ^2) * cos(2τ₁) +
        3ωₚ * ν * (b₃₁ * Aₓ^3 - b₃₂ * Aₓ * Aᵧ^2) * cos(3τ₁)
    )
    Ż = @. γ * (
        -ωₚ * ν * δₘ * Aᵧ * sin(τ₁) - 2ωₚ * ν * δₘ * d₂₁ * Aₓ * Aᵧ * sin(2τ₁) -
        3ωₚ * ν * δₘ * (d₃₂ * Aᵧ * Aₓ^2 - d₃₁ * Aᵧ^2) * sin(3τ₁)
    )

    return (; x = X, y = Y, z = Z, ẋ = Ẋ, ẏ = Ẏ, ż = Ż, Δt = T)

end



"""
Returns an analytical solution for a Halo orbit about `L`.

# Extended Help

## Arguments 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `L`: Lagrange point to orbit (L1 or L2).
- `Z`: Desired non-dimensional Z-amplitude for Halo orbit.
- `hemisphere`: Specifies northern or southern Halo orbit.
- `ϕ`: Desired Halo orbit phase.
- `steps`: Number of non-dimensional timepoints in returned state.

## Outputs
- Near-periodic initial condition `u`
- Halo orbit period `T`.
- Throws `ArgumentError` if L is not `1` or `2`.

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function richardson_halo(μ, L::Int; Z = 0.0, hemisphere = :northern, ϕ = 0.0, length = 10)
    if L == 1
        point = first(lagrange_point(μ, 1))
        γ = abs(one(μ) - μ - point)
        n = collect(1:4) .* one(μ)
        c = @. (μ + (-one(1))^n * (one(μ) - μ)γ^(n + 1)) / (γ^3 * (one(μ) - γ^(n + 1)))
    elseif L == 2
        point = first(lagrange_point(μ, 2))
        γ = abs(point - one(μ) + μ)
        n = collect(1:4) .* one(μ)
        c = @. ((-one(μ))^n * μ + (-one(μ))^n * (one(μ) - μ)γ^(n + 1)) /
           (γ^3 * (one(μ) + γ^(n + 1)))
    else
        throw(ArgumentError("Only Halo orbits about L1 or L2 are supported."))
    end

    if length < 2
        throw(
            ArgumentError(
                "The trajectory lenth must be two or greater. You provided: $length.",
            ),
        )
    end

    ωₚ = √((2 - c[2] + √((9c[2]^2 - 8c[2]))) / 2)
    k = (ωₚ^2 + 1 + 2c[2]) / (2ωₚ)

    d₁ = (3ωₚ^2 / k) * (k * (6ωₚ^2 - 1) - 2ωₚ)
    d₂ = (8ωₚ^2 / k) * (k * (11ωₚ^2 - 1) - 2ωₚ)
    a₂₁ = (3c[3] * (k^2 - 2)) / (4(1 + 2c[2]))
    a₂₂ = (3c[3]) / (4(1 + 2c[2]))
    a₂₃ = (-3c[3]ωₚ / (4k * d₁)) * (3k^3 * ωₚ - 6k * (k - ωₚ) + 4)
    a₂₄ = (-3c[3]ωₚ / (4k * d₁)) * (2 + 3k * ωₚ)
    b₂₁ = (-3c[3]ωₚ / (2d₁)) * (3k * ωₚ - 4)
    b₂₂ = -3c[3] * ωₚ / d₁
    d₂₁ = -c[3] / (2ωₚ^2)
    a₃₁ =
        (-9ωₚ / (4d₂)) * (4c[3] * (k * a₂₃ - b₂₁) + k * c[4] * (4 + k^2)) +
        ((9ωₚ^2 + 1 - c[2]) / (2d₂)) * (3c[3] * (2a₂₃ - k * b₂₁) + c[4] * (2 + 3k^2))
    a₃₂ =
        (-9ωₚ / (4d₂)) * (4c[3] * (3k * a₂₄ - b₂₂) + k * c[4]) -
        (3 / (2d₂)) * (9ωₚ^2 + 1 - c[2]) * (c[3] * (k * b₂₂ + d₂₁ - 2a₂₄) - c[4])
    b₃₁ =
        (3 / (8d₂)) * 8ωₚ * (3c[3] * (k * b₂₁ - 2a₂₃) - c[4] * (2 + 3k^2)) +
        (3 / (8d₂)) * (9ωₚ^2 + 1 + 2c[2]) * (4c[3] * (k * a₂₃ - b₂₁) + k * c[4] * (4 + k^2))
    b₃₂ =
        (9ωₚ / d₂) * (c[3] * (k * b₂₂ + d₂₁ - 2a₂₄) - c[4]) +
        (3(9ωₚ^2 + 1 + 2c[2]) / (8d₂) * (4c[3] * (k * a₂₄ - b₂₂) + k * c[4]))
    d₃₁ = (3 / (64ωₚ^2)) * (4c[3] * a₂₄ + c[4])
    d₃₂ = (3 / (64 + ωₚ^2)) * (4c[3] * (a₂₃ - d₂₁) + c[4] * (4 + k^2))

    s₁ =
        (1 / (2ωₚ * (ωₚ * (1 + k^2) - 2k))) * (
            3c[3] / 2 * (2a₂₁ * (k^2 - 2) - a₂₃ * (k^2 + 2) - 2k * b₂₁) -
            (3c[4] / 8) * (3k^4 - 8k^2 + 8)
        )
    s₂ =
        (1 / (2ωₚ * (ωₚ * (1 + k^2) - 2k))) * (
            3c[3] / 2 * (2a₂₂ * (k^2 - 2) + a₂₄ * (k^2 + 2) + 2k * b₂₂ + 5d₂₁) +
            (3c[4] / 8) * (12 - k^2)
        )
    l₁ = (-3c[3] / 2) * (2a₂₁ + a₂₃ + 5d₂₁) - (3c[4] / 8) * (12 - k^2) + 2ωₚ^2 * s₁
    l₂ = (3c[3] / 2) * (a₂₄ - 2a₂₂) + (9c[4] / 8) + 2ωₚ^2 * s₂
    Δ = ωₚ^2 - c[2]

    Aᵧ = Z / γ
    Aₓ = √((-l₂ * Aᵧ^2 - Δ) / l₁)

    ν = 1 + s₁ * Aₓ^2 + s₂ * Aᵧ^2
    T = 2π / (ωₚ * ν)
    τ = ν .* range(0, stop = T, length = length)

    if hemisphere == :northern
        m = 1
    elseif hemisphere == :southern
        m = 3
    else
        throw(ArgumentError("`hemisphere` must be `:northern` or `:southern`."))
    end

    δₘ = 2 - m
    τ₁ = @. ωₚ * τ + ϕ

    X = @. γ * (
        a₂₁ * Aₓ^2 + a₂₂ * Aᵧ^2 - Aₓ * cos(τ₁) +
        (a₂₃ * Aₓ^2 - a₂₄ * Aᵧ^2) * cos(2τ₁) +
        (a₃₁ * Aₓ^3 - a₃₂ * Aₓ * Aᵧ^2) * cos(3τ₁)
    ) + 1 - μ - (L == 1 ? γ : -γ)
    Y = @. γ * (
        k * Aₓ * sin(τ₁) +
        (b₂₁ * Aₓ^2 - b₂₂ * Aᵧ^2) * sin(2τ₁) +
        (b₃₁ * Aₓ^3 - b₃₂ * Aₓ * Aᵧ^2) * sin(3τ₁)
    )
    Z = @. γ * (
        δₘ * Aᵧ * cos(τ₁) +
        δₘ * d₂₁ * Aₓ * Aᵧ * (cos(2τ₁) - 3) +
        δₘ * (d₃₂ * Aᵧ * Aₓ^2 - d₃₁ * Aᵧ^3) * cos(3τ₁)
    )

    Ẋ = @. γ * (
        ωₚ * ν * Aₓ * sin(τ₁) - 2ωₚ * ν * (a₂₃ * Aₓ^2 - a₂₄ * Aᵧ^2) * sin(2τ₁) -
        3ωₚ * ν * (a₃₁ * Aₓ^3 - a₃₂ * Aₓ * Aᵧ^2) * sin(3τ₁)
    )
    Ẏ = @. γ * (
        ωₚ * ν * k * Aₓ * cos(τ₁) +
        2ωₚ * ν * (b₂₁ * Aₓ^2 - b₂₂ * Aᵧ^2) * cos(2τ₁) +
        3ωₚ * ν * (b₃₁ * Aₓ^3 - b₃₂ * Aₓ * Aᵧ^2) * cos(3τ₁)
    )
    Ż = @. γ * (
        -ωₚ * ν * δₘ * Aᵧ * sin(τ₁) - 2ωₚ * ν * δₘ * d₂₁ * Aₓ * Aᵧ * sin(2τ₁) -
        3ωₚ * ν * δₘ * (d₃₂ * Aᵧ * Aₓ^2 - d₃₁ * Aᵧ^2) * sin(3τ₁)
    )

    return [
        (; x = x, y = y, z = z, ẋ = ẋ, ẏ = ẏ, ż = ż, Δt = T) for
        (x, y, z, ẋ, ẏ, ż, Δt) in zip(X, Y, Z, Ẋ, Ẏ, Ż, T)
    ]

end


"""
Returns a `Vector` of `Matrix` values.
Each `Matrix` contains a 3-column nondimensional
position trajectory in the Synodic frame which
represents a Zero Velocity Curve.
"""
function zero_velocity_curves(
    r::AbstractVector,
    v::AbstractVector,
    μ::Real;
    nondimensional_range = range(-2; stop = 2, length = 1000),
)
    x = nondimensional_range
    y = nondimensional_range

    z = [
        2 * potential_energy([xs, ys, 0.0], μ) - jacobi_constant(r, v, μ) for xs in x,
        ys in y
    ]

    zvc_contour = contours(x, y, z, [0.0])

    curves = Vector{Matrix{Float64}}()

    for zvc_level in Contour.levels(zvc_contour)
        x = Float64[]
        y = Float64[]
        for zvc_line in Contour.lines(zvc_level)
            xs, ys = coordinates(zvc_line)
            if length(x) == length(y) == 0
                x = xs
                y = ys
            else
                x = vcat(x, xs)
                y = vcat(y, ys)
            end
        end

        if length(curves) == 0
            curves = [hcat(x, y)]
        else
            curves = vcat(curves, [hcat(x, y)])
        end
    end

    return curves
end

"""
Calculates the eigenvector associated with the stable manifold of a Monodromy matrix.
"""
function convergent_direction(stm::AbstractMatrix; atol = 1e-3)
    evals, evecs = eigen(stm)
    evals = filter(e -> isreal(e) && isposdef(e), evals) .|> real
    evecs =
        filter(
            x -> !isempty(x),
            map(vec -> filter(x -> all(isreal.(vec)), vec), eachcol(evecs)),
        ) .|> real

    imin = findmin(evals)[2]
    imax = findmax(evals)[2]

    if !isapprox(evals[imin] * evals[imax], 1, atol = atol)
        @warn "The dynamics appear to be ill-formed; the minimum and maximum real eigenvalues should be multiplicative inverses of one another. Product equals $(evals[imin] * evals[imax]), not 1."
    end

    return evecs[imin]
end

"""
Calculates the direction associated with the unstable manifold of a Monodromy matrix.
"""
function divergent_direction(stm::AbstractMatrix; atol = 1e-3)
    evals, evecs = eigen(stm)
    evals = filter(e -> isreal(e) && isposdef(e), evals) .|> real
    evecs =
        filter(
            x -> !isempty(x),
            map(vec -> filter(x -> all(isreal.(vec)), vec), eachcol(evecs)),
        ) .|> real

    imin = findmin(evals)[2]
    imax = findmax(evals)[2]

    if !isapprox(evals[imin] * evals[imax], 1, atol = atol)
        @warn "The dynamics appear to be ill-formed; the minimum and maximum real eigenvalues should be multiplicative inverses of one another. Product equals $(evals[imin] * evals[imax]), not 1."
    end

    return evecs[imax]
end

"""
Return the perturbation in Cartesian state space along a halo orbit onto the provided direction of the provided manifold.
"""
function perturbation(stm::AbstractMatrix, direction::AbstractVector; eps = 1e-8)
    return eps * normalize(stm * normalize(direction))
end

"""
Perturb a Cartesian state along a halo orbit onto a stable or unstable manifold.
"""
function perturb(
    u::AbstractVector,
    stm::AbstractMatrix,
    direction::AbstractVector;
    eps = 1e-8,
)
    p = similar(u)
    return perturb!(p, u, stm, direction; eps = eps)
end

"""
Perturb a Cartesian state in-place along a halo orbit onto a stable or unstable manifold.
"""
function perturb!(
    p::AbstractVector,
    u::AbstractVector,
    stm::AbstractMatrix,
    direction::AbstractVector;
    eps = 1e-8,
)
    return p .= @views(u[begin:begin+5]) + perturbation(stm, direction; eps = eps)
end

"""
Perturb halo orbit conditions onto the orbit's unstable manifold.
"""
function diverge(u::AbstractVector, stm::AbstractMatrix, Φ::AbstractMatrix; eps = 1e-8)
    p = similar(u)
    return diverge!(p, u, stm, Φ; eps = eps)
end

"""
Perturb halo orbit conditions in-place onto the orbit's unstable manifold.
"""
diverge!(
    p::AbstractVector,
    u::AbstractVector,
    stm::AbstractMatrix,
    Φ::AbstractMatrix;
    eps = 1e-8,
) = (p .= @views(u[begin:begin+5]) + perturbation(stm, divergent_direction(Φ); eps = eps))

"""
Perturb halo orbit conditions onto the orbit's unstable manifold.
"""
function converge(u::AbstractVector, stm::AbstractMatrix, Φ::AbstractMatrix; eps = 1e-8)
    p = similar(u)
    return converge!(p, u, stm, Φ; eps = eps)
end

"""
Perturb halo orbit conditions in-place onto the orbit's stable manifold.
"""
converge!(
    p::AbstractVector,
    u::AbstractVector,
    stm::AbstractMatrix,
    Φ::AbstractMatrix;
    eps = 1e-8,
) = (p .= @views(u[begin:begin+5]) + perturbation(stm, convergent_direction(Φ); eps = eps))

end