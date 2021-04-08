#
# Functions related to Halo orbits
#

"""
Returns the partial derivative matrix of potential `U`.

__Arguments:__
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `r`: Non-dimensional position vector for the spacecraft.

__Outputs:__
- Partial derivative matrix of potential `U`.

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/)
"""
potential_energy_hessian = let
    func = include("PotentialEnergyHessian.jl")
    (r,μ) -> func(r..., μ)
end

"""
Returns the derivative mapping of CR3BP state transition matrix, `F`.

__Arguments:__
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `r`: Non-dimensional position vector for the spacecraft.

__Outputs:__
- Linear mapping from Φ to Φ̇, `F`.

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/)
"""
function state_transition_dynamics(μ, r)

    return SMatrix{6,6}(vcat(
        hcat(zeros((3,3)), I(3)),
        hcat(potential_energy_hessian(r, μ), [0 2 0; -2 0 0; 0 0 0])
    ))

end

"""
Returns an analytical solution for a Halo orbit about `L`.

__Arguments:__ 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `Az`: Desired non-dimensional Z-amplitude for Halo orbit.
- `ϕ`: Desired Halo orbit phase.
- `steps`: Number of non-dimensional timepoints in returned state.
- `L`: Lagrange point to orbit (L1 or L2).
- `hemisphere`: Specifies northern or southern Halo orbit.

__Outputs:__
- Synodic position vector `r::Array{<:AbstractFloat}`
- Synodic velocity vector `v::Array{<:Abstractfloat}`.
- Halo orbit period `Τ`.
- Throws `ArgumentError` if L is not `1` or `2`.

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function analyticalhalo(μ; Az=0.00, ϕ=0.0, steps=1,
                       L=1, hemisphere=:northern)

    if L == 1
        point = first(lagrange(μ, 1))
        γ = abs(1 - μ - point)
        n = collect(1:4)
        c = @. (μ + (-1)^n * (1-μ)γ^(n+1)) / (γ^3 * (1 - γ^(n+1)))
    elseif L == 2
        point = first(lagrange(μ, 2))
        γ = abs(point - 1 + μ)
        n = collect(1:4)
        c = @. ((-1)^n * μ + (-1)^n * (1-μ)γ^(n+1)) / (γ^3 * (1 + γ^(n+1)))
    else
        throw(ArgumentError("Only Halo orbits about L1 or L2 are supported."))
    end

    ωₚ  = √((2 - c[2] + √((9c[2]^2 - 8c[2])))/2)
    k   = (ωₚ^2 + 1 + 2c[2]) / (2ωₚ)

    d₁  = (3ωₚ^2 / k) * (k*(6ωₚ^2 - 1) - 2ωₚ)
    d₂  = (8ωₚ^2 / k) * (k*(11ωₚ^2 - 1) - 2ωₚ)
    a₂₁ = (3c[3] * (k^2 - 2)) / (4(1 + 2c[2]))
    a₂₂ = (3c[3]) / (4(1 + 2c[2]))
    a₂₃ = (-3c[3]ωₚ / (4k*d₁)) * (3k^3 * ωₚ - 6k*(k-ωₚ) + 4)
    a₂₄ = (-3c[3]ωₚ / (4k*d₁)) * (2 + 3k*ωₚ)
    b₂₁ = (-3c[3]ωₚ / (2d₁)) * (3k*ωₚ - 4)
    b₂₂ = -3c[3]*ωₚ / d₁
    d₂₁ = -c[3] / (2ωₚ^2)
    a₃₁ = (-9ωₚ / (4d₂)) * (4c[3] * (k*a₂₃-b₂₁) + k*c[4]*(4+k^2)) + 
          ((9ωₚ^2 + 1 - c[2]) / (2d₂)) * (3c[3]*(2a₂₃-k*b₂₁) + c[4]*(2+3k^2))
    a₃₂ = (-9ωₚ / (4d₂)) * (4c[3] * (3k*a₂₄-b₂₂) + k*c[4]) -
          (3 / (2d₂)) * (9ωₚ^2 + 1 - c[2]) * (c[3]*(k*b₂₂+d₂₁-2a₂₄) - c[4])
    b₃₁ = (3 / (8d₂)) * 8ωₚ * (3c[3] * (k*b₂₁ - 2a₂₃) - c[4]*(2+3k^2)) + 
          (3/(8d₂)) * (9ωₚ^2 + 1 + 2c[2]) * (4c[3]*(k*a₂₃-b₂₁) + k*c[4]*(4+k^2))
    b₃₂ = (9ωₚ/d₂)*(c[3]*(k*b₂₂+d₂₁-2a₂₄)-c[4]) + 
          (3(9ωₚ^2 + 1 + 2c[2]) / (8d₂) * (4c[3]*(k*a₂₄-b₂₂)+k*c[4]))
    d₃₁ = (3 / (64ωₚ^2)) * (4c[3]*a₂₄ + c[4])
    d₃₂ = (3 / (64 + ωₚ^2)) * (4c[3]*(a₂₃ - d₂₁) + c[4]*(4+k^2))

    s₁  = (1 / (2ωₚ*(ωₚ*(1+k^2) - 2k))) * 
          (3c[3]/2 * (2a₂₁*(k^2 - 2) - a₂₃*(k^2 + 2) - 2k*b₂₁) - 
            (3c[4]/8) * (3k^4 - 8k^2 + 8))
    s₂  = (1 / (2ωₚ*(ωₚ*(1+k^2) - 2k))) * 
          (3c[3]/2 * (2a₂₂*(k^2-2) + a₂₄*(k^2 + 2) + 2k*b₂₂ + 5d₂₁) + 
            (3c[4]/8) * (12 - k^2))
    l₁  = (-3c[3] / 2) * (2a₂₁ + a₂₃ + 5d₂₁) - (3c[4]/8)*(12 - k^2) + 2ωₚ^2 * s₁
    l₂  = (3c[3]/2) * (a₂₄ - 2a₂₂) + (9c[4]/8) + 2ωₚ^2 * s₂
    Δ   = ωₚ^2 - c[2]

    Aᵧ  = Az / γ
    Aₓ  = √((-l₂*Aᵧ^2 - Δ) / l₁)

    ν   = 1 + s₁*Aₓ^2 + s₂*Aᵧ^2
    Τ   = 2π / (ωₚ*ν)
    τ   = ν .* (steps > 1 ? range(0, stop=Τ, length=steps) : range(0, stop=Τ, length=1000))

    if hemisphere == :northern
        m = 1.0
    elseif hemisphere == :southern
        m = 3.0
    else
        throw(ArgumentError("`hemisphere` must be `:northern` or `:southern`."))
    end

    δₘ = 2 - m
    τ₁ = @. ωₚ*τ + ϕ

    x = @. γ * (a₂₁*Aₓ^2 + a₂₂*Aᵧ^2 - Aₓ*cos(τ₁) + (a₂₃*Aₓ^2 - 
                    a₂₄*Aᵧ^2)*cos(2τ₁) + (a₃₁*Aₓ^3 - a₃₂*Aₓ*Aᵧ^2)*cos(3τ₁)) + 1 - μ - (L == 1 ? γ : -γ)
    y = @. γ * (k*Aₓ*sin(τ₁) + (b₂₁*Aₓ^2 - b₂₂*Aᵧ^2)*sin(2τ₁) + 
                    (b₃₁*Aₓ^3 - b₃₂*Aₓ*Aᵧ^2)*sin(3τ₁))
    z = @. γ * (δₘ*Aᵧ*cos(τ₁) + δₘ*d₂₁*Aₓ*Aᵧ*(cos(2τ₁)-3) + 
                    δₘ*(d₃₂*Aᵧ*Aₓ^2 - d₃₁*Aᵧ^3)*cos(3τ₁))

    ẋ = @. γ * (ωₚ*ν*Aₓ*sin(τ₁) - 2ωₚ*ν*(a₂₃*Aₓ^2 - a₂₄*Aᵧ^2)*sin(2τ₁) - 
                    3ωₚ*ν*(a₃₁*Aₓ^3 - a₃₂*Aₓ*Aᵧ^2)*sin(3τ₁))
    ẏ = @. γ * (ωₚ*ν*k*Aₓ*cos(τ₁) + 2ωₚ*ν*(b₂₁*Aₓ^2 - b₂₂*Aᵧ^2)*cos(2τ₁) + 
                    3ωₚ*ν*(b₃₁*Aₓ^3 - b₃₂*Aₓ*Aᵧ^2)*cos(3τ₁))
    ż = @. γ * (-ωₚ*ν*δₘ*Aᵧ*sin(τ₁) - 2ωₚ*ν*δₘ*d₂₁*Aₓ*Aᵧ*sin(2τ₁) - 
                    3ωₚ*ν*δₘ*(d₃₂*Aᵧ*Aₓ^2 - d₃₁*Aᵧ^2)*sin(3τ₁))

    return hcat(x, y, z)[1:steps, :], hcat(ẋ, ẏ, ż)[1:steps, :], Τ

end

"""
Returns a numerical solution for a Halo orbit about `L`.

__Arguments:__ 
- `μ`: Non-dimensional mass parameter for the CR3BP system.
- `Az`: Desired non-dimensional Z-amplitude for Halo orbit.
- `ϕ`: Desired Halo orbit phase.
- `L`: Lagrange point to orbit (L1 or L2).
- `hemisphere`: Specifies northern or southern Halo orbit.

__Outputs:__
- Tuple of initial states: `(r, v)` where `r::Vector{<:AbstractFloat}`, `v::Vector{<:Abstractfloat}`.
- Throws `ArgumentError` if L is not `:L1` or `:L2`

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function halo(μ; Az=0.0, L=1, hemisphere=:northern,
              tolerance=1e-8, max_iter=20,
              reltol=1e-14, abstol=1e-14)

    r₀, v₀, Τ = analyticalhalo(μ; Az=Az, ϕ=0.0, L=L, hemisphere=hemisphere)
    r₀ = r₀[1,:]
    v₀ = v₀[1,:]
    τ  = Τ/2
    
    Φ  = Matrix{promote_type(eltype(r₀), eltype(v₀), typeof(τ))}(undef, 6, 6)
    
    for i ∈ 1:max_iter
    
        problem = ODEProblem(
            RestrictedThreeBodySTMTic!,
            ComponentArray(rₛ  = r₀,
                           vₛ  = v₀,
                           Φ₁  = [1.0, 0, 0, 0, 0, 0],
                           Φ₂  = [0, 1.0, 0, 0, 0, 0],
                           Φ₃  = [0, 0, 1.0, 0, 0, 0],
                           Φ₄  = [0, 0, 0, 1.0, 0, 0],
                           Φ₅  = [0, 0, 0, 0, 1.0, 0],
                           Φ₆  = [0, 0, 0, 0, 0, 1.0]),
            (0.0, τ),
            ComponentArray(μ   =  μ)
        )    
    
        integrator = init(problem, Vern9(); reltol=reltol, abstol=abstol)
        solve!(integrator)
    
        rₛ = integrator.u.rₛ
        vₛ = integrator.u.vₛ
    
        Φ = hcat(integrator.u.Φ₁, integrator.u.Φ₂, integrator.u.Φ₃, integrator.u.Φ₄, integrator.u.Φ₅, integrator.u.Φ₆) |> transpose
    
        ∂vₛ = accel(rₛ, vₛ, μ)
    
        if Az ≉ 0
            F = @SMatrix [
                Φ[4,1] Φ[4,5] ∂vₛ[1];
                Φ[6,1] Φ[6,5] ∂vₛ[3];
                Φ[2,1] Φ[2,5]  vₛ[2]
            ]
        
            TERM1 = @SMatrix [r₀[1]; v₀[2]; τ] 
            TERM2 = - inv(F) * @SMatrix [vₛ[1]; vₛ[3]; rₛ[2]] 
            xᵪ = TERM1 + TERM2
        
            r₀[1] = xᵪ[1]
            v₀[2] = xᵪ[2]
            τ     = xᵪ[3]
        else
            F = @SMatrix [
                Φ[4,3] Φ[4,5] ∂vₛ[1];
                Φ[6,3] Φ[6,5] ∂vₛ[3];
                Φ[2,3] Φ[2,5]  vₛ[2]
            ]
        
            TERM1 = @SMatrix [r₀[3]; v₀[2]; τ] 
            TERM2 = - inv(F) * @SMatrix [vₛ[1]; vₛ[3]; rₛ[2]] 
            xᵪ = TERM1 + TERM2
        
            r₀[3] = xᵪ[1]
            v₀[2] = xᵪ[2]
            τ     = xᵪ[3]
        end

        if abs(integrator.u.vₛ[1]) ≤ tolerance && abs(integrator.u.vₛ[3]) ≤ tolerance
            break;
        elseif i == max_iter
            @warn "Desired tolerance was not reached, and iterations have hit the maximum number of iterations: $max_iter."
        end
    end

    return r₀, v₀, 2τ

end

"""
Returns a `NondimensionalThreeBodyState` type, instead of tuple `r,v,T`.
"""
function halo(μ, DU::T1, DT::T2; kwargs...) where T1 <: Unitful.Length where T2 <: Time
    r, v, T = halo(μ; kwargs...)
    return NondimensionalThreeBodyState(r, v, μ, T, DU, DT)
end

"""
Iterative halo solver, returns a `NondimensionalThreeBodyState` in-place.
"""
function halo!(state, μ; kwargs...) 
    r,v,T = halo(μ; kwargs...)
    state = NondimensionalThreeBodyState(r, v, μ, T)
    return nothing
end

"""
Returns dynamics tic for combined Halo iterative solver state vector.

__Arguments:__ 
- `∂u`: Derivative of state `u`.
- `u`: State vector: `[x, y, z, ẋ, ẏ, ż, Φ₁, Φ₂, Φ₃, Φ₄, Φ₅, Φ₆]`.
- `p`: Parameters (contains non-dimensional mass parameter `μ`, positions `x₁`, `x₂`, and configuration).
- `t`: Time in seconds.

__Outputs:__
- None (sets derivative `∂u` in-place).

__References:__
- [Rund, 2018](https://digitalcommons.calpoly.edu/theses/1853/).
"""
function RestrictedThreeBodySTMTic!(∂u, u, p, t)

    # Cartesian state
    ∂u.rₛ =  u.vₛ
    accel!(∂u.vₛ, u.rₛ, u.vₛ, p.μ)

    # State transition matrix
    ∂Φ  = state_transition_dynamics(p.μ, u.rₛ) * SMatrix{6,6}(transpose(hcat(u.Φ₁, u.Φ₂, u.Φ₃, u.Φ₄, u.Φ₅, u.Φ₆)))
    ∂u.Φ₁ = copy(∂Φ[1,:])[:]
    ∂u.Φ₂ = copy(∂Φ[2,:])[:]
    ∂u.Φ₃ = copy(∂Φ[3,:])[:]
    ∂u.Φ₄ = copy(∂Φ[4,:])[:]
    ∂u.Φ₅ = copy(∂Φ[5,:])[:]
    ∂u.Φ₆ = copy(∂Φ[6,:])[:]
    
end

"""
Returns the Monodromy Matrix for a Halo orbit.
"""
function monodromy(orbit::NondimensionalThreeBodyState; check_periodicity = true, reltol = 1e-14, abstol = 1e-14, atol = 1e-8)
    problem = ODEProblem(
        RestrictedThreeBodySTMTic!,
        ComponentArray(rₛ  = orbit.r,
                       vₛ  = orbit.v,
                       Φ₁  = [1.0, 0, 0, 0, 0, 0],
                       Φ₂  = [0, 1.0, 0, 0, 0, 0],
                       Φ₃  = [0, 0, 1.0, 0, 0, 0],
                       Φ₄  = [0, 0, 0, 1.0, 0, 0],
                       Φ₅  = [0, 0, 0, 0, 1.0, 0],
                       Φ₆  = [0, 0, 0, 0, 0, 1.0]),
        (0.0, orbit.Δt),
        ComponentArray(μ   =  orbit.μ)
    )

    u = solve(problem; reltol = reltol, abstol = abstol).u[end]
    
    if check_periodicity
        if !isapprox(orbit, NondimensionalThreeBodyState(u.rₛ, u.vₛ, orbit.μ, orbit.Δt, orbit.DU, orbit.DT); atol = atol)
            throw(ErrorException("Provided CR3BP system is NOT periodic!"))
        end
    end

    Matrix(transpose(hcat(u.Φ₁, u.Φ₂, u.Φ₃, u.Φ₄, u.Φ₅, u.Φ₆)))
end

"""
Returns true if a `RestrictedThreeBodySystem` is numerically periodic.
"""
function isperiodic(orbit::NondimensionalThreeBodyState; reltol = 1e-14, abstol = 1e-14, atol = 1e-8)
    problem = ODEProblem(
        RestrictedThreeBodySTMTic!,
        ComponentArray(rₛ  = orbit.r,
                       vₛ  = orbit.v,
                       Φ₁  = [1.0, 0, 0, 0, 0, 0],
                       Φ₂  = [0, 1.0, 0, 0, 0, 0],
                       Φ₃  = [0, 0, 1.0, 0, 0, 0],
                       Φ₄  = [0, 0, 0, 1.0, 0, 0],
                       Φ₅  = [0, 0, 0, 0, 1.0, 0],
                       Φ₆  = [0, 0, 0, 0, 0, 1.0]),
        (0.0, orbit.Δt),
        ComponentArray(μ   =  orbit.μ)
    )

    u = solve(problem; reltol = reltol, abstol = abstol).u[end]
    
    return isapprox(orbit, NondimensionalThreeBodyState(u.rₛ, u.vₛ, orbit.μ, orbit.Δt, orbit.DU, orbit.DT); atol = 1e-8)
end