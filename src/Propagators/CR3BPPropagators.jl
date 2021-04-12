#
# Propagators for CR3BP orbits.
#

"""
Dynamics for ideal Circular Restricted Three-body numerical integration.
Compatable with `DifferentialEquations` solvers!
"""
function CR3BPTic!(∂u, u, p, t=0)
    # Cartesian state
    ∂u[1:3] .=  u[4:6]
    accel!(∂u[4:6], u[1:3], u[4:6], p.μ)
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
function CR3BPSTMTic!(∂u, u, p, t)
    # Cartesian state
    ∂u.r .=  u.v
    accel!(∂u.v, u.r, u.v, p.μ)

    # State transition matrix
    ∂Φ     = state_transition_dynamics(p.μ, u.r) * SMatrix{6,6}(transpose(hcat(u.Φ₁, u.Φ₂, u.Φ₃, u.Φ₄, u.Φ₅, u.Φ₆)))
    ∂u.Φ₁ .= ∂Φ[1,:]
    ∂u.Φ₂ .= ∂Φ[2,:]
    ∂u.Φ₃ .= ∂Φ[3,:]
    ∂u.Φ₄ .= ∂Φ[4,:]
    ∂u.Φ₅ .= ∂Φ[5,:]
    ∂u.Φ₆ .= ∂Φ[6,:]

    return nothing
end

"""
Given a `CartesianOrbit`, return a `DifferentialEquations.ODEProblem` instance
that describes the Restricted Two-body Problem.
"""
function SciMLBase.ODEProblem(orbit::NormalizedSynodicCR3BPOrbit, Δt::Real) 
    # Initial conditions
    u = ComponentVector(vcat(position_vector(orbit.state), velocity_vector(orbit.state)), Axis(r=1:3, v=4:6))

    ts = (epoch(orbit.state), Δt)
    ts = (min(ts...), max(ts...))
    p  = let μ = normalized_mass_parameter(orbit.system)
        (μ=μ, x₁=-μ, x₂=1-μ)
    end

    return ODEProblem(CR3BPTic!, u, ts, p)
end


"""
Given a `CartesianOrbit`, return a `DifferentialEquations.ODEProblem` instance
that describes the Restricted Two-body Problem.
"""
function SciMLBase.ODEProblem(orbit::NormalizedSynodicSTMCR3BPOrbit, Δt::Real) 
    # Initial conditions
    u = ComponentVector(vcat(orbit.state.cart.r, orbit.state.cart.r, [row for row ∈ eachrow(I(6))]...),
                        Axis(r=1:3, v=4:6, Φ₁=7:12, Φ₂=13:18,
                             Φ₃=19:24, Φ₄=25:30, Φ₅=31:36, Φ₆=37:42))

    ts = (epoch(orbit.state), Δt)
    ts = (min(ts...), max(ts...))
    p  = let μ = normalized_mass_parameter(orbit.system)
        (μ=μ, x₁=-μ, x₂=1-μ)
    end

    return ODEProblem(CR3BPSTMTic!, u, ts, p)
end


"""
Given a `CartesianOrbit`, return a `DifferentialEquations.ODEProblem` instance
that describes the Restricted Two-body Problem.
"""
SciMLBase.ODEProblem(orbit::NormalizedSynodicCR3BPOrbit, Δt::Unitful.Time) = ODEProblem(orbit, ustrip(timeunit(orbit.state), Δt))

"""
Given a `CartesianOrbit`, return a `DifferentialEquations.ODEProblem` instance
that describes the Restricted Two-body Problem.
"""
SciMLBase.ODEProblem(orbit::NormalizedSynodicSTMCR3BPOrbit, Δt::Unitful.Time) = ODEProblem(orbit, ustrip(timeunit(orbit.state), Δt))

"""
Uses `DifferentialEquations` solvers to propagate normalized, 
Synodic CR3BP states Δt into the future. All keyword arguments are passed 
directly to `DifferentialEquations` solvers.
"""
function propagate(orbit::NormalizedSynodicCR3BPOrbit, Δt::Real; kwargs...)

    # Set default kwargs (modified from [3])
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    # Initial conditions
    problem = ODEProblem(orbit, Δt)

    # Numerically integrate!
    sols = solve(problem; options...)

    # Return Trajectory structure
    return [
        let 
            cart = CartesianState(sols.u[i][1:3], sols.u[i][4:6], sols.t[i], Synodic;
                                  lengthunit = NormalizedLengthUnit(), timeunit = NormalizedTimeUnit())
            CircularRestrictedThreeBodyOrbit(cart, orbit.system)
        end
        for i ∈ 1:length(sols.t)
    ]
end

"""
Uses `DifferentialEquations` solvers to propagate normalized, 
Synodic CR3BP states Δt into the future. All keyword arguments are passed 
directly to `DifferentialEquations` solvers.
"""
function propagate(orbit::NormalizedSynodicSTMCR3BPOrbit, Δt::Real; kwargs...)

    # Set default kwargs (modified from [3])
    defaults = (;  reltol=1e-14, abstol=1e-14)
    options = merge(defaults, kwargs)

    # Initial conditions
    problem = ODEProblem(orbit, Δt)

    # Numerically integrate!
    sols = solve(problem; options...)

    # Return Trajectory structure
    return [
        let 
            cart = CartesianState(sols.u[i][1:3], sols.u[i][4:6], sols.t[i], Synodic;
                                  lengthunit = NormalizedLengthUnit(), timeunit = NormalizedTimeUnit())
            state = SynodicCartesianSTMState(cart, SMatrix{6,6}(transpose(hcat(sols.u[i].Φ₁, sols.u[i].Φ₂, sols.u[i].Φ₃, sols.u[i].Φ₄, sols.u[i].Φ₅, sols.u[i].Φ₆))))
            CircularRestrictedThreeBodyOrbit(state, orbit.system)
        end
        for i ∈ 1:length(sols.t)
    ]
end

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
Returns the Monodromy Matrix for a Halo orbit.
"""
function monodromy(orbit::NormalizedSynodicSTMCR3BPOrbit, T; reltol = 1e-14, abstol = 1e-14, atol = 1e-8)
    problem = ODEProblem(orbit, T)

    integrator = init(problem, Vern9(); reltol=reltol, abstol=abstol)
    solve!(integrator)

    Matrix(integrator.u[7:end]...) |> transpose |> Matrix
end

"""
Returns the Monodromy Matrix for a Halo orbit.
"""
function monodromy(orbit::NormalizedSynodicCR3BPOrbit, T; reltol = 1e-14, abstol = 1e-14, atol = 1e-8)
    problem = ODEProblem(NormalizedSynodicSTMCR3BPOrbit(orbit), T)

    integrator = init(problem, Vern9(); reltol=reltol, abstol=abstol)
    solve!(integrator)

    SMatrix{6,6}(integrator.u[7:end]...) |> transpose |> Matrix
end

"""
Returns true if a `RestrictedThreeBodySystem` is numerically periodic.
"""
function isperiodic(orbit::CircularRestrictedThreeBodyOrbit, T; reltol = 1e-14, abstol = 1e-14, atol = 1e-8)
    problem = ODEProblem(orbit, T)

    integrator = init(problem, Vern9(); reltol=reltol, abstol=abstol)
    solve!(integrator)
        
    return isapprox.((orbit.r,orbit.v), (integrator.u.r, integrator.u.v); atol = atol)
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
            CR3BPSTMTic!,
            ComponentVector(vcat(r₀, v₀, [row for row ∈ eachrow(I(6))]...),
                            Axis(r=1:3, v=4:6, Φ₁=7:12, Φ₂=13:18,
                                 Φ₃=19:24, Φ₄=25:30, Φ₅=31:36, Φ₆=37:42)),
            (0.0, τ),
            (μ = μ,)
        )    

        integrator = init(problem, Vern9(); reltol=reltol, abstol=abstol)
        solve!(integrator)

        rₛ = integrator.u.r
        vₛ = integrator.u.v

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

        if abs(integrator.u.v[1]) ≤ tolerance && abs(integrator.u.v[3]) ≤ tolerance
            break;
        elseif i == max_iter
            @warn "Desired tolerance was not reached, and iterations have hit the maximum number of iterations: $max_iter."
        end

    end

    return r₀, v₀, 2τ

end

"""
A `halo` wrapper! Returns a `CircularRestrictedThreeBodyOrbit`.
Returns a tuple: `halo_orbit, halo_period`.
"""
function halo(sys::CircularRestrictedThreeBodySystem; kwargs...)
    r,v,T = halo(normalized_mass_parameter(sys); kwargs...)
    orbit = CircularRestrictedThreeBodyOrbit(redimensionalize_length.(r, normalized_length_unit(sys)), 
                                             redimensionalize_velocity.(v, normalized_length_unit(sys), normalized_time_unit(sys)),
                                             sys) |> normalize
    return orbit, T
end

"""
A wrapper for `monodromy`.
"""
function monodromy(orb::CircularRestrictedThreeBodyOrbit, period::Real)
    state = (NormalizedSynodicSTMCR3BPOrbit ∘ synodic ∘ normalize)(orb)
    final = propagate(orbit, period) |> last
    return final.state.stm
end


