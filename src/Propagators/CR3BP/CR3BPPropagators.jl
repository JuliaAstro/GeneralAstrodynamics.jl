#
# Propagators for CR3BP orbits.
#

"""
Dynamics for ideal Circular Restricted Three-body numerical integration.
Compatable with `DifferentialEquations` solvers!
"""
function CR3BPTic!(∂u, u, p, t=0)
    ∂u.r =  u.v
    accel!(∂u.v, u.r, u.v, p.μ)
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
potential_energy_hessian = include("PotentialEnergyHessian.jl")

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
        hcat(potential_energy_hessian(r..., μ), [0 2 0; -2 0 0; 0 0 0])
    ))

end