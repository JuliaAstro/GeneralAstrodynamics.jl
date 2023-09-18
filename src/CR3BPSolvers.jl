"""
Solvers specific to the Circular Restricted Three Body Problem.

# Extended Help

## Exports
$(EXPORTS)

## Imports
$(IMPORTS)
"""
module CR3BPSolvers

export halo, monodromy, converge, converge!, diverge, diverge!

using LinearAlgebra
using StaticArrays
using ModelingToolkit
using AstrodynamicalModels
using AstrodynamicalCalculations
using OrdinaryDiffEq
using DocStringExtensions

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

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)

                               $(DOCSTRING)
                               """

"""
Given a full state vector for CR3BP dynamics, including vertically concatenated
columns of the state transition matrix, return the differential correction term
for a planar periodic orbit.
"""
function planar_differential(state::AbstractVector, μ)

    f = CR3BPFunction()
    accel = f(state, (μ,), NaN)

    F = @views [
        state[22] state[34] accel[4]
        state[24] state[36] accel[6]
        state[20] state[32] state[5]
    ]

    δz, δẏ, δτ = -1 * (inv(F) * @SVector [state[4], state[6], state[2]])

    if isnan(δz) || isnan(δẏ) || isnan(δτ)
        error("The correction matrix does not have full rank, and thus cannot be applied. Try another initial condition, and if that does not work, please open an issue.")
    else
        return (; δz, δẏ, δτ)
    end

end

"""
Given a full state vector for CR3BP dynamics, including vertically concatenated
columns of the state transition matrix, return the differential correction term
for a periodic orbit.
"""
function differential(state::AbstractVector, μ)

    f = CR3BPFunction()
    accel = f(state, (μ,), NaN)

    F = @views [
        state[10] state[34] accel[4]
        state[12] state[36] accel[6]
        state[8] state[32] state[5]
    ]

    δx, δẏ, δτ = -1 * (inv(F) * @SVector [state[4], state[6], state[2]])

    if isnan(δx) || isnan(δẏ) || isnan(δτ)
        error("The correction matrix does not have full rank, and thus cannot be applied. Try another initial condition, and if that does not work, please open an issue.")
    else
        return (; δx, δẏ, δτ)
    end

end

"""
Iterate on an initial guess for Lyapunov orbit conditions.
"""
function lyapunov(x, ẏ, μ, T; reltol=1e-12, abstol=1e-12, maxiters=10)

    z = zero(x)
    y = zero(x)
    ẋ = zero(ẏ)
    ż = zero(ẏ)
    τ = T / 2

    ic = MVector{42}(x, y, z, ẋ, ẏ, ż, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1)
    p = (μ,)
    tspan = (zero(τ), τ)

    f = CR3BPFunction()
    problem = ODEProblem(CR3BPFunction(stm=true), ic, tspan, p)

    for _ in 1:maxiters

        solution = solve(problem, Vern9(), reltol=reltol, abstol=abstol, save_everystep=false)
        global fc = solution[end]

        if abs(fc[4]) <= abstol && abs(fc[6]) <= abstol
            return [fc[1], 0, 0, 0, fc[5], 0], 2 * solution.t[end]
        end

        correction = planar_differential(@views(solution[end]), μ)

        @debug "Differential Correction: $correction"

        z = @views solution[1][3] + correction.δz
        ẏ = @views solution[1][5] + correction.δẏ
        τ = @views solution.t[end] + correction.δτ

        ic = MVector{42}(x, y, z, ẋ, ẏ, ż, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1)

        tspan = (zero(τ), τ)
        _problem = remake(problem; u0=ic, tspan, p)
        problem = _problem

    end

    if abs(fc[4]) <= abstol && abs(fc[6]) <= abstol
        return (; x=problem.u0[1], ẏ=problem.u0[5], T=2solution.t[end])
    else
        error("Iterative solver failed to converge on a periodic orbit. Try another initial condition, or try relieving the provided tolerances.")
    end

end

function lyapunov(u::AbstractVector, μ, T; kwargs...)
    corrected = lyapunov(
        @views(u[begin]),
        @views(u[begin+4]),
        μ,
        T
    )

    return typeof(u)(@SVector [corrected.x, 0, 0, 0, corrected.ẏ, 0]), corrected.T
end

"""
Iterate on an initial guess for halo orbit conditions.
"""
function halo(x, z, ẏ, μ, T; reltol=1e-12, abstol=1e-12, maxiters=10)

    y = zero(x)
    ẋ = zero(ẏ)
    ż = zero(ẏ)
    τ = T / 2


    ic = MVector{42}(x, y, z, ẋ, ẏ, ż, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1)
    p = (μ,)
    tspan = (zero(τ), τ)

    f = CR3BPFunction()
    problem = ODEProblem(CR3BPFunction(stm=true), ic, tspan, p)

    for _ in 1:maxiters

        solution = solve(problem, Vern9(), reltol=reltol, abstol=abstol, save_everystep=false)
        global fc = solution[end]

        if abs(fc[4]) <= abstol && abs(fc[6]) <= abstol
            return [fc[1], 0, fc[3], 0, fc[5], 0], 2 * solution.t[end]
        end

        correction = differential(@views(solution[end]), μ)

        @debug "Differential Correction: $correction"

        x = @views solution[1][1] + correction.δx
        ẏ = @views solution[1][5] + correction.δẏ
        τ = @views solution.t[end] + correction.δτ

        ic = MVector{42}(x, y, z, ẋ, ẏ, ż, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1)

        tspan = (zero(τ), τ)
        _problem = remake(problem; u0=ic, tspan, p)
        problem = _problem

    end

    if abs(fc[4]) <= abstol && abs(fc[6]) <= abstol
        return (; x=problem.u0[1], z=problem.u0[3], ẏ=problem.u0[5], T=2solution.t[end])
    else
        error("Iterative solver failed to converge on a periodic orbit. Try another initial condition, or try relieving the provided tolerances.")
    end
end

function halo(u::AbstractVector, μ, T; kwargs...)
    corrected = halo(
        @views(u[begin]),
        @views(u[begin+2]),
        @views(u[begin+4]),
        μ,
        T
    )

    return typeof(u)(@SVector [corrected.x, 0, corrected.z, 0, corrected.ẏ, 0]), corrected.T
end

"""
Given a nondimensional mass parameter `μ`, and orbit characteristics, construct 
an initial guess using Richardson's analytical solution, and iterate on that
guess using a differential corrector. 
"""
function halo(μ, lagrange::Int; amplitude=0.0, phase=0.0, hemisphere=:northern, kwargs...)

    r, v, T = richardson_halo(μ, lagrange; Z=amplitude, hemisphere=hemisphere, ϕ=phase)

    x, y, z = r
    ẋ, ẏ, ż = v

    if z == zero(z)
        return lyapunov(x, ẏ, μ, T; kwargs...)
    else
        return halo(x, z, ẏ, μ, T; kwargs...)
    end

end

"""
Solve for the monodromy matrix of the periodic orbit.
"""
function monodromy(u::AbstractVector, μ, T; algorithm=Vern9(), reltol=1e-12, abstol=1e-12, save_everystep=false, kwargs...)
    problem = ODEProblem(CR3BPFunction(stm=true), MVector{42}(u[begin], u[begin+1], u[begin+2], u[begin+3], u[begin+4], u[begin+5], 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1), (zero(T), T), (μ,))
    solution = solve(problem, algorithm; reltol=reltol, abstol=abstol, save_everystep=save_everystep, kwargs...)

    if solution[begin][begin:begin+5] ≉ solution[end][begin:begin+5]
        @warn "The orbit does not appear to be periodic!"
    end

    return reshape((solution[end][begin+6:end]), 6, 6)
end

"""
Return conditions which diverge from the periodic orbit along its unstable manifold.
"""
function diverge(u::AbstractVector, stm::AbstractMatrix, μ, T; eps=1e-8, algorithm=Vern9(), reltol=1e-12, abstol=1e-12, save_everystep=false, kwargs...)

    Φ = monodromy(u, μ, T; algorithm=algorithm, reltol=reltol, abstol=abstol, save_everystep=save_everystep, kwargs...)
    return diverge(u, Φ; eps=eps)

end

function diverge(u::AbstractVector, stm::AbstractMatrix, Φ::AbstractMatrix; eps=1e-8)
    p = similar(u)
    return diverge!(p, u, stm, Φ; eps=eps)
end

"""
Perturb the orbit in-place along its unstable manifold.
"""
function diverge!(p::AbstractVector, u::AbstractVector, stm::AbstractMatrix, μ, T; eps=1e-8, algorithm=Vern9(), reltol=1e-12, abstol=1e-12, save_everystep=false, kwargs...)

    Φ = monodromy(u, μ, T; algorithm=algorithm, reltol=reltol, abstol=abstol, save_everystep=save_everystep, kwargs...)
    return diverge!(p, u, stm, Φ; eps=eps)

end

diverge!(p::AbstractVector, u::AbstractVector, stm::AbstractMatrix, Φ::AbstractMatrix; eps=1e-8) = (p .= @views(u[begin:begin+5]) + perturbation(Φ, divergent_direction(stm); eps=eps))


"""
Return conditions which converts to the periodic orbit along its stable manifold.
"""
function converge(u::AbstractVector, stm::AbstractMatrix, μ, T; eps=1e-8, algorithm=Vern9(), reltol=1e-12, abstol=1e-12, save_everystep=false, kwargs...)

    Φ = monodromy(u, μ, T; algorithm=algorithm, reltol=reltol, abstol=abstol, save_everystep=save_everystep, kwargs...)
    return converge(u, Φ; eps=eps)

end

function converge(u::AbstractVector, stm::AbstractMatrix, Φ::AbstractMatrix; eps=1e-8)
    p = similar(u)
    return converge!(p, u, stm, Φ; eps=eps)
end

"""
Perturb the orbit in-place along its stable manifold.
"""
function converge!(p::AbstractVector, u::AbstractVector, stm::AbstractMatrix, μ, T; eps=1e-8, algorithm=Vern9(), reltol=1e-12, abstol=1e-12, save_everystep=false, kwargs...)

    Φ = monodromy(u, μ, T; algorithm=algorithm, reltol=reltol, abstol=abstol, save_everystep=save_everystep, kwargs...)
    return converge!(p, u, Φ; eps=eps)

end

converge!(p::AbstractVector, u::AbstractVector, stm::AbstractMatrix, Φ::AbstractMatrix; eps=1e-8) = (p .= @views(u[begin:begin+5]) + perturbation(Φ, convergent_direction(Φ); eps=eps))

end # module