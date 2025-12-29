"""
Solvers specific to the Circular Restricted Three Body Problem.

# Extended Help

## Exports
$(EXPORTS)

## Imports
$(IMPORTS)
"""
module CR3BSolvers

export halo, lyapunov


using StaticArrays: @SMatrix, @SVector, SVector
using AstrodynamicalModels: CR3BFunction
using AstrodynamicalCalculations: richardson_ic
using OrdinaryDiffEqVerner: ODEProblem, Vern9, remake, solve
using DocStringExtensions: @template, DOCSTRING, EXPORTS, IMPORTS, LICENSE, SIGNATURES, TYPEDEF

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

    f = CR3BFunction()
    accel = f(state, SVector(μ), NaN)

    F = @SMatrix [
        state[22] state[34] accel[4]
        state[24] state[36] accel[6]
        state[20] state[32] state[5]
    ]

    δz, δẏ, δτ = -1 * (inv(F) * @SVector [state[4], state[6], state[2]])

    if isnan(δz) || isnan(δẏ) || isnan(δτ)
        error(
            "The correction matrix does not have full rank, and thus cannot be applied. Try another initial condition, and if that does not work, please open an issue.",
        )
    else
        return (; δz, δẏ, δτ)
    end

end

"""
Given a full state vector for CR3BP dynamics, including vertically concatenated
columns of the state transition matrix, return the differential correction term
for a periodic orbit.
"""
function extraplanar_differential(state::AbstractVector, μ)

    f = CR3BFunction()
    accel = f(state, SVector(μ), NaN)

    F = [
        state[10] state[34] accel[4]
        state[12] state[36] accel[6]
        state[8] state[32] state[5]
    ]

    δx, δẏ, δτ = -1 * (inv(F) * @SVector [state[4], state[6], state[2]])

    if isnan(δx) || isnan(δẏ) || isnan(δτ)
        error(
            "The correction matrix does not have full rank, and thus cannot be applied. Try another initial condition, and if that does not work, please open an issue.",
        )
    else
        return (; δx, δẏ, δτ)
    end

end

"""
Iterate on an initial guess for Lyapunov orbit conditions.
"""
function lyapunov(x, ẏ, μ, T; reltol = 1e-12, abstol = 1e-12, maxiters = 10)

    z = zero(x)
    y = zero(x)
    ẋ = zero(ẏ)
    ż = zero(ẏ)
    τ = T / 2

    ic = [
        x,
        y,
        z,
        ẋ,
        ẏ,
        ż,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
    ]
    p = @SVector [μ]
    tspan = (zero(τ), τ)

    problem = ODEProblem(CR3BFunction(stm = true), ic, tspan, p)

    for _ = 1:maxiters

        solution = solve(problem, Vern9();
            reltol,
            abstol,
            save_everystep = false,
        )
        global fc = solution.u[end]

        if abs(fc[4]) <= abstol && abs(fc[6]) <= abstol
            return (; x = x, ẏ = ẏ, Δt = 2τ)
        end

        correction = planar_differential(solution.u[end], μ)

        _z = z + correction.δz
        _ẏ = ẏ + correction.δẏ
        _τ = τ + correction.δτ

        z = _z
        ẏ = _ẏ
        τ = _τ

        ic = [
            x,
            y,
            z,
            ẋ,
            ẏ,
            ż,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
        ]

        tspan = (zero(τ), τ)
        _problem = remake(problem; u0 = ic, tspan, p)
        problem = _problem

    end

    if abs(fc[4]) <= abstol && abs(fc[6]) <= abstol
        return (; x = x, ẏ = ẏ, Δt = 2τ)
    else
        error(
            "Iterative solver failed to converge on a periodic orbit. Try another initial condition, or try relieving the provided tolerances.",
        )
    end

end

function lyapunov(u::AbstractVector, μ, T; kwargs...)
    corrected = lyapunov(u[begin], u[begin+4], μ, T)

    return typeof(u)(@SVector [corrected.x, 0, 0, 0, corrected.ẏ, 0]), corrected.T
end

"""
Iterate on an initial guess for halo orbit conditions.
"""
function halo(x, z, ẏ, μ, T; reltol = 1e-12, abstol = 1e-12, maxiters = 10)

    y = zero(x)
    ẋ = zero(ẏ)
    ż = zero(ẏ)
    τ = T / 2


    ic = [
        x,
        y,
        z,
        ẋ,
        ẏ,
        ż,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
    ]

    p = SVector(μ)
    tspan = (zero(τ), τ)

    f = CR3BFunction()
    problem = ODEProblem(CR3BFunction(stm = true), ic, tspan, p)

    for _ = 1:maxiters

        solution = solve(problem, Vern9();
            reltol,
            abstol,
            save_everystep = false,
        )
        global fc = solution.u[end]

        if abs(fc[4]) <= abstol && abs(fc[6]) <= abstol
            return (; x = x, z = z, ẏ = ẏ, Δt = 2τ)
        end

        correction = extraplanar_differential(@views(solution.u[end]), μ)
        _x = x + correction.δx
        _ẏ = ẏ + correction.δẏ
        _τ = τ + correction.δτ

        x = _x
        ẏ = _ẏ
        τ = _τ

        @debug """

               Differential Correction: (x = $x, ẏ = $ẏ, τ = $τ)

                New Initial Conditions: $correction


               """

        ic = [
            x,
            y,
            z,
            ẋ,
            ẏ,
            ż,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            1,
        ]

        tspan = (zero(τ), τ)
        _problem = remake(problem; u0 = ic, tspan, p)
        problem = _problem

    end

    if abs(fc[4]) <= abstol && abs(fc[6]) <= abstol
        return (; x = x, z = z, ẏ = ẏ, Δt = 2τ)
    else
        error(
            "Iterative solver failed to converge on a periodic orbit. Try another initial condition, or try relieving the provided tolerances.",
        )
    end
end

function halo(u::AbstractVector, μ, T; kwargs...)
    corrected = halo(u[begin], u[begin+2], u[begin+4], μ, T)

    return typeof(u)(@SVector [corrected.x, 0, corrected.z, 0, corrected.ẏ, 0]),
    corrected.Δt

end

"""
Given a nondimensional mass parameter `μ`, and orbit characteristics, construct
an initial guess using Richardson's analytical solution, and iterate on that
guess using a differential corrector.
"""
function halo(μ, lagrange::Int;
    amplitude = 0.0,
    phase = 0.0,
    hemisphere = :northern,
    kwargs...,
)

    u = richardson_ic(μ, lagrange; Z = amplitude, hemisphere, ϕ = phase)

    if u.z == 0
        return lyapunov(u.x, u.ẏ, μ, u.Δt; kwargs...)
    else
        return halo(u.x, u.z, u.ẏ, μ, u.Δt; kwargs...)
    end

end

end # module
