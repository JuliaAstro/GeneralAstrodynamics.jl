"""
Solvers specific to the Circular Restricted Three Body Problem.
"""
module CR3BPSolvers

export halo

using LinearAlgebra
using StaticArrays
using ModelingToolkit
using AstrodynamicalModels
using AstrodynamicalCalculations
using OrdinaryDiffEq
using DocStringExtensions

@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(SIGNATURES)

                                         $(DOCSTRING)
                                         """

@template (TYPES, CONSTANTS) = """
                               $(TYPEDEF)

                               $(DOCSTRING)
                               """

"""
Iterate on an initial guess for Lyapunov orbit conditions.
"""
function lyapunov_corrector(x, ẏ, μ, T; reltol=1e-12, abstol=1e-12, maxiters=10)

    z = zero(x)
    y = zero(x)
    ẋ = zero(ẏ)
    ż = zero(ẏ)
    τ = T / 2

    ic = MVector{42}(x, y, z, ẋ, ẏ, ż, reduce(vcat, eachcol(I(6)))...)
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

        fv = f(fc, p, NaN)

        F = @SMatrix [
            fc[22] fc[34] fv[4]
            fc[24] fc[36] fv[6]
            fc[20] fc[32] fc[5]
        ]

        if det(F) ≉ 0
            z, ẏ, τ = SVector(solution[1][3], solution[1][5], solution.t[end]) - inv(F) * SVector(fc[4], fc[6], fc[2])

            ic = MVector{42}(x, y, z, ẋ, ẏ, ż, reduce(vcat, eachcol(I(6)))...)

            tspan = (zero(τ), τ)
            _problem = remake(problem; u0=ic, tspan, p)
            problem = _problem
        else
            error("The correction matrix does not have full rank, and thus cannot be applied. Try another initial condition, and if that does not work, please open an issue.")
        end
    end

    if abs(fc[4]) <= abstol && abs(fc[6]) <= abstol
        return [fc[1], 0, 0, 0, fc[5], 0], 2 * solution.t[end]
    else
        error("Iterative solver failed to converge on a periodic orbit. Try another initial condition, or try relieving the provided tolerances.")
    end

end # lyapunov_corrector

"""
Iterate on an initial guess for halo orbit conditions.
"""
function halo_corrector(x, z, ẏ, μ, T; reltol=1e-12, abstol=1e-12, maxiters=10)

    y = zero(x)
    ẋ = zero(ẏ)
    ż = zero(ẏ)
    τ = T / 2


    ic = MVector{42}(x, y, z, ẋ, ẏ, ż, reduce(vcat, eachcol(I(6)))...)
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

        fv = f(fc, p, NaN)

        F = @SMatrix [
            fc[10] fc[34] fv[4]
            fc[12] fc[36] fv[6]
            fc[8] fc[32] fc[5]
        ]

        if det(F) ≉ 0
            x, ẏ, τ = SVector(solution[1][1], solution[1][5], solution.t[end]) - inv(F) * SVector(fc[4], fc[6], fc[2])

            ic = MVector{42}(x, y, z, ẋ, ẏ, ż, reduce(vcat, eachcol(I(6)))...)

            if any(isnan.((z, ẏ, τ)))
                error("Corrector returned NaN values.")
            end

            tspan = (zero(τ), τ)
            _problem = remake(problem; u0=ic, tspan, p)
            problem = _problem
        else
            error("Iterative solver failed to converge on a periodic orbit. Try another initial condition, or try relieving the provided tolerances.")
        end
    end

    if abs(fc[4]) <= abstol && abs(fc[6]) <= abstol
        return [fc[1], 0, fc[3], 0, fc[5], 0], 2 * solution.t[end]
    else
        error("Iterative solver failed to converge on a periodic orbit. Try another initial condition, or try relieving the provided tolerances.")
    end
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
        return lyapunov_corrector(x, ẏ, μ, T; kwargs...)
    else
        return halo_corrector(x, z, ẏ, μ, T; kwargs...)
    end

end

end # module