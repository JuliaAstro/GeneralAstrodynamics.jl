"""
Unit tests for various solvers.
"""
module AstrodynamicalSolversTests

using ParallelTestRunner: runtests, find_tests, parse_args

const init_code = quote
    using Test
    using AstrodynamicalSolvers:
        AstrodynamicalSolvers,
        halo,
        monodromy,
        propagate,
        propagate!
    using AstrodynamicalCalculations: converge, diverge, richardson_ic
    using OrdinaryDiffEqVerner: ODEProblem, ODESolution, Vern9, solve
    using AstrodynamicalModels: CR3BFunction, CR3BSystem, CartesianState, Orbit, R2BParameters
    using LinearAlgebra: I
    using ModelingToolkit: complete
end

args = parse_args(Base.ARGS)
testsuite = find_tests(@__DIR__)

runtests(AstrodynamicalSolvers, args; testsuite, init_code)

end # module
