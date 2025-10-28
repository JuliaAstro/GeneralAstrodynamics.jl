"""
Unit tests for various solvers.
"""
module AstrodynamicalSolversTests

using ParallelTestRunner: runtests, find_tests, parse_args
import AstrodynamicalSolvers

const init_code = quote
    using AstrodynamicalCalculations
    using AstrodynamicalModels
    using LinearAlgebra
    using OrdinaryDiffEqVerner
    using StaticArrays
end

args = parse_args(Base.ARGS)
testsuite = find_tests(@__DIR__)

runtests(AstrodynamicalSolvers, args; testsuite, init_code)

end # module
