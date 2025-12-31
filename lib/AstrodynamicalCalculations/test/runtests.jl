"""
Restricted Two-body Model tests.
"""
module AstrodynamicalCalculationsTests

using ParallelTestRunner: runtests, find_tests, parse_args
import AstrodynamicalCalculations

const init_code = quote
    # TODO: Use explicit imports in each test file after
    # # https://github.com/JuliaAstro/GeneralAstrodynamics.jl/pull/280
    # # is in
    using AstrodynamicalCalculations
    using LinearAlgebra
    using StaticArrays
end

args = parse_args(Base.ARGS)
testsuite = find_tests(@__DIR__)

runtests(AstrodynamicalCalculations, args; testsuite, init_code)

end # module
