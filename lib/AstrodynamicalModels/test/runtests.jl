"""
Unit tests for AstrodynamicalModels.jl
"""
module AstrodynamicalModelsTests

using ParallelTestRunner: runtests, find_tests, parse_args
import AstrodynamicalModels

const init_code = quote
    # TODO: use explicit imports in each test file after
    # https://github.com/JuliaAstro/GeneralAstrodynamics.jl/pull/280
    # is in
    using AstrodynamicalModels, ModelingToolkit, LinearAlgebra, Test
end

args = parse_args(Base.ARGS)
testsuite = find_tests(@__DIR__)

runtests(AstrodynamicalModels, args; testsuite, init_code)

end # module
