"""
Unit tests for AstrodynamicalModels.jl
"""
module AstrodynamicalModelsTests

using ParallelTestRunner: runtests, find_tests, parse_args
import AstrodynamicalModels

const init_code = quote
end

args = parse_args(Base.ARGS)
testsuite = find_tests(@__DIR__)

runtests(AstrodynamicalModels, args; testsuite, init_code)

end
