"""
Restricted Two-body Model tests.
"""
module AstrodynamicalCalculationsTests

using ParallelTestRunner: runtests, find_tests, parse_args
import AstrodynamicalCalculations

const init_code = quote
    using Test
    using AstrodynamicalCalculations:
        #CR3BP
        converge,
        diverge,
        jacobi_constant,
        richardson_ic,
        richardson_halo,

        # R2BP
        R2BPCalculations,
        argument_of_periapsis,
        cartesian_to_keplerian,
        cartesian_to_perifocal,
        conic,
        inclination,
        kepler,
        keplerian_to_cartesian,
        lambert,
        orbital_period,
        orbital_radius,
        right_ascension_ascending_node,
        specific_energy,
        specific_angular_momentum,
        specific_angular_momentum_vector,
        semimajor_axis
    using LinearAlgebra: norm
    using StaticArrays: SVector
end

args = parse_args(Base.ARGS)
testsuite = find_tests(@__DIR__)

runtests(AstrodynamicalCalculations, args; testsuite, init_code)

end # module
