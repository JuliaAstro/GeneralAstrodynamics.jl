"""
Unit tests for various solvers.
"""
module AstrodynamicalSolversTests

using Test, AstrodynamicalSolvers, AstrodynamicalCalculations
using StaticArrays

@testset "Lyapunov Orbit Correction" begin

    μ = 0.012150584395829193
    r, v, T = richardson_halo(μ, 1)
    x, _, _ = r
    _, ẏ, _ = v
    u, T = AstrodynamicalSolvers.CR3BPSolvers.lyapunov(x, ẏ, μ, T)

    @test u ≈ [0.8567678285004178, 0.0, 0.0, 0.0, -0.14693135696819282, 0.0]
    @test T ≈ 2.7536820160579087

end

@testset "Halo Orbit Correction" begin

    μ = 0.012150584395829193
    r, v, T = richardson_halo(μ, 2; Z=0.005)
    x, _, z = r
    _, ẏ, _ = v
    u, T = AstrodynamicalSolvers.CR3BPSolvers.halo(x, z, ẏ, μ, T)

    @test u ≈ [1.180859455641048, 0.0, -0.006335144846688764, 0.0, -0.15608881601817765, 0.0]
    @test T ≈ 3.415202902714686

end

@testset "Dynamic Orbit Correction" begin

    μ = 0.012150584395829193
    u, T = halo(μ, 1)

    @test u ≈ [0.8567678285004178, 0.0, 0.0, 0.0, -0.14693135696819282, 0.0]
    @test T ≈ 2.7536820160579087

    u, T = halo(μ, 2; amplitude=0.005)
    @test u ≈ [1.180859455641048, 0.0, -0.006335144846688764, 0.0, -0.15608881601817765, 0.0]
    @test T ≈ 3.415202902714686

end

end