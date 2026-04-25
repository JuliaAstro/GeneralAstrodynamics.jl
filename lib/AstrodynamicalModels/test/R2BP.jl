"""
Tests for R2BP dynamics.
"""
module R2BPTests

using Test
using AstrodynamicalModels:
    R2BFunction,
    R2BParameters,
    R2BState,
    R2BSystem,
    dynamics,
    system
using ModelingToolkit: System, ODEFunction, get_p, get_u0

@testset "R2BP Model Constructors" begin
    @test R2BSystem(; stm=false) isa System
    @test R2BSystem(; stm=true) isa System

    @test rand(R2BState) isa R2BState
    @test rand(R2BParameters) isa R2BParameters
    @test system(rand(R2BParameters)) isa System
    @test dynamics(rand(R2BParameters)) isa ODEFunction
end

@testset "R2BP Model Calculations" begin
    vectorfield = R2BFunction()
    sys = vectorfield.sys

    u0 = get_u0(sys, [
        [:x, :y, :z] .=> [-11e3, 0, 5e3]
        [:ẋ, :ẏ, :ż] .=> [0.0, 0.0, 0.0]
    ])

    p = get_p(sys, [:μ => 398600.4354360959])

    @test isapprox(
        vectorfield(u0, p, NaN),
        [0.0, 0.0, 0.0, 0.00248543, -0.0, -0.00112974];
        atol = 1e-8
    )
end

end # module
