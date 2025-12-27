"""
Tests for CR3BP dynamics.
"""
module CR3BPTests

using Test
using ModelingToolkit: ODEFunction, System, get_p, get_u0
using AstrodynamicalModels:
    CR3BSystem,
    CR3BState,
    CR3BParameters,
    CR3BFunction,
    system,
    dynamics

@testset "CR3BP Model Constructors" begin
    @test CR3BSystem(; stm=false) isa System
    @test CR3BSystem(; stm=true) isa System

    @test rand(CR3BState) isa CR3BState
    @test rand(CR3BParameters) isa CR3BParameters
    @test system(rand(CR3BParameters)) isa System
    @test dynamics(rand(CR3BParameters)) isa ODEFunction
end

@testset "CR3BP Model Calculations" begin
    vectorfield = CR3BFunction()
    sys = vectorfield.sys

    u0 = get_u0(sys, [
        [:x, :y, :z] .=> [0.9253021269565836, 0, 0]
        [:ẋ, :ẏ, :ż] .=> [0.0, 0.05852663414965813, 0.0]
    ])

    p = get_p(sys, [:μ => 0.0009536838895767625])

    @test isapprox(
        vectorfield(u0, p, NaN),
        [0.0, 0.05852663414965813, 0.0, 0.053265045303684255, 0.0, -0.0];
        atol = 1e-8
    )
end

end
