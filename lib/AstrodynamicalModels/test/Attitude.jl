"""
Tests for R2BP dynamics.
"""
module AttitudeTests

using AstrodynamicalModels, ModelingToolkit, LinearAlgebra, Test
using ModelingToolkit: get_p, get_u0

@testset "Attitude Model Constructors" begin
    model = AttitudeSystem()
    @test model isa System

    model = AttitudeSystem(stm=true)
    @test model isa System

    @test rand(AttitudeState) isa AttitudeState
    @test rand(AttitudeParameters) isa AttitudeParameters
    @test system(rand(AttitudeParameters)) isa ModelingToolkit.System
    @test dynamics(rand(AttitudeParameters)) isa ModelingToolkit.ODEFunction
end

@testset "Attitude Model Calculations" begin
    vectorfield = AttitudeFunction()
    sys = vectorfield.sys

    u0 = get_u0(sys, [
        :q => [0, 0, 0, 1]
        :ω => [0.1, 0.1, 0.1]
    ])

    p = get_p(sys, [
       :J => diagm([0.1, 0.2, 0.3])
       :L => [0, 0, 0]
       :f => [0, 0, 0]
    ])

    # see #270
    @test isapprox(
        vectorfield(u0, p, NaN),
        [0.05, 0.05, 0.05, -0.0, -0.009999999999999995, 0.01, -0.003333333333333335];
        atol=1e-8
    )
end

end
