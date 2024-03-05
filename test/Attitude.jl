"""
Tests for R2BP dynamics.
"""
module AttitudeTests

using AstrodynamicalModels, ModelingToolkit, LinearAlgebra, Test

@testset "Attitude Model Constructors" begin
    model = AttitudeSystem()
    @test model isa ODESystem

    model = AttitudeSystem(stm=true)
    @test model isa ODESystem

    @test rand(AttitudeState) isa AttitudeState
    @test rand(AttitudeParameters) isa AttitudeParameters
    @test system(rand(AttitudeParameters)) isa ModelingToolkit.ODESystem
    @test dynamics(rand(AttitudeParameters)) isa ModelingToolkit.ODEFunction
end

@testset "Attitude Model Calculations" begin
    vectorfield = AttitudeFunction()

    q = [0, 0, 0, 1]
    ω = [0.1, 0.1, 0.1]
    x = vcat(q, ω)

    J = diagm([0.1, 0.2, 0.3])
    L = [0, 0, 0]
    f = [0, 0, 0]
    p = vcat(vec(J), L, f)

    @test isapprox(
        vectorfield(x, p, NaN),
        [0.05, 0.05, 0.05, -0.0, -0.009999999999999995, 0.01, -0.003333333333333335];
        atol=1e-8
    )
end

end
