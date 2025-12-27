"""
Tests for R2BP dynamics.
"""
module AttitudeTests

using AstrodynamicalModels, ModelingToolkit, LinearAlgebra, Test

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

# See https://github.com/JuliaAstro/GeneralAstrodynamics.jl/issues/270
@testset "Attitude Model Calculations" begin
    vectorfield = AttitudeFunction()
    sys = vectorfield.sys

    op = [
        # u0
        :q => [0, 0, 0, 1],
        :ω => [0.1, 0.1, 0.1],
        # p
        :J => diagm([0.1, 0.2, 0.3]),
        :L => [0, 0, 0],
        :f => [0, 0, 0],
    ]

    prob =  ODEProblem(sys, op, NaN)

    @test isapprox(
        vectorfield(prob.u0, prob.p, NaN),
        [0.05, 0.05, 0.05, -0.0, -0.009999999999999995, 0.01, -0.003333333333333335];
        atol = 1e-8
    )

    ## TODO: Alternative approach. Wait to hear back from Joey
    """
    using ModelingToolkit: get_p

    let
       q = [0, 0, 0, 1]
       ω = [0.1, 0.1, 0.1]
       x = vcat(q, ω)

       p = get_p(sys, [
           :J => diagm([0.1, 0.2, 0.3]),
           :L => [0, 0, 0],
           :f => [0, 0, 0]
       ])

       vectorfield(x, p, NaN)
    end
    """
end

end
