"""
Tests for R2BP dynamics.
"""
module R2BPTests

using AstrodynamicalModels, ModelingToolkit, Test

@testset "R2BP Model Constructors" begin
    model = R2BSystem(; stm=false)
    model = R2BSystem(; stm=true)
    @test true
end

@testset "R2BP Model Calculations" begin
    vectorfield = R2BFunction()

    r = [-11e3, 0, 5e3]
    v = [0.0, 0.0, 0.0]
    μ = 398600.4354360959

    @test isapprox(
        vectorfield(vcat(r, v), [μ], NaN),
        [0.0, 0.0, 0.0, 0.00248543, -0.0, -0.00112974];
        atol=1e-8
    )
end

end
