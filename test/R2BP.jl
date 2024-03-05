"""
Tests for R2BP dynamics.
"""
module R2BPTests

using AstrodynamicalModels, ModelingToolkit, Test

@testset "R2BP Model Constructors" begin
    @test R2BSystem(; stm=false) isa ODESystem
    @test R2BSystem(; stm=true) isa ODESystem

    @test rand(R2BState) isa R2BState
    @test rand(R2BParameters) isa R2BParameters
    @test system(rand(R2BParameters)) isa ODESystem
    @test dynamics(rand(R2BParameters)) isa ODEFunction
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
